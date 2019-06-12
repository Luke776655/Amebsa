#include<cstdio>
#include<iostream>
#include<cstdlib>
#include<vector>
#include<cmath>

int N = 30;
double yhi = 0;
std::vector <double> y(N*N+1);
double temperature;
double tt;
/*
double func(std::vector <double> &x)
{
	double s = 0.0;
	for(int i=0; i<x.size(); i++)
	{
		s += x[i]*x[i];
	}
	return s;
}
*/
double func(std::vector <double> &A)
{
	double s = 0;
	for (int i = 0; i< N*N; i++)
	{
		if(i>=N) //ograniczenie górne
			s = s-cos(std::abs(A[i]-A[i-N]))+1;
		if(i<N*N-N) //ograniczenie dolne
			s = s-cos(std::abs(A[i]-A[i+N]))+1;
		if(i%N!=0) //ograniczenie lewe
			s = s-cos(std::abs(A[i]-A[i-1]))+1;
		if(i%N!=(N-1)) //ograniczenie prawe
			s = s-cos(std::abs(A[i]-A[i+1]))+1;
	}
	return s;
}

void print_arrow(double s)
{
	if(s<(M_PI/8+0*M_PI/4) || s>=(M_PI/8+7*M_PI/4))
		printf("↑ ");
	else if(s>=(M_PI/8+0*M_PI/4) && s<(M_PI/8+1*M_PI/4))
		printf("↗ ");
	else if(s>=(M_PI/8+1*M_PI/4) && s<(M_PI/8+2*M_PI/4))
		printf("→ ");
	else if(s>=(M_PI/8+2*M_PI/4) && s<(M_PI/8+3*M_PI/4))
		printf("↘ ");
	else if(s>=(M_PI/8+3*M_PI/4) && s<(M_PI/8+4*M_PI/4))
		printf("↓ ");
	else if(s>=(M_PI/8+4*M_PI/4) && s<(M_PI/8+5*M_PI/4))
		printf("↙ ");
	else if(s>=(M_PI/8+5*M_PI/4) && s<(M_PI/8+6*M_PI/4))
		printf("← ");
	else if(s>=(M_PI/8+6*M_PI/4) && s<(M_PI/8+7*M_PI/4))
		printf("↖ ");
	else
		printf("x ");
}

void print_arrow_matrix(std::vector <double> tab, int N)
{
	for(int i=0; i<N*N; i++)
	{
		if(i%N==0)
			printf("\n");
		print_arrow(std::fmod(std::fmod(tab[i], 2*M_PI)+2*M_PI, 2*M_PI));
	}
	printf("\n\n");
}

std::vector <double> create_matrix(int N)
{
	std::vector <double> tab(N*N);
	for(int i = 0; i<N*N; i++)
	{
		tab[i] = std::fmod(i, 2*M_PI);
	}
	return tab;
}

void print_matrix(std::vector <double> tab, int N)
{
	for(int i = 0; i<N*N; i++)
	{
		if(i%N==0)
		{
			printf("\n");
		}
			//if(tab[i]>=0)
			//	printf(" ");
			printf("%+.2f\t", tab[i]);
			
	}
	printf("\n");
}

std::vector <double> print_result(std::vector <std::vector <double> > &p, int ndim, int ilo, int N)
{
	double z = 0;
	std::vector <double> k = y;
	std::vector <std::vector <double> > q = p;

	z = k[0];
	k[0] = k[ilo];
	k[ilo] = z;
	std::vector <double> pmin(ndim);
	for(int i = 0; i<ndim; i++)
	{
		z = q[0][i];
		q[0][i] = q[ilo][i];
		q[ilo][i] = z;
		pmin[i] = q[0][i];
	}
	double fmin=k[0];
	print_matrix(pmin, N);
	print_arrow_matrix(pmin, N);
	printf("Energia układu, %f\n\n", fmin);
	return pmin;
}

std::vector <double> get_psum(std::vector <std::vector <double> > &p, int ndim, int mpts)
{
	std::vector <double> psum(ndim);
	for (int n=0; n<ndim; n++)
	{
		double sum=0.0;
		for(int m=0; m<mpts; m++)
		{
			sum += p[m][n];
		}
		psum[n]=sum;
	}
	return psum;
}

double amotsa(std::vector <std::vector <double> > &p, std::vector <double> &psum, const int ihi, int ndim, const double fac)
{
	std::vector <double> ptry(ndim);
	double fac1=(1.0-fac)/ndim;
	double fac2=fac1-fac;
	for(int j=0; j<ndim; j++)
	{
		ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	}
	double ytry=func(ptry);
	double yflu=ytry-tt*log(0.5*(double)rand()/RAND_MAX);
	if (yflu < yhi)
	{
		y[ihi]=ytry;
		yhi=yflu;
		for(int j=0; j<ndim; j++)
		{
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	return yflu;
}

int main()
{
std::vector<double> points(N*N);
std::vector<double> final_matrix;
for(int i = 0; i<4; i++)
{
if(i==0)
	points = create_matrix(N);
else
	points = final_matrix;

if(i==0)
	temperature = 1;
else if(i==1)
	temperature = 1e-4;
else if(i==2)
	temperature = 1e-12;
else if(i==3)
	temperature = 0;


tt = -temperature;		
print_matrix(points, N);
int ndim = points.size();
int delta = 1;
std::vector <int> delta_tab;
std::vector <std::vector <double> > p;
std::vector <std::vector <double> > pp;
int nfunc = 0;
double ftol = 1e-6;
srand (time(NULL));

for(int i=0; i<ndim; i++)
{
	delta_tab.push_back(delta);
}

for(int i=0; i<ndim+1; i++)
{
	std::vector<double> k;
	for(int j=0; j<ndim; j++)
	{
		k.push_back(points[j]);
	}
	pp.push_back(k);
	if(i!=0)
	{
		pp[i][i-1]+=delta_tab[i-1];
	}
}

int NMAX = 10000000;
double TINY = 1e-7;
int mpts = pp.size();
ndim = pp[0].size();
std::vector <double> psum;
p = pp;

psum = get_psum(p, ndim, mpts);


for(int i = 0; i<mpts; i++)
{
	std::vector <double> x;
	for(int j=0; j<ndim; j++)
	{
		x.push_back(p[i][j]);
	}
	y[i] = (func(x));
}

int it = 0;



while(true)
{
	//std::cout << "Iteracja" << it << std::endl;
	//for(int i=0; i<y.size(); i++)
	//	std::cout << y[i] << ' ';
	//std::cout<<std::endl;
	int ilo=0;
	double ylo=y[ilo]+tt*log(0.5*(double)rand()/RAND_MAX);
	int ihi=1;
	yhi=y[ihi]+tt*log(0.5*(double)rand()/RAND_MAX);
	int inhi;
	double ynhi;
	if(ylo>yhi)
	{
		inhi = 1;
		ynhi = yhi;
		ihi = 0;
		yhi = ylo;
	}
	else
	{
		inhi = 0;
		ynhi = ylo;
	}

	for (int i=0; i<mpts; i++)
	{
		double yi = y[i]+tt*log(0.5*(double)rand()/RAND_MAX);
		//std::cout<< "yi: " << yi << " ylo: "<< ylo << std::endl;
		if(yi<=ylo)
		{
			ilo = i;
			ylo = yi;
		}
		if(yi>yhi)
		{
			inhi = ihi;
			ynhi = yhi;
			ihi = i;
			yhi = yi;
		}
		else if(yi > ynhi && i != ihi)
		{
			inhi = i;
			ynhi = yi;
		}
	}
	double rtol=2.0*std::abs(yhi-ylo)/(std::abs(yhi)+std::abs(ylo)+TINY);
	//std::cout<<"abs "<<std::abs(yhi-ylo)<<std::endl;
	//std::cout<<"RTOL "<<rtol<<std::endl;

	if (rtol < ftol)
	{
		final_matrix = print_result(p, ndim, ilo, N);
		std::cout << "Iteracja " << it << std::endl;
		std::cout << "Temperatura " << temperature << std::endl;
		std::cout << "Sukces!\n";
		std::cout << std::endl;
		break;
	}

	if (nfunc >= NMAX)
	{
		final_matrix = print_result(p, ndim, ilo, N);
		std::cout << "Iteracja " << it << std::endl;
		std::cout << "Temperatura " << temperature << std::endl;
		std::cout << "MAX\n";
		std::cout << std::endl;
		break;
	}
	nfunc += 2;

	double ytry = amotsa(p, psum, ihi, ndim, -1.0);
	if(ytry <= ylo)
	{
		ytry=amotsa(p, psum, ihi, ndim, 2.0);
	}
	else if (ytry >= ynhi)
	{
		double ysave=yhi;
		ytry=amotsa(p, psum, ihi, ndim, 0.5);
		if(ytry >= ysave)
		{
			for (int i=0;i<mpts;i++) {
				if (i != ilo) {
					for (int j=0;j<ndim;j++)
					{
						psum[j]=0.5*(p[i][j]+p[ilo][j]);
						p[i][j]=psum[j];
					}
					y[i]=func(psum);
				}
			}
			nfunc += ndim;
			psum = get_psum(p, ndim, mpts);
		}
	} else nfunc = nfunc-1;

	if(it%1000 == 0)
	{
		print_result(p, ndim, ilo, N);
		std::cout << "Iteracja " << it << std::endl;
		std::cout << "Temperatura " << temperature << std::endl;
		//printf("%f\n", log(0.5*(double)rand()/RAND_MAX));
	}
	it++;
	
}
}
}

