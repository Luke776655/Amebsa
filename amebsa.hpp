#include<cstdio>
#include<iostream>
#include<cstdlib>
#include<vector>
#include<cmath>
#include<limits>

using namespace std;

template <class T> class Amebsa
{
	public:
		Amebsa(vector <T> &pt)
		{
			srand (time(NULL));
			point = pt;
            N = pt.size();
			ndim = N;
		};
		void minimize(double func(vector <double> &v));
        vector <double> temperature = {1, 0.00000000001, 0};
        double fmin = numeric_limits<double>::max();
        int NMAX = 1000000; //maximum allowed number of function evaluations
		int delta = 1; //displacement
		double ftol = 1e-6; //fractional convergence tolerance
		double TINY = 1e-7; //tiny value preventing from dividing by 0
        bool arrows_table_printing = false;
        bool show_iter_output = true;
        int iter_period = 1000;
	private:
		int N; //dimention of side of pseudo-square table
        void print_array(vector <double> tab);
		void print_table(vector <double> tab);
		void print_arrow(double s);
		void print_arrow_table(vector <double> tab);
		void print_result(vector <double> y, vector <vector <double> > &p, int ndim, int ilo);
		vector <double> get_psum(vector <vector <double> > &p, int ndim, int mpts);
		double amotsa(vector <vector <double> > &p, vector <double> &psum, vector <double> &y, int ihi, int ndim, double fac, double func(vector <double> &v));
		int amebsa_alg(double func(vector <double> &v));
		vector<double> point;
		int ndim;
		vector <vector <double> > p;
		vector <double> y;
		double yhi;
		int mpts;
		int iter;
		int nfunc;
        double tt;
		vector <double> psum;
};

template <class T> void Amebsa<T>::print_array(vector <double> tab)
{
	/*
	printing values of the array
	*/
	for(int i = 0; i<N; i++)
	{
		printf("%.2f  ", tab[i]);	
	}
	printf("\n");
}

template <class T> void Amebsa<T>::print_table(vector <double> tab)
{
	/*
	printing values of the table
	*/
	for(int i = 0; i<N; i++)
	{
		if(i%(int)sqrt(N)==0)
		{
			printf("\n");
		}
		printf("%+.2f\t", tab[i]);	
	}
	printf("\n");
}

template <class T> void Amebsa<T>::print_arrow(double s)
{
	/*
	printing value as an arrow with proper slope
	*/
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

template <class T> void Amebsa<T>::print_arrow_table(vector <double> tab)
{
	/*
	printing values of the table as arrows with proper slope
	*/
	for(int i=0; i<N; i++)
	{
		if(i%(int)sqrt(N)==0)
			printf("\n");
		print_arrow(fmod(fmod(tab[i], 2*M_PI)+2*M_PI, 2*M_PI));
	}
	printf("\n\n");
}

template <class T> void Amebsa<T>::print_result(vector <double> y, vector <vector <double> > &p, int ndim, int ilo)
{
	/*
	printing result of optimisation
	*/
	double z = 0;
	vector <double> k = y;
	vector <vector <double> > q = p;

	z = k[0];
	k[0] = k[ilo];
	k[ilo] = z;
	vector <double> pmin(ndim);
	for(int i = 0; i<ndim; i++)
	{
		z = q[0][i];
		q[0][i] = q[ilo][i];
		q[ilo][i] = z;
		pmin[i] = q[0][i];
	}
	fmin=k[0];
    printf("\nIteration %d\tTemperature %g\n", iter, -tt);
    if(arrows_table_printing == true)
    {
	    print_table(pmin);
	    print_arrow_table(pmin);
    }
    else
        print_array(pmin);
	printf("Energy of the system %g\n\n", fmin);
}

template <class T> vector <double> Amebsa<T>::get_psum(vector <vector <double> > &p, int ndim, int mpts)
{
	/*
	counting partial sum
	*/
	vector <double> psum(ndim);
	for (int j=0; j<ndim; j++)
	{
		double sum=0.0;
		for(int i=0; i<mpts; i++)
		{
			sum += p[i][j];
		}
		psum[j]=sum;
	}
	return psum;
}

template <class T> double Amebsa<T>::amotsa(vector <vector <double> > &p, vector <double> &psum, vector <double> &y, int ihi, int ndim, double fac, double func(vector <double> &v))
{
	/*
	extrapolation by a factor fac through the face of the simplex across from the high point
	replacing the high point if the new point is better
	*/
	vector <double> ptry(ndim);
	double fac1=(1.0-fac)/ndim;
	double fac2=fac1-fac;
	for(int j=0; j<ndim; j++)
	{
		ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	}
	double ytry=func(ptry);
	double yflu=ytry-tt*log(0.5*(double)rand()/RAND_MAX);
	if(yflu < yhi)
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

template <class T> void Amebsa<T>::minimize(double func(vector <double> &v))
{

	//creating delta values table
    vector <int> delta_tab;
	for(int i=0; i<ndim; i++)
	{
		delta_tab.push_back(delta);
	}

	//adding delta values to the point table with extended dimention as p simplex

	for(int i=0; i<ndim+1; i++)
	{
		vector<T> k;
		for(int j=0; j<ndim; j++)
		{
			k.push_back(point[j]);
		}
		p.push_back(k);
		if(i!=0)
		{
			p[i][i-1]+=delta_tab[i-1];
		}
	}
	//getting y table of solutions
	y.resize(ndim+1);
	yhi = 0;
	mpts = p.size(); //number of rows
	for(int i = 0; i<mpts; i++)
	{
		vector <T> x;
		for(int j=0; j<ndim; j++)
		{
			x.push_back(p[i][j]);
		}
		y[i] = (func(x));
	}
	//parameters for iterating
	nfunc = 0;
	psum = get_psum(p, ndim, mpts);
    for(int i=0; i<temperature.size(); i++)
    {
	    int w = -1;
        iter = 0;
        nfunc = 0;
        tt = -temperature[i];
	    while(w<0)
	    {
	    	w = amebsa_alg(func);
	    }
    }
}

template <class T> int Amebsa<T>::amebsa_alg(double func(vector <double> &v))
{
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

	double rtol=2.0*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo])+TINY);
	if (rtol < ftol)
	{
		print_result(y, p, ndim, ilo);
		printf("Optimisation succeeded\n");
		return 0;
	}

	if (nfunc >= NMAX)
	{
		print_result(y, p, ndim, ilo);
		printf("Maximal number of iterations exceeded\n");
		return 1;
	}
	nfunc += 2;

	double ytry = amotsa(p, psum, y, ihi, ndim, -1.0, func); //simplex reflection
	if(ytry <= ylo)
	{
		ytry = amotsa(p, psum, y, ihi, ndim, 2.0, func); //simplex extrapolation
	}
	else if (ytry >= ynhi)
	{
		double ysave = yhi;
		ytry = amotsa(p, psum, y, ihi, ndim, 0.5, func); //simplex contraction
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
    if(show_iter_output == true)
    {
	    if(iter%iter_period == 0)
	    {
		    print_result(y, p, ndim, ilo);
	    }
    }
	iter++;
	return -1;
}

