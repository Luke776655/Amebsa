/*
C++ implementation of the Nelder-Mead optimization algorithm with simulated annealing (Amebsa).
In this case algorithm search for the energetic minimum of the arrows system.
When all of the arrows are parallel, the system is in the global minimum.
GNU GPL license v3

Core algorithm code in line with
Numerical Receipes
The Art of Scientific Computing
Third Edition
W.H.Press, S.A.Teukolsky,
W.T Vetterling, B.P.Flannery

Łukasz Radziński
lukasz.radzinski _at_ gmail _dot_ com
*/

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
            pmin.resize(ndim);
		};
		void minimize(double func(vector <T> &v));
        vector <double> temperature = {1, 1e-10, 0};
        double fmin = numeric_limits<double>::max();
        vector <T> pmin;
        int NMAX = 1000000; //maximum allowed number of function evaluations
		int delta = 1; //displacement
		double ftol = 1e-16; //fractional convergence tolerance
		double TINY = 1e-7; //tiny value preventing from dividing by 0
        bool arrows_table_printing = false;
        bool show_iter_output = true;
        int iter_period = 1000;
	private:
		int N; //dimention of side of pseudo-square table
        void print_array(vector <T> tab);
		void print_table(vector <T> tab);
		void print_arrow(T s);
		void print_arrow_table(vector <T> tab);
		void print_result();
		vector <T> get_psum();
		double amotsa(double fac, double func(vector <T> &v));
		int amebsa_alg(double func(vector <T> &v));
		vector<T> point;
		int ndim;
		vector <vector <T> > p;
		vector <double> y;
        int ilo;
        int ihi;
        int inhi;
	    double ynhi;
		double yhi;
        double ylo;
        double ytry;
        double rtol;
		int mpts;
		int iter;
		int nfunc;
        double tt;
		vector <T> psum;
};

template <class T> void Amebsa<T>::print_array(vector <T> tab)
{
	/*
	printing values of the array
	*/
	for(int i = 0; i<N; i++)
	{
		cout<<tab[i]<<"  ";	
	}
	printf("\n");
}

template <class T> void Amebsa<T>::print_table(vector <T> tab)
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
		cout<<tab[i]<<"  ";	
	}
	printf("\n");
}

template <class T> void Amebsa<T>::print_arrow(T s)
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

template <class T> void Amebsa<T>::print_arrow_table(vector <T> tab)
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

template <class T> void Amebsa<T>::print_result()
{
	/*
	printing result of optimisation
	*/
	T z = 0;
	vector <double> k = y;
	vector <vector <T> > q = p;

	z = k[0];
	k[0] = k[ilo];
	k[ilo] = z;

    double g;
	for(int i = 0; i<ndim; i++)
	{
		g = q[0][i];
		q[0][i] = q[ilo][i];
		q[ilo][i] = g;
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

template <class T> vector <T> Amebsa<T>::get_psum()
{
	/*
	counting partial sum
	*/
	vector <T> psum(ndim);
	for (int j=0; j<ndim; j++)
	{
		T sum=0;
		for(int i=0; i<mpts; i++)
		{
			sum += p[i][j];
		}
		psum[j]=sum;
	}
	return psum;
}

template <class T> double Amebsa<T>::amotsa(double fac, double func(vector <T> &v))
{
	/*
	extrapolation by a factor fac through the face of the simplex across from the high point
	replacing the high point if the new point is better
	*/
	vector <T> ptry(ndim);
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

template <class T> void Amebsa<T>::minimize(double func(vector <T> &v))
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
		y[i] = func(x);
	}
	//parameters for iterating
	nfunc = 0;
	psum = get_psum();
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

template <class T> int Amebsa<T>::amebsa_alg(double func(vector <T> &v))
{
    ilo=0;
	ylo=y[ilo]+tt*log(0.5*(double)rand()/RAND_MAX);
	ihi=1;
	yhi=y[ihi]+tt*log(0.5*(double)rand()/RAND_MAX);

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

	rtol=2.0*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo])+TINY);
	if (rtol < ftol)
	{
		print_result();
		printf("Optimisation succeeded\n");
		return 0;
	}

	if (nfunc >= NMAX)
	{
		print_result();
		printf("Maximal number of iterations exceeded\n");
		return 1;
	}
	nfunc += 2;

	ytry = amotsa(-1.0, func); //simplex reflection
	if(ytry <= ylo)
	{
		ytry = amotsa(2.0, func); //simplex extrapolation
	}
	else if (ytry >= ynhi)
	{
		double ysave = yhi;
		ytry = amotsa(0.5, func); //simplex contraction
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
			psum = get_psum();
		}
	} else nfunc = nfunc-1;
    if(show_iter_output == true)
    {
	    if(iter%iter_period == 0)
	    {
		    print_result();
	    }
    }
	iter++;
	return -1;
}

