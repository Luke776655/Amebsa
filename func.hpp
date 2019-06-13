#include<cstdio>
#include<iostream>
#include<cstdlib>
#include<vector>
#include<cmath>

using namespace std;

template <typename T> double multi_parabole(vector <T> &A)
{
	//simple multidimentional parabola function for testing optimisation
	//global minimum in 0
	T s = 0;
	for(int i=0; i<A.size(); i++)
	{
		s += A[i]*A[i];
	}
	return s;
}


double ferromagnetic(vector <double> &A)
{
	/*
	potential energy function of the arrows system for testing optimisation
	when the neighbour arrow is parallel: E = 0
	when the neighbour arrow is orthogonal: E = 1
	when the neighbour arrow is antiparallel: E = 2
	global E is the sum of all arrows energy
	global energetic minimum is when all of the arrows are parallel, then global E = 0
	*/
	double s = 0;
	int N = A.size();
	for (int i = 0; i<N; i++)
	{
		int k = sqrt(N);
		if(i>=k) //top boundary
			s = s-cos(abs(A[i]-A[i-k]))+1;
		if(i<N-k) //bottom boundary
			s = s-cos(abs(A[i]-A[i+k]))+1;
		if(i%k!=0) //left boundary
			s = s-cos(abs(A[i]-A[i-1]))+1;
		if(i%k!=(k-1)) //right boundary
			s = s-cos(abs(A[i]-A[i+1]))+1;
	}
	return s;
}

vector <double> create_random_2PI_table(int n)
{
	/*
	creating a starting point of the system:
	n-dimentional table with random values
	values range: [0, 2*pi)
	*/
	vector <double> tab(n);
	for(int i = 0; i<n; i++)
	{
		tab[i] = fmod(rand(), 2*M_PI);
	}
	return tab;
}

template <typename T> vector <T> create_random_table(int n)
{
	/*
	creating a starting point of the system:
	n-dimentional table with random values
	*/
	vector <T> tab(n);
	for(int i = 0; i<n; i++)
	{
		tab[i] = rand()%1024;
	}
	return tab;
}
