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
#include"func.hpp"
#include"amebsa.hpp"

using namespace std;

int N = 1000;

int main()
{
    vector<int> v = create_random_table<int>(N);
    //vector<int> v = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	Amebsa<int> a(v);
    //a.show_iter_output = false;
    //cout<<a.fmin<<endl;
    a.temperature = {0};
	a.minimize(multi_parabole);
    //cout<<a.fmin<<endl;
	return 0;
}
