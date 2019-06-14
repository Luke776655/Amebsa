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

int N = 100;

int main()
{
    //vector<double> v = create_random_table<double>(N);
    vector<double> v = create_random_2PI_table(N);
    //vector<double> v = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	Amebsa<double> a(v);
    //a.show_iter_output = false;
    //a.temperature = {0};
    a.arrows_table_printing = true;
	a.minimize(ferromagnetic);
	return 0;
}
