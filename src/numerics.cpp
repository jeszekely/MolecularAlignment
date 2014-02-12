#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <new>
#include <fftw3.h>
#include <gsl/gsl_sf_coupling.h>
#include "array_structs.h"
#include "numerics.h"

double hermite(int n, double x)
{
	if (n == 0) 
	{
		return 1;
	}
	else if (n == 1) 
	{
		return 2*x;
	}
	else 
	{
		double poly = 2*x*hermite(n-1,x) - 2*(n-1)*hermite(n-2,x);
		return poly;
	}
}

double factorial(int n)
{
	double fact;
	if (n==0) 
	{
		fact = 1;
	}
	else 
	{
		fact = n*factorial(n-1);
	}
	return fact;
}

double ho_wvfxn(int n, double x, double xcen, double omega, double m, double hbar)
{ //nth harmonic oscillator solution evaluated at x, centered at xcen
	double coeff1 	= 1/sqrt(pow(2.0,n)*factorial(n));
	double coeff2 	= sqrt(sqrt(m*omega/M_PI/hbar));
	double coeff 	= coeff1 * coeff2;
	x 				-= xcen;
	double fxn 		= coeff * exp(-1*m*omega*x*x/2/hbar) * hermite(n,sqrt(m*omega/hbar)*x);
	return fxn;
}

int kron_delta(int x, int y)
{
	if (x == y) 
	{
		return 1;
	} 
	else 
	{
		return 0;
	}
}

/****************************************
Calculates < J K M | D^2_{QS}| J'K'M'>
"Field Matter Interaction Matrix Element"
*****************************************/

double FMIME (int J, int K, int M, int Q, int S, int j, int k, int m)
{
	double coeff1 = pow(-1.0,k+m); 
	double coeff2 = sqrt((2.0*J + 1.0)*(2.0*j+1.0)); 
	double J1 = gsl_sf_coupling_3j(2*J, 4, 2*j, 2*M, 2*Q, -2*m);
	double J2 = gsl_sf_coupling_3j(2*J, 4, 2*j, 2*K, 2*S, -2*k); 
	return coeff2*coeff1*J2*J1;  
}