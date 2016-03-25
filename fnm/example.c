#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gauss_legendre.h"

#ifndef PI
	#define PI 3.1415926535897932384626433832795028841971693993751
#endif
#ifndef FABS
	#define FABS(a) ((a)>=0?(a):-(a))
#endif


double f(double x, void* data)
{
	return sin(x);
}

int main(int argc, char* argv[])
{

	/* numerical approximation of integral */
	double approx;		

	/* true value of int(sin(x), x=0..Pi) = 2.0*/
	double exact = 2.0; 

	/* approximation error */
	double error;       

	int i;

	printf("Numerical Approximation of int(sin(x), x=0..Pi) by Gauss-Legendre Quadrature:\n");
	for (i=2;i<=128;i++)
	{
		approx = gauss_legendre(i,f,NULL,0,PI);
		error = approx-exact;
		printf("n = %4d: error = %.15g\n",i,FABS(error));
	}

}

