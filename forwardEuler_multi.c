
/*
Programmer: Kelechi Akwataghibe
Date: February 24, 2020
Purpose: To implement the euler forward-stepping algorithm to solve a simple linear ODE	System
*/


#include<stdlib.h>
#include<stdio.h>
#include<math.h>




struct Func FCN(double, double, double);
double formula(double);

//Double is used instead of float to allow for more precision
//Declared global variable to allow for easy access from different remote subroutines.

double x1n,x2n;
double arr[2];
FILE *file;
//Since C does not return multiple values, we create a struct in order to hold both F1 and F2 as they are returned from the FCN routine.
struct Func {
	double F1, F2;
};
int main() {
	double tEnd;
	double x10,x20;
	int t0, Nmax;
	file = fopen("output.txt", "w");
	printf("Enter the initial value of the first ODE(x10): ");
	scanf("%lf", &x10);
	printf("Enter the initial value of the second ODE(x20): ");
	scanf("%lf", &x20);
	printf("What is the maximum number of iterations to be performed: ");
	scanf("%d", &Nmax);
	printf("Enter an initial time: ");
	scanf("%d", &t0);
	printf("Enter the final time : ");
	scanf("%lf", &tEnd);

	EulerScheme1(t0, tEnd, Nmax,x10,x20);
	fclose(file);
	
	return 1;
}

EulerScheme1(int t0, double tEnd, int Nmax, double x10, double x20) {
	int n;

	x1n= x10;
	x2n = x20;
	double tn = t0;
	double dt = (tEnd - t0) / Nmax;
	double errN, errMax, y1Exact;
	errMax = 0.0;

	for (n = 0; n < Nmax+1; n++) {
		struct Func fcns = FCN(tn, x1n, x2n);
        if(x1n > 0 & x2n >0){
		x1n = x1n + dt * fcns.F1;
		x2n = x2n + dt * fcns.F2;
		tn = t0 + n*dt;
		
		fprintf(file,"%0.15e\t%0.15e\t%0.15e\n",tn, x1n, x2n);
        }
        else{
            fprintf(file,"%0.15e\t%0.15e\t%0.15e\n",tn, x1n, x2n);
            break;
        }
	}
	
	

}


struct Func FCN(double t, double x1, double x2) {

	double cStar = 7.52e-7;
	double x1Star = 0.05;
	double x2Star = 0.09;
	double gamma = 4.e-3;
	double mu = 1.e-3;
	double k = 5.e7;
	double 	DF1n = 0.0, DF2n = 0.0;
	double c = 0.0;

	c = 1.05*cStar + mu * (x1Star*x1Star *x1Star + x2Star*x2Star*x2Star) - mu * (x1*x1*x1 + x2*x2*x2);

	DF1n = k*(c-cStar * exp((gamma)/x1n)); 
	DF2n = k*(c - cStar * exp((gamma) / x2n));

	struct Func fcts = { DF1n,DF2n };
	return fcts;
	
}

