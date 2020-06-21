 
/*
Programmer: Kelechi Akwataghibe
Date: Feb 22, 2020
Purpose: To implement the euler forward-stepping algorithm to solve a simple linear ODE	
*/


#include<stdlib.h>
#include<stdio.h>
#include<math.h>

//Double is used instead of float to allow for more precision
//Declared global variable to allow for easy access from different remote subroutines.


double FCN(double, double);
double formula(double);




int main() {
	double xStar =  0.0949002896834991;
	double  tEnd;
	double t0;
	int Nmax;

	printf("What is the maximum number of iterations to be performed: ");
	scanf("%d", &Nmax);
	printf("Enter an initial time: ");
	scanf("%lf", &t0);
	printf("Enter the final time : ");
	scanf("%lf", &tEnd);

	EulerScheme(t0, tEnd, Nmax,xStar);

	
	return 1;
}

EulerScheme(double t0, double tEnd, int Nmax, double xZero) {
	FILE *f; //Creating File to Store Output
	f = fopen("values1.txt", "w"); //Creating/Opening the file to read values

	int n;
	double xn;
	double zed = 0;
	xn= xZero;
	double tn = t0;
	double dt = (tEnd - t0) /(double) Nmax;
	double errN, errMax, yExact;
	errMax = 0.0;

	for (n = 0; n < Nmax+1; n++) {
		if (xn > 0) {
			xn = xn + dt * FCN(tn, xn);
			tn = t0 + n*dt;

			fprintf(f, "%0.15e\t%0.15e\n", tn,xn); //Printing pair to the file
		}
	}
	
	fclose(f); //Close the file that was opened
	printf("\nYou output has been successful written to the file values.txt");//Inform user that output was written

	

}

double FCN(double t, double xn) {

	double cStar = 7.52e-7;
	double xStar = 0.0949002896834991;
	double gamma = 4.e-3;
	double mu = 1.e-3;
	double Fn = 0.0, DFn=0.0;
	double k = 5.e7;



	DFn = k*(1.05*cStar + mu*xStar*xStar*xStar - (mu)* (xn*xn*xn) - (cStar)*exp((gamma) / xn));
																						

	return DFn;
}

double formula(double tn) {

	return  exp(tn);
}