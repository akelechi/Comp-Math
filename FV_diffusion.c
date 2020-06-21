/*Programmer: Kelechi Akwataghibe
  Date: 10th April 2020
  Purpose: To simulate an explicit finite volume scheme 
           for diffusion
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double  dx,tend, dtout, factor, a,b,D, ERR, time;
double  *x, *U, *F, *uEXACT;
int i,MM,M;
FILE *file;

void MESH();
void INIT(double);
void FLUX();
void OUTPUT(double, int);
void COMPARE(double);
void PDE(double, double);

int main(){
double dtEXPL, dt, tout;
int Nend, nsteps;
//Create or open a file for reading the input values

file = fopen("input.txt","r");
while (!feof (file)) {   
//Read from the file
fscanf(file,"%lf,%lf,%lf,%d,%lf,%lf,%lf",&a,&b,&D,&MM,&tend,&dtout,&factor);
}
//Close the file
fclose(file);





dx = 1/(double)MM;
M = (b-a)*MM;


//Construct arrays to hold relevant values

x = (double *)malloc(sizeof(double)*(M + 2));
U = (double *)malloc(sizeof(double)*(M + 2));
F = (double *)malloc(sizeof(double)*(M + 2));


time = 0;
MESH();
INIT(time);

dtEXPL = dx*dx/(2*D);   //max timestep for stability
dt = factor * dtEXPL;

Nend = (int)(tend/dt)+1; //No of timesteps

  
nsteps = 0;
tout = fmax(dtout, dt);

//The first output run
OUTPUT(0, M+1);
printf("\n\n");


/*Timestepping*/ 

for (nsteps = 1; nsteps<= Nend; nsteps = nsteps+1){
time = nsteps * dt;
FLUX();
PDE(dt, time);

if(time >= tout){
    COMPARE(time);
    OUTPUT(time, nsteps);
    tout = tout + dtout;
}

}

if(time >= tend){
    printf("DONE at time = %lf , nsteps = %d", time, nsteps );
    printf("MAX ERROR = %lf", ERR);
}
else{
    printf("Out of timesteps: need a bigger Nend");
}

    return 0;
}

void MESH(){
x[0] = a;
x[1] = x[0] + dx/2;

x[M+1] = b;

	for (i = 2; i <= M; i = i + 1) {
		 x[i]=x[1]+(i-1)*dx ;
	} 

}

void INIT(double time) {
	for (i = 0; i <= M+1; i = i + 1) {
		 U[i] = 0.0;
	} 

}

void OUTPUT(double t, int n){
    int ntime;
    printf("Profile at time: %lf \t nsteps = %d\n", t,n);
	ntime = (int)(t);

//To print the profile at the first run

    if (ntime == 0) {
       for ( i = 0; i<=M+1;  i = i + 1 ){
       printf("%lf \t %lf\n", x[i], U[i]); 
    }

    }
    else {
         //To print first step , or desired timesteps
        for ( i = 0; i<=M;  i = i + ntime ){
            printf("%lf \t %lf\n", x[i], U[i]); 
            }

        printf("%lf \t %lf\n", x[M + 1], U[M + 1]);
    }
    



   

    
}

void FLUX(){
    for(i=1; i<=M+1; i = i+1){
        F[i] = 0 - D*(U[i]-U[i-1])/(x[i]-x[i-1]);
    }
}

void PDE(double dt, double time ){
    
	U[0] = 1.0;
	U[M + 1] = erfc(0.5*b / sqrt(D * time));
    for(i=1; i<=M; i=i+1){
		U[i] = U[i] + (dt / dx)*(F[i] - F[i + 1]);
    }
}

void COMPARE(double time){
    //Construct array for exact solution
   uEXACT = (double *)malloc(sizeof(double)*(M + 2));
   double arg, ERRi = 0.0;

    for(i=0; i<=M+1; i=i+1){
        
        arg = 0.5 * x[i]/ sqrt(D*time);
        uEXACT[i] = erfc(arg);
		
        ERRi = fabs(U[i]-uEXACT[i]);
        ERR = fmax(ERRi, ERR);
    }

}