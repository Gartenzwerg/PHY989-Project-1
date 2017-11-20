#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include<algorithm>
#include "lib.h"
#define EPS    3.0e-14
#define MAXIT  10
#define ZERO   1.0E-10
#define PI     3.14159265358979324E+00
#define HBARC	197.0		// hbarc = 197 MeV fm	    
using namespace std;

/*

   Function implementing a Kronecker Delta

*/

double kron (int a, int b) {
   if (a == b) return 1.0;
   return 0;
}

/*

   Function for the setting up the integration domain, mesh points, and mesh weights.

*/

void GaussLegendreQuadrature(double x1, double x2, double x[], double w[], int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
           ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
         p2 =0.0;

           /*
           ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

         for(j = 1; j <= n; j++) {
            p3 = p2;
            p2 = p1;
            p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
         }

           /*
           ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */

         pp = n * (z * p1 - p2)/(z * z - 1.0);
         z1 = z;
         z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /*
          ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
}

/*

   Function for the potential model, defined in momentum space. Takes in two 
   momenta k and k' (k1 and k2), a normalization constant V0, and a mass term mu.

*/

double V (double k1, double k2) {

	double mu = 137.0;	// Pion Mass in MeV
   double a = 1.0*mu;
   double Va = -10.463;
   double b = 4.0*mu;
   double Vb = -1650.6;
   double c = 7.0*mu;
   double Vc = 6484.3;

   double Vsum = 0.0;

   double k1k2p,k1k2m;
   double factor, num, denom;

   k1k2p = (k1 + k2)*(k1 + k2);
   k1k2m = (k1 - k2)*(k1 - k2);   
   factor = 1.0/(4.0*mu*k1*k2);

   Vsum += Va*factor*log((k1k2p + a*a)/(k1k2m + a*a));
   Vsum += Vb*factor*log((k1k2p + b*b)/(k1k2m + b*b));
   Vsum += Vc*factor*log((k1k2p + c*c)/(k1k2m + c*c));

   return Vsum;
}

/*

	Function for computing the inverse of a matrix

*/

void inverse(double **a, int n) {        
  int          i,j, *indx;
  double       d, *col, **y;

  // allocate space in memory
  indx = new int[n];
  col  = new double[n];
  y    = (double **) matrix(n, n, sizeof(double)); 
   
  ludcmp(a, n, indx, &d);   // LU decompose  a[][] 

  // find inverse of a[][] by columns 

  for(j = 0; j < n; j++) {

    // initialize right-side of linear equations 

    for(i = 0; i < n; i++) col[i] = 0.0;
    col[j] = 1.0;

    lubksb(a, n, indx, col);

    // save result in y[][] 

    for(i = 0; i < n; i++) y[i][j] = col[i];

  }   //j-loop over columns 
   
  // return the inverse matrix in a[][] 

  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) a[i][j] = y[i][j];
  } 
  free_matrix((void **) y);     // release local memory 
  delete [] col;
  delete []indx;

}

/*

	Function for multiplying two matricies

*/

void matmul (double** A, double** B, double** C, int N) {

	int i,j,k;
	double sum;

	for (i = 0; i<N; i++) {
		for (j = 0; j<N; j++) {
			sum = 0.0;
			for(k = 0; k<N; k++) {
				sum += A[i][k]*B[k][j];
			}
			C[i][j] = sum;
		}
	}

	return;
}

/*

   Main Method

*/

int main () {

      // Define any necessary variables

      int N = 100;

      double k0,k02;
      int EPoints = 70;
		double Emax = 350.0;
		double Emin = 0.0;
		double EStep = (Emax-Emin)/EPoints;

      int i,j,k;
      double temp1,temp2,temp_sum;

      double x1,x2;     // Variables for the integration in position space.
      double x[N];
      double w[N];

      double kMesh[N+1];      // Variable for the integration in momentum space.
      double kWeights[N];
		double kMax,kMin;

      double k1,k2;
      double const meshConst = HBARC; // Using MeV as the unit.
      double const mN = 939;

      double u[N+1];
      double** Vk;
      double** A;
		double** AInv;
		double** eye;
      double** R;

		A = (double **) matrix(N+1,N+1,sizeof(double));
		AInv = (double **) matrix(N+1,N+1,sizeof(double));
		eye = (double **) matrix(N+1,N+1,sizeof(double));
		Vk = (double **) matrix(N+1,N+1,sizeof(double));
		R = (double **) matrix(N+1,N+1,sizeof(double));

		double phase;
		double energy;

      bool printArrays = false;

      // Initialize arrays

      for (i = 0; i<N; i++) {
         x[i] = 0.0;
         w[i] = 0.0;
         kWeights[i] = 0.0;
      }
      for (i = 0; i<N+1; i++) {
         kMesh[i] = 0.0;
         u[i] = 0.0;
         for (j = 0; j<N+1; j++) {
            Vk[i][j] = 0.0;
            A[i][j] = 0.0;
            R[i][j] = 0.0;
         }
      }

		ofstream fout("data.txt");
		ofstream fvout("potential.txt");

      // Call GaussLegendreQuadrature to get x[] and w[]

      GaussLegendreQuadrature(-1.0,1.0,x,w,N);

      if (printArrays) {
         cout<<"Mesh points & weights in Position Space"<<"\n";
         for (i = 0; i<N; i++) {
            cout<<i<<"\t"<<x[i]<<"\t"<<w[i]<<"\n";
         }
         cout<<"\n";
      }

      // Convert mesh points to momentum space

      for (i = 0; i<N; i++) {
         temp1 = (PI/4.0)*(1.0 + x[i]);
			temp2 = 4.0*cos(temp1)*cos(temp1);
         kMesh[i] = meshConst*tan(temp1);
         kWeights[i] = meshConst*(PI)*(w[i]/temp2);
      }

      if (printArrays) {
         cout<<"Mesh points & weights in Momentum Space"<<"\n";
         for (i = 0; i<N; i++) {
            cout<<i<<"\t"<<kMesh[i]<<"\t"<<kWeights[i]<<"\n";
         }
         cout<<"\n";
      }

		kMin = *std::min_element(kMesh,kMesh+N);
		kMax = *std::max_element(kMesh,kMesh+N);

		cout<<"Smallest K Mesh Point: "<<kMin<<"\n";
		cout<<"Largest K Mesh Point: "<<kMax<<"\n";

      // Global loop over k0

      cout<<" E (MeV) "<<"\t"<<" Delta (Deg) "<<"\n";

      for (k = 0; k<EPoints; k++) {

         energy = (k+1)*EStep;
         k02 = energy*mN;
			k0 = sqrt(k02);
			kMesh[N] = k0;

         // Using these weights and mesh points, construct u[j] and Vk[i][j]

         for (i = 0; i<N; i++) {
            temp1 = kMesh[i]*kMesh[i];
				temp2 = (k02-temp1)/mN;
            u[i] = (2.0/PI)*((kWeights[i]*temp1)/temp2);
         }
         temp_sum = 0.0;
         for (i = 0; i<N; i++) {
				temp1 = kMesh[i]*kMesh[i];
				temp2 = (k02-temp1)/mN;
            temp_sum += (k02*kWeights[i])/temp2;
         }
         u[N] = (-2.0/PI)*temp_sum;

         if (printArrays) {
            cout<<"U Array"<<"\n";
            for (i = 0; i<N+1; i++) {
               cout<<i<<"\t"<<u[i]<<"\n";
            }
            cout<<"\n";
         }
      
         for (i = 0; i<N+1; i++) {
				k1 = kMesh[i];
            for (j = 0; j<N+1; j++) {
					k2 = kMesh[j];
               Vk[i][j] = V(k1,k2);
            }
         }

         if (printArrays) {
         	for (i = 0; i<N+1; i++) {
           		for (j = 0; j<N+1; j++) {
               	fvout<<i<<"\t"<<j<<"\t"<<Vk[i][j]<<"\n";
            	}
         	}
			}

         // Using u[], construct the A matrix

         for (i = 0; i<N+1; i++) {
            for (j = 0; j<N+1; j++) {
               A[i][j] = kron(i,j) - Vk[i][j]*u[j];
            }
         }

         if (printArrays) {
            cout<<"A Matrix"<<"\n";
            for (i = 0; i<N+1; i++) {
               for (j = 0; j<N+1; j++) {
                  cout<<i<<"\t"<<j<<"\t"<<A[i][j]<<"\n";
               }
            }
            cout<<"\n";
         }

			AInv = A;

         // Compute the R-matrix by R = A^-1 * V

			inverse(AInv,N+1);

         if (printArrays) {
            cout<<"A Inverse"<<"\n";
            for (i = 0; i<N+1; i++) {
               for (j = 0; j<N+1; j++) {
                  cout<<i<<"\t"<<j<<"\t"<<AInv[i][j]<<"\n";
               }
            }
            cout<<"\n";
         }

			matmul(AInv,A,eye,N+1);

         if (printArrays) {
            cout<<"A^-1 A = I"<<"\n";
            for (i = 0; i<N+1; i++) {
               for (j = 0; j<N+1; j++) {
                  cout<<i<<"\t"<<j<<"\t"<<eye[i][j]<<"\n";
               }
            }
            cout<<"\n";
         }

			matmul(AInv,Vk,R,N+1);

			//  Compute the phase shift for this particular k0 (or Energy)

			phase = (180/PI)*atan(-mN*k0*R[N][N]);

         // Print out the values

         cout<<" "<<energy<<"\t"<<phase<<"\n";
			fout<<" "<<energy<<"\t"<<phase<<"\n";


      }

return 0;
}
