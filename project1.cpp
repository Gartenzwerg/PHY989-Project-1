#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#define EPS    3.0e-14
#define MAXIT  10
#define ZERO   1.0E-10
#define PI     3.14159265358979324E+00    
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

double V (double k1, double k2, double V0, double mu) {

   double factor, num, denom;

   factor = V0/(4.0*k1*k2); 
   num = (k1 + k2)*(k1 + k2) + mu*mu;
   denom = (k1 - k2)*(k1 - k2) + mu*mu;

   return factor*log(num/denom);
}


/*

   Main Method

*/

int main () {

      // Define any necessary variables

      int N = 10;

      double k0,k02;
      int k0Points = 10;
      double k0Start = 1.0E-03;
      double k0Step = 1.0E-03;

      int i,j,k;
      double temp1,temp2;

      double x1,x2;     // Variables for the integration in position space.
      double x[N];
      double w[N];

      double kMesh[N];      // Variable for the integration in momentum space.
      double kWeights[N];

      double k1,k2;
      double const meshConst = 200.0; // Using MeV as the unit.
      double const V0 = 1.0;
      double const mN = 939.0;
      double const mu = mN/2.0;

      double u[N+1];
      double Vk[N+1][N+1];
      double A[N+1][N+1];
      double R[N+1][N+1];

      bool printArrays = false;

      // Initialize arrays

      for (i = 0; i<N; i++) {
         x[i] = 0.0;
         w[i] = 0.0;
         kMesh[i] = 0.0;
         kWeights[i] = 0.0;
      }
      for (i = 0; i<N+1; i++) {
         u[i] = 0.0;
         for (j = 0; j<N+1; j++) {
            Vk[i][j] = 0.0;
            A[i][j] = 0.0;
            R[i][j] = 0.012345;
         }
      }

      // Global loop over k0

      cout<<" k0 "<<"\t"<<" R(k0,k0) "<<"\n";

      for (k = 0; k<k0Points; k++) {

         k0 = k0Start + k0Step*k;
         k02 = k0*k0;
   
         // Call GaussLegendreQuadrature to get x[] and w[]

         GaussLegendreQuadrature(1.0,-1.0,x,w,N);

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
            kMesh[i] = meshConst*tan(temp1);
            kWeights[i] = meshConst*(PI/4.0)*(w[i]/(4.0*cos(temp1)*cos(temp1)));
         }

         if (printArrays) {
            cout<<"Mesh points & weights in Momentum Space"<<"\n";
            for (i = 0; i<N; i++) {
               cout<<i<<"\t"<<kMesh[i]<<"\t"<<kWeights[i]<<"\n";
            }
            cout<<"\n";
         }

         // Using these weights and mesh points, construct u[j] and Vk[i][j]

         for (i = 0; i<N; i++) {
            temp1 = kMesh[i]*kMesh[i];
            u[i] = (2.0/PI)*((kWeights[i]*temp1)/((k02-temp1)/mN));
         }
         temp2 = 0.0;
         for (i = 0; i<N; i++) {
            temp2 += (k02*kWeights[i])/((k02-temp1)/mN);
         }
         u[N] = (-2.0/PI)*temp2;

         if (printArrays) {
            cout<<"U Array"<<"\n";
            for (i = 0; i<N+1; i++) {
               cout<<i<<"\t"<<u[i]<<"\n";
            }
            cout<<"\n";
         }
      
         for (i = 0; i<N+1; i++) {
            if (i == N) {k1 = k0;} 
            else {k1 = kMesh[i];}

            for (j = 0; j<N+1; j++) {
               if (j == N) {k2 = k0;} 
               else {k2 = kMesh[j];}

               Vk[i][j] = V(k1,k2,V0,mu);

            }
         }
      
         if (printArrays) {
            cout<<"V Matrix"<<"\n";
            for (i = 0; i<N+1; i++) {
               for (j = 0; j<N+1; j++) {
                  cout<<i<<"\t"<<j<<"\t"<<Vk[i][j]<<"\n";
               }
            }
            cout<<"\n";
         }

         // Using u[], construct the A matrix

         for (i = 0; i<N+1; i++) {
            if (i == N) {k1 = k0;} 
            else {k1 = kMesh[i];}
            for (j = 0; j<N+1; j++) {
               if (i == N) {k1 = k0;} 
               else {k1 = kMesh[i];}
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

         // Compute the R-matrix by R = A^-1 * V









         // Print out the value of R(N+1,N+1) as R(k0,k0)

         cout<<" "<<k0<<"\t"<<R[N][N]<<"\n";


      }

return 0;
}
