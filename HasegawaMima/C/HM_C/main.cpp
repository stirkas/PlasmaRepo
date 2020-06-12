#include <iostream>
#include <blas/cblas.h>
#include <lapacke.h>
#include "HM_Definitions.hpp"
#include "HM_Modeler.hpp"

int main()
{
   HM::SimParams params {.nt = 4};
   HM::Modeler modeler(params);

   //Test BLAS works.
   double x[] = {1.0, 2.0, 3.0};
   double coeff = 4.323;
   int one = 1;
   int n = 3;
   cblas_dscal(n, coeff, x, one);
   for (int i = 0;  i < n; i++)
   {
      std::cout << " " << x[i];
   }
   std::cout << std::endl;

   //Test LAPACK works.
   const int N = 3, NRHS = 2, LDA = N, LDB = N;
   int ipiv[N] = {0};
   int info = 0;
   double a[LDA*N] = {
      6.80, -2.11,  5.66,
     -6.05, -3.30,  5.36,
     -0.45,  2.58, -2.70
   };
   double b[LDB*NRHS] = {
      4.02,  6.19, -8.22,
     -1.56,  4.00, -8.67
   };
   info = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, NRHS, a, LDA, ipiv, b, LDB);
   for (int i = 0; i < N; i++)
   {
      std::cout << " " << ipiv[i];
   }
   std::cout << std::endl;

   return 0;
}