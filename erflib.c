// erflib.c
// File imported from Perseo

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "erflib.h"



/**
 *  erflib.c
 *
 *     Fornisce una libraria di funzioni per determinare
 *  la densita' di probabilita' di una variabile casuale
 *  normale, la sua funzione cumulativa (erf) e la cumulativa
 *  complementare (erfc).
 *     Note queste fornisce dei metodi per valutare funzioni
 *  cumulative di variabili casuali con probabilita' relative
 *  uguali a quelle delle variabili normali ma con dominio
 *  chiuso e limitato.
 *     Sono messe a disposizione inoltre funzioni che determi-
 *  nano, per le varibili a dominio limitato, gli estremi degli
 *  intervalli nota la probabilita' di avere valore ivi contenuto.
 *
 *  doubleizzata da           Maurizio Mattia.
 *  Bibliografia       Numerical Recipies in C
 *  Versione              0.1, 18 agosto 1997.
 */


#define _TEST_ERFLIB  /* Se definita compila il main che testa le funzioni. */



/**
 *  Returns the complementary error function erfc(x) with
 *  fractional error everywhere less than 1.2 x 10^-7.
 */
double erflibClass::erfcc(double x)
{
   double t, z, ans;

   z = fabs(x);
   t = 1.0/(1.0+0.5*z);
   ans = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
         t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
         t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
   return x >= 0.0 ? ans : 2.0-ans;
}



/**
 *  Cumulative density function of normal distribution with
 *  variance sigma^2 and mean mu.
 */
double erflibClass::normalCumulative (double x, double mu, double sigma)
{
   return 1-0.5*erfcc((x-mu)/(sqrt(2.0)*sigma));
}


/**
 *  Normal distribution with variance sigma^2 and mean mu.
 */
double erflibClass::normal (double x, double mu, double sigma)
{
   return 1/(sqrt(3.1415927*2)*sigma)*exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma));
}


/**
 *  Cumulative density function for a random variable with
 *  relative probabilies equals to a normal variable with
 *  mean mu and variance sigma^2, and limited domain [a,b].
 */
double erflibClass::cutNormalCumulative (double x,
                          double mu, double sigma,
                          double a, double b)
{
   double n0;

   n0 = normalCumulative(a, mu, sigma);

   return (normalCumulative(x, mu, sigma)-n0)/(normalCumulative(b, mu, sigma)-n0);
}


/**
 *  Find the x value which gives cutNormalCumulative p.
 */
double erflibClass::findProbability (double p,
                      double mu, double sigma,
                      double a, double b)
{
   double        Error = 1e-6;
   int MaxIterations = 256;

   int  i=0;
   double x, xmin, xmax;

   if (p <= 0.0) return a;
   if (p >= 1.0) return b;

   xmin = a; xmax = b;
   while ((xmax-xmin) > Error && i<MaxIterations)
   {
      x = (xmax+xmin)*0.5;
      if (cutNormalCumulative(x, mu, sigma, a, b) < p)
         xmin = x;
      else
         xmax = x;
      i++;
   }

   if (i == MaxIterations)
      fprintf(stderr, "Iterazioni: %i. Errore: %g.\n", i, xmax-xmin);

   return (xmax+xmin)*0.5;
}


/**
 *   Builds or updates a generic look-up table (LUT) for the 
 *   off-line Monte Carlo. If needed allocates the memory for 
 *   the table. Table is a reference to an array that will 
 *   have IntervalNum elements. The gaussian distribution has 
 *   mean mu and standard deviation sigma, and the range of 
 *   the elements is [xmin,xmax].
 *   If Table is not a NULL pointer the LUT is not allocated
 *   and will be used the memory pointed by it, without any 
 *   control. For this reason the previous size of the LUT 
 *   have to match the new one in order to avoid any program 
 *   crash.
 *   Returns 0 if the function builds correctly the LUT, 
 *   1 otherwise.
 */

void erflibClass::makeGaussianLUT (double *Table, int IntervalNum,
                     double mu, double sigma,
                     double xmin, double xmax)
{
   double p, dp;
   double x0, x1;
   int i;

   /*** Allocates memory, if needed. ***/
   //if ((*Table) == NULL)
   //   if (!(*Table = malloc(sizeof(double)*IntervalNum))) return 1;

   /*** Table initialization. ***/
   if (sigma > 0.0) {

      dp = 1.0/(double)IntervalNum;
      x1 = xmin;
      for (p=dp, i=0; p<=1.0; p+=dp, i++)
      {
         x0 = x1;
         x1 = findProbability(p, mu, sigma, xmin, xmax);
         Table[i] = (x1+x0)*0.5;
      }
   } else

      for (i=0; i<IntervalNum; i++)
         Table[i] = mu;

}


uint32_t erflibClass::roundr2i (double r)
{
   static double i;

   i = floor(r);
   if (i+0.5 > r) 
      return (uint32_t)i;

   return (uint32_t)(i+1.0);
}


#ifdef TEST_ERFLIB


/**
 *  Main function of test program.
 */
int main ()
{
   double xmin=0.0, xmax=1.0;
   int  xsteps=256;
   double mu=1.0, sigma=0.25;
   double x0, x1;
   
   double x, dx;
   dx =(xmax-xmin)/xsteps;

/*
   for (x=xmin; x<=xmax; x+=dx)
      printf("%f %f\n", x,normalCumulative(x, mu,sigma));
*/
/*
   for (x=xmin; x<=xmax; x+=dx)
      printf("%f %f %f\n", x,
                           normalCumulative(x+dx, mu,sigma)-normalCumulative(x, mu,sigma),
                           normal(x+dx*0.5, mu, sigma)*dx);
*/
   /*** Deve venire fuori una gaussiana con media mu e d.s. sigma. ***/
   x1 = 0;
   for (x=xmin+dx, x0=xmin; x<=xmax; x+=dx)
   {
      x0 = x1;
      x1 = findProbability(x, mu, sigma, 0.0, 2.0);
      printf("%f %f %f\n", (x1+x0)*0.5, dx/(x1-x0), normal((x1+x0)*0.5, mu, sigma));
   }

   return 0;
}

#endif
