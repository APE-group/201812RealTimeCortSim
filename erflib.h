// erflib.h
// File imported from Perseo 

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

#include <stdint.h>

#ifndef erflibIncluded
#define erflibIncluded


class erflibClass {
 public:

/**
 *  Returns the complementary error function erfc(x) with
 *  fractional error everywhere less than 1.2 x 10^-7.
 */
double erfcc (double x);


/**
 *  Cumulative density function of normal distribution with
 *  variance sigma^2 and mean mu.
 */
double normalCumulative (double x, double mu, double sigma);


/**
 *  Normal distribution with variance sigma^2 and mean mu.
 */
double normal (double x, double mu, double sigma);


/**
 *  Cumulative density function for a random variable with
 *  relative probabilies equals to a normal variable with
 *  mean mu and variance sigma^2, and limited domain [a,b].
 */
double cutNormalCumulative (double x,
                          double mu, double sigma,
                          double a, double b);


/**
 *  Find the x value which gives cutNormalCumulative p.
 */
double findProbability (double p,
                      double mu, double sigma,
                      double a, double b);


/**
 *   Builds a generic look-up table (LUT) for the off-line 
 *   montecarlo. Allocates the memory for the table. Table is 
 *   a reference to an array that will have IntervalNum 
 *   elements. The gaussian distribution has mean mu and 
 *   standard deviation sigma, and the range of the elements is 
 *   [xmin,xmax].
 *   Returns 0 if the function builds correctly the LUT, 
 *   1 otherwise.
 */

void makeGaussianLUT (double *Table, int IntervalNum,
                     double mu, double sigma,
                     double xmin, double xmax);


/*
 *  Return the closest integer to <r>.
 *  This function has been taken from nalib.h file in PERSEO
 *
 */
 uint32_t roundr2i (double r);

};

#endif
