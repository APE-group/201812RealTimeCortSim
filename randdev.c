// randdev.c
// File imported from Perseo

/*
 *
 *   randdev.c
 *
 *     Library providing functions generating pseudo-random
 *   number from different probability distributions (RANDom
 *   DEViations).
 *
 *   Realized by Maurizio Mattia
 *   Started November 21st, 1993
 *
 */



#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "randdev.h"
#include "DPSNN_debug.h"


/*---------------------*
 *  LOCAL DEFINITIONS  *
 *---------------------*/

static int32_t    TimeSeed = 0;      /* Seme funzione del tempo del generatore. */
static int32_t rand49_idum = -77531; /* Varibile di appoggio per gli algoritmi. */
static uint32_t   localProc;
static uint64_t counter = 0;

/*--------------------*
 *  GLOBAL FUNCTIONS  *
 *--------------------*/

#ifdef RAND2

/*------------------------------------------------*
 *                                                *
 *   rand2 e Random                               *
 *                                                *
 *   Riporta  un float casuale  nell'intervallo   *
 *   [0,1[ con DISTRIBUZIONE UNIFORME. L'algo-    *
 *   ritmo usato e' da rand2 nel Numerical Reci-  *
 *   pice.                                        *
 *   NOTA: gli int devono essere a 32bit (Watcom  *
 *         e Unix, no Borland).                   *
 *------------------------------------------------*/

#define M 714025
#define IA 1366
#define IC 150889

double randdevClass::rand2 (int32_t * idum) {

	static int64_t iy,ir[98];
	static int32_t iff=0;
	int32_t j;

	if (*idum < 0  || iff == 0) {
		iff=1;
		if ((*idum=(IC-(*idum)) % M) < 0) *idum = -(*idum);
		for (j=1;j<=97;j++) {
			*idum=(IA*(*idum)+IC) % M;
			ir[j]=(*idum);
		}
		*idum=(IA*(*idum)+IC) % M;
		iy=(*idum);
	}
	j=(int32_t)(1 + 97.0*iy/M);
	if (j > 97 || j < 1) printf("RAN2: This cannot happen.");
	iy=ir[j];
	*idum=(IA*(*idum)+IC) % M;
	ir[j]=(*idum);
	return (double) iy/M;
}

#undef M
#undef IA
#undef IC

double randdevClass::Random (void)  { return rand2(&rand49_idum); }

#endif



#ifdef RAND3

/*------------------------------------------------------*
 *                                                      *
 *   rand3 e Random                                     *
 *                                                      *
 *   Riporta  un float casuale  nell'intervallo         *
 *   [0,1[ con DISTRIBUZIONE UNIFORME.                  *
 *   E' l'algoritmo preso da numerical recipice dovuto  *
 *   a knuth e nominato ran3. E' molto veloce ma poco   *
 *   studiato (non sanno bene che tipo di correlazioni  *
 *   puo' indurre).                                     *
 *   NOTA: gli int devono essere a 32bit (Watcom        *
 *         e Unix, no Borland). Altrimenti dove e' com- *
 *         mentato 'long' inserire 'long'.              *
 *------------------------------------------------------*/

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC ((1.0/(double)MBIG))


double randdevClass::rand3(int32_t * idum)
{
   static int32_t  inext, inextp;
   static int32_t  ma[56];        /*** long ***/
   static int32_t  iff = 0;
          int32_t  mj, mk;        /*** long ***/
          int32_t  i, ii, k;

/*** Da togliere in fase di debug.
   double out;
***/
	  
   /*** Percorso di inizializzazione. ***/
   if (*idum < 0 || iff == 0) {
      iff=1;
      mj=MSEED-(*idum < 0 ? -*idum : *idum);

      /*** Modifica mia. ***/
      mj = (mj < 0 ? -mj : mj); /*** Questa aggiunta l'ho fatta io, perche' ***
                                 *** se non ci fosse non vengono numeri tra ***
                                 *** 0.0 e 1.0. Io penso che in NR ci sia   ***
                                 *** un errore di interpretazione sul si-   ***
                                 *** gnificato di '%': pensano che riporti  ***
                                 *** sempre un numero positivo.             ***
                                 *** Se non ci credi prova con questo seed  ***
                                 *** 878104453 (da passare a setRanomSeed). ***/
      mj %= MBIG;
      ma[55]=mj;
      mk=1;
      for (i=1;i<=54;i++) {
         ii=(21*i) % 55;
         ma[ii]=mk;
         mk=mj-mk;
         if (mk < MZ) mk += MBIG;
         mj=ma[ii];
      }
      for (k=1;k<=4;k++)
         for (i=1;i<=55;i++) {
            ma[i] -= ma[1+(i+30) % 55];
            if (ma[i] < MZ) ma[i] += MBIG;
         }
      inext=0;
      inextp=31;
      *idum=1;
   }

   if (++inext == 56) inext=1;
   if (++inextp == 56) inextp=1;
   mj=ma[inext]-ma[inextp];
   if (mj < MZ) mj += MBIG;
   ma[inext]=mj;

   counter++;
/*** Da scommentare in fase di debug.
   out = mj*FAC;
   if (out < 0 || out >= 1) {
      fprintf(stderr, "ERROR: Bad random number.");
      abort();
   }
   return out;
***/
   DPSNNverboseStart(true,1,0);
   {
     double result;
     result = mj*FAC;
     if((result < 0.0) || (result >= 1.0)){
       printf("ERROR in random generation!!\n");
       fflush(stdout);exit(0);
     }
   }
   DPSNNverboseEnd();

   return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

double randdevClass::Random ()  { return rand3(&rand49_idum); }

uint64_t randdevClass::Statistics ()  { return counter; }
#endif /* of RAND3 */


#ifdef RAND4

/*------------------------------------------------------*
 *                                                      *
 *   rand4 e Random                                     *
 *                                                      *
 *   Riporta  un float casuale  nell'intervallo         *
 *   [0,1[ con DISTRIBUZIONE UNIFORME.                  *
 *   E' l'algoritmo preso da numerical recipice dovuto  *
 *   a knuth e nominato ran3. E' molto veloce ma poco   *
 *   studiato (non sanno bene che tipo di correlazioni  *
 *   puo' indurre).                                     *
 *   NOTA: gli int devono essere a 32bit (Watcom        *
 *         e Unix, no Borland). Altrimenti dove e' com- *
 *         mentato 'long' inserire 'long'.              *
 *------------------------------------------------------*/

//#define MBIG 1000000000
//#define MSEED 161803398
#define MZ 0
//#define FAC ((1.0/(double)MBIG))


double randdevClass::rand4(int32_t * idum)
{
   static int64_t  MBIG;  
   static double   FAC;
   static int64_t  MSEED;
   static int64_t  inext, inextp;
   static int64_t  ma[56];        /*** long ***/
   static int64_t  iff = 0;
          int64_t  mj, mk;        /*** long ***/
          int64_t  i, ii, k;
	  uint32_t numerator_uint32;
	  double   out_double;
   

/*** Da togliere in fase di debug.
   double out;
***/
	  
   /*** Percorso di inizializzazione. ***/
   if (*idum < 0 || iff == 0) {
     MBIG = (int64_t)1000000000<<32;
     printf("MBIG = %lld \n",MBIG);
     MSEED = ((int64_t)161803398 << 32) + 0xB76539E;
     printf("MSEED = %lld \n",MSEED);
     FAC = (double)1.0/(double)MBIG;
     printf("FAC = %e \n",FAC);
     iff=1;
      mj=MSEED-(*idum < 0 ? -*idum : *idum);

      /*** Modifica mia. ***/
      mj = (mj < 0 ? -mj : mj); /*** Questa aggiunta l'ho fatta io, perche' ***
                                 *** se non ci fosse non vengono numeri tra ***
                                 *** 0.0 e 1.0. Io penso che in NR ci sia   ***
                                 *** un errore di interpretazione sul si-   ***
                                 *** gnificato di '%': pensano che riporti  ***
                                 *** sempre un numero positivo.             ***
                                 *** Se non ci credi prova con questo seed  ***
                                 *** 878104453 (da passare a setRanomSeed). ***/
      mj %= MBIG;
      ma[55]=mj;
      mk=1;
      for (i=1;i<=54;i++) {
         ii=(21*i) % 55;
         ma[ii]=mk;
         mk=mj-mk;
         if (mk < MZ) mk += MBIG;
         mj=ma[ii];
      }
      for (k=1;k<=4;k++)
         for (i=1;i<=55;i++) {
            ma[i] -= ma[1+(i+30) % 55];
            if (ma[i] < MZ) ma[i] += MBIG;
         }
      inext=0;
      inextp=31;
      *idum=1;
   }

   if (++inext == 56) inext=1;
   if (++inextp == 56) inextp=1;
   mj=ma[inext]-ma[inextp];
   if (mj < MZ) mj += MBIG;
   ma[inext]=mj;

   numerator_uint32 = (uint32_t)((mj >> 16) && 0xFFFFFFFF);
   out_double = ((double)numerator_uint32/(double)((uint32_t)0xFFFFFFFF));

   DPSNNverboseStart(true,1,0);
     if((out_double < 0.0) || (out_double >= 1.0)){
       printf("ERROR in random generation!!\n");
       fflush(stdout);exit(0);
     }
   DPSNNverboseEnd();

   return out_double;

}

#undef MZ

double randdevClass::Random ()  { return rand4(&rand49_idum); }

#endif /* of RAND4 */


#ifdef RAND5

/*------------------------------------------------------*
 *                                                      *
 *   rand5 e Random                                     *
 *                                                      *
 *   Riporta  un float casuale  nell'intervallo         *
 *   [0,1[ con DISTRIBUZIONE UNIFORME.                  *
 *   E' l'algoritmo preso da numerical recipice dovuto  *
 *   a knuth e nominato ran3. E' molto veloce ma poco   *
 *   studiato (non sanno bene che tipo di correlazioni  *
 *   puo' indurre).                                     *
 *   NOTA: gli int devono essere a 32bit (Watcom        *
 *         e Unix, no Borland). Altrimenti dove e' com- *
 *         mentato 'long' inserire 'long'.              *
 *------------------------------------------------------*/

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC ((1.0/(double)MBIG))
#define MAX_RAND_H 8

double randdevClass::rand5(int32_t * idum)
{
   static int32_t  inext, inextp;
   static int32_t  ma[56];        /*** long ***/
   static int32_t  iff = 0;
          int32_t  mj, mk;        /*** long ***/
          int32_t  i, ii, k;
	  double   outRand[MAX_RAND_H];
	  
   /*** Percorso di inizializzazione. ***/
   if (*idum < 0 || iff == 0) {
      iff=1;
      mj=MSEED-(*idum < 0 ? -*idum : *idum);

      /*** Modifica mia. ***/
      mj = (mj < 0 ? -mj : mj); /*** Questa aggiunta l'ho fatta io, perche' ***
                                 *** se non ci fosse non vengono numeri tra ***
                                 *** 0.0 e 1.0. Io penso che in NR ci sia   ***
                                 *** un errore di interpretazione sul si-   ***
                                 *** gnificato di '%': pensano che riporti  ***
                                 *** sempre un numero positivo.             ***
                                 *** Se non ci credi prova con questo seed  ***
                                 *** 878104453 (da passare a setRanomSeed). ***/
      mj %= MBIG;
      ma[55]=mj;
      mk=1;
      for (i=1;i<=54;i++) {
         ii=(21*i) % 55;
         ma[ii]=mk;
         mk=mj-mk;
         if (mk < MZ) mk += MBIG;
         mj=ma[ii];
      }
      for (k=1;k<=4;k++)
         for (i=1;i<=55;i++) {
            ma[i] -= ma[1+(i+30) % 55];
            if (ma[i] < MZ) ma[i] += MBIG;
         }
      inext=0;
      inextp=31;
      *idum=1;
   }

   for(i=0;i<MAX_RAND_H;i++){
     if (++inext == 56) inext=1;
     if (++inextp == 56) inextp=1;
     mj=ma[inext]-ma[inextp];
     if (mj < MZ) mj += MBIG;
     ma[inext]=mj;
     outRand[i]=mj*FAC;
   }

   DPSNNverboseStart(true,1,0);
   {
     double result;
     result = outRand[localProc];
     if((result < 0.0) || (result >= 1.0)){
       printf("ERROR in random generation!!\n");
       fflush(stdout);exit(0);
     }
   }
   DPSNNverboseEnd();

   return outRand[localProc];
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

double randdevClass::Random ()  { return rand5(&rand49_idum); }

#endif /* of RAND5 */


#ifdef RAND1

/*------------------------------------------------*
 *                                                *
 *   UniDev.                                      *
 *                                                *
 *   Riporta  un  int  casuale  nell'intervallo   *
 *   [0,RAND_MAX] con DISTRIBUZIONE UNIFORME.     *
 *   Per inizializzare  il seme del  generatore   *
 *   utilizzare la seguente sintassi:             *
 *      UniDev( &Seed )                           *
 *   dove Seed deve essere una variabile intera   *
 *   contenente il seme iniziale. Per mantenere   *
 *   intatto il seme basta richiamare la funzio-  *
 *   ne con la seguente sintassi:                 *
 *      UniDev( NULL )                            *
 *   Questa procedura si basa sull'algoritmo di   *
 *   BAYS e DURHAM.                               *
 *                                                *
 *   NOTE: Il periodo del  generatore e' pratica- *
 *         infinito, ma i numeri forniti non so-  *
 *         no piu' di RAND_MAX+1.                 *
 *         Deve essere usata con i piedi di piom- *
 *         bo quando si utilizzano i bit meno si- *
 *         gnificativi.                           *
 *                                                *
 *------------------------------------------------*/

int32_t randdevClass::UniDev (int32_t * SeedPntr) /* Puntatore al seme del generatore. Se NULL il generatore *
                             * non viene inizializzato.                                */

{
   /*** Dichiarazione delle variabili locali. ***/
   #define SEQUENCE_DIM 98.0
   static int32_t y, Sequence[(int)SEQUENCE_DIM];
   static int32_t FirstCall = 1;         /* Prima chiamata della funzione? */
   static int32_t j;                     /* Contatore.                     */

   /*** Inizializzazione della sequenza. ***/
   if (SeedPntr != NULL || FirstCall)
     {
      FirstCall = 0;
      if (SeedPntr != NULL) srand(*SeedPntr);
      for (j=0; j<(int)SEQUENCE_DIM; j++) rand();
      for (j=0; j<(int)SEQUENCE_DIM; j++) Sequence[j] = rand();
      y = rand();
     }

   /*** Generazione del numero ps eudo-casuale. ***/
   j           = (SEQUENCE_DIM*y)/(RAND_MAX+1);
   y           = Sequence[j];
   Sequence[j] = rand();
   return y;
}

#endif /* of RAND1 */



/*----------------------------------------------*
 *                                              *
 *   SetRandomSeed (int Seed)                   *
 *                                              *
 *   Assegna un valore al seme del generatore   *
 *   numeri pseudo-casuali.                     *
 *                                              *
 *----------------------------------------------*/

void randdevClass::SetRandomSeed (int32_t Seed)

{
   TimeSeed = Seed;
   
#ifdef RAND2
   rand49_idum=-Seed;
   rand2(&rand49_idum);
#endif
#ifdef RAND3
   rand49_idum=-Seed;
   rand3(&rand49_idum);
#endif
#ifdef RAND4
   rand49_idum=-Seed;
   rand4(&rand49_idum);
#endif
#ifdef RAND5
   rand49_idum=-Seed;
   localProc = h;
   rand5(&rand49_idum);
#endif
#ifdef RAND1
   UniDev(&TimeSeed);
#endif
}



/*------------------------------------------------------*
 *   GetRandomSeed ()                                   *
 *                                                      *
 *   Riporta il seme del generatore di numeri pseudo-   *
 *   casuali TimeSeed.                                  *
 *------------------------------------------------------*/

int32_t randdevClass::GetRandomSeed (void)

{
   return TimeSeed;
}



/*----------------------------------------------*
 *                                              *
 *   Randomize.                                 *
 *                                              *
 *   Inizializza in modo casuale (attraverso    *
 *   la funzione time) il generatore di nume-   *
 *   ri pseudo-casuali UniDev.                  *
 *                                              *
 *----------------------------------------------*/

void randdevClass::Randomize ()

{
  SetRandomSeed((unsigned)time(NULL));
}



/*------------------------------------------------*
 *                                                *
 *   ExpDev.                                      *
 *                                                *
 *   Riporta un numero  casuale reale positivo    *
 *   (float) con  distribuzione esponenziale e    *
 *   media  unitaria, usando come  sorgente di    *
 *   numeri casuali con distribuzione uniforme    *
 *   Random.                                      *
 *     Fun. di distribuzione = e^-x con x>0.0    *
 *   Questa procedura si basa sull'algoritmo di   *
 *   Montecarlo.                                  *
 *                                                *
 *------------------------------------------------*/

double randdevClass::ExpDev ( void )

{
   /*** Dichiarazione delle variabili locali. ***/
   double Dum;      /* Variabile di appoggio. */

   /*** Ciclo di generazione del singolo numero. ***/
   do
      Dum = (double)Random();
   while (Dum == 0.0);
   return (double)-log(Dum);
}



/*------------------------------------------------*
 *                                                *
 *   NormDev.                                     *
 *                                                *
 *   Riporta un  numero  casuale reale (float)    *
 *   con distribuzione gaussiana,  media nulla    *
 *   e varianza unitaria. Come sorgente di nu-    *
 *   meri  casuali con  distribuzione uniforme    *
 *   viene usata Random().                        *
 *     Fun. di distr. = 1/(sqrt(2pi))*e^-x^2      *
 *   Questa procedura si basa sull'algoritmo di   *
 *   BOX-MULLER.                                  *
 *                                                *
 *------------------------------------------------*/

double randdevClass::NormDev ( void )

{
   /*** Dichiarazione delle variabili locali. ***/
   static int32_t   Set = 0;
   static double IIran;
   double Fac, r, v1, v2;

   /*** Decisione su quale set prendere il singolo numero. ***/
   if (Set == 0)
     {do
        {v1 = (double)(2.0 * Random() - 1.0);
         v2 = (double)(2.0 * Random() - 1.0);
         r  = v1 * v1 + v2 * v2;}
      while (r >= 1.0 || r==0.0);
      Fac   = (double)sqrt(-2.0 * log(r) / r);
      IIran = v1 * Fac;
      Set   = 1;
      return v2 * Fac;}
   else
     {Set = 0;
      return IIran;}
}



/*------------------------------------------------*
 *                                                *
 *   BernDev.                                     *
 *                                                *
 *   Riporta un  numero  casuale intero (int)     *
 *   con distribuzione Bernoulliana, con media    *
 *   N*p e varianza N*p*(1-p). Come sorgente di   *
 *   numeri casuali con distribuzione uniforme    *
 *   viene usata Random().                        *
 *                                                *
 *------------------------------------------------*/

int32_t randdevClass::BernDev (int32_t N, double p) 

{
   static double FDP;     /* Funzione di Distribuzione di Probabilita'. */
   static double FDC;     /* Funzione di Distribuzione Cumulativa.      */
   static double r;       /* Numero casuale in [0,1[.                   */
   static int32_t n;       /* Numero di eventi da riportare.             */
   static double q;       /* Prob. di NON emissione di uno spike.       */
   static double Precision = 0.00001;

   /* Iniz. delle Variab. utili a calcolare la distrib. di Bernoulli. */
   q = 1.0 - p;
   FDC = FDP = (double)(exp( N * log(q) ));
   r   = (double)Random();
   n   = 0;

   /*** Questa riga e' necessario poiche' nel ciclo successivo. ***
    *** FDC non tende a 1.0 bensi a qualcosa dell'ordine di     ***
    *** 1.0-Precision. Altrimenti si possono avere dei loop     ***
    *** infiniti.                                               ***/
   if ((1.0-r)<Precision) r = (double)(1.0-Precision);

   /* Estrazione del numero di eventi. */
   while (r>=FDC) {
      n++;
      FDP  *= (double)((N-n+1.)/n * p / q);
      FDC += FDP;
   }

   return n;
}
