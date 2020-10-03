#include "our_fft.h"

void mpreal_fft_sygrisi(complex_type *F, complex_type *DL_F, long int NN, int PRECISION, int DL_PRECISION)
{
	/* !!! Attention we must subtract 1 each time we use the matrix F
	in order to avoid the shift at the beggining and the end of
	the function !!!!! */
	/*  compl_type *F; */

	complex_type U, W, T;
	mpreal PI, endiamesi;

	complex_type DL_U, DL_W, DL_T;
	mpreal DL_PI, DL_endiamesi;

	long int I, IP, J, K, L, LE, LE1, NV2, NM1;
	long int p, kk, n;

	/*SET PRECISION*/
	PI = const_pi(PRECISION);
	endiamesi.setPrecision(PRECISION, MPFR_RNDN);

	U.real.setPrecision(PRECISION, MPFR_RNDN);
	U.imag.setPrecision(PRECISION, MPFR_RNDN);
	W.real.setPrecision(PRECISION, MPFR_RNDN);
	W.imag.setPrecision(PRECISION, MPFR_RNDN);
	T.real.setPrecision(PRECISION, MPFR_RNDN);
	T.imag.setPrecision(PRECISION, MPFR_RNDN);

	/*SET DOUBLE (DL) PREICISION*/
	DL_PI = const_pi(DL_PRECISION);
	DL_endiamesi.setPrecision(DL_PRECISION, MPFR_RNDN);

	DL_U.real.setPrecision(DL_PRECISION, MPFR_RNDN);
	DL_U.imag.setPrecision(DL_PRECISION, MPFR_RNDN);
	DL_W.real.setPrecision(DL_PRECISION, MPFR_RNDN);
	DL_W.imag.setPrecision(DL_PRECISION, MPFR_RNDN);
	DL_T.real.setPrecision(DL_PRECISION, MPFR_RNDN);
	DL_T.imag.setPrecision(DL_PRECISION, MPFR_RNDN);


	/*   for(index=NN;index>0;index--)
	{
	F[index].real=F[index-1].real;
	F[index].imag=F[index-1].imag;
	}
	*/
	/* Below Pointer Operation to avoid Shift  */
	/*F=(F1-1);*/

	/* THIS HOPEFULLY! CHECKS IF THE WINDOW IS A POWER OF 2 */
	kk = NN;
	n = 0;
	while (kk > 1)
	{
		if (kk % 2 != 0)
		{
			exit(0);  /* IN THE FUTURE: FILL IN WITH ZEROS!!!! */
		}
		kk = kk / 2;
		n++;
	}
	/*    END OF THIS HOPEFULLY! CHECKS IF THE WINDOW IS A POWER OF 2 */

	/*   NN=1;
	for(p=1;p<=n;p++)
	{
	NN=2*NN;
	}
	*/

	NV2 = NN / 2;
	NM1 = NN - 1;
	J = 1;
	for (I = 1; I <= NM1; I++)
	{
		if (I < J)
		{
			//Single Precision
			complexset(T, F[J - 1], PRECISION);
			complexset(F[J - 1], F[I - 1], PRECISION);
			complexset(F[I - 1], T, PRECISION);

			//Double Precision
			complexset(DL_T, DL_F[J - 1], DL_PRECISION);
			complexset(DL_F[J - 1], DL_F[I - 1], DL_PRECISION);
			complexset(DL_F[I - 1], DL_T, DL_PRECISION);
			
			//T = F[J - 1];
			//F[J - 1] = F[I - 1];
			//F[I - 1] = T;
		}
		K = NV2;
		while (K < J)
		{
			J = J - K;
			K = K / 2;
		}
		J = J + K;

	} /* FOR  */

	for (L = 1; L <= n; L++)
	{
		LE = 1;
		for (p = 1; p <= L; p++)
		{
			LE = 2 * LE;
		}
		LE1 = LE / 2;

		U.real = 1.0;
		U.imag = 0.0;
		DL_U.real = 1.0;
		DL_U.imag = 0.0;

		endiamesi = PI / (1.0*LE1);
		DL_endiamesi = DL_PI / (1.0*LE1);

		W.real = cos(endiamesi);
		W.imag = (-1.0)*sin(endiamesi);
		DL_W.real = cos(DL_endiamesi);
		DL_W.imag = (-1.0)*sin(DL_endiamesi);

		for (J = 1; J <= LE1; J++)
		{
			for (I = J; I <= NN; I = I + LE)
			{
				IP = I + LE1;

				complexmul(T, F[IP - 1], U, PRECISION);
				complexsub(F[IP - 1], F[I - 1], T, PRECISION);
				complexadd(F[I - 1], F[I - 1], T, PRECISION);

				complexmul(DL_T, F[IP - 1], DL_U, DL_PRECISION);
				complexsub(DL_F[IP - 1], DL_F[I - 1], DL_T, DL_PRECISION);
				complexadd(DL_F[I - 1], DL_F[I - 1], DL_T, DL_PRECISION);
			}
			complexmul(U, U, W, PRECISION);
			complexmul(DL_U, DL_U, DL_W, DL_PRECISION);
		}
	}

	/*  THIS IS NECESSARY ONLY FOR THE INVERSE FFT TRANSFORM
	for(I=1;I<=NN;I++)
	{
	F[I-1].real=F[I-1].real/NN;
	F[I-1].imag=F[I-1].imag/NN;
	}
	*/
	/*   for(index=0;index<NN;index++)
	{
	F[index].real=F[index+1].real;
	F[index].imag=F[index+1].imag;
	} */


}


