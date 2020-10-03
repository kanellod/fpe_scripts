#include "complex_operations_mpfr.h"
#include "our_fft.h"

#define PREC_DOUBLE 53 	//Precision
#define PREC_32 104		//Precision
#define PREC_64 235		//Precision

int main() 
{
	long int WINDOW_LENGTH;	
	complex_type *F,*DL_F;
	
	printf("Dwse Window Length: ");
	cin >> WINDOW_LENGTH;
	
	F    = new complex_type[WINDOW_LENGTH]; 
	DL_F = new complex_type[WINDOW_LENGTH];
	
	for (long int i=0; i<WINDOW_LENGTH ; i++)
	{
		F[i].real.setPrecision(PREC_DOUBLE,MPFR_RNDN);
		F[i].imag.setPrecision(PREC_DOUBLE,MPFR_RNDN);	
		DL_F[i].real.setPrecision(PREC_32,MPFR_RNDN);
		DL_F[i].imag.setPrecision(PREC_32,MPFR_RNDN);	
	}
	
	mpreal_fft_sygrisi(F,DL_F, WINDOW_LENGTH,PREC_DOUBLE, PREC_32);
	
}


