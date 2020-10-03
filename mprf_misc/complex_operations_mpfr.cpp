#include "complex_operations_mpfr.h"

void complexadd(complex_type &c, const complex_type &a, const complex_type &b, int PRECISION)
{

	//	c.real.setPrecision(PRECISION,MPFR_RNDN);
	//	c.imag.setPrecision(PRECISION,MPFR_RNDN);	

	c.real = a.real + b.real;
	c.imag = a.imag + b.imag;

}


void complexsub(complex_type &c, const complex_type &a, const complex_type &b, int PRECISION)
{

	//	c.real.setPrecision(PRECISION,MPFR_RNDN);
	//	c.imag.setPrecision(PRECISION,MPFR_RNDN);	

	c.real = a.real - b.real;
	c.imag = a.imag - b.imag;

}


void complexmul(complex_type &c, const complex_type &a, const complex_type &b, int PRECISION)
{

	//	c.real.setPrecision(PRECISION,MPFR_RNDN);
	//	c.imag.setPrecision(PRECISION,MPFR_RNDN);	

	c.real = (a.real*b.real) - (a.imag*b.imag);
	c.imag = (a.real*b.imag) + (a.imag*b.real);
}



void complexset(complex_type &c, const complex_type &a, int PRECISION)
{

	//	c.real.setPrecision(PRECISION,MPFR_RNDN);
	//	c.imag.setPrecision(PRECISION,MPFR_RNDN);	

	c.real = a.real;
	c.imag = a.imag;
}
