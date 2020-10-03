#ifndef __OUR_COMPLEX_OPERATIONS_H__
#define __OUR_COMPLEX_OPERATIONS_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmpxx.h>
#include <mpfr.h>
#include <mpreal.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace mpfr;
using namespace std;


#define PREC_DOUBLE 53 	//Precision
#define PREC_32 104		//Precision
#define PREC_64 235		//Precision


typedef struct 
{
	mpreal real;
	mpreal imag;
} complex_type;

/* COMPLEX TYPES */

  void complexadd(complex_type &c, const complex_type &a,const complex_type &b, int PRECISION);
  void complexsub(complex_type &c, const complex_type &a,const complex_type &b, int PRECISION);
  void complexmul(complex_type &c, const complex_type &a,const complex_type &b, int PRECISION);
  void complexset(complex_type &c, const complex_type &a, int PRECISION);

  #endif
