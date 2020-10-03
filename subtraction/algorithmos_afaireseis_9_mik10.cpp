/*
 * algorithmos_pollaplasiasmou_1.cpp
 *
 *  Created on: Sep 8, 2014
 *      Author: costas
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmpxx.h>
#include <mpfr.h>
#include <mpreal.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


using namespace std;
using namespace mpfr;

#define EPANAL 200000000


int strset_1(char  *num, char c, int length);
int mpfr_find_eddv2(mpfr_t Num_n, mpfr_t Num_nc);

#define SINGLE_PREC 16
#define DOUBLE_PREC 32


int main()
{
	mpfr_prec_t bit_precision , dl_bit_precision;
	mpfr_t fz,fx,fs;
	mpfr_t dz,dx,ds;
	mpfr_t arithmos,d_arithmos;
	char s[]="2.02"; //syntelesths
	long int ii,error,jj;

	bit_precision = digits2bits(SINGLE_PREC);
	dl_bit_precision = digits2bits(DOUBLE_PREC);

	mpfr_init2(fz,bit_precision);
	mpfr_init2(fx,bit_precision);
	mpfr_init2(fs,bit_precision);

	mpfr_init2(dz,dl_bit_precision);
	mpfr_init2(dx,dl_bit_precision);
	mpfr_init2(ds,dl_bit_precision);

	mpfr_init2(arithmos,bit_precision);
	mpfr_init2(d_arithmos,dl_bit_precision);


	mpfr_set_str(fx,"1.0333342432432",10,MPFR_RNDN);
	mpfr_set_str(dx,"1.0333342432432",10,MPFR_RNDN);

	mpfr_set_str(fs,s,10,MPFR_RNDN);
	mpfr_set_str(ds,s,10,MPFR_RNDN);

	mpfr_set(fz,fx,MPFR_RNDN);
	mpfr_set(dz,dx,MPFR_RNDN);

	for(ii=0 ; ii<EPANAL ; ii++)
	{


		mpfr_mul(fx,fz,fs,MPFR_RNDN);
		mpfr_mul(dx,dz,ds,MPFR_RNDN);




		if (mpfr_cmp_d(dx,10.0)<0)
		{
			mpfr_set(fz,fx,MPFR_RNDN);
			mpfr_set(dz,dx,MPFR_RNDN);
/*
			mpfr_fprintf(stderr,"%d %d :: \n",error,ii);
								mpfr_fprintf(stderr,"%.*Re\n",SINGLE_PREC-1,fz);
								mpfr_fprintf(stderr,"%.*Re\n",SINGLE_PREC-1,dz);
								getchar();
*/
			 if (mpfr_cmp_si(fz,10)>0)
							 {
								 mpfr_div_si(fz,fz,10,MPFR_RNDN);
								 mpfr_div_si(dz,dz,10,MPFR_RNDN);
							 }

		}
		else
		{

			//mpfr_sub_d(<metavliti oppou_tha_ekchorhthei_to_apotelesma>,<metablhth_mpfr>,<stathera>,MPFR_RNDN)
			//afairvw kai apo tis metavlites typou f,d (gia diaforetiko precision).

			//mpfr_div_si(fx,fz,4.9,MPFR_RNDN);
			//mpfr_div_si(dx,dz,4.9,MPFR_RNDN);

	mpfr_set_str(arithmos,"3.800000",10,MPFR_RNDN);
	mpfr_set_str(d_arithmos,"3.8000000",10,MPFR_RNDN);
			
mpfr_sub(fz,fz,arithmos,MPFR_RNDN);
mpfr_sub(dz,dz,d_arithmos,MPFR_RNDN);

			fprintf(stderr,"Afairesh!\n");
			//mpfr_fprintf(stderr,"error:%d ii:%d \n",error,ii);
			//mpfr_fprintf(stderr,"%.*Re\n",SINGLE_PREC-1,fz);
			//mpfr_fprintf(stderr,"%.*Re\n\n",SINGLE_PREC-1,dz);
			//getchar();

//			mpfr_add_si(fz,fz,8,MPFR_RNDN);
//			mpfr_add_si(dz,dz,8,MPFR_RNDN);

//			mpfr_mul(fz,fz,fs,MPFR_RNDN);
//			mpfr_mul(dz,dz,ds,MPFR_RNDN);
		
//			 if (mpfr_cmp_si(fz,10)>0)
//			 {
//					 mpfr_div_si(fz,fz,10,MPFR_RNDN);
//					 mpfr_div_si(dz,dz,10,MPFR_RNDN);
//			 }
//			mpfr_div_si(fz,fz,10,MPFR_RNDN);
//			mpfr_div_si(dz,dz,10,MPFR_RNDN);
//			mpfr_fprintf(stderr,"%d %d : \n",error,ii);
//									mpfr_fprintf(stderr,"%.*Re\n",SINGLE_PREC-1,fz);
//									mpfr_fprintf(stderr,"%.*Re\n",SINGLE_PREC-1,dz);
//									getchar();

		}

		error=mpfr_find_eddv2(fz,dz);


	//if (ii>=3000000 && ii<000021  )
	//	if (ii%100000==0)
	//	{
			mpfr_fprintf(stderr,"error:%d ii:%d \n",error,ii);
			mpfr_fprintf(stderr,"%.*Re\n",SINGLE_PREC-1,fz);
			mpfr_fprintf(stderr,"%.*Re\n\n",SINGLE_PREC-1,dz);
getchar();
//		}


	}

	mpfr_clear(fz);
	mpfr_clear(fx);
	mpfr_clear(fs);

	mpfr_clear(ds);
	mpfr_clear(dz);
	mpfr_clear(dx);

}

int mpfr_find_eddv2(mpfr_t Num_n, mpfr_t Num_nc)
{
	char * char_mantissa;
	long int size_man_array,result,lathos_psifia,length,vythisi,ekth_diaf,meg;
	mpfr_prec_t  bitprecision_n,bitprecision_nc;
	mpfr_t nf,nc,ncr;
	mpfr_t diafora,zero;
	mpfr_t mpf_ekth_diaf;
	char* z_pr;
	int precision;

	//eyresi precision arithnou se bits
	bitprecision_n=mpfr_get_prec(Num_n);
	bitprecision_nc=mpfr_get_prec(Num_nc);
	precision=bits2digits(bitprecision_n);

	z_pr=(char*)malloc(sizeof(char)*(3*precision));
	strset_1(z_pr,'\0',precision*3);

	//arxikopoiisi metavlitwn pros xrisi
	mpfr_init(mpf_ekth_diaf);
	mpfr_init2(nf,bitprecision_nc+16);
	mpfr_init2(ncr,bitprecision_nc+16);
	mpfr_init2(nc,bitprecision_nc);
	mpfr_init2(diafora,bitprecision_n);
	mpfr_init2(zero,bitprecision_n);
	mpfr_set_zero(zero,0);

	//extract ektheti twn dyo arithmwn
	mpfr_abs(nf,Num_n,MPFR_RNDN);
	mpfr_log10(nf,nf,MPFR_RNDN);
	mpfr_floor(nf,nf);

	mpfr_abs(nc,Num_nc,MPFR_RNDN);
	mpfr_log10(nc,nc,MPFR_RNDN);
	mpfr_floor(nc,nc);

	result=mpfr_cmp(nf,nc);


	if (result>=0)
	{
		meg=mpfr_get_si(nf,MPFR_RNDN);
	}
	else
	{
		meg=mpfr_get_si(nc,MPFR_RNDN);
	}

	size_man_array=precision*3;
	char_mantissa=(char*)malloc(sizeof(char)*size_man_array);
	strset_1(char_mantissa,'\0',size_man_array);

	mpfr_set(nf,Num_n, MPFR_RNDN);
	mpfr_set(nc,Num_nc, MPFR_RNDN);

	//round tou arithmou sthn epithymiti precision
	mpfr_sprintf(z_pr,"%.*Re",precision-1,nc);
	mpfr_set_str(ncr,z_pr,10,MPFR_RNDA);

	mpfr_sprintf(z_pr,"%.*Re",precision-1,nf);
	mpfr_set_str(nf,z_pr,10,MPFR_RNDA);

	//eyresi diaforas twn dyo arithmwn
	mpfr_sub(diafora,ncr,nf,MPFR_RNDN);
	mpfr_abs(diafora,diafora,MPFR_RNDN);

	result=mpfr_cmp(diafora,zero);

	if (result==0)
	{

		lathos_psifia=0;

	}
	else
	{
		mpfr_sprintf(char_mantissa,"%.*Re",precision-4,diafora);
		length=strlen(char_mantissa)-1;

		while (char_mantissa[length]!='e')
			length--;

		mpfr_set_str(mpf_ekth_diaf,&(char_mantissa[length+1]),10,MPFR_RNDN);
		ekth_diaf=mpfr_get_si(mpf_ekth_diaf,MPFR_RNDN);

		vythisi=meg-ekth_diaf;
		lathos_psifia=precision-vythisi;

		if (lathos_psifia<0)
			lathos_psifia=0;

	}
	free(char_mantissa);
	free(z_pr);

	mpfr_clear(mpf_ekth_diaf);
	mpfr_clear(nf);
	mpfr_clear(nc);
	mpfr_clear(ncr);
	mpfr_clear(diafora);
	mpfr_clear(zero);

	return lathos_psifia;
}









int strset_1(char  *num, char c, int length)
{
	int i;

	for ( i=0 ; i < length ; i++)
		num[i]=c;

	return 0;

}




