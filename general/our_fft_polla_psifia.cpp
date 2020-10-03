//Edw bgazoume statistika gia to FPE ston FFT. 
//Symplhrwnw WINDOW_LENGTH (megethos parathyrou poy tha paroume apo to kommati gia na ginei ypologismos FFT),
// SIZE_TMIMA (megethos parathyrou sto opoio tha ginei ypologismos FFT)
// MAX_SIZE (plithos deigmatwn pou tha paroume apo to kommati)
// arxh_song (apo pou thelw na ksekinhsw na pairnw deigmata apo to kommati)
// thelw_64 (to kanw 1 an thelw ypologismous me 64 bits kai ypologismos error 32-64. Alliws to bazw 0)
// myfile.open ("<onoma_arxeiou>.txt"); (to arxeio pou paragetai apo to Matlab (song_fft.m) ston Leghorn kai periexei ta deigmata tou kommatiou)

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
#include "complex_operations_mpfr.h"
#include "mpfr_fpe_functions.h"

using namespace mpfr;
using namespace std;

//#define PI_OUR 3.14159265358979
#define MAX_SIZE 10000000L
#define WINDOW_LENGTH 8192L*1024L
#define SIZE_TMIMA 8192L*1024L+10L	//Megethos tmhmatos pou dinoume ston fft
#define LATHOS_BITS 32

mpreal nums_read_16[MAX_SIZE],nums_read_32[MAX_SIZE],nums_read_64[MAX_SIZE]; //oi arithmoi pou tha diavastoun apo to arxeio
mpreal PI_OUR;
compl_type F_16[SIZE_TMIMA],F_32[SIZE_TMIMA],F_64[SIZE_TMIMA]; //Auto einai to tmhma pou dinoume ston fft


int main()
{	

	long int window_length, deiktis_read;	
	/* !!! Attention we must subtract 1 each time we use the matrix F
	   in order to avoid the shift at the beggining and the end of
	   the function !!!!! */
   	compl_type U_16,W_16,T_16,U_32,W_32,T_32,U_64,W_64,T_64;
   	long int I,IP,J,K,L,LE,LE1,NV2,NM1,NN;
   	long int p,kk,n,jj,i,ll,arxh_song,ii;
	mpreal our_pi_mpreal_16, our_pi_mpreal_32, our_pi_mpreal_64;	
	string string1;	
	int error_16_32_real, error_32_64_real, error_16_32_imag, error_32_64_imag, thelw_64;
	int plithos_errors_bits_real_16_32[LATHOS_BITS];
	int plithos_errors_bits_imag_16_32[LATHOS_BITS];
 	int plithos_errors_bits_real_32_64[LATHOS_BITS];
	int plithos_errors_bits_imag_32_64[LATHOS_BITS];
	mpreal pososto_16_32_real,our_10_eis_thn_meion_14;
		
	//arxh_song=500000L;
	arxh_song=0;
	
	thelw_64=0; //An thelw ypologismous me 64 bits kai ypologismos error 32-64	
	
	our_10_eis_thn_meion_14.setPrecision(PREC_DOUBLE,MPFR_RNDN);
	our_10_eis_thn_meion_14=1.;	
	for(ff=0;ff<14;ff++)
	{
		our_10_eis_thn_meion_14=our_10_eis_thn_meion_14/10.;
	}

//////////////////////////////////////////////////////////////////////
////////////////// PROSOXH!!!! PERIOXH DOKIMWN!!!!! //////////////////
//////////////////////////////////////////////////////////////////////
/*	mpreal dokimi_16,dokimi_32,dokimi_64;
	dokimi_16.setPrecision(digits2bits(16),MPFR_RNDN);
	dokimi_16=2.345347356745863;
	
	dokimi_32.setPrecision(digits2bits(30),MPFR_RNDN);
	//dokimi_32=2.345347356745863;

	dokimi_64.setPrecision(digits2bits(64),MPFR_RNDN);
	//dokimi_64=2.3453473567458635683645645724572476324673467247624672472552457487356835838248326824672457123561357247;
	
	
	dokimi_32=dokimi_16;
	dokimi_32=dokimi_32*dokimi_32;
	dokimi_64=dokimi_16;
	dokimi_64=dokimi_64*dokimi_64;
	dokimi_16=dokimi_16*dokimi_16;
	cout.precision(64);
	cout << dokimi_16 << endl << dokimi_32 << endl << dokimi_64 << endl;
	getchar();	
*/
//////////////////////////////////////////////////////////////////////
///////////////////// TELOS PERIOXHS DOKIMWN!!!!! ////////////////////
//////////////////////////////////////////////////////////////////////

	/*string pi_32,pi_64;	
	pi_32="3.14159265358979323846264338327950";
	pi_64="3.1415926535897932384626433832795028841971693993751058209749445923";
	mpreal str_pi_mpreal_32(pi_32,PREC_32,10,MPFR_RNDN);
	mpreal str_pi_mpreal_64(pi_64,PREC_64,10,MPFR_RNDN);*/	

	//cout << "///////////////////////////////////////\n";
	our_pi_mpreal_16.setPrecision(PREC_DOUBLE,MPFR_RNDN);
	our_pi_mpreal_16=const_pi(PREC_DOUBLE,MPFR_RNDN);
	//mpfr_printf("%.*Re\n",64,our_pi_mpreal_32.mpfr_ptr());	
	
	
	//cout << "///////////////////////////////////////\n";
	our_pi_mpreal_32.setPrecision(PREC_32,MPFR_RNDN);
	our_pi_mpreal_32=const_pi(PREC_32,MPFR_RNDN);
	//mpfr_printf("%.*Re\n",64,our_pi_mpreal_32.mpfr_ptr());	
	
	if(thelw_64==1)
	{
		our_pi_mpreal_64.setPrecision(PREC_64,MPFR_RNDN);
		our_pi_mpreal_64=const_pi(PREC_64,MPFR_RNDN);
		//mpfr_printf("%.*Re\n",64,our_pi_mpreal_64.mpfr_ptr());
	
		//cout << "///////////////////////////////////////\n";
		//mpfr_printf("%.*Re\n",64,str_pi_mpreal_32.mpfr_ptr());	
		//mpfr_printf("%.*Re\n",64,str_pi_mpreal_64.mpfr_ptr());
		//getchar();
	}
	
	/////////////////////////////////////////////////////////
	//ORISMOS PRECISION GIA TO PI KAI ORISMOS WINDOW LENGTH//
	/////////////////////////////////////////////////////////	
	//PI_OUR=our_pi_mpreal_32;	//Precision pi
	NN=WINDOW_LENGTH; //To megethos tou parathyrou gia to fft
	
	
	U_16.real.setPrecision(PREC_DOUBLE,MPFR_RNDN);
	U_16.imag.setPrecision(PREC_DOUBLE,MPFR_RNDN);	

	W_16.real.setPrecision(PREC_DOUBLE,MPFR_RNDN);
	W_16.imag.setPrecision(PREC_DOUBLE,MPFR_RNDN);

	T_16.real.setPrecision(PREC_DOUBLE,MPFR_RNDN);
	T_16.imag.setPrecision(PREC_DOUBLE,MPFR_RNDN);

	
	U_32.real.setPrecision(PREC_32,MPFR_RNDN);
	U_32.imag.setPrecision(PREC_32,MPFR_RNDN);	

	W_32.real.setPrecision(PREC_32,MPFR_RNDN);
	W_32.imag.setPrecision(PREC_32,MPFR_RNDN);

	T_32.real.setPrecision(PREC_32,MPFR_RNDN);
	T_32.imag.setPrecision(PREC_32,MPFR_RNDN);

	if(thelw_64==1)
	{
		U_64.real.setPrecision(PREC_64,MPFR_RNDN);
		U_64.imag.setPrecision(PREC_64,MPFR_RNDN);	

		W_64.real.setPrecision(PREC_64,MPFR_RNDN);
		W_64.imag.setPrecision(PREC_64,MPFR_RNDN);

		T_64.real.setPrecision(PREC_64,MPFR_RNDN);
		T_64.imag.setPrecision(PREC_64,MPFR_RNDN);
	}

	for(jj=0;jj<MAX_SIZE;jj++)
	{
		nums_read_16[jj].setPrecision(PREC_DOUBLE,MPFR_RNDN);
		//cout << "precision for " << jj << "element set!\n";
		nums_read_16[jj]=0.0;
		//cout << "Initialization complete!\n" << jj << "\n";

		nums_read_32[jj].setPrecision(PREC_32,MPFR_RNDN);
		//cout << "precision for " << jj << "element set!\n";
		nums_read_32[jj]=0.0;
		//cout << "Initialization complete!\n" << jj << "\n";
		if(thelw_64==1)
		{
			nums_read_64[jj].setPrecision(PREC_64,MPFR_RNDN);
			//cout << "precision for " << jj << "element set!\n";
			nums_read_64[jj]=0.0;
			//cout << "Initialization complete!\n" << jj << "\n";
		}
	}

	//cout << "Teleiwsa me set_precision-initialization!\n";
	
   //Anoigma arxeiou gia diavasma
	ifstream myfile;
	//myfile.open ("Simply_Red_Picture_Book.txt");
	//myfile.open ("Loreena_McKenitt_Marakesh_Night_Market.txt");
	//myfile.open ("Dulce_Pontes_Cancao_Do_Mar.txt");
	//myfile.open ("Guns_N_Roses_Paradise_City.txt");
	//myfile.open ("Donna_Summer_Hot_Stuff.txt");
	//myfile.open ("Emma_Shaplin_Spente_Le_Stelle.txt");	
	myfile.open ("ifft_Spente_Le_Stelle_WL8M_arx500K.txt");
	cerr << "File Opened\n";

	if(myfile.is_open()) //Elegxos gia to an anoikse to arxeio
	{
		deiktis_read=0;
		//while(!myfile.eof());
		while(getline(myfile,string1))
		{ /* while the file still has numbers to be read */
			istringstream costas_read(string1);
			//myfile >> nums_read[deiktis_read];
			costas_read >>	nums_read_16[deiktis_read];
			//mpfr_printf("%.*Re\n",6,nums_read_16[deiktis_read].mpfr_ptr());
			//costas_read >>	nums_read_32[deiktis_read];
			nums_read_32[deiktis_read]=nums_read_16[deiktis_read];
			//mpfr_printf("%.*Re\n",6,nums_read_32[deiktis_read].mpfr_ptr());
			if(thelw_64==1)
			{			
				//costas_read >>	nums_read_64[deiktis_read];
				nums_read_64[deiktis_read]=nums_read_16[deiktis_read];
				//mpfr_printf("%.*Re\n",6,nums_read_64[deiktis_read].mpfr_ptr());
				//cout << "Mexri stigmhs exw diavasei " << deiktis_read << "\n";		
				//getchar();
				//if(costas_read.fail())
				//{ /* something unexpected happened (read error/formatting error) */
				//	cout << "Error while reading...\n";
				//	break;
				//}
			}
			deiktis_read++;
			if(deiktis_read>=MAX_SIZE)
			{
				cerr << "KOUKOUXORITIKOTITA!!!!!\n";
				break;
			}
		} // while(!myfile.eof()) 
		myfile.close();
    } // if(myfile.is_open())
	else
		cout << "Unable to open file, sorry...\n";
	
	cerr << "End of reading!!!\n";

	//Proetoimasia pinaka pou exei diavastei - Orismos tmhmatos gia fft
 	/* Obtain the spectrum - fasma - of wind_sima=new_sima_wave */

			
	jj=0;
	for(ii=arxh_song;ii<(arxh_song+SIZE_TMIMA);ii++)
	{
		F_16[jj].real.setPrecision(PREC_DOUBLE,MPFR_RNDN);
		F_16[jj].imag.setPrecision(PREC_DOUBLE,MPFR_RNDN);

		F_32[jj].real.setPrecision(PREC_32,MPFR_RNDN);
		F_32[jj].imag.setPrecision(PREC_32,MPFR_RNDN);

		if(thelw_64==1)
		{
			F_64[jj].real.setPrecision(PREC_64,MPFR_RNDN);
			F_64[jj].imag.setPrecision(PREC_64,MPFR_RNDN);
		}
		jj++;
	}	


	//F[0].real=0.0;
	//F[0].imag=0.0;
	i=0;
	for(jj=arxh_song;jj<(arxh_song+SIZE_TMIMA);jj++)
	{
		//cerr << i << endl;
		//cerr << jj << endl << endl;
		F_16[i].real = nums_read_16[jj];
		F_16[i].imag = 0.0;

		F_32[i].real = nums_read_32[jj];
		F_32[i].imag = 0.0;

		if(thelw_64==1)
		{
			F_64[i].real = nums_read_64[jj];
			F_64[i].imag = 0.0;
		}
		i++;		
		
	}	
	//cout << "End of segment preparation for fft\n";
	//getchar();

	
	

	/////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	/////////////////////////////// FFT CODE ////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////

	cerr << "FFT Start!\n";

	/*   for(index=NN;index>0;index--)
		{
			F[index].real=F[index-1].real;
			F[index].imag=F[index-1].imag;
		}
 	*/
	/* Below Pointer Operation to avoid Shift  */
	/*F=(F1-1);*/

	 /* THIS HOPEFULLY! CHECKS IF THE WINDOW IS A POWER OF 2 */
	kk=NN;
	n=0;
   	while(kk>1)
	{
		if(kk==1)
		{
			break;
		}
   		if(kk%2!=0)
	  	{
			cout << "Window not power of 2...\n";
	  		return 0;  /* IN THE FUTURE: FILL IN WITH ZEROS!!!! */
	  	}
	 	kk=kk/2;
	 	n++;
	}
	/*    END OF THIS HOPEFULLY! CHECKS IF THE WINDOW IS A POWER OF 2 */
	/*   NN=1;
   	for(p=1;p<=n;p++)
	{
	 	NN=2*NN;
	}
	*/

   	NV2=NN/2;
   	NM1=NN-1;
   	J=1;
   	for(I=1;I<=NM1;I++)
	{
		if (I<J)
	  	{
	   		T_16=F_16[J-1];
	   		F_16[J-1]=F_16[I-1];
	   		F_16[I-1]=T_16;

			T_32=F_32[J-1];
	   		F_32[J-1]=F_32[I-1];
	   		F_32[I-1]=T_32;
			
			if(thelw_64==1)
			{
				T_64=F_64[J-1];
	   			F_64[J-1]=F_64[I-1];
	   			F_64[I-1]=T_64;
			}
	  	}
		K=NV2;
	 	while (K<J)
	  	{			
	   		J=J-K;
	   		K=K/2;
	  	}
		J=J+K;
	} /* FOR  */
   	for(L=1;L<=n;L++)
	{
	 	LE=1;
	 	for(p=1;p<=L;p++)
	  	{
	   		LE=2*LE;
	  	}
	 	LE1=LE/2;

	 	U_16.real=1.0;
	 	U_16.imag=0.0;
	 	W_16.real=cos(our_pi_mpreal_16/LE1);
	 	W_16.imag=(-1)*sin(our_pi_mpreal_16/LE1);

		U_32.real=1.0;
	 	U_32.imag=0.0;
	 	W_32.real=cos(our_pi_mpreal_32/LE1);
	 	W_32.imag=(-1)*sin(our_pi_mpreal_32/LE1);

		if(thelw_64==1)
		{
			U_64.real=1.0;
		 	U_64.imag=0.0;
		 	W_64.real=cos(our_pi_mpreal_64/LE1);
		 	W_64.imag=(-1)*sin(our_pi_mpreal_64/LE1);
		}
	 	
		for(J=1;J<=LE1;J++)
	  	{
	   		for(I=J;I<=NN;I=I+LE)
			{
	 			IP=I+LE1;
				
	 			T_16=complexmul_16(F_16[IP-1],U_16);
	 			F_16[IP-1]=complexsub_16(F_16[I-1],T_16);
	 			F_16[I-1]=complexadd_16(F_16[I-1],T_16);

				T_32=complexmul_32(F_32[IP-1],U_32);
	 			F_32[IP-1]=complexsub_32(F_32[I-1],T_32);
	 			F_32[I-1]=complexadd_32(F_32[I-1],T_32);

				if(thelw_64==1)
				{
					T_64=complexmul_64(F_64[IP-1],U_64);
		 			F_64[IP-1]=complexsub_64(F_64[I-1],T_64);
		 			F_64[I-1]=complexadd_64(F_64[I-1],T_64);
				}
			}

	   		U_16=complexmul_16(U_16,W_16);
			
			U_32=complexmul_32(U_32,W_32);
			
			if(thelw_64==1)
			{
				U_64=complexmul_64(U_64,W_64);
			}		
	  	}
	}


	for(ll=0;ll<LATHOS_BITS;ll++)
	{
 		plithos_errors_bits_real_16_32[ll]=0;
		plithos_errors_bits_imag_16_32[ll]=0;
		
		if(thelw_64==1)
		{
 			plithos_errors_bits_real_32_64[ll]=0;
			plithos_errors_bits_imag_32_64[ll]=0;
		}
	}

	cerr << "FFT End!\n\n";

	cerr << "Error Calculation Start!\n";
	for(ll=0;ll<NN;ll++)
	{
		if(ll%100000L==0)
		{
			fprintf(stderr,"%ld\n",ll);
			fflush(stdout);
		}
		error_16_32_real=0;
		error_16_32_imag=0;
		
		if(thelw_64==1)
		{
			error_32_64_real=0;
			error_32_64_imag=0;
		}

		//mpfr_printf("%.14Re  +  %.14Rei\n",F_16[ll].real.mpfr_ptr(),F_16[ll].imag.mpfr_ptr());
		//mpfr_printf("%.16Re + %.16Rei\n",F_32[ll].real.mpfr_ptr(),F_32[ll].imag.mpfr_ptr());
		//mpfr_printf("%.29Re  +  %.29Rei\n",F_32[ll].real.mpfr_ptr(),F_32[ll].imag.mpfr_ptr());
		//mpfr_printf("%.30Re + %.30Rei\n",F_64[ll].real.mpfr_ptr(),F_64[ll].imag.mpfr_ptr());
		//mpfr_printf("%.62Re + %.62Rei\n",F_64[ll].real.mpfr_ptr(),F_64[ll].imag.mpfr_ptr());

		error_16_32_real=mpfr_find_edd(F_16[ll].real.mpfr_ptr(), F_32[ll].real.mpfr_ptr(), 15);
		if(error_16_32_real<our_10_eis_thn_meion_14)
		{
			error_16_32_real=0;
		}
		plithos_errors_bits_real_16_32[error_16_32_real]++;
		
		error_16_32_imag=mpfr_find_edd(F_16[ll].imag.mpfr_ptr(), F_32[ll].imag.mpfr_ptr(), 15);
		if(error_16_32_imag<our_10_eis_thn_meion_14)
		{
			error_16_32_imag=0;
		}
		plithos_errors_bits_real_32_64[error_16_32_imag]++;

		if(thelw_64==1)
		{
			error_32_64_real=mpfr_find_edd(F_32[ll].real.mpfr_ptr(), F_64[ll].real.mpfr_ptr(), 30);
			if(error_32_64_real<our_10_eis_thn_meion_14)
			{
				error_32_64_real=0;
			}
			plithos_errors_bits_imag_16_32[error_32_64_real]++;

			error_32_64_imag=mpfr_find_edd(F_32[ll].imag.mpfr_ptr(), F_64[ll].imag.mpfr_ptr(), 30);
			if(error_32_64_imag<our_10_eis_thn_meion_14)
			{
				error_32_64_imag=0;
			}
			plithos_errors_bits_imag_32_64[error_32_64_imag]++;
		}

		//mpfr_printf("Error 16-32 real: %d, imaginary: %d\n", error_16_32_real, error_16_32_imag);
		//mpfr_printf("Error 32-64 real: %d, imaginary: %d\n\n", error_32_64_real, error_32_64_imag);
	}
	
	cerr << "Error Calculation End!\n Now Printing...\n";
	
	
	for(ll=0;ll<LATHOS_BITS;ll++)
	{
		printf("16-32 me error %ld digits real: %d (%le%%)\n",ll,plithos_errors_bits_real_16_32[ll],plithos_errors_bits_real_16_32[ll]/(deiktis_read*1.0)); 
		printf("16-32 me error %ld digits imag: %d (%le%%)\n",ll,plithos_errors_bits_imag_16_32[ll],plithos_errors_bits_imag_16_32[ll]/(deiktis_read*1.0)); 

		if(thelw_64==1)
		{
			printf("32-64 me error %ld digits real: %d",ll,plithos_errors_bits_real_32_64[ll]); 
			printf("(%le%%)\n",plithos_errors_bits_real_32_64[ll]/(deiktis_read*1.0));			
			printf("32-64 me error %ld digits imag: %d",ll,plithos_errors_bits_imag_32_64[ll]); 
			printf("(%le%%)\n",plithos_errors_bits_imag_32_64[ll]/(deiktis_read*1.0));
		}
	}

	mpfr_printf("Plithos arithmwn: %ld\n",deiktis_read);
	
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

	/////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	/////////////////////// END OF FFT CODE /////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////

	return 0;

} //end main


