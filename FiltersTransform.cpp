/*******************************************************************************
*       File Name: ProcessImage.cpp 
*		Purpose  : Program for performing follwing SPATIAL TRANSFORMATION
*					on an UNCOMPRESSED  Binary, Gray-Level, or Colored  Image
*					1. Notch Filter 
*					2. Lowpass Filter: Ideal,Butterworth and Gaussian
*					3. Highpass filter: Ideal,Butterworth and Gaussian
*					4. Laplacian Filter
*					5. Homomorphic filters
*		Input	 : Name of the INPUT and OUTPUT TIFF file and choice[0-9]
*		Output	 : Processed image stored in the OUTPUT TIFF File Format
*		Test File: lena50.tiff
*_______________________________________________________________________________
*	Reference : 1) TIFF, Revision 6.0, (June 1992). Adobe Developers Association.
*	               www.adobe.com/Support/TechNotes.html and
*                  ftp://ftp.adobe.com/pub/adobe/DeveloperSupport/TechNotes/PDFfiles
*               2) Digital Image Processing by Rafael C. Ganzalez and
*                  Richard E. Woods, Pearson Education, India, (for image
*                  processing operations)
*_______________________________________________________________________________
*		By	  : Maitrey G. Ranade
*
***********************************************************************************/

// Preprocessor Statements
#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>

// Some Constants
#define TRUE 1				// Boolean True
#define FALSE 0				// Boolean False
#define L 256				// Color Levels
#define MAX 400				// Maximum Dimension of Image
#define WSIZE 3				// Window Size

// legacy damage control
#define long int


// Global Shared variables to be used among all functions
unsigned long ImageWidthByte, ImageLength,ImageWidth, ImageSize;// Specify the dimensions of image
unsigned long file_position;									// set by read and used by write function
void *ImageBuffer,*OutBuffer;									// Stores the input and output raster image respectively
unsigned char *IB,*OB;											// Image processing pointer for In and out buffer
unsigned char pic_byte;											// Stores the picture byte for intermediate processing
double pi;



// Function Declarations
int screen();
void error_message(const char *error);
void read_tiff(char in_filename[14]);	// read to IB
void copy_tiff(char in_filename[14], char out_filename[14]);
void write_tiff(char out_filename[14]);	// Write from OB

void normalize(int M, int N, int out[MAX][MAX]);
void copy_image(int M, int N, double a[MAX][MAX],int out[MAX][MAX]);
void fourier(int M, int N, int fxy[MAX][MAX],double fuvreal[MAX][MAX],double fuvimg[MAX][MAX]);
void inverse_fourier(int M, int N, double fuvreal[MAX][MAX],double fuvimg[MAX][MAX],double gxyreal[MAX][MAX],double gxyimg[MAX][MAX]);
void transformcentre(int M, int N, int fxy[MAX][MAX], int fdashxy[MAX][MAX]);
void retransformcentre(int M, int N, double fdashxy[MAX][MAX], double fxy[MAX][MAX]);
void apply_filter(int M, int N, double H[MAX][MAX], double fuvreal[MAX][MAX], double fuvimg[MAX][MAX],double Guvreal[MAX][MAX],double Guvimg[MAX][MAX]);
void convolution(int M, int N, double f[MAX][MAX], double h[MAX][MAX], double fh[MAX][MAX]);
void convolution1(int M, int N, double f[MAX][MAX], double h[MAX][MAX], double fh[MAX][MAX]);
void MSRE(int M, int N, int fxy[MAX][MAX], int fxycouter[MAX][MAX]);

int main()
{

	// variable declaration

	char in_filename[14], out_filename[14];							// File names
	

	// For Fourier Transformations
	unsigned long x,y,u,v,M,N ;	// M Length N Width
	int fxy[MAX][MAX],out[MAX][MAX];					// input, output image 
	int fdashxy[MAX][MAX];								// preprocessing 
	int d,n=2;											// raduis							
	double fuvreal[MAX][MAX],fuvimg[MAX][MAX];			//output of fourier transform
	double Guvreal[MAX][MAX],Guvimg[MAX][MAX];			//after applying transform
	double gxyreal[MAX][MAX], gxyimg[MAX][MAX];			// reconstructed image after inverse fourier
	double gdashxy[MAX][MAX];							// post processing
	double D0,D0sq,Duv,Duvsq,temp,M_2,N_2,U,V;
	double YL,YH,c;										//constants for Homo Morphic filtering
	double H[MAX][MAX],Hreal[MAX][MAX],Himg[MAX][MAX];					// Filter
	double f[MAX][MAX],h[MAX][MAX],fh[MAX][MAX];		// Convolution
	int choice, FOURIER = FALSE,type;

	pi= acos(-1.0); 
	
	// Default intialization
	/*strcpy(in_filename,"image1.tif");
	strcpy(out_filename,"out.tif");
	read_tiff(in_filename);
	copy_tiff(in_filename,out_filename);
	
	// Intializing the output buffer for processing		
	OutBuffer =  malloc(ImageSize);
	if(OutBuffer == NULL)
	{
		printf("\n\nMemory Required: %d Bytes", ImageSize);
		error_message(" Memory Allocation Error- Insufficient Memory. ");
	}
	IB = (unsigned char *)ImageBuffer;
	OB = (unsigned char *)OutBuffer;		

	// copy image to array for processing
	M = ImageLength; 
	N = ImageWidth;
	for(x=0; x<M; x++)
	{
		for(y=0; y<N; y++)
		{
			fxy[x][y] = (*IB);
			IB++;

		}
	}*/



	// iterative choice taking loop
	do
	{
	   choice = screen();
	   printf("\n");
	   switch(choice)
	   {
		case  1 : // Fourier Transform
				// step 1: pre processing
				transformcentre(M,N,fxy,fdashxy);

				fourier(M, N, fdashxy,fuvreal,fuvimg);
				//copy_image(M, N, fuvreal,out);
				for(x=0; x<M; x++)
				{
					for(y=0; y<N; y++)
					{
				
						temp = sqrt(fuvreal[x][y]*fuvreal[x][y]+fuvimg[x][y]*fuvimg[x][y]); 
						out[x][y] = (int) (100*(log(1.0+temp)+0.5)); // round off

					}
				}

				normalize(M, N, out);
				FOURIER = TRUE;
				break;


		case  2 : // Inverse Fourier Transform
				if(FOURIER == TRUE)
				{
					inverse_fourier(M, N, fuvreal, fuvimg, gxyreal, gxyimg);
	
					//step5: post processing
					retransformcentre(M,N,gxyreal,gdashxy);	

					copy_image(M, N, gdashxy,out);
					normalize(M, N, out);
					MSRE(M,N,fxy,out);
				}
				else
				{
					printf("\nThe inverse Fourier not possible as the Fourier Transformation results not available");
				}

				break;


		case  3 : // Notch Filter
				//define filter
				for(u=0;u<M;u++)
					for(v=0;v<N;v++)
						H[u][v] = 1.0;
				H[M/2][N/2] = 0.0;
				
				
				// step 1: pre processing
				transformcentre(M,N,fxy,fdashxy);

				//step 2: applay fourier
				fourier(M, N, fdashxy,fuvreal, fuvimg);

				//step 3: apply filter
				apply_filter(M, N, H,fuvreal,fuvimg,Guvreal,Guvimg);

				// step 4: take invere transorm
				inverse_fourier(M, N, Guvreal, Guvimg, gxyreal, gxyimg);
				
				//step5: post processing
				retransformcentre(M,N,gxyreal,gdashxy);
				copy_image(M, N, gdashxy,out);
				normalize(M, N, out);
				MSRE(M,N,fxy,out);
				break;


		case  4 : // Low Pass Filters
				//define filter
				printf("\nEnter choice for Low Pass Filter");
				do{
					 printf("\nPress [1] for Ideal LPF [2] for Butterworth LPF [3] for Gaussian LPF : ");
					 scanf("%d",&type);
					}while((type<1)||(type>3));
				
				printf("\nEnter the value of D0: ");
				scanf("%f",&d);
				D0 = (double) d;


				M_2 = ((double) M)/2.0;
				N_2 = ((double) N)/2.0;

				if(type==1)	// Ideal low pas filter
				{
					for(u=0;u<M;u++)
						for(v=0;v<N;v++)
						{
							U = (double) u;
							V = (double) v;
							Duv = sqrt((U-M_2)*(U-M_2)+ (V-N_2)*(V-N_2));
							if(Duv<=D0)
								H[u][v] = 1.0;
							else
								H[u][v] = 0.0;
						}
				}
				if (type==2) // Butter Worth low pass filter
				{
					n = 2;// FOR PRACTICAL APPLICATIONS N IS TAKEN 2
					for(u=0;u<M;u++)
						for(v=0;v<N;v++)
						{
							U = (double) u;
							V = (double) v;
							Duv = sqrt((U-M_2)*(U-M_2)+ (V-N_2)*(V-N_2));
							
							H[u][v] = 1.0/(1.0+pow((Duv/D0),2*n));

						}

				}
				if (type==3) // Gaussian low pass filter
				{
					for(u=0;u<M;u++)
						for(v=0;v<N;v++)
						{
							U = (double) u;
							V = (double) v;
							Duv = sqrt((U-M_2)*(U-M_2)+ (V-N_2)*(V-N_2));
		
							H[u][v] = exp(-(Duv*Duv)/(2.0*D0*D0));

						}

				}
					
				
				// step 1: pre processing
				transformcentre(M,N,fxy,fdashxy);

				//step 2: applay fourier
				
				fourier(M, N, fdashxy,fuvreal, fuvimg);

				//step 3: apply filter
				apply_filter(M, N, H,fuvreal,fuvimg,Guvreal,Guvimg);

				// step 4: take invere transorm
				inverse_fourier(M, N, Guvreal, Guvimg, gxyreal, gxyimg);
				
				//step5: post processing
				retransformcentre(M,N,gxyreal,gdashxy);
				copy_image(M, N, gdashxy,out);
				normalize(M, N, out);
				MSRE(M,N,fxy,out);
				break;


		case  5 : // High Pass Filter
				printf("\nEnter choice for High Pass Filter");
				do{
					 printf("\nPress [1] for Ideal HPF [2] for Butterworth HPF [3] for Gaussian HPF : ");
					 scanf("%d",&type);
					}while((type<1)||(type>3));
				
				printf("\nEnter the value of D0: ");
				scanf("%f",&d);
				D0 = (double) d;


				M_2 = ((double) M)/2.0;
				N_2 = ((double) N)/2.0;

				if(type==1)	// Ideal High pass filter
				{
					for(u=0;u<M;u++)
						for(v=0;v<N;v++)
						{
							U = (double) u;
							V = (double) v;
							Duv = sqrt((U-M_2)*(U-M_2)+ (V-N_2)*(V-N_2));
							if(Duv<=D0)
								H[u][v] = 0.0;
							else
								H[u][v] = 1.0;
						}
				}
				if (type==2) // Butter Worth High pass filter
				{
					n = 2;// FOR PRACTICAL APPLICATIONS N IS TAKEN 2
					for(u=0;u<M;u++)
						for(v=0;v<N;v++)
						{
							U = (double) u;
							V = (double) v;
							Duv = sqrt((U-M_2)*(U-M_2)+ (V-N_2)*(V-N_2));
							
							H[u][v] = 1.0/(1.0+pow((D0/Duv),2*n));

						}

				}
				if (type==3) // Gaussian High pass filter
				{
					for(u=0;u<M;u++)
						for(v=0;v<N;v++)
						{
							U = (double) u;
							V = (double) v;
							Duv = sqrt((U-M_2)*(U-M_2)+ (V-N_2)*(V-N_2));
		
							H[u][v] = 1.0-exp(-(Duv*Duv)/(2.0*D0*D0));

						}

				}
					
				
				// step 1: pre processing
				transformcentre(M,N,fxy,fdashxy);

				//step 2: applay fourier
				
				fourier(M, N, fdashxy,fuvreal, fuvimg);

				//step 3: apply filter
				apply_filter(M, N, H,fuvreal,fuvimg,Guvreal,Guvimg);

				// step 4: take invere transorm
				inverse_fourier(M, N, Guvreal, Guvimg, gxyreal, gxyimg);
				
				//step5: post processing
				retransformcentre(M,N,gxyreal,gdashxy);
				copy_image(M, N, gdashxy,out);
				normalize(M, N, out);
				MSRE(M,N,fxy,out);
				break;


		case  6 : // Laplacian Filter
				//define filter
				for(u=0;u<M;u++)
					for(v=0;v<N;v++)
					{
						U = (double) u;
						V = (double) v;
						Duvsq = (U-M_2)*(U-M_2)+ (V-N_2)*(V-N_2);

						H[u][v] = -4.0*pi*pi*Duvsq;
					}
			
				
				
				// step 1: pre processing
				transformcentre(M,N,fxy,fdashxy);

				//step 2: applay fourier
				fourier(M, N, fdashxy,fuvreal, fuvimg);

				//step 3: apply filter
				apply_filter(M, N, H,fuvreal,fuvimg,Guvreal,Guvimg);

				// step 4: take invere transorm
				inverse_fourier(M, N, Guvreal, Guvimg, gxyreal, gxyimg);
				
				//step5: post processing
				retransformcentre(M,N,gxyreal,gdashxy);

				//obtain enhanced image
				copy_image(M, N, gdashxy,out);
				normalize(M, N, out);				
	
				for(x=0; x<M; x++)
					for(y=0; y<N; y++)				
						out[x][y] = fxy[x][y] - out[x][y];


				MSRE(M,N,fxy,out);
				break;


		case  7 : // Homo-morphic Filters
				//define filter
				printf("\nEnter the value of D0: ");
				scanf("%f",&d);
				D0 = (double) d;
				D0sq = d*d;

				YL= 0.01;YH=1.1;c=2.0;

				for(u=0;u<M;u++)
					for(v=0;v<N;v++)
					{
						U = (double) u;
						V = (double) v;
						Duvsq = (U-M_2)*(U-M_2)+ (V-N_2)*(V-N_2);

						H[u][v] = (YH-YL)*(1.0-exp(c*Duvsq/D0sq))+YL;

					}
			
				
				
				// step 1: pre processing
				transformcentre(M,N,fxy,fdashxy);

				//step 2: applay fourier
				fourier(M, N, fdashxy,fuvreal, fuvimg);

				//step 3: apply filter
				apply_filter(M, N, H,fuvreal,fuvimg,Guvreal,Guvimg);

				// step 4: take invere transorm
				inverse_fourier(M, N, Guvreal, Guvimg, gxyreal, gxyimg);
				
				//step5: post processing
				retransformcentre(M,N,gxyreal,gdashxy);

				//obtain enhanced image
				copy_image(M, N, gdashxy,out);
				normalize(M, N, out);				
				for(x=0; x<M; x++)
					for(y=0; y<N; y++)				
						out[x][y] = fxy[x][y] - out[x][y];

				MSRE(M,N,fxy,out);

				break;

		case  8 : // Testify Convolution operation
				
				//define signal h 
				printf("\nEnter choice for Low Pass Filter");
				do{
					 printf("\nPress [1] for Ideal LPF [2] for Butterworth LPF [3] for Gaussian LPF : ");
					 scanf("%d",&type);
					}while((type<1)||(type>3));
				
				printf("\nEnter the value of D0: ");
				scanf("%f",&d);
				D0 = (double) d;


				M_2 = ((double) M)/2.0;
				N_2 = ((double) N)/2.0;

				if(type==1)	// Ideal low pas filter
				{
					for(u=0;u<M;u++)
						for(v=0;v<N;v++)
						{
							U = (double) u;
							V = (double) v;
							
							Duv = sqrt((U-M_2)*(U-M_2)+ (V-N_2)*(V-N_2));
							if(Duv<=D0)
								H[u][v] = 1.0;
							else
								H[u][v] = 0.0;
						}
				}
				if (type==2) // Butter Worth low pass filter
				{
					n = 2;// FOR PRACTICAL APPLICATIONS N IS TAKEN 2
					for(u=0;u<M;u++)
						for(v=0;v<N;v++)
						{
							U = (double) u;
							V = (double) v;
							Duv = sqrt((U-M_2)*(U-M_2)+ (V-N_2)*(V-N_2));
							
							H[u][v] = 1.0/(1.0+pow((Duv/D0),2*n));

						}

				}
				if (type==3) // Gaussian low pass filter
				{
					for(u=0;u<M;u++)
						for(v=0;v<N;v++)
						{
							U = (double) u;
							V = (double) v;
							Duv = sqrt((U-M_2)*(U-M_2)+ (V-N_2)*(V-N_2));
		
							H[u][v] = exp(-(Duv*Duv)/(2.0*D0*D0));

						}

				}
					
			

				for(x=0;x<M;x++)
				{
					for(y=0;y<N;y++)
					{

						Himg[x][y] = 0.0;
						Hreal[x][y] = H[x][y];
					}
				}


				// step 1: take invere transorm of Filter to get its spatial domain
				inverse_fourier(M, N, Hreal, Himg, gxyreal, gxyimg);


			for(x=0;x<M;x++)
				{
					for(y=0;y<N;y++)
					{
						f[x][y] = (double) fxy[x][y];
						h[x][y] =  gxyreal[x][y];
						fh[x][y] = 0.0;

					}

				}

				//step 2: Convolute
				convolution(M,N,f,h,fh);

				copy_image(M, N, fh,out);
				normalize(M, N, out);				
				break;


		case  9 : // Getting Input File
				FOURIER = FALSE;
				free(ImageBuffer);
				free(OutBuffer);
		
				printf("Enter the name of Input TIFF File (*.tif): ");
				scanf( "%s", in_filename);
				read_tiff(in_filename);

				printf("\nEnter the name of Output TIFF File (*.tif): ");
				scanf( "%s", out_filename);
				copy_tiff(in_filename,out_filename);
		
				// Intializing the output buffer for processing		
				OutBuffer =  malloc(ImageSize);
				if(OutBuffer == NULL)
				{
					printf("\n\nMemory Required: %d Bytes", ImageSize);
					error_message(" Memory Allocation Error- Insufficient Memory. ");
				}
				IB = (unsigned char *)ImageBuffer;
				OB = (unsigned char *)OutBuffer;		
	
				// copy image to array for processing
				M = ImageLength; 
				N = ImageWidth;
				for(x=0; x<M; x++)
				{
					for(y=0; y<N; y++)
					{
						fxy[x][y] = (*IB);
						IB++;
					}
				}

				break;

		case  0 : //  Free Image Buffer and exit
					printf("\n\n Exiting..... \n");
					exit(0);
		default : printf("\nWrong choice  ");
			  break;


	    }//end switch

	if(choice != 9)
	{
		printf("\nWriting to output File...");
		OB = (unsigned char *)OutBuffer;	
		
		// save to out bufer
		for(x=0; x<M; x++)
		{
			for(y=0; y<N; y++)
			{
				(*OB) = out[x][y];
				OB++;
			}
		}


		write_tiff(out_filename);
		printf("\n\nProcessing SUCCESSFUL.....results in file: %s \n",out_filename);
	}


  }while(choice != 0);
  free(OutBuffer);
}//main program ends here
/*********************************************************************************************************
	Transform Centre multiply image by (-1)**(x+y) for preprocessing
**********************************************************************************************************/
void transformcentre(int M, int N, int fxy[MAX][MAX], int fdashxy[MAX][MAX])
{
	int x,y,sign;

	for(x=0; x<M; x++)
	{
		for(y=0; y<N; y++)
		{
			if((x+y)%2 == 0)
				sign = 1;
			else
				sign = -1;
			fdashxy[x][y] = sign*fxy[x][y];
		}
	}
}


/*********************************************************************************************************
	RE Transform Centre multiply image by (-1)**(x+y) for post processsing
**********************************************************************************************************/
void retransformcentre(int M, int N, double fdashxy[MAX][MAX], double fxy[MAX][MAX])
{
	int x,y,sign;

	for(x=0; x<M; x++)
	{
		for(y=0; y<N; y++)
		{
			if((x+y)%2 == 0)
				sign = 1;
			else
				sign = -1;
			fxy[x][y] = sign*fdashxy[x][y];
		}
	}
}
/*********************************************************************************************************
	Apply filter H to Fuv to obtain Guv
**********************************************************************************************************/
void apply_filter(int M, int N, double H[MAX][MAX], double fuvreal[MAX][MAX], double fuvimg[MAX][MAX],double Guvreal[MAX][MAX],double Guvimg[MAX][MAX])
{
	int u,v;

	printf("\nApplying Filter..........\n");
	for(u=0;u<M;u++)
		for(v=0;v<N;v++)
		{
			Guvreal[u][v] = H[u][v]*fuvreal[u][v];
			Guvimg[u][v]  = H[u][v]*fuvimg[u][v];
		}
}

					
/*********************************************************************************************************
	Normalize the values
**********************************************************************************************************/
void normalize(int M, int N, int out[MAX][MAX])
{
  int i,j;

  for(i=0;i<M;i++)
  {
    for(j=0;j<N;j++)
    {
      //  normalizing values
      if(out[i][j]<0)
			out[i][j] = 0;
      if(out[i][j]>(L-1))
			out[i][j] = L-1;
      }//end j
    }// end i
}

/*********************************************************************************************************
	Copy double array to int for image display
**********************************************************************************************************/
void copy_image(int M, int N, double a[MAX][MAX],int out[MAX][MAX])
{
	int i,j;

	for(i=0;i<M;i++)
		for(j=0;j<N;j++)
			out[i][j] = (int) (a[i][j]+0.5);		// round off
}

/*********************************************************************************************************
	Find the fourier transform
**********************************************************************************************************/
void fourier(int M, int N, int fxy[MAX][MAX],double fuvreal[MAX][MAX],double fuvimg[MAX][MAX])
{
	double theta,sumreal,sumimg;
	int u,v,x,y;
	
	//Apply Fourier Function
	printf("\nProcessing for Fourier Transform Please wait........");
	for(u=0; u<M; u++)
	{
		for(v=0; v<N; v++)
		{
			sumreal=0.0;
			sumimg=0.0;
			for(x=0;x<M;x++)
			{
				for(y=0;y<N;y++)
				{
					theta = 2.0 *pi*((double(u*x)/double (M))+(double (v*y)/double(N)));
					sumreal = sumreal+fxy[x][y]*cos(theta);
					sumimg  = sumimg +fxy[x][y]*sin(theta);	
				}
			}
			fuvreal[u][v]=sumreal/double(M*N);
			fuvimg[u][v]=-sumimg/double(M*N);
		}
	}	
	printf("\nFourier Transformation done successfully\n");
}


/*********************************************************************************************************
	Find the inverse fourier transform
**********************************************************************************************************/
void inverse_fourier(int M, int N, double fuvreal[MAX][MAX],double fuvimg[MAX][MAX],double gxyreal[MAX][MAX],double gxyimg[MAX][MAX])
						
{
	double theta,sumreal,sumimg;
	int u,v,x,y;

	//Apply Inverse Fourier Function
	printf("\nProcessing for Inverse Fourier Transform Please wait........");

	for(x=0; x<M; x++)
	{
		for(y=0; y<N; y++)
		{
			sumreal=0.0;
			sumimg=0.0;
			for(u=0;u<M;u++)
			{
				for(v=0;v<N;v++)
				{
					theta = 2.0 *pi*((double(u*x)/double (M))+(double (v*y)/double(N)));
					sumreal = sumreal+(double)fuvreal[u][v]*(double)cos(theta)-(double)fuvimg[u][v]*(double)sin(theta);
					sumimg  = sumimg +(double)fuvreal[u][v]*(double)sin(theta)+(double)fuvimg[u][v]*(double)cos(theta);	
				}
			}
			gxyreal[x][y]=sumreal;
			gxyimg[x][y] = sumimg;
		}
	}	
	printf("\nInverse Fourier Transformation done successfully\n");
}
	
/*********************************************************************************************************
	Perform convolution
**********************************************************************************************************/
void convolution(int M, int N, double f[MAX][MAX], double h[MAX][MAX], double fh[MAX][MAX])
{
	int x,y,s,t,xms,ymt;
	double sum,hxy,fst;
	


	//Apply Convolution
	printf("\nPerforming convolution........");

	for(x=0; x<M; x++)
	{
		for(y=0; y<N; y++)
		{
			sum = 0.0;
			for(s=-M/2+x;s<M/2+x;s++)
			{
				for(t=-N/2+y;t<N/2+y;t++)
				{
					xms = x-s;
					ymt = y-t;
			
					if((xms>=0)&&(xms<M)&&(ymt>=0)&&(ymt<N))
						hxy = h[xms][ymt];
					else
						hxy = 0.0;
					
					if((s>=0)&&(s<M)&&(t>=0)&&(t<N))
						fst = f[s][t];
					else
						fst = 0.0;					
					
					sum = sum+fst*hxy;

				}
			}
			fh[x][y] = sum;
		}
	}

	printf("\nConvolution done successfully\n");
}

/*********************************************************************************************************
	Perform convolution C++ code from internet
**********************************************************************************************************/
void convolution1(int M, int N, double f[MAX][MAX], double h[MAX][MAX], double fh[MAX][MAX])
{
	int x,y,s,t,ii,jj,ss,tt;
	double sum;
	int KcenterX = M/2;
	int KcenterY = N/2;


	//Apply Convolution
	printf("\nPerforming convolution........");

	for(x=0; x<M; x++)
	{
		for(y=0; y<N; y++)
		{
			sum = 0.0;
			for(s=0;s<M;s++)
			{
				ss = M-1-s;
				for(t=0;t<N;t++)
				{
					tt = N-1-t;
					ii = x+s-KcenterX;
					jj = y+t-KcenterY;
			
					if((ii>=0)&&(ii<M)&&(jj>=0)&&(jj<N))
						sum = sum+f[ii][jj]*h[ss][tt];

				}
			}
			fh[x][y] = sum;
		}
	}

	printf("\nConvolution done successfully\n");


}
/*********************************************************************************************************
	Calculation of Mean Square Reconstruction Error(MSRE)
**********************************************************************************************************/


void MSRE(int M, int N, int fxy[MAX][MAX], int fxycouter[MAX][MAX])
{
	int x, y, npixels=0;
	double sum1, sum2;
	double fxyij, fxycij;
	double msre=0.0;

	
	sum1=0.0; 
	sum2=0.0;
	// Initialize the value of ndiff
	for(x=0; x<M; x++)
		for(y=0; y<N; y++)
		{
			fxyij =(double) fxy[x][y];
			fxycij=(double) fxycouter[x][y];
			npixels++;
			sum1=sum1+fxyij*fxyij;
			sum2=sum2+(fxyij-fxycij)*(fxyij-fxycij);
		}  // y-loop
	if(npixels==0) 
		error_message("Total no. of pixels in MSRE calculation is zero");
	if(sum2==0)
		printf("\nOutput image is the same as the input image\n");

	if(sum1==0.0) 
		error_message("\nOriginal image is blank");
	else 
		msre=sum2/sum1;

	printf("\nMean Square Reconstruction Error(eps) : %10.5f\n", msre);


}



/*********************************************************************************************************
	Prints various choices
**********************************************************************************************************/
int screen()
{
   int ch;

   do
   {
       printf("\n\t��������������������������������������������ͻ");
       printf("\n\t�                MAIN  MENU                  �");
       printf("\n\t�      1.    Fourier Transform               �");
       printf("\n\t�      2.    Inverse Fourier Transform       �");
       printf("\n\t�      3.    Notch Filter                    �");
       printf("\n\t�      4.    Low Pass Filter                 �");
       printf("\n\t�      5.    High Pass Filter                �");
       printf("\n\t�      6.    Laplacican Filter               �");
       printf("\n\t�      7.    Homomorphic Filter              �");
       printf("\n\t�      8.    Convolution Theorem             �");
       printf("\n\t�      9.    Set Input/output Files          �");
       printf("\n\t�      0.    EXIT.                           �");
       printf("\n\t��������������������������������������������ͼ");
       printf("\n\t               Option.....");
       scanf("%d",&ch);
   }
   while(ch<0 || ch>9);
   return(ch);
}


/*********************************************************************************************************
	Read tiff file to an image buffer
**********************************************************************************************************/
void read_tiff(char in_filename[14])

{

// Variable Declarations for reading from File
	unsigned long  NoOfStrips, StripOffsets, RowsPerStrip,
				  StripsPerImage,   StripByteCounts;

	short int BitWhite, BitsPerSample;

	short unsigned int i, j;
	short unsigned int byte_order, TIFF_identifier;
	unsigned short int no_of_entries, tag, field_type;
	unsigned long offset_address, count, value, so, sbc, StripAddress,
		 StripSize;
//	unsigned char *IB;
	FILE *infile_ptr;



	// Open input file
	if((infile_ptr = fopen((const char *)in_filename, "rb+")) ==NULL)
		error_message("Input File Open Failure ");


	// Read and validate Byte Order Value
	fread((void *)&byte_order, (size_t)2, (size_t)1, infile_ptr);
	if(byte_order != 0x4949)
		error_message(" Invalid Value of Byte Order.");


	// Read and validate Tiff Indentifier
	fread((void *)&TIFF_identifier, (size_t)2, (size_t)1, infile_ptr);
	if(TIFF_identifier != 0x002A)
		error_message(" Invalid Tiff File Format.");


	// Read offset address and reposition file pointer at offset address
	fread((void *)&offset_address, (size_t)4, (size_t)1, infile_ptr);
	fseek(infile_ptr, offset_address, SEEK_SET);


	//Read the total number of entries and their contents
	fread((void *)&no_of_entries, (size_t)2, (size_t)1, infile_ptr);

	for(i=0; i<no_of_entries; i++)
	{
      fread((void *)&tag, (size_t)2, (size_t)1, infile_ptr);
      fread((void *)&field_type, (size_t)2, (size_t)1, infile_ptr);
	  fread((void *)&count, (size_t)4, (size_t)1, infile_ptr);
      fread((void *)&value, (size_t)4, (size_t)1, infile_ptr);


	  // Processing different type of relevant entries
      switch(tag)
      {
		case 256 :	// Number of pixels per row
					ImageWidth = value; 
					//printf("\nImage Width = %d",value);
					break;
		
		case 257 :	// Number of Pixel Rows
					ImageLength = value;
					//printf("\nImage Length = %d",value);
					break;
		
		case 258 :	// Number of bit for each color
					BitsPerSample = (short)value;
					//printf("\nBits per sample = %d",value);
					break;

		case 259 :  // Check if image is compressed
					if(value !=1)
						error_message(" Image is Compressed.");
					break;

		case 262 :	// Photochrometic Interpretation (for us it should be 0-3)
					if(value<0 || value >3)
						error_message("Image Type Cann't be processed by the demo....");

					//printf("\nPhotochrometic Level: %d", value);
					if(value ==1)	// Black and White image
						BitWhite =1;
					else 
						BitWhite =0;

					break;

		case 273 :	// Number of strips and their relative address w.r.t Start

					NoOfStrips = count;
					StripOffsets= value;
					//printf("\nNumber of strips = %d",count);
					break;

		case 278 :	//Number of rows per strip.
					RowsPerStrip = value;
					//printf("\nNumber of rows in a strip = %d",value);
					StripsPerImage = (ImageLength + RowsPerStrip -1)/RowsPerStrip;
					break;

		case 279 :	//Number of Bytes in a strip
					StripByteCounts = value;
					//printf("\nNumber of bytes in a strip = %d",value);
					break;
 
		default :   // Ignore other Tags
					break;

      }// End Switch
   }// End For


	// Calculation of Image size in bytes
	
	switch(BitsPerSample)
    {
		case 1 :	ImageWidthByte = (ImageWidth + 7)/8; 
					break;
		case 8 :	ImageWidthByte = ImageWidth;
					break;
		case 24:	ImageWidthByte = 3*ImageWidth;
					break;
		default:	break;
	}



	ImageSize = ImageWidthByte*ImageLength;


	//  Bringing Image in buffer for processing
	ImageBuffer =  malloc(ImageSize);
	if(ImageBuffer == NULL)
	{
		printf("\n\nMemory Required: %d Bytes", ImageSize);
		error_message(" Memory Allocation Error- Insufficient Memory. ");
	}


	// Intialising the buffer
	IB = (unsigned char *)ImageBuffer;
	for(i=0; i<ImageLength; i++)
		for(j=0; j<ImageWidthByte; j++)
		{
			(*IB) =0;
			IB++;
		}

	//  Read the image from strips in the TIFF format and copy to buffer

	so = StripOffsets;
	sbc = StripByteCounts;
	IB = (unsigned char *)ImageBuffer;

	// Repositioning input file pointer to Strip Offset
	fseek(infile_ptr, so, SEEK_SET);

	if(NoOfStrips ==1)
	{
		file_position = (unsigned long)ftell(infile_ptr);
		fread(ImageBuffer, (size_t)1, (size_t)ImageSize, infile_ptr);
	}

   else
		for(i=0; i<NoOfStrips; i++)
		{
			fread(&StripAddress, 4, 1, infile_ptr);
			so = so + 4;
			fseek(infile_ptr, sbc, SEEK_SET);
			sbc = sbc +4;
			fread(&StripSize, 4, 1, infile_ptr);
			fseek(infile_ptr, StripAddress, SEEK_SET);
			if (i==0)
					file_position = (unsigned long)ftell(infile_ptr);
			fread(IB, (size_t)StripSize, (size_t)1, infile_ptr);
			IB = (unsigned char *)IB + StripSize;
			fseek(infile_ptr, so, SEEK_SET);

		}
		fclose(infile_ptr);

}



/*******************************************************************************************************************
	Copy the input file to out put file
*******************************************************************************************************************/
void copy_tiff(char in_filename[14], char out_filename[14])
{

	// Variables for file pointers
	FILE *infile_ptr, *outfile_ptr;
	unsigned char data;


	// Open Input File

	if((infile_ptr = fopen((const char *)in_filename, "rb+")) ==NULL)
		error_message("Input File Open Failure ");

	// open Output File

	if((outfile_ptr = fopen((const char *)out_filename, "wb")) ==NULL)
	  error_message(" Output File Creation Failure.");



	// Copy input to Output
	while( fread((void *)&data, (size_t)1, (size_t)1, infile_ptr) !=0)
	  fwrite((void *)&data, (size_t)1, (size_t)1, outfile_ptr);



	// Closing the file pointer
	fclose(infile_ptr);
	fclose(outfile_ptr);
}



/****************************************************************************************************************
	Write an image to ouput tiff file to position specified by write_image_position
******************************************************************************************************************/

void write_tiff(char out_filename[14])

{

	// Variable Declaration

	FILE *outfile_ptr;
	unsigned char pic_byte;
	short unsigned int i, j;						// loop control variables
	unsigned long file_position_out = file_position;

	// Reopen the output file in read-write mode

	if((outfile_ptr = fopen((const char *)out_filename, "rb+")) ==NULL)
	error_message(" Output File reopening Failed.");


	//positioning the outfile pointer
	fseek(outfile_ptr, file_position_out, SEEK_SET);


	// Saving the Raster Image from outbuffer to Output File
	OB = (unsigned char *)OutBuffer;

	for(i=0; i<ImageLength; i++)
	{
		for(j=0; j<ImageWidthByte; j++)
		{
			pic_byte = (*OB);
			fwrite((const void *) &pic_byte, (size_t)sizeof(pic_byte), (size_t)1, outfile_ptr);
			file_position_out = (unsigned long)ftell(outfile_ptr);
			OB++;

		}
	}

	// Closing File pointers

	fclose(outfile_ptr);

}

/********************************************************************************************************************
	Function to exit grcefully when error is encountered
**********************************************************************************************************************/

void error_message(const char *error)
{
	printf("\n\n\n\aError: %s. Cann't Continue.....",error);
	printf("\nUNSUCCESFUL Termintation of Processing ....\n");
	//getch();
	//exit(0);
}
