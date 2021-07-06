/*******************************************************************************
*       File Name: ProcessImage.cpp 
*		Purpose  : Program for performing follwing SPATIAL TRANSFORMATION
*					on an UNCOMPRESSED  Binary, Gray-Level, or Colored  Image
*					1. Negation 
*					2. Logrithmic Transformation
*					3. Power Law Transformation
*					4. Gray Level Slicing
*					5. Histogram Equalisation
*					6. Histogram Matching
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

// Some Constants
#define TRUE 1				// Boolean True
#define FALSE 0				// Boolean False
#define L 256				// Color Levels
#define MAX 400				// Maximum Dimension of Image

#define long int


// Global Shared variables to be used among all functions
unsigned long ImageWidthByte, ImageLength,ImageWidth, ImageSize;// Specify the dimensions of image
unsigned long file_position;									// set by read and used by write function
void *ImageBuffer,*OutBuffer;									// Stores the input and output raster image respectively
unsigned char *IB,*OB;											// Image processing pointer for In and out buffer
unsigned char pic_byte;											// Stores the picture byte for intermediate processing


// Function Declarations
void error_message(const char *error);
void read_tiff(unsigned char in_filename[14]);
void copy_tiff(unsigned char in_filename[14], unsigned char out_filename[14]);
void write_tiff(unsigned char out_filename[14]);
int screen();

int main()
{

	// variable declaration
	short unsigned int i, j,k;								// loop control variables
	int choice,avinfile = FALSE,avoutfile = FALSE,processing = FALSE;	
	unsigned char in_filename[14], out_filename[14];		// File names
	float temp;
	
	
	// Required for histogram transformations
	float r[L]={0},s[L],cr[L]={0},cs[L];									
	int z[L];	
	float xt,dxt,st[L];

	// variables for Fourier Transform
	unsigned long x,y,u,v,M,N ;	// M Length N Width
	int f[MAX][MAX];							
	double fi[MAX][MAX],fr[MAX][MAX],fre[MAX][MAX], fim[MAX][MAX];
	double theta,real,img,PI;
	FILE *fp;

	fp = fopen("filedata.txt","w+");
	if(fp == NULL)
	{
		printf("Write file could not be opened");
		exit(1);
	}

	PI=acos(-1.0); // Courtsey Bawa Ji


	// iterative choice taking loop
	do
	{
		choice = screen();

		// Ensuring the presence of Input File Before Processing commence
		if ((avinfile == FALSE) && (choice !=0))
		{
				printf("\nInput File NOT available\n");
				printf("Enter the name of Input TIFF File (*.tif): ");
				scanf( "%s", in_filename);
				read_tiff(in_filename);
				avinfile = TRUE;

		}
		
		
		// Ensuring the presence of Output File Before Processing commence
		if ((avoutfile == FALSE) && (choice !=0))
		{
				printf("\nOutput File NOT available\n");
				printf("Enter the name of Output TIFF File (*.tif): ");
				scanf( "%s", out_filename);
				copy_tiff(in_filename,out_filename);
				avoutfile = TRUE;

		}


		// Intializing the output buffer for processing		
		processing = TRUE;
		OutBuffer =  malloc(ImageSize);
		if(OutBuffer == NULL)
		{
				printf("\n\nMemory Required: %d Bytes", ImageSize);
				error_message(" Memory Allocation Error- Insufficient Memory. ");
		}
		IB = (unsigned char *)ImageBuffer;
		OB = (unsigned char *)OutBuffer;		

		switch(choice)
		{


		case  1 :	// Negate the image
					printf("\nNegating the Image......");
						for(i=0; i<ImageLength; i++)
						{
							for(j=0; j<ImageWidthByte; j++)
							{
								pic_byte = (*IB);
								pic_byte = (unsigned char)(255) - pic_byte;
								(*OB)= pic_byte;
								OB++;
								IB++;
								fprintf(fp, "%d\t", pic_byte);
							}
								fprintf(fp, "\n");

						}



						break;


		case  2 : // Log Transform
					printf("\nApplying Log Transformation..... ");
						for(i=0; i<ImageLength; i++)
						{
							for(j=0; j<ImageWidthByte; j++)
							{
								pic_byte = (*IB);
								temp = (float) ((L-1)/8.0)*log(1+pic_byte);
								(*OB)= (unsigned char) temp;
								OB++;
								IB++;
							}
						}						
						break;


		case  3 : // Power Law Transform
						float n;
						printf("\nApplying Power Law Transformation..... ");
						printf("\nEnter the value of n [0-1] to apply the power function: ");
						scanf("%f",&n);
						for(i=0; i<ImageLength; i++)
						{
							for(j=0; j<ImageWidthByte; j++)
							{
								pic_byte = (*IB);
								temp = (float)pow(L-1,1-n)*(pow(pic_byte,n));
								(*OB)= (unsigned char) temp;
								OB++;
								IB++;
							}
						}
						break;	    

		
		case  4 : // Gray Level Slicing
					printf("\nApplying Gray Level Slicing..... ");
						for(i=0; i<ImageLength; i++)
						{
							for(j=0; j<ImageWidthByte; j++)
							{
								pic_byte = (*IB);
								if (pic_byte<=100)
									temp = float (50.0/100.0*pic_byte);
								else
									if (pic_byte<=200)
										temp = float(50.0+175/100.0*(pic_byte-100));
									else
										temp = float(225.0+30.0/55.0*(pic_byte-200.0));

								(*OB)= (unsigned char)temp;
								OB++;
								IB++;
							}
						}
		
						break;

		
		case  5 : // Histogram Equalisation 
						// generating histogram transformation Function
					printf("\nApplying Histogram Equalisation..... ");
						for(i=0; i<ImageLength; i++)
						{
							for(j=0; j<ImageWidthByte; j++)
							{
								pic_byte = (*IB);
								r[pic_byte]++;
								IB++;
							}
						}
						cr[0]=r[0];
						s[0]=ImageSize/L;
						cs[0]=s[0];
						for(i=1;i<L;i++)
						{
							s[i] = ImageSize/L;
							cr[i]=cr[i-1]+r[i];
							cs[i] = cs[i-1]+s[i];
							
						}

				
						//Applying transformtion
						IB = (unsigned char *)ImageBuffer;
						for(i=0; i<ImageLength; i++)
						{
							for(j=0; j<ImageWidthByte; j++)
							{
								pic_byte = (*IB);
								k=0;
								while(cr[pic_byte]>cs[k])
									k++;

								(*OB) = (unsigned char)k;
								OB++;
								IB++;
							}
						}			
						break;

		
		case  6 : // Histogram Matching
						// generating histogram transformation Function r to s
					printf("\nApplying Histogram Matching..... ");
						for(i=0; i<ImageLength; i++)
						{
							for(j=0; j<ImageWidthByte; j++)
							{
								pic_byte = (*IB);
								if ((pic_byte>=0)&& (pic_byte<L))
									r[pic_byte]++;

								IB++;
							}
						}
						r[0]=r[0]/ImageSize;
						s[0]=r[0];
						for(i=0;i<L-1;i++)
						{
							if(r[i] !=0)
								r[i]=r[i]/ImageSize;
							s[i+1]=s[i]+r[i];
						}
						r[L-1]=r[L-1]/ImageSize;

						// generating histogram transformation Function s to z
						st[0]=0;
						xt=0.0;
						dxt=1.0/(L-1);
						for(k=0;k<L;k++)
						{
							st[k]=(1-exp(xt))/1.7182;
							if(st[k]<0)
								st[k]=-st[k];
							if(k%10==0)
								xt=xt+dxt;
						}
						for(k=0;k<L/2;k++)
						{
							for(i=0;i<L/2;i++)
							{
								if((((s[i]-st[i])*100)<=0.01)||((s[i]-st[i])*100)>0.01)
									z[k]=i;
							}
							for(i=L/2;i<L;i++)
							{
								if((((s[i]-st[i])*10)<=0.1)||((s[i]-st[i])*10)>0.1)
									z[i]=i;
							}
						}

						//Applying transformtion
						IB = (unsigned char *)ImageBuffer;
						for(i=0; i<ImageLength; i++)
						{
							for(j=0; j<ImageWidthByte; j++)
							{
								pic_byte = (*IB);
								if ((pic_byte>=0)&& (pic_byte<L))
									temp =  z[pic_byte];
								(*OB)= (unsigned char)temp;
								OB++;
								IB++;
							}
						}				
						break;

		case  7 : // Fourier Transform
					printf("\nApplying Fourier Transformation..... ");
						i=0;
						M = ImageLength; 
						N = ImageWidth;
						for(x=0; x<M; x++)
						{
							for(y=0; y<N; y++)
							{
								f[x][y] = (*IB);
								IB++;
							}
						}
						
						//Apply Fourier Function
						printf("\n Procesing Please wait........");
						for(u=0; u<M; u++)
						{
							for(v=0; v<N; v++)
							{
								real=0.0;
								img=0.0;
								for(x=0;x<M;x++)
								{
									for(y=0;y<N;y++)
									{
										theta = 2.0 *PI*((double(u*x)/double (M))+(double (v*y)/double(N)));
										real = real+f[x][y]*cos(theta);
										img  = img +f[x][y]*sin(theta);	
									}
								}
								fr[u][v]=real/double(M*N);
								fi[u][v]=-img/double(M*N);

							}
						}	
						printf("\n Fourier Function Application Successful...");
						
						//Apply Inverse Fourier Function

						printf("\n Now applying Inverse Fourier Function please wait...");
						for(x=0; x<M; x++)
						{
							for(y=0; y<N; y++)
							{
								real=0.0;
								img=0.0;
								for(u=0;u<M;u++)
								{
									for(v=0;v<N;v++)
									{
										theta = 2.0 *PI*((double(u*x)/double (M))+(double (v*y)/double(N)));
										real = real+(double)fr[u][v]*(double)cos(theta)-(double)fi[u][v]*(double)sin(theta);
										img  = img +(double)fr[u][v]*(double)sin(theta)+(double)fi[u][v]*(double)cos(theta);	
									}
								}
								fre[x][y]=(float)real;
								fim[x][y]=(float)img;

								(*OB)=(unsigned char) fre[x][y];
								OB++;

							}
						}	
						printf("\n Inverse Fourier Function Successful...");


						break;

		case  8 : // Getting Input File
					if(ImageBuffer != NULL)
						free(ImageBuffer);
					printf("Enter the name of Input TIFF File (*.tif): ");
					scanf( "%s", in_filename);
					read_tiff(in_filename);
					avinfile= TRUE;
					if (avoutfile == TRUE)
						copy_tiff(in_filename,out_filename);
					

					processing = FALSE;

					break;

		case  9 : // Getting Output File
					printf("Enter the name of Output TIFF File (*.tif): ");
					scanf( "%s", out_filename);
					copy_tiff(in_filename,out_filename);
					avoutfile = TRUE;
					processing = FALSE;
					break;

		case  0 : //  Free Image Buffer and exit
					printf("\n\n Exiting..... \n");
					exit(0);
		default : printf("\nWrong choice  ");
			  break;


	    }
		if(processing == TRUE)
		{
			printf("\nWriting to output File...");
			write_tiff(out_filename);
			printf("\n\n Processing SUCCESSFUL..... \n");
			free(OutBuffer);
		}
  }while(choice != 0);

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
       printf("\n\t�      1.    Negate                          �");
       printf("\n\t�      2.    Log Transfom                    �");
       printf("\n\t�      3.    Power Law                       �");
       printf("\n\t�      4.    Gray Level Slicing              �");
       printf("\n\t�      5.    Histogram Equalisation          �");
       printf("\n\t�      6.    Histogram Matching              �");
       printf("\n\t�      7.    Fourier Transform               �");
       printf("\n\t�      8.    Select Input File               �");
	   printf("\n\t�      9.    Select Output File              �");
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
void read_tiff(unsigned char in_filename[14])

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
					printf("\nImage Width = %d",value);
					break;
		
		case 257 :	// Number of Pixel Rows
					ImageLength = value;
					printf("\nImage Length = %d",value);
					break;
		
		case 258 :	// Number of bit for each color
					BitsPerSample = (short)value;
					printf("\nBits per sample = %d",value);
					break;

		case 259 :  // Check if image is compressed
					if(value !=1)
						error_message(" Image is Compressed.");
					break;

		case 262 :	// Photochrometic Interpretation (for us it should be 0-3)
					if(value<0 || value >3)
						error_message("Image Type Cann't be processed by the demo....");

					printf("\nPhotochrometic Level: %d", value);
					if(value ==1)	// Black and White image
						BitWhite =1;
					else 
						BitWhite =0;

					break;

		case 273 :	// Number of strips and their relative address w.r.t Start

					NoOfStrips = count;
					StripOffsets= value;
					printf("\nNumber of strips = %d",count);
					break;

		case 278 :	//Number of rows per strip.
					RowsPerStrip = value;
					printf("\nNumber of rows in a strip = %d",value);
					StripsPerImage = (ImageLength + RowsPerStrip -1)/RowsPerStrip;
					break;

		case 279 :	//Number of Bytes in a strip
					StripByteCounts = value;
					printf("\nNumber of bytes in a strip = %d",value);
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
void copy_tiff(unsigned char in_filename[14], unsigned char out_filename[14])
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

void write_tiff(unsigned char out_filename[14])

{

	// Variable Declaration

	FILE *outfile_ptr;
	unsigned char pic_byte;
	short unsigned int i, j;						// loop control variables

	// Reopen the output file in read-write mode

	if((outfile_ptr = fopen((const char *)out_filename, "rb+")) ==NULL)
	error_message(" Output File reopening Failed.");


	//positioning the outfile pointer
	fseek(outfile_ptr, file_position, SEEK_SET);


	// Saving he Raster Image from outbuffer to Output File
	OB = (unsigned char *)OutBuffer;

	for(i=0; i<ImageLength; i++)
	{
		for(j=0; j<ImageWidthByte; j++)
		{
			pic_byte = (*OB);
			fwrite((const void *) &pic_byte, (size_t)sizeof(pic_byte), (size_t)1, outfile_ptr);
			file_position = (unsigned long)ftell(outfile_ptr);
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
	exit(0);
}
