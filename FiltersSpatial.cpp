/*******************************************************************************
*       File Name: ProcessImage.cpp 
*		Purpose  : Program for performing follwing SPATIAL TRANSFORMATION
*					on an UNCOMPRESSED  Binary, Gray-Level, or Colored  Image
*					1. Min, Max, Mean Filter 
*					2. Average Filters
*					3. Laplacian Filters
*					4. Robert's Cross Gradient Filter
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


// Function Declarations
int screen();
void error_message(const char *error);
void read_tiff(char in_filename[14]);	// read to IB
void copy_tiff(char in_filename[14], char out_filename[14]);
void write_tiff(char out_filename[14]);	// Write from OB

void print_matrix(int a[MAX][MAX], int size);
void print_op_matrix(int a[WSIZE][WSIZE], int wsize);
void normalize(int out[MAX][MAX], int size);
void apply_operator(int opw[WSIZE][WSIZE],int wsize, int size, int temp[MAX+2][MAX+2],int out[MAX][MAX]);
int  apply_operator_window(int w[WSIZE][WSIZE], int im[WSIZE][WSIZE], int wsize);


int  main()
{

	// variable declaration

	char in_filename[14], out_filename[14];		// File names

	int size,wsize = 3 ; // default image size, operator window size    
	int in[MAX][MAX];
    int ex [5][5] = {{1,3,7,2,6},{5,1,4,7,2},{1,3,1,2,4},{6,6,5,7,0},{1,0,3,4,1}};
    int temp[MAX+2][MAX+2];
    int out[MAX][MAX];
    int opw[WSIZE][WSIZE]; // operator window
    int a[WSIZE*WSIZE],min,max,mean;
    int i, j, k,m,s,t,sorted, choice, sum,chmask, sign, den;
	double outfloat,avg;	
	
	
	// Default intialization
	/*strcpy(in_filename,"image1.tif");
	strcpy(out_filename,"out1.tif");
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

	if(ImageLength==ImageWidthByte)
		size = ImageLength;
	else
		error_message("Program runs only for square images");

	// intialize in matrix
	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
			in[i][j] = (*IB);
			IB++;
		}
	}

	// intialize temp with zeros padded to input
    	for(i=0;i<size+2;i++)
       for(j=0;j<size+2;j++)
		temp[i][j] = 0;

   	//copy input to temp
   	for(i=0;i<size;i++)
       for(j=0;j<size;j++)
	  	temp[i+1][j+1] = in[i][j];*/


	// iterative choice taking loop
	do
	{
	   choice = screen();
	   printf("\n");
	   switch(choice)
	   {
		case  1 : // Min Filter
			for(i=0;i<size;i++)
			{
			  for(j=0;j<size;j++)
			  {
			    min = in[i][j];
			    for(s=0;s<wsize;s++)
			    for(t=0;t<wsize;t++)
				if(temp[i+s][j+t]<min)
				    min = temp[i+s][j+t];

			    out[i][j] = min;

			  } //end j
			}//end i
			printf("\nThe output of min filter is in %s: \n", out_filename);
			break;


		case  2 : // Max Filter
			for(i=0;i<size;i++)
			{
			  for(j=0;j<size;j++)
			  {
			    max = in[i][j];
			    for(s=0;s<wsize;s++)
			    for(t=0;t<wsize;t++)
			       if(temp[i+s][j+t]>max)
				    max = temp[i+s][j+t];
			    out[i][j] = max;
			} //end j
			}//end i
			printf("\nThe output of max filter is in %s \n",out_filename);
			break;


		case  3 : // Mean Filter
			for(i=0;i<size;i++)
			{
			  for(j=0;j<size;j++)
			  {
			    min = in[i][j];
			    max = in[i][j];
			    for(s=0;s<wsize;s++)
			    for(t=0;t<wsize;t++)
			    {
				   if(temp[i+s][j+t]<min)
				    min = temp[i+s][j+t];
				   if(temp[i+s][j+t]>max)
				    max = temp[i+s][j+t];
			    }
			    mean = (min+max)/2;
			    out[i][j] =mean;

			  } //end j
			}//end i
			printf("\nThe output of mean filter is in %s\n",out_filename);
			break;


		case  4 : // Median Filter
			for(i=0;i<size;i++)
			{
			  for(j=0;j<size;j++)
			  {
			    k=0;
			    for(s=0;s<wsize;s++)
			    for(t=0;t<wsize;t++)
			      a[k++] = temp[i+s][j+t];

			    // bubble sort
			    m = wsize*wsize;
			    sorted = FALSE;
			    while((m>1) && (sorted == FALSE))
			    {
				for(k=0;k<m-1;k++)
				{
				    sorted = TRUE;
				    if(a[k]>a[k+1])
				    {
					t = a[k];
					a[k] = a[k+1];
					a[k+1] = t;
					sorted = FALSE;
				    }
				}
				m--;
			     }//end while

			    out[i][j] = a[wsize*wsize/2];

			  } //end j
			}//end i
			printf("\nThe output of median filter is in %s \n",out_filename);
			break;


		case  5 : // Average Filter
			// define operator
			for(s=0;s<wsize;s++)
			  for(t=0;t<wsize;t++)
			    opw[s][t] = 1;

			den = wsize*wsize;

			apply_operator(opw,wsize,size,temp,out);


			for(i=0;i<size;i++)
			  for(j=0;j<size;j++)
			  {
			    outfloat = (double) out[i][j];
			    avg = outfloat/((double) den);
			    out[i][j] = (int) (avg+0.5);
			  }

			printf("\nThe output of average filter is in %s \n", out_filename);
			break;


		case  6 : // Laplacian
			//define operator
			printf("\nEnter choice for Laplacian Mask");
			do{
			 printf("\n[1] for gravity -4 [2] for gravit 4 [3] for gravit -8 [4] for gravity 8: ");
			 scanf("%d",&chmask);
			}while((chmask<1)||(chmask>4));
			if(chmask<3)	// gravity 4 or -4
			{
			  sign = (int) pow(-1,chmask);
			  opw[0][0] = opw[0][2] = opw[2][0] = opw[2][2] = 0;
			  opw[0][1] = opw[1][0] = opw[2][1] = opw[1][2] = -1*sign;
			  opw[1][1] = 4*sign;
			}else
			if(chmask>2)	// gravity 8 or -8
			{
			  sign = (int) pow(-1,chmask);
			  for(i=0;i<wsize;i++)
			  for(j=0;j<wsize;j++)
				opw[i][j] = -1*sign;
			  opw[1][1] = 8*sign;
			}
			printf("\nThe Laplacian operator is");
			print_op_matrix(opw,3);
			sum = 0;
			for(s=0;s<wsize;s++)
			  for(t=0;t<wsize;t++)
			    sum = sum+opw[s][t];
			if(sum != 0)
			 printf("\nError defining operator");
			else
			{
			  apply_operator(opw,wsize,size,temp,out);
			  normalize(out, size);
			  printf("\nThe output of Laplacian filter is in %s \n", out_filename);
			}
			break;

		case  7 : // Sobel Horz and vertical
			//define Horizontal operator
			do{
			 printf("\nPress [1] for Horizontal [2] for Vertical: ");
			 scanf("%d",&chmask);
			}while((chmask<1)||(chmask>2));
			if(chmask == 1)
			{
			    opw[0][0] = -1; opw[0][1] = -2; opw[0][2] = -1;
			    opw[1][0] =  0; opw[1][1] =  0; opw[1][2] =  0;
			    opw[2][0] =  1; opw[2][1] =  2; opw[2][2] =  1;
			    printf("\nThe Sobel's Horizontal operator is");
			}else
			if(chmask == 2)
			{
			    opw[0][0] = -1; opw[0][1] = 0; opw[0][2] = 1;
			    opw[1][0] = -2; opw[1][1] = 0; opw[1][2] = 2;
			    opw[2][0] = -1; opw[2][1] = 0; opw[2][2] = 1;
			    printf("\nThe Sobel's Vertical operator is");
			}
			print_op_matrix(opw,3);
			sum = 0;
			for(s=0;s<wsize;s++)
			  for(t=0;t<wsize;t++)
			    sum = sum+opw[s][t];
			if(sum != 0)
			 printf("\nError defining operator");
			else
			{
			  apply_operator(opw,wsize,size,temp,out);
			  normalize(out, size);
			  printf("\nThe output of Sobel's filter is in %s \n",out_filename);
			}
			break;

		case  8 : // Robert's Cross Gradient Operator
			do{
			 printf("\nPress [1] for Horizontal [2] for Vertical: ");
			 scanf("%d",&chmask);
			}while((chmask<1)||(chmask>2));
			if(chmask == 1)
			{
			    opw[0][0] = -1; opw[0][1] = 0;
			    opw[1][0] =  0; opw[1][1] = 1;
			    printf("\nThe Robert's Cross Gradient Horizontal operator is");
			}else
			if(chmask == 2)
			{
			    opw[0][0] = 0; opw[0][1] = -1;
			    opw[1][0] = 1; opw[1][1] = 0;
			    printf("\nThe Robert's Corss Gradient Vertical operator is");
			}
			print_op_matrix(opw,2);
			sum = 0;
			for(s=0;s<2;s++)
			  for(t=0;t<2;t++)
			    sum = sum+opw[s][t];
			if(sum != 0)
			 printf("\nError defining operator");
			else
			 {
			   apply_operator(opw,2,size,temp,out);
			   normalize(out, size);
			   printf("\nThe output of Robert's Cross Gradient filter is in %s \n",out_filename);
			 }
			break;

		case  9 : // Getting Input File
			free(ImageBuffer);
			free(OutBuffer);
			
			printf("Enter the name of Input TIFF File (*.tif): ");
			scanf( "%s", in_filename);
			read_tiff(in_filename);
			if(ImageLength==ImageWidthByte)
				size = ImageLength;
			else
				error_message("Program runs only for square images");

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

			// intialize in matrix
			for(i=0; i<size; i++)
			{
				for(j=0; j<size; j++)
				{
					in[i][j] = (*IB);
					IB++;
				}
			}
			// intialize temp with zeros padded to input
    		for(i=0;i<size+2;i++)
				for(j=0;j<size+2;j++)
					temp[i][j] = 0;

   			//copy input to temp
   			for(i=0;i<size;i++)
				for(j=0;j<size;j++)
	  				temp[i+1][j+1] = in[i][j];

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
		for(i=0; i<size; i++)
		{
			for(j=0; j<size; j++)
			{
				(*OB) = out[i][j];
				OB++;
			}
		}

		write_tiff(out_filename);
		printf("\n\n Processing SUCCESSFUL..... \n");
	}


  }while(choice != 0);
  free(OutBuffer);
}//main program ends here
/*********************************************************************************************************
	apply operator on entire input
**********************************************************************************************************/
void apply_operator(int opw[WSIZE][WSIZE],int wsize, int size, int temp[MAX+2][MAX+2],int out[MAX][MAX])
{
  int i,j,s,t;
  int imw[WSIZE][WSIZE];

  for(i=0;i<size;i++)
  {
    for(j=0;j<size;j++)
    {
      // extract image window
      for(s=0;s<wsize;s++)
      for(t=0;t<wsize;t++)
	imw[s][t] = temp[i+s][j+t];

      out[i][j] = apply_operator_window(opw,imw,wsize);

   } //end j
  }//end i
}

/*********************************************************************************************************
	apply operator on window
**********************************************************************************************************/
int apply_operator_window(int w[WSIZE][WSIZE], int im[WSIZE][WSIZE], int wsize)
{
   int i,j,fxy=0;
   for(i=0;i<wsize;i++)
     for(j=0;j<wsize;j++)
	fxy = fxy+w[i][j]*im[i][j];
   return(fxy);
}

/*********************************************************************************************************
	Normalize the values
**********************************************************************************************************/
void normalize(int out[MAX][MAX], int size)
{
  int i,j;

  for(i=0;i<size;i++)
  {
    for(j=0;j<size;j++)
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
	Prints the matrix
**********************************************************************************************************/
void print_matrix(int a[MAX][MAX], int size)
{
   int i,j;

   for(i=0;i<size;i++)
    {
	printf("\n\n");
	for(j=0;j<size;j++)
	 printf("%d\t",a[i][j]);
    }
}

/*********************************************************************************************************
	Prints the Operator matrix
**********************************************************************************************************/
void print_op_matrix(int a[WSIZE][WSIZE], int wsize)
{
   int i,j;

   for(i=0;i<wsize;i++)
    {
	printf("\n\n");
	for(j=0;j<wsize;j++)
	 printf("%d\t",a[i][j]);
    }
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
       printf("\n\t�      1.    Min Filter                      �");
       printf("\n\t�      2.    Max Filter                      �");
       printf("\n\t�      3.    Mean Filter                     �");
       printf("\n\t�      4.    Median Filter                   �");
       printf("\n\t�      5.    Average Filter                  �");
       printf("\n\t�      6.    Laplacican Filter               �");
       printf("\n\t�      7.    Sobel Horz & Vert Filter        �");
       printf("\n\t�      8.    Robert's Cross Gradient Filter  �");
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
					if(value !=1 )
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
	exit(0);
}
