/*
C. Moos
Template: fitsslpitbayer.c from Mischa Schirmer
Revise:  5  01.07.2020
cmoos@gmx.de
Interpolationsmethoden 
fitsiodemosaicbayer.c
Changes 1.3: black edge reduced
Changes 1.4: Edge smoothed, Error in PPG last interpolation in row and column was wrong
*/
/*
Changes 1.4 Rev3: Usage Text missleading: 2.21 1.00 1.51 1.00 = RGBg
Changes 1.5 Rev1: Added median filter
                  Added L[ab] export
                  200 exchanged with FILEMAX
changes 1.6 Rev1: make header interpolate.h
                  removed fitstools.h, not c99 conform
                  did not compile!! retun to source file without header interpolate.h
                  added colorsystems.h
changes 1.7     : -q 3 WCAPI added
                  diff, sum, median3x3, fmedian3x3, 
                  fill_intensities5x5[_2], fill_intensities7x7
                  filterm added
                  -t option added for generating test-patterns, not really a feature, mismatches other than RGGB, not complete
 Rev 2:
                   qfits_header_add for history information added
                   switch -q 3 in usage added
                   
 Rev 3: 
                   second try :declarations in demosaic.h
                   still doesn not compile, gave up
                   will deliver _demosaic.c, _demosaic.h, _myfitstools.c, _myfitstools.h for packaging by mischa
Rev 4:
                   fill_intensities5x5 and fill_intensities5x5_2 combined to fill_intensities5x5

open optimizations: treating the outer frame, 
                   fill_intensities_2 integration in fill_intensities
                   -m with pass > 2 is very slow
                   -d despeckle generates sometimes crosses
                  for fmedian3x3 maximum memory consumption is about 9ximage!!HEAVY didn't find a way to reduce. 
                  Maybe some functions are already in fitstools included, didn't look for	e.g fitslaplace, fitsmedian
                  sometimes functions works with reference to image, sometimes they don't (sum and diff versus median3x3 and fmedian3x3)
                  despeckle seems to become obsolete, related to success of artifact reducing bei fmedian,	then median3x3 is obsolete too
                  CIE-Lab has to proofe needs, maximum level is reduced to about 400 max, due to a 8 Bit model, should be watched further
                  -t for other than RGGB does not work properly


*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
//#include "eclipse.h"
#include "fitsio.h"

#include "fitsiotools.h"
//#include "myfitstools.h"//
//#include "demosaic.h"  doesnt compile because of fitstools.h
#include "colorsystems.h"   // CIEL*a*b

#define RED 1
#define GREEN1 2  // in blue line
#define BLUE 3
#define GREEN2 4  // in red line

// compiled with
// gcc -o fitsiodemosaicbayer fitsiodemosaicbayer.c colorsystems.c  -L.  -lcfitsio -lm //-lnsl

void usage(int i, char *argv[])
{
  if (i == 0) {
    fprintf(stderr,"\n");
    fprintf(stderr,"          Version 1.7 Rev 5 (2020-07-01) fitsio\n\n");
    fprintf(stderr,"  Author: Carsten Moos\n\n");
    fprintf(stderr,"  USAGE:  %s \n", argv[0]);
    fprintf(stderr,"           -i input_image \n");
    fprintf(stderr,"           -p pattern (available: RGGB, GRBG, GBRG, BGGR)\n");
    fprintf(stderr,"           -q 0    for bilinear Interpolation\n");
    fprintf(stderr,"           -q 1    for gradient Interpolation\n");
    fprintf(stderr,"           -q 2    for PPG Interpolation\n");
    fprintf(stderr,"           -q 3    for WACP Interpolation\n");
    fprintf(stderr,"           -c 2.21 1.00 1.51 1.006 for bayerfilter RGBG colorbalance\n");
    fprintf(stderr,"           -o 256  for Offsetcorrection before RGBG colorbalance\n");
    fprintf(stderr,"           -m 2 for 1-5 times artifact filter 3x3 of crominannce\n");
    fprintf(stderr,"           -d 1000 for despeckle above level 1000  filter 1x3 of crominannce\n");
    fprintf(stderr,"           -l with export of Luminance CIEL*a*b\n");
    fprintf(stderr,"           -t with generates RAW testimage: R=50 G=101[2] B=75\n");
    fprintf(stderr,"  PURPOSE: Demosaics a bayer matrix fits image into\n");
    fprintf(stderr,"  three monochrome fits images for RGB.\n\n");

    exit(1);
  }
}

// sollte nach library demosaic.c,
//laesst sich aber mit der fitstools.h compilieren ??
// in demosaic.c
float hue_transit(float l1,float l2,float l3,float v1,float v3)
{
// for PPG
//printf("hue_transit: l1 %5.1f l2 %5.1f l3 %5.1f v1 %5.1f v3 %5.1f\n",l1,l2,l3,v1,v3);
  if((l1<l2 && l2<l3) || (l1> l2 && l2 > l3))
     return(v1+(v3-v1) * (l2-l1)/(l3-l1));
  else
     return((v1+v3)/2.0+(l2*2.0-l1-l3)/4.0);
}

int direction(float N, float E, float W, float S)
{
// for PPG
if (N < E && W < S)
    if ( N < W)
        return(1);
    else
        return(3);
else
    if (E < S)
        return (2);
    else
        return (4);
}

void fill_intensities5x5(int k,int n, float *in,float array[6][6]){
// array is 6x6, but used for use 5x5 without null index
// k is actual pointer position in image
// n is images width
// array is local reference to filtermatrix
// input image pointer maybe usable with fill_intensities() as replacement
/*
FILL_INTENSITIES: fills an array of a givel center (k) for better access to surrounding
elements

I11 I12 I13 I14 I15
I21 I22 I23 I24 I25
I31 ..          I35
I41 ..          I45
I51 I52 I53 I54 I55

*/
int i,j;
for (j=1;j<=5;j++){ // rows
  for(i=1;i<=5;i++){  // cols
   array[j][i]=in[k + (i-3) + (j-3)*n];
   k=k;
  }
}

}


//void fill_intensities7x7(long k,int n,image_t *in,float array[8][8]){
void fill_intensities7x7(long k,int n,float *in,float array[8][8]){
// array is 8x8, but used for use 7x7 without null index
// k is actual pointer position in image
// n is images width
// array is local reference to filtermatrix
/*
FILL_INTENSITIES: fills an array of a givel center (k) for better access to surrounding
elements
I11 I12 I13 I14 I15 I16 I17
I21 I22 I23 I24 I25 I26 I27
I31 ..                  I37
I41 ..                  I47
I51 ..                  I57
I61 ..                  I67
I71 I72 I73 I74 I75 I76 I77

*/

int i,j;
for (j=1;j<=7;j++){ // rows
  for(i=1;i<=7;i++){  // cols
   //array[j][i]=in->data[k+ (i-4) + (j-4)*n];
   array[j][i]=in[k+ (i-4) + (j-4)*n];
  }
}


}

void despeckle (int m,int n,int filterdirection ,float level,float *in, float *out){
// m is image height
// n is image width
// filterdirection can be 1,2,3 or 4 for experimental
// level is a threshold, only pixelvalues above will be filtered
// in is a reference pointer to input_image
// out is a reference pointer to out_image
/*
DESPECKLE: median filter of 3 pixelvalues (3 cells)
to avoid hotpixel.
Filter is rather weak, but sometimes crosses appear.
*/
int k,filtersize=1;
int j,jt, i,it;
int cell;
long counter=0;

float filter[9];
//printf("despeckle filter, ");
 for (j=0; j<m; j++) { //row
      for (i=0; i<n; i++) {  // col
      if (filterdirection <1 || filterdirection >5){
         printf("no direction recognized, copying");
         out[i+n*j]=in[i+n*j];
      }
      else
      {
      k = 0; cell=0;
       if ( in[i+n*j]>level ){
        for (jt=j-filtersize;jt<=j+filtersize;jt++) {
	  for (it=i-filtersize;it<=i+filtersize;it++) {
	      if (it>=0 && jt>=0 && it<n && jt<m) {
                  cell=it-i+2+(jt-j+2-1)*3;
	          if (filterdirection==1 && (cell ==3 || cell==5 || cell==7) ){
	             filter[k] = in[it+n*jt];
	             k++;
	          }
	          if (filterdirection==2 && (cell ==1 || cell==5 || cell==9) ){
	             filter[k] = in[it+n*jt];
	             k++;
	          }
	          if (filterdirection==3 && (cell ==2 || cell==5 || cell==8) ){
	              filter[k] = in[it+n*jt];
	              k++;
                  }
	          if (filterdirection==4 && (cell ==4 || cell==5 || cell==6) ){
	              filter[k] = in[it+n*jt];
	              k++;
                  }
	      }
           }//it
         }//jt

       qsort(filter, k, sizeof(float), compare);
       if (k>1) { counter++;
	 if (k % 2 == 0)
	   out[i+n*j] = 0.5*(filter[k/2-1] + filter[k/2]);
	 else  
           out[i+n*j] = filter[k/2];    
       }
       else  // k=1
         {out[i+n*j]=in[i+n*j]; counter++;}
      }// below level
      else
         out[i+n*j]=in[i+n*j];
     }
    }// i
  }//j
}

void filter_m(int m, int n, float *in,float *out){
// m is image height
// n is image width
// in is reference to input image
// out  is reference to output image
/*
FILTER_M:
convolve input image with Laplacian operator LP
       1  9  1
1/11 * 9 -40 9
       1  9  1 
*/

int j,i,r,s;
float I[6][6],OP[6][6]={ {1, 9, 1},{ 9,-40, 9},{ 1, 9, 1} }, weight;
int T;
T=15; // experimental, suggestion of Lu & Tan

weight=0;
 for (j=0; j<m; j++) { //row
      for (i=0; i<n; i++) {  // col
        if (i>2&&i<n-2&&j>2&&j<m-2){
           //f=*(in[i-1+ n*(j-1)] +in[i+ n*(j-1)]*9 +in[i+1+ n*(j-1)]+ 
           // in[i-+n*j]*9 +in[i+ n*j]*40 + in[i+1+n*j]*9 + 
           //in[i-1+ n*(j+1)] + in[i+ n*(j+1)]*9 + in[i+1+ n*(j+1)]); 
            for (r=1;r<=3;r++){
                for(s=1;s<=3;s++){
                   weight+=I[r][s]*OP[r][s];
                }
            }
           weight=weight/11.0;
           if (weight>T) out[i+n*j]=1;
           else out[i+n*j]=0;
       }
      }
 }
}

void median3x3(int m, int n ,float *in ){
// m is image height
// n is image width
// in is refernce to input image, will be replaced by results
/*
MEDIAN3x3: does a median of surrounding 3x3 fields

*/
float *out;
int k,filtersize=1;
int j,jt, i,it;
out= (float*) calloc(m*n, sizeof(float));
float filter[9];
 for (j=0; j<m; j++) { //row
      for (i=0; i<n; i++) {  // col
      k = 0;
        for (jt=j-filtersize;jt<=j+filtersize;jt++) {
	  for (it=i-filtersize;it<=i+filtersize;it++) {
	      if (it>=0 && jt>=0 && it<n && jt<m) {
	        filter[k] = in[it+n*jt];
	        k++;
	      }
           }//it
         }//jt
       
       qsort(filter, k, sizeof(float), compare);
       if (k>1) {
	 if (k % 2 == 0)
	   out[i+n*j] = 0.5*(filter[k/2-1] + filter[k/2]);
	 else
           out[i+n*j] = filter[k/2];    
       }
       else // outside
         out[i+n*j]=in[i+n*j];
    }// i
//    if ((j%100)==0) printf(".");
  }//j
//copy out nach in
//printf("copy to out, ");
 for (j=0; j<m; j++) { //row
      for (i=0; i<n; i++) {  // col
         in[i+n*j]=out[i+n*j];
      }
 }
//printf("done\n");
free(out);
}

void fmedian3x3(int m, int n ,float *in ,float *mask ){
// m is image height
// n is image width
// in is refernce to input image, will be replaced by results
// mask is refernce to a pixelmask, same size as in
/*
FMEDIAN3x3: does a median of surrounding 3x3 fields if local mask is 1
is for artifact reduction, Lu & Tan
*/

float *out;
int k,filtersize=1;
int j,jt, i,it;
out= (float*) calloc(m*n, sizeof(float));
float filter[9];
 for (j=0; j<m; j++) { //row
      for (i=0; i<n; i++) {  // col
      if (mask){
      k = 0;
        for (jt=j-filtersize;jt<=j+filtersize;jt++) {
	  for (it=i-filtersize;it<=i+filtersize;it++) {
	      if (it>=0 && jt>=0 && it<n && jt<m) {
	        filter[k] = in[it+n*jt];
	        k++;
	      }
           }//it
         }//jt
       
       qsort(filter, k, sizeof(float), compare);
       if (k>1) {
	 if (k % 2 == 0)
	   out[i+n*j] = 0.5*(filter[k/2-1] + filter[k/2]);
	 else
           out[i+n*j] = filter[k/2];    
       }
       else // outside
         out[i+n*j]=in[i+n*j];
     }//masked 1
     else out[i+n*j]=in[i+n*j];
    }// i
  }//j
//copy out nach in

 for (j=0; j<m; j++) { //row
      for (i=0; i<n; i++) {  // col
         in[i+n*j]=out[i+n*j];
      }
 }

free(out);
}

void diff(int m,int n,float *in1,float *in2,float *out){
// m is image height
// n is image width
// in1 is refernce to input image1 
// in2 is refernce to input image2, will be subtracted from image1
// out is refernce to output

int i, j;
//printf("diff ");
 for (j=0; j<m; j++) { //row
      for (i=0; i<n; i++) {  // col
         out[i+n*j]=in1[i+n*j]-in2[i+n*j];
      }
 }
//printf("done\n");
}

void sum(int m,int n,float *in1,float *in2, float *out,float factor){
// m is image height
// n is image width
// in1 is refernce to input image1 
// in2 is refernce to input image2, will be subtracted from image1
// out is refernce to output
// factor is multiplicator to weight the sum, could mostly be 1.0

int i, j;
//printf("sum ");
 for (j=0; j<m; j++) { //row
      for (i=0; i<n; i++) {  // col
         out[i+j*n]=(in1[i+j*n]+in2[i+j*n])*factor;
      }
 }
}

//
// end for export to library demosaic.c
void printerror( int status);  // fitsio
void writeoutputfile(char *name, fitsfile *old , fitsfile *new, float *data,long *w_npixels, int *w_status); // fitsio
int main(int argc, char *argv[])
{

//******************************************************
// Demosaics a bayer matrix fits image (CFA)
//******************************************************

  int n, m, ns, ms, flag_q, flag_l, flag_t, flag_d,median, length;
  float R=1.0,G=1.0,g=1.0,B=1.0,o=0.0,level;  // Colorbalance of 4 Colorfilters
  int xoffset,yoffset; // for color determing
  long i, j, k1;
  float H,V;  // for Gradient horizontal, vertical
  float DN,DE,DW,DS,dne,dnw;   // for PPG
  float *channel1, *channel2, *channel3, *L=NULL;
  char input_image[FILEMAX], pattern[FILEMAX], *tmp6;
  char out1[FILEMAX], out2[FILEMAX], out3[FILEMAX], out4[FILEMAX];
  //image_t *image_in, *image_out1, *image_out2, *image_out3, *image_out4=NULL;
  //qfits_header *header;
// new for libchange
struct image_t {
  int lx;
  int ly;
  float *data;
};
int status=0, nkeys, keypos, nfound, anynull;  // fitsio
char card[FLEN_CARD];                          // fitsio
fitsfile *fptr,*ofptr1,*ofptr2,*ofptr3,*ofptr4;         // fitsio
long naxes[2];                                 // fitsio
long nbuffer, fpixel, npixels;                 // fitsio
float nullval;                                 // fitsio

  struct image_t *image_in, iimg1;               // fitsio
  char buffer[30]; // for History information 

  // print usage if no arguments were given

  // default: Linear interpolation
  flag_q = 0;
  flag_l = 0;
  flag_d = 0;
  flag_t = 0;
  level= 0;
  median = 0;
  k1=0;
  
  if (argc==1) usage(0, argv);
  
  for (i=1; i<argc; i++) {
      if (argv[i][0] == '-') {
	  switch(tolower((int)argv[i][1])) {
	      case 'i': strcpy(input_image,argv[++i]);
		  break;
	      case 'p': strcpy(pattern,argv[++i]);
		  break;
	      case 'q': flag_q = atoi(argv[++i]);
		  break;
	      case 'c': R = atof(argv[++i]);
	                G = atof(argv[++i]);
	                B = atof(argv[++i]);
	                g = atof(argv[++i]);
		  break;
	      case 'o': o = atof(argv[++i]);
                  break;
	      case 'm': median = atoi(argv[++i]);
                  break;
	      case 'd': level = atof(argv[++i]);flag_d=1;
                  break;
	      case 'l': flag_l = 1;
		  break;
	      case 't': flag_t = 1;
		  break;

	  }
      }
  }

  if (! (strcmp(pattern,"RGGB") == 0 || strcmp(pattern,"GRBG") == 0 || 
         strcmp(pattern,"GBRG") == 0 || strcmp(pattern,"BGGR") == 0)) {
      printf("\nBayer pattern not recognised. Nothing will be done.\n\n");
      exit (0);
  }

  // read the FITS header and the data block
//  checkfile(input_image);
if (fits_open_file(&fptr,input_image,READONLY, &status)){
	printerror( status );}
  //header = qfits_header_read(input_image);
if (fits_read_keys_lng(fptr,"NAXIS",1,2,naxes,&nfound, &status)){
	printerror( status );}
  //image_in = image_load(input_image);
  
npixels=naxes[0]*naxes[1];
fpixel=1;
nullval=0;
//image_in1->data[0] =(float *) malloc (npixels * sizeof (float));
iimg1.data =(float *) malloc (npixels * sizeof (float));
iimg1.data[npixels-1]=0;
if (fits_read_img(fptr,TFLOAT,fpixel,npixels,&nullval,iimg1.data, &anynull,&status)){
	printerror( status );}

  
// for new lib
iimg1.lx=naxes[0];
iimg1.ly=naxes[1];
// Uebergabe 
image_in=&iimg1;


  n = image_in->lx;
  m = image_in->ly;
  
  // chop the last row / column of pixels if the image dimensions are uneven
  if ( n % 2 != 0) n = n - 1;
  if ( m % 2 != 0) m = m - 1;
  
  ns = n / 2;
  ms = m / 2;
  
  // the file name prefix
  length = strlen(input_image) - 5;    // Oh har, hardcodiertes .fits
 // tmp6 = strndup(input_image, length);
  input_image[length]='\0';
  tmp6=input_image;
  strcpy(out1, tmp6);
  strcpy(out2, tmp6);
  strcpy(out3, tmp6);
  strcpy(out4, tmp6);

// cut all pattern to RGGB
 if (strcmp(pattern,"RGGB") == 0){
    xoffset=0;
    yoffset=0;
 }
 if (strcmp(pattern,"GRBG") == 0){
    xoffset=1;
    yoffset=0;
    
 }
 if (strcmp(pattern,"GBRG") == 0){
    xoffset=0;
    yoffset=1;
 }
 if (strcmp(pattern,"BGGR") == 0){
    xoffset=1;
    yoffset=1;
 }
/*
RGRGRGRGR
gBgBgBgBg
RGRGRGRGR
gBgBgBgBg

*/
/*
// only for testing and distinguishing the output files
if (median)  {
 strcat(out1,".F");
 strcat(out2,".F");
 strcat(out3,".F");
 strcat(out4,".F");
}
if (flag_d)  {
 strcat(out1,".D");
 strcat(out2,".D");
 strcat(out3,".D");
 strcat(out4,".D");
}
*/

 strcat(out1,".R.fits");
 strcat(out2,".G.fits");
 strcat(out3,".B.fits");
 strcat(out4,".L.fits");
 
 channel1 = (float*) calloc(n*m, sizeof(float));   //red
 channel2 = (float*) calloc(n*m, sizeof(float));   //green
 channel3 = (float*) calloc(n*m, sizeof(float));   //blue
sprintf(buffer,"Bayerpattern used: %s",pattern);
//qfits_header_add(header,"HISTORY",buffer,"","");



 // correct Offset &  correct colors
if (!(R==1.0 && G==1.0 && B==1.0 && o==0.0))
{
 printf("subtracts offset:%4.2f\ncolorbalance with sensormultipliers of R:G:B:G = %4.2f:%4.2f:%4.2f:%4.2f\n",o,R,G,B,g);

sprintf(buffer,"Correct offset: %4.2f ",o);
//qfits_header_add(header,"HISTORY",buffer,"","");
sprintf(buffer,"Correct colorbalance R:G:B:G = %4.2f:%4.2f:%4.2f:%4.2f",R,G,B,g);
//qfits_header_add(header,"HISTORY",buffer,"","");


 for (j=0; j<m; j++) { //Zeilen
      for (i=0; i<n; i++) { // Spalten

           if(j==0||j>=m-1||i==0|| i>=n-1) {  // Rand oben unten links rechts
              // don't care here
              k1++;
           }
           else{
           if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // Rotes Feld
              image_in->data[(i+n*j)] = R*(image_in->data[(i+n*j)]-o);
              k1++;
           }
           if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
              image_in->data[(i+n*j)]= G*(image_in->data[(i+n*j)]-o);
              k1++;
           }
          if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
              image_in->data[(i+n*j)]= g*(image_in->data[(i+n*j)]-o);
              k1++;;
          }
           if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
             image_in->data[(i+n*j)]  = B*(image_in->data[(i+n*j)]-o);
              k1++;
           }
          }
      }
  }
}
k1=0;
if (flag_t){
printf("creating testmap.fits ...");
strcpy(out1, "testmap.fits");
//image_out1 = image_new(600,400);
npixels=600*400;
for (j=yoffset; j<(400+yoffset); j++){ //rows
    for(i=xoffset;i<(600+xoffset);i++){ // cols
      if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0) channel1[k1]=50; // rotes Feld
      if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0) channel1[k1]=100; // grunes Feld
      if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1) channel1[k1]=101; // grunes Feld
      if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1) channel1[k1]=75; // blaues Feld
      k1++;
    }

}
//image_save_fits_hdrdump(image_out1, out1, header, BPP_IEEE_FLOAT);
writeoutputfile(out1, fptr , ofptr1, channel1, &npixels, &status);
printf("done!\n");
exit (0);
}


  k1 = 0;  // red
if (flag_q == 0){
//qfits_header_add(header,"HISTORY","Interpolation: bilinear","","");
printf("Using Bilinear\n");
  for (j=0; j<m; j++) { //Zeilen
      for (i=0; i<n; i++) { // Spalten

           if(j==0||j>=m-1||i==0|| i>=n-1) {  // Rand oben unten links rechts
              /*channel1[k1]=image_in->data[(i+n*j)];
              channel2[k1]=image_in->data[(i+n*j)];
              channel3[k1]=image_in->data[(i+n*j)];*/
              if(!xoffset && !yoffset){  // RGGB
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i+1+n*j)];
                    channel3[k1] = image_in->data[(i+1+n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i-n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i-1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i-1-n+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }
              if(xoffset && !yoffset){ //GRBG
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i-1+n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i-n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i+1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i+1-n+n*j)];
                    channel2[k1] = image_in->data[(i+1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }
              if(!xoffset && yoffset){ //GBRG
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i-n+n*j)];
                    channel3[k1] = image_in->data[(i+1-n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i+n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i-1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i-1+n+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }
              if(xoffset && yoffset){ //BGGR
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i-1-n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i+n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i+1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i+1+n+n*j)];
                    channel2[k1] = image_in->data[(i+1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }
              k1++;
           }
           else{
           if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
              channel1[k1] = image_in->data[(i+n*j)];
              channel2[k1]=(image_in->data[(i-1+n*j)]+ image_in->data[(i+1+n*j)]
              +image_in->data[(i-n+n*j)] + image_in->data[(i+n+n*j)]) / 4.0;
              channel3[k1] = (image_in->data[(i-n-1+n*j)] + image_in->data[(i-n+1+n*j)]+
              image_in->data[(i+n-1+n*j)] + image_in->data[(i+n+1+n*j)]) /4.0;
              k1++;
           }
           if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
              channel1[k1] = (image_in->data[(i-n+n*j)]+image_in->data[(i+n+n*j)])/2.0;
              channel2[k1] = image_in->data[(i+n*j)];
              channel3[k1] = (image_in->data[(i-1+n*j)]+image_in->data[(i+1+n*j)])/2.0;
              k1++;
           }
          if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
              channel1[k1] = (image_in->data[(i-1+n*j)]+image_in->data[(i+1+n*j)])/2.0;
              channel2[k1] = image_in->data[(i+n*j)];
              channel3[k1] = (image_in->data[(i-n+n*j)]+image_in->data[(i+n+n*j)])/2.0;
              k1++;;
          }
           if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
              channel1[k1] = (image_in->data[(i-n-1+n*j)] + image_in->data[(i-n+1+n*j)]+
              image_in->data[(i+n-1+n*j)] + image_in->data[(i+n+1+n*j)]) /4.0;
              channel2[k1]=(image_in->data[(i-1+n*j)]+ image_in->data[(i+1+n*j)]
              +image_in->data[(i-n+n*j)] + image_in->data[(i+n+n*j)]) / 4.0;
              channel3[k1] = image_in->data[(i+n*j)];
              k1++;
           }
          }
      }
  }
}

// Gradienten Mode;  refers to: grafics.cs.msa.ru/en/publications/text/gc2004lk.pdf (A. Lukin
// and D. Kubasov)
if (flag_q == 1){
printf("Using Gradient\n");
//qfits_header_add(header,"HISTORY","Interpolation: gradient","","");

  for (j=0; j<m; j++) { //Zeilen
      for (i=0; i<n; i++) { // Spalten

           if(j==0||j>=m-1||i==0|| i>=n-1) {  // Rand oben unten links rechts
              /*channel1[k1]=image_in->data[(i+n*j)];
              channel2[k1]=image_in->data[(i+n*j)];
              channel3[k1]=image_in->data[(i+n*j)];*/
              if(!xoffset && !yoffset){  // RGGB
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i+1+n*j)];
                    channel3[k1] = image_in->data[(i+1+n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i-n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i-1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i-1-n+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }
              if(xoffset && !yoffset){ //GRBG
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i-1+n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i-n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i+1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i+1-n+n*j)];
                    channel2[k1] = image_in->data[(i+1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }
              if(!xoffset && yoffset){ //GBRG
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i-n+n*j)];
                    channel3[k1] = image_in->data[(i+1-n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i+n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i-1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i-1+n+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }
              if(xoffset && yoffset){ //BGGR
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i-1-n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i+n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i+1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i+1+n+n*j)];
                    channel2[k1] = image_in->data[(i+1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }


              k1++;
           }
           else{
	   if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld 
              channel1[k1] = image_in->data[(i+n*j)];
		// Gradients determing for green
	      H=abs(image_in->data[(i-1+n*j)] - image_in->data[(i+1+n*j)]);
	      V=abs(image_in->data[(i-n+n*j)] - image_in->data[(i+n+n*j)]);
              if ( H>=V )
                 channel2[k1]=(image_in->data[(i-n+n*j)] + image_in->data[(i+n+n*j)]) / 2.0;
	      else
	         channel2[k1]=(image_in->data[(i-1+n*j)] + image_in->data[(i+1+n*j)]) / 2.0;
		// Gradients determing for blue H=L (H=high to Right down)
	      H=abs(image_in->data[(i-n-1+n*j)] - image_in->data[(i+1+n*j)]);
	      V=abs(image_in->data[(i+n-1+n*j)] - image_in->data[(i-n+1+n*j)]);
              if (H>=V)
              channel3[k1] = (image_in->data[(i+n-1+n*j)] + image_in->data[(i-n+1+n*j)]) /2.0;
	      else
              channel3[k1]=(image_in->data[(i-n-1+n*j)] + image_in->data[(i+n+1+n*j)]) /2.0;
              k1++;
           }
           if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
              channel1[k1] = (image_in->data[(i-n+n*j)]+image_in->data[(i+n+n*j)])/2.0;
	      channel2[k1] = image_in->data[(i+n*j)];
              channel3[k1] = (image_in->data[(i-1+n*j)]+image_in->data[(i+1+n*j)])/2.0;
              k1++;
           }
          if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
              channel1[k1] = (image_in->data[(i-1+n*j)]+image_in->data[(i+1+n*j)])/2.0;
	      channel2[k1] = image_in->data[(i+n*j)];
              channel3[k1] = (image_in->data[(i-n+n*j)]+image_in->data[(i+n+n*j)])/2.0;
              k1++;;
	  }
           if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
		// Gradients determing for Red H=L (H=high to Right down)
	      H=abs(image_in->data[(i-n-1+n*j)] - image_in->data[(i+1+n*j)]);
	      V=abs(image_in->data[(i+n-1+n*j)] - image_in->data[(i-n+1+n*j)]);
              if (H>=V)
              channel1[k1] = (image_in->data[(i+n-1+n*j)] + image_in->data[(i-n+1+n*j)]) /2.0;
	      else
              channel1[k1]=(image_in->data[(i-n-1+n*j)] + image_in->data[(i+n+1+n*j)]) /2.0;
		// Gradients determing for green
	      H=abs(image_in->data[(i-1+n*j)] - image_in->data[(i+1+n*j)]);
	      V=abs(image_in->data[(i-n+n*j)] - image_in->data[(i+n+n*j)]);
	      if ( H>=V )
                 channel2[k1]=(image_in->data[(i-n+n*j)] + image_in->data[(i+n+n*j)]) / 2.0;
	      else
	         channel2[k1]=(image_in->data[(i-1+n*j)] + image_in->data[(i+1+n*j)]) / 2.0;
              channel3[k1] = image_in->data[(i+n*j)];
              k1++;
           }
          }
      }
  }
}

// PPG, web.cecs.pdx.edu/~cklin/demosaic/  (Chuan-Kai Lin)
if (flag_q == 2){
printf("Using Pixel Grouping\n");
//qfits_header_add(header,"HISTORY","Interpolation: PPG","","");

// Calculation the green values at red and blue pixels
  for (j=0; j<m; j++) { //Zeilen
      for (i=0; i<n; i++) { // Spalten

           //if(j<=2||j>=m-2||i<=2|| i>=n-2) {         // 3 Reihen Rand oben unten links rechts
           if(j<=2||j>=m-3||i<=2|| i>=n-3) {
              channel1[k1]=image_in->data[(i+n*j)];  // alle gleich gemacht mit aktueller Farbe
              channel2[k1]=image_in->data[(i+n*j)];
              channel3[k1]=image_in->data[(i+n*j)];
              

              k1++;
           }
           else{
//Gradients calculation  for green values at red or blue pixels
DN=abs(image_in->data[(i-2*n+n*j)]-image_in->data[(i+n*j)])*2.0+abs(image_in->data[(i-n+n*j)]-image_in->data[(i+n+n*j)]);
DE=abs(image_in->data[(i+n*j)]-image_in->data[(i+2+n*j)])*2.0+abs(image_in->data[(i-1+n*j)]-image_in->data[(i+1+n*j)]);
DW=abs(image_in->data[(i-2+n*j)]-image_in->data[(i+n*j)])*2.0+abs(image_in->data[(i-1+n*j)]-image_in->data[(i+1+n*j)]);
DS=abs(image_in->data[(i+n*j)]-image_in->data[(i+2*n+n*j)])*2.0+abs(image_in->data[(i-n+n*j)]-image_in->data[(i+n+n*j)]);

	  if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld 
              channel1[k1] = image_in->data[(i+n*j)];
	      switch (direction(DN,DE,DW,DS)){
	      case 1: channel2[k1]=(image_in->data[(i-n+n*j)]*3.0+image_in->data[(i+n*j)]
              +image_in->data[(i+n+n*j)] - image_in->data[(i-2*n+n*j)]) / 4.0;break;
	      case 2: channel2[k1]=(image_in->data[(i+1+n*j)]*3.0+image_in->data[(i+n*j)]
              +image_in->data[(i-1+n*j)] - image_in->data[(i+2+n*j)]) / 4.0;break;
	      case 3: channel2[k1]=(image_in->data[(i-1+n*j)]*3.0+image_in->data[(i+n*j)]
              +image_in->data[(i+1+n*j)] - image_in->data[(i-2+n*j)]) / 4.0;break;
	      case 4: channel2[k1]=(image_in->data[(i+n+n*j)]*3.0+image_in->data[(i+n*j)]
              +image_in->data[(i-n+n*j)] - image_in->data[(i+2*n+n*j)]) / 4.0;break;
              }
              //channel3[k1] = (image_in->data[(i-n-1+n*j)] + image_in->data[(i-n+1+n*j)]+
              //image_in->data[(i+n-1+n*j)] + image_in->data[(i+n+1+n*j)]) /4.0;
              k1++;
           }
           if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
              //channel1[k1] = (image_in->data[(i-n-1+n*j)] + image_in->data[(i-n+1+n*j)]+
              //image_in->data[(i+n-1+n*j)] + image_in->data[(i+n+1+n*j)]) /4.0;
	      switch (direction(DN,DE,DW,DS)){
	      case 1: channel2[k1]=(image_in->data[(i-n+n*j)]*3.0+image_in->data[(i+n*j)]
              +image_in->data[(i+n+n*j)] - image_in->data[(i-2*n+n*j)]) / 4.0;break;
	      case 2: channel2[k1]=(image_in->data[(i+1+n*j)]*3.0+image_in->data[(i+n*j)]
              +image_in->data[(i-1+n*j)] - image_in->data[(i+2+n*j)]) / 4.0;break;
	      case 3: channel2[k1]=(image_in->data[(i-1+n*j)]*3.0+image_in->data[(i+n*j)]
              +image_in->data[(i+1+n*j)] - image_in->data[(i-2+n*j)]) / 4.0;break;
	      case 4: channel2[k1]=(image_in->data[(i+n+n*j)]*3.0+image_in->data[(i+n*j)]
              +image_in->data[(i-n+n*j)] - image_in->data[(i+2*n+n*j)]) / 4.0;break;
              }
              channel3[k1] = image_in->data[(i+n*j)];
              k1++;
           }
           if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
              //channel1[k1] = (image_in->data[(i-n+n*j)]+image_in->data[(i+n+n*j)])/2.0;
	      channel2[k1] = image_in->data[(i+n*j)];
              //channel3[k1] = (image_in->data[(i-1+n*j)]+image_in->data[(i+1+n*j)])/2.0;
              k1++;
           }
          if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
              //channel1[k1] = (image_in->data[(i-1+n*j)]+image_in->data[(i+1+n*j)])/2.0;
	      channel2[k1] = image_in->data[(i+n*j)];
              //channel3[k1] = (image_in->data[(i-n+n*j)]+image_in->data[(i+n+n*j)])/2.0;
              k1++;;
	  }
          }
      }
  }

k1 = 0; 
// Calculating blue and red at green pixels
// Calculating blue or red at red or blue pixels
 for (j=0; j<m; j++) { //Zeilen
      for (i=0; i<n; i++) { // Spalten

           if(j<=2||j>=m-3||i<=2|| i>=n-3) {  // drei Reihen Rand oben unten links rechts
              /*channel1[k1]=image_in->data[(i+n*j)];  // alle gleich gemacht mit aktueller Farbe
              channel2[k1]=image_in->data[(i+n*j)];
              channel3[k1]=image_in->data[(i+n*j)];*/
              if(!xoffset && !yoffset){  // RGGB
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i+1+n*j)];
                    channel3[k1] = image_in->data[(i+1+n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i-n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i-1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i-1-n+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }
              if(xoffset && !yoffset){ //GRBG
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i-1+n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i-n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i+1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i+1-n+n*j)];
                    channel2[k1] = image_in->data[(i+1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }
              if(!xoffset && yoffset){ //GBRG
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i-n+n*j)];
                    channel3[k1] = image_in->data[(i+1-n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i+n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i-1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i-1+n+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }
              if(xoffset && yoffset){ //BGGR
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i-1-n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i+n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i+1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i+1+n+n*j)];
                    channel2[k1] = image_in->data[(i+1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }

	
              k1++;
           }
           else{
// diagonal Gradients
dne=abs(image_in->data[(i-n+1+n*j)]-image_in->data[(i+n-1+n*j)])+abs(image_in->data[(i-2*n+2+n*j)]-image_in->data[(i+n*j)])+abs(image_in->data[(i+n*j)]-image_in->data[(i+2*n-2+n*j)])+abs(channel2[(i-n+1+n*j)]-channel2[(i+n*j)])+abs(channel2[(i+n*j)]-channel2[(i+n-1+n*j)]);
dnw=abs(image_in->data[(i-n-1+n*j)]-image_in->data[(i+n+1+n*j)])+abs(image_in->data[(i-2-2*n+n*j)]-image_in->data[(i+n*j)])+abs(image_in->data[(i+n*j)]-image_in->data[(i+2+2*n+n*j)])+abs(channel2[(i-n-1+n*j)]-channel2[(i+n*j)])+abs(channel2[(i+n*j)]-channel2[(i+n+1+n*j)]);
	   if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld 
              if (dne<=dnw)
		//printf("hue_transit Value: %5.1f\n",hue_transit(channel2[i-n+1+n*j],channel2[i+n*j],channel2[i+n-1+n*j], image_in->data[(i-n+1+n*j)], image_in->data[(i+n-1+n*j)]));
              channel3[k1]=hue_transit(channel2[i-n+1+n*j],channel2[i+n*j],channel2[i+n-1+n*j], image_in->data[(i-n+1+n*j)], image_in->data[(i+n-1+n*j)]);
	      else
	      channel3[k1]=hue_transit(channel2[i-n-1+n*j],channel2[i+n*j],channel2[i+n+1+n*j], image_in->data[(i-n-1+n*j)], image_in->data[(i+n+1+n*j)]);
              k1++;
           }
           if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
              channel1[k1] = hue_transit(channel2[(i-n+n*j)],image_in->data[(i+n*j)],channel2[(i+n+n*j)],image_in->data[(i-n+n*j)],image_in->data[(i+n+n*j)]);
              channel3[k1] = hue_transit(channel2[(i-1+n*j)],image_in->data[(i+n*j)],channel2[(i+1+n*j)],image_in->data[(i-1+n*j)],image_in->data[(i+1+n*j)]);
              k1++;
           }
          if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
              channel1[k1] = hue_transit(channel2[(i-1+n*j)],image_in->data[(i+n*j)],channel2[(i+1+n*j)],image_in->data[(i-1+n*j)],image_in->data[(i+1+n*j)]);
              channel3[k1] = hue_transit(channel2[(i-n+n*j)],image_in->data[(i+n*j)],channel2[(i+n+n*j)],image_in->data[(i-n+n*j)],image_in->data[(i+n+n*j)]);
              k1++;;
	  }
           if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
              if (dne<=dnw)
              channel1[k1]=hue_transit(channel2[(i-n+1+n*j)],channel2[(i+n*j)],channel2[(i+n-1+n*j)], image_in->data[(i-n+1+n*j)], image_in->data[(i+n-1+n*j)]);
	      else
	      channel1[k1]=hue_transit(channel2[(i-n-1+n*j)],channel2[(i+n*j)],channel2[(i+n+1+n*j)], image_in->data[(i-n-1+n*j)], image_in->data[(i+n+1+n*j)]);
              k1++;
           }
          }
     }

 }

}

//LAPLACE=WACPI, Lu & Tan, 2003, IEEE
// color filter array demosaicking new method, Wenmiao Lu and Yap-Peng Tan
// IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL 12, NO. 10, OCTOBER 2003
// added in Version 1.7
/*
Interpolation of Bayer mask, 
takes advantages of spectral correlation and selection of suitable schemes
== hybrid demosaicking method
3 steps:
1. interpolation all green value
2. interpolating all red values at blue and all blue values at red pixel
3. interpolating all red and blue values at green pixel

calculates 4 different alpha-weighting factors (a[][])
calculates 4 different direction vectors dir[][]

expects that the local the ratios between red [blue] and green are highly similar

result must be cominance filtered by fmedian3x3 to reduce zipper and color artifacts 
( -m 1 option is rather adviseable)
*/
if (flag_q == 3){
//qfits_header_add(header,"HISTORY","Interpolation: WACPI","","");

float I[6][6], I2[8][8],G[6][6],dir[8][8];
int color;
float a[6][6];
k1=0;
printf("Using Weighted Adaptive Color Plane\n");
// Calculation the green values at red and blue pixels
printf("green ... ");
  for (j=0; j<m; j++) { //rows
      for (i=0; i<n; i++) { // cols

       // 5 Reihen Rand oben unten links rechts
           if(j<=4||j>=m-5||i<=4|| i>=n-5) {   // Rand
              if(!xoffset && !yoffset){  // RGGB
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i+1+n*j)];
                    channel3[k1] = image_in->data[(i+1+n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i-n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i-1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i-1-n+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }
              if(xoffset && !yoffset){ //GRBG
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i-1+n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i-n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i+1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i+1-n+n*j)];
                    channel2[k1] = image_in->data[(i+1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }
              if(!xoffset && yoffset){ //GBRG
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i-n+n*j)];
                    channel3[k1] = image_in->data[(i+1-n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i+n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i+1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i-1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i-1+n+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }
              if(xoffset && yoffset){ //BGGR
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0){ // rotes Feld
                    channel1[k1] = image_in->data[(i+n*j)];
                    channel2[k1] = image_in->data[(i-1+n*j)];
                    channel3[k1] = image_in->data[(i-1-n+n*j)]; 
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0){ // gruenes Feld1, rot oben
                    channel1[k1] = image_in->data[(i+n+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-1+n*j)];
                 }
                 if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1){ // gruenes Feld 2 blau oben
                    channel1[k1] = image_in->data[(i+1+n*j)];
                    channel2[k1] = image_in->data[(i+n*j)];
                    channel3[k1] = image_in->data[(i-n+n*j)];
                 }
                 if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1){ // blaues Feld
                    channel1[k1] = image_in->data[(i+1+n+n*j)];
                    channel2[k1] = image_in->data[(i+1+n*j)];
                    channel3[k1] = image_in->data[(i+n*j)];
                 }
              }

	
              k1++;
           }
           else{

//defining color of current pixel
          color=0;
          if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0) color=RED;
          if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0) color=GREEN1; // gruenes Feld1, rot oben
          if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1) color=GREEN2; // gruenes Feld 2 blau oben
          if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1) color=BLUE;
//estimations  for green values at red or blue pixels
// defining filter array to not get confused
          fill_intensities7x7((i+n*j),n,&image_in->data[0],I2);
// calculate weights 

          a[3][4]=1.0/(1. + abs(I2[5][4]-I2[3][4]) + abs(I2[3][4]-I2[1][4]) + abs(I2[4][4]-I2[2][4]) + abs((I2[4][3]-I2[2][3])/2.)  + abs((I2[4][5]-I2[2][5])/2.));
          a[4][3]=1.0/(1. + abs(I2[4][5]-I2[4][3]) + abs(I2[4][3]-I2[4][1]) + abs(I2[4][4]-I2[4][2]) + abs((I2[3][4]-I2[3][2])/2.)  + abs((I2[5][4]-I2[5][2])/2.));
          a[4][5]=1.0/(1. + abs(I2[4][3]-I2[4][5]) + abs(I2[4][5]-I2[4][7]) + abs(I2[4][4]-I2[4][6]) + abs((I2[3][4]-I2[3][6])/2.)  + abs((I2[5][4]-I2[5][6])/2.));
          a[5][4]=1.0/(1. + abs(I2[3][4]-I2[5][4]) + abs(I2[5][4]-I2[7][4]) + abs(I2[4][4]-I2[6][4]) + abs((I2[4][3]-I2[6][3])/2.)  + abs((I2[4][5]-I2[6][5])/2.));

	  if (color==RED||color==BLUE){ // red or blue field same
              dir[3][4]=I2[3][4]+(I2[4][4] -I2[2][4] )/2.;
              dir[4][3]=I2[4][3]+(I2[4][4] -I2[4][2] )/2.;
              dir[4][5]=I2[4][5]+(I2[4][4] -I2[4][6] )/2.;
              dir[5][4]=I2[5][4]+(I2[4][4] -I2[6][4] )/2.;
              channel2[k1] = (a[3][4]*dir[3][4] + a[4][3]*dir[4][3] + a[4][5]*dir[4][5] +a[5][4]*dir[5][4])/(a[3][4]+a[4][3]+a[4][5]+a[5][4]);
              k1++;
          }
          else{  // Green
              channel2[k1]=image_in->data[(i+n*j)];
              k1++;
          }
    }
  }
}
printf("done\nred and blue at each other ...");
k1=0;
// Calculation the red and blue values at green pixel
  for (j=0; j<m; j++) { //Zeilen
      for (i=0; i<n; i++) { // Spalten

           // 5 Reihen Rand oben unten links rechts
          if(j<=4||j>=m-5||i<=4|| i>=n-5) {      //innerhalb des Randes
          k1++; }
          else
          {
          color=0;
          if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0) color=RED;
          if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0) color=GREEN1; // gruenes Feld1, rot oben
          if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1) color=GREEN2; // gruenes Feld 2 blau oben
          if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1) color=BLUE;
//estimations  for green values at red or blue pixels
// defining filter array to not get confused
          fill_intensities5x5((i+n*j),n,channel2,G);  // already estimated complete green plane
// calculate weights
          a[2][2]=1.0/(1. + abs(G[1][1]-G[2][2]) + abs(G[2][2]-G[3][3]) +  abs((G[2][1]-G[3][2])/2.)  + abs((G[1][2]-G[2][3])/2.));
          a[2][4]=1.0/(1. + abs(G[1][5]-G[2][4]) + abs(G[2][4]-G[3][3]) +  abs((G[1][4]-G[2][3])/2.)  + abs((G[2][5]-G[3][4])/2.));
          a[4][2]=1.0/(1. + abs(G[5][1]-G[4][2]) + abs(G[4][2]-G[3][3]) +  abs((G[4][1]-G[3][2])/2.)  + abs((G[5][2]-G[4][3])/2.));
          a[4][4]=1.0/(1. + abs(G[5][5]-G[4][4]) + abs(G[4][4]-G[3][3]) +  abs((G[5][4]-G[4][3])/2.)  + abs((G[4][5]-G[3][4])/2.));

          //fill_intensities5x5((i+n*j),n,image_in,I);  // I2 only partly used!!
          fill_intensities5x5((i+n*j),n,&image_in->data[0],I);  // I2 only partly used!!
	  if (color==BLUE){ // red at blue field       // I still green plane
              channel3[k1]=image_in->data[(i+n*j)];
              dir[2][2]=I[2][2]+(G[3][3] -G[2][2] );
              dir[2][4]=I[2][4]+(G[3][3] -G[2][4] );
              dir[4][2]=I[4][2]+(G[3][3] -G[4][2] );
              dir[4][4]=I[4][4]+(G[3][3] -G[4][4] );

              channel1[k1] = (a[2][2]*dir[2][2] + a[2][4]*dir[2][4] + a[4][2]*dir[4][2] +a[4][4]*dir[4][4])/(a[2][2]+a[2][4]+a[4][2]+a[4][4]);
              k1++;
          }

	  if (color==RED){ // blue at red field
              channel1[k1]=image_in->data[(i+n*j)];
              dir[2][2]=I[2][2]+(G[3][3] -G[2][2] );
              dir[2][4]=I[2][4]+(G[3][3] -G[2][4] );
              dir[4][2]=I[4][2]+(G[3][3] -G[4][2] );
              dir[4][4]=I[4][4]+(G[3][3] -G[4][4] );

              channel3[k1] = (a[2][2]*dir[2][2] + a[2][4]*dir[2][4] + a[4][2]*dir[4][2] +a[4][4]*dir[4][4])/(a[2][2]+a[2][4]+a[4][2]+a[4][4]);
              k1++;
              }
          if (color==GREEN1||color==GREEN2) k1++;
       }
    }
} 
printf("done\nred and blue at green ... ");
k1=0;
// Calculation the red and blue values at green pixel
  for (j=0; j<m; j++) { //Zeilen
      for (i=0; i<n; i++) { // Spalten
           // 5 Reihen Rand oben unten links rechts
          if(j<=4||j>=m-5||i<=4|| i>=n-5) {      //innerhalb des Randes
          k1++; }
          else
          {
          color=0;
          if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==0) color=RED;
          if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==0) color=GREEN1; // gruenes Feld1, rot oben
          if ((j+yoffset)%2 ==0 && (i+xoffset)%2 ==1) color=GREEN2; // gruenes Feld 2 blau oben
          if ((j+yoffset)%2 ==1 && (i+xoffset)%2 ==1) color=BLUE;
// defining filter array to not get confused
          fill_intensities5x5((i+n*j),n,channel2,G);  // already estimated complete green 
//          fill_intensities5x5((i+n*j),n,channel2,I);  // already estimated complete green channel

          a[2][3]=1.0/(1. + abs(G[3][3]-G[2][3]) + abs(G[2][3]-G[1][3]) +  abs((G[3][2]-G[2][2])/2.)  + abs((G[2][2]-G[1][2])/2.) + abs((G[3][4]-G[2][4])/2.)  + abs((G[2][4]-G[1][4])/2.));
          a[3][2]=1.0/(1. + abs(G[3][3]-G[3][2]) + abs(G[3][2]-G[3][1]) +  abs((G[2][3]-G[2][2])/2.)  + abs((G[2][2]-G[2][1])/2.) + abs((G[4][3]-G[4][2])/2.)  + abs((G[4][2]-G[4][1])/2.)) ;
          a[3][4]=1.0/(1. + abs(G[3][3]-G[3][4]) + abs(G[3][4]-G[3][5]) +  abs((G[2][3]-G[2][4])/2.)  + abs((G[2][4]-G[2][5])/2.) + abs((G[4][3]-G[4][4])/2.)  + abs((G[4][4]-G[4][5])/2.) );
          a[4][3]=1.0/(1. + abs(G[3][3]-G[4][3]) + abs(G[4][3]-G[5][3]) +  abs((G[3][2]-G[4][2])/2.)  + abs((G[4][2]-G[5][2])/2.) + abs((G[3][4]-G[4][4])/2.)  + abs((G[4][4]-G[5][4])/2.));

	  if (color==GREEN1||color==GREEN2){ // red at green field G1=blue row G2=red row
              fill_intensities5x5((i+n*j),n,channel1,I);  // red already partly estimated
              dir[2][3]=I[2][3]+(G[3][3] -G[2][3] );
              dir[3][2]=I[3][2]+(G[3][3] -G[3][2] );
              dir[3][4]=I[3][4]+(G[3][3] -G[3][4] );
              dir[4][3]=I[4][3]+(G[3][3] -G[4][3] );

              channel1[k1] = (a[2][3]*dir[2][3] + a[3][2]*dir[3][2] + a[3][4]*dir[3][4] +a[4][3]*dir[4][3])/(a[2][3]+a[3][2]+a[3][4]+a[4][3]);

              fill_intensities5x5((i+n*j),n,channel3,I);  // blue already partly estimated
              dir[2][3]=I[2][3]+(G[3][3] -G[2][3] );
              dir[3][2]=I[3][2]+(G[3][3] -G[3][2] );
              dir[3][4]=I[3][4]+(G[3][3] -G[3][4] );
              dir[4][3]=I[4][3]+(G[3][3] -G[4][3] );

              channel3[k1] = (a[2][3]*dir[2][3] + a[3][2]*dir[3][2] + a[3][4]*dir[3][4] +a[4][3]*dir[4][3])/(a[2][3]+a[3][2]+a[3][4]+a[4][3]);
              k1++;
          }
          else{
              k1++;
          }
        }
      }
   }
printf("done\n");
}

//DESPECKLE
/*
based on AHD (ADAPTIVE HOMEGENEITY-DIRECTED )interpolation by Hirakawa and Thomas W. Parks
IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL 14, NO.3, MARCH 2005
fast=1 uses special MEDIAN of 3 elements in configurable directions : / \ - |

*/
if (flag_d ){
sprintf(buffer,"Filter: *despeckle* with level %6.2f",level);
//qfits_header_add(header,"HISTORY",buffer,"","");

// hot pixel reduction with special threshold median
int fast=1;  // fast mode: attention forms the PSF in unlinear matter
float M;
printf("removing hot pixel ... ");
float *channel_op,*channel_des, *channel_opA,*channel_opB;
float *channel_op2;
channel_op=  (float*) calloc(m*n, sizeof(float));
channel_des= (float*) calloc(m*n, sizeof(float));
if (fast) channel_op2= (float*) calloc(m*n, sizeof(float));

//  method
  M=1.0;
//Red
   despeckle (m,n,3,level,channel2,channel_des);  // G'
   diff(m,n,channel1,channel_des,channel_op);     //R-G'
   printf(" R-G'");
   if (!fast){
     median3x3(m,n,channel_op);
     sum(m,n,channel_op,channel_des,channel1,M);    //R+G'
   }
   else{
     despeckle (m,n,4,level,channel_op,channel_op2);
     sum(m,n,channel_op2,channel_des,channel1,M);    //R+G'
   }
//Blue
//   despeckle (m,n,3,channel2,channel_opA);       // G'
   diff(m,n,channel3,channel_des,channel_op);    // B-G'
   printf(" B-G'");
   if (!fast){
     median3x3(m,n,channel_op);
     sum(m,n,channel_op,channel_des,channel3,M);   // B+G'
   }
   else{
     despeckle (m,n,4,level,channel_op,channel_op2);
     sum(m,n,channel_op2,channel_des,channel3,M);    //R+G'
   }
//GreenA
channel_opA= (float*) calloc(m*n, sizeof(float));
   despeckle (m,n,3,level,channel1,channel_des);  // R'
   diff(m,n,channel2,channel_des,channel_op);     //G-R'
   printf(" G-R'");
   if (!fast){
     median3x3(m,n,channel_op);
     sum(m,n,channel_op,channel_des,channel_opA,M); //G+R'
   }
   else{
     despeckle (m,n,4,level,channel_op,channel_op2);
     sum(m,n,channel_op2,channel_des,channel_opA,M);    //R+G'
   }

channel_opB= (float*) calloc(m*n, sizeof(float));
   despeckle (m,n,3,level,channel3,channel_des);  // B'
   diff(m,n,channel2,channel_des,channel_op);     // G-B'
   printf(" G-B'");
   if(!fast){
     median3x3(m,n,channel_op);
     sum(m,n,channel_op,channel_des,channel_opB,M); // G+B' 
   }
   else{
     despeckle (m,n,4,level,channel_op,channel_op2);
     sum(m,n,channel_op2,channel_des,channel_opB,M);    //R+G'
   }

M=0.5;
   sum(m,n,channel_opA,channel_opB,channel2,M);

free (channel_opB);
free (channel_opA);
free (channel_op);
if(fast) free (channel_op2);
free (channel_des);

printf(" done\n");

}

//ARTIFACT REDUCTION
/*
should be done after interpolation with WACPI, Lu & Tan
uses a threshold mask for local crominance to detect edges and zippers
the laplacian operater convolves with the interpolated green channel
to not desaturate the colors as simple median does.

*/
if (median>0){
sprintf(buffer,"Filter: artifact reduction %i x",median);
//qfits_header_add(header,"HISTORY",buffer,"","");

// interpolation artifacts reduction ,part of WACPI, Lu and Tan
int i;
float M;
float *channel_op,*channel_opA,*channel_opB,*mask;
channel_op= (float*) calloc(m*n, sizeof(float));
mask= (float*) calloc(m*n, sizeof(float));
for (i=1;i<=median;i++){
printf("\ninterpolating artifacts, pass %i of %i ",i,median);
M=1.0;
//Green
filter_m(m,n,channel2,mask);
channel_opA= (float*) calloc(m*n, sizeof(float));
   diff(m,n,channel1,channel2,channel_op);
   printf(" R'-G");
   fmedian3x3(m,n,channel_op,mask);
   diff(m,n,channel1,channel_op,channel_opA);
channel_opB= (float*) calloc(m*n, sizeof(float));
   diff(m,n,channel3,channel2,channel_op);
   printf(" B'-G");
   fmedian3x3(m,n,channel_op,mask);
   diff(m,n,channel3,channel_op,channel_opB);
M=0.5;
   sum(m,n,channel_opA,channel_opB,channel2,M);
M=1.0;
//Red
   diff(m,n,channel1,channel2,channel_op);
   printf(" R-G'");
   fmedian3x3(m,n,channel_op,mask);
   sum(m,n,channel_op,channel2,channel1,M);
//Blue
   printf(" B-G'");
   diff(m,n,channel3,channel2,channel_op);
   fmedian3x3(m,n,channel_op,mask);
   sum(m,n,channel_op,channel2,channel3,M);
}
free (channel_opA);
free (channel_opB);
free (channel_op);
printf(" done\n");
}

//L CHANNEL
/*
generates a true Luminance channel of the 3 interpolated RGB channels
for experiments with an artifical luminance 
maximum level is reduced to about 400 max, due to a 8 Bit model
should be watched further
*/
if ( flag_l ){
//qfits_header_add(header,"HISTORY","CIE-L as artifical luminance ALUM","","");
printf("generating CIE-L artifical luminance ... ");
// needs colorsystems.h
float P[3],l[3];
L=(float*) calloc(m*n, sizeof(float));
  for (j=0; j<m; j++) {
      for (i=0; i<n; i++) {
         RGBtoXYZ(channel1[i+n*j], channel2[i+n*j], channel3[i+n*j],P);
	 XYZtoLab(P[0],P[1],P[2],l);
         L[i+n*j]=l[0];
      }
  }
printf("done\n");
}

// write the demosaiced images as fits each
 // image_out1 = image_new(n,m);
 // image_out2 = image_new(n,m);
 // image_out3 = image_new(n,m);
//if (flag_l) image_out4 = image_new(n,m);

writeoutputfile(out1, fptr , ofptr1, channel1, &npixels, &status);
writeoutputfile(out2, fptr , ofptr2, channel2, &npixels, &status);
writeoutputfile(out3, fptr , ofptr3, channel3, &npixels, &status);
if (flag_l) writeoutputfile(out4, fptr , ofptr4, L, &npixels, &status);


// free the memory
  free(channel1);
  free(channel2);
  free(channel3);
//  free(L);

 free(iimg1.data);

  exit (0);
}

void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/


    if (status)
    {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
    return;
}

void writeoutputfile(char *name, fitsfile *old , fitsfile *new, float *data, long *w_npixels, int *w_status)
{
int type;

remove(name);
//printf("remove done\n");
*w_status=0;
if (fits_create_file(&new,name, w_status))
   printerror( *w_status );
//printf("Create  new done\n");
if (fits_copy_header(old,new,w_status))
   printerror( *w_status );
//printf("copy header done\n");
//if (fits_create_img(ofptr1, FLOAT_IMG, 2, naxes,w_status))
//   printerror( *w_status ) ;
//printf("create image done\n");
type=FLOAT_IMG;
if (fits_update_key(new,TLONG,"BITPIX",&type,0,w_status))
   printerror( *w_status );
//printf("copy header done\n");


if (fits_write_img(new, TFLOAT, 1, *w_npixels ,data ,w_status))
   printerror( *w_status ) ;
//printf("write image done\n");
if (fits_write_history(new,"demosaicing by fitsiodemosaicbayer" ,w_status))
   printerror( *w_status ) ;
//printf("Keyupdate done\n");
if (fits_close_file(new,w_status))
   printerror( *w_status ) ;
//printf("close done\n");

}
