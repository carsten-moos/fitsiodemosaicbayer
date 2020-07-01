#ifndef FITSTOOLS_H
#define FITSTOOLS_H
#define true 1;
#define false 0;

typedef int bool;

int FILEMAX = 4096;
double const PI = 3.1415926535;
double const RAD = 3.1415926535/180.;


//****************************************************************************
double get_posangle(double cd11, double cd12, double cd21, double cd22)
{
  double pscale1, pscale2;
  double pa, nfac1, nfac2;

  // the pixel scale
  pscale1 = sqrt(cd11 * cd11 + cd21 * cd21);
  pscale2 = sqrt(cd12 * cd12 + cd22 * cd22);

  // take out the pixel scale
  cd11 /= pscale1;
  cd12 /= pscale2;
  cd21 /= pscale1;
  cd22 /= pscale2;

  // sqrt(CD elements) is very close to one, but not perfectly
  // as coordinate system is not necessarily orthogonal (shear + contraction)
  nfac1 = sqrt(cd11 * cd11 + cd21 * cd21);
  nfac2 = sqrt(cd12 * cd12 + cd22 * cd22);

  // make sure it is one, so that we can do the inverse trig functions
  cd11 /= nfac1;
  cd21 /= nfac1;
  cd12 /= nfac2;
  cd22 /= nfac2;

  // due to flipping the rotation angle is ambiguous.

  // possi bly, the following could be simplified with det(CD), 
  // but at the moment i don't see how to identify the additional 2*PI easily

  // normal
  if      (cd11 <  0 && cd12 <= 0 && cd21 <= 0 && cd22 >  0) pa = acos(-cd11);       //   0 <= phi <  90
  else if (cd11 >= 0 && cd12 <  0 && cd21 <  0 && cd22 <= 0) pa = acos(-cd11);       //  90 <= phi < 180
  else if (cd11 >  0 && cd12 >= 0 && cd21 >= 0 && cd22 <  0) pa = 2.*PI-acos(-cd11); // 180 <= phi < 270
  else if (cd11 <= 0 && cd12 >  0 && cd21 >  0 && cd22 >= 0) pa = 2.*PI-acos(-cd11); // 270 <= phi < 360

  // flipped
  else if (cd11 >= 0 && cd12 >= 0 && cd21 <= 0 && cd22 >= 0) pa = acos(cd11);       //   0 <= phi <  90
  else if (cd11 <  0 && cd12 >  0 && cd21 <  0 && cd22 <  0) pa = acos(cd11);       //  90 <= phi < 180
  else if (cd11 <  0 && cd12 <= 0 && cd21 >= 0 && cd22 <  0) pa = 2.*PI-acos(cd11); // 180 <= phi < 270
  else if (cd11 >= 0 && cd12 <  0 && cd21 >  0 && cd22 >= 0) pa = 2.*PI-acos(cd11); // 270 <= phi < 360
  else {
    // we are very likely close to 0, 90, 180 or 270 degrees, and the CD matrix is slightly non-orthogonal.
    // In this case, lock onto 0, 90, 180 or 270 degrees. Otherwise, exit with an error.
    double cd11cd12 = fabs(cd11/cd12);
    double cd22cd21 = fabs(cd22/cd21);

    if (cd11cd12 > 20. && cd22cd21 > 20.) {
      if (cd11 > 0. && cd22 > 0.) pa = 0.*PI/2.;
      if (cd11 < 0. && cd22 > 0.) pa = 0.*PI/2.;
      if (cd11 > 0. && cd22 < 0.) pa = 2.*PI/2.;
      if (cd11 < 0. && cd22 < 0.) pa = 2.*PI/2.;
    }    

    else if (cd11cd12 < 0.05 && cd22cd21 < 0.05) {
      if (cd12 > 0. && cd21 < 0.) pa = 1.*PI/2.;
      if (cd12 < 0. && cd21 < 0.) pa = 1.*PI/2.;
      if (cd12 > 0. && cd21 > 0.) pa = 3.*PI/2.;
      if (cd12 < 0. && cd21 > 0.) pa = 3.*PI/2.;
    }

    else {
      fprintf(stderr, "ERROR: Could not determine position angle from CD matrix!\n");
      exit (1);
    }
  }

  return (pa);
}


//********************************************************
void matrix_mult(double a11, double a12, double a21, double a22, 
		 double *b11_ptr, double *b12_ptr, double *b21_ptr, 
		 double *b22_ptr)
{
  double c11, c12, c21, c22;

  c11 = a11 * *b11_ptr + a12 * *b21_ptr;
  c12 = a11 * *b12_ptr + a12 * *b22_ptr;
  c21 = a21 * *b11_ptr + a22 * *b21_ptr;
  c22 = a21 * *b12_ptr + a22 * *b22_ptr;

  *b11_ptr = c11;
  *b12_ptr = c12;
  *b21_ptr = c21;
  *b22_ptr = c22;
}


//*****************************************************************
double get_exptime(char *filename)
{
  double exptime;
// cut

  return (exptime);
}


//*****************************************************************
void has_table(char *filename, int extension)
{
// cut
}


//*****************************************************************
void checkfile(char *filename)
{
//cut
}


//*****************************************************************
int compare(const void *a, const void *b)
{
  float a1 = *((float*)a);
  float b1 = *((float*)b);

  if (a1<b1) 
      return -1;
  else if (a1>b1) 
      return 1;
  else 
      return 0;    
}


//*****************************************************************
int compare_double(const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  
  return (*da > *db) - (*da < *db);
}


//*****************************************************************
long max_of_long_array(long *data, long dim)
{
  long i;
  long max;
  
  max = data[0];
  
  for (i = 0; i < dim; i++)  {
    if (data[i] > max) max = data[i];
  }
  
  return (max);
}

//*****************************************************************
int max_of_int_array(int *data, long dim)
{
  long i;
  int max;
  
  max = data[0];
  
  for (i = 0; i < dim; i++)  {
    if (data[i] > max) max = data[i];
  }
  
  return (max);
}

//*****************************************************************

float max_of_float_array(float *data, long dim)
{
  long i;
  float max;
  
  max = data[0];
  
  for (i = 0; i < dim; i++)  {
    if (data[i] > max) max = data[i];
  }
  
  return (max);
}

//*****************************************************************

double max_of_double_array(double *data, long dim)
{
  long i;
  double max;
  
  max = data[0];
  
  for (i = 0; i < dim; i++)  {
    if (data[i] > max) max = data[i];
  }
  
  return (max);
}


//*****************************************************************

long maxindex_of_float_array(float *data, long dim)
{
  long i, maxindex;
  float max;
  
  max = data[0];
  maxindex = 0;
  
  for (i = 0; i < dim; i++)  {
    if (data[i] > max) {
      max = data[i];
      maxindex = i;
    }
  }
  
  return (maxindex);
}


//*****************************************************************

float min_of_float_array(float *data, long dim)
{
  long i;
  float min;
  
  min = data[0];
  
  for (i = 0; i < dim; i++)  {
    if (data[i] < min) min = data[i];
  }
  
  return (min);
}

//*****************************************************************

double min_of_double_array(double *data, long dim)
{
  long i;
  double min;
  
  min = data[0];
  
  for (i = 0; i < dim; i++)  {
    if (data[i] < min) min = data[i];
  }
  
  return (min);
}

//*****************************************************************

float min_of_long_array(long *data, long dim)
{
  long i;
  long min;
  
  min = data[0];
  
  for (i = 0; i < dim; i++)  {
    if (data[i] < min) min = data[i];
  }
  
  return (min);
}

//*****************************************************************

long minindex_of_float_array(float *data, long dim)
{
  long i, minindex;
  float min;
  
  min = data[0];
  minindex = 0;
  
  for (i = 0; i < dim; i++)  {
    if (data[i] < min) {
      min = data[i];
      minindex = i;
    }
  }
  
  return (minindex);
}

//************************************************************
float mean_of_array(float *data, long dim)
{
  long i, n;
  float avg;
  
  avg = 0.;
  n = 0;
  
  for (i = 0; i < dim; i++) {
    avg += data[i];
    n++;
  }
  
  if (n != 0.) {
    return (avg / (float) n);
  }
  else return (0.);
}


// calculate the self-weighted mean
//************************************************************
float weighted_mean_of_array(float *data, long dim)
{
  long i, n;
  float avg, w;
  
  avg = 0.;
  n = 0;
  w = 0.;
  
  for (i = 0; i < dim; i++) {
    avg += data[i]*data[i];
    n++;
    w += data[i];
  }
  
  if (n != 0.) {
    return (avg / w);
  }
  else return (0.);
}


//************************************************************
double mean_of_double_array(double *data, long dim)
{
  long i, n;
  double avg;
  
  avg = 0.;
  n = 0;
  
  for (i = 0; i < dim; i++) {
    avg = avg + data[i];
    n++;
  }
  
  if (n != 0.) {
    return (avg / (float) n);
  }
  else return (0.);
}


//************************************************************
float variance(float *data, long dim)
{
  long i, n;
  float avg, sdev;
  
  n = 0;
  avg = mean_of_array(data, dim);
  
  sdev = 0.;
  for (i = 0; i < dim; i++)  {
    sdev = sdev + (avg - data[i]) * (avg - data[i]);
    n++;
  }
  
  if (n > 1) {
    return (sdev / (float) (n - 1));
  }
  else return (0.);
}


//************************************************************
float rms_of_array(float *data, long dim)
{
  long i, n;
  float avg, sdev;
  
  n = 0;
  avg = mean_of_array(data, dim);
  
  sdev = 0.;
  for (i = 0; i < dim; i++)  {
    sdev = sdev + (avg - data[i]) * (avg - data[i]);
    n++;
  }
  
  if (n > 1) {
    return sqrt(sdev / (float) (n - 1));
  }
  else return (0.);
}


//************************************************************
double rms_of_double_array(double *data, long dim)
{
  long i, n;
  double avg, sdev;
  
  n = 0;
  avg = mean_of_double_array(data, dim);
  
  sdev = 0.;
  for (i = 0; i < dim; i++)  {
    sdev = sdev + (avg - data[i]) * (avg - data[i]);
    n++;
  }
  
  if (n > 1) {
    return sqrt(sdev / (float) (n - 1));
  }
  else return (0.);
}


//************************************************************
float crossvariance(float *data1, float *data2, long dim)
{
  long i, n;
  float mean1, mean2, crossvar;
  
  n = 0;
  mean1 = mean_of_array(data1, dim);
  mean2 = mean_of_array(data2, dim);
  
  crossvar = 0.;
  for (i = 0; i < dim; i++)  {
    crossvar = crossvar + (mean1 - data1[i]) * (mean2 - data2[i]);
    n++;
  }
  
  if (n > 1) {
    return (crossvar / (float) (n - 1));
  }
  else return (0.);
}


//************************************************************
void linfit(float *data, long dim, float *result)
{
  long i;
  float sx, sxy;
  float *abscissa;
  
  // fill in the abscissa
  abscissa = (float*) calloc(dim, sizeof(float));
  for (i=0; i<dim; i++) {
    abscissa[i] = (float) (i+1);
  }    
  
  // do the linfit
  sx  = variance(abscissa, dim);
  sxy = crossvariance(abscissa, data, dim);
  // the slope
  result[0] = sxy / sx;
  result[1] = mean_of_array(data, dim) - 
    result[0] * mean_of_array(abscissa, dim);
}


//************************************************************
void get_zscale(float *data, int n, int m, float *zscale, 
		float shrink, float contrast)
{
  long l, t, j;
  long dim_small, dim, i, idim;
  float *sub, *fit, min, max, *intensity;
  
  // define some odd number
  // then select every l-th pixel in the array
  // that is always smaller than n for n > 9
  // the sqrt ensures that the array is probed in a quasi-random pattern
  
  l = (long) (2/3*n + sqrt(n)); 
  
  // the number of elements we test
  dim = n * m;
  dim_small = dim / l;
  
  sub = (float*) calloc(dim_small, sizeof(float));
  
  // select the small array
  t = 0;
  j = 0;
  for (i=0; i<dim; i++) {
    if (t==0 && j < dim_small) {
      sub[j] = data[i];
      j++;
    }
    if (t<l) t++;
    if (t == l) t = 0;
  }
  
  idim = j - 1;
  
  // sort the array
  qsort(sub, idim, sizeof(float), compare);
  
  intensity = (float*) calloc(idim, sizeof(float));
  for (i=0; i<j-1; i++) {
    intensity[i] = sub[i];
  }
  
  free(sub);
  
  min = intensity[0];
  max = intensity[idim-1];
  
  // the linear fit
  fit = (float*) calloc(2, sizeof(float));
  linfit(intensity, idim, fit);
  // replace the data by the fit
  for (i=0; i<idim; i++) {
    intensity[i] = fit[1] + fit[0] * ((float) (i - idim/2));
  }
  
  zscale[0] = intensity[idim/2] + fit[0] / contrast * (float) (1 - idim/2);
  zscale[1] = intensity[idim/2] + fit[0] / contrast * (float) (idim/2);
  
  // don't let the range get larger / smaller than the actual data extremes
  if (zscale[0] < min) zscale[0] = min;
  if (zscale[1] > max) zscale[1] = max;
  
  // shrink the dynamic range if requested (increase the lower threshold)
  if (shrink != 0.) {
    zscale[0] = zscale[0] + shrink * (zscale[1] - zscale[0]);
  }
  
  free(fit);
  free(intensity);
}


//*****************************************************************

double t_mean_of_array(double *data, long dim, float threshold)
{
  long i, n;
  double avg;
  
  avg = 0.;
  n = 0;
  
  for (i = 0; i < dim; i++)  {
    if (data[i] > threshold) {
      avg = avg + data[i];
      n++;
    }
  }
  
  if (n != 0.) {
    return (avg / (double) n);
  }
  else avg = 0.;
  
  return (avg);
}

//*****************************************************************

double t_meanadd_of_array(double *data, long dim, float threshold, int add)
{
  long i, n;
  double avg;
  
  avg = 0.;
  n = 0;
  
  for (i = 0; i < dim; i++)  {
    if (data[i] > threshold) {
      avg = avg + data[i];
      n++;
    }
  }
  
  if (n != 0.) {
    if (add == 0) return (avg / (double) n);
    else return (avg);
  }
  else return (0.);
}

//*****************************************************************

double t_rms_of_array(double *data, long dim, float threshold)
{
  long i, n;
  double avg, sdev;
  
  avg = t_mean_of_array(data, dim, threshold);
  sdev = 0.;
  n = 0;
  
  for (i = 0; i < dim; i++)  {
    if (data[i] > threshold) {
      sdev = sdev + (avg - data[i]) * (avg - data[i]);
      n++;
    }
  }
  
  if (sdev != 0. && n++ > 0) {
    sdev = sqrt(sdev /( (double) n - 1.));
  }
  else sdev = 0.;
  
  return (sdev);
}

//*****************************************************************
float estimate(double *data, long n, long m, float threshold)
{
  long l, t, j;
  long dim_small, dim, dim90, i;
  float *sub, global_rms;
  double *arr90;
  
  // define some odd number
  // then select every l-th pixel in the array
  // that is always smaller than n for n > 9
  // the sqrt ensures that the array is probed in a quasi-random pattern
  
  l = 2/3*n + (long) sqrt(n); 
  
  // the number of elements we test
  dim = n * m;
  dim_small = dim / l;
  
  sub = (float*) calloc(dim_small, sizeof(float));
  
  // select the small array
  t = 0;
  j = 0;
  for (i=0; i<dim; i++) {
    if (t==0 && j < dim_small) {
      sub[j] = data[i];
      j++;
    }
    if (t<l) t++;
    if (t == l) t = 0;
  }
  
  // sort the array
  qsort(&sub[0], j-1, sizeof(float), compare);
  
  // get the rms from the middle 90%    
  float tmin, tmax;
  long indlow, indhigh;
  
  for (i=0; i<j; i++) {
    printf("%f\n", sub[i]);
  }
  
  indlow = (long) (0.05*j);
  indhigh = (long) (0.95*j);
  tmin = sub[indlow];
  tmax = sub[indhigh];
  printf("%lf %lf %ld %ld %ld\n", tmin, tmax, j, indlow, indhigh);
  dim90 = 0;
  for (j=0; j<m; j++) {
    for (i=0; i<n; i++) {
      if (data[i+n*j] >= tmin && data[i+n*j] <= tmax) dim90++;
    }
  }
  
  arr90 = (double*) calloc(dim90, sizeof(double));
  dim90 = 0;
  for (j=0; j<m; j++) {
    for (i=0; i<n; i++) {
      if (data[i+n*j] >= tmin && data[i+n*j] <= tmax) {
	arr90[dim90] = data[i+n*j];
	dim90++;
      }
    }
  }
  
  free(sub);
  
  global_rms = t_rms_of_array(arr90, dim90, threshold);
  
  free(arr90);
  
  return(global_rms);
}

//***************************************************************
float median(float *data, long n)
{
  float median;

  qsort(data, n, sizeof(float), compare);

  if (n % 2) 
    median = data[n / 2];
  else {
    if (n != 1) 
      median = 0.5 * (data[n / 2] + data[n / 2 - 1]);
    if (n == 1)
      median = data[0];
  }

  return(median);
}

//***************************************************************
float medianconserve(float *data, long n)
{
  long i;
  float median;
  float *data2;

  data2 = (float*) calloc(n, sizeof(float));

  for (i=0; i<n; i++) {
    data2[i] = data[i];
  }

  qsort(data2, n, sizeof(float), compare);

  if (n % 2) 
    median = data2[n / 2];
  else 
    median = 0.5 * (data2[n / 2] + data2[n / 2 - 1]);

  free(data2);

  return(median);
}

//***************************************************************
double median_double(double *data, long n)
{
  double median;

  qsort(data, n, sizeof(double), compare);

  if (n % 2) 
    median = data[n / 2];
  else 
    median = 0.5 * (data[n / 2] + data[n / 2 - 1]);

  return(median);
}

//***************************************************************
float medianfast(float *data, long n, long m)
{
  long l, t, j;
  long dim_small, dim, i;
  float *sub, median;
  
  // define some odd number
  // then select every l-th pixel in the array
  // that is always smaller than n for n > 9
  // the sqrt ensures that the array is probed in a quasi-random pattern
  
  l = 2/3*n + sqrt(n); 
  
  // the number of elements we test
  dim = n * m;
  dim_small = dim / l;
  
  sub = (float*) calloc(dim_small, sizeof(float));
  
  // select the small array
  t = 0;
  j = 0;
  for (i=0; i<dim; i++) {
    if (t==0 && j < dim_small) {
      sub[j] = data[i];
      j++;
    }
    if (t<l) t++;
    if (t == l) t = 0;
  }
  
  // sort the array
  qsort(&sub[0], j-1, sizeof(float), compare);
  
  // get the median
  median = sub[(long)(0.5*j)];
  
  free(sub);
  
  return median;
}

//*******************************************************
char *RemoveChars( char *src, char *key )
{
  char *dest;
  size_t len_src;
  size_t len_key;
  int found;
  int i;
  int j;
  int k;
  
  i = 0;
  j = 0;
  k = 0;
  len_src = 0;
  len_key = 0;
  dest = NULL;
  
  len_src = strlen( src );
  len_key = strlen( key );
  
  dest = (char*) malloc(sizeof(char) * len_src + 1 );
  
  if ( NULL == dest ) {
    printf("Unable to allocate memory\n");
    exit (1);
  }
  
  memset( dest, 0x00, sizeof( char ) * len_src + 1 );
  
  for ( i = 0; i < len_src; i++ ) {
    found = false;
    for ( j = 0; j < len_key; j++ ) {
      if ( src[i] == key[j] ) found = true;
    }
    
    if (found == 0) {
      dest[k] = src[i];
      k++;
    }
  }
  
  return (dest);
}

//************************************************
// remove leading and trailing whitespace
char *trimwhitespace(char *str)
{
  char *end;

  // Trim leading space
  while(isspace(*str)) str++;

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;

  // Write new null terminator
  *(end+1) = 0;

  return str;
}


//************************************************
bool isInt(char* str, int test)
{
  int length = strlen(str), i;
  char sth;

  for (i=0; i<length; i++) {
    sth = str[i];
    // check ASCII codes that can make up integers

    // test for positive integers only!
    if (test > 0) {
      // string starts with a zero
      if (i==0 && sth == 48) return false;
      if (! (sth >= 48 && sth <= 57)) return false;
    }
    // positive or zero integer
    if (test == 0) {
      if (! (sth >= 48 && sth <= 57)) return false;		
    }
    // all integers
    if (test < 0) {
      if (! ((sth >= 48 && sth <= 57) || sth == 45) ) return false;		
    }
  }
  return true;
}


//************************************************
bool isFloat(char* str)
{
  int length = strlen(str), i;
  char sth;
  for (i=0; i<length; i++) {
    sth = str[i];
    // check ASCII codes that can make up float numbers
    if(! ((sth >= 48 && sth <= 57) || sth == 45 || sth == 46) )
      return false;		
  }
  return true;
}

/*
Polygon tester (pnpoly)

Copyright (c) 1970-2003, Wm. Randolph Franklin

Permission is hereby granted, free of charge, to any person obtaining a copy 
of this software and associated documentation files (the "Software"), to deal 
in the Software without restriction, including without limitation the rights 
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions: 

(1) Redistributions of source code must retain the above copyright notice, this 
list of conditions and the following disclaimers.
(2) Redistributions in binary form must reproduce the above copyright notice in 
the documentation and/or other materials provided with the distribution.
(3) The name of W. Randolph Franklin may not be used to endorse or promote 
products derived from this Software without specific prior written permission. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
SOFTWARE. 
*/


int pnpoly(long nvert, float *vertx, float *verty, float testx, float testy)
{
  long i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}


//**********************************************************************
void apply_poly(int *data, long nvert, float *vertx, float *verty, 
		long n, long m, int maskvalue)
{
  int polytest;
  long i, j;
  float x, y;

  // apply the polygon mask
  for (j=0; j<m; j++) {
    // we have to add +1 to compensate for arrays in C starting with 0
    y = (float) j + 1;
    for (i=0; i<n; i++) {
      x = (float) i + 1;
      // test if the pixel is inside or outside the polygon
      polytest = pnpoly(nvert, vertx, verty, x, y);

      // mark the pixels inside the polygon
      if (polytest == 1) data[i+n*j] += maskvalue;
    }
  }
}


//**********************************************************************
void apply_polycircle(int *data, float x, float y, float r, 
		      long n, long m, int maskvalue)
{
  long i, j;
  float ii, jj, d;

  // apply the polygon mask
  for (j=0; j<m; j++) {
    jj = (float) j;
    for (i=0; i<n; i++) {
      ii = (float) i;
      // mark the pixels if inside the circle
      d = (ii-x) * (ii-x) + (jj-y) * (jj-y);
      if (d <= r*r) data[i+n*j] += maskvalue;
    }
  }
}


//*****************************************************************
// the number of vertices in a polygon
//*****************************************************************
long polygon_nvert(char *polystring)
{
  size_t ii, length;
  long nvert;
  char polystring_vert[FILEMAX];

  length = strlen(polystring);
  strcpy(polystring_vert,"");

  // remove the 'polygon(' string and the closing ')'
  for (ii=8; ii<length-1 && polystring[ii] != '\0'; ii++)
    polystring_vert[ii-8] = polystring[ii];
  for ( ; ii < length-1; ii++)
    polystring_vert[ii] = '\0';

  length = ii-8;

  // how many vertices do we have
  nvert = 0;
  for (ii=0; ii<length; ii++) {
    if (polystring_vert[ii] == ',') nvert++;
  }
  
  return (nvert + 1) / 2;
}


//****************************************************************
// split a ds9 region POLYGON string into x- and y vertice arrays
//****************************************************************
void polygon2vertices(char *polystring, float *vertx, float *verty)
{
  size_t ii;
  long j, k, length;
  char *polystring_vert, *current;

  length = strlen(polystring);

  polystring_vert = polystring;

  // remove the 'polygon(' string
  for (ii=8; ii<length-1 && polystring[ii] != '\0'; ii++)
    polystring_vert[ii-8] = polystring[ii];
  // cut off the closing ')'
  polystring_vert[length-9] = '\0';

  length = ii-8;

  // split the polygon string into its vertices
  j = 0;
  k = 1;
  current = strtok(polystring_vert, ",");
  vertx[j] = atof(current);
  while ( (current = strtok(NULL, ",")) != NULL) {
    if (k%2 == 0) vertx[j] = atof(current);
    else {
      verty[j] = atof(current);
      j++;
    }
    k++;
  }
}


//****************************************************************
// split a ds9 region POLYGON string into x- and y vertice arrays
//****************************************************************
void polygon2circle(char *polystring, float *x_ptr, float *y_ptr, float *r_ptr)
{
  size_t ii;
  long j, length;
  char *polystring_vert, *current;

  length = strlen(polystring);

  polystring_vert = polystring;

  // remove the 'circle(' string
  for (ii=7; ii<length-1 && polystring[ii] != '\0'; ii++)
    polystring_vert[ii-7] = polystring[ii];
  // cut off the closing ')'
  polystring_vert[length-8] = '\0';

  length = ii-7;

  // extract the middle point and the radius
  j = 0;
  current = strtok(polystring_vert, ",");
  *x_ptr = atof(current);
  while ( (current = strtok(NULL, ",")) != NULL) {
    if (j==0) *y_ptr = atof(current);
    if (j==1) *r_ptr = atof(current);
    j++;
  }
}


//****************************************************************
// mask a ds9 polygon region
//****************************************************************
void mask_ds9_polygon(int *mask, char *polystring, long n, long m)
{
  size_t ii, length;
  long j, k, nvert;
  char *polystring_vert, *current;
  float *vertx, *verty;

  length = strlen(polystring);
  polystring_vert = polystring;

  // remove the 'polygon(' string
  for (ii=8; ii<length-1 && polystring[ii] != '\0'; ii++)
    polystring_vert[ii-8] = polystring[ii];
  // cut off the closing ')'
  polystring_vert[length-9] = '\0';

  length = ii-8;

  // how many vertices do we have
  nvert = 0;
  for (ii=0; ii<length; ii++) {
    if (polystring_vert[ii] == ',') nvert++;
  }
  nvert = (nvert + 1) / 2;

  vertx = (float*) calloc(nvert, sizeof(float));
  verty = (float*) calloc(nvert, sizeof(float));

  // split the polygon string into its vertices
  j = 0;
  k = 1;
  current = strtok(polystring_vert, ",");
  vertx[j] = atof(current);
  while ( (current = strtok(NULL, ",")) != NULL) {
    if (k%2 == 0) vertx[j] = atof(current);
    else {
      verty[j] = atof(current);
      j++;
    }
    k++;
  }

  apply_poly(mask, nvert, vertx, verty, n, m, 1);
  
  free(vertx);
  free(verty);
}

//****************************************************************
// mask a ds9 circle region
//****************************************************************
void mask_ds9_circle(int *mask, char *polystring, long n, long m)
{
  float xcen, ycen, radius;
  
  xcen = 0.;
  ycen = 0.;
  radius = 0.;

  // extract the circle midpoint and radius
  polygon2circle(polystring, &xcen, &ycen, &radius);

  // apply the mask
  apply_polycircle(mask, xcen, ycen, radius, n, m, 1);
}

//****************************************************************
// Translate a ds9 polygon mask into a mask image;
// The mask image is set to 1 whereever a pixel is found inside 
// ds9-type polygon or circle; unless there is a header line:
//      # Sense: out
// In this case the mask will be inverted
//****************************************************************
void make_ds9_maskimage(int *mask, char *ds9reg, long n, long m, 
			char *area_arg, char *combine_arg)
{
  FILE *_file_;
  size_t ii;
  long length, i;
  char dummy[FILEMAX], line_copy[FILEMAX], line[FILEMAX], 
    area[FILEMAX], combine[FILEMAX];

  // if the mask orientation is determined from the ds9 region file:
  if (strcmp(area_arg,"") == 0) strcpy(area, "in");
  else strcpy(area, area_arg);

  // if the combine type is determined from the ds9 region file:
  if (strcmp(combine_arg,"") == 0) strcpy(combine, "or");
  else strcpy(combine, combine_arg);

  // read the polygon
  if ( (_file_ = fopen(ds9reg, "r")) == NULL) {
    printf("\tError: Could not read from %s!\n", ds9reg);
    exit (1);
  }
  while ( fgets(dummy, FILEMAX, _file_) != NULL) {
    strcpy(line, trimwhitespace(dummy));

    // obtain the 'sense', i.e. do we keep pixels inside or outside
    if (strncmp(line, "# Sense:", 8) == 0 && strcmp(area_arg,"") == 0) {
      strcpy(line_copy,line);

      length = strlen(line_copy);
      // remove the '# Sense:' string
      for (ii=9; ii<length; ii++)
	area[ii-9] = line_copy[ii];
    }

    // obtain the 'combine type', i.e. do we keep pixels in all regions or only in the overlap
    if (strncmp(line, "# Combine:", 10) == 0 && strcmp(combine_arg,"") == 0) {
      strcpy(line_copy,line);

      length = strlen(line_copy);
      // remove the '# Combine:' string
      for (ii=11; ii<length; ii++)
	combine[ii-11] = line_copy[ii];
    }

    // mask a ds9 polygon
    if (strncmp(line, "polygon(", 8) == 0 || 
	strncmp(line, "POLYGON(", 8) == 0) {
      strcpy(line_copy,line); 
      mask_ds9_polygon(mask, line_copy, n, m);
    }

    // mask a ds9 circle
    if (strncmp(line, "circle(", 7) == 0 || 
	strncmp(line, "CIRCLE(", 7) == 0) {
      strcpy(line_copy,line);
      mask_ds9_circle(mask, line_copy, n, m);
    }
  }
  fclose(_file_);

  // if requested, keep only pixels common to all ds9 region elements
  if (strcmp(combine,"and") == 0) {
    float maskmax = max_of_int_array(mask, n*m);
    for (i=0; i<n*m; i++) {
      if (mask[i] < maskmax) mask[i] = 0;
      else mask[i] = 1;
    }
  }

  // invert the mask if area = "out"
  if (strcmp(area,"out") == 0) {
    for (i=0; i<n*m; i++) {
      if (mask[i] > 0) mask[i] = 0;
      else mask[i] = 1;
    }
  }
}

#endif
