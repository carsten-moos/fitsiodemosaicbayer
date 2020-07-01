/*
C. Moos

Revise:  1  11.4.2012
cmoos@gmx.de
transform RGB to XYZ to L*a*b to XYZ to RGB
colorsystems.h
changes:
*/

#include <math.h>
#include "colorsystems.h"

void RGBtoXYZ(float R, float G, float B, float P[3]){
float X,Y,Z;

float T[3][3];

T[0][0]=0.412453;
T[0][1]=0.357580;
T[0][2]=0.180423;

T[1][0]=0.212671;
T[1][1]=0.715160;
T[1][2]=0.072169;

T[2][0]=0.019334;
T[2][1]=0.119193;
T[2][2]=0.950227;

X= T[0][0]*R + T[0][1]*G + T[0][2]*B;
Y= T[1][0]*R + T[1][1]*G + T[1][2]*B;
Z= T[2][0]*R + T[2][1]*G + T[2][2]*B;

P[0]=X;
P[1]=Y;
P[2]=Z;

}

void XYZtoLab(float X, float Y, float Z, float P[3]){
float L,a=0,b=0;

if (X/Xn > LAB_Segment_limit && Y/Yn > LAB_Segment_limit && Z/Zn > LAB_Segment_limit  )
   {
     L = 116 * pow((Y/Yn),1.0/3) - 16;
   }
   else
   {
     L = 903.3 * Y/Yn ;
   }

if (X/Xn > LAB_Segment_limit && Y/Yn > LAB_Segment_limit )
   {
      a = 500 * ( pow(X/Xn,1.0/3) - pow(Y/Yn,1.0/3) );
   }
if (X/Xn <= LAB_Segment_limit && Y/Yn > LAB_Segment_limit )
   {
     a = 500 * ( 7.787 * X/Xn + 16/116 - pow(Y/Yn,1.0/3) ) ;
   }
if (X/Xn > LAB_Segment_limit && Y/Yn <= LAB_Segment_limit )
   {
     a = 500 * ( pow(X/Xn,1.0/3) - 7.787 * Y/Yn + 16/116  ) ;
   }
if (X/Xn <= LAB_Segment_limit && Y/Yn <= LAB_Segment_limit )
   {
     a = 500 * ( 7.787 * Y/Yn - 7.787 * Y/Yn  ) ;
   }

if (Y/Yn > LAB_Segment_limit && Z/Zn > LAB_Segment_limit )
   {
      b = 200 * ( pow(Y/Yn,1.0/3) - pow(Z/Zn,1.0/3) );
   }
if (Y/Yn <= LAB_Segment_limit && Z/Zn > LAB_Segment_limit )
   {
     b = 200 * ( 7.787 * Y/Yn + 16/116 - pow(Z/Zn,1.0/3) ) ;
   }
if (Y/Yn > LAB_Segment_limit && Z/Zn <= LAB_Segment_limit )
   {
     b = 200 * ( pow(Y/Yn,1.0/3) - 7.787 * Z/Zn + 16/116  ) ;
   }
if (Y/Yn <= LAB_Segment_limit && Z/Zn <= LAB_Segment_limit )
   {
     b = 200 * ( 7.787 * Y/Yn - 7.787 * Z/Zn  ) ;
   }

P[0]=L;
P[1]=a;
P[2]=b;
}
void LabtoXYZ(float L, float a, float b, float P[3]){
float X,Y,Z,K;

K= (L + 16) / 116;
X = Xn * pow(( K + a / 500 ), 3);
Y = Yn * pow(K, 3);
Z = Zn * pow(( K - b / 200 ), 3);

P[0]=X;
P[1]=Y;
P[2]=Z;

}
void XYZtoRGB(float X, float Y, float Z, float P[3]){
float R,G,B;

float T[3][3];

T[0][0]=3.240479;
T[0][1]=-1.537150;
T[0][2]=-0.498535;

T[1][0]=-0.969256;
T[1][1]=1.875992;
T[1][2]=0.041556;

T[2][0]=0.055648;
T[2][1]=-0.204043;
T[2][2]=1.057311;

R= T[0][0]*X + T[0][1]*Y + T[0][2]*Z;
G= T[1][0]*X + T[1][1]*Y + T[1][2]*Z;
B= T[2][0]*X + T[2][1]*Y + T[2][2]*Z;


P[0]=R;
P[1]=G;
P[2]=B;
}
