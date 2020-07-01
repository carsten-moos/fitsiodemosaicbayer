#ifndef _COLORSYSTEMS_H
#define _COLORSYSTEMS_H
/*
C. Moos
Template: 
Revise:  1  11.4.2012
cmoos@gmx.de
transform RGB to XYZ to L*a*b to XYZ to RGB
colorsystems.h
changes:
*/

//#define D65[3]={95.05,100,108.9}
#define Xn 95.05
#define Yn 100.0
#define Zn 108.9

#define LAB_Segment_limit 0.008856

void RGBtoXYZ(float R, float G, float B, float P[3]);
void XYZtoLab(float X, float Y, float Z, float P[3]);
void LabtoXYZ(float L, float a, float b, float P[3]);
void XYZtoRGB(float X, float Y, float Z, float P[3]);

#endif
