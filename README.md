# fitsiodemosaicbayer

Demosaicing of a fitsfile for astronomical data.
fitsiodemosaicbayer is a command line tool for the MS console. 

  :red_circle: :white_check_mark: :red_circle: :white_check_mark: :red_circle: :white_check_mark:
  
  :white_check_mark:  :large_blue_diamond: :white_check_mark:  :large_blue_diamond: :white_check_mark:  :large_blue_diamond:
# Releases

Download latest binary releases for MS Windows on x86_64 bit processors

Here are some sample datasets to play with:fitmosaic-canon.fit


# References

the core of demosaicing fitsfiles is used in these projects:

theli V2,  astronomical data pipeline,      see https://www.astro.uni-bonn.de/theli/

theli V3,  astronomical data pipeline,      see https://github.com/schirmermischa/THELI

nightlight, astronomical image processing,  see https://github.com/mlnoga/nightlight


# Capabilities

    Reads FITS files, demosaics with 4 qualities ( Bilinear, gradient, PPG and WACP )
    Stores the 3 colorchannels to a single FITS files,

# Limitations

    Does not support RAW file input from regular digital cameras, only FITS
    
    Not all switches can be combined and may lead to trash
    
    uses libcfitsio.a (Version 3.24) precompiled at mingw gcc @ windows 10, 64, the library may be outdated or not suitable for every windows version.
    

# Usage via commandline

Usage of fitsiodemosaicbayer is only for console commandline.

          Version 1.7 Rev 4 (2012-04-21) fitsio

  Author: Carsten Moos

  USAGE:  fitsiodemosaicbayer
  
           -i input_image
           
           -p pattern (available: RGGB, GRBG, GBRG, BGGR)
           
           -q 0    for bilinear Interpolation
           -q 1    for gradient Interpolation
           -q 2    for PPG Interpolation
           -q 3    for WACP Interpolation
           -c 2.21 1.00 1.51 1.006 for bayerfilter RGBG colorbalance
           -o 256  for Offsetcorrection before RGBG colorbalance
           -m 2 for 1-5 times artifact filter 3x3 of crominannce
           -d 1000 for despeckle above level 1000  filter 1x3 of crominannce
           -l with export of Luminance CIEL*a*b
           -t with generates RAW testimage: R=50 G=101[2] B=75
           
  PURPOSE: Demosaics a bayer matrix osc fits image into three monochrome fits images for RGB.
