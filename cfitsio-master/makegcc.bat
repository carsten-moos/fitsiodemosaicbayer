rem:  this batch file builds the cfitsio library 
rem:  author Carsten Moos
rem:  using gcc with DEVcpp 4.9xx or mingw
rem:
gcc -c buffers.c
gcc -c cfileio.c
gcc -c checksum.c
gcc -c compress.c
gcc -c drvrfile.c
gcc -c drvrmem.c
gcc -c drvrnet.c
gcc -c drvrsmem.c 
gcc -c drvrgsiftp.c
gcc -c editcol.c
gcc -c edithdu.c
gcc -c eval_l.c
gcc -c eval_y.c
gcc -c eval_f.c
gcc -c fitscore.c
gcc -c getcol.c
gcc -c getcolb.c
gcc -c getcolsb.c
gcc -c getcoli.c
gcc -c getcolj.c
gcc -c getcolui.c
gcc -c getcoluj.c
gcc -c getcoluk.c
gcc -c getcolk.c
gcc -c getcole.c
gcc -c getcold.c
gcc -c getcoll.c
gcc -c getcols.c
gcc -c getkey.c
gcc -c group.c
gcc -c grparser.c
gcc -c histo.c
gcc -c iraffits.c
gcc -c modkey.c
gcc -c putcol.c
gcc -c putcolb.c
gcc -c putcolsb.c
gcc -c putcoli.c
gcc -c putcolj.c
gcc -c putcolui.c
gcc -c putcoluj.c
gcc -c putcoluk.c
gcc -c putcolk.c
gcc -c putcole.c
gcc -c putcold.c
gcc -c putcols.c
gcc -c putcoll.c
gcc -c putcolu.c
gcc -c putkey.c
gcc -c region.c
gcc -c scalnull.c
gcc -c swapproc.c
gcc -c wcsutil.c
gcc -c wcssub.c
gcc -c imcompress.c
gcc -c quantize.c
gcc -c ricecomp.c
gcc -c pliocomp.c
gcc -c fits_hcompress.c
gcc -c fits_hdecompress.c
gcc -c f77_wrap1.c
gcc -c f77_wrap2.c
gcc -c f77_wrap3.c
gcc -c f77_wrap4.c
del libcfitsio.a
pause
ar r libcfitsio.a buffers.o cfileio.o checksum.o compress.o drvrfile.o drvrmem.o drvrnet.o drvrsmem.o drvrgsiftp.o editcol.o edithdu.o eval_l.o eval_y.o eval_f.o fitscore.o getcol.o getcolb.o getcold.o getcole.o getcoli.o getcolj.o getcolk.o getcoll.o getcols.o getcolsb.o getcoluk.o getcolui.o getcoluj.o getkey.o group.o grparser.o histo.o iraffits.o modkey.o putcol.o putcolb.o putcold.o putcole.o putcoli.o putcolj.o putcolk.o putcoluk.o putcoll.o putcols.o putcolsb.o putcolu.o putcolui.o putcoluj.o putkey.o region.o scalnull.o swapproc.o wcssub.o wcsutil.o imcompress.o quantize.o ricecomp.o pliocomp.o fits_hcompress.o fits_hdecompress.o f77_wrap1.o f77_wrap2.o f77_wrap3.o f77_wrap4.o
ranlib libcfitsio.a
pause
gcc -o testprog.exe testprog.c -L. -lcfitsio
rem gcc -f cookbook.c cfitsio.lib

