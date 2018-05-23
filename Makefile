#### This is a Makefile for the jet package. 
#### !!!! WARNING: USE ONLY GNU MAKE !!!!
####
#### Use:
####
####    gmake -f Makefile [VPATH=<path>] [PDF=pdflib] <target name>
####
#### where the entries in square brackets are optional. VPATH needs
#### to be specified when the source (fortran) files sit in a directory
#### different from the working one. When PDF=pdflib is used, the package
#### uses the CERN library of parton densities, instead of the library
#### which is provided with this package. The user may need to adjust
#### the link command to the libraries, depending on the implementation
#### on his machine; look for the string "pdflib" in this file to easily
#### find the place where modifications may be necessary. Finally, 
####   <target name> = MAINPH  MAINHD  MAINHD_NUC
#### are the name of the executables relevant to the pointlike component,
#### the hadronic component, and the hadronic component for nuclear
#### collisions. These executables are referred to in the documentation
#### file as POINTLIKE, HADRONIC, and HI respectively. Notice that in
#### the case of hadronic component no electron PDFs is linked here
#### (such as ELPDF_GRV mentioned in the documentation file). The 
#### modification necessary to include it is trivial
####
ifeq ($(shell uname),AIX)
F77=xlf -qextname -qflttrap=overflow:zerodivide:invalid:enable -O3 -qstrict \
#       -qautodbl=dblpad
SYSOBJ=aix.o
#AUTODBL=-qautodbl=dblpad
endif
ifeq ($(shell uname),SunOS)
F77= gfortran -fnonstd
SYSOBJ=sun.o
endif
ifeq ($(shell uname),Linux)
CC=gcc
F77= gfortran -Wall -fno-automatic # -fno-globals  -fdebug-kludge
SYSOBJ=linux.o trapfpe.o
endif
ifeq ($(shell uname),HP-UX)
F77= gfortran -Wall
SYSOBJ=aix.o
endif

FORTHD=hdyjetdiff.for hdyjetcrs.for
FORTHD_NUC=hdyjetdiff_nuc.for hdyjetcrs.for
FORTPH=phyjetdiff.for phyjetcrs.for
FORTUTI=jetint.for jetuti.for
FORTNPDF= jetpdf.for dummy.for
FORTNPDF_NUC= jetpdf.for jetpdf_nuc.for dummy.for
FORTPDFLIB=jetpdflib.for dummy.for
FORTPDFLIB_NUC=jetpdflib_nuc.for

SYSDEP=avi.for vax.for sun.for aix.for linux.for trapfpe.c

OBJBASEHD = $(patsubst %.for,%.o,$(FORTHD)) $(SYSOBJ)
OBJBASEHD_NUC = $(patsubst %.for,%.o,$(FORTHD_NUC)) $(SYSOBJ)
OBJBASEPH = $(patsubst %.for,%.o,$(FORTPH)) $(SYSOBJ)
OBJUTI = $(patsubst %.for,%.o,$(FORTUTI)) $(SYSOBJ)
OBJNPDF= $(patsubst %.for,%.o,$(FORTNPDF))
OBJNPDF_NUC= $(patsubst %.for,%.o,$(FORTNPDF_NUC))
OBJPDFLIB = $(patsubst %.for,%.o,$(FORTPDFLIB))
OBJPDFLIB_NUC = $(patsubst %.for,%.o,$(FORTPDFLIB_NUC))

ifeq ($(PDF),pdflib)
  OBJHD=$(OBJBASEHD) $(OBJUTI) $(OBJPDFLIB)
  OBJHD_NUC=$(OBJBASEHD_NUC) $(OBJUTI) $(OBJPDFLIB_NUC)
  OBJPH=$(OBJBASEPH) $(OBJUTI) $(OBJPDFLIB)
  LIBS= -lc `cernlib pdflib genlib`
  LIBS804= -lc `cernlib pdflib804 genlib`
else
  OBJHD=$(OBJBASEHD) $(OBJUTI) $(OBJNPDF)
  OBJHD_NUC=$(OBJBASEHD_NUC) $(OBJUTI) $(OBJNPDF_NUC)
  OBJPH=$(OBJBASEPH) $(OBJUTI) $(OBJNPDF)
  LIBS=
  LIBS804=
endif

FOR = $(F77) $(DEBUG)

vpath %.for ../
vpath %.c ../
vpath %.h ../

%.o: %.for
	$(FOR) $(DEBUG) $(AUTODBL) -c $^
%.o: %.c
	$(CC) $(DEBUG) -c $^

MAINHD: $(OBJHD) yjetuser.o
	$(F77) $(DEBUG)  -o $@ $^ $(LIBS)

MAINHD_NUC: $(OBJHD_NUC) yjetuser.o
	$(F77) $(DEBUG)  -o $@ $^ $(LIBS804)

MAINPH: $(OBJPH) yjetuser.o
	$(F77) $(DEBUG)  -o $@ $^ $(LIBS)
