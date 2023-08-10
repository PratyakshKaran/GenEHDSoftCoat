# Makefile for compiling dlss numerical solution
#-----------------------------------------------------------------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------------------------------------------------------------
# Setting compiler and linking libraries
compile 	= 	gfortran
link 		= 	-llapack -lf77blas -lcblas -latlas -lfftw3
#-----------------------------------------------------------------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------------------------------------------------------------
# Declaring object and module files
object 			=	datain.o allocinit.o meshgen.o march.o finout.o constparam.o
subobject		=	transparam.o linelast.o viscelast.o instout.o
module 			= 	varinit.o derivative.o
#-----------------------------------------------------------------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------------------------------------------------------------
# Defining compilation for wrapper
ttsi.out:	$(libdmumps) $(object) $(subobject) $(module)
		$(compile) -o ttsi.out ttsi.f90 $(object) $(module) $(subobject) $(link)
#-----------------------------------------------------------------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------------------------------------------------------------
# Defining compilation for individual subroutines
derivative.o:			derivative.f90
						$(compile) -c -o derivative.o derivative.f90
varinit.o: 				varinit.f90
						$(compile) -c -o varinit.o varinit.f90
datain.o:				datain.f90 $(module)
						$(compile) -c -o datain.o datain.f90
allocinit.o:			allocinit.f90 $(module)
						$(compile) -c -o allocinit.o allocinit.f90
constparam.o:			constparam.f90 $(module)
						$(compile) -c -o constparam.o constparam.f90
meshgen.o:				meshgen.f90 $(module)
						$(compile) -c -o meshgen.o meshgen.f90
march.o:				march.f90 $(module) $(subobject)
						$(compile) -c -o march.o march.f90
transparam.o:			transparam.f90 $(module)
						$(compile) -c -o transparam.o transparam.f90
linelast.o:				linelast.f90 $(module)
						$(compile) -c -o linelast.o linelast.f90 $(link)
viscelast.o:			viscelast.f90 $(module)
						$(compile) -c -o viscelast.o viscelast.f90 $(link)
instout.o:				instout.f90 $(module)
						$(compile) -c -o instout.o instout.f90
finout.o:				finout.f90 $(module)
						$(compile) -c -o finout.o finout.f90				
#-----------------------------------------------------------------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------------------------------------------------------------
# Cleaning operation
clean:
	rm -f *.o *.mod core *.out
	rm -f *.ini diss* Wxyz* Vk* esp* *.stackdump
	rm -f fort.* *.dat
#-----------------------------------------------------------------------------------------------------------------------------------
