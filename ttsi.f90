!-----------------------------------------------------------------------------------------------------------------------------------
! wrapper - 	analytical solution of oscillatory loading of a spherical-profile probe over a deformable substrate coating
!
!				employs Hankel-transformation and Fourier-Transformation solutions to solve for the axisymmetric vertical loading on
!				substrate of arbitrary thickness (which can be extended to semi-infinite) solving for the deformation employing the 
!				Cauchy-Green stress tensor with small-strain models of compressible linear-elastic and visco-elastic constitutive 
!				models; the latter is specified in a hereditarty integral formulation, and, its solution methodology poses the 
!				restrictions of one-way pressure-deflection coupling and constant characteristic system scales, implying, 
!				small-oscillation-amplitude, small-oscillation-frequency and small-substrate-deflection 
!				Hankel transformation (forward and reverse) kernel is based on the Baddour and Chouinard, JOSAA, 2015, Fourier
!				transformation uses the FFTW3 subroutine, LAPACK's DGESV subroutine used for solution of sets of linear algebraic
!				equations wherever needed

! author -		pratyaksh karan (pratyakshkaran@iitkgp.ac.in)
! date -        2019-05-24
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
program ttsi

!	variable declaration
	use varinit
	use derivative
	implicit none

!	calling subroutines

!	importing input parameters from in-files
	call datain
!	calculating pertinent constant parameter values
	call constparam
!	allocating size for and initializing matrices and arrays
	call allocinit
!	generating mesh
	call meshgen
!	time-marching
	call march
!	writing output data
	call finout

end program ttsi
!-----------------------------------------------------------------------------------------------------------------------------------
