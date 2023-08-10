!-----------------------------------------------------------------------------------------------------------------------------------
! variable declaration

! all variables are globally declared here except variables required for calling external functions/subrouties/libraries
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
module varinit
	
!	variables to be loaded from input files
!	system properties
!	fluid/solvent
	double precision, parameter :: 			pi = 3.141592653589793
!											density, dynamic viscosity, hard-core radius
	double precision ::						rho, mu, sigma
!	electrolyte
!											number of electroytic species
	integer ::								nspecies
!											electroneutral number density, relative permittivity
	double precision ::						n0, relperm
!											diffusion coefficients, valencies
	double precision, allocatable ::		diffcoeff(:), valence(:)
!	solid
!											substrate rigididity modulus 1, substrate rigididity modulus 2, Hamaker's constant, 
!											sol pressure amplitude, sol pressure phase
	double precision ::						subsmodul1, subsmodul2, capAsfw, capA, phase
!											switch for solid model (0 = linear-elastic solid, 1 =  linear elastic limit of VE model 
!											2 = Standard Linear Solid, 3 = Kelvin-Voigt), what is the set of 
!											elasticity parameters provided (0 = Young's modulus and Poisson's ratio, 1 = Lame's
!											first parameter and shear modulus, 2 = bulk modulus and shear modulus), mode of 
!											viscoelastic material compressibility (1 = sum of LE and Newtonian, 
!											2 = complex Poisson ratio)
	integer ::								model, iscapLame
!											surface potential (is ideally dependent on the electrolyte properties as well)
	double precision, allocatable ::		psis(:)
!	universal constants
!											universal charge, Boltzmann constant, permittivity of free space
	double precision ::						elementaryq, kcapB, permvac
!	geometry and dynamic parameters
!											probe oscillation amplitude, mean separation of probe from origin, 
!											undeformed substrate thickness, probe radius, probe oscillation frequency, temperature
	double precision ::						h0, capD, capL, capR, omega, capT
!											is the sphere approaching/receding or oscillating (0 = oscillating; 1 = approaching; 
!											2 = receding), whether to force harmonic response with single sinusoidal term for the
!											pressure (and resultantly force) and deflection response of the system (works in
!											viscoelastic solution mode only, hence the ending with 've'), whether a spring is
!											attached above the sphere (0 = no, 1 = yes)
	integer :: 								iappr, iforceharmonicve, ispring
!											spring constant 
	double precision :: 					capKspring
!	simulation parameters
!											consider adaptive grid for y, t, tau (0 = no, 1 = yes)
	integer ::								ady, adt, adtau
!											m and n (skew parameters of	original and uniform grids) for y, t, tau
	double precision ::						mady, nady, madt, nadt, madtau, nadtau
!											scheme to be used for two-way coupling solution (-1 = brute-force using envelope iter.,  
!											0 = Newton-Raphson), whether symmetric boundary condition is imposed for deflection at 
!											centerline (0 = false, 1 = true), guess (0 = zero, 1 = both, 2 = only phd, 3 = only l), 
!											start oscillations from the bottom (0 = false 1 = true), whether to consider limiting 
!											expressions for X matrix
	integer ::								scheme, symbc, guess, revtime, considerlim
!											reduce output to file {real-space and transform-space grid, probe profile, constant and 
!											transient parameters, and, force response always written to file}
!											(-1 = write all variables and don't write interface origin variables and 
!											quarter-phase-radial variables in separate files 0 = write all variables 2 = write all
!											interface variables and deformation field variables 3 = write all interface variables 
!											and flow field variables 4 = write only interface variables 5 = write only interface 
!											origin variables and quarter-phase-radial variables), 
!											show progress of simulation in terminal (0 = false 1 = true)
	integer ::								minout, termout
!											error tolerance for iterative solution, maximum error tolerance to allow break if
!											iterative scheme starts diverging, threshold for small amplitude (maximum allowed 
!	 										value of alpha), threshold for one-way FSI- coupling by deflection mode (maximum 
!											allowed theta/epsi0), threshold for one-way FSI- coupling by deflection time-rate mode 
!											relaxation parameter for root-finding scheme (maximum allowed 
!											(theta/(alph*epsi0))*n*max(abs(dndt,dhdt,1))), relaxation of Newton-Raphson solver, 
!											threshold for adhesion (maximum Gamma to not consider adhesion), threshold for push-in 
!											(maximum Gamma to not consider push-in cavity), threshold of applied pressure scale to 
!											compliance ratio for linear strain assumption (relegatable to a posteriori check), 
!											threshold value for smallness of Wommersley number, number of steps to include substrate 
!											velocity in Reynolds equation in brute-force solution, threshold large value of 
!											L/sqrt((D-h0)*R) for substrate to be considered thin, threshold small value of 
!											L/sqrt((D-h0)*R) for substrate to be considered semi-infinite, threshold large value of 
!											nu for substrate material to be considered perfectly compressible, threshold small value 
!											of nu for substrate material to be considered perfectly incompressible, dummy Poisson 
!											ratio to adjust expressions of X matrix and theta for incompressible substrate material, 
!											threshold for minimum ratio of characteristic hydrodynamic pressure to characteristic 
!											total pressure to consider pressure to be hydrodynamic pressure dominated, threshold for
!											condition 2 for reduction of viscoelastic equations to Fourier-space form
	double precision ::						tol, tolbreak, thressmallamp, thresoneway1, thresoneway2, relax, thresadhesion
	double precision :: 					threscavity, threslin, thresstokes
	integer ::								nsubsvel
	double precision :: 					thresthin, thressemi, threspcompr, thresincompr, nudummy, thresdomhd, thresreduce
!											number of nodes in time grid, radial grid, substrate vertical grid, 
!											solvent vertical grid, Hankel grid, Fourier grid, and VE-history grid 
	integer ::								nt, nr, ny, nz, ntau
!											start and end time (non-dimensional, i.e. omega*t), and, infinity values for r-grid, 
!											x-grid (Hankel), and tau-grid (VE-history)
	double precision ::						tstart, tend, rinf, tauinf

!	system variables

!	real-space grid											
!											time grid, radial grid, vertical grid in substrate, vertical grid in solvent
	double precision, allocatable ::		t(:), r(:), y(:)
	double precision, allocatable ::		z(:,:)
!	transformed-space grid
!	 										Hankel grid, Fourier grid, VE-history grid
	double precision, allocatable ::		x(:), fr(:), tau(:)
!	field variables
!											substrate deflection, total pressure
!											linear-elastic solution
	double precision, allocatable ::		l(:), p(:)
!											visco-elastic solution
	complex*16, allocatable ::				lve(:,:), pve(:,:), pmve(:,:), phdve(:,:)
!											hd pressure, EDL disjoining pressure, van der Waals pressure, solvation pressure, 
	double precision, allocatable ::		phd(:), pdl(:), pvdw(:), ps(:)
!											velocity components
	double precision, allocatable ::		vr(:,:), vz(:,:)
!											deformation components
!											linear-elastic solution
	double precision, allocatable ::		ur(:,:), uy(:,:)
!											visco-elastic solution
	complex*16, allocatable ::				urve(:,:,:), uyve(:,:,:)

!	force response
	double precision, allocatable ::		capF(:), capF_d(:)

!	probe profile and spring compression
	double precision, allocatable ::		capH(:), lcapS(:)


!	transformed co-ordinate field variables and support functions (as per requirement)
	double precision, allocatable ::		capX(:)
	complex*16, allocatable ::				capGbar(:), capKbulkbar(:)
	double precision, allocatable ::		lcapH0(:), pcapH0(:)
	complex*16, allocatable ::				lvecapH0(:,:), pvecapH0(:,:)
	complex*16, allocatable ::				lmvecapH0(:), pmvecapH0(:)
	complex*16, allocatable ::				probevel(:)
	complex*16, allocatable ::				phdfreq(:), phdfreqprev(:), lfreq(:), lfreqprev(:)
	double precision, allocatable ::		urcapH1(:,:), uycapH0(:,:)
	complex*16, allocatable ::				urvecapH1(:,:), uyvecapH0(:,:)
	double precision, allocatable ::		capAcoeff(:), capBcoeff(:), capCcoeff(:), capDcoeff(:)
	complex*16, allocatable ::				capAcoeffve(:), capBcoeffve(:), capCcoeffve(:), capDcoeffve(:)
	double precision, allocatable ::		gcapH(:,:)
	complex*16, allocatable ::				gcapHve(:,:)

!	miscellaneous
!	system parameters (to be calculated)
!	constant non-dimensional parameters
	double precision ::						epsi, alpha, beta, theta, thetacapS, delta, zeta, nu, capGamma, capGammacapS
!	transient non-dimensional parameters
	double precision,allocatable ::			h(:), n(:), ncapS(:), d(:), s(:), m(:), w(:)
!	characteristic total pressure, Debye length, substrate rigidity modulus magnitudes, relaxation frequency
	double precision ::						pic, capK, capEy, nupois, capG, lambda, capKbulk, frrelax
!	number of terms in the prony series
	integer ::								nprony
!	moduli amplitudes and relaxation times
	double precision, allocatable ::		capKbulkprony(:), capGprony(:), taukprony(:), taugprony(:)
!	viscoelasticity moduli definition variables
	double precision :: 					capKbulkc, capGc, capKbulkcnd, capGcnd
	double precision, allocatable :: 		capKbulkpronynd(:), capGpronynd(:), kappaprony(:), gammaprony(:), kprony(:), gprony(:)
!	loss and storage modulus from sinusoidal response
	complex*16 ::							capGdlosstorsin, capGlosstorsin
!	scheme and regime switches
!											small-amplitude, simplified viscoelastic formulation, whether material is thin or 
!											semi-infinite or none, whether material is perfectlyincompressible or perfectly 
!											compressible or none, location in the limit matrix
	integer ::								smallamp, reduce, isthin, isincompr, loclimmat
!											one-way coupling
	integer, allocatable ::					oneway(:)
!	last time step and last iteration containers (system variables)
	double precision, allocatable ::		llast(:), llast2(:), lprev(:), phdprev(:)
	double precision :: 					lcapSprev
	double precision, allocatable ::		tmparr1(:), tmparr2(:), tmparr3(:), tmparr4(:)
	double precision, allocatable ::		tmparx1(:), tmparx2(:), tmparx3(:), tmparx4(:)
!	index iterators
	integer ::								itemp, itemp1, iter, itersubsvel, iph, iprony
	integer ::								ir, iy, iz, it, ix, ifr, itau
	integer ::								ir1, iy1, iz1, it1, ix1, ifr1, itau1 
	integer ::								ir2, iy2, iz2, it2, ix2, ifr2, itau2
	integer ::								jr, jy, jz, jt, jx, jfr, jtau
!	temporary value holders
	double precision ::						containers(1:1000)
	double precision ::						atemp(1:20,-2:2,-2:2)
	integer ::								nt1
	double precision ::						off = 0.0
!	simulation variables
	integer ::								nx, nfr
	integer ::								nodeph(1:4)
	double precision ::						tph(1:4), tphd(1:4)
	double precision ::						rmin, rmax, xmin, xmax, tend1, ymin, ymax, zmin, zmax, xinf
	double precision ::						tstep, rstep, frstep, ystep, zstep, tstep1, xstep
	double precision, allocatable ::		capRes, capReslast
	double precision, allocatable ::		dndtbyn(:), dncapSdtbyncapS(:), dhdtbyh(:)
	double precision, allocatable ::		dpdl(:), dpdldl(:), dpvdwdl(:), dpsdl(:)
	double precision, allocatable ::		dgcapHdy(:,:), d2gcapHdy2(:,:)
	complex*16, allocatable ::				dgcapHvedy(:,:), d2gcapHvedy2(:,:)
	double precision, allocatable ::		dphddr(:),d2phddr2(:), dldr(:)
	double precision, allocatable ::		Hnkkrnlf(:,:), Hnkkrnli(:,:), Hnkkrnli1(:,:)
	double precision, allocatable ::		Bslroot(:), Bslroot1(:), r1(:)
	double precision, allocatable :: 		jacobian(:,:), jacobian_in(:,:), jacobian1(:,:), jacobian1_in(:,:)
	double precision, allocatable :: 		residual(:), residual_in(:), residual1(:), residual1_in(:)
	complex*16, allocatable ::				jacobianfr(:,:), residualfr(:), jacobianfr_in(:,:), residualfr_in(:)
	double precision, allocatable :: 		solution(:)
	complex*16, allocatable :: 				solutionfr(:)
	complex*16, allocatable :: 				pmve_holder(:,:), pve_holder(:,:)

!	environment generation variables
	integer*4 :: 							today(3), now(3)
	character (len=8) :: 					datestamp
	character (len=15) :: 					siminstance
	character (len=15) :: 					probtype
	character (len=2) :: 					dirtr
	logical :: 								exists(20)
	character :: 							chartemp
!	external library/function/submodule support variables
	integer, allocatable :: 				ipiv1(:)
	integer, allocatable :: 				ipiv(:)
	integer :: 								info
	integer :: 								mpiswitch=0
	
end module varinit
!-----------------------------------------------------------------------------------------------------------------------------------
