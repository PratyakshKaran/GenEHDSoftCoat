!-----------------------------------------------------------------------------------------------------------------------------------
! constant parameter calculation subroutine

! obtains constant parameter values (non-dimensional parameters and characteristic values) and passes warnings and cautions
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
subroutine constparam

	use varinit
	implicit none

!	opening logfiles
	open(			unit=93,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/mumps.log', &
					status='new',action='write')
	open(			unit=94,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/main.log', &
					status='new',action='write')
	open(			unit=95,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/caution.log', &
					status='new',action='write')
	open(			unit=96,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/fourier.log', &
					status='new',action='write')
	open(			unit=97,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/hankel.log', &
					status='new',action='write')

!	obtaining the three combinations of elasticity moduli (Ey and nu, and, lambda and G, and, K and G)
	if (iscapLame .eq. 0) then
		capEy =							subsmodul1 	
		nupois =						subsmodul2
		lambda = 						(capEy*nupois)/((1.0+nupois)*(1.0-2.0*nupois))
		capG = 							capEy/(2.0*(1.0+nupois))
		capKbulk = 						capEy/(3.0*(1.0-2.0*nupois))
	elseif (iscapLame .eq. 1) then
		lambda =						subsmodul1 	
		capG =							subsmodul2	
		capEy = 						(capG*(3.0*lambda+2*capG))/(lambda+capG)
		nupois =						lambda/(2.0*(lambda+capG))
		capKbulk = 						lambda+((2.0*capG)/3.0)
	else
		capKbulk =						subsmodul1 	
		capG =							subsmodul2	
		capEy = 						(9.0*capKbulk*capG)/(3.0*capKbulk+capG)
		nupois =						(3.0*capKbulk-2*capG)/(2.0*(3.0*capKbulk+capG))
		lambda = 						capKbulk-((2.0*capG)/3.0)
	end if
!	obtaining Lame parameters and setting compressibility and thinness switches
	if (considerlim .eq. 1) then
!		limits-considered mode
!		determining the limits of incompressibility and thinness
		containers(3) = 		capL/sqrt(capR*(capD-h0))
		if (nupois .ge. thresincompr) then
			isincompr = 				1
		elseif (nupois .le. threspcompr) then
			isincompr = 				-1
		else
			isincompr = 				0
		end if
		if (containers(3) .le. thresthin) then
			isthin = 					1
		elseif (containers(3) .ge. thressemi) then
			isthin = 					-1
		end if
!		determining where in the limiting matrix (matrix for thinness and compressiblity limiting expressions of X matrix and theta) 
!		does the system belong
		if 		(	(isthin .eq. 1	) .and. (isincompr .eq. -1	)	) then
			loclimmat =			1
		elseif 	(	(isthin .eq. 0	) .and. (isincompr .eq. -1	)	) then
			loclimmat =			2
		elseif 	(	(isthin .eq. -1) .and. (isincompr .eq. -1	)	) then
			loclimmat =			3
		elseif 	(	(isthin .eq. 1	) .and. (isincompr .eq. 0	)	) then
			loclimmat =			4
		elseif 	(	(isthin .eq. 0	) .and. (isincompr .eq. 0	)	) then
			loclimmat =			5
		elseif 	(	(isthin .eq. -1) .and. (isincompr .eq. 0	)	) then
			loclimmat =			6
		elseif 	(	(isthin .eq. 1	) .and. (isincompr .eq. 1	)	) then
			loclimmat =			7
		elseif 	(	(isthin .eq. 0	) .and. (isincompr .eq. 1	)	) then
			loclimmat =			8
		elseif 	(	(isthin .eq. -1) .and. (isincompr .eq. 1	)	) then
			loclimmat =			9
		end if
!		setting Lame parameters
		if (isincompr .eq. 1) then
			lambda =					(capEy*nudummy)/((1.0+nudummy)*(1.0-2.0*nudummy))
			capG = 						capEy/3.0
			capKbulk = 					capEy/(3.0*(1.0-2.0*nudummy))
		elseif (isincompr .eq. -1) then
			lambda =					0.0
			capG = 						capEy/2.0
			capKbulk = 					0.0
		end if
	else
!		limits-not-considered mode
		isincompr = 			0
		if (nupois .ge. thresincompr) then
			isincompr = 				1
			lambda =					(capEy*nudummy)/((1.0+nudummy)*(1.0-2.0*nudummy))
			capG = 						capEy/(2.0*(1.0+nupois))
			capKbulk = 					capEy/(3.0*(1.0-2.0*nudummy))
		end if
	end if
	
!	setting the Prony series constants
	if (model .eq. 1) then
		capKbulkprony =						0.0
		capGprony =							0.0
	end if
	allocate(capKbulkpronynd				(0:nprony))
	allocate(capGpronynd					(0:nprony))
	allocate(kappaprony						(1:nprony))
	allocate(gammaprony						(1:nprony))
	allocate(kprony							(0:nprony))
	allocate(gprony							(0:nprony))
	capKbulkpronynd =						0.0
	capGpronynd =							0.0
	kappaprony =							0.0
	gammaprony =							0.0
	kprony =								0.0
	gprony =								0.0
	capKbulkc = 							capKbulk
	capGc = 								capG
	do iprony = 1,nprony
		capKbulkc = 						capKbulkc + capKbulkprony(iprony)
		capGc = 							capGc + capGprony(iprony)
	end do
	capKbulkpronynd(0) = 					capKbulk/(3*capKbulkc+4*capGc)
	capGpronynd(0) = 						capG/(3*capKbulkc+4*capGc)
	do iprony = 1,nprony
		capKbulkpronynd(iprony) =			capKbulkprony(iprony)/(3*capKbulkc+4*capGc)
		capGpronynd(iprony) =				capGprony(iprony)/(3*capKbulkc+4*capGc)
	end do
	capKbulkcnd = 							0.0
	capGcnd = 								0.0
	do iprony = 0,nprony
		capKbulkcnd = 						capKbulkcnd + capKbulkpronynd(iprony)
		capGcnd = 							capGcnd + capGpronynd(iprony)
	end do
	do iprony = 1,nprony
		kappaprony(iprony) =		 		1.0/(2*pi*omega*taukprony(iprony))
		gammaprony(iprony) = 				1.0/(2*pi*omega*taugprony(iprony))
	end do
	do iprony = 0,nprony
		kprony(iprony) = 					capKbulkpronynd(iprony)/capKbulkcnd
		gprony(iprony) = 					capGpronynd(iprony)/capGcnd
	end do
!	setting same frequency-grid size as time-grid size and same Hankel-grid size as radial-grid size
	nfr =	nt
	nx =	nr

!	obtaining far-end value of x for kernel-mode x-axis adaptation
	xinf =								Bslroot(nr-1)/rinf
	
!	calculating constant non-dimensional parameters
	epsi =								capD/capR
	alpha =								h0/capD
	beta =								capL/capR

!	setting of characteristic total pressure, characteristic substrate rigidity modulus and Debye length
	containers(1) =						((abs(psis(1))+abs(psis(2)))/(abs(psis(1))+abs(psis(2))+1e-15))
	pic =								((mu*omega*alpha)/(epsi*(1.0-alpha)**2)) + 2*n0*kcapB*capT*containers(1) + &
										(capAsfw/(6*pi*epsi**3*capR**3*(1.0-alpha)**3)) + abs(capA)
	capK =								sqrt((elementaryq**2*n0*(abs(valence(1)*valence(2)**2)+abs(valence(1)**2*valence(2))))/ &
										(relperm*permvac*kcapB*capT))

!	calculating constant non-dimensional parameters emerging from scaling principles
	containers(2) =						sqrt(epsi)
	delta =								min(beta,containers(2))
	if (model .eq. 0) then
		theta =							((delta*pic)/(lambda+2*capG))
	else
		theta =							((3*delta*pic)/(3*capKbulkc+4*capGc))	
	end if
	zeta =								beta/delta
	nu =								beta/containers(2)
	capGamma =							theta/epsi

!	condition checks (linear strain and small Wommersley number)
	if (pic/(lambda+2*capG) .ge. threslin) then
		if (termout .eq. 1)	then
			write(*,*)	'CAUTION: Linear-strain assumption violated (due to compliance): '
			write(*,*) 	'Max Strain Expected = (pic/(lambda+2*capG)) = '  , (pic/(lambda+2*capG))
			write(*,*) 	'Threshold set                                 = ', threslin
		end if
		write(95,*)		'CAUTION: Linear-strain assumption violated (due to compliance): '
		write(95,*) 	'Max Strain Expected = (pic/(lambda+2*capG)) = '  , (pic/(lambda+2*capG))
		write(95,*) 	'Threshold set                                 = ', threslin
		stop
	end if
	if (capD/sqrt(mu/(rho*omega)) .ge. thresstokes) then
		if (termout .eq. 1)	then
			write(*,*)	'CAUTION: Condition for Stokes flow not satisfied'
			write(*,*)	'D/sqrt(mu/(rho*omega)) = ', D/sqrt(mu/(rho*omega))
			write(*,*)	'Threshold set as =       ', thresstokes
		end if
		write(95,*)		'CAUTION: Condition for Stokes flow not satisfied'
		write(95,*)		'D/sqrt(mu/(rho*omega)) = ', D/sqrt(mu/(rho*omega))
		write(95,*)		'Threshold set as =       ', thresstokes
		stop
	end if

!	determining system dynamics regimes (i.e. small-amplitude or not, Fourier-space-reduced viscoelastic formulation or not)
	if (alpha .le. thressmallamp) then
		smallamp =		1
	else
		smallamp =		0
		if (model .ne. 0) then
			if (iforceharmonicve .ne. 0) then
				iforceharmonicve = 	0
				if (termout .eq. 1)	then
					write(*,*)		'WARNING: oscillation amplitude is large, overriding forced harmonic'
				end if
				write(95,*)			'WARNING: oscillation amplitude is large, overriding forced harmonic'			
			end if
		end if
	end if
	if (	(iforceharmonicve .ne. 0) .and. (.not. &
			(((n0 .eq. 0) .or. ((psis(1)+psis(2)) .eq. 0)) .and. (capAsfw .eq. 0) .and. (capA .eq. 0)))	.and.	&
			((((mu*omega*alpha)/(epsi*(1.0-alpha)**2))/pic) .le. thresdomhd)	) then
			iforceharmonicve = 	0
		if (termout .eq. 1)	then
			write(*,*)		'WARNING: non-hydrodynamic forces are present and not weak, overriding forced harmonic'
		end if
		write(95,*)			'WARNING: non-hydrodynamic forces are present and not weak, overriding forced harmonic'			
	end if
	if (model .ne. 0) then
		reduce =			0
		if (capGamma .le. thresoneway1) then
			if (smallamp .eq. 1) then
				if ((theta/sqrt(epsi)) .le. thresreduce) then
					if ((capGamma/alpha) .ge. thresoneway2) then
						reduce =	1
						if (termout .eq. 1) then
							write(*,*)	'WARNING: High Interfacial Velocity Two-Way Coupling - bulk field output is incorrect'
						end if
						write(95,*) 	'WARNING: High Interfacial Velocity Two-Way Coupling - bulk field output is incorrect'
					else
						reduce =	2
					end if
				else
					if (termout .eq. 1) then
						write(*,*)	'CAUTION: The second condition for reduction to Fourier-space viscoelastic not met'
						write(*,*)	'(theta/sqrt(epsi)) = ',(theta/sqrt(epsi)), ' is higher than ', thresreduce
					end if
					write(95,*) 	'CAUTION: The second condition for reduction to Fourier-space viscoelastic not met'
					write(95,*)		'(theta/sqrt(epsi)) = ',(theta/sqrt(epsi)), ' is higher than ', thresreduce
				end if
			else
				if (termout .eq. 1) then
					write(*,*)		'CAUTION: Oscillation amplitude or range of motion is too high'
				end if
				write(95,*)			'CAUTION: Oscillation amplitude or range of motion is too high'
			end if
		else
			if (termout .eq. 1) then
				write(*,*)			'CAUTION: Viscoelastic substrate but two-way coupling imminent - killing simulation'
				write(*,*)			'capGamma = ', capGamma, ' is higher than ', thresoneway1
			end if
			write(95,*)				'CAUTION: Viscoelastic substrate but two-way coupling imminent - killing simulation'		
			write(95,*)				'capGamma = ', capGamma, ' is higher than ', thresoneway1
		end if
	end if			
!	stopping simulation if reduced formulation is not obtained
	if ((model .ne. 0) .and. (reduce .eq. 0)) then
		stop
		if (termout .eq. 1) then
			write(*,*)			'CAUTION: Reduced formulation not obtained - killing simulation'
		end if
		write(95,*)				'CAUTION: Reduced formulation not obtained - killing simulation'
	end if
!	overriding different features for viscoelastic solution with brute-force
	if ((model .ne. 0) .and. (scheme .ne. 0)) then
		if (termout .eq. 1) then
			write(*,*)				'WARNING: Brute-force iteration mode not built-in for viscoelastic - overriding'
		end if
		write(95,*)					'WARNING: Brute-force iteration mode not built-in for viscoelastic - overriding'
		scheme = 0
	end if
	if ((model .ne. 0) .and. (reduce .ne. 0) .and. (adt .ne. 0)) then
		if (termout .eq. 1) then
			write(*,*)				'WARNING: Adaptive time grid not built-in for viscoelastic - overriding'
		end if
		write(95,*)					'WARNING: Adaptive time grid not built-in for viscoelastic - overriding'
		adt = 0
	end if
	if ((model .ne. 0) .and. (reduce .ne. 0) .and. (ispring .ne. 0)) then
		if (termout .eq. 1) then
			write(*,*)				'WARNING: Spring connection not built-in for viscoelastic - overriding'
		end if
		write(95,*)					'WARNING: Spring connection not built-in for viscoelastic - overriding'
		ispring = 0
	end if

!	for problem of approach, setting start time as zero and limiting end time to undeformed contact 
	if (iappr .ne. 0) then
		if (termout .eq. 1) then
			write(*,*)		'Note: approach mode solution - shifting time to start at zero'
		end if
		write(94,*)			'Note: approach mode solution - shifting time to start at zero'
		tend =				tend-tstart
		tstart =			0.0
		if ((ispring .eq. 0) .and. (tend .gt. 2.0)) then
			if (termout .eq. 1) then
				write(*,*)	'Warning: no spring attached, prescribed time duration too high, manually overriding'
			end if
			write(94,*)		'Warning: no spring attached, prescribed time duration too high, manually overriding'
			tend = 			2.0
		end if
	end if
	
!	for oscillatory loading with Fourier-mode reduced formulation for viscoelastic material, setting time from 0 to 2*pi 
	if ((model .ne. 0) .and. (reduce .ne. 0) .and. (iappr .eq. 0)) then
		if (termout .eq. 1) then
			write(*,*)		'Note: oscillatory mode solution with Fourier-reduced formulation - shifting time to start at zero'
		end if
		write(94,*)			'Note: oscillatory mode solution with Fourier-reduced formulation - shifting time to start at zero'
		tend =				tend-tstart
		tstart =			0.0
		if (tend .gt. 2.0*pi) then
			if (termout .eq. 1) then
				write(*,*)	'Warning: prescribed time duration too high, manually overriding'
			end if
			write(94,*)		'Warning: prescribed time duration too high, manually overriding'
			tend = 			2.0*pi
		end if
	end if
	
!	overriding adaptive time grid if working with Fourier-mode reduced formulation for viscoelastic material 
	if ((model .ne. 0) .and. (reduce .ne. 0) .and. (adt .ne. 0)) then
		adt = 				0
		if (termout .eq. 1) then
			write(*,*)	'Warning: adaptive time grid cannot be used in Fourier mode, manually overriding'
		end if
		write(94,*)		'Warning: adaptive time grid cannot be used in Fourier mode, manually overriding'
	end if

!	calculating constant non-dimensional parameters emerging from scaling principles - for spring
	if (ispring .ne. 0) then
		thetacapS =						(epsi*pic*capR)/capKspring
	else
		thetacapS =						0.0
	end if
	capGammacapS =						thetacapS/epsi
	
end subroutine constparam
!-----------------------------------------------------------------------------------------------------------------------------------
