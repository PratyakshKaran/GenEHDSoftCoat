!-----------------------------------------------------------------------------------------------------------------------------------
! matrix space allocation and initialization subroutine

! allocates size to matrices and arrays and initializes with zero or default values
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
subroutine allocinit

	use varinit
	implicit none

! 	allocating dimension to matrices/arrays

!	common matrices/arrays
	allocate(t						(1:nt))
	allocate(r						(1:nr))
	allocate(x						(1:nx))
	if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
		allocate(y					(1:ny))
	end if
	allocate(capH					(1:nr))
	allocate(lcapS					(1:nt))
	allocate(phd					(1:nr))
	allocate(pdl					(1:nr))
	allocate(pvdw					(1:nr))
	allocate(ps						(1:nr))
	allocate(capF					(1:nt))
	allocate(capF_d					(1:nt))
	allocate(capX					(1:nx))
	if (.not. ((minout .eq. 2) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
		allocate(vr					(1:nr,1:nz))
		allocate(vz					(1:nr,1:nz))
		allocate(dphddr				(1:nr))
		allocate(d2phddr2			(1:nr))
		allocate(dldr				(1:nr))
		allocate(z					(1:nr,1:nz))
	end if
	allocate(Hnkkrnlf				(1:nx,1:nr))
	allocate(Hnkkrnli				(1:nr,1:nx))

!	reduced viscoelastic solution matrices/arrays
	if ((model .ne. 0) .and. (reduce .ne. 0)) then
		allocate(lve				(1:nr,1:nt))
		allocate(pve				(1:nr,1:nt))
		allocate(pve_holder			(1:nr,1:nfr))
		if (reduce .eq. 1) then
			allocate(pmve			(1:nr,1:nt))
			allocate(pmve_holder	(1:nr,1:nfr))
			allocate(phdve			(1:nr,1:nt))
			allocate(pmvecapH0		(1:nx))
			allocate(lmvecapH0		(1:nx))
			allocate(probevel		(1:nt))
			allocate(phdfreq		(1:nr))
			allocate(phdfreqprev	(1:nr))
			allocate(lfreq			(1:nr))
			allocate(lfreqprev		(1:nr))
		else
			allocate(lvecapH0		(1:nx,1:nfr))
			allocate(pvecapH0		(1:nx,1:nfr))
		end if
		allocate(capGbar			(1:nfr))
		allocate(capKbulkbar		(1:nfr))
		allocate(fr					(1:nfr))
		if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
			allocate(urve			(1:nr,1:ny,1:nt))
			allocate(uyve			(1:nr,1:ny,1:nt))
			allocate(urvecapH1		(1:nx,1:ny))
			allocate(uyvecapH0		(1:nx,1:ny))
			allocate(gcapHve		(1:nx,1:ny))
			allocate(dgcapHvedy		(1:nx,1:ny))
			allocate(d2gcapHvedy2	(1:nx,1:ny))
			allocate(capAcoeffve	(1:nx))
			allocate(capBcoeffve	(1:nx))
			allocate(capCcoeffve	(1:nx))
			allocate(capDcoeffve	(1:nx))
			allocate(Hnkkrnli		(1:nr,1:nx))
			allocate(r1				(1:nr))
		end if
	end if

!	linear elastic solution and full viscoelastic solution matrices/arrays
	if (.not. ((model .ne. 0) .and. (reduce .ne. 0))) then
		allocate(p					(1:nr))
		allocate(l					(1:nr))
		allocate(lcapH0				(1:nx))
		allocate(pcapH0				(1:nx))
		if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
			allocate(ur				(1:nr,1:ny))
			allocate(uy				(1:nr,1:ny))
			allocate(urcapH1		(1:nx,1:ny))
			allocate(uycapH0		(1:nx,1:ny))
			allocate(gcapH			(1:nx,1:ny))
			allocate(dgcapHdy		(1:nx,1:ny))
			allocate(d2gcapHdy2		(1:nx,1:ny))
			allocate(capAcoeff		(1:nx))
			allocate(capBcoeff		(1:nx))
			allocate(capCcoeff		(1:nx))
			allocate(capDcoeff		(1:nx))
			allocate(Hnkkrnli		(1:nx,1:nr))
			allocate(r1				(1:nr))
		end if
		if (model .ne. 0) then
			allocate(tau			(1:ntau))
		end if
		if (smallamp .eq. 0) then
			allocate(h				(1:nt))
			allocate(n				(1:nt))
			allocate(ncapS			(1:nt))
			allocate(d				(1:nt))
			allocate(s				(1:nt))
			allocate(m				(1:nt))
			allocate(w				(1:nt))
			allocate(dndtbyn		(1:nt))
			allocate(dncapSdtbyncapS(1:nt))
			allocate(dhdtbyh		(1:nt))
		end if
		allocate(dpdl				(1:nr))
		allocate(dpDLdl				(1:nr))
		allocate(dpvdWdl			(1:nr))
		allocate(dpSdl				(1:nr))		
		allocate(oneway				(1:nt))
		allocate(llast				(1:nr))
		allocate(llast2				(1:nr))
		allocate(lprev				(1:nt))
		allocate(phdprev			(1:nt))
		allocate(tmparr1			(1:nr))
		allocate(tmparr2			(1:nr))
		allocate(tmparr3			(1:nr))
		allocate(tmparr4			(1:nr))
		allocate(tmparx1			(1:nx))
		allocate(tmparx2			(1:nx))
		allocate(tmparx3			(1:nx))
		allocate(tmparx4			(1:nx))
	end if
	
!	allocating jacobian and residual matrices
	if (.not. ((model .ne. 0) .and. (reduce .ne. 0))) then 
		if (scheme .eq. -1) then
			allocate(jacobian1		(1:nr,1:nr))
			allocate(residual1		(1:nr))
			allocate(jacobian1_in	(1:nr,1:nr))
			allocate(residual1_in	(1:nr))
			allocate(ipiv1			(1:nr))
		else
			allocate(jacobian		(1:nr+nx+1,1:nr+nx+1))
			allocate(residual		(1:nr+nx+1))
			allocate(jacobian_in	(1:nr+nx+1,1:nr+nx+1))
			allocate(residual_in	(1:nr+nx+1))
			allocate(ipiv			(1:nr+nx+1))
			allocate(solution		(1:nr+nx+1))
		end if
	end if
	if ((model .ne. 0) .and. (reduce .eq. 1)) then 
		allocate(jacobianfr			(1:nr+nx,1:nr+nx))
		allocate(residualfr			(1:nr+nx))
		allocate(jacobianfr_in		(1:nr+nx,1:nr+nx))
		allocate(residualfr_in		(1:nr+nx))
		allocate(ipiv				(1:nr+nx))
		allocate(solutionfr			(1:nr+nx))
	end if

!	initializing matrices and arrays
!	common matrices/arrays
	t =								0.0
	r =								0.0
	x =								0.0
	if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then 
		y =							0.0
	end if
	capH =							0.0
	lcapS =							0.0
	phd =							0.0
	pdl =							0.0
	pvdw =							0.0
	ps =							0.0
	capF =							0.0
	capF_d =						0.0
	capX =							0.0
	if (.not. ((minout .eq. 2) .or. (minout .eq. 4) .or. (minout .eq. 5))) then 
		vr =						0.0
		vz =						0.0
		dphddr =					0.0
		d2phddr2 =					0.0
		dldr = 						0.0
		z =							0.0
	end if
	Hnkkrnlf =						0.0
	Hnkkrnli =						0.0

!	reduced viscoelastic solution matrices/arrays
	if ((model .ne. 0) .and. (reduce .ne. 0)) then 
		lve =						cmplx(0.0,0.0)
		pve =						cmplx(0.0,0.0)
		pve_holder =				cmplx(0.0,0.0)
		if (reduce .eq. 1) then
			pmve =					cmplx(0.0,0.0)
			pmve_holder =			cmplx(0.0,0.0)
			phdve =					cmplx(0.0,0.0)
			pmvecapH0 =				cmplx(0.0,0.0)
			lmvecapH0 =				cmplx(0.0,0.0)
			probevel =				cmplx(0.0,0.0)
			phdfreq =				cmplx(0.0,0.0)
			phdfreqprev =			cmplx(0.0,0.0)
			lfreq =					cmplx(0.0,0.0)
			lfreqprev =				cmplx(0.0,0.0)
		else
			lvecapH0 =				cmplx(0.0,0.0)
			pvecapH0 =				cmplx(0.0,0.0)
		end if
		capGbar =					cmplx(0.0,0.0)
		capKbulkbar =				cmplx(0.0,0.0)
		fr =						0.0
		if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then 
			urve =					cmplx(0.0,0.0)
			uyve =					cmplx(0.0,0.0)
			urvecapH1 =				cmplx(0.0,0.0)
			uyvecapH0 =				cmplx(0.0,0.0)
			gcapHve =				cmplx(0.0,0.0)
			dgcapHvedy =			cmplx(0.0,0.0)
			d2gcapHvedy2 =			cmplx(0.0,0.0)
			capAcoeffve =			cmplx(0.0,0.0)
			capBcoeffve =			cmplx(0.0,0.0)
			capCcoeffve =			cmplx(0.0,0.0)
			capDcoeffve =			cmplx(0.0,0.0)
			Hnkkrnli1 =				0.0
			r1 =					0.0
		end if
	end if

!	linear elastic solution and full viscoelastic solution matrices/arrays
	if (.not. ((model .ne. 0) .and. (reduce .ne. 0))) then 
		p =							0.0
		l =							0.0
		pcapH0 =					0.0
		lcapH0 =					0.0
		if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then 
			ur =					0.0
			uy =					0.0
			urcapH1 =				0.0
			uycapH0 =				0.0
			gcapH =					0.0
			dgcapHdy =				0.0
			d2gcapHdy2 =			0.0
			capAcoeff =				0.0
			capBcoeff =				0.0
			capCcoeff =				0.0
			capDcoeff =				0.0
			Hnkkrnli1 =				0.0
			r1 = 					0.0
		end if	
		if (model .ne. 0) then 
			tau =					0.0
		end if
		if (smallamp .eq. 0) then
			h =						1.0
			n =						1.0
			ncapS =					1.0
			d =						1.0
			s =						1.0
			m =						1.0
			w =						1.0		
			dndtbyn = 				0.0
			dncapSdtbyncapS = 		0.0
			dhdtbyh = 				0.0
		end if
		dpdl =						0.0
		dpdldl =					0.0
		dpvdwdl =					0.0
		dpsdl =						0.0
		if ((smallamp .eq. 1) .and. &
			((capGamma .ge. thresoneway1) .or. (capGammacapS .ge. thresoneway1) .or. &
			((capGamma/alpha) .ge. thresoneway2) .or. ((capGammacapS/alpha) .ge. thresoneway2))) then
			oneway =				0
		else
			oneway =				1
		end if
		llast =						0.0
		llast2 =					0.0
		lprev =						0.0
		phdprev =					0.0
		tmparr1 =					0.0
		tmparr2 =					0.0
		tmparr3 =					0.0
		tmparr4 =					0.0
		tmparx1 =					0.0
		tmparx2 =					0.0
		tmparx3 =					0.0
		tmparx4 =					0.0
	end if
	if (.not. ((model .ne. 0) .and. (reduce .ne. 0))) then 
		if (scheme .eq. -1) then
			jacobian1 =				0.0
			residual1 =				0.0
			jacobian1_in =			0.0
			residual1_in =			0.0
			ipiv1 =					0
		else
			jacobian_in =			0.0
			residual_in =			0.0
			ipiv = 					0
			solution =				0.0
		end if
	end if
	if ((model .ne. 0) .and. (reduce .eq. 1)) then 
		jacobianfr =			cmplx(0.0,0.0)
		residualfr =			cmplx(0.0,0.0)
		jacobianfr_in =			cmplx(0.0,0.0)
		residualfr_in =			cmplx(0.0,0.0)
		ipiv = 					0
		solutionfr =			cmplx(0.0,0.0)
	end if

!	creating files for dumping output

!	real-space and transform space interfacial-grid variables
	open(		unit=11,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/t.dat', &
				status='new',action='write')
	open(		unit=12,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/r.dat', &
				status='new',action='write')
	open(		unit=13,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/x.dat', &
				status='new',action='write')
	if ((model .ne. 0) .and. (reduce .ne. 0)) then
		open(	unit=14,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/fr.dat', &
				status='new',action='write')
	end if
!	**unit=15 has been used up for tau**
	open(		unit=16,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/capX.dat', &
				status='new',action='write')
	open(		unit=51,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/t_d.dat', &
				status='new',action='write')
	open(		unit=52,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/r_d.dat', &
				status='new',action='write')
	if ((model .ne. 0) .and. (reduce .eq. 0)) then
		open(	unit=15,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/tau.dat', &
				status='new',action='write')
		open(	unit=55,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/tau_d.dat', &
				status='new',action='write')
	end if

!	interfacial-field variables
	if (minout .ne. 5) then
		open(	unit=21,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/l.dat', &
				status='new',action='write')
		open(	unit=22,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/p.dat', &
				status='new',action='write')
		open(	unit=23,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/phd.dat', &
				status='new',action='write')
		open(	unit=24,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/pdl.dat', &
				status='new',action='write')
		open(	unit=25,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/pvdw.dat', &
				status='new',action='write')
		open(	unit=26,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/ps.dat', &
				status='new',action='write')
		open(	unit=61,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/l_d.dat', &
				status='new',action='write')
		open(	unit=62,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/p_d.dat', &
				status='new',action='write')
		open(	unit=63,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/phd_d.dat', &
				status='new',action='write')
		open(	unit=64,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/pdl_d.dat', &
				status='new',action='write')
		open(	unit=65,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/pvdw_d.dat', &
				status='new',action='write')
		open(	unit=66,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/ps_d.dat', &
				status='new',action='write')
	end if

!	flow-field variables
	if (.not. ((minout .eq. 2) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
		open(	unit=31,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/vr.dat', &
				status='new',action='write')
		open(	unit=32,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/vz.dat', &
				status='new',action='write')
		open(	unit=33,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/z.dat', &
				status='new',action='write')
		open(	unit=71,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/vr_d.dat', &
				status='new',action='write')
		open(	unit=72,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/vz_d.dat', &
				status='new',action='write')
		open(	unit=73,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/z_d.dat', &
				status='new',action='write')
	end if

!	deformation-field variables
	if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
		open(	unit=41,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/ur.dat', &
				status='new',action='write')
		open(	unit=42,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/uy.dat', &
				status='new',action='write')
		open(	unit=43,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/y.dat', &
				status='new',action='write')
		open(	unit=44,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/r1.dat', &
				status='new',action='write')
		open(	unit=81,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/ur_d.dat', &
				status='new',action='write')
		open(	unit=82,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/uy_d.dat', &
				status='new',action='write')
		open(	unit=83,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/y_d.dat', &
				status='new',action='write')
		open(	unit=84,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/r1_d.dat', &
				status='new',action='write')
	end if

!	interfacial-field origin and quarter-phase-radial variables
	if (minout .ne. -1) then
		open(	unit=111,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/t_q.dat', &
				status='new',action='write')
		open(	unit=112,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/t_q_d.dat', &
				status='new',action='write')
		if (smallamp .eq. 0) then
			open(	unit=113,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/r_q_d.dat', &
					status='new',action='write')
		end if
		open(	unit=121,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/l_o.dat', &
				status='new',action='write')
		open(	unit=122,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/p_o.dat', &
				status='new',action='write')
		open(	unit=123,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/phd_o.dat', &
				status='new',action='write')
		open(	unit=124,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/pdl_o.dat', &
				status='new',action='write')
		open(	unit=125,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/pvdw_o.dat', &
				status='new',action='write')
		open(	unit=126,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/ps_o.dat', &
				status='new',action='write')
		open(	unit=161,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/l_o_d.dat', &
				status='new',action='write')
		open(	unit=162,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/p_o_d.dat', &
				status='new',action='write')
		open(	unit=163,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/phd_o_d.dat', &
				status='new',action='write')
		open(	unit=164,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/pdl_o_d.dat', &
				status='new',action='write')
		open(	unit=165,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/pvdw_o_d.dat', &
				status='new',action='write')
		open(	unit=166,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/ps_o_d.dat', &
				status='new',action='write')
		open(	unit=131,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/l_q.dat', &
				status='new',action='write')
		open(	unit=132,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/p_q.dat', &
				status='new',action='write')
		open(	unit=133,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/phd_q.dat', &
				status='new',action='write')
		open(	unit=134,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/pdl_q.dat', &
				status='new',action='write')
		open(	unit=135,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/pvdw_q.dat', &
				status='new',action='write')
		open(	unit=136,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/ps_q.dat', &
				status='new',action='write')
		open(	unit=171,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/l_q_d.dat', &
				status='new',action='write')
		open(	unit=172,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/p_q_d.dat', &
				status='new',action='write')
		open(	unit=173,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/phd_q_d.dat', &
				status='new',action='write')
		open(	unit=174,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/pdl_q_d.dat', &
				status='new',action='write')
		open(	unit=175,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/pvdw_q_d.dat', &
				status='new',action='write')
		open(	unit=176,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/ps_q_d.dat', &
				status='new',action='write')
	end if

!	interfacial-field variables for sinusoidal-mode solution
	if ((model .ne. 0) .and. (reduce .eq. 1) .and. (iforceharmonicve .ne. 0)) then
		open(	unit=141,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/l_sine.dat', &
				status='new',action='write')
		open(	unit=142,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/phd_sine.dat', &
				status='new',action='write')		
		open(	unit=143,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/G_sine.dat', &
				status='new',action='write')		
		open(	unit=181,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/l_sine_d.dat', &
				status='new',action='write')
		open(	unit=182,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/phd_sine_d.dat', &
				status='new',action='write')
		open(	unit=183,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/G_sine_d.dat', &
				status='new',action='write')
	end if

!	force-responce
	open(		unit=211,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/F.dat', &
				status='new',action='write')
	open(		unit=251,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/F_d.dat', &
				status='new',action='write')

!	probe profile, spring profile and spring compression
	open(		unit=221,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/H.dat', &
				status='new',action='write')
	open(		unit=222,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/HS.dat', &
				status='new',action='write')
	open(		unit=223,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/lS.dat', &
				status='new',action='write')
	open(		unit=261,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/H_d.dat', &
				status='new',action='write')
	open(		unit=262,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/HS_d.dat', &
				status='new',action='write')
	open(		unit=263,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/lS_d.dat', &
				status='new',action='write')

!	constant and transient non-dimensional parameters
	open(		unit=91,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/constparam.dat', &
				status='new',action='write')
	if (smallamp .eq. 0) then
		open(	unit=92,file='../../../../out/ttsi/'//datestamp//'/'//trim(siminstance)//'/transparam.dat', &
				status='new',action='write')
	end if


!	writing case parameters to output file

    write(91,*)  'solvent: '
    write(91,*)  'rho, mu, sigma'
    write(91,*)  rho, mu, sigma
    write(91,*)  'electrolyte: '
    write(91,*)  'nspecies, n0, relperm, (diffcoeff(itemp), itemp=1,nspecies)'//&
				'(valence(itemp), itemp=1,nspecies)'
    write(91,*) nspecies, n0, relperm, (diffcoeff(itemp), itemp=1,nspecies), &
				(valence(itemp), itemp=1,nspecies)
    write(91,*) 'solid: '
    write(91,*) 'capEy, nupois, lambda, capG, capKbulk, capAsfw, capA, phase, model, ' // &
				'iscapLame, psis(1),psis(2), nprony'
    write(91,*) capEy, nupois, lambda, capG, capKbulk, capAsfw, capA, phase, model, &
				iscapLame, psis(1), psis(2), nprony
 	if (model .ne. 0) then
		write(91,*) 'Prony series variables: capKprony, capGprony, taukprony, taugprony'
		do iprony = 1,nprony
			write(91,*)	capKbulkprony(iprony), capGprony(iprony), taukprony(iprony), taugprony(iprony)
		end do
	end if
	write(91,*) 'univ:'
	write(91,*)	'elementaryq, kcapB, permvac'
	write(91,*)	elementaryq, kcapB, permvac
    write(91,*) 'geometry: '
    write(91,*) 'h0, capD, capL, capR, omega, capT, iappr, iforceharmonicve, ispring, capKspring'
    write(91,*) h0, capD, capL, capR, omega, capT, iappr, iforceharmonicve, ispring, capKspring
    write(91,*) 'grid: '
	write(91,*)	'ady, adt, adtau'
	write(91,*)	ady, adt, adtau
	write(91,*)	'mady, nady, madt, nadt, madtau, nadtau'
	write(91,*)	mady, nady, madt, nadt, madtau, nadtau
	write(91,*)	'scheme, symbc, guess, revtime, considerlim'
	write(91,*)	scheme, symbc, guess, revtime, considerlim
	write(91,*)	'minout, termout'
	write(91,*)	minout, termout
	write(91,*)	'tol, tolbreak, thressmallamp, thresoneway1, thresoneway2, relax, thresadhesion, threscavity, threslin, ' // &
				'thresstokes, nsubsvel, thresthin, thressemi, threspcompr, thresincompr, nudummy, thresdomhd, thresreduce'
	write(91,*)	tol, tolbreak, thressmallamp, thresoneway1, thresoneway2, relax, thresadhesion, threscavity, threslin, &
				thresstokes, nsubsvel, thresthin, thressemi, threspcompr, thresincompr, nudummy, thresdomhd, thresreduce
	write(91,*)	'nt, nr, ny, nz, ntau'
	write(91,*)	nt, nr, ny, nz, ntau
	write(91,*)	'tstart, tend, rinf, tauinf'
	write(91,*)	tstart, tend, rinf, tauinf
	write(91,*)	'============================================================================================'
    write(91,*) 'probtype: '
    write(91,*) probtype
    write(91,*) 'parameters: '
    write(91,*) 'epsi, alpha, beta, pic, capK, theta, delta, zeta, nu ,capGamma, capGammacapS, smallamp, reduce, isincompr'
    write(91,*) epsi, alpha, beta, pic, capK, theta, delta, zeta, nu ,capGamma, capGammacapS, smallamp, reduce, isincompr

	close(91)

end subroutine allocinit
!-----------------------------------------------------------------------------------------------------------------------------------
