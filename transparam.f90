!-----------------------------------------------------------------------------------------------------------------------------------
! transient parameter calculation subroutine
!
! calculates pertinent transient parameter values
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine transparam

	use varinit
	use derivative
	implicit none

!	calculating transient non-dimensional parameters
	if (iappr .eq. 0) then
		h(it) =			1.0+alpha*cos(t(it))
	else
		if (iappr .eq. 1) then
			h(it) =		1.0+alpha*(1-t(it))
		else
			h(it) =		1.0-alpha*(1-t(it))
		end if
	end if
	containers(1) =		sqrt(h(it))
	d(it) = 			min(containers(1)*containers(10),beta)/delta
	s(it) =				1.0/d(it)
	m(it) =				(s(it)*d(it))/containers(1)
	n(it) =				d(it)
	ncapS(it) =			h(it)
	w(it) =				n(it)/h(it)

!	a posteriori checks (deflection doesn't overpower z-scaling and simplification of traction balance conditions stay valid)
	ir = 1
	if (l(ir) .ge. 0) then
		if (w(it)*capGamma .ge. threscavity) then
			if (termout .eq. 1)	then
				write(*,*)	'CAUTION: Cavity into substrate large enough to dominate z-scaling'
				write(*,*) 	'w(it)*capGamma = ', w(it)*capGamma
				write(*,*) 	'Threshold set  = ', threscavity
			end if
			write(95,*)		'CAUTION: Cavity into substrate large enough to dominate z-scaling'
			write(95,*) 	'w(it)*capGamma = ', w(it)*capGamma
			write(95,*)	 	'Threshold set  = ', threscavity
			stop
		end if
	else
		if (w(it)*capGamma .ge. thresadhesion) then
			if (termout .eq. 1)	then
				write(*,*)	'CAUTION: Adhesion of substrate and sphere proximate enough to dominate z-scaling'
				write(*,*) 	'w(it)*capGamma = ', w(it)*capGamma
				write(*,*) 	'Threshold set  = ', thresadhesion
			end if
			write(95,*)		'CAUTION: Adhesion of substrate and sphere proximate enough to dominate z-scaling'
			write(95,*) 	'w(it)*capGamma = ', w(it)*capGamma
			write(95,*) 	'Threshold set  = ', thresadhesion
			stop
		end if
	end if
	
!	checking one-way coupling or two-way coupling for the upcoming time-step and computing the required derivative if latter
	if (smallamp .eq. 0) then
		if (iappr .eq. 0) then
			dhdtbyh(it) =		(-alpha*sin(t(it)))/(1.0+alpha*cos(t(it)))
		else
			if (iappr .eq. 1) then
				dhdtbyh(it) =	-alpha/(1.0+alpha*(1-t(it)))
			else
				dhdtbyh(it) =	alpha/(1.0-alpha*(1-t(it)))
			end if
		end if
	end if
	dncapSdtbyncapS(it) = 	dhdtbyh(it)
	if (h(it)*epsi .le. beta**2) then
		dndtbyn(it) = 		0.5*dhdtbyh(it)
	else
		dndtbyn(it) = 		0.0
	end if
	if (	(w(it)*capGamma .ge. thresoneway1) .or. &
			(n(it)*containers(15)*max(max(abs(dndtbyn(it)),abs(dhdtbyh(it))),1.0) .ge. thresoneway2) .or. &
			(capGammacapS .ge. thresoneway1) .or. &
			(ncapS(it)*containers(31)*max(abs(dncapSdtbyncapS(it)),1.0) .ge. thresoneway2)	) then
		oneway(it) = 0
	end if

end subroutine transparam
!-----------------------------------------------------------------------------------------------------------------------------------
