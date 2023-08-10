!-----------------------------------------------------------------------------------------------------------------------------------
! time marching subroutine
!
! marches across time steps and writes time step output to file for linear-elastic substrate
! obtains the aggregate solution for one-way-FSI small-amplitude solution for visco-elastic substrate
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine march

	use varinit
	use derivative
	implicit none

	if (model .eq. 0) then
!		constant containers
		containers(10) = 			sqrt(epsi)
		containers(11) = 			lambda+capG
		containers(12) =			lambda+2*capG
		containers(13) =			lambda+3*capG
		containers(14) =			((zeta*delta*pic)/theta)/capG
		containers(15) =			capGamma/alpha
		containers(31) = 			capGammacapS/alpha
		containers(16) =			32*(tanh((elementaryq*(0.5*(psis(1)+psis(2))))/(4*kcapB*capT)))**2
		containers(17) = 			epsi*capR
		containers(18) =			(mu*omega*alpha)/(epsi*pic)
		containers(19) =			(2*n0*kcapB*capT)/pic
		containers(20) = 			capAsfw/(6*pi*pic*epsi**3*capR**3)
		containers(111) = 			capA/pic
		containers(112) = 			capGamma*epsi*capR
		containers(113) = 			theta/(alpha*epsi)
		containers(114) =			epsi*pic*capR**2
		containers(115) =			(pic*capR)/capKspring
	if (considerlim .eq. 0) then
			if (isincompr .eq. 0) then
				containers(21) =	containers(13)*containers(12)
				containers(22) =	2*containers(11)*containers(12)
				containers(23) =	containers(13)*containers(11)
				containers(24) =	2*(containers(12)**2+capG**2)
				containers(25) =	containers(11)**2
			else
				containers(21) =	1.0
				containers(22) =	2.0
				containers(23) =	1.0
				containers(24) =	2.0
				containers(25) =	1.0
			end if
		else
			if 	(loclimmat .eq. 2) then
				containers(21) =	3.0
				containers(22) =	4.0
				containers(23) =	3.0
				containers(24) =	10.0
				containers(25) =	4.0			
			elseif 	(loclimmat .eq. 5) then
				containers(21) =	1.0*containers(13)*containers(12)
				containers(22) =	4.0*containers(12)*containers(11)
				containers(23) =	2.0*containers(13)*containers(11)
				containers(24) =	4.0*(containers(12)**2+capG**2)
				containers(25) =	8.0*(containers(11)**2)
			elseif 	(loclimmat .eq. 8) then
				containers(21) =	1.0
				containers(22) = 	4.0
				containers(23) =	2.0
				containers(24) =	4.0
				containers(25) = 	8.0			
			elseif 	(loclimmat .eq. 6) then
				containers(21) =	(containers(12)**2)/(2*capG*containers(11))
				containers(22) = 	-(lambda**2+3*lambda*capG+2*capG**2)/(2*lambda**2+9*lambda*capG+11*capG**2)
			elseif 	(loclimmat .eq. 9) then
				containers(21) = 	containers(12)/(2*capG)
			end if
		end if
		containers(26) = 			containers(17)*capK
		containers(27) =			containers(17)/sigma
		containers(28) = 			containers(112)*capK
		containers(29) = 			containers(112)/sigma
		containers(30) = 			containers(12)/capG
		if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
			containers(301) = 		(zeta**2*delta*pic)/(nu**2*theta)
			containers(302) =		zeta/nu
			containers(304) = 		containers(11)**2
			containers(305) = 		2*lambda*containers(11)
			containers(306) = 		capG*containers(13)
			containers(307) = 		lambda*containers(13)
			containers(308) = 		2*containers(304)*containers(11)
			containers(309) = 		4*capG*containers(11)*containers(12)
			containers(310) = 		containers(304)*containers(13)
			containers(311) = 		lambda*containers(13)
			containers(312) = 		delta/containers(10)
			containers(313) = 		((lambda+2*capG)/capG)*((delta**2)/epsi)
		end if
		containers(314) = 			containers(10)*capR
		containers(315) = 			theta*capR
		containers(323) = 			thetacapS*capR
		containers(316) = 			epsi*capR
		if (.not. ((minout .eq. 2) .or. (minout .eq. 4) .or. (minout .eq. 5))) then		
			containers(317) =		alpha*epsi*omega*capR
			containers(318) =		containers(317)/containers(10)
		end if
		if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
			containers(319) = 		delta*capR
		end if
		containers(320) =			1.0/(nsubsvel-1)
		containers(321) =			(capR**2)/Bslroot(nr)
		containers(322) =			1.0/containers(321)
!		transient containers for small-amplitude case (behave like constant containers) 
		if (smallamp .eq. 1)  then
			containers(41) = 		2*nu
			containers(43) = 		containers(26)
			containers(44) = 		containers(27)
			containers(45) = 		containers(18)
			containers(46) = 		containers(19)
			containers(47) = 		containers(20)
			containers(48) = 		containers(111)
			containers(49) = 		1.0
			containers(50) = 		1.0
			containers(141) = 		capGamma
			containers(149) = 		capGammacapS
			containers(142) = 		containers(28)
			containers(143) = 		3*capGamma
			containers(150) = 		3*capGammacapS
			containers(144) = 		containers(29)
			containers(145) = 		1.0
			containers(148) = 		containers(114)
			containers(51) = 		2*containers(41)
			containers(52) =		containers(41)**3
			containers(53) = 		containers(14)
			containers(54) = 		2*pi*containers(44)
			containers(55) = 		12*containers(145)*containers(113)	
			if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
				containers(351) = 	containers(301)
				containers(352) =	containers(302)
				containers(353) = 	nu
				containers(354) = 	containers(353)**2
				containers(355) = 	-containers(312)
				containers(356) = 	containers(313)
			end if
			if (.not. ((minout .eq. 2) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
				containers(361) = 	capGamma
			end if
			containers(362) = 		containers(314)
			if (smallamp .eq. 0) then
				containers(362) = 	containers(314)
			end if
			containers(363) = 		containers(315)
			containers(375) = 		containers(323)
			containers(364) =		pic
			if ((minout .ne. 5) .or. (minout .ne. -1)) then			
				containers(365) =	containers(45)*containers(364)
				containers(366) =	containers(46)*containers(364)
				containers(367) =	containers(47)*containers(364)
				containers(368) =	containers(48)*containers(364)
			end if
			containers(369) = 		containers(316)
			if (.not. ((minout .eq. 2) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
				containers(370) =	containers(318)
			end if
			if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then				
				containers(371) =	containers(363)
				containers(372) = 	containers(319)
			end if
			containers(373) = 		containers(15)
			containers(374) = 		containers(31)
		end if
!		obtaining the Hankel space flexibility function
		if (smallamp .eq. 1) then
			if (considerlim .eq. 0) then
				do ix = 1,nx
					if (x(ix) .eq. 0) then
						capX(ix) =	-(containers(22)/(2*containers(23)+containers(24)))*containers(53)
					else
						capX(ix) =	((containers(21)*(1.0-exp(-containers(51)*x(ix)))- &
									containers(22)*containers(41)*x(ix)*exp(-containers(41)*x(ix)))/ &
									(containers(23)*containers(41)*x(ix)*(1.0+exp(-containers(51)*x(ix)))+ &
									(containers(24)*containers(41)*x(ix)+(containers(25)*containers(52)*x(ix)**3))* &
									exp(-containers(41)*x(ix))))*containers(53)
					end if
				end do
			else
				if ((loclimmat .eq. 1) .or. (loclimmat .eq. 4)) then
					do ix = 1,nx
						capX(ix) =	1.0
					end do
				elseif (loclimmat .eq. 7) then
					do ix = 1,nx
						capX(ix) =	0.0
					end do
				elseif ((loclimmat .eq. 2) .or. (loclimmat .eq. 5) .or. (loclimmat .eq. 8)) then
					do ix = 1,nx
						if (x(ix) .eq. 0) then
							capX(ix) =	&
									-(containers(22)/(2*containers(23)+containers(24)))*containers(53)
						else
							capX(ix) =	&
									((containers(21)*(1.0-exp(-containers(51)*x(ix)))- &
									containers(22)*0.5*containers(41)*x(ix)*exp(-containers(41)*x(ix)))/ &
									(containers(23)*0.5*containers(41)*x(ix)*(1.0+exp(-containers(51)*x(ix)))+ &
									((containers(24)*0.5*containers(41)*x(ix)+containers(25)*0.125*containers(41)**3*x(ix)**3))* &
									exp(-containers(41)*x(ix))))*containers(53)
						end if
					end do
				elseif (loclimmat .eq. 3) then
					do ix = 1,nx
						if (x(ix) .eq. 0) then
							capX(ix) =	&
									-0.25
						else
							capX(ix) =	&
									2.0/x(ix)
						end if
					end do
				elseif (loclimmat .eq. 6) then
					do ix = 1,nx
						if (x(ix) .eq. 0) then
							capX(ix) =	&
									containers(22)
						else
							capX(ix) =	&
									containers(21)/x(ix)
						end if
					end do
				elseif (loclimmat .eq. 9) then
					do ix = 1,nx
						if (x(ix) .eq. 0) then
							capX(ix) =	&
									-0.5
						else
							capX(ix) =	&
									containers(21)/x(ix)
						end if
					end do					
				end if
			end if
		end if
!		initializing jacobian and residual matrices and setting constant values for jacobian and residual matrices as applicable
!		filling spring compression jacobian for spring compression expression			
		if (scheme .eq. 0) then
			jacobian(nr+nx+1,nr+nx+1) = 				1.0/pi
		end if
		if ((scheme .eq. 0) .and. (.not.((smallamp .eq. 1) .and. (capGamma .le. thresoneway1)))) then
			ir = 1
			atemp(1,:,:) =								a('x','fd',ir,ir,r,r);
!			filling hydrodynamic pressure jacobian for hydrodynamic pressure centerline symmetry condition			
			ir = 1
			do jr = 0,2
				jacobian(ir,ir+jr) =					atemp(1,jr,0)
			end do
!			filling hydrodynamic pressure jacobian for hydrodynamic pressure far-end zero condition			
			ir = nr
			jacobian(ir,ir) =							1.0
!			filling deflection jacobian for deflection far-end zero condition
			ir = nr
			jacobian(nr+ir,nr+ir) =						1.0
!			filling deflection jacobian for deflection centerline symmetry condition	
			if (symbc .eq. 1) then
				ir = 1
				do jr = 0,2		
					jacobian(nr+ir,nr+ir+jr) =			atemp(1,jr,0)
				end do
			end if
!			filling hydrodynamic pressure jacobian for Hankel-space deflection-pressure relation
			if ((smallamp .eq. 1) .and. (scheme .eq. 0)) then
				do ix = 2,nx-1
					do ir = 1,nr
						jacobian(nr+ix,ir) =			-Hnkkrnlf(ix,ir)*containers(45)*capX(ix)
					end do
				end do
				if (symbc .eq. 0) then
					ix = 1
					do ir = 1,nr
						jacobian(nr+ix,ir) =			-Hnkkrnlf(ix,ir)*containers(45)*capX(ix)
					end do
				end if
			end if
		end if
		
!	obtaining solution for linear-elastic solid
!		starting time-marching
		do it = 1,nt
!			calculating transient non-dimensionalized parameters
			if (smallamp .eq. 0) then
				call transparam
			end if
!			transient non-iterable containers
			if (it .eq. 1) then
				containers(39) = 		0.0
			else
				if (it .eq. 2) then
					containers(39) = 	1.0/(t(it)-t(it-1))
				else
					atemp(20,:,:) =		a('x','bd',it,it,t,t)
					containers(39) = 	atemp(20,0,0)
				end if
			end if
			if (smallamp .eq. 0) then
				containers(38) = 		-(4*h(it)*containers(39))/alpha
				containers(37) = 		-(4*h(it)*dndtbyn(it))/alpha
				containers(36) = 		-(4*h(it)*dncapSdtbyncapS(it))/alpha
				containers(35) = 		(2*h(it)*dhdtbyh(it))/alpha
			else
				containers(38) = 		-(4*containers(39))/alpha
				containers(37) = 		0
				containers(36) = 		0
				containers(35) = 		0
			end if
			if (iappr .eq. 0) then
				containers(40) = 		-sin(t(it))
			elseif (iappr .eq. 1) then
				containers(40) = 		-1
			else
				containers(40) = 		1
			end if
			containers(42) = 			-3*containers(40)
			containers(147) = 			sin(t(it))
			if (smallamp .eq. 0) then
				containers(41) = 		2*m(it)*nu
				containers(43) = 		containers(26)*h(it)
				containers(44) = 		containers(27)*h(it)
				containers(45) = 		(1.0/(h(it)**2))*containers(18)
				containers(46) = 		containers(19)
				containers(47) = 		(1.0/(h(it)**3))*containers(20)
				containers(48) = 		containers(111)
				if (it .ne. 1) then
					containers(49) = 	n(it)/n(it-1)
					containers(50) = 	h(it-1)/h(it)
				end if
				containers(141) = 		w(it)*capGamma
				containers(149) = 		1.0*capGammacapS
				containers(142) = 		w(it)*h(it)*containers(28)
				containers(143) = 		3*w(it)*capGamma
				containers(150) = 		3*1.0*capGammacapS
				containers(144) = 		w(it)*h(it)*containers(29)
				containers(145) = 		n(it)
				containers(148) = 		containers(114)*h(it)
				containers(51) = 		2*containers(41)
				containers(52) =		containers(41)**3
				containers(53) = 		s(it)*containers(14)
				containers(54) = 		2*pi*containers(44)
				containers(55) = 		12*containers(145)*containers(113)
				if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
					containers(351) = 	((d(it)*s(it)**2)/(n(it)*m(it)**2))*containers(301)
					containers(352) =	(s(it)/m(it))*containers(302)
					containers(353) = 	m(it)*nu
					containers(354) = 	containers(353)**2
					containers(355) = 	-containers(312)*(d(it)/sqrt(h(it)))
					containers(356) = 	containers(313)*((d(it)**2)/h(it))
				end if
				if (.not. ((minout .eq. 2) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
					containers(361) = 	w(it)*capGamma
				end if
				containers(362) = 		containers(314)
				if (smallamp .eq. 0) then
					containers(362) = 	sqrt(h(it))*containers(314)
				end if
				containers(363) = 		n(it)*containers(315)
				containers(375) = 		ncapS(it)*containers(323)
				containers(364) =		pic
				if ((minout .ne. 5) .or. (minout .ne. -1)) then
					containers(365) =	containers(45)*containers(364)
					containers(366) =	containers(46)*containers(364)
					containers(367) =	containers(47)*containers(364)
					containers(368) =	containers(48)*containers(364)
				end if
				containers(369) = 		h(it)*containers(316)
				if (.not. ((minout .eq. 2) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
					containers(370) =	containers(318)/sqrt(h(it))
				end if
				if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
					containers(371) =	containers(363)
					containers(372) = 	d(it)*containers(319)
				end if
				containers(373) = 		n(it)*containers(15)
				containers(374) = 		ncapS(it)*containers(31)
			end if
!			obtaining probe profile for small amplitude case
			if (smallamp .eq. 1) then
				if (iappr .eq. 0) then
					capH(1) =				1.0+alpha*cos(t(it))
				else
					if (iappr .eq. 1) then
						capH(1) =			1.0+alpha*(1.0-t(it))
					else
						capH(1) =			1.0-alpha*(1.0-t(it))
					end if
				end if
				do ir = 2,nr
					capH(ir) =				capH(1) + 0.5*r(ir)**2
				end do
			end if
!			setting guess as zero if zero-switch is on
			select case (guess)
				case (0)
					do ir = 1,nr
						l(ir) =			0.0
						phd(ir) =		0.0
					end do
				case (2)
					do ir = 1,nr
						l(ir) =			0.0
					end do				
				case (3)
					do ir = 1,nr
						phd(ir) =		0.0
					end do					
			end select
!			calling linear-elastic solver
			call linelast
!			terminal output of progress
			if (termout .eq. 1) then 
				write(*,*)	't =', (t(it)),',',(it),'of',(nt)
			end if
!			calling instance file-write
			call instout
!			saving last time-step value for reference
			llast2 = 		llast
			llast = 		l		
		end do

	else
	
!		calling solver
		call viscelast

	end if

end subroutine march
!-----------------------------------------------------------------------------------------------------------------------------------s
