!-----------------------------------------------------------------------------------------------------------------------------------
! time-step solution output write subroutine (for linelast [AND VISCELASTFULL])
!
! writes individual time step solution to output files for linear-elastic substrate solution
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine instout

	use varinit
	use derivative
	implicit none

!	writing instance-wise output to files

	if (smallamp .eq. 0) then
		write(52,*)			(containers(362)*r(ir), ir=1,nr)
		if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then 
			write(84,*)		(containers(362)*r1(ir), ir=1,nr)
		end if
	end if
	if (smallamp .eq. 0) then
		write(16,*)		(capX(ix), ix=1,nx)
	end if
	if ((model .ne. 0) .and. (reduce .eq. 0)) then
		write(15,*)			(tau(it), it=1,nt)
		write(55,*)			(tau(it)/omega, it=1,nt)
	end if

	if (minout .ne. 5) then
		write(21,*)			(l(ir), ir=1,nr)
		write(61,*)			(containers(363)*l(ir), ir=1,nr)
		write(22,*)			(p(ir), ir=1,nr)
		write(62,*)			(containers(364)*p(ir), ir=1,nr)
		write(23,*)			(phd(ir), ir=1,nr)
		write(63,*)			(containers(365)*phd(ir), ir=1,nr)
		if (.not. ((n0 .eq. 0.0) .or. (abs(psis(1))+abs(psis(2)) .eq. 0))) then
			write(24,*)		(pdl(ir), ir=1,nr)
			write(64,*)		(containers(366)*pdl(ir), ir=1,nr)
		end if
		if (capAsfw .ne. 0) then
			write(25,*)		(pvdw(ir), ir=1,nr)
			write(65,*)		(containers(367)*pvdw(ir), ir=1,nr)
		end if
		if (capA .ne. 0) then
			write(26,*)		(ps(ir), ir=1,nr)
			write(66,*)		(containers(368)*ps(ir), ir=1,nr)
		end if
	end if

	if (smallamp .eq. 1) then
		write(221,*)		(capH(ir), ir=1,nr)
	end if
	write(261,*)			(containers(369)*capH(ir), ir=1,nr)
	if (ispring .ne. 0) then
		if (smallamp .eq. 0) then
			write(222,*)	(capH(ir)+ncapS(it)*capGammacapS*lcapS(it), ir=1,nr)
		else
			write(222,*)	(capH(ir)+capGammacapS*lcapS(it), ir=1,nr)
		end if
		write(262,*)		(containers(369)*capH(ir)+containers(375)*lcapS(it), ir=1,nr)
	end if

	if (.not. ((minout .eq. 2) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
		do iz = 1,nz
			write(31,*)		(vr(ir,iz), ir=1,nr)
			write(32,*)		(vz(ir,iz), ir=1,nr)
			write(33,*)		(z(ir,iz), ir=1,nr)
			write(71,*)		(containers(370)*vr(ir,iz), ir=1,nr)
			write(72,*)		(containers(317)*vz(ir,iz), ir=1,nr)
			write(73,*)		(containers(369)*z(ir,iz), ir=1,nr)
		end do
	end if

	if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
		do iy = 1,ny
			write(41,*)		(ur(ir,iy), ir=1,nr)
			write(42,*)		(uy(ir,iy), ir=1,nr)
		end do
		if (smallamp .eq. 0) then
			write(43,*)		(y(iy), iy=1,ny)
		end if
		do iy = 1,ny
			write(81,*)		(containers(371)*ur(ir,iy), ir=1,nr)
			write(82,*)		(containers(363)*uy(ir,iy), ir=1,nr)
		end do
		if (smallamp .eq. 0) then
			write(83,*)		(containers(372)*y(iy), iy=1,ny)
		end if
	end if

	if (minout .ne. -1) then
		write(121,*)		l(1)
		write(161,*)		containers(363)*l(1)
		write(122,*)		p(1)
		write(162,*)		containers(364)*p(1)
		write(123,*)		phd(1)
		write(163,*)		containers(365)*phd(1)
		if (.not. ((n0 .eq. 0.0) .or. (abs(psis(1))+abs(psis(2)) .eq. 0))) then
			write(124,*)	pdl(1)
			write(164,*)	containers(366)*pdl(1)
		end if
		if (capAsfw .ne. 0) then
			write(125,*)	pvdw(1)
			write(165,*)	containers(367)*pvdw(1)
		end if
		if (capA .ne. 0) then
			write(126,*)	ps(1)
			write(166,*)	containers(368)*ps(1)
		end if
		if (it .eq. nodeph(iph)) then
			write(131,*)		(l(ir), ir=1,nr)
			write(171,*)		(containers(363)*l(ir), ir=1,nr)
			write(132,*)		(p(ir), ir=1,nr)
			write(172,*)		(containers(364)*p(ir), ir=1,nr)
			write(133,*)		(phd(ir), ir=1,nr)
			write(173,*)		(containers(365)*phd(ir), ir=1,nr)
			if (.not. ((n0 .eq. 0.0) .or. (abs(psis(1))+abs(psis(2)) .eq. 0))) then
				write(134,*)	(pdl(ir), ir=1,nr)
				write(174,*)	(containers(366)*pdl(ir), ir=1,nr)
			end if
			if (capAsfw .ne. 0) then
				write(135,*)	(pvdw(ir), ir=1,nr)
				write(175,*)	(containers(367)*pvdw(ir), ir=1,nr)
			end if
			if (capA .ne. 0) then
				write(136,*)	(ps(ir), ir=1,nr)
				write(176,*)	(containers(368)*ps(ir), ir=1,nr)
			end if
			if (smallamp .eq. 0) then
				write(113,*)	(containers(362)*r(ir), ir=1,nr)
			end if
			iph = 			iph+1
		end if
	end if

end subroutine instout
!-----------------------------------------------------------------------------------------------------------------------------------
