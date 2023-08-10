!-----------------------------------------------------------------------------------------------------------------------------------
! final output subroutine
!
! writes all pertinent variables to output files that aren't written with each time-step
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
subroutine finout

	use varinit
	implicit none

!	writing output to files
	if (minout .ne. -1) then 
		write(111,*)	(tph(iph), 			iph=1,4)
		write(112,*)	(tph(iph)/omega,	iph=1,4)
	end if

	write(11,*)			(t(it), 			it=1,nt)
	write(12,*)			(r(ir),				ir=1,nr)
	if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then 
		write(44,*)		(r1(ir),			ir=1,nr)
	end if
	write(13,*)			(x(ix),				ix=1,nx)
	if (smallamp .eq. 1) then
		write(16,*)		(capX(ix), 			ix=1,nx)
	end if
	if ((model .ne. 0) .and. (reduce .ne. 0)) then 
		write(14,*)		(fr(ifr),			ifr=1,nfr)
	end if
	write(51,*)			(t(it)/omega,		it=1,nt)
	if (smallamp .eq. 1) then 
		write(52,*)		(containers(314)*r(ir), ir=1,nr)
		if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then 
			write(84,*)	(containers(314)*r1(ir), ir=1,nr)
		end if
	end if
	if (smallamp .eq. 0) then
		write(221,*)	(capH(ir), ir=1,nr)
	end if
	if (ispring .ne. 0) then
		write(223,*)		(lcapS(it), it=1,nt)
	end if
	if (ispring .ne. 0) then
		if (smallamp .eq. 0) then
			write(263,*)	(ncapS(it)*containers(323)*lcapS(it), it=1,nt)
		else
			write(263,*)	(containers(323)*lcapS(it), it=1,nt)		
		end if
	end if

	if (smallamp .eq. 0) then
		write(92,*)		(h(it), it=1,nt)
		write(92,*)		(d(it), it=1,nt)
		write(92,*)		(s(it), it=1,nt)
		write(92,*)		(m(it), it=1,nt)
		write(92,*)		(n(it), it=1,nt)
		write(92,*)		(ncapS(it), it=1,nt)
		write(92,*)		(w(it), it=1,nt)
		write(92,*)		(oneway(it), it=1,nt)
		write(92,*)		(dhdtbyh(it), it=1,nt)
		write(92,*)		(dndtbyn(it), it=1,nt)
		write(92,*)		(dncapSdtbyncapS(it), it=1,nt)
	end if

	if ((.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) .and. (smallamp .eq. 1)) then
		write(43,*)		(y(iy), iy=1,ny)
		write(83,*)		(containers(316)*y(iy), iy=1,ny)
	end if

	write(211,*)		(capF(it), it=1,nt)
	write(251,*)		(capF_d(it), it=1,nt)

	if ((model .ne. 0) .and. (reduce .eq. 1) .and. (iforceharmonicve .ne. 0)) then
		write(141,*)	(real(lfreq(ir)), ir=1,nr)
		write(141,*)	(imag(lfreq(ir)), ir=1,nr)
		write(142,*)	(containers(365)*real(phdfreq(ir)), ir=1,nr)
		write(142,*)	(containers(365)*imag(phdfreq(ir)), ir=1,nr)
		write(143,*)	real(capGlosstorsin), imag(capGlosstorsin), sqrt(real(capGlosstorsin)**2+imag(capGlosstorsin)**2), &
						atan(imag(capGlosstorsin)/real(capGlosstorsin)), imag(capGlosstorsin)/real(capGlosstorsin)
		write(181,*)	(real(lfreq(ir)), ir=1,nr)
		write(181,*)	(imag(lfreq(ir)), ir=1,nr)
		write(182,*)	(containers(363)*real(phdfreq(ir)), ir=1,nr)
		write(182,*)	(containers(363)*imag(phdfreq(ir)), ir=1,nr)
		write(183,*)	real(capGdlosstorsin), imag(capGdlosstorsin), sqrt(real(capGdlosstorsin)**2+imag(capGdlosstorsin)**2), &
						atan(imag(capGdlosstorsin)/real(capGdlosstorsin)), imag(capGdlosstorsin)/real(capGdlosstorsin)
	end if

!	closing all output files
	close(93)
	close(94)
	close(95)
	close(96)
	close(97)
	close(11)
	close(12)
	close(13)
	if ((model .ne. 0) .and. (reduce .ne. 0)) then
		close(14)
	end if
	close(16)
	close(51)
	close(52)
	if ((model .ne. 0) .and. (reduce .eq. 0)) then
		close(15)
		close(55)
	end if
	if (minout .ne. 5) then
		close(21)
		close(22)
		close(23)
		close(24)
		close(25)
		close(26)
		close(61)
		close(62)
		close(63)
		close(64)
		close(65)
		close(66)
	end if
	if (.not. ((minout .eq. 2) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
		close(31)
		close(32)
		close(33)
		close(71)
		close(72)
		close(73)
	end if
	if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
		close(41)
		close(42)
		close(43)
		close(81)
		close(82)
		close(83)
	end if
	if (minout .ne. -1) then
		close(111)
		close(112)
		if (smallamp .eq. 0) then
			close(113)
		end if
		close(121)
		close(122)
		close(123)
		close(124)
		close(125)
		close(126)
		close(161)
		close(162)
		close(163)
		close(164)
		close(165)
		close(166)
		close(131)
		close(132)
		close(133)
		close(134)
		close(135)
		close(136)
		if ((model .ne. 0) .and. (reduce .eq. 1) .and. (iforceharmonicve .ne. 0)) then
			close(141)
			close(142)
			close(143)
		end if
		close(171)
		close(172)
		close(173)
		close(174)
		close(175)
		close(176)
		if ((model .ne. 0) .and. (reduce .eq. 1) .and. (iforceharmonicve .ne. 0)) then
			close(181)
			close(182)
			close(183)
		end if
	end if
	close(211)
	close(251)
	close(221)
	close(222)
	close(223)
	close(261)
	close(262)
	close(263)
	if (smallamp .eq. 0) then
		close(92)
	end if

end subroutine finout
!-----------------------------------------------------------------------------------------------------------------------------------
