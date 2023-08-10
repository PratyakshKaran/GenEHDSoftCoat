!-----------------------------------------------------------------------------------------------------------------------------------
! mesh generation subroutine

! generates the grid do t, r, x, and fr
! generates the array do capH do large-amplitude case
! generates the grid do y do small amplitude case 
! generates mesh do z do small amplitude and one-way coupling case
! identifies the four quarter=phase time-instants
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
subroutine meshgen

	use varinit
	implicit none

!	generating time grid
	if (revtime .eq. 0) then
		if (adt .eq. 0) then
			tstep = ((tend-tstart)/(nt-1))
			t(1) =						tstart
			do it = 2,nt
				t(it) =					t(it-1)+tstep
			end do
		else
			nt1 = (nt-1)/2+1
			tend1 =						tstart+((tend-tstart)/2.0)
			tstep =						(1.0/(nt1-1))
			t(1) =						0.0			
			do it1 = 2,nt1
				t(it1) =				t(it1-1)+tstep
			end do
			containers(2) =				(madt*nadt*(tend1-tstart)-nadt*tend1)/ &
										((tstart-madt*(tend1-tstart))*tend1+madt*nadt*(tend1-tstart)**2)
			containers(3) =				((madt-nadt)*(tend1-tstart)-tstart)/ &
										((tstart-madt*(tend1-tstart))*tend1+madt*nadt*(tend1-tstart)**2)
			containers(1) =				-containers(2)*tstart
			do it1 = 1,nt1
				t(it1) =				(containers(1)-t(it1))/(t(it1)*containers(3)-containers(2))
				t(nt-it1+1) =			tend-t(it1)
			end	do
		end if
	else
		if (adt .eq. 0) then
			tstep =						((tend-tstart)/(nt-1))
			t(1) =						omega*pi+tstart
			do it = 2,nt
				t(it) =					t(it-1)+tstep
			end do
		else
			nt1 =						(nt-1)/2+1
			tend1 =						tstart+((tend-tstart)/2.0)
			tstep1 =					(1.0/(nt1-1))
			t(1) =						0.0
			do it1 = 2,nt1
				t(it1) =				t(it1-1)+tstep
			end do
			containers(2) =				((1.0-madt)*(1.0-nadt)*(tend1-tstart)-(1.0-nadt)*tend1)/ &
										((tstart-(1.0-madt)*(tend1-tstart))*tend1+(1.0-madt)*(1.0-nadt)*(tend1-tstart)**2)
			containers(3) =				(((1.0-madt)-(1.0-nadt))*(tend1-tstart)-tstart)/ &
										((tstart-(1.0-madt)*(tend1-tstart))*tend1+(1.0-madt)*(1.0-nadt)*(tend1-tstart)**2)
			containers(1) =				-containers(2)*tstart
			do it1 = 1,nt1
				t(it1) =				omega*pi+((containers(1)-t(it1))/(t(it1)*containers(3)-containers(2)))
				t(nt-it1+1) =			omega*2*pi+tend-t(it1)
			end do
		end if
	end if

!	generating radial grid
	rmin =								0.0
	rmax =								rinf
	containers(1) =						rmax/Bslroot(nr-1)
	ir =								1
	r(ir) =								0.0
	do ir = 2,nr
		r(ir) =							Bslroot(ir-1)*containers(1)
	end do


!	generating Hankel grid
	xmin =								0.0
	xmax =								xinf
	containers(1) =						xmax/Bslroot(nr-1)
	ix =								1
	x(ix) =								0.0
	do ix = 2,nx
		x(ix) =							Bslroot(ix-1)*containers(1)
	end do
	
!	generating y-grid
	if (smallamp .eq. 1) then 
		if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
			ymin =						0.0
			ymax =						beta/delta
			if (ady .eq. 0) then
				ystep =					((ymax-ymin)/(ny-1))
				y(1) =					ymin
				do iy = 2,ny
					y(iy) =				y(iy-1)+ystep
				end do
			else
				ystep =					(1.0/(ny-1))
				y(1) =					0.0
				do iy = 2,ny
					y(iy) =				y(iy-1)+ystep
				end do
				containers(2) =			(mady*nady*(ymax-ymin)-nady*ymax)/((ymin-mady*(ymax-ymin))*ymax+mady*nady*(ymax-ymin)**2)
				containers(3) =			((mady-nady)*(ymax-ymin)-ymin)/((ymin-mady*(ymax-ymin))*ymax+mady*nady*(ymax-ymin)**2)
				containers(1) =			-containers(2)*ymin
				do iy = 1,ny
					y(iy) =				(containers(1)-y(iy))/(y(iy)*containers(3)-containers(2))
				end do
			end if
		end if
	end if

!	generating probe profile (do large-dedomation case)
	if (smallamp .eq. 0) then 
		do ir = 1,nr
			capH(ir) =					1.0 + 0.5*r(ir)**2
		end do
	end if

!	obtaining frequency grid do small-amplitude one-way FSI viscoelastic substrates
	if ((model .ne. 0) .and. (reduce .ne. 0)) then
		frstep =						(1.0/tend)/(nfr-1)
		fr(1) = 0.0
		do ifr = 2,nfr
			fr(ifr) =					fr(ifr-1)+frstep
		end do
	end if

!	identifying quarter-phase time nodes
	containers(1) =						(tend-tstart)/4
	tph(1) =							tstart
	do iph = 2,4
		tph(iph) =						tph(iph-1)+containers(1)
	end do
	iph = 1
	do it = 1,nt
		if ((iph .eq. 1) .and. (t(it) > tph(iph))) then
			nodeph(iph) =				it-1
			tph(iph) =					t(it-1) 
			iph = 2
		end if
		if ((iph .eq. 2) .and. (t(it) > tph(iph))) then
			nodeph(iph) =				it-1
			tph(iph) =					t(it-1) 
			iph = 3
		end if
		if ((iph .eq. 3) .and. (t(it) > tph(iph))) then
			nodeph(iph) =				it
			tph(iph) =					t(it) 
			iph = 4
		end if
		if ((iph .eq. 4) .and. (t(it) >= tph(iph))) then
			nodeph(iph) =				it
			tph(iph) =					t(it) 
			iph = 4
		end if
	end do
	iph = 1

!	generating Hankel transdomation Kernels
	containers(1) =				2.0/Bslroot(nr)
	do ix = 1,nx-1
		do ir = 1,nr-1
			Hnkkrnlf(ix,ir) =	containers(1)*(bessel_jn(0,(Bslroot(ir)*Bslroot(ix))/Bslroot(nr))/ &
								((bessel_jn(1,Bslroot(ir)))**2))
		end do
	end do
	do ix = 1,nx-1
		do ir = 1,nr-1
			Hnkkrnli(ir,ix) =	Hnkkrnlf(ir,ix)
		end do
	end do	

	if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then	
		containers(1) =			rmax/Bslroot1(nr-1)
		ir =					1
		r1(ir) =				0.0
		do ir = 2,nr
			r1(ir) =			Bslroot(ir-1)*containers(1)
		end do
		containers(1) =			2.0/Bslroot1(nr)
		do ix = 1,nx-1
			do ir = 1,nr-1
				Hnkkrnli1(ix,ir) =	&
								containers(1)*(bessel_j1((Bslroot1(ir)*Bslroot1(ix))/Bslroot1(nr))/ &
								(bessel_jn(2,Bslroot1(ir))**2))
			end do
		end do
	end if
	
end subroutine meshgen
!-----------------------------------------------------------------------------------------------------------------------------------
