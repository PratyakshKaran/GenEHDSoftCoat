!-----------------------------------------------------------------------------------------------------------------------------------
! time step computation subroutine for linear-elastic substrate
!
! obtains solution for both small-deflection and large-deflection cases
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine linelast

	use varinit
	use derivative
	implicit none

	if (smallamp .eq. 0) then
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

!	setting non-iterative values for jacobian and residual matrices as applicable
	if ((scheme .ne. -1) .and. (.not.((smallamp .eq. 1) .and. (capGamma .le. thresoneway1)))) then
!		filling hydrodynamic pressure jacobian for Hankel-space deflection-pressure relation
		if ((smallamp .ne. 1) .and. (scheme .eq. 0)) then
			do ix = 2,nx-1
				do ir = 1,nr
					jacobian(nr+ix,ir) =		-Hnkkrnlf(ix,ir)*containers(45)*capX(ix)
				end do
			end do
			if (symbc .eq. 0) then
				ix = 1
				do ir = 1,nr
					jacobian(nr+ix,ir) =		-Hnkkrnlf(ix,ir)*containers(45)*capX(ix)
				end do
			end if
		end if
	end if

	if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
		do ix = 1,nx
			containers(401) = 		containers(353)*x(ix)
			containers(411) = 		exp(containers(401))
			containers(412) = 		exp(-containers(401))
			containers(402) = 		containers(401)**2
			containers(403) = 		(2*containers(402)+1.0)
			containers(404) = 		containers(308)*containers(403)+containers(309)+containers(310)* &
									(exp(containers(41)*x(ix))+exp(-containers(41)*x(ix)))
			containers(405) = 		(containers(304)*containers(403)+containers(305)*containers(401)+containers(306))* &
									containers(411)+containers(311)*containers(412)
			containers(406) = 		containers(11)*(2*containers(401)+1)*containers(411)+containers(13)*containers(412)
			containers(407) = 		(containers(304)*containers(403)-containers(305)*containers(401)+containers(306))* &
									containers(412)+containers(311)*containers(411)
			containers(408) = 		containers(11)*(2*containers(401)-1)*containers(412)-containers(13)*containers(411)
		end do
	end if

!	obtaining solution for one-way coupling
	if (oneway(it) .eq. 1) then
!	 	obtaining node-wise pressure components
		do ir = 1,nr
			phd(ir) =				containers(42)/(capH(ir)**2)
			pdl(ir) =				containers(16)*exp(-containers(43)*capH(ir))
			pvdw(ir) =				-1.0/(capH(ir)**3)
			ps(ir) = 				exp(-containers(44)*capH(ir))*cos(containers(54)*capH(ir)+phase)
			p(ir) = 				containers(45)*phd(ir)+containers(46)*pdl(ir)+containers(47)*pvdw(ir)+containers(48)*ps(ir)
		end do
!		obtaining solution for deflection in Hankel space
		do ix = 1,nx
			pcapH0(ix) =			0.0
			do ir = 1,nr
				pcapH0(ix) =		pcapH0(ix)+Hnkkrnlf(ix,ir)*p(ir)
			end do
		end do
!		obtaining Hankel-space deflection from Hankel-space pressure
		do ix = 1,nx
			lcapH0(ix) =			capX(ix)*pcapH0(ix)
		end do
!		inverting solution of deflection in Hankel-space to get solution of deflection in real-space
		do ir = 1,nr
			l(ir) =					0.0
			do ix = 1,nx
				l(ir) =				l(ir)+Hnkkrnli(ir,ix)*lcapH0(ix)
			end do
		end do
		if (termout .eq. 1) then
			write(*,*)				't = ', (t(it)), ' time step ', (it), ' of ', (nt), ' one-way'
		end if
		write(94,*) 				't = ', (t(it)), ' time step ', (it), ' of ', (nt), ' one-way'

	else

!		obtaining solution using brute-force method of envelope iterations (scheme = -1)
		if (scheme .eq. -1) then
!			envelope iterations for substrate velocity inclusion
			containers(603) =	0.0
			do itersubsvel = 1,nsubsvel-1
				containers(603) =			containers(603)+containers(320)
				capRes = tol*2.0
				iter =	0
				do while (capRes .ge. tol)
					iter =					iter+1
					lprev =					l
					phdprev =				phd
					containers(520) = 		containers(149)*lcapS(it)
!					numerically solving Reynolds equation for hydrodynamic pressure
					jacobian1 = 			0.0
					residual1 = 			0.0
					ir = 1
					atemp(1,:,:) =			a('x','fd',ir,ir,r,r)
					do jr = ir,ir+2
						jacobian1(ir,jr) =	atemp(1,jr-ir,0)
					end do
					residual1(ir) =			0.0
					do ir = 2,nr-1
						atemp(1,:,:) =						a('xx','cd',ir,ir,r,r)
						atemp(2,:,:) =						a('x','cd',ir,ir,r,r)
						if (it .eq. 1) then
							containers(501) = 				0.0
							containers(512) = 				0.0
						else
							if (it .eq. 2) then
								containers(501) =			(l(ir)-llast(ir))*containers(146)	
								containers(502) =			(lcapS(it)-lcapS(it-1))*containers(146)	
							else
								containers(501) =			atemp(20,0,0)*l(ir)+atemp(20,-1,0)*llast(ir)+ &
															atemp(20,-2,0)*llast2(ir)
								containers(502) =			atemp(20,0,0)*lcapS(it)+atemp(20,-1,0)*lcapS(it-1)+ &
															atemp(20,-2,0)*lcapS(it-2)
							end if
						end if
						if (smallamp .eq. 0) then
							containers(501) = 				containers(373)* &
															(containers(501)+l(ir)*dndtbyn(it)- &
															0.5*r(ir)*differential1d(l,'x','cd',ir,ir,'x',r,r)*dhdtbyh(it))
							containers(502) = 				containers(374)*(containers(502)+lcapS(it)*dncapSdtbyncapS(it))
						else
							containers(501) = 				containers(373)*containers(501)
							containers(502) = 				containers(374)*containers(502)
						end if
						do jr = -1,1
							jacobian1(ir,ir+jr) = &
												r(ir)*(capH(ir)+containers(520)+containers(141)*l(ir))**3*atemp(1,jr,0)+ &
												((capH(ir)+containers(520)+containers(141)*l(ir))**3+3*r(ir)* &
												(capH(ir)+containers(520)+containers(141)*l(ir))**2* &
												(r(ir)+containers(141)*differential1d(l,'x','cd',ir,ir,'x',r,r)))* &
												atemp(2,jr,0)
						end do
						residual1(ir) =			12*r(ir)*(containers(603)*(containers(501)+containers(502))+containers(40))			
					end do
					ir = nr
					jacobian1(ir,ir) =			1.0
					residual1(ir) =				0.0
					jacobian1_in = 				jacobian1
					residual1_in =				residual1
					call dgesv(nr,1,jacobian1_in,nr,ipiv1,residual1_in,nr,info)
					if (info .ne. 0) then
						if (termout .eq. 1)	then
							write(*,*)			'CAUTION: Brute force Reynolds equation solver - linear solution not obtained'
						end if
						write(95,*) 			'CAUTION: Brute force Reynolds equation solver - linear solution not obtained'
					end if
					phd = 						residual1_in
!					obtaining node-wise non-hydrodynamic pressure components and total pressure  
					do ir = 1,nr
						pdl(ir) =				containers(16)*exp(-containers(43)*(capH(ir)+containers(520)+containers(141)*l(ir)))
						pvdw(ir) =				-1.0/((capH(ir)+containers(520)+containers(141)*l(ir))**3)
						ps(ir) =				exp(-containers(44)*(capH(ir)+containers(520)+containers(141)*l(ir)))* &
												cos(containers(54)*(capH(ir)+containers(520)+containers(141)*l(ir))+phase)
						p(ir) =					containers(45)*phd(ir)+containers(46)*pdl(ir)+containers(47)*pvdw(ir)+ &
												containers(48)*ps(ir)
					end do
!					obtaining pressure in Hankel space from pressure in real-space
					do ix = 1,nx
						pcapH0(ix) =			0.0
						do ir = 1,nr
							pcapH0(ix) =		pcapH0(ix)+Hnkkrnlf(ix,ir)*p(ir)
						end do
					end do
!					obtaining Hankel-space deflection from Hankel-space pressure
					do ix = 1,nx
						lcapH0(ix) =			capX(ix)*pcapH0(ix)
					end do
!					inverting solution of deflection in Hankel-space to get solution of deflection in real-space
					do ir = 1,nr
						l(ir) =					0.0
						do ix = 1,nx
							l(ir) =				l(ir)+Hnkkrnli(ir,ix)*lcapH0(ix)
						end do
					end do
!					obtaining spring deflection
					if (ispring .ne. 0) then
						lcapSprev = 			lcapS(it)
						lcapS(it) = 			0.0
						ir = 1
						lcapS(it) = 			lcapS(it)+r(ir)*p(ir)*(r(ir+1)-r(ir))
						do ir = 2,nr-1
							lcapS(it) = 		lcapS(it)+r(ir)*p(ir)*(r(ir+1)-r(ir-1))
						end do
						ir = nr
						lcapS(it) = 			lcapS(it)+r(ir)*p(ir)*(r(ir)-r(ir-1))
						lcapS(it) = 			lcapS(it)*pi
					end if
!					obtaining residual error
					if (maxval(abs(lprev)) .eq. 0) then
						capRes = sum(abs(l-lprev))
					else
						capRes = sum(abs(l-lprev))/maxval(abs(lprev))
					end if
					if (maxval(abs(phdprev)) .eq. 0) then
						capRes = capRes + sum(abs(phd-phdprev))
						capRes = capRes/(2*(nr-1))
					else
						capRes = capRes + sum(abs(phd-phdprev))/maxval(abs(phdprev))
					end if
					if (ispring .ne. 0) then
						if (abs(lcapSprev) .eq. 0) then
							capRes = capRes + abs(lcapS(it)-lcapSprev)
						else
							capRes = capRes + abs(lcapS(it)-lcapSprev)/abs(lcapSprev)
						end if
					end if
					if (termout .eq. 1) then
						write(*,*)	't =', (t(it)), ',', (it), 'of', (nt), 'iter. no',(iter), &
									'subsvel. iter.',(itersubsvel), 'error',(capRes),'b.f.'
					end if
					write(94,*)		't = ', (t(it)), ',', (it), ' of ', (nt), ' iter. no ',(iter), &
									' subsvel. iter. ',(itersubsvel), ' error ',(capRes),' b.f.'
				end	do
			end do

		else

!			obtaining solution using scheme and method as specified in the variable initiation module
			capRes =	tol*2.0
			capReslast=	1.0e20
			iter =		0
			do while (capRes .ge. tol)
				iter =							iter+1
				lprev =							l
				phdprev =						phd
				containers(520) = 				containers(149)*lcapS(it)
!				obtaining node-wise non-hydrodynamic pressure components and total pressure  
				do ir = 1,nr
					pdl(ir) =					containers(16)*exp(-containers(43)*(capH(ir)+containers(520)+containers(141)*l(ir)))
					pvdw(ir) =					-1.0/((capH(ir)+containers(520)+containers(141)*l(ir))**3)
					ps(ir) =					exp(-containers(44)*(capH(ir)+containers(520)+containers(141)*l(ir)))* &
												cos(containers(54)*(capH(ir)+containers(520)+containers(141)*l(ir))+phase)
					p(ir) =						containers(45)*phd(ir)+containers(46)*pdl(ir)+containers(47)*pvdw(ir)+ &
												containers(48)*ps(ir)
				end do
!				obtaining derivative with deflection for node-wise non-hydrodynamic pressure components and total pressure  
				do ir = 1,nr
					dpdldl(ir) =	-containers(142)*pdl(ir)
					dpvdwdl(ir) =	-(containers(143)/(capH(ir)+containers(520)+containers(141)*l(ir)))*pvdw(ir)
					dpsdl(ir) =		-containers(144)* &
									(1.0+2*pi*tan(containers(54)*(capH(ir)+containers(520)+containers(141)*l(ir))+phase))*ps(ir)
					dpdl(ir) =		containers(46)*dpdldl(ir)+containers(47)*dpvdwdl(ir)+containers(48)*dpsdl(ir)
				end do
!				filling residual for Reynolds equation
				do ir = 2,nr-1
					if (it .eq. 1) then
						containers(501) = 		0.0
						containers(502) = 		0.0
					else
						if (it .eq. 2) then
							containers(501) =	(l(ir)-llast(ir))*containers(39)
							containers(502) =	(lcapS(it)-lcapS(it-1))*containers(39)
						else							
							containers(501) =	atemp(20,0,0)*l(ir)+atemp(20,-1,0)*llast(ir)+atemp(20,-2,0)*llast2(ir)
							containers(502) =	atemp(20,0,0)*lcapS(it)+atemp(20,-1,0)*lcapS(it-1)+atemp(20,-2,0)*lcapS(it-2)
						end if
					end if
					if (smallamp .eq. 0) then
						containers(501) = 		containers(373)* &
												(containers(501)+l(ir)*dndtbyn(it)- &
												0.5*r(ir)*differential1d(l,'x','cd',ir,ir,'x',r,r)*dhdtbyh(it))
						containers(502) = 		containers(374)*(containers(502)+lcapS(it)*dncapSdtbyncapS(it))
					else
						containers(501) = 		containers(373)*containers(501)
						containers(502) = 		containers(374)*containers(502)
					end if
					residual(ir) =				r(ir)*(capH(ir)+containers(520)+containers(141)*l(ir))**3* &
												differential1d(phd,'xx','cd',ir,ir,'x',r,r)+ &
												((capH(ir)+containers(520)+containers(141)*l(ir))**3+ &
												3*r(ir)*(capH(ir)+containers(520)+containers(141)*l(ir))**2* &
												(r(ir)+containers(141)*differential1d(l,'x','cd',ir,ir,'x',r,r)))* &
												differential1d(phd,'x','cd',ir,ir,'x',r,r)- &
												12*r(ir)*(containers(501)+containers(502)+containers(40))
				end do
!				filling jacobian for Reynolds equation
				do ir = 2,nr-1
					do ir1 = 1,nr+nx+1
						jacobian(ir,ir1) = 		0.0
					end do
				end do
				do ir = 2,nr-1
					atemp(1,:,:) =				a('','cd',ir,ir,r,r)
					atemp(2,:,:) =				a('x','cd',ir,ir,r,r)
					atemp(3,:,:) =				a('xx','cd',ir,ir,r,r)
					containers(503) =			capH(ir)+containers(520)+containers(141)*l(ir)
					containers(504) =			containers(503)**2
					containers(505) =			containers(503)**3
					containers(506) = 			r(ir)+containers(141)*differential1d(l,'x','cd',ir,ir,'x',r,r)
					containers(507) = 			r(ir)*containers(504)
					containers(508) = 			r(ir)*containers(505)
					containers(509) = 			containers(505)+3*containers(507)*containers(506)
					containers(510) = 			containers(507)*differential1d(phd,'xx','cd',ir,ir,'x',r,r)+ &
												(containers(504)+2*r(ir)*containers(503)*containers(506))* &
												differential1d(phd,'x','cd',ir,ir,'x',r,r)+r(ir)*containers(38)
					containers(511) = 			containers(143)*(containers(510)+r(ir)*containers(37))
					containers(512) = 			containers(143)* &
												(containers(507)*differential1d(phd,'x','cd',ir,ir,'x',r,r)+r(ir)**2*containers(35))
					do jr = -1,1
						jacobian(ir,ir+jr) =	containers(508)*atemp(3,jr,0)+containers(509)*atemp(2,jr,0)
						jacobian(ir,nr+ir+jr) =	containers(511)*atemp(1,jr,0)+containers(512)*atemp(2,jr,0)+containers(514)
					end do
					if (ispring .ne. 0) then
						jacobian(ir,nr+nx+1) =	containers(150)*(containers(510)+r(ir)*containers(36))
					end if
				end do
!				filling residual for hydrodynamic pressure centerline symmetry and far-end atmospheric conditions 
				ir = 1
				residual(ir) =								differential1d(phd,'x','fd',ir,ir,'x',r,r)
				ir = nr
				residual(ir) =								phd(ir)
!				filling residual for deflection-pressure Hankel-space relation
				do ix = 2,nx-1
					residual(nr+ix) =						0.0
					do ir = 1,nr
						residual(nr+ix) =					residual(nr+ix)+Hnkkrnlf(ix,ir)*(l(ir)-capX(ix)*p(ir))
					end do
				end do
				if (symbc .eq. 0) then
					ix = 1
					residual(nr+ix) =						0.0
					do ir = 1,nr
						residual(nr+ix) =					residual(nr+ix)+Hnkkrnlf(ix,ir)*(l(ir)-capX(ix)*p(ir))
					end do
				else
!				filling residual for symmetric deflection at centerline (if applied)
					ix = 1
					ir = ix
					residual(nr+ix) =						differential1d(l,'x','fd',ir,ir,'x',r,r)
				end if
!				filling jacobian for deflection for deflection-pressure Hankel-space relation
				do ix = 2,nx-1
					do ir = 1,nr
						jacobian(nr+ix,nr+ir) =				Hnkkrnlf(ix,ir)*(1.0-capX(ix)*dpdl(ir))
					end do
				end do
				if (symbc .eq. 0) then
					ix = 1
					do ir = 1,nr
						jacobian(nr+ix,nr+ir) =				Hnkkrnlf(ix,ir)*(1.0-capX(ix)*dpdl(ir))
					end do
				end if
!				filling residual for spring compression expression
				residual(nr+nx+1) = 						lcapS(it)/pi
				if (ispring .ne. 0) then
					do ir = 1,nr-1
						residual(nr+nx+1) = 				residual(nr+nx+1)-(p(ir+1)*r(ir+1)+p(ir)*r(ir))*(r(ir+1)-r(ir))
					end do
				end if
!				filling jacobian for hydrodynamic pressure and deflection for spring compression expression
				if (ispring .ne. 0) then
					ir = 1
					jacobian(nr+nx+1,ir) = 					r(ir)*(r(ir+1)-r(ir))
					jacobian(nr+nx+1,nr+ir) = 				jacobian(nr+nx+1,ir)
					do ir = 2,nr-1
						jacobian(nr+nx+1,ir) = 				r(ir)*(r(ir+1)-r(ir-1))
						jacobian(nr+nx+1,nr+ir) = 			jacobian(nr+nx+1,ir)
					end do
					ir = nr
					jacobian(nr+nx+1,ir) = 					r(ir)*(r(ir)-r(ir-1))
					jacobian(nr+nx+1,nr+ir) = 				jacobian(nr+nx+1,ir)
					do ir = 1,nr
						jacobian(nr+nx+1,ir) = 				-containers(45)*jacobian(nr+nx+1,ir)
						jacobian(nr+nx+1,nr+ir) = 			-dpdl(ir)*jacobian(nr+nx+1,nr+ir)
					end do
				end if
!				obtaining iterative solution
				jacobian_in = jacobian
				residual_in = residual
				call dgesv(nr+nx+1,1,jacobian_in,nr+nx+1,ipiv,residual_in,nr+nx+1,info)
				solution = 									residual_in
				do ir = 1,nr
					phd(ir) =								phd(ir)-relax*solution(ir)
					l(ir) =									l(ir)-relax*solution(nr+ir)
				end do
				if (ispring .ne. 0) then
					lcapS(it) = 							lcapS(it)-relax*solution(nr+nx+1)
				end if
				capRes = 0.0
				do ir = 1,nr+nx+1
					if (abs(residual(ir)) .ge. capRes) then 
						capRes =	abs(residual(ir))
					end if
				end do
				if (termout .eq. 1) then
					write(*,*)	't =', (t(it)), ',', (it), 'of', (nt), 'iter. no',(iter), &
								'error',(capRes),(scheme)
				end if
				write(94,*)		't =', (t(it)), ',', (it), 'of', (nt), 'iter. no',(iter), &
								'error',(capRes),(scheme)
				if ((capRes .gt. capReslast) .and. (capReslast .le. tolbreak)) then
					do ir = 1,nr
						phd(ir) =							phd(ir)+relax*solution(ir)
						l(ir) =								l(ir)+relax*solution(nr+ir)
					end do	
					if (ispring .ne. 0) then
						lcapS(it) = 						lcapS(it)+relax*solution(nr+nx+1)
					end if
					capRes = 0.0
					if (termout .eq. 1) then
						write(*,*)	'Warning: solution diverging but error in allowed range -> breaking out of loop' 
					end if
					write(94,*)		'Warning: solution diverging but error in allowed range -> breaking out of loop'
				end if
				capReslast = 	capRes
			end do
		end if
	end if

!	obtaining other field variables (deformation field)
	if (.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) then 
		if (oneway(it) .eq. 0) then 
			do ix = 1,nx
				pcapH0(ix) =			0.0
				do ir = 1,nr
					pcapH0(ix) =		pcapH0(ix)+Hnkkrnlf(ix,ir)*p(ir)
				end do
			end do
		end if
!		generating y-grid
		if (smallamp .eq. 0) then
			ymin = 		0.0
			ymax =						(1.0/d(it))*(beta/delta)
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
		do ix = 1,nx
			containers(409) =			containers(351)*(pcapH0(ix)/(2*x(ix)**2))
			containers(410) =			(containers(409)*containers(352))/x(ix)
			capAcoeff(ix) =				(containers(405)/containers(404))*containers(412)*containers(410)
			capBcoeff(ix) =				-(containers(406)/containers(404))*containers(412)*containers(409)
			capCcoeff(ix) =				-(containers(407)/containers(404))*containers(411)*containers(410)
			capDcoeff(ix) =				(containers(408)/containers(404))*containers(411)*containers(409)
		end do
		do ix = 1,nx
			do iy = 1,ny
				gcapH(ix,iy) =			(capAcoeff(ix)+capBcoeff(ix)*y(iy))*exp(containers(401)*y(iy)) + &
										(capCcoeff(ix)+capDcoeff(ix)*y(iy))*exp(-containers(401)*y(iy))
			end do
		end do
		do ix = 1,nx
			iy = 1
			dgcapHdy(ix,iy) = 			differential2d (gcapH,'y','fd',ix,iy,'xy',x,y)		
			do iy = 2,ny-1
				dgcapHdy(ix,iy) =		differential2d (gcapH,'y','cd',ix,iy,'xy',x,y)
			end do
			iy = ny
			dgcapHdy(ix,iy) = 			differential2d (gcapH,'y','bd',ix,iy,'xy',x,y)
		end do
		do ix = 1,nr
			iy = 1
			d2gcapHdy2(ix,iy) = 		differential2d (dgcapHdy,'y','fd',ix,iy,'xy',x,y)
			do iy = 2,ny-1
				d2gcapHdy2(ix,iy)=		differential2d (dgcapHdy,'y','cd',ix,iy,'xy',x,y)
			end do
			iy = ny
			d2gcapHdy2(ix,iy) = 		differential2d (dgcapHdy,'y','bd',ix,iy,'xy',x,y)
		end do
		do ix = 1,nx
			do iy = 1,ny
				urcapH1(ix,iy) =		containers(355)*x(ix)*dgcapHdy(ix,iy)
				uycapH0(ix,iy) =		d2gcapHdy2(ix,iy)-containers(356)*x(ix)**2*gcapH(ix,iy)
			end do
		end do
		do iy = 1,ny
			ur(ir,iy) = 			0.0
			uy(ir,iy) =				0.0
			do ix = 1,nx
				ur(ir,iy) =			ur(ir,iy)+Hnkkrnli1(ir,ix)*urcapH1(ix,iy)
				uy(ir,iy) = 		uy(ir,iy)+Hnkkrnli(ir,ix)*uycapH0(ix,iy)
			end do
			ur =					((rmax**2)/Bslroot1(nr))*ur
			uy =					((rmax**2)/Bslroot1(nr))*ur
		end do
	end if

!	obtaining other field variables (flow field)
	if (.not. ((minout .eq. 2) .or. (minout .eq. 4) .or. (minout .eq. 5))) then
		do ir = 1,nr
			zmax =					capH(ir)
			z(ir,1) =				-containers(361)*l(ir)
			zstep =					(capH(ir)-z(ir,1))/(nz-1)
			do iz = 2,nz
				z(ir,iz) =			z(ir,iz-1)+zstep
			end do
		end do
		ir = 1
		dldr(ir) = 					differential1d(l,'x','fd',ir,ir,'x',r,r)
		do ir = 2,nr-1	
			dldr(ir) = 				differential1d(l,'x','cd',ir,ir,'x',r,r)
		end do
		ir = nr
		dldr(ir) = 					differential1d(l,'x','bd',ir,ir,'x',r,r)
		ir = 1
		dphddr(ir) = 				differential1d(phd,'x','fd',ir,ir,'x',r,r)
		do ir = 2,nr-1
			dphddr(ir) = 			differential1d(phd,'x','cd',ir,ir,'x',r,r)
		end do
		ir = nr
		dphddr(ir) = 				differential1d(phd,'x','bd',ir,ir,'x',r,r)
		ir = 1
		d2phddr2(ir) = 				differential1d(dphddr,'x','fd',ir,ir,'x',r,r)
		do ir = 2,nr-1
			d2phddr2(ir) = 			differential1d(phd,'xx','cd',ir,ir,'x',r,r)
		end do
		ir = nr
		d2phddr2(ir) = 				differential1d(dphddr,'x','bd',ir,ir,'x',r,r)
		if (it .eq. 1) then
			containers(501) = 		0.0
			containers(512) = 		0.0
		else
			if (it .eq. 2) then
				containers(502) =	(lcapS(it)-lcapS(it-1))*containers(146)	
			else
				containers(502) =	atemp(20,0,0)*lcapS(it)+atemp(20,-1,0)*lcapS(it-1)+ &
									atemp(20,-2,0)*lcapS(it-2)
			end if
		end if
		if (smallamp .eq. 0) then
			containers(502) = 		containers(374)*(containers(502)+lcapS(it)*dncapSdtbyncapS(it))
		else
			containers(502) = 		containers(374)*containers(502)
		end if
		do ir = 1,nr
			do iz = 1,nz
				vr(ir,iz) =			0.5*dphddr(ir)* &
									(z(ir,iz)**2-(capH(ir)+containers(149)*lcapS(it)-containers(141)*l(ir))*z(ir,iz)- &
									containers(141)*l(ir)*(capH(ir)+containers(149)*lcapS(it)))
				vz(ir,iz) =			(containers(502)+containers(40))-(1.0/12.0)*(d2phddr2(ir)+(1.0/r(ir))*dphddr(ir))* &
									(2.0*(z(ir,iz)**3-(capH(ir)+containers(149)*lcapS(it))**3)- &
									3.0*(capH(ir)+containers(149)*lcapS(it)-containers(141)*l(ir))* &
									(z(ir,iz)**2-(capH(ir)+containers(149)*lcapS(it))**2)- &
									6.0*containers(141)*l(ir)*(capH(ir)+containers(149)*lcapS(it))* &
									(z(ir,iz)-(capH(ir)+containers(149)*lcapS(it)))) - &
									(1.0/4.0)*dphddr(ir)* &
									(r(ir)*(capH(ir)+containers(149)*lcapS(it))**2.0+ &
									2.0*containers(141)*r(ir)*l(ir)*(capH(ir)+containers(149)*lcapS(it))- &
									r(ir)*z(ir,iz)*(z(ir,iz)+2.0*containers(141)*l(ir))- &
									2*containers(141)*z(ir,iz)*(capH(ir)+containers(149)*lcapS(it))*dldr(ir)+ &
									containers(141)*z(ir,iz)**2.0*dldr(ir))
			end do
		end do
		ir = 1
		do iz = 1,nz
			vz(ir,iz) =				(containers(502)+containers(40))-(1.0/12.0)*d2phddr2(ir)* &
									(2.0*(z(ir,iz)**3-(capH(ir)+containers(149)*lcapS(it))**3)- &
									3.0*(capH(ir)+containers(149)*lcapS(it)-containers(141)*l(ir))* &
									(z(ir,iz)**2-(capH(ir)+containers(149)*lcapS(it))**2)- &
									6.0*containers(141)*l(ir)*(capH(ir)+containers(149)*lcapS(it))* &
									(z(ir,iz)-(capH(ir)+containers(149)*lcapS(it))))
		end do
	end if

!	obtaining force response
	capF(it) =						0.0
	do ir = 2,nr
		capF(it) =					capF(it)+ &
									0.5*(r(ir-1)*p(ir-1)+r(ir)*p(ir))*(r(ir)-r(ir-1))
	end do
	capF(it) =						2*pi*capF(it)
	capF_d(it) =					containers(148)*capF(it)

end subroutine linelast
!-----------------------------------------------------------------------------------------------------------------------------------
