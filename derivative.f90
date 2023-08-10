!-----------------------------------------------------------------------------------------------------------------------------------
! derivative functions module

! provides functions for obtaining derivatives as well as obtaining co-efficients for discretized FD equations (non-uniform grid)
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
module derivative

    implicit none

    contains

!		function defining coefficients of the appropriate boundaries for different derivatives
!		arguments are: diff - differentiation type 'x', 'xx', 'y', 'yy', 'xy'; dirr - direction of derivative 'fd', 'cd', 'bd', 
!		with only 'cd' being applicable for 2nd order derivatives i.e. 'xx', 'yy', 'xy'; ix - index of first coordinate x where 
!		derivative is to be obtained; iy - index of second coordinate y where derivative is to be obtained; x - x-grid; y - y-grid
        function a(diff,dir,ix,iy,x,y)
            double precision :: a(-2:2,-2:2)
            character (len=*), intent(in) :: diff, dir
            integer, intent(in) :: ix, iy
            double precision, intent(in) :: x(:), y(:)
            double precision :: iidx, idx, dx, dxi, dxii
            double precision :: jjdy, jdy, dy, dyj, dyjj
            a = 0.0
            select case (trim(diff))
            case ('')
                a(0,0) = 1.0
            case ('x')
                select case (dir)
                case ('fd')
                    dx =        x(ix+1)-x(ix)
                    dxi =       x(ix+2)-x(ix+1)
                    a(2,0) =    -(dx**2.0)/(dxi*dx*(dx+dxi))
                    a(1,0) =    ((dx+dxi)**2.0)/(dxi*dx*(dx+dxi))
                    a(0,0) =    -((dx+dxi)**2.0-dx**2.0)/(dxi*dx*(dx+dxi))
                case ('bd')
                    idx =       x(ix)-x(ix-1)
                    iidx =      x(ix-1)-x(ix-2)
                    a(0,0) =    ((idx+iidx)**2.0-idx**2.0)/(iidx*idx*(idx+iidx))
                    a(-1,0) =   -((idx+iidx)**2.0)/(iidx*idx*(idx+iidx))
                    a(-2,0) =   (idx**2.0)/(iidx*idx*(idx+iidx))
                case ('cd')
                    dx =        x(ix+1)-x(ix)
                    idx =       x(ix)-x(ix-1)
                    a(1,0) =    idx/(2.0*dx*idx)
                    a(0,0) =    (dx-idx)/(2.0*dx*idx)
                    a(-1,0) =   -dx/(2.0*dx*idx)
                end select
            case ('y')
                select case (dir)
                case ('fd')
                    dy =        y(iy+1)-y(iy)
                    dyj =       y(iy+2)-y(iy+1)
                    a(0,2) =    -(dy**2.0)/(dyj*dy*(dy+dyj))
                    a(0,1) =    ((dy+dyj)**2.0)/(dyj*dy*(dy+dyj))
                    a(0,0) =    -((dy+dyj)**2.0-dy**2.0)/(dyj*dy*(dy+dyj))
                case ('bd')
                    jdy =       y(iy)-y(iy-1)
                    jjdy =      y(iy-1)-y(iy-2)
                    a(0,0) =    ((jdy+jjdy)**2.0-jdy**2.0)/(jjdy*jdy*(jdy+jjdy))
                    a(0,-1) =   -((jdy+jjdy)**2.0)/(jjdy*jdy*(jdy+jjdy))
                    a(0,-2) =   (jdy**2.0)/(jjdy*jdy*(jdy+jjdy))
                case ('cd')
                    dy =        y(iy+1)-y(iy)
                    jdy =       y(iy)-y(iy-1)
                    a(0,1) =    jdy/(2.0*dy*jdy)
                    a(0,0) =    (dy-jdy)/(2.0*dy*jdy)
                    a(0,-1) =   -dy/(2.0*dy*jdy)
                end select
            case ('xx')
                a(-1,0) = 		2.0/((x(ix)-x(ix-1))*(x(ix+1)-x(ix-1)))
                a(1,0) = 		2.0/((x(ix+1)-x(ix))*(x(ix+1)-x(ix-1)))
                a(0,0) = 		-2.0/((x(ix+1)-x(ix))*(x(ix)-x(ix-1)))
            case ('yy')
                a(0,-1) = 		2.0/((y(iy)-y(iy-1))*(y(iy+1)-y(iy-1)))
                a(0,1) = 		2.0/((y(iy+1)-y(iy))*(y(iy+1)-y(iy-1)))
                a(0,0) = 		-2.0/((y(iy+1)-y(iy))*(y(iy)-y(iy-1)))
            case ('xy')
                dx =        	x(ix+1)-x(ix)
                idx =      		x(ix)-x(ix-1)
                dy =        	y(iy+1)-y(iy)
                jdy =       	y(iy)-y(iy-1)
                a(1,1) =    	(idx*jdy)/(4.0*idx*dx*jdy*dy)
                a(1,0) =    	(idx*(dy-jdy))/(4.0*idx*dx*jdy*dy)
                a(0,1) =    	((dx-idx)*jdy)/(4.0*idx*dx*jdy*dy)
                a(1,-1) =   	(-idx*dy)/(4.0*idx*dx*jdy*dy)
                a(0,0) =    	((dx-idx)*(dy-jdy))/(4.0*idx*dx*jdy*dy)
                a(-1,1) =   	(-dx*jdy)/(4.0*idx*dx*jdy*dy)
                a(0,-1) =   	(-(dx-idx)*dy)/(4.0*idx*dx*jdy*dy)
                a(-1,0) =   	(-dx*(dy-jdy))/(4.0*idx*dx*jdy*dy)
                a(-1,-1) =  	(dx*dy)/(4.0*idx*dx*jdy*dy)
            end select
        end function a

!		function computing derivative of variable defined over only one coordinate
!		arguments are: var - variable of which derivative is computed; diff - differentiation type 'x', 'xx', 'y', 'yy'; dirr - 
!		direction of derivative 'fd', 'cd', 'bd', with only 'cd' being applicable for 2nd order derivatives i.e. 'xx', 'yy'; ix - 
!		index of first coordinate x where derivative is to be obtained; iy - index of second coordinate y where derivative is to be 
!		obtained; dimn - coordinate for which the variable is defined; x - x-grid; y - y-grid
        function differential1d (var,diff,dir,ix,iy,dimn,x,y)
            double precision :: differential1d
            double precision, dimension(-2:2,-2:2) :: acoeff
            integer :: jx, jy
            character(len=*), intent(in) :: diff, dir, dimn
            double precision, intent(in) :: x(:), y(:)
            integer, intent(in) :: ix,iy
            double precision, intent(in) :: var(:)
            differential1d = 0.0
            if ((((diff=='x').or.(diff=='xx')).and.(dimn=='x')).or.(((diff=='y').or.(diff=='yy')).and.(dimn=='y'))) then
                acoeff = a(diff,dir,ix,iy,x,y)
                if (dimn == 'x') then
                    do jx = -2,2
                        jy = 0
                        if (acoeff(jx,jy) .ne. 0.0) then
                            differential1d = differential1d + acoeff(jx,jy)*var(ix+jx)
                        end if
                    end do
                else
                    jx = 0
                    do jy = -2,2
                        if (acoeff(jx,jy) .ne. 0.0) then
                            differential1d = differential1d + acoeff(jx,jy)*var(iy+jy)
                        end if
                    end do
                end if
            end if
        end function differential1d

!		function computing derivative of complex variable defined over only one coordinate (same parameters as for real)
        function differential1dcmplx (var,diff,dir,ix,iy,dimn,x,y)
           complex*16 :: differential1dcmplx
            double precision, dimension(-2:2,-2:2) :: acoeff
            integer :: jx, jy
            character(len=*), intent(in) :: diff, dir, dimn
            double precision, intent(in) :: x(:), y(:)
            integer, intent(in) :: ix,iy
            complex*16, intent(in) :: var(:)
            differential1dcmplx = 0.0
            if ((((diff=='x').or.(diff=='xx')).and.(dimn=='x')).or.(((diff=='y').or.(diff=='yy')).and.(dimn=='y'))) then
                acoeff = a(diff,dir,ix,iy,x,y)
                if (dimn == 'x') then
                    do jx = -2,2
                        jy = 0
                        if (acoeff(jx,jy) .ne. 0.0) then
                            differential1dcmplx = differential1dcmplx + acoeff(jx,jy)*var(ix+jx)
                        end if
                    end do
                else
                    jx = 0
                    do jy = -2,2
                        if (acoeff(jx,jy) .ne. 0.0) then
                            differential1dcmplx = differential1dcmplx + acoeff(jx,jy)*var(iy+jy)
                        end if
                    end do
                end if
            end if
        end function differential1dcmplx

!		function computing derivative of variable defined over both coordinates
!		arguments are: var - variable of which derivative is computed; diff - differentiation type 'x', 'xx', 'y', 'yy', 'xy'; 
!		dirr - direction of derivative 'fd', 'cd', 'bd', with only 'cd' being applicable for 2nd order derivatives i.e. 'xx', 'yy'
!		and 'xy'; ix - index of first coordinate x where derivative is to be obtained; iy - index of second coordinate y where 
!		derivative is to be obtained; dimn - dummy variable with 'xy' the only applicable option; x - x-grid; y - y-grid
        function differential2d (var,diff,dir,ix,iy,dimn,x,y)
            double precision :: differential2d
            double precision, dimension(-2:2,-2:2) :: acoeff
            integer :: jx, jy
            character(len=*), intent(in) :: diff, dir, dimn
            double precision, intent(in) :: x(:), y(:)
            integer, intent(in) :: ix,iy
            double precision, intent(in) :: var(:,:)
            differential2d = 0.0
            acoeff = a(diff,dir,ix,iy,x,y)
            do jx = -2,2
                do jy = -2,2
                    if (acoeff(jx,jy) .ne. 0.0) then
                        differential2d = differential2d + acoeff(jx,jy)*var(ix+jx,iy+jy)
                    end if
                end do
            end do
        end function differential2d

!		function computing derivative of complex variable defined over both coordinates (same parameters as for real)
        function differential2dcmplx (var,diff,dir,ix,iy,dimn,x,y)
            complex*16 :: differential2dcmplx
            double precision, dimension(-2:2,-2:2) :: acoeff
            integer :: jx, jy
            character(len=*), intent(in) :: diff, dir, dimn
            double precision, intent(in) :: x(:), y(:)
            integer, intent(in) :: ix,iy
            complex*16, intent(in) :: var(:,:)
            differential2dcmplx = 0.0
            acoeff = a(diff,dir,ix,iy,x,y)
            do jx = -2,2
                do jy = -2,2
                    if (acoeff(jx,jy) .ne. 0.0) then
                        differential2dcmplx = differential2dcmplx + acoeff(jx,jy)*var(ix+jx,iy+jy)
                    end if
                end do
            end do
        end function differential2dcmplx

!		DESCRIPTION TBD
        function itoa(i) result(res)
            character(:),allocatable :: res
            integer,intent(in) :: i
            character(range(i)+2) :: tmp
            write(tmp,'(i0)') i
            res = trim(tmp)
        end function itoa
		
!		DESCRIPTION TBD
		function jrange(ix,iy,nx,ny,diffrange)		
			integer, intent (in) :: ix,iy,nx,ny,diffrange
			integer :: jrange(1:2,1:2)
			if (diffrange .eq. 3) then
				if (ix .eq. 1) then
					jrange(1,1) = 0
					jrange(1,2) = 2
				else
					if (ix .eq. nx) then
						jrange(1,1) = -2
						jrange(1,2) = 0
					else
						jrange(1,1) = -1
						jrange(1,2) = 1
					end if
				end if
				if (iy .eq. 1) then
					jrange(2,1) = 0
					jrange(2,2) = 2
				else
					if (iy .eq. ny) then
						jrange(2,1) = -2
						jrange(2,2) = 0
					else
						jrange(2,1) = -1
						jrange(2,2) = 1
					end if
				end if			
			end if
			if (diffrange .ne. 3) then
				write(*,*) "	schematics yet to be drawn for other-ranged derivatives, &
								program as per three-ranged derivatives"
				jrange = 0
			end if			
		end function jrange

end module derivative
!-----------------------------------------------------------------------------------------------------------------------------------
