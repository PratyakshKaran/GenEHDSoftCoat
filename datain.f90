!-----------------------------------------------------------------------------------------------------------------------------------
! data import subroutine
!
! imports variable values from input file (see path convention in the subroutine for the environment)
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
subroutine datain

	use varinit
	implicit none

!	importing solvent properties
	open(unit=3,file='../../../../in/ttsi/solvent.csv',status='old',action='read')
	read(3,*)
	read(3,*)	rho, mu, sigma
	close(3)

!	importing electrolyte properties
	open(unit=1,file='../../../../in/ttsi/electrolyte.csv',status='old',action='read')
	read(1,*)
	read(1,*)	nspecies
	close(1)
	allocate(diffcoeff(1:nspecies))
	allocate(valence(1:nspecies))
	open(unit=1,file='../../../../in/ttsi/electrolyte.csv',status='old',action='read')
	read(1,*)
	read(1,*)	nspecies, n0, relperm, (diffcoeff(itemp), itemp=1,nspecies), &
				(valence(itemp), itemp=1,nspecies)
	close(1)

!	importing solid properties
	open(unit=2,file='../../../../in/ttsi/solid.csv',status='old',action='read')
	allocate(psis(2))
	read(2,*)
	read(2,*)	subsmodul1, subsmodul2, capAsfw, capA, phase, model, iscapLame, psis(1), psis(2)
	close(2)
	open(unit=2,file='../../../../in/ttsi/prony.csv',status='old',action='read')
	read(2,*)
	nprony = 0
	do
		read(2,*,iostat=itemp)
		if (itemp .ne. 0) exit
		nprony = nprony + 1
	end do
	close(2)	
	allocate(capKbulkprony(1:nprony))
	allocate(capGprony(1:nprony))
	allocate(taukprony(1:nprony))
	allocate(taugprony(1:nprony))
	open(unit=2,file='../../../../in/ttsi/prony.csv',status='old',action='read')
	read(2,*)
	do iprony = 1,nprony
		read(2,*)	capKbulkprony(iprony), capGprony(iprony), taukprony(iprony), taugprony(iprony)
	end do
	close(2)	

!	importing universal constants
	open(unit=7,file='../../../../in/ttsi/univ.csv',status='old',action='read')
	read(7,*)
	read(7,*)	elementaryq, kcapB, permvac
	close(7)

!	importing geometry parameteres
	open(unit=4,file='../../../../in/ttsi/geometry.csv',status='old',action='read')
	read(4,*)
	read(4,*)	h0, capD, capL, capR, omega, capT, iappr, iforceharmonicve, ispring, capKspring
	close(4)

!	importing grid parameteres
	open(unit=5,file='../../../../in/ttsi/grid.csv',status='old',action='read')
	read(5,*)
	read(5,*)	ady, adt, adtau
	read(5,*)
	read(5,*)	mady, nady, madt, nadt, madtau, nadtau
	read(5,*)
	read(5,*)	scheme, symbc, guess, revtime, considerlim
	read(5,*)
	read(5,*)	minout, termout
	read(5,*)
	read(5,*)
	read(5,*)	tol, tolbreak, thressmallamp, thresoneway1, thresoneway2, relax, thresadhesion, threscavity, threslin
	read(5,*)	thresstokes, nsubsvel, thresthin, thressemi, threspcompr, thresincompr, nudummy, thresdomhd, thresreduce
	read(5,*)
	read(5,*)	nt, nr, ny, nz, ntau
	read(5,*)
	read(5,*)	tstart, tend, rinf, tauinf 		
	close(5)

! 	obtaining Bessel function roots
	allocate(Bslroot(1:nr))
	open(unit=8,file='../../../supp/besselroots_1_0.dat',status='old',action='read')
	read(8,*)	(Bslroot(ir), ir=1,nr)
	close(8)
	if ((.not. ((minout .eq. 3) .or. (minout .eq. 4) .or. (minout .eq. 5))) .and. (scheme .eq. 0)) then	
		allocate(Bslroot1(1:nr))
		open(unit=9,file='../../../supp/besselroots_1_1.dat',status='old',action='read')
		read(9,*)	(Bslroot1(ir), ir=1,nr)
		close(9)
	end if
		

!	importing problem prob
	open(unit=60,file='../../../../in/ttsi/type.csv',status='old',action='read')
	read(60,'(A)')	probtype
	close(60)

!	getting timestamp
	call idate(today)
	call itime(now)
	write(datestamp,10)			today(3),today(2),today(1)
	write(siminstance,20)		now(1),now(2),now(3)
	siminstance = trim(probtype)//"_"//siminstance

	call execute_command_line(	"mkdir -p ../../../../out/ttsi/"//datestamp//"/"//trim(siminstance),exitstat=itemp)

	10 format (I4,I2.2,I2.2)
	20 format (I2.2,I2.2,I2.2)

end subroutine datain
!-----------------------------------------------------------------------------------------------------------------------------------
