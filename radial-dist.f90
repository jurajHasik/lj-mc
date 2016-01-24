program radialDist
implicit none

real*8, parameter :: Pi     = 3.14159265359
integer, parameter :: conf_u = 30
integer, parameter :: out_u = 30
character(len=*), parameter :: out_F = 'rho_r.dat'

integer :: N, bins, smpl, bind
real*8 :: overflow
! In LJ units
real*8 :: rL, rDens, rMax, binRes, rBinD, rBinU
real*8, allocatable, dimension(:,:) :: r
real*8, allocatable, dimension(:) :: rho

! ***** Arguments *****
integer :: numArg, IARGC
character(len=20) :: argConf, argSamples, argRMax, argNBins, argN

integer :: i,j,a1,a2
real*8, dimension(3) :: dr
real*8 :: d

! ***** read in variables *****
numArg = IARGC()
if(numArg .ne. 5) then
    ! Unique id for output & input files
    write(*,'("Invalid Number of Arguments")')
    call EXIT(0)
endif
call GETARG(1, argConf)     ! File with configurations
call GETARG(2, argN)        ! #atoms
call GETARG(3, argSamples)  ! Number of samples
call GETARG(4, argRMax) 	! Maximal value for arg of pair dist f
call GETARG(5, argNBins)    ! Number of Bins
read(argN, *) N
read(argSamples, *) smpl
read(argRMax, *) rMax
read(argNBins, *) bins

allocate(r(3,N),rho(bins))
binRes = rMax/real(bins)


overflow = 0
rho = 0.0d0
open(unit=conf_u, file=trim(argConf), form='formatted',&
        status='old', action='read')
do i=1, smpl
	write(*,'("Processing sample: ",I10)') i
	call readConf(conf_u, N, r, rL)
	! Loop over pairs
	do a1=1, N-1
	do a2=a1+1, N
		dr = r(:,a2)-r(:,a1)
		dr = dr - rL*anint(dr/rL)
		d = sqrt(DOT_PRODUCT(dr,dr))
		if(d .lt. rMax) then
			bind = 1 + floor(d/binRes)
			rho(bind)=rho(bind)+2.0d0
		else
			overflow = overflow+1.0d0
		endif
	end do
	end do
enddo
close(conf_u)

! Take average
rho = rho/real(N*smpl)
overflow = overflow/real(N*smpl)
write(*,'("Overflow: ",1f10.5)') overflow
! Normalize with respect to id.gas 
! In the case of NVT ensemble rL is constant
rDens = N/(rL*rL*rL)
open(unit=out_u, file=out_F, form='formatted',&
        status='new', action='write')
do i=1, bins
	rBinD = (i-1)*binRes
	rBinU = i*binRes
	rho(i)=rho(i)/(4.0d0/3.0d0*Pi*rDens*(rBinU**3.0d0-rBinD**3.0d0))
	write(out_u,'(2f10.5)') (rBinD+(rBinU-rBinD)/2.0d0), rho(i)
enddo
close(out_u)

end program radialDist

subroutine readConf(in_unit, N, r, rL)
    implicit none

    Integer, intent(in) :: in_unit
    Integer, intent(in) :: N
    real*8, intent(out) :: rL
    Real*8, dimension(3,N), intent(out) :: r

    Integer :: i, nn

    ! Configuration is stored in .xyz formatted file
    ! Read simulation box dimensions from comment
    read(in_unit,*) nn
    if(nn .ne. N) then
        write(*,'("ERROR - Unexpected #Atoms: ",I10)') nn
        stop
    endif
    read(in_unit,'(4X,1F10.5)') rL
    do i=1, N
        read(in_unit,'(2X,3(1X,F20.10))') r(1,i), r(2,i), r(3,i)
    enddo
end subroutine readConf