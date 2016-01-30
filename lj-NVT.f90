program ljNVT

use mersenne_twister

implicit none

! ##### VARIABLES #####################################################
! ***** Sim variables *****
real*8, parameter :: Pi     = 3.14159265359
real*8, parameter :: kToEV  = 8.617342301d-5 ! [eV.K^-1]
! John A. White "Lennard-Jones as a model for argon 
! and test of extended renormalization group calculations",
! Journal of Chemical Physics 111 pp. 9352-9356 (1999)
! epsilon/k_b = 125.7  [K]
! sigma       = 0.3345 [nm]
real*8, parameter :: sig = 1.0d0 ! [Ang]
real*8, parameter :: eps = 1.0d0 ! [K]
integer :: N, A
real*8 :: temp, rT, rB, rD, rDens
real*8 :: L, rL, rV, rCut, rCut6

integer :: accSAmv ! Acceptance of single atom moves
real*8, parameter :: accRtSAmv = 0.35d0 ! Desired acceptance rate
! Frequency of updating SAmv delta in [MC steps]
integer, parameter :: SAmvAdjstRt = 100
integer :: eqS, prS
integer :: smplRt

! ***** Atom coords *****
real*8, allocatable, dimension(:,:) :: r
real*8 :: e6, e12 ! holds contribution from r^-6 & r^-12 part
! Long range corrections to energy and virial
real*8 :: eLRC6, eLRC12, vLRC6, vLRC12
real*8 :: eTot, vir, pres
real*8 :: rate

! ***** PRNG variables *****
integer :: seed
real, dimension(4) :: rnd4
real, dimension(10) :: rnd10
real, allocatable, dimension(:) :: rndN

! ***** File output *****
character(len=20), parameter :: initSc = "init-sc.xyz"
character(len=*), parameter :: conf   = "conf.xyz"
integer, parameter :: genU  = 20
integer, parameter :: confU = 21

! ***** Arguments *****
integer :: numArg, IARGC
character(len=20) :: argUuid, argN, argT, argL, argEqS, argPrS
character(len=20) :: argSmplRt, argSeed

! ***** Utility variables *****
integer :: i, k
real*8 :: t0, t1

! ***** Functions *****
real*8 :: virial

! ##### SETUP SIMULATION ##############################################
! ***** read in variables *****
numArg = IARGC()
if(numArg .ne. 8) then
    ! Unique id for output & input files
    write(*,'("Invalid Number of Arguments")')
    call EXIT(0)
endif
call GETARG(1, argUuid)
call GETARG(2, argN)
call GETARG(3, argT) ! read in temperature in Kelvins
call GETARG(4, argL) ! read in position in Angstrom 
call GETARG(5, argEqS)
call GETARG(6, argPrS)
call GETARG(7, argSmplRt)
call GETARG(8, argSeed)
read(argN, *) N
read(argT, *) temp
read(argL, *) L
read(argEqS, *) eqS
read(argPrS, *) prS
read(argSmplRt, *) smplRt
read(argSeed, *) seed

allocate(r(3,N), rndN(N))
rL = L/sig
rD = rL/10.0d0 ! Set initial delta to rL/10
! cutoff fixed to rL/2
rCut = rL/2.0d0
! we actually only care about a rCut^2, rCut^6
rCut6 = rCut**6.0d0
rV = rL*rL*rL
rT = temp/eps ! Get temp to reduced units
rB = 1.0d0/rT
rDens = N/rV

write(*,'("N",19X,I10)') N
!write(*,'("T [K]",15X,1f10.5)') temp/kToEV
!write(*,'("T [eV]",15X,1f10.5)') temp
!write(*,'("Beta [1/eV]",9X,1f10.5)') beta
!write(*,'("L [Ang]",13X,1f10.5)') L
!write(*,'("delta [Ang]",9X,1f10.5)') delta
!write(*,'("sigma [Ang]",9X,1f10.5)') sig
!write(*,'("epsilon [eV]",8X,1f10.5)') eps
write(*,'("REDUCED UNITS")')
write(*,'("rL",18X,1f10.5)') rL
write(*,'("rCut",16X,1f10.5)') rCut
write(*,'("rD",18X,1f10.5)') rD
write(*,'("rT",18X,1f10.5)') rT 
write(*,'("rBeta",15X,1f10.5)') rB
write(*,'("rDens",15X,1f10.5)') rDens

! ***** Setup initial sc arrangement *****
call setup_sc(N, r, rL)
!call print_InitR(N, r, initSc)
write(*,'("INITIAL SETUP INTO SC LATTICE DONE")')
! compute initial energy contributions e6 & e12
call ljTot(N, r, rL, rCut6, e6, e12)
! compute energy LRC
eLRC6  = -8.0d0*Pi*N*(N/rV)*(1.0d0/(3.0d0*(rCut**3.0d0)))
eLRC12 = 8.0d0*Pi*N*(N/rV)*(1.0d0/(9.0d0*(rCut**9.0d0)))
vLRC6  = 2.0d0*eLRC6
vLRC12 = 4.0d0*eLRC12
eTot = e12+e6 + (eLRC12+eLRC6) ! No double counting coming from ljTot
write(*,'("init E",14X,1f10.5)') eTot
write(*,'("E_lrc: ",1f10.5," V_rlc: ",1f10.5)') eLRC12+eLRC6, vLRC12+vLRC6

! Initialize MT PRNG
call random_setseed(seed)
write(*,'("Seed",16X,I10)') seed
write(*,'("-sample output-")')
call random_number(rnd10)
do i=1, 10
    write(*,*) rnd10(i)
enddo

accSAmv = 0
call CPU_TIME(t0)
do i=1, eqS
    ! Perform MC Step - attempt to Move all @N atoms
    call random_number(rndN)
    do k=1, N
        call random_number(rnd4)
        !Always select random atom instead of keeping fixed order
        !call SAmv(N, k, r, rD, rL, rB, lj(1,k), lj(2,k), rnd4, accSAmv)
        A = min(1 + floor(N*rndN(k)),N) 
        ! floor(N*rnd) maps to 0, N-1 -> in exceptionally rare cases to 0-N
        call SAmv(N, A, r, rD, rL, rB, rCut6, e6, e12, rnd4, accSAmv)
    enddo
    ! Take EQ values
    if(mod(i,smplRt) .eq. 0) then
        call CPU_TIME(t1)
        eTot = e12+e6 + (eLRC12+eLRC6)
        vir  = virial(N, r, rL, rCut)
        pres = (dble(N)/rB+vir+vLRC12+vLRC6)/rV !Eq of state + correction
        write(*,'("EQ STEP: ",I10," E: ",1f20.10," P: ",1f20.10,"&
            & time: ",1f10.5)') i, eTot, pres, t1-t0
        call CPU_TIME(t0)
    endif
    ! For gas with low density, nearly all moves 
    ! are accepted regardless of size of delta due to PBC
    if(mod(i,SAmvAdjstRt) .eq. 0) then
        rate = dble(accSAmv)/dble(N*SAmvAdjstRt)
        if(rate .lt. accRtSAmv) then
            rD = rD*0.95d0
        else
            rD = rD*1.05d0
        endif
        if(rD .gt. rL/4.0d0) then
            rD = rD/1.05d0
        endif

        accSAmv = 0
        write(*,'("EQ STEP: ",I10," rD: ",1f20.10)') i, rD
    endif
enddo

open(unit=confU, file=trim(argUuid)//conf, form='formatted',&
        status='new', action='write')
do i=1, prS
    ! Perform MC Step - attempt to Move all @N atoms
    call random_number(rndN)
    do k=1, N
        A = min(1 + floor(N*rndN(k)),N) 
        ! floor(N*rnd) maps to 0, N-1 -> in exceptionally rare cases to 0-N
        call random_number(rnd4)
        call SAmv(N, A, r, rD, rL, rB, rCut6, e6, e12, rnd4, accSAmv)
    enddo
    ! Take PROD values
    if(mod(i,smplRt) .eq. 0) then
        call CPU_TIME(t1)
        eTot = e12+e6 + (eLRC12+eLRC6)
        vir  = virial(N, r, rL, rCut)
        pres = (dble(N)/rB+vir+vLRC12+vLRC6)/rV
        write(*,'("PR STEP: ",I10," E: ",1f10.5," P: ",1f10.5,"&
            & vir: ",1f10.5," time: ",1f10.5)') i, eTot, pres, vir, t1-t0
        call append_r(N, r, confU, rL, eTot, pres)
        call CPU_TIME(t0)
    endif    
enddo
close(confU)

end program ljNVT

subroutine SApot(N, A, r, rL, rCut6, e6, e12)
    ! Compute potential of single atom within the gas
    implicit none

    integer, intent(in) :: N, A
    real*8, dimension(3,N), intent(in) :: r
    real*8, intent(in) :: rL, rCut6
    real*8, intent(out) :: e6, e12

    integer :: i
    real*8, dimension(3) :: dr
    real*8 :: dr6

    e6  = 0.0d0
    e12 = 0.0d0
    do i=1, N
        if(i .ne. A) then
            dr  = r(:,i) - r(:,A)
            dr  = dr - rL*anint(dr/rL)
            dr6 = DOT_PRODUCT(dr,dr)**3.0d0
            ! cutoff
            if(dr6 .lt. rCut6) then
                e6  = e6  - 1.0d0/dr6
                e12 = e12 + 1.0d0/(dr6*dr6)
            endif
        endif
    enddo
    e6  = 4.0d0*e6
    e12 = 4.0d0*e12
end subroutine SApot

subroutine ljTot(N, r, rL, rCut6, e6, e12)
    ! Compute potential of single atom within the gas
    implicit none

    integer, intent(in) :: N
    real*8, dimension(3,N), intent(in) :: r
    real*8, intent(in) :: rL, rCut6
    real*8, intent(out) :: e6, e12

    integer :: i, j
    real*8, dimension(3) :: dr
    real*8 :: dr6
    
    e6  = 0.0d0
    e12 = 0.0d0
    do i=1, N-1
        do j=i+1, N
            dr  = r(:,j) - r(:,i)
            dr  = dr - rL*anint(dr/rL)
            dr6 = DOT_PRODUCT(dr,dr)**3.0d0
            ! cutoff - enforces single image
            if(dr6 .lt. rCut6) then
                e6  = e6  - 1.0d0/dr6
                e12 = e12 + 1.0d0/(dr6*dr6)
            endif
        enddo
    enddo
    e6  = 4.0d0*e6
    e12 = 4.0d0*e12
end subroutine ljTot

subroutine SAmv(N, A, r, dlt, rL, rB, rCut6, e6, e12, rnd4, acc)
    ! Perform single atom move on atom @A with
    ! magnitude of displacement given by @dlt  
    implicit none

    integer, intent(in) :: N, A
    real*8, dimension(3,N), intent(inout) :: r
    real, dimension(4), intent(in) :: rnd4
    real*8, intent(in) :: rL, rB, dlt, rCut6
    real*8, intent(inout) :: e6, e12
    integer, intent(inout) :: acc

    real*8, dimension(3) :: old_R
    real*8 :: o_e12, o_e6, n_e12, n_e6, dE6, dE12, dE

    call SApot(N, A, r, rL, rCut6, o_e6, o_e12)
    old_R = r(:,A)

    ! Displace atom @A
    r(1,A) = r(1,A) + dlt*(rnd4(1)-0.5d0)
    r(2,A) = r(2,A) + dlt*(rnd4(2)-0.5d0)
    r(3,A) = r(3,A) + dlt*(rnd4(3)-0.5d0)
    ! Apply PBC wrt @rL
    r(1,A) = r(1,A) - rL*anint(r(1,A)/rL)
    r(2,A) = r(2,A) - rL*anint(r(2,A)/rL)
    r(3,A) = r(3,A) - rL*anint(r(3,A)/rL)

    call SApot(N, A, r, rL, rCut6, n_e6, n_e12)
    dE6  = n_e6  - o_e6
    dE12 = n_e12 - o_e12

    if(exp(-rB*(dE12+dE6)) .gt. rnd4(4)) then   
        e12 = e12 + dE12
        e6  = e6  + dE6
        acc = acc + 1
    else
        r(:,A) = old_R
    endif
end subroutine SAmv

subroutine setup_sc(N, r, rL)
    ! Initialize postions @r of @N atoms into SC
    ! lattice and apply PBC wrt @rL
    implicit none

    integer, intent(in) :: N
    real*8, intent(in) :: rL
    real*8, dimension(3,N), intent(inout) :: r

    real*8 :: del, a_3
    real*8, allocatable, dimension(:,:) :: sc
    integer :: a, id, i, j, k

    write(*,'("CREATING SC LATTICE")')
    a_3 = dble(N)**(1.0d0/3.0d0)
    write(*,'("N^(1/3)",13X,1f10.5)') a_3
    if(nint(a_3) - aint(a_3) .gt. 0.5d0) then
        a = nint(a_3)
    else
        a = nint(a_3)+1
    endif
    write(*,'("a of SC",13X,I10)') a
    del = rL/dble(a)
    allocate(sc(3,a*a*a))
    do i=1, a
        do j=1, a
            do k=1, a
                id = k + (j-1)*a + (i-1)*a*a
                sc(1,id) = -rL/2.0d0 + del*k
                sc(2,id) = -rL/2.0d0 + del*j
                sc(3,id) = -rL/2.0d0 + del*i
            enddo
        enddo
    enddo
    ! arrange atoms in sc and apply PBC
    do i=1, N
        r(:,i) = sc(:,i)
        r(:,i) = r(:,i) - rL*anint(r(:,i)/rL)
    enddo
end subroutine setup_sc

subroutine print_InitR(N, r, name)
    ! Print formatted positions @r of all @N atoms to file @name
    ! in .xyz format 
    implicit none

    integer, intent(in) :: N
    character(len=20), intent(in) :: name
    real*8, dimension(3,N), intent(in) :: r

    integer, parameter :: outU = 20
    integer :: i

    open(unit=outU, file=trim(name), form='formatted',&
    status='new', action='write')
    write(outU, *) N
    write(outU, '("initial sc arrangement")')
    do i=1, N
        write(outU, '("Ar ",1F20.10," ",1F20.10," ",1F20.10)') r(:,i) 
    enddo
    close(outU)
end subroutine print_InitR

subroutine append_r(N, r, unt, rL, E, pres)
    ! Print formatted positions @r of all @N atoms to file @name
    ! in .xyz format 
    implicit none

    integer, intent(in) :: N, unt
    real*8, intent(in) :: rL, E, pres
    real*8, dimension(3,N), intent(in) :: r

    integer :: i

    write(unt, *) N
    write(unt, '("rL= ",1f10.5," E= ",1f10.5," P= ",1f10.5)') rL, E, pres
    do i=1, N
        write(unt, '("Ar ",1F20.10," ",1F20.10," ",1F20.10)') r(:,i) 
    enddo
end subroutine append_r

function virial(N, r, rL, rCut) result(vir)
    ! Computes virial 
    implicit none

    integer, intent(in) :: N
    real*8, intent(in) :: rL, rCut
    real*8, dimension(3,N), intent(in) :: r

    real*8, dimension(3) :: rij, fij
    integer :: i, j
    real*8 :: vir, f, dr

    vir = 0.0d0
    do i=1, N-1
        do j=i+1, N
            rij = r(:,j)-r(:,i)
            rij = rij - rL*anint(rij/rL)
            dr  = sqrt(DOT_PRODUCT(rij,rij))
            if(dr .lt. rCut) then
                f = 4.0d0*(12.0d0/(dr**13.00)-6.0d0/(dr**7.0d0))
                vir = vir + f*dr
            endif
        enddo
    end do
    vir = vir/3.0d0
end function virial
