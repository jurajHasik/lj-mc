program ljNPT
! NPT simulation of LJ Gas with cutoff rescaling

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
real*8 :: eps = 1.0d0 ! [eV]
integer :: N, A
real*8 :: temp, rT, rB, rD, rDens, pres, rP
real*8 :: rL, rV, rCut, rCut6
real*8 :: rDV

integer :: accSAmv, accSclMv ! Acceptance of single atom moves
real*8, parameter :: accRtSAmv = 0.35d0 ! Desired acceptance rate
real*8, parameter :: accRtSclMv = 0.35d0 ! Desired acceptance rate
! Frequency of updating SAmv delta in [MC steps]
integer, parameter :: SAmvAdjstRt = 100
integer, parameter :: sclMvAdjstRt = 100
integer :: eqS, prS
integer :: smplRt
integer, parameter :: sclMvRt = 1

! ***** Atom coords *****
real*8, allocatable, dimension(:,:) :: r
real*8 :: e6, e12 ! holds contribution from r^-6 & r^-12 part
! Long range corrections to energy and virial
real*8 :: eLRC6, eLRC12, vLRC6, vLRC12
real*8 :: eTot, vir, vir6, vir12
real*8 :: rate

! ***** PRNG variables *****
integer :: seed
real, dimension(2) :: rnd2
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
character(len=20) :: argUuid, argN, argT, argrP, argEqS, argPrS
character(len=20) :: argSmplRt, argSeed

! ***** Utility variables *****
integer :: i, k
real*8 :: t0, t1

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
call GETARG(4, argrP) ! read in pressure in LJUnits
call GETARG(5, argEqS)
call GETARG(6, argPrS)
call GETARG(7, argSmplRt)
call GETARG(8, argSeed)
read(argN, *) N
read(argT, *) temp
read(argrP, *) rP
read(argEqS, *) eqS
read(argPrS, *) prS
read(argSmplRt, *) smplRt
read(argSeed, *) seed

allocate(r(3,N), rndN(N))

rT = temp/eps
rB = 1.0d0/rT

rV = N*rT/rP ! Assume init volume as the one for ideal gas
rDens = N/rV 
rL = rV**(1.0d0/3.0d0)
! cutoff fixed to L/2
rCut = rL/2.0d0
! we actually only care about a rCut^2
rCut6 = rCut**6.0d0
rD = rL/10.0d0 ! Set initial delta to rL/10 
rDV = rV/1000.0d0 ! Set initial deltaVol to rV/1000


write(*,'("N",19X,I10)') N
write(*,'("REDUCED UNITS")')
write(*,'("rP",18X,1f10.5)') rP
write(*,'("rL",18X,1f10.5)') rL
write(*,'("rCut",16X,1f10.5)') rCut
write(*,'("rD",18X,1f10.5)') rD
write(*,'("rT",18X,1f10.5)') rT 
write(*,'("rBeta",15X,1f10.5)') rB
write(*,'("rDens",15X,1f10.5)') rDens

! ***** Setup initial sc arrangement *****
call setup_sc(N, r, rL)
! Optional check of init conf of sc-lattice
!call print_InitR(N, r, initSc)
write(*,'("INITIAL SETUP INTO SC LATTICE DONE")')
! compute initial energy
call ljTot(N, r, rL, rCut6, e6, e12)
! compute energy LRC
eLRC6  = -8.0d0*Pi*N*rDens*(1.0d0/(3.0d0*(rCut**3.0d0)))
eLRC12 = 8.0d0*Pi*N*rDens*(1.0d0/(9.0d0*(rCut**9.0d0)))
vLRC6  = 2.0d0*eLRC6
vLRC12 = 4.0d0*eLRC12
eTot = e12+e6 + (eLRC12+eLRC6) ! No double counting coming from ljTot
write(*,'("init E",14X,1f10.5)') eTot
write(*,'("E_lrc: ",1f10.5," V_rlc: ",1f10.5)') eLRC12+eLRC6, vLRC12+vLRC6

! Initialize MT PRNG
call random_setseed(seed)
write(*,'("Seed",16X,I10)') seed
! Check the output of PRNG
write(*,'("-sample output-")')
call random_number(rnd10)
do i=1, 10
    write(*,*) rnd10(i)
enddo

accSAmv = 0
accSclMv = 0
call CPU_TIME(t0)
do i=1, eqS
    ! Perform MC Step - attempt to Move all @N atoms
    call random_number(rndN)
    do k=1, N
        call random_number(rnd4)
        !Always select random atom instead of keeping fixed order
!call SAmv(N, k, r, rD, rL, rB, lj(1,k), lj(2,k), rnd4, accSAmv)
        A = min(1 + floor(N*rndN(k)),N)
        ! floor(N*rnd) maps to 0, N - since rnd is elem of [0,1]
        call SAmv(N, A, r, rD, rL, rB, rCut6, e6, e12, rnd4, accSAmv)
    enddo
    if(mod(i,sclMvRt) .eq. 0) then
        call random_number(rnd2)
        call sclMv(N, rP, rDV, rnd2, r, e6, e12, vir6, vir12, eLRC6, &
            &eLRC12, vLRC6, vLRC12, rB, rD, rL, rV, rCut, rCut6, rDens,&
            &accSclMv)
    endif
    ! Take EQ values
    if(mod(i,smplRt) .eq. 0) then
        call CPU_TIME(t1)
        eTot = e12+e6 + (eLRC12+eLRC6)
        call virial(vir6, vir12, N, r, rL, rCut)
        !Eq of state + correction
        vir = vir6+vir12+vLRC12+vLRC6
        !Only contribution from virial since p_ext is fixed
        pres = vir/rV
        write(*,'("EQ STEP: ",I10," E: ",1f10.5," P: ",1f10.5,"&
            & rDens: ",1f10.5," time: ",1f10.5)') i, eTot, pres, rDens,&
            & t1-t0
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
        write(*,'("EQ STEP: ",I10," rD [sigma]: ",1f20.10)') i, rD
    endif
    if(mod(i,sclMvAdjstRt) .eq. 0) then
        rate = dble(accSclMv)/dble(sclMvAdjstRt/sclMvRt)
        write(*,'("sclMv rate: ",1f10.5)') rate
        if(rate .lt. accRtSclMv) then
            rDV = rDV*0.95d0
        else
            rDV = rDV*1.05d0
        endif
        accSclMv = 0
        write(*,'("EQ STEP: ",I10," rDV: ",1f20.10)') i, rDV
    endif
enddo

open(unit=confU, file=trim(argUuid)//conf, form='formatted',&
        status='new', action='write')
do i=1, prS
    ! Perform MC Step - attempt to Move all @N atoms
    call random_number(rndN)
    do k=1, N
        A = min(1 + floor(N*rndN(k)),N)
        call random_number(rnd4)
        call SAmv(N, A, r, rD, rL, rB, rCut6, e6, e12, rnd4, accSAmv)
    enddo
    if(mod(i,sclMvRt) .eq. 0) then
        call random_number(rnd2)
        call sclMv(N, rP, rDV, rnd2, r, e6, e12, vir6, vir12, eLRC6, &
            &eLRC12, vLRC6, vLRC12, rB, rD, rL, rV, rCut, rCut6, rDens,&
            &accSclMv)
    endif
    ! Take PROD values
    if(mod(i,smplRt) .eq. 0) then
        call CPU_TIME(t1)
        eTot = e12+e6 + (eLRC12+eLRC6)
        call virial(vir6, vir12, N, r, rL, rCut)
        !Eq of state + correction
        vir = vir6+vir12+vLRC12+vLRC6
        !Only contribution from virial since p_ext is fixed
        pres = vir/rV
        write(*,'("PR STEP: ",I10," E: ",1f10.5," P: ",1f10.5,"&
            & rDens: ",1f10.5," vir: ",1f10.5," time: ",1f10.5)') i, eTot,&
            & pres, rDens, vir, t1-t0
        call append_r(N, r, confU, rL, eTot, pres)
        accSAmv  = 0
        accSclMv = 0
        call CPU_TIME(t0)
    endif    
enddo
close(confU)

end program ljNPT

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
    real*8, intent(inout) :: e12, e6
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

subroutine sclMv(N, ext_p, rDV, rnd2, r, e6, e12, v6, v12, eLRC6, &
    &eLRC12, vLRC6, vLRC12, rB, rD, rL, rV, rCut, rCut6, rDens, acc)
    implicit none

    integer, intent(in) :: N
    real, dimension(2), intent(in) :: rnd2
    real*8, dimension(3,N), intent(inout) :: r
    real*8, intent(in)    :: rDV, ext_p, rB
    real*8, intent(inout) :: rD, rL, rV, rCut, rCut6, rDens
    real*8, intent(inout) :: e6, e12, v6, v12
    real*8, intent(inout) :: eLRC6, eLRC12, vLRC6, vLRC12
    integer, intent(inout) :: acc

    real*8 :: scl, iscl3, iscl6, delH
    real*8 :: n_rV, n_eLRC6, n_eLRC12

    n_rV = rV + rDV*(rnd2(1)-0.5d0)
    scl = (n_rV/rV)**(1.0d0/3.0d0)

    iscl6 = scl**(-6.0d0)
    iscl3 = scl**(-3.0d0)
    ! the analytic expression for new eLRC
    ! eLRC6 = -8.0d0*Pi*N*rDens*(1.0d0/(3.0d0*(rCut**3.0d0)))
    ! eLRC12 = 8.0d0*Pi*N*rDens*(1.0d0/(9.0d0*(rCut**9.0d0)))
    ! iscl^3 for rescaled rDens, iscl^3 for rescaled cutoff
    n_eLRC6  = eLRC6*iscl6
    ! iscl^3 for rescaled rDens, iscl^9 for rescaled cutoff
    n_eLRC12 = eLRC12*iscl6*iscl6

    ! NPT ensamble is weighted as exp(-beta*(PV+U)+NlogV)
    delH = (e12*iscl6 + e6)*iscl6 - (e12+e6) + &
        &(n_eLRC12+n_eLRC6) - (eLRC12+eLRC6)+&
        &ext_p*rV*(scl*scl*scl-1.0d0) - real(3*N)*log(scl)/rB
    if(exp(-rB*delH) .gt. rnd2(2)) then
        r=r*scl
        rD=rD*scl
        rL=rL*scl
        rV=n_rV
        rDens=N/n_rV
        rCut=rCut*scl
        rCut6=rCut6*(scl**6.0d0)
        e6=e6*iscl6
        v6=v6*iscl6
        e12=e12*iscl6*iscl6
        v12=v12*iscl6*iscl6
        eLRC6  = n_eLRC6
        eLRC12 = n_eLRC12
        vLRC6  = 2.0d0*n_eLRC6
        vLRC12 = 4.0d0*n_eLRC12
        acc=acc+1
    endif

end subroutine sclMv

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

subroutine virial(v6, v12, N, r, rL, rCut)
    ! Computes virial 
    implicit none

    integer, intent(in) :: N
    real*8, intent(in) :: rL, rCut
    real*8, dimension(3,N), intent(in) :: r
    real*8, intent(out) :: v6, v12

    real*8, dimension(3) :: rij, fij
    integer :: i, j
    real*8 :: f, dr

    v6  = 0.0d0
    v12 = 0.0d0
    do i=1, N-1
        do j=i+1, N
            rij = r(:,j)-r(:,i)
            rij = rij - rL*anint(rij/rL)
            dr  = sqrt(DOT_PRODUCT(rij,rij))
            if(dr .lt. rCut) then
                v6  = v6  - 4.0d0*6.0d0/(dr**7.0d0)*dr/3.0d0
                v12 = v12 + 4.0d0*12.0d0/(dr**13.00)*dr/3.0d0
            endif
        enddo
    end do
end subroutine virial 
