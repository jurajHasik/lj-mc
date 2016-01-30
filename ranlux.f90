MODULE Ran_Lux_Mod
! Subtract-and-borrow random number generator proposed by Marsaglia and Zaman, implemented by F. James with the name RCARRY in 1991,
! and later improved by Martin Luescher in 1993 to produce "Luxury Pseudorandom Numbers". Fortran 77 coded by F. James, 1993.
! Converted to Fortran 90 [and Lahey Elf90 subset] by Loren Meissner, 1995.
!
! References: M. Luscher, Computer Physics Communications 79 (1994) 100; F. James, Computer Physics Communications 79 (1994) 111
!
! LUXURY LEVELS. -- The available luxury levels are:
! level 0 (p = 24) : equivalent to the original RCARRY of Marsaglia and Zaman, very long period, but fails many tests.
! level 1 (p = 48) : considerable improvement in quality over level 0, now passes the gap test, but still fails spectral test.
! level 2 (p = 97) : passes all known tests, but theoretically still defective.
! level 3 (p = 223) : DEFAULT VALUE. Any theoretically possible correlations have very small chance of being observed.
! level 4 (p = 389) : highest possible luxury, all 24 bits chaotic.
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
! Calling sequences for RanLux:
! call RanLux (RVec)
!   Returns a vector RVec of Len(RVec) 32-bit random floating point numbers X, such that 0.0 < X < 1.0 .
! call RLuxGo (Lux, Int, K1, K2)
!   Initializes the generator from one 32-bit integer INT and sets Luxury Level LUX which is integer between zero and MaxLev, or
!   if Lux > 24, it sets p = Lux directly. K1 and K2 should be set to zero unless restarting at a break point given by output of
!   RLuxAt (see RLuxAt).
! call RLuxAt (Lux, Int, K1, K2)
!   Gets the values of four integers which can be used to restart the RanLux generator at the current point by calling RLuxGo.
!   K1 and K2 specify how many numbers were generated since the initialization with Lux and Int. The restarting skips over
!   K1 + K2 * 1E9 numbers, so it can be long. A more efficient but less convenient way of restarting is by:
!     call RLuxIn (ISVec) ! Restart the generator from vector ISVec of 25 32-bit integers (see RLXUt)
!     call RLuxUt (ISVec) ! Output the current values of the 25 32-bit integer Seeds, to be used for restarting.
! The array argument to RLuxIn or RLuxUt must be dimensioned 25 in the calling program
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
!                                     default
!       Luxury Level    0     1     2   * 3 *   4
!                       0    24    73   199   365
!   corresponds to p = 24    48    97   223   389
!           time factor 1     2     3     6    10 on slow workstation
!                       1   1.5     2     3     5 on fast mainframe
!
! NotYet is .TRUE. if no initialization has been performed yet.

  IMPLICIT NONE

    private
    public RanLux, RLuxIn, RLuxUt, RLuxAt, RLuxGo 

  INTEGER, PARAMETER :: NSeeds = 25, MaxLev = 4, LxDflt = 3
  REAL, PARAMETER :: TwoP12 = 4096.0
  INTEGER, PARAMETER :: IGiga = 1000000000, JSDFlt = 314159265, ITwo24 = 2 ** 24, ICons = 2147483563
  INTEGER :: I
  INTEGER, PARAMETER :: Next(NSeeds - 1) = (/ NSeeds - 1, (I, I = 1, NSeeds - 2) /)  ! Table look-up (faster than Mod function).

  INTEGER :: I24 = 24, J24 = 10, In24 = 0, Kount = 0, LuxLev = LxDflt, MKount = 0 ! Initialized variables are automatically saved.
  INTEGER, DIMENSION(0: MaxLev) :: NDSkip = (/ 0, 24, 73, 199, 365 /) ! Initialized variables are automatically saved.
  INTEGER, SAVE :: NSkip, InSeed
  REAL :: Carry = 0.0 ! Initialized variables are automatically saved.
  REAL, SAVE :: Seeds(NSeeds - 1), TwoM24, TwoM12
  LOGICAL, SAVE :: NotYet = .TRUE.

  REAL :: Uni

  PRIVATE :: RCarry

CONTAINS

  SUBROUTINE RanLux (RVec)
    ! Default Initialization by Multiplicative Congruential
    REAL, INTENT(out) :: RVec(:)
    INTEGER :: ISeeds(NSeeds - 1), I, IVec, JSeed, K, LEnv, LP

! start subroutine RanLux
    LEnv = SIZE (RVec)
    IF (NotYet) THEN
      NotYet = .FALSE.
      JSeed = JSDFlt
      InSeed = JSeed
      PRINT *, " RanLux default initialization: ", JSeed
      LuxLev = LxDflt
      NSkip = NDSkip(LuxLev)
      LP = NSkip + NSeeds - 1
      In24 = 0
      Kount = 0
      MKount = 0
      PRINT *, " RanLux default luxury level = ", LuxLev, " p = ", LP
      TwoM24 = 1.0
      DO I = 1, NSeeds - 1
        TwoM24 = TwoM24 * 0.5
        K = JSeed / 53668
        JSeed = 40014 * (JSeed - K * 53668) - K * 12211
        IF (JSeed < 0) JSeed = JSeed + ICons
        ISeeds(I) = MOD (JSeed, ITwo24)
      END DO
      TwoM12 = TwoM24 * 4096.0
      Seeds = REAL (ISeeds) * TwoM24
      I24 = NSeeds - 1
      J24 = 10
      Carry = MERGE (TwoM24, 0.0, Seeds(NSeeds - 1) == 0.0)
    END IF

    DO IVec = 1, LEnv
      RVec(IVec) = RCarry (1)
    ! Skipping to Luxury. As proposed by Martin Luscher.
      In24 = In24 + 1
      IF (In24 == NSeeds - 1) THEN
        In24 = 0
        Kount = Kount + NSkip
        Uni = RCarry (NSkip)
      END IF
   END DO
   ! "Pad" small numbers (with less than 12 "significant" bits) and eliminate zero values (in case someone takes a logarithm)
    WHERE (RVec < TwoM12) RVec = RVec + TwoM24 * Seeds(J24)
    WHERE (Rvec == 0.0) RVec = TwoM24 * TwoM24
    Kount = Kount + LEnv
    IF (Kount >= IGiga) THEN
      MKount = MKount + 1
      Kount = Kount - IGiga
    END IF
    RETURN
  END SUBROUTINE RanLux

! Input and float integer Seeds from previous run
  SUBROUTINE RLuxIn (ISDext)
    INTEGER, INTENT(in) :: ISDext(:)
    INTEGER :: I, ISD
! start subroutine RLuxIn
    IF (SIZE(ISDext) /= NSeeds) THEN
      PRINT *, " Array size for RLuxIn must be ", NSeeds
      RETURN
    END IF
    ! The following IF block added by Phillip Helbig, based on conversation with Fred James;
    ! an equivalent correction has been published by James.
    IF (NotYet) THEN
      PRINT *, " Proper results only with initialisation from 25 integers obtained with RLuxUt"
      NotYet = .FALSE.
    END IF
    TwoM24 = 1.0
    DO I = 1, NSeeds - 1
      TwoM24 = TwoM24 * 0.5
    END DO
    TwoM12 = TwoM24 * 4096.0
    PRINT *, " Full initialization of RanLux with 25 integers:"
    PRINT *, ISDext
    Seeds = REAL (ISDext(: NSeeds - 1)) * TwoM24
    Carry = 0.0
    IF (ISDext(NSeeds) < 0) Carry = TwoM24
    ISD = ABS (ISDext(NSeeds))
    I24 = MOD (ISD, 100)
    ISD = ISD / 100
    J24 = MOD (ISD, 100)
    ISD = ISD / 100
    In24 = MOD (ISD, 100)
    ISD = ISD / 100
    LuxLev = ISD
    IF (LuxLev <= MaxLev) THEN
      NSkip = NDSkip(LuxLev)
      PRINT *, " RanLux luxury level set by RLuxIn to: ", LuxLev
    ELSE IF (LuxLev >= NSeeds - 1) THEN
      NSkip = LuxLev - NSeeds + 1
      PRINT *, " RanLux p-value set by RLuxIn to:", LuxLev
    ELSE
      NSkip = NDSkip(MaxLev)
      PRINT *, " RanLux illegal luxury RLuxIn: ", LuxLev
      LuxLev = MaxLev
    END IF
    InSeed = - 1
    RETURN
  END SUBROUTINE RLuxIn

! Ouput Seeds as integers
  SUBROUTINE RLuxUt (ISDext)
    INTEGER, INTENT(out) :: ISDext(:)
! start subroutine RLuxUt
    IF (SIZE(ISDext) /= NSeeds) THEN
      ISDext = 0
      PRINT *, " Array size for RLuxUt must be ", NSeeds
      RETURN
    END IF
    ISDext(: NSeeds - 1) = INT (Seeds * TwoP12 * TwoP12)
    ISDext(NSeeds) = MERGE (-ISDext(NSeeds), I24 + 100 * J24 + 10000 * In24 + 1000000 * LuxLev, Carry > 0.0)
    RETURN
  END SUBROUTINE RLuxUt

! Output the "convenient" restart point
  SUBROUTINE RLuxAt (LOut, InOut, K1, K2)
    INTEGER, INTENT(out) :: LOut, InOut, K1, K2
! start subroutine RLuxAt
    LOut = LuxLev
    InOut = InSeed
    K1 = Kount
    K2 = MKount
    RETURN
  END SUBROUTINE RLuxAt

! Initialize from one or three integers
  SUBROUTINE RLuxGo (Lux, Int, K1, K2)
    INTEGER, INTENT(in) :: Lux, Int, K1, K2
    INTEGER :: ISeeds(NSeeds - 1), ILx, I, IOuter, IZip, IZip2, JSeed, K
! start subroutine RLuxGo
    IF (Lux < 0) THEN
      LuxLev = LxDflt
    ELSE IF (Lux <= MaxLev) THEN
      LuxLev = Lux
    ELSE IF (Lux < NSeeds - 1 .OR. Lux > 2000) THEN
      LuxLev = MaxLev
      PRINT *, " RanLux illegal luxury level in RLuxGo: ", Lux
    ELSE
      LuxLev = Lux
      DO ILx = 0, MaxLev
        IF (Lux == NDSkip(ILx) + NSeeds - 1) THEN
          LuxLev = ILx
        END IF
      END DO
    END IF
    IF (LuxLev <= MaxLev) THEN
      NSkip = NDSkip(LuxLev)
      PRINT *, " RanLux luxury level set by RLuxGo :", LuxLev, " p = ", NSkip + NSeeds - 1
    ELSE
      NSkip = LuxLev - 24
      PRINT *, " RanLux p-value set by RLuxGo to:", LuxLev
    END IF
    In24 = 0
    IF (Int < 0) THEN
      PRINT *, " Illegal initialization by RLuxGo, negative input seed"
    ELSE IF (Int > 0) THEN
      JSeed = Int
      PRINT *, " RanLux initialized by RLuxGo from Seeds", JSeed, K1, K2
    ELSE
      JSeed = JSDFlt
      PRINT *, " RanLux initialized by RLuxGo from default seed"
    END IF
    InSeed = JSeed
    NotYet = .FALSE.
    TwoM24 = 1.0
    DO I = 1, NSeeds - 1
      TwoM24 = TwoM24 * 0.5
      K = JSeed / 53668
      JSeed = 40014 * (JSeed - K * 53668) - K * 12211
      IF (JSeed < 0) JSeed = JSeed + ICons
      ISeeds(I) = MOD (JSeed, ITwo24)
    END DO
    TwoM12 = TwoM24 * 4096.0
    Seeds = REAL (ISeeds) * TwoM24
    I24 = NSeeds - 1
    J24 = 10
    Carry = MERGE (TwoM24, 0.0, Seeds(NSeeds - 1) == 0.0)

    ! If restarting at a break point, skip K1 + IGIGA * K2
    ! Note that this is the number of numbers delivered to the user PLUS the number skipped (if Luxury > 0) .
    Kount = ABS (K1)
    MKount = ABS (K2)
    IF (Kount + MKount /= 0) THEN
      DO IOuter = 1, MKount + 1
        Uni = RCarry (MERGE (Kount, IGiga, IOuter == MKount + 1))
      END DO
      ! Get the right value of IN24 by direct calculation
      In24 = MOD (Kount, NSkip + NSeeds - 1)
      IF (MKount > 0) THEN
        IZip = MOD (IGiga, NSkip + NSeeds - 1)
        IZip2 = MKount * IZip + In24
        In24 = MOD (IZip2, NSkip + NSeeds - 1)
      END IF
      ! Now IN24 had better be between zero and 23 inclusive
      IF ((In24 < 1) .OR. (In24 >= NSeeds - 1)) THEN
        PRINT *, " Error in restarting with RLuxGo: the values", Int, K1, K2, " cannot occur at luxury level", LuxLev
        In24 = 0
      END IF
    END IF
    RETURN
  END SUBROUTINE RLuxGo

  FUNCTION RCarry (N) RESULT (Uni)  ! Private (in module); generates a sequence of N uniform random numbers; returns the last one.
    REAL :: Uni
    INTEGER, INTENT(in) :: N
    INTEGER :: Many
! start function RCarry
    DO Many = 1, N
    ! The Generator proper: "Subtract-with-borrow", as proposed by Marsaglia and Zaman, Florida State University, March, 1989
      Uni = Seeds(J24) - Seeds(I24) - Carry
      IF (Uni < 0.0) THEN
        Uni = Uni + 1.0
        Carry = TwoM24
      ELSE
        Carry = 0.0
      END IF
      Seeds(I24) = Uni
      I24 = Next(I24)
      J24 = Next(J24)
    END DO
    RETURN
  END FUNCTION RCarry

END MODULE Ran_Lux_Mod

!!$! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
!!$  PROGRAM Lux_Tst
!!$  ! Exercise for the RANLUX Pseudorandom number generator.
!!$
!!$  USE Ran_Lux_Mod
!!$  IMPLICIT NONE
!!$
!!$  REAL, ALLOCATABLE :: RVec(:)
!!$  INTEGER :: ISDext(NSeeds), I1, I2, I3, I4, Li
!!$
!!$  ALLOCATE (RVec(100))
!!$  PRINT *, " Vector length is now 100 "
!!$
!!$! Check that we get the right numbers (machine-indep.)
!!$  PRINT *, " call RanLux(RVec) "
!!$  CALL RanLux (RVec)
!!$  PRINT *, " RanLux default numbers 1 - 5:"
!!$  PRINT *, RVec(: 5)
!!$  CALL RanLux (RVec)
!!$  PRINT *, " RanLux default numbers 101 - 105:"
!!$  PRINT *, RVec(: 5)
!!$
!!$  PRINT *, " call RLuxGo (0, 0, 0, 0) "
!!$  CALL RLuxGo (0, 0, 0, 0)
!!$  CALL RanLux (RVec)
!!$  PRINT *, " RanLux luxury level 0, 1 - 5:"
!!$  PRINT *, RVec(: 5)
!!$  CALL RanLux (RVec)
!!$  PRINT *, " RanLux luxury level 0, 101 - 105:"
!!$  PRINT *, RVec(: 5)
!!$
!!$  PRINT *, " call RLuxGo (389, 1, 0, 0) "
!!$  CALL RLuxGo (389, 1, 0, 0)
!!$  CALL RanLux (RVec)
!!$  PRINT *, " RanLux luxury p = 389, 1 - 5:"
!!$  PRINT *, RVec(: 5)
!!$  CALL RanLux (RVec)
!!$  PRINT *, " RanLux luxury p = 389, 101 - 105:"
!!$  PRINT *, RVec(: 5)
!!$
!!$  PRINT *, " call RLuxGo (75, 0, 0, 0) "
!!$  CALL RLuxGo (75, 0, 0, 0)
!!$  CALL RanLux (RVec)
!!$  PRINT *, " RanLux luxury p = 75, 1 - 5:"
!!$  PRINT *, RVec(: 5)
!!$  CALL RanLux (RVec)
!!$  PRINT *, " RanLux luxury p = 75, 101 - 105:"
!!$  PRINT *, RVec(: 5)
!!$!  if (RVec(1) > 0.0) stop "OK!"
!!$
!!$  PRINT *, " Test restarting from the full vector"
!!$  CALL RLuxUt (ISDext)
!!$  PRINT *, " Current RanLux status saved:"
!!$  PRINT *, ISDext
!!$  CALL RanLux (RVec)
!!$  PRINT *, " RanLux numbers 1-5:"
!!$  PRINT *, RVec(: 5)
!!$  CALL RanLux (RVec)
!!$  PRINT *, " RanLux numbers 101 - 105:"
!!$  PRINT *, RVec(: 5)
!!$
!!$  PRINT *, " Previous RanLux status will be restored"
!!$  CALL RLuxIn (ISDext)
!!$  CALL RanLux (RVec)
!!$  PRINT *, " RanLux numbers 1-5:"
!!$  PRINT *, RVec(: 5)
!!$  CALL RanLux (RVec)
!!$  PRINT *, " RanLux numbers 101 - 105:"
!!$  PRINT *, RVec(: 5)
!!$
!!$  PRINT *, " Test the restarting by skipping"
!!$  CALL RLuxGo (4, 7674985, 0, 0)
!!$  CALL RLuxAt (i1, i2, i3, i4)
!!$  PRINT *, " RLuxAt values = ", i1, i2, i3, i4
!!$
!!$  DEALLOCATE (RVec)
!!$  ALLOCATE (RVec(1000))
!!$  PRINT *, " Vector length is now 1000 "
!!$
!!$  DO Li = 1, 10
!!$    CALL RanLux (RVec)
!!$  END DO
!!$  CALL RLuxAt (i1, i2, i3, i4)
!!$  PRINT *, " RLuxAt values = ", i1, i2, i3, i4
!!$
!!$  DEALLOCATE (RVec)
!!$  ALLOCATE (RVec(200))
!!$  PRINT *, " Vector length is now 200 "
!!$
!!$  CALL RanLux (RVec)
!!$  PRINT *, " Next and 200th numbers are:", RVec(1), RVec(200)
!!$  CALL RLuxGo (i1, i2, i3, i4)
!!$  CALL RanLux (RVec)
!!$  PRINT *, " Next and 200th numbers are:", RVec(1), RVec(200)
!!$
!!$  PRINT *, " The following should provoke an error message"
!!$  CALL RLuxGo (4, 11111, 31, 0)
!!$  STOP
!!$
!!$END PROGRAM Lux_Tst
!!$
!!$! Results from NAG F90 on IBM RS-6000 [L. Meissner, 22 Aug 1995]
!!$!   Vector length is now 100
!!$!   call RanLux(RVec)
!!$!   RanLux default initialization:  314159265
!!$!   RanLux default luxury level =  3  p =  223
!!$!   RanLux default numbers 1 - 5:
!!$!    0.5398182   0.7615504   6.0299397E-02   0.7960026   0.3063122
!!$!   RanLux default numbers 101 - 105:
!!$!    0.4315674   3.7744164E-02   0.2489711   1.4778376E-03   0.9027445
!!$!   call RLuxGo (0, 0, 0, 0)
!!$!   RanLux luxury level set by RLuxGo : 0  p =  24
!!$!   RanLux initialized by RLuxGo from default seed
!!$!   RanLux luxury level 0, 1 - 5:
!!$!    0.5398182   0.7615504   6.0299397E-02   0.7960026   0.3063122
!!$!   RanLux luxury level 0, 101 - 105:
!!$!    0.4153877   5.3309321E-02   0.5819531   0.9139745   0.6703444
!!$!   call RLuxGo (389, 1, 0, 0)
!!$!   RanLux luxury level set by RLuxGo : 4  p =  389
!!$!   RanLux initialized by RLuxGo from Seeds 1 0 0
!!$!   RanLux luxury p = 389, 1 - 5:
!!$!    0.9458949   0.4734785   0.9515279   0.4297197   9.1273844E-02
!!$!   RanLux luxury p = 389, 101 - 105:
!!$!    2.6182652E-02   3.7753463E-02   0.9727478   0.1330217   0.4312606
!!$!   call RLuxGo (75, 0, 0, 0)
!!$!   RanLux p-value set by RLuxGo to: 75
!!$!   RanLux initialized by RLuxGo from default seed
!!$!   RanLux luxury p = 75, 1 - 5:
!!$!    0.5398182   0.7615504   6.0299397E-02   0.7960026   0.3063122
!!$!   RanLux luxury p = 75, 101 - 105:
!!$!    0.2560073   0.2344321   0.5916438   0.5903584   7.0114136E-02
!!$!   Test restarting from the full vector
!!$!   Current RanLux status saved:
!!$!  16156027 16534309 15243811 2751687 6002207 7979506 1301976 4567313
!!$!  4305996 5872599 12003090 2146823 12606367 4111505 5979640 12739666
!!$!  10489318 14036909 11729352 8061448 7832659 6069758 3197719 1832730
!!$!  75080216
!!$!   RanLux numbers 1-5:
!!$!    0.2261783   0.6065599   0.8641744   0.4392008   0.2338251
!!$!   RanLux numbers 101 - 105:
!!$!    8.1071973E-02   0.2146685   0.8485673   0.9407805   0.8562623
!!$!   Previous RanLux status will be restored
!!$!   Full initialization of RanLux with 25 integers:
!!$!  16156027 16534309 15243811 2751687 6002207 7979506 1301976 4567313
!!$!  4305996 5872599 12003090 2146823 12606367 4111505 5979640 12739666
!!$!  10489318 14036909 11729352 8061448 7832659 6069758 3197719 1832730
!!$!  75080216
!!$!   RanLux p-value set by RLuxIn to: 75
!!$!   RanLux numbers 1-5:
!!$!    0.2261783   0.6065599   0.8641744   0.4392008   0.2338251
!!$!   RanLux numbers 101 - 105:
!!$!    8.1071973E-02   0.2146685   0.8485673   0.9407805   0.8562623
!!$!   Test the restarting by skipping
!!$!   RanLux luxury level set by RLuxGo : 4  p =  389
!!$!   RanLux initialized by RLuxGo from Seeds 7674985 0 0
!!$!   RLuxAt values =  4 7674985 0 0
!!$!   Vector length is now 1000
!!$!   RLuxAt values =  4 7674985 161840 0
!!$!   Vector length is now 200
!!$!   Next and 200th numbers are:   1.9647598E-02   0.5905859
!!$!   RanLux luxury level set by RLuxGo : 4  p =  389
!!$!   RanLux initialized by RLuxGo from Seeds 7674985 161840 0
!!$!   Next and 200th numbers are:   1.9647598E-02   0.5905859
!!$!   The following should provoke an error message
!!$!   RanLux luxury level set by RLuxGo : 4  p =  389
!!$!   RanLux initialized by RLuxGo from Seeds 11111 31 0
!!$!   Error in restarting with RLuxGo:
!!$! the values 11111 31 0  cannot occur at luxury level 4
!!$
