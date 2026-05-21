!=======================================================================
!  UMAT variant A: Linear elasticity + Norton-Bailey creep ONLY
!  Plasticity removed.
!
!  Based on the user's original coupled DP plasticity + Norton-Bailey UMAT.
!
!  PROPS kept compatible with the original input deck:
!    PROPS(1) = E
!    PROPS(2) = NU
!    PROPS(6) = A_CR
!    PROPS(7) = N_CR
!    PROPS(8) = M_CR
!    PROPS(9) = ECR0
!
!  STATEV kept compatible with the original output:
!    STATEV(1) = EPBAR_P   = 0 in this variant
!    STATEV(2) = EPBAR_CR  accumulated equivalent creep strain
!    STATEV(3) = Q         von Mises equivalent stress after creep
!    STATEV(4) = P         hydrostatic stress after creep
!    STATEV(5) = DECR_BAR  equivalent creep increment of last increment
!=======================================================================

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      INTEGER NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT
      INTEGER JSTEP(4),KINC
      DOUBLE PRECISION
     1 STRESS(NTENS),STATEV(NSTATV),
     2 DDSDDE(NTENS,NTENS),SSE,SPD,SCD,RPL,
     3 DDSDDT(NTENS),DRPLDE(NTENS),DRPLDT,
     4 STRAN(NTENS),DSTRAN(NTENS),TIME(2),DTIME,TEMP,DTEMP,
     5 PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),
     6 PNEWDT,CELENT,DFGRD0(3,3),DFGRD1(3,3)

      DOUBLE PRECISION E,NU,A_CR,N_CR,M_CR,ECR0
      DOUBLE PRECISION LAMBDA,MU,KBULK
      DOUBLE PRECISION STRESS_TR(6),STRESS_CR(6),S_TR(6),S_CR(6)
      DOUBLE PRECISION P_TR,Q_TR,P_CR,Q_CR
      DOUBLE PRECISION EPBAR_CR,DECR_BAR
      DOUBLE PRECISION Q_MID,ECR_MID,ECR_EFF_MID
      DOUBLE PRECISION RES,DRES,DECR_IT,ECRDOT
      DOUBLE PRECISION CR_FACT
      DOUBLE PRECISION DDSDDE_E(6,6)
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,HALF,ONETHIRD,TWO3
      DOUBLE PRECISION SQRT32,TOL_CR,TOL_NR,ECR_MIN
      INTEGER I,J,ITER,MAXITER_CR
      LOGICAL CR_ACTIVE

      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0)
      PARAMETER (HALF=0.5D0, ONETHIRD=1.D0/3.D0, TWO3=2.D0/3.D0)
      PARAMETER (SQRT32=1.224744871391589D0)
      PARAMETER (TOL_CR=1.D-10)
      PARAMETER (TOL_NR=1.D-12)
      PARAMETER (MAXITER_CR=100)
      PARAMETER (ECR_MIN=1.D-15)

!=======================================================================
!  1. MATERIAL PARAMETERS
!=======================================================================
      E    = PROPS(1)
      NU   = PROPS(2)

!     Keep original PROPS numbering for creep
      A_CR = ZERO
      N_CR = ZERO
      M_CR = ZERO
      ECR0 = ECR_MIN

      IF (NPROPS .GE. 6) A_CR = PROPS(6)
      IF (NPROPS .GE. 7) N_CR = PROPS(7)
      IF (NPROPS .GE. 8) M_CR = PROPS(8)
      IF (NPROPS .GE. 9) ECR0 = PROPS(9)
      IF (ECR0 .LT. ECR_MIN) ECR0 = ECR_MIN

      LAMBDA = E*NU / ((ONE+NU)*(ONE-TWO*NU))
      MU     = E / (TWO*(ONE+NU))
      KBULK  = E / (THREE*(ONE-TWO*NU))

!=======================================================================
!  2. READ STATE VARIABLES
!=======================================================================
      EPBAR_CR = ZERO
      IF (NSTATV .GE. 2) EPBAR_CR = STATEV(2)
      IF (EPBAR_CR .LT. ZERO) EPBAR_CR = ZERO

!=======================================================================
!  3. ELASTIC STIFFNESS
!=======================================================================
      DO I = 1, NTENS
        DO J = 1, NTENS
          DDSDDE_E(I,J) = ZERO
        END DO
      END DO

      DO I = 1, NDI
        DO J = 1, NDI
          DDSDDE_E(I,J) = LAMBDA
        END DO
        DDSDDE_E(I,I) = LAMBDA + TWO*MU
      END DO

      DO I = NDI+1, NTENS
        DDSDDE_E(I,I) = MU
      END DO

!=======================================================================
!  4. ELASTIC TRIAL STRESS
!=======================================================================
      DO I = 1, NTENS
        STRESS_TR(I) = STRESS(I)
        DO J = 1, NTENS
          STRESS_TR(I) = STRESS_TR(I) + DDSDDE_E(I,J)*DSTRAN(J)
        END DO
      END DO

      P_TR = ZERO
      DO I = 1, NDI
        P_TR = P_TR + ONETHIRD*STRESS_TR(I)
      END DO

      DO I = 1, NDI
        S_TR(I) = STRESS_TR(I) - P_TR
      END DO
      DO I = NDI+1, NTENS
        S_TR(I) = STRESS_TR(I)
      END DO

      Q_TR = ZERO
      DO I = 1, NDI
        Q_TR = Q_TR + S_TR(I)*S_TR(I)
      END DO
      DO I = NDI+1, NTENS
        Q_TR = Q_TR + TWO*S_TR(I)*S_TR(I)
      END DO
      Q_TR = SQRT32*DSQRT(Q_TR)

!=======================================================================
!  5. NORTON-BAILEY CREEP CORRECTION ONLY
!     DECR_BAR = A * q_mid^n * (EPBAR_CR + 0.5*DECR_BAR + ECR0)^m * DTIME
!=======================================================================
      DECR_BAR = ZERO

      CR_ACTIVE = (DTIME .GT. ZERO) .AND.
     1            (A_CR .GT. ZERO) .AND.
     2            (Q_TR .GT. TOL_CR)

      IF (CR_ACTIVE) THEN

        Q_MID = Q_TR
        ECR_MID = EPBAR_CR
        ECR_EFF_MID = ECR_MID + ECR0
        IF (ECR_EFF_MID .LT. ECR_MIN) ECR_EFF_MID = ECR_MIN

        DECR_BAR = A_CR*(Q_MID**N_CR)*(ECR_EFF_MID**M_CR)*DTIME

        IF (DECR_BAR .GT. Q_TR/(THREE*MU+TOL_CR)) THEN
          DECR_BAR = Q_TR/(THREE*MU+TOL_CR)*0.9D0
        END IF

        RES = ZERO

        DO ITER = 1, MAXITER_CR

          Q_MID = Q_TR - THREE*MU*DECR_BAR
          ECR_MID = EPBAR_CR + HALF*DECR_BAR
          ECR_EFF_MID = ECR_MID + ECR0

          IF (Q_MID .LT. TOL_CR) Q_MID = TOL_CR
          IF (ECR_EFF_MID .LT. ECR_MIN) ECR_EFF_MID = ECR_MIN

          ECRDOT = A_CR*(Q_MID**N_CR)*(ECR_EFF_MID**M_CR)
          RES = DECR_BAR - ECRDOT*DTIME

          IF (DABS(RES) .LT. TOL_NR*(ONE + DECR_BAR)) EXIT

          DRES = ONE - DTIME*(
     1      A_CR*N_CR*(Q_MID**(N_CR-ONE))*(-THREE*MU)
     2           *(ECR_EFF_MID**M_CR)
     3    + A_CR*(Q_MID**N_CR)*M_CR
     4           *(ECR_EFF_MID**(M_CR-ONE))*HALF )

          IF (DABS(DRES) .LT. TOL_CR) DRES = ONE

          DECR_IT = -RES/DRES

          IF (DECR_IT .GT. HALF*DECR_BAR+TOL_CR)
     1      DECR_IT = HALF*(DECR_BAR+TOL_CR)
          IF (DECR_IT .LT. -HALF*DECR_BAR)
     1      DECR_IT = -HALF*DECR_BAR

          DECR_BAR = DECR_BAR + DECR_IT

          IF (DECR_BAR .LT. ZERO) DECR_BAR = ZERO
          IF (DECR_BAR .GT. Q_TR/(THREE*MU+TOL_CR))
     1      DECR_BAR = Q_TR/(THREE*MU+TOL_CR)*0.99D0

        END DO

        IF (DABS(RES) .GT. TOL_NR*(ONE + DECR_BAR)) THEN
          PNEWDT = 0.25D0
          RETURN
        END IF

      END IF

!=======================================================================
!  6. STRESS AFTER CREEP
!=======================================================================
      DO I = 1, NTENS
        STRESS_CR(I) = STRESS_TR(I)
      END DO

      IF (DECR_BAR .GT. ZERO .AND. Q_TR .GT. TOL_CR) THEN
        CR_FACT = THREE*MU*DECR_BAR/Q_TR
        DO I = 1, NTENS
          STRESS_CR(I) = STRESS_TR(I) - CR_FACT*S_TR(I)
        END DO
      END IF

      P_CR = ZERO
      DO I = 1, NDI
        P_CR = P_CR + ONETHIRD*STRESS_CR(I)
      END DO

      DO I = 1, NDI
        S_CR(I) = STRESS_CR(I) - P_CR
      END DO
      DO I = NDI+1, NTENS
        S_CR(I) = STRESS_CR(I)
      END DO

      Q_CR = ZERO
      DO I = 1, NDI
        Q_CR = Q_CR + S_CR(I)*S_CR(I)
      END DO
      DO I = NDI+1, NTENS
        Q_CR = Q_CR + TWO*S_CR(I)*S_CR(I)
      END DO
      Q_CR = SQRT32*DSQRT(Q_CR)

      DO I = 1, NTENS
        STRESS(I) = STRESS_CR(I)
      END DO

!=======================================================================
!  7. UPDATE STATE VARIABLES
!=======================================================================
      IF (NSTATV .GE. 1) STATEV(1) = ZERO
      IF (NSTATV .GE. 2) STATEV(2) = EPBAR_CR + DECR_BAR
      IF (NSTATV .GE. 3) STATEV(3) = Q_CR
      IF (NSTATV .GE. 4) STATEV(4) = P_CR
      IF (NSTATV .GE. 5) STATEV(5) = DECR_BAR

!=======================================================================
!  8. CREEP SECANT TANGENT
!=======================================================================
      IF (Q_TR .GT. TOL_CR .AND. DECR_BAR .GT. ZERO) THEN
        CR_FACT = ONE - THREE*MU*DECR_BAR/Q_TR
        IF (CR_FACT .LT. ZERO) CR_FACT = ZERO
      ELSE
        CR_FACT = ONE
      END IF

      DO I = 1, NTENS
        DO J = 1, NTENS
          DDSDDE(I,J) = ZERO
        END DO
      END DO

      DO I = 1, NDI
        DO J = 1, NDI
          DDSDDE(I,J) = KBULK
        END DO
      END DO

      DO I = 1, NDI
        DDSDDE(I,I) = DDSDDE(I,I) + TWO*MU*CR_FACT
        DO J = 1, NDI
          DDSDDE(I,J) = DDSDDE(I,J)
     1                 - TWO*MU*CR_FACT*ONETHIRD
        END DO
      END DO

      DO I = NDI+1, NTENS
        DDSDDE(I,I) = MU*CR_FACT
      END DO

      RETURN
      END
