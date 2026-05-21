!=======================================================================
!  UMAT variant B: Linear elasticity + Drucker-Prager plasticity ONLY
!  Creep removed.
!
!  Based on the user's original coupled DP plasticity + Norton-Bailey UMAT.
!
!  PROPS kept compatible with the original input deck:
!    PROPS(1) = E
!    PROPS(2) = NU
!    PROPS(3) = ALPHA_DP
!    PROPS(4) = K_COH
!    PROPS(5) = PSI_DP
!    PROPS(6:9) ignored in this variant
!
!  STATEV kept compatible with the original output:
!    STATEV(1) = EPBAR_P   accumulated equivalent plastic strain
!    STATEV(2) = EPBAR_CR  = 0 in this variant
!    STATEV(3) = Q         von Mises equivalent stress before/after DP check
!    STATEV(4) = P         hydrostatic stress before/after DP check
!    STATEV(5) = DECR_BAR  = 0 in this variant
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

      DOUBLE PRECISION E,NU,ALPHA_DP,K_COH,PSI_DP
      DOUBLE PRECISION LAMBDA,MU,KBULK
      DOUBLE PRECISION STRESS_TR(6),S_TR(6),S_N(6)
      DOUBLE PRECISION P_TR,Q_TR,P_N,Q_N
      DOUBLE PRECISION EPBAR_P,DEPBAR_P,DGAMMA,F_TR,BETA_PL
      DOUBLE PRECISION DDSDDE_E(6,6)
      DOUBLE PRECISION THETA_PL
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,ONETHIRD,TWO3
      DOUBLE PRECISION SQRT32,TOL_PL
      INTEGER I,J
      LOGICAL PLASTIC

      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0)
      PARAMETER (ONETHIRD=1.D0/3.D0, TWO3=2.D0/3.D0)
      PARAMETER (SQRT32=1.224744871391589D0)
      PARAMETER (TOL_PL=1.D-10)

!=======================================================================
!  1. MATERIAL PARAMETERS
!=======================================================================
      E        = PROPS(1)
      NU       = PROPS(2)
      ALPHA_DP = PROPS(3)
      K_COH    = PROPS(4)
      PSI_DP   = PROPS(5)

!     Same convention as original code:
!     if PSI_DP = 0, use associated flow, PSI_DP = ALPHA_DP
      IF (DABS(PSI_DP) .LT. TOL_PL) PSI_DP = ALPHA_DP

      LAMBDA = E*NU / ((ONE+NU)*(ONE-TWO*NU))
      MU     = E / (TWO*(ONE+NU))
      KBULK  = E / (THREE*(ONE-TWO*NU))

!=======================================================================
!  2. READ STATE VARIABLES
!=======================================================================
      EPBAR_P = ZERO
      IF (NSTATV .GE. 1) EPBAR_P = STATEV(1)
      IF (EPBAR_P .LT. ZERO) EPBAR_P = ZERO

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
!  5. DRUCKER-PRAGER YIELD CHECK
!     Same sign convention as the original code:
!       F = q + ALPHA_DP*p - K_COH
!=======================================================================
      F_TR = Q_TR + ALPHA_DP*P_TR - K_COH

      IF (F_TR .LE. TOL_PL) THEN
        PLASTIC = .FALSE.
      ELSE
        PLASTIC = .TRUE.
      END IF

!=======================================================================
!  6. PLASTIC RETURN MAPPING ONLY
!=======================================================================
      IF (.NOT. PLASTIC) THEN

        DO I = 1, NTENS
          STRESS(I) = STRESS_TR(I)
        END DO
        DEPBAR_P = ZERO
        DGAMMA   = ZERO
        P_N      = P_TR
        Q_N      = Q_TR

      ELSE

        IF (Q_TR .LT. TOL_PL) THEN

          IF (DABS(ALPHA_DP) .GT. TOL_PL .AND.
     1        DABS(PSI_DP) .GT. TOL_PL) THEN
            DGAMMA = -(P_TR + K_COH/ALPHA_DP)/(PSI_DP*KBULK)
          ELSE
            DGAMMA = ZERO
          END IF

          IF (DGAMMA .LT. ZERO) DGAMMA = ZERO

          DO I = 1, NDI
            STRESS(I) = STRESS_TR(I) - KBULK*PSI_DP*DGAMMA
          END DO
          DO I = NDI+1, NTENS
            STRESS(I) = ZERO
          END DO

          P_N = ZERO
          DO I = 1, NDI
            P_N = P_N + ONETHIRD*STRESS(I)
          END DO
          Q_N = ZERO

        ELSE

          DGAMMA = F_TR/(THREE*MU + ALPHA_DP*PSI_DP*KBULK)

          IF (DGAMMA .LT. ZERO) DGAMMA = ZERO

          BETA_PL = ONE - THREE*MU*DGAMMA/Q_TR
          IF (BETA_PL .LT. ZERO) BETA_PL = ZERO

          P_N = P_TR - KBULK*PSI_DP*DGAMMA

          DO I = 1, NDI
            S_N(I) = BETA_PL*S_TR(I)
            STRESS(I) = S_N(I) + P_N
          END DO
          DO I = NDI+1, NTENS
            S_N(I) = BETA_PL*S_TR(I)
            STRESS(I) = S_N(I)
          END DO

          Q_N = Q_TR - THREE*MU*DGAMMA
          IF (Q_N .LT. ZERO) Q_N = ZERO

        END IF

        DEPBAR_P = DGAMMA

      END IF

!=======================================================================
!  7. UPDATE STATE VARIABLES
!=======================================================================
      IF (NSTATV .GE. 1) STATEV(1) = EPBAR_P + DEPBAR_P
      IF (NSTATV .GE. 2) STATEV(2) = ZERO
      IF (NSTATV .GE. 3) STATEV(3) = Q_N
      IF (NSTATV .GE. 4) STATEV(4) = P_N
      IF (NSTATV .GE. 5) STATEV(5) = ZERO

!=======================================================================
!  8. ELASTO-PLASTIC SECANT / APPROXIMATE TANGENT
!=======================================================================
      DO I = 1, NTENS
        DO J = 1, NTENS
          DDSDDE(I,J) = DDSDDE_E(I,J)
        END DO
      END DO

      IF (PLASTIC .AND. Q_TR .GT. TOL_PL .AND. DGAMMA .GT. ZERO) THEN

        THETA_PL = ONE/(THREE*MU + ALPHA_DP*PSI_DP*KBULK)

!       Start from elastoplastic volumetric/deviatoric approximation
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

        BETA_PL = ONE - THREE*MU*DGAMMA/Q_TR
        IF (BETA_PL .LT. ZERO) BETA_PL = ZERO

        DO I = 1, NDI
          DDSDDE(I,I) = DDSDDE(I,I) + TWO*MU*BETA_PL
          DO J = 1, NDI
            DDSDDE(I,J) = DDSDDE(I,J)
     1                   - TWO*MU*BETA_PL*ONETHIRD
          END DO
        END DO

        DO I = NDI+1, NTENS
          DDSDDE(I,I) = MU*BETA_PL
        END DO

!       Rank-one correction in the plastic loading direction
        DO I = 1, NTENS
          DO J = 1, NTENS
            DDSDDE(I,J) = DDSDDE(I,J)
     1        - (THREE*MU)**2*THETA_PL
     2        * S_TR(I)*S_TR(J)*TWO3/(Q_TR**2)
          END DO
        END DO

!       Volumetric-deviatoric coupling from friction/dilatancy
        DO I = 1, NDI
          DO J = 1, NTENS
            DDSDDE(I,J) = DDSDDE(I,J)
     1        - KBULK*PSI_DP*THREE*MU*THETA_PL
     2        * S_TR(J)*TWO3/Q_TR
          END DO

          DO J = 1, NDI
            DDSDDE(I,J) = DDSDDE(I,J)
     1        - THREE*MU*THETA_PL*KBULK*ALPHA_DP
     2        * S_TR(I)*TWO3/Q_TR
          END DO
        END DO

      END IF

      RETURN
      END
