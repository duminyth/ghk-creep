!=======================================================================
!  UMAT COPY: Drucker-Prager plasticity + Norton-Bailey creep
!             with simple internal time substepping
!
!  Purpose:
!    Diagnostic version for time-discretization sensitivity.
!
!  Main change compared with the original operator split:
!    The Abaqus increment DSTRAN, DTIME is divided into NSUB local
!    subincrements. On each subincrement the same sequence is applied:
!
!       elastic trial -> creep correction -> DP plastic correction
!
!  This does not make the algorithm fully coupled, but it strongly
!  reduces the operator-splitting error when creep and plasticity are
!  both active.
!
!  PROPS:
!    1  E
!    2  NU
!    3  ALPHA_DP
!    4  K_COH
!    5  PSI_DP
!    6  A_CR
!    7  N_CR
!    8  M_CR
!    9  ECR0
!   10  optional: DTIME_MAX_LOCAL  [s]
!   11  optional: DECR_MAX_SUB     [-]
!
!  STATEV diagnostics:
!    1  PEEQ
!    2  CEEQ
!    3  q after creep, last substep
!    4  p after creep, last substep
!    5  dCEEQ total over Abaqus increment
!    6  DTIME Abaqus increment
!    7  NSUB used
!    8  q before creep, last substep
!    9  q after creep, last substep
!   10  q after plasticity, last substep
!   11  p before creep, last substep
!   12  p after creep, last substep
!   13  p after plasticity, last substep
!   14  yield F before plasticity, last substep
!   15  yield F after plasticity, last substep
!   16  dCEEQ last substep
!   17  dPEEQ last substep
!   18  creep active flag, last substep
!   19  plastic active flag, last substep
!   20  first creep onset time, if NSTATV >= 20
!   21  first plastic onset time, if NSTATV >= 21
!=======================================================================

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
      RPL,DDSDDT,DRPLDE,DRPLDT, &
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
      NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT, &
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      INTEGER NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT
      INTEGER JSTEP(4),KINC

      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV)
      DOUBLE PRECISION DDSDDE(NTENS,NTENS),SSE,SPD,SCD,RPL
      DOUBLE PRECISION DDSDDT(NTENS),DRPLDE(NTENS),DRPLDT
      DOUBLE PRECISION STRAN(NTENS),DSTRAN(NTENS),TIME(2),DTIME
      DOUBLE PRECISION TEMP,DTEMP,PREDEF(1),DPRED(1),PROPS(NPROPS)
      DOUBLE PRECISION COORDS(3),DROT(3,3),PNEWDT,CELENT
      DOUBLE PRECISION DFGRD0(3,3),DFGRD1(3,3)

      DOUBLE PRECISION E,NU,ALPHA_DP,K_COH,PSI_DP
      DOUBLE PRECISION A_COMP,N_COMP,M_COMP
      DOUBLE PRECISION A_SHEAR,N_SHEAR,M_SHEAR
      DOUBLE PRECISION ECR0
      DOUBLE PRECISION THETA_SWITCH_DEG,THETA_SHEAR_DEG
      DOUBLE PRECISION THETA_PQ,W_SHEAR,P_COMP,PI
      DOUBLE PRECISION ECRDOT_COMP,ECRDOT_SHEAR
      DOUBLE PRECISION B_POW
      DOUBLE PRECISION LAMBDA,MU,KBULK
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,HALF,ONETHIRD,TWO3
      DOUBLE PRECISION SQRT32,TOL_CR,TOL_PL,TOL_NR,ECR_MIN
      DOUBLE PRECISION DTIME_MAX_LOCAL,DECR_MAX_SUB
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0)
      PARAMETER (HALF=0.5D0,ONETHIRD=1.D0/3.D0,TWO3=2.D0/3.D0)
      PARAMETER (SQRT32=1.224744871391589D0)
      PARAMETER (TOL_CR=1.D-12,TOL_PL=1.D-10,TOL_NR=1.D-12)
      PARAMETER (ECR_MIN=1.D-15)
      PARAMETER (PI=3.141592653589793D0)

      INTEGER I,J,ISUB,NSUB,ITER,MAXITER_CR
      PARAMETER (MAXITER_CR=80)

      DOUBLE PRECISION DDSDDE_E(6,6),DELTA(6)
      DOUBLE PRECISION DSTRAN_SUB(6),STRESS_TR(6),STRESS_CR(6)
      DOUBLE PRECISION S_TR(6),S_CR(6),S_N(6)
      DOUBLE PRECISION P_TR,Q_TR,P_CR,Q_CR,P_N,Q_N
      DOUBLE PRECISION P_END,Q_END
      DOUBLE PRECISION EPBAR_P,EPBAR_CR,DEPBAR_P,DECR_BAR
      DOUBLE PRECISION DEPBAR_P_TOT,DECR_BAR_TOT
      DOUBLE PRECISION DGAMMA,F_TR_CR,F_END,BETA_PL
      DOUBLE PRECISION Q_MID,ECR_MID,ECR_EFF_MID,ECRDOT
      DOUBLE PRECISION RES,DRES,DECR_IT,CR_FACT
      DOUBLE PRECISION DTSUB
      DOUBLE PRECISION CR_FACT_LAST,THETA_PL
      LOGICAL CR_ACTIVE,PLASTIC

!=======================================================================
!  1. Material parameters
!=======================================================================
      E        = PROPS(1)
      NU       = PROPS(2)
      ALPHA_DP = PROPS(3)
      K_COH    = PROPS(4)
      PSI_DP   = PROPS(5)

      A_COMP  = PROPS(6)
      N_COMP  = PROPS(7)
      M_COMP  = PROPS(8)

      A_SHEAR = PROPS(9)
      N_SHEAR = PROPS(10)
      M_SHEAR = PROPS(11)

      B_POW   = PROPS(12)

      ECR0    = PROPS(13)

      THETA_SWITCH_DEG = 71.57D0
      THETA_SHEAR_DEG  = 90.D0


      IF (ECR0 .LT. ECR_MIN) ECR0          = ECR_MIN
      IF (DABS(PSI_DP) .LT. TOL_PL) PSI_DP = ALPHA_DP

      DTIME_MAX_LOCAL = DTIME
      DECR_MAX_SUB    = 1.D-12


      IF (NPROPS .GE. 14) THEN
        IF (PROPS(14) .GT. ZERO) DTIME_MAX_LOCAL = PROPS(14)
      END IF
      IF (NPROPS .GE. 15) THEN
        IF (PROPS(15) .GT. ZERO) DECR_MAX_SUB = PROPS(15)
      END IF

      LAMBDA = E*NU / ((ONE+NU)*(ONE-TWO*NU))
      MU     = E / (TWO*(ONE+NU))
      KBULK  = E / (THREE*(ONE-TWO*NU))

      DO I=1,6
        DELTA(I)=ZERO
      END DO
      DO I=1,NDI
        DELTA(I)=ONE
      END DO

!=======================================================================
!  2. Elastic stiffness
!=======================================================================
      DO I=1,6
        DO J=1,6
          DDSDDE_E(I,J)=ZERO
        END DO
      END DO

      DO I=1,NDI
        DO J=1,NDI
          DDSDDE_E(I,J)=LAMBDA
        END DO
        DDSDDE_E(I,I)=LAMBDA+TWO*MU
      END DO
      DO I=NDI+1,NTENS
        DDSDDE_E(I,I)=MU
      END DO

!=======================================================================
!  3. Initialize state and substepping
!=======================================================================
      EPBAR_P  = STATEV(1)
      EPBAR_CR = STATEV(2)
      IF (EPBAR_P  .LT. ZERO) EPBAR_P  = ZERO
      IF (EPBAR_CR .LT. ZERO) EPBAR_CR = ZERO

      DEPBAR_P_TOT = ZERO
      DECR_BAR_TOT = ZERO

      NSUB = 1
      IF (DTIME .GT. ZERO .AND. DTIME_MAX_LOCAL .GT. ZERO) THEN
        NSUB = INT(DTIME/DTIME_MAX_LOCAL)
        IF (DBLE(NSUB)*DTIME_MAX_LOCAL .LT. DTIME) NSUB=NSUB+1
        IF (NSUB .LT. 1) NSUB=1
      END IF



!     First estimate based on explicit creep at beginning of increment.
!     This is deliberately conservative.
      DO I=1,NTENS
        STRESS_TR(I)=STRESS(I)
        DO J=1,NTENS
          STRESS_TR(I)=STRESS_TR(I)+DDSDDE_E(I,J)*DSTRAN(J)
        END DO
      END DO
      CALL GET_PQ(NDI,NTENS,STRESS_TR,P_TR,Q_TR,S_TR)

      P_COMP = -P_TR

      IF (P_COMP .LE. TOL_CR) THEN
      THETA_PQ = 90.D0
      ELSE
      THETA_PQ = DATAN2(Q_TR,P_COMP) * 180.D0 / PI
      END IF

      IF (THETA_PQ .LE. THETA_SWITCH_DEG) THEN
      W_SHEAR = ZERO
      ELSE
      W_SHEAR = ((THETA_PQ - THETA_SWITCH_DEG) / &
              (THETA_SHEAR_DEG - THETA_SWITCH_DEG))**B_POW
      IF (W_SHEAR .LT. ZERO) W_SHEAR = ZERO
      IF (W_SHEAR .GT. ONE)  W_SHEAR = ONE
      END IF



        IF ((A_COMP .GT. ZERO .OR. A_SHEAR .GT. ZERO) .AND. &
            DTIME .GT. ZERO .AND. Q_TR .GT. TOL_CR) THEN
        ECR_EFF_MID = EPBAR_CR + ECR0
        IF (ECR_EFF_MID .LT. ECR_MIN) ECR_EFF_MID=ECR_MIN
        
        P_COMP = -P_TR

        IF (P_COMP .LE. TOL_CR) THEN
        THETA_PQ = 90.D0
        ELSE
        THETA_PQ = DATAN2(Q_TR,P_COMP) * 180.D0 / PI
        END IF

        IF (THETA_PQ .LE. THETA_SWITCH_DEG) THEN
        W_SHEAR = ZERO
        ELSE
        W_SHEAR = ((THETA_PQ - THETA_SWITCH_DEG) / &
                (THETA_SHEAR_DEG - THETA_SWITCH_DEG))**B_POW
        IF (W_SHEAR .LT. ZERO) W_SHEAR = ZERO
        IF (W_SHEAR .GT. ONE)  W_SHEAR = ONE
        END IF

        ECRDOT_COMP = &
             A_COMP*(Q_TR**N_COMP)*(ECR_EFF_MID**M_COMP) 

        ECRDOT_SHEAR = &
              A_SHEAR*(Q_TR**N_SHEAR)*(ECR_EFF_MID**M_SHEAR)

        DECR_BAR = &
              ((ONE-W_SHEAR)*ECRDOT_COMP &
              + W_SHEAR*ECRDOT_SHEAR)*DTIME


        IF (DECR_BAR .GT. DECR_MAX_SUB) THEN
          NSUB = MAX(NSUB,INT(DECR_BAR/DECR_MAX_SUB)+1)
        END IF
      END IF

      DTSUB = DTIME/DBLE(NSUB)
      DO I=1,NTENS
        DSTRAN_SUB(I)=DSTRAN(I)/DBLE(NSUB)
      END DO

!=======================================================================
!  4. Local substepping loop
!=======================================================================
      DO ISUB=1,NSUB

!       Elastic trial for this substep
        DO I=1,NTENS
          STRESS_TR(I)=STRESS(I)
          DO J=1,NTENS
            STRESS_TR(I)=STRESS_TR(I)+DDSDDE_E(I,J)*DSTRAN_SUB(J)
          END DO
        END DO

        CALL GET_PQ(NDI,NTENS,STRESS_TR,P_TR,Q_TR,S_TR)

        P_COMP = -P_TR

        IF (P_COMP .LE. TOL_CR) THEN
        THETA_PQ = 90.D0
        ELSE
        THETA_PQ = DATAN2(Q_TR,P_COMP) * 180.D0 / PI
        END IF

        IF (THETA_PQ .LE. THETA_SWITCH_DEG) THEN
        W_SHEAR = ZERO
        ELSE
        W_SHEAR = ((THETA_PQ - THETA_SWITCH_DEG) / &
               (THETA_SHEAR_DEG - THETA_SWITCH_DEG))**B_POW
        IF (W_SHEAR .LT. ZERO) W_SHEAR = ZERO
        IF (W_SHEAR .GT. ONE)  W_SHEAR = ONE
        END IF

!       Creep correction
        CR_ACTIVE = (DTSUB .GT. ZERO) .AND. &
                    (A_COMP .GT. ZERO .OR. A_SHEAR .GT. ZERO) .AND. &
                    (Q_TR .GT. TOL_CR)
        DECR_BAR = ZERO

        IF (CR_ACTIVE) THEN
          Q_MID = Q_TR
          ECR_EFF_MID = EPBAR_CR + ECR0
          IF (ECR_EFF_MID .LT. ECR_MIN) ECR_EFF_MID=ECR_MIN


          ECRDOT_COMP = A_COMP*(Q_MID**N_COMP)*(ECR_EFF_MID**M_COMP)

          ECRDOT_SHEAR = A_SHEAR*(Q_MID**N_SHEAR)*(ECR_EFF_MID**M_SHEAR)

          ECRDOT = (ONE-W_SHEAR)*ECRDOT_COMP + W_SHEAR*ECRDOT_SHEAR

          DECR_BAR = ECRDOT*DTSUB

          IF (DECR_BAR .GT. Q_TR/(THREE*MU+TOL_CR)) THEN
            DECR_BAR = Q_TR/(THREE*MU+TOL_CR)*0.9D0
          END IF

          DO ITER=1,MAXITER_CR
            Q_MID = Q_TR - THREE*MU*DECR_BAR
            IF (Q_MID .LT. TOL_CR) Q_MID = TOL_CR

            ECR_MID = EPBAR_CR + HALF*DECR_BAR
            ECR_EFF_MID = ECR_MID + ECR0
            IF (ECR_EFF_MID .LT. ECR_MIN) ECR_EFF_MID=ECR_MIN

            ECRDOT_COMP = &
                  A_COMP*(Q_MID**N_COMP)*(ECR_EFF_MID**M_COMP)

            ECRDOT_SHEAR = &
                  A_SHEAR*(Q_MID**N_SHEAR)*(ECR_EFF_MID**M_SHEAR)

            ECRDOT = (ONE-W_SHEAR)*ECRDOT_COMP &
                    + W_SHEAR*ECRDOT_SHEAR


            RES = DECR_BAR - ECRDOT*DTSUB

            IF (DABS(RES) .LT. TOL_NR*(ONE+DECR_BAR)) GOTO 100

            DRES = ONE - DTSUB*( &
                  (ONE-W_SHEAR)* &
                  ( &
                    A_COMP*N_COMP*(Q_MID**(N_COMP-ONE))*(-THREE*MU) &
                    *(ECR_EFF_MID**M_COMP) &
                  + A_COMP*(Q_MID**N_COMP)*M_COMP &
                    *(ECR_EFF_MID**(M_COMP-ONE))*HALF &
                  ) &
                  + W_SHEAR* &
                  ( &
                    A_SHEAR*N_SHEAR*(Q_MID**(N_SHEAR-ONE))*(-THREE*MU) &
                    *(ECR_EFF_MID**M_SHEAR) &
                  + A_SHEAR*(Q_MID**N_SHEAR)*M_SHEAR &
                    *(ECR_EFF_MID**(M_SHEAR-ONE))*HALF &
                  ) &
                  )

            IF (DABS(DRES) .LT. TOL_CR) DRES=ONE
            DECR_IT = -RES/DRES

            IF (DECR_IT .GT. HALF*(DECR_BAR+TOL_CR)) THEN
              DECR_IT = HALF*(DECR_BAR+TOL_CR)
            END IF
            IF (DECR_IT .LT. -HALF*DECR_BAR) THEN
              DECR_IT = -HALF*DECR_BAR
            END IF

            DECR_BAR = DECR_BAR + DECR_IT
            IF (DECR_BAR .LT. ZERO) DECR_BAR=ZERO
            IF (DECR_BAR .GT. Q_TR/(THREE*MU+TOL_CR)) THEN
              DECR_BAR = Q_TR/(THREE*MU+TOL_CR)*0.99D0
            END IF
          END DO

!         If local Newton still fails, ask Abaqus for smaller increment.
          PNEWDT = MIN(PNEWDT,0.25D0)
          RETURN
 100      CONTINUE
        END IF

!       Stress after creep correction
        DO I=1,NTENS
          STRESS_CR(I)=STRESS_TR(I)
        END DO
        IF (DECR_BAR .GT. ZERO .AND. Q_TR .GT. TOL_CR) THEN
          CR_FACT = THREE*MU*DECR_BAR/Q_TR
          DO I=1,NTENS
            STRESS_CR(I)=STRESS_TR(I)-CR_FACT*S_TR(I)
          END DO
        END IF

        CALL GET_PQ(NDI,NTENS,STRESS_CR,P_CR,Q_CR,S_CR)

!       Drucker-Prager yield check after creep correction
        F_TR_CR = Q_CR + ALPHA_DP*P_CR - K_COH
        PLASTIC = F_TR_CR .GT. TOL_PL
        DGAMMA = ZERO
        DEPBAR_P = ZERO

        IF (.NOT. PLASTIC) THEN
          DO I=1,NTENS
            STRESS(I)=STRESS_CR(I)
          END DO
        ELSE
          IF (Q_CR .LT. TOL_PL) THEN
            DGAMMA = -(P_CR + K_COH/ALPHA_DP)/(PSI_DP*KBULK)
            IF (DGAMMA .LT. ZERO) DGAMMA=ZERO
            DO I=1,NDI
              STRESS(I)=STRESS_CR(I)-KBULK*PSI_DP*DGAMMA
            END DO
            DO I=NDI+1,NTENS
              STRESS(I)=ZERO
            END DO
          ELSE
            DGAMMA = F_TR_CR/(THREE*MU+ALPHA_DP*PSI_DP*KBULK)
            IF (DGAMMA .LT. ZERO) DGAMMA=ZERO
            BETA_PL = ONE - THREE*MU*DGAMMA/Q_CR
            P_N = P_CR - KBULK*PSI_DP*DGAMMA
            DO I=1,NDI
              S_N(I)=BETA_PL*S_CR(I)
              STRESS(I)=S_N(I)+P_N
            END DO
            DO I=NDI+1,NTENS
              S_N(I)=BETA_PL*S_CR(I)
              STRESS(I)=S_N(I)
            END DO
          END IF
          DEPBAR_P = DGAMMA
        END IF

!       Update state locally after this substep
        EPBAR_CR = EPBAR_CR + DECR_BAR
        EPBAR_P  = EPBAR_P  + DEPBAR_P
        DECR_BAR_TOT  = DECR_BAR_TOT  + DECR_BAR
        DEPBAR_P_TOT = DEPBAR_P_TOT + DEPBAR_P

        CALL GET_PQ(NDI,NTENS,STRESS,P_END,Q_END,S_N)
        F_END = Q_END + ALPHA_DP*P_END - K_COH

!       Store onset time once.
        IF (NSTATV .GE. 20) THEN
          IF (STATEV(20) .EQ. ZERO .AND. DECR_BAR .GT. TOL_CR) THEN
            STATEV(20)=TIME(2)
          END IF
        END IF
        IF (NSTATV .GE. 21) THEN
          IF (STATEV(21) .EQ. ZERO .AND. DEPBAR_P .GT. TOL_PL) THEN
            STATEV(21)=TIME(2)
          END IF
        END IF

      END DO

!=======================================================================
!  5. Final state variables and diagnostics
!=======================================================================
      STATEV(1)=EPBAR_P
      STATEV(2)=EPBAR_CR
      STATEV(3)=Q_CR
      STATEV(4)=P_CR
      STATEV(5)=DECR_BAR_TOT
      IF (NSTATV .GE. 6)  STATEV(6)=DTIME
      IF (NSTATV .GE. 7)  STATEV(7)=DBLE(NSUB)
      IF (NSTATV .GE. 8)  STATEV(8)=Q_TR
      IF (NSTATV .GE. 9)  STATEV(9)=Q_CR
      IF (NSTATV .GE. 10) STATEV(10)=Q_END
      IF (NSTATV .GE. 11) STATEV(11)=P_TR
      IF (NSTATV .GE. 12) STATEV(12)=P_CR
      IF (NSTATV .GE. 13) STATEV(13)=P_END
      IF (NSTATV .GE. 14) STATEV(14)=F_TR_CR
      IF (NSTATV .GE. 15) STATEV(15)=F_END
      IF (NSTATV .GE. 16) STATEV(16)=DECR_BAR
      IF (NSTATV .GE. 17) STATEV(17)=DEPBAR_P
      IF (NSTATV .GE. 18) THEN
        IF (DECR_BAR .GT. TOL_CR) THEN
          STATEV(18)=ONE
        ELSE
          STATEV(18)=ZERO
        END IF
      END IF
      IF (NSTATV .GE. 19) THEN
        IF (DEPBAR_P .GT. TOL_PL) THEN
          STATEV(19)=ONE
        ELSE
          STATEV(19)=ZERO
        END IF
      END IF

!=======================================================================
!  6. Approximate tangent
!     For this diagnostic version, use a safe elastic/secant-like tangent.
!     This is not the fully consistent coupled tangent.
!=======================================================================
      CR_FACT_LAST = ONE
      IF (Q_TR .GT. TOL_CR .AND. DECR_BAR .GT. ZERO) THEN
        CR_FACT_LAST = ONE - THREE*MU*DECR_BAR/Q_TR
        IF (CR_FACT_LAST .LT. ZERO) CR_FACT_LAST=ZERO
      END IF

      DO I=1,NTENS
        DO J=1,NTENS
          DDSDDE(I,J)=ZERO
        END DO
      END DO

      DO I=1,NDI
        DO J=1,NDI
          DDSDDE(I,J)=KBULK
        END DO
      END DO

      DO I=1,NDI
        DDSDDE(I,I)=DDSDDE(I,I)+TWO*MU*CR_FACT_LAST
        DO J=1,NDI
          DDSDDE(I,J)=DDSDDE(I,J)-TWO*MU*CR_FACT_LAST*ONETHIRD
        END DO
      END DO
      DO I=NDI+1,NTENS
        DDSDDE(I,I)=MU*CR_FACT_LAST
      END DO

      RETURN
      END

!=======================================================================
!  Utility: compute p, q and deviatoric stress in Abaqus Voigt notation
!=======================================================================
      SUBROUTINE GET_PQ(NDI,NTENS,SIG,P,Q,SDEV)
      INCLUDE 'ABA_PARAM.INC'
      INTEGER NDI,NTENS,I
      DOUBLE PRECISION SIG(NTENS),SDEV(6),P,Q
      DOUBLE PRECISION ZERO,TWO,THIRD,SQRT32
      PARAMETER (ZERO=0.D0,TWO=2.D0,THIRD=1.D0/3.D0)
      PARAMETER (SQRT32=1.224744871391589D0)

      P=ZERO
      DO I=1,NDI
        P=P+THIRD*SIG(I)
      END DO

      DO I=1,NDI
        SDEV(I)=SIG(I)-P
      END DO
      DO I=NDI+1,NTENS
        SDEV(I)=SIG(I)
      END DO

      Q=ZERO
      DO I=1,NDI
        Q=Q+SDEV(I)*SDEV(I)
      END DO
      DO I=NDI+1,NTENS
        Q=Q+TWO*SDEV(I)*SDEV(I)
      END DO
      Q=SQRT32*DSQRT(Q)

      RETURN
      END
