!=======================================================================
!  UMAT: Drucker-Prager Plastizitaet + Norton-Bailey Kriechen
!        (Strain Hardening, von Mises Vergleichsspannung)
!
!  Autor:   generiert via Claude
!  Version: 2.0
!
! ----------------------------------------------------------------------
!  PHYSIKALISCHES MODELL
! ----------------------------------------------------------------------
!
!  Gesamtdehnung: eps = eps_e + eps_p + eps_cr
!
!  [A] ELASTIZITAET
!      sigma = Ce : eps_e
!
!  [B] DRUCKER-PRAGER PLASTIZITAET (nur ausserhalb der Yield-Flaeche)
!      Yield-Funktion:   F  = q - ALPHA*p - K  = 0
!      Plast. Potential: G  = q - PSI*p
!      Flieβregel:       deps_p = dGamma * dG/dsigma
!
!  [C] NORTON-BAILEY KRIECHEN (strain hardening, aktiv in GESAMTER Region)
!      Kriechgesetz (Potenzansatz):
!        d(ecr)/dt = A * q^N * ecr_bar^(-M/(1+M))   [strain hardening]
!
!      Schreibweise als Inkrement:
!        Decr_bar = A * q_eff^N * ecr_bar_old^(-M/(1+M)) * Dt
!
!      Dabei: q_eff = von-Mises-Vergleichsspannung am Inkrementanfang
!             Kriechrichtung = von-Mises-Richtung (deviatorisch)
!
!      Strain-Hardening-Form (ABAQUS-Konvention):
!        ecr_dot = A * q^N * ecr_bar^M       [M <= 0 fuer Verfestigung]
!
! ----------------------------------------------------------------------
!  MATERIALPARAMETER (PROPS)
! ----------------------------------------------------------------------
!    PROPS(1)  = E       [Pa]   Elastizitaetsmodul
!    PROPS(2)  = NU      [-]    Querkontraktionszahl
!    PROPS(3)  = ALPHA   [-]    DP Reibungsparameter
!    PROPS(4)  = K       [Pa]   DP Kohaesion
!    PROPS(5)  = PSI     [-]    DP Dilatanzparameter (0 => assoziiert)
!    PROPS(6)  = A_CR    [1/Pa^N/s]  Norton-Bailey Koeffizient A
!    PROPS(7)  = N_CR    [-]    Norton-Bailey Spannungsexponent n
!    PROPS(8)  = M_CR    [-]    Norton-Bailey Verfestigungsexponent m
!                               (m < 0: Verfestigung, m=0: primaer->sek.)
!
! ----------------------------------------------------------------------
!  ZUSTANDSVARIABLEN (STATEV)
! ----------------------------------------------------------------------
!    STATEV(1)  = EPBAR_P   Aequivalente plastische Dehnung               // Accumulated equivalent plastic strain
!    STATEV(2)  = EPBAR_CR  Aequivalente Kriechdehnung (akkumuliert)      // Accumulated equivalent creep strain
!    STATEV(3)  = Q_EFF     Von-Mises-Vergleichsspannung (Ausgabe)        // von Mises equivalent stress (for output)
!    STATEV(4)  = P_MECH    Hydrostatischer Druck (Ausgabe)               // hydrostatic pressure (for output)
!    STATEV(5)  = DECR_BAR  Kriechdehnungsinkrement (letztes Inkrement)   // creep strain increment from the last increment
!
! ----------------------------------------------------------------------
!  ALGORITHMUS
! ----------------------------------------------------------------------
!  Operator-Split (sequentiell):
!
!  1. Elastischer Versuch (Trial):
!       sigma_TR = sigma_n + Ce : Deps
!
!  2. Kriechkorrektur (implizit, Newton-Raphson):
!       Kriechinkrement Decr_bar bestimmen, so dass:
!       R(Decr_bar) = Decr_bar - A*(q_mid)^N * ecr_bar_mid^M * Dt = 0
!       sigma_cr    = sigma_TR - 2*MU * (3/2) * (s/q) * Decr_bar
!       => mid-point rule fuer q und ecr_bar (verbesserte Stabilitaet)
!
!  3. Plastizitaetspruefung und ggf. Return-Mapping (wie bisher):
!       Auf sigma_cr anwenden
!
!  HINWEIS zur Tangente:
!    Die konsistente algorithmische Tangente wird fuer den
!    gekoppelten Fall (Kriechen + Plastizitaet) als
!    "viskoplastische Sekantentangente" approximiert.
!    Dies ist konservativ aber robust und konvergiert linear.
!    Fuer quadratische Konvergenz waere eine vollgekoppelte
!    algorithmische Tangente noetig (erheblich aufwaendiger).
!=======================================================================

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

!      IMPLICIT NONE
      INCLUDE 'ABA_PARAM.INC'

! --- Uebergabeparameter ---
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

! --- Materialparameter ---
      DOUBLE PRECISION E, NU, ALPHA_DP, K_COH, PSI_DP
      DOUBLE PRECISION A_CR, N_CR, M_CR
      DOUBLE PRECISION ECR0
!     PSI_DP is the dilatency parameter in the DP plastic potential
!     K is the coehsive-related yield stress 


! --- Elastische Moduln ---
      DOUBLE PRECISION LAMBDA, MU, KBULK

! --- Spannungsgroessen ---
      DOUBLE PRECISION STRESS_TR(6)   ! Trial-Spannung
      DOUBLE PRECISION STRESS_CR(6)   ! Spannung nach Kriechkorrektur
      DOUBLE PRECISION S_TR(6)        ! Trial-Deviator
      DOUBLE PRECISION S_CR(6)        ! Deviator nach Kriechkorrektur
      DOUBLE PRECISION P_TR, Q_TR     ! Trial p und q
      DOUBLE PRECISION P_CR, Q_CR     ! p und q nach Kriechkorrektur
      DOUBLE PRECISION S_N(6)         ! Deviator am Inkrementende (DP)
      DOUBLE PRECISION P_N, Q_N       ! p,q am Inkrementende (DP)

! --- Plastizitaetsgroessen ---
      DOUBLE PRECISION EPBAR_P, DEPBAR_P, DGAMMA, F_TR_CR, BETA_PL
!     EPBAR_P  is the accumulated equivalent plastic strain
!     DEPBAR_P is the increment of EPBAR_P in the step
!     DGAMMA   is Delta gamma, the plastic multiplier (from DP flow rule)
!     F_TR_CR  is the yield function evaluated at the creep-corrected trial stress (from DP yield function)
!     BETA_PL  is a scaling factor

! --- Kriechgroessen ---
      DOUBLE PRECISION EPBAR_CR, DECR_BAR
      DOUBLE PRECISION ECR_EFF_MID
!     EPBAR_CR  is the accumulated equivalent creep strain
!     DECR_BAR is the increment of EPBAR_CR in the step

      DOUBLE PRECISION Q_MID, ECR_MID           ! Midpoint-Groessen
!     Q_MID is the VM stress at mid point, ECR_MID is the creep strain at mid point -> for improving the Newton-Raphson stability

      DOUBLE PRECISION RES, DRES, DECR_IT        ! Newton-Raphson
!     RES si the Newton-Raphson residual, DRES is the derivative of RES, DECR_IT is the correction, i.e. the update per iteration

      DOUBLE PRECISION ECRDOT_OLD                ! Kriechrate alt
!     ECRDOT_OLD is the creep rate from the previous increment

      DOUBLE PRECISION CR_FACT                   ! skalarer Kriechfaktor
!     A scalar creep factor

! --- Hilfsvariablen ---
      DOUBLE PRECISION DELTA(6)       ! Kronecker-Delta (Voigt)
      DOUBLE PRECISION DDSDDE_E(6,6)  ! Elastische Steifigkeit
      DOUBLE PRECISION THETA_PL, A1   ! Tangenten-Parameter (DP)
!     Tangeant parameter for DP, related to the return mapping scalling

      DOUBLE PRECISION ECR_MIN        ! Schutz gegen ecr=0
      DOUBLE PRECISION NORM_Q, NORM_S
!     Norm of the von Mises stress and devatoric stress, respectively

! --- Numerische Konstanten ---
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,HALF,ONETHIRD,TWO3
      DOUBLE PRECISION SQRT32, TOL_CR, TOL_PL, TOL_NR
      INTEGER I, J, ITER, MAXITER_CR
      LOGICAL PLASTIC, CR_ACTIVE

      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0)
      PARAMETER (HALF=0.5D0, ONETHIRD=1.D0/3.D0, TWO3=2.D0/3.D0)
      PARAMETER (SQRT32=1.224744871391589D0)   ! sqrt(3/2)
      PARAMETER (TOL_CR=1.D-10)   ! Toleranz Kriechiteration
      PARAMETER (TOL_PL=1.D-10)   ! Toleranz Yield-Check
      PARAMETER (TOL_NR=1.D-12)   ! Newton-Raphson Toleranz
      PARAMETER (MAXITER_CR=100)  ! Max. Kriechiterationen
      PARAMETER (ECR_MIN=1.D-15)  ! Mindestwert ecr_bar (Num. Schutz)

!=======================================================================
!  1. MATERIALPARAMETER
!=======================================================================
      E        = PROPS(1)
      NU       = PROPS(2)
      ALPHA_DP = PROPS(3)
      K_COH    = PROPS(4)
      PSI_DP   = PROPS(5)
      A_CR     = PROPS(6)
      N_CR     = PROPS(7)
      M_CR     = PROPS(8)   ! Erwartet M_CR <= 0 fuer Verfestigung
      ECR0     = PROPS(9)
      IF (ECR0 .LT. ECR_MIN) ECR0 = ECR_MIN
!     Assoziierte Flieβregel falls PSI=0
      IF (DABS(PSI_DP) .LT. TOL_PL) PSI_DP = ALPHA_DP

!     Lame-Konstanten
      LAMBDA = E*NU / ((ONE+NU)*(ONE-TWO*NU))
      MU     = E / (TWO*(ONE+NU))
      KBULK  = E / (THREE*(ONE-TWO*NU))

!     Kronecker-Delta (Voigt-Notation)
      DO I = 1, NTENS
        DELTA(I) = ZERO
      END DO
      DO I = 1, NDI
        DELTA(I) = ONE
      END DO

!=======================================================================
!  2. ZUSTANDSVARIABLEN LESEN
!=======================================================================
      EPBAR_P  = STATEV(1)
      EPBAR_CR = STATEV(2)

!     Schutz: Kriechdehnung darf nicht negativ oder null sein
!     change from EPBAR_CR = ECR_MIN to EPBAR_CR = 0
      IF (EPBAR_CR .LT. ZERO) EPBAR_CR = ZERO 


!=======================================================================
!  3. ELASTISCHE STEIFIGKEITSMATRIX
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
!  4. ELASTISCHER TRIAL-STRESS
!     sigma_TR = sigma_n + Ce : Deps
!     with STRESS the stress at the beg. of the increment,
!     DSTRAN = Deps
!=======================================================================
      DO I = 1, NTENS
        STRESS_TR(I) = STRESS(I)
        DO J = 1, NTENS
          STRESS_TR(I) = STRESS_TR(I) + DDSDDE_E(I,J)*DSTRAN(J)
        END DO
      END DO

!     Trial: p (Zugspannung positiv) und Deviator s
      P_TR = ZERO
      DO I = 1, NDI
        P_TR = P_TR + ONETHIRD * STRESS_TR(I)
      END DO
      DO I = 1, NDI
        S_TR(I) = STRESS_TR(I) - P_TR
      END DO
      DO I = NDI+1, NTENS
        S_TR(I) = STRESS_TR(I)
      END DO

!     Trial: von-Mises-Vergleichsspannung
!     q = sqrt(3/2 * s:s),  Voigt: s:s = sum_NDI(si^2) + 2*sum_NSHR(si^2)
      Q_TR = ZERO
      DO I = 1, NDI
        Q_TR = Q_TR + S_TR(I)*S_TR(I)
      END DO
      DO I = NDI+1, NTENS
        Q_TR = Q_TR + TWO*S_TR(I)*S_TR(I)
      END DO
      Q_TR = SQRT32 * DSQRT(Q_TR)

!=======================================================================
!  5. NORTON-BAILEY KRIECHKORREKTUR (implizit, mid-point)
!
!  Norton-Bailey (strain hardening):
!    ecr_dot = A * q^n * ecr_bar^m
!
!  Inkrementelle Form (Backward-Euler / mid-point):
!    Decr_bar = A * q_mid^n * ecr_bar_mid^m * Dt
!
!  Midpoint-Approximation fuer numerische Stabilitaet:
!    q_mid       = q_TR   - 3*MU * Decr_bar        (linear in Decr_bar)
!    ecr_bar_mid = ecr_bar_old + HALF * Decr_bar    (mean value)
!
!  Residuum:
!    R(Decr_bar) = Decr_bar - A * q_mid^n * ecr_bar_mid^m * Dt = 0
!
!  Loesung mit Newton-Raphson.
!=======================================================================

!     Nur Kriechen wenn DTIME > 0 und Kriechparameter gesetzt
!     CR_Active is TRUE if time increment is >0, A_CR>0, deviatoric creep >tol, i.e. without there is nothing to drive creep
      CR_ACTIVE = (DTIME .GT. ZERO) .AND.
     1            (A_CR .GT. ZERO) .AND.
     2            (Q_TR .GT. TOL_CR)

      DECR_BAR = ZERO

      IF (CR_ACTIVE) THEN

!       Startwert: explizites Euler-Inkrement (Vorschaetzung)
        Q_MID   = Q_TR
        ECR_MID = EPBAR_CR
        ECR_EFF_MID = ECR_MID + ECR0

!       Make sure the creep strain is strictely positive
        IF (ECR_EFF_MID .LT. ECR_MIN) ECR_EFF_MID  = ECR_MIN

!       Define the increment of accumulated equivalent creep strain in the step
        DECR_BAR = A_CR * (Q_MID**N_CR) * (ECR_EFF_MID **M_CR) * DTIME

!       Schutz: Decr_bar nicht groesser als q/(3*MU) (verhindert q<0)
!       Safety: prevent the vM stress to be negative by reducing the strain -> prevent overshooting 
        IF (DECR_BAR .GT. Q_TR/(THREE*MU+TOL_CR)) THEN
          DECR_BAR = Q_TR / (THREE*MU+TOL_CR) * 0.9D0
        END IF

!       Newton-Raphson fuer implizites Kriechinkrement
        DO ITER = 1, MAXITER_CR

!         Midpoint-Groessen
          Q_MID   = Q_TR - THREE*MU*DECR_BAR
          ECR_MID = EPBAR_CR + HALF*DECR_BAR
          ECR_EFF_MID = ECR_MID + ECR0

!         Schutz gegen negative Werte
          IF (Q_MID   .LT. TOL_CR)  Q_MID   = TOL_CR
          IF (ECR_EFF_MID  .LT. ECR_MIN) ECR_EFF_MID  = ECR_MIN

!         Kriechrate am Midpoint
          ECRDOT_OLD = A_CR * (Q_MID**N_CR) * (ECR_EFF_MID **M_CR)

!         Residuum: R = Decr_bar - ecr_dot * Dt
          RES = DECR_BAR - ECRDOT_OLD * DTIME

          IF (DABS(RES) .LT. TOL_NR*(ONE + DECR_BAR)) EXIT

!         Ableitung dR/d(Decr_bar) fuer Newton-Schritt:
!         dR/d(Decr_bar) = 1
!                        - Dt * d(ecr_dot)/d(Decr_bar)
!         d(ecr_dot)/d(Decr_bar) =
!           A*(Q_MID^(N-1))*(-3*MU)*N*ECR_MID^M
!         + A*(Q_MID^N)*M*ECR_MID^(M-1)*(HALF)

!         Compute the derivative of the residual with respect to DECR_BAR -> Jacobian of the NR algorithm
          DRES = ONE - DTIME * (
     1      A_CR * N_CR * (Q_MID**(N_CR-ONE)) * (-THREE*MU)
     2                 * (ECR_EFF_MID**M_CR)
     3    + A_CR * (Q_MID**N_CR)
     4           * M_CR * (ECR_EFF_MID**(M_CR-ONE)) * HALF )


!         Numerischer Schutz: Jacobian nicht zu klein
          IF (DABS(DRES) .LT. TOL_CR) DRES = ONE


!         Newton increment af current NR step 
          DECR_IT = -RES / DRES
!         Daempfung: maximal 50% Aenderung pro Iteration
!         if the update is more than half of DECR_BAR, fix the update to half DECR_BAR
          IF (DECR_IT .GT. HALF*DECR_BAR+TOL_CR)
     1      DECR_IT = HALF*(DECR_BAR+TOL_CR)

!         if the update is more than half of DECR_BAR, fix the update to half DECR_BAR (in the negative side)
          IF (DECR_IT .LT. -HALF*DECR_BAR)
     1      DECR_IT = -HALF*DECR_BAR

!         after damping, actualize DECR_BAR
          DECR_BAR = DECR_BAR + DECR_IT

!         Physikalische Schranke: Decr_bar >= 0 und q nach Korrektur > 0
          IF (DECR_BAR .LT. ZERO) DECR_BAR = ZERO
          IF (DECR_BAR .GT. Q_TR/(THREE*MU+TOL_CR))
     1      DECR_BAR = Q_TR/(THREE*MU+TOL_CR)*0.99D0

        END DO ! Newton-Raphson

!       In case of non convergence at the end of the MAXITER_CR iteration, restart the ABQ increment, but with a smaller time step (1/4 smaller)
        IF (DABS(RES) .GT. TOL_NR*(ONE + DECR_BAR)) THEN
          PNEWDT = 0.25D0
          RETURN
        END IF

      END IF ! CR_ACTIVE

!=======================================================================
!  6. SPANNUNGSUPDATE NACH KRIECHKORREKTUR
!
!  Kriechrichtung = Richtung des deviatorischen Spannungstensors
!  (von Mises assoziiert):
!    n_cr = (3/2) * s / q
!
!  Spannungskorrektur:
!    sigma_CR = sigma_TR - 2*MU * n_cr * Decr_bar
!             = sigma_TR - 3*MU * (s_TR/q_TR) * Decr_bar
!
!  Hydrostatischer Anteil unveraendert (vol. Kriechen = 0 fuer vMises)
!=======================================================================

!     Creep stress     
      DO I = 1, NTENS
        STRESS_CR(I) = STRESS_TR(I)
      END DO

      IF (DECR_BAR .GT. ZERO .AND. Q_TR .GT. TOL_CR) THEN
        CR_FACT = THREE*MU*DECR_BAR / Q_TR
        DO I = 1, NDI
          STRESS_CR(I) = STRESS_TR(I) - CR_FACT * S_TR(I)
        END DO
        DO I = NDI+1, NTENS
          STRESS_CR(I) = STRESS_TR(I) - CR_FACT * S_TR(I)
        END DO
      END IF

!     Hydrostatic creep stress
      P_CR = ZERO
      DO I = 1, NDI
        P_CR = P_CR + ONETHIRD * STRESS_CR(I)
      END DO

!     Deviatoric creep stress
      DO I = 1, NDI
        S_CR(I) = STRESS_CR(I) - P_CR
      END DO
      DO I = NDI+1, NTENS
        S_CR(I) = STRESS_CR(I)
      END DO

!     vM creep stress
      Q_CR = ZERO
      DO I = 1, NDI
        Q_CR = Q_CR + S_CR(I)*S_CR(I)
      END DO
      DO I = NDI+1, NTENS
        Q_CR = Q_CR + TWO*S_CR(I)*S_CR(I)
      END DO
      Q_CR = SQRT32 * DSQRT(Q_CR)

!=======================================================================
!  7. DRUCKER-PRAGER YIELD-CHECK (auf kriechkorrigierte Spannung)
!
!  F = q_CR - ALPHA*p_CR - K
!=======================================================================

!     Definition of the yield function with the creep corrected stress
      F_TR_CR = Q_CR + ALPHA_DP*P_CR - K_COH

      IF (F_TR_CR .LE. TOL_PL) THEN
        PLASTIC = .FALSE.
      ELSE
        PLASTIC = .TRUE.
      END IF

!=======================================================================
!  8. PLASTISCHE KORREKTUR (Return-Mapping, analytisch)
!
!  Analytische Loesung fuer linearen Drucker-Prager:
!    DGAMMA = F_TR_CR / (3*MU + ALPHA*PSI*KBULK)
!=======================================================================
      IF (.NOT. PLASTIC) THEN
!       Keine Plastizitaet: kriechkorrigierte Spannung ist Endzustand
        DO I = 1, NTENS
          STRESS(I) = STRESS_CR(I)
        END DO
        DEPBAR_P = ZERO
        DGAMMA   = ZERO

      ELSE
!       Plastisches Return-Mapping auf sigma_CR

!       if vM creep stress (equivalently the deviatorix stress) is near zero --> the stress is almost only hydrostatic
!       Then the stress update is purely hydrostatic and the shear components are killed
        IF (Q_CR .LT. TOL_PL) THEN
!         Rein hydrostatischer Fall nach Kriechkorrektur
          DGAMMA = -(P_CR + K_COH/ALPHA_DP) / (PSI_DP * KBULK)
          IF (DGAMMA .LT. ZERO) DGAMMA = ZERO
          DO I = 1, NDI
            STRESS(I) = STRESS_CR(I) - KBULK*PSI_DP*DGAMMA
          END DO
          DO I = NDI+1, NTENS
            STRESS(I) = ZERO
          END DO



!       Normal case for the return mapping 
        ELSE
!         Define DGAMMA by assuming that return maping requires solving for DGAMMA such that F (the yield function) = 0 after correction
          DGAMMA  = F_TR_CR / (THREE*MU + ALPHA_DP*PSI_DP*KBULK)
          
!         Deviatoric scalling 
          BETA_PL = ONE - THREE*MU*DGAMMA/Q_CR

!         Pressure correction
          P_N     = P_CR - KBULK*PSI_DP*DGAMMA

          DO I = 1, NDI
            S_N(I)  = BETA_PL * S_CR(I)
            STRESS(I) = S_N(I) + P_N
          END DO
          DO I = NDI+1, NTENS
            S_N(I)  = BETA_PL * S_CR(I)
            STRESS(I) = S_N(I)
          END DO
        END IF

        DEPBAR_P = DGAMMA  ! aequiv. plast. Dehnung (DP-Konvention)

      END IF

!=======================================================================
!  9. ZUSTANDSVARIABLEN UPDATEN
!=======================================================================
      STATEV(1) = EPBAR_P  + DEPBAR_P
      STATEV(2) = EPBAR_CR + DECR_BAR
      STATEV(3) = Q_CR                  ! von-Mises q (Ausgabe)
      STATEV(4) = P_CR                  ! hydrostatischer Druck (Ausgabe)
      STATEV(5) = DECR_BAR              ! Kriechinkrement (Ausgabe)

!=======================================================================
!  10. TANGENTENSTEIFIGKEIT (algorithmisch)
!
!  Strategie: "viskoplastische Sekantentangente"
!
!  Der Gesamtoperator ist: Ce_cr + Ce_pl
!
!  [A] Rein kriechend (kein Plastizitaet):
!      Ceff = Ce * (1 - 3*MU*Decr_bar/Q_TR * ...)  (reduzierter Schub)
!
!  [B] Kriechend + plastisch:
!      Ceff = Ce_cr weiter reduziert durch Plastizitaet
!
!  Fuer beide Faelle wird eine konsistente Sekantentangente gebaut:
!
!    DDSDDE = KBULK * (delta x delta)
!           + 2*MU * theta_cr * IIdev
!           - [Rang-1 Term aus DP-Return-Mapping, falls plastisch]
!
!  mit theta_cr = 1 - 3*MU*Decr_bar/Q_TR  (Kriech-Abminderung)
!=======================================================================

!     Kriech-Reduktionsfaktor fuer die Schertangente
      IF (Q_TR .GT. TOL_CR .AND. DECR_BAR .GT. ZERO) THEN
        CR_FACT = ONE - THREE*MU*DECR_BAR/Q_TR
        IF (CR_FACT .LT. ZERO) CR_FACT = ZERO
      ELSE
        CR_FACT = ONE
      END IF

!     Initialisierung
      DO I = 1, NTENS
        DO J = 1, NTENS
          DDSDDE(I,J) = ZERO
        END DO
      END DO

!     Volumenanteil: KBULK * delta_i * delta_j
      DO I = 1, NDI
        DO J = 1, NDI
          DDSDDE(I,J) = KBULK
        END DO
      END DO

!     Deviatorischer Anteil (kriechreduziert)
      DO I = 1, NDI
        DDSDDE(I,I) = DDSDDE(I,I) + TWO*MU*CR_FACT
        DO J = 1, NDI
          DDSDDE(I,J) = DDSDDE(I,J) - TWO*MU*CR_FACT*ONETHIRD
        END DO
      END DO
      DO I = NDI+1, NTENS
        DDSDDE(I,I) = MU*CR_FACT
      END DO

!     Plastische Rang-1-Korrektur (falls Plastizitaet aktiv)
      IF (PLASTIC .AND. Q_CR .GT. TOL_PL .AND. DGAMMA .GT. ZERO) THEN

        THETA_PL = ONE / (THREE*MU + ALPHA_DP*PSI_DP*KBULK)

!       Rang-1: -(3MU)^2 * THETA * (s_CR x s_CR) / q_CR^2
        DO I = 1, NTENS
          DO J = 1, NTENS
            DDSDDE(I,J) = DDSDDE(I,J)
     1        - (THREE*MU)**2 * THETA_PL
     2        * S_CR(I)*S_CR(J) * TWO3 / Q_CR**2
          END DO
        END DO

!       Kopplung Volumen-Deviator (Reibung/Dilatanz):
        DO I = 1, NDI
          DO J = 1, NTENS
            DDSDDE(I,J) = DDSDDE(I,J)
     1        - KBULK*PSI_DP * THREE*MU*THETA_PL
     2        * S_CR(J) * TWO3 / Q_CR
          END DO
          DO J = 1, NDI
            DDSDDE(I,J) = DDSDDE(I,J)
     1        - THREE*MU*THETA_PL * KBULK*ALPHA_DP
     2        * S_CR(I) * TWO3 / Q_CR
          END DO
        END DO

!       Zusaetzlicher deviatorischer Term (Plastizitaets-Abminderung)
        Q_N = Q_CR - THREE*MU*DGAMMA
        IF (Q_N .LT. TOL_PL) Q_N = TOL_PL
        DO I = 1, NDI
          DDSDDE(I,I) = DDSDDE(I,I)
     1      + TWO*MU*(ONE - THREE*MU*DGAMMA/Q_CR - CR_FACT)
          DO J = 1, NDI
            DDSDDE(I,J) = DDSDDE(I,J)
     1        - TWO*MU*(ONE - THREE*MU*DGAMMA/Q_CR - CR_FACT)*ONETHIRD
          END DO
        END DO
        DO I = NDI+1, NTENS
          DDSDDE(I,I) = DDSDDE(I,I)
     1      + MU*(ONE - THREE*MU*DGAMMA/Q_CR - CR_FACT)
        END DO

      END IF

      RETURN
      END SUBROUTINE UMAT

!=======================================================================
!  ANMERKUNGEN ZUR KALIBRIERUNG
! -----------------------------------------------------------------------
!
!  Norton-Bailey Strain-Hardening (ABAQUS-Konvention):
!    ecr_dot = A * q^n * ecr_bar^m
!
!  Typische Werte (Gestein/Boden):
!    A: 1e-20 bis 1e-10  (stark materialabhaengig)
!    n: 1 bis 5          (Spannungsexponent)
!    m: -0.9 bis -0.1    (negativ = Verfestigung)
!         m = 0  => primaeres Kriechen (strain rate ~ const)
!         m < 0  => strain hardening (Kriechrate nimmt ab)
!
!  ABAQUS-Inputdeck Beispiel:
!  *Material, name=DP_NortonBailey
!  *User Material, constants=8
!   210000., 0.3, 0.3, 100., 0.0, 1.0e-16, 3.0, -0.5
!  *Depvar
!   5
!  *Initial Conditions, type=SOLUTION
!   <instance>, 0., 1.e-15, 0., 0., 0.
!
!  ZUSTANDSVARIABLEN-Ausgabe (Field Output: SDV):
!   SDV1 = aequiv. plast. Dehnung
!   SDV2 = aequiv. Kriechdehnung
!   SDV3 = von-Mises q
!   SDV4 = hydrostatischer Druck p
!   SDV5 = Kriechdehnungsinkrement (letztes Dt)
!=======================================================================
