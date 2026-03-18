C ================= Same as in ORIENT_SPOTS & HKL_GEN =========================

C ================= Generate Main HKLs ========================================
C	SUBROUTINE GEN_HEMISPHERE(HKLS,NHKLS, UB,ICEN,D_MIN,D_MAX)
C ================= Check for HKL special conditions ==========================
C	SUBROUTINE CHECK_SPECIALS(LALLOW, IH,IK,IL)
C ================= Add satellite/modulation HKLs =============================
C	SUBROUTINE ADD_MOD_SAT_HKLS(NHKLS, UB,D_MIN,D_MAX)
C	SUBROUTINE ADD_MOD_SAT_HKL(IHKL_MOD, IMULT,IHKL,IMOD)
C ================= Prune HKL list to X,Y, and wavelength limits ==============
C	SUBROUTINE PRUNE_XY_WAV(NHKLS, WAV_MIN,WAV_MAX,X_MIN,X_MAX,Y_MIN,Y_MAX)
C ================= Prune HKL list by merge/remove harmonic multiples =========
C	SUBROUTINE PRUNE_MULTIPLES(NHKLS)
C	SUBROUTINE PRUNE_XY_WAV(NHKLS, WAV_MIN,WAV_MAX,X_MIN,X_MAX,Y_MIN,Y_MAX)
C =============================================================================

	SUBROUTINE GEN_HEMISPHERE(HKLS,NHKLS, UB,ICEN,D_MIN,D_MAX)
C
C Calculate a list of HKLs from UB matrix within d-spacing limits, and
C obeying the centering given by ICEN = 0(P) 1(A) 2(B) 3(C) 4(I) 5(F) 6(R).
C
C HKLs are also checked against special conditions in /SPECIALS_COM/
C The list is half of all reflections as it excludes the Friedel pair
C NB: Since we ignore X,Y,wav limits we don't need to know the PHI value
C
	PARAMETER (NHKLS_MAX=300000)
	REAL HKLS(3,NHKLS_MAX),UB(3,3)
C
	LOGICAL LALLOW
	INTEGER H,H_MAX
C
C The max & min values of |UB.HKL|^2
C
	SSQ_MAX=1/D_MIN**2
	SSQ_MIN=1/D_MAX**2
C
C Dot products needed to calculate the limits in H,K,L
C
	DOT11=UB(1,1)*UB(1,1)+UB(2,1)*UB(2,1)+UB(3,1)*UB(3,1)
	DOT12=UB(1,1)*UB(1,2)+UB(2,1)*UB(2,2)+UB(3,1)*UB(3,2)
	DOT22=UB(1,2)*UB(1,2)+UB(2,2)*UB(2,2)+UB(3,2)*UB(3,2)
	DOT13=UB(1,1)*UB(1,3)+UB(2,1)*UB(2,3)+UB(3,1)*UB(3,3)
	DOT23=UB(1,2)*UB(1,3)+UB(2,2)*UB(2,3)+UB(3,2)*UB(3,3)
	DOT33=UB(1,3)*UB(1,3)+UB(2,3)*UB(2,3)+UB(3,3)*UB(3,3)
C
C Calculate the limits in H to satisfy SSQ_MAX, for any K or L
C
	TEMP=(DOT11*DOT33-DOT13**2) - (DOT33*DOT12-DOT13*DOT23)**2 / (DOT22*DOT33-DOT23**2)
	H_MAX=FLOOR(SQRT( SSQ_MAX*DOT33/TEMP ))
C
C Zero the list of hkl's
C
	NHKLS=0
C
C Loop over H
C
	DO H=0,H_MAX
C Calculate the limits in K to satisfy SSQ_MAX, for any L
	  TEMP1=(DOT22-DOT23**2/DOT33)*SSQ_MAX
	  TEMP2=(DOT11-DOT13**2/DOT33)*(DOT22-DOT23**2/DOT33) - (DOT12-DOT13*DOT23/DOT33)**2
	  TEMP=TEMP1-TEMP2*H**2
	  IF(TEMP .LT. 0.0) CALL QUIT('BUG: TEMP < 0')
	  RK1=( (DOT13*DOT23/DOT33-DOT12)*H - SQRT(TEMP) )/(DOT22-DOT23**2/DOT33)
	  RK2=( (DOT13*DOT23/DOT33-DOT12)*H + SQRT(TEMP) )/(DOT22-DOT23**2/DOT33)
	  K_LO=CEILING( MIN(RK1,RK2) )
	  K_HI=FLOOR  ( MAX(RK1,RK2) )
	  K_INC=1
C If H=0, make K_LO > -1 so we don't create Friedel pairs
	  IF(H .EQ. 0) K_LO=MAX(0,K_LO)
C For C & F centering (H+K=2N), we advance K in steps of 2 starting
C from H+K is even
	  IF(ICEN.EQ.3 .OR. ICEN.EQ.5) THEN
	    K_INC=2
	    IF((H+K_LO)/2*2 .NE. H+K_LO) K_LO=K_LO+1
	  ENDIF
C
C Loop over K
C
	  DO K=K_LO,K_HI,K_INC
C Calculate the limits in L that satisfy SSQ_MAX
	    TEMP=DOT33*(SSQ_MAX-DOT11*H**2-2*DOT12*H*K-DOT22*K**2) + (DOT13*H+DOT23*K)**2
	    IF(TEMP .LT. 0.0) CALL QUIT('BUG: TEMP(L) < 0')
	    RL1=( -(DOT13*H+DOT23*K) - SQRT(TEMP) )/DOT33
	    RL2=( -(DOT13*H+DOT23*K) + SQRT(TEMP) )/DOT33
	    L_LO=CEILING( MIN(RL1,RL2) )
	    L_HI=FLOOR  ( MAX(RL1,RL2) )
	    L_INC=1
C If H&K=0, make L_LO > 0 so we don't create Friedel pairs or 0,0,0
	  IF(H.EQ.0 .AND. K.EQ.0) L_LO=MAX(1,L_LO)
C For A & F centering (K+L=2N), we advance L in steps of 2 starting
C from K+L is even
	    IF(ICEN.EQ.1 .OR. ICEN.EQ.5) THEN
	      L_INC=2
	      IF((K+L_LO)/2*2 .NE. K+L_LO) L_LO=L_LO+1
	    ENDIF
C For B centering (H+L=2N), we advance L in steps of 2 starting
C from H+L is even
	    IF(ICEN .EQ. 2) THEN
	      L_INC=2
	      IF((H+L_LO)/2*2 .NE. H+L_LO) L_LO=L_LO+1
	    ENDIF
C For I centering (H+K+L=2N), we advance L in steps of 2 starting
C from H+K+L is even
	    IF(ICEN .EQ. 4) THEN
	      L_INC=2
	      IF((H+K+L_LO)/2*2 .NE. H+K+L_LO) L_LO=L_LO+1
	    ENDIF
C If R centering (rhom in hex setting) (-H+K+L=3N), we advance L in steps of 3
C starting from -H+K+L is a multiple of three
	    IF(ICEN .EQ. 6) THEN
	      L_INC=3
	      IF((-H+K+L_LO)/3*3 .NE. -H+K+L_LO) L_LO=L_LO+1
	      IF((-H+K+L_LO)/3*3 .NE. -H+K+L_LO) L_LO=L_LO+1
	    ENDIF
C
C Loop over L, and add suitable H,K,L to HKLS()
C
	    DO L=L_LO,L_HI,L_INC
C Ignore H,K,L with |UB.HKL|^2 < SSQ_MIN
	      SSQ= DOT11*H*H + DOT12*H*K + DOT13*H*L + 
	1	   DOT12*K*H + DOT22*K*K + DOT23*K*L + 
	2	   DOT13*L*H + DOT23*L*K + DOT33*L*L
	      IF(SSQ .LT. SSQ_MIN) CYCLE
C Ignore H,K,L if not allowed by special conditions
	      CALL CHECK_SPECIALS(LALLOW, H,K,L)
	      IF( .NOT.LALLOW ) CYCLE
C Add H,K,L to end of list
	      NHKLS=NHKLS+1
	      IF(NHKLS .GT. NHKLS_MAX) CALL QUIT('ERROR: Max. number of HKLs exceeded')
	      HKLS(1,NHKLS)=H
	      HKLS(2,NHKLS)=K
	      HKLS(3,NHKLS)=L
C
C End of loops over H, K & L
C
	    ENDDO
	  ENDDO
	ENDDO
C
	RETURN
	END


C ================= Check for HKL special conditions =============

	SUBROUTINE CHECK_SPECIALS(LALLOW, IH,IK,IL)
C
	LOGICAL LALLOW
C
	COMMON /SPECIALS_COM/ SPECS(3,100),NSPECS
C
C SPECS() contains 3-digit codes for linear and planar special conditions
C   0 means HKL must also be zero to match this condition
C   1,2,3,4,6 means divide HKL by this number and sum up values
C   If HKL matches all zeros, a non-integer sum means the HKL is forbidden
C
C 3-digit codes for linear/planar special conditions where h=0:
C   00l,l=2n  00l,l=3n  00l,l=4n  00l,l=6n  0kl,k=2n  0kl,l=2n  0kl,k+l=2n  0kl,k+l=4n
C     002       003	  004	    006	      021       012       022         044
C
C Quick checks for allowed reflections
C	No special conditions, or all HKL are non-zero
C
	LALLOW=.TRUE.
	IF(NSPECS .EQ. 0) RETURN
	IF(MIN( ABS(IH), ABS(IK), ABS(IL) ) .NE. 0) RETURN
C
C Check that HKL obeys all matching special conditions
C
	LALLOW=.FALSE.
	DO I=1,NSPECS
C If zero indices don't match, skip to next special condition
	  IF( (SPECS(1,I) .EQ. 0.0) .AND. (IH .NE. 0) ) CYCLE
	  IF( (SPECS(2,I) .EQ. 0.0) .AND. (IK .NE. 0) ) CYCLE
	  IF( (SPECS(3,I) .EQ. 0.0) .AND. (IL .NE. 0) ) CYCLE
C Create sum of HKL/SPECS(), with trick to avoid 0/0
	  SUM=IH/MAX(1E-6,SPECS(1,I)) + IK/MAX(1E-6,SPECS(2,I)) + IL/MAX(1E-6,SPECS(3,I))
C If sum is non-integer, return with HKL not allowed
	  IF(ABS( SUM - NINT(SUM) ) .GT. 1E-3) RETURN
C
	ENDDO
C
	LALLOW=.TRUE.
	RETURN
	END


C ================= Add satellite HKLs =============

	SUBROUTINE ADD_SAT_HKLS(HKLS,IMODS,NHKLS, UB,D_MIN,D_MAX)
C
C Adds satellite spots to main spots (or their Friedels) with d-spacing > DMIN_MOD
C Only keep satellites with positive wavelengths and d-spacing of D_MIN to D_MAX
C
	INTEGER IMODS(*)
	REAL HKLS(3,*),UB(3,3)
C
	COMMON /MODULATE_COM/ NHKLS_MOD,DMIN_MOD,IMOD_STRONG,RHKLS_MOD(3,100)
C
	PARAMETER (NHKLS_MAX=300000)
	LOGICAL LVALID
	REAL HVECS(3,NHKLS_MAX),HVECS_MOD(3,100),HVEC_SAT(3)
C
C Calculate reciprocal lattive vectors (HVEC) for the main hkls
C
	CALL CALC_HKLS_TO_HVECS(HVECS, UB,HKLS,NHKLS)
C
C Calculate reciprocal lattive vectors (HVEC) for the fraction hkls of the modulation
C
	CALL CALC_HKLS_TO_HVECS(HVECS_MOD, UB,RHKLS_MOD,NHKLS_MOD)
C
C Try making satellites for all main peaks and add acceptable ones to the
C end of the reflection list
C
	NHKLS_MAIN=NHKLS
	DO IHKL_MAIN=1,NHKLS_MAIN
	  DSPACE=1.0/SQRT(HVECS(IHKL_MAIN,1)**2 + HVECS(IHKL_MAIN,2)**2 + HVECS(IHKL_MAIN,3)**2 )
C Ignore main HKLs if d-spacing less than "modulation" limit
	  IF(DSPACE .LT. DMIN_MOD) CYCLE
C
C Check then add valid satellite HKLs created from the main HKLs and their Friedel pairs
	  DO IMOD=1,NHKLS_MOD
	    DO MULT=-1,1,2
	      HVEC_SAT(1)=MULT*HVECS(1,IHKL_MAIN) + HVECS_MOD(1,IMOD)
	      HVEC_SAT(2)=MULT*HVECS(2,IHKL_MAIN) + HVECS_MOD(2,IMOD)
	      HVEC_SAT(3)=MULT*HVECS(3,IHKL_MAIN) + HVECS_MOD(3,IMOD)
	      CALL CHECK_SAT_HKL(LVALID, HVEC_SAT, D_MIN,D_MAX)
C If valid, add satellite to the end of the HKL list
	      IF( LVALID ) THEN
	        NHKLS=NHKLS+1
	        IF(NHKLS .GT. NHKLS_MAX) CALL QUIT('ERROR: Too many satellites')
C Store main hkl (or Friedel) and modulation number
	        HKLS(1,NHKLS)=MULT*HKLS(1,IHKL_MAIN)
	        HKLS(2,NHKLS)=MULT*HKLS(2,IHKL_MAIN)
	        HKLS(3,NHKLS)=MULT*HKLS(3,IHKL_MAIN)
	        IMODS( NHKLS)=IMOD
	      ENDIF
	    ENDDO
	  ENDDO
C
C Iterate on IHKL_MAIN
	ENDDO
C
C We are finished if "modulation" d-spacing limit > 0
C
	IF(DMIN_MOD .GT. 0.0) RETURN
C
C Add modulation for main HKL = (0,0,0)
C
	DO IMOD=1,NHKLS_MOD
	  CALL CHECK_SAT_HKL(LVALID, HVECS_MOD(1,IMOD), D_MIN,D_MAX)
C If valid, add satellite to the end of the HKL list
	  IF( LVALID ) THEN
	    NHKLS=NHKLS+1
	    IF(NHKLS .GT. NHKLS_MAX) CALL QUIT('ERROR: Too many satellites')
C Store modulation number with main hkl = (0,0,0)
	    HKLS(1,NHKLS)=0.0
	    HKLS(2,NHKLS)=0.0
	    HKLS(3,NHKLS)=0.0
	    IMODS( NHKLS)=IMOD
	  ENDIF
	ENDDO
C
	RETURN
	END


	SUBROUTINE CHECK_SAT_HKL(LVALID, HVEC, D_MIN,D_MAX)
C
	LOGICAL LVALID
	REAL HVEC(3)
C
C Set status to failed
C
	LVALID=.FALSE.
C
C Skip if d-spacing outside general limits
C
	DSPACE=1.0/MAX(1E-6, SQRT(HVEC(1)**2 + HVEC(2)**2 + HVEC(3)**2) )
	IF(DSPACE.LT.D_MIN .OR. DSPACE.GT.D_MAX) RETURN
C
C Skip if negative wavelength
C
	CALL CALC_HVEC_TO_SVEC(S1X,S1Y,S1Z,WAV, HVEC(1),HVEC(2),HVEC(3))
	IF(WAV .LT. 0.0) RETURN
C
C Set status to success, and return
C
	LVALID=.TRUE.
	RETURN
	END

C ================= Prune HKL list to X, Y, and wavelength limits =============

	SUBROUTINE PRUNE_XY_WAV(HKLS,IHKLS,IMODS,NHKLS, UB, X_MIN,X_MAX,Y_MIN,Y_MAX, WAV_MIN,WAV_MAX)
C
C Prune HKLs with pixel X,Y or wavelength not between limits
C However, first converts any HKLS with negative wavelength to -H,-K,-L
C
	INTEGER IHKLS(3,*),IMODS(*)
	REAL HKLS(3,*),UB(3,3)
C
	PARAMETER (NHKLS_MAX=300000)
	REAL PIXWAVS(3,NHKLS_MAX)
C
	CALL CALC_HKLS_TO_PIXWAVS(PIXWAVS, UB,HKLS,NHKLS)
C
C Loop through HKLs, only keeping valid X,Y,WAV
C
	N=0
	DO I=1,NHKLS
C Skip if outside X,Y limits
	  IF(PIXWAVS(1,I).LT.X_MIN .OR. PIXWAVS(1,I).GT.X_MAX) CYCLE
	  IF(PIXWAVS(2,I).LT.Y_MIN .OR. PIXWAVS(2,I).GT.Y_MAX) CYCLE
C Skip if outside wavelength limits
	  WAV=ABS( PIXWAVS(3,I) )
	  IF(WAV.LT.WAV_MIN .OR. WAV.GT.WAV_MAX) CYCLE
C Flip the HKL to the Friedel if the wavelength is negative
	  IF(PIXWAVS(3,I) .LE. 0.0) THEN
	    HKLS(1,I)=-HKLS(1,I)
	    HKLS(2,I)=-HKLS(2,I)
	    HKLS(3,I)=-HKLS(3,I)
	    IHKLS(1,I)=-IHKLS(1,I)
	    IHKLS(2,I)=-IHKLS(2,I)
	    IHKLS(3,I)=-IHKLS(3,I)
	  ENDIF
C Shuffle HKLS,IHKLS,IMODS arrays (skipping invalid spots)
	  N=N+1
	  HKLS(1,N)=HKLS(1,I)
	  HKLS(2,N)=HKLS(2,I)
	  HKLS(3,N)=HKLS(3,I)
	  IHKLS(1,N)=IHKLS(1,I)
	  IHKLS(2,N)=IHKLS(2,I)
	  IHKLS(3,N)=IHKLS(3,I)
	  IMODS(N)=IMODS(I)
	ENDDO
	NHKLS=N
C
	RETURN
	END


C ================= Prune HKL list by merge/remove harmonic multiples =============

	SUBROUTINE PRUNE_MULTIPLES(HKLS,IHKLS,IMODS,IMULTS,NHKLS, UB)
C
C Find repeated reflections (probably harmonic multiples) with X & Y within
C one pixel of each other. Only keep the repeat with the largest wavelength.
C The number of repeats (mulitiplicity) is returned in IMULTS().
C If LMOD_WEAK is false the satellites are treated the same as main spots,
C otherwise the satellites are treated as much weak and they are ignored if
C they overlap a main spot and will not increase the multiplicity.
C
	INTEGER IHKLS(3,*),IMODS(*),IMULTS(*)
	REAL HKLS(3,*),UB(3,3)
C
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	LOGICAL LMOD_WEAK
	COMMON /MODULATE_COM/ NHKLS_MOD,DMIN_MOD,LMOD_WEAK,RHKLS_MOD(3,100)
C
	PARAMETER (NHKLS_MAX=300000)
	REAL PIXWAVS(3,NHKLS_MAX)
C
	INTEGER IXY_HKLS(8000,4000),IXY_MULT(8000,4000)
C
	IF(NUMX.GT.8000 .OR. NUMY.GT.4000) CALL QUIT('Detector size > 8000 x 4000 pixels')
C
C Zero IXY_HKLS() & IXY_MULT(), the XY maps for index to HKLS() and for multiplicity
C
	DO IY=1,NUMY
	  DO IX=1,NUMX
	    IXY_HKLS(IX,IY)=0
	    IXY_MULT(IX,IY)=0
	  ENDDO
	ENDDO
C
C Calculate pixel x,y and wavelengths
C
	CALL CALC_HKLS_TO_PIXWAVS(PIXWAVS, UB,HKLS,NHKLS)
C
C Set IXY_HKLS(IX,IY) to the index of HKLS() matching IX & IY. If more than
C one hkl at a pixel, update the index to the largest wavelength and store
C multiplicity in IXY_MULT().
C If LMOD_WEAK is set, satellite spots overlapping main spots are ignored.
C NB: It is assumed that main spots precede satellites.
C
	DO I=1,NHKLS
	  IX=NINT(PIXWAVS(1,I))
	  IY=NINT(PIXWAVS(2,I))
	  I2=IXY_HKLS(IX,IY)
C
C If (IX,IY) is empty, add index to hkl and set multiplicity to 1
	  IF(I2 .EQ. 0) THEN
	    IXY_HKLS(IX,IY)=I
	    IXY_MULT(IX,IY)=1
C If (IX,IY) is not empty:
	  ELSE
C If LMOD_WEAK is set and (IX,IY) is a main spot, ignore any satellite spots
	    IF(LMOD_WEAK .AND. (IMODS(I2).EQ.0) .AND. (IMODS(I).NE.0)) CYCLE
C Increment multiplicity, and update index if larger wavelength
	    IXY_MULT(IX,IY)=IXY_MULT(IX,IY)+1
	    IF(PIXWAVS(3,I) .GT. PIXWAVS(3,I2)) IXY_HKLS(IX,IY)=I
	  ENDIF
C
	ENDDO
C
C Shuffle hkls keeping only those that match IXY_HKLS()
C
	N=0
	DO I=1,NHKLS
	  IX=NINT(PIXWAVS(1,I))
	  IY=NINT(PIXWAVS(2,I))
C If index matches IXY_HKLS(), then it is the one to keep
	  IF(I .EQ. IXY_HKLS(IX,IY)) THEN
C Shuffle HKLS,IHKLS,IMODS arrays
	    N=N+1
	    HKLS(1,N)=HKLS(1,I)
	    HKLS(2,N)=HKLS(2,I)
	    HKLS(3,N)=HKLS(3,I)
	    IHKLS(1,N)=IHKLS(1,I)
	    IHKLS(2,N)=IHKLS(2,I)
	    IHKLS(3,N)=IHKLS(3,I)
	    IMODS(N)=IMODS(I)
C Record multiplicity in IMULTS
	    IMULTS(N)=IXY_MULT(IX,IY)
	  ENDIF
	ENDDO
	NHKLS=N
C
	RETURN
	END
