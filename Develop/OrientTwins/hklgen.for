C ================= Common to ORIENT_SPOTS & HKL_GEN ============

	SUBROUTINE GEN_HEMISPHERE(NHKLS, UB,ICEN,D_MIN,D_MAX)
C
C Calculate a list of HKLs using the UB matrix,the minimum d-spacing D_MIN,
C and the centering given by ICEN = 0(P) 1(A) 2(B) 3(C) 4(I) 5(F) 6(R).
C
C The list is for half of all reflections as it excludes Friedel pairs.
C The hemisphere is defined by (H>0,K,L) & (0,K>0,L) & (0,0,L>0)
C
	REAL UB(3,3)
C
	INTEGER HKLS
	COMMON /HKLGEN_COM/ HKLS(3,100000),PIXWAVS(3,100000),IMULTS(100000)
C
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
	TEMP=(DOT11*DOT33-DOT13**2) -
	1	(DOT33*DOT12-DOT13*DOT23)**2 / (DOT22*DOT33-DOT23**2)
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
	  TEMP2=(DOT11-DOT13**2/DOT33)*(DOT22-DOT23**2/DOT33)
	1		-(DOT12-DOT13*DOT23/DOT33)**2
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
	    TEMP=DOT33*(SSQ_MAX-DOT11*H**2-2*DOT12*H*K-DOT22*K**2) +
	1			(DOT13*H+DOT23*K)**2
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
C Loop over L and store any hkl with |UB.HKL|^2 > SSQ_MIN
	    DO L=L_LO,L_HI,L_INC
	      SSQ= DOT11*H*H + DOT12*H*K + DOT13*H*L + 
	1	   DOT12*K*H + DOT22*K*K + DOT23*K*L + 
	2	   DOT13*L*H + DOT23*L*K + DOT33*L*L
	      IF(SSQ .GE. SSQ_MIN) THEN
	        IF(NHKLS .EQ. 100000) RETURN
	        NHKLS=NHKLS+1
	        HKLS(1,NHKLS)=H
	        HKLS(2,NHKLS)=K
	        HKLS(3,NHKLS)=L
	      ENDIF
C
C End of loops over H, K & L
C
	    ENDDO
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE PRUNE_MULTIPLES(NHKLS)
C
	EXTERNAL XY_COMPARE
C
C Need PIXWAV in common so that XY_COMPARE() can access it
C
	INTEGER HKLS
	COMMON /HKLGEN_COM/ HKLS(3,100000),PIXWAVS(3,100000),IMULTS(100000)
C
	INTEGER ITAGS(100000)
C
C Zero the multiplicity values for all hkls
C
	DO I=1,NHKLS
	  IMULTS(I)=0
	ENDDO
C
C Do a tagged sort of hkls based on differences in X & Y
C
	CALL TAG_SORT(XY_COMPARE,ITAGS,NHKLS)
C
C Find groups of hkls with the "same" X & Y and set IMULTS()
C for the hkl with the largest wavelength in the group
C
	IFIRST=1
	DO I=2,NHKLS
C If ITAGS(I) and ITAGS(IFIRST) have the "different" X & Y,
C then hkls IFIRST to I-1 are a group of multiples hkls
C NB: XY_COMPARE() return zero if spot separation < 0.5 pixels
	  IF(XY_COMPARE(ITAGS(I),ITAGS(IFIRST)) .NE. 0.0) THEN
	    ILAST=I-1
C Find the related hkl in the group with the longest wavelength
	    ILONG=IFIRST
	    DO I2=IFIRST+1,ILAST
	      IF(PIXWAVS(3,ITAGS(I2)) .GT. PIXWAVS(3,ITAGS(ILONG))) ILONG=I2
	    ENDDO
C Set IMULTS() for the longest wavelength hkl in the group
	    IMULTS(ITAGS(ILONG))=ILAST-IFIRST+1
C Update IFIRST to the "different" hkl
	    IFIRST=I
	  ENDIF
	ENDDO
C
C Set IMULTS() for the final group
C
	ILONG=IFIRST
	DO I2=IFIRST+1,NHKLS
	  IF(PIXWAVS(3,ITAGS(I2)) .GT. PIXWAVS(3,ITAGS(ILONG))) ILONG=I2
	ENDDO
	IMULTS(ITAGS(ILONG))=NHKLS-IFIRST+1
C
C Remove from the arrays any hkls with zero multiplicity
C
	N=0
	DO I=1,NHKLS
	  IF(IMULTS(I) .NE. 0) THEN
	    N=N+1
	    DO I2=1,3
	      HKLS(I2,N)=HKLS(I2,I)
	      PIXWAVS(I2,N)=PIXWAVS(I2,I)
	    ENDDO
	    IMULTS(N)=IMULTS(I)
	  ENDIF
	ENDDO
	NHKLS=N
C
	RETURN
	END


	SUBROUTINE PRUNE_DETECTED(NHKLS, UB,WAV_MIN,WAV_MAX,X_MIN,X_MAX,Y_MIN,Y_MAX)
C
	REAL UB(3,3)
C
	INTEGER HKLS
	COMMON /HKLGEN_COM/ HKLS(3,100000),PIXWAVS(3,100000),IMULTS(100000)
C
	CALL CALC_INT_HKLS_TO_PIXWAVS(PIXWAVS, UB,HKLS,NHKLS)
C
C Loop through HKLs, only keeping valid X,Y,WAV and
C inverting any H,K,L,WAV with a negative WAV value
C
	N=0
	DO I=1,NHKLS
	  IF(PIXWAVS(1,I).LT.X_MIN .OR. PIXWAVS(1,I).GT.X_MAX .OR.
	1     PIXWAVS(2,I).LT.Y_MIN .OR. PIXWAVS(2,I).GT.Y_MAX .OR.
	2     ABS(PIXWAVS(3,I)).LT.WAV_MIN .OR.
	3     ABS(PIXWAVS(3,I)).GT.WAV_MAX                    ) CYCLE
C
	  N=N+1
	  IF(PIXWAVS(3,I) .GT. 0.0) THEN
	    HKLS(1,N)=HKLS(1,I)
	    HKLS(2,N)=HKLS(2,I)
	    HKLS(3,N)=HKLS(3,I)
	    PIXWAVS(1,N)=PIXWAVS(1,I)
	    PIXWAVS(2,N)=PIXWAVS(2,I)
	    PIXWAVS(3,N)=PIXWAVS(3,I)
	  ELSE
	    HKLS(1,N)=-HKLS(1,I)
	    HKLS(2,N)=-HKLS(2,I)
	    HKLS(3,N)=-HKLS(3,I)
	    PIXWAVS(1,N)=PIXWAVS(1,I)
	    PIXWAVS(2,N)=PIXWAVS(2,I)
	    PIXWAVS(3,N)=-PIXWAVS(3,I)
	  ENDIF
	ENDDO
C
	NHKLS=N
	RETURN
	END


	SUBROUTINE CALC_INT_HKLS_TO_PIXWAVS(PIXWAVS, UB,IHKLS,NHKLS)
C
	INTEGER IHKLS(3,100000)
	REAL PIXWAVS(3,100000),UB(3,3)
C
	REAL RHKLS(3,100000)
C
	DO I=1,NHKLS
	   RHKLS(1,I)=IHKLS(1,I)
	   RHKLS(2,I)=IHKLS(2,I)
	   RHKLS(3,I)=IHKLS(3,I)
	ENDDO
C
	CALL CALC_HKLS_TO_PIXWAVS(PIXWAVS, UB,RHKLS,NHKLS)
C
	RETURN
	END


	FUNCTION XY_COMPARE(I1,I2)
C
C Return the difference in X, or Y if the X values agree, or zero
C if both X & Y agree. Agreement in X or Y is set to 0.5 pixels.
C
	INTEGER HKLS
	COMMON /HKLGEN_COM/ HKLS(3,100000),PIXWAVS(3,100000),IMULT(100000)
C
C Set agreement cutoff to 0.5 pixels
C
	CUTOFF=0.5
C
	DX=PIXWAVS(1,I1)-PIXWAVS(1,I2)
	IF(ABS(DX) .GE. CUTOFF) THEN
	  XY_COMPARE=DX
	ELSE
	  DY=PIXWAVS(2,I1)-PIXWAVS(2,I2)
	  IF(ABS(DY) .GE. CUTOFF) THEN
	    XY_COMPARE=DY
	  ELSE
	    XY_COMPARE=0.0
	  ENDIF
	ENDIF
C
	RETURN
	END