	SUBROUTINE SOLN_PRUNE_REPEATS(UB_SOLN,NSOLN, NPEAKS)
C
C Calculate the average difference between peak and hkl unit scattering
C vectors and remove solutions that have the same average difference.
C
	PARAMETER (NSOLN_MAX=100000)
	REAL UB_SOLN(3,3,NSOLN_MAX)
C
	EXTERNAL DIFF_ANGLE_COMPARE
C
	COMMON /DIFF_ANGLE_COM/ DIFF(NSOLN_MAX)
C
	LOGICAL LREPEAT(NSOLN_MAX)
	INTEGER ITAGS(NSOLN_MAX)
C
C Calculate the average angle difference between peak and hkl
C scattering vectors for all triplet-solutions and put in DIFF().
C
	CALL SOLN_CALC_DIFFS(DIFF, UB_SOLN,NSOLN, NPEAKS)
C
C Tagged sort on DIFF().
C
	CALL TAG_SORT(DIFF_ANGLE_COMPARE,ITAGS,NSOLN)
C
C Set LREPEAT() if DIFF is the same as the previous value
C
	TRIP_LAST=-1.0
	DO I=1,NSOLN
	  TRIP_NOW=DIFF(ITAGS(I))
	  LREPEAT(I)=(ABS(TRIP_NOW-TRIP_LAST) .LE. 1E-5)
	  TRIP_LAST=TRIP_NOW
	ENDDO
C
C Rebuild U_SOLN() with only the unique (non repeated) solutions
C
	N=0
	DO I=1,NSOLN
	  IF( .NOT.LREPEAT(I) ) THEN
	    N=N+1
	    CALL MX3_COPY(UB_SOLN(1,1,N),UB_SOLN(1,1,I))
	  ENDIF
	ENDDO
	NSOLN=N
C
C Output some info
C
	PRINT '(1X,I5,A)',NSOLN,' different solutions found'
	RETURN
	END



	SUBROUTINE SOLN_PRUNE_MATCHES(DXY_SOLN,U_SOLN,NSOLN,NTARGET, NPEAKS)
C
C Remove solutions with less spot matches until reach NTARGET, or less, solutions
C
	PARAMETER (NSOLN_MAX=100000)
	REAL DXY_SOLN(2,NSOLN_MAX),U_SOLN(3,3,NSOLN_MAX)
C
	PARAMETER (NSPOTS_MAX=10000)
	COMMON /SPOTS_COM/ NSPOTS,XSPOT(NSPOTS_MAX),YSPOT(NSPOTS_MAX),ESDS(NSPOTS_MAX)
C
	COMMON /HKL_BEST_COM/ IHKL_MAX1,IHKL_MAX2
C
	INTEGER NSUM(NSOLN_MAX)
	REAL HKLS(3,NSPOTS_MAX),HVECS(3,NSPOTS_MAX)
C
C Setup LSQ for NPEAKS spots
C
	CALL SETUP_LSQ(XSPOT,YSPOT,NPEAKS)
C
C Calculate the scattering vectors of spots from observed X,Y
C
	DO I=1,NPEAKS
	  CALL CALC_PIXEL_TO_HVEC(HVECS(1,I), XSPOT(I),YSPOT(I))
	ENDDO
C
C STEP 1: Calculate NSUM (how many matched spots) for each solution
C
C Try to refine UB orientation and DXY for all solutions
C
C Use an initial distance cutoff of 20 pixels for a match
	CUTOFF=20.0
100	DO I=1,NSOLN
C Start with a zero XY offset
	  DXY_SOLN(1,I)=0.0
	  DXY_SOLN(2,I)=0.0
C Calculate the most likely HKLs for this solution's UB
	  CALL CALC_BEST_HKLS(HKLS, U_SOLN(1,1,I),HVECS,NPEAKS,IHKL_MAX1)
C Refine UB rotation and DXY with matches within CUTOFF pixels
	  CALL ITER_LSQ(U_SOLN(1,1,I),DXY_SOLN(1,I),RMS1,NSUM(I),
	1				XSPOT,YSPOT,HKLS,CUTOFF,NPEAKS)
	  IF(NSUM(I) .GT. 0) THEN
	    DO ITER=2,4
	      CALL ITER_LSQ(U_SOLN(1,1,I),DXY_SOLN(1,I),RMS2,NSUM(I),
	1				XSPOT,YSPOT,HKLS,CUTOFF,NPEAKS)
	      IF(ABS(RMS1-RMS2) .LT. 0.1) EXIT
	      RMS1=RMS2
	    ENDDO
	  ENDIF
C Loop back for the next solution
	ENDDO
C
C Find the smallest ICUTOFF where we have less than NTARGET
C occurrences of NSUM() > ICUTOFF
C
	DO ICUTOFF=5,100
	  NTOT=0
	  DO I=1,NSOLN
	    IF(NSUM(I) .GT. ICUTOFF) NTOT=NTOT+1
	  ENDDO
	  IF(NTOT .LT. NTARGET) EXIT
	ENDDO
C
C If NTOT is much less than NTARGET, repeat with larger CUTOFF
C
	IF(NTOT .LT. NTARGET/2) THEN
	  CUTOFF=CUTOFF+5.0
	  IF(CUTOFF .LE. 50.1) GOTO 100
	ENDIF
C
C STEP 2: Remove solutions with NSUM() < ICUTOFF
C
	N=0
	DO I=1,NSOLN
	  IF(NSUM(I) .LT. ICUTOFF) CYCLE
	  IF(N .GT. NTARGET*3/2) EXIT
	  N=N+1
	  NSUM(N)=NSUM(I)
	  DXY_SOLN(1,N)=DXY_SOLN(1,I)
	  DXY_SOLN(2,N)=DXY_SOLN(2,I)
	  CALL MX3_COPY(U_SOLN(1,1,N), U_SOLN(1,1,I))
	ENDDO
	NSOLN=N
	PRINT '(1X,I5,A,I2,A)',NSOLN,' solutions with ',ICUTOFF,
	1				' or more matched spots'
C
	RETURN
	END



	SUBROUTINE SOLN_PRUNE_EQUIVS(DXY_SOLN,U_SOLN,NSOLN, CUTOFF,NPEAKS)
C
	PARAMETER (NSOLN_MAX=100000)
	REAL DXY_SOLN(2,NSOLN_MAX),U_SOLN(3,3,NSOLN_MAX)
C
	PARAMETER (NSPOTS_MAX=10000)
	COMMON /SPOTS_COM/ NSPOTS,XSPOT(NSPOTS_MAX),YSPOT(NSPOTS_MAX),ESDS(NSPOTS_MAX)
C
	COMMON /HKL_BEST_COM/ IHKL_MAX1,IHKL_MAX2
C
	PARAMETER (NHKLS_MAX=1000)
	LOGICAL LINCLUDE(NHKLS_MAX)
	INTEGER NSUM(NHKLS_MAX),IHASH(NHKLS_MAX)
	REAL HKLS(3,NHKLS_MAX),HVECS(3,NHKLS_MAX),RMS(NHKLS_MAX)
C
	REAL HVEC(3)
C
C Calculate the scattering vectors for all observed spots
C
	DO I=1,NPEAKS
	  CALL CALC_PIXEL_TO_HVEC(HVECS(1,I), XSPOT(I),YSPOT(I))
	ENDDO
C
C For each solution create a hash number corresponding to which
C observed spot is matched within CUTOFF pixels
C NB: I use a 32 bit hash, so only the first 31 spots are used
C
	DO I=1,NSOLN
	  CALL CALC_BEST_HKLS(HKLS, U_SOLN(1,1,I),HVECS,NPEAKS,IHKL_MAX1)
	  IHASH(I)=0
	  IBIT=1
	  DO I2=1,NPEAKS
	    CALL VC3_MULT(HVEC, U_SOLN(1,1,I),HKLS(1,I2))
	    CALL CALC_HVEC_TO_PIXEL(X0,Y0, HVEC)
	    DSQ=(XSPOT(I2)-X0-DXY_SOLN(1,I))**2 + (YSPOT(I2)-Y0-DXY_SOLN(2,I))**2
	    LINCLUDE(I2)=(DSQ .LE. CUTOFF**2)
	    IF(IBIT .GT. 0) THEN
	      IF( LINCLUDE(I2) ) IHASH(I)=IHASH(I)+IBIT
	      IBIT=IBIT*2
	    ENDIF
	  ENDDO
	ENDDO
C
C Set HASH=0 for repeated solutions with the same HASH,NSUM and RMS values
C
	DO I=1,NSOLN
	  IF(IHASH(I) .EQ. 0) CYCLE
	  DO I2=I+1,NSOLN
	    IF( IHASH(I2).EQ.IHASH(I) .AND.
	1        NSUM(I2).EQ.NSUM(I)   .AND.
	2        ABS(RMS(I2)-RMS(I)).LT.0.05 ) IHASH(I2)=0
	  ENDDO
	ENDDO
C
C Remove any solutions with HASH=0
C
	N=0
	DO I=1,NSOLN
	  IF(IHASH(I) .EQ. 0) CYCLE
	  N=N+1
	  NSUM(N)=NSUM(I)
	  RMS(N)=RMS(I)
	  DXY_SOLN(1,N)=DXY_SOLN(1,I)
	  DXY_SOLN(2,N)=DXY_SOLN(2,I)
	  CALL MX3_COPY(U_SOLN(1,1,N), U_SOLN(1,1,I))
	ENDDO
	NSOLN=N
C
	PRINT '(1X,I5,A)',NSOLN,' solutions after removing equivalents'
C
	RETURN
	END



	SUBROUTINE SOLN_PRUNE_MERIT(DXY_SOLN, UB_SOLN,NSOLN, NPEAKS)
C
	PARAMETER (NSOLN_MAX=100000)
	REAL DXY_SOLN(2,NSOLN_MAX),UB_SOLN(3,3,NSOLN_MAX)
C
	PARAMETER (NSPOTS_MAX=10000)
	COMMON /SPOTS_COM/ NSPOTS,XSPOT(NSPOTS_MAX),YSPOT(NSPOTS_MAX),ESDS(NSPOTS_MAX)
C
	PARAMETER (NHKLS_MAX=1000)
	REAL HVECS(3,NHKLS_MAX)
C
	PARAMETER (NSOLN2_MAX=1000)
	INTEGER ITAGS(NSOLN2_MAX)
	REAL RMERIT(NSOLN2_MAX)
C
C Calculate the scattering vectors for all observed spots
C
	DO I=1,NPEAKS
	  CALL CALC_PIXEL_TO_HVEC(HVECS(1,I), XSPOT(I),YSPOT(I))
	ENDDO
C
C Loop through solutions, refining DXY_SOLN & UB_SOLN, then using these
C values to calculate the merit (0 to 1) of each solution
C
	DO I=1,NSOLN
	  CALL SOLN_ROBUST_REFINE(DXY_SOLN(1,I),UB_SOLN(1,1,I), HVECS,NPEAKS)
	  CALL SOLN_CALC_MERIT(RMERIT(I),DXY_SOLN(1,I),UB_SOLN(1,1,I), HVECS,NPEAKS)
	ENDDO
C
C Do an "increasing" tag sort on the MERIT
C
	CALL TAG_SORT_REALS(ITAGS,RMERIT,NSOLN)
C
C Copy the sorted list (only the 10 largest merits) to
C array elements (NSOLN+1 .. NSOLN+10)
C
	DO I=1,MIN(10,NSOLN)
	  I2=ITAGS(NSOLN-I+1)
	  DXY_SOLN(1,NSOLN+I)=DXY_SOLN(1,I2)
	  DXY_SOLN(2,NSOLN+I)=DXY_SOLN(2,I2)
	  CALL MX3_COPY(UB_SOLN(1,1,NSOLN+I), UB_SOLN(1,1,I2))
	  RMERIT(NSOLN+I)=RMERIT(I2)
	ENDDO
C
C Now copy the sorted values to the start of the arrays
C
	DO I=1,MIN(10,NSOLN)
	  DXY_SOLN(1,I)=DXY_SOLN(1,NSOLN+I)
	  DXY_SOLN(2,I)=DXY_SOLN(2,NSOLN+I)
	  CALL MX3_COPY(UB_SOLN(1,1,I), UB_SOLN(1,1,NSOLN+I))
	  RMERIT(I)=RMERIT(NSOLN+I)
	ENDDO
C
	IF(NSOLN .LE. 10) THEN
	  PRINT '(/,1X,A,I2,A)','The merit of the final',NSOLN,' solutions: (higher is better)'
	ELSE
	  PRINT '(/,1X,A)','The merit of the top 10 solutions: (higher is better)'
	  NSOLN=10
	ENDIF
	PRINT '(3X,10I4)',(NINT(100.0*RMERIT(K)),K=1,NSOLN)
C
	RETURN
	END
