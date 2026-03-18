	SUBROUTINE FIND_SOLUTIONS(UB,XYCEN, B,NPEAKS_MATCH,NPEAKS_SOLN)
C
	REAL UB(3,3),XYCEN(2),B(3,3)
C
	COMMON /INSTRUM_COM/ ITYPE,NUMXY(2),XY_CEN(2),XY_SIZE(2),DRUM_RAD
C
	COMMON /SOLN_SELECT_COM/ ISOLN_SELECT
C
	PARAMETER (NSOLN_MAX=100000)
	REAL UB_SOLN(3,3,NSOLN_MAX),DXY_SOLN(2,NSOLN_MAX)
C
C Make "solutions" from the stored triplets
C The UB matrix for each solution is returned in UB_SOLN()
C
	CALL MAKE_SOLUTIONS(UB_SOLN,NSOLN, B)
	IF(NSOLN .LT. 1) GOTO 1000
C
C Remove solutions with identical angle differences between obs & calc
C
	CALL SOLN_PRUNE_REPEATS(UB_SOLN,NSOLN, NPEAKS_MATCH)
	IF(NSOLN .LT. 1) GOTO 1000
C
C Keep ~NTARGET solutions with the most matched obs & calc spots
C NB: Returns at most 1.5*NTARGET solutions
C
	NTARGET=100
	CALL SOLN_PRUNE_MATCHES(DXY_SOLN,UB_SOLN,NSOLN,NTARGET, NPEAKS_SOLN)
C
	CUTOFF=20.0
	CALL SOLN_PRUNE_EQUIVS(DXY_SOLN,UB_SOLN,NSOLN, CUTOFF,NPEAKS_SOLN)
C
C
	CALL SOLN_PRUNE_MERIT(DXY_SOLN, UB_SOLN,NSOLN, NPEAKS_SOLN)
	IF(NSOLN .LT. 1) GOTO 1000
C
	ISOLN_SELECT=MIN(NSOLN,ISOLN_SELECT)
	IF(ISOLN_SELECT .EQ. 1) THEN
	  PRINT *,'The highest merit solution is written to the output file'
	ELSE
	  PRINT '(1X,A,I2,A)','Number',ISOLN_SELECT,
	1	' highest merit solution is written to the output file'
	ENDIF
	CALL MX3_COPY(UB,UB_SOLN(1,1,ISOLN_SELECT))
	XYCEN(1)=XY_CEN(1)+DXY_SOLN(1,ISOLN_SELECT)
	XYCEN(2)=XY_CEN(2)+DXY_SOLN(2,ISOLN_SELECT)
	RETURN
C
C Complain and give up if we have no valid solutions!
C
1000	CALL QUIT('UNABLE TO FIND A VALID SOLUTION')
	END



	SUBROUTINE MAKE_SOLUTIONS(UB_SOLN,NSOLN, B)
C
	REAL B(3,3)
	PARAMETER (NSOLN_MAX=100000)
	REAL UB_SOLN(3,3,NSOLN_MAX)
C
	PARAMETER (NTRIP_MAX=100000)
	COMMON /TRIPLETS_COM/ NTRIP,IPEAKS_TRIP(3,NTRIP_MAX),
	1			IHKLS_TRIP(3,NTRIP_MAX)
C
	REAL U(3,3)
C
C Calculate U from a triplet and store in U_SOLN() if it is valid
C
	NSOLN=0
	DO I=1,NTRIP
	  CALL CALC_U_FROM_TRIPLET(U,IHKLS_TRIP(1,I),IPEAKS_TRIP(1,I), ISTATUS)
	  IF(ISTATUS .EQ. 1)THEN
	    NSOLN=NSOLN+1
	    CALL MX3_MULT(UB_SOLN(1,1,NSOLN), U,B)
	  ENDIF
	ENDDO
C
	IF(NSOLN .NE. NTRIP) PRINT '(1X,I5,A)',
	1	NTRIP-NSOLN,' triplets rejected as coplanar'
	RETURN
	END



	SUBROUTINE SOLN_ROBUST_REFINE(DXY,UB, HVECS,NPEAKS)
C
C Refine the UB rotation & DXY for a set of scattering vectors
C The method calculates a distance cutoff for matching spots
C As well as the cutoff, all possible combinations of two extra
C spots are also rejected
C The best refinement based on rms error is returned in DXY & UB
C
	REAL DXY(2),UB(3,3),HVECS(3,NPEAKS)
C
	PARAMETER (NSPOTS_MAX=10000)
	COMMON /SPOTS_COM/ NSPOTS,XSPOT(NSPOTS_MAX),YSPOT(NSPOTS_MAX),ESDS(NSPOTS_MAX)
C
	COMMON /HKL_BEST_COM/ IHKL_MAX1,IHKL_MAX2
C
	PARAMETER (NHKLS_MAX=1000)
	INTEGER ITAGS(NHKLS_MAX)
	REAL HKLS(3,NHKLS_MAX),DIST(NHKLS_MAX),DIST2(NHKLS_MAX)
	REAL XSPOT2(NHKLS_MAX),YSPOT2(NHKLS_MAX),HKLS2(3,NHKLS_MAX)
C
	REAL HVEC(3),UB_TRY(3,3),DXY_TRY(2),UB_BEST(3,3),DXY_BEST(2)
C
C Calculate the best HKL for the list of scattering vectors HVECS
C
	CALL CALC_BEST_HKLS(HKLS, UB,HVECS,NPEAKS,IHKL_MAX2)
C
C Calculate the obs. to calc. distance for all peaks
C
	DO I=1,NPEAKS
	  CALL VC3_MULT(HVEC, UB,HKLS(1,I))
	  CALL CALC_HVEC_TO_PIXEL(X0,Y0, HVEC)
	  DSQ=(XSPOT(I)-X0-DXY(1))**2 + (YSPOT(I)-Y0-DXY(2))**2
	  DIST(I)=SQRT(DSQ)
	ENDDO
C
C Do an "increasing" tag sort on the distance
C
	CALL TAG_SORT_REALS(ITAGS,DIST,NPEAKS)
C
C Calculate DIST2 = "(dy/dx)/y" of the sorted distance
C
	N23=NPEAKS*2/3
	DO I=1,N23
	  DIST2(I)=DIST(ITAGS(I+1))/DIST(ITAGS(I))
	ENDDO
C
C Find ICUT, the point of maximum DIST2
C
	ICUT=4
	DO I=5,N23
	  IF(DIST2(I) .GT. DIST2(ICUT)) ICUT=I
	ENDDO
C
C Create new arrays with only the ICUT+1 best spots
C
	NPEAKS2=ICUT
	DO I=1,NPEAKS2
	  I2=ITAGS(I)
	  XSPOT2(I)=XSPOT(I2)
	  YSPOT2(I)=YSPOT(I2)
	  HKLS2(1,I)=HKLS(1,I2)
	  HKLS2(2,I)=HKLS(2,I2)
	  HKLS2(3,I)=HKLS(3,I2)
	ENDDO
C
C Do a LSQ iteration on all of the ICUT spots (CUTOFF=1E5)
C
	CALL SETUP_LSQ(XSPOT2,YSPOT2,NPEAKS2)
	CALL ITER_LSQ(UB,DXY,RMS2,NSUM2, XSPOT2,YSPOT2,HKLS2,1.0E5,NPEAKS2)
C
C Try refinements with every combination of two spots removed
C The spots are "removed" by setting XSPOT to more than the CUTOFF
C
	RMS3B=1E7
	DO I1=1,NPEAKS2
	  X1=XSPOT2(I1)
	  XSPOT2(I1)=1.0E6
	  DO I2=1,I1-1
	    X2=XSPOT2(I2)
	    XSPOT2(I2)=1.0E6
C Restore the initial UB & DXY
	    CALL MX3_COPY(UB_TRY, UB)
	    DXY_TRY(1)=DXY(1)
	    DXY_TRY(2)=DXY(2)
	    CALL ITER_LSQ(UB_TRY,DXY_TRY,RMS3,NSUM3, XSPOT2,YSPOT2,HKLS2,1.0E5,NPEAKS2)
	    CALL ITER_LSQ(UB_TRY,DXY_TRY,RMS4,NSUM3, XSPOT2,YSPOT2,HKLS2,1.0E5,NPEAKS2)
C Save the best DXY & UB so far
	    IF(RMS4 .LT. RMS3B) THEN
	      RMS3B=RMS4
	      CALL MX3_COPY(UB_BEST, UB_TRY)
	      DXY_BEST(1)=DXY_TRY(1)
	      DXY_BEST(2)=DXY_TRY(2)
	    ENDIF
C Restore the original XSPOTS2 values
	    XSPOT2(I2)=X2
	  ENDDO
	  XSPOT2(I1)=X1
	ENDDO
C
C Copy the best refined values to UB & DXY
C
	CALL MX3_COPY(UB, UB_BEST)
	DXY(1)=DXY_BEST(1)
	DXY(2)=DXY_BEST(2)
C
	RETURN
	END



	SUBROUTINE SOLN_CALC_MERIT(RMERIT,DXY,UB, HVECS,NPEAKS)
C
C Calculates the merit of a solution, RDIST, which is the fraction
C of the NPEAKS that are matched to within 10 pixels using the
C given DXY & UB values
C
	REAL HVECS(3,NPEAKS),DXY(2),UB(3,3)
C
	PARAMETER (NSPOTS_MAX=10000)
	COMMON /SPOTS_COM/ NSPOTS,XSPOT(NSPOTS_MAX),YSPOT(NSPOTS_MAX),ESDS(NSPOTS_MAX)
C
	COMMON /HKL_BEST_COM/ IHKL_MAX1,IHKL_MAX2
C
	PARAMETER (NHKLS_MAX=1000)
	INTEGER ITAGS(NHKLS_MAX)
	REAL HKLS(3,NHKLS_MAX)
	REAL DIST(NHKLS_MAX),DIST2(NHKLS_MAX)
C
	REAL HVEC(3)
C
C Calculate the best hkls for the UB & HVECS
C
	CALL CALC_BEST_HKLS(HKLS, UB,HVECS,NPEAKS,IHKL_MAX2)
C
C Calculate the obs. to calc. distance of all peaks
C
	DO I=1,NPEAKS
	  CALL VC3_MULT(HVEC, UB,HKLS(1,I))
	  CALL CALC_HVEC_TO_PIXEL(X0,Y0, HVEC)
	  DSQ=(XSPOT(I)-X0-DXY(1))**2 + (YSPOT(I)-Y0-DXY(2))**2
	  DIST(I)=SQRT(DSQ)
	ENDDO
C
C Do an "increasing" sort on the distance
C
	CALL TAG_SORT_REALS(ITAGS,DIST,NPEAKS)
	DO I=1,NPEAKS
	  DIST2(I)=DIST(ITAGS(I))
	ENDDO
C
C Calculate when the sorted distance exceeds 10 pixels
C Use linear interpolation to get a fractional value
C
	RDIST=0.0
	DO I=1,NPEAKS-1
	  IF(DIST2(I).LT.10.0 .AND. DIST2(I+1).GE.10.0) THEN
	    RDIST=I+(10.0-DIST2(I))/(DIST2(I+1)-DIST2(I))
	    EXIT
	  ENDIF
	ENDDO
	IF(RDIST.EQ.0.0 .AND. DIST2(1).LT.10.0) RDIST=NPEAKS
C
C Finally, divide by NPEAKS to give the merit
C
	RMERIT=RDIST/NPEAKS
	RETURN
	END



	SUBROUTINE SOLN_CALC_DIFFS(DIFF, UB_SOLN,NSOLN, NPEAKS)
C
C Calculate the ave. angle difference between obs & calc scattering vectors
C
	PARAMETER (NSOLN_MAX=100000)
	REAL DIFF(NSOLN_MAX),UB_SOLN(3,3,NSOLN_MAX)
C
	COMMON /HKL_BEST_COM/ IHKL_MAX1,IHKL_MAX2
C
	COMMON /HKL_VECTORS_COM/ NHKLS,H_HKLS(3,1000)
C
	REAL H_PEAKS
	COMMON /PEAKS_COM/ H_PEAKS(3,10000)
C
	REAL HKLS(3,100),VEC(3)
C
	IF(NPEAKS .GT. 100) CALL QUIT('BUG(soln_calc_diffs) NPEAKS > 100')
C
C Calculate the average angle difference between peak and hkl
C scattering vectors for each triplet-solution and put in DIFF().
C
	DO I=1,NSOLN
C Determine best Miller indices HKLS() to match the scattering vectors
	  CALL CALC_BEST_HKLS(HKLS, UB_SOLN(1,1,I),H_PEAKS,NPEAKS,IHKL_MAX1)
C Calculate DIFF(), the average angle difference between obs & calc scattering vectors
	  DEV=0.0
	  DO I2=1,NPEAKS
	    CALL VC3_MULT(VEC, UB_SOLN(1,1,I),HKLS(1,I2))
	    CALL VC3_UNIT(VEC)
	    CALL VC3_SUB(VEC,H_PEAKS(1,I2),VEC)
	    DEV=DEV+VC3_SIZE(VEC)*57.3
	  ENDDO
	  DIFF(I)=DEV/NPEAKS
	ENDDO
C
	RETURN
	END



	SUBROUTINE CALC_BEST_HKLS(HKLS, UB,H_PEAKS,NPEAKS, IHKL_MAX)
C
C Given the UB matrix calculate the best matching Miller indices for the
C first NPEAKS peaks. The indices are stored as reals in HKLS().
C The absolute value of the indices will not exceed IHKL_MAX.
C
	REAL UB(3,3),HKLS(3,NPEAKS),H_PEAKS(3,NPEAKS)
C
	REAL UB_INV(3,3),HKL(3)
C
	CALL MX3_INVERT(UB_INV,DET, UB)
C DET is usually about 1e-3
	IF(ABS(DET) .LT. 1E-6) CALL QUIT('BUG(calc_best_hkls): DET < 1e-6')
C
C Loop through peaks, calculate best HKL and store in HKLS().
C The absolute value of indices will not exceed IHKL_MAX.
C
	DO I=1,NPEAKS
C Calculate the non-integer hkl that best matches the vector H_PEAK
	  CALL VC3_MULT(HKL, UB_INV,H_PEAKS(1,I))
C Scale hkl so its largest magnitude index is +/- 1
	  SCALE=MAX(ABS(HKL(1)),ABS(HKL(2)),ABS(HKL(3)))
	  HKL(1)=HKL(1)/SCALE
	  HKL(2)=HKL(2)/SCALE
	  HKL(3)=HKL(3)/SCALE
C Multiply hkl by 1 to IHKL_MAX and save the best multiplier
	  DIFF_BEST=1E6
	  DO IMUL=1,IHKL_MAX
	    DIFF=ABS(IMUL*HKL(1)-NINT(IMUL*HKL(1))) +
	1		ABS(IMUL*HKL(2)-NINT(IMUL*HKL(2))) +
	1		ABS(IMUL*HKL(3)-NINT(IMUL*HKL(3)))
	    IF(DIFF .LT. DIFF_BEST) THEN
	      DIFF_BEST=DIFF
	      IBEST=IMUL
	    ENDIF
	  ENDDO
C Save the "best integer" values for hkl
	  HKLS(1,I)=NINT(IBEST*HKL(1))
	  HKLS(2,I)=NINT(IBEST*HKL(2))
	  HKLS(3,I)=NINT(IBEST*HKL(3))
C Loop back for the next peak
	ENDDO
C
	RETURN
	END
