C========== All routines specific to the initial spot orientation ==========
C	SUBROUTINE ORIENT_SPOTS_LEVEL1(RMS, UB, LCEN_FIXED)
C	SUBROUTINE MATCH_LEVEL1_TEST1(RMS,N_BEST, UB,XY_CEN, LCEN_FIXED)
C	SUBROUTINE MATCH_LEVEL1_TEST2(RMS,N_BEST, UB,XY_CEN, LCEN_FIXED)
C	SUBROUTINE MATCH_LEVEL1_TEST3(RMS,N_BEST, UB,XY_CEN, LCEN_FIXED)
C	SUBROUTINE MATCH_LEVEL1_TEST4(RMS,NMATCH2, UB,XY_CEN, LCEN_FIXED)
C	INTEGER FUNCTION MAKE_INITIAL_MATCH(UB,HSQ_MAX,DIST_MAX)
C	SUBROUTINE REFINE_ROT_UB(UB_OUT, UB)
C	SUBROUTINE REFINE_XY(XY, UB_IN,XY_IN)
C	SUBROUTINE REFINE_PIX_CEN(UB)
C ==========================================================================

	SUBROUTINE ORIENT_SPOTS_LEVEL1(RMS, UB, LCEN_FIXED)
C
	LOGICAL LCEN_FIXED
	REAL UB(3,3)
C
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
C
	INTEGER NMATCHS(4)
	REAL PIX_CEN0(2),UBS(3,3,4),PIX_CENS(2,4),RMSS(4)
C
C Save a copy of PIX_CEN in /PIXEL_COM/
C
	PIX_CEN0(1)=PIX_CEN(1)
	PIX_CEN0(2)=PIX_CEN(2)
C
C Make four copies of PIX_CEN & UB for each test
C
	DO I=1,4
	  CALL MX3_COPY(UBS(1,1,I), UB)
	  PIX_CENS(1,I)=PIX_CEN(1)
	  PIX_CENS(2,I)=PIX_CEN(2)
	ENDDO
C
C Run all four matching tests
C
	CALL MATCH_LEVEL1_TEST1(RMSS(1),NMATCHS(1), UBS(1,1,1),PIX_CENS(1,1), LCEN_FIXED)
	CALL MATCH_LEVEL1_TEST2(RMSS(2),NMATCHS(2), UBS(1,1,2),PIX_CENS(1,2), LCEN_FIXED)
	CALL MATCH_LEVEL1_TEST3(RMSS(3),NMATCHS(3), UBS(1,1,3),PIX_CENS(1,3), LCEN_FIXED)
	CALL MATCH_LEVEL1_TEST4(RMSS(4),NMATCHS(4), UBS(1,1,4),PIX_CENS(1,4), LCEN_FIXED)
C
C Determine the best matches to use
C
	NBEST=MAX(NMATCHS(1),NMATCHS(2),NMATCHS(3),NMATCHS(4))
	RBEST=MIN(RMSS(1),RMSS(2),RMSS(3),RMSS(4))
C Best NMATCHS if within 3 pixels of best RMSS
	DO I=1,4
	  IF(NMATCHS(I).EQ.NBEST .AND. RMSS(I).LE.RBEST+3.0) GOTO 100
	ENDDO
C First match within 2 of best NMATCHS (and > 9) and 3 pixels of best RMSS
	DO I=1,4
	  IF(NMATCHS(I).GE.MAX(10,NBEST-2) .AND. RMSS(I).LE.RBEST+3.0) GOTO 100
	ENDDO
C First match within 75% of best NMATCHS (and > 9) and 3 pixels of the best RMSS
	DO I=1,4
	  IF(NMATCHS(I).GE.MAX(10,NBEST*3/4) .AND. RMSS(I).LE.RBEST+3.0) GOTO 100
	ENDDO
C First match with >9 NMATCHS and RMSS < 5
	DO I=1,4
	  IF(NMATCHS(I).GT.9 .AND. RMSS(I).LE.5.0) GOTO 100
	ENDDO
C First match with >5 NMATCHS and RMSS < 10
	DO I=1,4
	  IF(NMATCHS(I).GT.5 .AND. RMSS(I).LE.10.0) GOTO 100
	ENDDO
C First match with >3 NMATCHS and RMSS < 20
	DO I=1,4
	  IF(NMATCHS(I).GT.3 .AND. RMSS(I).LE.20.0) GOTO 100
	ENDDO
C First match with >2 NMATCHS
	DO I=1,4
	  IF(NMATCHS(I) .GT. 2) GOTO 100
	ENDDO
C Not even 3 matches, so give up.
	STOP 'ERROR: Unable to match 3 spots in initial phase'
C
C Copy optimised UB & PIX_CEN from selected test results
C
100	PIX_CEN(1)=PIX_CENS(1,I)
	PIX_CEN(2)=PIX_CENS(2,I)
	CALL MX3_COPY(UB,UBS(1,1,I))
	RMS=RMSS(I)
	PRINT '(3X,A,I2)','Using matches from Test',I
C
	RETURN
	END


	SUBROUTINE MATCH_LEVEL1_TEST1(RMS,N_BEST, UB,XY_CEN, LCEN_FIXED)
C
C Match low index hkls with medium mismatch distance while
C optimising PIX_CEN.
C This test should work for non-pathological problems.
C
	LOGICAL LCEN_FIXED
	REAL UB(3,3),XY_CEN(2)
C
	INTEGER MAKE_INITIAL_MATCH
C
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
C
	REAL UB_BEST(3,3),XY_BEST(2)
C
	REAL DX(21),DY(21)
	DATA DX/   0.0,  5.0, -5.0, 10.0,-10.0, 15.0,-15.0,
	1          0.0,  5.0, -5.0, 10.0,-10.0, 15.0,-15.0,
	2          0.0,  5.0, -5.0, 10.0,-10.0, 15.0,-15.0/
	DATA DY/   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
	1          5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,
	2         -5.0, -5.0, -5.0, -5.0, -5.0, -5.0, -5.0/
C
C Load PIX_CEN from XY_CEN
C
	PIX_CEN(1)=XY_CEN(1)
	PIX_CEN(2)=XY_CEN(2)
C
C Loop through shifts in PIX_CEN to find a good match
C at the next level.
C
	RMS=99.9
	N_BEST=-1
	ILOOP=21
	IF( LCEN_FIXED ) ILOOP=1
	DO I=1,ILOOP
	  PIX_CEN(1)=XY_CEN(1)+DX(I)
	  PIX_CEN(2)=XY_CEN(2)+DY(I)
C Try |h|^2 <= 6, matching error of 20 pixels
	  NMATCH1=MAKE_INITIAL_MATCH(UB,6.1,20.0)
C Ignore if < 3 matches
	  IF(NMATCH1 .LT. 3) CYCLE
C Optimise PIX_CEN, but only up to 5 pixels
	  IF( .NOT.LCEN_FIXED ) CALL REFINE_PIX_CEN(UB)
C Repeat matching spots, and calculate its RMS
	  NMATCH=MAKE_INITIAL_MATCH(UB,6.1,20.0)
	  RMS=CALC_RMS_ERR(UB)
C If the rms is bad, ignore it
	  IF(RMS .GT. 5.0) CYCLE
C Record the best match so far
	  IF(NMATCH .GT. N_BEST) THEN
	    N_BEST=NMATCH
	    RMS_BEST=RMS
	    CALL MX3_COPY(UB_BEST,UB)
	    XY_BEST(1)=PIX_CEN(1)
	    XY_BEST(2)=PIX_CEN(2)
C If matched 30, jump out of loop
	    IF(N_BEST .GE. 30) EXIT
	  ENDIF
	ENDDO
C
C Copy the best values back to the arguments
C
	CALL MX3_COPY(UB,UB_BEST)
	XY_CEN(1)=XY_BEST(1)
	XY_CEN(2)=XY_BEST(2)
C
C Output the results
C
	IF( N_BEST .GT. 0) THEN
	  PRINT '(3X,A,I3,A,F6.3,A)','Test 1:',N_BEST,' matches, rms=',
	1			RMS_BEST,' pixels for HSQ= 6 DIST=  20'
	ELSE
	  PRINT '(3X,A)','Test 1: Failed'
	ENDIF
C
	RETURN
	END


	SUBROUTINE MATCH_LEVEL1_TEST2(RMS,N_BEST, UB,XY_CEN, LCEN_FIXED)
C
C Match very low index hkls with increasing mismatch distance
C while optimising both PIX_CEN and rotations of UB.
C This test is designed for large jumps in crystal rotations.
C
	LOGICAL LCEN_FIXED
	REAL UB(3,3),XY_CEN(2)
C
	INTEGER MAKE_INITIAL_MATCH
C
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
C
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH000
C
	REAL UB_BEST(3,3),XY_BEST(2),UB2(3,3)
C
C Try |h|^2 <= 3 and a mismatch distance of 30 to 200 until
C we get a good match
C
	RMS=99.9
	N_BEST=-1
	DO DIST_MAX=30.0,200.1,10.0
C Restore PIX_CEN and do an initial match
	  PIX_CEN(1)=XY_CEN(1)
	  PIX_CEN(2)=XY_CEN(2)
	  NMATCH=MAKE_INITIAL_MATCH(UB,3.1,DIST_MAX)
C Ignore if < 3 matches
	  IF(NMATCH .LT. 3) CYCLE
C Refine UB rotation & PIX_CEN
	  CALL REFINE_ROT_UB(UB2, UB)
	  IF( .NOT.LCEN_FIXED ) CALL REFINE_PIX_CEN(UB2)
C Repeat matching spots with smaller mismatch, and calculate its RMS
	  NMATCH=MAKE_INITIAL_MATCH(UB,3.1,DIST_MAX/2.0)
	  RMS=CALC_RMS_ERR(UB)
C If the rms is bad, ignore it
	  IF(RMS .GT. 20.0) CYCLE
C Record the best match so far
	  IF(NMATCH .GT. N_BEST) THEN
	    N_BEST=NMATCH
	    RMS_BEST=RMS
	    DIST_BEST=DIST_MAX
	    CALL MX3_COPY(UB_BEST,UB)
	    XY_BEST(1)=PIX_CEN(1)
	    XY_BEST(2)=PIX_CEN(2)
C If matched 30, jump out of loop
	    IF(N_BEST .GE. 30) EXIT
	  ENDIF
	ENDDO
C
C Copy the best values back to the arguments
C
	CALL MX3_COPY(UB,UB_BEST)
	XY_CEN(1)=XY_BEST(1)
	XY_CEN(2)=XY_BEST(2)
C
C Print out summary of test
C
	IF( N_BEST .GT. 0) THEN
	  PRINT '(3X,A,I3,A,F6.3,A,I4)','Test 2:',N_BEST,
	1	' matches, rms=',RMS_BEST,' pixels for HSQ= 3 DIST=',NINT(DIST_BEST)
	ELSE
	  PRINT '(3X,A)','Test 2: Failed'
	ENDIF
C
	RETURN
	END



	SUBROUTINE MATCH_LEVEL1_TEST3(RMS,N_BEST, UB,XY_CEN, LCEN_FIXED)
C
C Match increasingly higher order hkls using a medium mismatch
C distance while optimising both PIX_CEN and rotations of UB.
C This test is designed to deal with weak supercells.
C
	LOGICAL LCEN_FIXED
	REAL UB(3,3),XY_CEN(2)
C
	INTEGER MAKE_INITIAL_MATCH
C
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
C
	REAL UB_BEST(3,3),XY_BEST(2),UB2(3,3)
C
C Try |h|^2 <= 9 to 14 with a mismatch distance of 20 until
C we get a good match
C
	RMS=99.9
	N_BEST=-1
	DO HSQ_MAX=9.0,14.1,1.0
C Restore PIX_CEN then do a level 1 match
	  PIX_CEN(1)=XY_CEN(1)
	  PIX_CEN(2)=XY_CEN(2)
	  NMATCH=MAKE_INITIAL_MATCH(UB,HSQ_MAX,20.0)
C Ignore if < 3 matches
	  IF(NMATCH .LT. 3) CYCLE
C Refine UB rotation & PIX_CEN
	  CALL REFINE_ROT_UB(UB2, UB)
	  IF( .NOT.LCEN_FIXED ) CALL REFINE_PIX_CEN(UB2)
C Repeat matching spots, and calculate its RMS
	  NMATCH=MAKE_INITIAL_MATCH(UB,HSQ_MAX,20.0)
C Ignore if < 3 matches
	  IF(NMATCH .LT. 3) CYCLE
C Ignore if the rms is bad
	  RMS=CALC_RMS_ERR(UB)
	  IF(RMS .GT. 10.0) CYCLE
C Record the best match so far
	  IF(NMATCH .GT. N_BEST) THEN
	    N_BEST=NMATCH
	    RMS_BEST=RMS
	    HSQ_BEST=HSQ_MAX
	    CALL MX3_COPY(UB_BEST,UB)
	    XY_BEST(1)=PIX_CEN(1)
	    XY_BEST(2)=PIX_CEN(2)
C If matched 30, jump out of loop
	    IF(N_BEST .GE. 30) EXIT
	  ENDIF
	ENDDO
C
C Copy the best values back to the arguments
C
	CALL MX3_COPY(UB,UB_BEST)
	XY_CEN(1)=XY_BEST(1)
	XY_CEN(2)=XY_BEST(2)
C
C Print out summary of test
C
	IF( N_BEST .GT. 0) THEN
	  PRINT '(3X,A,I3,A,F6.3,A,I4)','Test 3:',N_BEST,
	1	' matches, rms=',RMS_BEST,' pixels for HSQ=',NINT(HSQ_BEST)
	ELSE
	  PRINT '(3X,A)','Test 3: Failed'
	ENDIF
C
	RETURN
	END


	SUBROUTINE MATCH_LEVEL1_TEST4(RMS,NMATCH2, UB,XY_CEN, LCEN_FIXED)
C
C This test is completely different to the others.
C It trys to deduce the integer HKL for the 50 strongest spots then
C returns UB & XY_CEN for matches with 2 sigma outliers in x,y removed.
C
	LOGICAL LCEN_FIXED
	REAL UB(3,3),XY_CEN(2)
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
C
	INTEGER ITAGS(100)
	REAL PIXWAV(3,100),DIFF_BEST(100),HKL(3),RHKL(3),TEMP(3),UB2(3,3)
C
C Load PIX_CEN from XY_CEN
C
	PIX_CEN(1)=XY_CEN(1)
	PIX_CEN(2)=XY_CEN(2)
C
C Try to match half the spots, but not <3 or >30, from
C the 50 strongest spots
C
	NMATCH=MAX(3,MIN(30, NSPOTS/2 ))
	NSTRONG=MIN(50,NSPOTS)
C
C Calculate best integer HKLs
C
	WAV_MIN=1.0
	DO I=1,NSTRONG
C Calculate HKL assuming the minimum wavelength possible
	  PIXWAV(1,I)=SPOTS(1,I)
	  PIXWAV(2,I)=SPOTS(2,I)
	  PIXWAV(3,I)=WAV_MIN
	  CALL CALC_PIXWAVS_TO_HKLS(HKL, UB,PIXWAV(1,I),1)
C Find best integer HKL with same ratios but small indices
	  RMAX=MAX(ABS(HKL(1)),ABS(HKL(2)),ABS(HKL(3)))
	  DO I2=1,3
	    RHKL(I2)=HKL(I2)/RMAX
	  ENDDO
	  DIFF_BEST(I)=1E6
	  DO ITRY=1,IFIX(RMAX)
	    DO I2=1,3
	      HKL(I2)=NINT( RHKL(I2)*ITRY )
	    ENDDO
	    CALL CALC_HKLS_TO_PIXWAVS(TEMP, UB,HKL,1)
	    DIFF=SQRT( (PIXWAV(1,I)-TEMP(1))**2 +  (PIXWAV(2,I)-TEMP(2))**2 )
	    IF(DIFF .LT. DIFF_BEST(I)) THEN
	      DIFF_BEST(I)=DIFF
	      IBEST=ITRY
	      PIXWAV(3,I)=TEMP(3)
	    ENDIF
	  ENDDO
C Copy best integer HKL into HKL_MATCH()
	  DO I2=1,3
	    HKL_MATCH(I2,I)=NINT( RHKL(I2)*ITRY )
	  ENDDO
	ENDDO
C
C Do tagged bubble-sort on calculated wavelengths
C
	DO I=1,NSTRONG
	  ITAGS(I)=I
	ENDDO
C
	DO I1=1,NSTRONG
	  DO I2=1,NSTRONG
	    IF(PIXWAV(3,ITAGS(I1)) .GT. PIXWAV(3,ITAGS(I2))) THEN
	      ITEMP=ITAGS(I1)
	      ITAGS(I1)=ITAGS(I2)
	      ITAGS(I2)=ITEMP
	    ENDIF
	  ENDDO
	ENDDO
C
C Copy NMATCH spots with largest wavelengths into the matching lists
C
	DO I=1,NMATCH
	  I2=ITAGS(I)
	  XY_MATCH(1,I)=SPOTS(1,I2)
	  XY_MATCH(2,I)=SPOTS(2,I2)
C Recalculate HKLs from UB and load integer values into HKL_MATCH
	  CALL CALC_PIXWAVS_TO_HKLS(HKL,UB,PIXWAV(1,I2),1)
	  HKL_MATCH(1,I)=NINT(HKL(1))
	  HKL_MATCH(2,I)=NINT(HKL(2))
	  HKL_MATCH(3,I)=NINT(HKL(3))
	ENDDO
C
C Recalculate x,y,wav using integer HKLs
C
	CALL CALC_HKLS_TO_PIXWAVS(PIXWAV, UB,HKL_MATCH,NMATCH)
C
C Get RMS for x & y differences
C
	SUMX2=0.0
	SUMY2=0.0
	DO I=1,NMATCH
	  DX=PIXWAV(1,I)-XY_MATCH(1,I)
	  DY=PIXWAV(2,I)-XY_MATCH(2,I)
	  SUMX2=SUMX2+DX**2
	  SUMY2=SUMY2+DY**2
	ENDDO
	RMS=SQRT((SUMX2+SUMY2)/NMATCH)
	RMSX=SQRT(SUMX2/NMATCH)
	RMSY=SQRT(SUMY2/NMATCH)
C
C Remove any outliers in x & y, but not if less than 20 pixels
C
	DX=MAX(20.0,RMSX*2.0)
	DY=MAX(20.0,RMSY*2.0)
C
	N=0
	DO I=1,NMATCH
	  IF(ABS(PIXWAV(1,I)-XY_MATCH(1,I)) .GT. DX) CYCLE
	  IF(ABS(PIXWAV(2,I)-XY_MATCH(2,I)) .GT. DY) CYCLE
	  N=N+1
	  XY_MATCH(1,N)=XY_MATCH(1,I)
	  XY_MATCH(2,N)=XY_MATCH(2,I)
	  HKL_MATCH(1,N)=HKL_MATCH(1,I)
	  HKL_MATCH(2,N)=HKL_MATCH(2,I)
	  HKL_MATCH(3,N)=HKL_MATCH(3,I)
	ENDDO
	NMATCH=N
C
C Refine UB rotation & PIX_CEN
C
	IF(NMATCH .GE. 3) THEN
	  CALL REFINE_ROT_UB(UB2, UB)
	  CALL MX3_COPY(UB,UB2)
	  IF( .NOT.LCEN_FIXED ) CALL REFINE_PIX_CEN(UB)
	ENDIF
C
C Calculate RMS of matches
C
	RMS=CALC_RMS_ERR(UB)
C
C Copy the final values back to the arguments
C
	XY_CEN(1)=PIX_CEN(1)
	XY_CEN(2)=PIX_CEN(2)
	NMATCH2=NMATCH
C
C Print out summary of test
C
	PRINT '(3X,A,I3,A,F6.3,A)','Test 4:',NMATCH,' matches, rms=',RMS,' pixels'
C
	RETURN
	END


	INTEGER FUNCTION MAKE_INITIAL_MATCH(UB,HSQ_MAX,DIST_MAX)
C
C Trys to return ~10 obs & calc spot matches
C
	REAL UB(3,3)
C
	REAL HKLS(3,10000),PIXWAV(3,10000)
C
C Calculate x,y,wav for the lowest index hkls
C
	CALL GEN_HKLS_SIMPLE(HKLS,PIXWAV,NHKLS, UB,HSQ_MAX,2.0)
C If < 20 hkls, try again with a smaller wav-min
	IF(NHKLS .LT. 20) THEN
	  CALL GEN_HKLS_SIMPLE(HKLS,PIXWAV,NHKLS, UB,HSQ_MAX,1.3)
	ENDIF
C
C Give up if less than 3 hkls
C
	IF(NHKLS .LT. 3) THEN
	  MAKE_INITIAL_MATCH=0
	  RETURN
	ENDIF
C
C Try to match calculated spots to an increasing number of the
C strongest observed spots. If we have 10 matches (or half the
C calculated spots are matched), jump out of the loop.
C
	MULT=MAX(3,100/NHKLS)
	DO NHKLS2=NHKLS,NHKLS*MULT,NHKLS
	  CALL FIND_MATCH(NMATCH,NOVER, HKLS,PIXWAV,NHKLS, NHKLS2,DIST_MAX)
	  IF(NMATCH .GE. MIN(NHKLS/2,10)) EXIT
	ENDDO
C
	MAKE_INITIAL_MATCH=NMATCH
	RETURN
	END


	SUBROUTINE REFINE_ROT_UB(UB_OUT, UB)
C
	REAL UB(3,3),UB_OUT(3,3)
C
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
C
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	REAL PIXWAV(3),U_ROT(3,3)
	REAL H_HKL(3,10000),H_SPOT(3,10000)
	REAL UB_ROT(3,3),ROT_PHI(3,3)
C
C Create the UB adjusted for the PHI rotation
C
	CALL CALC_YROT_MX(ROT_PHI, -PHI)
	CALL MX3_MULT(UB_ROT, ROT_PHI,UB)
C
C ALIGN_VECTORS can only handle 100 matches
C
	NMATCH2=MIN(100,NMATCH)
C
C Calculate reciprocal (unit) vectors for the hkls
C
	DO I=1,NMATCH2
	  CALL VC3_MULT(H_HKL(1,I), UB_ROT,HKL_MATCH(1,I))
	  CALL VC3_UNIT(H_HKL(1,I))
	ENDDO
C
C Calculate reciprocal (unit) vectors for the found spots
C
	PIXWAV(3)=1.0
	DO I=1,NMATCH2
	  PIXWAV(1)=XY_MATCH(1,I)
	  PIXWAV(2)=XY_MATCH(2,I)
	  CALL CALC_PIXWAVS_TO_HVECS(H_SPOT(1,I), PIXWAV,1)
	  CALL VC3_UNIT(H_SPOT(1,I))
	ENDDO
C
C Rotate UB to align the H_SPOT & H_HKL vectors
C
	CALL ALIGN_VECTORS(U_ROT, H_SPOT,H_HKL,NMATCH2)
C
C Rotate UB (not UB_ROT)
C
	CALL MX3_MULT(UB_OUT, U_ROT,UB)
C
	RETURN
	END



	SUBROUTINE REFINE_XY(XY, UB_IN,XY_IN)
C
	REAL XY(2),UB_IN(3,3),XY_IN(2)
C
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
C
	REAL PIXWAV(3)
C
	SUM_DX=0.0
	SUM_DY=0.0
	DO I=1,NMATCH
	  CALL CALC_HKLS_TO_PIXWAVS(PIXWAV, UB_IN,HKL_MATCH(1,I),1)
	  SUM_DX=SUM_DX+XY_MATCH(1,I)-PIXWAV(1)
	  SUM_DY=SUM_DY+XY_MATCH(2,I)-PIXWAV(2)
	ENDDO
C
C Limit shift in XY to +/-20 pixels
C
	XY(1)=XY_IN(1)+MAX(-20.0,MIN(+20.0, SUM_DX/NMATCH ))
	XY(2)=XY_IN(2)+MAX(-20.0,MIN(+20.0, SUM_DY/NMATCH ))
C
	RETURN
	END



	SUBROUTINE REFINE_PIX_CEN(UB)
C
	REAL UB(3,3)
C
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
C
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
C
	REAL PIXWAV(3)
C
	SUM_DX=0.0
	SUM_DY=0.0
	DO I=1,NMATCH
	  CALL CALC_HKLS_TO_PIXWAVS(PIXWAV, UB,HKL_MATCH(1,I),1)
	  SUM_DX=SUM_DX+XY_MATCH(1,I)-PIXWAV(1)
	  SUM_DY=SUM_DY+XY_MATCH(2,I)-PIXWAV(2)
	ENDDO
	AVE_DX=SUM_DX/MAX(1,NMATCH)
	AVE_DY=SUM_DY/MAX(1,NMATCH)
C
C Limit shift in PIX_CEN to +/-5 pixels
C
	PIX_CEN(1)=PIX_CEN(1) + MAX(-5.0,MIN(+5.0, AVE_DX ))
	PIX_CEN(2)=PIX_CEN(2) + MAX(-5.0,MIN(+5.0, AVE_DY ))
C
	RETURN
	END
