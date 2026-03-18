C ===============================================================================
C	SUBROUTINE MAKE_FULL_MATCH()
C	SUBROUTINE SAVE_MATCH_LIST()
C	SUBROUTINE RESTORE_MATCH_LIST()
C	SUBROUTINE MAKE_PARTIAL_MATCH(RMS,NMATCH,NTRY)
C	INTEGER FUNCTION MAKE_HKL_MATCH(NSTRONG,DIST_MAX)
C	INTEGER FUNCTION MAKE_SPOT_MATCH(NSPOTS2,NCALC,DIST_MAX)
C	SUBROUTINE FIND_MATCH(NMATCH2,NOVER, HKLS,PIXWAV,NHKLS, NSPOTS2,DIST_MAX)
C	FUNCTION CALC_RMS_ERR(UB)
C ===============================================================================


	SUBROUTINE MAKE_FULL_MATCH()
C
C Try to get 50 obs & calc spot matches
C
	INTEGER MAKE_SPOT_MATCH
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
C
	REAL UB(3,3)
C
C Calculate UB matrix from the last refinement
C
	CALL CALC_UB_MATRIX(UB)
C
C Check how many fairly good spots are in the list
C
	NGOOD=2
	DO I=3,NSPOTS
	  IF(SPOTS(3,I) .GT. 2.5) NGOOD =I
	ENDDO
C
C Empirical equations for no. of spots with SIG < CUTOFF
C    NCUTOFF=NGOOD * (2.5/CUTOFF)^POWER
C    where POWER is approximated by (2300-NGOOD)/1300
C If the empirical equations are valid then good targets for
C NOBS & NCALC are: NOBS=23+0.6*NGOOD; NCALC=(15+0.02*NGOOD)^2
C
	NOBS=23.0+0.6*NGOOD
	NCALC=(15.0+0.02*NGOOD)**2
C
	NMATCH=MAKE_SPOT_MATCH(NOBS,NCALC,10.0)
C
C Special case when few matches, try more NOBS
C
	IF(NMATCH .LT. 100) THEN
C Save the current match list
	  CALL SAVE_MATCH_LIST()
C Save the efficiency for matching spots
	  EFF=NMATCH/NOBS
C Try to make larger match lists
	  DO FRAC=1.25,2.51,0.25
	    NTRY=FRAC*NOBS
	    NMATCH2=MAKE_SPOT_MATCH(NTRY,NCALC,10.0)
	    EFF2=(NMATCH2-NMATCH)/(0.25*NOBS)
C If the efficiency significantly drops, give up
	    IF(EFF2 .LT. 0.33*EFF) EXIT
C Seems ok, so overwrite the saved match list
	    CALL SAVE_MATCH_LIST()
	  ENDDO
C Restore the last saved match list
	  CALL RESTORE_MATCH_LIST()
	ENDIF
C
C If still few matches, try the alternative HKL method
C
	IF(NMATCH .LT. 100) THEN
C Save the current match list
	  CALL SAVE_MATCH_LIST()
	  NBEST=NMATCH
C Do HKL method with a small match to get an idea of the fit
	  NMATCH=MAKE_HKL_MATCH(MIN(30,NSPOTS),5.0)
C If better than the original save this match
	  IF(NMATCH .GT. NBEST) CALL SAVE_MATCH_LIST()
	  NBEST=NMATCH		! SAVE_MATCH_LIST updates NMATCH
C Try HKL method with more spots until RMS > 2
	  DO ISTRONG=50,NSPOTS,25
	    NMATCH=MAKE_HKL_MATCH(MIN(ISTRONG,NSPOTS),5.0)
	    RMS=CALC_RMS_ERR(UB)
	    IF(RMS .GT. 2.0) EXIT
	    IF(NMATCH .GT. NBEST) CALL SAVE_MATCH_LIST()
	    NBEST=NMATCH		! SAVE_MATCH_LIST updates NMATCH
	  ENDDO
C Restore the last saved match list
	  CALL RESTORE_MATCH_LIST()
	ENDIF
C
C Calculate fit and output results
C
	RMS=CALC_RMS_ERR(UB)
	PRINT '(1X,A,I5,A,F6.3,A)','Matching',NMATCH,' spots, initial rms=',RMS,' mm'
C
	RETURN
	END



	SUBROUTINE SAVE_MATCH_LIST()
C
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
	COMMON /MATCH_SAVE/ XY_SAVE(2,10000),HKL_SAVE(3,10000),NSAVE
C
	DO I=1,NMATCH
	  XY_SAVE(1,I)=XY_MATCH(1,I)
	  XY_SAVE(2,I)=XY_MATCH(2,I)
	  HKL_SAVE(1,I)=HKL_MATCH(1,I)
	  HKL_SAVE(2,I)=HKL_MATCH(2,I)
	  HKL_SAVE(3,I)=HKL_MATCH(3,I)
	ENDDO
	NSAVE=NMATCH
C
	RETURN
	END



	SUBROUTINE RESTORE_MATCH_LIST()
C
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
	COMMON /MATCH_SAVE/ XY_SAVE(2,10000),HKL_SAVE(3,10000),NSAVE
C
	NMATCH=NSAVE
	DO I=1,NMATCH
	  XY_MATCH(1,I)=XY_SAVE(1,I)
	  XY_MATCH(2,I)=XY_SAVE(2,I)
	  HKL_MATCH(1,I)=HKL_SAVE(1,I)
	  HKL_MATCH(2,I)=HKL_SAVE(2,I)
	  HKL_MATCH(3,I)=HKL_SAVE(3,I)
	ENDDO
C
	RETURN
	END


	SUBROUTINE MAKE_PARTIAL_MATCH(RMS,NMATCH,NTRY)
C
C Try to get 50 obs & calc spot matches
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS
C
	REAL UB(3,3)
C
	CALL CALC_UB_MATRIX(UB)
C
C If less than 100 spots, try to match all of them
C
	IF(NSPOTS .LT. 100) THEN
	  NMATCH=MAKE_SPOT_MATCH(NSPOTS,100,10.0)
	  RMS=CALC_RMS_ERR(UB)
C
C Try for 50 matches from the strongest 10% of found spots,
C if it fails then the strongest 20%,...,50% of spots.
	ELSE
	  DO FRAC=0.1,0.51,0.1
	    NTRY=INT(NSPOTS*FRAC)
	    NMATCH=MAKE_SPOT_MATCH(NTRY,100,10.0)
	    IF(NMATCH .GE. 50) EXIT
C Give up if spot merit below 3.5
	    IF(SPOTS(3,NTRY) .LT. 3.5) EXIT
	  ENDDO
	ENDIF
C
C If very few matches, try the alternative HKL method
C
	IF(NMATCH .LT. 10) NMATCH=MAKE_HKL_MATCH(MIN(30,NSPOTS),3.0)
	IF(NMATCH .LT. 10) NMATCH=MAKE_HKL_MATCH(MIN(50,NSPOTS),5.0)
C
	RMS=CALC_RMS_ERR(UB)
	PRINT '(1X,A,I5,A,F6.3,A)','Matching',NMATCH,' spots, initial rms=',RMS,' mm'
	RETURN
	END



	INTEGER FUNCTION MAKE_HKL_MATCH(NSTRONG,DIST_MAX)
C
C Try to match the first NSTRONG observed spots by calculating
C integer HKLs from X,Y and keeping those within DIST_MAX pixels
C of the observed.
C NB: This routine is probably only for triclinic, so I ignore centering
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS
C
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
C
	REAL UB(3,3),DIFF_BEST(10000),HKL(3),RHKL(3),PIXWAV(3),TEMP(3)
C
C Estimate WAV_MIN(=2*D_MIN) that gives NCALC spots on the detector
C
	CALL CALC_UB_MATRIX(UB)
C
C Calculate best integer HKLs for observed spots
C
	WAV_MIN=1.0
	DO I=1,NSTRONG
C Calculate HKL assuming the minimum wavelength possible
	  PIXWAV(1)=SPOTS(1,I)
	  PIXWAV(2)=SPOTS(2,I)
	  PIXWAV(3)=WAV_MIN
	  CALL CALC_PIXWAVS_TO_HKLS(HKL, UB,PIXWAV,1)
C Find best integer HKL with same ratios but small indices
	  RMAX=MAX(ABS(HKL(1)),ABS(HKL(2)),ABS(HKL(3)))
	  DO I2=1,3
	    RHKL(I2)=HKL(I2)/RMAX
	  ENDDO
C Record scaled HKL with best Obs-Calc distance
	  DIFF_BEST(I)=1E6
	  DO ITRY=1,IFIX(RMAX)
	    DO I2=1,3
	      HKL(I2)=NINT( RHKL(I2)*ITRY )
	    ENDDO
	    CALL CALC_HKLS_TO_PIXWAVS(TEMP, UB,HKL,1)
	    DIFF=SQRT( (PIXWAV(1)-TEMP(1))**2 +  (PIXWAV(2)-TEMP(2))**2 )
	    IF(DIFF .LT. DIFF_BEST(I)) THEN
	      DIFF_BEST(I)=DIFF
	      DO I2=1,3
	        HKL_MATCH(I2,I)=HKL(I2)
	      ENDDO
	    ENDIF
	  ENDDO
C Copy best integer HKL into HKL_MATCH()
	ENDDO
C
C Reject any matches worse than DIST_MAX
C
	N=0
	DO I=1,NSTRONG
	  IF(DIFF_BEST(I) .GT. DIST_MAX**2) CYCLE
	  N=N+1
	  XY_MATCH(1,N)=SPOTS(1,I)
	  XY_MATCH(2,N)=SPOTS(2,I)
	  HKL_MATCH(1,N)=HKL_MATCH(1,I)
	  HKL_MATCH(2,N)=HKL_MATCH(2,I)
	  HKL_MATCH(3,N)=HKL_MATCH(3,I)
	ENDDO
C
	MAKE_HKL_MATCH=N
	RETURN
	END



	INTEGER FUNCTION MAKE_SPOT_MATCH(NSPOTS2,NCALC,DIST_MAX)
C
C Try to match the first NSPOTS2 observed spots with NCALC
C generated HKLs within a distance of DIST_MAX pixels
C
	COMMON /CELL_COM/ CELL(6),ILATT,ICEN,CELL_PROD
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS
C
	REAL UB(3,3),DUMMY(3,3)
	REAL HKLS(3,10000),PIXWAV(3,10000)
C
C Estimate WAV_MIN(=2*D_MIN) that gives NCALC spots on the detector
C
	CALL CALC_UB_MATRIX(UB)
	CALL MX3_INVERT(DUMMY,DET_UB, UB)
	VCELL=1.0/MAX(1E-6,ABS(DET_UB))
      WAV_MIN = ( 1.75*VCELL/NCALC )**0.3333
C
C Do correction for lattice centering
C
	IF(ICEN.EQ. 5) THEN
	  WAV_MIN=WAV_MIN * 0.63	! FCC
	ELSE IF(ICEN .EQ. 6) THEN
	  WAV_MIN=WAV_MIN * 0.69	! R (rhom in hex)
	ELSE IF(ICEN .NE. 0) THEN
	  WAV_MIN=WAV_MIN * 0.79	! A,B,C,BCC
	ENDIF
C
C Generate hkl & pixwav lists (limit WAV_MIN to > 0.9)
C
	CALL GEN_HKLS_LSQ(HKLS,PIXWAV,NHKLS, MAX(0.9,WAV_MIN),WAV_MIN/2)
C
C Match spots using the NOBS strongest spots
C
	CALL FIND_MATCH(NMATCH,NOVER, HKLS,PIXWAV,NHKLS, NSPOTS2,DIST_MAX)
C
	MAKE_SPOT_MATCH=NMATCH
	RETURN
	END



	SUBROUTINE FIND_MATCH(NMATCH2,NOVER, HKLS,PIXWAV,NHKLS, NSPOTS2,DIST_MAX)
C
	REAL PIXWAV(3,10000),HKLS(3,10000)
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
C
C Loop over calc. spots adding to the list of MATCH
C
	NMATCH=0
	NOVER=0
	DO ICALC=1,NHKLS
	  NFIND=0
	  DO IOBS=1,MIN(NSPOTS,NSPOTS2)
C Find observed spots within DIST_MAX pixels of the ICALC spot
	    DSQ=(PIXWAV(1,ICALC)-SPOTS(1,IOBS))**2 +
	1			(PIXWAV(2,ICALC)-SPOTS(2,IOBS))**2
	    IF(DSQ < DIST_MAX**2) THEN
		    IFIND=IOBS
	      NFIND=NFIND+1
	    ENDIF
	  ENDDO
C If a unique match, add it to	MATCH
	  IF(NFIND .EQ. 1) THEN
	    NMATCH=NMATCH+1
	    XY_MATCH(1,NMATCH)=SPOTS(1,IFIND)
	    XY_MATCH(2,NMATCH)=SPOTS(2,IFIND)
	    HKL_MATCH(1,NMATCH)=HKLS(1,ICALC)
	    HKL_MATCH(2,NMATCH)=HKLS(2,ICALC)
	    HKL_MATCH(3,NMATCH)=HKLS(3,ICALC)
C If more than one match, increment the overlap counter
	  ELSE IF(NFIND .GT. 1) THEN
	    NOVER=NOVER+1
	  ENDIF
	ENDDO
C
	NMATCH2=NMATCH
	RETURN
	END


	FUNCTION CALC_RMS_ERR(UB)
C
	REAL UB(3,3)
C
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
C
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
C
	REAL PIXWAV(3,10000)
C
	IF(NMATCH .EQ. 0) THEN
	  CALC_RMS_ERR=0.0
	  RETURN
	ENDIF
C
	CALL CALC_HKLS_TO_PIXWAVS(PIXWAV, UB,HKL_MATCH,NMATCH)
C
	SUMX2=0.0
	SUMY2=0.0
	DO I=1,NMATCH
	  SUMX2=SUMX2 + (XY_MATCH(1,I)-PIXWAV(1,I))**2
	  SUMY2=SUMY2 + (XY_MATCH(2,I)-PIXWAV(2,I))**2
	ENDDO
C
	SSQ=SUMX2*PIX_SIZE(1)**2 + SUMY2*PIX_SIZE(2)**2
	CALC_RMS_ERR=SQRT(SSQ/NMATCH)
	RETURN
	END
