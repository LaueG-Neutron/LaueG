C----- Manage the non-parametric smoothing of wavelength spectra --------------
c	SUBROUTINE CLEAR_WAV_NONPAR
c	FUNCTION CALC_WAV_NONPAR(WAV0)
c	SUBROUTINE FIT_WAV_NONPAR
c	SUBROUTINE LOAD_WAV_NONPAR(NBINS_NPAR,GAMMA)
c	SUBROUTINE MERGE_WAV_RATIOS(RAT,SRAT,RAT2, WAV_MIN,WAV_FRAC, NWAVS)
c	SUBROUTINE CALC_TAUT_SMOOTH(SMOOTH,TRACE, OBS,SIG,NDATA, GAMMA)
c	SUBROUTINE SOLVE_5DIAG_EQN_TRACE(TRACE, A_BAND,NDATA)
c	SUBROUTINE SOLVE_5DIAG_EQN(X_VEC, A_MAT,Y_VEC,NDATA)
C------------------------------------------------------------------------------


C----- Manage the non-parametric smoothing of wavelength spectra --------------

	SUBROUTINE CLEAR_WAV_NONPAR
C
C Set NWAV_NPAR=0 to signal CALC_WAV_NONPAR to return a value of 1
C
	PARAMETER (NBINS_MAX=5000)
	COMMON /WAV_NPAR_COM/ NWAV_NPAR,WAV_NPAR_MIN,WAV_NPAR_STEP,DOF,RAT(NBINS_MAX)
C
	NWAV_NPAR=0
	RETURN
	END


	FUNCTION CALC_WAV_NONPAR(WAV0)
C
C Calculate the smoothed ratio IOBS/ICALC and its esd for a given wavelength
C
	PARAMETER (NBINS_MAX=5000)
	COMMON /WAV_NPAR_COM/ NWAV_NPAR,WAV_NPAR_MIN,WAV_NPAR_STEP,DOF,RAT(NBINS_MAX)
C
C NWAV_NPAR=0 flags that the non-parameteric correction has not been "fitted".
C In this case just return 1.
C
	IF(NWAV_NPAR .EQ. 0) THEN
	  CALC_WAV_NONPAR=1.0
	  RETURN
	ENDIF
C
C Calculate the correction via linear interpolation of the table values
C
	CALC_WAV_NONPAR=CALC_LIN_SPLINE(WAV0-WAV_NPAR_MIN, WAV_NPAR_STEP,RAT,NWAV_NPAR)
C
      RETURN
	END


	SUBROUTINE FIT_WAV_NONPAR
C
C Use taut splines to "fit" the non-param distribution, if the option is selected.
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	COMMON /WAV_CORR_COM/ IWAV_OPT,NWAV,WAV_GAMMA,WAV_MIN,WAV_STEP,WAV_PAR(20)
C
	PARAMETER (NBINS_MAX=5000)
	COMMON /WAV_NPAR_COM/ NWAV_NPAR,WAV_NPAR_MIN,WAV_NPAR_STEP,DOF,RAT(NBINS_MAX)
C
C Set DOF and return if not wavelength correction with non-parametric mode
C
	DOF=0.0
	IF( .NOT.(LWAV_CORR .AND. (IWAV_OPT .EQ. 2)) ) RETURN
C
C Output statistics before the correction for comparison
C
	PRINT '(1X,A,/)','Calculating non-parametric wavelength correction'
	WRITE(10,'(/,A,/)') '====== Non-parametric wavelength correction ======'
C
	WRITE(10,'(A)') 'Merge statistics before correction:'
	CALL CALC_ALL_COUNTS_SEQ
	CALL LOG_LSQ_MERGE
C
C Load tables of the ratio of IOBS/ICALC
C Recalculate all ICALC with new wavelength distribution by "refining"
C the COUNTS_SEQ() parameters, all other LSQ parameters are fixed.
C
	CALL LOAD_WAV_NONPAR(NBINS_NPAR,WAV_GAMMA)
C
	WRITE(10,'(/,A,F6.1,A,I4,A)') 'Smoothing algorithm for gamma =',
	1					WAV_GAMMA,' used',NINT(DOF),' degrees of freedom'
C
	WRITE(10,'(/,A,I2)') 'Merge statistics after correction:'
	CALL CALC_ALL_COUNTS_SEQ
	CALL LOG_LSQ_MERGE
C
	RETURN
	END


	SUBROUTINE LOAD_WAV_NONPAR(NBINS_NPAR,GAMMA)
C
C Smooth the ratio IOBS/ICALC in terms of wavelength using a non-parameteric
C smoother. Then create tables of the smoothed ratio suitable for linear
C interpolation by the routine CALC_WAV_NONPAR.
C The degrees of freedom used by the smoothing algorithm saved in /WAV_NPAR_COM/
C
	PARAMETER (NBINS_MAX=5000)
	COMMON /WAV_NPAR_COM/ NWAV_NPAR,WAV_NPAR_MIN,WAV_NPAR_STEP,DOF,RAT(NBINS_MAX)
C
	REAL SMOOTH(NBINS_MAX),CURRENT(NBINS_MAX),RAT2(NBINS_MAX)
	REAL TRACE(NBINS_MAX),RAT_SQ(NBINS_MAX),SRAT(NBINS_MAX)
C
C Sanity checks
C
	IF(NWAV_NPAR .NE. 0) STOP 'BUG(load_wav_nonpar): NWAV_NPAR != 0'
	IF(NBINS_NPAR+20 .GT. NBINS_MAX) STOP 'BUG(load_wav_nonpar): NBINS_NPAR too large'
C
C Merge the ratio IOBS/ICALC into NWAV_NPAR wavelength bins which are in
C the geometric sequence WAV_MIN*WAV_FRAC**(IBIN-1). We use logarithmically
C spaced bins as peaks/dips are narrower at shorter wavelengths.
C
	CALL MERGE_WAV_RATIOS(RAT,SRAT,RAT_SQ, WAV_MIN,WAV_FRAC, NBINS_NPAR)
C
C Extend bins by 10 points at both ends using fake ratios of 1.
C This helps the ends after smoothing from "curling up".
C
	DO I=NBINS_NPAR,1,-1
	  RAT   (I+10)=RAT   (I)
	  RAT_SQ(I+10)=RAT_SQ(I)
	  SRAT  (I+10)=SRAT  (I)
      ENDDO
C
      NBINS=NBINS_NPAR+20
	WAV_MIN=WAV_MIN*WAV_FRAC**(-10)
C
      DO I=1,10
	  RAT   (I)=1.0
	  RAT_SQ(I)=0.0
	  SRAT  (I)=0.1
	  RAT   (NBINS-I+1)=1.0
	  RAT_SQ(NBINS-I+1)=0.0
	  SRAT  (NBINS-I+1)=0.1
      ENDDO
C
C We have to smooth the total wavelength correction, which includes
C the current correction from the calibration file and the LSQ spline
C correction, and the IOBS/ICALC in RAT(). We can't do this directly.
C Instead we modify RAT() as if ICALC() was calculated using a smoothed
C version of the calibration and spline correction. We then smooth this
C new ratio and multiply it by the smoothed version of the calibration
C and spline correction to get a smoothed version of the total
C correction.
C
C Calculate the current wavelength correction (calibration file and
C LSQ spline) for each bin.
C
      WAV=WAV_MIN
      DO I=1,NBINS
	  CURRENT(I)=CALC_WAV_SCALE(WAV)
        WAV=WAV*WAV_FRAC
      ENDDO
C
C Smooth the correction (strongly oversmooth using 10% fractional errors)
C
	CALL CALC_TAUT_SMOOTH(SMOOTH,TRACE, CURRENT,CURRENT,NBINS, GAMMA*0.1)
C
C Adjust RAT to be relative to the smoothed wavelength correction, then
C subtract 1 from RAT() so nominal value is 0 instead of 1
C
	DO I=1,NBINS
	  RAT(I)=RAT(I) * CURRENT(I) / SMOOTH(I) - 1.0
	ENDDO
C
C Smooth the ratio values and put in RAT2()
C
	CALL CALC_TAUT_SMOOTH(RAT2,TRACE, RAT,SRAT,NBINS, GAMMA)
C
C Copy the smoothed value for the last "real" data points at
C both ends over the extended parts of RAT()
C
C
      DO I=1,10
	  RAT(I)=RAT2(11)
	  RAT(NBINS-I+1)=RAT2(NBINS-10)
	ENDDO
C
C Smooth the ratios again, now with the smoothed end points
C
	CALL CALC_TAUT_SMOOTH(RAT2,TRACE, RAT,SRAT,NBINS, GAMMA)
C
C Sum the diagonal to get "degrees of freedom"
C
	DOF=0.0
	DO I=1,NBINS
	  DOF=DOF+TRACE(I)
      ENDDO
C
C Convert RAT2() to a ratio relative to the current wavelength correction
C
	DO I=1,NBINS
	  RAT2(I)=(RAT2(I)+1.0) * SMOOTH(I) / CURRENT(I)
      ENDDO
C
C Change from the logarithmic binning of RAT2() to constant 0.001
C binning for RAT(). Use the same wavelength min and max.
C
      WAV_NPAR_MIN=WAV_MIN
      WAV_MAX=WAV_NPAR_MIN*WAV_FRAC**(NBINS-1)
	WAV_NPAR_STEP=0.001
	NWAV_NPAR=1+NINT( (WAV_MAX-WAV_NPAR_MIN) / WAV_NPAR_STEP )
C
	WAV=WAV_NPAR_MIN
	DO I=1,NWAV_NPAR
	  X=ALOG(WAV/WAV_NPAR_MIN)/ALOG(WAV_FRAC)
	  RAT(I)=CALC_LIN_SPLINE(X, 1.0,RAT2,NBINS)
	  WAV=WAV+WAV_NPAR_STEP
	ENDDO	
C
	RETURN
	END


	SUBROUTINE MERGE_WAV_RATIOS(RAT,SRAT,RAT2, WAV_MIN,WAV_FRAC, NWAVS)
C
C Perform weighted averages of IOBS(I)/ICALC(I) for a set of NWAVS bins in wavelength.
C The binned ratio, its esds, and the ratio-squared are stored in RAT(), SRAT() & RAT2().
C The wavelength of the bins is a geometric sequence [WAV_MIN, WAV_MIN*WAV_FRAC, ...]
C
	REAL RAT(NWAVS),SRAT(NWAVS),RAT2(NWAVS)
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
C Find the wavelength limits of the data, and calculate WAV_FRAC from NWAVS
C
	WAV_MIN=WAVS(1)
	WAV_MAX=WAVS(1)
	DO I=1,ILAST(NSEQ)
	  WAV_MIN=MIN(WAV_MIN,WAVS(I))
	  WAV_MAX=MAX(WAV_MAX,WAVS(I))
	ENDDO
	WAV_FRAC=EXP( ALOG(WAV_MAX/WAV_MIN)/(NWAVS-1) )
C
C Do the sums over wavelength bins
C
	DO I=1,NWAVS
	  RAT(I)=0.0
	  RAT2(I)=0.0
	  SRAT(I)=0.0
	ENDDO
C
	DO I=1,NDATA
C Exclude really poor data
	  IF(ICALC(I).GT.0.0 .AND. ICALC(I).GT.ISIG(I)) THEN
	    IWAV=1 + NINT( ALOG(WAVS(I)/WAV_MIN) / ALOG(WAV_FRAC) )
	    IWAV=MAX(1,MIN(NWAVS, IWAV ))
	    W=( ICALC(I)/ISIG(I) )**2
	    SRAT(IWAV)=SRAT(IWAV)+W
	    RAT(IWAV)=RAT(IWAV)+ W*IOBS(I)/ICALC(I)
	    RAT2(IWAV)=RAT2(IWAV)+ W*(IOBS(I)/ICALC(I))**2
	  ENDIF
	ENDDO
C
C Convert the sums to weighted averages and esds
C
	DO I=1,NWAVS
	  IF(SRAT(I) .EQ. 0.0) THEN
C If no data, assume the default ratio with a 10% uncertainty
	    RAT(I)=1.0
	    RAT2(I)=0.0
	    SRAT(I)=0.1
	  ELSE
C Put average of obs/calc into RAT() and its esd into SRAT()
	    RAT(I)=RAT(I)/SRAT(I)
	    RAT2(I)=RAT2(I)/SRAT(I)
	    SRAT(I)=1.0/SQRT(SRAT(I))
	  ENDIF
	ENDDO
C
C
	RETURN
	END


	SUBROUTINE CALC_TAUT_SMOOTH(SMOOTH,TRACE, OBS,SIG,NDATA, GAMMA)
C
C Create a smoothed version, SMOOTH, of data OBS with esds SIG.
C The smoothing penalises the sum of the square of second-derivatives
C of SMOOTH, the Lagrange multiplier being GAMMA**2.
C Also return the diagonal of the smoothing matrix in TRACE.
C Data is assumed to be equispaced and of size NDATA.
C
	REAL SMOOTH(NDATA),TRACE(NDATA),OBS(NDATA),SIG(NDATA)
C
	PARAMETER (NBINS_MAX=5000)
	REAL A_BAND(-2:2,NBINS_MAX)
C
	IF(NDATA .GT. NBINS_MAX) STOP 'BUG(calc_taut_smooth): NDATA too large'
C
C Create the 5 diagonal banded matrix of the problem
C
	DO I=3,NDATA-2
	  A_BAND(+2,I-2)=     (GAMMA*SIG(I))**2
	  A_BAND(+1,I-1)=-4.0*(GAMMA*SIG(I))**2
	  A_BAND( 0,I  )= 6.0*(GAMMA*SIG(I))**2 + 1.0
	  A_BAND(-1,I+1)=-4.0*(GAMMA*SIG(I))**2
	  A_BAND(-2,I+2)=     (GAMMA*SIG(I))**2
	ENDDO
C For last 2 end-points turn off smoothing
	DO I1=1,2
	  DO I2=-2,2
	    IF(I1 .GT. I2) A_BAND(I2,I1-I2)=0.0
	    IF(I1 .GT. -I2) A_BAND(I2,NDATA+1-I1-I2)=0.0
	  ENDDO
	  A_BAND(0,I1)=1.0
	  A_BAND(0,NDATA+1-I1)=1.0
	ENDDO
C
C Solve the linear equations OBS = A_BAND * SMOOTH
C
	CALL SOLVE_5DIAG_EQN_TRACE(TRACE, A_BAND,NDATA)
	CALL SOLVE_5DIAG_EQN(SMOOTH, A_BAND,OBS,NDATA)
C
	RETURN
	END


	SUBROUTINE SOLVE_5DIAG_EQN_TRACE(TRACE, A_BAND,NDATA)
C
C Calculates the degrees of freedom from the Trace of A_BAND-inverse.
C NB: Assumes 10 dummy points have been added to both ends of the data
C
	REAL A_BAND(-2:2,NDATA),TRACE(NDATA)
C
	PARAMETER (NBINS_MAX=5000)
	INTEGER IPIVOTS(NBINS_MAX)
	REAL*8 B_VEC(NBINS_MAX,NBINS_MAX)
	REAL*8 A_MAT(NBINS_MAX,NBINS_MAX)
C
	IF(NDATA .GT. NBINS_MAX) STOP 'BUG(solve_5diag_eqn): NDATA too large'
C
C Setup the matrices as DGBSV wants it (NB: change to double precision)
C
	DO I1=1,NDATA
	  DO I2=1,NDATA
	    A_MAT(I2,I1)=0.0
	    IF(ABS(I2-I1) .LE. 2) A_MAT(I2,I1)=A_BAND(I2-I1,I1)
	    B_VEC(I1,I2)=0.0
	  ENDDO
	  B_VEC(I1,I1)=1.0
	ENDDO
C
C Solve the linear equations and complain and die if it unexpectedly fails
C
	NB=NDATA
	CALL DGESV(NDATA, NB, A_MAT, NBINS_MAX, IPIVOTS, B_VEC, NBINS_MAX, ISTATUS)
	IF(ISTATUS .NE. 0) STOP	'BUG(solve_5diag_eqn): Failed'
C
C Copy the diagonal values into TRACE()
C
	DO I=11,NDATA-10
	  TRACE(I)=B_VEC(I,I)
	ENDDO
C
	RETURN
	END


	SUBROUTINE SOLVE_5DIAG_EQN(X_VEC, A_MAT,Y_VEC,NDATA)
C
C Solves X_VEC for the linear equation A_MAT * X_VEC = Y_VEC
C where A_MAT is a 5-diagonal matrix of order NDATA.
C Uses LAPACK to solve the linear equations
C
	REAL X_VEC(NDATA),Y_VEC(NDATA)
	REAL A_MAT(-2:2,NDATA)
C
	PARAMETER (NBINS_MAX=5000)
	INTEGER IPIVOTS(NBINS_MAX)
	REAL*8 A_BAND(7,NBINS_MAX),B_VEC(NBINS_MAX,1)
C
	IF(NDATA .GT. NBINS_MAX) STOP 'BUG(solve_5diag_eqn): NDATA too large'
C
C Setup the matrices as DGBSV wants it (NB: change to double precision)
C
	DO I=1,NDATA
	  DO IDIAG=MAX(-2,1-I),MIN(2,NDATA-I)
	    A_BAND(IDIAG+5,I)=A_MAT(IDIAG,I)
	  ENDDO
	  B_VEC(I,1)=Y_VEC(I)
	ENDDO
C
C Solve the linear equations and complain and die if it unexpectedly fails
C
	CALL DGBSV(NDATA, 2,2,1, A_BAND,7, IPIVOTS, B_VEC, NBINS_MAX, ISTATUS)
	IF(ISTATUS .NE. 0) STOP	'BUG(solve_5diag_eqn): Failed'
C
C Copy the result back to the single precision vector
C
	DO I=1,NDATA
	  X_VEC(I)=B_VEC(I,1)
	ENDDO
C
	RETURN
	END
