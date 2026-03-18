C----- Refine the data set for various scenarios -------------------------
c	SUBROUTINE PRELIM_REFINE
c	SUBROUTINE FULL_REFINE
c	SUBROUTINE FINAL_REFINE
c	SUBROUTINE INITIAL_FILE_SCALES
C-------------------------------------------------------------------------


C----- Refine the data set for various scenarios -------------------------

	SUBROUTINE PRELIM_REFINE
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	PRINT '(/,1X,A)','Preliminary Parameter Refinement'
	WRITE(10,'(/,A,/)') '====== Preliminary parameter refinement ======'
C
C Load COUNTS & DCOUNTS for refinement
C
	CALL LOAD_LSQ_DATA
C
C Increase esds using empirical formulae to emphasize spots with wav ~ 1.25
C
	CALL ADD_ISIG_INITIAL
C
C Load variables and other initialisation tasks for the parameters
C NB: Don't initialise EXTI value until COUNTS_SEQ is calculated
C
	CALL INIT_FILE_PARS
	CALL INIT_TWINS_PARS
	CALL INIT_WAV_PARS
	CALL INIT_EFF_PARS
	CALL INIT_ABS_PARS
C
C Do a fast estimate of scale factors using statistics instead of LSQ
C
	IF( LFILE_CORR ) CALL INITIAL_FILE_SCALES
C
C Calculate COUNTS_SEQ (always with extinction off, i.e. EXTI=0)
C
	CALL CALC_ALL_COUNTS_SEQ
C
C Initialise EXTI value now that COUNTS_SEQ is calculated
C
	CALL INIT_EXTI_PARS
C
C Do initial refinement
C
C Set refinement for file, exti, absorb, and possibly twin ratio
	CALL REFINE_FILE_PARS(LFILE_CORR)
cccccccccccccccccccccccc	CALL REFINE_TWIN_PARS
	CALL REFINE_WAV_PARS(.FALSE.)
	CALL REFINE_EFF_PARS(.FALSE.)
	CALL REFINE_EXTI_PARS(LEXTI_CORR)
	CALL REFINE_ABS_PARS(LABS_CORR)
C Run the LSQ refinement and output summary of corrections
	CALL RUN_LSQ
C
C
	RETURN
	END



	SUBROUTINE FULL_REFINE
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	PRINT *,'Full Parameter Refinement'
	WRITE(10,'(/,A,/)') '====== Full refinement of parameters ======'
C
C Re-load COUNTS & DCOUNTS to wipe out preliminary weighting scheme
C
	CALL LOAD_LSQ_DATA
C
C Modify refinement weights for various model uncertainties
C
	CALL ADD_ISIG_FINAL
C
C Refine all parameters as given in EDIT_OPTIONS, and possibly twin ratio
C
	CALL REFINE_FILE_PARS(LFILE_CORR)
	CALL REFINE_TWIN_PARS
	CALL REFINE_WAV_PARS (LWAV_CORR)
	CALL REFINE_EFF_PARS (LEFF_CORR)
	CALL REFINE_ABS_PARS (LABS_CORR)
	CALL REFINE_EXTI_PARS(LEXTI_CORR)
C
C Run the actual LSQ refinement and give some refinement information
C
	CALL RUN_LSQ
C
	RETURN
	END


	SUBROUTINE FINAL_REFINE
C
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
	PARAMETER (NPARS_MAX=500)
	COMMON /LSQ_PARS_COM/ PARS(NPARS_MAX),ESDS(NPARS_MAX),COVAR(NPARS_MAX,NPARS_MAX),NPARS
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	COMMON /EXTI_CORR_COM/ IEXTI_OPT,EXTI
C
C Only need DOF, the degrees of freedom used by the FIT_WAV_NONPAR smoother
	PARAMETER (NBINS_MAX=5000)
	COMMON /WAV_NPAR_COM/ NWAV_NPAR,WAV_NPAR_MIN,WAV_NPAR_STEP,DOF,RAT(NBINS_MAX)
C
	PRINT *,'Final Parameter Refinement'
	WRITE(10,'(/,A,/)') '====== Final refinement of parameters ======'
C
C Re-load COUNTS & DCOUNTS to wipe out initial weighting scheme
	CALL LOAD_LSQ_DATA
C Output merge statistics before we fiddle with the esds
	WRITE(10,'(A)') 'Merge statistics using initial weights'
 	CALL CALC_ALL_COUNTS_SEQ
	CALL LOG_LSQ_MERGE
C
C Modify refinement weights for various model uncertainties
C
	CALL ADD_ISIG_FINAL
C
C Output merge statistics for final weights
C
	WRITE(10,'(A)') 'Merge statistics using final weights:'
	CALL CALC_ALL_COUNTS_SEQ
	CALL LOG_LSQ_MERGE
C
C Run the LSQ refinement according to parameters in FLAG_OPTS_COM, and possibly twin ratio
C
	CALL REFINE_FILE_PARS(LFILE_CORR)
	CALL REFINE_TWIN_PARS
	CALL REFINE_WAV_PARS (LWAV_CORR)
	CALL REFINE_EFF_PARS (LEFF_CORR)
	CALL REFINE_ABS_PARS (LABS_CORR)
	CALL REFINE_EXTI_PARS(LEXTI_CORR)
	CALL RUN_LSQ
C
C Output merge statistics to the console
C
	PRINT '(1X,A)','Merge statistics for refinement intensities:'
	CALL PRINT_LSQ_MERGE
C
C Output merge statistics then summary of corrections to the log file
C
	WRITE(10,'(/,A)') 'Merge statistics for final refinement:'
	CALL LOG_LSQ_MERGE
	WRITE(10,'(/,A)') 'Summary of refineable correction factors:'
	CALL LOG_CORRECT_SUMMARY
C
C Complain and die if extinction parameter is negative
C
	IF( LEXTI_CORR .AND. EXTI.LT.0.0) THEN
	  WRITE(10,*) 'ERROR: Negative extinction parameter detected'
	  PRINT '(/,2(1X,A,/))','ERROR: Negative extinction parameter detected',
	1				'       Disable extinction correction and re-run'
	  STOP
	ENDIF
C
C Output info on number of data and degrees of freedom to log file
C
	NTOTAL=NDATA-NPARS-NINT(DOF)-NSEQ
	WRITE(10,'(/,A,/,(A,I5))')
	1		'Degrees of freedom:',
	2		'  Number of observations        ',NDATA,
	2		'  Least squares parameters      ',-NPARS,
	3		'  Wavelength non-param smoother ',-NINT(DOF),
	4		'  Merged intensities            ',-NSEQ,
	5		'Remaining degrees of freedom =  ',NTOTAL 
C
	RETURN
	END



	SUBROUTINE INITIAL_FILE_SCALES
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_CORR_COM/ FSCALE(NFILES_MAX)
C
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	REAL SUM(NFILES_MAX),SUMW(NFILES_MAX),SCALES(NFILES_MAX)
C
	WRITE(10,'(A)') 'Calculating initial image file scale factors'
C
C Calculate all COUNTS_SEQ() & IOBS()
C
	CALL CALC_ALL_COUNTS_SEQ
C
C Zero the arrays used for weighted average of file scale factors
C
	DO I=1,NFILES
	  SUM(I)=0.0
	  SUMW(I)=0.0
	ENDDO
C
C Do weighted sum (for > 3 esd spots) of ICALC()/IOBS()
C
	DO ISEQ=1,NSEQ
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    WEIGHT=MAX(0.0,IOBS(I))/ISIG(I)
	    SCALE=ICALC(I)/MAX(0.1,IOBS(I))
	    IF(WEIGHT .GT. 3.0) THEN
	      SUM(IFILES(I))=SUM(IFILES(I))+SCALE*WEIGHT
	      SUMW(IFILES(I))=SUMW(IFILES(I))+WEIGHT
	    ENDIF
	  ENDDO
	ENDDO
C
C Put weighted average of ICALC()/IOBS() into SCALES()
C
	DO I=1,NFILES
	  SCALES(I)=MAX(0.1,SUM(I))/MAX(0.1,SUMW(I))
	ENDDO
C
C Calculate the modal values of SCALES()
C NB: The technique also destroys SCALES()
C
	DO IREJ=1,NFILES/2
C Average all SCALES() > 0   (=0 means deleted)
	  AVE=0.0
	  DO I=1,NFILES
	    IF(SCALES(I) .GT. 0.0) AVE=AVE+SCALES(I)
	  ENDDO
	  AVE=AVE/(NFILES-IREJ+1)
C Find SCALES() furthest from AVE and delete it (value = 0)
	  DIFF=-1.0
	  DO I=1,NFILES
	    IF(SCALES(I) .GT. 0.0) THEN
		  IF(ABS(SCALES(I)-AVE) .GT. DIFF) THEN
	        DIFF=ABS(SCALES(I)-AVE)
	        IDIFF=I
	      ENDIF
	    ENDIF
	  ENDDO
	  SCALES(IDIFF)=0.0
	ENDDO
	RMODAL=AVE
C
C Put weighted average of ICALC()/IOBS() into SCALES()
C Then divide RMODAL by SCALES() and the value in the refinable scale-factors
C
	DO I=1,NFILES
	  SCALES(I)=MAX(0.1,SUM(I))/MAX(0.1,SUMW(I))
	  FSCALE(I)=RMODAL/SCALES(I)
	ENDDO
C
	RETURN
	END
