C----- Adjust LSQ weighting (using ISIG) for initial/final/correct runs -------
c	SUBROUTINE ADD_ISIG_INITIAL
c	SUBROUTINE ADD_ISIG_FINAL
c	SUBROUTINE ADD_ISIG_CORRECT
C----- Components used in above routines --------------------------------------
c	SUBROUTINE ADD_MULT_ISIG
c	SUBROUTINE ADD_ISIG_REL
c	SUBROUTINE ADD_ISIG_EXTI(EXTI_ERR)
c	SUBROUTINE ADD_ISIG_WAV
c	FUNCTION CALC_WAV_VARY(WAV)
C------------------------------------------------------------------------------


C----- Adjust LSQ weighting (using ISIG) for initial/final/correct runs -------

	SUBROUTINE ADD_ISIG_INITIAL
C
C This empirical weighting scheme works well for the preliminary refinement
C
C Add 1% error to all intensities, plus a term that stresses the importance
C of wavelengths around 1.25
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
	REL_ERR=0.01
	FACTOR=0.5
C
	WRITE(10,'(A,F5.1,A)') 'Modifying weights for preliminary refinements'
	NOBS=0
	DO ISEQ=1,NSEQ
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    ISIG(I)=ISIG(I) + ABS( IOBS(I)*REL_ERR )
	    ISIG(I)=ISIG(I) + ABS( IOBS(I)*FACTOR*(WAVS(I)-1.25)**2 )
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE ADD_ISIG_FINAL
C
C Adjust the weights via ISIG() suitable for the final parameter refinements
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	PARAMETER (NDATA_MAX=2000000)
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
C Multiply the weak sigmas by a factor
	CALL ADD_MULT_ISIG
C Adding fixed relative error to refinement esds
	CALL ADD_ISIG_REL
C Add estimated error due to uncertain wavelength correction
	CALL ADD_ISIG_WAV
C Add the extinction correction uncertainty, but only at the 10% level
	IF( LEXTI_CORR ) CALL ADD_ISIG_EXTI(0.1)
C
C Ensure ISIG() are > 1
C
	DO I=1,ILAST(NSEQ)
	  ISIG(I)=MAX(1.0, ISIG(I) )
	ENDDO
C
	RETURN
	END


	SUBROUTINE ADD_ISIG_CORRECT
C
C Adjust the weights via ISIG() suitable for the data correction
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	PARAMETER (NDATA_MAX=2000000)
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	REAL ISIG_SAVE(NDATA_MAX)
C
C Save current ISIG() values
C
	DO I=1,ILAST(NSEQ)
	  ISIG_SAVE(I)=ISIG(I)
	ENDDO
C
C Multiply the weak sigmas by a factor
	CALL ADD_MULT_ISIG
C Adding fixed relative error to refinement esds
	CALL ADD_ISIG_REL
C Add estimated error due to uncertain wavelength correction
	CALL ADD_ISIG_WAV
C Add the extinction correction uncertainty using EXTI_ERR 
	IF( LEXTI_CORR ) CALL ADD_ISIG_EXTI(EXTI_ERR)
C
C Ensure all flagged outliers (ISIG_SAVE() < 0) have negative ISIG()
C but all others have a positive ISIG()
C
	DO I=1,ILAST(NSEQ)
	  IF(ISIG_SAVE(I) .LT. 0.0) THEN
	    ISIG(I)=-ABS(ISIG(I))
	  ELSE
	    ISIG(I)=MAX(1.0, ISIG(I) )
	  ENDIF
	ENDDO
C
	RETURN
	END


C----- Components used in above routines --------------------------------------

	SUBROUTINE ADD_MULT_ISIG
C
C Multiply all sigmas by ALL_SIG_MULT and weak sigmas by WEAK_SIG_MULT as well
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
C Save the "background" sigma value, =0 means not initialised yet
C
	REAL ISIG_BKG
	DATA ISIG_BKG /0.0/
	SAVE ISIG_BKG
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
C If both factors = 1, nothing to do
C
	IF(ALL_SIG_MULT.EQ.1.0 .AND. WEAK_SIG_MULT.EQ.1.0) RETURN
C
C If not initialised, get the value of ISIG_BKG
C
	IF(ISIG_BKG .EQ. 0.0) THEN
C Average ISIG for all counts of 1 to 3 sigma
	  NSUM=0
	  SUM=0.0
	  DO ISEQ=1,NSEQ
	    DO I=IFIRST(ISEQ),ILAST(ISEQ)
	      IF(IOBS(I).GT.ISIG(I) .AND. IOBS(I).LT.3.0*ISIG(I)) THEN
	        SUM=SUM+ISIG(I)
	        NSUM=NSUM+1
	      ENDIF
	    ENDDO
	  ENDDO
	  ISIG_BKG=SUM/NSUM
	ENDIF
C
C Adjust ISIG such that weakest are scaled by WEAK_SIG_MULT, and
C all ISIG by ALL_SIG_MULT
C 
	DO ISEQ=1,NSEQ
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    IF(ISIG(I) .LT. ISIG_BKG) THEN
	      ISIG(I)=WEAK_SIG_MULT*ISIG(I)
	    ELSE
	      ISIG(I)=SQRT( ISIG(I)**2 - (1.0-WEAK_SIG_MULT**2)*ISIG_BKG**2 )
	    ENDIF
	    ISIG(I)=ALL_SIG_MULT*ISIG(I)
	  ENDDO
	ENDDO
C
C Log any corrections made to the esds
C
	IF(WEAK_SIG_MULT .NE. 1.0) WRITE(10,'(A,I4,A,F5.2)')
	1				'Multiplying errors for weak reflections (esd ~<',
	1					NINT(ISIG_BKG),') by',WEAK_SIG_MULT  
	IF(ALL_SIG_MULT .NE. 1.0) WRITE(10,'(A,F5.2)')
	1		'Multiplying errors for all reflections by',ALL_SIG_MULT
C
	RETURN
	END


	SUBROUTINE ADD_ISIG_REL
C
C Add ALL_REL_ERR times the observed intensity to ISIG
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
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
	WRITE(10,'(A,2P,F5.1,A)') 'Adding',ALL_REL_ERR,
	1							'% relative error to all intensities'
C
C
	ALL_REL_ERR_SQ=ALL_REL_ERR**2
	IF(ALL_REL_ERR .LT. 0.0) ALL_REL_ERR_SQ=-ALL_REL_ERR_SQ
C
	DO ISEQ=1,NSEQ
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
C Do RMS addition, as if uncorrelated random errors
	    ISIG(I)=SQRT(MAX(10.0, ISIG(I)**2 + ALL_REL_ERR_SQ*IOBS(I)**2 ))
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE ADD_ISIG_EXTI(EXTI_ERR)
C
C Add EXTI_ERR times the extinction correction to ISIG
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
	DO ISEQ=1,NSEQ
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    X=CALC_EXTI_FACTOR(COUNTS_SEQ(ISEQ),WAVS(I),TTHS(I))
	    Y=1.0/SQRT(MAX(1E-4,1.0+X))
C Do RMS addition, as if uncorrelated random errors
	    ISIG(I)=SQRT( ISIG(I)**2 + ( IOBS(I)*ABS(1.0-Y)*EXTI_ERR )**2 )
	  ENDDO
	ENDDO
C
	WRITE(10,'(A,I3,A)') 'Adding ',NINT(EXTI_ERR*100.0),
	1					'% error to the extinction correction factor'
C
	RETURN
	END


	SUBROUTINE ADD_ISIG_WAV
C
C Add to ISIG the uncertainty in the wavelength correction
C
C Increase ISIG() due to the uncertainty in calculating ICALC() due
C to the uncertainty in the wavelength correction, times WAV_ERR1 and
C WAV_ERR2/TAN(TH). The first accounts for a relative error in the
C wavelength, the second for an absolute error in the TH angle.
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
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
C Loop through all reflections, increasing ISIG() 
C
	DO I=1,ILAST(NSEQ)
	  VAR_ADD=( WAV_ERR1**2 + (WAV_ERR2/TAND(TTHS(I)/2.0))**2 ) *
	1								( IOBS(I) * CALC_WAV_VARY(WAVS(I)) )**2
	  ISIG(I)=SQRT( ISIG(I)**2 + VAR_ADD )
	ENDDO
C
	WRITE(10,'(A,F5.2,A)') 'Adding error due to uncertainty in wavelength and theta'
C
	RETURN
	END


	FUNCTION CALC_WAV_VARY(WAV)
C
C Calculates the relative error in calculating ICALC() due to the uncertainty
C in the wavelength assumed to lie in the range WAV +/- 0.01
C
C The refined wavelength scales factors and the wavelength dependence of the
C efficiency, absorption and extinction corrections is ignored. Only the
C refined wavelength distribution and the Wav^4/Sin(TH)^2 term is considered.
C Wav and Sin(TH) are 100% correlated, so the later reduces to a Wav^2 term.
C
C Subdivide wavelength range into 5 steps and measure the variation
C
	FMIN=+1E6
	FMAX=-1E6
	DO RWAV=-0.01,0.0101,0.005
	  WAV2=WAV + RWAV
	  F=CALC_WAV_SCALE(WAV2) * WAV2**2
	  FMIN=MIN(FMIN, F)
	  FMAX=MAX(FMAX, F)
	ENDDO
C
C Calculate and return half the relative variation from the min. value
C
	CALC_WAV_VARY=(FMAX-FMIN) / FMIN
C
	RETURN
	END
