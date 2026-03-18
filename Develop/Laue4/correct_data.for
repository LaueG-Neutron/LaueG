C----- Correct the data set using refined parameters --------------------------
c	SUBROUTINE CORRECT_DATA
C----- Output summaries about the intensity corrections -----------------------
c	SUBROUTINE PRINT_CORRECT_SUMMARY
c	SUBROUTINE LOG_CORRECT_SUMMARY
c	SUBROUTINE WRITE_CORRECT_SUMMARY(IUNIT,ISPACE)
c	SUBROUTINE CALC_CORRECT_SUMMARY(
c	1				FILE_MIN,FILE_MAX,FILE_AVE,FILE_ESD,
c	2				WAV_MIN,WAV_MAX,WAV_AVE,WAV_ESD,
c	3				EFF_MIN,EFF_MAX,EFF_AVE,EFF_ESD,
c	4				ABS_MIN,ABS_MAX,ABS_AVE,ABS_ESD,
c	5				EXTI_MIN,EXTI_MAX,EXTI_AVE,EXTI_ESD		)
C------------------------------------------------------------------------------


C----- Correct the data set using refined parameters --------------------------

	SUBROUTINE CORRECT_DATA
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	WRITE(10,'(/,A,/)') '====== Apply refined corrections ======'
C
C Load IOBS & ISIG from COUNTS & DCOUNTS
C
	CALL LOAD_LSQ_DATA
	CALL OUTPUT_DEBUG_DATA
	CALL OUTPUT_DEBUG_CALC
C
C Output merge statistics before we fiddle with the esds
C
 	CALL CALC_ALL_COUNTS_SEQ
	CALL OUTPUT_DEBUG_SEQ
	WRITE(10,'(A)') 'Merge statistics using initial weights'
	CALL LOG_LSQ_MERGE
C
C Adjust ISIG using various terms
C
	CALL ADD_ISIG_CORRECT
C
C Recalculate merged intensities using the final weights
C
	CALL CALC_ALL_COUNTS_SEQ
	CALL OUTPUT_DEBUG_DATA
	CALL OUTPUT_DEBUG_CALC
	CALL OUTPUT_DEBUG_SEQ
C
C Output correction factors and merge statistics to the console
C
	PRINT '(/,1X,A)','Summary of refineable correction factors:'
	CALL PRINT_CORRECT_SUMMARY
	PRINT '(/,1X,A)','Merge statistics for corrected intensities:'
	CALL PRINT_LSQ_MERGE
C
C Output correction factors and merge statistics to the log file
C
	WRITE(10,'(/,A)') 'Merge statistics using final weights:'
	CALL LOG_LSQ_MERGE
	WRITE(10,'(/,A)') 'Summary of refineable correction factors:'
	CALL LOG_CORRECT_SUMMARY
C
C If required, do harmonic correction to intensities
C
	IF(LHARM_CORR) THEN
	  CALL HARMONIC_CORR
	  CALL OUTPUT_DEBUG_DATA
	  CALL OUTPUT_DEBUG_CALC
	  CALL OUTPUT_DEBUG_SEQ
	ENDIF
C
	RETURN
	END
	

C----- Output summaries about the intensity corrections -----------------------

	SUBROUTINE PRINT_CORRECT_SUMMARY
C
C Output a summary of the size of the corrections being used
C
	CALL WRITE_CORRECT_SUMMARY(6,3)
	RETURN
	END



	SUBROUTINE LOG_CORRECT_SUMMARY
C
C Output a summary of the size of the corrections being used
C
	CALL WRITE_CORRECT_SUMMARY(10,2)
	RETURN
	END



	SUBROUTINE WRITE_CORRECT_SUMMARY(IUNIT,ISPACE)
C
C Output a summary of the size of the corrections being used
C
	CALL CALC_CORRECT_SUMMARY(
	1				FILE_MIN,FILE_MAX,FILE_AVE,FILE_ESD,
	2				WAV_MIN,WAV_MAX,WAV_AVE,WAV_ESD,
	3				EFF_MIN,EFF_MAX,EFF_AVE,EFF_ESD,
	4				ABS_MIN,ABS_MAX,ABS_AVE,ABS_ESD,
	5				EXTI_MIN,EXTI_MAX,EXTI_AVE,EXTI_ESD		)
C
	WRITE(IUNIT,'(<ISPACE>X,2A,2P,4(F5.1,A))',IOSTAT=IDUM)
	1		'File scale factor ',
	2		'   min=',FILE_MIN-1.0,'%  max=',FILE_MAX-1.0,
     3		'%  ave=',FILE_AVE-1.0,'%  st.dev.=',FILE_ESD,'%'
C
	WRITE(IUNIT,'(<ISPACE>X,2A,2P,4(F5.1,A))',IOSTAT=IDUM)
	1		'Wavelength spectra',
	2		'   min=',WAV_MIN-1.0,'%  max=',WAV_MAX-1.0,
     3		'%  ave=',WAV_AVE-1.0,'%  st.dev.=',WAV_ESD,'%'
C
	WRITE(IUNIT,'(<ISPACE>X,2A,2P,4(F5.1,A))',IOSTAT=IDUM)
	1		'Refined efficiency',
	2		'   min=',EFF_MIN-1.0,'%  max=',EFF_MAX-1.0,
     3		'%  ave=',EFF_AVE-1.0,'%  st.dev.=',EFF_ESD,'%'
C
	WRITE(IUNIT,'(<ISPACE>X,2A,2P,4(F5.1,A))',IOSTAT=IDUM)
	1		'Sample absorption ',
	2		'   min=',ABS_MIN-1.0,'%  max=',ABS_MAX-1.0,
     3		'%  ave=',ABS_AVE-1.0,'%  st.dev.=',ABS_ESD,'%'
C
	WRITE(IUNIT,'(<ISPACE>X,2A,2P,4(F5.1,A))',IOSTAT=IDUM)
	1		'Sample extinction ',
	2		'   min=',EXTI_MIN-1.0,'%  max=',EXTI_MAX-1.0,
     3		'%  ave=',EXTI_AVE-1.0,'%  st.dev.=',EXTI_ESD,'%'
C
	RETURN
	END



	SUBROUTINE CALC_CORRECT_SUMMARY(
	1				FILE_MIN,FILE_MAX,FILE_AVE,FILE_ESD,
	2				WAV_MIN,WAV_MAX,WAV_AVE,WAV_ESD,
	3				EFF_MIN,EFF_MAX,EFF_AVE,EFF_ESD,
	4				ABS_MIN,ABS_MAX,ABS_AVE,ABS_ESD,
	5				EXTI_MIN,EXTI_MAX,EXTI_AVE,EXTI_ESD		)
C
C Output a summary of the size of the corrections being used
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
	FILE_SUM=0.0
	FILE_SUM2=0.0
	FILE_MIN=+1E6
	FILE_MAX=-1E6
C
	WAV_SUM=0.0
	WAV_SUM2=0.0
	WAV_MIN=+1E6
	WAV_MAX=-1E6
C
	EFF_SUM=0.0
	EFF_SUM2=0.0
	EFF_MIN=+1E6
	EFF_MAX=-1E6
C
	ABS_SUM=0.0
	ABS_SUM2=0.0
	ABS_MIN=+1E6
	ABS_MAX=-1E6
C
	EXTI_SUM=0.0
	EXTI_SUM2=0.0
	EXTI_MIN=+1E6
	EXTI_MAX=-1E6
C
	DO ISEQ=1,NSEQ
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
C
	    FILE=1.0/CALC_FILE_SCALE(IFILES(I))
	    FILE_SUM=FILE_SUM+FILE
	    FILE_SUM2=FILE_SUM2+FILE**2
	    FILE_MIN=MIN(FILE_MIN,FILE)
	    FILE_MAX=MAX(FILE_MAX,FILE)
C
	    WAV=CALC_WAV_ROUGH(WAVS(I))/CALC_WAV_SCALE(WAVS(I))
	    WAV_SUM=WAV_SUM+WAV
	    WAV_SUM2=WAV_SUM2+WAV**2
	    WAV_MIN=MIN(WAV_MIN,WAV)
	    WAV_MAX=MAX(WAV_MAX,WAV)
C
	    EFF=1.0/CALC_EFFIC_REFINE(XY_PIX(2,I))
	    EFF_SUM=EFF_SUM+EFF
	    EFF_SUM2=EFF_SUM2+EFF**2
	    EFF_MIN=MIN(EFF_MIN,EFF)
	    EFF_MAX=MAX(EFF_MAX,EFF)
C
	    ABCOR=1.0/CALC_ABS_FACTOR(WAVS(I),TTHS(I))
	    ABS_SUM=ABS_SUM+ABCOR
	    ABS_SUM2=ABS_SUM2+ABCOR**2
	    ABS_MIN=MIN(ABS_MIN,ABCOR)
	    ABS_MAX=MAX(ABS_MAX,ABCOR)
C
	    EXTI_X=CALC_EXTI_FACTOR(COUNTS_SEQ(ISEQ),WAVS(I),TTHS(I))
	    EXTI_Y=SQRT(MAX(1E-4,1.0+EXTI_X))	! reciprocal of factor
	    EXTI_SUM=EXTI_SUM+EXTI_Y
	    EXTI_SUM2=EXTI_SUM2+EXTI_Y**2
	    EXTI_MIN=MIN(EXTI_MIN,EXTI_Y)
	    EXTI_MAX=MAX(EXTI_MAX,EXTI_Y)
C
	  ENDDO
	ENDDO
C
	FILE_AVE=FILE_SUM/NDATA
	FILE_ESD=SQRT(MAX(0.0, FILE_SUM2/NDATA - FILE_AVE**2 ))
C
	WAV_AVE=WAV_SUM/NDATA
	WAV_ESD=SQRT(MAX(0.0, WAV_SUM2/NDATA - WAV_AVE**2 ))
C
	EFF_AVE=EFF_SUM/NDATA
	EFF_ESD=SQRT(MAX(0.0, EFF_SUM2/NDATA - EFF_AVE**2 ))
C
	ABS_AVE=ABS_SUM/NDATA
	ABS_ESD=SQRT(MAX(0.0, ABS_SUM2/NDATA - ABS_AVE**2 ))
C
	EXTI_AVE=EXTI_SUM/NDATA
	EXTI_ESD=SQRT(MAX(0.0, EXTI_SUM2/NDATA - EXTI_AVE**2 ))
C
	RETURN
	END
