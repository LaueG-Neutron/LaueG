C----- Routines to output final data and list files ---------------------------
c	SUBROUTINE OUTPUT_FINAL_RESULTS
c	SUBROUTINE CALC_BEST_SCALES(SCALE_F,SCALE_C)
C----- Output any systematic absences -----------------------------------------
c	SUBROUTINE OUTPUT_SYSTEM_ABSENCES
C----- Output the CIF file ----------------------------------------------------
c	SUBROUTINE OUTPUT_CIF_FILE
C----- Output the merged and un-merged hkl intensity files --------------------
c	SUBROUTINE OUTPUT_MERGED_INTS(FILE_NAME,SCALE_C)
c	SUBROUTINE OUTPUT_UNMERGED_INTS(FILE_NAME,SCALE_C)
C----- Output the wav. calibration file + related info. -----------------------
c	SUBROUTINE OUTPUT_WAV_DATA(FILE_NAME)
C----- Output details of individual merge & corrections -----------------------
c	SUBROUTINE OUTPUT_MERGE_INFO(FILE_NAME,SCALE_F,SCALE_C)
c	SUBROUTINE CALC_EQUIV_MOD(HKLM_EQV, HKLM)
C----- Output the extended information intensity file -------------------------
c	SUBROUTINE OUTPUT_VERBOSE_INTS(FILE_NAME)
C----- Output the efficiency correction file ----------------------------------
c	SUBROUTINE OUTPUT_EFFIC_FILE(FILE_NAME)
C----- Output the parts of the CIF file ---------------------------------------
c	SUBROUTINE WRITE_CIF_CELL
c	SUBROUTINE WRITE_CIF_INSTRUM
c	SUBROUTINE WRITE_CIF_EXPERI
c	SUBROUTINE WRITE_CIF_LAUEG
c	SUBROUTINE WRITE_CIF_SOFTWARE
c	SUBROUTINE WRITE_CIF_REFERENCES
c	SUBROUTINE GET_EXPERI_SUMMARY(SUMMARY)
C------------------------------------------------------------------------------


C----- Routines to output final data and list files ---------------------------

	SUBROUTINE OUTPUT_FINAL_RESULTS
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	CHARACTER BASE_NAME*80
C
	WRITE(10,'(/,A)') '====== Output of final results ======'
C
C Determine a reasonable factor-of-ten scale factor for IOBS() & ICALC().
C
	CALL CALC_BEST_SCALES(SCALE_F,SCALE_C)
C
C Output the merged & unmerged *.HKL files using the stored base name
C
	CALL MAKE_FILE_NAME(BASE_NAME, -1,'')
C
	CALL OUTPUT_MERGED_INTS(TRIM(BASE_NAME)//'_mrg.hkl',SCALE_C)
	CALL OUTPUT_UNMERGED_INTS(TRIM(BASE_NAME)//'_wav.hkl',SCALE_C)
C
C Output the new wavelength distribution file and related info.
C 
	CALL OUTPUT_WAV_DATA('laue4_wav.new')
C
C In verbose mode we make some extra files
C
	IF(IVERBOSE .NE. 0) THEN
	  CALL OUTPUT_MERGE_INFO('laue4_mrg.lis',SCALE_F,SCALE_C)
	  CALL OUTPUT_VERBOSE_INTS('laue4_all.int')
	  CALL OUTPUT_EFFIC_FILE('laue4_eff.lis')
	ENDIF
C
C Output Jana information file(for Option 1), and cif file
C
	CALL OUTPUT_CIF_FILE
	WRITE(10,'(/,2A)') 'laue4.cif file created'
C
C And if used, output the final debug info
C
	CALL OUTPUT_DEBUG_CALC
	CALL OUTPUT_DEBUG_SEQ
C
	RETURN
	END



	SUBROUTINE CALC_BEST_SCALES(SCALE_F,SCALE_C)
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	PARAMETER (NDATA_MAX=2000000)
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	REAL DCOUNTS_SEQ(NSEQ_MAX)
C
C Calculate the esds in COUNTS_SEQ() from the LSQ uncertainty in the fit
C
	CALL CALC_ALL_DCOUNTS_SEQ(DCOUNTS_SEQ)
C
C Determine a reasonable factor-of-ten scale factor for IOBS() & ICALC().
C
	FMAX=-1E9
	FMIN=+1E9
	DO I=1,ILAST(NSEQ)
C Ignore outliers which are marked with ISIG() < 0
	  IF(ISIG(I) .GE. 0.0) THEN
	    FMAX=MAX(FMAX, IOBS(I),ICALC(I))
	    FMIN=MIN(FMIN, ISIG(I))
	  ENDIF
	ENDDO
C
	SCALE_F=1.0
C Increase SCALE_F so the smallest ISIG is > 0.3
	IF(FMIN*SCALE_F .LT. 0.3) SCALE_F=SCALE_F*10.0
	IF(FMIN*SCALE_F .LT. 0.3) SCALE_F=SCALE_F*10.0
	IF(FMIN*SCALE_F .LT. 0.3) SCALE_F=SCALE_F*10.0
C Decrease SCALE_F so the largest IOBS and ICALC is < 1e5
	IF(FMAX*SCALE_F .GE. 1E5) SCALE_F=SCALE_F/10.0
	IF(FMAX*SCALE_F .GE. 1E5) SCALE_F=SCALE_F/10.0
	IF(FMAX*SCALE_F .GE. 1E5) SCALE_F=SCALE_F/10.0
C
C Determine a reasonable factor-of-ten scale factor for COUNTS_SEQ()
C
	CMAX=-1E9
	CMIN=+1E9
	DO ISEQ=1,NSEQ
C Ignore any sequences containing only outliers
	  IF(COUNTS_SEQ(ISEQ) .GE. 0.0) THEN
	    CMAX=MAX(CMAX, COUNTS_SEQ(ISEQ))
	    CMIN=MIN(CMIN, DCOUNTS_SEQ(ISEQ))
	  ENDIF
	ENDDO
C
	SCALE_C=1.0
C Increase SCALE_C so the smallest DCOUNTS_SEQ is > 0.3
	IF(CMIN*SCALE_C .LT. 0.3) SCALE_C=SCALE_C*10.0
	IF(CMIN*SCALE_C .LT. 0.3) SCALE_C=SCALE_C*10.0
	IF(CMIN*SCALE_C .LT. 0.3) SCALE_C=SCALE_C*10.0
C Decrease SCALE_C so the largest COUNTS_SEQ is < 1e5
	IF(CMAX*SCALE_C .GE. 1E5) SCALE_C=SCALE_C/10.0
	IF(CMAX*SCALE_C .GE. 1E5) SCALE_C=SCALE_C/10.0
	IF(CMAX*SCALE_C .GE. 1E5) SCALE_C=SCALE_C/10.0
C
	RETURN
	END


C-------------- Output any systematic absences -----------

	SUBROUTINE OUTPUT_SYSTEM_ABSENCES
C
C Output info on odd/even systematic absences and suggest possible centerings
C
	PARAMETER (NSEQ_MAX=100000)
	REAL DCOUNTS_SEQ(NSEQ_MAX)
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	CHARACTER STRING*80
	INTEGER NSUM(8)
	REAL SUM(8),I_CEN
C
	WRITE(10,'(/,A,/)') '====== Search for odd/even systematic absences ======'
C
C Calculate the esds in COUNTS_SEQ() from the LSQ uncertainty in the fit
C
	CALL CALC_ALL_DCOUNTS_SEQ(DCOUNTS_SEQ)
C
C Zero the sum arrays
C
	DO I=1,8
	  SUM(I)=0.0
	  NSUM(I)=0
	ENDDO
C
C Sum various odd/even arrangements according to:
C   Absences  eee eeo eoe eoo oee oeo ooe ooo
C   Number     1   2   3   4   5   6   7   8
C
	DO ISEQ=1,NSEQ
C Ignore "rejected" sequences
	  IF( DCOUNTS_SEQ(ISEQ) .LE. 0.0) CYCLE
C Ignore h*H + k*k +l*l > 100
	    I=IFIRST(ISEQ)
	    IF(HKLMS(1,I)**2 + HKLMS(2,I)**2 + HKLMS(3,I)**2 .GT. 100) CYCLE
C Create the Number of the odd/even arrangement
	    IHMOD2=ABS(HKLMS(1,I) - HKLMS(1,I)/2*2)
	    IKMOD2=ABS(HKLMS(2,I) - HKLMS(2,I)/2*2)
	    ILMOD2=ABS(HKLMS(3,I) - HKLMS(3,I)/2*2)
	    I123=4*IHMOD2 +2*IKMOD2 +ILMOD2 +1
C Add to the appropriate sums
	  SUM(I123)=SUM(I123)+COUNTS_SEQ(ISEQ)/DCOUNTS_SEQ(ISEQ)
	  NSUM(I123)=NSUM(I123)+1
	ENDDO
C
C Convert sums to averages, but avoid any divide by zeroes
C
	DO I=1,8
	  SUM(I)=SUM(I)/MAX(1,NSUM(I))
	ENDDO
C
C Output to log file the results for various odd/evens
C
	WRITE(10,'(2X,A)') 'hkl        eee   eeo   eoe   eoo   oee   oeo   ooe   ooo'
	WRITE(10,'(2X,A,8F6.1,A)') 'Ave Ints',SUM,' sig'
	WRITE(10,'(2X,A,8I6)') 'Num Ave',NSUM
C
C Get largest "absence" for various odd/even centering types
C
C A centered:	2,3,6,7
	A_CEN=MAX( SUM(2) , SUM(3) , SUM(6) , SUM(7) )
C B centered:	2,4,5,7
	B_CEN=MAX( SUM(2) , SUM(4) , SUM(5) , SUM(7) )
C C centered:	3,4,5,6
	C_CEN=MAX( SUM(3) , SUM(4) , SUM(5) , SUM(6) )
C I centered:	2,3,5,8
	I_CEN=MAX( SUM(2) , SUM(3) , SUM(5) , SUM(8) )
C F centered:	2,3,4,5,6,7
	F_CEN=MAX( SUM(2) , SUM(3) , SUM(4) , SUM(5) , SUM(6) , SUM(7) )
C A  doubled: 5,6,7,8
	A_DBL=MAX( SUM(5) , SUM(6) , SUM(7) , SUM(8) )
C B  doubled: 3,4,7,8
	B_DBL=MAX( SUM(3) , SUM(4) , SUM(7) , SUM(8) )
C C  doubled: 2,4,6,8
	C_DBL=MAX( SUM(2) , SUM(4) , SUM(6) , SUM(8) )
C
C Get max and min of possible absences
C
	SYS_MIN=MIN(A_CEN,B_CEN,C_CEN,I_CEN,F_CEN,A_DBL,B_DBL,C_DBL)
	SYS_MAX=MAX(A_CEN,B_CEN,C_CEN,I_CEN,F_CEN,A_DBL,B_DBL,C_DBL)
C
C If no significant absences, output to logfile and return
C
	IF(SYS_MIN .GT. MAX(2.0,SYS_MAX/10.0) ) THEN
	  WRITE(10,'(/,2X,A)') 'No systematic absences detected'
	  RETURN
	ENDIF
C
C If HKLs generated under centering, output to logfile and return
C
	IF( F_CEN .EQ. 0.0) THEN
	  WRITE(10,'(/,A)') '  HKLs generated using F centering'
	ELSEIF( I_CEN .EQ. 0.0) THEN
	  WRITE(10,'(/,A)') '  HKLs generated using I centering'
	ELSEIF( A_CEN .EQ. 0.0) THEN
	  WRITE(10,'(/,A)') '  HKLs generated using A centering'
	ELSEIF( B_CEN .EQ. 0.0) THEN
	  WRITE(10,'(/,A)') '  HKLs generated using B centering'
	ELSEIF( C_CEN .EQ. 0.0) THEN
	  WRITE(10,'(/,A)') '  HKLs generated using C centering'
	ENDIF
	IF(SYS_MIN .EQ. 0.0) THEN
	  WRITE(10,'(/,A)') '  No tests performed for systematic absences'
	  RETURN
	ENDIF
C
C Must have significant absences, so load STRING with the best explanation
C
	CUTOFF=1.5*SYS_MIN +0.1
C If ooo is absent: possible I centering or Cell doubling
	IF(SUM(8) .LT. CUTOFF) THEN
	  IF(I_CEN .LT. CUTOFF) THEN
	    STRING='I centered'
	  ELSEIF(MAX(A_DBL,B_DBL,C_DBL) .LT. CUTOFF) THEN
	    STRING='doubled along a,b,c'
	  ELSEIF(MAX(A_DBL,B_DBL) .LT. CUTOFF) THEN
	    STRING='doubled along a,b'
	  ELSEIF(MAX(B_DBL,C_DBL) .LT. CUTOFF) THEN
	    STRING='doubled along b,c'
	  ELSEIF(MAX(A_DBL,C_DBL) .LT. CUTOFF) THEN
	    STRING='doubled along a,c'
	  ELSEIF(A_DBL .LT. CUTOFF) THEN
	    STRING='doubled along a'
	  ELSEIF(B_DBL .LT. CUTOFF) THEN
	    STRING='doubled along b'
	  ELSEIF(C_DBL .LT. CUTOFF) THEN
	    STRING='doubled along c'
	  ELSE
	    STRING='cell doubled twins (if Tetra or Cubic)'
	  ENDIF
C ELSE: possible A,B,C,F centering
	ELSE
	  IF(F_CEN .LT. CUTOFF) THEN
	    STRING='F centered'
	  ELSEIF(MAX(A_CEN,B_CEN) .LT. CUTOFF) THEN
	    STRING='A & B centered'
	  ELSEIF(MAX(A_CEN,C_CEN) .LT. CUTOFF) THEN
	    STRING='A & C centered'
	  ELSEIF(MAX(B_CEN,C_CEN) .LT. CUTOFF) THEN
	    STRING='B & C centered'
	  ELSEIF(A_CEN .LT. CUTOFF) THEN
	    STRING='A centered'
	  ELSEIF(B_CEN .LT. CUTOFF) THEN
	    STRING='B centered'
	  ELSEIF(C_CEN .LT. CUTOFF) THEN
	    STRING='C centered'
C No idea what it is, so just return
	  ELSE
	    RETURN
	  ENDIF
	ENDIF
C
C Output the explanation of absences to log file and console
C
	IF(SYS_MIN .LT. 1.0) THEN
	  STRING='Systematic absences are consistent with '//TRIM(STRING)
	ELSE
	  STRING='Systematic absences are possibly consistent with '//TRIM(STRING)
	ENDIF
	WRITE(10,'(/,2X,A)') TRIM(STRING)
	PRINT '(1X,A,/)',TRIM(STRING)
C
	RETURN
	END


C-------------- Output the CIF file -----------

	SUBROUTINE OUTPUT_CIF_FILE
C
C Open the CIF file on UNIT 2
C
	OPEN(UNIT=2,FILE='laue4.cif',STATUS='UNKNOWN')
C
C Audit records (needed?)
C
	WRITE(2,'(A,//,(A))') "data_1","_audit_block_code       LaueG",
	1			"_audit_creation_method	'LaueG software'"
C
C Cell and lattice records
C
	WRITE(2,'(/,1H#,79(1H=),/,A,/)') "# CRYSTAL"
	CALL WRITE_CIF_CELL
C
C Instrument records
C
	WRITE(2,'(/,1H#,79(1H=),/,A,/)') "# INSTRUMENT"
	CALL WRITE_CIF_INSTRUM
C
C Experiment records
C
	WRITE(2,'(/,1H#,79(1H=),/,A,/)') "# EXPERIMENT"
	CALL WRITE_CIF_EXPERI
C
C LaueG specific records
C
	WRITE(2,'(/,1H#,79(1H=),/,A,/)') "# LAUEG SPECIFIC"
	CALL WRITE_CIF_LAUEG
C
C Software records
C
	WRITE(2,'(/,1H#,79(1H=),/,A,/)') "# SOFTWARE"
	CALL WRITE_CIF_SOFTWARE
C
C References records
C
	WRITE(2,'(/,1H#,79(1H=),/,A,/)') "# REFERENCES"
	CALL WRITE_CIF_REFERENCES
C
	CLOSE(UNIT=2)
C
	RETURN
	END


C-------------- Output the merged and un-merged hkl intensity files -----------

	SUBROUTINE OUTPUT_MERGED_INTS(FILE_NAME,SCALE_C)
C
	CHARACTER*(*) FILE_NAME
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	PARAMETER (NSEQ_MAX=100000)
	REAL DCOUNTS_SEQ(NSEQ_MAX)
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	REAL MOD_VEC
	COMMON /MODULATE_COM/ NMOD_VEC,NMOD_IDX,MOD_VEC(3,100),MOD_IDX(100,100)
C
	CHARACTER STEMP*100
	INTEGER HKL_EQV(3),HKLM_EQV(4)
C
C Calculate the esds in COUNTS_SEQ() from the LSQ uncertainty in the fit
C
	CALL CALC_ALL_DCOUNTS_SEQ(DCOUNTS_SEQ)
C
C Output the corrected-merged intensities
C In modulation mode, use a *.hklq file with HKL and q-vector indices
C
	IF(IDATA_MODE .NE. 3) THEN
	  OPEN(UNIT=3,FILE=FILE_NAME,STATUS='UNKNOWN')
	ELSE
	  OPEN(UNIT=3,FILE=FILE_NAME//'q',STATUS='UNKNOWN')
	ENDIF
C
	I_BN=1
	NOUT=0
	DO ISEQ=1,NSEQ
C Ignore sequences that only contain outliers
	  IF(DCOUNTS_SEQ(ISEQ) .LT. 0.0) CYCLE
C Output equivalent HKL, sequence intensity & esd, and batch number
	  NOUT=NOUT+1
	  I=IFIRST(ISEQ)
      IF(IDATA_MODE .LE. 2) THEN
        IF(IDATA_MODE .EQ. 2) I_BN=HKLMS(4,I)
        CALL CALC_EQUIV_HKL(HKL_EQV,HKLMS(1,I))
	    WRITE(3,'(3I4,2F8.1,I4)') HKL_EQV,SCALE_C*COUNTS_SEQ(ISEQ),
	1							SCALE_C*DCOUNTS_SEQ(ISEQ),I_BN
	  ELSE
	    CALL CALC_EQUIV_MOD(HKLM_EQV, HKLMS(1,I))
        IMOD_EQV=HKLM_EQV(4)
        IF(IMOD_EQV .NE. 0) THEN
	      WRITE(STEMP,'(<NMOD_VEC>I4)') (MOD_IDX(IMOD_EQV,K),K=1,NMOD_VEC)
		ELSE
	      WRITE(STEMP,'(<NMOD_VEC>I4)') (0,K=1,NMOD_VEC)
		ENDIF
	    WRITE(3,'(3I4,A,2F8.1)') (HKLM_EQV(K),K=1,3),TRIM(STEMP),          
	1				SCALE_C*COUNTS_SEQ(ISEQ),SCALE_C*DCOUNTS_SEQ(ISEQ)
	  ENDIF
C Loop for more sequences
	ENDDO
C
	CLOSE(UNIT=3)
C
C Output summary of intensity files written
C
	IF(IDATA_MODE .LE. 2) THEN
	  PRINT '(1X,I6,2A)',NOUT,' merged intensities written to ',FILE_NAME
	  WRITE(10,'(I6,2A)') NOUT,' merged intensities written to ',FILE_NAME
	  WRITE(10,'(2X,2A)') 'File in SHELX HKLF 5 format'
	  IF(IDATA_MODE .EQ. 2) WRITE(10,'(2X,A)') 'Final integer is the Twin number'
	ELSE
	  PRINT '(1X,I6,2A)',NOUT,' merged intensities written to ',FILE_NAME//'q'
	  WRITE(10,'(I6,2A)') NOUT,' merged intensities written to ',FILE_NAME//'q'
	  WRITE(10,'(2X,2A)') 'File uses HKL and wavevector indices'
	ENDIF
	WRITE(10,'(2X,A,F6.2)') 'Intensities multiplied by',SCALE_C
C
	RETURN
	END


	SUBROUTINE OUTPUT_UNMERGED_INTS(FILE_NAME,SCALE_C)
C
	CHARACTER*(*) FILE_NAME
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
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
	COMMON /EXTI_CORR_COM/ IEXTI_OPT,EXTI
C
	REAL MOD_VEC
	COMMON /MODULATE_COM/ NMOD_VEC,NMOD_IDX,MOD_VEC(3,100),MOD_IDX(100,100)
C
	CHARACTER STEMP*100
C
C Output reflections corrected for everything except extinction
C In modulation mode, use a *.hklq file with HKL and q-vector indices
C
	IF(IDATA_MODE .NE. 3) THEN
	  OPEN(UNIT=3,FILE=FILE_NAME,STATUS='UNKNOWN')
	ELSE
	  OPEN(UNIT=3,FILE=FILE_NAME//'q',STATUS='UNKNOWN')
	ENDIF
C
	I_BN=1
	NOUT=0
	DO I=1,ILAST(NSEQ)
C Ignore outliers which are marked with ISIG() < 0
	  IF(ISIG(I) .LT. 0.0) CYCLE
C Use the Batch Number (BN) for the Twin number
	  IF(IDATA_MODE .EQ. 2) I_BN=HKLMS(4,I)
C Remove extinction correction by calculating ICALC/COUNTS_SEQ (=DERIV)
	  CALL DERIV_COUNTS_SEQ(DUMMY,DERIV, 0.0,I)
	  R_OBS=SCALE_C*IOBS(I)/DERIV
	  R_SIG=SCALE_C*ISIG(I)/DERIV
C Output HKL, uncorrected-intensity & esd, batch number, and wavelength
	  NOUT=NOUT+1
	  IF(IDATA_MODE .NE. 3) THEN
	    WRITE(3,'(3I4,2F8.1,I4,F8.3)') (HKLMS(K,I),K=1,3),R_OBS,R_SIG,I_BN,WAVS(I)
	  ELSE
		IMOD=HKLMS(4,I)
	    IF(IMOD .NE. 0) THEN
	      WRITE(STEMP,'(<NMOD_VEC>I4)') (MOD_IDX(IMOD,K),K=1,NMOD_VEC)
		ELSE
	      WRITE(STEMP,'(<NMOD_VEC>I4)') (0,K=1,NMOD_VEC)
		ENDIF
	    WRITE(3,'(3I4,A,2F8.1,I4,F8.3)') (HKLMS(K,I),K=1,3),TRIM(STEMP),R_OBS,R_SIG,1,WAVS(I)
	  ENDIF
C
	ENDDO
C
	CLOSE(UNIT=3)
C
C Output summary of intensity files written
C
	IF(IDATA_MODE .NE. 3) THEN
	  PRINT '(1X,I6,2A)',NOUT,' unmerged intensities written to ',FILE_NAME
	  WRITE(10,'(I6,2A)') NOUT,' unmerged intensities with wavelength written to ',FILE_NAME
	  WRITE(10,'(2X,2A)') 'File in SHELX HKLF 2 format'
	  IF(IDATA_MODE .EQ. 2) WRITE(10,'(2X,A)')
	1						'Integer before wavelength is the Twin number'
	ELSE
	  PRINT '(1X,I6,2A)',NOUT,' unmerged intensities written to ',FILE_NAME//'q'
	  WRITE(10,'(I6,2A)') NOUT,' unmerged intensities with wavelength written to ',FILE_NAME//'q'
	  WRITE(10,'(2X,2A)') 'File uses HKL and wavevector indices'
	ENDIF
	WRITE(10,'(2X,A)') 'Intensities are not extinction corrected'
	WRITE(10,'(2X,A,F6.2)') 'Intensities multiplied by',SCALE_C
C
	RETURN
	END


C --------------- Output the wav. calibration file + related info. ------------

	SUBROUTINE OUTPUT_WAV_DATA(FILE_NAME)
C
	CHARACTER FILE_NAME*(*)
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	COMMON /WAV_CORR_COM/ IWAV_OPT,NWAV,WAV_GAMMA,WAV_MIN,WAV_STEP,WAV_PAR(20)
	COMMON /WAV_FILE_COM/ NWAV_FILE,WAV_FILE_MIN,WAV_FILE_STEP,WAV_FILE(10000)
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	CHARACTER HOST_NAME*10,DATE*10
	REAL WAV_FILE2(10000)
C
C If no wavelength refinement, delete the new wavelength file and return
C
	IF( .NOT.LWAV_CORR ) THEN
	  OPEN(UNIT=3,FILE=FILE_NAME,STATUS='UNKNOWN')
	  CLOSE(UNIT=3,STATUS='DELETE')
	  RETURN
	ENDIF
C
C Write something to the log file
C
	WRITE(10,'(/,2A)') 'Writing new wavelength correction file to ',FILE_NAME
C
C Read the wavelength input file to get the nominal distribution which
C may be present as the third column of the file
C
	OPEN(UNIT=1,STATUS='OLD',FILE='laue4_wav.dat')
C Skip header and info lines
	READ(1,*)
	READ(1,*)
	READ(1,*)
C Try to read data in 3 column format
	DO I=1,NWAV_FILE
	  READ(1,*,IOSTAT=IERR) WAV,DUMMY,WAV_FILE2(I)
	ENDDO
	CLOSE(UNIT=1)
C If the third column is missing, warn in the log file and copy the
C input distribution to the nominal
	IF(IERR .NE. 0) THEN
	  WRITE(10,'(A)') 'WARNING: Nominal wavelength distribution data is missing'
	  DO I=1,NWAV_FILE
	    WAV_FILE2(I)=WAV_FILE(I)
	  ENDDO
	ENDIF
C
C Open file and write a short header
C
	OPEN(UNIT=3,FILE=FILE_NAME,STATUS='UNKNOWN')
C
	WRITE(3,'(A)') '===LAUE4 WAVELENGTH CALIBRATION FILE==='
C
	CALL GET_HOST_DATE(HOST_NAME,DATE)
	WRITE(3,'(4A)') 'Instrument: ',HOST_NAME(1:9),' Calibration date: ',DATE
C
	IF(IWAV_OPT .EQ. 1) THEN
	  WRITE(3,'(A,F6.3)') 'Comments: Generated by LAUE4 using quad-spline method'
	ELSE
	  WRITE(3,'(A,F6.3)') 'Comments: Generated by LAUE4 using spline & non-param method'
	ENDIF
C
C Output the refined and nominal distributions and close the file
C
	DO I=1,NWAV_FILE
	  WAV=WAV_FILE_MIN+WAV_FILE_STEP*(I-1)
	  WRITE(3,'(F6.3,2F9.6)') WAV,MAX(1E-4,CALC_WAV_SCALE(WAV)),WAV_FILE2(I)
	ENDDO
C
	CLOSE(UNIT=3)
C
C Calculate the weighted mean for the initial, final and nominal wavelength
C distributions. The 1/WAV**2	weighting downscales the long wavelengths.
C
	SNUM0=0.0
	SNUM1=0.0
	SNUM2=0.0
	SDEN0=0.0
	SDEN1=0.0
	SDEN2=0.0
	DO I=1,NWAV_FILE
	  WAV=WAV_FILE_MIN+WAV_FILE_STEP*(I-1)
	  CALC0=WAV_FILE2(I)/WAV**2
	  CALC1=WAV_FILE(I)/WAV**2
	  CALC2=CALC_WAV_SCALE(WAV)/WAV**2
	  SNUM0=SNUM0+CALC0*WAV
	  SNUM1=SNUM1+CALC1*WAV
	  SNUM2=SNUM2+CALC2*WAV
	  SDEN0=SDEN0+CALC0
	  SDEN1=SDEN1+CALC1
	  SDEN2=SDEN2+CALC2
	ENDDO
C
C Output info on wavelength shifts to log file
C
	RAT_21=( SNUM2/SDEN2 ) / ( SNUM1/SDEN1 )
	RAT_20=( SNUM2/SDEN2 ) / ( SNUM0/SDEN0 )
	WRITE(10,'(/,A,SP,2P,F6.1,A)') 'Wav. distribution shifted by',
	1			RAT_21-1.0,'% (this refinement)'
	WRITE(10,'(A,SP,2P,F6.1,A)') 'Wav. distribution shifted by',
	1			RAT_20-1.0,'% (compared to original)'
	WRITE(10,'(A,F7.3)') 'Possible cell length multiplier =',
	1			CELL_MULT/RAT_21
C
	RETURN
	END


C------------- Output details of individual merge & corrections ----------

	SUBROUTINE OUTPUT_MERGE_INFO(FILE_NAME,SCALE_F,SCALE_C)
C
C Output details on the intensity merge to the laue4_mrg.lis file
C
	CHARACTER FILE_NAME*(*)
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	CHARACTER STR_HKL_M*15
	INTEGER HKL_EQV(3),HKLM_EQV(4)
	REAL DCOUNTS_SEQ(NSEQ_MAX)
C
C Calculate the esds in COUNTS_SEQ() from the LSQ uncertainty in the fit
C
	CALL CALC_ALL_DCOUNTS_SEQ(DCOUNTS_SEQ)
C
C Open and write header for the merging list file
C
	WRITE(10,'(/,A)') 'Details of intensity merge written to '//FILE_NAME
	OPEN(UNIT=3,FILE=FILE_NAME,STATUS='UNKNOWN')
C Output header (including Twin or Modulation number if used)
	IF(IDATA_MODE .EQ. 1) THEN
	  WRITE(3,'(2A,/)') '   h   k   l   Calc   Obs   Sig   Difference',
	1		'   File  Wav  TTH   Xpix Ypix  Effic  Abs  Exti'
	ELSE
	  WRITE(3,'(2A,/)') '   h   k   l  m   Calc   Obs   Sig   Difference',
	1		'   File  Wav  TTH   Xpix Ypix  Effic  Abs  Exti'
	ENDIF
C
C Loop through sequences of reflections
C
	NWRITE=1
	DO ISEQ=1,NSEQ
C Add a divider between any non-empty sequences
	  IF(NWRITE .GT. 0) WRITE(3,*)
C
C Output details of individual reflections in this sequence
	  NWRITE=0
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
C Ignore outliers which are marked with ISIG() < 0
	    IF(ISIG(I) .LT. 0.0) CYCLE
C Calculate the information to output (and set upper/lower limits)
	    IOBS2 =MAX(-999999,MIN(9999999, NINT(IOBS(I) *SCALE_F) ))
	    ICALC2=MAX(-999999,MIN(9999999, NINT(ICALC(I) *SCALE_F) ))
	    ISIG2 =MIN(99999, NINT(ISIG(I) *SCALE_F) )
	    REL100=100.0 * (IOBS(I)-ICALC(I)) / MAX(0.1,ICALC(I))
	    IREL100=MAX(-99,MIN(99, NINT(REL100) ))
	    ESD100=100.0 * ISIG(I) / MAX(0.1,ICALC(I))
	    IESD100=MIN(99, NINT(ESD100) )
	    RSIGMA=MAX(-9.9,MIN(9.9, (IOBS(I)-ICALC(I))/ISIG(I) ))
C Calculate specific intensity factors
	    EFF_FACTOR=CALC_EFFIC_FACTOR(I)
	    ABS_FACTOR=CALC_ABS_FACTOR(WAVS(I),TTHS(I))
	    EXTI_X=CALC_EXTI_FACTOR(COUNTS_SEQ(ISEQ),WAVS(I),TTHS(I))
	    EXTI_FACTOR=1.0/SQRT(MAX(1E-4,1.0+EXTI_X))
C Make string with HKL and an optional twin/modulation number
	    IF(IDATA_MODE .EQ. 1) THEN
	      WRITE(STR_HKL_M,'(3I4)') (HKLMS(K,I),K=1,3)
	    ELSE
	      WRITE(STR_HKL_M,'(3I4,I3)') (HKLMS(K,I),K=1,4)
	    ENDIF
C Output the information for one reflection
	    NWRITE=NWRITE+1
	    WRITE(3,'(A,2I7,I5,I4,A,I2,A,F5.1, I5,F6.3,F6.1,1X,2I5, F7.3,2F6.3)',
	1														IOSTAT=IDUMMY)
	2			TRIM(STR_HKL_M),ICALC2,IOBS2,ISIG2,IREL100,'(',IESD100,')%',
	3			RSIGMA,IFILE_NUM(IFILES(I)),WAVS(I),TTHS(I),
	4			(NINT(XY_PIX(K,I)),K=1,2),EFF_FACTOR,ABS_FACTOR,EXTI_FACTOR
	  ENDDO
C
C Output summary information for merge of the sequence
C
C Make string with HKL and twin/modulation number for equivalent HKL
	  I=IFIRST(ISEQ)
	  IF(IDATA_MODE .EQ. 1) THEN
	    CALL CALC_EQUIV_HKL(HKL_EQV, HKLMS(1,I) )
	    WRITE(STR_HKL_M,'(3I4)') HKL_EQV
	  ELSEIF(IDATA_MODE .EQ. 2) THEN
	    CALL CALC_EQUIV_HKL(HKL_EQV, HKLMS(1,I) )
	    WRITE(STR_HKL_M,'(3I4,I3)') HKL_EQV,HKLMS(4,I)
        ELSE
	    CALL CALC_EQUIV_MOD(HKLM_EQV, HKLMS(1,I))
	    WRITE(STR_HKL_M,'(3I4,I3)') HKLM_EQV
	  ENDIF
C If just one reflection, output the corrected intensity
	  IF(NWRITE .EQ. 1) THEN
	    WRITE(3,'(A,2X,A,I6,I4)',IOSTAT=IDUMMY) TRIM(STR_HKL_M),'Corrected Intensity =',
	1		NINT(COUNTS_SEQ(ISEQ)*SCALE_C),	NINT(DCOUNTS_SEQ(ISEQ)*SCALE_C)
C If more than one reflection, output merge statistics as well
	  ELSE
	    CALL CALC_LSQ_MERGE_SEQ(R1,R1_ESD,R2,R2_ESD,GOOF,NSUM,
	1					ISEQ,ISEQ,-1.0)
	    IR2=MIN(99,NINT( R2 ))
	    IR2_ESD=MIN(99,NINT( R2_ESD ))
	    IR1=MIN(99,NINT( R1 ))
	    IR1_ESD=MIN(99,NINT( R1_ESD ))
	    WRITE(3,'(A,2X,A,I6,I4,2X,2(A,I3,A,I2,A,2X),A,F4.1)',IOSTAT=IDUMMY)
	1		TRIM(STR_HKL_M),'Corrected Intensity =',NINT(COUNTS_SEQ(ISEQ)*SCALE_C),
	2		NINT(DCOUNTS_SEQ(ISEQ)*SCALE_C),'wR2 =',IR2,' (',IR2_ESD,')%',
	3		'R =',IR1,' (',IR1_ESD,')%','GOOF =',GOOF
	  ENDIF
C
C Loop back for next sequence
	ENDDO
C
C Close output file and return
C
	CLOSE(UNIT=3)
C
	WRITE(10,'(2X,A,F6.2)') 'Raw intensities multiplied by',SCALE_F
	WRITE(10,'(2X,A,F6.2)') 'Corrected intensities multiplied by',SCALE_C
C
	RETURN
	END


	SUBROUTINE CALC_EQUIV_MOD(HKLM_EQV, HKLM)
C
C Returns the HKL indices and modulation number which are equivalent
C to the given values. This routine may not give a unique answer if
C the modulation list does not have matching positive and negative
C modulations. In this case it may just return the given values.
C As this routine is not used in the sorting of sequences this
C limitation shouldn't cause a problem.
C
	INTEGER HKLM(4),HKLM_EQV(4)
C
	REAL MOD_VEC
	COMMON /MODULATE_COM/ NMOD_VEC,NMOD_IDX,MOD_VEC(3,100),MOD_IDX(100,100)
C
C Calculate equivalent to main reflection (ignores modulation)
C
	CALL CALC_EQUIV_HKL(HKLM_EQV,HKLM)
C
C If no modulation, return with equivalent HKL and M=0
C
	IF(HKLM(4) .EQ. 0) THEN
		HKLM_EQV(4)=0
		RETURN
	ENDIF
C
C Find a modulation for HKLM_EQV main reflection that is equivalent to HKLM()
C NB: May not be possible if any modulation does not have a negative
C
	DO I=1,NMOD_IDX
		HKLM_EQV(4)=I
		IF(ORDER_SEQ(HKLM_EQV, HKLM) .EQ. 0.0) RETURN
	ENDDO
C
C Give up and use original HKLM
C
	DO I=1,4
		HKLM_EQV(I)=HKLM(I)
	ENDDO
C
	RETURN
	END


C --------------- Output the extended information intensity file ------------

	SUBROUTINE OUTPUT_VERBOSE_INTS(FILE_NAME)
C
	CHARACTER*(*) FILE_NAME
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
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
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
	REAL LOR_CORR
C
C Open output file, write header and update log file
C
	WRITE(10,'(A)') 'Extended intensity information written to '//FILE_NAME
	OPEN(UNIT=3,FILE=FILE_NAME,STATUS='UNKNOWN')
C Output header (including Twin or Modulation number if used)
	IF(IDATA_MODE .EQ. 1) THEN
	  WRITE(3,'(2A,/)') '   H   K   L   Counts & esd    Wav    TTH    X      ',
	1		'Y   File  Factor  Loren  Scale  Wav   Eff   Ext   Abs'
	ELSE
	  WRITE(3,'(2A,/)') '   H   K   L  m    Counts & esd    Wav    TTH    X      ',
	1		'Y   File  Factor  Loren  Scale  Wav   Eff   Ext   Abs'
	ENDIF
C
C
	DO ISEQ=1,NSEQ
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    LOR_CORR=WAVS(I)**4 /( SIND(TTHS(I)/2) )**2
	    FIL_CORR=CALC_FILE_SCALE(IFILES(I))
	    WAV_CORR=CALC_WAV_SCALE(WAVS(I))
	    EFF_CORR=CALC_EFFIC_FACTOR(I)
	    EXT_FAC=CALC_EXTI_FACTOR(COUNTS_SEQ(ISEQ),WAVS(I),TTHS(I))
	    EXT_CORR=1.0 / SQRT( MAX(0.5,1.0+EXT_FAC) )
	    ABS_CORR=CALC_ABS_FACTOR(WAVS(I),TTHS(I))
	    FACTOR=LOR_CORR*FIL_CORR*WAV_CORR*EFF_CORR*EXT_CORR*ABS_CORR
C Output factors (including Twin or Modulation number if used)
	    IF(IDATA_MODE .EQ. 1) THEN
	      WRITE(3,'(3I4,2I7, F8.3,F7.2,2F7.1,I4, F9.3,F7.1,7F6.3)',IOSTAT=IDUM)
	1				(HKLMS(K,I),K=1,3),NINT(IOBS(I)),NINT(ISIG(I)),
	2				WAVS(I),TTHS(I),(XY_PIX(K,I),K=1,2),IFILE_NUM(IFILES(I)),
	3				FACTOR,LOR_CORR,FIL_CORR,WAV_CORR,EFF_CORR,EXT_CORR,ABS_CORR
	    ELSE
	      WRITE(3,'(3I4,I3,2I7, F8.3,F7.2,2F7.1,I4, F9.3,F7.1,7F6.3)',IOSTAT=IDUM)
	1				(HKLMS(K,I),K=1,4),NINT(IOBS(I)),NINT(ISIG(I)),
	2				WAVS(I),TTHS(I),(XY_PIX(K,I),K=1,2),IFILE_NUM(IFILES(I)),
	3				FACTOR,LOR_CORR,FIL_CORR,WAV_CORR,EFF_CORR,EXT_CORR,ABS_CORR
	    ENDIF
	  ENDDO
	ENDDO
	CLOSE(UNIT=3)
C
	RETURN
	END


C --------------- Output the efficiency correction file ------------

	SUBROUTINE OUTPUT_EFFIC_FILE(FILE_NAME)
C
	CHARACTER FILE_NAME*(*)
C
	COMMON /EFF_CORR_COM/ IEFF_OPT,EFF_POLY(10)
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	REAL EFF(0:20)
C
C Don't output the file unless we are in verbose mode
C
	IF(IVERBOSE .EQ. 0) RETURN
C
	WRITE(10,'(2A)'),'Efficiency corrections written to ',FILE_NAME
	OPEN(UNIT=3,FILE=FILE_NAME,STATUS='UNKNOWN')
C
	CALL GET_NUMXY(NUMX,NUMY)
	WRITE(3,'(A,/)') 'Total efficiency correction versus Y & wavelength'
	WRITE(3,'(6X,21I8)') (IY,IY=0,NUMY,NUMY/20)
	X=NUMX/2.0
	DO WAV=0.80,1.801,0.05
	  DO IY=0,20
	    Y=IY*NUMY/20
	    EFF(IY)=CALC_EFFIC_FIXED(WAV,Y) * CALC_EFFIC_REFINE(Y)
	  ENDDO
	  WRITE(3,'(F6.2,21F8.4)',IOSTAT=IDUMMY) WAV,EFF
	ENDDO
C
	WRITE(3,'(/,A,/)') 'Refined efficiency correction versus Y & wavelength'
	WRITE(3,'(6X,21I8)') (IY,IY=0,NUMY,NUMY/20)
	DO WAV=0.80,1.801,0.05
	  DO IY=0,20
	    EFF(IY)=CALC_EFFIC_REFINE(FLOAT(IY*NUMY/20))
	  ENDDO
	  WRITE(3,'(F6.2,21F8.4)',IOSTAT=IDUMMY) WAV,EFF
	ENDDO
C
	IF(IEFF_OPT .NE. 0) THEN
	  WRITE(3,'(/,A,/)') 'Polynomial efficiency correction versus Y'
	  WRITE(3,'(21I8)') (IY,IY=0,NUMY,NUMY/20)
	  DO IY=0,20
	    EFF(IY)=CALC_EFFIC_POLY(FLOAT(IY*NUMY/20), EFF_POLY,IEFF_OPT)
	  ENDDO
	  WRITE(3,'(21F8.4)',IOSTAT=IDUMMY) EFF
	ENDIF
C
	WRITE(3,'(/,A,/)') 'Fixed efficiency correction versus Y & wavelength'
	WRITE(3,'(6X,21I8)') (IY,IY=0,NUMY,NUMY/20)
	DO WAV=0.80,1.801,0.05
	  DO IY=0,20
	    EFF(IY)=CALC_EFFIC_FIXED(WAV,FLOAT(IY*NUMY/20))
	  ENDDO
	  WRITE(3,'(F6.2,21F8.4)',IOSTAT=IDUMMY) WAV,EFF
	ENDDO
C
	WRITE(3,'(/,A,/)') 'Unnormalized IP efficiency versus Y & wavelength'
	WRITE(3,'(6X,21I8)') (IY,IY=0,NUMY,NUMY/20)
	DO WAV=0.80,1.801,0.05
	  DO IY=0,20
	    EFF(IY)=CALC_IP_EFFIC(WAV,FLOAT(IY*NUMY/20))
	  ENDDO
	  WRITE(3,'(F6.2,21F8.4)',IOSTAT=IDUMMY) WAV,EFF
	ENDDO
C
	CLOSE(UNIT=3)
C
	RETURN
	END


C ------------------- Output the parts of the CIF file ----------------

	SUBROUTINE WRITE_CIF_CELL
C
	COMMON /CELL_TYPE_COM/ ITYPE,IFRIEDEL,ITYPE_ORIG,ICEN,CELL(6)
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	CHARACTER*1 CEN(7)
	DATA CEN /"P","A","B","C","I","F","R"/
	CHARACTER*12 LATT(12)
	DATA LATT/"cubic       ","tetragonal  ","tetragonal  ","tetragonal  ",
	1	  "orthorhombic","hexagonal   ","trigonal    ","trigonal    ",
	2	  "monoclinic  ","monoclinic  ","monoclinic  ","triclinic   "/
C
C Crystal records
C
	WRITE(2,'(2A,/,2A,//,A)')
	1	"_symmetry_cell_setting        ",TRIM(LATT(ITYPE)),
	2	"_space_group.centring_type    ",CEN(ICEN),
	3	"#  Change the following if you don't use the neutron cell dimensions"
	WRITE(2,'(A,F7.2)')
	1	"_cell_length_a              ",CELL(1)*CELL_MULT,
	2	"_cell_length_b              ",CELL(2)*CELL_MULT,
	2	"_cell_length_c              ",CELL(3)*CELL_MULT
	WRITE(2,'(A,F7.1)')
	1	"_cell_angle_alpha           ",CELL(4),
	2	"_cell_angle_beta            ",CELL(5),
	3	"_cell_angle_gamma           ",CELL(6)
C
	
	FACTOR=1.0 - COSD(CELL(4))**2 - COSD(CELL(5))**2 - COSD(CELL(6))**2 
	2			+ 2.0*COSD(CELL(4))*COSD(CELL(5))*COSD(CELL(6))
	VOL=CELL(1)*CELL(2)*CELL(3)*SQRT(FACTOR)*CELL_MULT**3
	WRITE(2,'(A,F8.1)')
	1	"_cell_volume                ",VOL
C
	WRITE(2,'(A)')
	1	"_cell_measurement_radiation  'neutron'",
	2	"_cell_special_details",
	3	"'Cell obtained by the Laue method is inaccurate."//
	3				" Treat as indicative only.'"
C
	RETURN
	END


	SUBROUTINE WRITE_CIF_INSTRUM
C
	CHARACTER HOST_NAME*10, DATE*10
C
	WRITE(2,'(A)')
	1	"_diffrn_radiation_probe          'neutron'",
	2	"_diffrn_radiation_type           'neutron'",
	3	"_diffrn_source                   'nuclear reactor'",
	4	"_diffrn_source_details           'thermal neutrons, supermirror guides'"
C
	CALL GET_HOST_DATE(HOST_NAME,DATE)
	IF(HOST_NAME .EQ. 'KOALA') THEN
	  WRITE(2,'(A)')
	1	"_diffrn_source_type              'OPAL reactor, ANSTO, Lucas Heights, Australia'",
	2	"_diffrn_measurement_device_type  'KOALA'"
	ELSEIF(HOST_NAME .EQ. 'VIVALDI') THEN
	  WRITE(2,'(A)')
	1	"_diffrn_source_type              'High-Flux Reactor, ILL, Grenoble, France'",
	2	"_diffrn_measurement_device_type  'VIVALDI'"
	ELSEIF(HOST_NAME .EQ. 'IMAGINE') THEN
	  WRITE(2,'(A)')
	1	"_diffrn_source_type              'High Flux Isotope Reactor, ORNL, Oak Ridge, USA'",
	2	"_diffrn_measurement_device_type  'IMAGINE'"
	ELSEIF(HOST_NAME .EQ. 'CYCLOPS') THEN
	  WRITE(2,'(A)')
	1	"_diffrn_source_type              'High-Flux Reactor, ILL, Grenoble, France'",
	2	"_diffrn_measurement_device_type  'CYCLOPS'"
	ENDIF
C
	WRITE(2,'(A)')
	1	"_diffrn_measurement_device       'Laue image-plate diffractometer'",
	2	"_diffrn_measurement_method       'Laue'",
	3	"_diffrn_reflns_reduction_process 'LaueG (Piltz, 2018a)'"
C
	RETURN
	END


	SUBROUTINE WRITE_CIF_EXPERI
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	COMMON /CELL_TYPE_COM/ ITYPE,IFRIEDEL,ITYPE_ORIG,ICEN,CELL(6)
C
	DATA ICUB,ITETA,ITETB,ITETC,IORTH,IHEX,ITRIG,IRHOM,IMONA,
CCC	1		IMONB,IMONC,IDENT/1,2,3,4,5,6,7,8,9,10,11,12/
	1		IMONB,IMONC/1,2,3,4,5,6,7,8,9,10,11/
C
	CHARACTER SUMMARY*80
	INTEGER H_MIN,H_MAX
	REAL DCOUNTS_SEQ(NSEQ_MAX)
C
C Output a summary e.g. "12 images, exposure time of 10 to 60 minutes'
C
	CALL GET_EXPERI_SUMMARY(SUMMARY)
	WRITE(2,"(A,1H',A,1H')")
	1	"_diffrn_measurement_details  ",TRIM(SUMMARY)
C
C Calculate and output limits in HKL, WAV & TTH
C
	H_MIN=+999
	H_MAX=-999
	K_MIN=+999
	K_MAX=-999
	L_MIN=+999
	L_MAX=-999
	WAV_MIN=1E6
	TH_MAX=0.0
	D_MIN=1E6
	DO I=1,ILAST(NSEQ)
	  H_MIN=MIN(H_MIN,HKLMS(1,I))
	  H_MAX=MAX(H_MAX,HKLMS(1,I))
	  K_MIN=MIN(K_MIN,HKLMS(2,I))
	  K_MAX=MAX(K_MAX,HKLMS(2,I))
	  L_MIN=MIN(L_MIN,HKLMS(3,I))
	  L_MAX=MAX(L_MAX,HKLMS(3,I))
	  WAV_MIN=MIN(WAV_MIN,WAVS(I))
	  TH_MAX=MAX(TH_MAX,TTHS(I)/2.0)
	  D_MIN=MIN(D_MIN, 0.5*WAVS(I)/SIND(TTHS(I)/2.0) )
	ENDDO
C
	WRITE(2,'(A,I6)')
	1	"_diffrn_reflns_limit_h_min          ",H_MIN,
	2	"_diffrn_reflns_limit_h_max          ",H_MAX,
	3	"_diffrn_reflns_limit_k_min          ",K_MIN,
	4	"_diffrn_reflns_limit_k_max          ",K_MAX,
	5	"_diffrn_reflns_limit_l_min          ",L_MIN,
	6	"_diffrn_reflns_limit_l_max          ",L_MAX
C
	WRITE(2,'(A,F5.2)')
	1	"_diffrn_radiation_wavelength         ",WAV_MIN,
	2	"_diffrn_reflns_theta_full            ",TH_MAX-0.1,
	3	"_diffrn_reflns_theta_max             ",TH_MAX,
	4	"_diffrn_reflns_resolution_max        ",1.0/D_MIN
C
C
	VOL=CELL(1)*CELL(2)*CELL(3)*
	1		SQRT( 1.0 - COSD(CELL(4))**2	- COSD(CELL(5))**2 - COSD(CELL(6))**2 
	2					+ 2.0*COSD(CELL(4))*COSD(CELL(5))*COSD(CELL(6)) )
	FULL=4.0*3.14159/3.0*VOL/D_MIN**3
	IF(ICEN .EQ. 2) FULL=FULL/2.0
	IF(ICEN .EQ. 3) FULL=FULL/2.0
	IF(ICEN .EQ. 4) FULL=FULL/2.0
	IF(ICEN .EQ. 5) FULL=FULL/4.0
	IF(ICEN .EQ. 6) FULL=FULL/2.0
	IF(ICEN .EQ. 7) FULL=FULL/3.0
C
	NEQV=1
	IF(ITYPE.EQ.IMONA .OR. ITYPE.EQ.IMONB .OR. ITYPE.EQ.IMONC) NEQV=2
	IF(ITYPE .EQ. IORTH) NEQV=4
	IF(ITYPE.EQ.ITETA .OR. ITYPE.EQ.ITETB .OR. ITYPE.EQ.ITETC) NEQV=8
	IF(ITYPE.EQ.ITRIG .OR. ITYPE.EQ.IRHOM) NEQV=3
	IF(ITYPE .EQ. IHEX) NEQV=6
	IF(ITYPE .EQ. ICUB) NEQV=12
	IF(IFRIEDEL .EQ. 1) NEQV=2*NEQV
C
	FRAC=NSEQ*NEQV/FULL
	WRITE(2,'(A,F6.3)')
	1	"_diffrn_measured_fraction_theta_full",FRAC,
	2	"_diffrn_measured_fraction_theta_max ",FRAC
C
C Calculate various statistics from the data set
C
	CALL CALC_ALL_DCOUNTS_SEQ(DCOUNTS_SEQ)
C
	NTOTAL=0
	NSEQ_GT=0
C
	SUM_I1=0.0
	SUM_D1=0.0
	SUM_S1=0.0
C
	DO ISEQ=1,NSEQ
	  IF(COUNTS_SEQ(ISEQ) .GT. 2.0*DCOUNTS_SEQ(ISEQ)) NSEQ_GT=NSEQ_GT+1
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
C Ignore any outlier reflns
	    IF(ISIG(I) .LE. 0.0) CYCLE
	    NTOTAL=NTOTAL+1
	    SUM_I1=SUM_I1+IOBS(I)
	    SUM_D1=SUM_D1+ABS(IOBS(I)-ICALC(I))
	    SUM_S1=SUM_S1+ABS(ISIG(I))
	  ENDDO
	ENDDO
C
	R_MERG=SUM_D1/SUM_I1
	R_ESD=SUM_S1/SUM_I1
C
	WRITE(2,'(2(A,F7.4,/),3(A,I6,/),A,/,A)')
	1	"_diffrn_reflns_av_R_equivalents     ",R_MERG,
	2	"_diffrn_reflns_av_unetI/netI        ",R_ESD,
	3	"_diffrn_reflns_number               ",NTOTAL,
	4	"_reflns_number_total                ",NSEQ,
	5	"_reflns_number_gt                   ",NSEQ_GT,
	6	"_reflns_threshold_expression      '>2sigma(I)'",
	7	"_diffrn_standards_number                 0",
	8	"_exptl_absorpt_coefficient_mu          0.0"
C
C Output absorption correction method, or not
C
	IF( LABS_CORR ) THEN
	  WRITE(2,'(A)')
	1	"_exptl_absorpt_correction_type				'empirical'",
	2	"_exptl_absorpt_process_details",
	3	"'multiwavelength empirical method, LaueG (Piltz, 2018a)'"
	ELSE
	  WRITE(2,'(A)') "_exptl_absorpt_correction_type       'none'"
	ENDIF
C
	RETURN
	END


	SUBROUTINE WRITE_CIF_LAUEG
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	COMMON /WAV_CORR_COM/ IWAV_OPT,NWAV,WAV_GAMMA,WAV_MIN,WAV_STEP,WAV_PAR(20)
	COMMON /EFF_CORR_COM/ IEFF_OPT,EFF_POLY(10)
	COMMON /EXTI_CORR_COM/ IEXTI_OPT,EXTI
	COMMON /ABS_CORR_COM/ IABS_OPT,UR0,UR_LIN
C
	COMMON /CELL_TYPE_COM/ ITYPE,IFRIEDEL,ITYPE_ORIG,ICEN,CELL(6)
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
	COMMON /HARM_INFO_COM/ HARM_MAX,NHKLS_1ESD
C
C Begin the multiline string "_reflns_special_details"
C
	WRITE(2,'(A,/,1H;)') "_reflns_special_details"
	WRITE(2,'(A)') "LAUEG SPECIFIC INFORMATION:"
C
C Calculate and output limits in wav & d-space of the data
C
	WAV_MIN=1E6
	WAV_MAX=0.0
	D_MIN=1E6
	DO I=1,ILAST(NSEQ)
	  WAV_MIN=MIN(WAV_MIN,WAVS(I))
	  WAV_MAX=MAX(WAV_MAX,WAVS(I))
	  D_MIN=MIN(D_MIN, 0.5*WAVS(I)/SIND(TTHS(I)/2.0) )
	ENDDO
C
	WRITE(2,'(2X,3(A,F6.3))') "HKL generation: wav =",WAV_MIN,
	1			" to",WAV_MAX,", d-space >",D_MIN
C
C Output the corrections applied to the data
C
	IF( LFILE_CORR ) THEN
	  WRITE(2,'(2X,A)') 'Image scale factors: refined'
	ELSE
	  WRITE(2,'(2X,A)') 'Image scale factors: fixed to 1'
	ENDIF
C
	IF( .NOT.LWAV_CORR ) THEN
	  WRITE(2,'(2X,A)') 'Wavelength spectra model: not refined'
	ELSE IF(IWAV_OPT .EQ. 1) THEN
	  WRITE(2,'(2X,A,I3,A)') 'Wavelength spectra model: ',NWAV,' point spline'
	ELSE IF(IWAV_OPT .EQ. 2) THEN
	  WRITE(2,'(2X,A,I3,A,F6.1)') 'Wavelength spectra model: ',NWAV,
	1			' point spline + non-param, gamma =',WAV_GAMMA
	ELSE
	  STOP 'BUG(log_options): Invalid IWAV_OPT'
	ENDIF
C
	IF( LEFF_CORR ) THEN
	  IF(IEFF_OPT .EQ. 1) WRITE(2,'(2X,A)') 'Efficiency model: Linear in Y'
	  IF(IEFF_OPT .EQ. 2) WRITE(2,'(2X,A)') 'Efficiency model: Quadratic in Y'
	  IF(IEFF_OPT .GT. 2) WRITE(2,'(2X,A,I2,A)')
	1			'Efficiency model: Polynomial(',IEFF_OPT,') in Y'
	ELSE
	  WRITE(2,'(2X,A)') 'Efficiency model: not refined'
	ENDIF
C
	IF( LEXTI_CORR ) THEN
	  IF(IEXTI_OPT .EQ. 1) WRITE(2,'(2X,A)') 'Extinction correction: Zach. Type 1 / SHELX'
	  IF(IEXTI_OPT .EQ. 2) WRITE(2,'(2X,A)') 'Extinction correction: Zach. Type 2'
	  IF(IEXTI_OPT .EQ. 3) WRITE(2,'(2X,A)') 'Extinction correction: Zach. Type 1 (+ B&C)'
	  IF(IEXTI_OPT .EQ. 4) WRITE(2,'(2X,A)') 'Extinction correction: B&C Type 1,G'
	  IF(IEXTI_OPT .EQ. 5) WRITE(2,'(2X,A)') 'Extinction correction: B&C Type 1,L'
	  IF(IEXTI_OPT .EQ. 6) WRITE(2,'(2X,A)') 'Extinction correction: B&C Type 2'
	ELSE
	  WRITE(2,'(2X,A)') 'Extinction correction: off'
	ENDIF
C
	IF( LHARM_CORR ) THEN
	  WRITE(2,'(2X,A)') 'Harmonics (wav/2) correction: on'
	ELSE
	  WRITE(2,'(2X,A)') 'Harmonics (wav/2) correction: off'
	ENDIF
C
	IF( LABS_CORR ) THEN
	  WRITE(2,'(2X,2A)') 'Absorption correction: ',
	1		'Isotropic model (near-spherical sample)'
	  IF(IABS_OPT .EQ. 1) THEN
	    WRITE(2,'(4X,A)') 'Absorption model: non-hydrogenous'
	  ELSE IF(IABS_OPT .EQ. 2) THEN
	    WRITE(2,'(4X,A)') 'Absorption model: hydrogenous'
	  ELSE IF(IABS_OPT .EQ. 3) THEN
	    WRITE(2,'(4X,A)') 'Absorption model: general case'
	  ELSE
	    STOP 'BUG(log_options): Invalid IABS_OPT'
	  ENDIF
	ELSE
	  WRITE(2,'(2X,A)') 'Absorption correction: off'
	ENDIF
C
C Output if Friedel pairs are averaged, or not
C
	IF(IFRIEDEL .EQ. 1)  THEN
	  WRITE(2,'(2X,A)') "Friedel pairs are merged"
	ELSE
	  WRITE(2,'(2X,A)') "Friedel pairs not merged"
	ENDIF
C
C Output the effect of the corrections
C
	WRITE(2,'(2X,A)') 'Merge statistics using final weights:'
	CALL WRITE_LSQ_MERGE(2,4)
	IF( LHARM_CORR ) WRITE(2,'(2X,A,F5.1,A,I7,A)')
	1	'Max. harmonic correction =',MIN(999.9,HARM_MAX),' esds,',
	2	NHKLS_1ESD,' hkls corrected by > 1 esd'
	WRITE(2,'(2X,A)') 'Summary of refineable correction factors:'
	CALL WRITE_CORRECT_SUMMARY(2,4)
C
C Close off the multiline string "_reflns_special_details"
C
	WRITE(2,'(1H;)')
C
	RETURN
	END


	SUBROUTINE WRITE_CIF_SOFTWARE
C
	CHARACTER HOST_NAME*10,DATE*10
C
	CALL GET_HOST_DATE(HOST_NAME,DATE)
	IF(HOST_NAME .EQ. 'KOALA') THEN
		WRITE(2,'(A)')
	1		"_computing_data_collection  'MAATEL/ANSTO control program'"
	ELSE
		WRITE(2,'(A)')
	1		"_computing_data_collection  'MAATEL control program'"
	ENDIF
C
	WRITE(2,'(A)')
	1	"_computing_cell_refinement  'LaueG (Piltz, 2018b)'",
	2	"_computing_data_reduction",
	3	"'argonne_boxes (Wilkinson et al., 1988) & LaueG (Piltz, 2018b)'"
C
	RETURN
	END


	SUBROUTINE WRITE_CIF_REFERENCES
C
C References records
C
	WRITE(2,'(A,/,1H;,/,2(A,/),/,3(A,/),1H;)')
	1	"_publ_section_references",
	2	"Laue4:",
	3	"R.O. Piltz, J. Appl. Cryst. (2018). 51, 635-645.",
	4	"https://doi.org/10.1107/S1600576718005058",
	5	"LaueG:",
	6	"R.O. Piltz, J. Appl. Cryst. (2018). 51, 963-965.",
	7	"https://doi.org/10.1107/S1600576718005046",
	8	"argonne_boxes:",
	9	"C. Wilkinson, H. W. Khamis, R. F. D. Stansfield and G. J. McIntyre (1988),",
	1	"J. Appl. Cryst., 21, 471."
C
	RETURN
	END


	SUBROUTINE GET_EXPERI_SUMMARY(SUMMARY)
C
	CHARACTER SUMMARY*(*)
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
	CHARACTER*2000 COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,ADC_INFO
	COMMON /TIFFINFO_COM/ COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,
	1			IOFF_DATA,NUMX,NUMY,IXYRES,ISTARTX,ISTARTY,ISPEED,
	2			EXPOSE_TIME,EXPOSE_PHI,TEMP_BEGIN,TEMP_END,TEMP_MIN,TEMP_MAX,
	3			IADC,ADC_ADD,ADC_DIV,ADC_INFO,IAPER_1,IAPER_2,IKAPPA
C
	CHARACTER TIFF_NAME*132
C
C Load defaults values to return on errors
C
	WRITE(SUMMARY,"(I3,A)") NFILES," images with unknown exposure times"
C
C Read the TIF headers for info to make up the summary
C
	EXPO_MIN=1E6
	EXPO_MAX=0.0
	DO IFILE=1,NFILES
	  CALL MAKE_FILE_NAME(TIFF_NAME, IFILE_NUM(IFILE),'.tif')
	  CALL READ_TIFF_INFO(ISTATUS, TRIM(TIFF_NAME),1)
C If file missing or invalid, jump to error messages
	  IF(ISTATUS .EQ. -1) GOTO 900
	  IF(ISTATUS .EQ. 1) GOTO 910
C Update info using header values
	  EXPO_MIN=MIN(EXPO_MIN,EXPOSE_TIME)
	  EXPO_MAX=MAX(EXPO_MAX,EXPOSE_TIME)
	ENDDO
C
	IEXPO_MIN=NINT(EXPO_MIN/60.0)
	IEXPO_MAX=NINT(EXPO_MAX/60.0)
	IF(IEXPO_MIN .EQ. IEXPO_MAX) THEN
		WRITE(SUMMARY,"(I3,A,I4,A)") NFILES,
	1		" images, exposure time",IEXPO_MIN," minutes"
	ELSE
		WRITE(SUMMARY,"(I3,A,I4,A,I4,A)") NFILES,
	1		" images, exposure time",IEXPO_MIN," to",IEXPO_MAX," minutes"
	ENDIF
	SUMMARY=ADJUSTL(SUMMARY)
	RETURN
C
C On error, complain and return with the dummy values
C
900	WRITE(10,'(/,3A)') 'Unable to read ',TRIM(TIFF_NAME),
	1			  ', the CIF file will contain dummy information'
	RETURN
910	WRITE(10,'(/,3A)') TRIM(TIFF_NAME),' appears invalid',
	1			  ', the CIF file will contain dummy information'
	RETURN
	END
