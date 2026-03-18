C =========== Read and prepare data for refinement or correction ===========

C----- Top level routines -----------------------------------------------
c	SUBROUTINE PREP_REFINE_DATA
c	SUBROUTINE PREP_CORRECT_DATA
C----- Helper routines for above ----------------------------------------
c	SUBROUTINE PREP_TWIN_REFS
c	SUBROUTINE INFO_MOD_REFS
C----- Find sequences of equivalent HKLs --------------------------------
c	SUBROUTINE SORT_HKLS
c	SUBROUTINE FIND_EQUIV_SEQ
C----- Pruning the data lists -------------------------------------------
c	SUBROUTINE PRUNE_EQUIV_SEQ(NSEQ_MIN)
c	SUBROUTINE PRUNE_WAV_REFS(WAV_MIN,WAV_MAX)
c	SUBROUTINE PRUNE_XY_REFS(XLO,XHI,YLO,YHI)
c	SUBROUTINE PRUNE_TWIN_REFS
c	SUBROUTINE PRUNE_MOD_REFS
c	SUBROUTINE PRUNE_WEAK_SEQ(FRAC)
C----- Shuffling reflection data arrays ---------------------------------
c	SUBROUTINE APPEND_DATA_ITEM(IOUT, IIN)
c	SUBROUTINE REORDER_DATA_ITEMS(ITAGS)
C------------------------------------------------------------------------


C----- Top level routines -----------------------------------------------

	SUBROUTINE PREP_REFINE_DATA
C
C  Read and prepare the reduced data set for refinement of parameters
C
      COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	WRITE(10,'(/,A,/)') '====== Read data files and prepare refinement data ======'
	CALL READ_DATA_FILES(WAV_LO,CELL_MULT)
	CALL OUTPUT_DEBUG_DATA
C
C Removes any spots with wavelengths outside WAV_LO ... WAV_HI
C
	CALL PRUNE_WAV_REFS(WAV_LO,WAV_HI)
	CALL OUTPUT_DEBUG_DATA
C
C Remove any spots with X,Y outside of limits
C
	CALL PRUNE_XY_REFS(X_LO,X_HI, Y_LO,Y_HI)
	CALL OUTPUT_DEBUG_DATA
C
C Log twin information, and prune unused twin data
C
	CALL PREP_TWIN_REFS
	CALL PRUNE_TWIN_REFS
	CALL OUTPUT_DEBUG_DATA
C
C Log modulation information, and remove all satellite reflections
C
	CALL INFO_MOD_REFS
	CALL PRUNE_MOD_REFS
	CALL OUTPUT_DEBUG_DATA
C
C Sort reflections into "sequences" of consecutive HKLs that are equivalent.
C Uses Twin/Modulation numbers if present.
C
	CALL SORT_HKLS
	CALL OUTPUT_DEBUG_DATA
C
C Identify the sequences of equivalents
C
      CALL FIND_EQUIV_SEQ
	CALL OUTPUT_DEBUG_SEQ
C
C Remove sequences (and reflection data) if less than NSEQ_MIN items
C
	CALL PRUNE_EQUIV_SEQ( MAX(2,NSEQ_MIN) )
	CALL OUTPUT_DEBUG_DATA
	CALL OUTPUT_DEBUG_SEQ
C
C Remove a sequence if its average intensity corrected for the wavelength
C distribution is ranked in the lower SEQ_FRAC proportion of sequences.
C
	CALL PRUNE_WEAK_SEQ(SEQ_FRAC)
	CALL OUTPUT_DEBUG_DATA
	CALL OUTPUT_DEBUG_SEQ
C
	RETURN
      END


	SUBROUTINE PREP_CORRECT_DATA
C
C  Read and prepare the full data set for correction using refined parameters
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
      COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	PRINT '(/,1X,A)','Reread data files for final correction of intensities'
	WRITE(10,'(/,A,/)') '====== Reread data files and prepare for corrections ======'
C
	CALL READ_DATA_FILES(WAV_LO,CELL_MULT)
	CALL OUTPUT_DEBUG_DATA
C
C Removes any spots with wavelengths < WAV_LO, or > WAV_HI
C
	CALL PRUNE_WAV_REFS(WAV_LO,WAV_HI)
	CALL OUTPUT_DEBUG_DATA
C
C Remove any spots with X,Y outside of limits
C
	CALL PRUNE_XY_REFS(X_LO,X_HI, Y_LO,Y_HI)
	CALL OUTPUT_DEBUG_DATA
C
C Log twin information, and prune unused twin data
C
	CALL PREP_TWIN_REFS
	CALL PRUNE_TWIN_REFS
	CALL OUTPUT_DEBUG_DATA
C
C Log modulation information
C
	CALL INFO_MOD_REFS
C
C Sort reflections into "sequences" of consecutive HKLs that are equivalent.
C Uses Twin/Modulation numbers if present.
C
	CALL SORT_HKLS
	CALL OUTPUT_DEBUG_DATA
C
C Identify the sequences of equivalents
C
      CALL FIND_EQUIV_SEQ
	CALL OUTPUT_DEBUG_SEQ
C
C Remove sequences (and reflection data) if less than NSEQ_MIN items
C
	CALL PRUNE_EQUIV_SEQ(NSEQ_MIN)
	CALL OUTPUT_DEBUG_DATA
	CALL OUTPUT_DEBUG_SEQ
C
C Output how many reflections left to be corrected
C
	PRINT '(I7,A)',NDATA,' reflections suitable for correction'
C
	RETURN
      END


C----- Helper routines for above ----------------------------------------

	SUBROUTINE PREP_TWIN_REFS
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
      LOGICAL LTWIN1,LTWIN2,LTWIN11,LTWIN12,LRATIO
      COMMON /TWINS_CORR_COM/ ITWIN_OPT,TWIN_RATIO,LTWIN1,LTWIN2,LTWIN11,LTWIN12,LRATIO
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	IF(IDATA_MODE .NE. 2) RETURN
C
C Count types of twin reflections, and remove any singletons
C
	NSUM1=0
	NSUM2=0
      NSUM3=0
      NSUM11=0
      NSUM12=0
      N=0
      DO I=1,NDATA
	  IF(HKLMS(4,I) .EQ. 1) THEN
          NSUM1=NSUM1+1
	    CALL APPEND_DATA_ITEM(N,I)
	  ELSEIF(HKLMS(4,I) .EQ. 2) THEN
          NSUM2=NSUM2+1
	    CALL APPEND_DATA_ITEM(N,I)
	  ELSEIF( (HKLMS(4,I) .EQ. 11) .AND. (HKLMS(4,I+1) .EQ. 12) ) THEN
          NSUM3=NSUM3+1
	    CALL APPEND_DATA_ITEM(N,I)
	    CALL APPEND_DATA_ITEM(N,I+1)
          HKLMS(4,I+1)=-1
        ELSEIF(HKLMS(4,I) .EQ. 11) THEN
          NSUM11=NSUM11+1
        ELSEIF(HKLMS(4,I) .EQ. 12) THEN
          NSUM12=NSUM12+1
        ENDIF
      ENDDO
      NDATA=N
C
C Output tally to log file
C
      WRITE(10,'(/,A)') 'Intensity data for TWINNED reflections'
      WRITE(10,'(I7,A)') NSUM1,' separable intensities for Twin 1'
      WRITE(10,'(I7,A)') NSUM2,' separable intensities for Twin 2'
      WRITE(10,'(I7,A)') 2*NSUM3,' overlapped intensities as Twin pairs'
      WRITE(10,'(I7,A)') NSUM11+NSUM12,' incomplete Twin pairs (rejected)'
C
      RETURN
	END


	SUBROUTINE INFO_MOD_REFS
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
      IF(IDATA_MODE .NE. 3) RETURN
C
C Output modulation info to log file
C
	NSUM1=0
	DO I=1,NDATA
	  IF(HKLMS(4,I) .EQ. 1) NSUM1=NSUM1+1
      ENDDO
      WRITE(10,'(/,A)') 'Intensity data for MODULATED reflections'
      WRITE(10,'(I7,A)') NSUM1,' main reflections'
      WRITE(10,'(I7,A)') NDATA-NSUM1,' satellite reflections'
C
	RETURN
      END


C----- Find sequences of equivalent HKLs --------------------------------

	SUBROUTINE SORT_HKLS
C
C Sort reflections in order of importance: sequence number;
C		twin number, if applicable; actual HKL; wavelength (reversed)
C The sort on actual HKL and wavelength is to make the printouts prettier.
C
C This sort creates "sequences" of consecutive equivalent and repeated HKLs. 
C
	EXTERNAL ORDER_REF
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
      INTEGER ITAGS(NDATA_MAX)
C
C Sort the reflections using a tagged sort
C
	CALL TAG_SORT(ORDER_REF,ITAGS,NDATA)
C
C Reorder the reflections according to ITAGS()
C
	CALL REORDER_DATA_ITEMS(ITAGS)
C
	RETURN
	END


	SUBROUTINE FIND_EQUIV_SEQ
C
C Find repeats of the same equivalent hkl and put into /SEQ_COM/
C Assumes SORT_HKLS() has already been used to sort data
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
C Loop through reflection data finding all sequences of equivalents
C
	NSEQ=1
	IFIRST(1)=1
	ILAST(1)=1
C
	IREF_LAST=1
	DO IREF=1,NDATA
C If reflection in current sequence, update ILAST() of the sequence
	  IF(ORDER_SEQ( HKLMS(1,IREF) , HKLMS(1,IREF_LAST) ) .EQ. 0.0) THEN
	    ILAST(NSEQ)=IREF
	  ELSE
C Else, start a new sequence
	    NSEQ=NSEQ+1
	    IF(NSEQ .GT. NSEQ_MAX) STOP 'ERROR: Too many groups of equivalents'
	    IFIRST(NSEQ)=IREF
	    ILAST(NSEQ)=IREF
	  ENDIF
C Update IREF_LAST
		IREF_LAST=IREF
	ENDDO
C
C Zero COUNTS_SEQ (for OUTPUT_DEBUG_SEQ)
C
	DO I=1,NSEQ
	  COUNTS_SEQ(I)=0.0
	ENDDO
C
	RETURN
	END



C----- Pruning the data lists -------------------------------------------

	SUBROUTINE PRUNE_EQUIV_SEQ(NSEQ_MIN)
C
C Remove sequences with less than NSEQ_MIN items, and their reflection data
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
C Remove sequences with less than NSEQ_MIN members
C
	I2=0
	DO I1=1,NSEQ
	  IF(ILAST(I1)-IFIRST(I1)+1 .GE. NSEQ_MIN) THEN
	    I2=I2+1
	    IFIRST(I2)=IFIRST(I1)
	    ILAST(I2)=ILAST(I1)
	  ENDIF
	ENDDO
	NSEQ_REJ=NSEQ-I2
	NSEQ=I2
C
C Prune the reflection lists to only those in a sequence
C
	NDATA2=0
	NMULT=0
	DO ISEQ=1,NSEQ
	  ISTART=NDATA2+1
C Copy lists from IFIRST(ISEQ):ILAST(ISEQ) to ISTART++
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    CALL APPEND_DATA_ITEM(NDATA2, I)
	  ENDDO
C Update indices to lists
	  IFIRST(ISEQ)=ISTART
	  ILAST(ISEQ)=NDATA2
C Count sequences with multiple equivalents
	  IF(ILAST(ISEQ) .GT. IFIRST(ISEQ)) NMULT=NMULT+1
	ENDDO
	NREJECT=NDATA-NDATA2
	NDATA=NDATA2
C
C Output info about equivalents
C
	WRITE(10,'(/,A,I5,A,I2,A)') 'HKLs sorted into',NSEQ,
	1		' groups of equivalents with at least',NSEQ_MIN,' members'
	WRITE(10,'(I7,1X,A)') NREJECT,'reflections rejected as not in a group'
	WRITE(10,'(I5,A)') NMULT,
	1	' groups have > 1 member and will be used for merge statistics'
C
	RETURN
	END


	SUBROUTINE PRUNE_WAV_REFS(WAV_MIN,WAV_MAX)
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	N=0
	DO I=1,NDATA
C
C If we accept the hkl, move its to the new list position.
C Otherwise, it will be overwritten in the list.
C
	  IF(WAVS(I).GE.WAV_MIN .AND. WAVS(I).LE.WAV_MAX) THEN
	    CALL APPEND_DATA_ITEM(N, I)
	  ENDIF
C
	ENDDO
C
	WRITE(10,'(I7,A)') NDATA-N,' rejected due to wavelength limits'
	NDATA=N
C
	RETURN
	END


	SUBROUTINE PRUNE_XY_REFS(XLO,XHI,YLO,YHI)
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	N=0
	DO I=1,NDATA
	  IF(XY_PIX(1,I).GT.XLO .AND. XY_PIX(1,I).LT.XHI .AND.
	1		XY_PIX(2,I).GT.YLO .AND. XY_PIX(2,I).LT.YHI   ) THEN
	    CALL APPEND_DATA_ITEM(N,I)
	  ENDIF
	ENDDO
C
	WRITE(10,'(I7,A)') NDATA-N,' rejected due to XY limits'
	NDATA=N
C
	RETURN
      END


	SUBROUTINE PRUNE_TWIN_REFS
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
      LOGICAL LTWIN1,LTWIN2,LTWIN11,LTWIN12,LRATIO
      COMMON /TWINS_CORR_COM/ ITWIN_OPT,TWIN_RATIO,LTWIN1,LTWIN2,LTWIN11,LTWIN12,LRATIO
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	IF(IDATA_MODE .NE. 2) RETURN
C
      IF( .NOT.LTWIN1 ) WRITE(10,'(2X,A)') 'Excluding separable intensities for Twin 1'
      IF( .NOT.LTWIN2 ) WRITE(10,'(2X,A)') 'Excluding separable intensities for Twin 2'
      IF( .NOT.LTWIN11 ) WRITE(10,'(2X,A)') 'Excluding the intensity sum of overlapped pairs'
      IF( .NOT.LTWIN12 ) WRITE(10,'(2X,A)') 'Excluding the intensity difference of overlapped pairs'
C
	N=0
      DO I=1,NDATA
        ITWIN=HKLMS(4,I)
	  IF( (ITWIN .EQ. 1) .AND. LTWIN1 ) CALL APPEND_DATA_ITEM(N, I)
	  IF( (ITWIN .EQ. 2) .AND. LTWIN2 ) CALL APPEND_DATA_ITEM(N, I)
	  IF( (ITWIN .EQ. 11) .AND. LTWIN11 ) CALL APPEND_DATA_ITEM(N, I)
	  IF( (ITWIN .EQ. 12) .AND. LTWIN12 ) CALL APPEND_DATA_ITEM(N, I)
      ENDDO
C
      WRITE(10,'(I7,A)') NDATA-N,' twin intensities excluded'
	NDATA=N
C
      RETURN
	END


	SUBROUTINE PRUNE_MOD_REFS
C
C Remove all satellites (for parameter refinement)
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
      IF(IDATA_MODE .NE. 3) RETURN
C
	N=0
	DO I=1,NDATA
	  IF(HKLMS(4,I) .EQ. 1) THEN
	    CALL APPEND_DATA_ITEM(N, I)
	  ENDIF
	ENDDO
	WRITE(10,'(A)') 'All modulation satellites removed'
	NDATA=N
C
	RETURN
      END


	SUBROUTINE PRUNE_WEAK_SEQ(FRAC)
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	EXTERNAL SORT_WEAK_FUNC
C
	COMMON /SORT_WEAK_COM/ CSEQ(NSEQ_MAX)
C
	INTEGER ITAGS(NSEQ_MAX)
C
C Make weighted averages of a rough estimate for COUNTS_SEQ()
C for each sequence. The estimates only correct intensities
C for the rough wavelength distribution and Lorentz factor.
C This should be sufficient to rank the relative intensities.
C
	DO ISEQ=1,NSEQ
	  SUM_C=0.0
	  SUM_D=0.0
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    FACTOR=CALC_WAV_ROUGH(WAVS(I)) * WAVS(I)**4 / ( SIND(TTHS(I)/2) )**2
	    C=COUNTS(I)/FACTOR
	    D=DCOUNTS(I)/FACTOR
	    SUM_C=SUM_C+C/D
	    SUM_D=SUM_D+1.0/D
	  ENDDO
	  CSEQ(ISEQ)=SUM_C/SUM_D
	ENDDO
C
C Do tagged sort of CSEQ()
C
	CALL TAG_SORT(SORT_WEAK_FUNC,ITAGS,NSEQ)
C
C Get CSEQ() value where FRAC of the sequences have a lower CSEQ()
C
	IMED=1+NINT(FRAC*(NSEQ-1))
	CSEQ_CUTOFF=CSEQ(ITAGS(IMED))
C
C Remove the sequences with a lower CSEQ() than CSEQ_CUTOFF
C Also identify which remaining sequence is closest to the cutoff.
C
	NSEQ2=0
	CSEQ_MIN=1E10
	DO ISEQ=1,NSEQ
	  IF(CSEQ(ISEQ) .GT. CSEQ_CUTOFF) THEN
	    NSEQ2=NSEQ2+1
	    IFIRST(NSEQ2)=IFIRST(ISEQ)
	    ILAST(NSEQ2)=ILAST(ISEQ)
	    IF(CSEQ(ISEQ) .LT. CSEQ_MIN) THEN
	      CSEQ_MIN=CSEQ(ISEQ)
	      ISEQ_MIN=NSEQ2
	    ENDIF
	  ENDIF
	ENDDO
C
C Prune the list to only those in a sequence
C
	N=0
	DO ISEQ=1,NSEQ2
	  ISTART=N+1
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    CALL APPEND_DATA_ITEM(N, I)
	  ENDDO
	  IFIRST(ISEQ)=ISTART
	  ILAST(ISEQ)=N
	ENDDO
C
C Log results about pruning the groups
C
	WRITE(10,'(A,I3,A)') 'Removing weakest',NINT(100.0*FRAC),
	1						'% of equivalent HKL groups'
	WRITE(10,'(2X,2(I6,A,I5,A))') NDATA,' refl. in',NSEQ,
	1				' groups reduced to',N,' refl. in',NSEQ2,' groups'
C
C Output specific information about the group just above the cutoff
C
	RSIG_MIN=+1E6
	RSIG_MAX=-1E6
	RSIG_AVE=0.0
	N1=0
	N3=0
	N10=0
	DO I=IFIRST(ISEQ_MIN),ILAST(ISEQ_MIN)
	  RSIG=COUNTS(I)/DCOUNTS(I)
	  RSIG_MIN=MIN(RSIG_MIN,RSIG)
	  RSIG_MAX=MAX(RSIG_MAX,RSIG)
	  RSIG_AVE=RSIG_AVE+RSIG
	  IF(RSIG .GE. 1.0) N1=N1+1
	  IF(RSIG .GE. 3.0) N3=N3+1
	  IF(RSIG .GE. 10.0) N10=N10+1
	ENDDO
	NGROUP=ILAST(ISEQ_MIN)-IFIRST(ISEQ_MIN)+1
	RSIG_AVE=RSIG_AVE/NGROUP
C
C Update the NSEQ & NDATA and output data used for refinement
C
	NSEQ=NSEQ2
	NDATA=N
	PRINT '(1X,2(I6,A))',NDATA,' reflections (',NSEQ,
	1		' equivalents) used for parameter refinement'
C
C Output to console and log file the weakest group of equivalents used
C
	WRITE(6,'(/,1X,A)') 'Weakest group of equivalents used in refinement:'
	WRITE(10,'(/,A)') 'Weakest group of equivalents used in refinement:'
	WRITE(6,'(6X,A,4(I3,A))') 'Equiv reflns:',NGROUP,' total,',
	1							N1,' >1 sig,',N3,' >3 sig,',N10,' >10 sig'
	WRITE(10,'(5X,A,4(I3,A))') 'Equiv reflns:',NGROUP,' total,',
	1							N1,' >1 sig,',N3,' >3 sig,',N10,' >10 sig'
C
	RETURN
	END


C----- Shuffling reflection data arrays ---------------------------------

	SUBROUTINE APPEND_DATA_ITEM(IOUT, IIN)
C
C Copy all arrays in DATA_COM & TWINS_COM from index IIN to index IOUT+1
C The incremented IOUT is returned
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
      IOUT=IOUT+1
C
	HKLMS(1,IOUT)=HKLMS(1,IIN)
	HKLMS(2,IOUT)=HKLMS(2,IIN)
	HKLMS(3,IOUT)=HKLMS(3,IIN)
	HKLMS(4,IOUT)=HKLMS(4,IIN)
	IFILES(IOUT)=IFILES(IIN)
	WAVS(IOUT)=WAVS(IIN)
	COUNTS(IOUT)=COUNTS(IIN)
	DCOUNTS(IOUT)=DCOUNTS(IIN)
	TTHS(IOUT)=TTHS(IIN)
	XY_PIX(1,IOUT)=XY_PIX(1,IIN)
	XY_PIX(2,IOUT)=XY_PIX(2,IIN)
C
	RETURN
	END


	SUBROUTINE REORDER_DATA_ITEMS(ITAGS)
C
C Reorder all data items according to ITAGS()
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	INTEGER ITAGS(NDATA_MAX),ISAVE(NDATA_MAX)
	REAL RSAVE(NDATA_MAX)
C HKLs and modulation/twin numbers
	DO I2=1,4
	  DO I=1,NDATA
	    ISAVE(I)=HKLMS(I2,ITAGS(I))
	  ENDDO
	  DO I=1,NDATA
	    HKLMS(I2,I)=ISAVE(I)
	  ENDDO
	ENDDO
C IFILES
	DO I=1,NDATA
	  ISAVE(I)=IFILES(ITAGS(I))
	ENDDO
	DO I=1,NDATA
	  IFILES(I)=ISAVE(I)
	ENDDO
C WAVS
	DO I=1,NDATA
	  RSAVE(I)=WAVS(ITAGS(I))
	ENDDO
	DO I=1,NDATA
	  WAVS(I)=RSAVE(I)
	ENDDO
C COUNTS
	DO I=1,NDATA
	  RSAVE(I)=COUNTS(ITAGS(I))
	ENDDO
	DO I=1,NDATA
	  COUNTS(I)=RSAVE(I)
	ENDDO
C DCOUNTS
	DO I=1,NDATA
	  RSAVE(I)=DCOUNTS(ITAGS(I))
	ENDDO
	DO I=1,NDATA
	  DCOUNTS(I)=RSAVE(I)
	ENDDO
C TTHS
	DO I=1,NDATA
	  RSAVE(I)=TTHS(ITAGS(I))
	ENDDO
	DO I=1,NDATA
	  TTHS(I)=RSAVE(I)
	ENDDO
C XY_PIX
	DO I2=1,2
	  DO I=1,NDATA
	    RSAVE(I)=XY_PIX(I2,ITAGS(I))
	  ENDDO
	  DO I=1,NDATA
	    XY_PIX(I2,I)=RSAVE(I)
	  ENDDO
	ENDDO
C
	RETURN
	END
