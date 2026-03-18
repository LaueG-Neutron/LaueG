C----- Detect if HKL(M) is in the rejects list --------------------------------
c	LOGICAL FUNCTION IS_REJECT_HKL(IHKL,ITWIN)
C----- Detect if a debug HKL for file number IFILE ----------------------------
c	LOGICAL FUNCTION IS_DEBUG_HKL(IFILE,IHKL)
C----- Read reject HKLs for "BASE_NAME" data file -----------------------------
c	SUBROUTINE READ_REJECT_HKLS(NREJ2, BASE_NAME,REJECT_NAME)
C----- Load debug HKLs (and file-numbers) -------------------------------------
c	SUBROUTINE LOAD_DEBUG_HKLS
C----- Read a line of "file-name plus hkl" for rejects/debug file -------------
c	SUBROUTINE READ_NAME_HKLM(NAME,HKLM,ISTATUS, NHKLM,IUNIT)
C----- Routines to output debug info ------------------------------------------
c	SUBROUTINE OUTPUT_DEBUG_DATA
c	SUBROUTINE OUTPUT_DEBUG_CALC
c	SUBROUTINE OUTPUT_DEBUG_SEQ
C------------------------------------------------------------------------------

C WARNING: DEBUG mode may not work correctly


C----- Detect if HKL(M) is in the rejects list --------------------------------

	LOGICAL FUNCTION IS_REJECT_HKL(IHKL,ITWIN)
C
C Used by READ_GEASC_FILE & READ_INT_FILE to ignore hkls on the rejects list
C NB: Rejects list is loaded for individual data files
C
	INTEGER IHKL(3)
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	PARAMETER (NREJ_MAX=10000)
	INTEGER HKLM_REJECT
	COMMON /REJECT_COM/ NREJ,HKLM_REJECT(4,NREJ_MAX)
C
	DO I=1,NREJ
		IF(IHKL(1) .NE. HKLM_REJECT(1,I)) CYCLE
		IF(IHKL(2) .NE. HKLM_REJECT(2,I)) CYCLE
		IF(IHKL(3) .NE. HKLM_REJECT(3,I)) CYCLE
		IF( (IDATA_MODE .GT. 1) .AND. 
	1				(ITWIN .NE. HKLM_REJECT(4,I)) ) CYCLE
C Everything matches, so return TRUE
	  IS_REJECT_HKL=.TRUE.
	  RETURN
	ENDDO
C
C No match, so return FALSE
	IS_REJECT_HKL=.FALSE.
	END

C ============ Detect if a debug HKL for file number IFILE  ===========

	LOGICAL FUNCTION IS_DEBUG_HKL(IFILE,IHKL)
C
C Return true if reflection is in the debug list
C NB: Twin/satellite numbers are ignored
C
	INTEGER IHKL(3)
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
	PARAMETER (NDEBUG_MAX=100)
	CHARACTER*80 NAMES_DEBUG(NDEBUG_MAX)
	INTEGER HKL_DEBUG(3,NDEBUG_MAX)
	COMMON /DEBUG_COM/ NAMES_DEBUG,HKL_DEBUG,NDEBUG
C
	CHARACTER*80 BASE_NAME
C
	DO I=1,NDEBUG
C Check if matches IHKL
	  IF(IHKL(1) .NE. HKL_DEBUG(1,I)) CYCLE
	  IF(IHKL(2) .NE. HKL_DEBUG(2,I)) CYCLE
	  IF(IHKL(3) .NE. HKL_DEBUG(3,I)) CYCLE
C Check if matches the base-name for IFILE
	  INUM=IFILE_NUM(IFILE)
	  CALL MAKE_FILE_NAME(BASE_NAME, INUM,'')
	  CALL UP2LOW_CASE(BASE_NAME, BASE_NAME)
	  IF(BASE_NAME .NE. NAMES_DEBUG(I)) CYCLE
C Everything matches, so return with TRUE
	  IS_DEBUG_HKL=.TRUE.
	  RETURN
	ENDDO
C
C Nothing matches, so return FALSE
C
	IS_DEBUG_HKL=.FALSE.
	RETURN
	END


C ============ Read reject HKLs for "BASE_NAME" data file  ===========

	SUBROUTINE READ_REJECT_HKLS(NREJ2, BASE_NAME,REJECT_NAME)
C
C Reads the reject file (REJECT_NAME) and creates a reject list of
C HKLM_REJECT(1..3/4,1..NREJ2) of hkl(m)s for the data file BASE_NAME.
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	CHARACTER BASE_NAME*(*),REJECT_NAME*(*)
C
	PARAMETER (NREJ_MAX=10000)
	INTEGER HKLM_REJECT
	COMMON /REJECT_COM/ NREJ,HKLM_REJECT(4,NREJ_MAX)
C
	CHARACTER*255 NAME,TEMP
	INTEGER HKLM(4)
C
C Try to open the rejection list file on UNIT=3
C
	IUNIT=3
	OPEN(UNIT=IUNIT,FILE=REJECT_NAME,STATUS='OLD',ERR=900)
C
	TEMP=REPEAT(' ',LEN(TEMP))
	READ(IUNIT,'(Q,A,//)') ILEN,TEMP(1:ILEN)
C
	IF(TEMP(1:26) .NE. 'LaueG Rejects File Version') THEN
	  WRITE(10,'(/,2A)') 'ERROR: Invalid header in ',REJECT_NAME
	  STOP 'ERROR: Invalid header in reject file'
	ELSE IF(TEMP(29:36) .NE. ', Option') THEN
	  WRITE(10,'(/,2A)') 'ERROR: Invalid header in ',REJECT_NAME
	  STOP 'ERROR: Invalid header in reject file'
	ENDIF
C
	READ(TEMP(27:38),'(I2,8X,I2)') IVER,IOPT
	IF(IVER.LT.1 .OR. IVER.GT.2) THEN
	  WRITE(10,'(/,A,I2,2A)') 'ERROR: Invalid header version=',
	1				IVER,' in ',REJECT_NAME
	  STOP 'ERROR: Invalid header version in reject file'
	ELSE IF(IOPT.LT.1 .OR. IOPT.GT.IVER) THEN
	  WRITE(10,'(/,A,I2,2A)') 'ERROR: Invalid header option=',
	1				IOPT,' in ',REJECT_NAME
	  STOP 'ERROR: Invalid header option in reject file'
	ENDIF
C
C Check that data & rejects file have the same data type
C NB: Assumes GET_DATA_MODE has already been run
C
	IF(MIN(2,IDATA_MODE) .NE. IOPT) THEN
	  IF(IOPT .EQ. 2) THEN
	    WRITE(10,'(/,A,I2,2A)') 'ERROR: Rejects file contains twin/satellite numbers'
	    STOP 'ERROR: Rejects file contains twins/satellite numbers'
	  ELSE
	    WRITE(10,'(/,A,I2,2A)') 'ERROR: Rejects file does not contain twin/satellite numbers'
	    STOP 'ERROR: Rejects file does not contain twins/satellite numbers'
	  ENDIF
	ENDIF
C
C Set NHKLM=3,4 for hkl or hklm indices
C
	NHKLM=IOPT+2
C
C Loop to read lines of "file-name plus hkl(m)"
C
	NREJ=0
100	CALL READ_NAME_HKLM(NAME,HKLM,ISTATUS, NHKLM,IUNIT)
C If EOF, jump out of loop
	IF(ISTATUS .EQ. -1) GOTO 1000
C If an error, complain and die
	IF(ISTATUS .NE. 0) STOP 'ERROR: Invalid data in reject file'
C
C Convert NAME and BASE_NAME to lowercase (put BASE_NAME in TEMP)
	CALL UP2LOW_CASE(NAME, NAME)
	CALL UP2LOW_CASE(TEMP, BASE_NAME)
C
C If NAME and BASE_NAME don't match, jump back to read new line
	IF(TEMP .NE. TRIM(NAME)) GOTO 100
C
C Everything is OK, so append to HKLM_REJECT() and jump back to read again
	NREJ=NREJ+1
	IF(NREJ .GT. NREJ_MAX) STOP 'ERROR: Rejection list is too large'
	DO I=1,4
		HKLM_REJECT(I,NREJ)=HKLM(I)
	ENDDO
	GOTO 100
C
C Successful completion
C
1000	CLOSE(UNIT=IUNIT)
	NREJ2=NREJ
	RETURN
C
C Failed
900	NREJ2=-1
	RETURN
	END


C ============ Load debug HKLs (and file-numbers) ===========

	SUBROUTINE LOAD_DEBUG_HKLS
C
C Reads the debug file and loads the list of base-names and hkls
C
	PARAMETER (NDEBUG_MAX=100)
	CHARACTER*80 NAMES_DEBUG(NDEBUG_MAX)
	INTEGER HKL_DEBUG(3,NDEBUG_MAX)
	COMMON /DEBUG_COM/ NAMES_DEBUG,HKL_DEBUG,NDEBUG
C
	CHARACTER*80 BASE_NAME
	INTEGER HKLM(4)
C
C Try to open the debug file on UNIT=3
C
	NDEBUG=0
	IUNIT=3
	OPEN(UNIT=IUNIT,FILE='laue4_debug.dat',STATUS='OLD',ERR=900)
C
	WRITE(10,'(A,/)') "DEBUG: Loading laue4_debug.dat"
C
C Loop to read lines of "base-name plus hkl"
C
	DO I=1,NDEBUG_MAX+1
C Read a line of "base-name plus hkl"
	  CALL READ_NAME_HKLM(BASE_NAME,HKLM,ISTATUS, 3,IUNIT)
C Check for EOF and any errors
	  IF(ISTATUS .EQ. -1) EXIT
	  IF(ISTATUS .NE. 0) STOP 'ERROR: Invalid data in debug file'
	  IF(I .GE. NDEBUG_MAX) STOP 'ERROR: Too many debug HKLs'
C Store base-name (in lower case) and hkls
	  CALL UP2LOW_CASE(NAMES_DEBUG(I),BASE_NAME)
	  HKL_DEBUG(1,I)=HKLM(1)
	  HKL_DEBUG(2,I)=HKLM(2)
	  HKL_DEBUG(3,I)=HKLM(3)
	ENDDO
	NDEBUG=I-1
C
900	RETURN
	END


C ======= Read a line of "file-name plus hkl" for rejects/debug file =======

	SUBROUTINE READ_NAME_HKLM(NAME,HKLM,ISTATUS, NHKLM,IUNIT)
C
C Reads a line of the rejects, or debug, file on unit IUNIT and
C trys to decode the line as NAME followed by HKL.
C NHKLM=3 for normal "hkl", =4 for twinned and modulated "hklm"
C ISTATUS = 0 Success
C          -1 EOF encountered
C           1 Invalid format
C
C Assumes "%-40s%4i%4i  ! %s\n" for HKL, "%-40s%4i%4i%4i  ! %s\n" for HKLM
C
	CHARACTER NAME*(*)
	INTEGER HKLM(4)
C
	CHARACTER LINE*132
C
C Read a line from file
C
1	READ(IUNIT,'(A)',ERR=999,END=1000) LINE
C
C If an empty or blank line, skip it and read the next line
C
	IF(LEN_TRIM(LINE) .EQ. 0) GOTO 1
C
C Find beginning of comments, error if not found
C
	ICOMM=INDEX(LINE,'  ! ')
	IF(ICOMM .EQ. 0) GOTO 999
C Calculate end of file name, error if too small
	IEND=ICOMM-4*NHKLM-1
	IF(IEND .LT. 1) GOTO 999
C Copy file name, and decode HLKM
	NAME=LINE(1:IEND)
	HKLM(4)=0
	READ(LINE((IEND+1):(ICOMM-1)),'(<NHKLM>I4)',ERR=999) (HKLM(K),K=1,NHKLM)
C Success, return with ISTATUS=0
	ISTATUS=0
	RETURN
C
C ERROR: print out bad line and return with ISTATUS=1
C
999	IF(NHKLM .EQ. 3) THEN
		PRINT '(1X,A)','Expected "file-name h k l  ! ", instead of:'
	ELSE
		PRINT '(1X,A)','Expected "file-name h k l m  ! ", instead of:'
	ENDIF
	PRINT '(1X,A)',LINE(1:LEN_LINE)
	ISTATUS=1
	RETURN
C
C Reached EOF, return with ISTATUS=-1
C
1000	ISTATUS=-1
	RETURN
C
	END


C ============ Routines to output debug info  ===========

	SUBROUTINE OUTPUT_DEBUG_DATA
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
	PARAMETER (NDEBUG_MAX=100)
	CHARACTER*80 NAMES_DEBUG(NDEBUG_MAX)
	INTEGER HKL_DEBUG(3,NDEBUG_MAX)
	COMMON /DEBUG_COM/ NAMES_DEBUG,HKL_DEBUG,NDEBUG
C
	LOGICAL IS_DEBUG_HKL,LMATCH
C
	IF(NDEBUG .EQ. 0) RETURN
C
	LMATCH=.FALSE.
	DO I=1,NDATA
	  IF( .NOT.IS_DEBUG_HKL(IFILES(I),HKLMS(1,I)) ) CYCLE
C
	  LMATCH=.TRUE.
	  WRITE(10,'(A,3I4,A,I3,A,I6)') 'DEBUG:',(HKLMS(K,I),K=1,3),
	1					' (File',IFILES(I),') in global array',I
	  WRITE(10,'(A,F7.3,2F7.1,F8.2,F8.1,F7.1)') 'DEBUG: WAV,XY,TTH,COUNT,DCOUNT=',
	1				WAVS(I),XY_PIX(1,I),XY_PIX(2,I),TTHS(I),COUNTS(I),DCOUNTS(I)
	  IF(IDATA_MODE .GT. 1) WRITE(10,'(A,I6)') 'DEBUG: ITWIN=',HKLMS(4,I)
C
	ENDDO
C
	IF( .NOT.LMATCH ) WRITE(10,'(A)') 'DEBUG: No matching HKLs for OUTPUT_DEBUG_DATA'
C
	RETURN
	END


	SUBROUTINE OUTPUT_DEBUG_CALC
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	PARAMETER (NDEBUG_MAX=100)
	CHARACTER*80 NAMES_DEBUG(NDEBUG_MAX)
	INTEGER HKL_DEBUG(3,NDEBUG_MAX)
	COMMON /DEBUG_COM/ NAMES_DEBUG,HKL_DEBUG,NDEBUG
C
	LOGICAL IS_DEBUG_HKL,LMATCH
C
	IF(NDEBUG .EQ. 0) RETURN
C
	LMATCH=.FALSE.
	DO I=1,ILAST(NSEQ)
	  IF( .NOT.IS_DEBUG_HKL(IFILES(I),HKLMS(1,I)) ) CYCLE
C
	  LMATCH=.TRUE.
	  IF(IDATA_MODE .EQ. 1) THEN
	    WRITE(10,'(A,3I4,A,I3,A,3F10.1)') 'DEBUG:',(HKLMS(K,I),K=1,3),
	1					' (File',IFILES(I),') ICALC,IOBS,ISIG=',ICALC(I),IOBS(I),ISIG(I)
	  ELSE
	    WRITE(10,'(A,3I4,A,I3,A,I3,3F10.1)') 'DEBUG:',(HKLMS(K,I),K=1,3),
	1					' (File',IFILES(I),') ITWIN,ICALC,IOBS,ISIG=',
	2					HKLMS(4,I),ICALC(I),IOBS(I),ISIG(I)
	  ENDIF
C
	ENDDO
C
	IF( .NOT.LMATCH ) WRITE(10,'(A)') 'DEBUG: No matching HKLs for OUTPUT_DEBUG_CALC'
C
	RETURN
	END


	SUBROUTINE OUTPUT_DEBUG_SEQ
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	PARAMETER (NDEBUG_MAX=100)
	CHARACTER*80 NAMES_DEBUG(NDEBUG_MAX)
	INTEGER HKL_DEBUG(3,NDEBUG_MAX)
	COMMON /DEBUG_COM/ NAMES_DEBUG,HKL_DEBUG,NDEBUG
C
	LOGICAL IS_DEBUG_HKL,LMATCH
C
	IF(NDEBUG .EQ. 0) RETURN
C
	LMATCH=.FALSE.
	DO ISEQ=1,NSEQ
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    IF( .NOT.IS_DEBUG_HKL(IFILES(I),HKLMS(1,I)) ) CYCLE
C
	    LMATCH=.TRUE.
	    IF(IDATA_MODE .EQ. 1) THEN
	      WRITE(10,'(A,3I4,A,I3,A,I5,F10.1)') 'DEBUG:',(HKLMS(K,I),K=1,3),
	1					' (File',IFILES(I),') ISEQ,COUNTS_SEQ=',ISEQ,COUNTS_SEQ(ISEQ)
	    ELSE
	      WRITE(10,'(A,3I4,A,I3,A,I2,I5,F10.1)') 'DEBUG:',(HKLMS(K,I),K=1,3),
	1					' (File',IFILES(I),') ITWIN,ISEQ,COUNTS_SEQ=',
	2					HKLMS(4,I),ISEQ,COUNTS_SEQ(ISEQ)
	    ENDIF
C
	  ENDDO
	ENDDO
C
	IF( .NOT.LMATCH ) WRITE(10,'(A)') 'DEBUG: No matching HKLs for OUTPUT_DEBUG_SEQ'
C
	RETURN
	END
