C----- Routine to read/set parameters according to the run mode ---------------
c	SUBROUTINE READ_PARAM_FILES
C----- Routine to read parameters from the LaueG input file	-------------------
c	SUBROUTINE READ_LAUEG_INPUT_FILE(IDATA_MODE, IFILE_NUM,NFILES)
c	SUBROUTINE PRINT_TWINS_INPUT
c	SUBROUTINE PRINT_MODUL_INPUT
C----- Routine to read a *.ldm file and extract cell type	---------------------
c	SUBROUTINE READ_LDM_FILE(FILE_NAME)
C----- Get useful info. from the image file -----------------------------------
c	SUBROUTINE GET_IMAGE_INFO(IFILE)
C----- Routine to read wavelength calibration file ----------------------------
c	SUBROUTINE READ_WAV_FILE
C----- Utility routines to help with our file naming scheme	-------------------
c	SUBROUTINE GET_FILE_ARGS(IFILE1,IFILE2)
c	SUBROUTINE STORE_BASE_NAME(BASE_NAME0)
c	SUBROUTINE CHECK_FILE_NAMES(IFILE,EXT)
c	SUBROUTINE MAKE_FILE_NAME(FILE_NAME, IFILE,EXT)
C----- Routines to read the argonne_boxes data files --------------------------
c	SUBROUTINE READ_DATA_FILES(WAV_LO,CELL_MULT)
c	SUBROUTINE OPEN_INT_FILE(IOPT, FILE_NAME)
c	SUBROUTINE READ_INT_FILE(HKLS,XY_PIX,WAVS,COUNTS,DCOUNTS,
c	1						TTHS,ITWINS,NDATA,NMULT2,NMULT3,NFAIL,NREJ,
c	2                                            WAV_LO,FILE_NAME,IFILE)
c	SUBROUTINE READ_GEASC_FILE(HKLS,XY_PIX,WAVS,COUNTS,DCOUNTS,TTHS,
c	1					NDATA,NMULT2,NMULT3,NFAIL,NREJ, WAV_LO,FILE_NAME)
C------------------------------------------------------------------------------


C----- Routine to read/set parameters according to the run mode ---------------

	SUBROUTINE READ_PARAM_FILES
C
C Sets IFILE_MODE & IDATA_MODE in /MODES_COM/
C	IFILE_MODE=1 for LaueG mode, =2 for interactive Lauegen mode
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
	CHARACTER FILE_NAME*100,HOST_NAME*10,DATE*10
C
C In LaueG mode, all information is in the input file
C
	IF(IFILE_MODE .EQ. 1) THEN
		CALL READ_LAUEG_INPUT_FILE(IDATA_MODE, IFILE_NUM,NFILES)
C Check if files are numbered  as "1" or "001"
	  CALL CHECK_FILE_NAMES(IFILE_NUM(1),'.int')
C
C In the interactive Lauegen mode, information is in *.ldm & *.tif files
C
	ELSE
	  CALL GET_FILE_ARGS(IFILE1,IFILE2)
	  NFILES=IFILE2-IFILE1+1
	  IF(NFILES .GT. NFILES_MAX) THEN
	    WRITE(10,'(/,A,I4)') 'ERROR: Number of input files exceeds',NFILES_MAX
	    STOP 'ERROR: Too many input files'
	  ENDIF
C Check if files are numbered  as "1" or "001"
	  CALL CHECK_FILE_NAMES(IFILE1,'.ldm')
C Convert first & last files numbers to an array of file numbers
	  DO I=1,NFILES
	    IFILE_NUM(I)=I+IFILE1-1
	  ENDDO
C Read first *.ldm file to get cell/symmetry information
	  CALL MAKE_FILE_NAME(FILE_NAME, IFILE1,'.ldm')
	  CALL READ_LDM_FILE(TRIM(FILE_NAME))
C Get information from the header of the first image file
	  CALL GET_IMAGE_INFO(IFILE1)
C Laugen mode doesn't have twins or modulations
	  IDATA_MODE=1
	ENDIF
C
C Output the instrument used and the experiment date
C
	CALL GET_HOST_DATE(HOST_NAME,DATE)
	WRITE(10,'(/,4A)') 'Data from ',TRIM(HOST_NAME),' instrument, measured on ',TRIM(DATE)
C
	RETURN
	END


C----- Routine to read parameters from the LaueG input file -------------------

	SUBROUTINE READ_LAUEG_INPUT_FILE(IDATA_MODE, IFILE_NUM,NFILES)
C
	INTEGER IFILE_NUM(*)
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_CORR_COM/ FSCALE(NFILES_MAX)
	COMMON /WAV_CORR_COM/ IWAV_OPT,NWAV,WAV_GAMMA,WAV_MIN,WAV_STEP,WAV_PAR(20)
	COMMON /EFF_CORR_COM/ IEFF_OPT,EFF_POLY(10)
	COMMON /EXTI_CORR_COM/ IEXTI_OPT,EXTI
	COMMON /ABS_CORR_COM/ IABS_OPT,UR0,UR_LIN
C
	COMMON /CELL_TYPE_COM/ ITYPE,IFRIEDEL,ITYPE_ORIG,ICEN,CELL(6)
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	CHARACTER*10 HOST_NAME,DATE
	COMMON /FILE_INFO_COM/ NUMX,NUMY,HOST_NAME,DATE
C
	REAL MOD_VEC
	COMMON /MODULATE_COM/ NMOD_VEC,NMOD_IDX,MOD_VEC(3,100),MOD_IDX(100,100)
C
	CHARACTER STRING*100
C
	OPEN(UNIT=1,FILE='___laueg_laue4.in',STATUS='OLD',ERR=900)
C
	NLINE=1
	READ(1,*,ERR=910,END=920) IVERSION
	IF(IVERSION .NE. 6) THEN
	  WRITE(10,*) 'ERROR: LaueG input file not Version 6'
	  STOP 'ERROR: LaueG input file not Version 6'
	ENDIF
C
C The base name and the indices of the data files
C
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) NFILES
	NLINE=NLINE+1
	READ(1,'(A)',ERR=910,END=920) STRING
	CALL STORE_BASE_NAME(TRIM(STRING))
C 
C If NFILES=0, read the first and last file number
C Otherwise, read the NFILES file numbers
C Complain and die if too many files
C
	IF(NFILES .EQ. 0) THEN
	  NLINE=NLINE+1
	  READ(1,*,ERR=910,END=920) IFILE1,IFILE2
	  NFILES=IFILE2-IFILE1+1
	  IF(NFILES .GT. NFILES_MAX) THEN
	    WRITE(10,'(/,A,I4)') 'ERROR: Number of input files exceeds',NFILES_MAX
	    STOP 'ERROR: Too many input files'
	  ENDIF
	  DO I=1,NFILES
	    IFILE_NUM(I)=I+IFILE1-1
	  ENDDO
	ELSE
	  IF(NFILES .GT. NFILES_MAX) THEN
	    WRITE(10,'(/,A,I4)') 'ERROR: Number of input files exceeds',NFILES_MAX
	    STOP 'ERROR: Too many input files'
	  ENDIF
	  DO I=1,NFILES
	    NLINE=NLINE+1
	    READ(1,*,ERR=910,END=920) IFILE_NUM(I)
	  ENDDO
	ENDIF
C
C Read image size
C
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) NUMX,NUMY
C
C Read host-name, time-date and convert to IHOST & DATE
C
	NLINE=NLINE+1
	READ(1,'(A)',ERR=910,END=920) HOST_NAME						! host name
	NLINE=NLINE+1
	READ(1,'(A)',ERR=910,END=920) STRING							! date
	DATE='23/06/1808'	! Weird default date (based on Koalas!)
	I=LEN_TRIM(STRING)
	IF(I .GE. 8) DATE=STRING
C
C Flags to refine & turn on various corrections
C
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) LFILE_CORR
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) LWAV_CORR
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) LEFF_CORR
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) LEXTI_CORR
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) LHARM_CORR
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) LABS_CORR
C
C Read in the merge rules
C
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) ITYPE,IFRIEDEL,ITYPE_ORIG
C
C Read cell information for *.cif & *.inf files
C
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) ICEN,CELL
C
C Parameters for correction models
C
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) IWAV_OPT,NWAV,WAV_GAMMA,NBINS_NPAR
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) IEFF_OPT
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) IEXTI_OPT,EXTI_CORR_MAX
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) IABS_OPT,UR0,UR_LIN
C
C Misc. "internal parameters"
C
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) NSEQ_MIN,SEQ_FRAC
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) X_LO,X_HI,Y_LO,Y_HI
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) WAV_LO,WAV_HI
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) CELL_MULT
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) ALL_SIG_MULT,ALL_REL_ERR,WEAK_SIG_MULT
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) EFF_WIDTH,EFF_THICK
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) WAV_ERR1,WAV_ERR2,EXTI_ERR
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) IVERBOSE,HARM_CUTOFF_MULT
C
C Read flags for modulated and twinned lattices
C
	NLINE=NLINE+1
	READ(1,*,ERR=910,END=920) ITWINS,IMODUL
	IDATA_MODE=1
	IF(ITWINS .NE. 0) IDATA_MODE=2
	IF(IMODUL .NE. 0) IDATA_MODE=3
C
C Finish up with modulation vectors and indices
C
	IF(IDATA_MODE .EQ. 3) THEN
C Read vectors and indices
		NLINE=NLINE+1
		READ(1,*,ERR=910,END=920) NMOD_VEC,NMOD_IDX
		DO I=1,NMOD_VEC
		  NLINE=NLINE+1
		  READ(1,*,ERR=910,END=920) (MOD_VEC(K,I),K=1,3)
		ENDDO
		DO I=1,NMOD_VEC
		  NLINE=NLINE+1
		  READ(1,*,ERR=910,END=920) (MOD_IDX(K,I),K=1,NMOD_IDX)
		ENDDO
C Output vectors and indices to log file
		CALL PRINT_MODUL_INPUT
	ENDIF
C
C Close input file
C
	CLOSE(UNIT=1)
C
	RETURN
C
C Error messages
C
900	WRITE(10,*) 'ERROR: Unable to open LaueG input file'
	STOP 'ERROR: Unable to open input file'
910	WRITE(10,'(A,I3,A)') 'ERROR: Reading line',
	1		NLINE,' of the LaueG input file'
	STOP 'ERROR: Error reading the LaueG input file'
920	WRITE(10,'(A,I3)') 'ERROR: LaueG input file truncated at line',NLINE
	STOP 'ERROR: LaueG input file is truncated'
C
	END


	SUBROUTINE PRINT_MODUL_INPUT
C
	REAL MOD_VEC
	COMMON /MODULATE_COM/ NMOD_VEC,NMOD_IDX,MOD_VEC(3,100),MOD_IDX(100,100)
C
	REAL QVEC(3)
C
	WRITE(10,'(/,A)') 'Modulation Information:'
C
	WRITE(10,'(9X,A)') 'Wavevectors'
	DO I=1,NMOD_VEC
	  WRITE(10,'(2X,A,I1,1X,3F8.3)') 'q',I,(MOD_VEC(K,I),K=1,3)
	ENDDO
C
	WRITE(10,'(A)') '  Modulation Numbers (m) as Vector Indices and HKL Displacement'
	WRITE(10,'(3X,A,<NMOD_VEC>(A,I1),A)') 'm ',('   q',K,K=1,NMOD_VEC),'       dH     dK     dL'
C Also output the main reflection as m=0
	WRITE(10,'(I4,1X,<NMOD_VEC>I5,4X,3F7.3)') 0,(0,K=1,NMOD_VEC),0.0,0.0,0.0
	DO I=1,NMOD_IDX
	  CALL CALC_MOD_QVEC(QVEC,I)
	  WRITE(10,'(I4,1X,<NMOD_VEC>I5,4X,3F7.3)') I,(MOD_IDX(I,K),K=1,NMOD_VEC),QVEC
	ENDDO
C
	RETURN
	END


C ========= Routine to read a *.ldm file and extract cell type  =========

	SUBROUTINE READ_LDM_FILE(FILE_NAME)
C
C Read the *.ldm file FILE_NAME to get the cell/symmetry data
C
	CHARACTER FILE_NAME*(*)
C
	CHARACTER TAG*20,CTYPE*4
C
	COMMON /CELL_TYPE_COM/ ITYPE,IFRIEDEL,ITYPE_ORIG,ICEN,CELL(6)
C
	DATA ICUB,ITETA,ITETB,ITETC,IORTH,IHEX,ITRIG,IRHOM,IMONA,
	1	IMONB,IMONC,IDENT/1,2,3,4,5,6,7,8,9,10,11,12/
C
C Set CTYPE & CELL so we can check if they have been updated
C
	CTYPE=REPEAT(' ',LEN(CTYPE))
	CELL(1)=0.0
C
C Open the *.ldm file FILE_NAME and read cell constants and type
C
	OPEN(UNIT=1,FILE=FILE_NAME,STATUS='OLD',ERR=900)
100	TAG=REPEAT(' ',LEN(TAG))
	READ(1,*,END=110) TAG
	CALL UP2LOW_CASE(TAG,TAG)
	IF(TAG .EQ. 'system') THEN
	  BACKSPACE(1)
	  READ(1,*,ERR=901) TAG,CTYPE
	  CALL UP2LOW_CASE(CTYPE, CTYPE)
	ELSE IF(TAG .EQ. 'a ') THEN
	  BACKSPACE(1)
	  READ(1,*,ERR=902) (TAG,CELL(K),K=1,6)
	  CALL UP2LOW_CASE(TAG,TAG)
	  IF(TAG .NE. 'gamma') GOTO 902
	ENDIF
	GOTO 100
110	CLOSE(UNIT=1)
C
C If CTYPE or CELL are the initial values print out an error message
C
	IF(CTYPE .EQ. ' ') GOTO 901
	IF(CELL(1) .EQ. 0.0) GOTO 902
C
	ITYPE=0
	IF     (CTYPE .EQ. 'cubi') THEN
	  ITYPE=ICUB
	ELSE IF(CTYPE .EQ. 'tetr') THEN
	  IF(ABS(CELL(2)-CELL(3)) .LT. 1E-5) ITYPE=ITETA
	  IF(ABS(CELL(1)-CELL(3)) .LT. 1E-5) ITYPE=ITETB
	  IF(ABS(CELL(1)-CELL(2)) .LT. 1E-5) ITYPE=ITETC
	  IF(ITYPE .EQ. 0) GOTO 910
	ELSE IF(CTYPE .EQ. 'orth') THEN
	  ITYPE=IORTH
	ELSE IF(CTYPE .EQ. 'hexa') THEN
	  ITYPE=IHEX					
	ELSE IF(CTYPE .EQ. 'rhom') THEN
	  ITYPE=IRHOM
	ELSE IF(CTYPE .EQ. 'mono') THEN
	  IF(ABS(CELL(4)-90.0) .GT. 1E-4) ITYPE=IMONA
	  IF(ABS(CELL(5)-90.0) .GT. 1E-4) ITYPE=IMONB
	  IF(ABS(CELL(6)-90.0) .GT. 1E-4) ITYPE=IMONC
	  IF(ITYPE .EQ. 0) GOTO 910
	ELSE IF(CTYPE .EQ. 'tric') THEN
	  ITYPE=IDENT
	ELSE
	  GOTO 911
	ENDIF
C
C Save the original cell type
C
	ITYPE_ORIG=ITYPE
C
C Lauegen doesn't distinguish trigonal from hexagonal, so safest
C to assume it is trigonal. Can be reset in OPTIONS menu.
C
	IF(ITYPE .EQ. IHEX) ITYPE=ITRIG					
C
C Default is merge Friedels (set even if implied by cell type)
C
	IFRIEDEL=1
C
	RETURN
C
C Print out error messages
C
900	STOP 'BUG(read_ldm_file): Cannot open LDM file'
C
901	PRINT *,'SYSTEM line is invalid or missing'
	STOP 'ERROR: First *.ldm file is invalid'
C
902	PRINT *,'A,B,C,ALPHA,BETA,GAMMA line is invalid or missing'
	STOP 'ERROR: First *.ldm file is invalid'
C
910	PRINT '(1X,A,3F8.4,3F8.3,2A)','Inconsistent cell=',CELL,' for ',CTYPE
	STOP 'ERROR: Invalid first *.ldm file'
C
911	PRINT '(1X,2A)','Unknown SYSTEM type =',CTYPE
	STOP 'ERROR: Invalid first *.ldm file'
	END


C ===== Get useful info. from the image file =====

	SUBROUTINE GET_IMAGE_INFO(IFILE)
C
C In lauegen mode we try to get needed info from first TIF header
C
	CHARACTER*10 HOST_NAME,DATE
	COMMON /FILE_INFO_COM/ NUMX,NUMY,HOST_NAME,DATE
C
	CHARACTER*2000 COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,ADC_INFO
	COMMON /TIFFINFO_COM/ COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,
	1			IOFF_DATA,NUMX2,NUMY2,IXYRES,ISTARTX,ISTARTY,ISPEED,
	2			EXPOSE_TIME,EXPOSE_PHI,TEMP_BEGIN,TEMP_END,TEMP_MIN,TEMP_MAX,
	3			IADC,ADC_ADD,ADC_DIV,ADC_INFO,IAPER_1,IAPER_2,IKAPPA
C
	CHARACTER TIFF_NAME*132,STRING*80
C
C Try reading the host name and date/time from the TIF file with index IFILE
C
	CALL MAKE_FILE_NAME(TIFF_NAME, IFILE,'.tif')
	CALL READ_TIFF_INFO(ISTATUS, TRIM(TIFF_NAME),1)
C
C Handle the cases where the tif header could not be read
C
C If TIF header OK, overwrite the default values
C
	IF(ISTATUS .EQ. 0) THEN
	  CALL GET_TIFF_INFO(STRING,VALUE,NUMX, 'NUMX')
	  CALL GET_TIFF_INFO(STRING,VALUE,NUMY, 'NUMY')
C Convert the date & time string to date only
	  CALL GET_TIFF_INFO(STRING,VALUE,IVALUE, 'DATE')
	  I=LEN_TRIM(STRING)
	  IF(I .GT. 16) DATE=STRING(I-9:I)
C Copy the host name
	  CALL GET_TIFF_INFO(HOST_NAME,VALUE,IVALUE, 'HOST')
	ELSE
C Can't read header, so set some silly default values and complain
		NUMX=4000
		NUMY=2000
		HOST_NAME='unknown'
		DATE='23/06/1808'	! Weird default date (based on Koalas!)
C
		IF(ISTATUS .EQ. -1) THEN
			WRITE(10,'(/,3A)') 'Unable to find ',TRIM(TIFF_NAME),
	1			  ', cannot determine all experiment information'
		ELSE
			WRITE(10,'(/,3A)') TRIM(TIFF_NAME),' appears invalid',
	1				  ', cannot determine all experiment information'
	  ENDIF
C
	  WRITE(10,'(2A)') 'WARNING: Assuming a Laue IP instrument ',
	1			'with 4000 x 2000 pixels'
	ENDIF
C
	RETURN
	END


C ============ Routine to read wavelength calibration file  ===========

	SUBROUTINE READ_WAV_FILE
C
	COMMON /WAV_FILE_COM/ NWAV_FILE,WAV_FILE_MIN,WAV_FILE_STEP,WAV_FILE(10000)
C
	CHARACTER STRING*80
C
C Open the wavelength calibration file
C
	WRITE(10,'(/,A)') 'Reading wavelength calibration file: laue4_wav.dat'
	OPEN(UNIT=1,STATUS='OLD',FILE='laue4_wav.dat',ERR=900)
C
C Read header line from file, complain and die if it is incorrect
C
	STRING=REPEAT(' ',LEN(STRING))
	READ(1,'(Q,A)') ILEN,STRING(1:ILEN)
	IF(STRING(1:39) .NE. '===LAUE4 WAVELENGTH CALIBRATION FILE===') THEN
	  PRINT *,'Invalid header in wavelength calibration file laue4_wav.dat'
	  STOP 'ERROR: Invalid (or old version) wavelength calibration file'
	ENDIF
C
C Copy next two lines to the log file
C
	STRING=REPEAT(' ',LEN(STRING))
	READ(1,'(Q,A)') ILEN,STRING(1:ILEN)
	WRITE(10,'(2X,A)') TRIM(STRING)
	STRING=REPEAT(' ',LEN(STRING))
	READ(1,'(Q,A)') ILEN,STRING(1:ILEN)
	WRITE(10,'(2X,A)') TRIM(STRING)
C
C Read the first two lines of numeric data to get WAV_FILE_STEP
C
	READ(1,*) WAV_FILE_MIN,WAV_FILE(1)
	READ(1,*) WAV,WAV_FILE(2)
	WAV_FILE_STEP=WAV-WAV_FILE_MIN
C
C Read remaining numeric data and complain and die if it
C doesn't agree with WAV_FILE_STEP
C
	DO I=3,10000
	  READ(1,*,END=100) WAV,WAV_FILE(I)
	  IF(ABS((WAV-WAV_FILE_MIN)/(I-1) - WAV_FILE_STEP) .GT. 0.001)
	1	STOP 'ERROR: Wavelength calibration file not regular spaced'
	ENDDO
	STOP 'ERROR: Wavelength file > 10,000 points'
100	CLOSE(UNIT=1)
C
C Recalculate WAV_FILE_STEP now we have all the data
C
	WAV_MAX=WAV
	NWAV_FILE=I-1
	WAV_FILE_STEP=(WAV_MAX-WAV_FILE_MIN)/(NWAV_FILE-1)
C
C Output the wavelength range of the calibration file, and return
C
	WRITE(10,'(2X,3(A,F6.3))') 'Calibration data from',WAV_FILE_MIN,' to',
	1				WAV_MAX,' in steps of',WAV_FILE_STEP
	RETURN
C
900	PRINT *,'Cannot find the wavelength calibration file laue4_wav.dat'
	STOP 'ERROR: Unable to find calibration file'
	END



C ============ Utility routines to help with our file naming scheme  ===========

	SUBROUTINE GET_FILE_ARGS(IFILE1,IFILE2)
C
	INTEGER*2  I2ARG,I2LEN
	CHARACTER*80 NAME1,NAME2
C
C Set boolean if we have 2 command line arguments
C
	IF(NARGS() .LT. 2) THEN
	  PRINT *,'First and last file names must follow command'
	  GOTO 1000
	ENDIF
C
C Load NAME1 & 2 from the command line arguments
C
	NAME1=REPEAT(' ',LEN(NAME1))
	I2ARG=1
	CALL GETARG(I2ARG, NAME1, I2LEN)
C
	NAME2=REPEAT(' ',LEN(NAME2))
	I2ARG=2
	CALL GETARG(I2ARG, NAME2, I2LEN)
C
C Get length of names (trimmed of blanks)
C
	LEN1=LEN_TRIM(NAME1)
	LEN2=LEN_TRIM(NAME2)
C
C If NAME* doesn't end in a numeral, try removing a .*** extension
C
	IF(NAME1(LEN1:LEN1).LT.'0' .OR. NAME1(LEN1:LEN1).GT.'9') THEN
	  IF(LEN1 .GT. 4) THEN
	    IF(NAME1(LEN1-3:LEN1-3) .EQ. '.') LEN1=LEN1-4
	  ENDIF
	ENDIF
C
	IF(NAME2(LEN2:LEN2).LT.'0' .OR. NAME2(LEN2:LEN2).GT.'9') THEN
	  IF(LEN2 .GT. 4) THEN
	    IF(NAME2(LEN2-3:LEN2-3) .EQ. '.') LEN2=LEN2-4
	  ENDIF
	ENDIF
C
C Find the position of the last non-numeric in NAME*
C
	INON1=0
	DO I=1,LEN1
	  IF(NAME1(I:I).LT.'0' .OR. NAME1(I:I).GT.'9') INON1=I
	ENDDO
C
	INON2=0
	DO I=1,LEN2
	  IF(NAME2(I:I).LT.'0' .OR. NAME2(I:I).GT.'9') INON2=I
	ENDDO
C
C Do various checks, on errors complain and die 
C
	IF(INON1.EQ.LEN1 .OR. INON2.EQ.LEN2) THEN
	  PRINT *,'Numeric part of one or both filenames is missing'
	  GOTO 1000
	ENDIF
C
	IF(INON1 .NE. INON2) THEN
	  PRINT *,'The base part of the file names do not match' 
	  GOTO 1000
	ENDIF
	IF(INON1 .GT. 0) THEN
	  IF(NAME1(1:INON1) .NE. NAME1(1:INON2)) THEN
	    PRINT *,'The base part of the file names do not match' 
	    GOTO 1000
	  ENDIF
	ENDIF
C
C Decode the numeric parts
C
	READ(NAME1(INON1+1:LEN1),'(I20)') IFILE1
	READ(NAME2(INON2+1:LEN2),'(I20)') IFILE2
	IF(IFILE1 .EQ. IFILE2) THEN
	  PRINT *,'The two file names cannot be identical'
	  GOTO 1000
	ELSE IF(IFILE1 .GT. IFILE2) THEN
	  PRINT *,'Numeric part of second file smaller than first'
	  GOTO 1000
	ENDIF
C
C Store the base name for later use
C
	CALL STORE_BASE_NAME(NAME1(1:INON1))
C
C SUCCESS, so return
C
	RETURN
C
C Print out syntax info, then die
C
1000	PRINT *
	PRINT *,'Please try again using a command similar to:'
	PRINT *,'  laue4 Si_7 Si_24'
	PRINT *
	STOP 'ERROR: Invalid command syntax'
C
	END



	SUBROUTINE STORE_BASE_NAME(BASE_NAME0)
C
	CHARACTER BASE_NAME0*(*)
C
	LOGICAL LSTRIP
	CHARACTER BASE_NAME*255
	COMMON /FILE_NAME_COM/ BASE_NAME,LSTRIP
C
C Store BASE_NAME0 into the COMMON variable BASE_NAME
C
      ILEN=LEN_TRIM(BASE_NAME0)
	IF(ILEN .GT. 254) STOP 'ERROR: Base name of file too long'
	BASE_NAME=BASE_NAME0(1:ILEN)
C
C Make sure last character is a "_"
C
      IF(BASE_NAME(ILEN:ILEN) .NE. '_') BASE_NAME=TRIM(BASE_NAME)//'_'
C
	RETURN
	END


	SUBROUTINE CHECK_FILE_NAMES(IFILE,EXT)
C
C Checks for the file "base_name"//<numeric>//EXT, where <numeric> is
C the number IFILE expressed as "001" or "1". If neither or both are
C found we complain and die. Otherwise we set LSTRIP to TRUE for "1"
C and use this in MAKE_FILENAME.
C
	CHARACTER EXT*(*)
C
	LOGICAL LSTRIP
	CHARACTER BASE_NAME*255
	COMMON /FILE_NAME_COM/ BASE_NAME,LSTRIP
C
	LOGICAL LFOUND1,LFOUND2
	CHARACTER SNUM*10,FILE_NAME1*100,FILE_NAME2*100
C
C Make FILE_NAME1 using a numeric part stripped of zeroes (eg. "1")
C
	FILE_NAME1=REPEAT(' ',100)
	WRITE(SNUM,'(I0)',IOSTAT=IDUMMY) IFILE
	FILE_NAME1=TRIM(BASE_NAME)//TRIM(SNUM)//TRIM(EXT)
C
C Make FILE_NAME2 using a 3 digit numeric part (eg. "001")
C
	FILE_NAME2=REPEAT(' ',100)
	WRITE(SNUM,'(I3.3)',IOSTAT=IDUMMY) IFILE
	FILE_NAME2=TRIM(BASE_NAME)//SNUM(1:3)//TRIM(EXT)
C
C Try to find both types of files
C
	INQUIRE(FILE=FILE_NAME1,EXIST=LFOUND1)
	INQUIRE(FILE=FILE_NAME2,EXIST=LFOUND2)
C
C Do sanity check on results, complain and die on errors
C
	IF( .NOT.(LFOUND1 .OR. LFOUND2) ) THEN
	  WRITE(10,'(4A)') 'ERROR: Unable to find ',TRIM(FILE_NAME1),
	1					' or ',TRIM(FILE_NAME2)
	  STOP 'ERROR: Cannot find input files'
	ELSEIF( LFOUND1 .AND. LFOUND2 .AND. (FILE_NAME1 .NE. FILE_NAME2) ) THEN
	  WRITE(10,'(4A)') 'ERROR: Found both ',TRIM(FILE_NAME1),
	1					' and ',TRIM(FILE_NAME2)
	  STOP 'ERROR: Ambiguous input file names'
	ENDIF
C
	LSTRIP=LFOUND1
	RETURN
	END


	SUBROUTINE MAKE_FILE_NAME(FILE_NAME, IFILE,EXT)
C
C Returns FILE_NAME made from: stored BASE_NAME, integer IFILE as a string,
C the string EXT, and right padded with spaces.
C If IFILE<0, ignore the numeric part
C The IFILE string is zero padded if LSTRIP in /FILE_NAME_COM/ is FALSE
C
	CHARACTER FILE_NAME*(*),EXT*(*)
C
	LOGICAL LSTRIP
	CHARACTER BASE_NAME*255
	COMMON /FILE_NAME_COM/ BASE_NAME,LSTRIP
C
C Do case of no numeric part (strip trailing '_' in BASE_NAME)
C
      IF(IFILE .LE. 0) THEN
        ILEN=LEN_TRIM(BASE_NAME)
        IF(BASE_NAME(ILEN:ILEN) .EQ. '_') ILEN=ILEN-1
        FILE_NAME=BASE_NAME(1:ILEN)//TRIM(EXT)
        RETURN
	ENDIF
C
C Append IFILE as a string (with/without leading zeroes) to BASE_NAME
C
	IF( LSTRIP ) THEN
	  WRITE(FILE_NAME,'(A,I0)') TRIM(BASE_NAME),IFILE
	ELSE
	  WRITE(FILE_NAME,'(A,I3.3)') TRIM(BASE_NAME),IFILE
	ENDIF
C
C Add EXT string
C
      IF(EXT .NE. '') FILE_NAME=TRIM(FILE_NAME)//TRIM(EXT)
C
	RETURN
	END


C ========= Routines to read the argonne_boxes data files ========

	SUBROUTINE READ_DATA_FILES(WAV_LO,CELL_MULT)
C
C IFILE_MODE=1, LaueG mode with *.int files
C            2, Lauegen mode with *.geasc files
C
C For Lauegen mode this routine sets IDATA_MODE=1 (off), for
C LaueG mode the routine READ_INT_FILE returns IDATA_MODE.
C Assumes GET_DATA_MODE has already been called
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	LOGICAL IS_DEBUG_HKL
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
	CHARACTER FILE_NAME*132,BASE_NAME*132
C
	PARAMETER (NREAD_MAX=20000)
	INTEGER HKLS0(3,NREAD_MAX),ITWINS0(NREAD_MAX)
	REAL XY_PIX0(2,NREAD_MAX),TTHS0(NREAD_MAX),WAVS0(NREAD_MAX)
	REAL COUNTS0(NREAD_MAX),DCOUNTS0(NREAD_MAX)
C
C Load debug HKLs and base-names
C
	CALL LOAD_DEBUG_HKLS
C
C Tally the number of different invalid spots across all files
C
	NNREJ1=0
	NNREJ2=0
	NNMULT2=0
	NNMULT3=0
	NNFAIL=0
C
C Loop through the requested file numbers
C   NDATA is the total number of valid spots entered into arrays
C   NFILES is the number of files we expect to read
C   NFILES2 is the number of valid files that have been read
C   IFILE is the number of files we have attempted to read
C   INUM is the file number appended to the base name
C   IFILE_NUM() is a list of file numbers where INUM=IFILE_NUM(IFILE)
C   HKLMS(4,*) is the twin or modulation number for IDATA_MODE=2,3
C
	NDATA=0
	NFILES2=0
	DO IFILE=1,NFILES
C Create the base-name and the mode specific file-name
	  INUM=IFILE_NUM(IFILE)
	  CALL MAKE_FILE_NAME(BASE_NAME, INUM,'')
	  FILE_NAME=TRIM(BASE_NAME)//'.int'
	  IF(IFILE_MODE .NE. 1) FILE_NAME=TRIM(BASE_NAME)//'.geasc'
C Log which data file we are trying to read
	  WRITE(10,'(A,I3,2A)') 'Reading data file',NFILES2+1,' : ',TRIM(FILE_NAME)
C Load reject HKLs for this particular base-name
	  CALL READ_REJECT_HKLS(NREJ1, TRIM(BASE_NAME),'laue4_rej.dat')
C Attempt to open and read the '.geasc' or '.int' input file
	  IF(IFILE_MODE .EQ. 1) THEN
C Open file, check header, and keep unit=1 open
          CALL OPEN_INT_FILE(IOPT, TRIM(FILE_NAME))
C Complain and die if this file has a different twin mode
	    IF(IOPT .NE. IDATA_MODE) THEN
	      WRITE(10,'(/,A)') 'ERROR: *.int files have mixed normal/twin/modulate options'
	      STOP 'ERROR: *.int files have mixed normal/twin/modulate options'
	    ENDIF
C Read the rest of the file and close unit=1
          CALL READ_INT_FILE(HKLS0,XY_PIX0,WAVS0,COUNTS0,DCOUNTS0,
	1						TTHS0,ITWINS0,NHKLS,NMULT2,
	2						NMULT3,NFAIL,NREJ2, WAV_LO,TRIM(FILE_NAME),IFILE)
	  ELSE
	    CALL READ_GEASC_FILE(HKLS0,XY_PIX0,WAVS0,COUNTS0,DCOUNTS0,TTHS0,
	1				NHKLS,NMULT2,NMULT3,NFAIL,NREJ2, WAV_LO,TRIM(FILE_NAME))
	  ENDIF
C Skip the file if it is missing (NHKLS<0)
	  IF(NHKLS .LT. 0) THEN
	    PRINT '(1X,2A)','WARNING: Skipped the missing file ',TRIM(FILE_NAME)
	    WRITE(10,'(1X,2A)') 'WARNING: Skipped the missing file ',TRIM(FILE_NAME)
	    CYCLE
	  ENDIF
C Increment NFILES2 and update IFILE_NUM() in case we skipped a file,
C then log the file name and number
	  NFILES2=NFILES2+1
	  IFILE_NUM(NFILES2)=INUM
C Update the tally of different invalid spots for all files
	  NNMULT2=NNMULT2+NMULT2
	  NNMULT3=NNMULT3+NMULT3
	  NNFAIL=NNFAIL+NFAIL
	  NNREJ1=NNREJ1+NREJ1
	  NNREJ2=NNREJ2+NREJ2
C Loop through all data in file, copying from local to COMMON storage
	  DO I=1,NHKLS
	    NDATA=NDATA+1
	    IF(NDATA .GT. NDATA_MAX) THEN
	      WRITE(10,'(/,A,I8)') 'ERROR: Total number of valid spots exceeds',NDATA_MAX
	      STOP 'ERROR: Too many spots in total'
	    ENDIF
	    HKLMS(1,NDATA)=HKLS0(1,I)
	    HKLMS(2,NDATA)=HKLS0(2,I)
	    HKLMS(3,NDATA)=HKLS0(3,I)
	    HKLMS(4,NDATA)=ITWINS0(I)
	    IFILES(NDATA)=NFILES2
	    WAVS(NDATA)=WAVS0(I)*CELL_MULT
	    COUNTS(NDATA)=COUNTS0(I)
C The +1.0 is to compensate for integer rounding of COUNTS() in data files
	    DCOUNTS(NDATA)=MAX(0.0,DCOUNTS0(I))+1.0
	    TTHS(NDATA)=TTHS0(I)
	    XY_PIX(1,NDATA)=XY_PIX0(1,I)
	    XY_PIX(2,NDATA)=XY_PIX0(2,I)
C If a debug hkl, output where it is copied to
	    IF( .NOT.IS_DEBUG_HKL(IFILE,HKLS0(1,I)) ) CYCLE
	    IF(IDATA_MODE .EQ. 1) THEN
	      WRITE(10,'(A,3I4,A,I3,A,I6)') 'DEBUG:',(HKLMS(K,NDATA),K=1,3),
	1					' (File',IFILE,') moved to global array',NDATA
	    ELSE
	      WRITE(10,'(A,4I4,A,I3,A,I6)') 'DEBUG:',(HKLMS(K,NDATA),K=1,4),
	1					' (File',IFILE,') moved to global array',NDATA
          ENDIF
C End of loop over data in file
        ENDDO
C End of loop over data files
	ENDDO
C Set NFILES to NFILES2 in case we skipped some files
	NFILES=NFILES2
C
C Output to console and log file the total number of files read
C
	NTOTAL=NDATA+NNMULT3+NNMULT2+NNFAIL+NNREJ2
	PRINT '(I7,A)',NTOTAL,' reflections read from data files'
	WRITE(10,'(I7,A)') NTOTAL,' reflections read from data files'
C
	IF(ABS(CELL_MULT-1.0) .GT. 1E-6) THEN
	  WRITE(10,*) 'Data file wavelengths have been scaled to simulate'
	  WRITE(10,'(3X,A,F7.3)') 'the cell lengths multiplied by',CELL_MULT
	ENDIF
C
C Output rejected spots tally to log file (and 1 warning to console)
C
	IF(NNREJ1 .LT. 0) THEN
	  WRITE(10,'(/,A)') 'Rejection list file "laue4_rej.dat" not found'
	ELSE IF(NNREJ1 .EQ. NNREJ2) THEN
	  WRITE(10,'(/,I7,A)') NNREJ2,' rejected as entries in "laue4_rej.dat"'
	ELSE
	  WRITE(10,'(/,I7,A,I5,A)') NNREJ2,' rejected (WARNING: ',NNREJ1,
	1								' entries in "laue4_rej.dat")'
	  PRINT '(I7,A,I5,A)',NNREJ2,' rejected (WARNING: ',NNREJ1,
	1								' entries in "laue4_rej.dat")'
	ENDIF
	WRITE(10,'(I7,A)') NNFAIL,' rejected due to integration failure'
	WRITE(10,'(I7,A)') NNMULT3,' rejected due to multiplicity > 2'
	WRITE(10,'(I7,A)') NNMULT2,' rejected due to multiplicity = 2'
C
      RETURN
      END


      SUBROUTINE OPEN_INT_FILE(IOPT, FILE_NAME)
C
C Open a LaueG intensity file, read/check the header, and
C returns the Twin mode, and keep the file open on Unit=1
C
C IOPT = 1 (normal mode), 2 (twin mode), 3 (modulate mode)
C
	CHARACTER FILE_NAME*(*)
C
	CHARACTER LINE*100
C
	OPEN(UNIT=1,FILE=FILE_NAME,STATUS='OLD',ERR=900)
C
	READ(1,'(Q,A,//)') ILEN,LINE
C
C Check the fixed text for the header line
C
	IF(LINE(1:30) .NE. 'LaueG Intensities File Version') THEN
	  WRITE(10,'(/,2A)') 'ERROR: Invalid header in ',FILE_NAME
	  STOP 'ERROR: Invalid data file header'
	ELSE IF(LINE(33:40) .NE. ', Option') THEN
	  WRITE(10,'(/,2A)') 'ERROR: Invalid header in ',FILE_NAME
	  STOP 'ERROR: Invalid data file header'
	ENDIF
C
C Complain and die if the version and option numbers are invalid
C
	READ(LINE(31:42),'(I2,8X,I2)') IVER,IOPT
	IF(IVER .NE. 1) THEN
	  WRITE(10,'(/,A,I2,2A)') 'ERROR: Invalid header version=',IVER,' in ',FILE_NAME
	  STOP 'ERROR: Invalid data file header'
	ENDIF
C
	IF(IOPT.LT.1 .OR. IOPT.GT.3) THEN
	  WRITE(10,'(/,A,I2,2A)') 'ERROR: Invalid header option=',IOPT,' in ',FILE_NAME
	  STOP 'ERROR: Invalid data file header'
	ENDIF
C
	RETURN
C
C Complain and die if something goes wrong
C
900	WRITE(10,'(/,2A)') 'ERROR: Cannot open file ',FILE_NAME
	STOP 'ERROR: Missing data file'
C
901	WRITE(10,'(/,2A)') 'ERROR: Cannot find "REFLECTION DATA" line in ',FILE_NAME
	STOP 'ERROR: Invalid data file'
	RETURN
      END

      
	SUBROUTINE READ_INT_FILE(HKLS,XY_PIX,WAVS,COUNTS,DCOUNTS,
	1						TTHS,ITWINS,NDATA,NMULT2,NMULT3,NFAIL,NREJ,
     2                                            WAV_LO,FILE_NAME,IFILE)
C
	CHARACTER FILE_NAME*(*)
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	PARAMETER (NREAD_MAX=20000)
	INTEGER HKLS(3,NREAD_MAX),ITWINS(NREAD_MAX)
	REAL XY_PIX(2,NREAD_MAX),TTHS(NREAD_MAX),WAVS(NREAD_MAX)
	REAL COUNTS(NREAD_MAX),DCOUNTS(NREAD_MAX)
C
	LOGICAL IS_REJECT_HKL,IS_DEBUG_HKL,LDEBUG
	INTEGER IHKL(3)
C
C MAIN LOOP: Read in data and keep a tally of types of spots
C
	NDATA=0
	NREJ=0
	NFAIL=0
	NMULT2=0
	NMULT3=0
	ITWIN=0
	DO IREAD=1,NREAD_MAX+1
C
C IHKL() is the Miller indices
C ITWIN is the optional twin or modulation number
C X & Y are the pixel coords of the spot
C WAV is the wavelength of the lowest multiple of the spot
C IMULT is the wavelength multiplicity of the spot
C TTH is the two-theta
C PERP is the angle of the beam to the IP perpendicular
C ICNTS & IESDS the integrated intensity and its e.s.d.
C
	  IF(IDATA_MODE .EQ. 1) THEN
	    READ(1,'(3I4,F7.3,I3,F8.1,F7.1,F8.2,I7,I6)',END=1000) IHKL,WAV,IMULT,X,Y,TTH,ICNTS,IESDS
	  ELSE
	    READ(1,'(4I4,F7.3,I3,F8.1,F7.1,F8.2,I7,I6)',END=1000) IHKL,ITWIN,WAV,IMULT,X,Y,TTH,ICNTS,IESDS
	  ENDIF
C
C Check if a debug hkl, if so output the data read in
	  LDEBUG=IS_DEBUG_HKL(IFILE,IHKL)
	  IF( LDEBUG ) THEN
	    WRITE(10,'(A,3I4,A,I3,A)') 'DEBUG:',IHKL,' (File',IFILE,') read from '//TRIM(FILE_NAME)
	    IF(IDATA_MODE .EQ. 1) THEN
	      WRITE(10,'(A,F7.3,I3,F8.1,F7.1,F8.2,I7,I6)') 'DEBUG: WAV,IMULT,X,Y,TTH,COUNTS+ESD=',
	1				WAV,IMULT,X,Y,TTH,ICNTS,IESDS
	    ELSE
	      WRITE(10,'(A,F7.3,I3,F8.1,F7.1,F8.2,I7,I6)') 'DEBUG: ITWIN,WAV,IMULT,X,Y,TTH,COUNTS+ESD=',
	1				ITWIN,WAV,IMULT,X,Y,TTH,ICNTS,IESDS
	    ENDIF
	  ENDIF
C
C If spot is valid add to the data list, else update rejection statistics
	  IF( IS_REJECT_HKL(IHKL,ITWIN) ) THEN		! Reject spot if on the rejection list
	    NREJ=NREJ+1
	    IF( LDEBUG ) WRITE(10,'(A)') 'DEBUG: Rejected'
	  ELSE IF(IESDS .LT. 0) THEN		! Reject spot if integration failure
	    NFAIL=NFAIL+1
	    IF( LDEBUG ) WRITE(10,'(A)') 'DEBUG: Integration failure'
	  ELSE IF(IMULT .GE. 3) THEN			! Reject spot if 3 or more multiples
	    NMULT3=NMULT3+1
	    IF( LDEBUG ) WRITE(10,'(A)') 'DEBUG: 3 or more multiples'
	  ELSE IF(IMULT.EQ.2 .AND.			! Reject if 2 multiples and
	1				WAV/2.0.LT.WAV_LO) THEN		! wav/2 is not within limits
	    NMULT2=NMULT2+1
	    IF( LDEBUG ) WRITE(10,'(A)') 'DEBUG: 2 multiples and not within wav/2 limits'
	  ELSE								! Add to the data lists
	    NDATA=NDATA+1
	    IF(NDATA .GT. NREAD_MAX) GOTO 900
	    HKLS(1,NDATA)=IHKL(1)
	    HKLS(2,NDATA)=IHKL(2)
	    HKLS(3,NDATA)=IHKL(3)
	    ITWINS(NDATA)=ITWIN
	    WAVS(NDATA)=WAV
	    COUNTS(NDATA)=ICNTS
	    DCOUNTS(NDATA)=IESDS
	    TTHS(NDATA)=TTH
	    XY_PIX(1,NDATA)=X
	    XY_PIX(2,NDATA)=Y
	    IF( LDEBUG ) WRITE(10,'(A,I5)') 'DEBUG: Added to local array',NDATA
	  ENDIF
C END OF MAIN LOOP
	ENDDO
C
1000  CLOSE(UNIT=1)
      RETURN
C
900	WRITE(10,'(/,A,I6,2A)') 'ERROR: More than',NREAD_MAX,' valid spots in ',FILE_NAME
	STOP 'ERROR: Too many spots in a single data file'
	END


	SUBROUTINE READ_GEASC_FILE(HKLS,XY_PIX,WAVS,COUNTS,DCOUNTS,TTHS,
	1					NDATA,NMULT2,NMULT3,NFAIL,NREJ, WAV_LO,FILE_NAME)
C
	CHARACTER FILE_NAME*(*)
C
	PARAMETER (NREAD_MAX=20000)
	INTEGER HKLS(3,NREAD_MAX)
	REAL XY_PIX(2,NREAD_MAX),TTHS(NREAD_MAX),WAVS(NREAD_MAX)
	REAL COUNTS(NREAD_MAX),DCOUNTS(NREAD_MAX)
C
	LOGICAL IS_REJECT_HKL
C
	CHARACTER LINE*100
	INTEGER FLAGS(3),NOD(3),IHKL(3),CNTS1(6),ESDS1(6),CNTS2(6),ESDS2(6)
C
	OPEN(UNIT=1,FILE=FILE_NAME,STATUS='OLD',ERR=900)
C
C Read in various parameters from the header
C
100	LINE=REPEAT(' ',LEN(LINE))
	READ(1,'(Q,A)',END=901) ILEN,LINE(1:ILEN)
C All values except YSCA are in mm
	IF(LINE(1:5) .EQ. 'CTOF ') READ(LINE,'(5X,F10.3)') CTOF ! drum radius
	IF(LINE(1:5) .EQ. 'RAST ') READ(LINE,'(5X,F10.3)') RAST ! pixel size
	IF(LINE(1:5) .EQ. 'YSCA ') READ(LINE,'(5X,F12.5)') YSCA ! Y/X pixel ratio
	IF(LINE(1:5) .EQ. 'XC_S ') READ(LINE,'(5X,F12.5)') XC_S ! X center offset
	IF(LINE(1:5) .EQ. 'YC_S ') READ(LINE,'(5X,F12.5)') YC_S ! Y center offset
	IF(LINE(1:5) .EQ. 'FIDX ') READ(LINE,'(5X,F12.5)') FIDX ! nom. X center
	IF(LINE(1:5) .EQ. 'FIDY ') READ(LINE,'(5X,F12.5)') FIDY ! nom. Ycenter
C If we haven't reached the data section, continue reading the header
	IF(LINE(1:15) .NE. 'REFLECTION DATA') GOTO 100
C
C Read the number of reflections, and die and complain if the wrong format
C
	IF(LINE(23:34) .NE. ' REFLECTIONS') THEN
	  WRITE(10,'(/,2A)') 'ERROR: Invalid data file header in ',FILE_NAME
	  STOP 'ERROR: Invalid data file header'
	ENDIF
C
	READ(LINE(16:22),'(I7)') NREAD
	IF(NREAD .GT. NREAD_MAX) THEN
	  WRITE(10,'(/,A,I6,2A)') 'ERROR: More than',NREAD_MAX,' valid spots in ',FILE_NAME
	  STOP 'ERROR: Too many spots in a single data file'
	ENDIF
C
C Extra two lines for Jacqui, made by LADIABS
C
ccc	READ(1,*)
ccc	READ(1,*)
C
C MAIN LOOP: Read in data and keep a tally of valid, > 2 multiples, and failed spots
C
	NDATA=0
	NREJ=0
	NFAIL=0
	NMULT2=0
	NMULT3=0
	DO IREAD=1,NREAD
C
C IHKL is the Miller indices
C X & Y are the spot positions in cm relative to the image center
C WAV is the wavelength of the spot (the lowest value if IMULT>1)
C IMULT is the wavelength multiplicity of the spot
C FLAGS, NOD & NIDX are not used by argonne_boxes. I don't know what they are!
	  READ(1,'(3I5,2F12.2,F8.3,I7,3I2,4I5)') IHKL,X,Y,WAV,IMULT,FLAGS,NOD,NIDX
C If IMULT=2, WAV2 is the second wavelength.
C IPOINT & DMINTH are not used by argonne_boxes. I don't know what they are!
	  READ(1,'(E12.4,I6,E12.4)') WAV2,IPOINT,DMINTH
C CNTS2(1) & ESDS2(1) contain the spot intensity and esd (we only use the first one)
	  READ(1,'(6I7,6I6)') CNTS1,ESDS1
	  READ(1,'(6I7,6I6)',IOSTAT=IERR) CNTS2,ESDS2
C On an error try to reread the line in case the first intensity is I8 not I7
C Garry's absorption correction program can do this to the files
	  IF(IERR .NE. 0) THEN
	    BACKSPACE(1)
	    READ(1,'(I8,5I7,6I6)',IOSTAT=IERR) CNTS2,ESDS2
	  ENDIF
C If spot is valid add to the data list, else update rejection statistics
	  IF( IS_REJECT_HKL(IHKL,0) ) THEN		! Reject spot if on the rejection list
	    NREJ=NREJ+1
	  ELSE IF(ESDS2(1) .LT. 0.0) THEN		! Reject spot if integration failure
	    NFAIL=NFAIL+1
	  ELSE IF(IMULT .GE. 3) THEN			! Reject spot if 3 or more multiples
	    NMULT3=NMULT3+1
	  ELSE IF(IMULT.EQ.2 .AND.			! Reject if 2 multiples and
	1		WAV/2.0.LT.WAV_LO) THEN		! wav/2 is not within limits
	    NMULT2=NMULT2+1
	  ELSE						! Add to the data lists
	    NDATA=NDATA+1
	    HKLS(1,NDATA)=IHKL(1)
	    HKLS(2,NDATA)=IHKL(2)
	    HKLS(3,NDATA)=IHKL(3)
	    WAVS(NDATA)=WAV
	    COUNTS(NDATA)=CNTS2(1)
	    DCOUNTS(NDATA)=ESDS2(1)
C Calculate the TTH from the spot position
	    RADIUS=CTOF
	    PHI=X/RADIUS
	    X2=RADIUS*SIN(PHI)
	    Y2=RADIUS*COS(PHI)
	    Z2=Y
	    SIZE=SQRT(X2**2 + Y2**2 + Z2**2)
	    DOT=Y2/SIZE
	    TTHS(NDATA)=ACOSD(DOT)
C Convert X & Y in mm to pixel values
	    XY_PIX(1,NDATA)=(X+FIDX+XC_S)/RAST
	    XY_PIX(2,NDATA)=(Y+FIDY+YC_S)/RAST*YSCA
C
	  ENDIF
C END OF MAIN LOOP
	ENDDO
C
	CLOSE(UNIT=1)
      RETURN
C
C Complain and die if something goes wrong
C
900	WRITE(10,'(/,2A)') 'ERROR: Cannot open file ',FILE_NAME
	STOP 'ERROR: Missing data file'
C
901	WRITE(10,'(/,2A)') 'ERROR: Cannot find "REFLECTION DATA" line in ',FILE_NAME
	STOP 'ERROR: Invalid data file'
	END
