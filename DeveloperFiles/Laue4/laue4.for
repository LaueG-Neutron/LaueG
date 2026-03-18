C ====================== Routines in file =====================
C	PROGRAM LAUE4
C	SUBROUTINE PRINT_BANNERS(VERSION)
C	SUBROUTINE GET_RUN_MODE
C	SUBROUTINE GET_DATA_MODE(LTWINS,LMODUL)
C	SUBROUTINE FINISH_FILES
C	SUBROUTINE DELETE_FILE(FILE_NAME)
C	SUBROUTINE QUIT(TEXT)
C =============================================================


	PROGRAM LAUE4
C
C Output banners with version date to the console and log file
C
	CALL PRINT_BANNERS('(build: 5/9/2024)')
C=========================== TO DO ===========================
C
C Twins are in the test phase (I have actually given up on them)
C For modulations, outliers and maybe other statistics are incorrect
C
C==================== Main array sizes =======================
c	PARAMETER (NFILES_MAX=400)			! Max no of data files
c	PARAMETER (NREAD_MAX=20000)			! Max no of hkls in a data file
c	PARAMETER (NDATA_MAX=2000000)		! Max no of all hkls (before pruning)
c	PARAMETER (NSEQ_MAX=100000)			! Max no of unique hkls (sequences of equivalents)
c	PARAMETER (NPARS_MAX=500)			! Max no parameters refined by LSQ
c	PARAMETER (IRW=20000000)			! Size of work array used by LSQ
c	PARAMETER (NREJ_MAX=10000)			! Max no reject spots
c	PARAMETER (NDEBUG_MAX=100)			! Max no debug spots
c	PARAMETER (NBINS_MAX=5000)			! Max no bins (-20) for non-param wav fit
C=========================== I/O units =======================
c	1   Temporary input file
c  2-3   Temporary output files
c   10  Main log file 'laue4.lis'
c   11  CIF output file 'laue4.cif'
C========================== Some Commons =====================
C Turn on/off intensity corrections
c	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C True if parameters are refined in an intensity corrections
c	COMMON /LSQ_LPAR_COM/ LFILE,LWAV,LEFF,LABS,LEXTI,LTWIN_RATIO
C Defines number and position of refined parameters in LSQ parameter list
c	COMMON /PARS_NUM_COM/ NFILE,NWAV,NEFF,NABS,NEXTI,NTWIN
C Options and parameters (many refineable) for individual intensity corrections
c	COMMON /ABS_CORR_COM/ IABS_OPT,UR0,UR_LIN
c	COMMON /EFF_CORR_COM/ IEFF_OPT,EFF_POLY(10)
c	COMMON /EXTI_CORR_COM/ IEXTI_OPT,EXTI
c	COMMON /FILE_CORR_COM/ FSCALE(NFILES_MAX)
c	COMMON /TWINS_CORR_COM/ ITWIN_OPT,TWIN_RATIO,LTWIN1,LTWIN2,LTWIN11,LTWIN12,LRATIO
c	COMMON /WAV_CORR_COM/ IWAV_OPT,NWAV,WAV_GAMMA,WAV_MIN,WAV_STEP,WAV_PAR(20)
C Relates sequential numbering of files to actual number used in file name
c	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C Main data arrays of individual spots
c	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C Defines first/last indices for spots in a "sequence" and the merged/corrected intensity
c	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C Miscellaneous parameters
c	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI, ...
C=============================================================
C Arrays of packed data and parameters used by LSQ_FUNC() & NLSCON()
c	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
c	COMMON /LSQ_PARS_COM/ PARS(NPARS_MAX),ESDS(NPARS_MAX),COVAR(NPARS_MAX,NPARS_MAX),NPARS

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Types of data files (LaueG/lauegen) and intensity type (HKL/Twin/Modulated)
c	COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C Symmetry and cell parameters
c	COMMON /CELL_TYPE_COM/ ITYPE,IFRIEDEL,ITYPE_ORIG,ICEN,CELL(6)
C Modulation vectors and indices
c	COMMON /MODULATE_COM/ NMOD_VEC,NMOD_IDX,MOD_VEC(3,100),MOD_IDX(100,100)
C Defines spot for which verbose debug information is output
c	COMMON /DEBUG_COM/ NAMES_DEBUG,HKL_DEBUG,NDEBUG
C Save info from the header of the first TIF file
c	COMMON /FILE_INFO_COM/ NUMX,NUMY,HOST_NAME,DATE
C For base file-name of input files
c	COMMON /FILE_NAME_COM/ BASE_NAME,LSTRIP
c	COMMON /HARM_INFO_COM/ HARM_MAX,NHKLS_1ESD
c	COMMON /READ_TIFF_COM/ BUFFER1,IREC,IUNIT
c	COMMON /REJECT_COM/ NREJ,HKLM_REJECT(4,NREJ_MAX)
c	COMMON /SORT_BIN_COM/ BIN_DATA(NDATA_MAX)
c	COMMON /SORT_WEAK_COM/ CSEQ(NSEQ_MAX)
c	COMMON /TIFFINFO_COM/ COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM, ...
c	COMMON /WAV_FILE_COM/ NWAV_FILE,WAV_FILE_MIN,WAV_FILE_STEP,WAV_FILE(10000)
c	COMMON /WAV_NPAR_COM/ NWAV_NPAR,WAV_NPAR_MIN,WAV_NPAR_STEP,DOF,RAT(NBINS_MAX)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C Determine file modes: LaueG or Laugen
C
	CALL GET_RUN_MODE
C
C Set all options to their defaults
C
	CALL SET_DEFAULT_OPTIONS
C
C Read files to get needed parameters and file names
C
	CALL READ_PARAM_FILES
C
C Determine integrated intensity type: normal, twinned, modulated
C
      CALL GET_DATA_MODE
C
C Read default wavelength distribution file and clear the non-param correction
C
	CALL READ_WAV_FILE
	CALL CLEAR_WAV_NONPAR
C
C Allow manual editting of the options and log the options values used
C
	CALL EDIT_OPTIONS
	CALL LOG_OPTIONS
C
C Refine parameters on the strongest groups of equivalents
C
	CALL PREP_REFINE_DATA
	CALL PRELIM_REFINE
	CALL FULL_REFINE
	CALL FIT_WAV_NONPAR
	CALL FINAL_REFINE
C
C Apply refined parameters to correct all data
C
	CALL PREP_CORRECT_DATA
	CALL CORRECT_DATA
C
C Output useful information and statistics
C
	CALL OUTPUT_HKL_REDUNDANCY
	CALL OUTPUT_INTENSITY_STATS
	CALL OUTPUT_DSPACE_STATS
C
C Find and possibly reject outlier intensities
C
	CALL PROCESS_OUTLIERS
C
C Search for systematic absences and output results
C
	CALL OUTPUT_SYSTEM_ABSENCES
C
C Output corrected intensities plus various information files
C
	CALL OUTPUT_FINAL_RESULTS
C
C Output a final farewell and (potentially) delete the input file
C
	CALL FINISH_FILES
C
	END


	SUBROUTINE PRINT_BANNERS(VERSION)
C
	CHARACTER VERSION*(*)
C
C Output short banner to the console
C
	PRINT '(1X,2A,/)','LAUE4  Laue normalization software for neutrons ',
	1			VERSION
C
C Open log file and output a longer banner to it
C
	OPEN(UNIT=10,FILE='laue4.lis',STATUS='unknown')
	WRITE(10,'(A)')
	1	'==========================================================================',
	2	'  LAUE4  Laue normalization software for for neutrons '//VERSION,
	3	'               Method and software by Ross Piltz (ANSTO)',
     4    '       Please cite: R.O. Piltz, J. Appl. Cryst. (2018). 51, 635-645',
	5	'=========================================================================='
C
	RETURN
	END


	SUBROUTINE GET_RUN_MODE
C
C Saves the type of intensity files to /MODES_COM/
C
C LaueG mode (IFILE_MODE=1) indicated if the input file exists,
C otherwise Lauegen mode (IFILE_MODE=2) if any command line arguments,
C otherwise complain and die.
C
	COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
      LOGICAL LFOUND
C
	INQUIRE(FILE='___laueg_laue4.in',EXIST=LFOUND)
	IF(LFOUND) THEN
	  IFILE_MODE=1
	  WRITE(10,'(/,A)') 'Running in LaueG mode (LaueG input file found)'
	ELSE IF(NARGS() .GT. 1) THEN
	  IFILE_MODE=2
	  WRITE(10,'(/,A)') 'Running in Interactive Lauegen mode (line arguments found)'
	ELSE
	  WRITE(10,'(/,A)') 'ERROR: No parameter files or command line arguments'
	  STOP 'ERROR: No parameter files or command line arguments'
	ENDIF
C
	RETURN
	END


	SUBROUTINE GET_DATA_MODE
C
C Saves the type of data in intensity files to /MODES_COM/
C 	IDATA_MODE = 1 (normal mode), 2 (twin mode), 3 (modulated mode)
C Also, checks if data mode consistent with initial value from parameter file
C
C For IFILE_MODE=1, read the IDATA_MODE from the header of the first int file
C For IFILE_MODE=2, Lauegen mode which is always IDATA_MODE=1
C
C NB: In Laueg mode, must run before READ_DATA_FILES
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
      PARAMETER (NFILES_MAX=400)
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
	CHARACTER FILE_NAME*132
C
C DATA_MODE always =1 in Lauegen mode, so nothing more to do
C
	IF(IFILE_MODE .EQ. 2) RETURN
C
C In LauG mode, read data mode from first *.int file
C
	CALL MAKE_FILE_NAME(FILE_NAME, IFILE_NUM(1),'.int')
      CALL OPEN_INT_FILE(IDATA_MODE2,TRIM(FILE_NAME))
      CLOSE(UNIT=1)
C
C Output to console which data mode we are using
C
	IF(IDATA_MODE2 .EQ. 1) THEN
	  PRINT('(1X,A)'),'Integrated intensities use hkl indices of a reciprocal lattice'
	ELSEIF(IDATA_MODE2 .EQ. 2) THEN
	  PRINT('(1X,A)'),'Integrated intensities use hkl,N indices of a TWINNED lattice'
	ELSEIF(IDATA_MODE2 .EQ. 3) THEN
	  PRINT('(1X,A)'),'Integrated intensities use hkl,M indices of a MODULATED lattice'
      ENDIF
C
C Complain and die if parameter file uses a different data type
C
	IF(IDATA_MODE .NE. IDATA_MODE2) THEN
		CALL QUIT('ERROR: Inconsistent intensity types (i.e. twins or modulations)')
	ENDIF
C
	RETURN
      END


	SUBROUTINE FINISH_FILES
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	CHARACTER BASE_NAME*132,BUFF*132
C
C Output the final info
C
	PRINT '(/,1X,A)','See "laue4.lis" for more details'
      CALL MAKE_FILE_NAME(BASE_NAME, -1,'')
	PRINT '(1X,A)','laue4.lis & laue4.cif copied to '//TRIM(BASE_NAME)//'.*'
	WRITE(10,'(A)') 'laue4.lis & laue4.cif copied to '//TRIM(BASE_NAME)//'.*'
C
C This must be the last line of the log file so SCILAB knows it was successful
C
	WRITE(10,'(/,A)') '====== LAUE4 finished ======'
	CLOSE(UNIT=10)
C
C Copy laue4.lis to <BASE_NAME>.lis (for Max)
C
      OPEN(UNIT=10,FILE='laue4.lis',STATUS='OLD',ERR=200)
	OPEN(UNIT=11,FILE=TRIM(BASE_NAME)//'.lis',STATUS='unknown',ERR=200)
110   READ(10,'(A)',END=120) BUFF
      WRITE(11,'(A)') TRIM(BUFF)
      GOTO 110
120   CLOSE(UNIT=10)
      CLOSE(UNIT=11)
C
C Copy laue4.cif to <BASE_NAME>.cif (for Max)
C
200   OPEN(UNIT=10,FILE='laue4.cif',STATUS='OLD',ERR=300)
	OPEN(UNIT=11,FILE=TRIM(BASE_NAME)//'.cif',STATUS='unknown',ERR=300)
210   READ(10,'(A)',END=220) BUFF
      WRITE(11,'(A)') TRIM(BUFF)
      GOTO 210
220   CLOSE(UNIT=10)
      CLOSE(UNIT=11)
C
C If LaueG mode, (potentially) delete the input file
C
300	IF(IFILE_MODE .EQ. 1) CALL DELETE_FILE('___laueg_laue4.in')
C
	RETURN
	END


	SUBROUTINE DELETE_FILE(FILE_NAME)
C
C Delete the file if the "delete files" is missing or contains a TRUE
C
	CHARACTER FILE_NAME*(*)
C
	LOGICAL LDELETE
C
	OPEN(UNIT=17,FILE='___laueg_delete_files.in',STATUS='OLD',ERR=100)
	READ(17,*,ERR=100,END=100) LDELETE
	IF( .NOT.LDELETE ) RETURN
100	CLOSE(UNIT=17,IOSTAT=IDUMMY)
C
	OPEN(UNIT=17,FILE=FILE_NAME,STATUS='OLD',IOSTAT=IDUMMY)
	CLOSE(UNIT=17,STATUS='DELETE',IOSTAT=IDUMMY)
C
	RETURN
	END


	SUBROUTINE QUIT(TEXT)
C
C Workaround as SCILAB doesn't return STOP messages as they
C go to stderr not stdout
C
	CHARACTER*(*) TEXT
C
	PRINT *,TEXT
	CALL EXIT()
C
	END
