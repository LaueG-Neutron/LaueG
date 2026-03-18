C----- Main routine to edit the data processing options -------------------
c	SUBROUTINE EDIT_OPTIONS
C----- Secondary routines to edit specific options ------------------------
c	SUBROUTINE ASK_TWIN_MODEL(ITWIN_OPT)
c	SUBROUTINE ASK_WAV_MODEL(LWAV_CORR,IWAV_OPT,NWAV,WAV_GAMMA)
c	SUBROUTINE ASK_EFF_MODEL(LEFF_CORR,IEFF_OPT)
c	SUBROUTINE ASK_ABSORB_MODEL(LABS_CORR,IABS_OPT,UR0,UR_LIN)
c	SUBROUTINE ASK_EXTI_MODEL(LEXTI_CORR,IEXTI_OPT)
c	SUBROUTINE ASK_PARAM_OPTS
C----- Routines to Initialise and Log the Options -------------------------
c	SUBROUTINE SET_DEFAULT_OPTIONS
c	SUBROUTINE LOG_OPTIONS
C----- Routines to Read/Write the Options File ----------------------------
c	SUBROUTINE CREATE_OPTIONS_FILE
c	SUBROUTINE READ_OPTIONS_FILE
c	SUBROUTINE READ_OPTIONS_FILE_V2
C--------------------------------------------------------------------------


C----- Main routine to edit the data processing options -------------------

	SUBROUTINE EDIT_OPTIONS
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	COMMON /WAV_CORR_COM/ IWAV_OPT,NWAV,WAV_GAMMA,WAV_MIN,WAV_STEP,WAV_PAR(20)
	COMMON /EFF_CORR_COM/ IEFF_OPT,EFF_POLY(10)
	COMMON /EXTI_CORR_COM/ IEXTI_OPT,EXTI
	COMMON /ABS_CORR_COM/ IABS_OPT,UR0,UR_LIN
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
      LOGICAL LTWIN1,LTWIN2,LTWIN11,LTWIN12,LRATIO
      COMMON /TWINS_CORR_COM/ ITWIN_OPT,TWIN_RATIO,LTWIN1,LTWIN2,LTWIN11,LTWIN12,LRATIO
C
	CHARACTER CELL_TYPE*40,STR*80
C
C List the intensity correction factors that can be turned on or off
C
1	PRINT '(/,1X,A)','Intensity correction factors:'
C
	IF( LFILE_CORR ) THEN
	  PRINT *,'  1  Image scale factors: Refined'
	ELSE
	  PRINT *,'  1  Image scale factors: Off'
	ENDIF
C
	IF( LWAV_CORR ) THEN
	  PRINT *,'  2  Wavelength spectra : Refined'
	ELSE
	  PRINT *,'  2  Wavelength spectra : Off'
	ENDIF
C
	IF( LEFF_CORR ) THEN
	  PRINT *,'  3  Efficiency         : Refined'
	ELSE
	  PRINT *,'  3  Efficiency         : Off'
	ENDIF
C
	IF( LEXTI_CORR ) THEN
	  PRINT *,'  4  Extinction         : Refined'
	ELSE
	  PRINT *,'  4  Extinction         : Off'
	ENDIF
C
	IF( LHARM_CORR ) THEN
	  PRINT *,'  5  Harmonics (wav/2)  : Corrected'
	ELSE
	  PRINT *,'  5  Harmonics (wav/2)  : Off'
	ENDIF
C
	IF( LABS_CORR ) THEN
	  PRINT *,'  6  Sample absorption  : Refined'
	ELSE
	  PRINT *,'  6  Sample absorption  : Off'
	ENDIF
C
C Second, list various models needed for merging or corrections
C
	PRINT *,'Models used for data processing:'
C
	CALL GET_CELL_TYPE(CELL_TYPE)
	PRINT   *,'  7  HKL merging rules  : '//TRIM(CELL_TYPE)
C
	STR='Off'
	IF( LWAV_CORR ) THEN
	  IF(IWAV_OPT .EQ. 1) WRITE(STR,'(I2,A)') NWAV,' point spline'
	  IF(IWAV_OPT .EQ. 2) WRITE(STR,'(I2,A,I5,A)') NWAV,
	1	' point spline + non-parametric (',NINT(WAV_GAMMA),')'
	ENDIF
	PRINT *,'  8  Wavelength model   : '//TRIM(STR)
C
	STR='Off'
	IF( LEFF_CORR ) THEN
	  IF(IEFF_OPT .EQ. 1) STR='Linear in Y'
	  IF(IEFF_OPT .EQ. 2) STR='Quadratic in Y'
	  IF(IEFF_OPT .EQ. 3) STR='Cubic in Y'
	  IF(IEFF_OPT .GT. 3) WRITE(STR,'(A,I2,A)') 'Polynomial(',IEFF_OPT,') in Y'
	ENDIF
	PRINT *,'  9  Efficiency model   : '//TRIM(STR)
C
	STR='Off'
	IF( LEXTI_CORR ) THEN
	  IF(IEXTI_OPT .EQ. 1) STR='Zach. Type 1 / SHELX'
	  IF(IEXTI_OPT .EQ. 2) STR='Zach. Type 2'
	  IF(IEXTI_OPT .EQ. 3) STR='Zach. Type 1 (+ B&C)'
	  IF(IEXTI_OPT .EQ. 4) STR='B&C Type 1,G'
	  IF(IEXTI_OPT .EQ. 5) STR='B&C Type 1,L'
	  IF(IEXTI_OPT .EQ. 6) STR='B&C Type 2'
	ENDIF
	PRINT *,' 10  Extinction model   : '//TRIM(STR)
C
	STR='Off'
	IF( LABS_CORR ) THEN
	  IF(IABS_OPT .EQ. 1) STR='Non-hydrogenous'
	  IF(IABS_OPT .EQ. 2) STR='Hydrogenous'
	  IF(IABS_OPT .EQ. 3) STR='General case'
	ENDIF
	PRINT *,' 11  Absorption model   : '//TRIM(STR)
C
C Third, provide access to internal parameters
C
	PRINT *,'Miscellaneous:'
	PRINT *,' 12  Edit internal parameters'
	PRINT *,' 13  Save options to "laue4.opt"'
	PRINT *,' 14  Read options from "laue4.opt"'
C
C Ask which parameters to change.
C
	PRINT '(1X,A,$)','Input 1 - 14 to change options [RETURN = skip]: '
	IOPT=0
	READ(*,'(I10)',ERR=1000) IOPT
	IF(IOPT .EQ. 0) GOTO 1000
C
C First handle the boolean values, just switch their values
C
C File scale factors
	IF(IOPT .EQ. 1) THEN
	  LFILE_CORR=.NOT.LFILE_CORR
C Wavelength correction refinement
	ELSE IF(IOPT .EQ. 2) THEN
	  LWAV_CORR=.NOT.LWAV_CORR
C Efficiency(Y & wav) corrections
	ELSE IF(IOPT .EQ. 3) THEN
	  LEFF_CORR=.NOT.LEFF_CORR
C Extinction correction refinement
	ELSE IF(IOPT .EQ. 4) THEN
	  LEXTI_CORR=.NOT.LEXTI_CORR
C WAV/2 correction
	ELSE IF(IOPT .EQ. 5) THEN
	  LHARM_CORR=.NOT.LHARM_CORR
C Sample absorption refinement
	ELSE IF(IOPT .EQ. 6) THEN
	  LABS_CORR=.NOT.LABS_CORR
C
C Next handle the more complex model parameters by calling procedures
C
C Merging: cell & Friedels
	ELSE IF(IOPT .EQ. 7) THEN
	  CALL ASK_CELL_TYPE
C Wavelength correction model
	ELSE IF(IOPT .EQ. 8) THEN
	  CALL ASK_WAV_MODEL(LWAV_CORR,IWAV_OPT,NWAV,WAV_GAMMA)
C Efficiency correction model
	ELSE IF(IOPT .EQ. 9) THEN
	  CALL ASK_EFF_MODEL(LEFF_CORR,IEFF_OPT)
C Extinction model
	ELSE IF(IOPT .EQ. 10) THEN
	  CALL ASK_EXTI_MODEL(LEXTI_CORR,IEXTI_OPT)
C Sample absorption model
	ELSE IF(IOPT .EQ. 11) THEN
	  CALL ASK_ABSORB_MODEL(LABS_CORR,IABS_OPT,UR0,UR_LIN)
C
C Finally the miscellaneous "options"
C
C Edit the internal parameters
	ELSE IF(IOPT .EQ. 12) THEN
	  CALL ASK_PARAM_OPTS
C Save options to a file
	ELSE IF(IOPT .EQ. 13) THEN
	  CALL CREATE_OPTIONS_FILE
C Read options from a file
	ELSE IF(IOPT .EQ. 14) THEN
	  CALL READ_OPTIONS_FILE
C
	ENDIF
C
C Handled all the valid options, so jump back to change more parameters
C
	GOTO 1
C
C Jump point to get out of "option asking" loop
C
1000  PRINT *
C
C Get options for processing twins (if twinning detected)
C ITWIN_OPT=0 if no twin data
C
      ITWIN_OPT=0
	IF(IDATA_MODE .EQ. 2) CALL ASK_TWIN_MODEL(ITWIN_OPT)
C Set which twin data types (=1,2,11,12) we will process
      LTWIN1 =(ITWIN_OPT .NE. 2) .AND. (ITWIN_OPT .NE. 6)
	LTWIN2 =(ITWIN_OPT .NE. 1) .AND. (ITWIN_OPT .NE. 6)
	LTWIN11=(ITWIN_OPT .GE. 4)
	LTWIN12=(ITWIN_OPT .GE. 5)
C Set twin ratio and if we need to refine it
      TWIN_RATIO=1.0
	LRATIO =(ITWIN_OPT .GE. 3)
C
	RETURN
	END


C ========== Secondary routines to edit specific options ===========

	SUBROUTINE ASK_TWIN_MODEL(ITWIN_OPT)
C
	PRINT *,'TWINNING options:'
	PRINT *,' 1  Use only Twin 1 separable spots'
	PRINT *,' 2  Use only Twin 2 separable spots'
	PRINT *,' 3  Use all separable spots'
	PRINT *,' 4  Use all separable spots and total of overlapped pairs'
	PRINT *,' 5  Use all separable and overlapped pairs'
	PRINT *,' 6  Use only overlapped pairs (not recommended)'
100	PRINT '(1X,A,$)','Select option 1 - 6: '
C
	READ(*,'(I10)',ERR=100) ITWIN_OPT
	IF(ITWIN_OPT.LT.1 .OR. ITWIN_OPT.GT.6) GOTO 100
      PRINT *
C
	RETURN
	END


	SUBROUTINE ASK_WAV_MODEL(LWAV_CORR,IWAV_OPT,NWAV,WAV_GAMMA)
C
	LOGICAL LWAV_CORR
C
C Ask for the model to use
C
	PRINT '(/,1X,A)','Wavelength distribution correction model'
	PRINT *,' 1  quadratic spline'
	PRINT *,' 2  quadratic spline + non-parametric correction'
	PRINT '(1X,A,$)','Input 1 - 2 to change model type: '
	READ(*,'(I1)',IOSTAT=IERR) IOPT
C
C Only change IWAV_OPT if a valid answer
C
	IF(IERR.NE.0 .OR. IOPT.LT.1 .OR. IOPT.GT.2) RETURN
	IWAV_OPT=IOPT
C
C Ask for model options
C
	PRINT '(1X,A,$)','   Number of points in spline [3-20]: '
	READ(*,'(I10)',IOSTAT=IERR) NWAV
	IF(IERR .NE. 0) NWAV=8
	NWAV=MAX(3,MIN(20, NWAV ))
C
	IF(IWAV_OPT .EQ. 2) THEN
	  PRINT '(1X,A,$)','   Gamma value for smoothing [~100-1000]: '
	  READ(*,'(F10.0)',IOSTAT=IERR) WAV_GAMMA
	  IF(IERR .NE. 0) WAV_GAMMA=300
	  WAV_GAMMA=MAX(1.0,MIN(9.99E4, WAV_GAMMA ))
	ENDIF
C
C Make sure the correction is turned on
C
	LWAV_CORR=.TRUE.
C
	RETURN
	END



	SUBROUTINE ASK_EFF_MODEL(LEFF_CORR,IEFF_OPT)
C
	LOGICAL LEFF_CORR
C
C Ask for the model to use
C
	PRINT '(/,1X,A)','Model for fitting Y efficiency variation'
	PRINT *,' 0  Constant'
	PRINT *,' 1  Linear'
	PRINT *,'>1  Polynomial of order 2 - 10'
	PRINT '(1X,A,$)','Input 0 - 10 to change model type: '
	READ(*,'(I)',IOSTAT=IERR) IOPT
C If an invalid answer, return without doing anything
	IF(IERR.NE.0 .OR. IOPT.LT.0 .OR. IOPT.GT.10) RETURN
	IEFF_OPT=IOPT
C
	LEFF_CORR= (IEFF_OPT .NE. 0)
C
	RETURN
	END



	SUBROUTINE ASK_ABSORB_MODEL(LABS_CORR,IABS_OPT,UR0,UR_LIN)
C
	LOGICAL LABS_CORR
C
C Ask for which model to use
C
	PRINT '(/,1X,A)','Isotropic absorption model (near-spherical sample)'
	PRINT *,'Select dominant atom type used for absorption calculations:'
	PRINT *,' 1  non-hydrogenous'
	PRINT *,' 2  hydrogenous'
	PRINT *,' 3  general case'
	PRINT '(1X,A,$)','Input 1 - 3 to change sample type: '
	READ(*,'(I1)',IOSTAT=IERR) IOPT
C If an invalid answer, just return without doing anything
	IF(IERR.NE.0 .OR. IOPT.LT.1 .OR. IOPT.GT.3) RETURN
C
C Set starting values for UR0 & UR_LIN
C
	IABS_OPT=IOPT
	UR0=0.1
	UR_LIN=0.1
	IF(IABS_OPT .EQ. 1) THEN
	  UR0=0.0
	ELSE IF(IABS_OPT .EQ. 2) THEN
	  UR_LIN=UR0*0.93
	ENDIF
C
C Make sure the correction is turned on
C
	LABS_CORR=.TRUE.
C
	RETURN
	END



	SUBROUTINE ASK_EXTI_MODEL(LEXTI_CORR,IEXTI_OPT)
C
	LOGICAL LEXTI_CORR
C
C Ask for which model to use
C
	PRINT '(/,1X,A)','Select model to use for isotropic extinction:'
	PRINT *,' 1  Zach. Type 1 / SHELX'
	PRINT *,' 2  Zach. Type 2'
	PRINT *,' 3  Zach. Type 1 (+ B&C)'
	PRINT *,' 4  B&C Type 1,G'
	PRINT *,' 5  B&C Type 1,L'
	PRINT *,' 6  B&C Type 2'
	PRINT '(1X,A,$)','Input 1 - 6 to change sample type: '
	READ(*,'(I1)',IOSTAT=IERR) IOPT
C If an invalid answer, just return without doing anything
	IF(IERR.NE.0 .OR. IOPT.LT.1 .OR. IOPT.GT.6) RETURN
C
C Turn on correction, store model type
C
	IEXTI_OPT=IOPT
	LEXTI_CORR=.TRUE.
C
	RETURN
	END



	SUBROUTINE ASK_PARAM_OPTS
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	COMMON /EFF_CORR_COM/ IEFF_OPT,EFF_POLY(10)
C
1	PRINT '(/,1X,A)','Internal parameters:'
	PRINT '(A,I3,F6.2)'   ,'  1  NSEQ_MIN,SEQ_FRAC=',NSEQ_MIN,SEQ_FRAC
	PRINT '(A,4F7.1)'     ,'  2  X_LO,X_HI,Y_LO,Y_HI=',X_LO,X_HI,Y_LO,Y_HI
	PRINT '(A,F6.3,F6.2)' ,'  3  WAV_LO,WAV_HI=',WAV_LO,WAV_HI
	PRINT '(A,F6.3)'      ,'  4  CELL_MULT=',CELL_MULT
	PRINT '(A,F5.2,F6.3)' ,'  5  ALL_SIG_MULT,ALL_REL_ERR=',
	2				ALL_SIG_MULT,ALL_REL_ERR
	PRINT '(A,2F7.1,F6.3)','  6  EFF_WIDTH,EFF_THICK=',
	1				EFF_WIDTH,EFF_THICK
	PRINT '(A,3F5.2)'     ,'  7  WAV_ERR1,WAV_ERR2,EXTI_ERR=',
	1				WAV_ERR1,WAV_ERR2,EXTI_ERR
	PRINT '(A,I5,F5.2)'   ,'  8  NBINS_NPAR,EXTI_CORR_MAX=',
	1				NBINS_NPAR,EXTI_CORR_MAX
	PRINT '(A,F8.3,F5.2)' ,'  9  HARM_CUTOFF_MULT,WEAK_SIG_MULT=',
	1				HARM_CUTOFF_MULT,WEAK_SIG_MULT
	PRINT '(A,I5)'        ,' 10  IVERBOSE=',IVERBOSE
C
C Ask which parameters to change.
C
	PRINT '(1X,A,$)','Input 1 - 10 to change options [RETURN = skip]: '
	IOPT=0
	READ(*,'(I11)',ERR=1000) IOPT
	IF(IOPT .EQ. 0) GOTO 1000
C
	PRINT '(1X,A,$)','Input parameter values: '
	IF     (IOPT .EQ. 1) THEN
	  READ(*,*) NSEQ_MIN,SEQ_FRAC
	ELSE IF(IOPT .EQ. 2) THEN
	  READ(*,*) X_LO,X_HI,Y_LO,Y_HI
	ELSE IF(IOPT .EQ. 3) THEN
	  READ(*,*) WAV_LO,WAV_HI
	ELSE IF(IOPT .EQ. 4) THEN
	  READ(*,*) CELL_MULT
	ELSE IF(IOPT .EQ. 5) THEN
	  READ(*,*) ALL_SIG_MULT,ALL_REL_ERR
	ELSE IF(IOPT .EQ. 6) THEN
	  READ(*,*) EFF_WIDTH,EFF_THICK
	ELSE IF(IOPT .EQ. 7) THEN
	  READ(*,*) WAV_ERR1,WAV_ERR2,EXTI_ERR
	ELSE IF(IOPT .EQ. 8) THEN
	  READ(*,*) NBINS_NPAR,EXTI_CORR_MAX
	ELSE IF(IOPT .EQ. 9) THEN
	  READ(*,*) HARM_CUTOFF_MULT,WEAK_SIG_MULT
	ELSE IF(IOPT .EQ. 10) THEN
	  READ(*,*) IVERBOSE
	ELSE
	  PRINT *,' Invalid option'
	  GOTO 1
	ENDIF
	PRINT *,'WARNING: Entered values are not checked for validity'
C
	GOTO 1
C
C
1000	RETURN
	END


C ================ Routines to Initialise and Log the Options ===================

	SUBROUTINE SET_DEFAULT_OPTIONS
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	COMMON /WAV_CORR_COM/ IWAV_OPT,NWAV,WAV_GAMMA,WAV_MIN,WAV_STEP,WAV_PAR(20)
	COMMON /EFF_CORR_COM/ IEFF_OPT,EFF_POLY(10)
	COMMON /EXTI_CORR_COM/ IEXTI_OPT,EXTI
	COMMON /ABS_CORR_COM/ IABS_OPT,UR0,UR_LIN
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
C Flags to enable various corrections
C
	LFILE_CORR=.TRUE.
	LWAV_CORR =.TRUE.
	LEFF_CORR =.TRUE.
	LEXTI_CORR=.FALSE.
	LABS_CORR =.FALSE.
	LHARM_CORR=.FALSE.
C
C Options parameters for correction models
C
	IWAV_OPT=2
	NWAV=10
	WAV_GAMMA=300.0
C
	IEFF_OPT=2
	DO I=1,10
	  EFF_POLY(I)=0.005
	ENDDO
C
	IEXTI_OPT=1
	EXTI=0.0				! must start as zero
C
	IABS_OPT=1
	UR0=0.0
	UR_LIN=0.1
C
C The "internal parameters"
C
	NSEQ_MIN=2
	SEQ_FRAC=0.75
C
	X_LO=3+50
	X_HI=3981-50
	Y_LO=1+50
	Y_HI=1976-50-20
C
	WAV_LO=0.85
	WAV_HI=1.70
C
	CELL_MULT=1.0
C
	ALL_REL_ERR=0.005
	ALL_SIG_MULT=1.0
C
	EFF_THICK=1.0	! amount of Gd relative to Koala IPs
	EFF_WIDTH=796.0	! half-radius of drum in pixels
C
	WAV_ERR1=0.10
	WAV_ERR2=0.10
	EXTI_ERR=0.20
C
	NBINS_NPAR=700
	EXTI_CORR_MAX=0.10
	WEAK_SIG_MULT=1.0
	HARM_CUTOFF_MULT=1.0
C
	IVERBOSE=0
C
	RETURN
	END


	SUBROUTINE LOG_OPTIONS
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	COMMON /WAV_CORR_COM/ IWAV_OPT,NWAV,WAV_GAMMA,WAV_MIN,WAV_STEP,WAV_PAR(20)
	COMMON /EFF_CORR_COM/ IEFF_OPT,EFF_POLY(10)
	COMMON /EXTI_CORR_COM/ IEXTI_OPT,EXTI
	COMMON /ABS_CORR_COM/ IABS_OPT,UR0,UR_LIN
C
      LOGICAL LTWIN1,LTWIN2,LTWIN11,LTWIN12,LRATIO
      COMMON /TWINS_CORR_COM/ ITWIN_OPT,TWIN_RATIO,LTWIN1,LTWIN2,LTWIN11,LTWIN12,LRATIO
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	CHARACTER CELL_TYPE*40,STEMP*30
C
	WRITE(10,'(/,A,/)') '====== Data processing options ======'
C
	CALL GET_CELL_TYPE(CELL_TYPE)
	WRITE(10,'(2A,/)') 'HKL merging rules: ',TRIM(CELL_TYPE)
C
C
	IF( LFILE_CORR ) THEN
	  WRITE(10,'(A)') 'Image scale factors: refined'
	ELSE
	  WRITE(10,'(A)') 'Image scale factors: fixed to 1'
	ENDIF
C
C
	IF( .NOT.LWAV_CORR ) THEN
	  WRITE(10,'(A)') 'Wavelength spectra model: not refined'
	ELSE IF(IWAV_OPT .EQ. 1) THEN
	  WRITE(10,'(A,I3,A)') 'Wavelength spectra model: ',NWAV,' point spline'
	ELSE IF(IWAV_OPT .EQ. 2) THEN
	  WRITE(10,'(A,I3,A,F6.1)') 'Wavelength spectra model: ',NWAV,
	1			' point spline + non-param, gamma =',WAV_GAMMA
	ELSE
	  STOP 'BUG(log_options): Invalid IWAV_OPT'
	ENDIF
C
C
	IF( LEFF_CORR ) THEN
	  IF(IEFF_OPT .EQ. 1) WRITE(10,'(A)') 'Efficiency model: Linear in Y'
	  IF(IEFF_OPT .EQ. 2) WRITE(10,'(A)') 'Efficiency model: Quadratic in Y'
	  IF(IEFF_OPT .GT. 2) WRITE(10,'(A,I2,A)') 'Efficiency model: Polynomial(',IEFF_OPT,') in Y'
	ELSE
	  WRITE(10,'(A)') 'Efficiency model: not refined'
	ENDIF
C
C
	IF( LEXTI_CORR ) THEN
	  IF(IEXTI_OPT .EQ. 1) WRITE(10,'(A)') 'Extinction correction: Zach. Type 1 / SHELX'
	  IF(IEXTI_OPT .EQ. 2) WRITE(10,'(A)') 'Extinction correction: Zach. Type 2'
	  IF(IEXTI_OPT .EQ. 3) WRITE(10,'(A)') 'Extinction correction: Zach. Type 1 (+ B&C)'
	  IF(IEXTI_OPT .EQ. 4) WRITE(10,'(A)') 'Extinction correction: B&C Type 1,G'
	  IF(IEXTI_OPT .EQ. 5) WRITE(10,'(A)') 'Extinction correction: B&C Type 1,L'
	  IF(IEXTI_OPT .EQ. 6) WRITE(10,'(A)') 'Extinction correction: B&C Type 2'
	  WRITE(10,'(2X,A,2P,F5.1,A)') 'Expected largest correction = ',EXTI_CORR_MAX,'%'
	ELSE
	  WRITE(10,'(A)') 'Extinction correction: off'
	ENDIF
C
C
	IF( LHARM_CORR ) THEN
	  WRITE(10,'(A)') 'Harmonics (wav/2) correction: on'
	ELSE
	  WRITE(10,'(A)') 'Harmonics (wav/2) correction: off'
	ENDIF
C
C
	IF( LABS_CORR ) THEN
	  WRITE(10,'(2A)') 'Absorption correction: ',
	1			'Isotropic model (near-spherical sample)'
	  IF(IABS_OPT .EQ. 1) THEN
	    WRITE(10,'(2X,A)') 'Absorption model: non-hydrogenous'
	    WRITE(10,'(5X,A,F6.2,A)') 'uR =',UR_LIN,' * wavelength(A)'
	  ELSE IF(IABS_OPT .EQ. 2) THEN
	    WRITE(10,'(2X,A)') 'Absorption model: hydrogenous'
	    WRITE(10,'(5X,A,F6.2,A)') 'uR =',UR0,' * ( 1 + 0.93 * wavelength(A) )'
	  ELSE IF(IABS_OPT .EQ. 3) THEN
	    WRITE(10,'(2X,A)') 'Absorption model: general case'
	    WRITE(10,'(5X,A,2(F6.2,A))') 'uR =',UR0,' +',UR_LIN,' * wavelength(A)'
	  ELSE
	    STOP 'BUG(log_options): Invalid IABS_OPT'
	  ENDIF
	ELSE
	  WRITE(10,'(A)') 'Absorption correction: off'
	ENDIF
C
C Info about processing twins
	IF(ITWIN_OPT .GT. 0) THEN
C Create list of N values to be used
		STEMP=' (N='
	  IF( LTWIN1 ) STEMP=TRIM(STEMP)//' 1'
	  IF( LTWIN2 ) STEMP=TRIM(STEMP)//' 2'
	  IF( LTWIN11 ) STEMP=TRIM(STEMP)//' 11'
	  IF( LTWIN12 ) STEMP=TRIM(STEMP)//' 12'
	  STEMP=TRIM(STEMP)//')'
C Output what spots are processed and if twin ratio is refined
		WRITE(10,'(/,A)') 'TWINNED options:'
		IF(ITWIN_OPT .EQ. 1) WRITE(10,'(2X,2A)') 'Processing twin 1 separable spots',TRIM(STEMP)
		IF(ITWIN_OPT .EQ. 2) WRITE(10,'(2X,2A)') 'Processing twin 2 separable spots',TRIM(STEMP)
		IF(ITWIN_OPT .EQ. 3) WRITE(10,'(2X,2A)') 'Processing twin 1 & 2 separable spots',TRIM(STEMP)
		IF(ITWIN_OPT .EQ. 4) WRITE(10,'(2X,2A)') 'Processing separable spots, and sum of overlapped pairs',TRIM(STEMP)
		IF(ITWIN_OPT .EQ. 5) WRITE(10,'(2X,2A)') 'Processing separable spots, and overlapped spots pairs',TRIM(STEMP)
		IF(ITWIN_OPT .EQ. 6) WRITE(10,'(2X,2A)') 'Processing overlapped intensities',TRIM(STEMP)
		IF( LRATIO ) WRITE(10,'(2X,A)') 'Refining twinned pair ratio'
	ENDIF
C
C List the miscellaneous "internal paramaeters"
	WRITE(10,'(/,A)') 'Internal Parameters:'
	WRITE(10,'(2X,A,I3,F6.2)') 'NSEQ_MIN,SEQ_FRAC=',NSEQ_MIN,SEQ_FRAC
	WRITE(10,'(2X,A,4F7.1)') 'X_LO,X_HI,Y_LO,Y_HI=',X_LO,X_HI,Y_LO,Y_HI
	WRITE(10,'(2X,A,F6.3,F6.2)') 'WAV_LO,WAV_HI=',WAV_LO,WAV_HI
	WRITE(10,'(2X,A,F6.3)') 'CELL_MULT=',CELL_MULT
	WRITE(10,'(2X,A,F5.2,F6.3)') 'ALL_SIG_MULT,ALL_REL_ERR=',
	1				ALL_SIG_MULT,ALL_REL_ERR
	WRITE(10,'(2X,A,F7.1,F6.3)') 'EFF_WIDTH,EFF_THICK=',
	1				EFF_WIDTH,EFF_THICK
	WRITE(10,'(2X,A,3F5.2)') 'WAV_ERR1,WAV_ERR2,EXTI_ERR=',
	1				WAV_ERR1,WAV_ERR2,EXTI_ERR
	WRITE(10,'(2X,A,I5,F5.2)') 'NBINS_NPAR,EXTI_CORR_MAX=',
	1				NBINS_NPAR,EXTI_CORR_MAX
	WRITE(10,'(2X,A,F8.3,F5.2)') 'HARM_CUTOFF_MULT,WEAK_SIG_MULT=',
	1				HARM_CUTOFF_MULT,WEAK_SIG_MULT
	WRITE(10,'(2X,A,I5)') 'IVERBOSE=',IVERBOSE
C
	IF(IVERBOSE .NE. 0) WRITE(10,'(/,A,/)') 'Running in VERBOSE mode'
C
	RETURN
	END


C ================ Routines to Read/Write the Options File ===================

	SUBROUTINE CREATE_OPTIONS_FILE
C
C Update the version number when we change the format of the laue4.opt file
C
	DATA IVERSION/4/
C
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
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	OPEN(UNIT=2,FILE='laue4.opt',STATUS='UNKNOWN')
C
	WRITE(2,'(A,/,I8)') 'LAUE4 OPTION FILE VERSION:',IVERSION
C
C Flags to refine & turn on various corrections
C
	WRITE(2,'(A,/,L8)') 'LFILE_CORR:',LFILE_CORR
	WRITE(2,'(A,/,L8)') 'LWAV_CORR:',LWAV_CORR
	WRITE(2,'(A,/,L8)') 'LEFF_CORR:',LEFF_CORR
	WRITE(2,'(A,/,L8)') 'LEXTI_CORR:',LEXTI_CORR
	WRITE(2,'(A,/,L8)') 'LHARM_CORR:',LHARM_CORR
	WRITE(2,'(A,/,L8)') 'LABS_CORR:',LABS_CORR
C
C Merge rules
C
	WRITE(2,'(A,/,3I8)') 'ITYPE,IFRIEDEL,ITYPE_ORIG:',ITYPE,IFRIEDEL,ITYPE_ORIG
C
C Options parameters for correction models
C
	WRITE(2,'(A,/,2I8,F12.5)') 'IWAV_OPT,NWAV,WAV_GAMMA:',IWAV_OPT,NWAV,WAV_GAMMA
	WRITE(2,'(A,/,I8)') 'IEFF_OPT:',IEFF_OPT
	WRITE(2,'(A,/,I8)') 'IEXTI_OPT:',IEXTI_OPT
	WRITE(2,'(A,/,I8,2F12.5)') 'IABS_OPT,UR0,UR_LIN:',IABS_OPT,UR0,UR_LIN
C
C The "internal parameters"
C
	WRITE(2,'(A,/,I8,F12.5)') 'NSEQ_MIN,SEQ_FRAC:',NSEQ_MIN,SEQ_FRAC
	WRITE(2,'(A,/,4F12.5)') 'X_LO,X_HI,Y_LO,Y_HI:',X_LO,X_HI,Y_LO,Y_HI
	WRITE(2,'(A,/,2F12.5)') 'WAV_LO,WAV_HI:',WAV_LO,WAV_HI
	WRITE(2,'(A,/,F12.5)') 'CELL_MULT:',CELL_MULT
	WRITE(2,'(A,/,3F12.5)') 'ALL_SIG_MULT,ALL_REL_ERR:',
	1				ALL_SIG_MULT,ALL_REL_ERR
	WRITE(2,'(A,/,3F12.5)') 'dummy,EFF_WIDTH,EFF_THICK:',
	1				1000.0,EFF_WIDTH,EFF_THICK
	WRITE(2,'(A,/,3F12.5)') 'WAV_ERR1,WAV_ERR2,EXTI_ERR',
	1				WAV_ERR1,WAV_ERR2,EXTI_ERR
	WRITE(2,'(A,/,I8,F12.5)') 'NBINS_NPAR,EXTI_CORR_MAX:',
	1				NBINS_NPAR,EXTI_CORR_MAX
	WRITE(2,'(A,/,2F12.5)') 'HARM_CUTOFF_MULT,WEAK_SIG_MULT:',
	1				HARM_CUTOFF_MULT,WEAK_SIG_MULT
	WRITE(2,'(A,/,I8)') 'IVERBOSE:',IVERBOSE
C
	CLOSE(UNIT=2)
	RETURN
	END


	SUBROUTINE READ_OPTIONS_FILE
C
C Read the options file laue4.opt
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
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	CHARACTER STRING*50
C
C Open options file, read/check the header and version number
C
	OPEN(UNIT=1,FILE='laue4.opt',STATUS='OLD',ERR=900)
	STRING=REPEAT(' ',LEN(STRING))
	READ(1,'(Q,A)') ILEN,STRING(1:ILEN)
	IF(STRING(1:26) .NE. 'LAUE4 OPTION FILE VERSION:') GOTO 910
	READ(1,'(I8)',ERR=920) IVERSION
	IF(IVERSION.LT.1 .OR. IVERSION.GT.4) GOTO 910
C
C Give warning then handle the case of and older version
C
	IF(IVERSION .NE. 4) THEN
	  PRINT '(/,1X,A)','WARNING: The options file is an old version'
	  PRINT '(1X,A,/)','  Some incompatibilities may occur, check laue4.lis'
	  WRITE(10,'(/,A)') 'WARNING: Reading an old version options file'
	  WRITE(10,'(A,/)') 'WARNING: Some incompatibilities may arise'
	  IF(IVERSION .LE. 2) THEN
	    CLOSE(UNIT=1)
	    CALL READ_OPTIONS_FILE_V2
	    RETURN
	  ENDIF
	ENDIF
C
C Flags to refine & turn on various corrections
C
	READ(1,'(/,L8)',ERR=920) LFILE_CORR
	READ(1,'(/,L8)',ERR=920) LWAV_CORR
	READ(1,'(/,L8)',ERR=920) LEFF_CORR
	READ(1,'(/,L8)',ERR=920) LEXTI_CORR
	READ(1,'(/,L8)',ERR=920) LHARM_CORR
	READ(1,'(/,L8)',ERR=920) LABS_CORR
C
C Read in the merge rules
C Accept values if original cell types match (ignoring Friedels)
C
	READ(1,'(/,3I8)',ERR=920) ITYPE2,IFRIEDEL2,ITYPE_ORIG2
	IF(ITYPE_ORIG/2 .EQ. ITYPE_ORIG2/2) THEN
	  ITYPE=ITYPE2
	  IFRIEDEL=IFRIEDEL2
	  ITYPE_ORIG=ITYPE_ORIG2
	ELSE
	  PRINT *,'Ignoring inconsistent cell type in options file'
	ENDIF
C
C Options parameters for correction models
C
	READ(1,'(/,2I8,F12.5)',ERR=920) IWAV_OPT,NWAV,WAV_GAMMA
C Relative-shift wavelength correction removed for version 4
	IF(IVERSION .EQ. 3) THEN
	  IWAV_OPT=IWAV_OPT-1
	  IF(IWAV_OPT .EQ. 0) THEN
	    IWAV_OPT=1
	    NWAV=3
	  ENDIF
	ENDIF
	READ(1,'(/,I8)',ERR=920) IEFF_OPT
	READ(1,'(/,I8)',ERR=920) IEXTI_OPT
	READ(1,'(/,I8,2F12.5)',ERR=920) IABS_OPT,UR0,UR_LIN
C
C Any remaining "internal parameters"
C
	READ(1,'(/,I8,F12.5)',ERR=920) NSEQ_MIN,SEQ_FRAC
	READ(1,'(/,4F12.5)',ERR=920) X_LO,X_HI,Y_LO,Y_HI
	READ(1,'(/,2F12.5)',ERR=920) WAV_LO,WAV_HI
	READ(1,'(/,F12.5)',ERR=920) CELL_MULT
	READ(1,'(/,2F12.5)',ERR=920) ALL_SIG_MULT,ALL_REL_ERR
C Fixed EFF_YCEN for versions 3 & 4
	READ(1,'(/,3F12.5)',ERR=920) dummy,EFF_WIDTH,EFF_THICK
	READ(1,'(/,3F12.5)',ERR=920) WAV_ERR1,WAV_ERR2,EXTI_ERR
	READ(1,'(/,I8,F12.5)',ERR=920) NBINS_NPAR,EXTI_CORR_MAX
	READ(1,'(/,2F12.5)',ERR=920) HARM_CUTOFF_MULT,WEAK_SIG_MULT
	READ(1,'(/,I8)',ERR=920) IVERBOSE
C
	CLOSE(UNIT=1)
C
C Due to a BUG, width & thick may have be swapped around
C
	IF(EFF_WIDTH .LT. EFF_THICK) THEN
	  TMP=EFF_WIDTH
	  EFF_WIDTH=EFF_THICK
	  EFF_THICK=TMP
	ENDIF
C
	RETURN
C
C Error messages
C
900	STOP 'ERROR: Unable to find options file "laue4.opt"'
910	STOP 'ERROR: Options file "laue4.opt" is invalid or incorrect version'
920	STOP 'ERROR: Error reading options file "laue4.opt"'
C
	END



	SUBROUTINE READ_OPTIONS_FILE_V2
C
C Read a Version 2 options file laue4.opt as best as we can
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
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
	CHARACTER STRING*50
C
	OPEN(UNIT=1,FILE='laue4.opt',STATUS='OLD',ERR=900)
	STRING=REPEAT(' ',LEN(STRING))
	READ(1,'(Q,A)') ILEN,STRING(1:ILEN)
	IF(STRING(1:26) .NE. 'LAUE4 OPTION FILE VERSION:') GOTO 910
C
	READ(1,'(I8)',ERR=920) IVERSION
	IF(IVERSION.NE.1 .AND. IVERSION.NE.2) GOTO 910
C
C Flags to refine & turn on various corrections
C
	READ(1,'(/,L8)',ERR=920) LFILE_CORR
	READ(1,'(/,L8)',ERR=920) LWAV_CORR
	READ(1,'(/,L8)',ERR=920) LEFF_CORR
	READ(1,'(/,L8)',ERR=920) LABS_CORR
	READ(1,'(/,L8)',ERR=920) LEXTI_CORR
	READ(1,'(/,L8)',ERR=920) LHARM_CORR
C
C Read in the merge rules
C For IVERSION =1, ignore the values
C              =2, accept values if original cell types match
C
	IF(IVERSION .EQ. 1) THEN
	  READ(1,'(/,2I8)',ERR=920) ITYPE2,IFRIEDEL2
	ELSE
	  READ(1,'(/,3I8)',ERR=920) ITYPE2,IFRIEDEL2,ITYPE_ORIG2
	  IF(ITYPE_ORIG/2 .EQ. ITYPE_ORIG2/2) THEN
	    ITYPE=ITYPE2
	    IFRIEDEL=IFRIEDEL2
	    ITYPE_ORIG=ITYPE_ORIG2
	  ELSE
	    PRINT *,'Ignoring inconsistent cell type in options file'
	  ENDIF
	ENDIF
C
C Options parameters for correction models
C
	READ(1,'(/,I8)',ERR=920) IEFF_OPT
	READ(1,'(/,I8,2F12.5)',ERR=920) IABS_OPT,UR0,UR_LIN
	READ(1,'(/,2I8)',ERR=920) IWAV_OPT,NWAV
	READ(1,'(/,I8,F12.5)',ERR=920) IEXTI_OPT,EXTI_CORR_MAX
C
C The "internal parameters"
C
	READ(1,'(/,I8,F12.5)',ERR=920) NSEQ_MIN,SEQ_FRAC
	READ(1,'(/,4F12.5)',ERR=920) X_LO,X_HI,Y_LO,Y_HI
	READ(1,'(/,2F12.5)',ERR=920) WAV_LO,WAV_HI
	READ(1,'(/,3F12.5)',ERR=920) ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT
	READ(1,'(/,4F12.5)',ERR=920) EXTI_ERR,DUMMY,WAV_ERR1,WAV_ERR2
	READ(1,'(/,F12.5)',ERR=920) CELL_MULT
	READ(1,'(/,3F12.5)',ERR=920) dummy,EFF_WIDTH,EFF_THICK
	READ(1,'(/,I8,F12.5)',ERR=920) NBINS_NPAR,WAV_GAMMA
C
	CLOSE(UNIT=1)
	RETURN
C
C Error messages
C
900	STOP 'ERROR: Unable to find options file "laue4.opt"'
910	STOP 'ERROR: Options file "laue4.opt" is invalid or incorrect version'
920	STOP 'ERROR: Error reading options file "laue4.opt"'
C
	END
