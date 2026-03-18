	PROGRAM ARGONNE_BOXES
C
	COMMON /DEBUG_COM/ IDEBUG
C
	IDEBUG=-1									! if >0, output debug info for refln IDEBUG
C
	CALL MAIN(  '(build: 12/9/2024)'  )
C
C The following is used by LaueG to check the program didn't crash
C
	PRINT *,'SUCCESSFUL COMPLETION'
C
	CALL DELETE_FILE('___laueg_argonne_boxes.in')
C
	END

C ========================== Routines in this file ============================
C	SUBROUTINE MAIN(VERSION_STRING0)
C
C	SUBROUTINE DELETE_FILE(FILE_NAME)
C	SUBROUTINE QUIT(TEXT)
C =============================================================================


	SUBROUTINE MAIN(VERSION_STRING0)
C
C This version abandons the OLD mode option
C
C Unit numbers 1,2,3 for temporary I/O
C Unit numbers 30,40,50 for *.lis,*.ell,*.mdl files
C
	CHARACTER VERSION_STRING0*(*)
C
	CHARACTER VERSION_STRING*80
	COMMON /VERSION_STRING_COM/ VERSION_STRING
C
	COMMON /IZONES_COM/ IZONES(20000)
	COMMON /LIB_COM/ PLIB(1000,6),ILIB(1000),NLIB
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /SPOT2_COM/ IHKLM(4,20000),WAV(20000),TTH(20000),MULT(20000)
	COMMON /BGLOCAL_COM/ BKG(20000),BKG_ESD(20000),BKG_VAR(20000)
	COMMON /INTINT_COM/ RINTINT(20000),RSIGMA(20000),IOPTION(20000),IOVER(20000)
	COMMON /GEOM_COM/ PIX_CEN(2),PIX_SIZE(2),DRUM_RAD
	COMMON /MODELS_PFRAC_COM/ PFRAC_TARGET,PFRAC_ERROR,DPFRAC_ZERO,DPFRAC_GRAD
C
	COMMON /OPTIONS_COM/ IOPT_HKL,NHKLM
	COMMON /GAIN_COM/ GAIN,ESD_FACTOR,CROSSTALK,CROSS_ELLIPSE,CROSS_SMOOTH
C
	LOGICAL LSPACER,LCENTER
	CHARACTER BASE_NAME*80
	REAL RAD_MULT(3),ZCUTOFFS(5)
C
C
C Copy the version number string to COMMON, and print it out
C
	VERSION_STRING=VERSION_STRING0
	PRINT '(1X,A,5X,A)','ARGONNE_BOXES for LaueG ',TRIM(VERSION_STRING)
	PRINT *
C
C Delete the LaueG output file, it is exists
C
	OPEN(UNIT=1,FILE='___laueg_argonne_boxes.out',STATUS='OLD',IOSTAT=IDUMMY)
	CLOSE(UNIT=1,STATUS='DELETE',IOSTAT=IDUMMY)
C
C Read and setup various parameters
C
	CALL ARGBOX_SETUP_PARS(BASE_NAME, RAD_MULT,CONTOUR,AREA_FAC,TOLER_NBOUR,
	1			FFILL,LCENTER, NMODEL_CEN,NMODEL_OUT, SCALE,A_MAX,R_SMOOTH, IPIXCEN)
C
	CALL ARGBOX_SETUP_DATA(ZCUTOFFS, BASE_NAME, NUMX,NUMY,NREFS, IOPT_HKL, IMAGE_ZERO, IPIXCEN)
C
C Do spot centering, if requested
C
	IF( LCENTER ) CALL REFINE_SPOTS_CENTERS(RAD_MULT,ZCUTOFFS,TOLER_NBOUR,A_MAX,NREFS)
C
C Loop over spots trying to find good ones to add to the model library
C
100	WRITE(6,'(1X,A)') 'Finding model spots .....'
	WRITE(30,'(//,A,/)') '================== Finding Model Spots ====================='
C
C----------------------------------------------------
C Start of loop over reflections to find model spots
C
C Zero the counter for models in the model library
	NLIB=0
C Set LSPACER to prevent an initial print out of the spacer line
	LSPACER=.TRUE.
C Loop through spots adding any that meet the "model" criterion
	DO IR=1,NREFS
	  CALL ADD_MODEL_SPOT(A_MAX,FFILL,TOLER_NBOUR,CONTOUR,RAD_MULT,ZCUTOFFS, IR)
	ENDDO
C
C End of loop over reflections to find model spots
C----------------------------------------------------
C
C Explain which reflections rejected without any output
C
	IF(IOPT_HKL .NE. 3) THEN
	  WRITE(30,'(/,A)') 'Unlisted spots are below the zone intensity cutoff'
	ELSE
	  WRITE(30,'(/,A)') 'Unlisted spots are satellites or are below the zone intensity cutoffs'
	ENDIF
C
C Scale the peak ellipses to give the required peak fraction of PFRAC_TARGET,
C which may remove some spots. The library is then pruned of weaker models to
C give approximately NMODEL_CEN and NMODEL_OUT models per zone.
C
	WRITE(30,'(//,A)') '=============== Optimising Model Parameters ================'
	CALL SET_PEAK_FRAC(RAD_MULT,PFRAC_TARGET,TOLER_NBOUR)
	CALL PRUNE_MODELS(NMODEL_CEN,NMODEL_OUT)
C
C Give summary, complain and die if not enough models
C
	WRITE(30,'(/,A,I5,A)') 'Total of',NLIB,' model spots remain'
	IF(NLIB .LT. 6) CALL QUIT('ERROR: Less than 6 model spots were found')
C
C Smooth or fit the ellipse parameters for each model
C
	CALL FIT_MODEL_PARAMS(RAD_MULT,TOLER_NBOUR)
C
C Write the final model parameters to the *.mdl file
C
	OPEN(UNIT=50,FILE=TRIM(BASE_NAME)//'.mdl',STATUS='UNKNOWN')
	WRITE(50,'(A)') 'Model Spot   xcen ycen      (e,f,h * 100)    Pfrac'
	DO I=1,NLIB
	  WRITE(50,'(I4,I6,2X,2I5,2X,3F6.2,F8.3)') I,ILIB(I),
	1	(NINT(PLIB(I,K)),K=1,2),(PLIB(I,K)*100.0,K=3,5),PLIB(I,6)
	ENDDO
	CLOSE(50)
C
C Write banners that we are starting to integrate spots
C
	WRITE(6,*) 'Integrating spots .....'
	WRITE(30,'(//,A,/)') '================== Integrating Spots ====================='
C
C Open *.ell file on unit 40 to store elllipse information from each integration
C
	OPEN(UNIT=40,STATUS='UNKNOWN',FILE=TRIM(BASE_NAME)//'.ell')
	IF(IOPT_HKL .EQ. 1) THEN
	  WRITE(40,'(A,//,A)') 'LaueG Ellipse File Version 1, Option 1',
	1			'   h   k   l  xcen   ycen         Ellipse parameters,           Areas, Type'
	ELSE
	  WRITE(40,'(A,//,A)') 'LaueG Ellipse File Version 1, Option 2',
	1			'   h   k   l  xcen   ycen         Ellipse parameters,           Areas,  Type,m'
	ENDIF
C
C Do the actual integrations
C
	CALL INTEGRATE_ALL(NREFS, RAD_MULT,TOLER_NBOUR)
C
C Output the results
C
	CALL OUTPUT_RESULTS(BASE_NAME,SCALE,NREFS)
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
	LOGICAL LOPEN
C
	INQUIRE(UNIT=30,OPENED=LOPEN)
	IF( LOPEN ) WRITE(30,'(A)') TEXT
C
	PRINT *,TEXT
	CALL EXIT()
C
	END
