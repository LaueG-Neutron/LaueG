C ================ Routines in this file =========================
C	SUBROUTINE READ_PARS_FILE(BASE_NAME, CONTOUR,AREA_FAC,TOLER_NBOUR,FFILL,
C	1                                 LCENTER, NMODEL_CEN,NMODEL_OUT, IPIXCEN)
C	SUBROUTINE READ_HKLS_FILE(NLINES_HKLS)
C	SUBROUTINE READ_IMAGE_FILE(FILE_NAME)
C	SUBROUTINE GET_IMAGE_NUMXY(NUMX,NUMY)
C ================================================================


	SUBROUTINE READ_PARS_FILE(BASE_NAME, CONTOUR,AREA_FAC,TOLER_NBOUR,FFILL,
	1                                 LCENTER, NMODEL_CEN,NMODEL_OUT, IPIXCEN)
C
C Read the parameter section of the file "___laueg_argonne_boxes.in"
C
	LOGICAL LCENTER
	CHARACTER BASE_NAME*(*)
C
	COMMON /OPTIONS_COM/ IOPT_HKL,NHKLM
	COMMON /GEOM_COM/ PIX_CEN(2),PIX_SIZE(2),DRUM_RAD
	COMMON /GAIN_COM/ GAIN,ESD_FACTOR,CROSSTALK,CROSS_ELLIPSE,CROSS_SMOOTH
	COMMON /MODELS_PFRAC_COM/ PFRAC_TARGET,PFRAC_ERROR,DPFRAC_ZERO,DPFRAC_GRAD
	COMMON /EXCLUDE_COM/ NCIRC,ICIRC(3,999),NRECT,IRECT(4,999)
C
C Open LaueG input file
C
	OPEN(UNIT=1,FILE='___laueg_argonne_boxes.in',STATUS='OLD',ERR=900)
C
C Read the initial line which defines options and section sizes
C
	READ(1,*) IOPT_HKL,NLINES_PARS,NLINES_HKLS
	IF(IOPT_HKL.LT.1 .OR. IOPT_HKL.GT.3) GOTO 901
C
C Set number of HKLM indices to use for twin or incommensurate options
C
	NHKLM=3
	IF(IOPT_HKL .NE. 1) NHKLM=4
C
C Read the base part of the data file name
C
	READ(1,'(A)') BASE_NAME
C
	CALL OPEN_LIST_FILE(TRIM(BASE_NAME)//'.lis')
C
	READ(1,*) NMODEL_CEN,NMODEL_OUT,PFRAC_TARGET,PFRAC_ERROR,ICEN
	LCENTER=(ICEN .NE. 0)
C
	READ(1,*) CONTOUR,AREA_FAC
	READ(1,*) TOLER_NBOUR,FFILL
	READ(1,*) PIX_CEN,DRUM_PIX
C IPIXCEN=0 if pixel centre is unstable (VIVALDI & KOALA 1)
	READ(1,*) PIX_SIZE,IPIXCEN
	READ(1,*) CROSSTALK,CROSS_ELLIPSE,CROSS_SMOOTH
C Convert drum radius from multiples of X pixel size to mm
	DRUM_RAD=DRUM_PIX*PIX_SIZE(1)
C
	READ(1,*) NCIRC
	IF(NCIRC .GT. 999) GOTO 902
	DO I=1,NCIRC
	  READ(1,*) (ICIRC(K,I),K=1,3)
	ENDDO
	READ(1,*) NRECT
	IF(NRECT .GT. 999) GOTO 903
	DO I=1,NRECT
	  READ(1,*) (IRECT(K,I),K=1,4)
	ENDDO
C
C Run Sanity check on size of parameters section
C <<< Must update here if you add more parameters lines >>>
C
	IF(NLINES_PARS .NE. 9+NCIRC+NRECT) GOTO 904
C
	CLOSE(UNIT=1)
	RETURN
C
C Error messages:
C
900	CALL QUIT('ERROR: Can''t open LaueG input file')
901	CALL QUIT('ERROR: Unknown option in LaueG input file')
902	CALL QUIT('ERROR: Too many exclusion circles')
903	CALL QUIT('ERROR: Too many exclusion rectangles')
904	CALL QUIT('ERROR: Inconsistent LaueG input file')
	END


	SUBROUTINE READ_HKLS_FILE(NLINES_HKLS)
C
C Read the hkls section of the file "___laueg_argonne_boxes.in"
C
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /SPOT2_COM/ IHKLM(4,20000),WAV(20000),TTH(20000),MULT(20000)
	COMMON /OPTIONS_COM/ IOPT_HKL,NHKLM
	COMMON /DEBUG_COM/ IDEBUG
C
	WRITE(30,'(/,1A)') 'Reading reflection list from LaueG input file'
	OPEN(UNIT=1,FILE='___laueg_argonne_boxes.in',STATUS='OLD',ERR=900)
C
C Read header line and do sanity check on IOPT_HKL & NLINES_PARS values
C
	READ(1,*,ERR=901) IOPT_HKL,NLINES_PARS,NLINES_HKLS
	IF(IOPT_HKL.LT.1 .OR. IOPT_HKL.GT.3) GOTO 901
	IF(NLINES_HKLS .GT. 20000) CALL QUIT('ERROR: Number of reflections exceeds 20,000')
C
C Set flag for four indices to include twin or satellite number)
C
	IF(IOPT_HKL .EQ. 2) THEN
	  WRITE(6,'(1X,A)') 'Running in TWIN mode'
	  WRITE(30,'(A,/)') 'Running in TWIN mode'
	ELSEIF(IOPT_HKL .EQ. 3) THEN
	  WRITE(6,'(1X,A)') 'Running in SATELLITE mode'
	  WRITE(30,'(A)') 'Running in SATELLITE mode'
	  WRITE(30,'(A)') 'Reflections stored as H,K,L,M where M is the satellite number'
	  WRITE(30,'(A,/)') '  Non-satellite reflections denoted by M = 0'
	ENDIF
C
C Skip the parameter section
C
	DO I=1,NLINES_PARS
	  READ(1,*)
	ENDDO
C
C Read in the hkl list, and save limits in wavelength & d-spacing
C
	WAV_MIN=1E6
	WAV_MAX=0.0
	D_MIN=1E6
	DO I=1,NLINES_HKLS
	  READ(1,*,ERR=901) (IHKLM(K,I),K=1,NHKLM),X(I),Y(I),WAV(I),MULT(I),TTH(I)
	  WAV_MIN=MIN(WAV_MIN,WAV(I))
	  WAV_MAX=MAX(WAV_MAX,WAV(I))
	  D_MIN=MIN(D_MIN,WAV(I)/(2.0*SIND(TTH(I)/2.0)))
	ENDDO
C
	CLOSE(UNIT=1)
C
C Output info on the spots read and wavelength/d-spacing limits
C
	WRITE(6,'(I6,2A)') NLINES_HKLS,' spots read from LaueG input file'
	WRITE(30,'(2X,A,I5,A,F5.2,A,F6.2,A,F5.2,A)') 'Total of',NLINES_HKLS,
	1	' reflections  (wav =',WAV_MIN,' to',WAV_MAX,', d-min =',D_MIN,')'
C
C For "normal" reflections, zero the M index
C
	IF(IOPT_HKL .EQ. 1) THEN
	  DO I=1,NLINES_HKLS
	    IHKLM(4,I)=0
	  ENDDO
	ENDIF
C
C For twins, output number of separate and paired twins
C
	IF(IOPT_HKL .EQ. 2) THEN
C Count number of m=1,2,11,12 spots
	  M1=0
	  M2=0
	  M11=0
	  M12=0
	  DO I=1,NLINES_HKLS
	    IF(IHKLM(4,I) .EQ. 1) M1=M1+1
	    IF(IHKLM(4,I) .EQ. 2) M2=M2+1
	    IF(IHKLM(4,I) .EQ. 11) M11=M11+1
	    IF(IHKLM(4,I) .EQ. 12) M12=M12+1
	  ENDDO
C Output numbers of twins and possible overlaps
	  WRITE(30,'(5X,I5,A)') M1,' separate spots from twin 1'
	  WRITE(30,'(5X,I5,A)') M2,' separate spots from twin 2'
	  WRITE(30,'(5X,I5,A)') M11+M12,' spots in paired twins'
CCC	  IF(NLINES_HKLS .NE. M1+M2+M11+M12) GOTO 904
CCC	  IF(M11 .NE. M12) GOTO 905
	ENDIF
C
C For modulations, output number of main and satellite peaks
C
	IF(IOPT_HKL .EQ. 3) THEN
C Count number of main and satellite peaks
	  M0=0
	  M_MAX=0
	  DO I=1,NLINES_HKLS
	    IF(IHKLM(4,I) .EQ. 0) M0=M0+1
	    M_MAX=MAX(M_MAX, IHKLM(4,I) )
	    IF(IHKLM(4,I).LT.0 .OR. IHKLM(4,I).GT.99) CALL QUIT('ERROR: Satellite M indices not in range 0 to 99')
	  ENDDO
C Output numbers of normal and satellite peaks
	  WRITE(30,'(5X,I5,A)') M0,' main peaks'
	  WRITE(30,'(5X,I5,A,I3)') NLINES_HKLS-M0,' satellite peaks with M up to',M_MAX
C
	ENDIF
C
C Output extra info on debug reflection
C
	IF(IDEBUG .GT. 0) THEN
		WRITE(30,'(/,A,I5,A,I6,A)') '*** DEBUG: Refln ',IDEBUG, ' of ',NLINES_HKLS,' reflections'
		WRITE(30,'(A,4I5,2F7.1)') '*** DEBUG: IHKLM,X,Y=',(IHKLM(K,IDEBUG),K=1,4),X(IDEBUG),Y(IDEBUG)
		WRITE(30,'(A,F7.3,I4,F7.2)') '*** DEBUG: WAV,MULT,TTH=',WAV(IDEBUG),MULT(IDEBUG),TTH(IDEBUG)
		WRITE(30,'(A,2F8.1)') '*** DEBUG: Estimated BKG + ESD=',GET_GLOBAL_BKG(X(IDEBUG),Y(IDEBUG)),
	1							GET_GLOBAL_BKG_ESD(X(IDEBUG),Y(IDEBUG))
	ENDIF
C
	RETURN
C
C Error messages:
C
900	CALL QUIT('ERROR: Can''t open LaueG input file')
901	CALL QUIT('ERROR: Invalid or corrupted LaueG input file')
C
	END


	SUBROUTINE READ_IMAGE_FILE(IMAGE_ZERO, FILE_NAME)
C
C NB: This routine can be a hog, so put IY in outer loops
C
	CHARACTER FILE_NAME*(*)
C
	COMMON /GEOM_COM/ PIX_CEN(2),PIX_SIZE(2),DRUM_RAD
C
	COMMON /NUMXY_COM/ NUMX0,NUMY0
	COMMON /RIMAGE_COM/ RIMAGE(8000,2500)
C
	CHARACTER*2000 COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,ADC_INFO
	COMMON /TIFFINFO_COM/ COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,
	1			IOFF_DATA,NUMX,NUMY,IXYRES,ISTARTX,ISTARTY,ISPEED,
	2			EXPOSE_TIME,EXPOSE_PHI,TEMP_BEGIN,TEMP_END,TEMP_MIN,TEMP_MAX,
	3			IADC,ADC_ADD,ADC_DIV,ADC_INFO,IAPER_1,IAPER_2,IKAPPA
C
	LOGICAL LCYCLOPS
	INTEGER*2 IMAGE2_READ(8000*2500)
C
C Check the TIF file is valid and get some useful information from the header
C
	WRITE(6,'(1X,2A)') 'Reading image data from ',FILE_NAME
	WRITE(30,'(/,2A)') 'Reading image data from ',FILE_NAME
C
C Open the TIF file, load header info to /TIFFINFO_COM/, and complain
C and die if it is missing or invalid
C NB: The TIFF routines use UNIT=1 to read the file
C
	CALL CHECK_LAUE_TIFF(NUMX,NUMY, FILE_NAME)
	IF(NUMX .LT. 0) THEN
	  WRITE(30,'(/,A)') 'ERROR: Image data file is missing or corrupted'
	  CALL QUIT('ERROR: Image data file is missing or corrupted')
	ENDIF
	NUMX0=NUMX
	NUMY0=NUMY
C
C Output some of the header info and give up if the image is too large
C
	LCYCLOPS=(TRIM(HOST) .EQ. 'CYCLOPS')
	WRITE(30,'(2X,3A)') TRIM(HOST),' data created on ',TRIM(DATETIME)
	WRITE(30,'(2X,A,F7.1)') 'Exposure PHI value:',EXPOSE_PHI
	WRITE(30,'(2X,A,2(I5,A))') 'Image size:',NUMX,' x',NUMY,' pixels'
	IF(NUMX .GT. 8000) CALL QUIT('ERROR: Image > 8000 pixels wide')
	IF(NUMY .GT. 2500) CALL QUIT('ERROR: Image > 2500 pixels high')
C
C Do a fast read of the image data into IMAGE2_READ
C
	CALL READ_LAUE_TIFF_QUIET(IMAGE2_READ,NUMX,NUMY, FILE_NAME)
C
C Convert to REAL intensities and copy to /RIMAGE_COM/
C
	DO IY=1,NUMY
		DO IX=1,NUMX
			RIMAGE(IX,IY)=ZEXT( IMAGE2_READ( IX+(IY-1)*NUMX ) )
		ENDDO
	ENDDO
C
C If not CYCLOPS, determine the intensity zero level by averaging
C intensities over IX0,IY0 +/-5 pixels. If CYCLOPS, use a value of 0.
C
	IF( LCYCLOPS ) THEN
	  IMAGE_ZERO=0
	ELSE
C Search for clear spot in exit hole
	  IBEST=999999
	  DO IX=NINT(PIX_CEN(1))-50,NINT(PIX_CEN(1))+50,5
	    DO IY=NINT(PIX_CEN(2))-50,NINT(PIX_CEN(2))+50,5
	      IAVE=( RIMAGE(IX,IY)+RIMAGE(IX-10,IY)+RIMAGE(IX+10,IY)+
	1							RIMAGE(IX,IY-10)+RIMAGE(IX,IY+10) ) /5
	      IF(IAVE .LT. IBEST) THEN
			    IX5=IX
	        IY5=IY
	        IBEST=IAVE
		  ENDIF
	    ENDDO
	  ENDDO
C Make average over X & Y
	  SUM=0.0
	  DO IX=IX5-5,IX5+5
	    DO IY=IY5-5,IY5+5
	      SUM=SUM+RIMAGE(IX,IY)
	    ENDDO
	  ENDDO
C Deliberately underestimate value by 1 (see intensity > 0 below)
	  IMAGE_ZERO=IFIX( SUM / 121.0 ) - 1
	ENDIF
C
C Subtract zero level from intensities (later added back for *.img file)
C Make all intensities > 0 or they will get confused with "edge" pixels
C
	WRITE(30,'(A,I4,A)') 'Subtracting zero level =',IMAGE_ZERO,' from all intensities'
	DO IY=1,NUMY
	  DO IX=1,NUMX
	    RIMAGE(IX,IY)=MAX(1.0, RIMAGE(IX,IY)-IMAGE_ZERO )	! ensure intensity > 0
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE GET_IMAGE_NUMXY(NUMX,NUMY)
C
C Return the size of the image
C
	COMMON /NUMXY_COM/ NUMX0,NUMY0
C
	NUMX=NUMX0
	NUMY=NUMY0
	RETURN
	END
