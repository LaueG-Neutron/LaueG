	PROGRAM IMAGE_INFO
C
C
	CHARACTER*2000 COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,ADC_INFO
	COMMON /TIFFINFO_COM/ COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,
	1			IOFF_DATA,NUMX,NUMY,IXYRES,ISTARTX,ISTARTY,ISPEED,
	2			EXPOSE_TIME,EXPOSE_PHI,TEMP_BEGIN,TEMP_END,TEMP_MIN,TEMP_MAX,
	3			IADC,ADC_ADD,ADC_DIV,ADC_INFO,IAPER_1,IAPER_2,IKAPPA
C
	CHARACTER TIFF_NAME*80
C
C Output a simple banner
C
	PRINT '(1X,A)','Image info program for LaueG (Ross Piltz, 11/12/2023)'
	PRINT *
C
C Delete the output files if they still exist
C
	OPEN(UNIT=1,FILE='___laueg_image_info.out',STATUS='OLD',IOSTAT=IDUMMY)
	CLOSE(UNIT=1,STATUS='DELETE',IOSTAT=IDUMMY)
C
C Read the TIF file name and IOPTION from the input file
C If IOPTION=1, we just want the header information
C If IOPTION=2, also add central hole and base intensity from image data
C
	OPEN(UNIT=1,FILE='___laueg_image_info.in',STATUS='OLD',IOSTAT=IERR)
	IF(IERR .NE. 0) CALL QUIT('ERROR: Unable to open input file')
	TIFF_NAME=REPEAT(' ',LEN(TIFF_NAME))
	READ(1,'(Q,A)') ILEN,TIFF_NAME(1:ILEN)
	READ(1,*) IOPTION
	CLOSE(1)
C
C Read the tiff header (IOPTION=1), or the full file (IOPTION=2)
C
	PRINT *,'Reading image file: '//TRIM(TIFF_NAME)//'.tif'
	IF(IOPTION .EQ. 1) THEN
	  CALL CHECK_LAUE_TIFF(NUMX2,NUMY2,TRIM(TIFF_NAME)//'.tif')
	  IF(NUMX2 .LT. 0) THEN
	    IF(NUMY2 .LT. 0) CALL QUIT('ERROR: Unable to read image file')
	    CALL QUIT('ERROR: Image file header is invalid')
	  ENDIF
	ELSEIF(IOPTION .EQ. 2) THEN
	  CALL READ_IMAGE(TRIM(TIFF_NAME)//'.tif')
	ELSE
	  CALL QUIT('ERROR: Invalid IOPTION in input file')
	ENDIF
C
C Set some default parameters
C
	IBASE=0
	IXHOLE=NUMX/2
	IYHOLE=NUMY/2
	IRAD=10800/IXYRES
	IF(TRIM(HOST) .EQ. 'CYCLOPS_LAUEG') THEN
	  IXHOLE=2160
	  IYHOLE=600
	  IRAD=10000/IXYRES
	ENDIF
C Up the compress ratio for larger images
	ICOMPRESS=4
	IF(NUMX .GT. 4000) ICOMPRESS=1 + (NUMX-1)/1000
	IF(TRIM(HOST) .EQ. 'KOALA2') ICOMPRESS=5
C
C If IOPTION=2, use image data to get base intensity (and possibly central hole)
C
	IF(IOPTION .EQ. 2) THEN
C
C For both KOALA & VIVALDI, find the central hole as it moves about
	  IF( (TRIM(HOST) .EQ. 'KOALA') .OR. (TRIM(HOST) .EQ. 'VIVALDI') ) THEN
C Fit the position of the central hole
          CALL FIT_HOLE_APPROX(IXHOLE,IYHOLE,IRAD, NUMX/8,NUMY/8)
          DO ITER=1,2
	      CALL FIT_HOLE_XCEN(IXHOLE,IYHOLE,IRAD,10)
	      CALL FIT_HOLE_YCEN(IXHOLE,IYHOLE,10,IRAD)
	    ENDDO
	  ENDIF
C
C Get the "base intensity" from intensity in the exit hole
	  CALL GET_BASE_COUNTS(IBASE, IXHOLE,IYHOLE,IRAD)
C
	ENDIF
C
C Write to output file
C
	OPEN(UNIT=2,FILE='___laueg_image_info.out',STATUS='NEW')
	WRITE(2,'(A)'   ) TRIM(COMMENT)
	WRITE(2,'(A)'   ) TRIM(SAMPLE)
	WRITE(2,'(A)'   ) TRIM(USER)
	WRITE(2,'(A)'   ) TRIM(HOST)
	WRITE(2,'(A)'   ) TRIM(DATETIME)
	WRITE(2,'(F8.3)') MAX(-360.0,MIN(+360.0, EXPOSE_PHI ))
	WRITE(2,'(I8)'  ) NINT(EXPOSE_TIME)
	WRITE(2,'(2I8)' ) NUMX,NUMY
	WRITE(2,'(I8)'  ) IOFF_DATA
	WRITE(2,'(2I8)' ) ISTARTX,ISTARTY
	WRITE(2,'(I8)'  ) IXYRES
	WRITE(2,'(3I8)' ) IXHOLE,IYHOLE,IRAD
	WRITE(2,'(I8)'  ) IBASE
	WRITE(2,'(I8)'  ) ICOMPRESS
C
C Additional parameters for KOALA2
C
	IF(TRIM(HOST) .EQ. 'KOALA2') WRITE(2,'(4I8)') IADC,IAPER_1,IAPER_2,IKAPPA
C
	CLOSE(UNIT=2)
C
C Delete the input file (optionally)
C
	CALL DELETE_FILE('___laueg_image_info.in')
C
	PRINT *,'SUCCESSFUL COMPLETION'
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



      SUBROUTINE FIT_HOLE_APPROX(IXCEN,IYCEN,IRCEN, IDX,IDY)
C
C Find approximate hole center in image at IXCEN+/-IDX, IYCEN+/-IDY using
C a very approximate radius of IRCEN.
C New hole center is returned in IXCEN,IYCEN.
C
      BEST=1E9
      DO IX=IXCEN-IDX,IXCEN+IDX
	  DO IY=IYCEN-IDY,IYCEN+IDY
          VAL=COUNTS(IX,IY)*(IRCEN+ABS(IX-IXCEN)+ABS(IY-IYCEN))
          IF(VAL .LT. BEST) THEN
            BEST=VAL
            IXCEN2=IX
            IYCEN2=IY
          ENDIF
        ENDDO
      ENDDO
      IF(BEST .GT. 0.9E9) STOP 'BUG: fit_hole_approx #1'
C
      IXCEN=IXCEN2
      IYCEN=IYCEN2
C
      RETURN
      END



	SUBROUTINE FIT_HOLE_XCEN(IXCEN,IYCEN,IDX,IDY)
C
C
	REAL XCOUNTS(10000),SMOOTH(10000)
C
C Across the detector (in X), take the minimum counts in two
C strips IDY high (in Y) just above and below the center IYCEN.
C Then take the maximum of both strips for each value of X and
C store in XCOUNTS.
C
	DO IX=IXCEN-2*IDX-2,IXCEN+2*IDX
C
	  XCOUNT1=1E6
	  DO IY=IYCEN-IDY,IYCEN-1
	    XCOUNT1=MIN(XCOUNT1,COUNTS(IX,IY))
	  ENDDO
C
	  XCOUNT2=1E6
	  DO IY=IYCEN+1,IYCEN+IDY
	    XCOUNT2=MIN(XCOUNT2,COUNTS(IX,IY))
	  ENDDO
C
	  XCOUNTS(IX)=MAX(XCOUNT1,XCOUNT2)
	ENDDO
C
C Smooth XCOUNTS over 2*IDX+1 points
C
	SMOOTH(IXCEN-2*IDX-2)=XCOUNTS(IXCEN-2*IDX-2)
	DO I=IXCEN-2*IDX-1,IXCEN+2*IDX
	  SMOOTH(I)=SMOOTH(I-1)+XCOUNTS(I)
	ENDDO
C
	DO I=IXCEN-IDX,IXCEN+IDX
	  XCOUNTS(I)=(SMOOTH(I+IDX)-SMOOTH(I-IDX-1))/(2*IDX+1)
	ENDDO
C
C Search for minimum in the smoothed XCOUNTS and return as IXCEN
C
	IMIN=IXCEN
	DO I=IXCEN-IDX,IXCEN+IDX
	  IF(XCOUNTS(I) .LT. XCOUNTS(IMIN)) IMIN=I
	ENDDO
	IXCEN=IMIN
C
	RETURN
	END



	SUBROUTINE FIT_HOLE_YCEN(IXCEN,IYCEN,IDX,IDY)
C
C
	REAL YCOUNTS(10000),SMOOTH(10000)
C
C Across the detector (in Y), take the minimum counts in two
C strips IDX high (in X) just above and below the center IXCEN.
C Then take the maximum of both strips for each value of Y and
C store in YCOUNTS.
C
	DO IY=IYCEN-2*IDY-2,IYCEN+2*IDY
C
	  YCOUNT1=1E6
	  DO IX=IXCEN-IDX,IXCEN-1
	    YCOUNT1=MIN(YCOUNT1,COUNTS(IX,IY))
	  ENDDO
C
	  YCOUNT2=1E6
	  DO IX=IXCEN+1,IXCEN+IDX
	    YCOUNT2=MIN(YCOUNT2,COUNTS(IX,IY))
	  ENDDO
C
	  YCOUNTS(IY)=MAX(YCOUNT1,YCOUNT2)
	ENDDO
C
C Smooth YCOUNTS over 2*IDY+1 points
C
	SMOOTH(IYCEN-2*IDY-2)=YCOUNTS(IYCEN-2*IDY-2)
	DO I=IYCEN-2*IDY-1,IYCEN+2*IDY
	  SMOOTH(I)=SMOOTH(I-1)+YCOUNTS(I)
	ENDDO
C
	DO I=IYCEN-IDY,IYCEN+IDY
	  YCOUNTS(I)=(SMOOTH(I+IDY)-SMOOTH(I-IDY-1))/(2*IDY+1)
	ENDDO
C
C Search for minimum in the smoothed YCOUNTS and return as IYCEN
C
	IMIN=IYCEN
	DO I=IYCEN-IDY,IYCEN+IDY
	  IF(YCOUNTS(I) .LT. YCOUNTS(IMIN)) IMIN=I
	ENDDO
	IYCEN=IMIN
C
	RETURN
	END



	SUBROUTINE GET_BASE_COUNTS(ICOUNTS, IXCEN,IYCEN,IRAD0)
C
C
	IRAD=IRAD0/4
C
	NSUM=0
	SUM=0
	DO IX=IXCEN-IRAD,IXCEN+IRAD
	  DO IY=IYCEN-IRAD,IYCEN+IRAD
	    IF((IXCEN-IX)**2 + (IYCEN-IY)**2 .LE. IRAD**2) THEN
	      NSUM=NSUM+1
		  SUM=SUM+COUNTS(IX,IY)
	    ENDIF
	  ENDDO
	ENDDO
C
	ICOUNTS=NINT(SUM/MAX(1,NSUM))
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

