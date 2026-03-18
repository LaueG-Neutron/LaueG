	PROGRAM STRIP_BACK
C
	COMMON /EXCLUDE_COM/ NCIRC,ICIRC(3,999),NRECT,IRECT(4,999)
	COMMON /BACK_GEOM_COM/ NUMX,NUMY,IXCEN,IYCEN,XSIZE,YSIZE,DRUM_RAD,IPIXCEN
C
	CHARACTER TIFF_NAME*80
	REAL BACK(8000,2500),RIMAGE(8000,2500)
C
C Output a simple banner
C
	PRINT '(1X,A)','Background stripping program for LaueG (Ross Piltz, 10/1/2025)'
	PRINT *
C
C Delete any output files that still exist
C
	OPEN(UNIT=1,FILE='___laueg_strip_back.raw',STATUS='UNKNOWN',IOSTAT=IDUMMY)
	CLOSE(UNIT=1,STATUS='DELETE',IOSTAT=IDUMMY)
C
C Read the input parameter file
C
	OPEN(UNIT=1,FILE='___laueg_strip_back.in',STATUS='OLD',IOSTAT=IERR)
	IF(IERR .NE. 0) CALL QUIT('ERROR: Unable to open input file')
C Read image file name and image dimensions
	READ(1,'(A)') TIFF_NAME
	READ(1,*) NUMX,NUMY
	IF(NUMX.GT.8000 .OR. NUMY.GT.2500) CALL QUIT('ERROR: Image > 8000 x 2500 pixels')
C Read base counts (intensity in exit hole)
	READ(1,*) IBASE
C Read pixel centre
	READ(1,*) IXCEN,IYCEN
C Read pixel sizes and make sure they are positive
	READ(1,*) XSIZE,YSIZE
	XSIZE=ABS(XSIZE)
	YSIZE=ABS(YSIZE)
C Read drum size in mm
	READ(1,*) DRUM_RAD
C Read pixcen stability: 0=unstable, 1=stable
	READ(1,*) IPIXCEN
C Read image margins, circular and rectangular exclusion regions
	READ(1,*) IXLO,IXHI,IYLO,IYHI
	READ(1,*) NCIRC
	IF(NCIRC .GT. 999) CALL QUIT('ERROR: NCIRC > 999')
	IF(NCIRC .GT. 0) READ(1,*) ((ICIRC(K1,K2),K1=1,3),K2=1,NCIRC)
	READ(1,*) NRECT
	IF(NRECT .GT. 999) CALL QUIT('ERROR: NRECT > 999')
	IF(NRECT .GT. 0) READ(1,*) ((IRECT(K1,K2),K1=1,4),K2=1,NRECT)
C
	CLOSE(1)
C
C Read in the image file
C NB: Base intensity (IBASE) is subtracted from RIMAGE()
C
	CALL READ_IMAGE(RIMAGE, IBASE,NUMX2,NUMY2, TRIM(TIFF_NAME)//'.tif')
	IF(NUMX2.NE.NUMX .OR. NUMY2.NE.NUMY) CALL QUIT('ERROR: Inconsistent image size')
C
C Add the image margins as four exclusion rectangles
C
	CALL ADD_MARGINS(IXLO,IXHI,IYLO,IYHI,NUMX,NUMY)
C
C Extract the background component from RIMAGE() and copy to BACK().
C
	CALL CALC_BACK_ARRAY(BACK, RIMAGE)
C
C Zero RIMAGE() in the exit hole (and entrance hole if it exists)
C NB: This has already happened to the BACK()
C
	CALL ZERO_HOLE_COUNTS(RIMAGE)
C
C Subtract background from RIMAGE(), and add base-counts
C
	DO IY=1,NUMY
	  DO IX=1,NUMX
	    RIMAGE(IX,IY)=RIMAGE(IX,IY)-BACK(IX,IY) +IBASE
	  ENDDO
	ENDDO
C
C Output the stripped image file
C
	CALL WRITE_STRIP_FILE('___laueg_strip_back.raw',RIMAGE,NUMX,NUMY)
C
C Delete the input file
C
	CALL DELETE_FILE('___laueg_strip_back.in')
C
	PRINT *,'SUCCESSFUL COMPLETION'
	END


	SUBROUTINE ADD_MARGINS(IXLO,IXHI,IYLO,IYHI,NUMX,NUMY)
C
C Add the IP margins as four exclusion rectangles
C
	COMMON /EXCLUDE_COM/ NCIRC,ICIRC(3,999),NRECT,IRECT(4,999)
C
	IRECT(1,NRECT+1)=1
	IRECT(2,NRECT+1)=NUMX
	IRECT(3,NRECT+1)=1
	IRECT(4,NRECT+1)=IYLO-1
C
	IRECT(1,NRECT+2)=1
	IRECT(2,NRECT+2)=NUMX
	IRECT(3,NRECT+2)=IYHI+1
	IRECT(4,NRECT+2)=NUMY
C
	IRECT(1,NRECT+3)=1
	IRECT(2,NRECT+3)=IXLO-1
	IRECT(3,NRECT+3)=1
	IRECT(4,NRECT+3)=NUMY
C
	IRECT(1,NRECT+4)=IXHI+1
	IRECT(2,NRECT+4)=NUMX
	IRECT(3,NRECT+4)=1
	IRECT(4,NRECT+4)=NUMY
C
	NRECT=NRECT+4
C
	RETURN
	END


	SUBROUTINE READ_IMAGE(RIMAGE, IBASE_COUNTS,NUMX,NUMY,FILE_NAME)
C
	CHARACTER FILE_NAME*(*)
	REAL RIMAGE(8000,2500)
C
	INTEGER*2 IMAGE(8000*2500)
C
	PRINT '(/,1X,2A)','Reading data file: ',FILE_NAME
C
	CALL CHECK_LAUE_TIFF(NUMX,NUMY, FILE_NAME)
	IF(NUMX*NUMY .GT. 8000*2500) CALL QUIT('ERROR: Image too large')
C
	CALL READ_LAUE_TIFF(IMAGE,NUMX,NUMY, FILE_NAME)
	IF(NUMX .LT. 0) CALL QUIT('ERROR: Unable to read data file')
C
C Copy the unsigned I*16 1D array IMAGE() to the real 2D array RIMAGE().
C Also subtract IBASE_COUNTS from RIMAGE().
C
	PRINT '(1X,A,I4)','Using base counts of',IBASE_COUNTS
	I=0
	DO IY=1,NUMY
	  DO IX=1,NUMX
	    I=I+1
	    RIMAGE(IX,IY)=ZEXT(IMAGE(I)) -IBASE_COUNTS
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE WRITE_STRIP_FILE(FILE_NAME,RIMAGE,NUMX,NUMY)
C
	CHARACTER FILE_NAME*(*)
	REAL RIMAGE(8000,2500),RIMAGE_OUT(8000*2500)
C
	IOUT=0
	DO IY=1,NUMY
	  DO IX=1,NUMX
C
	    IOUT=IOUT+1
	    RIMAGE_OUT(IOUT)=RIMAGE(IX,IY)
C
	  ENDDO
	ENDDO
C
	PRINT *,'Writing background stripped image to ',TRIM(FILE_NAME)
	CALL WRITE_R4_FILE(FILE_NAME,RIMAGE_OUT,NUMX*NUMY)
C
	RETURN
	END



	SUBROUTINE WRITE_R4_FILE(FILE_NAME,R4_OUT,NOUT)
C
	CHARACTER FILE_NAME*(*)
	REAL R4_OUT(NOUT)
C
	OPEN(UNIT=2,FILE=FILE_NAME,STATUS='UNKNOWN',FORM='BINARY',
	1					ACCESS='SEQUENTIAL')
	WRITE(2) R4_OUT
	CLOSE(UNIT=2)
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
