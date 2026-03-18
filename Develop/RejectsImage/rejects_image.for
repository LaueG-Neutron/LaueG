	PROGRAM REJECTS_IMAGE
C
	CHARACTER*80 FILE_NAMES(100)
	INTEGER*2 IMAGE2(1000*500)
	REAL XCEN(100),YCEN(100)
C
C Output a simple banner
C
	PRINT '(1X,A)','Make rejects image for LaueG (Ross Piltz, 5/9/2024)'
	PRINT *
C
C Delete the image file, if it exists
C
	OPEN(UNIT=2,FILE='___laueg_rejects_image.img',STATUS='OLD',IOSTAT=IDUMMY)
	CLOSE(2,STATUS='DELETE',IOSTAT=IDUMMY)
C
C Read the input file for its list of files and XY
C
	OPEN(UNIT=1,FILE='___laueg_rejects_image.in',STATUS='OLD',IOSTAT=IERR)
	IF(IERR .NE. 0) CALL QUIT('ERROR: Unable to open input file')
	READ(1,*) NFILES
	IF(NFILES .LT. 1) CALL QUIT('ERROR: Invalid NFILES < 1')
	IF(NFILES .GT. 100) CALL QUIT('ERROR: Invalid NFILES > 100')
	DO I=1,NFILES
	  FILE_NAMES(I)=REPEAT(' ',LEN(FILE_NAMES(I)))
	  READ(1,'(Q,A)') ILEN,FILE_NAMES(I)(1:ILEN)
	  READ(1,*) XCEN(I),YCEN(I)
	ENDDO
	CLOSE(1)
C
C Zero the array IMAGE2()
C
	DO I=1,1000*500
	  IMAGE2(I)=0
	ENDDO
C
C Read the tiles of the images into IMAGE2()
C
	DO I=1,NFILES
	  CALL READ_TILE_IMAGE(IMAGE2,
	1	TRIM(FILE_NAMES(I)),I,XCEN(I),YCEN(I))
	ENDDO
C
C Output the image file
C
	CALL WRITE_I2_FILE('___laueg_rejects_image.img',IMAGE2,1000*500)
C
C Delete the input file
C
	CALL DELETE_FILE('___laueg_rejects_image.in')
C
	PRINT *,'SUCCESSFUL COMPLETION'
	END


	SUBROUTINE READ_TILE_IMAGE(IMAGE2, FILE_NAME,ITILE,XCEN,YCEN)
C
C Copy a 99x99 image centered on XCEN,YCEN read from FILE_NAME onto
C tile number ITILE of IMAGE2().
C Tiles  1 to  9 start at (1,  2) (101,  2) ... (901,  2)
C Tiles 11 to 19 start at (1,102) (101,102) ... (901,102)
C
	CHARACTER FILE_NAME*(*)
	INTEGER*2 IMAGE2(1000,500)
C
	INTEGER*2 IMAGE2_SQUARE(99,99)
C
C Read a 99x99 square of pixels centered on XCEN,YXCEN into IMAGE2_SQUARE()
C
	IXCEN=NINT(XCEN)
	IYCEN=NINT(YCEN)-1		! not sure why, but it works 
	CALL READ_IMAGE_SQUARE(IMAGE2_SQUARE, FILE_NAME,IXCEN,IYCEN,99)
C
C Copy IMAGE2_SQUARE() onto tile ITILE of IMAGE2()
C NB: Invert Y at this stage
C
	IY_TILE=1+(ITILE-1)/10
	IX_TILE=ITILE-10*(IY_TILE-1)
	IX0=(IX_TILE-1)*100
	IY0=(IY_TILE-1)*100
C
	DO IY=1,99
	  DO IX=1,99
	    IMAGE2(IX0+IX,500-IY0-IY)=IMAGE2_SQUARE(IX,IY)
	  ENDDO
	ENDDO
C
	RETURN
	END



	SUBROUTINE READ_IMAGE_SQUARE(IMAGE2, FILE_NAME,IXCEN,IYCEN,IWID)
C
C Read a IWID x IWID square of the image centered on IXCEN,IYCEN from the
C TIF file FILE_NAME into the INTEGER*2 array IMAGE2
C
	CHARACTER FILE_NAME*(*)
	INTEGER*2 IMAGE2(IWID,IWID)
C
	LOGICAL LEXISTS
	INTEGER*2 IMAGE2_READ(8000*2500)
C
C Read the complete image of FILE_NAME into IMAGE2_READ
C
	INQUIRE(FILE=FILE_NAME,EXIST=LEXISTS)
	IF( .NOT.LEXISTS ) CALL QUIT('ERROR: Unable to find file: '//FILE_NAME)
	CALL READ_LAUE_TIFF_QUIET(IMAGE2_READ,NUMX,NUMY, FILE_NAME)
	IF(NUMX .LT. 0) CALL QUIT('ERROR: Unable to read: '//FILE_NAME)
C
C Zero the output image
C
	DO IY=1,IWID
	  DO IX=1,IWID
	    IMAGE2(IX,IY)=0
	  ENDDO
	ENDDO
C
C Copy an IWID * IWID square of IMAGE2_READ into IMAGE2, taking care
C of the indices going past the edges of IMAGE2_READ
C
	IXLO=IXCEN-IWID/2
	IYLO=IYCEN-IWID/2
	DO IY=MAX(1,IYLO),MIN(NUMY,IYLO+IWID-1)
	  DO IX=MAX(1,IXLO),MIN(NUMX,IXLO+IWID-1)
	    IMAGE2(IX-IXLO+1,IY-IYLO+1)=IMAGE2_READ(IX+NUMX*(IY-1))
	  ENDDO
	ENDDO
C
	RETURN
	END



	SUBROUTINE WRITE_I2_FILE(FILE_NAME,I2_OUT,NOUT)
C
	CHARACTER FILE_NAME*(*)
	INTEGER*2 I2_OUT(NOUT)
C
	OPEN(UNIT=2,FILE=FILE_NAME,STATUS='UNKNOWN',FORM='BINARY',
	1					ACCESS='SEQUENTIAL')
	WRITE(2) I2_OUT
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
