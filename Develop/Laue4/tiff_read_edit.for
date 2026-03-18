C KNOWN ISSUES:
C Assumes the TIF file is a 16bit contiguous image created on VIVALDI, etc.
C Unchecked for big-endian computers (especially double & short strings)
C Currently set for the Intel compiler with default options 
C
C Reading and editting TIF files should use the following routines:
C	SUBROUTINE WRITE_LAUE_TIFF(FILE_NAME, IMAGE2,NUMX,NUMY,
C	1			HOST,USER,SAMPLE,DATETIME,COMMENTS,PHI)
C	SUBROUTINE EDIT_LAUE_TIFF(IMAGE, NUMX,NUMY, FILE_NAME_WRITE,FILE_NAME_READ)
C	SUBROUTINE CHECK_LAUE_TIFF(NUMX,NUMY, FILE_NAME)
C	SUBROUTINE READ_LAUE_TIFF(IMAGE,NUMX,NUMY, FILE_NAME)
C	SUBROUTINE READ_LAUE_TIFF_QUIET(IMAGE,NUMX,NUMY, FILE_NAME)
C
C After calling any of the above, except EDIT_LAUE_TIFF, the header information
C of the TIF file is accessible via the common area:
C	CHARACTER*2000 COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,ADC_INFO
C	COMMON /TIFFINFO_COM/ COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,
C	1			IOFF_DATA,NUMX,NUMY,IXYRES,ISTARTX,ISTARTY,ISPEED,
C	2			EXPOSE_TIME,EXPOSE_PHI,TEMP_BEGIN,TEMP_END,TEMP_MIN,TEMP_MAX,
C	3			IADC,ADC_ADD,ADC_DIV,ADC_INFO,IAPER_1,IAPER_2,IKAPPA
C or by calling
C	SUBROUTINE GET_TIFF_INFO(STRING,VALUE,IVALUE, NAME)

C NB: Needs a QUIT() routine. See bottom of file for example.


C======================================================================
C  User level routines to read and edit Laue TIF files
C======================================================================

	SUBROUTINE WRITE_LAUE_TIFF(FILE_NAME, IMAGE2,NUMX,NUMY,
	1			HOST,USER,SAMPLE,DATETIME,COMMENTS,PHI)
C
C Create a Laue tif file which is non-standard as it uses
C user defined TAGs 1002 & 1004 for the sample & PHI value
C
	CHARACTER*(*) FILE_NAME,HOST,COMMENTS,DATETIME,USER,SAMPLE
	INTEGER*2 IMAGE2(NUMX*NUMY)
C
	OPEN(UNIT=2,FILE=FILE_NAME,STATUS='UNKNOWN',
	1			FORM='UNFORMATTED',ACCESS='STREAM')
C
C Write the magic numbers
C
	CALL TIFF_WRITE_I2(257*ICHAR('I'))
	CALL TIFF_WRITE_I2(42)
C
C Write the offset for the trailing IFD
C
	IOFFSET_DATA=8
	NBYTES=2*NUMX*NUMY
	IOFFSET_IFD=IOFFSET_DATA+NBYTES
	CALL TIFF_WRITE_I4(IOFFSET_IFD)
C
C Write the image as one block
C
	WRITE(2) IMAGE2
C
C Write TAGS to the trailing IFD
C
	NTAGS=14
	CALL TIFF_WRITE_I2(NTAGS)
C Write the I2 tags containing values
	NBITS=16
	IPHOTO=1
	IUP=1
	NROWS=NUMY
	CALL TIFF_WRITE_I2_TAG(256,NUMX)
	CALL TIFF_WRITE_I2_TAG(257,NUMY)
	CALL TIFF_WRITE_I2_TAG(258,NBITS)
	CALL TIFF_WRITE_I2_TAG(259,IPHOTO)
	CALL TIFF_WRITE_I2_TAG(262,IUP)
	CALL TIFF_WRITE_I2_TAG(278,NROWS)
C Write the I4 tags containing values (or offset to image data)
	CALL TIFF_WRITE_I4_TAG(279,NBYTES)
	CALL TIFF_WRITE_I4_TAG(273,IOFFSET_DATA)
C Write the tags with offsets to values after the IFD
	IOFFSET=IOFFSET_IFD+2+12*NTAGS
	CALL TIFF_WRITE_R8_TAG(1004,IOFFSET)
	IOFFSET=IOFFSET+8
	CALL TIFF_WRITE_STRING_TAG(316,IOFFSET,HOST)
	IOFFSET=IOFFSET+LEN_TRIM(HOST)+1
	CALL TIFF_WRITE_STRING_TAG(270,IOFFSET,COMMENTS)
	IOFFSET=IOFFSET+LEN_TRIM(COMMENTS)+1
	CALL TIFF_WRITE_STRING_TAG(306,IOFFSET_DATETIME,DATETIME)
	IOFFSET=IOFFSET+LEN_TRIM(DATETIME)+1
	CALL TIFF_WRITE_STRING_TAG(315,IOFFSET,USER)
	IOFFSET=IOFFSET+LEN_TRIM(USER)+1
	CALL TIFF_WRITE_STRING_TAG(1002,IOFFSET,SAMPLE)
C Write the values after the IFD
	CALL TIFF_WRITE_R8(PHI)
	CALL TIFF_WRITE_STRING(HOST)
	CALL TIFF_WRITE_STRING(COMMENTS)
	CALL TIFF_WRITE_STRING(DATETIME)
	CALL TIFF_WRITE_STRING(USER)
	CALL TIFF_WRITE_STRING(SAMPLE)
C
	CLOSE(UNIT=2)
	RETURN
	END



	SUBROUTINE EDIT_TIFF_IMAGE(FILE_NAME_WRITE, IMAGE2,FILE_NAME_READ)
C
C Read the TIF file with the name FILE_NAME_READ and copy all header information
C to a new TIF file with the name FILE_NAME_WRITE, but substitute the image data
C with the data in IMAGE(). The NUMX,NUMY of the original file must be given as
C input to the routine, and it is assumed that it matches the input file, but if
C it doesn't, no warning is given and the output file will probably be useless.
C
C NB: The routine overwrites IMAGE()
C NB: The routine uses units 1 & 2 to read and write the files
C
	CHARACTER*(*) FILE_NAME_WRITE,FILE_NAME_READ
	INTEGER*2 IMAGE2(*)
C
C Check the input file and read the header into common
C
	CALL CHECK_LAUE_TIFF(NUMX,NUMY, FILE_NAME_READ)
	IF(NUMX .LT. 0) CALL QUIT('ERROR: Input *.tif file is missing or invalid')
C
C Open input and output files to read as binary streams
C
	OPEN(UNIT=1,FILE=FILE_NAME_READ,STATUS='OLD',
	1				FORM='UNFORMATTED',ACCESS='STREAM')
	OPEN(UNIT=2,FILE=FILE_NAME_WRITE,STATUS='UNKNOWN',
	1				FORM='UNFORMATTED',ACCESS='STREAM')
C
C Use helper routine to read/write so we can do IMAGE in one block
C
	CALL EDIT_TIFF_IMAGE_IO(IMAGE2,NUMX*NUMY)
C
	CLOSE(UNIT=1)
	CLOSE(UNIT=2)
C
	RETURN
	END


	SUBROUTINE CHECK_LAUE_TIFF(NUMX,NUMY, FILE_NAME)
C
C Reads and checks the header of the Laue TIF file FILE_NAME and returns the size
C of the image in NUMX NUMY. On an error, NUMX=-1 and ISTATUS is returned in NUMY.
C NB: The routine uses UNIT=1 to read the file.
C
	CHARACTER FILE_NAME*(*)
C
	INTEGER*2 IMAGE(1)
C
	IUNIT=1
	CALL READ_TIFF_IMAGE(IMAGE,NUMX,NUMY,ISTATUS, FILE_NAME,IUNIT,.FALSE.,.FALSE.)
C
	IF(ISTATUS .NE. 0) THEN
	  NUMX=-1
	  NUMY=ISTATUS
	ENDIF
C
	RETURN
	END


	SUBROUTINE READ_LAUE_TIFF(IMAGE,NUMX,NUMY, FILE_NAME)
C
C Reads the Laue TIF file FILE_NAME into the 2byte array IMAGE2 and returns
C the size of the image in NUMX NUMY. On an error, NUMX=NUMY=-1 is returned.
C This routine outputs information on the experiment and IO errors.
C NB: The routine uses UNIT=1 to read the file.
C
	CHARACTER FILE_NAME*(*)
	INTEGER*2 IMAGE(*)
C
	IUNIT=1
	CALL READ_TIFF_IMAGE(IMAGE,NUMX,NUMY,ISTATUS, FILE_NAME,IUNIT,.TRUE.,.TRUE.)
C
	IF(ISTATUS .NE. 0) THEN
	  NUMX=-1
	  NUMY=-1
	ENDIF
C
	RETURN
	END


	SUBROUTINE READ_LAUE_TIFF_QUIET(IMAGE,NUMX,NUMY, FILE_NAME)
C
C Reads the Laue TIF file FILE_NAME into the 2byte array IMAGE2 and return
C the size of the image in NUMX NUMY. On an error, NUMX=NUMY=-1 is returned.
C This routine outputs NO information on the experiment or IO errors.
C NB: The routine uses UNIT=1 to read the file.
C
	CHARACTER FILE_NAME*(*)
	INTEGER*2 IMAGE(*)
C
	IUNIT=1
	CALL READ_TIFF_IMAGE(IMAGE,NUMX,NUMY,ISTATUS, FILE_NAME,IUNIT,.TRUE.,.FALSE.)
C
	IF(ISTATUS .NE. 0) THEN
	  NUMX=-1
	  NUMY=-1
	ENDIF
C
	RETURN
	END


	SUBROUTINE GET_TIFF_INFO(SVALUE,RVALUE,IVALUE, NAME)
C
C Utility to return various values read from the header
C
	CHARACTER SVALUE*(*),NAME*(*)
C
	CHARACTER*2000 COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,ADC_INFO
	COMMON /TIFFINFO_COM/ COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,
	1			IOFF_DATA,NUMX,NUMY,IXYRES,ISTARTX,ISTARTY,ISPEED,
	2			EXPOSE_TIME,EXPOSE_PHI,TEMP_BEGIN,TEMP_END,TEMP_MIN,TEMP_MAX,
	3			IADC,ADC_ADD,ADC_DIV,ADC_INFO,IAPER_1,IAPER_2,IKAPPA
C
C Copy string parameters
C
	IF(NAME .EQ. 'COMMENT') THEN
	  SVALUE=COMMENT
	ELSEIF(NAME .EQ. 'SAMPLE') THEN
	  SVALUE=SAMPLE
	ELSEIF(NAME .EQ. 'USER') THEN
	  SVALUE=USER
	ELSEIF(NAME .EQ. 'DATETIME') THEN
	  SVALUE=DATETIME
	ELSEIF(NAME .EQ. 'HOST') THEN
	  SVALUE=HOST
	ELSEIF(NAME .EQ. 'BEAM') THEN
	  SVALUE=BEAM
C
C Copy integer parameters
C
	ELSEIF(NAME .EQ. 'NUMX') THEN
	  IVALUE=NUMX
	ELSEIF(NAME .EQ. 'NUMY') THEN
	  IVALUE=NUMY
	ELSEIF(NAME .EQ. 'IXYRES') THEN
	  IVALUE=IXYRES
	ELSEIF(NAME .EQ. 'ISTARTX') THEN
	  IVALUE=ISTARTX
	ELSEIF(NAME .EQ. 'ISTARTY') THEN
	  IVALUE=ISTARTY
	ELSEIF(NAME .EQ. 'IOFFSET') THEN
	  IVALUE=IOFF_DATA
C
C Copy real parameters
C
	ELSEIF(NAME .EQ. 'EXPOSE_TIME') THEN
	  RVALUE=EXPOSE_TIME
	ELSEIF(NAME .EQ. 'EXPOSE_PHI') THEN
	  RVALUE=EXPOSE_PHI
	ELSEIF(NAME .EQ. 'TEMP_BEGIN') THEN
	  RVALUE=TEMP_BEGIN
	ELSEIF(NAME .EQ. 'TEMP_END') THEN
	  RVALUE=TEMP_END
	ELSEIF(NAME .EQ. 'TEMP_MIN') THEN
	  RVALUE=TEMP_MIN
	ELSEIF(NAME .EQ. 'TEMP_MAX') THEN
	  RVALUE=TEMP_MAX
C
C Copy new Koala 2 parameters
C
	ELSEIF(NAME .EQ. 'IADC') THEN
	  IVALUE=IADC
	ELSEIF(NAME .EQ. 'ADC_ADD') THEN
	  RVALUE=ADC_ADD
	ELSEIF(NAME .EQ. 'ADC_DIV') THEN
	  RVALUE=ADC_DIV
	ELSEIF(NAME .EQ. 'ADC_INFO') THEN
	  SVALUE=ADC_INFO
	ELSEIF(NAME .EQ. 'IAPER_1') THEN
	  IVALUE=IAPER_1
	ELSEIF(NAME .EQ. 'IAPER_2') THEN
	  IVALUE=IAPER_2
	ELSEIF(NAME .EQ. 'IKAPPA') THEN
	  IVALUE=IKAPPA
C
C None of the above, so complain and die
C
	ELSE
	  PRINT *,'NAME=',TRIM(NAME)
	  CALL QUIT('BUG(get_tiff_info ) Unknown NAME='//TRIM(NAME))
	ENDIF
C
	RETURN
	END


C======================================================================
C  Medium level routines to open the file and read the info and image
C======================================================================

	SUBROUTINE EDIT_TIFF_IMAGE_IO(IMAGE,ISIZE)
C
	INTEGER*2 IMAGE(ISIZE)
C
	CHARACTER*2000 COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,ADC_INFO
	COMMON /TIFFINFO_COM/ COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,
	1			IOFF_DATA,NUMX,NUMY,IXYRES,ISTARTX,ISTARTY,ISPEED,
	2			EXPOSE_TIME,EXPOSE_PHI,TEMP_BEGIN,TEMP_END,TEMP_MIN,TEMP_MAX,
	3			IADC,ADC_ADD,ADC_DIV,ADC_INFO,IAPER_1,IAPER_2,IKAPPA
C
	INTEGER*1 I1_BUFF(4000)
C
C Copy the IOFFSET bytes at the start of the file
C
	IF(IOFF_DATA .GT. 4000) CALL QUIT('BUG(edit_laue_tiff): IOFF_DATA > 4000')
	READ (1) (I1_BUFF(K),K=1,IOFF_DATA)
	WRITE(2) (I1_BUFF(K),K=1,IOFF_DATA)
	NBYTES=IOFF_DATA
C
C Write IMAGE to the output file
C
	WRITE(2) IMAGE
	NBYTES=NBYTES+2*NUMX*NUMY
C
C Read the same amount from the input file to skip to the IFD section
C NB: This overwrites IMAGE() but it works really fast
C
	READ(1) IMAGE
C
C Copy any trailing data byte by byte until we hit an EOF
C
	DO I=1,100000
	  READ (1,END=100) I1_BUFF(1)
	  WRITE(2) I1_BUFF(1)
	  NBYTES=NBYTES+1
	ENDDO
C
c100	PRINT '(1X,A,I9)','Total number of bytes copied =',NBYTES
100	RETURN
	END


	SUBROUTINE READ_TIFF_IMAGE(IMAGE,NUMX0,NUMY0,ISTATUS,
	1				FILE_NAME,IUNIT,LREAD,LVERBOSE)
C
C Read the 16bit image of the Laue TIF file FILE_NAME using unit IUNIT
C Return ISTATUS=0 for success, -1 for file not found, +1 otherwise
C The header information is accessible via /TIFFINFO_COM/
C If LREAD is false the header is checked but the image is not read
C If LVERBOSE is true the header information is printed on the console
C The file unit is closed before return from this routine
C
	LOGICAL LREAD,LVERBOSE
	CHARACTER FILE_NAME*(*)
	INTEGER*2 IMAGE(*)
C
	CHARACTER*2000 COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,ADC_INFO
	COMMON /TIFFINFO_COM/ COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,
	1			IOFF_DATA,NUMX,NUMY,IXYRES,ISTARTX,ISTARTY,ISPEED,
	2			EXPOSE_TIME,EXPOSE_PHI,TEMP_BEGIN,TEMP_END,TEMP_MIN,TEMP_MAX,
	3			IADC,ADC_ADD,ADC_DIV,ADC_INFO,IAPER_1,IAPER_2,IKAPPA
C
	CHARACTER*10 SADC,SKAPPA
C
C Open the Laue TIF file, check it, copy header info to TIFFINFO_COM,
C then close the unit and return with the status
C
	CALL READ_TIFF_INFO(ISTATUS, FILE_NAME,IUNIT)
	NUMX0=NUMX
	NUMY0=NUMY
C
C On error, output info if LVERBOSE is set, then return
C
	IF(LVERBOSE) THEN
	  IF(ISTATUS .EQ. -1) THEN
			PRINT *,'ERROR: Laue TIF file not found'
	  ELSEIF(ISTATUS .GE. 100) THEN
			PRINT *,'WARNING: Laue TIF file may be invalid'
	  ELSEIF(ISTATUS .NE. 0) THEN
			PRINT *,'ERROR: Laue TIF file is invalid'
		ENDIF
	ENDIF
	IF( (ISTATUS .NE. 0) .AND. (ISTATUS .LT. 100) ) RETURN
C
C If LVERBOSE is true, print out experimental information
C
	IF(LVERBOSE) THEN
	  IF(TRIM(COMMENT) .NE. '') PRINT *,TRIM(COMMENT)
C
C Do 'ancient' CYCLOPS instrument
	  IF(TRIM(HOST) .EQ. 'CYCLOPS') THEN
	    PRINT '(1X,A,2X,I5,A,I5)','CYCLOPS',NUMX,' x',NUMY
	    PRINT '(1X,A,F6.1)','PHI = ',EXPOSE_PHI
C
C Do all other LADI-related instruments
	  ELSE
	    PRINT '(1X,7A)',TRIM(SAMPLE),', ',TRIM(USER),', ',TRIM(HOST),
	1				', ',TRIM(DATETIME)
	    PRINT '(1X,A,I6,A,F6.1)','Exposure =',NINT(EXPOSE_TIME),
	1				' s    PHI =',EXPOSE_PHI
	    PRINT '(1X,A,4F8.3,A)','Temp (K) =',TEMP_BEGIN,TEMP_END,
	1				TEMP_MIN,TEMP_MAX,' (start,end,min,max)'
	    PRINT '(1X,2(I5,A),I4,A,2I4)',NUMX,' x',NUMY,' at',IXYRES,
	1				' microns    Start coords =',ISTARTX,ISTARTY
	    PRINT '(1X,2A)','Neutron beam = ',TRIM(BEAM)
	  ENDIF
C
C Extra info for KOALA2
	  IF(TRIM(HOST) .EQ. 'KOALA2') THEN
	    SADC='Gain = 0?'
	    IF(IADC .EQ. 1) SADC='Low Gain'
	    IF(IADC .EQ. 2) SADC='High Gain'
	    SKAPPA='OFF'
	    IF(IKAPPA .EQ. 1) SKAPPA='LOW'
	    IF(IKAPPA .EQ. 2) SKAPPA='HIGH'
	    PRINT '(1X,2A,2I3,2A)',TRIM(SADC),' File, Guide Apertures =',
	1                            IAPER_1,IAPER_2,', KAPPA = ',TRIM(SKAPPA)
	  ENDIF
C
	ENDIF
C
C Do a fast read of the image data into IMAGE (if LREAD is TRUE)
C
	IF(LREAD) CALL FAST_TIFF_READ(IMAGE, NUMX,NUMY,IOFF_DATA, FILE_NAME)
C
	RETURN
	END


	SUBROUTINE FAST_TIFF_READ(IMAGE, NUMX,NUMY,IOFF_DATA, FILE_NAME)
C
C Simple routine to read IMAGE from a Laue TIFF as a single block
C NB: No error checking is done, it is assumed the file is valid
C
	CHARACTER FILE_NAME*(*)
	INTEGER*2 IMAGE(NUMX*NUMY)
C
	INTEGER*1 BUFFER1(4000)
	COMMON /READ_TIF_COM/ BUFFER1,IREC,IUNIT,ISWAP
C
	INTEGER*1 DUMMY(1000)
C
C Open file as a binary stream
C
	OPEN(UNIT=IUNIT,FILE=FILE_NAME,STATUS='OLD',
	1			FORM='UNFORMATTED',ACCESS='STREAM')
C
C Skip the header bytes before the data block
C
	READ(IUNIT) DUMMY(1:IOFF_DATA)
C
C Read the image as a single block
C
	READ(IUNIT) IMAGE
C
C If required, swap the byte ordering
C
	IF(ISWAP .NE. 0) THEN
		DO I=1,NUMX*NUMY
			IMAGE(I)=ISHFTC(IMAGE(I),8)
		ENDDO
	ENDIF
C
C All finished, so close file and return
C
	CLOSE(IUNIT)
C
	RETURN
	END


	SUBROUTINE READ_TIFF_INFO(ISTATUS, FILE_NAME,IUNIT0)
C
C Read the header part of the Laue TIF file FILE_NAME using unit IUNIT
C Status is returned in ISTATUS
C	ISTATUS= 0	Success
C	ISTATUS=-1	Can't open file
C	ISTATUS= 1	Not a TIF file (Invalid magic numbers)
C	ISTATUS= 2	Invalid IFD offset
C	ISTATUS= 3	Invalid IFD NTAGS
C	ISTATUS=11	Invalid Laue image (apparently an ImageJ copy)
C	ISTATUS=12	Not a laue image (big-endian TIF)
C	ISTATUS=100+n	Unknown user defined tags found (n = how many)
C The header information is accessible via /TIFFINFO_COM/
C The file unit is closed before return from this routine
C
	CHARACTER*(*) FILE_NAME
C
	CHARACTER*2000 COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,ADC_INFO
	COMMON /TIFFINFO_COM/ COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,
	1			IOFF_DATA,NUMX,NUMY,IXYRES,ISTARTX,ISTARTY,ISPEED,
	2			EXPOSE_TIME,EXPOSE_PHI,TEMP_BEGIN,TEMP_END,TEMP_MIN,TEMP_MAX,
	3			IADC,ADC_ADD,ADC_DIV,ADC_INFO,IAPER_1,IAPER_2,IKAPPA
C
	INTEGER*1 BUFFER1(4000)
	COMMON /READ_TIF_COM/ BUFFER1,IREC,IUNIT,ISWAP
C
	CHARACTER STRING*2000
C
C Set strings to empty (else they are filled with NUL)
C
	COMMENT=''
	SAMPLE=''
	USER=''
	DATETIME=''
	HOST=''
	BEAM=''
	ADC_INFO=''
C
C Clear numeric parameters
C
	NUMX=0
	NUMY=0
	IXYRES=0
	ISTARTX=0
	ISTARTY=0
	EXPOSE_TIME=0.0
	EXPOSE_PHI=0.0
	TEMP_BEGIN=0.0
	TEMP_END=0.0
	TEMP_MIN=0.0
	TEMP_MAX=0.0
C
C Initialise parameters used by TIFF_READ_INT1() which actually reads the file
C
	IREC=-1
	IUNIT=IUNIT0
C
C Open file by direct access and 4000 byte records (1000 words for Intel)
C
	OPEN(UNIT=IUNIT,FILE=FILE_NAME,STATUS='OLD',
	1       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1000,ERR=900)
C
C Check if a little-endian TIF
C
	CALL TIFF_READ_INT2(IMAGIC,0)
	ISWAP=0
	IF(IMAGIC .EQ. 257*ICHAR('M')) THEN
	  ISWAP=1
	ELSEIF(IMAGIC .NE. 257*ICHAR('I')) THEN
	  GOTO 901
	ENDIF
C
	CALL TIFF_READ_INT2(I42,2)
	IF(I42 .NE. 42) GOTO 901
C
C Get file offset to the first IFD
C
	CALL TIFF_READ_INT4(IOFF_IFD,4)
	IF(IOFF_IFD .LE. 0) GOTO 902
C
C Get number of tags in IFD
C
	CALL TIFF_READ_INT2(NTAGS,IOFF_IFD)
	IF(NTAGS .LE. 0) GOTO 903
	IOFF_IFD=IOFF_IFD+2
C
C Loop through all records in the IFD
C
	NEXTRAS=0
	DO ITAG=1,NTAGS
C
C Read and unpack the next record in the IFD
	  CALL TIFF_READ_INT2(ITAG_ID,IOFF_IFD)
	  CALL TIFF_READ_INT2(ITAG_TYPE,IOFF_IFD+2)
	  CALL TIFF_READ_INT4(ITAG_COUNT,IOFF_IFD+4)
	  CALL TIFF_READ_INT2(ITAG_VAL2,IOFF_IFD+8)
	  CALL TIFF_READ_INT4(ITAG_VAL4,IOFF_IFD+8)
	  ITAG_OFFSET=ITAG_VAL4		! the last I4 may also be an offset
C Create the value (STRING, VALUE) depending on ITAG_TYPE
	  IF(ITAG_TYPE .EQ. 2) THEN
	    CALL TIFF_GET_STRING(STRING,ITAG_OFFSET,ITAG_COUNT)
	    READ(STRING,'(F20.0)',IOSTAT=IDUMMY) VALUE
	  ELSE IF(ITAG_TYPE .EQ. 3) THEN
	    VALUE=ITAG_VAL2
	  ELSE IF(ITAG_TYPE .EQ. 4) THEN
	    VALUE=ITAG_VAL4
	  ELSE IF(ITAG_TYPE .EQ. 5) THEN
	    CALL TIFF_GET_RATIONAL(VALUE,ITAG_OFFSET)
	  ELSE IF(ITAG_TYPE .EQ. 12) THEN
	    CALL TIFF_GET_DOUBLE(VALUE,ITAG_OFFSET)
	  ENDIF
C
C Save some of the general TIFF tags
	  IF     (ITAG_ID .EQ. 256) THEN
	    NUMX=NINT(VALUE)
	  ELSE IF(ITAG_ID .EQ. 257) THEN
	    NUMY=NINT(VALUE)
	  ELSE IF(ITAG_ID .EQ. 273) THEN
C If count=1, the data offset is IVAL4, else stored as an array and have to read value
	    IF(ITAG_COUNT .EQ. 1) THEN
	      IOFF_DATA=ITAG_OFFSET
	    ELSE
	      CALL TIFF_READ_INT4(IOFF_DATA,ITAG_OFFSET)
	    ENDIF
	  ELSE IF(ITAG_ID .EQ. 270) THEN
	    COMMENT=TRIM(STRING)
	  ELSE IF(ITAG_ID.EQ.282) THEN
	    IXYRES=NINT(1E4*VALUE)
	  ELSE IF(ITAG_ID.EQ.286) THEN
	    ISTARTX=NINT(VALUE)
	  ELSE IF(ITAG_ID.EQ.287) THEN
	    ISTARTY=NINT(VALUE)
	  ELSE IF(ITAG_ID .EQ. 306) THEN
C Convert date to a 24 hour Australian format
	    CALL CONVERT_DATE_TIME(STRING,ITAG_COUNT)
	    DATETIME=TRIM(STRING)
	  ELSE IF(ITAG_ID .EQ. 315) THEN
	    USER=TRIM(STRING)
	  ELSE IF(ITAG_ID .EQ. 316) THEN
	    HOST=TRIM(STRING)
C
C Save known user defined tags
	  ELSE IF(ITAG_ID .EQ. 1000) THEN
	    TEMP_BEGIN=VALUE
	  ELSE IF(ITAG_ID .EQ. 1001) THEN
	    TEMP_END=VALUE
	  ELSE IF(ITAG_ID .EQ. 1002) THEN
	    SAMPLE=TRIM(STRING)
	  ELSE IF(ITAG_ID .EQ. 1003) THEN
	    EXPOSE_TIME=VALUE
	  ELSE IF(ITAG_ID .EQ. 1004) THEN
	    EXPOSE_PHI=VALUE
	  ELSE IF(ITAG_ID .EQ. 1005) THEN
	    ISTARTX=NINT(VALUE)
	  ELSE IF(ITAG_ID .EQ. 1006) THEN
	    ISTARTY=NINT(VALUE)
	  ELSE IF(ITAG_ID .EQ. 1007) THEN
	    ISPEED=NINT(VALUE)
	  ELSE IF(ITAG_ID .EQ. 1008) THEN
	    TEMP_MIN=VALUE
	  ELSE IF(ITAG_ID .EQ. 1009) THEN
	    TEMP_MAX=VALUE
	  ELSE IF(ITAG_ID .EQ. 1010) THEN
	    BEAM=TRIM(STRING)
C Extra user defined tags for Koala 2
	  ELSE IF(ITAG_ID .EQ. 1011) THEN
	    IADC=NINT(VALUE)
	  ELSE IF(ITAG_ID .EQ. 1012) THEN
	    ADC_ADD=VALUE
	  ELSE IF(ITAG_ID .EQ. 1013) THEN
	    ADC_DIV=VALUE
	  ELSE IF(ITAG_ID .EQ. 1014) THEN
	    IAPER_1=NINT(VALUE)
	  ELSE IF(ITAG_ID .EQ. 1015) THEN
	    IAPER_2=NINT(VALUE)
	  ELSE IF(ITAG_ID .EQ. 1016) THEN
	    ADC_INFO=TRIM(STRING)
	  ELSE IF(ITAG_ID .EQ. 1017) THEN
	    IKAPPA=NINT(VALUE)
C
C Increment counter for unknown user defined tags
	  ELSE IF(ITAG_ID .GE. 1000) THEN
	    NEXTRAS=NEXTRAS+1
	  ENDIF
C
C Advance the offset to the next IFD record, then loop back
	  IOFF_IFD=IOFF_IFD+12
	ENDDO
C
	CLOSE(UNIT=IUNIT)
C
C Test for ImageJ copies and big-endian TIF files
C
	IF( (COMMENT(1:7) .EQ. 'ImageJ=') .AND. (HOST .EQ. '') ) GOTO 910
	IF(ISWAP .NE. 0) GOTO 911
C
C Set BEAM, if an empty string
C
	IF(BEAM .EQ. '') BEAM='White Beam'
C
C Decode the instrument being used
C
	IF( (IXYRES .EQ. 0.0) .AND. (IADC .NE. 0) ) THEN
		HOST='KOALA2'
		IXYRES=200
		ISTARTX=0
		ISTARTY=0
		ISPEED=600
	ENDIF
C If HOST not set and COMMENTS indicates CYCLOPS, set HOST and clear COMMENT
C NB: laug_prep.exe will write the HOST tag as CYLOPS_LAUEG
	    IF(HOST .EQ. '') THEN
	      IF(INDEX(STRING,'_CCD_Cyclops_') .NE. 0) THEN
	        HOST='CYCLOPS'
	        COMMENT=''
	      ENDIF
	    ENDIF
C
C Return warning about extra tags, else return as successful
C
	ISTATUS=0
	IF(NEXTRAS .GT. 0) ISTATUS=100+NEXTRAS
	RETURN
C
900	ISTATUS=-1	! Can't open file
	RETURN
901	ISTATUS=1		! Not a TIF file (Invalid magic numbers)
	CLOSE(UNIT=IUNIT)
	RETURN
902	ISTATUS=2   ! Invalid IFD offset
	CLOSE(UNIT=IUNIT)
	RETURN
903	ISTATUS=3   ! Invalid IFD NTAGS
	CLOSE(UNIT=IUNIT)
	RETURN
910	ISTATUS=11		! Invalid Laue image (apparently an ImageJ copy)
	RETURN
911	ISTATUS=12		! Not a laue image (big-endian TIF)
	RETURN
C
	END


	SUBROUTINE CONVERT_DATE_TIME(STRING,LSTRING)
C
C Convert the date & time string STRING(1:LSTRING) to
C a 24 hour Australian format "11:24:38 26/07/2003"
C
	CHARACTER STRING*60
C
	CHARACTER TEMP*19
C
	IF(STRING(5:5) .EQ. ':') THEN
C KOALA/VIVALDI style, "2002:06:23 17:19:14"
	  TEMP=STRING(12:19)//' '//
	1		STRING(9:10)//'/'//STRING(6:7)//'/'//STRING(1:4)
	  STRING=TEMP
	ELSE
C BIO-C style, "7/26/2013 8:23:24 PM'
	  I1=1
	  IF(STRING(3:3) .EQ. '/') I1=2
	  I2=1
	  IF(STRING(I1+4:I1+4) .EQ. '/') I2=2
	  I3=1
	  IF(STRING(I1+I2+10:I1+I2+10) .EQ. ':') I3=2
	  I4=1
	  IF(STRING(I1+I2+I3+11:I1+I2+I3+11) .EQ. ':') I4=2
	  I5=1
	  IF(STRING(I1+I2+I3+I4+12:I1+I2+I3+I4+12) .EQ. ' ') I5=2
	  READ(STRING,'(I<I1>,1X,I<I2>,1X,I4,1X,I<I3>,1X,I<I4>,1X,I<I5>)')
	1				IMON,IDAY,IYEAR,IHOUR,IMIN,ISEC
	  IF(STRING(LSTRING-2:LSTRING-2) .EQ. 'P') IHOUR=IHOUR+12
	  WRITE(STRING,'(5(I2.2,A),I4)')
	1	IHOUR,':',IMIN,':',ISEC,' ',IDAY,'/',IMON,'/',IYEAR
	ENDIF
C
	LSTRING=19
C
	RETURN
	END


C===============================================================
C  Low level routines that get TIF data from the file
C  All call TIFF_READ_INT1 to actually read from the file
C===============================================================

	SUBROUTINE TIFF_GET_STRING(STRING,IOFFSET,N_VALUES)
C
	CHARACTER STRING*(*)
C
	STRING=''
C
C If a short string the characters are in the bytes of IOFFSET
C
	IF(N_VALUES .LE. 4) THEN
C This method ensures IOFFSET is unpacked as READ_TIFF_INT4() packed it
	  I4=IOFFSET
	  DO I=1,N_VALUES
	    STRING(I:I)=CHAR(IAND(I4,255))
	    I4=I4/256
	  ENDDO
C
C If a long string, get it from the bytes according to the IOFFSET
C
	ELSE
	  DO I=1,N_VALUES
	    CALL TIFF_READ_INT1(IVALUE,IOFFSET+I-1)
	    STRING(I:I)=CHAR(IVALUE)
	  ENDDO
	ENDIF
C
C Convert any trailing NUL to SPACE
C
	IF(N_VALUES .GT. 0) THEN
		IF(STRING(N_VALUES:N_VALUES) .EQ. CHAR(0)) STRING(N_VALUES:N_VALUES)=' '
	ENDIF
C
	RETURN
	END


	SUBROUTINE TIFF_GET_RATIONAL(RATIO,IOFFSET)
C
	CALL TIFF_READ_INT4(IDENOM,IOFFSET)
	CALL TIFF_READ_INT4(INUMER,IOFFSET+4)
	IF(IDENOM .EQ. 0) THEN
	  RATIO=0.0
	ELSE
	  RATIO=FLOAT(INUMER)/FLOAT(IDENOM)
	ENDIF
C
	RETURN
	END


	SUBROUTINE TIFF_GET_DOUBLE(VALUE,IOFFSET)
C
	REAL*8 VAL_R8(1)
	INTEGER*4 VAL_I4(2)
	EQUIVALENCE (VAL_I4,VAL_R8)
C
	CALL TIFF_READ_INT4(VAL_I4(1),IOFFSET)
	CALL TIFF_READ_INT4(VAL_I4(2),IOFFSET+4)
C
	VALUE=VAL_R8(1)
	RETURN
	END


	SUBROUTINE TIFF_READ_INT4(IVALUE,IOFFSET)
C
	INTEGER*1 BUFFER1(4000)
	COMMON /READ_TIF_COM/ BUFFER1,IREC,IUNIT,ISWAP
C
	CALL TIFF_READ_INT2(IVAL1,IOFFSET)
	CALL TIFF_READ_INT2(IVAL2,IOFFSET+2)
C
	IF(ISWAP .EQ. 0) THEN
	  IVALUE=IVAL1+ISHFT(IVAL2,16)
	ELSE
		IVALUE=IVAL2+ISHFT(IVAL1,16)
	ENDIF
C
	RETURN
	END


	SUBROUTINE TIFF_READ_INT2(IVALUE,IOFFSET)
C
	INTEGER*1 BUFFER1(4000)
	COMMON /READ_TIF_COM/ BUFFER1,IREC,IUNIT,ISWAP
C
	CALL TIFF_READ_INT1(IVAL1,IOFFSET)
	CALL TIFF_READ_INT1(IVAL2,IOFFSET+1)
C
	IF(ISWAP .EQ. 0) THEN
	  IVALUE=IVAL1+ISHFT(IVAL2,8)
	ELSE
	  IVALUE=IVAL2+ISHFT(IVAL1,8)
	ENDIF
C
	RETURN
	END


	SUBROUTINE TIFF_READ_INT1(IVALUE,IOFFSET)
C
	INTEGER*1 BUFFER1(4000)
	COMMON /READ_TIF_COM/ BUFFER1,IREC,IUNIT,ISWAP
C
C Read the required record if not in the buffer
C
	IREC2=IOFFSET/4000+1
	IF(IREC2 .NE. IREC) THEN
	  IREC=IREC2
C To read the last record of the file we ignore all I/O errors
	  READ(IUNIT,REC=IREC,IOSTAT=IDUMMY) BUFFER1
	ENDIF
C
C Extract the required byte from the buffer
C
	ISTART=IOFFSET-(IREC-1)*4000
	IVALUE=IAND(BUFFER1(ISTART+1),255)
C
	RETURN
	END


C===============================================================
C  Low level routines to write tags for new TIF files
C===============================================================

C ============= BUNCH OF ROUTINES TO HELP CREATE NEW TIF FILES =============

	SUBROUTINE TIFF_WRITE_I2_TAG(ITAG,IVALUE)
C Write a tag for a single INTEGER*2 value
	CALL TIFF_WRITE_I2(ITAG)
	CALL TIFF_WRITE_I2(3)
	CALL TIFF_WRITE_I4(1)
	CALL TIFF_WRITE_I2(IVALUE)
C Tags must be 12 bytes long, so have to pad an extra 2 bytes
	CALL TIFF_WRITE_I2(0)
C
	RETURN
	END


	SUBROUTINE TIFF_WRITE_I4_TAG(ITAG,IVALUE)
C Write a tag for a single INTEGER*4 value
	CALL TIFF_WRITE_I2(ITAG)
	CALL TIFF_WRITE_I2(4)
	CALL TIFF_WRITE_I4(1)
	CALL TIFF_WRITE_I4(IVALUE)
C
	RETURN
	END


	SUBROUTINE TIFF_WRITE_R8_TAG(ITAG,IOFFSET)
C Write am offset tag for a single REAL*8
	CALL TIFF_WRITE_I2(ITAG)
	CALL TIFF_WRITE_I2(12)
	CALL TIFF_WRITE_I4(1)
	CALL TIFF_WRITE_I4(IOFFSET)
C
	RETURN
	END


	SUBROUTINE TIFF_WRITE_STRING_TAG(ITAG,IOFFSET,STRING)
C Write an offset tag for a string
	CHARACTER STRING*(*)
C
	CALL TIFF_WRITE_I2(ITAG)
	CALL TIFF_WRITE_I2(2)
	CALL TIFF_WRITE_I4(LEN(STRING)+1)
	CALL TIFF_WRITE_I4(IOFFSET)
C
	RETURN
	END


	SUBROUTINE TIFF_WRITE_I2(IVALUE)
C
	INTEGER*2 I2
C
	I2=IVALUE
	WRITE(2) I2
	RETURN
	END


	SUBROUTINE TIFF_WRITE_I4(IVALUE)
C
	WRITE(2) IVALUE
	RETURN
	END


	SUBROUTINE TIFF_WRITE_STRING(STRING)
C
	CHARACTER STRING*(*)
C
	INTEGER*1 IBUFF
C
	DO I=1,LEN(STRING)
	  IBUFF=ICHAR(STRING(I:I))
	  WRITE(2) IBUFF
	ENDDO
C
	IBUFF=0
	WRITE(2) IBUFF
C
	RETURN
	END


	SUBROUTINE TIFF_WRITE_R8(VALUE)
C
	REAL*8 R8
C
	R8=VALUE
	WRITE(2) R8
C
	RETURN
	END
