	PROGRAM LAUEG_PREP
C
	CHARACTER*80 BASE_NAME,DATE_TIME,USER,SAMPLE
C
C Output a simple banner
C
	PRINT *,'CYCLOPS preparation program for LaueG (Ross Piltz, 21/7/2023)'
C
C Try to open the input file, if no file give an example of what we want
C
	OPEN(UNIT=11,FILE='laueg_cyclops.dat',STATUS='OLD',ERR=900)
C
C Get user name, sample, and date for experiment
C
	PRINT '(1X,A,$)','User name : '
	USER=REPEAT(' ',LEN(USER))
	READ(*,'(Q,A)') ILEN,USER(1:ILEN)
C
	PRINT '(1X,A,$)','Sample : '
	SAMPLE=REPEAT(' ',LEN(SAMPLE))
	READ(*,'(Q,A)') ILEN,SAMPLE(1:ILEN)
C
10	PRINT '(1X,A,$)','The experiment date as numbers "day month year" : '
	READ(*,*,IOSTAT=IERR) IDAY,IMONTH,IYEAR
	IF(IERR .NE. 0) THEN
	  PRINT *,'Enter the date as three numbers, such as: 15 3 2013'
	  GOTO 10
	ENDIF
	WRITE(DATE_TIME,'(I4,A,I2,A,I2,A)') IYEAR,':',IMONTH,':',IDAY,' 00:00:00'
	PRINT *,DATE_TIME
C
C Read the input file line by line processing files as we go
C
100	BASE_NAME=REPEAT(' ',80)
	READ(11,*,END=200,ERR=901) BASE_NAME,NFILES,INDEX,PHI
	PRINT '(1X,A,I4.4,A,I3.3,A,F6.1,A)',TRIM(BASE_NAME)//'0001 to ',NFILES,
	1			'.tif ==> laueg_prep_',INDEX,'.tif  (Phi =',PHI,')'
	CALL PROCESS_IMAGES(BASE_NAME,INDEX,DATE_TIME,USER,SAMPLE,PHI,NFILES)
	GOTO 100
C
200	PRINT *,'SUCCESSFUL COMPLETION'
	GOTO 1000
C
C Output instructions if we can't read or find the input file
C
900	PRINT *,'ERROR: Input file "laueg_cyclops.dat" is missing'
	GOTO 950
901	PRINT *,'ERROR: Input file "laueg_cyclops.dat" is invalid'
C
950	PRINT *
	PRINT *,'Use an editor to create an input file with 4 parameters per line:'
	PRINT *,'1) Base-name of the CYCLOPS files (exclude the 4 digit number and ".tif")'
	PRINT *,'2) Number of files for this base-name'
	PRINT *,'3) The index number used for the "laueg_*.tif" output file'
	PRINT *,'4) The OMEGA angle for these files (NB: PHI = -OMEGA)'
	PRINT *
	PRINT *,'An example input file is:'
	PRINT *,' ano_300K_w10G20_ 3 11 10'
	PRINT *,' ano_300K_w15G20_ 3 12 15'
	PRINT *,' ano_300K_w20G20_ 3 13 20'
	PRINT *
	PRINT *,'The first line instructs processing of files "ano_300K_w10G20_0001.tif",'
	PRINT *,'"ano_300K_w10G20_0002.tif", and "ano_300K_w10G20_0003.tif" with the results'
	PRINT *,'written to "laueg_011.tif". The Phi angle of these exposures, 10 degrees,'
	PRINT *,'will be stored in the header of the output file.'
	PRINT *,'The second line converts "ano_300K_w15G20_000*.tif" to "laueg_012.tif",'
	PRINT *,'and similar for the remaining lines'
C
1000	END


	SUBROUTINE PROCESS_IMAGES(BASE_NAME,INDEX,DATETIME,USER,SAMPLE,PHI,NFILES)
C
	CHARACTER*(*) BASE_NAME,USER,SAMPLE,DATETIME
C
	CHARACTER*2000 COMMENT,SAMPLE2,USER2,DATETIME2,HOST,BEAM,ADC_INFO
	COMMON /TIFFINFO_COM/ COMMENT,SAMPLE2,USER2,DATETIME2,HOST,BEAM,
	1			IOFF_DATA,NUMX,NUMY,IXYRES,ISTARTX,ISTARTY,ISPEED,
	2			EXPOSE_TIME,EXPOSE_PHI,TEMP_BEGIN,TEMP_END,TEMP_MIN,TEMP_MAX,
	3			IADC,ADC_ADD,ADC_DIV,ADC_INFO,IAPER_1,IAPER_2,IKAPPA
C
	CHARACTER FILE_NAME*80,COMMENTS*240
	INTEGER*2 IMAGE2(20000000),IMAGE2B(5000000)
	INTEGER*4 IMAGE4(20000000),IMAGE4B(5000000)
C
	DO IFILE=1,NFILES
C
	  WRITE(FILE_NAME,'(A,I4.4,A)') TRIM(BASE_NAME),IFILE,'.tif'
	  PRINT '(1X,2A)','Reading data file: ',TRIM(FILE_NAME)
	  CALL CHECK_LAUE_TIFF(NUMX,NUMY, FILE_NAME)
	  IF(NUMX .LT. 0)CALL QUIT('ERROR: Unable to read data file')
	  IF(NUMX*NUMY .GT. 20000000) CALL QUIT('ERROR: Image too large (> 20M pixels)')
	  IF(TRIM(HOST) .NE. 'CYCLOPS') CALL QUIT('ERROR: Image file is not a raw CYCLOPS file')
C
	  CALL READ_LAUE_TIFF_QUIET(IMAGE2,NUMX,NUMY, FILE_NAME)
C
	  IF(IFILE .EQ. 1) THEN
	    DO I=1,NUMX*NUMY
	      IMAGE4(I)=ZEXT(IMAGE2(I))
	    ENDDO
	  ELSE
	    DO I=1,NUMX*NUMY
	      IMAGE4(I)=IMAGE4(I)+ZEXT(IMAGE2(I))
	    ENDDO
	  ENDIF
C
	ENDDO
C
C First remove any hot pixels
C
	PREC=0.02
	CUTOFF=100
	DO ITER=1,4
	NBAD=0
	DO IX=4,NUMX-3
	  DO IY=4,NUMY-3
	    I=IX+NUMX*(IY-1)
	    AVEX=( IMAGE4(I-3)+IMAGE4(I-2)+IMAGE4(I+2)+IMAGE4(I+3) )/4.0
	    AVEY=( IMAGE4(I-3*NUMX)+IMAGE4(I-2*NUMX)+
	1			IMAGE4(I+2*NUMX)+IMAGE4(I+3*NUMX) )/4.0
	    SIGX=SQRT(AVEX)+AVEX*PREC
	    SIGY=SQRT(AVEY)+AVEY*PREC
	    IF(IMAGE4(I) .GT. AVEX+CUTOFF*SIGX) THEN
	      IF(IMAGE4(I) .GT. AVEY+CUTOFF*SIGY) THEN
	        IMAGE4(I)=MIN(AVEX,AVEY)
	        NBAD=NBAD+1
	      ENDIF
	    ENDIF
	  ENDDO
	ENDDO
	PRINT '(I7,A)',NBAD,' hot pixels removed'
	CUTOFF=2.0+20.0/ITER
	ENDDO
C
C Second remove any cold pixels
C
	PREC=0.02
	CUTOFF=10.0
	NBAD=0
	DO IX=4,NUMX-3
	  DO IY=4,NUMY-3
	    I=IX+NUMX*(IY-1)
	    AVEX=( IMAGE4(I-3)+IMAGE4(I-2)+IMAGE4(I+2)+IMAGE4(I+3) )/4.0
	    AVEY=( IMAGE4(I-3*NUMX)+IMAGE4(I-2*NUMX)+
	1			IMAGE4(I+2*NUMX)+IMAGE4(I+3*NUMX) )/4.0
	    SIGX=SQRT(AVEX)+AVEX*PREC
	    SIGY=SQRT(AVEY)+AVEY*PREC
	    IF(IMAGE4(I) .LT. AVEX-CUTOFF*SIGX) THEN
	      IF(IMAGE4(I) .LT. AVEY-CUTOFF*SIGY) THEN
	        IMAGE4(I)=MAX(AVEX,AVEY)
	        NBAD=NBAD+1
	      ENDIF
	    ENDIF
	  ENDDO
	ENDDO
	PRINT '(I7,A)',NBAD,' cold pixels removed'
C
C Zero the border intensities
C
	NBORDER=10
C
	DO IY=1,NBORDER
	  DO IX=1,NUMX
	    IMAGE4(IX+NUMX*(IY-1))=0
	    IMAGE4(IX+NUMX*(NUMY-IY))=0
	  ENDDO
	ENDDO
C
	DO IX=1,NBORDER
	  DO IY=1,NUMY
	    IMAGE4(IX+NUMX*(IY-1))=0
	    IMAGE4(NUMX-IX+1+NUMX*(IY-1))=0
	  ENDDO
	ENDDO
C
	I=0
	DO IY=1,NUMY,2
	  DO IX=1,NUMX,2
	    IXY=IX+NUMX*(IY-1)
	    ISUM=IMAGE4(IXY)+IMAGE4(IXY+1)+IMAGE4(IXY+NUMX)+IMAGE4(IXY+1+NUMX)
	    IMAX=MAX(IMAGE4(IXY),IMAGE4(IXY+1),IMAGE4(IXY+NUMX),IMAGE4(IXY+1+NUMX))
	    IMIN=MIN(IMAGE4(IXY),IMAGE4(IXY+1),IMAGE4(IXY+NUMX),IMAGE4(IXY+1+NUMX))
	    I=I+1
	    IMAGE4B(I)=(ISUM-IMAX-IMIN)/2
	  ENDDO
	ENDDO
C
	NUMX2=NUMX/2
	NUMY2=NUMY/2
	DO I=1,NUMX2*NUMY2
	  ITEMP=IMAGE4B(I)/NFILES
	  IF(ITEMP .GE. 128*256) ITEMP=ITEMP-256*256
	  IMAGE2B(I)=ITEMP
	ENDDO
	PRINT '(1X,A,2I5)','Reducing image size to',NUMX2,NUMY2
C
	COMMENTS=REPEAT(' ',LEN(COMMENTS))
	WRITE(COMMENTS,'(4A,I4.4,A)') 'Created by laueg_prep.exe from files ',
	1	TRIM(BASE_NAME),'0001.tif to ',TRIM(BASE_NAME),NFILES,'.tif'
C
	WRITE(FILE_NAME,'(A,I3.3,A)') 'laueg_prep_',INDEX,'.tif'
	PRINT *,'Creating output file "',TRIM(FILE_NAME),'"'
	CALL WRITE_LAUE_TIFF(TRIM(FILE_NAME), IMAGE2B,NUMX2,NUMY2,
	1			'CYCLOPS_LAUEG',TRIM(USER),TRIM(SAMPLE),
	2			TRIM(DATETIME),TRIM(COMMENTS),PHI_SET)
	PRINT *
C
	RETURN
	END



	SUBROUTINE GET_NUM_XY(NUMX0,NUMY0)
C
	REAL IMAGE_RAW(8000,2500),IMAGE(8000,2500)
	COMMON /IMAGE_COM/ IMAGE_RAW,IMAGE,NUMX,NUMY
C
	NUMX0=NUMX
	NUMY0=NUMY
C
	RETURN
	END




	FUNCTION COUNTS_RAW(IX,IY)
C
	REAL IMAGE_RAW(8000,2500),IMAGE(8000,2500)
	COMMON /IMAGE_COM/ IMAGE_RAW,IMAGE,NUMX,NUMY
C
	IF(IX.GE.1 .AND. IX.LE.NUMX .AND. IY.GE.1 .AND. IY.LE.NUMY) THEN
	  COUNTS_RAW=IMAGE_RAW(IX,IY)
	ELSE
	  COUNTS_RAW=0.0
	ENDIF
C
	RETURN
	END
