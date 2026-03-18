	PROGRAM REJECTS_DATA
C
	CHARACTER*80 STRING1,STRING2,COMMENT
	CHARACTER*250 REJ_NAME,ELL_NAME
	INTEGER HKLS(3)
	REAL XYCEN(2),EFH(3),AMULT(3)
C
C Output a simple banner
C
	PRINT '(1X,A)','Rejects file processing for LaueG (Ross Piltz, 19/8/2019)'
	PRINT *
C
C Delete output file, if it still exists
C
	OPEN(UNIT=1,FILE='___laueg_rejects_data.out',STATUS='OLD',IOSTAT=IDUMMY)
	CLOSE(1,STATUS='DELETE',IOSTAT=IDUMMY)
C
C Read reject file name from input file and delete input file
C
	OPEN(UNIT=1,FILE='___laueg_rejects_data.in',STATUS='OLD',IOSTAT=IERR)
	IF(IERR .NE. 0) CALL QUIT('ERROR: Unable to open input file')
	REJ_NAME=REPEAT(' ',LEN(REJ_NAME))
	READ(1,'(Q,A)',ERR=900) ILEN,REJ_NAME(1:ILEN)
	CLOSE(1)
C
C Open reject file
C
	PRINT *,'Reading the reject file = ',TRIM(REJ_NAME)
	OPEN(UNIT=1,STATUS='OLD',FILE=TRIM(REJ_NAME),ERR=910)
C
C Read header with version number and options and complain and die if invalid
C
	READ(1,'(A26,I2,A8,I2)',ERR=920) STRING1,IVER,STRING2,IOPT
	IF(STRING1 .NE. 'LaueG Rejects File Version') CALL QUIT('ERROR: Invalid file header')
	IF(IVER.LT.1 .OR. IVER.GT.2) CALL QUIT('ERROR: Invalid file version')
	IF(STRING2 .NE. ', Option') CALL QUIT('ERROR: Invalid file header')
	IF(IOPT.LT.1 .OR. IOPT.GT.2) CALL QUIT('ERROR: Invalid file option')
C
	ISIZE=20
	IF(IVER .GE. 2) ISIZE=40
C
	IF(IOPT .NE. 1) THEN
	  PRINT *,'Reject spots are twins'
	ENDIF
C
C Skip 2 comment lines
C
	READ(1,*)
	READ(1,*)
C
C Open the output file on unit 3
C
	OPEN(UNIT=3,FILE='___laueg_rejects_data.out',STATUS='UNKNOWN')
C
C Start the loop to read single reject spots
C If twins are not set, use ITWIN=0
	ITWIN=0
100	ELL_NAME=REPEAT(' ',LEN(ELL_NAME))
	IF(IOPT .EQ. 1) THEN
	  READ(1,'(A<ISIZE>,3I4,A)',END=200,ERR=930) ELL_NAME,HKLS,COMMENT
	ELSE
	  READ(1,'(A<ISIZE>,4I4,A)',END=200,ERR=930) ELL_NAME,HKLS,ITWIN,COMMENT
	ENDIF
C
C Try to find the ellipse that matches the spot
C Returns ITYPE=-1 if no match
C
	PRINT *,'Reading ellipse file = ',TRIM(ELL_NAME)//'.ell'
	CALL FIND_ELLIPSE(XYCEN,EFH,AMULT,ITYPE, TRIM(ELL_NAME)//'.ell',HKLS,ITWIN)
C
C Write the results to the output file
C
	WRITE(3,'(A<ISIZE>,3I4,I3,2F7.1,3P,3F8.3,0P,3F4.1,I3,A)')
	1	ELL_NAME,HKLS,ITWIN,XYCEN,EFH,AMULT,ITYPE,TRIM(COMMENT)
C
C Loop back for more spots
C
	GOTO 100
C
C Several error messages
C
900	CALL QUIT('ERROR: Unable to read input file')
910	CALL QUIT('ERROR: Unable to open reject file')
920	CALL QUIT('ERROR: Unable to read reject file header')
930	CALL QUIT('ERROR: Unable to read reject file data')
C
C Finish up
C
200	CLOSE(UNIT=1)
	CLOSE(UNIT=3)
C
	CALL DELETE_FILE('___laueg_rejects_data.in')
C
	PRINT *,'SUCCESSFUL COMPLETION'
	END



	SUBROUTINE FIND_ELLIPSE(XYCEN,EFH,AMULT,ITYPE, ELL_NAME,HKLS,ITWIN)
C
C Returns XYCEN,EFH,AMULT & ITYPE values for the first matching ellipse.
C Returns a negative ITYPE is the ellipse is not found.
C
	CHARACTER ELL_NAME*(*)
	INTEGER HKLS(3)
	REAL XYCEN(2),EFH(3),AMULT(3)
C
	CHARACTER*80 STRING1,STRING2
C
C Set ITYPE to indicate can't find spot in *.ell file
C
	ITYPE=-1
C
C Try to open the file
C
	OPEN(UNIT=2,STATUS='OLD',FILE=ELL_NAME,ERR=900)
C
C Read header with version number/option, complain and die if invalid
C
	READ(2,'(A26,I2,A8,I2)') STRING1,IVER,STRING2,IOPT
	IF(STRING1 .NE. 'LaueG Ellipse File Version') CALL QUIT('ERROR: Invalid file header')
	IF(IVER .NE. 1) CALL QUIT('ERROR: Invalid file version')
	IF(STRING2 .NE. ', Option') CALL QUIT('ERROR: Invalid file header')
	IF(IOPT.LT.1 .OR. IOPT.GT.2) CALL QUIT('ERROR: Invalid file option')
	IF(ITWIN.EQ.0 .AND. IOPT.NE.1) CALL QUIT('ERROR: File has twin information')
	IF(ITWIN.NE.0 .AND. IOPT.NE.2) CALL QUIT('ERROR: File has no twin information')
C
C Skip 2 comment lines
C
	READ(2,*)
	READ(2,*)
C
C Start the loop to read individual ellipses
C
100	IF(IOPT .EQ. 1) THEN
	  READ(2,'(3I4,2F7.1,3E11.3,1X,3F4.1,I3)',END=200,ERR=910)
 	1				IH,IK,IL,XYCEN,EFH,AMULT,ITYPE
	ELSE
	  READ(2,'(3I4,2F7.1,3E11.3,1X,3F4.1,2I3)',END=200,ERR=910)
 	1				IH,IK,IL,XYCEN,EFH,AMULT,ITYPE,ITWIN2
	ENDIF
C
C If the HKL (+ twin) do not match, loop back to read another ellipse
C
	IF(HKLS(1).NE.IH .OR. HKLS(2).NE.IK .OR. HKLS(3).NE.IL) GOTO 100
	IF(IOPT.EQ.2 .AND. ITWIN2.NE.ITWIN) GOTO 100
C
C Finished, either found a match or hit the EOF
C
200	CLOSE(UNIT=2)
	RETURN
C
C Some more error messages
C
900	CALL QUIT('ERROR: Cannot find or open file')
910	CALL QUIT('ERROR: Invalid file data')
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
