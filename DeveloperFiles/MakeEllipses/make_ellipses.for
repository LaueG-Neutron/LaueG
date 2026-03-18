	PROGRAM DRAW_ELLIPSES
C
C Calculates the individual pixels that will be drawn given a set of ellipses.
C The input file consists of groups of similar ellipses, the output file is a
C list of pixel X & Y for the pixels to draw, followed by an index table so
C individual ellipses can be found.
C
	INTEGER IXY(2,10000)
	INTEGER IPOINT(10000)
C
C Output a simple banner
C
	PRINT '(1X,A)','Ellipse drawing program for LaueG (Ross Piltz, 15/12/2014)'
	PRINT *
C
	OPEN(UNIT=1,FILE='___laueg_make_ellipses.in',STATUS='OLD',IOSTAT=IERR)
	IF(IERR .NE. 0) CALL QUIT('ERROR: Unable to open input file')
	OPEN(UNIT=2,FILE='___laueg_make_ellipses.out',STATUS='UNKNOWN',IOSTAT=IERR)
	IF(IERR .NE. 0) CALL QUIT('ERROR: Unable to open output file')
C
C Loop for reading groups of ellipses from the input file
C
C Read a header line for a group of NREAD ellipses
	IGROUP=0
	IWRITE=0
100	READ(1,*,END=200) NREAD,ELLIM,IXLO,IXHI,IYLO,IYHI
	IF(NREAD .EQ. 0) GOTO 200
C Read in the ellipses in this group and output the pixels to draw
	DO I=1,NREAD
	  READ(1,*) XCEN,YCEN,E,F,H
	  CALL CALC_ELLIPSE_PIXELS(IXY,NXY, XCEN,YCEN,E,F,H,ELLIM,
	1					IXLO,IXHI,IYLO,IYHI)
	  WRITE(2,'(I5,1X,I6)') (IXY(1,K),IXY(2,K),K=1,NXY)
	  IWRITE=IWRITE+NXY
	ENDDO
C Save the number of the last pixel written to the file
	IGROUP=IGROUP+1
	IPOINT(IGROUP)=IWRITE
C Jump back for the next group of ellipses
	GOTO 100
C Add a list of (99999,IPOINT) "XY" values at the end of the file
C which contains in "Y" the number of the last pixel in the group
200	DO I=1,IGROUP
	  WRITE(2,*) 99999,IPOINT(I)
	ENDDO
C
	CLOSE(UNIT=1)
	CLOSE(UNIT=2)
C
	CALL DELETE_FILE('___laueg_make_ellipses.in')
C
	PRINT *,'SUCCESSFUL COMPLETION'
	END



	SUBROUTINE CALC_ELLIPSE_PIXELS(IXY,NXY, XCEN,YCEN,E,F,H,ELLIM,
	1					IXLO,IXHI,IYLO,IYHI)
C
	INTEGER IXY(2,1000000)
C
C Start with no pixels
C
	NXY=0
C
C Nothing to do if ellipse is completely outside the X,Y limits
C
	XWIDTH=SQRT(ELLIM*F/(E*F-H**2))
	IF( XCEN+XWIDTH.LT.FLOAT(IXLO) .OR.
	1    XCEN-XWIDTH.GT.FLOAT(IXHI)     ) RETURN
	YWIDTH=SQRT(ELLIM*E/(E*F-H**2))
	IF( YCEN+YWIDTH.LT.FLOAT(IYLO) .OR.
	1    YCEN-YWIDTH.GT.FLOAT(IYHI)     ) RETURN
C
C Loop through all possible X values finding Y where X,Y is
C inside the ellipse but within +/-1 in X or Y of being outside
C
	DO IX=CEILING(XCEN-XWIDTH),FLOOR(XCEN+XWIDTH)
	  DY=SQRT( (1.0 - ((IX-XCEN)/XWIDTH)**2 )*ELLIM/F )
	  Y0=YCEN-(IX-XCEN)*H/F
C Set 1 or 2 loops on IY with different limits
	  IF( IX.LE.CEILING(XCEN-XWIDTH)+1 .OR.
	1		IX.GE.FLOOR  (XCEN+XWIDTH)-1       ) THEN
C If close to the IX limits, test the complete set of IY
	    IY_LO=CEILING(Y0-DY)
	    IY_HI=FLOOR(Y0+DY)
	    NLOOP=1
	  ELSE
C First test 3 IY near the minimum, then loop for the maximum
	    IY_LO=CEILING(Y0-DY)
	    IY_HI=IY_LO+2
	    NLOOP=2
	  ENDIF
C Find pixels +/-1 in X or Y from the outside of the ellipse
	  DO I=1,NLOOP
	    DO IY=IY_LO,IY_HI
	      DELX=IX-XCEN
	      DELY=IY-YCEN 
	      DSQ0=E*DELX**2+F*DELY**2+2.0*H*DELX*DELY
	      DSQ_MAX=DSQ0 + MAX( E+2.0*ABS(E*DELX+H*DELY) ,
	1				F+2.0*ABS(F*DELY+H*DELX) )
	      IF(DSQ_MAX .GT. ELLIM) THEN
	        IF( IX.GE.IXLO .AND. IX.LE.IXHI .AND.
	1            IY.GE.IYLO .AND. IY.LE.IYHI      ) THEN
	          NXY=NXY+1
	          IXY(1,NXY)=IX
	          IXY(2,NXY)=IY
	        ENDIF
	      ENDIF
	    ENDDO
C If NLOOP=2, then loop back with these IY limits
	    IY_HI=FLOOR(Y0+DY)
	    IY_LO=IY_HI-2
	  ENDDO
C
	ENDDO
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
