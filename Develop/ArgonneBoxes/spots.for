C ========================== Routines in this file ============================
C	SUBROUTINE REFINE_SPOTS_CENTERS(RAD_MULT,ZCUTOFFS,TOLER_NBOUR,A_MAX,NREFS)
C	SUBROUTINE RECENTER_SPOTS(X,Y, ITWIN,NREFS)
C =============================================================================

	SUBROUTINE REFINE_SPOTS_CENTERS(RAD_MULT,ZCUTOFFS,TOLER_NBOUR,A_MAX,NREFS)
C
	REAL RAD_MULT(3),ZCUTOFFS(5)
C
	COMMON /IZONES_COM/ IZONES(20000)
	COMMON /LIB_COM/ PLIB(1000,6),ILIB(1000),NLIB
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /SPOT2_COM/ IHKLM(4,20000),WAV(20000),TTH(20000),MULT(20000)
	COMMON /BGLOCAL_COM/ BKG(20000),BKG_ESD(20000),BKG_VAR(20000)
C
	COMMON /OPTIONS_COM/ IOPT_HKL,NHKLM
C
	INTEGER NSUM(3)
	REAL CSUM(3)
C
C Output some info to the log files
C
	WRITE(6,'(1X,A)') 'Refining spot centres .....'
	WRITE(30,'(//,A,/)') '================== Refining spot centres ====================='
	WRITE(30,'(A)') 'Observed spot centres'
C
C----------------------------------------------------
C Start of loop building a table of good centered spots
C
	NLIB=0
	DO IR=1,NREFS
C
C Reasons to skip a spot:-
C Too weak:
	  PKBG=(RMAX(IR)-BKG(IR))/SQRT(MAX(10.0, BKG(IR) ))
	  IF(PKBG .LT. ZCUTOFFS(IZONES(IR))) CYCLE
C Within 20 pixels to edges, holes, etc:
	  CALL REJECT_RADIUS(IOUT, X(IR),Y(IR),20.0)
	  IF(IOUT .EQ. 1) CYCLE
C Second spot of closely overlapped twin:
	  IF(IOPT_HKL.EQ.2 .AND. IHKLM(4,IR).EQ.52) CYCLE
C Partially overlapped twin:
	  IF(IOPT_HKL.EQ.2 .AND. IHKLM(4,IR).GT.60) CYCLE
C A satellite reflection:
	  IF(IOPT_HKL.EQ.3 .AND. IHKLM(4,IR).NE.0) CYCLE
C
C Couldn't make contour ellipse:
	  CALL CALC_CONTOUR_CENTER(XSPOT,YSPOT, E,F,H,FFRAC, A_MAX, IR,ISTATUS)
	  IF(ISTATUS .EQ. 1) CYCLE
C Pixel intensity overload:
	  IF(GET_ELLIPSE_MAXINT(XSPOT,YSPOT,E,F,H) .GT. 6.0E4) CYCLE
C Peak ellipse includes edge, holes, etc:
	  CALL REJECT_ELLIPSE(ISTATUS, XSPOT,YSPOT,E,F,H,RAD_MULT(2))
	  IF(ISTATUS .NE. 0) CYCLE
C
C Calculate degree of spot overlap
	  CALL CALC_NEIGH_OVERLAP(ROVER,IOVER2, IR, E,F,H, 0.8,RAD_MULT, .FALSE.)
C
C Core-core overlap:
	  IF(IOVER2 .EQ. 4) CYCLE
C Strong neighbour overlap:
	  IF(ROVER .GT. 0.5*TOLER_NBOUR) CYCLE
C
C Calculate core and peak intensities
	  CALL SUM_ELLIPSE_INTS(CSUM,NSUM, XSPOT,YSPOT,E,F,H,RAD_MULT,2)
	  CORE=CSUM(1)-NSUM(1)*BKG(IR)
	  PEAK=CSUM(2)-NSUM(2)*BKG(IR) - CORE
	  DCORE=SQRT( CSUM(1) + (BKG_ESD(IR)*NSUM(1))**2 )
	  DPEAK=SQRT( ABS(CSUM(2)-CSUM(1)) +
	1				(BKG_ESD(IR)*(NSUM(2)-NSUM(1)))**2 )

C Peak fraction or esd are ridiculous:
	  PFRAC=CORE/(CORE+PEAK)
	  DPFRAC=DPEAK/MAX(DPEAK/4.0,CORE+PEAK)		! approximation
	  IF(PFRAC.LE.0.5 .OR. PFRAC.GE.1.1 .OR. DPFRAC.GE.0.3) CYCLE
C
C If there is no space in model library, exit loop
	  NLIB=NLIB+1
	  IF(NLIB .GE. 1000) EXIT
C
C Output information on spot
	  IM=IHKLM(4,IR)-(IHKLM(4,IR)/10)*10
C Non-twin reflections, including main reflections for modulations
	  IF(IOPT_HKL .NE. 2) THEN
	    IM=0
	    WRITE(30,'(A,I5,3X,A,2I5,3X,A,F6.1,3X,A,F5.2,3X,A,2F5.1)',IOSTAT=IDUMMY)
	1				'Spot',IR,'XY',NINT(X(IR)),NINT(Y(IR)),'Peak/esd',PKBG,
	2				'Fill',FFRAC,'dX,dY',XSPOT-X(IR),YSPOT-Y(IR)
C Closely overlapped twins
	  ELSEIF(IHKLM(4,IR) .EQ. 51) THEN
	    IM=3
	    WRITE(30,'(A,I5,3X,A,2I5,3X,A,F6.1,3X,A,F5.2,3X,A,2F5.1)',IOSTAT=IDUMMY)
	1				'Spot',IR,'Twin 1+2  XY',NINT(X(IR)),NINT(Y(IR)),
	2				'Peak/esd',PKBG,'Fill',FFRAC,'dX,dY',XSPOT-X(IR),YSPOT-Y(IR)
C Non-overlapped twin 1
	  ELSEIF(IM .EQ. 1) THEN
	    WRITE(30,'(A,I5,3X,A,2I5,3X,A,F6.1,3X,A,F5.2,3X,A,2F5.1)',IOSTAT=IDUMMY)
	1				'Spot',IR,'Twin 1    XY',NINT(X(IR)),NINT(Y(IR)),
	2				'Peak/esd',PKBG,'Fill',FFRAC,'dX,dY',XSPOT-X(IR),YSPOT-Y(IR)
C Non-overlapped twin 2
	  ELSEIF(IM .EQ. 2) THEN
	    WRITE(30,'(A,I5,3X,A,2I5,3X,A,F6.1,3X,A,F5.2,3X,A,2F5.1)',IOSTAT=IDUMMY)
	1				'Spot',IR,'Twin 2    XY',NINT(X(IR)),NINT(Y(IR)),
	2				'Peak/esd',PKBG,'Fill',FFRAC,'dX,dY',XSPOT-X(IR),YSPOT-Y(IR)
C Debug message
	  ELSE
	    CALL QUIT("Unexpected IM in MAIN at 296")
	  ENDIF
C
C Add spot center info to model library
	  ILIB(NLIB)  =IM
	  PLIB(NLIB,1)=X(IR)
	  PLIB(NLIB,2)=Y(IR)
	  PLIB(NLIB,3)=PKBG
	  PLIB(NLIB,4)=XSPOT-X(IR)
	  PLIB(NLIB,5)=YSPOT-Y(IR)
C
	ENDDO
C
C End of loop building a table of good centered spots
C----------------------------------------------------
C
C Recenter all spots using the center offsets in PLIB()
C
	IF(IOPT_HKL .NE. 2) THEN
	  CALL RECENTER_SPOTS(X,Y, 0,NREFS)
	ELSE
	  WRITE(30,'(/,A)') 'TWIN 1:'
	  CALL RECENTER_SPOTS(X,Y, 1,NREFS)
	  WRITE(30,'(/,A)') 'TWIN 2:'
	  CALL RECENTER_SPOTS(X,Y, 2,NREFS)
	ENDIF
	WRITE(30,'(/,A)') 'X and Y offsets added to all spot positions'
C
	RETURN
	END


	SUBROUTINE RECENTER_SPOTS(X,Y, ITWIN,NREFS)
C
	REAL X(20000),Y(20000)
C
	COMMON /SPOT2_COM/ IHKLM(4,20000),WAV(20000),TTH(20000),MULT(20000)
	COMMON /LIB_COM/ PLIB(1000,6),ILIB(1000),NLIB
C
	REAL XOFFSET(0:40,0:20),YOFFSET(0:40,0:20)
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
	ISTEPX=NUMX/40
	ISTEPY=NUMY/20
	DIST0=200
C
C Count how many spots of type ITWIN, and exit if < 10
C
	NUM=0
	DO I=1,NLIB
	  IF(ILIB(I).EQ.ITWIN .OR. ILIB(I).EQ.3) NUM=NUM+1
	ENDDO
	IF(NUM .LT. 10) THEN
	  WRITE(30,'(/,A)') 'WARNING: Less than 10 good spots, X and Y unchanged'
	  RETURN
	ENDIF
C
C Make grids of 41 x 21 for X and Y spot center offsets from the
C table of spots in PLIB(). DIST0 is the approximate distance
C over which results from similar "quality" spots are averaged.
C
	DO IX=0,40
	  X0=IX*ISTEPX
	  DO IY=0,20
	    Y0=IY*ISTEPY
	    SUMW=0.0
	    SUMDX=0.0
	    SUMDY=0.0
	    DO I=1,NLIB
	      IF(ILIB(I).EQ.ITWIN .OR. ILIB(I).EQ.3) THEN
	        DIST=SQRT( (PLIB(I,1)-X0)**2 + (PLIB(I,2)-Y0)**2 )
	        WGHT=EXP(-DIST/DIST0) * SQRT(PLIB(I,3))
	        SUMW=SUMW+WGHT
	        SUMDX=SUMDX+ WGHT*PLIB(I,4)
	        SUMDY=SUMDY+ WGHT*PLIB(I,5)
	      ENDIF
	    ENDDO
	    XOFFSET(IX,IY)=SUMDX/SUMW
	    YOFFSET(IX,IY)=SUMDY/SUMW
	  ENDDO
	ENDDO
C
	WRITE(30,'(/,A)') 'Smoothed spot X offsets for X(horiz) by Y(vert)'
	WRITE(30,'(I10,1X,10I6)',IOSTAT=IDUMMY) (K*ISTEPX,K=0,40,4)
	DO IY=0,20,2
	  WRITE(30,'(I5,21F6.1)',IOSTAT=IDUMMY)
	1				IY*ISTEPY,(XOFFSET(K,IY),K=0,40,4)
	ENDDO
C
	WRITE(30,'(/,A)') 'Smoothed spot Y offsets for X(horiz) by Y(vert)'
	WRITE(30,'(I10,1X,10I6)',IOSTAT=IDUMMY) (K*ISTEPX,K=0,40,4)
	DO IY=0,20,2
	  WRITE(30,'(I5,21F6.1)',IOSTAT=IDUMMY)
	1				IY*ISTEPY,(YOFFSET(K,IY),K=0,40,4)
	ENDDO
C
C Correct X & Y centers using bilinear interpolation of grid values
C
	DO IR=1,NREFS
	  IM=IHKLM(4,IR)-(IHKLM(4,IR)/10)*10
	  IF(IM .NE. ITWIN) CYCLE
C
	  XSTEPS=X(IR)/ISTEPX
	  YSTEPS=Y(IR)/ISTEPY
	  IX=MAX(0,MIN(39, IFIX(XSTEPS) ))
	  IY=MAX(0,MIN(19, IFIX(YSTEPS) ))
	  XSTEPS=XSTEPS-IX
	  YSTEPS=YSTEPS-IY
C
	  XGRAD=XOFFSET(IX+1,IY)-XOFFSET(IX,IY)
	  YGRAD=XOFFSET(IX,IY+1)-XOFFSET(IX,IY)
	  XYGRAD=XOFFSET(IX+1,IY)+XOFFSET(IX,IY)-XOFFSET(IX+1,IY)-XOFFSET(IX,IY+1)
	  XVAL=XOFFSET(IX,IY) + XGRAD*XSTEPS + YGRAD*YSTEPS + XYGRAD*XSTEPS*YSTEPS
	  X(IR)=X(IR)+XVAL
C
	  XGRAD=YOFFSET(IX+1,IY)-YOFFSET(IX,IY)
	  YGRAD=YOFFSET(IX,IY+1)-YOFFSET(IX,IY)
	  XYGRAD=YOFFSET(IX+1,IY)+YOFFSET(IX,IY)-YOFFSET(IX+1,IY)-YOFFSET(IX,IY+1)
	  YVAL=YOFFSET(IX,IY) + XGRAD*XSTEPS + YGRAD*YSTEPS + XYGRAD*XSTEPS*YSTEPS
	  Y(IR)=Y(IR)+YVAL
C
	ENDDO
C
	RETURN
	END
