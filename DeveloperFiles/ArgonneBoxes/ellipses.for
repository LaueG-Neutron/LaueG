C ===================== Routines in this file ==============================
C	SUBROUTINE CALC_CONTOUR_CENTER(XCEN,YCEN, E,F,H,FFRAC, A_MAX, IR,ISTATUS)
C	SUBROUTINE CALC_CONTOUR_SHAPE(E,F,H,FFRAC,A,B,ANGLE, DXCEN,DYCEN,
C	1										SHKLM, CONTOUR,FFILL,A_MAX, IR,ISTATUS)
C	SUBROUTINE CALC_ELL_MOMENTS(DXCEN,DYCEN,A,B,ANGLE, XCEN,YCEN, E,F,H,
C	1								RLEVEL_LO,RLEVEL_HI, LFIXED)
C	FUNCTION GET_ELLIPSE_FILL(XCEN,YCEN,E,F,H, RLEVEL_LO)
C	SUBROUTINE SUM_ELLIPSE_INTS(CSUM,NSUM, XCEN,YCEN,E,F,H, RAD,NRAD)
C	FUNCTION GET_ELLIPSE_MAXINT(XCEN,YCEN,E,F,H)
C	SUBROUTINE AVE_ELLIPSE_BKGS(BKG_AVE, XCEN,YCEN,E,F,H, RAD,NRAD)
C ==================== Some short ellipse handling routines ==============
C	SUBROUTINE CALC_COG2AXES(A,B,ANGLE, T11,T22,T12)
C	SUBROUTINE CALC_EFH2AXES(A,B,ANGLE, E,F,H)
C	SUBROUTINE FIX_ELLI_AXES(A,B,ANGLE)
C	SUBROUTINE CALC_AXES2EFH(E,F,H, A0,B0,ANGLE)
C	FUNCTION CALC_ELLIPSE_DSQ(DX,DY,E,F,H)
C ===================== Unused test & Obsolete routines ====================
C	SUBROUTINE TEST_ELLIPSES()
C	SUBROUTINE GET_ELLIPSE_TEST(XCEN,YCEN,E,F,H,RMAX,CONTOUR)
C	SUBROUTINE CALC_CONTOUR_ELLIPSE_OLD(E,F,H,CG,FFRAC,A,B,ANGLE,
C	1							XCEN,YCEN,RMAX, CONTOUR, LFIXED, ISTATUS)
C	SUBROUTINE CALC_CONTOUR_ELLIPSE_OLDER(E,F,H,CG,FFRAC,A,B,ANGLE,
C	1							XCEN,YCEN,RMAX, CONTOUR, LFIXED, ISTATUS)
C	SUBROUTINE CALC_IMAGE_ELLIPSE(XCEN,YCEN,E,F,H,RAD_MULT)
C =========================================================================

C Based on ellipse formula: e*x^2 + f*y^2 + 2*h*x*y = 1

	SUBROUTINE CALC_CONTOUR_CENTER(XCEN,YCEN, E,F,H,FFRAC, A_MAX, IR,ISTATUS)
C
C Version used for centering spots
C
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /BACK_COM/ BACK(8000,2500)
C
C ISTATUS=1 if the routine fails.
C
	ISTATUS=1
C
	XCEN=X(IR)
	YCEN=Y(IR)
C
C RMAX is the maximum of the smoothed image within 4 pixels of XCEN,YCEN
	PEAK0=RMAX(IR)-BACK(NINT(XCEN),NINT(YCEN))
	IF(PEAK0 .LE. 0.0) RETURN
C
C If edges, scratches, etc. within a 20 pixel radius, return with ISTATUS=1
C
	CALL REJECT_RADIUS(IOUT, XCEN,YCEN,20.0)
	IF(IOUT .NE. 0) RETURN
C
C Use radius 5 for initial search area, with moveable center, and weights
C from 0 to 100% of peak height. Exit if silly ellipse shape or off-center.
C
	CALL CALC_AXES2EFH(E,F,H, 5.0,5.0,0.0)
	CALL CALC_ELL_MOMENTS(DXCEN,DYCEN, A,B,ANGLE,
	1					 XCEN,YCEN, E,F,H, 0.0,PEAK0, .FALSE.)
	IF(A.EQ.0.0 .OR. A.GT.50.0 .OR. B.LT.1.0) RETURN
	IF(DXCEN**2 + DYCEN**2 .GT. 25.0) RETURN
C
	XCEN=XCEN+DXCEN
	YCEN=YCEN+DYCEN
	A=MIN(10.0,2.0*A)
	B=MAX( 4.0,2.0*B)
	CALL CALC_AXES2EFH(E,F,H, A,B,ANGLE)
C
	CALL EDGE_ELLIPSE(IOUT, XCEN,YCEN,E,F,H,1.5)
	IF(IOUT .NE. 0) RETURN
C
	DO ITER=1,5
	  A0=A
C
	  CALL CALC_ELL_MOMENTS(DXCEN,DYCEN, A,B,ANGLE,
	1						XCEN,YCEN, E/2.0,F/2.0,H/2.0,
	2						0.2*PEAK0,0.24*PEAK0, .FALSE.)
	  IF(A.EQ.0.0 .OR. A.GT.50.0 .OR. B.LT.1.0) RETURN
C
	  XCEN=XCEN+DXCEN
	  YCEN=YCEN+DYCEN
	  IF((XCEN-X(IR))**2 + (YCEN-Y(IR))**2 .GT. 25.0) RETURN
C
	  AREA=A*B
	  A=MIN(A_MAX,A0*1.5, A)
	  B=MAX(2.0,A/5.0, B)
	  CALL CALC_AXES2EFH(E,F,H, A,B,ANGLE)
C
CCC	  FACT=2.0
	  FACT=2.0*A*B/AREA
	  FFRAC=GET_ELLIPSE_FILL(XCEN,YCEN,E/2.0,F/2.0,H/2.0, 0.22*PEAK0)*FACT
	  IF(FFRAC.LT.0.5 .OR. FFRAC .GT. 1.2) RETURN
C
	  IF(ABS(A-A0) .LT. 0.01) EXIT
C
	ENDDO
C
	ISTATUS=0
	RETURN
	END


	SUBROUTINE CALC_CONTOUR_SHAPE(E,F,H,FFRAC,A,B,ANGLE, DXCEN,DYCEN,
	1										SHKLM, CONTOUR,FFILL,A_MAX, IR,ISTATUS)
C
	CHARACTER SHKLM*(*)
C
C Version used for getting the ellipse for a centered spot
C
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /BACK_COM/ BACK(8000,2500)
	COMMON /DEBUG_COM/ IDEBUG
C
	XCEN=X(IR)
	YCEN=Y(IR)
C
C Set ISTATUS=1 for failure
C
	ISTATUS=1
C
C If edges, scratches, etc. within a 20 pixel radius, return with ISTATUS=1
C
	CALL REJECT_RADIUS(IOUT, XCEN,YCEN,20.0)
	IF(IOUT .NE. 0) THEN
	  WRITE(30,'(2A)') SHKLM,' rejected as ellipse overlaps edges, holes, etc.'
	  RETURN
	ENDIF
C
C RMAX is the maximum of the smoothed image within 4 pixels of XCEN,YCEN
C
	PEAK0=RMAX(IR)-BACK(NINT(XCEN),NINT(YCEN))
	IF(IDEBUG .EQ. IR) WRITE(30,'(/,A,I6,2F8.1)'),'*** DEBUG: IR,RMAX,PEAK0 = ',IR,RMAX(IR),PEAK0
	IF(PEAK0 .LE. 0.0) RETURN
C
C Use radius 5 circle for initial summation area
C
	CALL CALC_AXES2EFH(ESUM,FSUM,HSUM, 5.0,5.0,0.0)
C
C Calc. momments with fixed center, and weights from 0 to 100% of peak height
C
	CALL CALC_ELL_MOMENTS(DXCEN,DYCEN, A,B,ANGLE, XCEN,YCEN, ESUM,FSUM,HSUM, 0.0,PEAK0, .FALSE.)
	IF(IDEBUG .EQ. IR) WRITE(30,'(A,4F8.3,F8.0)'),'*** DEBUG: DXCEN,DYCEN,A,B,ANGLE = ',DXCEN,DYCEN,A,B,ANGLE*57.3
C  Return for bad ellipse shape or far off center
	IF(A.LT.1.0 .OR. B.LT.1.0 .OR. A.GT.50.0) THEN
	  WRITE(30,'(2A,2F5.1)',IOSTAT=IDUMMY) SHKLM,' rejected as  A,B =',A,B
	  RETURN
	ENDIF
	SHIFT=SQRT( DXCEN**2 + DYCEN**2 )
	IF(SHIFT .GT. 5.0) THEN
	  WRITE(30,'(2A,F5.1,A)',IOSTAT=IDUMMY) SHKLM,' rejected as off centre by',SHIFT,' pixels'
	  RETURN
	ENDIF
C
C Create new summation ellipse approx. twice the fitted ellipse area
C
	ASUM=MIN(10.0, A)*1.4
	BSUM=MAX(2.0,A/5.0, B)*1.4
	IF(IDEBUG .EQ. IR) WRITE(30,'(A,3F8.3)'),'*** DEBUG: ASUM,BSUM = ',ASUM,BSUM
	CALL CALC_AXES2EFH(ESUM,FSUM,HSUM, ASUM,BSUM,ANGLE)
C Return if summation ellipse crosses image boundaries
	CALL EDGE_ELLIPSE(IOUT, XCEN,YCEN,ESUM,FSUM,HSUM,1.0)
	IF(IOUT .NE. 0) RETURN
C
C Set cutoff for smoothed image peak
C
	CUTOFF=CONTOUR*PEAK0
C
	DO ITER=1,5
	  A_LAST=A
C
C Calculate moments inside the summation ellipse using 20% cutoff window and fixed spot centers
	  CALL CALC_ELL_MOMENTS(DXCEN,DYCEN, A,B,ANGLE, XCEN,YCEN, ESUM,FSUM,HSUM, CUTOFF*0.9,CUTOFF*1.1, .TRUE.)
	  IF(IDEBUG .EQ. IR) WRITE(30,'(A,I3,2F7.1)'),'*** DEBUG: ITER,CUTOFF = ',ITER,CUTOFF*0.9,CUTOFF*1.1
	  IF(IDEBUG .EQ. IR) WRITE(30,'(A,4F8.3,F8.0)'),'*** DEBUG: DXCEN,DYCEN,A,B,ANGLE = ',DXCEN,DYCEN,A,B,ANGLE*57.3
C Reject ridiculous values of A,B or a large peak shift
	  IF(A.LT.1.0 .OR. B.LT.1.0 .OR. A.GT.50.0) THEN
	    WRITE(30,'(2A,2F5.1)',IOSTAT=IDUMMY) SHKLM,' rejected as  A,B =',A,B
	    RETURN
	  ENDIF
	  SHIFT=SQRT( DXCEN**2 + DYCEN**2 )
	  IF(SHIFT .GT. 5.0) THEN
	    WRITE(30,'(2A,F5.1,A)',IOSTAT=IDUMMY) SHKLM,' rejected as off centre by',SHIFT,' pixels'
	    RETURN
	  ENDIF
C
C In debug mode, output current fractional fill value
	  IF(IDEBUG .EQ. IR) THEN
		  CALL CALC_AXES2EFH(E,F,H, A,B,ANGLE)
	    FFRAC=GET_ELLIPSE_FILL(XCEN,YCEN,E,F,H, CUTOFF)
			WRITE(30,'(A,2F8.3)'),'*** DEBUG: FFRAC = ',FFRAC
	  ENDIF
C
C Create new summation ellipse approx. twice the area of fitted ellipse
	  ASUM=MIN(A_MAX, A)*1.414
	  BSUM=MAX(2.0,A/5.0, B)*1.414
	  CALL CALC_AXES2EFH(ESUM,FSUM,HSUM, ASUM,BSUM,ANGLE)
	  IF(IDEBUG .EQ. IR) WRITE(30,'(A,2F8.3)'),'*** DEBUG: ASUM,BSUM = ',ASUM,BSUM
C
C Exit do loop if A value has converged
	  IF(ABS(A-A_LAST) .LT. 0.01) EXIT
	ENDDO
C
C Calculate final values of E,F,H,FFRAC
C
	CALL CALC_AXES2EFH(E,F,H, A,B,ANGLE)
	FFRAC=GET_ELLIPSE_FILL(XCEN,YCEN,E,F,H, CUTOFF)
C Reject if a bad fill fraction
	IF(FFRAC .LT. FFILL) THEN
	  WRITE(30,'(2A,F5.2)') SHKLM,' rejected as fill fraction =',FFRAC
	  RETURN
	ENDIF

	ISTATUS=0
	RETURN
	END



	SUBROUTINE CALC_ELL_MOMENTS(DXCEN,DYCEN,A,B,ANGLE, XCEN,YCEN, E,F,H,
	1								RLEVEL_LO,RLEVEL_HI, LFIXED)
C
	LOGICAL LFIXED
C
	COMMON /SMOOTH_COM/ SMOOTH_IMAGE(8000,2500)
C
C Do weighted sums of 1st and 2nd moments
C
	RMASS=0.0
	DXCEN=0.0
	DYCEN=0.0
	T11=0.0
	T22=0.0
	T12=0.0
C
	XWID=SQRT(F/(E*F-H**2))
	DO IX2=CEILING(XCEN-XWID),FLOOR(XCEN+XWID)
	  YWID=SQRT(MAX(0.0, (1.0 - ((IX2-XCEN)/XWID)**2 )/F ))
	  Y0=YCEN-(IX2-XCEN)*H/F
	  DO IY2=CEILING(Y0-YWID),FLOOR(Y0+YWID)
C
C Use the background-subtracted then smoothed image for the intensities
	    COUNTS=SMOOTH_IMAGE(IX2,IY2)
C
C Calculate the weight 0 to 1 for COUNTS = RLEVEL to RLEVEL+WINDOW
	    IF(COUNTS .LT. RLEVEL_LO) CYCLE
	    WEIGHT=MIN(1.0, (COUNTS-RLEVEL_LO)/(RLEVEL_HI-RLEVEL_LO) )
C Do weighted sums for X,Y relative to XCEN,YCEN
	    RMASS=RMASS+WEIGHT
	    DX=IX2-XCEN
	    DY=IY2-YCEN
	    DXCEN=DXCEN+WEIGHT*DX
	    DYCEN=DYCEN+WEIGHT*DY
	    T11=T11+WEIGHT*DX**2
	    T22=T22+WEIGHT*DY**2
	    T12=T12+WEIGHT*DX*DY
C
	  ENDDO
	ENDDO
C
C Calculate moments
C
	RMASS=MAX(0.1,RMASS)
	DXCEN=DXCEN/RMASS
	DYCEN=DYCEN/RMASS
	T11=T11/RMASS
	T22=T22/RMASS
	T12=T12/RMASS
C
C If spot center is not fixed, correct ellipse for X/Y offsets
C
	IF( .NOT.LFIXED ) THEN
	  T11=T11-DXCEN**2
	  T22=T22-DYCEN**2
	  T12=T12-DXCEN*DYCEN
	ENDIF
C
C Calculate the major and minor axes and tilt of the ellipse
C
	CALL CALC_COG2AXES(A,B,ANGLE, T11,T22,T12)
	CALL FIX_ELLI_AXES(A,B,ANGLE)
C
	RETURN
	END


	FUNCTION GET_ELLIPSE_FILL(XCEN,YCEN,E,F,H, CUTOFF)
C
C Return the maximum intensity in the specified ellipse
C Returns -1 if the ellipse crosses any edges
C
	COMMON /SMOOTH_COM/ SMOOTH_IMAGE(8000,2500)
C
C Reject ellipses with non-physical parameters
C
	GET_ELLIPSE_FILL=-1.0
	IF(E.LE.0.0 .OR. F.LE.0.0 .OR. E*F.LE.H**2) RETURN
C
C Calculate tangents of ellipse
C
	DX=SQRT( F/(E*F-H**2) )
	DY=SQRT( E/(E*F-H**2) )
C
C Return with VALUE OF -1 if edge is crossed
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
	IF(NINT(XCEN-DX).LT.1 .OR. NINT(XCEN+DX).GT.NUMX) RETURN
	IF(NINT(YCEN-DY).LT.1 .OR. NINT(YCEN+DY).GT.NUMY) RETURN
C
C Count pixels in ellipse (ISUM) and how many above CUTOFF (IFFILL)
C Use the background-subtracted then smoothed image for the intensities
C
	ISUM=0
	IFILL=0
	DO IX=CEILING(XCEN-DX),FLOOR(XCEN+DX)
	  DY=SQRT(MAX(0.0, (1.0 - ((IX-XCEN)/DX)**2 )/F ))
	  Y0=YCEN-(IX-XCEN)*H/F
	  DO IY=CEILING(Y0-DY),FLOOR(Y0+DY)
	    COUNTS=SMOOTH_IMAGE(IX,IY)
	    ISUM=ISUM+1
	    IF(COUNTS .GT. CUTOFF) IFILL=IFILL+1
	  ENDDO
	ENDDO
C
	GET_ELLIPSE_FILL=IFILL/(ISUM+0.01)
	RETURN
      END


	SUBROUTINE SUM_ELLIPSE_INTS(CSUM,NSUM, XCEN,YCEN,E,F,H, RAD,NRAD)
C
C Sum the counts within NRAD ellipses and radius multipliers RAD().
C Return the total counts in CSUM() and the number of pixels in NSUM().
C
C The sums are cumulative, RAD(n) contains all of RAD(1 to n-1),
C and it is assumed that the values in RAD() increase with index.
C
C No check is made for the validity of pixels used in the sums.
C
	INTEGER NSUM(NRAD)
	REAL CSUM(NRAD),RAD(NRAD)
C
	COMMON /RIMAGE_COM/ RIMAGE(8000,2500)
C
C Quick sanity check
C
	DO I=1,NRAD-1
	  IF(RAD(I) .GT. RAD(I+1)) CALL QUIT("BUG(sum_ellipse_ints): Invalid RAD")
	ENDDO
C
C Calculate tangents of ellipse
C
	TX=SQRT(F/(E*F-H**2))
	TY=SQRT(E/(E*F-H**2))
C
C Zero the sum arrays
C
	DO I=1,NRAD
	  NSUM(I)=0
	  CSUM(I)=0.0
	ENDDO
C
C Loop over pixels for last ellipse.
C Sum values in the first ellipse that contains the pixel.
C
	DX=TX*RAD(NRAD)
	DO IX=CEILING(XCEN-DX),FLOOR(XCEN+DX)
	  DY=SQRT(MAX(0.0, (1.0 - ((IX-XCEN)/DX)**2 )/F ))*RAD(NRAD)
	  Y0=YCEN-(IX-XCEN)*H/F
	  DO IY=CEILING(Y0-DY),FLOOR(Y0+DY)
C
	    RAD2=E*(IX-XCEN)**2 + F*(IY-YCEN)**2 + 2*H*(IX-XCEN)*(IY-YCEN)
	    DO IRAD=1,NRAD-1
	      IF(RAD2 .LT. RAD(IRAD)**2) EXIT
	    ENDDO
	    NSUM(IRAD)=NSUM(IRAD)+1
	    CSUM(IRAD)=CSUM(IRAD)+RIMAGE(IX,IY)
C
	  ENDDO
	ENDDO
C
C Convert to cumulative sums
C
	DO I=2,NRAD
	  NSUM(I)=NSUM(I)+NSUM(I-1)
	  CSUM(I)=CSUM(I)+CSUM(I-1)
	ENDDO
C
	RETURN
	END


	FUNCTION GET_ELLIPSE_MAXINT(XCEN,YCEN,E,F,H)
C
C Return the maximum intensity in the specified ellipse
C No check is made for the validity of pixels
C
	COMMON /RIMAGE_COM/ RIMAGE(8000,2500)
C
	RMAXINT=0.0
	DX=SQRT(F/(E*F-H**2))
	DO IX=CEILING(XCEN-DX),FLOOR(XCEN+DX)
	  DY=SQRT(MAX(0.0, (1.0 - ((IX-XCEN)/DX)**2 )/F ))
	  Y0=YCEN-(IX-XCEN)*H/F
	  DO IY=CEILING(Y0-DY),FLOOR(Y0+DY)
	    RMAXINT=MAX(RMAXINT,RIMAGE(IX,IY))
	  ENDDO
	ENDDO
C
	GET_ELLIPSE_MAXINT=RMAXINT
	RETURN
	END


	SUBROUTINE AVE_ELLIPSE_BKGS(BKG_AVE, XCEN,YCEN,E,F,H, RAD,NRAD)
C
C Averages the global background within NRAD ellipses with radius
C multipliers RAD(), exactly the same as for SUM_ELLIPSE_INTS().
C The averages are cumulative, RAD(n) contains all of RAD(1 to n-1).
C
	REAL BKG_AVE(NRAD),RAD(NRAD)
C
	COMMON /BACK_COM/ BACK(8000,2500)
C
	INTEGER NSUM(10)
	REAL CSUM(10)
C
C Sanity checks
C
	IF(NRAD .GT. 10) CALL QUIT("BUG(ave_ellipse_bkgs): NRAD > 10")
	DO I=1,NRAD-1
	  IF(RAD(I) .GT. RAD(I+1)) CALL QUIT("BUG(ave_ellipse_bkgs): Invalid RAD")
	ENDDO
C
C Calculate tangents of ellipse
C
	TX=SQRT(F/(E*F-H**2))
	TY=SQRT(E/(E*F-H**2))
C
C Zero the sum arrays
C
	DO I=1,NRAD
	  NSUM(I)=0
	  CSUM(I)=0.0
	ENDDO
C
C Loop over pixels for last ellipse.
C Sum values in the first ellipse that contains the pixel.
C
	DX=TX*RAD(NRAD)
	DO IX=CEILING(XCEN-DX),FLOOR(XCEN+DX)
	  DY=SQRT(MAX(0.0, (1.0 - ((IX-XCEN)/DX)**2 )/F ))*RAD(NRAD)
	  Y0=YCEN-(IX-XCEN)*H/F
	  DO IY=CEILING(Y0-DY),FLOOR(Y0+DY)
C
	    RAD2=E*(IX-XCEN)**2 + F*(IY-YCEN)**2 + 2*H*(IX-XCEN)*(IY-YCEN)
	    DO IRAD=1,NRAD-1
	      IF(RAD2 .LT. RAD(IRAD)**2) EXIT
	    ENDDO
	    NSUM(IRAD)=NSUM(IRAD)+1
	    CSUM(IRAD)=CSUM(IRAD)+BACK(IX,IY)
C
	  ENDDO
	ENDDO
C
C Convert to cumulative sums
C
	DO I=2,NRAD
	  NSUM(I)=NSUM(I)+NSUM(I-1)
	  CSUM(I)=CSUM(I)+CSUM(I-1)
	ENDDO
C
C Calculate the averages
C
	DO I=1,NRAD
	  BKG_AVE(I)=CSUM(I)/MAX(1,NSUM(I))
	ENDDO
C
	RETURN
	END


C ======================= Some short ellipse handling routines ======================

	SUBROUTINE CALC_COG2AXES(A,B,ANGLE, T11,T22,T12)
C
CC EQUATION CHECKED NUMERICALLY:
	TEMP=SQRT( (T11-T22)**2 + 4.0*T12**2 )
	A=2*SQRT(MAX(1E-9, 2.0*(T11*T22-T12**2) / MAX(1E-9, T11+T22 - TEMP) ))
	B=2*SQRT(MAX(1E-9, 2.0*(T11*T22-T12**2) / MAX(1E-9, T11+T22 + TEMP) ))
	ANGLE=0.5*ATAN2( 2.0*T12 , T11-T22 )
C
	RETURN
	END


	SUBROUTINE CALC_EFH2AXES(A,B,ANGLE, E,F,H)
C
C Calculate the major and minor axes and tilt of the ellipse
C from the E,F,H ellipse parameters. "A" is returned as the
C major axis, and "ANGLE" is the rotation of the major axes
C anticlockwise from the +X direction, measured in radians.
C
CC EQUATION CHECKED NUMERICALLY:
	TEMP=SQRT( (E-F)**2 + 4.0*H**2 )
	A=1.0/SQRT(0.5*MAX(1E-9, E+F - TEMP ))
	B=1.0/SQRT(0.5*MAX(1E-9, E+F + TEMP ))
	ANGLE=0.5*ATAN2( -2.0*H , F-E )
C
	RETURN
	END


	SUBROUTINE FIX_ELLI_AXES(A,B,ANGLE)
C
	IF(A .LT. B) THEN
	  TEMP=A
	  A=B
	  B=TEMP
	  IF(ANGLE .GT. 0.0) THEN
	    ANGLE=ANGLE-1.5708
	  ELSE
	    ANGLE=ANGLE+1.5708
	  ENDIF
	ENDIF
C
	RETURN
	END


	SUBROUTINE CALC_AXES2EFH(E,F,H, A0,B0,ANGLE)
C
C Calculate the E,F,H ellipse parameters from the major,
C minor axes and rotation of the ellipse.
C Angle is the rotation of the "A" axes anticlockwise
C from the +X direction, measured in radians.
C
	A=MAX(1E-3,A0)
	B=MAX(1E-3,B0)
C
CC EQUATION CHECKED NUMERICALLY:
	E=(COS(ANGLE)/A)**2+(SIN(ANGLE)/B)**2
	F=(SIN(ANGLE)/A)**2+(COS(ANGLE)/B)**2
	H=((1.0/A)**2-(1.0/B)**2)*SIN(ANGLE)*COS(ANGLE)
C
C Prevent non-positive ellipses
C
	IF(E*F .LT. H**2) H=H*0.999999
C
	RETURN
	END


	FUNCTION CALC_ELLIPSE_DSQ(DX,DY,E,F,H)
C
	CALC_ELLIPSE_DSQ=E*DX**2 +F*DY**2 +2.0*H*DX*DY
	RETURN
      END


C ===================== UNUSED TEST ROUTINES ====================

	SUBROUTINE TEST_ELLIPSES()
C
	REAL CG(2)
C
	PI=4.0*ATAN(1.0)
	RADS=180.0/PI
C
C Check CALC_ELLIPSE_DSQ() agrees with CALC_AXES2EFH(), and that ANGLE is the
C rotation in radians of the "A" axes anticlockwise from the +X direction.
C
	IF( .FALSE. ) THEN
	PRINT '(/,1X,A)','Start Test 0'
	A0=2.0
	B0=1.0
	ANGLE0=110
C
	CALL CALC_AXES2EFH(E,F,H, A0,B0,ANGLE0/RADS)
C
	DSQ_MAX=-1E6
	DSQ_MIN=1E6
	DO X=-3.0,3.0,0.001
	  DO Y=-3.0,3.0,0.001
	    VAL=CALC_ELLIPSE_DSQ(X,Y,E,F,H)
	    DSQ=X**2 + Y**2
	    IF(VAL .LE. 1.0) THEN
	      IF(DSQ .GT. DSQ_MAX) THEN
	        DSQ_MAX=DSQ
	        X_MAX=X
			Y_MAX=Y
	      ENDIF
	    ELSE
	      IF(DSQ .LT. DSQ_MIN) THEN
	        DSQ_MIN=DSQ
	        X_MIN=X
			Y_MIN=Y
	      ENDIF
		ENDIF
	  ENDDO
	ENDDO
C
	PRINT *,A0,B0,ANGLE0
	PRINT *,SQRT(DSQ_MAX),SQRT(DSQ_MIN),ATAN2D(Y_MAX,X_MAX)
	PRINT *,'End Test 0'
	ENDIF
C
C Check CALC_COG2AXES() gives the correct A,B,ANGLE from the moments summation
C
	IF( .FALSE. ) THEN
	PRINT '(/,1X,A)','Start Test 1'
	A0=2.0
	B0=1.0
	ANGLE0=45
	XCEN=0.3
	YCEN=0.1
C
	CALL CALC_AXES2EFH(E,F,H, A0,B0,ANGLE0/RADS)
C
	RMASS=0.0
	CG(1)=0.0
	CG(2)=0.0
	T11=0.0
	T22=0.0
	T12=0.0
	DO X=-3.0,3.001,0.01
	  DO Y=-3.0,3.001,0.01
	    VAL=CALC_ELLIPSE_DSQ(X,Y,E,F,H)
	    IF(VAL .LE. 1.0) THEN
C Do weighted sums of X,Y relative to XCEN,YCEN
C Use unit weights for the test
	      WEIGHT=1.0
	      RMASS=RMASS+WEIGHT
	      DX=X-XCEN
	      DY=Y-YCEN
	      CG(1)=CG(1)+WEIGHT*DX
	      CG(2)=CG(2)+WEIGHT*DY
	      T11=T11+WEIGHT*DX**2
	      T22=T22+WEIGHT*DY**2
	      T12=T12+WEIGHT*DX*DY
	    ENDIF
	  ENDDO
	ENDDO
C Calculate moments
	CG(1)=CG(1)/RMASS
	CG(2)=CG(2)/RMASS
	T11=T11/RMASS
	T22=T22/RMASS
	T12=T12/RMASS
C If spot center is not fixed, correct ellipse for X/Y offsets
	PRINT '(2F7.2,2X,2F7.2)',XCEN,YCEN,CG
	T11=T11-CG(1)**2
	T22=T22-CG(2)**2
	T12=T12-CG(1)*CG(2)
C Calculate the major and minor axes and tilt of the ellipse
	CALL CALC_COG2AXES(A,B,ANGLE, T11,T22,T12)
	CALL FIX_ELLI_AXES(A,B,ANGLE)
C
	PRINT '(3F7.2,2X,3F7.3,2X,3F7.2)',A0,B0,ANGLE0,E,F,H,A,B,ANGLE*RADS
	PRINT *,'End Test 1'
	ENDIF
C
C Check CALC_AXES2EFH() and CALC_EFH2AXES() agree for A > B ellipses
C
	IF( .FALSE. ) THEN
	PRINT '(/,1X,A)','Start Test 2'
	A0=3.25
	B0=1.7956
	DO ANGLE0=0.0,180.0,10.0
	  CALL CALC_AXES2EFH(E,F,H, A0,B0,ANGLE0/RADS)
	  CALL CALC_EFH2AXES(A,B,ANGLE, E,F,H)
	  ANGLE=ANGLE*RADS
	  PRINT '(3F7.2,2X,3F7.3,2X,3F7.2)',A0,B0,ANGLE0,E,F,H,A,B,ANGLE
	ENDDO
	PRINT *,'End Test 2'
	ENDIF
C
C Check CALC_AXES2EFH() and CALC_EFH2AXES() agree for A < B ellipses
C
	IF( .FALSE. ) THEN
	PRINT '(/,1X,A)','Start Test 3'
	A0=1.7956
	B0=3.25
	DO ANGLE0=0.0,180.0,10.0
	  CALL CALC_AXES2EFH(E,F,H, A0,B0,ANGLE0/RADS)
	  CALL CALC_EFH2AXES(A,B,ANGLE, E,F,H)
	  PRINT '(3F7.2,2X,3F7.3,2X,3F7.2)',A0,B0,ANGLE0,E,F,H,A,B,ANGLE*RADS
	ENDDO
	PRINT *,'End Test 3'
	ENDIF
C
C Check FIX_ELLI_AXES() is correct
C
	IF( .FALSE. ) THEN
	PRINT '(/,1X,A)','Start Test 4'
	A1=1.7956
	B1=3.25
	DO ANGLE1=0.0,180.0,10.0
	  A0=A1
	  B0=B1
	  ANGLE0=ANGLE1/RADS
	  CALL CALC_AXES2EFH(E,F,H, A0,B0,ANGLE0)
	  CALL FIX_ELLI_AXES(A0,B0,ANGLE0)
	  CALL CALC_AXES2EFH(E2,F2,H2, A0,B0,ANGLE0)
	  PRINT '(3F7.2,3X,3F7.2,3X,3F8.4)',A1,B1,ANGLE1,A0,B0,ANGLE0*RADS,E-E2,F-F2,H-H2
	ENDDO
	PRINT *,'End Test 4'
	ENDIF
C
	STOP 'End of TEST_ELLIPSES'
      END


	SUBROUTINE GET_ELLIPSE_TEST(XCEN,YCEN,E,F,H,RMAX,CONTOUR)
C
	COMMON /BACK_COM/ BACK(8000,2500)
	COMMON /SMOOTH_COM/ SMOOTH_IMAGE(8000,2500)
C
	IX=NINT(XCEN)
	IY=NINT(YCEN)
C
	BG=BACK(IX,IY)
cc	RLEVEL_LO=CONTOUR*(RMAX-BG) + BG
	RLEVEL_LO=CONTOUR*(RMAX-BG)
	RLEVEL_HI=RLEVEL_LO+SQRT(BG)
c
	do iy2=iy-9,iy+9
	  print '(19i4)',(nint(smooth_image(k,iy2)),k=ix-9,ix+9)
	enddo
C
	ISUM=0
	IFILL=0
	DX=SQRT(F/(E*F-H**2))
	DO IX=CEILING(XCEN-DX),FLOOR(XCEN+DX)
	  DY=SQRT(MAX(0.0, (1.0 - ((IX-XCEN)/DX)**2 )/F ))
	  Y0=YCEN-(IX-XCEN)*H/F
	  DO IY=CEILING(Y0-DY),FLOOR(Y0+DY)
CCC Use the smoothed and background-subtracted image for the intensities
cc	    COUNTS=IMAGE(IX,IY)
CC	    PRINT '(2I5,3F10.1)',IX,IY,COUNTS-bg,RLEVEL_LO-bg,RLEVEL_HI-bg
	    COUNTS=SMOOTH_IMAGE(IX,IY)
	    PRINT '(2I5,3F10.1)',IX,IY,COUNTS,RLEVEL_LO,RLEVEL_HI
	    ISUM=ISUM+1
	    IF(COUNTS .GT. RLEVEL_LO) IFILL=IFILL+1
	  ENDDO
	ENDDO
	print *,ifill,isum
C
	RETURN
	END


	SUBROUTINE CALC_CONTOUR_ELLIPSE_OLD(E,F,H,CG,FFRAC,A,B,ANGLE,
	1							XCEN,YCEN,RMAX, CONTOUR, LFIXED, ISTATUS)
C
C Calculate the ellipse shape for the peak IR using the CONTOUR level.
C The routine iterates the total area used to search for the contour.
C This modification allows ellipses to be correctly found for large
C streaky spots.
C
C ISTATUS=1 if the routine fails.
C LFIXED is true when the ellipse are calculated assuming
C the spot center is fixed to the (XCEN,YCEN).
C
	LOGICAL LFIXED
C
	COMMON /BACK_COM/ BACK(8000,2500)
	COMMON /SMOOTH_COM/ SMOOTH_IMAGE(8000,2500)
C
	REAL CG(2)
C
	IX=NINT(XCEN)
	IY=NINT(YCEN)
C
C Setup a weighting scheme for pixels within 1 esd (no gain correction)
C above the nominal cutoff level
C	
	BG=BACK(IX,IY)
	RLEVEL_LO=CONTOUR*(RMAX-BG)
	RLEVEL_HI=RLEVEL_LO+SQRT(BG)
C
C Do 3 iterations first summing pixels within a radius of 9, later with
C a radius determined by the observed ellipse size
C
	IRAD=9
	DO ITER=1,3
C If problems integrating over radius=IRAD, give up and signal failure
	  CALL REJECT_RADIUS(IOUT, XCEN,YCEN,FLOAT(IRAD))
	  IF(IOUT .EQ. 1) THEN
	    ISTATUS=1
	    RETURN
	  ENDIF
C
	  RMASS=0.0
	  CG(1)=0.0
	  CG(2)=0.0
	  T11=0.0
	  T22=0.0
	  T12=0.0
	  DO IDY=-IRAD,IRAD
	    IRADX=IFIX(SQRT(FLOAT( IRAD**2 - IDY**2 )))
	    DO IDX=-IRADX,IRADX
C Use the smoothed and background-subtracted image for the intensities
	      COUNTS=SMOOTH_IMAGE(IX+IDX,IY+IDY)
C Calculate the weight 0 to 1 for COUNTS = RLEVEL to RLEVEL+WINDOW
	      IF(COUNTS .LT. RLEVEL_LO) CYCLE
	      WEIGHT=MIN(1.0, (COUNTS-RLEVEL_LO)/(RLEVEL_HI-RLEVEL_LO) )
C Do weighted sums of X,Y relative to XCEN,YCEN
	      RMASS=RMASS+WEIGHT
	      DX=IX+IDX-XCEN
	      DY=IY+IDY-YCEN
	      CG(1)=CG(1)+WEIGHT*DX
	      CG(2)=CG(2)+WEIGHT*DY
	      T11=T11+WEIGHT*DX**2
	      T22=T22+WEIGHT*DY**2
	      T12=T12+WEIGHT*DX*DY
	    ENDDO
	  ENDDO
C If nothing in ellipse(!?!), double IRAD and cycle loop
	  IF(RMASS .EQ. 0.0) THEN
	    IRAD=2*IRAD
	    CYCLE
	  ENDIF
C Calculate moments
	  CG(1)=CG(1)/RMASS
	  CG(2)=CG(2)/RMASS
	  T11=T11/RMASS
	  T22=T22/RMASS
	  T12=T12/RMASS
C If spot center is not fixed, correct ellipse for X/Y offsets
	  IF( .NOT.LFIXED ) THEN
	    T11=T11-CG(1)**2
	    T22=T22-CG(2)**2
	    T12=T12-CG(1)*CG(2)
	  ENDIF
C Calculate the major and minor axes and tilt of the ellipse
	  CALL CALC_COG2AXES(A,B,ANGLE, T11,T22,T12)
	  CALL FIX_ELLI_AXES(A,B,ANGLE)
C
C Increase IRAD to x1.5 the max. ellipse axis for the next iteration
C If IRAD is not increased, exit the loop early
	  IRAD2=NINT(1.5*A)
	  IF(IRAD2 .LE. IRAD) EXIT
	  IRAD=IRAD2
C
	ENDDO
C
C Calculate the E,F,H ellipse parameters from the major,
C minor axes and tilt of the ellipse.
C
	CALL CALC_AXES2EFH(E,F,H, A,B,ANGLE)
C
C Calculate contour filling fraction
C
	FFRAC=GET_ELLIPSE_FILL(XCEN,YCEN,E,F,H, RLEVEL_LO)
cc	FFRAC=GET_ELLIPSE_FILL(XCEN,YCEN,E,F,H, BG+RLEVEL_LO)
ccc	ffrac=(1.0+ffrac)/2.0
C
	ISTATUS=0
	RETURN
	END


	SUBROUTINE CALC_CONTOUR_ELLIPSE_OLDER(E,F,H,CG,FFRAC,A,B,ANGLE,
	1							XCEN,YCEN,RMAX, CONTOUR, LFIXED, ISTATUS)
C
C Calculate the ellipse shape for the peak IR using the CONTOUR level.
C The routine iterates the total area used to search for the contour.
C This modification allows ellipses to be correctly found for large
C streaky spots.
C
C LFIXED is true when the ellipse are calculated assuming
C the spot center is fixed to the (XCEN,YCEN).
C
C XCEN & YCEN are no longer modified by this routine.
C
	LOGICAL LFIXED
C
	COMMON /BACK_COM/ BACK(8000,2500)
	COMMON /SMOOTH_COM/ SMOOTH_IMAGE(8000,2500)
C
	REAL CG(2)
C
	IX=NINT(XCEN)
	IY=NINT(YCEN)
C
C Setup a weighting scheme for pixels within 1 esd (no gain correction)
C above the nominal cutoff level
C	
	BG=BACK(IX,IY)
	RLEVEL_LO=CONTOUR*(RMAX-BG)
	RLEVEL_HI=RLEVEL_LO+SQRT(BG)
C
C Do 3 iterations first summing pixels within a radius of 9, later with
C a radius determined by the observed ellipse size
C
	IRAD=9
	DO ITER=1,3
C If problems integrating over radius=IRAD, give up and signal failure
	  CALL REJECT_RADIUS(IOUT, XCEN,YCEN,FLOAT(IRAD))
	  IF(IOUT .EQ. 1) THEN
	    ISTATUS=1
	    RETURN
	  ENDIF
C
	  RMASS=0.0
	  CG(1)=0.0
	  CG(2)=0.0
	  T11=0.0
	  T22=0.0
	  T12=0.0
	  DO IDY=-IRAD,IRAD
	    IRADX=IFIX(SQRT(FLOAT( IRAD**2 - IDY**2 )))
	    DO IDX=-IRADX,IRADX
C Use the smoothed image for the intensities
	      COUNTS=SMOOTH_IMAGE(IX+IDX,IY+IDY)
C Calculate the weight 0 to 1 for COUNTS = RLEVEL to RLEVEL+WINDOW
	      IF(COUNTS .LT. RLEVEL_LO) CYCLE
	      WEIGHT=MIN(1.0, (COUNTS-RLEVEL_LO)/(RLEVEL_HI-RLEVEL_LO) )
C Do weighted sums of X,Y relative to XCEN,YCEN
	      RMASS=RMASS+WEIGHT
	      DX=IX+IDX-XCEN
	      DY=IY+IDY-YCEN
	      CG(1)=CG(1)+WEIGHT*DX
	      CG(2)=CG(2)+WEIGHT*DY
	      T11=T11+WEIGHT*DX**2
	      T22=T22+WEIGHT*DY**2
	      T12=T12+WEIGHT*DX*DY
	    ENDDO
	  ENDDO
C If nothing in ellipse(!?!), double IRAD and cycle loop
	  IF(RMASS .LE. 0.0) THEN
	    IRAD=2*IRAD
	    CYCLE
	  ENDIF
C Calculate moments
	  CG(1)=CG(1)/RMASS
	  CG(2)=CG(2)/RMASS
	  T11=T11/RMASS
	  T22=T22/RMASS
	  T12=T12/RMASS
C If spot center is not fixed, correct ellipse for X/Y offsets
	  IF( .NOT.LFIXED ) THEN
	    T11=T11-CG(1)**2
	    T22=T22-CG(2)**2
	    T12=T12-CG(1)*CG(2)
	  ENDIF
C Calculate the major and minor axes and tilt of the ellipse
	  CALL CALC_COG2AXES(A,B,ANGLE, T11,T22,T12)
	  CALL FIX_ELLI_AXES(A,B,ANGLE)
C
C Increase IRAD to x1.5 the max. ellipse axis for the next iteration
C If IRAD is not increased, exit the loop early
	  IRAD2=NINT(1.5*A)
	  IF(IRAD2 .LE. IRAD) EXIT
	  IRAD=IRAD2
C
	ENDDO
C
C Calculate the E,F,H ellipse parameters from the major,
C minor axes and tilt of the ellipse.
C
	CALL CALC_AXES2EFH(E,F,H, A,B,ANGLE)
C
C Calculate contour filling fraction
C
	FFRAC=RMASS/(3.1415926*A*B)
C
	ISTATUS=0
	RETURN
	END


	SUBROUTINE CALC_IMAGE_ELLIPSE(XCEN,YCEN,E,F,H,RAD_MULT)
C
C Store in /IFLY/ data used to later draw the ellipses
C
	REAL RAD_MULT(3)
C
	COMMON/IFLY/IFLYX(1000),IFLYY(1000),ICONT(1000),NFLY
C
	INTEGER IY_LO(2),IY_HI(2)
C
C Zero the list of points used to draw ellipses on the image file
C
	NFLY=0
C
C Loop over the three ellipses sizes
C
	TY=SQRT(E/(E*F-H**2))		! tangent lines (pixels)
	TX=SQRT(F/(E*F-H**2))		! tangent lines
	DO IRAD=1,3
C
	  ELLIM=RAD_MULT(IRAD)
	  DX=SQRT(ELLIM)*TX
	  DO IX=CEILING(XCEN-DX),FLOOR(XCEN+DX)
	    DY=SQRT(MAX(0.0, (1.0 - ((IX-XCEN)/DX)**2 )*ELLIM/F ))
	    Y0=YCEN-(IX-XCEN)*H/F
C Choose 1 or 2 sets of IY limits to search over
	    IF(IX.LE.CEILING(XCEN-DX)+1 .OR. IX.GE.FLOOR(XCEN+DX)-1) THEN
C If close to the IX limits, use 1 set of the complete range of IY
	      NSET=1
	      IY_LO(1)=CEILING(Y0-DY)
	      IY_HI(1)=FLOOR(Y0+DY)
	    ELSE
C Set 1 is 3 IY near the minimum, set 2 is 3 IY near the maximum
	      NSET=2
	      IY_LO(1)=CEILING(Y0-DY)
	      IY_HI(1)=IY_LO(1)+2
	      IY_HI(2)=FLOOR(Y0+DY)
	      IY_LO(2)=IY_HI(1)-2
	    ENDIF
C Loop IY over the 1 or 2 sets
	    DO I=1,NSET
	      DO IY=IY_LO(I),IY_HI(I)
	        DX2=IX-XCEN
	        DY2=IY-YCEN
	        DSQ0=E*DX2**2 +F*DY2**2 +2.0*H*DX2*DY2
	        DSQX=DSQ0+E+2.0*ABS(E*DX2+H*DY2)
	        DSQY=DSQ0+F+2.0*ABS(H*DX2+F*DY2)
C If within 1 in X or Y of the ellipse edge, store pixel values
	        IF(MAX(DSQX,DSQY) .GE. ELLIM) THEN
	          IF(NFLY .LT. 1000) NFLY=NFLY+1
	          IFLYX(NFLY)=IX
	          IFLYY(NFLY)=IY
	          ICONT(NFLY)=IRAD
	        ENDIF
	      ENDDO
	    ENDDO
C
	  ENDDO
C
	ENDDO

  	RETURN
	END
