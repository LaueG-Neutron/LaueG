	SUBROUTINE CALC_PEAK_VECTORS(NPEAKS)
C
C Calculate the scattering vectors for all observed spots.
C The list of vectors is stored in /PEAKS_COM/.
C
	PARAMETER (NSPOTS_MAX=10000)
	COMMON /SPOTS_COM/ NSPOTS,XSPOT(NSPOTS_MAX),YSPOT(NSPOTS_MAX),ESDS(NSPOTS_MAX)
C
	PARAMETER (NPEAKS_MAX=10000)
	REAL H_PEAKS
	COMMON /PEAKS_COM/ H_PEAKS(3,NPEAKS_MAX)
C
C Calculate unit scattering vectors from XSPOT,YSPOT values
C NB: Calculate vectors for all spots as we need them to validate solutions
C
	DO I=1,NPEAKS
	  CALL CALC_PIXEL_TO_HVEC(H_PEAKS(1,I), XSPOT(I),YSPOT(I))
	  CALL CALC_HVEC_TO_PIXEL(X,Y, H_PEAKS(1,I))
	ENDDO
C
	RETURN
	END



	SUBROUTINE CALC_PIXEL_TO_HVEC(HSPOT, XSPOT,YSPOT)
C
C Calculate the "unit" scattering vector for the given pixel X & Y
C
	REAL HSPOT(3)
C
C To begin with use hardwired info, later read it from the input file
C
	COMMON /INSTRUM_COM/ ITYPE,NUMXY(2),XY_CEN(2),XY_SIZE(2),DRUM_RAD
C
	REAL S0(3),SH(3),ROT_MX(3,3)
C
C Cylindrical IP
	IF(ITYPE .NE. 3) THEN
	  ANGLE=(XSPOT-XY_CEN(1))*XY_SIZE(1)/DRUM_RAD
	  X=DRUM_RAD*SIN(ANGLE)
	  Y=(YSPOT-XY_CEN(2))*XY_SIZE(2)
	  Z=DRUM_RAD*COS(ANGLE)
	  CALL VC3_SET(SH, X,Y,Z)
C Octagonal CCD
	ELSE
	  IPLATE=1+INT( XSPOT/(NUMXY(1)/8.0) )
	  YROT=(IPLATE-5)*45.0
	  X0=XY_CEN(1)+(IPLATE-5)*NUMXY(1)/8.0
	  DX=(XSPOT-X0)*XY_SIZE(1)
	  DY=(YSPOT-XY_CEN(2))*XY_SIZE(2)
	  CALL VC3_SET(SH, DX,DY,DRUM_RAD)
	  CALL CALC_YROT_MX(ROT_MX, -YROT)
	  CALL VC3_MULT(SH, ROT_MX,SH)
	ENDIF
C
	CALL VC3_UNIT(SH)
	CALL VC3_SET(S0, 0.0,0.0,1.0)
	CALL VC3_SUB(HSPOT, SH,S0)
	CALL VC3_UNIT(HSPOT)
C
	RETURN
	END



	SUBROUTINE CALC_HVEC_TO_PIXEL(XSPOT,YSPOT, HSPOT)
C
C Calculate the pixel X & Y for the given scattering vector
C
	REAL HSPOT(3)
C
C To begin with use hardwired info, later read it from the input file
C
	COMMON /INSTRUM_COM/ ITYPE,NUMXY(2),XY_CEN(2),XY_SIZE(2),DRUM_RAD
C
C
	REAL SH(3),ROT_MX(3,3)
C
	WAV = -(HSPOT(1)**2 +  HSPOT(2)**2 + HSPOT(3)**2) / (2.0*HSPOT(3))
C
	SH(1)=HSPOT(1)/WAV
	SH(2)=HSPOT(2)/WAV
	SH(3)=HSPOT(3)/WAV+1.0
C
C
C
C Cylindrical IP
	IF(ITYPE .NE. 3) THEN
	  ANGLE=ATAN2(SH(1),SH(3))
	  DX=ANGLE*DRUM_RAD
	  DY=SH(2)/SQRT(SH(1)**2+SH(3)**2)*DRUM_RAD
	  XSPOT=XY_CEN(1)+DX/XY_SIZE(1)
	  YSPOT=XY_CEN(2)+DY/XY_SIZE(2)
C Octagonal CCD
	ELSE
	  YROT=ATAN2D(SH(1),SH(3))
	  IPLATE=5+NINT( YROT/45.0 )
	  IF(IPLATE .EQ. 9) IPLATE=1
C
	  YROT=(IPLATE-5)*45.0
	  CALL CALC_YROT_MX(ROT_MX, YROT)
	  CALL VC3_MULT(SH, ROT_MX,SH)
C
	  DX=SH(1)/SH(3)*DRUM_RAD
	  DY=SH(2)/SH(3)*DRUM_RAD
	  X0=XY_CEN(1)+(IPLATE-5)*NUMXY(1)/8.0
	  XSPOT=X0+DX/XY_SIZE(1)
	  YSPOT=XY_CEN(2)+DY/XY_SIZE(2)
	ENDIF
C
	RETURN
	END
