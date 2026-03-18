C ================ Routines in this file =========================
C	SUBROUTINE MAKE_GLOBAL_BKG(PIX_CEN,PIX_SIZE,DRUM_RAD0,IPIXCEN0)
C	FUNCTION GET_GLOBAL_BKG(X,Y)
C	FUNCTION GET_GLOBAL_BKG_ESD(X,Y)
C	SUBROUTINE CALC_SPOTS_BKG(NREFS)
C ================================================================

	SUBROUTINE MAKE_GLOBAL_BKG(PIX_CEN,PIX_SIZE,DRUM_RAD0,IPIXCEN0)
C
C Global background means the background at an arbitrary pixel position
C
	REAL PIX_CEN(2),PIX_SIZE(2)
C
C Calculate background array in /BACK_COM/ from image in /RIMAGE_COM/
C
	COMMON /RIMAGE_COM/ RIMAGE(8000,2500)
	COMMON /BACK_COM/ BACK(8000,2500)
	COMMON /BACK_GEOM_COM/ NUMX,NUMY,IXCEN,IYCEN,XSIZE,YSIZE,DRUM_RAD,IPIXCEN
C
	INTEGER IGRIDX(10),IGRIDY(10)
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
C
	IXCEN=NINT(PIX_CEN(1))
	IYCEN=NINT(PIX_CEN(2))
	XSIZE=PIX_SIZE(1)
	YSIZE=PIX_SIZE(2)
	IPIXCEN=IPIXCEN0
	DRUM_RAD=DRUM_RAD0
C
C Calculate background while redirecting PRINT output to unit=30
C
	CALL CALC_BACK_ARRAY2(BACK, RIMAGE, 30)
C
C Output global background calculated on a grid of 10 x 10 across X & Y
C
	DO I=1,10
		IGRIDX(I)=NINT( NUMX*(0.1*I-0.05) )
		IGRIDY(I)=NINT( NUMY*(0.1*I-0.05) )
	ENDDO
C
	WRITE(30,'(/,A)') 'Calculated global background for X(horiz) by Y(vert)'
	WRITE(30,'(5X,10I6)',IOSTAT=IDUMMY) (IGRIDX(KX),KX=1,10)
	DO IY=1,10
	  WRITE(30,'(I5,21I6)',IOSTAT=IDUMMY) IGRIDY(IY),
	1							(NINT(BACK(IGRIDX(KX),IGRIDY(IY))), KX=1,10)
	ENDDO
C
	RETURN
	END


	FUNCTION GET_GLOBAL_BKG(X,Y)
C
C Global background means the background at an arbitrary pixel position
C
	COMMON /BACK_COM/ BACK(8000,2500)
C
	GET_GLOBAL_BKG=BACK(NINT(X),NINT(Y))
	RETURN
	END


	FUNCTION GET_GLOBAL_BKG_ESD(X,Y)
C
C Global background means the background at an arbitrary pixel position
C
	COMMON /BACK_COM/ BACK(8000,2500)
	COMMON /GAIN_COM/ GAIN,ESD_FACTOR,CROSSTALK,CROSS_ELLIPSE,CROSS_SMOOTH
C
C Factor to reduce "raw" sigma
C Estimated using BACK_ESDS.EXE simulation
C
	FACTOR=70.0
C
	BKG=BACK(NINT(X),NINT(Y))
	GET_GLOBAL_BKG_ESD=ESD_FACTOR*SQRT(MAX(1.0,BKG))/FACTOR
	RETURN
	END


	SUBROUTINE CALC_SPOTS_BKG(NREFS)
C
C Calculate, and store in COMMON, the background at the centre of each spot
C
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /BGLOCAL_COM/ BKG(20000),BKG_ESD(20000),BKG_VAR(20000)
C
	DO IR=1,NREFS
	  BKG(IR)=GET_GLOBAL_BKG(X(IR),Y(IR))
	  BKG_ESD(IR)=GET_GLOBAL_BKG_ESD(X(IR),Y(IR))
	ENDDO
C
	RETURN
	END
