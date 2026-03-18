C ========================== Routines in this file ============================
C	SUBROUTINE ARGBOX_SETUP_PARS(BASE_NAME, RAD_MULT,CONTOUR,AREA_FAC,TOLER_NBOUR,
C	1			FFILL,LCENTER, NMODEL_CEN,NMODEL_OUT, SCALE,A_MAX,R_SMOOTH, IPIXCEN)
C	SUBROUTINE ARGBOX_SETUP_DATA(ZCUTOFFS, BASE_NAME, NUMX,NUMY,NREFS, IOPT_HKL, IMAGE_ZERO)
C	SUBROUTINE SETUP_GAIN()
C	SUBROUTINE SETUP_BOXES_REFS(NREFS)
C	SUBROUTINE SETUP_SPOT_AVE(NREFS)
C ============= (Image - Background) Smoothing ==================================
C	SUBROUTINE CALC_SMOOTH_IMAGE()
C	SUBROUTINE SMOOTH_IMAGE()
C	SUBROUTINE SMOOTH_IMAGE2(RIMAGE, NUMX,NUMY,ISMOOTH)
C ============ Setup model spots in zones =====================
C	SUBROUTINE LOAD_IZONES(NREFS)
C	SUBROUTINE CALC_ZONE_CUTOFFS(ZCUTOFFS,NWANT,NREFS)
C =============================================================================

	SUBROUTINE ARGBOX_SETUP_PARS(BASE_NAME, RAD_MULT,CONTOUR,AREA_FAC,TOLER_NBOUR,
	1				FFILL,LCENTER, NMODEL_CEN,NMODEL_OUT, SCALE,A_MAX,R_SMOOTH, IPIXCEN)
C
	LOGICAL LCENTER
	CHARACTER BASE_NAME*80
	REAL RAD_MULT(3)
C
	COMMON /OPTIONS_COM/ IOPT_HKL,NHKLM
	COMMON /GAIN_COM/ GAIN,ESD_FACTOR,CROSSTALK,CROSS_ELLIPSE,CROSS_SMOOTH
	COMMON /MODELS_PFRAC_COM/ PFRAC_TARGET,PFRAC_ERROR,DPFRAC_ZERO,DPFRAC_GRAD
C
C Some hardwired parameters
C
	SCALE=10.0				! downscale for output intensities
	A_MAX=10.0
	R_SMOOTH=10
C
C DATA AND PARAMETER INPUT
C
C Read the parameter part of the input file ___laueg_argonne_boxes.in
C NB: Opens the *.lis file and writes header, acknowledgements, version
C
	CALL READ_PARS_FILE(BASE_NAME, CONTOUR,AREA_FAC,TOLER_NBOUR,FFILL,
	1				LCENTER, NMODEL_CEN,NMODEL_OUT, IPIXCEN)
C
C Log the integration parameters read from the parameter file
C
	WRITE(30,'(//,A)') '================ Parameter Input ==============='
	WRITE(30,'(/,A)') 'Reading parameters from the LaueG input file'
	WRITE(30,'(/,A)') 'NEW background/model/integration algorithms:'
C
	WRITE(30,'(2X,A,I3)') 'Model spots for zone 1:',NMODEL_CEN
	WRITE(30,'(2X,A,I3)') 'Model spots for zones 2-5:',NMODEL_OUT
	WRITE(30,'(2X,A,F6.2)') 'Contour level for ellipse calculation:',CONTOUR
	WRITE(30,'(2X,A,F6.2)') 'Peak fraction target:',PFRAC_TARGET
	WRITE(30,'(2X,A,F6.1,A)') 'Peak fraction esd factor:',100.0*PFRAC_ERROR,' %'
	WRITE(30,'(2X,A,F5.1)') 'Area multiplier for peak ellipse:',AREA_FAC
	WRITE(30,'(2X,A,F6.1,A)') 'Rejection of overlapping neighbours:',TOLER_NBOUR,' %'
	WRITE(30,'(2X,A,F6.2)') 'Min allowed filling fraction:',FFILL
	IF( LCENTER ) THEN
	  WRITE(30,'(2X,A,F6.2)') 'Recenter all spots: YES'
	ELSE
	  WRITE(30,'(2X,A,F6.2)') 'Recenter all spots: NO'
	ENDIF
C
C Output crosstalk parameters
C
	WRITE(30,'(/,A,F5.2,2(2X,A,F5.2))') 'Crosstalk Parameters: value',
	1	CROSSTALK,'ellipse',CROSS_ELLIPSE,'smoother',CROSS_SMOOTH
C
C Store relative sizes of ellipses compared to core ellipse
C
	RAD_MULT(1)=1.0
	RAD_MULT(2)=SQRT(AREA_FAC)
	RAD_MULT(3)=SQRT(2.0*AREA_FAC)
C
	RETURN
	END


	SUBROUTINE ARGBOX_SETUP_DATA(ZCUTOFFS, BASE_NAME, NUMX,NUMY,NREFS, IOPT_HKL, IMAGE_ZERO, IPIXCEN)
C
	CHARACTER BASE_NAME*80
	REAL ZCUTOFFS(5)
C
	COMMON /GEOM_COM/ PIX_CEN(2),PIX_SIZE(2),DRUM_RAD
	COMMON /GAIN_COM/ GAIN,ESD_FACTOR,CROSSTALK,CROSS_ELLIPSE,CROSS_SMOOTH
	COMMON /MODELS_PFRAC_COM/ PFRAC_TARGET,PFRAC_ERROR,DPFRAC_ZERO,DPFRAC_GRAD
	COMMON /DEBUG_COM/ IDEBUG
C
C
	WRITE(30,'(//,A)') '=============== Data Input and Setup =============='
C
C Read image data (stripped of zero intensity value), and load in /RIMAGE_COM/
C
	CALL READ_IMAGE_FILE(IMAGE_ZERO, TRIM(BASE_NAME)//'.tif')
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
C
C Calculate global background and create smoothed image with bkg subtracted
C
	WRITE(6,'(1X,A)') 'Calculating background .....'
	CALL MAKE_GLOBAL_BKG(PIX_CEN,PIX_SIZE,DRUM_RAD,IPIXCEN)
	CALL CALC_SMOOTH_IMAGE()
C
C Output bkg and smoothed images if debug is on
C
	IF(IDEBUG .GT. 0) CALL WRITE_BACK_SMOOTH_IMAGES()
C
C Calculate the detector gain and output the effect on the counting statistics.
C
	CALL SETUP_GAIN
	ESD_FACTOR=SQRT(GAIN)
	WRITE(30,'(2X,A,F6.2)') 'Counting statistics are multiplied by =',ESD_FACTOR
C
C Apply the exclusion areas given in the parameter file
C
	CALL APPLY_REJECT_AREAS()
C
C Read reflections from the LaueG file
C
	CALL READ_HKLS_FILE(NREFS)
C
C If twinned, do initial categorisation of twins
C NB: This may reorder the reflection list
C
	IF(IOPT_HKL .EQ. 2) CALL TWIN_OVERLAP_ANALYSIS(NREFS)
C
C Calculate the background at the centre of each spot
C
	CALL CALC_SPOTS_BKG(NREFS)
C
C Create a table of all spots in a particular box
C
	CALL SETUP_BOXES_REFS(NREFS)
C
C Create a table of the maximum, or average, spot intensity with
C background intensity included
C
	CALL SETUP_SPOT_AVE(NREFS)
C
C Output the zones used for model spots
C
	WRITE(30,'(/,A)') 'Zones used for model spots'
	WRITE(30,'(2X,A,4(I5,A))') '1 : X =',NUMX/3,' -',NUMX*2/3,',  Y =',NUMY/3,' -',NUMY*2/3
	WRITE(30,'(2X,A,2(I5,A))') '2 : X <=',NUMX/2,',  Y <=',NUMY/2,' (excluding Zone 1)'
	WRITE(30,'(2X,A,2(I5,A))') '3 : X > ',NUMX/2,',  Y <=',NUMY/2,' (excluding Zone 1)'
	WRITE(30,'(2X,A,2(I5,A))') '4 : X <=',NUMX/2,',  Y > ',NUMY/2,' (excluding Zone 1)'
	WRITE(30,'(2X,A,2(I5,A))') '5 : X > ',NUMX/2,',  Y > ',NUMY/2,' (excluding Zone 1)'
C
C Load into IZONES_COM the zone numbers of all spots and calculate cutoffs
C per zone to get NWANT potential models
C
	NWANT=100
	CALL LOAD_IZONES(NREFS)
	CALL CALC_ZONE_CUTOFFS(ZCUTOFFS,NWANT,NREFS)
C
	RETURN
	END


	SUBROUTINE SETUP_GAIN()
C
C Calculate detector gain using the global background values.
C GAIN value (not corrected for crosstalk) is stored in /GAIN_COM/
C NB: Must be performed after global background calculation
C NB: This routine can be a hog, so put IY in outer loops
C
	COMMON /RIMAGE_COM/ RIMAGE(8000,2500)
	COMMON /BACK_COM/ BACK(8000,2500)
C
	COMMON /GAIN_COM/ GAIN,ESD_FACTOR,CROSSTALK,CROSS_ELLIPSE,CROSS_SMOOTH
C
	LOGICAL LREJECT(8000,2500)
	REAL SMOOTH(8000,2500),RDIFF2(8000,2500)
C
C Apply a wide margin to edges as they often contain artefacts
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
	IXLO=NUMX/10
	IXHI=NUMX-IXLO
	IYLO=NUMY/10
	IYHI=NUMY-IYLO
C
C Create map LREJECT() of pixels within 7 pixels of being
C rejected by holes, edges, etc. Only do every 5x5 pixel.
C
	  DO IY=IYLO,IYHI,5
	    DO IX=IXLO,IXHI,5
	      CALL REJECT_RADIUS(IOUT, FLOAT(IX),FLOAT(IY),7.0)
	      LREJECT(IX,IY)=(IOUT .EQ. 1)
	    ENDDO
	  ENDDO
C
C Average IMAGE() within +/-2 and put in SMOOTH(), put (diff/esd)^2 in RDIFF2()
C NB: To reduce computation only calculate for every 5x5 pixels
C
	DO IY=IYLO,IYHI,5
	  DO IX=IXLO,IXHI,5
	    SUM=0.0
	    DO IDY=-2,2
	      DO IDX=-2,2
	        SUM=SUM+RIMAGE(IX+IDX,IY+IDY)
	      ENDDO
	    ENDDO
	    SMOOTH(IX,IY)=SUM/25.0
	    RDIFF2(IX,IY)=25.0*(SMOOTH(IX,IY)-BACK(IX,IY))**2 / SMOOTH(IX,IY)
	  ENDDO
	ENDDO
C
C Iterate twice since algorithm needs an initial GAIN value
C
	GAIN=1.0
	DO ITER=1,2
	  CUTOFF2=(4.0*GAIN)**2
C
C Calculate the statistics for differences between original, background,
C and smoothed images. Avoid any pixels close to an edge, hole or any
C region where the difference between background and smoothed images is
C above a cutoff.
	  SUM_VAR1=0.0
	  SUM_BG1=0.0
	  SUM_VAR2=0.0
	  SUM_BG2=0.0
	  DO IY=IYLO+5,IYHI-5,5
	    DO IX=IXLO+5,IXHI-5,5
	      IF( LREJECT(IX,IY) ) CYCLE
		  RMAX2=MAX(RDIFF2(IX,IY),RDIFF2(IX-5,IY),RDIFF2(IX+5,IY),RDIFF2(IX,IY-5),RDIFF2(IX,IY+5))
	      IF(RMAX2 .GT. CUTOFF2) CYCLE
	      SUM_VAR1=SUM_VAR1+(RIMAGE(IX,IY)-BACK(IX,IY))**2
	      SUM_BG1=SUM_BG1+BACK(IX,IY)
	      SUM_VAR2=SUM_VAR2+(RIMAGE(IX,IY)-SMOOTH(IX,IY))**2
	      SUM_BG2=SUM_BG2+SMOOTH(IX,IY)
	    ENDDO
	  ENDDO
C
C Calculate gain from both statistics, then average the values
	  GAIN_AVE1=SUM_VAR1/SUM_BG1
	  GAIN_AVE2=SUM_VAR2/SUM_BG2*CROSS_SMOOTH**2
	  GAIN_AVE2=GAIN_AVE2*0.9		! empirical correction
	  GAIN=(GAIN_AVE1+GAIN_AVE2)/2.0
C
	ENDDO
C
C Output results
C
	GAIN=GAIN*CROSSTALK
	WRITE(30,'(/,A,3(/,2X,A,F6.2))')
	1	'Detector Gain calculated from background intensities:',
	2	'Using global background =',GAIN_AVE1,
	3	'Using smoothed image + correction =',GAIN_AVE2,
	4	'Average, corrected for crosstalk =',GAIN
C
	RETURN
	END


	SUBROUTINE SETUP_BOXES_REFS(NREFS)
C
C Store in NBOX(IX,IY,1..IBOX(IX,IY)) the spots numbers for
C all spots in the box (IX,IY).
C
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /NEIGH_BOX_COM/ IBOX(160,50,100),NBOX(160,50),NUMX_BOX,NUMY_BOX
C
C Calculate the limits of the boxes from the image size
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
	NUMX_BOX=1+(NUMX-1)/50
	NUMY_BOX=1+(NUMY-1)/50
	IF(NUMX_BOX.GT.160 .OR. NUMY_BOX.GT.50)
	1			CALL QUIT('ERROR: Image larger than 160*50 x 50*50 pixels')
C
C Zero arrays
C
	DO IX=1,NUMX_BOX
	  DO IY=1,NUMY_BOX
	    NBOX(IX,IY)=0
	    DO I=1,100
	      IBOX(IX,IY,I)=0
	    ENDDO
	  ENDDO
	ENDDO
C
C Put spot numbers into NBOX(), also find the limits of unempty boxes
C
	DO IR=1,NREFS
	  NX=MAX(1,MIN(NUMX_BOX, INT(X(IR)-1.0)/50+1 ))	! get box numbers
	  NY=MAX(1,MIN(NUMY_BOX, INT(Y(IR)-1.0)/50+1 ))
	  NBOX(NX,NY)=NBOX(NX,NY)+1
	  IF(NBOX(NX,NY) .GE. 100)
	1			CALL QUIT('ERROR: >100 reflns in a region of 50 x 50 pixels')
	  IBOX(NX,NY,NBOX(NX,NY))=IR
	ENDDO
C
	RETURN
C

C
	END


	SUBROUTINE SETUP_SPOT_AVE(NREFS)
C
C Find maximum in smoothed image within 4 pixels of each spot's nominal center
C
	COMMON /RIMAGE_COM/ RIMAGE(8000,2500)
	COMMON /BACK_COM/ BACK(8000,2500)
	COMMON /SMOOTH_COM/ SMOOTH_IMAGE(8000,2500)
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
C
	DO IR=1,NREFS
	  IXCEN=NINT(X(IR))
	  IYCEN=NINT(Y(IR))
C Get maximum values of RIMAGE() and SMOOTH_IMAGE()
	  RMAX0=0.0
	  RMAX(IR)=0.0
	  DO IX=MAX(1,IXCEN-4),MIN(NUMX,IXCEN+4)
	    DO IY=MAX(1,IYCEN-4),MIN(NUMY,IYCEN+4)
	      IF((IX-IXCEN)**2 + (IY-IYCEN)**2 .GT. 4**2) CYCLE
C Maximum smoothed intensity minus background
	      RMAX(IR)=MAX(RMAX(IR), SMOOTH_IMAGE(IX,IY) )
C Maximum raw pixel value
	      RMAX0=MAX(RMAX0, RIMAGE(IX,IY) )
	    ENDDO
	  ENDDO
C Add background to RMAX()
	  RMAX(IR)=RMAX(IR)+BACK(IXCEN,IYCEN)
C If any pixel was overloaded, make RMAX() overloaded
	  IF(RMAX0 .GT. 6E5) RMAX(IR)=RMAX0
C
	ENDDO
C
	RETURN
	END


C ============= (Image - Background) Smoothing ================

	SUBROUTINE CALC_SMOOTH_IMAGE()
C
C Subtract background from image, then smooth it
C NB: This routine can be a hog, so put IY in outer loops
c
c
	COMMON /RIMAGE_COM/ RIMAGE(8000,2500)
	COMMON /BACK_COM/ BACK(8000,2500)
	COMMON /SMOOTH_COM/ SMOOTH_IMAGE(8000,2500)
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
C
C Sanity checks
C
	IF(NUMX .GT. 8000) CALL QUIT('BUG(smooth_image): NUMX>8000')
	IF(NUMY .GT. 2500) CALL QUIT('BUG(smooth_image): NUMX>2500')
C
	DO IY=1,NUMY
	  DO IX=1,NUMX
		SMOOTH_IMAGE(IX,IY)=RIMAGE(IX,IY) - BACK(IX,IY)
	  ENDDO
	ENDDO
C
C Do two muffin tin smooths with radiuses of 1 and 2
C
	CALL SMOOTH_IMAGE2(SMOOTH_IMAGE, NUMX,NUMY,1)
	CALL SMOOTH_IMAGE2(SMOOTH_IMAGE, NUMX,NUMY,2)
C
	RETURN
	END


	SUBROUTINE SMOOTH_IMAGE2(RIMAGE, NUMX,NUMY,ISMOOTH)
C
C Optimised muffin-tin smooth of radius ISMOOTH pixels.
C CIMAGE() is used for both input and output.
C
	REAL RIMAGE(8000,2500)
C
	INTEGER IDX(-100:100)
	REAL RIMAGE2(8000,2500),RLINE(8001)
C
C Create list of X limits versus Y that define a circle of
C radius ISMOOTH. Also, count number of pixels in the circle.
C
	NXY=0
	DO IY=-ISMOOTH,ISMOOTH
	  IDX(IY)=IFIX( SQRT( ISMOOTH**2 - IY**2 + 0.1 ) )
	  NXY=NXY + 2*IDX(IY) + 1
	ENDDO
C
C Clear the output array
C
	DO IY=1,NUMY
	  DO IX=1,NUMX
	    RIMAGE2(IX,IY)=0.0
	  ENDDO
	ENDDO
C
C Average RIMAGE() according to the above list, and put values into RIMAGE2()
C
	DO IY2=1,NUMY
C
C Make cumulative sum of the line IMAGE(...,IY2)
	  RLINE(1)=0.0
	  DO IX=1,NUMX
	    RLINE(IX+1)=RLINE(IX)+RIMAGE(IX,IY2)
	  ENDDO
C
C Do all smoothing that uses the line IMAGE(...,IY2)
	  DO IDY=-ISMOOTH,ISMOOTH
	    IY=IY2+IDY
		IF(IY.LT.1 .OR. IY.GT.NUMY) CYCLE
C
C Do sum along X using the cumulative sum in RLINE()
	    DO IX=MAX( 1, 1+IDX(IDY) ), MIN( NUMX, NUMX-IDX(IDY) )
	      IX1=IX-IDX(IDY)
	      IX2=IX+IDX(IDY)
	      RIMAGE2(IX,IY)=RIMAGE2(IX,IY) + RLINE(IX2+1) - RLINE(IX1)
	    ENDDO
C
	  ENDDO
	ENDDO
C
C Copy from RIMAGE2() to RIMAGE() and normalize by NXY
C
	DO IY=1,NUMY
	  DO IX=1,NUMX
	    RIMAGE(IX,IY)=RIMAGE2(IX,IY)/NXY
	  ENDDO
	ENDDO
C
	RETURN
	END


C ============ Setup model spots in zones =====================

	SUBROUTINE LOAD_IZONES(NREFS)
C
	COMMON /IZONES_COM/ IZONES(20000)
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
C
C Determine zones (1-5) for every spot
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
	DO IR=1,NREFS
	  IX=NINT(X(IR))
	  IY=NINT(Y(IR))
	  IF(IX.GE.NUMX/3 .AND. IX.LE.NUMX*2/3 .AND.
	1		IY.GE.NUMY/3 .AND. IY.LE.NUMY*2/3) THEN
	    IZONES(IR)=1
	  ELSEIF(IX.LE.NUMX/2 .AND. IY.LE.NUMY/2) THEN
	    IZONES(IR)=2
	  ELSEIF(IY .LE. NUMY/2) THEN
	    IZONES(IR)=3
	  ELSEIF(IX .LE. NUMX/2) THEN
	    IZONES(IR)=4
	  ELSE
	    IZONES(IR)=5
	  ENDIF
	ENDDO
C
	RETURN
	END


	SUBROUTINE CALC_ZONE_CUTOFFS(ZCUTOFFS,NWANT,NREFS)
C
	REAL ZCUTOFFS(5)
C
	COMMON /IZONES_COM/ IZONES(20000)
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
C
	INTEGER NUM_ZONE(1000,5)
C
C Make histogram of "peak size" for each zone
C
	DO I1=1,1000
	  DO I2=1,5
	    NUM_ZONE(I1,I2)=0
	  ENDDO
	ENDDO
C
	DO IR=1,NREFS
	  BKG=GET_GLOBAL_BKG(X(IR),Y(IR))
	  PKBG=(RMAX(IR)-BKG)/SQRT(MAX(10.0, BKG ))
	  IPKBG=MAX(1,MIN(1000, NINT(10.0*PKBG) ))
	  NUM_ZONE(IPKBG,IZONES(IR))=NUM_ZONE(IPKBG,IZONES(IR))+1
	ENDDO
C
C Determine peak-size cutoff for each zone to get NWANT potential models
C
	NTOTAL=0
	DO IZ=1,5
	  N=0
	  DO I=1000,1,-1
	    N=N+NUM_ZONE(I,IZ)
	    IF(N .GE. NWANT) EXIT
	  ENDDO
	  ZCUTOFFS(IZ)=I*0.1
	  NTOTAL=NTOTAL+N
	ENDDO
C
C Output info on cutoffs for strong peaks in each zone
C
	WRITE(30,'(/,A,I4,A)') 'Initial Peak/esd(Bkg) cutoffs for',NWANT,' strong spots per zone'
	WRITE(30,'(2X,A,5I6)') 'Zone  ',(K,K=1,5)
	WRITE(30,'(2X,A,5F6.1)') 'Cutoff',ZCUTOFFS
	WRITE(30,'(A,I4,A)') 'Total of',NTOTAL,' spots considered for model peaks'
C
	RETURN
	END