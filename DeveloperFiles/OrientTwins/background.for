C =======================================================================
C Common code for: argonne_boxes, find_spots, orient_twins, strip_back
C =======================================================================
C	SUBROUTINE CALC_BACK_ARRAY(BACK, RIMAGE)
C	SUBROUTINE CALC_BACK_ARRAY2(BACK, RIMAGE, IUNIT)
C ============ Helper routines for CALC_BACK_ARRAY() ===============
C	SUBROUTINE MOVE_CENTRAL_HOLE(RIMAGE)
C	SUBROUTINE GET_INT_OFFSET(SUMX,SUMY, IRHOLE,IXCEN,IYCEN, RIMAGE)
C	SUBROUTINE CREATE_BACK_MASK(LIMAGE, RIMAGE)
C	SUBROUTINE FLATTEN_BKG(CIMAGE,VIMAGE, RIMAGE,IRHOLE)
C	SUBROUTINE MAKE_4X4_ARRAYS(CAVE,SCAVE,CMIN,CMAX,LIMAGE4, CIMAGE,VIMAGE,LIMAGE)
C	SUBROUTINE REMOVE_4X4_PEAKS(BACK4, CAVE,SCAVE,CMIN,CMAX,LIMAGE4)
C	SUBROUTINE FILL_4X4_BACK(BACK4, LIMAGE4)
C	SUBROUTINE FILL_4X4_RING(BACK4, IRAD,LIMAGE4)
C	SUBROUTINE SMOOTH_4X4_IMAGE(BACK4, 5)
C	SUBROUTINE EXPAND_4X4_IMAGE(BACK, BACK4)
C	SUBROUTINE UNFLATTEN_BKG(BACK,IRHOLE)
C	SUBROUTINE ZERO_HOLE_COUNTS(BACK)
C ================== Set RIMAGE() within exclusions (Optional) ==========
C	SUBROUTINE IMAGE_EXCLUDE(RIMAGE,VALUE)
C =======================================================================

C ============ User entry routines ===============

	SUBROUTINE CALC_BACK_ARRAY(BACK, RIMAGE)
C
C Stub to run CALC_BACK_ARRAY2() with output to UNIT=6
C NB: Input and output arrays must be sized 8000x2500.
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
	REAL BACK(MAX_NUMX,MAX_NUMY),RIMAGE(MAX_NUMX,MAX_NUMY)
C
	CALL CALC_BACK_ARRAY2(BACK, RIMAGE,6)
C
	RETURN
	END


	SUBROUTINE CALC_BACK_ARRAY2(BACK, RIMAGE,IUNIT0)
C
C Calculate the background component of RIMAGE() and return it in BACK().
C Parameters are input using /EXCLUDE_COM/ & /BACK_GEOM_COM/.
C NB: Input and output arrays must be sized 8000x2500.
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
C
	REAL BACK(MAX_NUMX,MAX_NUMY),RIMAGE(MAX_NUMX,MAX_NUMY)
C
	COMMON /EXCLUDE_COM/ NCIRC,ICIRC(3,999),NRECT,IRECT(4,999)
	COMMON /BACK_GEOM_COM/ NUMX,NUMY,IXCEN,IYCEN,XSIZE,YSIZE,DRUM_RAD,IPIXCEN
C
	COMMON /BACK_MISC_COM/ NUMX4,NUMY4,IUNIT
C
	LOGICAL LIMAGE(MAX_NUMX,MAX_NUMY),LIMAGE4(MAX_NUMX/4,MAX_NUMY/4)
	REAL CIMAGE(MAX_NUMX,MAX_NUMY),VIMAGE(MAX_NUMX,MAX_NUMY)
C
	REAL BACK4(MAX_NUMX/4,MAX_NUMY/4)
	REAL CAVE(MAX_NUMX/4,MAX_NUMY/4),CMIN(MAX_NUMX/4,MAX_NUMY/4)
	REAL CMAX(MAX_NUMX/4,MAX_NUMY/4),SCAVE(MAX_NUMX/4,MAX_NUMY/4)
C
	IF(NUMX.GT.MAX_NUMX .OR. NUMY.GT.MAX_NUMY) CALL
	1	QUIT('BUG(calc_back_array2): Image exceeds maximum size')
C
	IUNIT=IUNIT0
	WRITE(IUNIT,'(/,1X,A)') 'Calculating global background using maxima removal'
C
C Reposition the central hole and the pixel centers if pixcen is not stable
C
	IF(IPIXCEN .EQ. 0) CALL MOVE_CENTRAL_HOLE(RIMAGE)
C
	IXHOLE=ICIRC(1,1)
	IYHOLE=ICIRC(2,1)
	IRHOLE=ICIRC(3,1)
C
C Create mask where LIMAGE() is true for pixels not in the central hole or
C excluded areas or very weak due to a deep scratch, edge, etc.
C
	CALL CREATE_BACK_MASK(LIMAGE, RIMAGE)
C
C Subtract and divide by components to approximately flatten background.
C CIMAGE() contains the flattened image and VIMAGE() its variance.
C
	CALL FLATTEN_BKG(CIMAGE,VIMAGE, RIMAGE,LIMAGE,IRHOLE)
C
C Create reduced images of 4x4 blocks of pixels with CAVE(),SCAVE(),CMIN(),
C being the CMAX(),LIMAGE4() are the average/esd/min/max/validity of the block,
C
	NUMX4=NUMX/4
	NUMY4=NUMY/4
	WRITE(IUNIT,'(2X,A,I7,A)') 'Reducing image to',NUMX4*NUMY4,' 4x4 blocks'
	CALL MAKE_4X4_ARRAYS(CAVE,SCAVE,CMIN,CMAX,LIMAGE4, CIMAGE,VIMAGE,LIMAGE)
C
C Copy CAVE() values to BACK4(), except if close to an intensity maximum then
C mark CAVE() with -1e20 to indicate this block is not a valid background.
C The radius of marked blocks increases with the intensity of the maxima.
C
	CALL REMOVE_4X4_PEAKS(BACK4, CAVE,SCAVE,CMIN,CMAX,LIMAGE4)
C
C Fill in the marked blocks using the nearby background in the unmarked blocks
C
	CALL FILL_4X4_BACK(BACK4, LIMAGE4)
C
C Smooth reduced background twice with a Top Hat smooth of 5 blocks radius.
C NB: This is similar to bilinear smoothing
C
	WRITE(IUNIT,'(2X,A)') 'Linear smoothing of background'
	CALL SMOOTH_4X4_IMAGE(BACK4, 5)
	CALL SMOOTH_4X4_IMAGE(BACK4, 5)
C
C Expand smoothed background to full size
C
	CALL EXPAND_4X4_IMAGE(BACK, BACK4)
C
C Reverse the background flattening correction on the smoothed background.
C
	CALL UNFLATTEN_BKG(BACK,IRHOLE)
C
C Zero BACK() in the exit hole (and entrance hole if it exists)
C
	CALL ZERO_HOLE_COUNTS(BACK)
C
	RETURN
	END


C ============ Helper routines for CALC_BACK_ARRAY() ===============
C NB: These routines can be hogs, so put IY before IX in long DO loops

	SUBROUTINE MOVE_CENTRAL_HOLE(RIMAGE)
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
C
	REAL RIMAGE(MAX_NUMX,MAX_NUMY)
C
	COMMON /EXCLUDE_COM/ NCIRC,ICIRC(3,999),NRECT,IRECT(4,999)
	COMMON /BACK_GEOM_COM/ NUMX,NUMY,IXCEN,IYCEN,XSIZE,YSIZE,DRUM_RAD,IPIXCEN
	COMMON /BACK_MISC_COM/ NUMX4,NUMY4,IUNIT
C
C Copy values from common
C
	IXHOLE=ICIRC(1,1)
	IYHOLE=ICIRC(2,1)
	IRHOLE=ICIRC(3,1)
C
C Recenter the central hole in X
C
	SUM0=0.0
	SUMX=0.0
	DO IX=IXHOLE-2*IRHOLE,IXHOLE+2*IRHOLE
	  DO IY=IYHOLE-4,IYHOLE+4
	    SUM0=SUM0+1.0/MAX(20.0,RIMAGE(IX,IY))
	    SUMX=SUMX+IX/MAX(20.0,RIMAGE(IX,IY))
	  ENDDO
	ENDDO
	IXHOLE=SUMX/SUM0
C
C Recenter the central hole in Y
C
	SUM0=0.0
	SUMY=0.0
	DO IY=IYHOLE-2*IRHOLE,IYHOLE+2*IRHOLE
	  DO IX=IXHOLE-4,IXHOLE+4
	    SUM0=SUM0+1.0/MAX(20.0,RIMAGE(IX,IY))
	    SUMY=SUMY+IY/MAX(20.0,RIMAGE(IX,IY))
	  ENDDO
	ENDDO
	IYHOLE=SUMY/SUM0
C
C Copy new central hole position to ICIRC(*,1)
C
	ICIRC(1,1)=IXHOLE
	ICIRC(2,1)=IYHOLE
C
C Output the new central hole position
C
	WRITE(IUNIT,'(2X,A,2(I4,A))') 'Central hole moved to (',
	1				IXHOLE,',',IYHOLE,')'
C
C Iterate grid search for IXCEN & IYCEN with min. SUMX & SUMY
C
	ISTEP=9
	DO ITER=1,3
C Grid search on IXCEN
	  IXBEST=0
	  BEST=1E9
	  DO IX=IXCEN-3*ISTEP,IXCEN+3*ISTEP,ISTEP
	    CALL GET_INT_OFFSET(SUMX,SUMY, IRHOLE,IX,IYCEN, RIMAGE)
	    IF(ABS(SUMX) .LT. BEST) THEN
	      IXBEST=IX
	      BEST=ABS(SUMX)
	    ENDIF
	  ENDDO
	  IXCEN=IXBEST
C Grid search on IYCEN
	  BEST=1E9
	  DO IY=IYCEN-3*ISTEP,IYCEN+3*ISTEP,ISTEP
	    CALL GET_INT_OFFSET(SUMX,SUMY, IRHOLE,IXCEN,IY, RIMAGE)
	    IF(ABS(SUMY) .LT. BEST) THEN
	      IYBEST=IY
	      BEST=ABS(SUMY)
	   ENDIF
	  ENDDO
	  IYCEN=IYBEST
C Reduce grid size by 3
	  ISTEP=ISTEP/3
	ENDDO
C
C Output the new pixel center (values returned in /BACK_GEOM_COM/)
C
	WRITE(IUNIT,'(2X,A,2(I4,A))') 'Pixel centre moved to (',
	1					IXCEN,',',IYCEN,')'
C
	RETURN
	END


	SUBROUTINE GET_INT_OFFSET(SUMX,SUMY, IRHOLE,IXCEN,IYCEN, RIMAGE)
C
C Used by MOVE_CENTRAL_HOLE()
C
C Returns the "intensity offset" around the central hole.
C SUMX is the sum for X > IXCEN minus the sum for X < IXCEN of all pixels
C within 2 and 4 * IRHOLE pixels from (IXCEN,IYCEN). SUMY is similar for Y.
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
C
	REAL RIMAGE(MAX_NUMX,MAX_NUMY)
C
	SUMX=0.0
	SUMY=0.0
	DO IY=IYCEN-4*IRHOLE,IYCEN+4*IRHOLE
	  DO IX=IXCEN-4*IRHOLE,IXCEN+4*IRHOLE
	    IRAD2=(IX-IXCEN)**2 + (IY-IYCEN)**2
		IF(IRAD2 .LE. (2*IRHOLE)**2) CYCLE
	    IF(IX .LT. IXCEN) SUMX=SUMX-RIMAGE(IX,IY)
	    IF(IX .GT. IXCEN) SUMX=SUMX+RIMAGE(IX,IY)
	    IF(IY .LT. IYCEN) SUMY=SUMY-RIMAGE(IX,IY)
	    IF(IY .GT. IYCEN) SUMY=SUMY+RIMAGE(IX,IY)
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE CREATE_BACK_MASK(LIMAGE, RIMAGE)
C
C Create mask for pixels not in the central hole or excluded areas.
C LIMAGE() is true if the X,Y pixel is valid.
C

	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
C
	LOGICAL LIMAGE(MAX_NUMX,MAX_NUMY)
	REAL RIMAGE(MAX_NUMX,MAX_NUMY)
C
	COMMON /EXCLUDE_COM/ NCIRC,ICIRC(3,999),NRECT,IRECT(4,999)
	COMMON /BACK_GEOM_COM/ NUMX,NUMY,IXCEN,IYCEN,XSIZE,YSIZE,DRUM_RAD,IPIXCEN
	COMMON /BACK_MISC_COM/ NUMX4,NUMY4,IUNIT
C
C Sanity check
C
	IF(NCIRC .LT. 1) CALL QUIT('BUG(create_back_mask): NCIRC  <1')
C
C Rough estimates of minimum, mean, and median intensities from 1/64 of the pixels
C
	RMIN=RIMAGE(1,1)
	SUM=0.0
	DO IY=1,NUMY,8
	  DO IX=1,NUMX,8
	    RMIN=MIN(RMIN, RIMAGE(IX,IY) )
	    SUM=SUM+RIMAGE(IX,IY)
	  ENDDO
	ENDDO
	RMEAN=SUM/(NUMX*NUMY/64)
C
	SUM=0.0
	DO IY=1,NUMY,8
	  DO IX=1,NUMX,8
	    SUM=SUM + MIN( MAX(10.0,2.0*RMEAN), RIMAGE(IX,IY) )
	  ENDDO
	ENDDO
	RMEDIAN=SUM/(NUMX*NUMY/64)
C
C Set LIMAGE() as TRUE, except for very low intensities which are
C probably edges, deep scratches, etc.
C
	CUTOFF=RMIN+0.01*(RMEDIAN-RMIN)
	DO IY=1,NUMY
	  DO IX=1,NUMX
	    LIMAGE(IX,IY)=( RIMAGE(IX,IY) .GE. CUTOFF)
	  ENDDO
	ENDDO
C
C Add circular exclusions as FALSE
C
	DO I=1,NCIRC
	  DO IDY=-ICIRC(3,I),ICIRC(3,I)
	    DO IDX=-ICIRC(3,I),ICIRC(3,I)
	      IF(IDX**2 + IDY**2 .GT. ICIRC(3,I)**2) CYCLE
	      IX=ICIRC(1,I)+IDX
	      IY=ICIRC(2,I)+IDY
	      IF(IX.GE.1 .AND. IX.LE.NUMX .AND.
	1          IY.GE.1 .AND. IY.LE.NUMY) LIMAGE(IX,IY)=.FALSE.
	    ENDDO
	  ENDDO
	ENDDO
C
C Add rectangular exclusions as FALSE
C NB: The margins are given by the first (or last) four rectangles
C
	DO I=1,NRECT
	  DO IY=MAX(1,IRECT(3,I)),MIN(NUMY,IRECT(4,I))
	    DO IX=MAX(1,IRECT(1,I)),MIN(NUMX,IRECT(2,I))
	      LIMAGE(IX,IY)=.FALSE.
	    ENDDO
	  ENDDO
	ENDDO
C
C Count and output the proportion of rejected or insensitive pixels
C
	NZERO=0
	DO IY=1,NUMY
	  DO IX=1,NUMX
	    IF( .NOT.LIMAGE(IX,IY) ) NZERO=NZERO+1
	  ENDDO
	ENDDO
C
	RZERO=FLOAT(NZERO)/(NUMX*NUMY)
	WRITE(IUNIT,'(F7.1,A)') RZERO*100.0,'% of pixels neglected as possibly inaccurate'
C
	RETURN
	END


	SUBROUTINE FLATTEN_BKG(CIMAGE,VIMAGE, RIMAGE,LIMAGE,IRHOLE)
C
C Flatten background for Y dependence and scatter around exit and entrance holes.
C CIMAGE() contains the flattened image and VIMAGE() its variance.
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
C
	LOGICAL LIMAGE(MAX_NUMX,MAX_NUMY)
	REAL CIMAGE(MAX_NUMX,MAX_NUMY),VIMAGE(MAX_NUMX,MAX_NUMY),RIMAGE(MAX_NUMX,MAX_NUMY)
C
	COMMON /BACK_GEOM_COM/ NUMX,NUMY,IXCEN,IYCEN,XSIZE,YSIZE,DRUM_RAD,IPIXCEN
C Keep RIMAGE_AVE for UNFLATTEN_BKG()
	COMMON /BACK_FLAT_COM/ RIMAGE_AVE
C
C Average RIMAGE() in 30x30 blocks around Y=Ycen and X for TTH=+/-90 degrees.
C
	IX_TTH90=NINT(DRUM_RAD/XSIZE*1.5708)
	NSUM=0
	SUM=0.0
	DO IDX=-30,30
	  DO IDY=-30,30
	    IF( LIMAGE(IXCEN-IX_TTH90+IDX,IYCEN+IDY) ) THEN
	      NSUM=NSUM+1
	      SUM=SUM+RIMAGE(IXCEN-IX_TTH90+IDX,IYCEN+IDY)
	    ENDIF
	    IF( LIMAGE(IXCEN+IX_TTH90+IDX,IYCEN+IDY) ) THEN
	      NSUM=NSUM+1
	      SUM=SUM+RIMAGE(IXCEN+IX_TTH90+IDX,IYCEN+IDY)
	    ENDIF
	  ENDDO
	ENDDO
	RIMAGE_AVE=SUM/NSUM
C
C Subtract Y dependent background, then multiply by a TTH dependant factor
C to flatten the scattering around exit/entrance holes.
C Factor is sqrt(|sin(TTH)|) except within the exit/entrance holes where
C it is set to the value used just outside the holes.
C
	DX_180=DRUM_RAD/XSIZE*3.14159
	SIN_MIN=SIN( (0.5*IRHOLE)*XSIZE/DRUM_RAD )
	DO IY=1,NUMY
	  DO IX=1,NUMX
C TTH_X is relative to exit or entrance holes
			DX=ABS(IX-IXCEN)
			IF(DX .GT. 0.5*DX_180) DX=ABS(DX_180-DX)
			TTH_X=DX*XSIZE/DRUM_RAD
C TTH_Y for vertical component of TTH
			TTH_Y=(IY-IYCEN)*YSIZE/DRUM_RAD
			TTH = SQRT( TTH_X**2 + TTH_Y**2 )
C Subtract RIMAGE_AVE corrected for 1/R^2 variation in Y
	    RIMAGE2=RIMAGE(IX,IY) - RIMAGE_AVE/(1.0 + TTH_Y**2 )
C Multiply remaining intensity by TTH factor to give CIMAGE()
			FACTOR = SQRT(MAX(SIN_MIN, ABS(SIN(TTH)) ))
	    CIMAGE(IX,IY)=RIMAGE2*FACTOR
C Calculate variance for CIMAGE() and put in VIMAGE()
	    VIMAGE(IX,IY)=MAX(3.0, RIMAGE(IX,IY) )*FACTOR**2
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE MAKE_4X4_ARRAYS(CAVE,SCAVE,CMIN,CMAX,LIMAGE4, CIMAGE,VIMAGE,LIMAGE)
C
C Make arrays for blocks of 4x4 pixels with CAVE(),SCAVE(),CMIN(),CMAX()
C being the average pixel value, its esd, minimum and maximum pixel values.
C Also, LIMAGE4() which is true if all 4x4 pixels are valid.
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
C
	LOGICAL LIMAGE4(MAX_NUMX/4,MAX_NUMY/4)
	REAL CAVE(MAX_NUMX/4,MAX_NUMY/4),CMIN(MAX_NUMX/4,MAX_NUMY/4)
	REAL CMAX(MAX_NUMX/4,MAX_NUMY/4),SCAVE(MAX_NUMX/4,MAX_NUMY/4)
C
	LOGICAL LIMAGE(MAX_NUMX,MAX_NUMY)
	REAL CIMAGE(MAX_NUMX,MAX_NUMY),VIMAGE(MAX_NUMX,MAX_NUMY)
C
	COMMON /BACK_MISC_COM/ NUMX4,NUMY4,IUNIT
C
C Create arrays for blocks of 4x4 pixels
C
	DO IY4=1,NUMY4
	  DO IX4=1,NUMX4
C
	    CAVE(IX4,IY4)=0.0
	    SCAVE(IX4,IY4)=0.0
	    CMIN(IX4,IY4)=+1E10
	    CMAX(IX4,IY4)=-1E10
C
          LIMAGE4(IX4,IY4)=.TRUE.
	    DO IX=4*IX4-3,4*IX4
	      DO IY=4*IY4-3,4*IY4
	        VAR=VIMAGE(IX,IY)
	        COUNTS=CIMAGE(IX,IY)
	        CAVE(IX4,IY4)=    CAVE(IX4,IY4) + COUNTS
	        SCAVE(IX4,IY4)=  SCAVE(IX4,IY4) + VAR
	        CMIN(IX4,IY4)=MIN(CMIN(IX4,IY4),  COUNTS)
	        CMAX(IX4,IY4)=MAX(CMAX(IX4,IY4),  COUNTS)
	        IF( .NOT.LIMAGE(IX,IY) ) LIMAGE4(IX4,IY4)=.FALSE.
	      ENDDO
	    ENDDO
C
	    CAVE(IX4,IY4)=CAVE(IX4,IY4)/16.0
	    SCAVE(IX4,IY4)=SQRT( SCAVE(IX4,IY4)/16.0 )
C
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE REMOVE_4X4_PEAKS(BACK4, CAVE,SCAVE,CMIN,CMAX,LIMAGE4)
C
C Uses an "erosion" method to remove peaks.
C
C Identifies maxima in the image and marks blocks using BACK()=-2E10 when the
C block is possibly contaminated by a peak, otherwise copies CAVE() to BACK().
C Maxima are ranked in size with the largest maxima causing a 4 block radius to be
C marked, while smaller maxima mark a radius of 3 to 0 blocks.
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
C
	LOGICAL LIMAGE4(MAX_NUMX/4,MAX_NUMY/4)
	REAL BACK4(MAX_NUMX/4,MAX_NUMY/4)
	REAL CAVE(MAX_NUMX/4,MAX_NUMY/4),CMIN(MAX_NUMX/4,MAX_NUMY/4)
	REAL CMAX(MAX_NUMX/4,MAX_NUMY/4),SCAVE(MAX_NUMX/4,MAX_NUMY/4)
C
	COMMON /BACK_MISC_COM/ NUMX4,NUMY4,IUNIT
C
	INTEGER IBINS(0:4)
	REAL RPEAK(MAX_NUMX/4,MAX_NUMY/4)
C
C Factor used to estimate the width of spots versus their intensity. An intensity increase
C of SHAPE_FACTOR corresponds to an increase in the spot radius by one 4x4 block of pixels.
C
	SHAPE_FACTOR=3.0
C
C Calculate RPEAK(), an estimate of peak height divided by esd
C
	DO IX4=1,NUMX4
	  DO IY4=1,NUMY4
		  IF(CMIN(IX4,IY4) .LT. 20.0) THEN
C Zero RPEAK() for pixels on edges or with deep narrow scratches
	      RPEAK(IX4,IY4)=0.0
	    ELSE
C For normal case, calculate height/esd
	      RPEAK(IX4,IY4)=( CMAX(IX4,IY4) - CAVE(IX4,IY4) )/MAX(3.0, SCAVE(IX4,IY4) )
	    ENDIF
	  ENDDO
	ENDDO
C
C Zero RPEAK() for any invalid blocks
C
	DO IX4=1,NUMX4
	  DO IY4=1,NUMY4
	    IF( .NOT.LIMAGE4(IX4,IY4) ) RPEAK(IX4,IY4)=0.0
	  ENDDO
	ENDDO
C
C Iterate reducing RPEAK() until marked blocks are less than 25%
C
	DO ITER=1,20
C
C Count how many pixels will be marked for radius of 0 to 4 blocks
	  DO I=0,4
	    IBINS=0
	  ENDDO
C
	  DO IX4=1,NUMX4
	    DO IY4=1,NUMY4
C
	      RPEAK(IX4,IY4)=RPEAK(IX4,IY4)*0.8
			  DO I=0,4
	        IF(RPEAK(IX4,IY4) .GT. SHAPE_FACTOR**I) IBINS(I)=IBINS(I)+1
	      ENDDO
C
	    ENDDO
	  ENDDO
C
C Rough estimate of the fraction of pixels that would be removed
	  RFILL=FLOAT(IBINS(0))/(NUMX4*NUMY4)
	  IF(RFILL .LT. 0.25) EXIT
C
	ENDDO
C
C Output some information
C
	WRITE(IUNIT,'(2X,A)') 'Number of maxima and assumed exclusion width:'
	WRITE(IUNIT,'(6X,I9,A20)') IBINS(4),'36 pixels',IBINS(3),'28 pixels',IBINS(2),'20 pixels',
	1														IBINS(1),'12 pixels',IBINS(0),'4 pixels'
C
C Expand size of maxima by 4,3,2,1 adjacent blocks
C
	DO ITER=4,1,-1
	  RCUTOFF=SHAPE_FACTOR**ITER
C
	  DO IX4=2,NUMX4-1
	    DO IY4=2,NUMY4-1
C Ignore blocks with RPEAK() below cutoff
	      IF(RPEAK(IX4,IY4) .LT. RCUTOFF) CYCLE
C Ensure adjacent blocks are above cutoff for the next iteration
	      RPEAK(IX4-1,IY4-1)=MAX(RPEAK(IX4-1,IY4-1), RCUTOFF/2.0)
	      RPEAK(IX4-1,IY4  )=MAX(RPEAK(IX4-1,IY4  ), RCUTOFF/2.0)
	      RPEAK(IX4-1,IY4+1)=MAX(RPEAK(IX4-1,IY4+1), RCUTOFF/2.0)
	      RPEAK(IX4  ,IY4-1)=MAX(RPEAK(IX4  ,IY4-1), RCUTOFF/2.0)
	      RPEAK(IX4  ,IY4  )=MAX(RPEAK(IX4  ,IY4  ), RCUTOFF/2.0)
	      RPEAK(IX4  ,IY4+1)=MAX(RPEAK(IX4  ,IY4+1), RCUTOFF/2.0)
	      RPEAK(IX4+1,IY4-1)=MAX(RPEAK(IX4+1,IY4-1), RCUTOFF/2.0)
	      RPEAK(IX4+1,IY4  )=MAX(RPEAK(IX4+1,IY4  ), RCUTOFF/2.0)
	      RPEAK(IX4+1,IY4+1)=MAX(RPEAK(IX4+1,IY4+1), RCUTOFF/2.0)
	    ENDDO
	  ENDDO
C
	ENDDO
C
C Now create BACK4() using CAVE() for all unmarked blocks (RPEAK() < 1), or
C mark block using BACK4()=-2E10 indicating possible peak contamination.
C
	DO IX4=1,NUMX4
	  DO IY4=1,NUMY4
	    IF(RPEAK(IX4,IY4) .LE. 1.0) THEN
	      BACK4(IX4,IY4)=CAVE(IX4,IY4)
	    ELSE
	      BACK4(IX4,IY4)=-2.0E10
	    ENDIF
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE FILL_4X4_BACK(BACK4, LIMAGE4)
C
C Fill the blocks marked as possibly contaminated by peaks with the
C average from surrounding unmarked blocks
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
C
	LOGICAL LIMAGE4(MAX_NUMX/4,MAX_NUMY/4)
	REAL BACK4(MAX_NUMX/4,MAX_NUMY/4)
C
	COMMON /BACK_MISC_COM/ NUMX4,NUMY4,IUNIT
C
C Count the number of marked pixels
C
	NMARK=0
	DO IX4=1,NUMX4
	  DO IY4=1,NUMY4
	    IF(BACK4(IX4,IY4) .LT. -1E10) NMARK=NMARK+1
	  ENDDO
	ENDDO
C
C Output number of marked blocks
C
	RMARK=FLOAT(NMARK)/(NUMX4*NUMY4)
	WRITE(IUNIT,'(2X,A,I7)') 'Total 4x4 blocks removed',NMARK
C
C
C Fill in marked pixels with averages over rings of radius 6 to 48
C NB: Exit loop early if we have removed all marked blocks.
C
	DO ITER=1,10
C
C Fill in blocks marked as contaminated by the average of surrounding unmarked blocks
C within a ring of thickness 2 and outer radius IRAD
	  IRAD=MIN(48,6*ITER)
	  CALL FILL_4X4_RING(BACK4, IRAD,LIMAGE4)
C
C Count the remaining marked blocks
	  NMARK=0
	  DO IX=1,NUMX4
	    DO IY=1,NUMY4
	      IF(BACK4(IX,IY) .LT. -1E10) NMARK=NMARK+1
	    ENDDO
	  ENDDO
C
C Exit loop if no marks left, else output how many to go
	  IF(NMARK .EQ. 0) EXIT
	  WRITE(IUNIT,'(2X,A,I7,A,I2,A)') 'Remaining blocks to fill',NMARK,' (Iter',ITER,')'
C
	ENDDO
C
C Warn if any blocks weren't filled, then fill them with zeroes
C
	IF(NMARK .NE. 0) THEN
	  WRITE(IUNIT,*) 'WARNING: Remaining excluded filled with zero'
	  DO IX=1,NUMX4
	    DO IY=1,NUMY4
	      IF(BACK4(IX,IY) .LT. -1E10) BACK4(IX,IY)=0.0
	    ENDDO
	  ENDDO
	ENDIF
C
	RETURN
	END


	SUBROUTINE FILL_4X4_RING(BACK4, IRAD,LIMAGE4)
C
C Fill in blocks marked as contaminated by the average of surrounding unmarked blocks
C within a ring of thickness 2 and outer radius IRAD
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
C
	LOGICAL LIMAGE4(MAX_NUMX/4,MAX_NUMY/4)
	REAL BACK4(MAX_NUMX/4,MAX_NUMY/4)
C
	COMMON /BACK_MISC_COM/ NUMX4,NUMY4,IUNIT
C
	INTEGER IXY(2,1000)
	REAL BACK4_SAVE(MAX_NUMX/4,MAX_NUMY/4)
C
C Create offsets in X,Y for the ring
C
	NXY=0
	DO IX=-IRAD,IRAD
	  DO IY=-IRAD,IRAD
	    IF( (IX**2+IY**2 .GT. (IRAD-2)**2) .AND. (IX**2+IY**2 .LE. IRAD**2) ) THEN
	      NXY=NXY+1
	      IXY(1,NXY)=IX
	      IXY(2,NXY)=IY
	    ENDIF
	  ENDDO
	ENDDO
C
C Copy BACK4() to BACK4_SAVE()
C
	DO IX=1,NUMX4
	  DO IY=1,NUMY4
	    BACK4_SAVE(IX,IY)=BACK4(IX,IY)
	  ENDDO
	ENDDO
C
C For blocks marked as contaminated by a peak, replace BACK() value by the
C average of unmarked blocks in the ring
	DO IX=1,NUMX4
	  DO IY=1,NUMY4
C Ignore block if not marked as contaminated by a peak
	    IF(BACK4_SAVE(IX,IY) .GT. -1E10) CYCLE
C
C Do summation over valid unmarked blocks in the ring
	    NVALID=0
	    NSUM=0
	    SUM=0.0
	    DO I=1,NXY
C Ignore blocks in the ring outside the image or inside excluded areas
	      IX2=IX+IXY(1,I)
	      IY2=IY+IXY(2,I)
	      IF(IX2.LT.1 .OR. IX2.GT.NUMX4 .OR. IY2.LT.1 .OR. IY2.GT.NUMY4) CYCLE
	      IF( .NOT.LIMAGE4(IX2,IY2) ) CYCLE
C Count number of valid blocks in the ring
	      NVALID=NVALID+1
C Count number of unmarked blocks in ring, and sum BACK() values
	      IF(BACK4_SAVE(IX2,IY2) .GT. -1E10) THEN
	        NSUM=NSUM+1
	        SUM=SUM+BACK4_SAVE(IX2,IY2)
	      ENDIF
	    ENDDO
C If significant number of summed blocks, overwrite BACK4(IX,IY) by the ring average
	    IF(NSUM .GT. NVALID/4) BACK4(IX,IY)=SUM/MAX(1,NSUM)
C
C Loop back for next block (IX,IY)
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE SMOOTH_4X4_IMAGE(CIMAGE, ISMOOTH)
C
C Peform a Top Hat smooth of radius ISMOOTH on CIMAGE()
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
C
	REAL CIMAGE(MAX_NUMX/4,MAX_NUMY/4)
C
	COMMON /BACK_MISC_COM/ NUMX4,NUMY4,IUNIT
C
	INTEGER IXY(2,1000)
	REAL CIMAGE2(MAX_NUMX/4,MAX_NUMY/4)
C
C Create an array of (X,Y) offsets to use for smoothing
C
	NXY=0
	DO IY=-ISMOOTH,ISMOOTH
	  DO IX=-ISMOOTH,ISMOOTH
	    IF(IX**2 + IY**2 .LE. ISMOOTH**2) THEN
	      NXY=NXY+1
	      IF(NXY .GT. 1000) CALL
	1			QUIT('BUG(smooth_background): NXY>1000')
	      IXY(1,NXY)=IX
	      IXY(2,NXY)=IY
	    ENDIF
	  ENDDO
	ENDDO
C
C Average CIMAGE() according to the above table, and put values into CIMAGE2()
C
	DO IY=1,NUMY4
	  DO IX=1,NUMX4
C
	    SUM=0.0
	    DO I=1,NXY
	      IX2=MAX(1,MIN(NUMX4, IX+IXY(1,I) ))
	      IY2=MAX(1,MIN(NUMY4, IY+IXY(2,I) ))
	      SUM=SUM+CIMAGE(IX2,IY2)
	    ENDDO
	    CIMAGE2(IX,IY)=SUM/NXY
C
	  ENDDO
	ENDDO
C
C Average CIMAGE() according to the above table, and put values into CIMAGE2()
C
	DO IY=1,NUMY4
	  DO IX=1,NUMX4
	    CIMAGE(IX,IY)=CIMAGE2(IX,IY)
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE EXPAND_4X4_IMAGE(CIMAGE, CAVE)
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
C
	REAL CIMAGE(MAX_NUMX,MAX_NUMY),CAVE(MAX_NUMX/4,MAX_NUMY/4)
C
	COMMON /BACK_MISC_COM/ NUMX4,NUMY4,IUNIT
C
	REAL CIMAGE2(MAX_NUMX,MAX_NUMY)
C
C Expand CAVE() to CIMAGE() using a bilinear interpolation
C
	NUMX=4*NUMX4
	NUMY=4*NUMY4
C
	DO IY=1,NUMY
	  DO IX=1,NUMX
	    IX4=1+(IX-3)/4
	    IY4=1+(IY-3)/4
	    IX4=MAX(1,MIN(NUMX4-1, IX4 ))
	    IY4=MAX(1,MIN(NUMY4-1, IY4 ))
	    RX=(4*IX4-IX+0.5)/4.0
	    RY=(4*IY4-IY+0.5)/4.0
	    CIMAGE(IX,IY)=CAVE(IX4,IY4)*(0.5+RX)*(0.5+RY) +
	1	              CAVE(IX4+1,IY4)*(0.5-RX)*(0.5+RY) +
	2	              CAVE(IX4,IY4+1)*(0.5+RX)*(0.5-RY) +
	3	              CAVE(IX4+1,IY4+1)*(0.5-RX)*(0.5-RY)
	  ENDDO
	ENDDO
C
C Do a 3x3 box convolution of CIMAGE()
C
	NWID=1
C Integrate CIMAGE() along Y, put in CIMAGE2()
	DO IX=1,NUMX
	  SUM=0.0
	  DO IY=1,NUMY
	    SUM=SUM+CIMAGE(IX,IY)
	    CIMAGE2(IX,IY)=SUM
	  ENDDO
	ENDDO
C Differentiate CIMAGE2() along Y, put in CIMAGE()
	DO IX=1,NUMX
	  SUM=0
	  DO IY=1+NWID,NUMY-NWID
	    CIMAGE(IX,IY)=( CIMAGE2(IX,IY+NWID) - SUM )/(1+2*NWID)
	    SUM=CIMAGE2(IX,IY-NWID)
	  ENDDO
	ENDDO
C Integrate CIMAGE() along X, put in CIMAGE2()
	DO IY=1,NUMY
	  SUM=0.0
	  DO IX=1,NUMX
	    SUM=SUM+CIMAGE(IX,IY)
	    CIMAGE2(IX,IY)=SUM
	  ENDDO
	ENDDO
C Differentiate CIMAGE2() along X, put in CIMAGE()
	DO IY=1,NUMY
	  SUM=0.0
	  DO IX=1+NWID,NUMX-NWID
	    CIMAGE(IX,IY)=( CIMAGE2(IX+NWID,IY) - SUM )/(1+2*NWID)
	    SUM=CIMAGE2(IX-NWID,IY)
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE UNFLATTEN_BKG(BACK, IRHOLE)
C
C Reverse the background flattening correction on the smoothed background.
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
C
	REAL BACK(MAX_NUMX,MAX_NUMY)
C
	COMMON /BACK_GEOM_COM/ NUMX,NUMY,IXCEN,IYCEN,XSIZE,YSIZE,DRUM_RAD,IPIXCEN
	COMMON /BACK_FLAT_COM/ RIMAGE_AVE
C
C Divide smoothed background by FACTOR used in FLATTEN_BKG(), then
C add Y dependent background also used in FLATTEN_BKG()
C
	DX_180=DRUM_RAD/XSIZE*3.14159
	SIN_MIN=SIN( (0.5*IRHOLE)*XSIZE/DRUM_RAD )
	DO IY=1,NUMY
	  DO IX=1,NUMX
C TTH_X is relative to exit or entrance holes
			DX=ABS(IX-IXCEN)
			IF(DX .GT. 0.5*DX_180) DX=ABS(DX_180-DX)
			TTH_X=DX*XSIZE/DRUM_RAD
C TTH_Y for vertical component of TTH
			TTH_Y=(IY-IYCEN)*YSIZE/DRUM_RAD
			TTH = SQRT( TTH_X**2 + TTH_Y**2 )
C Divide BACK() by FACTOR to give corrected intensity CIMAGE()
			FACTOR = SQRT(MAX(SIN_MIN, ABS(SIN(TTH)) ))
	    BACK(IX,IY)=BACK(IX,IY)/FACTOR
C Add 1/R^2 variation in Y component
	    BACK(IX,IY)=BACK(IX,IY) + RIMAGE_AVE/(1.0 + TTH_Y**2 )
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE ZERO_HOLE_COUNTS(BACK)
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
C
	REAL BACK(MAX_NUMX,MAX_NUMY)
C
	COMMON /EXCLUDE_COM/ NCIRC,ICIRC(3,999),NRECT,IRECT(4,999)
	COMMON /BACK_GEOM_COM/ NUMX,NUMY,IXCEN,IYCEN,XSIZE,YSIZE,DRUM_RAD,IPIXCEN
C
C Zero BACK() in the exit hole (and entrance hole if it exists)
C
	IRHOLE=ICIRC(3,1)+1
	IDRUM=NINT(3.14159*DRUM_RAD/XSIZE)
C For all pixels within IRHOLE of X,Y=ICIRC(1:2,1)
	DO IDX=-IRHOLE,IRHOLE
	  DO IDY=-IRHOLE,IRHOLE
		  IF(IDX**2+IDY**2 .LE. IRHOLE**2) THEN
C Set intensity to zero within exit hole
	      IX=ICIRC(1,1)+IDX
	      IY=ICIRC(2,1)+IDY
		    BACK(IX,IY)=0.0
C Set intensity to zero within entrance hole (if it exists)
				IF(IX+IDRUM .LE. NUMX) BACK(IX+IDRUM,IY)=0.0
				IF(IX-IDRUM .GE.		1) BACK(IX-IDRUM,IY)=0.0
	    ENDIF
	  ENDDO
	ENDDO
C
	RETURN
	END


C ================== Set RIMAGE() within exclusions (Optional) ============

	SUBROUTINE IMAGE_EXCLUDE(RIMAGE,VALUE)
C
C Optional routine to set RIMAGE() to VALUE within all exclusions
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500)
C
	REAL RIMAGE(MAX_NUMX,MAX_NUMY)
C
	COMMON /EXCLUDE_COM/ NCIRC,ICIRC(3,999),NRECT,IRECT(4,999)
	COMMON /BACK_GEOM_COM/ NUMX,NUMY,IXCEN,IYCEN,XSIZE,YSIZE,DRUM_RAD,IPIXCEN
C
C Zero RIMAGE() within all circular exclusions
C
	DO I=1,NCIRC
	  DO IDY=-ICIRC(3,I),ICIRC(3,I)
	    DO IDX=-ICIRC(3,I),ICIRC(3,I)
	      IF(IDX**2 + IDY**2 .GT. ICIRC(3,I)**2) CYCLE
	      IX=ICIRC(1,I)+IDX
	      IY=ICIRC(2,I)+IDY
	      IF(IX.GE.1 .AND. IX.LE.NUMX .AND.
	1                IY.GE.1 .AND. IY.LE.NUMY) RIMAGE(IX,IY)=VALUE
	    ENDDO
	  ENDDO
	ENDDO
C
C Zero RIMAGE() within all rectangular exclusions
C
	DO I=1,NRECT
	  DO IY=MAX(1,IRECT(3,I)),MIN(NUMY,IRECT(4,I))
	    DO IX=MAX(1,IRECT(1,I)),MIN(NUMX,IRECT(2,I))
	      RIMAGE(IX,IY)=VALUE
	    ENDDO
	  ENDDO
	ENDDO
C
	RETURN
	END
