	SUBROUTINE FIND_MERGED_SPOTS(DMERGE)
C
	COMMON /SPOTS_COM/ NSPOTS,XSPOT(10000),YSPOT(10000),CSPOT(10000)
C
C NSPOTS_MAX is the maximum number of spots we want to find
C
	NSPOTS_MAX=2000
C
C Find NSPOTS potential spots with approximate positions XSPOT() & YSPOT()
C
	PRINT '(/,1X,A,I5,A)','Searching for up to',NSPOTS_MAX,' spots'
	IRADIUS=DMERGE
	CALL FIND_COUNT_MAXIMA(XSPOT,YSPOT,NSPOTS, NSPOTS_MAX,IRADIUS)
	PRINT '(3X,I4,A)',NSPOTS,' local maxima found'
C
C Uses intensity weighted averages to center the spots
C
	AVE_CUTOFF=0.1
	IRADIUS=DMERGE/3.5
	CALL CENTER_SPOTS(IRADIUS,AVE_CUTOFF)
	IRADIUS=DMERGE/2.0
	CALL CENTER_SPOTS(IRADIUS,AVE_CUTOFF)
	PRINT '(1X,I6,A)',NSPOTS,' spots after initial rejection test'
C
C Merge any spots within IRADIUS pixels of each other
C
	IRADIUS=DMERGE/1.5
	CALL MERGE_SPOTS(XSPOT,YSPOT,CSPOT,NSPOTS, IRADIUS)
	PRINT '(1X,I6,A,I3,A)',NSPOTS,' spots after merging within',
	1					IRADIUS,' pixel radius'
C
C Recalculate new centers
C
	AVE_CUTOFF=1.0
	IRADIUS=DMERGE/3.5
	CALL CENTER_SPOTS(IRADIUS,AVE_CUTOFF)
	CALL CENTER_SPOTS(IRADIUS,AVE_CUTOFF)
	PRINT '(1X,I6,A)',NSPOTS,' spots after final rejection test'
C
C Merge any spots within IRADIUS pixels of each other
C
	IRADIUS=DMERGE
	CALL MERGE_SPOTS(XSPOT,YSPOT,CSPOT,NSPOTS, IRADIUS)
	PRINT '(1X,I6,A,I3,A)',NSPOTS,' spots after merging within',
	1					IRADIUS,' pixel radius'
C
	RETURN
	END


C =========== Functions to return stripped-image & back values ==============

	FUNCTION GET_STRIP(IX,IY)
C
	COMMON /IMAGE_BACK_COM/ RIMAGE(8000,2500),BACK(8000,2500),NUMX,NUMY
C
	IF(IX.GE.1 .AND. IX.LE.NUMX .AND. IY.GE.1 .AND. IY.LE.NUMY) THEN
	  GET_STRIP=RIMAGE(IX,IY) - BACK(IX,IY)
	ELSE
	  GET_STRIP=0.0
	ENDIF
C
	RETURN
      END


	FUNCTION GET_BACK(IX,IY)
C
	COMMON /IMAGE_BACK_COM/ RIMAGE(8000,2500),BACK(8000,2500),NUMX,NUMY
C
	IF(IX.GE.1 .AND. IX.LE.NUMX .AND. IY.GE.1 .AND. IY.LE.NUMY) THEN
	  GET_BACK=BACK(IX,IY)
	ELSE
	  GET_BACK=0.0
	ENDIF
C
	RETURN
	END


	SUBROUTINE GET_NUM_XY(NUMX0,NUMY0)
C
	COMMON /IMAGE_BACK_COM/ RIMAGE(8000,2500),BACK(8000,2500),NUMX,NUMY
C
	NUMX0=NUMX
	NUMY0=NUMY
C
	RETURN
	END


C========================================================================

  	SUBROUTINE FIND_COUNT_MAXIMA(XMAX,YMAX,NMAX, NMAX_MAX,IRADIUS)
C
C Create a list of the local maxima in intensity separated by at least
C IRADIUS pixels from other maxima.
C
C The routine returns the NMAX approximate positions in XMAX() & YMAX().
C The routine adjusts the signicance threshold for peaks so that
C NMAX is close to, but never exceeds, NMAX_MAX.
C
C
	REAL XMAX(NMAX_MAX),YMAX(NMAX_MAX)
C
	INTEGER IXY4MAX(2,100000)
	REAL CSIG(2000,625),CMAX(2000,625),CAVE(2000,625)
C
C If NMAX_MAX is too big we must also increase XMAX(), YMAX() & ITAGS()
C
	IF(NMAX_MAX .GT. 10000) CALL QUIT('BUG(find_count_maxima): NMAX_MAX > 10,000')
C
	CALL GET_NUM_XY(NUMX,NUMY)
C
C Create arrays of the maximum and average intensities for blocks of 4 x 4
C pixels and store in CMAX() & CAVE()
C
	NUMX4=NUMX/4
	NUMY4=NUMY/4
	DO IX4=1,NUMX4
	  DO IY4=1,NUMY4
	    CAVE(IX4,IY4)=0.0
	    CMAX(IX4,IY4)=0.0
	    DO IX=4*IX4-3,4*IX4
	      DO IY=4*IY4-3,4*IY4
	        CAVE(IX4,IY4)=CAVE(IX4,IY4)+GET_STRIP(IX,IY)
	        CMAX(IX4,IY4)=MAX(CMAX(IX4,IY4), GET_STRIP(IX,IY) )
	      ENDDO
	    ENDDO
	    CAVE(IX4,IY4)=CAVE(IX4,IY4)/16.0
	  ENDDO
      ENDDO
C
C Search for local maxima in CMAX (but avoiding edges), and store their
C significance in esd's in CSIG(). If not a valid maximum, store CSIG()=0.
C
	DO IX4=1,NUMX4
	  DO IY4=1,NUMY4
          CSIG(IX4,IY4)=0.0
	  ENDDO
	ENDDO
	DO IX4=4,NUMX4-3
	  DO IY4=4,NUMY4-3
          CSIG(IX4,IY4)=GET_MAXIMA_ESDS(CMAX,CAVE,IX4,IY4)
	  ENDDO
	ENDDO
C
C Merge any maxima within a radius of IRADIUS pixels into a single maxima which
C is added to IXMAX,IYMAX arrays.
C The routine limits NMAX to NMAX_MAX by returning the maxima with the
C greatest significance.
C
	IRAD4=IRADIUS/4
	CALL MERGE_MAXIMA(IXY4MAX,NMAX, CSIG,NMAX_MAX,IRAD4)
C
C Convert X & Y from 4 x 4 blocks to actual pixel values
C
      DO I=1,NMAX
        XMAX(I)=4*IXY4MAX(1,I)-1.5
        YMAX(I)=4*IXY4MAX(2,I)-1.5
      ENDDO
C
	RETURN
	END


	FUNCTION GET_MAXIMA_ESDS(CMAX,CAVE,IX,IY)
C
C Check if IX,IY is a local maxima at least CUTOFF esds above background.
C
	REAL CMAX(2000,625),CAVE(2000,625)
C
	LOGICAL VALID_PIXEL_CIRCLE
C
C Set return value for failure
C
      GET_MAXIMA_ESDS=0.0      
C
C Check IX,IY is a local maximum in CMAX, if not return
C
	IF(CMAX(IX,IY) .LE. CMAX(IX-1,IY)) RETURN
	IF(CMAX(IX,IY) .LE. CMAX(IX+1,IY)) RETURN
	TEMP=MAX( CMAX(IX-1,IY+1),CMAX(IX,IY+1),CMAX(IX+1,IY+1),
	1	      CMAX(IX-1,IY-1),CMAX(IX,IY-1),CMAX(IX+1,IY-1) )
	IF(CMAX(IX,IY) .LE. TEMP) RETURN
C
C Check if all pixels within radius of 8 are valid, if not return
C
	IF( .NOT.VALID_PIXEL_CIRCLE(4*IX-2,4*IY-2,8) ) RETURN
C
C Use the minimum of CAVE +/- 2 in X or Y as the background intensity
C
	BACK=MIN(CAVE(IX-2,IY),CAVE(IX+2,IY),CAVE(IX,IY-2),CAVE(IX,IY+2))
C
C Calculate background intensity including what was removed by STRIP_BACKGROUND
C
	STRIP=BACK+GET_BACK(4*IX,4*IY)
C
C Return the significance in esd's of the maxima
C
      GET_MAXIMA_ESDS=(CMAX(IX,IY)-BACK)/SQRT(MAX(10.0,STRIP))
C
	RETURN
	END



	SUBROUTINE MERGE_MAXIMA(IXY4MAX,NMAX, CSIG,NMAX_MAX,IRAD4)
C
C Accept close to, but never exceeding, NMAX_MAX largest values in CSIG()
C as possible maxima. Maxima within IRAD4, or closer, pixel-groups are
C merged into one.
C Position of maxima returned in IXMAX(1..NMAX),IYMAX(1..NMAX).
C
	INTEGER IXY4MAX(2,NMAX_MAX)
	REAL CSIG(2000,625)
C
	LOGICAL LSPOT(2000,625)
	INTEGER HIST(1000)
C
	CALL GET_NUM_XY(NUMX,NUMY)
C
C Create a frequency histogram of log(CSIG())
C
	DO IX=1,NUMX/4
	  DO IY=1,NUMY/4
	    ISIG=IFIX( 100.0*ALOG10(MAX(1.0, CSIG(IX,IY)*10.0 )) )
	    ISIG=MAX(1,MIN(1000, ISIG ))
	    HIST(ISIG)=HIST(ISIG)+1
	  ENDDO
	ENDDO
C
C Find the largest index with >NMAX_MAX points above it,
C then convert to a cutoff in CSIG()
C
	NMAX=0
	DO ISIG=1000,1,-1
	  NMAX=NMAX+HIST(ISIG)
	  IF(NMAX .GT. NMAX_MAX) EXIT
	ENDDO
	CUTOFF=10.0**(0.01*ISIG-0.99)
	PRINT '(3X,A,F7.2)','Using a cutoff for maxima of',CUTOFF
C
C Set all CSIG() values less than CUTOFF to zero
C
	DO IX=1,NUMX/4
	  DO IY=1,NUMY/4
	    IF(CSIG(IX,IY) .LT. CUTOFF) CSIG(IX,IY)=0.0
	  ENDDO
	ENDDO
C
C Set LSPOT if maxima is above CUTOFF and it is not within
C 3 pixel-groups of a more significant maxima
C
	DO IX=IRAD4+1,NUMX/4-IRAD4
	  DO IY=IRAD4+1,NUMY/4-IRAD4
	    IF(CSIG(IX,IY) .EQ. 0.0) THEN
	      LSPOT(IX,IY)=.FALSE.
	    ELSE
	      LSPOT(IX,IY)=.TRUE.
	      DO IX2=IX-IRAD4,IX+IRAD4
	        DO IY2=IY-IRAD4,IY+IRAD4
	          IF((IX-IX2)**2 + (IY-IY2)**2 .GT. IRAD4**2) CYCLE
	          IF(CSIG(IX2,IY2) .GT. CSIG(IX,IY)) LSPOT(IX,IY)=.FALSE.
	        ENDDO
	      ENDDO
	    ENDIF
	  ENDDO
	ENDDO
C
C Search LSPOT to make a new list of maxima
C
	NMAX=0
	DO IX=4,NUMX/4-3
	  DO IY=4,NUMY/4-3
	    IF( LSPOT(IX,IY) ) THEN
	      IF(NMAX .EQ. NMAX_MAX) EXIT
	      NMAX=NMAX+1
	      IXY4MAX(1,NMAX)=IX
	      IXY4MAX(2,NMAX)=IY
	    ENDIF
	  ENDDO
	ENDDO
C
	RETURN
	END



	SUBROUTINE CENTER_SPOTS(IRAD,CUTOFF)
C
C Use intensity weighted averages to calculate new centers of spots in
C the tables XCEN,YCEN(1..NCEN). The weighted averages use a circle of
C of radius IRAD surrounding the original position. If the average
C intensity is less than CUTOFF the spot is removed from the lists.
C NB: It is assumed that STRIP_BACKGROUND() has already been run
C
	COMMON /SPOTS_COM/ NSPOTS,XSPOT(10000),YSPOT(10000),CSPOT(10000)
C
C Loop through the lists of spots overwriting these lists with new
C X,Y values and skipping spots with low intensities
C
	NSPOTS2=0
	DO ISPOT=1,NSPOTS
C
	  IXCEN=NINT(XSPOT(ISPOT))
	  IYCEN=NINT(YSPOT(ISPOT))
	  NSUM=0
	  NBACK=0
C Calculate the intensity weighted average in a circle of IRAD radius
	  SUMC=0.0
	  SUMX=0.0
	  SUMY=0.0
	  SUMB=0.0
	  DO IX=IXCEN-IRAD,IXCEN+IRAD
	    DO IY=IYCEN-IRAD,IYCEN+IRAD
	      IDXY2=(IX-IXCEN)**2+(IY-IYCEN)**2
	      IF(IDXY2 .LE. IRAD**2) THEN
	        NSUM=NSUM+1
	        SUMC=SUMC+GET_STRIP(IX,IY)
	        SUMX=SUMX+IX*GET_STRIP(IX,IY)
	        SUMY=SUMY+IY*GET_STRIP(IX,IY)
	      ENDIF
	    ENDDO
	  ENDDO
	  XNEW=SUMX/MAX(1.0,SUMC)
	  YNEW=SUMY/MAX(1.0,SUMC)
C
C If the average intensity is above the cutoff:
	  IF(SUMC .GT. CUTOFF*NSUM) THEN
C Limit X,Y shifts to +/-IRAD
	    XNEW=MAX(XSPOT(ISPOT)-IRAD,MIN(XSPOT(ISPOT)+IRAD, XNEW ))
	    YNEW=MAX(YSPOT(ISPOT)-IRAD,MIN(YSPOT(ISPOT)+IRAD, YNEW ))
	    NSPOTS2=NSPOTS2+1
	    XSPOT(NSPOTS2)=XNEW
	    YSPOT(NSPOTS2)=YNEW
	    CSPOT(NSPOTS2)=SUMC/NSUM
	  ENDIF
C
	ENDDO
C
C Update the new size of the lists
C
	NSPOTS=NSPOTS2
C
	RETURN
	END



	SUBROUTINE MERGE_SPOTS(XSPOT,YSPOT,CSPOT,NSPOTS, IRAD)
C
C Given a list of NSPOT spot positions in XSPOT,YSPOT merge any spots
C within IRAD pixels of each other and make a new list with only
C the merged positions.
C
	REAL XSPOT(10000),YSPOT(10000),CSPOT(10000)
C
	LOGICAL VALID_PIXEL_CIRCLE
C
	LOGICAL LSPOT(10000)
C
C Loop through all pairs of spots and set LSPOT to true only if there
C is no better spot within RAD pixels of it
C
	DO I1=1,NSPOTS
	  LSPOT(I1)=.TRUE.
	  DO I2=1,NSPOTS
	    IF(I2 .EQ. I1) CYCLE
C Skip I2 if not within IRAD pixels of I1 spot
	    IF((XSPOT(I2)-XSPOT(I1))**2 + (YSPOT(I2)-YSPOT(I1))**2
	1					.GT. FLOAT(IRAD)**2) CYCLE
C Skip I2 if smaller than the I1 spot
	    IF(CSPOT(I2) .LE. CSPOT(I1)) CYCLE
C I1 is close to a better spot, so mark it as merged
	    LSPOT(I1)=.FALSE.
	  ENDDO
	ENDDO

C
C Only accept spots if clear of edges, holes, etc.
C
	DO I=1,NSPOTS
	  IF( LSPOT(I) ) LSPOT(I)=VALID_PIXEL_CIRCLE(NINT(XSPOT(I)),NINT(YSPOT(I)),IRAD)
	ENDDO

C
C Only keep elements where LSPOT() = TRUE
C
	NSPOTS2=0
	DO I=1,NSPOTS
	  IF( LSPOT(I) ) THEN
	    NSPOTS2=NSPOTS2+1
	    XSPOT(NSPOTS2)=XSPOT(I)
	    YSPOT(NSPOTS2)=YSPOT(I)
	    CSPOT(NSPOTS2)=CSPOT(I)
	  ENDIF
	ENDDO
	NSPOTS=NSPOTS2
C
	RETURN
	END


	SUBROUTINE GET_PEAK_ESDS(ESDS, IX,IY)
C
C Return the number of esds for the peak intensity
C
	IRADIUS=4
C
	CALL GET_NUM_XY(NUMX,NUMY)
C
C Sum intensities with/without background in a radius=4 circle
C
      IRADIUS=4
      SUM1=0.0
      SUM2=0.0
	DO IX2=MAX(1,IX-IRADIUS),MIN(NUMX,IX+IRADIUS)
	  DO IY2=MAX(1,IY-IRADIUS),MIN(NUMY,IY+IRADIUS)
	    IF((IX-IX2)**2 + (IY-IY2)**2 .GT. IRADIUS**2) CONTINUE
C
          SUM1=SUM1+GET_STRIP(IX2,IY2)
	    SUM2=SUM2+MAX(1.0, GET_STRIP(IX2,IY2)+GET_BACK(IX2,IY2) )
C
        ENDDO
	ENDDO
C
C Calculate obs/esd including gain factor for KOALA
C
	GAIN_FACTOR=2.2
	ESDS=SUM1/( GAIN_FACTOR*SQRT(SUM2) )
C
	RETURN
	END



	LOGICAL FUNCTION VALID_PIXEL_CIRCLE(IX,IY,IRAD)
C
	COMMON /EXCLUDE_COM/ NCIRC,ICIRC(3,999),NRECT,IRECT(4,999)
C
	VALID_PIXEL_CIRCLE=.FALSE.
C
	DO I=1,NCIRC
	  IF((IX-ICIRC(1,I))**2 + (IY-ICIRC(2,I))**2 .LE.
	1			(ICIRC(3,I)+IRAD)**2) RETURN
	ENDDO
C
	DO I=1,NRECT
C If we are within a rectangle expanded by 2*IRAD in X & Y, then we
C might also be within IRAD of the true size rectangle
	  IF( IX+IRAD.GE.IRECT(1,I) .AND. IX-IRAD.LE.IRECT(2,I) .AND. 
	1		IY+IRAD.GE.IRECT(3,I) .AND. IY-IRAD.LE.IRECT(4,I) ) THEN
C Test if any points in the circle are in the rectangle
	    DO IX2=IX-IRAD,IX+IRAD
	      DO IY2=IY-IRAD,IY+IRAD
	        IF((IX-IX2)**2 + (IY-IY2)**2 .GT. IRAD**2) CYCLE
	          IF( IX2.GE.IRECT(1,I) .AND. IX2.LE.IRECT(2,I) .AND. 
	1		      IY2.GE.IRECT(3,I) .AND. IY2.LE.IRECT(4,I)  ) RETURN
	      ENDDO
	    ENDDO
	  ENDIF
	ENDDO
C
	VALID_PIXEL_CIRCLE=.TRUE.
	RETURN
	END
