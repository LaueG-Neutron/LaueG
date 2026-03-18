C ===================== Routines in this file ==============================
C	SUBROUTINE CALC_NEIGH_OVERLAP(ROVER,IOVER, IR, E,F,H,PEAK_FRAC,RAD_MULT, LVERBOSE)
C	SUBROUTINE CALC_NEIGH_INTS(RSIGMA,IWORST, IR, DIST,NEIGH,N_NEIGH,
C	1							E,F,H,PEAK_FRAC,RAD_MULT, LVERBOSE)
C	SUBROUTINE FIND_NEIGHBOURS(NEIGH,DIST,N_NEIGH, IR,E,F,H,RAD)
C ==========================================================================
	SUBROUTINE CALC_NEIGH_OVERLAP(ROVER,IOVER, IR, E,F,H,PEAK_FRAC,RAD_MULT, LVERBOSE)
C
	LOGICAL LVERBOSE
	REAL RAD_MULT(*)
C
C Calculate neighbours of peak IR and the effect of overlapping ellipses.
C Returns ROVER (ratio of overlap contribution to counting stats of core average)
C Returns IOVER (the worst overlap type: 1=none, 2=P-P, 3=C-P, 4=C-C)
C
	COMMON /SPOT2_COM/ IHKLM(4,20000),WAV(20000),TTH(20000),MULT(20000)
	COMMON /OPTIONS_COM/ IOPT_HKL,NHKLM
C
	INTEGER NEIGH(100)
	REAL DIST(100)
C
C Find neighbours where the peak and background ellipses touch
C
	CALL FIND_NEIGHBOURS(NEIGH,DIST,N_NEIGH, IR,E,F,H,RAD_MULT)
C
C If we filled the list then something is major wrong, so return as failure
C
	IF(N_NEIGH .GE. 100) THEN
	  ROVER=0.0
	  IOVER=-1
	  RETURN
	ENDIF
C
C Remove any peak-back overlaps
C
	N=0
	DO I=1,N_NEIGH
	  IF(DIST(I) .LE. RAD_MULT(2)+RAD_MULT(2)) THEN
	    N=N+1
	    NEIGH(N)=NEIGH(I)
	    DIST(N)=DIST(I)
	  ENDIF
	ENDDO
	N_NEIGH=N
C
C If an overlapping pair, remove other twin of pair
C
	IF(IOPT_HKL.EQ.2 .AND. IHKLM(4,IR).GT.10) THEN
	  IREJ=IR+1
	  IF(MODULO(IHKLM(4,IR),10) .EQ. 2) IREJ=IR-1
	  N=0
	  DO I=1,N_NEIGH
	    IF(NEIGH(I) .NE. IREJ) THEN
	      N=N+1
	      NEIGH(N)=NEIGH(I)
	      DIST(N)=DIST(I)
	    ENDIF
	  ENDDO
	  N_NEIGH=N
	ENDIF
C
C Estimate how much effect the remaining neighbours have on
C the intensities and optionally output some information.
C Update the worst overlap type
C
	CALL CALC_NEIGH_INTS(ROVER,IWORST, IR, DIST,NEIGH,N_NEIGH,
	1							E,F,H,PEAK_FRAC,RAD_MULT, LVERBOSE)
	IOVER=5-IWORST
C
	RETURN
	END


	SUBROUTINE CALC_NEIGH_INTS(RSIGMA,IWORST, IR, DIST,NEIGH,N_NEIGH,
	1							E,F,H,PEAK_FRAC,RAD_MULT, LVERBOSE)
C
C Returns RSIGMA the ratio of overlap contribution per pixel compared to
C the average core intensity.
C Also returns IWORST (the worst overlap type: 1=C-C, 2=C-P, 3=P-P, 4=none).
C
C
C DIST(I) is the distance between the target and neighbour spot divided by
C the size, RCORE, of the core ellipse along the line connecting the two spots.
C Along this line it is assumed the intensity is a gaussian given by
C   I(DIST) = I0 * EXP( - GAMMA * DIST**2 ), where GAMMA = - ALOG( 1 - PEAK_FRAC )
C and I0 is best calculated from the average intensity in the core elllipse
C   CSUM(1)/NSUM(1) - BKG = I0 * PEAK_FRAC / GAMMA
C The average, not integrated, intensity allows for the shape of the ellipse.
C
C If DIST is large the intensity varies exponentially not as a Gaussian. We
C will use a Gaussian with exponential tails that smoothly join the Gaussian
C when it's value is 0.25.
C For the Gaussian exp(-K*x^2 ), the function with tails is:
C   F(x) = exp(-K*x*x)        ; x < x0 = sqrt{ln(4)/K}
C        = 4*exp(-2*K*x0*x)   ; x > x0
C
C For C-P overlap, I assume the only contribution is from the pixel closest
C to the neighbour and any pixels within 1 pixel of this minimum distance.
C All these pixels will lie on the edge of the target's core ellipse closest
C to the neighbour, and the number is given by
C    NPIXELS = 2 * SQRT{ 2( Dpixels - Rpixels ) + 1 }
C where Dpixels & Rpixels are DIST & RCORE measured in terms of pixels.
C The average overlap contribution to each pixel in the core is then
C    NPIXELS * I(DIST-1) / NSUM(1)
C To be conservative we will actually calculate I() for 1 pixel closer
C to the neighbour spot.
C
C
	LOGICAL LVERBOSE
	INTEGER NEIGH(100)
	REAL DIST(100),RAD_MULT(3)
C
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
C
	CHARACTER*3 NAMES(3)
	INTEGER IDIST(100),NSUM(3)
	REAL RCORE(100),ROVER(100),CSUM(3)
C
C If no neighbours, set IWORST & RSIGMA for zero overlap and return
C
	IF(N_NEIGH .EQ. 0) THEN
	  IWORST=4
	  RSIGMA=0.0
	  RETURN
	ENDIF
C
C Store in IDIST() the degree of ellipse overlap: 1=C-C, 2=C-P, 3=P-P
C
	DO I=1,N_NEIGH
	  IF(DIST(I) .LT. 2.0*RAD_MULT(1)) THEN
	    IDIST(I)=1
	  ELSEIF(DIST(I) .LT. RAD_MULT(1)+RAD_MULT(2)) THEN
	    IDIST(I)=2
	  ELSE
	    IDIST(I)=3
	  ENDIF
	ENDDO
C
C Calculate the ave. core - background intensity of the target peak
C
	CALL SUM_ELLIPSE_INTS(CSUM,NSUM, X(IR),Y(IR),E,F,H,RAD_MULT,2)
	RCORE0=CSUM(1)/NSUM(1) - GET_GLOBAL_BKG(X(IR),Y(IR))
	RCORE0_ESD=SQRT(CSUM(1)) / NSUM(1)
C
C Calculate the ave. core - background intensity of the neighbours
C
	DO I=1,N_NEIGH
	  I2=NEIGH(I)
	  CALL SUM_ELLIPSE_INTS(CSUM,NSUM, X(I2),Y(I2),E,F,H,RAD_MULT,2)
	  BKG_NEIGH=GET_GLOBAL_BKG(X(I2),Y(I2))
	  RCORE(I)=CSUM(1)/NSUM(1)-BKG_NEIGH
	ENDDO
C
	GAMMA=-ALOG(1.0-PEAK_FRAC)
C
C Do an approximate correction to RCORE() due to double counted intensities
C in overlapping cores (ignore neighbour-neighbour overlaps)
C
	DO I=1,N_NEIGH
	  RCORE_CORE=EXP(-DIST(I)**2) *(0.40+0.75*DIST(I))*MAX(0.0,2.0-DIST(I))
	  RCORE0=RCORE0-RCORE(I)*RCORE_CORE
	ENDDO
C
	DO I=1,N_NEIGH
	  RCORE_CORE=EXP(-DIST(I)**2) *(0.40+0.75*DIST(I))*MAX(0.0,2.0-DIST(I))
	  RCORE(I)=RCORE(I)-RCORE0*RCORE_CORE
	ENDDO
C
C Estimate the intensity increase in the average core intensity due to
C overlap with neighbours. Use different algorithms for C-C & C-P overlaps.
C
	DO I=1,N_NEIGH
	  ROVER(I)=0.0
C Add the core-peak overlap, if spots are not too close together
	  IF(DIST(I) .GE. 1) THEN
	    I2=NEIGH(I)
	    DIST_PIXELS=SQRT( (X(IR)-X(I2))**2 + (Y(IR)-Y(I2))**2 )
	    CORE_PIXELS=DIST_PIXELS/DIST(I)
C Calculate the number of pixels within 1 pixel to the closest distance
	    PIXELS=2.0*SQRT( 2.0*(DIST_PIXELS-CORE_PIXELS) + 1.0 )
C Calculate the intensity for 1 pixel closer than the closest distance
	    RDIST=DIST(I)-1.0-1.0/CORE_PIXELS
C Calculate the Gaussian with exponential tails
	    RDIST0=SQRT( ALOG(4.0)/GAMMA )
	    IF(RDIST .LE. RDIST0) THEN
	      GAUSS=EXP(-GAMMA*RDIST**2)
	    ELSE
	      GAUSS=4.0*EXP(-2.0*GAMMA*RDIST0*RDIST)
	    ENDIF
C Finally, calculate the overlap contribution
	    ROVER(I)=(RCORE(I)/PEAK_FRAC*GAMMA)*GAUSS * PIXELS/NSUM(1)
	  ENDIF
C Add the core-core overlap, if the cores touch
	  IF(DIST(I) .LT. 2) THEN
	    ROVER(I)=(0.40+0.75*DIST(I)) * (2.0-DIST(I)) *
	1				RCORE(I) * EXP(-DIST(I)**2)
	  ENDIF
	  ROVER(I)=MAX(0.0, ROVER(I) )
	ENDDO
C
C Save the IDIST() of the worst overlap in IWORST (but favour a lower IDIST)
C Also calculate the max & total overlap contribution
C
	ROVER_MAX=0.0
	SUM_OVER=0.0
	IWORST=4
	DO I=1,N_NEIGH
	  IF(ROVER(I) .GE. ROVER_MAX+0.05) THEN
	    IWORST=IDIST(I)
	  ELSE IF(ROVER(I) .GE. ROVER_MAX-0.05) THEN
	    IWORST=MIN(IWORST,IDIST(I))
	  ENDIF
	  ROVER_MAX=MAX(ROVER_MAX,ROVER(I))
	  SUM_OVER=SUM_OVER+ROVER(I)
	ENDDO
C
C Optionally output information on the overlaps
C
	IF (LVERBOSE) THEN
	  WRITE(30,'(A)') 'Overlapping neighbours:'
	  NAMES(1)='C-C'
	  NAMES(2)='C-P'
	  NAMES(3)='P-P'
	  DO I=1,N_NEIGH
	    I2=NEIGH(I)
	    WRITE(30,'(2X,A,I5,A,I4,A,I4,A,F4.1,1X,2A,I6,A,F6.1)',IOSTAT=IDUMMY)
	1		'Ref',I2,' (',NINT(X(I2)),',',NINT(Y(I2)),') Rmult',DIST(I),
	2		NAMES(IDIST(I)),' Icore',NINT(RCORE(I)),' Iover',ROVER(I)
	  ENDDO
	  WRITE(30,'(2X,A,F7.1,1X,A)',IOSTAT=IDUMMY)
	1		'Estimated intensity overlap:',SUM_OVER,'per pixel'
	ENDIF
C
C Calculate degree of total overlap by the ratio of overlap contribution
C to counting stats of core average
C
	RSIGMA=SUM_OVER/MAX(0.1,RCORE0_ESD)
C
	RETURN
	END


	SUBROUTINE FIND_NEIGHBOURS(NEIGH,DIST,N_NEIGH, IR,E,F,H,RAD)
C
C Find neighbours to peak IR where the peak and background ellipses touch
C Also calculate inter-spot distances as ratios of the core ellipse size
C NB: This list does not include the target peak as a "neighbour"
C
	INTEGER NEIGH(100)
	REAL DIST(100),RAD(3)
C
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
C
	COMMON /NEIGH_BOX_COM/ IBOX(160,50,100),NBOX(160,50),NUMX_BOX,NUMY_BOX
C
C Box number of target peak
C
	NX=MAX(1,MIN(NUMX_BOX, INT(X(IR)-1.0)/50+1 ))
	NY=MAX(1,MIN(NUMY_BOX, INT(Y(IR)-1.0)/50+1 ))
C
C Limits in X & Y pixel positions which could be in background ellipse (or box)
C NB: IXLO .. IYHI can exceed physical limits of data (correct this later on)
C
	TY=SQRT(E/(E*F-H**2))
	TX=SQRT(F/(E*F-H**2))
	IXLO=INT(X(IR)-RAD(3)*TX)
	IXHI=INT(X(IR)+RAD(3)*TX)
	IYLO=INT(Y(IR)-RAD(3)*TY)
	IYHI=INT(Y(IR)+RAD(3)*TY)
C
C Make list of possibly overlapping neighbours in NEIGH(1..N_NEIGH)
C NB: List also contains the target peak 
C Calculate distance between target centre and neighbour ellipses, store in DIST(,).
C Units are multiples of target core ellipse size along a line joining centres.
C
	N_NEIGH=0
	NCORES=0
	DO NX2=MAX(1,NX-1),MIN(NUMX_BOX,NX+1) ! search over 3 x 3 boxes surrounding peak
	  DO NY2=MAX(1,NY-1),MIN(NUMY_BOX,NY+1)
	    DO I=1,NBOX(NX2,NY2)
	      I2=IBOX(NX2,NY2,I)
C Ignore neighbour if it's the target peak
	      IF(I2 .EQ. IR) CYCLE
	      DSQ=CALC_ELLIPSE_DSQ(X(I2)-X(IR),Y(I2)-Y(IR),E,F,H)
C Ignore neighbour if ellipses don't touch
	      IF(DSQ .GE. (RAD(2)+RAD(3))**2) CYCLE
C Add neighbour to the list, unless list is full
	      IF(N_NEIGH .EQ. 100) RETURN
	      N_NEIGH=N_NEIGH+1
	      NEIGH(N_NEIGH)=I2		! reflection number of neighbour
	      DIST(N_NEIGH)=SQRT(DSQ)
	    ENDDO
	  ENDDO
	ENDDO
C
	RETURN
	END
