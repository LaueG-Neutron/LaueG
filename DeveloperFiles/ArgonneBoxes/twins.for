C ==================== Routines in this file =========================
C	SUBROUTINE TWIN_OVERLAP_ANALYSIS(NREFS)
C	SUBROUTINE INTEGRATE_TWIN_CHECK(LOVERLAP, E,F,H,  IM,ITWIN,IR)
C	SUBROUTINE CALC_TWIN_INTINTS(IR, E,F,H,PEAK_FRAC,DPEAK_FRAC)
C	SUBROUTINE SUM_ELLIPSE_PAIR(CSUM,BSUM,NSUM, X1,Y1, X2,Y2, E,F,H)
C ====================================================================

	SUBROUTINE TWIN_OVERLAP_ANALYSIS(NREFS)
C
C Analyse the twinned pairs and change the M index, IHKLM(4,*)
C Output a summary of types of twinned pairs
C
C Assumes the spots are in three blocks of ascending HKLM
C	(1) twin M=1, (2) twin M=2, (3) pairs of twins M=11 & 12
C
C The M indices for twinned pairs are replaced with:
C  (11,12)  weak spots
C  (21,22)  strong twin 1
C  (31,32)  strong twin 2
C  (41,42)  separable spots
C  (51,52)  overlapping spots < 2 pixels separation
C  (61,62)  overlapping spots > 2 pixels separation
C
	COMMON /RIMAGE_COM/ RIMAGE(8000,2500)
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /SPOT2_COM/ IHKLM(4,20000),WAV(20000),TTH(20000),MULT(20000)
	COMMON /BGLOCAL_COM/ BKG(20000),BKG_ESD(20000),BKG_VAR(20000)
	COMMON /GAIN_COM/ GAIN,ESD_FACTOR,CROSSTALK,CROSS_ELLIPSE,CROSS_SMOOTH
C
C Count numbers of M=1,2 separate spots and M=11,12 pairs of spots
C
	M1=0
	DO I=1,NREFS
	  IF(IHKLM(4,I) .NE. 1) EXIT
        M1=I
	ENDDO
C
	M2=0
	DO I=M1+1,NREFS
	  IF(IHKLM(4,I) .NE. 2) EXIT
        M2=I-M1
      ENDDO
C
      DO I=M1+M2+1,NREFS,2
        IF( (IHKLM(1,I) .NE. IHKLM(1,I+1)) .OR.
     1      (IHKLM(2,I) .NE. IHKLM(2,I+1)) .OR. (IHKLM(3,I) .NE. IHKLM(3,I+1)) ) THEN
          CALL QUIT('BUG(twin_overlap_analysis): Invalid pair HKL')
        ENDIF
        IF( (IHKLM(4,I).NE.11) .OR. (IHKLM(4,I+1).NE.12) ) THEN
          CALL QUIT('BUG(twin_overlap_analysis): Invalid pair M')
        ENDIF
      ENDDO
C
      NTWIN1=M1
	NTWIN2=M2
	NPAIRS=(NREFS-M1-M2)/2
C
C ================================================
C Analyse the twin shape and update M index value
C
C Counter for pixel separation and twin shapes
C
	N2=0
	N5=0
	N10=0
	N20=0
	NFAR=0
C
	NCLOSE=0
	NSTR1=0
	NSTR2=0
	NOVER=0
C
C Loop through the twinned pairs
C
	DO IPAIR=1,NPAIRS
	  IR=NTWIN1+NTWIN2+2*IPAIR-1
C
	  DXY=SQRT( (X(IR)-X(IR+1))**2 + (Y(IR)-Y(IR+1))**2 )
C Update counters on the spot separation
	  IF(DXY .LT. 2.0) THEN
	    N2=N2+1
	  ELSEIF(DXY .LT. 5.0) THEN
	    N5=N5+1
	  ELSEIF(DXY .LT. 10.0) THEN
	    N10=N10+1
	  ELSEIF(DXY .LT. 20.0) THEN
	    N20=N20+1
	  ELSE
	    NFAR=NFAR+1
        ENDIF
C
C Get information on the intensities of the twins and the mid-point
	  IX1=NINT( X(IR) )
	  IY1=NINT( Y(IR) )
	  IX2=NINT( X(IR+1) )
	  IY2=NINT( Y(IR+1) )
	  IX12=NINT( (X(IR)+X(IR+1))/2.0 )
	  IY12=NINT( (Y(IR)+Y(IR+1))/2.0 )
	  CNT1=(	RIMAGE(IX1,IY1) + RIMAGE(IX1-1,IY1) + RIMAGE(IX1+1,IY1) + 
	1						  RIMAGE(IX1,IY1-1) + RIMAGE(IX1,IY1+1) )/5.0
	  CNT2=( RIMAGE(IX2,IY2) + RIMAGE(IX2-1,IY2) + RIMAGE(IX2+1,IY2) + 
	1						  RIMAGE(IX2,IY2-1) + RIMAGE(IX2,IY2+1) )/5.0
	  CNT12=( RIMAGE(IX12,IY12) + RIMAGE(IX12-1,IY12) + RIMAGE(IX12+1,IY12) + 
	1						  RIMAGE(IX12,IY12-1) + RIMAGE(IX12,IY12+1) )/5.0
	  ESD=SQRT(CNT12/5.0)*ESD_FACTOR
	  CNT1=CNT1 - GET_GLOBAL_BKG(FLOAT(IX1),FLOAT(IY1))
	  CNT2=CNT2 - GET_GLOBAL_BKG(FLOAT(IX2),FLOAT(IY2))
	  CNT12=CNT12 - GET_GLOBAL_BKG(FLOAT(IX12),FLOAT(IY12))
C
C Categorise the twin shape and update the M index and counters
C
C Closely overlapping twins separated by less than cutoff (51,52)
	  DXY_CUTOFF=5.0
	  IF(DXY .LT. DXY_CUTOFF) THEN
	      IHKLM(4,IR)  =51
	      IHKLM(4,IR+1)=52
	      NCLOSE=NCLOSE+1
	      CYCLE
	  ENDIF
C
C Cases for 1 or 2 strong twins
	  IF(MAX(CNT1,CNT2) .GT. 5.0*ESD) THEN
C Twin 1 stronger (21,22)
	    IF(CNT1 .GT. 5.0*CNT2) THEN
	      IHKLM(4,IR)  =21
	      IHKLM(4,IR+1)=22
	      NSTR1=NSTR1+1
	      CYCLE
C Twin 2 stronger (31,32)
	    ELSEIF(CNT2 .GT. 5.0*CNT1) THEN
	      IHKLM(4,IR)  =31
	      IHKLM(4,IR+1)=32
	      NSTR2=NSTR2+1
	      CYCLE
	    ENDIF
	  ENDIF
C Everything else, treat as possibly overlapping twins (61,62)
	  IHKLM(4,IR)  =61
	  IHKLM(4,IR+1)=62
	  NOVER=NOVER+1
C
C (41,42) are reserved for possibly overlapping twins that are
C integrated as a pair but then deemed to be independent
C
	ENDDO
C
C ================================================
C Output statistics on the twinned pairs
C
	WRITE(30,'(/,A)') 'Initial analysis of twinned spots'
C
	WRITE(30,'(2X,A,I5,A)') 'Spot separation:',N2,' < 2 pixels'
	WRITE(30,'(18X,I5,A)') N5,  '   2 -  5 pixels'
	WRITE(30,'(18X,I5,A)') N10, '   5 - 10 pixels'
	WRITE(30,'(18X,I5,A)') N20, '  10 - 20 pixels'
	WRITE(30,'(18X,I5,A)') NFAR,'  > 20 pixels'
C
	WRITE(30,'(2X,A)') 'Twin shapes:'
	WRITE(30,'(14X,I5,A)') NSTR1, '  strong twin 1 spot'
	WRITE(30,'(14X,I5,A)') NSTR2, '  strong twin 2 spot'
	WRITE(30,'(14X,I5,A)') NCLOSE,'  close overlapping spots, <2 pixels separation'
	WRITE(30,'(14X,I5,A)') NOVER, '  other overlapping spots, >2 pixels separation'
	WRITE(30,'(/,A)') 'Overlapping (>5 pixels) twins will be excluded from the model search'
C
	RETURN
	END


	SUBROUTINE INTEGRATE_TWIN_CHECK(LOVERLAP, E,F,H,  IM,ITWIN,IR)
C
	LOGICAL LOVERLAP
C
	COMMON /RIMAGE_COM/ RIMAGE(8000,2500)
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /SPOT2_COM/ IHKLM(4,20000),WAV(20000),TTH(20000),MULT(20000)
	COMMON /INTINT_COM/ RINTINT(20000),RSIGMA(20000),IOPTION(20000),IOVER(20000)
C
C Add info about the nominally twinned spot
C
	IF(IM .LT. 10) THEN
	  WRITE(30,'(A,I2,A)') 'Twin',ITWIN,', not part of a twinned pair'
	ELSEIF(IM .LT. 20) THEN
	  CALL QUIT('BUG(main): invalid IHKLM(4,IR) value')
	ELSEIF(IM .LT. 30) THEN
	  WRITE(30,'(A,I2,A)') 'Twin',ITWIN,' of a pair where Twin 1 is strong'
	ELSEIF(IM .LT. 40) THEN
	  WRITE(30,'(A,I2,A)') 'Twin',ITWIN,' of a pair where Twin 2 is strong'
	ELSEIF(IM .LT. 50) THEN
	  WRITE(30,'(A,I2,A)') 'Twin',ITWIN,' of a pair, but deemed independent'
	ELSEIF(IM .LT. 60) THEN
	  WRITE(30,'(A,I2,A)') 'Twin',ITWIN,' of a closely overlapped pair'
	ELSE
	  WRITE(30,'(A,I2,A)') 'Twin',ITWIN,' of an overlapped pair'
	ENDIF
C
C If not overlapping pairs, clear overlap flag and return
C
	IF(IM .LT. 40) THEN
		LOVERLAP=.FALSE.
		RETURN
	ENDIF
C
C If twin 2, just return
C
	IF(ITWIN .EQ. 2) RETURN
C
C Output separation and peak/midpoint intensities
C
	DX=X(IR+1)-X(IR)
	DY=Y(IR+1)-Y(IR)
	DRAD=CALC_ELLIPSE_DSQ(DX,DY, E,F,H)
	WRITE(30,'(A,F5.1,A,F5.1,A,F5.1,A)',IOSTAT=IDUMMY)
	1	'Twin Separation: (',DX,',',DY,') =',DRAD,' ellipse radii'
C
	INT1=RIMAGE( NINT(X(IR)) , NINT(Y(IR)) )
	INT2=RIMAGE( NINT(X(IR+1)) , NINT(Y(IR+1)) )
	INT12=RIMAGE( NINT(X(IR)+DX/2.0) , NINT(Y(IR)+DY/2.0) )
	WRITE(30,'(A,I6,2(2X,A,I6))',IOSTAT=IDUMMY)
	1		'Intensities: Twin 1',INT1,'Twin 2',INT2,'Midpoint',INT12
C
C If large separation, recategorise twin pair as separable spots
C
	IF(DRAD .GE. 4) THEN
	  WRITE(30,'(A)') 'Deemed not overlapped as twin separation > 4 radii'
	  IHKLM(4,IR)=1
	  IHKLM(4,IR+1)=2
	  LOVERLAP=.FALSE.
	ENDIF
C
	RETURN
	END


	SUBROUTINE CALC_TWIN_INTINTS(IR, E,F,H,PEAK_FRAC,DPEAK_FRAC)
C
C ISTATUS = 5	Success: overlapped twins integration
C
	COMMON /INTINT_COM/ RINTINT(20000),RSIGMA(20000),IOPTION(20000),IOVER(20000)
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /SPOT2_COM/ IHKLM(4,20000),WAV(20000),TTH(20000),MULT(20000)
	COMMON /BGLOCAL_COM/ BKG(20000),BKG_ESD(20000),BKG_VAR(20000)
	COMMON /GAIN_COM/ GAIN,ESD_FACTOR,CROSSTALK,CROSS_ELLIPSE,CROSS_SMOOTH
C
	INTEGER NSUM(3)
	REAL CSUM(3),BSUM(3)
C
	WRITE(30,'(A)') 'Overlapped twins integration using sig(I)/I'
C
C Sum counts and background in the core ellipses for pixels in 3 areas.
C 1 & 2 are for exclusively twins 1 & 2, 3 for the overlapping area.
C
	CALL SUM_ELLIPSE_PAIR(CSUM,BSUM,NSUM,
	1				X(IR),Y(IR),X(IR+1),Y(IR+1), E,F,H)
C
	WRITE(30,'(A)') 'Intensity & background summations over twinned spots:'
	WRITE(30,'(1X,A,I8,I7,A,I5,A)') 'Twin 1 ',NINT(CSUM(1)),NINT(BSUM(1)),
	1											'   (',NSUM(1),' pixels)'
	WRITE(30,'(1X,A,I8,I7,A,I5,A)') 'Twin 2 ',NINT(CSUM(2)),NINT(BSUM(2)),
	1											'   (',NSUM(2),' pixels)'
	WRITE(30,'(1X,A,I8,I7,A,I5,A)') 'Overlap',NINT(CSUM(3)),NINT(BSUM(3)),
	1											'   (',NSUM(3),' pixels)'
C
C Calc. integrated intensity and esd over all 3 areas.
C NB: Ignore the crosstalk effects for an ellipse
C
	RINT_12=CSUM(1)+CSUM(2)+CSUM(3) - ( BSUM(1)+BSUM(2)+BSUM(3) )
	RSIG_12=SQRT( (CSUM(1)+CSUM(2)+CSUM(3))*ESD_FACTOR**2 +
	1			( BKG_ESD(IR) * (NSUM(1)+NSUM(2)+NSUM(3)) )**2 )
C
C Calc. integrated intensity for twins 1 & 2 including half overlapped region
C
	ROVERLAP=(CSUM(3)-BSUM(3))/2.0
	RINT_1=CSUM(1)-BSUM(1) + ROVERLAP
	RINT_2=CSUM(2)-BSUM(2) + ROVERLAP
C
C Calculate individual uncertainties in RINT_1 & RINT_2.
C |RINT_1-RINT_2| can be in error up to +/-(CSUM(3)-BSUM(3))
C due to the half-intensity assumption. We assume the esd
C is on average half this uncertainty. Also assume the ratio
C of individual uncertainties in RINT_1,2 is equal to the
C ratio of RINT_1,2. This gives us the uncertainties:
C
	RSIG_1=ROVERLAP
	RSIG_2=ROVERLAP
C
C Add the counting statistics part of the esd's
C
	RSIG_1=SQRT( RSIG_1**2 + (CSUM(1)+CSUM(3)/2.0)*ESD_FACTOR**2 +
	1			( BKG_ESD(IR) * (NSUM(1)+NSUM(3)/2.0) )**2 )
	RSIG_2=SQRT( RSIG_2**2 + (CSUM(2)+CSUM(3)/2.0)*ESD_FACTOR**2 +
	1			( BKG_ESD(IR) * (NSUM(2)+NSUM(3)/2.0) )**2 )
C
C Add the uncertainty in peak fraction
C
	RSIG_1 =SQRT(RSIG_1**2  + ( MAX(0.0,RINT_1 )*DPEAK_FRAC/PEAK_FRAC )**2 )
	RSIG_2 =SQRT(RSIG_2**2  + ( MAX(0.0,RINT_2 )*DPEAK_FRAC/PEAK_FRAC )**2 )
	RSIG_12=SQRT(RSIG_12**2 + ( MAX(0.0,RINT_12)*DPEAK_FRAC/PEAK_FRAC )**2 )
C
C Do the peak fraction scaling
C
	RINT_1 =RINT_1 /PEAK_FRAC
	RINT_2 =RINT_2 /PEAK_FRAC
	RINT_12=RINT_12/PEAK_FRAC
	RSIG_1 =RSIG_1 /PEAK_FRAC
	RSIG_2 =RSIG_2 /PEAK_FRAC
	RSIG_12=RSIG_12/PEAK_FRAC
C
C Output integration information
C
	WRITE(30,'(A,I8,A,I6,A)',IOSTAT=IDUMMY) 'Integrated Intensity:',NINT(RINT_1),
	1								' +/-',NINT(RSIG_1),'   (Twin 1)'
	WRITE(30,'(A,I8,A,I6,A)',IOSTAT=IDUMMY) 'Integrated Intensity:',NINT(RINT_2),
	1								' +/-',NINT(RSIG_2),'   (Twin 2)'
	WRITE(30,'(A,I8,A,I6,A)',IOSTAT=IDUMMY) 'Integrated Intensity:',NINT(RINT_12),
	1								' +/-',NINT(RSIG_12),'   (Twins 1+2)'
C
C For the case of ROVERLAP < RSIG_12, the overlap contribution to twins 1 & 2 are
C deemed insignificant and the intensities are output as separable spots
C
	IF(ROVERLAP .LT. RSIG_12) THEN
	  WRITE(30,'(A)') 'Twin 1 & 2 treated as independent due to low overlap intensity'
	  WRITE(30,'(A)') 'The intensity for this reflection is the Twin 1 value'
	  RINTINT(IR)=RINT_1
	  RSIGMA(IR)=RSIG_1
	  IHKLM(4,IR)=41
	  IHKLM(4,IR+1)=42
	ELSE
	  WRITE(30,'(A)') 'The intensity for this reflection is the Twins 1+2 value'
	  RINTINT(IR)=RINT_12
	  RSIGMA(IR)=RSIG_12
	ENDIF
C
	WRITE(30,'(A)') 'The intensity for the next reflection is the Twin 2 value'
	RINTINT(IR+1)=RINT_2
	RSIGMA(IR+1)=RSIG_2
C
	RETURN
	END



	SUBROUTINE SUM_ELLIPSE_PAIR(CSUM,BSUM,NSUM, X1,Y1, X2,Y2, E,F,H)
C
C Return the total counts in CSUM() and the number of pixels in NSUM().
C
	INTEGER NSUM(3)
	REAL CSUM(3),BSUM(3)
C
	COMMON /RIMAGE_COM/ RIMAGE(8000,2500)
C
C Calculate tangents of ellipse
C
	TY=SQRT(E/(E*F-H**2))
	TX=SQRT(F/(E*F-H**2))
C
C Calculate extremes of both integration areas
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
	IXHI=MIN(NUMX, NINT(MAX(X1,X2) + TX) )
	IXLO=MAX(1, NINT(MIN(X1,X2) - TX) )
	IYHI=MIN(NUMY, NINT(MAX(Y1,Y2) + TY) )
	IYLO=MAX(1, NINT(MIN(Y1,Y2) - TY) )
C
C Zero the sums
C
	DO I=1,3
	  NSUM(I)=0
	  CSUM(I)=0.0
	  BSUM(I)=0.0
	ENDDO
C
C Sum intensity and background over both integration areas
C
	DO IX=IXLO,IXHI
	  DO IY=IYLO,IYHI
	    DSQ1=E*(IX-X1)**2 + F*(IY-Y1)**2 + 2.0*H*(IX-X1)*(IY-Y1)
	    DSQ2=E*(IX-X2)**2 + F*(IY-Y2)**2 + 2.0*H*(IX-X2)*(IY-Y2)
	    IF(DSQ1.LE.1.0 .AND. DSQ2.LE.1.0) THEN
	      NSUM(3)=NSUM(3)+1
	      CSUM(3)=CSUM(3)+RIMAGE(IX,IY)
	      BSUM(3)=BSUM(3)+GET_GLOBAL_BKG(FLOAT(IX),FLOAT(IY))
	    ELSEIF(DSQ1 .LE. 1.0) THEN
	      NSUM(1)=NSUM(1)+1
	      CSUM(1)=CSUM(1)+RIMAGE(IX,IY)
	      BSUM(1)=BSUM(1)+GET_GLOBAL_BKG(FLOAT(IX),FLOAT(IY))
	    ELSEIF(DSQ2 .LE. 1.0) THEN
	      NSUM(2)=NSUM(2)+1
	      CSUM(2)=CSUM(2)+RIMAGE(IX,IY)
	      BSUM(2)=BSUM(2)+GET_GLOBAL_BKG(FLOAT(IX),FLOAT(IY))
	    ENDIF
	  ENDDO
	ENDDO
C
	RETURN
	END
