C ===================== Routines in this file ==============================
C	SUBROUTINE INTEGRATE_ALL(NREFS, RAD_MULT,TOLER_NBOUR)
C	SUBROUTINE INTEGRATE_CLEAR(RINTINT,RSIGMA, NREFS)
C	SUBROUTINE INTEGRATE_SPOT(IR, RAD_MULT,TOLER_NBOUR, DPFRAC,PFRAC_OBS,PFRAC_SIG)
C	SUBROUTINE CALC_INT_INT(RINTINT,RSIGMA,ROVER,IOVER,ISTATUS,
C	1                        PFRAC_OBS,PFRAC_SIG, X,Y,BKG_ESD,
C	2											  E,F,H,RAD_MULT, PFRAC,DPFRAC,  TOLER_NBOUR)
C	SUBROUTINE CALC_INT_INT_STRONG(RINTINT,RSIGMA, CSUM,NSUM, BKG,BKG_ESD)
C	SUBROUTINE CALC_INT_INT_WEAK(RINTINT,RSIGMA, CSUM,NSUM, BKG,BKG_ESD, PFRAC,DPFRAC)
C ==========================================================================

	SUBROUTINE INTEGRATE_ALL(NREFS, RAD_MULT,TOLER_NBOUR)
C
	REAL RAD_MULT(3)
C
	COMMON /OPTIONS_COM/ IOPT_HKL,NHKLM
	COMMON /LIB_COM/ PLIB(1000,6),ILIB(1000),NLIB
C
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /SPOT2_COM/ IHKLM(4,20000),WAV(20000),TTH(20000),MULT(20000)
	COMMON /BGLOCAL_COM/ BKG(20000),BKG_ESD(20000),BKG_VAR(20000)
	COMMON /INTINT_COM/ RINTINT(20000),RSIGMA(20000),IOPTION(20000),IOVER(20000)
C
	COMMON /PFRAC_CHECK_COM/ PFRAC_LIST(3,20000),PFRAC_TARGET,NPFRAC
C
	COMMON /INT_COUNTERS_COM/ NSUCCESS,NFAIL,NFULL,NWEAK1,NWEAK2,NWEAK3,
	1			NWEAK4, NPAIR,NEDGE,NNOBACK,NPIXEL,NNEIGH,NFAR,NELLIP
C
	CHARACTER SHKLM*34
C
C Zero integration counters and set arrays as if integration failed
C
	CALL INTEGRATE_CLEAR(RINTINT,RSIGMA, NREFS)
C
C Zero number of results in PFRAC_LIST()
C
	NPFRAC=0
C
C Start loop over all reflections
C
	DO IR=1,NREFS
C
C Output a banner and basic info about the spot
	  WRITE(30,'(A)') '**************************************'
C
C Output hkl (or hklm), X, Y, Wav
	  CALL MAKE_HKLM_IDENT(SHKLM, IR,IHKLM(1,IR))
	  WRITE(30,'(A,2(A,I4),A,F5.2)',IOSTAT=IDUMMY) TRIM(SHKLM),
	1          '  XY (',NINT(X(IR)),',',NINT(Y(IR)),')  Wav',WAV(IR)
C
C If a single spot, try to integrate spot intensity
C If an overlapping twin pairs, integrate both for Twin 1, and nothing for Twin 2
	  CALL INTEGRATE_SPOT(IR, RAD_MULT,TOLER_NBOUR, DPFRAC,PFRAC_OBS,PFRAC_SIG)
C
C For success, increment counter and log the integrated intensity
 	  IF(IOPTION(IR) .GT. 0) THEN
	    NSUCCESS=NSUCCESS+1
	    WRITE(30,'(A,I8,A,I6)',IOSTAT=IDUMMY) 'Integrated Intensity:',
	1			NINT(RINTINT(IR)),' +/-',NINT(RSIGMA(IR))
C If no overlaps and PFRAC_SIG reasonable, update PFRAC_LIST()
	    IF(IOVER(IR).EQ.1 .AND. PFRAC_SIG.LE.0.5) THEN
	      NPFRAC=NPFRAC+1
	      PFRAC_LIST(1,NPFRAC)=PFRAC_OBS
	      PFRAC_LIST(2,NPFRAC)=PFRAC_SIG
	      PFRAC_LIST(3,NPFRAC)=DPFRAC
	    ENDIF
	  ELSE
C
C For failure, increment counter and zero IOVER
	    NFAIL=NFAIL+1
	    IOVER(IR)=0
	  ENDIF
C	
	ENDDO
C
	RETURN
	END


	SUBROUTINE INTEGRATE_CLEAR(RINTINT,RSIGMA, NREFS)
C
	REAL RINTINT(NREFS),RSIGMA(NREFS)
C
	COMMON /INT_COUNTERS_COM/ NSUCCESS,NFAIL,NFULL,NWEAK1,NWEAK2,NWEAK3,
	1			NWEAK4, NPAIR,NEDGE,NNOBACK,NPIXEL,NNEIGH,NFAR,NELLIP
C
C Total number of successful and failed integrations
C
	NSUCCESS=0
	NFAIL=0
C
C Number of successful integrations in different modes
C
	NFULL=0		! full peak mode
	NWEAK1=0	! sig(i)/i mode with no overlap
	NWEAK2=0	! sig(i)/i mode with peak-peak overlap
	NWEAK3=0	! sig(i)/i mode with core-peak overlap
	NWEAK4=0	! sig(i)/i mode with core-core overlap
	NPAIR=0		! overlapped twin pair integration
C
C Number of various integration errors
C
	NEDGE=0		! edges/holes
	NNOBACK=0	! background points
	NPIXEL=0	! pixel overload
	NNEIGH=0	! strong neighbour overlap
	NFAR=0		! spot off-center
	NELLIP=0	! bad calculated ellipse shape
C
C Make all integrations look like failures
C
	DO IR=1,NREFS
	  RINTINT(IR)=-9999.0
	  RSIGMA(IR)=-9999.0
	ENDDO
C
	RETURN
	END


	SUBROUTINE INTEGRATE_SPOT(IR, RAD_MULT,TOLER_NBOUR, DPFRAC,PFRAC_OBS,PFRAC_SIG)
C
	REAL RAD_MULT(3)
C
C Performs integrations of twinned or non-twinned spots
C
C Sets IOPTION(IR) & IOVER(IR) to show the type of success or failure
C Successes:
C      IOPTION=3  IOVER=1-2	Full integration, strong peak
C      IOPTION=2  IOVER=1-4	sigi/i method, due to spot overlap
C      IOPTION=2  IOVER=5  	overlapping spot method for twin pair
C      IOPTION=1  IOVER=1-4	sigi/i method, due to weak intensity
C      IOPTION=1  IOVER=1
C      IOPTION=1  IOVER=1
C Failures:
C      IOPTION=-1	 Spot overlap of an edge, hole etc.
C      IOPTION=-2	 Pixel intensity overload
C      IOPTION=-3	 Significant overlap with neighbour
C      IOPTION=-4	 Invalid ellipse shape or background
C IOVER contains the worst overlap type: 1=none, 2=P-P, 3=C-P, 4=C-C
C
	COMMON /OPTIONS_COM/ IOPT_HKL,NHKLM
C
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /SPOT2_COM/ IHKLM(4,20000),WAV(20000),TTH(20000),MULT(20000)
	COMMON /BGLOCAL_COM/ BKG(20000),BKG_ESD(20000),BKG_VAR(20000)
	COMMON /INTINT_COM/ RINTINT(20000),RSIGMA(20000),IOPTION(20000),IOVER(20000)
C
	COMMON /MODELS_PFRAC_COM/ PFRAC_TARGET,PFRAC_ERROR,DPFRAC_ZERO,DPFRAC_GRAD
	COMMON /GEOM_COM/ PIX_CEN(2),PIX_SIZE(2),DRUM_RAD
	COMMON /LIB_COM/ PLIB(1000,6),ILIB(1000),NLIB
C
C Counters for success/fail integrations, success modes and fail types
C
	COMMON /INT_COUNTERS_COM/ NSUCCESS,NFAIL,NFULL,NWEAK1,NWEAK2,NWEAK3,
	1			NWEAK4, NPAIR,NEDGE,NNOBACK,NPIXEL,NNEIGH,NFAR,NELLIP
C
C Flag for an overlapped twin pair
C
	LOGICAL LOVERLAP
C
C Calculate fitted ellipse parameters and peak fraction for (X,Y)
C
	CALL CALC_FITTED_EFH(E,F,H, X(IR),Y(IR))
	PFRAC=PFRAC_TARGET
	DPFRAC=DPFRAC_ZERO+DPFRAC_GRAD*SQRT( (X(IR)-PIX_CEN(1))**2 + (Y(IR)-PIX_CEN(2))**2 )
C
C Output more information on spot
C
	CALL CALC_EFH2AXES(A,B,ANGLE, E,F,H)
	WRITE(30,'(A,2P,3F6.2,0P,2X,A,2F5.1,I4)',IOSTAT=IDUMMY)
	1		'EFH(*100)',E,F,H,'A,B,Angle',A,B,NINT(ANGLE*57.3)
	WRITE(30,'(A,F7.3,A,F6.3,3X,A,F8.2,A,F6.2)',IOSTAT=IDUMMY)
	1	'Peak Frac',PFRAC,' +/-', DPFRAC,
	2	'Background',BKG(IR),' +/-',BKG_ESD(IR)
C
C Spot rejected if it includes edges, scratches, etc.
C
	CALL REJECT_ELLIPSE(IREJECT, X(IR),Y(IR),E,F,H,RAD_MULT(2))
C
C Do special cases for twins
C
	LOVERLAP=.FALSE.
	IF(IOPT_HKL .EQ. 2) THEN
	  IM=IHKLM(4,IR)
	  ITWIN=IM - (IM/10)*10
	  LOVERLAP=(IM .GE. 40)
C If twin 1 of an overlapping twin, also test twin 2 for edges, etc.
	  IF(LOVERLAP .AND. IREJECT.EQ.0) CALL REJECT_ELLIPSE(IREJECT, X(IR+1),Y(IR+1),E,F,H,RAD_MULT(2))
C If still not rejected, do more twin tests
		IF(IREJECT .EQ. 0) THEN
		  CALL INTEGRATE_TWIN_CHECK(LOVERLAP, E,F,H,  IM,ITWIN,IR)
C If twin 2 of an overlapped pair, skip pair integration and return
	    IF(LOVERLAP .AND. ITWIN.EQ.2) THEN
			  WRITE(30,'(A)') 'Results shown in previous integration'
				RETURN
			ENDIF
		ENDIF
	ENDIF
C
C If rejected, log reason, set flag, update *.ell file, then return
C Increment the counters twice if an overlapping twin
C
	IR2=IR
	IF( LOVERLAP ) IR2=IR+1
	IF(IREJECT .EQ. -1) THEN
	  WRITE(30,'(A)') 'Integration failed: Invalid ellipse shape'
	  DO I=IR,IR2
	    NELLIP=NELLIP+1
	    IOPTION(I)=-4
	    CALL OUTPUT_ELL_RECORD(I, E,F,H,RAD_MULT)
	  ENDDO
	  RETURN
	ELSEIF(IREJECT .EQ. 1) THEN
	  WRITE(30,'(A)') 'Integration failed: Spot overlaps edge, hole etc'
	  DO I=IR,IR2
	    NEDGE=NEDGE+1
	    IOPTION(I)=-1
	    CALL OUTPUT_ELL_RECORD(I, E,F,H,RAD_MULT)
	  ENDDO
	  RETURN
	ENDIF
C
C Estimate spot overlap, ROVER & IOVER, needed for CALC_INTINTS().
C Also load BKG(IR) and BKG_ESD(IR).
C NB:  Overlapped twin is ignored if M>50 (LOVERLAP is true).
C
	CALL CALC_NEIGH_OVERLAP(ROVER,IOVER(IR), IR, E,F,H, PFRAC,RAD_MULT, .TRUE.)
	IF(BKG_ESD(IR) .LE. 0.0) THEN
	  DO I=IR,IR2
	    NNOBACK=NNOBACK+1
	    IOPTION(I)=-4
	    CALL OUTPUT_ELL_RECORD(I, E,F,H,RAD_MULT)
	  ENDDO
	  RETURN
	ENDIF
C
C Calculate the integrated intensity for overlapping pairs of twins
C
	IF( LOVERLAP ) THEN
	  CALL CALC_TWIN_INTINTS(IR, E,F,H,PFRAC,DPFRAC)
	  DO I=IR,IR2
	    NPAIR=NPAIR+1
	    IOPTION(I)=2
	    IOVER(I)=5
	    CALL OUTPUT_ELL_RECORD(I, E,F,H,RAD_MULT)
	  ENDDO
	  RETURN
	ENDIF
C
C Calculate the integrated intensity and output results to the log file
C
	CALL CALC_INT_INT(RINTINT(IR),RSIGMA(IR),ROVER,IOVER(IR),ISTATUS,
	1                  PFRAC_OBS,PFRAC_SIG, X(IR),Y(IR),BKG_ESD(IR),
	2                  E,F,H,RAD_MULT, PFRAC,DPFRAC,  TOLER_NBOUR)
C
C Decode ISTATUS & IOVER(IR) and update counters & IOPTION(IR)
C
	IF(ISTATUS .EQ. 0) THEN
	  IOPTION(IR)=3
	  NFULL=NFULL+1
	ELSEIF(ISTATUS .GT. 0) THEN
	  IOPTION(IR)=ISTATUS
	  IF(IOVER(IR) .EQ. 1) NWEAK1=NWEAK1+1
	  IF(IOVER(IR) .EQ. 2) NWEAK2=NWEAK2+1
	  IF(IOVER(IR) .EQ. 3) NWEAK3=NWEAK3+1
	  IF(IOVER(IR) .EQ. 4) NWEAK4=NWEAK4+1
	ELSEIF(ISTATUS .EQ. -1) THEN
	  IOPTION(IR)=-1
	  NEDGE=NEDGE+1
	  WRITE(30,'(A)') 'Integration failed: Spot overlaps edge, hole etc'
	ELSEIF(ISTATUS .EQ. -2) THEN
	  IOPTION(IR)=-3
	  NNEIGH=NNEIGH+1
	  WRITE(30,'(A)') 'Integration failed: Significant overlap with neighbour'
	ELSEIF(ISTATUS .EQ. -3) THEN
	  IOPTION(IR)=-2
	  NPIXEL=NPIXEL+1
	  WRITE(30,'(A)') 'Integration failed: Pixel intensity overload'
	ENDIF
C
C Add spot shape info for *.ell file
C
	CALL OUTPUT_ELL_RECORD(IR, E,F,H,RAD_MULT)
C
	RETURN
	END


	SUBROUTINE CALC_INT_INT(RINTINT,RSIGMA,ROVER,IOVER,ISTATUS,
	1                        PFRAC_OBS,PFRAC_SIG, X,Y,BKG_ESD,
	2											  E,F,H,RAD_MULT, PFRAC,DPFRAC,  TOLER_NBOUR)
C
	REAL RAD_MULT(3)
C
C Performs full peak integration or sigi/i method to integrate single spots
C depending on the relative sizes of peak esd's and spot overlaps
C
C ISTATUS = 0	Success: full integration
C		    1	Success: sigi/i due to weakness
C		    2	Success: sigi/i due to overlap
C          -1	Failure: ellipse includes an edge, hole etc.
C          -2	Failure: overlap with strong neighbour
C          -3	Failure: Pixel overload
C IOVER returns the worst overlap type: 1=none, 2=P-P, 3=C-P, 4=C-C
C
	COMMON /GAIN_COM/ GAIN,ESD_FACTOR,CROSSTALK,CROSS_ELLIPSE,CROSS_SMOOTH
C
	LOGICAL LSTRONG
	INTEGER NSUM(2)
	REAL CSUM(2),BKG_AVE(2)
C
C Skip integration if peak ellipse includes edge, holes, etc.
C
	CALL REJECT_ELLIPSE(IREJ, X,Y,E,F,H,RAD_MULT(2))
	IF(IREJ .NE. 0) THEN
	  ISTATUS=-1
	  RETURN
	ENDIF
C
C Check for pixel intensity overload, which occurs much less than 65K
C
	IF(GET_ELLIPSE_MAXINT(X,Y,E,F,H) .GT. 6.0E4) THEN
	  ISTATUS=-3
	  RETURN
	ENDIF
C
C Skip integration if core-core or core-peak overlap, and the intensity
C overlap from the neighbours is above the cutoff
C
	IF(IOVER .GE. 3) THEN
	  IF(ROVER .GT. TOLER_NBOUR) THEN
	    ISTATUS=-2
	    RETURN
	  ENDIF
	ENDIF
C
C Sum counts in the core/peak ellipses
C
	CALL SUM_ELLIPSE_INTS(CSUM,NSUM, X,Y,E,F,H, RAD_MULT,2)
C
C Get average background counts for the weak and strong integration areas
C
	CALL AVE_ELLIPSE_BKGS(BKG_AVE,X,Y,E,F,H, RAD_MULT,2)
C
C Calculate int.ints using both the weak and strong methods
C
	CALL CALC_INT_INT_WEAK(RINT1,RSIG1, CSUM(1),NSUM(1), BKG_AVE(1),BKG_ESD, PFRAC,DPFRAC)
	CALL CALC_INT_INT_STRONG(RINT2,RSIG2, CSUM(2),NSUM(2), BKG_AVE(2),BKG_ESD)
C
	WRITE(30,'(A,I6,A,I6,A)',IOSTAT=IDUMMY)
     1				'Int. Int. = ',NINT(RINT1),' (core),',
	2				NINT(RINT2),' (core+peak)'
	WRITE(30,'(A,I7,A,I7,A)',IOSTAT=IDUMMY)
     1				'Esds for Int. Int. = ',NINT(RSIG1),' (core),',
	2				NINT(RSIG2),' (core+peak)'
C Try the strong integration if its esd estimate is less
	LSTRONG=(RSIG2 .LT. RSIG1)
C
C If strong with non, or minor, overlap, do the full integration
C
	IF(LSTRONG .AND. IOVER.LE.2) THEN
C Use the strong (or full integration) method
	  RINTINT=RINT2
	  RSIGMA=RSIG2
C Set STATUS=0 for full integration
	  ISTATUS=0
C Output integration mode and any overlap
	  IF(IOVER .EQ. 1) THEN
	    WRITE(30,'(A)') 'Full integration, strong peak + no overlap'
	  ELSE
	    WRITE(30,'(A)') 'Full integration, strong peak + peak-peak overlap'
	  ENDIF
C
C Otherwise, use the sigI/I integration method
C
	ELSE
C Use the weak (or sigI/I integration) method
	  RINTINT=RINT1
	  RSIGMA=RSIG1
C Set STATUS=1 (or 2 if downgraded from STRONG) for sigI/I integration
	  ISTATUS=1
	  IF( LSTRONG ) THEN
	    ISTATUS=2
	    WRITE(30,'(2X,A)') 'Downgrade from full integration due to overlap'
	  ENDIF
C Output integration type and any overlap
	  IF(IOVER .EQ. 1) WRITE(30,'(A)') 'sig(I)/I integration, no overlap'
	  IF(IOVER .EQ. 2) WRITE(30,'(A)') 'sig(I)/I integration, peak-peak overlap'
	  IF(IOVER .EQ. 3) WRITE(30,'(A)') 'sig(I)/I integration, core-peak overlap + weak neighbour'
	  IF(IOVER .EQ. 4) WRITE(30,'(A)') 'sig(I)/I integration, core-core overlap + weak neighbour'
C
	ENDIF
C
	PFRAC_OBS=0.0
	PFRAC_SIG=1E3
	IF(RINT1.GT.0.0 .AND. RINT2.GT.0.0) THEN
	  PFRAC_OBS=PFRAC*RINT1/RINT2
	  PFRAC_SIG=PFRAC_OBS*SQRT( (RSIG1/RINT1)**2 + (RSIG2/RINT2)**2 )
	ENDIF
C
C Output ellipse sums and the observed peak-fraction
C
	  WRITE(30,'(2X,2(A,I8,A,I4),A,F5.2)',IOSTAT=IDUMMY)
	1		'Counts(N): Core',NINT(CSUM(1)),'(',NSUM(1),') Full',
	2		NINT(CSUM(2)),'(',NSUM(2),')  Pfrac',PFRAC_OBS
C
	RETURN
	END


	SUBROUTINE CALC_INT_INT_STRONG(RINTINT,RSIGMA, CSUM,NSUM, BKG,BKG_ESD)
C
	COMMON /GAIN_COM/ GAIN,ESD_FACTOR,CROSSTALK,CROSS_ELLIPSE,CROSS_SMOOTH
C
C Calculate integrated intensity and esd for the peak ellipse
C
C The summation of the "peak area" will be affected by the crosstalk term, but only
C near the edges which we assume have the same intensity as the background level.
C The correction to the esd of CSUM(2) will be the same value as if the peak
C area has the same average intensity as the background, ie. TOTI(2)=0.
C The counting statistics are multiplied by the detector gain factor ESD_FACTOR.
C
	CSUM_BKG=NSUM*BKG
	RINTINT=CSUM-CSUM_BKG
C
	RAD_PEAK=SQRT(NSUM/3.14159)
	ESD_CORR=(1.0-CROSS_ELLIPSE/RAD_PEAK)
	CSUM_CORR=MAX(0.0,RINTINT) + CSUM_BKG*ESD_CORR**2
	RSIGMA=SQRT( CSUM_CORR*ESD_FACTOR**2 + (NSUM*BKG_ESD)**2 )
C
	RETURN
	END


	SUBROUTINE CALC_INT_INT_WEAK(RINTINT,RSIGMA, CSUM,NSUM, BKG,BKG_ESD, PFRAC,DPFRAC)
C
	COMMON /GAIN_COM/ GAIN,ESD_FACTOR,CROSSTALK,CROSS_ELLIPSE,CROSS_SMOOTH
C
C Calculate the integrated intensity from the core ellipse
C
	CSUM_BKG=NSUM*BKG
	RINTINT=CSUM-CSUM_BKG
C
C Calculate the esd of PPK from counting statistics of the peak and the error
C in background contribution. The counting statistics are multiplied by the
C detector gain factor ESD_FACTOR.
C
C The esd of the "core area" will be affected by crosstalk effects, but only near
C the edges which we assume have the same intensity as the background level.
C Any intensity above the background level in CSUM(1) is assumed to be from
C the peak and will not be affected by crosstalk.
C
	RAD_PEAK=MAX(1.0, SQRT(NSUM/3.14159) )
	ESD_CORR=1.0 - CROSS_ELLIPSE/RAD_PEAK
	CSUM_CORR=MAX(0.0,RINTINT) + CSUM_BKG*ESD_CORR**2
	PFRAC_CORR=MAX(0.0,RINTINT)*DPFRAC/PFRAC
	RSIGMA=SQRT( CSUM_CORR*ESD_FACTOR**2 + (NSUM*BKG_ESD)**2 + PFRAC_CORR**2 )
C
C Scale RINTINT & RSIGMA BY peak-fraction
C
	RINTINT=RINTINT/PFRAC
	RSIGMA =RSIGMA /PFRAC
C
	RETURN
	END
