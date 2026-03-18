C ========================== Routines in this file ===============================
C ---------  Routines to add eligible model ellipses	--------
C	SUBROUTINE ADD_MODEL_SPOT(A_MAX,FFILL,TOLER_NBOUR,CONTOUR,RAD_MULT,ZCUTOFFS, IR)
C	SUBROUTINE MAKE_HKLM_IDENT(SHKLM, IR,HKLM)
C	SUBROUTINE PRINT_SPACER(LSPACER)
C --------- Routines for massaging model ellipses --------
C	SUBROUTINE FIT_MODEL_PARAMS(RAD_MULT,TOLER_NBOUR)
C	SUBROUTINE FIT_PEAK_FRAC(EFH_SCALE,PFRAC1,DPFRAC1,PFRAC2,DPFRAC2,
C	1						RAD_MULT,TOLER_NBOUR,PFRAC_TARGET, IMODEL)
C	SUBROUTINE SET_PEAK_FRAC(RAD_MULT,PFRAC_TARGET,TOLER_NBOUR)
C	SUBROUTINE MEASURE_PEAK_FRAC(PEAK_FRAC,DPEAK_FRAC, IR,E,F,H,RAD,TOLER_NBOUR)
C	SUBROUTINE CALC_EFH_SCALE(EFH_SCALE, XCEN,YCEN,E,F,H,RAD_MULT, PFRAC_TARGET)
C ---------- Reduce number of model spots in zones ---------
C	SUBROUTINE PRUNE_MODELS(NWANT_CEN,NWANT_OUT)
C ===============================================================================

C ---------  Routines to add eligible model ellipses	--------
	SUBROUTINE ADD_MODEL_SPOT(A_MAX,FFILL,TOLER_NBOUR,CONTOUR,RAD_MULT,ZCUTOFFS, IR)
C
	REAL RAD_MULT(3),ZCUTOFFS(5)
C
	COMMON /RIMAGE_COM/ RIMAGE(8000,2500)
	COMMON /IZONES_COM/ IZONES(20000)
	COMMON /LIB_COM/ PLIB(1000,6),ILIB(1000),NLIB
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /SPOT2_COM/ IHKLM(4,20000),WAV(20000),TTH(20000),MULT(20000)
	COMMON /BGLOCAL_COM/ BKG(20000),BKG_ESD(20000),BKG_VAR(20000)
C
	COMMON /OPTIONS_COM/ IOPT_HKL,NHKLM
C
	LOGICAL LSPACER
	SAVE LSPACER
C
	CHARACTER SHKLM*34
	INTEGER NSUM(3)
	REAL CSUM(3)
C
C Set LSPACER to prevent an initial print out of the spacer line
C
	IF(IR .EQ. 1) LSPACER=.TRUE.
C
C If spot too weak, skip to next spot without	any output
C
	PKBG=(RMAX(IR)-BKG(IR))/SQRT(MAX(10.0, BKG(IR) ))
	IF(PKBG .LT. ZCUTOFFS(IZONES(IR))) RETURN
C If satellite reflection, skip to next spot without any output
	IF(IOPT_HKL.EQ.3 .AND. IHKLM(4,IR).NE.0) RETURN
C
C Create identifier string for spot with optional twin or modulation number
C
	CALL MAKE_HKLM_IDENT(SHKLM, IR,IHKLM(1,IR))
C
C Reject various twin pairs that are likely to have a distorted shape
C
	IF(IOPT_HKL.EQ.2 .AND. IHKLM(4,IR).EQ.52) THEN
	  WRITE(30,'(2A)') TRIM(SHKLM),' rejected as second of an overlapped twin'
	  LSPACER=.FALSE.
	  RETURN
	ELSEIF(IOPT_HKL.EQ.2 .AND. IHKLM(4,IR).GT.60) THEN
	  WRITE(30,'(2A)') TRIM(SHKLM),' rejected as a partially overlapped twin'
	  LSPACER=.FALSE.
	  RETURN
	ENDIF
C
C Fit a contour ellipse to the spot
C
	CALL CALC_CONTOUR_SHAPE(E,F,H, FFRAC,A,B,ANGLE, DXCEN,DYCEN,
	1							SHKLM, CONTOUR,FFILL,A_MAX, IR,ISTATUS)
C Reject on failure
	IF(ISTATUS .EQ. 1) THEN
	  LSPACER=.FALSE.
	  RETURN
	ENDIF
C
C Reject if pixel intensity overload (occurs << 65K)
C
	IF(GET_ELLIPSE_MAXINT(X(IR),Y(IR),E,F,H) .GT. 6.0E4) THEN
	  WRITE(30,'(2A)') TRIM(SHKLM),' rejected due to intensity overload'
	  LSPACER=.FALSE.
	  RETURN
	ENDIF
C
C The spot may be OK, so output details to the log file
C
	IX=NINT(X(IR))
	IY=NINT(Y(IR))
C Output a spacer line if we don't already have one
	IF( .NOT.LSPACER ) WRITE(30,'(38(1H*))')
C Output spot identifier plus XY & wav
	WRITE(30,'(2A,2(I4,A),F5.2)',IOSTAT=IDUMMY) TRIM(SHKLM),
	1	'  XY (',IX,',',IY,')  Wav',WAV(IR)
C Output if the first spot of a closely overlapped twin
	  IF(IOPT_HKL.EQ.2 .AND. IHKLM(4,IR).EQ.51)
	1	WRITE(30,'(A)') 'Spot accepted as twin 1 of a closely overlapped pair'
C Output intensities and background
	  WRITE(30,'(A,I6,A,I6,A,I5,A,F6.1,A)',IOSTAT=IDUM)
	1	'Intensities',NINT(RIMAGE(IX,IY)),'(cen)',
	2	NINT(RMAX(IR)),'(ave)',NINT(BKG(IR)),'(bkg)',PKBG,'(Peak/esd-bkg)'
C Output the ellipse information
	WRITE(30,'(A,2F5.1,2X,A,2F5.1,I4,2X,A,F6.3)',IOSTAT=IDUM)
	1          'dX,dY',DXCEN,DYCEN,'Ell. axes & angle',
	2					A,B,NINT(ANGLE*57.296),'Ffrac',FFRAC
	WRITE(30,'(A,3F7.3)') 'Core ellipse parameters(*100)',
	1				  E*100.0,F*100.0,H*100.0
C
C Reject and print spacer if peak ellipse includes edge, holes, etc.
C
	CALL REJECT_ELLIPSE(ISTATUS, X(IR),Y(IR),E,F,H,RAD_MULT(2))
	IF(ISTATUS .NE. 0) THEN
	  WRITE(30,'(A)') 'Spot rejected as ellipse overlaps edges, holes, etc.'
	  CALL PRINT_SPACER(LSPACER)
	  RETURN
	ENDIF
C
C Calculate degree of spot overlap
C
	CALL CALC_NEIGH_OVERLAP(ROVER,IOVER2, IR, E,F,H, 0.8,RAD_MULT, .TRUE.)
C Reject and print spacer if any core-core overlap
	IF(IOVER2 .EQ. 4) THEN
	  WRITE(30,'(A)') 'Spot rejected due to core-core overlap'
	  CALL PRINT_SPACER(LSPACER)
	  RETURN
	ENDIF
C Reject and print spacer if any strong neighbour overlap
	IF(ROVER .GT. 0.5*TOLER_NBOUR) THEN
	  WRITE(30,'(A)') 'Spot rejected due to strong neighbour'
	  CALL PRINT_SPACER(LSPACER)
	  RETURN
	ENDIF
C
C Calculate and check peak fraction
C Sum intensities in core and peak ellipses
C
	CALL SUM_ELLIPSE_INTS(CSUM,NSUM, X(IR),Y(IR),E,F,H,RAD_MULT,2)
C Calculate and output core and peak intensities
	CORE=CSUM(1)-NSUM(1)*BKG(IR)
	PEAK=CSUM(2)-NSUM(2)*BKG(IR) - CORE
	DCORE=SQRT( CSUM(1) + (BKG_ESD(IR)*NSUM(1))**2 )
	DPEAK=SQRT( ABS(CSUM(2)-CSUM(1)) + (BKG_ESD(IR)*(NSUM(2)-NSUM(1)))**2 )
	WRITE(30,'(2(A,I8,I5,2X),A,F7.1,1X,A,F5.1)')
	1         'Core',NINT(CSUM(1)),NSUM(1),
	2         'Full peak',NINT(CSUM(2)),NSUM(2),
	3         'bkg',BKG(IR),'+/-',BKG_ESD(IR)
C Calculate and output peak fraction
	PFRAC=CORE/(CORE+PEAK)
	DPFRAC=DPEAK/MAX(DPEAK/4.0,CORE+PEAK)		! approximation
	WRITE(30,'(A,F6.3,1X,A,F6.3)',IOSTAT=IDUM) 'Peak frac',PFRAC,'+/-',DPFRAC
C Reject and print spacer if a ridiculous PFRAC value
	IF(PFRAC.LE.0.5 .OR. PFRAC.GE.1.1) THEN
	  WRITE(30,'(A)') 'Spot rejected as peak fraction not within 0.5 to 1.1'
	  CALL PRINT_SPACER(LSPACER)
	  RETURN
	ELSEIF(DPFRAC .GE. 0.3) THEN
	  WRITE(30,'(A)') 'Spot rejected as peak fraction esd > 0.3'
	  CALL PRINT_SPACER(LSPACER)
	  RETURN
	ENDIF
C
C If there is space, add model to model library, else output a warning
C
	IF(NLIB .GE. 1000) THEN
	  WRITE(30,'(A)') 'WARNING: Spot not stored, model list full'
	ELSE
	  NLIB=NLIB+1
	  ILIB(NLIB)  =IR
	  PLIB(NLIB,1)=X(IR)
	  PLIB(NLIB,2)=Y(IR)
	  PLIB(NLIB,3)=E
	  PLIB(NLIB,4)=F
	  PLIB(NLIB,5)=H
	  PLIB(NLIB,6)=PFRAC
	  WRITE(30,'(A,I4)') 'Spot stored as potential model',NLIB
	ENDIF
C
	CALL PRINT_SPACER(LSPACER)
C
	RETURN
	END


	SUBROUTINE MAKE_HKLM_IDENT(SHKLM, IR,HKLM)
C
C Create an identifier string containing IR, HKL, and an optional
C twin or modulation number
C
	CHARACTER SHKLM*(*)
	INTEGER HKLM(4)
C
	COMMON /OPTIONS_COM/ IOPT_HKL,NHKLM
C
C Create identifier string for spot with HKL and optional twin or modulation index
	WRITE(SHKLM,'(A4,I5,A6,3I4)',IOSTAT=IDUMMY) 'Spot',IR,' HKL',(HKLM(K),K=1,3)
	IF(IOPT_HKL .EQ. 2) THEN
	  IM=HKLM(4)-(HKLM(4)/10)*10
	  WRITE(SHKLM,'(A,A5,I2)',IOSTAT=IDUMMY) TRIM(SHKLM),' Twin',IM
	ELSEIF(IOPT_HKL .EQ. 3) THEN
	  WRITE(SHKLM,'(A,A5,I2)',IOSTAT=IDUMMY) TRIM(SHKLM),' Sat.',HKLM(4)
	ENDIF
C
	RETURN
	END


	SUBROUTINE PRINT_SPACER(LSPACER)
C
	LOGICAL LSPACER
C
	WRITE(30,'(38(1H*))')
	LSPACER=.TRUE.
C
	RETURN
	END


C --------- Routines for massaging model ellipses --------

	SUBROUTINE FIT_MODEL_PARAMS(RAD_MULT,TOLER_NBOUR)
C
	REAL RAD_MULT(3)
C
	COMMON /LIB_COM/ PLIB(1000,6),ILIB(1000),NLIB
C
	COMMON /MODELS_PFRAC_COM/ PFRAC_TARGET,PFRAC_ERROR,DPFRAC_ZERO,DPFRAC_GRAD
C
	COMMON /GEOM_COM/ PIX_CEN(2),PIX_SIZE(2),DRUM_RAD
C
	REAL EFHP(4,1000)
	REAL ELL(3,0:10,0:10)
C
	NITERS=3
	DO ITER=1,NITERS
C
	WRITE(30,'(/,A,I2)') '=== Optimising ellipse shape: Iteration',ITER
C
C Sanity check on library size
C
	IF(NLIB .LT. 10) CALL QUIT('ERROR: Less than 10 model spots in table')
C
C Save the E,F,H,Pfrac for all models
C
	DO I=1,NLIB
	  EFHP(1,I)=PLIB(I,3)
	  EFHP(2,I)=PLIB(I,4)
	  EFHP(3,I)=PLIB(I,5)
	  EFHP(4,I)=PLIB(I,6)
	ENDDO
C
C Fit model ellipses to sample shape and mosaicity. The fitted
C values are stored in /FIT_EFH_COM/.
C
	CALL FIT_EFH(PLIB,NLIB)
C
	WRITE(30,'(/,A)') 'Model   X    Y      EFH*100(initial)      EFH*100(fitted)'
	DO I=1,NLIB
	  WRITE(30,'(I4,I6,I5,2P,2(1X,3F7.2),0P)',IOSTAT=IDUMMY) I,
	1				(NINT(PLIB(I,K)),K=1,2),(EFHP(K,I),K=1,3),(PLIB(I,K),K=3,5)
	ENDDO
C
C Sanity check for unphysical ellipses
C
	DO I=1,NLIB
	  E=PLIB(I,3)
	  F=PLIB(I,4)
	  H=PLIB(I,5)
	  IF(MIN(E,F) .LE. 0.0) CALL QUIT('BUG(fit_all_models): Invalid E,F')
	  IF(E*F .LE. H**2) CALL QUIT('BUG(fit_all_models): Invalid H')
	ENDDO
C
C Loop over all models, recalculating Pfrac values in PLIB()
C
	WRITE(30,'(/,A,F6.2)') 'Scaling ellipse size to achieve peak fraction =',PFRAC_TARGET
	WRITE(30,'(/,A)') 'Model  X    Y     P-Frac(initial,fitted)  EFH scale factor'
C
	DO I=1,NLIB
	  CALL FIT_PEAK_FRAC(EFH_SCALE,PFRAC1,DPFRAC1,PFRAC2,DPFRAC2,
	1						RAD_MULT,TOLER_NBOUR,PFRAC_TARGET, I)
C If successful output results and scale EFH
C Otherwise, output a warning and set ILIB=0 to indicate failure
	  IF(EFH_SCALE.GE.0.5 .AND. EFH_SCALE.LE.2.0) THEN
	    PLIB(I,3)=PLIB(I,3)*EFH_SCALE
	    PLIB(I,4)=PLIB(I,4)*EFH_SCALE
	    PLIB(I,5)=PLIB(I,5)*EFH_SCALE
	    PLIB(I,6)=PFRAC2
	    WRITE(30,'(I4,I6,I5,8X,2F8.3,F11.2)') I,(NINT(PLIB(I,K)),K=1,2),
	1										PFRAC1,PFRAC2,EFH_SCALE
	  ELSE
	    WRITE(30,'(I4,I6,I5,11X,A)') I,(NINT(PLIB(I,K)),K=1,2),
	1				'........................ Scaling failed'
	    ILIB(I)=0
	  ENDIF
	ENDDO
C
C Update library by skipping any ILIB=0
C
	N=0
	DO I=1,NLIB
	  IF(ILIB(I) .EQ. 0) CYCLE
	  N=N+1
	  ILIB(N)=ILIB(I)
	  DO I2=1,6
	    PLIB(N,I2)=PLIB(I,I2)
	  ENDDO
	ENDDO
	NLIB=N
C
	ENDDO
C
C Loop over all models, recalculating Pfrac values in PLIB()
C
	WRITE(30,'(/,A)') '=== Optimising ellipse shape: Final Results'
	WRITE(30,'(/,A)') 'Peak fraction versus radial distance in pixels'
	WRITE(30,'(A)') 'Model  X    Y     P-Frac(final,  esd)   Radial Distance'
	NSUM=0
	SUMX=0.0
	SUMY=0.0
	SUMXY=0.0
	SUMX2=0.0
C
	DO I=1,NLIB
C
	  CALL MEASURE_PEAK_FRAC(PFRAC2,DPFRAC2, ILIB(I),
	1		PLIB(I,3),PLIB(I,4),PLIB(I,5), RAD_MULT,TOLER_NBOUR)
	  IF(ABS(PLIB(I,6)-PFRAC2) .GT. 0.01)
	1			CALL QUIT('BUG(fit_model_params): FRAC2 invalid')
	  RDIST=SQRT( (PLIB(I,1)-PIX_CEN(1))**2 + (PLIB(I,2)-PIX_CEN(2))**2 )
	  WRITE(30,'(I4,I6,I5,F15.3,F7.3,I10)') I,(NINT(PLIB(I,K)),K=1,2),
	1											PFRAC2,DPFRAC2,NINT(RDIST)
	  DIFF=ABS(PFRAC1-PFRAC_TARGET)
	  NSUM=NSUM+1
	  SUMX=SUMX+RDIST
	  SUMY=SUMY+DIFF
	  SUMXY=SUMXY+RDIST*DIFF
	  SUMX2=SUMX2+RDIST**2
	ENDDO
C
	AVEX=SUMX/NSUM
	AVEY=SUMY/NSUM
	DPFRAC_GRAD=SUMXY/SUMX2
	DPFRAC_ZERO=AVEY-DPFRAC_GRAD*AVEX
	WRITE(30,'(/,A)') 'Statistical analysis of P-Frac error versus Radial Distance'
	WRITE(30,'(A,2(F5.2,A))') ' | P-Frac -',PFRAC_TARGET,' | =',AVEY,' (average)'
	WRITE(30,'(19X,A,2(F6.3,A))') '=',DPFRAC_ZERO,' +',DPFRAC_GRAD*1000.0,' * R-DIST / 1000 (fit)'
	IF(DPFRAC_GRAD .LT. 0) THEN
	  DPFRAC_GRAD=0.0
	  DPFRAC_ZERO=AVEX
	ENDIF
C
	DPFRAC_GRAD=DPFRAC_GRAD*PFRAC_ERROR
	DPFRAC_ZERO=MAX(0.01, DPFRAC_ZERO*PFRAC_ERROR )
C
	WRITE(30,'(/,A,F6.3)')    'Setting peak fraction =',PFRAC_TARGET
	WRITE(30,'(A,2(F6.3,A))') '             with esd =',DPFRAC_ZERO,' +',
	1									DPFRAC_GRAD*1000.0,' * R-DIST / 1000'
C
C Make grids of 11 x 11 for A,B,ANGLE calculated for across X,Y.
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
	ISTEPX=NUMX/10
	ISTEPY=NUMY/10
C
	DO IX=0,10
	  X0=IX*ISTEPX
	  DO IY=0,10
	    Y0=IY*ISTEPY
	    CALL CALC_FITTED_EFH(E0,F0,H0, X0,Y0)
	    CALL CALC_EFH2AXES(ELL(1,IX,IY),ELL(2,IX,IY),ELL(3,IX,IY), E0,F0,H0)
	  ENDDO
	ENDDO
C
C Output calculated values in tables to the log file
C
	WRITE(30,'(/,A)') 'Calculated spot ellipse A for X(horiz) by Y(vert)'
	WRITE(30,'(I10,1X,10I6)',IOSTAT=IDUMMY) (K*ISTEPX,K=0,10)
	DO IY=0,10
	  WRITE(30,'(I5,21F6.1)',IOSTAT=IDUMMY)
	1				IY*ISTEPY,(ELL(1,K,IY),K=0,10)
	ENDDO
C
	WRITE(30,'(/,A)') 'Calculated spot ellipse B for X(horiz) by Y(vert)'
	WRITE(30,'(I10,1X,10I6)',IOSTAT=IDUMMY) (K*ISTEPX,K=0,10)
	DO IY=0,10
	  WRITE(30,'(I5,21F6.1)',IOSTAT=IDUMMY)
	1				IY*ISTEPY,(ELL(2,K,IY),K=0,10)
	ENDDO
C
	WRITE(30,'(/,A)') 'Calculated spot ellipse ANGLE for X(horiz) by Y(vert)'
	WRITE(30,'(I10,1X,10I6)',IOSTAT=IDUMMY) (K*ISTEPX,K=0,10)
	DO IY=0,10
	  WRITE(30,'(I5,21F6.0)',IOSTAT=IDUMMY)
	1				IY*ISTEPY,(ELL(3,K,IY)*57.296,K=0,10)
	ENDDO
C
C
	RETURN
	END


	SUBROUTINE FIT_PEAK_FRAC(EFH_SCALE,PFRAC1,DPFRAC1,PFRAC2,DPFRAC2,
	1						RAD_MULT,TOLER_NBOUR,PFRAC_TARGET, IMODEL)
C
	REAL RAD_MULT(3)
C
	COMMON /LIB_COM/ PLIB(1000,6),ILIB(1000),NLIB
C
C Load spot info for model IMODEL in PLIB()
C
	IR=ILIB(IMODEL)
	X=PLIB(IMODEL,1)
	Y=PLIB(IMODEL,2)
	E=PLIB(IMODEL,3)
	F=PLIB(IMODEL,4)
	H=PLIB(IMODEL,5)
C
C Measure the peak fraction for the initial E,F,H
C
	CALL MEASURE_PEAK_FRAC(PFRAC1,DPFRAC1, IR,
	1						E,F,H, RAD_MULT,TOLER_NBOUR)
C
C Calculate the scale factor for E,F,H to get to the target peak fraction
C
	CALL CALC_EFH_SCALE(EFH_SCALE, X,Y,E,F,H,RAD_MULT, PFRAC_TARGET)
C
C If routine fails, return with EFH_SCALE < 0
C
	IF(EFH_SCALE .LE. 0.0) RETURN
C
C Measure the peak fraction for the scaled E,F,H
C
	E=E*EFH_SCALE
	F=F*EFH_SCALE
	H=H*EFH_SCALE
	CALL MEASURE_PEAK_FRAC(PFRAC2,DPFRAC2, IR,
	1						E,F,H, RAD_MULT,TOLER_NBOUR)
C
	RETURN
	END


	SUBROUTINE MEASURE_PEAK_FRAC(PEAK_FRAC,DPEAK_FRAC, IR,E,F,H,RAD,TOLER_NBOUR)
C
	REAL RAD(3)
C
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /BGLOCAL_COM/ BKG(20000),BKG_ESD(20000),BKG_VAR(20000)
C
	INTEGER NSUM(3)
	REAL CSUM(3)
C
C Set return values to indicate failure (at this point)
C
	PEAK_FRAC=-1.0
	DPEAK_FRAC=-1.0
C
C If peak ellipse contains edges, etc., return as a failure
C
	CALL REJECT_ELLIPSE(IOUT,X(IR),Y(IR),E,F,H,RAD(2))
	IF(IOUT .EQ. 1) RETURN
C
C Calculate ROVER and IOVER as estimates of spot overlap
C Also load BKG(IR) & BKG_ESD(IR) with background and esd
C
	CALL CALC_NEIGH_OVERLAP(ROVER,IOVER, IR, E,F,H,0.8,RAD, .FALSE.)
C
C Return on failure, core-core or significant overlap
C
	IF(  (BKG_ESD(IR) .LE. 0.0) .OR. (IOVER .EQ. 4) .OR.
	1					(ROVER .GT. 0.5*TOLER_NBOUR)   ) RETURN
C
C Calculate the intensities in core & peak ellipses
C 
	CALL SUM_ELLIPSE_INTS(CSUM,NSUM, X(IR),Y(IR),E,F,H,RAD,2)
C
C Calculate core and peak intensities and esd's
C
	CORE=CSUM(1)-NSUM(1)*BKG(IR)
	PEAK=CSUM(2)-NSUM(2)*BKG(IR)
	DPEAK_ONLY=SQRT( ABS(CSUM(2)-CSUM(1)) + ( BKG_ESD(IR)*(NSUM(2)-NSUM(1)) )**2 )
C
C Calculate the peak fraction, PEAK_FRAC, and its esd, DPEAK_FRAC
C
	PEAK_FRAC=CORE/PEAK
	DPEAK_FRAC=DPEAK_ONLY/MAX(DPEAK_ONLY/4.0,PEAK)		! a fudge but good enough
C
	RETURN
	END


	SUBROUTINE SET_PEAK_FRAC(RAD_MULT,PFRAC_TARGET,TOLER_NBOUR)
C
C Scale E,F,H of model peaks to get the target peak fraction. Checks if the
C spot is valid in terms of exclusion areas and spot overlap, and removes
C any invalid spots from the model list.
C Also load esd in peak fraction into DPFRAC_LIB() for use with zones.
C
	REAL RAD_MULT(3)
C
	COMMON /IZONES_COM/ IZONES(20000),DPFRAC_LIB(1000)
	COMMON /LIB_COM/ PLIB(1000,6),ILIB(1000),NLIB
	COMMON /BGLOCAL_COM/ BKG(20000),BKG_ESD(20000),BKG_VAR(20000)
C
	INTEGER NSUM(2)
	REAL CSUM(2)
C
	WRITE(30,'(/,A,F6.2)') 'Scaling ellipse size to achieve peak fraction =',PFRAC_TARGET
	WRITE(30,'(/,A)') 'Model Zone  Xpix Ypix   PFrac  esd'
C
C Loop through all model spots
C
	N=0
	DO I=1,NLIB
C Unload ellipse info from models table
	  IR=ILIB(I)
	  X0=PLIB(I,1)
	  Y0=PLIB(I,2)
	  E=PLIB(I,3)
	  F=PLIB(I,4)
	  H=PLIB(I,5)
C Get scale factor E,F,H needed to get to the target peak fraction
	  CALL CALC_EFH_SCALE(EFH_SCALE, X0,Y0,E,F,H,RAD_MULT, PFRAC_TARGET)
C If routine failed, skip to next model peak
	  IF(EFH_SCALE .LE. 0.0) THEN
		WRITE(30,'(I4,I5,I7,I5,3X,A)') I,IZONES(IR),NINT(X0),NINT(Y0),
	1									'........... Scaling failed'
	    CYCLE
	  ENDIF
C Scale factor E,F,H to get to the target peak fraction
	  E=E*EFH_SCALE
	  F=F*EFH_SCALE
	  H=H*EFH_SCALE
C If peak ellipse contains edge, holes, etc., skip to next model spot
	  CALL REJECT_ELLIPSE(IOUT,X0,Y0,E,F,H,RAD_MULT(2))
	  IF(IOUT.EQ.1) THEN
	    WRITE(30,'(I4,I5,I7,I5,3X,A)') I,IZONES(IR),NINT(X0),NINT(Y0),
	1									'........... Edge, hole, etc.'
	    CYCLE
	  ENDIF
C If significant spot overlap, skip to next model spot
	  CALL CALC_NEIGH_OVERLAP(ROVER,IOVER, ILIB(I),	E,F,H,PFRAC_TARGET,RAD_MULT, .FALSE.)
C NB: This criterion is 2 times stricter than for peak integration
	  IF(ROVER .GT. 0.5*TOLER_NBOUR) THEN
	    WRITE(30,'(I4,I5,I7,I5,3X,A)') I,IZONES(IR),NINT(X0),NINT(Y0),
	1									'........... Strong neighbour'
	    CYCLE
	  ENDIF
C Don't allow any core-core overlap
	  IF(IOVER .EQ. 4) THEN
		WRITE(30,'(I4,I5,I7,I5,3X,A)') I,IZONES(IR),NINT(X0),NINT(Y0),
	1									'........... Core-core overlap'
	    CYCLE
	  ENDIF
C Recalculate the peak fraction
	  CALL SUM_ELLIPSE_INTS(CSUM,NSUM, X0,Y0,E,F,H,RAD_MULT,2)
	  PFRAC2=(CSUM(1)-NSUM(1)*BKG(IR))/(CSUM(2)-NSUM(2)*BKG(IR))
	  DPFRAC2=SQRT( CSUM(2)-CSUM(1) + ( (NSUM(2)-NSUM(1))*BKG_ESD(IR) )**2 ) *
	1				(CSUM(1)-NSUM(1)*BKG(IR)) / (CSUM(2)-NSUM(2)*BKG(IR))**2

C Reject if peak fraction is not within 50% of (1 - target value )
	  IF(ABS(PFRAC2-PFRAC_TARGET) .GT. 0.5*(1.0-PFRAC_TARGET)) THEN
		WRITE(30,'(I4,I5,I7,I5,F8.3,F6.3,A)') I,IZONES(IR),NINT(X0),NINT(Y0),
	1									PFRAC2,DPFRAC2,' Invalid peak fraction'
	    CYCLE
	  ENDIF
C Copy model to new location in list
	  N=N+1
	  ILIB(N)=IR
	  PLIB(N,1)=X0
	  PLIB(N,2)=Y0
	  PLIB(N,3)=E
	  PLIB(N,4)=F
	  PLIB(N,5)=H
	  PLIB(N,6)=PFRAC2
C
	  DPFRAC_LIB(N)=DPFRAC2
	  WRITE(30,'(I4,I5,I7,I5,F8.3,F6.3)') I,IZONES(IR),NINT(X0),NINT(Y0),PFRAC2,DPFRAC2
C
	ENDDO
	NLIB=N
C
	RETURN
	END


	SUBROUTINE CALC_EFH_SCALE(EFH_SCALE, XCEN,YCEN,E,F,H,RAD_MULT, PFRAC_TARGET)
C
C Calculate the scale factor for E,F,H to get the desired fraction, FRAC
C On failure returns EFH_SCALE < 0
C
	REAL RAD_MULT(2)
C
	COMMON /GAIN_COM/ GAIN,ESD_FACTOR,CROSSTALK,CROSS_ELLIPSE,CROSS_SMOOTH
C
	INTEGER NSUM(2)
	REAL CSUM(2),RMULT(2)
C
	REAL VAL1(7),VAL2(7),SVAL1(7),SVAL2(7)
C
C Inverse Hessian (*1512) for solving cubic polynomial for x=1 to 7
C
	REAL RMAT(4,4)
	DATA RMAT/ 12744., -12276.,  3240.,  -252.,
	1		  -12276.,  12973., -3588.,   287.,
	2		    3240.,  -3588.,  1026.,   -84.,
	3		    -252.,    287.,   -84.,     7. /
C
	REAL VEC1(4),VEC2(4),PARS1(4),PARS2(4)
C
C Set EFH_SCALE to indicate failure
C
	EFH_SCALE=-1.0
C
C Get the background intensity
C
	BKG=GET_GLOBAL_BKG(XCEN,YCEN)
C
C Iterate three times (if unable to interpolate RFRAC)
C
	RFRAC0=1.0
	DO ITER=1,3
C
C Calculate core & peaks intensities for seven RFRAC = 0.75 - 1.33
	  RFRAC=RFRAC0*0.75
	  DO I=1,7
	    RMULT(1)=RAD_MULT(1)*RFRAC
	    RMULT(2)=RAD_MULT(2)*RFRAC
	    CALL SUM_ELLIPSE_INTS(CSUM,NSUM, XCEN,YCEN,E,F,H,RMULT,2)
	    IF(CSUM(1).LE.0.0 .OR. CSUM(2).LE.0.0) GOTO 900
	    VAL1(I)=CSUM(1)-NSUM(1)*BKG
	    VAL2(I)=CSUM(2)-NSUM(2)*BKG
	    SVAL1(I)=0.02*VAL1(I) + SQRT(CSUM(1))*ESD_FACTOR
	    SVAL2(I)=0.02*VAL2(I) + SQRT(CSUM(2))*ESD_FACTOR
	    RFRAC=RFRAC*1.1
	  ENDDO
C
C Solve LSQ for core & peak intensities using a cubic polynomial
C as a function of I=1 - 7
	  DO I1=1,4
	    VEC1(I1)=0.0
	    VEC2(I1)=0.0
	    DO I2=1,7
	      VEC1(I1)=VEC1(I1)+VAL1(I2)*I2**(I1-1)
	      VEC2(I1)=VEC2(I1)+VAL2(I2)*I2**(I1-1)
	    ENDDO
	  ENDDO
C
	  DO I1=1,4
	    PARS1(I1)=0.0
	    PARS2(I1)=0.0
	    DO I2=1,4
	      PARS1(I1)=PARS1(I1)+RMAT(I1,I2)*VEC1(I2)/1512.0
	      PARS2(I1)=PARS2(I1)+RMAT(I1,I2)*VEC2(I2)/1512.0
	    ENDDO
	  ENDDO
C
C Sanity check 1 : Does Val1() increase from step 1 to 7 by > 1 esd?
	  DIFF1=(VAL1(7)-VAL1(1))/SVAL1(7)
	  DIFF2=(VAL2(7)-VAL2(1))/SVAL2(7)
	  IF(diff1 .LE. 1.0) GOTO 900
C
C Sanity check 2 : Does Val1() decrease by < 1 esd per step?
	  WORST=1E6
	  DO I=1,6
	    DIFF1=(VAL1(I+1)-VAL1(I))/SQRT( SVAL1(I+1)**2 + SVAL1(I)**2 )
	    DIFF2=(VAL2(I+1)-VAL2(I))/SQRT( SVAL2(I+1)**2 + SVAL2(I)**2 )
	    WORST=MIN(WORST, DIFF1 )
	  ENDDO
	  IF(WORST .LE. -1.0) GOTO 900
C
C Sanity check 3 : Is the fit within 2 esds of Val()?
	  WORST=0.0
	  DO I=1,7
	    DIFF1=PARS1(1)+I*(PARS1(2)+I*(PARS1(3)+I*PARS1(4))) - VAL1(I)
	    DIFF2=PARS2(1)+I*(PARS2(2)+I*(PARS2(3)+I*PARS2(4))) - VAL2(I)
	    WORST=MAX(WORST, ABS(DIFF1/SVAL1(I)), ABS(DIFF2/SVAL2(I)) )
	  ENDDO
	  IF(WORST .GT. 2.0) GOTO 900
C
C Find where "fitted" peak-frac first exceeds target value
	  DO I=1,7
	    F1=PARS1(1)+I*(PARS1(2)+I*(PARS1(3)+I*PARS1(4)))
	    F2=PARS2(1)+I*(PARS2(2)+I*(PARS2(3)+I*PARS2(4)))
	    PFRAC=F1/F2
	    IF(PFRAC .GT. PFRAC_TARGET) EXIT
	    PLAST=PFRAC
	  ENDDO
C
C Interpolate RFAC, but used capped values instead of extrapolation
C NB: Clamp RFRAC max. to 1.2 even though I=7 is not the extrema
	  IF(I .LE. 1) THEN
	    RFRAC=0.75
	  ELSEIF(I .GE. 7) THEN
	    RFRAC=1.20
	  ELSE
	    RFRAC=0.75*1.1**(I-1 - (PFRAC-PFRAC_TARGET)/(PFRAC-PLAST) )
	  ENDIF
C
C Update cumulative-product value of RFRAC
	  RFRAC0=RFRAC0*RFRAC
C
	ENDDO
C
C Convert from RFRAC0 to scale of E,F,H
C
	EFH_SCALE=1.0/RFRAC0**2
	RETURN
C
900	RETURN
	END


C ---------- Reduce number of model spots in zones ---------

	SUBROUTINE PRUNE_MODELS(NWANT_CEN,NWANT_OUT)
C
	COMMON /IZONES_COM/ IZONES(20000),DPFRAC_LIB(1000)
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /LIB_COM/ PLIB(1000,6),ILIB(1000),NLIB
C
	INTEGER NUM_ZONE(500,5),IZCUTOFFS(5),IBEFORE(5),IAFTER(5)
C
	WRITE(30,'(/,A,2(I3,A))') 'Pruning models to',NWANT_CEN,
	1					' for Zone 1, and',NWANT_OUT,' for Zones 2-5'
C
C Count how many models per zone
C
	DO IZ=1,5
	  IBEFORE(IZ)=0
	  DO I=1,NLIB
	    IF(IZONES(ILIB(I)) .EQ. IZ) IBEFORE(IZ)=IBEFORE(IZ)+1
	  ENDDO
	ENDDO
C
C Make histogram of DPFRAC_LIB (esd of the peak fraction) for models in each zone
C
	DO I1=1,500
	  DO I2=1,5
	    NUM_ZONE(I1,I2)=0
	  ENDDO
	ENDDO
C
	DO I=1,NLIB
	  IR=ILIB(I)
	  IESD=MAX(1,MIN(500, NINT(DPFRAC_LIB(I)*500.0) ))
	  NUM_ZONE(IESD,IZONES(IR))=NUM_ZONE(IESD,IZONES(IR))+1
	ENDDO
C
C Determine cutoff to get NWANT models per zone
C
	DO IZ=1,5
	  N=0
	  DO I=1,499
	    N=N+NUM_ZONE(I,IZ)
	    IF(IZ .EQ. 1) THEN
	      IF(N .GE. NWANT_CEN) EXIT
	    ELSE
	      IF(N .GE. NWANT_OUT) EXIT
	    ENDIF
	  ENDDO
	  IZCUTOFFS(IZ)=I
	ENDDO
C
C Prune model list to the new cutoffs
C
	N=0
	DO I=1,NLIB
	  IR=ILIB(I)
	  IESD=MAX(1,MIN(500, NINT(DPFRAC_LIB(I)*500.0) ))
	  IF(IESD .LE. IZCUTOFFS(IZONES(IR))) THEN
	    N=N+1
	    ILIB(N)=ILIB(I)
	    DO I2=1,6
	      PLIB(N,I2)=PLIB(I,I2)
	    ENDDO
	  ENDIF
	ENDDO
	NLIB=N
C
C Output info on strongest peaks in each zone
C
	NPRUNE=0
	DO IZ=1,5
	  IAFTER(IZ)=0
	  DO I=1,NLIB
	    IF(IZONES(ILIB(I)) .EQ. IZ) IAFTER(IZ)=IAFTER(IZ)+1
	  ENDDO
	  NPRUNE=NPRUNE+IBEFORE(IZ)-IAFTER(IZ)
	ENDDO
C
C Output summary info on cutoffs per zone
C
	WRITE(30,'(/,A)') 'Model numbers and cutoff per zone:'
	WRITE(30,'(2X,A,5I8)') 'Zone  ',(K,K=1,5)
	WRITE(30,'(2X,A,5I8)') 'Before',IBEFORE
	WRITE(30,'(2X,A,5F8.3)') 'Cutoff',IZCUTOFFS*0.001
	WRITE(30,'(2X,A,5I8)') 'After ',IAFTER
	WRITE(30,'(/,A,I5,A)') 'Total of',NLIB,' model spots remain'
C
C Complain and die if too few models
C
	IF(NLIB .LT. 6) CALL QUIT('ERROR: Less than 6 model spots were found')
C
	RETURN
	END
