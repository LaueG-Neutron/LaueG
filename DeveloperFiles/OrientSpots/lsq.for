C ============= Similar to other LSQ.FOR files used to wrap NLSCON =============
C	SUBROUTINE RUN_LSQ()
C	SUBROUTINE SET_LSQ_REFINE(IPCEN,ICELL,IROT,ISIZE,IXTAL,IBEAM)
C	SUBROUTINE INIT_CELL_REFINE()
C	SUBROUTINE LSQ_FUNC(FCALC,NOBS, PARS,NPARS)
C ============ Routines to handle variable (un)packing =========
C	SUBROUTINE PACK_LSQ_PARS(PARS,NPARS)
C	SUBROUTINE UNPACK_LSQ_PARS(PARS,ESDS,LPRINT)
C	SUBROUTINE PRINT_LSQ_ESDS
C	FUNCTION IESD_CONV(ESD,IFIGS)
C	SUBROUTINE GET_PAR_NAME(NAME, IPAR)
C----------- Calculate the parameter esds from the covariance matrix ---------
C	SUBROUTINE GET_PAR_ESDS(ESDS, COV,NPAR)
C---------- ROUTINES CONTAINING NLSCON SPECIFICS -----------
C	SUBROUTINE START_NLSCON(COVAR,FSIG, PARS,NPARS, FOBS,NOBS, ISTATUS,LPRINT)
C	SUBROUTINE SET_NLSCON_OPTIONS(IOPT,IW,RW, IIW,IRW)
C	SUBROUTINE NLSCON_DUMMY_FUNC(N,M,MCON,X,DFX,IFAIL)
C	SUBROUTINE NLSCON_FUNC(NPARS, NOBS,NCON, PARS,FCALC, IFAIL)
C ==============================================================================

	SUBROUTINE RUN_LSQ()
C
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
	COMMON /HKL_MODE_COM/ IMODE5,CELL0(6)
C
	COMMON /LSQ_PARS_COM/ PARS(50),ESDS(50),COVAR(50,50)
C
	REAL FCALC(20000),FSIG(20000)
	REAL UB(3,3)
C
C Load LSQ parameters into PARS, and return NPARS
C
	CALL PACK_LSQ_PARS(PARS,NPARS)
	IF(NPARS .GT. 50) CALL QUIT('BUG(run_lsq): NPARS > 50')
C
C Some info before we start
C
	WRITE(6,'(1X,A,I4,A,I6,A)') 'Refinement started with',NPARS,
	1				' parameters using',NMATCH,' reflections'
C
C Setup number of observations for LSQ and weights to use
C
	NOBS=2*NMATCH
	DO I=1,NOBS
	  FSIG(I)=1.0
	ENDDO
C In HKL mode, add cell dimensions as mild-restraints to the refinement
	IF(IMODE5 .EQ. 5) THEN
	  XY_MATCH(1,NMATCH+1)=CELL0(1)
	  XY_MATCH(2,NMATCH+1)=CELL0(2)
	  XY_MATCH(1,NMATCH+2)=CELL0(3)
	  XY_MATCH(2,NMATCH+2)=CELL0(4)
	  XY_MATCH(1,NMATCH+3)=CELL0(5)
	  XY_MATCH(2,NMATCH+3)=CELL0(6)
	  DO I=1,3
	    FSIG(NOBS+I)=1.0
	    FSIG(NOBS+I+3)=10.0
	  ENDDO
	  NOBS=NOBS+6
	ENDIF
C
C Call routine to setup and start NLSCON refinement routines
C
	CALL START_NLSCON(COVAR,FSIG, PARS,NPARS, 
	1		XY_MATCH,NOBS, ISTATUS,.FALSE.)
C
C Calculate FCALC() using the final parameters.
C
	CALL LSQ_FUNC(FCALC,NOBS, PARS,NPARS)
C
C Calculate the parameter esds and output any bad correlations
C
	CALL GET_PAR_ESDS(ESDS, COVAR,NPARS)
C
C Unpack the refined parameters and output the refined result and esds
C
	CALL UNPACK_LSQ_PARS(PARS,ESDS,.TRUE.)
C
C If outright failure, give up and complain
C
	IF(ISTATUS .EQ. -1) CALL QUIT('ERROR: Least squares refinement failed')
C
	CALL CALC_UB_MATRIX(UB)
	RMS=MIN(99.999, CALC_RMS_ERR(UB) )
	PRINT '(I5,A,F6.3,A)',NMATCH,' spots, final rms=',RMS,' mm'
C
	RETURN
	END


	SUBROUTINE SET_LSQ_REFINE(IPCEN,ICELL,IROT,ISIZE,IXTAL,IBEAM)
C
C Set the parameter groups to refine
C
	COMMON /CELL_COM/ CELL(6),ILATT,ICEN,CELL_PROD
C
	LOGICAL LCELL,LROT,LPCEN,LSIZE,LXTAL,LBEAM
	COMMON /LSQ_REF_COM/ NPARS2,IREF_CELL(6),LCELL,LROT,LPCEN,
	1					LSIZE,LXTAL,LBEAM
C
	LCELL=(ICELL .EQ. 1)
	LROT =(IROT  .EQ. 1)
	LPCEN=(IPCEN .EQ. 1)
	LSIZE=(ISIZE .EQ. 1)
	LXTAL=(IXTAL .EQ. 1)
	LBEAM=(IBEAM .EQ. 1)
C
C No cell dimensions to refine if its cubic!
C
	IF(ILATT .EQ. 5) LCELL=.FALSE.
C
C Output which parameters are being refined
C
	IF( LCELL ) PRINT *,'Refining cell parameters'
	IF( LROT  ) PRINT *,'Refining UB rotation'
	IF( LPCEN ) PRINT *,'Refining pixel centers'
	IF( LSIZE ) PRINT *,'Refining pixel size'
	IF( LXTAL ) PRINT *,'Refining sample offset'
	IF( LBEAM ) PRINT *,'Refining vertical beam angle'
C
	RETURN
	END


	SUBROUTINE INIT_CELL_REFINE()
C
C Setup refineable and dependent parameters according to cell type
C
	COMMON /CELL_COM/ CELL(6),ILATT,ICEN,CELL_PROD
C
	LOGICAL LCELL,LROT,LPCEN,LSIZE,LXTAL,LBEAM
	COMMON /LSQ_REF_COM/ NPARS2,IREF_CELL(6),LCELL,LROT,LPCEN,
	1					LSIZE,LXTAL,LBEAM
C
C IREF contains the parameter refinement links for different cells
C ILATT=1(Tric),2(Mono),3(Ortho),4(Tetr),5(Cubic),6(Rhom),7(Trig & Hex)
C
	INTEGER IREF(6,7)
	DATA IREF/ 0,9,9,9,9,9,  0,9,9,0,9,0,  0,9,9,0,0,0,  0,1,9,0,0,0,
	1		   0,1,1,0,0,0,  0,1,1,9,4,4,  0,1,9,0,0,0/
C
C Copy the appropriate cell refinement rules for ILATT
C
	DO I=1,6
	  IREF_CELL(I)=IREF(I,ILATT)
	ENDDO
C
	RETURN
	END


	SUBROUTINE LSQ_FUNC(FCALC,NOBS, PARS,NPARS)
C
	REAL PARS(NPARS),FCALC(NOBS)
C
	COMMON /CELL_COM/ CELL(6),ILATT,ICEN,CELL_PROD
	COMMON /HKL_MODE_COM/ IMODE5,CELL0(6)
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
C
	REAL UB(3,3),PIXWAV(3,10000)
C
	IF(IMODE5 .EQ. 5) THEN
	  IF(NOBS .NE. NMATCH*2+6) CALL QUIT('BUG(lsq_func): NOBS invalid')
	ELSE
	  IF(NOBS .NE. NMATCH*2) CALL QUIT('BUG(lsq_func): NOBS invalid')
	ENDIF
C
C Unpack parameters values from PARS(1..NPARS) into COMMONs
C PARS() in the second argument is just a dummy instead of the esds array
C
	CALL UNPACK_LSQ_PARS(PARS,PARS,.FALSE.)
C
	CALL CALC_UB_MATRIX(UB)
C
	CALL CALC_HKLS_TO_PIXWAVS(PIXWAV, UB,HKL_MATCH,NMATCH)
C
C Copy pixel positions to FCALC()
C
	DO I=1,NMATCH
	  FCALC(2*I-1)=PIXWAV(1,I)
	  FCALC(2*I  )=PIXWAV(2,I)
	ENDDO
C
C For HKL mode, append current cell constants to calculated values
C
	IF(IMODE5 .EQ. 5) THEN
	  DO I=1,6
	    FCALC(2*NMATCH+I)=CELL(I)
	  ENDDO
	ENDIF
C
	RETURN
	END


C ============ Routines to handle variable (un)packing =========

	SUBROUTINE PACK_LSQ_PARS(PARS,NPARS)
C
	REAL PARS(*)
C
	COMMON /CELL_COM/ CELL(6),ILATT,ICEN,CELL_PROD
	COMMON /U_COM/ U_BASE(3,3),RXYZ(3)
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
C
	LOGICAL LCELL,LROT,LPCEN,LSIZE,LXTAL,LBEAM
	COMMON /LSQ_REF_COM/ NPARS2,IREF_CELL(6),LCELL,LROT,LPCEN,
	1					LSIZE,LXTAL,LBEAM
C
	NPARS=0
C
C Pack cell refinement parameters according to IREF_CELL rules
C
	IF( LCELL ) THEN
	  DO I=1,6
	    IF(IREF_CELL(I) .EQ. 9) THEN
	      NPARS=NPARS+1
	      PARS(NPARS)=CELL(I)
	    ENDIF
	  ENDDO
	ENDIF
C
C Pack U matrix rotation & pixel centers refinement parameters
C Do corrections to RXYZ and centers to reduce parameter correlations
C
	IF( LROT ) THEN
	  PARS(NPARS+1)=RXYZ(1)+0.5*BEAM_VERT
	  PARS(NPARS+2)=RXYZ(2)+0.5*ATAND(XTAL_OFF(1)/159.16)
	  PARS(NPARS+3)=RXYZ(3)
	  NPARS=NPARS+3
	ENDIF
C
	IF( LPCEN ) THEN
	  PARS(NPARS+1)=PIX_CEN(1)+2.0*XTAL_OFF(1)/0.2
	  PARS(NPARS+2)=PIX_CEN(2)+TAND(BEAM_VERT)*159.16/0.2
	  NPARS=NPARS+2
	ENDIF
C
C Pack pixel size and skewness refinement parameters
C
	IF( LSIZE ) THEN
	  PARS(NPARS+1)=PIX_SIZE(1)
	  PARS(NPARS+2)=PIX_SIZE(2)
	  PARS(NPARS+3)=PIX_SKEW
	  NPARS=NPARS+3
	ENDIF
C
C Pack offsets and beam angle refinement parameters
C
	IF( LXTAL ) THEN
	  PARS(NPARS+1)=XTAL_OFF(1)
	  PARS(NPARS+2)=XTAL_OFF(3)
	  NPARS=NPARS+2
	ENDIF
C
	IF( LBEAM ) THEN
	  PARS(NPARS+1)=BEAM_VERT
	  NPARS=NPARS+1
	ENDIF
C
	NPARS2=NPARS
	RETURN
	END


	SUBROUTINE UNPACK_LSQ_PARS(PARS,ESDS,LPRINT)
C
	LOGICAL LPRINT
	REAL PARS(*),ESDS(*)
C
	COMMON /CELL_COM/ CELL(6),ILATT,ICEN,CELL_PROD
	COMMON /U_COM/ U_BASE(3,3),RXYZ(3)
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
C
	COMMON /ESDS_COM/ ESDS_RXYZ(3),ESDS_CEN(2),ESDS_SIZE(2),ESDS_SKEW,
	1				ESDS_CELL(6),ESDS_XTAL(3),ESDS_BEAM
C
	LOGICAL LCELL,LROT,LPCEN,LSIZE,LXTAL,LBEAM
	COMMON /LSQ_REF_COM/ NPARS2,IREF_CELL(6),LCELL,LROT,LPCEN,
	1					LSIZE,LXTAL,LBEAM
C
	NPARS=0
C
C Unpack cell parameters according to IREF_CELL rules
C
	IF( LCELL ) THEN
	  DO I=1,6
	    IF(IREF_CELL(I) .EQ. 9) THEN
	      NPARS=NPARS+1
	      CELL(I)=PARS(NPARS)
	      ESDS_CELL(I)=ESDS(NPARS)
	    ELSE IF(IREF_CELL(I) .NE. 0) THEN
	      CELL(I)=CELL(IREF_CELL(I))
	    ENDIF
	  ENDDO
	ENDIF
C
C Unpack RXYZ angles and pixel centers
C
	IF( LROT ) THEN
	  RXYZ(1)=PARS(NPARS+1)
	  RXYZ(2)=PARS(NPARS+2)
	  RXYZ(3)=PARS(NPARS+3)
	  ESDS_RXYZ(1)=ESDS(NPARS+1)
	  ESDS_RXYZ(2)=ESDS(NPARS+2)
	  ESDS_RXYZ(3)=ESDS(NPARS+3)
	  NPARS=NPARS+3
	ENDIF
C
	IF( LPCEN ) THEN
	  PIX_CEN(1)=PARS(NPARS+1)
	  PIX_CEN(2)=PARS(NPARS+2)
	  ESDS_CEN(1)=ESDS(NPARS+1)
	  ESDS_CEN(2)=ESDS(NPARS+2)
	  NPARS=NPARS+2
	ENDIF
C
C Unpack pixel size and shape
C
	IF( LSIZE ) THEN
	  PIX_SIZE(1)=PARS(NPARS+1)
	  PIX_SIZE(2)=PARS(NPARS+2)
	  PIX_SKEW=PARS(NPARS+3)
	  ESDS_SIZE(1)=ESDS(NPARS+1)
	  ESDS_SIZE(2)=ESDS(NPARS+2)
	  ESDS_SKEW=ESDS(NPARS+3)
	  NPARS=NPARS+3
	ENDIF
C
C Unpack offsets and beam angle
C
	IF( LXTAL ) THEN
	  XTAL_OFF(1)=PARS(NPARS+1)
	  XTAL_OFF(3)=PARS(NPARS+2)
	  ESDS_XTAL(1)=ESDS(NPARS+1)
	  ESDS_XTAL(3)=ESDS(NPARS+2)
	  NPARS=NPARS+2
	ENDIF
C
	IF( LBEAM ) THEN
	  BEAM_VERT=PARS(NPARS+1)
	  ESDS_BEAM=ESDS(NPARS+1)
	  NPARS=NPARS+1
	ENDIF
C
C Reverse the corrections made in PACK_LSQ_PARS
C
	IF( LROT ) THEN
	  RXYZ(1)=RXYZ(1)-0.5*BEAM_VERT
	  RXYZ(2)=RXYZ(2)-0.5*ATAND(XTAL_OFF(1)/159.16)
	ENDIF
C
	IF( LPCEN ) THEN
	  PIX_CEN(1)=PIX_CEN(1)-2.0*XTAL_OFF(1)/0.2
	  PIX_CEN(2)=PIX_CEN(2)-TAND(BEAM_VERT)*159.16/0.2
	ENDIF
C
C
	IF(NPARS .NE. NPARS2) CALL QUIT('BUG(unpack_lsq_pars): NPARS incorrect')
C
	IF(LPRINT) CALL PRINT_LSQ_ESDS()
C
	RETURN
	END


	SUBROUTINE PRINT_LSQ_ESDS
C
	COMMON /CELL_COM/ CELL(6),ILATT,ICEN,CELL_PROD
	COMMON /U_COM/ U_BASE(3,3),RXYZ(3)
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
C
	COMMON /ESDS_COM/ ESDS_RXYZ(3),ESDS_CEN(2),ESDS_SIZE(2),ESDS_SKEW,
	1				ESDS_CELL(6),ESDS_XTAL(3),ESDS_BEAM
C
	LOGICAL LCELL,LROT,LPCEN,LSIZE,LXTAL,LBEAM
	COMMON /LSQ_REF_COM/ NPARS2,IREF_CELL(6),LCELL,LROT,LPCEN,
	1					LSIZE,LXTAL,LBEAM
C
	IF( LCELL ) THEN
	  PRINT '(1X,A,F8.4)','A =',CELL(1)
	    IF(IREF_CELL(2) .EQ. 9) PRINT '(1X,A,F8.4,A,I3,A)',
	1				'B =',CELL(2),'(',IESD_CONV(ESDS_CELL(2),4),')'
	    IF(IREF_CELL(3) .EQ. 9) PRINT '(1X,A,F8.4,A,I3,A)',
	1				'C =',CELL(3),'(',IESD_CONV(ESDS_CELL(3),4),')'
	    IF(IREF_CELL(4) .EQ. 9) PRINT '(1X,A,F8.3,A,I3,A)',
	1				'alpha=',CELL(4),'(',IESD_CONV(ESDS_CELL(4),3),')'
	    IF(IREF_CELL(5) .EQ. 9) PRINT '(1X,A,F8.3,A,I3,A)',
	1				'beta =',CELL(5),'(',IESD_CONV(ESDS_CELL(5),3),')'
	    IF(IREF_CELL(6) .EQ. 9) PRINT '(1X,A,F8.3,A,I3,A)',
	1				'gamma=',CELL(6),'(',IESD_CONV(ESDS_CELL(6),3),')'
	ENDIF
C
	IF( LROT ) WRITE(6,'(1X,A,3(F8.3,A,I3,A),A)',IOSTAT=IDUM)
	1	'ROT(xyz)=',(RXYZ(K),
     2	'(',IESD_CONV(ESDS_RXYZ(K),3),')',K=1,3),' degrees'
C
	IF( LPCEN ) WRITE(6,'(1X,A,2(F8.2,A,I3,A),A)',IOSTAT=IDUM)
	1	'Center(xy)=',(PIX_CEN(K),
     2	'(',IESD_CONV(ESDS_CEN(K),2),')',K=1,2),' pixels'
C
	IF( LSIZE ) WRITE(6,'(1X,A,2(F8.2,A,I3,A),3X,A,F8.5,A,I3,A)',
	1	IOSTAT=IDUM) 'PIXEL(xy)=',(1E3*PIX_SIZE(K),
     2	'(',IESD_CONV(ESDS_SIZE(K),5),')',K=1,2),
     3	'Skew=',PIX_SKEW,'(',IESD_CONV(ESDS_SKEW,5),')'
C
	IF( LXTAL ) THEN
	  WRITE(6,'(1X,A,2(F8.3,A,I3,A))',IOSTAT=IDUM )
	1		'Xtal offset(xz)=',(XTAL_OFF(K),
	2		'(',IESD_CONV(ESDS_XTAL(K),3),')',K=1,3,2)
	ENDIF
C
	IF( LBEAM ) THEN
	  WRITE(6,'(1X,A,F8.3,A,I3,A)',IOSTAT=IDUM)
	1		'Beam(v)=',BEAM_VERT,'(',IESD_CONV(ESDS_BEAM,3),')'
	ENDIF
C
	RETURN
	END


	FUNCTION IESD_CONV(ESD,IFIGS)
C
	ISCALE=10**IFIGS
	IESD_CONV=NINT(MIN(999.0, ISCALE*ESD ))
	RETURN
	END


	SUBROUTINE GET_PAR_NAME(NAME, IPAR)
C
	CHARACTER NAME*(*)
C
	LOGICAL LCELL,LROT,LPCEN,LSIZE,LXTAL,LBEAM
	COMMON /LSQ_REF_COM/ NPARS2,IREF_CELL(6),LCELL,LROT,LPCEN,
	1					LSIZE,LXTAL,LBEAM
C
	LOGICAL LPARS(5)
	CHARACTER NAMES(5)*5,NAMES2(5)*6
	INTEGER NPARS(5)
	DATA NAMES/'B    ','C    ','alpha','beta ','gamma'/
	DATA NAMES2/'Rot   ','Center','Pixel ','Xtal  ','Beam  '/
	DATA NPARS /    3   ,   2    ,   3    ,   2    ,   1    /
C
	N=IPAR
	NAME=REPEAT(' ',LEN(NAME))
C
	IF( LCELL ) THEN
	  DO I=1,5
	    IF(IREF_CELL(I+1) .EQ. 9) THEN
	      IF(N .EQ. 1) NAME=NAMES(I)
	      N=N-1
	    ENDIF
	  ENDDO
	ENDIF
C
	LPARS(1)=LROT
	LPARS(2)=LPCEN
	LPARS(3)=LSIZE
	LPARS(4)=LXTAL
	LPARS(5)=LBEAM
C
	DO I=1,5
	  IF( LPARS(I) ) THEN
	    IF(N.GE.1 .AND. N.LE.NPARS(I))
	1	WRITE(NAME,'(2A,I1,A)') TRIM(NAMES2(I)),'(',N,')'
	    N=N-NPARS(I)
	  ENDIF
	ENDDO
C
	IF(N .GT. 0) CALL QUIT('BUG(get_par_name): Invalid IPAR')
C
	RETURN
	END


C----------- Calculate the parameter esds from the covariance matrix ---------

	SUBROUTINE GET_PAR_ESDS(ESDS, COV,NPAR)
C
	REAL ESDS(NPAR),COV(NPAR,NPAR)
C
	CHARACTER*20 NAME1,NAME2
C
C Extract parameter esds from the covariance matrix
C
	DO I=1,NPAR
	  ESDS(I)=SQRT(COV(I,I))
	ENDDO
C
C Convert the covariance matrix to a correlation matrix
C
	DO I1=1,NPAR
	  DO I2=1,NPAR
	    IF(ESDS(I1)*ESDS(I2) .EQ. 0.0) THEN
	      COV(I1,I2)=0.0
	    ELSE
	      COV(I1,I2)=COV(I1,I2)/(ESDS(I1)*ESDS(I2))
	    ENDIF
	  ENDDO
	ENDDO
C
C Output any cases of correlation more than +/-90%
C
	DO I1=1,NPAR
	  DO I2=I1+1,NPAR
	    IF(ABS(COV(I1,I2)) .GT. 0.90) THEN
	      CALL GET_PAR_NAME(NAME1, I1)
	      CALL GET_PAR_NAME(NAME2, I2)
	      PRINT '(1X,A,I4,A)','High correlation =',NINT(100.0*COV(I1,I2)),
	1		  '% between parameters '//TRIM(NAME1)//' & '//TRIM(NAME2)
	    ENDIF
	  ENDDO
	ENDDO
C
	RETURN
	END


C---------- ROUTINES CONTAINING NLSCON SPECIFICS -----------

	SUBROUTINE START_NLSCON(COVAR,FSIG, PARS,NPARS, FOBS,NOBS, ISTATUS,LPRINT)
C
C Suppress warning and informational output if LPRINT is false.
C Return ISTATUS=1 for success, =0 for minor problem, =-1 for failure.
C
C NB: COVAR() is the returned covariance matrix
C NB: FSIG() is the uncertainty in FOBS() used in a weighted refinement
C
	LOGICAL LPRINT
	REAL COVAR(NPARS*NPARS),FSIG(NOBS),PARS(NPARS),FOBS(NOBS)
C
C Size of work arrays, hopefully big enough!
C
	PARAMETER (IRW=100000)
	PARAMETER (IIW=100)
	REAL RW(IRW)
	INTEGER IW(IIW)
C
	INTEGER IOPT(50)
C
C Two more arrays hopefully big enough, will check later.
C
	PARAMETER (NPAR_MAX=30)
	PARAMETER (NOBS_MAX=10000)
	REAL PAR_SCALE(NPAR_MAX)
C
	CHARACTER*20 NAME1,NAME2
C
	EXTERNAL NLSCON_FUNC,NLSCON_DUMMY_FUNC
C
C Complain and die if we haven't compiled in double-precision
C
	IF(1.0 .EQ. 1.0+1E-10) CALL QUIT('BUG: RECOMPILE USING DOUBLE PRECISION')
C
C Sanity check on array sizes.
C
	IF(NPARS .GT. NPAR_MAX) CALL QUIT('BUG(start_nlscon): Too many parameters')
	IF(NOBS .GT. NOBS_MAX) CALL QUIT('BUG(start_nlscon): Too many reflections')
	IF(IIW .LT. NPARS+52) CALL QUIT('BUG(start_nlscon): IIW too small')
	IF(IRW .LT. (2*NOBS+NPARS)*NPARS + 8*NOBS + 10*NPARS +
	1	MAX(NOBS,NPARS) + 61) CALL QUIT('BUG(start_nlscon): IRW too small')
C
C Setup usual NLSCON options.
C
	CALL SET_NLSCON_OPTIONS(IOPT,IW,RW,IIW,IRW)
C
C Set lower threshold for parameter scaling to << 1
C
	DO I=1,NPARS
		PAR_SCALE(I)=MAX(1.0,ABS(PARS(I)))*1E-8
	ENDDO
C
C Set required parameter precision.
C
      EPS = 1.0D-5
C
	IF(LPRINT) PRINT *
      ITER=0
100	ITER=ITER+1
C Have edited NLSCON so RW(50...) returns esd's of parameters.
	IF(LPRINT) PRINT '(1H+,A,I3)','Iteration',ITER
	CALL NLSCON(NPARS, NOBS,NOBS, NLSCON_FUNC,NLSCON_DUMMY_FUNC,
	1		PARS,PAR_SCALE, FOBS,FSIG, EPS, IOPT, IERR, IIW,IW,IRW,RW)
	IF (IERR.EQ.-1) GOTO 100
	IF(LPRINT) WRITE(6,'(1X,A,I4,A)') 'Total of',IW(5)+IW(8),
	1				' function calls by NLSCON'
	IF(RW(32) .GT. 0.5) THEN
	  PRINT *,'WARNING: Refinement may be highly non-linear'
	  WRITE(6,'(1X,A,F7.2,A)') 'WARNING(start_nlscon): SKAP =',RW(32),
	1			' > 0.5 indicates possible linearity problem'
	ENDIF
C
C Copy covariance matrix into COVAR().
C NB: The NLSCON documentation says RW(...) contains the
C     correlation matrix, this is incorrect.
C
	DO I=1,NPARS*NPARS
		COVAR(I)=RW(50+2*NPARS+I)
	ENDDO
C
C If IERR=0, then success and print a message to the console.
C If IERR=1,2,3, then minor failure, so warn and continue.
C If IERR>3, then major failure, warn and stop program.
C
	IF(IERR .EQ. 0) THEN		! Success, output some info
	  ISTATUS=1
	  IF(LPRINT) WRITE(6,'(1X,A,I3,A)') 'Least-squares successful in',ITER,' iterations'
	ELSE IF(IERR .EQ. 1) THEN	! Stationary point, so output some info
	  ISTATUS=0			! treat this as a warning, not an error
	  IF( .NOT.LPRINT ) RETURN
C
	  PRINT '(/,1X,A)','Refinement detected a stationary point'
C
	  PRINT '(2X,A)','The following parameters have a zero Jacobian:'
	  DO I=1,NPARS
	    VAR=COVAR( I + (I-1)*NPARS )
	    CALL GET_PAR_NAME(NAME1, I)
	    IF(VAR .EQ. 0.0) THEN
		  PRINT '(4X,A)',TRIM(NAME1)
	    ENDIF
	  ENDDO
C
	  PRINT '(2X,A)','The following parameters have extreme uncertainties:'
	  DO I=1,NPARS
	    VAR=SQRT( COVAR( I + (I-1)*NPARS ) )
	    CALL GET_PAR_NAME(NAME1, I)
	    IF(VAR .GE. 1000.0) THEN
		  PRINT '(4X,A,1P,E16.2)',TRIM(NAME1),VAR
	    ENDIF
	  ENDDO
C
	  PRINT '(2X,A)','The following parameters have extreme correlations:'
	  DO I1=1,NPARS
	    DO I2=I1+1,NPARS
	      VAR=COVAR(I1+(I2-1)*NPARS)/SQRT(COVAR(I1+(I1-1)*NPARS)*COVAR(I2+(I2-1)*NPARS))
	      IF(ABS(VAR) .GT. 0.999) THEN
	        CALL GET_PAR_NAME(NAME1, I1)
	        CALL GET_PAR_NAME(NAME2, I2)
		  PRINT '(4X,3A,F10.4,A)',TRIM(NAME1),', ',TRIM(NAME2),VAR*100.0,'%'
	      ENDIF
	    ENDDO
	  ENDDO
C
	ELSE IF(IERR .EQ. 2) THEN		! Reached maximum number of iterations
	  ISTATUS=0
	  IF(LPRINT) PRINT '(/,1X,A,/)','WARNING: Refinement not fully converged'
	ELSE IF(IERR .EQ. 3) THEN		! Damping factor too small
	  ISTATUS=0
	  IF(LPRINT) PRINT '(/,1X,A,/)','WARNING: Refinement has poor convergence'
	ELSE		! Trap all remaining errors and incorrect parameter values
C Give explicit complaints if NLSCON is given the wrong values, then stop
	  IF(IERR .EQ. 10) CALL QUIT('BUG(start_nlscon): IWK or RWK too small')
	  IF(IERR .EQ. 20) CALL QUIT('BUG(start_nlscon): Invalid N,M,MFIT')
	  IF(IERR .EQ. 21) CALL QUIT('BUG(start_nlscon): RTOL < 0')
	  IF(IERR .EQ. 22) CALL QUIT('BUG(start_nlscon): XSCAL() < 0')
	  IF(IERR .EQ. 30) CALL QUIT('BUG(start_nlscon): Invalid IOPT()')
C Not sure what is wrong, so print out some information and stop
	  PRINT '(/,1X,A,I4,A)','Total of',IW(5)+IW(8),' function calls by NLSCON'
	  PRINT *,'Error',IERR,' returned from NLSCON'
	  WRITE(6,'(/,1X,A,I4,A)') 'Total of',IW(5)+IW(8),' function calls by NLSCON'
	  WRITE(6,'(1X,A,I4,A)') 'Error',IERR,' returned from NLSCON'
	  CALL QUIT('ERROR: Least-squares failed')
	ENDIF
C
	RETURN
      END


	SUBROUTINE SET_NLSCON_OPTIONS(IOPT,IW,RW, IIW,IRW)
C
	INTEGER IOPT(50),IW(IIW)
	REAL RW(IRW)
C
C Begin by zeroing arrays to give default values.
C
	DO I=1,50
		IOPT(I)=0
	ENDDO
	DO I=1,IIW
		IW(I)=0
	ENDDO
	DO I=1,IRW
		RW(I)=0.0D0
	ENDDO
C
C Now override some defaults using non-zero values.
C
C (=1) Execution mode: Stepwise mode
      IOPT(2)=1
C (=3) Jacobian: computed by numerical differentation (with feedback)
      IOPT(3)=3
C (=1) A posteriori statistical analysis: yes
      IOPT(21)=1
C (=0) Broyden updates: inhibit
      IOPT(32)=0
C (=2) Problem classification: mildly nonlinear
C (=3) Problem classification: highly nonlinear
      IOPT(31)=3
C (=0) Automatic row scaling: allowed
      IOPT(35)=0
C (=2) Output error and warning messages to unit 6
C	IOPT(11)=6
C	IOPT(12)=6
C
C
C Override maximum allowed number of iterations:
C     IW(31)       NITMAX      Maximum number of permitted iteration (default: 50)
      IW(31)=20
C     Override initial pseudo-rank:
C     IW(32)=N
C
C
C     Override starting damping factor:
      RW(21)=1.0D+1
C     Override minimal allowed damping factor:
      RW(22)=1.0D-6
C     Override rank1-decision parameter SIGMA:
C     RW(23)=2.0D0
C     Override maximum permitted subcondition for DECCON:
      RW(25)= 1.0D+16
C
	RETURN
	END


      SUBROUTINE NLSCON_DUMMY_FUNC(N,M,MCON,X,DFX,IFAIL)
C Dummy function for Jacobian evaluation.
	PRINT *,N,M,MCON,X,DFX,IFAIL
	CALL QUIT('BUG(nlscon_dummy_func): Function should never be called!')
      END


	SUBROUTINE NLSCON_FUNC(NPARS, NOBS,NCON, PARS,FCALC, IFAIL)
C
	REAL PARS(*),FCALC(*)
C
C A sanity checks (There should be no constraints)
C
	IF(NCON .NE. 0) CALL QUIT('BUG(f): NCON .NE. 0')
C
	CALL LSQ_FUNC(FCALC,NOBS, PARS,NPARS)
C
	IFAIL=0
	RETURN
	END