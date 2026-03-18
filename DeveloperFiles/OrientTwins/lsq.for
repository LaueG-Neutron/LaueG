	SUBROUTINE RUN_LSQ(X1,X2,E1,E2, FOBS,NOBS)
C
	REAL FOBS(NOBS)
C
	COMMON /FIT_PEAKS_COM/ WAXAX(5)
	COMMON /ESDS_COM/ ESDS_WAXAX(5)
C
	COMMON /LSQ_PARS_COM/ PARS(50),ESDS(50),COVAR(50,50)
C
	REAL FCALC(2000),FSIG(2000)
C
C Load LSQ parameters into PARS, and return NPARS
C
	CALL PACK_LSQ_PARS(PARS,NPARS)
	IF(NPARS .GT. 50) CALL QUIT('BUG(run_lsq): NPARS > 50')
C
C Some info before we start
C
CC	WRITE(6,'(1X,A,I4,A,I6,A)') 'Fitting of',NOBS,
CC	1							' points using',NPARS,' parameters'
C
C Load observed data into COMMON
C
	CALL LOAD_LSQ_DATA(FOBS,NOBS)
C
C Use fixed weights
C
	DO I=1,NOBS
	  FSIG(I)=1.0
	ENDDO
C
C Call routine to setup and start NLSCON refinement routines
C
	CALL START_NLSCON(COVAR,FSIG, PARS,NPARS, FOBS,NOBS,
	1										ISTATUS,.FALSE.)
C
C Calculate FCALC() using the final parameters.
C
	CALL LSQ_FUNC(FCALC,NOBS, PARS,NPARS)
C
C Calculate the parameter esds (optionally, output bad correlations)
C
	CALL GET_PAR_ESDS(ESDS, COVAR,NPARS,.FALSE.)
C
C Unpack the refined parameters and output the refined result and esds
C
	CALL UNPACK_LSQ_PARS(PARS,ESDS,.FALSE.)
C
	X1 =WAXAX(3)
	X2 =WAXAX(5)
	E1=ESDS_WAXAX(3)
	E2=ESDS_WAXAX(5)
C
C If outright failure, give up and complain
C
	IF(ISTATUS .EQ. -1) CALL QUIT('ERROR: Least squares refinement failed')
C
	RETURN
	END


	SUBROUTINE LOAD_LSQ_DATA(FOBS0,NOBS0)
C
	REAL FOBS0(NOBS0)
C
	COMMON /LSQ_DATA_COM/ FOBS(1000),NOBS
C
	NOBS=NOBS0
	DO I=1,NOBS
		FOBS(I)=FOBS0(I)
	ENDDO
C
	RETURN
	END


	SUBROUTINE LSQ_FUNC(FCALC,NOBS, PARS,NPARS)
C
	REAL FCALC(NOBS),PARS(NPARS)
C
	COMMON /FIT_PEAKS_COM/ WAXAX(5)
C
C Unpack parameters values from PARS(1..NPARS) into /FIT_PEAKS_COM/
C PARS() in the second argument is just a dummy instead of the esds array
C
	CALL UNPACK_LSQ_PARS(PARS,PARS,.FALSE.)
C
	DO I=1,NOBS
		FCALC(I)=2.0**( WAXAX(2) - ((I-WAXAX(3))/WAXAX(1))**2 )
	  IF(NPARS .GT. 3) FCALC(I)=FCALC(I)+2.0**( WAXAX(4) - ((I-WAXAX(5))/WAXAX(1))**2 )
	ENDDO
C
	RETURN
	END



C ============ Routines to handle variable (un)packing =========

	SUBROUTINE PACK_LSQ_PARS(PARS,NPARS2)
C
	REAL PARS(*)
C
	COMMON /FIT_PEAKS_COM/ WAXAX(5)
	COMMON /ESDS_COM/ ESDS_WAXAX(5)
C
	COMMON /LSQ_REF_COM/ NPARS
C
	NPARS=4
	NPARS2=NPARS
C
C Pack parameters for 1 or 2 peaks
C
	DO I=1,NPARS
		PARS(I)=WAXAX(I+1)
	ENDDO
C
	RETURN
	END


	SUBROUTINE UNPACK_LSQ_PARS(PARS,ESDS,LPRINT)
C
	LOGICAL LPRINT
	REAL PARS(*),ESDS(*)
C
	COMMON /FIT_PEAKS_COM/ WAXAX(5)
	COMMON /ESDS_COM/ ESDS_WAXAX(5)
C
	COMMON /LSQ_REF_COM/ NPARS
C
	NPARS2=4
	IF(NPARS .NE. NPARS2) CALL QUIT('BUG(unpack_lsq_pars): NPARS incorrect')
C
C Unpack parameters for 1 or 2 peaks
C
	DO I=1,NPARS
		WAXAX(I+1)=PARS(I)
		ESDS_WAXAX(I+1)=ESDS(I)
	ENDDO
C
	IF(LPRINT) CALL PRINT_LSQ_ESDS()
C
	RETURN
	END


	SUBROUTINE PRINT_LSQ_ESDS
C
	COMMON /FIT_PEAKS_COM/ WAXAX(5)
C
	COMMON /ESDS_COM/ ESDS_WAXAX(5)
C
	COMMON /LSQ_REF_COM/ NPARS
C
	PRINT '(1X,A,F8.1,A,I3,A)','X1 =',WAXAX(3),'(',IESD_CONV(ESDS_WAXAX(3),1),')'
	PRINT '(1X,A,F8.1,A,I3,A)','X2 =',WAXAX(5),'(',IESD_CONV(ESDS_WAXAX(5),1),')'
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
	CHARACTER NAMES(5)*5
	DATA NAMES/'Width','A1   ','X1   ','A2   ','X2   '/
C
	NAME=NAMES(IPAR)
C
	RETURN
	END



C----------- Calculate the parameter esds from the covariance matrix ---------

	SUBROUTINE GET_PAR_ESDS(ESDS, COV,NPAR,LPRINT)
C
	LOGICAL LPRINT
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
C All finished if not printing anything
C
	IF( .NOT. LPRINT ) RETURN
C
C Output any cases of correlation more than +/-90%
C
	DO I1=1,NPAR
	  DO I2=I1+1,NPAR
	    IF(ABS(COV(I1,I2)) .GT. 0.90) THEN
	      CALL GET_PAR_NAME(NAME1, I1)
	      CALL GET_PAR_NAME(NAME2, I2)
	      PRINT '(1X,A,I4,A)','High correlation =',NINT(100.0*COV(I1,I2)),
	1							'% between '//TRIM(NAME1)//' & '//TRIM(NAME2)
	    ENDIF
	  ENDDO
	ENDDO
C
	RETURN
	END



C---------- ROUTINES CONTAINING NLSCON SPECIFICS -----------

	SUBROUTINE START_NLSCON(COVAR,FSIG, PARS,NPARS, FOBS,NOBS,
	1						ISTATUS,LPRINT)
C
C Suppress warning and informational output if LPRINT is false.
C Return ISTATUS=1 for success, =0 for minor problem, =-1 for failure.
C
C NB: COVAR() is the returned covariance matrix
C NB: FSIG() is the uncertainty in FOBS() used in a weighted refinement
C
C
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
      EPS = 1.0D-2
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
      IW(31)=10
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
