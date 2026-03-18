C----- Perform the LSQ refinement -------------------------------------
c	SUBROUTINE RUN_LSQ
C----- The function called by LSQ refinement --------------------------
c	SUBROUTINE LSQ_FUNC(ICALC0,NDATA, PARS,NPARS)
C----- Load data used for LSQ refinement ------------------------------
c	SUBROUTINE LOAD_LSQ_DATA
C----- Routines for information on refinement results -----------------
c	SUBROUTINE PRINT_LSQ_MERGE
c	SUBROUTINE LOG_LSQ_MERGE
c	SUBROUTINE WRITE_LSQ_MERGE(IUNIT,ISPACE)
c	SUBROUTINE CALC_LSQ_MERGE_SEQ(R1,R1_ESD,R2,R2_ESD,GOOF,NSUM,
c	1												ISEQ1,ISEQ2,RSIG)
C----- Calculate the parameter esds from the covariance matrix --------
c	SUBROUTINE GET_PAR_ESDS(ESDS, COV,NPARS)
C----- NLSCON SPECIFIC ROUTINES ---------------------------------------
c	SUBROUTINE START_NLSCON(COVAR,ISIG, PAR,NPARS, IOBS,NOBS)
c	SUBROUTINE SET_NLSCON_OPTIONS(IOPT,IW,RW, IIW,IRW)
c	SUBROUTINE NLSCON_DUMMY_FUNC(N,M,MCON,X,DFX,IFAIL)
c	SUBROUTINE NLSCON_FUNC(NPARS, NOBS,NCON, PAR,ICALC, IFAIL)
C----------------------------------------------------------------------


C----- Perform the LSQ refinement -------------------------------------
	SUBROUTINE RUN_LSQ
C
C Start the LSQ refinements of NDATA reflections. The number of refined
C parameters is returned in NPARS.
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
      REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	PARAMETER (NPARS_MAX=500)
	COMMON /LSQ_PARS_COM/ PARS(NPARS_MAX),ESDS(NPARS_MAX),COVAR(NPARS_MAX,NPARS_MAX),NPARS
C
C Sanity check on data size
C
	IF(NDATA .GT. NDATA_MAX) STOP 'BUG(run_lsq): NDATA > NDATA_MAX'
C
C Load LSQ parameters into PARS, and check on the total number
C Also, save the position of the various refined parameters
C
	CALL PACK_LSQ_PARS(PARS,NPARS)
	IF(NPARS .GT. NPARS_MAX) STOP 'BUG(run_lsq): NPARS > NPARS_MAX'
C
C Starting from a value of 100, calculate COUNTS_SEQ() given the
C current normalisation parameters
C
	DO ISEQ=1,NSEQ
	  COUNTS_SEQ(ISEQ)=100.0
	ENDDO
	CALL CALC_ALL_COUNTS_SEQ
C
C Some info before we start
C
	WRITE(10,'(A,I4,A,I6,A)') 'Refinement started with',NPARS,
	1						' parameters using',NDATA,' reflections'
C
C Call routine to setup and start NLSCON refinement routines
C
	CALL START_NLSCON(COVAR,ISIG, PARS,NPARS, IOBS,NDATA)
C
C Calculate ICALC() with the final parameters.
C
	CALL LSQ_FUNC(ICALC,NDATA, PARS,NPARS)
C
C Calculate the parameter esds and output any bad correlations
C
	CALL GET_PAR_ESDS(ESDS, COVAR,NPARS)
C
C Unpack the refined parameters and output the refined result and esds
C
	CALL UNPACK_LSQ_PARS(PARS,ESDS,.TRUE.)
C
	PRINT *
	RETURN
	END


C----- The function called by LSQ refinement --------------------------

	SUBROUTINE LSQ_FUNC(ICALC0,NDATA, PARS,NPARS)
C
	REAL ICALC0(NDATA),PARS(NPARS)
C
	PARAMETER (NDATA_MAX=2000000)
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
C Unpack PARS() so CALC_ALL_COUNTS_SEQ uses the new parameter values
C PARS() in the second argument is just a dummy instead of an esds array
C
	CALL UNPACK_LSQ_PARS(PARS,PARS,.FALSE.)
C
C Calculate COUNTS_SEQ() for all sequences, which in turn calculates ICALC()
C for all spots using the current normalisation parameters
C
	CALL CALC_ALL_COUNTS_SEQ
C
C Copy the ICALC() made by CALC_SEQ_COUNTS to the argument ICALC0().
C ICALC() & ICALC0() should be the same array, but copy it just in case!
C
	DO I=1,NDATA
	  ICALC0(I)=ICALC(I)
	ENDDO
C
	RETURN
	END


C----- Load data used for LSQ refinement ------------------------------

	SUBROUTINE LOAD_LSQ_DATA
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	WRITE(10,'(A,I6,A)') 'Loading data and refinement weights for',
	1											NDATA,' reflections'
C
C Load data into IOBS() & ISIG()
C
	NDATA2=0
	DO ISEQ=1,NSEQ
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    NDATA2=NDATA2+1
	    IF(NDATA2 .GT. NDATA_MAX) STOP 'BUG(load_lsq_data): NDATA2 > NDATA_MAX'
	    IF(NDATA2 .NE. I) STOP 'BUG(load_lsq_data): NDATA2 .NE. I'
	    IF(DCOUNTS(I) .LT. 0.0) STOP 'BUG(load_lsq_data): DCOUNTS() < 0'
	    IOBS(NDATA2)=COUNTS(I)
	    ISIG(NDATA2)=DCOUNTS(I)
	  ENDDO
	ENDDO
C
	IF(NDATA2 .NE. NDATA) STOP 'Bug(LOAD_LSQ_DATA): NDATA2 .NE. NDATA'
C
C Load starting values of 100 for COUNTS_SEQ()
C NB: Call CALC_ALL_COUNTS_SEQ later on to refine better values
C
	DO ISEQ=1,NSEQ
	  COUNTS_SEQ(ISEQ)=100.0
	ENDDO
C
	RETURN
	END


C----- Routines for information on refinement results -----------------

	SUBROUTINE PRINT_LSQ_MERGE
C
C Write R-factors, GOOF, etc. to the console
C
	CALL WRITE_LSQ_MERGE(6,3)
	RETURN
	END



	SUBROUTINE LOG_LSQ_MERGE
C
C Write R-factors, GOOF, etc. to the log file
C
	CALL WRITE_LSQ_MERGE(10,2)
	RETURN
	END



	SUBROUTINE WRITE_LSQ_MERGE(IUNIT,ISPACE)
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	CALL CALC_LSQ_MERGE_SEQ(R1,R1_ESD,R2,R2_ESD,GOOF,NSUM, 1,NSEQ,4.0)
	ISIZ=CEILING(ALOG10(REAL(NSUM)))
	WRITE(IUNIT,'(<ISPACE>X,2(A,F4.1,A,F4.1,A,2X),A,F4.1,2X,A,I<ISIZ>,A)',IOSTAT=IDUM)
	1	'R =',R1,' (',R1_ESD,')%',
	2	'wR2 =',R2,' (',R2_ESD,')%',
	3	'GOOF =',GOOF,'(',NSUM,' reflns > 4 sig)'
C
	CALL CALC_LSQ_MERGE_SEQ(R1,R1_ESD,R2,R2_ESD,GOOF,NSUM, 1,NSEQ,-1.0)
	ISIZ=CEILING(ALOG10(REAL(NSUM)))
	WRITE(IUNIT,'(<ISPACE>X,2(A,F4.1,A,F4.1,A,2X),A,F4.1,2X,A,I<ISIZ>,A)',IOSTAT=IDUM)
	1	'R =',R1,' (',R1_ESD,')%',
	2	'wR2 =',R2,' (',R2_ESD,')%',
	3	'GOOF =',GOOF,'(all ',NSUM,' reflns)'
C
	RETURN
	END



	SUBROUTINE CALC_LSQ_MERGE_SEQ(R1,R1_ESD,R2,R2_ESD,GOOF,NSUM,
	1												ISEQ1,ISEQ2,RSIG)
C
C Do the sums to calculate R_MERGE for intensities > RSIG * esd
C Ignore any reflections marked as outliers (ISIG < 0)
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
C Perform the usual sums to determine the statistics.
C
	DOF=0.0
	NSUM=0
	SUM_I1=0.0
	SUM_D1=0.0
	SUM_S1=0.0
	SUM_I2=0.0
	SUM_D2=0.0
	DO ISEQ=ISEQ1,ISEQ2
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
C Is the refln acceptable in terms of Int/esd?
	    IF(RSIG.LT.0.0 .OR. IOBS(I).GT.RSIG*ISIG(I)) THEN
C Ignore any outlier reflns and sequences of just 1 refln
	      IF(ISIG(I).GT.0.0 .AND. IFIRST(ISEQ).NE.ILAST(ISEQ)) THEN
C Calculating the degrees of freedom, DOF, is a bit tricky as we are not
C summing over complete sequences. For each refln, we add 1 to the d.o.f.
C and subtract 1/(number in sequence).
	        DOF=DOF +1.0 -1.0/(ILAST(ISEQ)-IFIRST(ISEQ)+1)
	        NSUM=NSUM+1
	        SUM_I1=SUM_I1+IOBS(I)
	        SUM_D1=SUM_D1+ABS(IOBS(I)-ICALC(I))
	        SUM_S1=SUM_S1+ABS(ISIG(I))
	        SUM_I2=SUM_I2+IOBS(I)**2/ISIG(I)**2
	        SUM_D2=SUM_D2+(IOBS(I)-ICALC(I))**2/ISIG(I)**2
C
	      ENDIF
	    ENDIF
C
	  ENDDO
	ENDDO
C
C Ensure no divide by zero errors
C
	NSUM=MAX(1,NSUM)
	DOF=MAX(1.0,DOF)
	SUM_I1=MAX(1E-6,SUM_I1)
	SUM_I2=MAX(1E-6,SUM_I2)
C
C R1 is the unweighted R = SUM{|Ihkl -Imerge|} / SUM{Ihkl}
C R2 is the weighted-rms R = SQRT(  SUM{w * |Ihkl -Imerge|^2} / SUM{w * Ihkl^2}  )
C GOOF is the goodness-of-fit = SQRT(  SUM{w * |Ihkl -Imerge|^2} / DOF  )
C                   where w is the weight = 1.0 / variance(Ihkl)
C
	R1=        SUM_D1/SUM_I1
	R2=  SQRT( SUM_D2/SUM_I2 )
	GOOF=SQRT( SUM_D2/DOF )
C The factor (pi/2)^1/2 is the ratio of the "mean absolute deviation" to the
C standard deviation for a Normal distribution.
	R1_ESD=SUM_S1/SUM_I1  /SQRT(3.141593/2.0) *SQRT(DOF/NSUM)
	R2_ESD=R2/MAX(0.1,GOOF)
C
C Convert R1 & R2 to percentages, and limits the range of all values
C
	R1    =MIN(99.9,100.0*R1)
	R1_ESD=MIN(99.9,100.0*R1_ESD)
	R2    =MIN(99.9,100.0*R2)
	R2_ESD=MIN(99.9,100.0*R2_ESD)
	GOOF  =MIN(9.9, GOOF )
C
	RETURN
	END


C----- Calculate the parameter esds from the covariance matrix --------

	SUBROUTINE GET_PAR_ESDS(ESDS, COV,NPARS)
C
	REAL ESDS(NPARS),COV(NPARS,NPARS)
C
	CHARACTER NAME1*8,NAME2*8
C
C Extract parameter esds from the covariance matrix
C
	DO I=1,NPARS
	  ESDS(I)=SQRT(COV(I,I))
	ENDDO
C
C Convert the covariance matrix to a correlation matrix
C
	DO I1=1,NPARS
	  DO I2=1,NPARS
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
	DO I1=1,NPARS
	  DO I2=I1+1,NPARS
	    IF(ABS(COV(I1,I2)) .GT. 0.90) THEN
	      CALL GET_PAR_NAME(NAME1, I1)
	      CALL GET_PAR_NAME(NAME2, I2)
	      PRINT '(1X,A,I4,A)','High correlation =',NINT(100.0*COV(I1,I2)),
	1		  '% between parameters '//TRIM(NAME1)//' & '//TRIM(NAME2)
	      WRITE(10,'(A,I4,A)') 'High correlation =',NINT(100.0*COV(I1,I2)),
	1		  '% between parameters '//TRIM(NAME1)//' & '//TRIM(NAME2)
	    ENDIF
	  ENDDO
	ENDDO
C
	RETURN
	END



C----- NLSCON SPECIFIC ROUTINES ---------------------------------------

	SUBROUTINE START_NLSCON(COVAR,ISIG, PAR,NPARS, IOBS,NOBS)
C
C NB: COVAR() is the returned covariance matrix
C NB: ISIG() is the uncertainty in IOBS() used in a weighted refinement
C
	REAL COVAR(NPARS*NPARS),ISIG(NOBS),PAR(NPARS),IOBS(NOBS)
C
C Size of work arrays, should be big enough!
C
	PARAMETER (IRW=20000000)
	PARAMETER (IIW=600)
	REAL RW(IRW)
	INTEGER IW(IIW)
C
	INTEGER IOPT(50)
C
C Two more arrays hopefully big enough, will check later.
C
	PARAMETER (NPARS_MAX=500)
	REAL PAR_SCALE(NPARS_MAX)
C
	CHARACTER NAME1*8,NAME2*8
C
	EXTERNAL NLSCON_FUNC,NLSCON_DUMMY_FUNC
C
C Complain and die if we haven't compiled in double-precision
C
	IF(1.0 .EQ. 1.0+1E-10) STOP 'BUG: RECOMPILE USING DOUBLE PRECISION'
C
C Setup usual NLSCON options.
C
	CALL SET_NLSCON_OPTIONS(IOPT,IW,RW,IIW,IRW)
C
C Set lower threshold for parameter scaling to << 1
C
	DO I=1,NPARS
		PAR_SCALE(I)=MAX(1.0,ABS(PAR(I)))*1E-8
	ENDDO
C
C Set required parameter precision.
C
      EPS = 1.0D-5
C
	PRINT *
      ITER=0
100	ITER=ITER+1
C Have editted NLSCON so RW(50...) returns esd's of parameters.
	PRINT '(1H+,A,I3)','Iteration',ITER
	CALL NLSCON(NPARS, NOBS,NOBS, NLSCON_FUNC,NLSCON_DUMMY_FUNC,
	1		PAR,PAR_SCALE, IOBS,ISIG, EPS, IOPT, IERR, IIW,IW,IRW,RW)
	IF (IERR .EQ. -1) GOTO 100
C Specific error for work arrays too small
	IF (IERR .EQ. 10) GOTO 999
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
	IF(IERR .EQ. 0) THEN			! Success, so give minimal info
	  WRITE(10,'(A,I3,A)') 'Least-squares successful in',ITER,' iterations'
	ELSE IF(IERR .EQ. 1) THEN		! Stationary point, so output some info
	  PRINT '(/,1X,A)','Refinement detected a stationary point'
	  WRITE(10,'(/,A)') 'Refinement detected a stationary point'
C
	  PRINT '(2X,A)','The following parameters have a zero Jacobian:'
	  WRITE(10,'(1X,A)') 'The following parameters have a zero Jacobian:'
	  DO I=1,NPARS
	    VAR=COVAR( I + (I-1)*NPARS )
	    CALL GET_PAR_NAME(NAME1, I)
	    IF(VAR .EQ. 0.0) THEN
		  PRINT '(4X,A)',NAME1
		  WRITE(10,'(3X,A)') NAME1
	    ENDIF
	  ENDDO
C
	  PRINT '(2X,A)','The following parameters have extreme uncertainties:'
	  WRITE(10,'(1X,A)') 'The following parameters have extreme uncertainties:'
	  DO I=1,NPARS
	    VAR=SQRT( COVAR( I + (I-1)*NPARS ) )
	    CALL GET_PAR_NAME(NAME1, I)
	    IF(VAR .GE. 1000.0) THEN
		  PRINT '(4X,A,1P,E16.2)',NAME1,VAR
		  WRITE(10,'(3X,A,1P,E16.2)') NAME1,VAR
	    ENDIF
	  ENDDO
C
	  PRINT '(2X,A)','The following parameters have extreme correlations:'
	  WRITE(10,'(1X,A)') 'The following parameters have extreme correlations:'
	  DO I1=1,NPARS
	    DO I2=I1+1,NPARS
	      VAR=COVAR(I1+(I2-1)*NPARS)/SQRT(COVAR(I1+(I1-1)*NPARS)*COVAR(I2+(I2-1)*NPARS))
	      IF(ABS(VAR) .GT. 0.999) THEN
	        CALL GET_PAR_NAME(NAME1, I1)
	        CALL GET_PAR_NAME(NAME2, I2)
		    PRINT '(4X,3A,F10.4,A)',TRIM(NAME1),', ',TRIM(NAME2),VAR*100.0,'%'
		    WRITE(10,'(3X,3A,F10.4,A)') TRIM(NAME1),', ',TRIM(NAME2),VAR*100.0,'%'
	      ENDIF
	    ENDDO
	  ENDDO
	  STOP 'ERROR: Least-squares failed at stationary point'
	ELSE IF(IERR .EQ. 2) THEN		! Reached maximum number of iterations
	  PRINT '(/,1X,A,/)','WARNING: Refinement not fully converged'
	  WRITE(10,'(/,A,/)') 'WARNING: Refinement not fully converged'
	ELSE IF(IERR .EQ. 3) THEN		! Damping factor too small
	  PRINT '(/,1X,A,/)','WARNING: Refinement has poor convergence'
	  WRITE(10,'(/,A,/)') 'WARNING: Refinement has poor convergence'
	ELSE		! Trap all remaining errors and incorrect parameter values
C Give explicit complaints if NLSCON is given the wrong values, then stop
	  IF(IERR .EQ. 20) STOP 'BUG(start_nlscon): Invalid N,M,MFIT'
	  IF(IERR .EQ. 21) STOP 'BUG(start_nlscon): RTOL < 0'
	  IF(IERR .EQ. 22) STOP 'BUG(start_nlscon): XSCAL() < 0'
	  IF(IERR .EQ. 30) STOP 'BUG(start_nlscon): Invalid IOPT()'
C Not sure what is wrong, so print out some information and stop
	  PRINT '(/,1X,A,I4,A)','Total of',IW(5)+IW(8),' function calls by NLSCON'
	  PRINT *,'Error',IERR,' returned from NLSCON'
	  WRITE(10,'(/,A,I4,A)') 'Total of',IW(5)+IW(8),' function calls by NLSCON'
	  WRITE(10,'(A,I4,A)') 'Error',IERR,' returned from NLSCON'
	  STOP 'ERROR: Least-squares failed' 
	ENDIF
C
	RETURN
C
C Specific error message for work arrays too small
C
999	WRITE(10,'(A,I4,A)') "BUG(start_nlscon): Work arrays too small"
	WRITE(10,'(A,2I9,A,2I9)') 'Using',IIW,IRW,'  Needs', IW(18),IW(19)
	STOP "BUG(start_nlscon): Work arrays too small"
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
C     RW(21)=1.0D0
C     Override minimal allowed damping factor:
C     RW(22)=1.0D-3
C     Override rank1-decision parameter SIGMA:
C     RW(23)=2.0D0
C     Override maximum permitted subcondition for DECCON:
      RW(25)= 1.0D+16
C
	RETURN
	END



      SUBROUTINE NLSCON_DUMMY_FUNC(N,M,MCON,X,DFX,IFAIL)
C Dummy function for Jacobian evaluation.
	PRINT *,'N,M,MCON,X,DFX,IFAIL=',N,M,MCON,X,DFX,IFAIL
	STOP 'BUG(nlscon_dummy_func): Function should never be called!'
      END



	SUBROUTINE NLSCON_FUNC(NPARS, NOBS,NCON, PAR,ICALC, IFAIL)
C
	REAL PAR(*),ICALC(*)
C
C A sanity checks (There should be no constraints)
C
	IF(NCON .NE. 0) STOP 'BUG(f): NCON .NE. 0'
C
	CALL LSQ_FUNC(ICALC,NOBS, PAR,NPARS)
C
	IFAIL=0
	RETURN
	END
