C----- Bunch of miscellaneous routines ----------------------------------
c	SUBROUTINE CALC_MOD_QVEC(QVEC,IMOD)
c	SUBROUTINE GET_HOST_DATE(HOST_NAME2,DATE2)
c	SUBROUTINE GET_NUMXY(NUMX2,NUMY2)
c	SUBROUTINE MAT_SOLVE_EIGEN(V_OUT, AMAT,V_IN,NDIM,NDIM_MAX)
c	SUBROUTINE MAT_INVERT_EIGEN(AMAT_INV, AMAT,NDIM)
c	SUBROUTINE MAT_SOLVE(V_OUT, MAT,V_IN,NDIM,NDIM_MAX)
c	SUBROUTINE FIT_POLY_ORTHO(COEFFS,NORDER, YDATA,NDATA)
c	FUNCTION POLY_ORTHO(X,IORDER)
c	SUBROUTINE ASK_YES_NO(IANSWER,STRING)
c	FUNCTION CALC_LIN_SPLINE(X, XSTEP,TABLE,NTABLE)
c	FUNCTION CALC_QUAD_SPLINE(X, XSTEP,Y,NY)
c	SUBROUTINE UP2LOW_CASE(OUT, IN)
c	FUNCTION CUMUL_INV(P)
C------------------------------------------------------------------------


C----- Bunch of miscellaneous routines ----------------------------------

	SUBROUTINE CALC_MOD_QVEC(QVEC,IMOD)
C
C Calculate the Q-vector (satellite offset) for modulation number IMOD
C IMOD < 1 returns (0,0,0)
C
	REAL QVEC(3)
C
	REAL MOD_VEC
	COMMON /MODULATE_COM/ NMOD_VEC,NMOD_IDX,MOD_VEC(3,100),MOD_IDX(100,100)
C
	QVEC(1)=0.0
	QVEC(2)=0.0
	QVEC(3)=0.0
	IF(IMOD .GE. 1) THEN
	  DO I1=1,3
	    QVEC(I1)=0.0
	    DO I2=1,NMOD_VEC
	      QVEC(I1)=QVEC(I1) + MOD_IDX(IMOD,I2)*MOD_VEC(I1,I2)
	    ENDDO
	  ENDDO
	ENDIF
C
	RETURN
	END


	SUBROUTINE GET_HOST_DATE(HOST_NAME2,DATE2)
C
	CHARACTER*10 HOST_NAME2,DATE2
C
	CHARACTER*10 HOST_NAME,DATE
	COMMON /FILE_INFO_COM/ NUMX,NUMY,HOST_NAME,DATE
C
	HOST_NAME2=HOST_NAME
	DATE2=DATE
C
	RETURN
	END


	SUBROUTINE GET_NUMXY(NUMX2,NUMY2)
C
	CHARACTER*10 HOST_NAME,DATE
	COMMON /FILE_INFO_COM/ NUMX,NUMY,HOST_NAME,DATE
C
	NUMX2=NUMX
	NUMY2=NUMY
C
	RETURN
	END


	SUBROUTINE MAT_SOLVE_EIGEN(V_OUT, AMAT,V_IN,NDIM,NDIM_MAX)
C
	REAL AMAT(NDIM_MAX,NDIM),V_IN(NDIM),V_OUT(NDIM)
C
	REAL*8 AMAT2(NDIM,NDIM),AINV(NDIM,NDIM)
C
	DO I1=1,NDIM
	  DO I2=1,NDIM
	    AMAT2(I1,I2)=AMAT(I1,I2)
	  ENDDO
	ENDDO
C
	CALL MAT_INVERT_EIGEN(AINV, AMAT2,NDIM)
C
	DO I1=1,NDIM
	  V_OUT(I1)=0.0
	  DO I2=1,NDIM
	    V_OUT(I1)=V_OUT(I1)+AINV(I1,I2)*V_IN(I2)
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE MAT_INVERT_EIGEN(AMAT_INV, AMAT,NDIM)
C
	REAL AMAT(NDIM,NDIM),AMAT_INV(NDIM,NDIM)
C
	REAL*8 EVEC(NDIM,NDIM),EVAL(NDIM),WORK(NDIM*(NDIM+2))
C
	DO I1=1,NDIM
	  DO I2=1,NDIM
	    EVEC(I1,I2)=AMAT(I1,I2)
	  ENDDO
	ENDDO
C
	NWORK=NDIM*(NDIM+2)
	CALL DSYEV('V','U', NDIM,EVEC,NDIM,EVAL,WORK,NWORK,INFO)
	IF(INFO .NE. 0) STOP 'DSYEV failed'
C
C	PRINT '(1X,A,1P,10G10.2)','Eigenvalues:',(EVAL(K),K=1,NDIM)
C
	DO I1=1,NDIM
	  DO I2=1,NDIM
	    AMAT_INV(I1,I2)=0.0
	    DO I3=1,NDIM
	      AMAT_INV(I1,I2)=AMAT_INV(I1,I2)+EVEC(I1,I3)*EVEC(I2,I3)/EVAL(I3)
	    ENDDO
	  ENDDO
	ENDDO	
C
	RETURN
	END


	SUBROUTINE MAT_SOLVE(V_OUT, MAT,V_IN,NDIM,NDIM_MAX)
C
	REAL MAT(NDIM_MAX,NDIM_MAX),V_IN(NDIM_MAX),V_OUT(NDIM_MAX)
C
	INTEGER IPIVOTS(NDIM)
	REAL A_MAT(NDIM,NDIM),B_VEC(NDIM)
C
C Setup the matrices as DGBSV wants it (NB: change to double precision)
C
	DO I1=1,NDIM
	  DO I2=1,NDIM
	    A_MAT(I1,I2)=MAT(I1,I2)
	  ENDDO
	  B_VEC(I1)=V_IN(I1)
C	  IF(ABS(A_MAT(I1,I1)) .LT. 1E-6) A_MAT(I1,I1)=1.0
	ENDDO
C
C Solve the linear equations and complain and die if it unexpectedly fails
C
	CALL DGESV(NDIM, 1, A_MAT, NDIM, IPIVOTS, B_VEC, NDIM, ISTATUS)
	IF(ISTATUS .NE. 0) STOP	'BUG: DGESV failed'
C
C Copy the diagonal values into TRACE()
C
	DO I=1,NDIM
	  V_OUT(I)=B_VEC(I)
	ENDDO
C
	RETURN
	END


	SUBROUTINE FIT_POLY_ORTHO(COEFFS,NORDER, YDATA,NDATA)
C
C Use Legendre orthogonal polynomials of order NORDER to fit the Y values
C in YDATA(1..NDATA) versus X values of -1 to +1, corresponding to the
C first and last elements of YDATA().
c NB: Assumes data is regular
C
	REAL YDATA(NDATA),COEFFS(0:NORDER)
C
C Sanity check
C
	IF(NORDER .GT. 8) STOP 'BUG(fit_poly_ortho): NORDER > 8'
C
C Zero the polynomial coefficients
C
	DO I=0,NORDER
	  COEFFS(I)=0.0
	ENDDO
C
C Create sums to approximate correlations of Y values with polynomials
C|
	DX=2.0/(NDATA-1)
	DO I=1,NDATA
	  X=-1.0+DX*(I-1)
	  DO ISUM=0,NORDER
	    COEFFS(ISUM)=COEFFS(ISUM)+POLY_ORTHO(X,ISUM)*YDATA(I)
	  ENDDO
	ENDDO
C
C Calculate polynomial coefficients from orthogonality rules
C
	DO I=0,NORDER
	  COEFFS(I)=(2*I+1)*COEFFS(I)/NDATA
	ENDDO
C
	RETURN
	END


	FUNCTION POLY_ORTHO(X,IORDER)
C
C Orthogonal (Legendre) polynomials for order 0 to 10
C
	IF(IORDER .EQ. 0) THEN
	  Y=1.0
	ELSE IF(IORDER .EQ. 1) THEN
	  Y=X
	ELSE IF(IORDER .EQ. 2) THEN
	  Y=(3.0*X**2-1.0)/2.0
	ELSE IF(IORDER .EQ. 3) THEN
	  Y=((5.0*X**2-3.0)*X)/2.0
	ELSE IF(IORDER .EQ. 4) THEN
	  Y=((35.0*X**2-30.0)*X**2+3.0)/8.0
	ELSE IF(IORDER .EQ. 5) THEN
	  Y=(((63.0*X**2-70.0)*X**2+15.0)*X)/8.0
	ELSE IF(IORDER .EQ. 6) THEN
	  Y=(((231.0*X**2-315.0)*X**2+105.0)*X**2-5.0)/16.0
	ELSE IF(IORDER .EQ. 7) THEN
	  Y=((((429.0*X**2-693.0)*X**2+315.0)*X**2-35.0)*X)/16.0
	ELSE IF(IORDER .EQ. 8) THEN
	  Y=((((6435.0*X**2-12012.0)*X**2+6930.0)*X**2-1260.0)*X**2+35.0)/128.0
	ELSE IF(IORDER .EQ. 9) THEN
	  Y=((((12155.0*X**2-25740.0)*X**2+18018.0)*X**2-4620.0)*X**2+315.0)*X/128.0
	ELSE IF(IORDER .EQ. 10) THEN
	  Y=(((((46189.0*X**2-109395.0)*X**2+90090.0)*X**2-30030.0)*X**2+3465.0)*X**2-63.0)/256.0
	ELSE
	  STOP 'BUG(poly_ortho): IORDER outside range of 0 - 10'
	ENDIF
C
	POLY_ORTHO=Y
	RETURN
	END


	SUBROUTINE ASK_YES_NO(IANSWER,STRING)
C
C Ask a yes/no question with STRING as the query
C Returns IANSWER=1 if yes, =0 if no, =-1 anything else
C
	CHARACTER*(*) STRING
C
	CHARACTER*1 KEY
C
C
	PRINT '(1X,2A,$)',TRIM(STRING),' [y/n] '
	READ(*,'(A1)',ERR=1000) KEY
	IF(KEY.EQ.'Y' .OR. KEY.EQ.'y') IANSWER=1
	IF(KEY.EQ.'N' .OR. KEY.EQ.'n') IANSWER=0
C
1000	RETURN
	END


	FUNCTION CALC_LIN_SPLINE(X, XSTEP,TABLE,NTABLE)
C
C For a table of data TABLE(1..NTABLE) linearly interpolate the
C values for the variable X. X=(N-1)*XSTEP refers to the Nth value
C in the table, that is the first value corresponds to X=0.
C Beyond the table limits the values are just the end values.
C
	REAL TABLE(NTABLE)
C
C NB: IX starts from 0 not 1
	IX=INT(X/XSTEP)
	IF(IX .LT. 0) THEN
	  CALC_LIN_SPLINE=TABLE(1)
	ELSE IF(IX .GT. NTABLE-2) THEN
	  CALC_LIN_SPLINE=TABLE(NTABLE)
	ELSE
	  FREL=(X-IX*XSTEP)/XSTEP
C NB: IX starts from 0 not 1
	  CALC_LIN_SPLINE=(1.0-FREL)*TABLE(IX+1)+FREL*TABLE(IX+2)
	ENDIF
C
	RETURN
	END


	FUNCTION CALC_QUAD_SPLINE(X, XSTEP,Y,NY)
C
C For a table of data, Y(1..NY), return a piecewise quadratic spline
C that approximates the Y() values for values of X. The approximation
C improves with the smoothness of the data. X=(N-1)*XSTEP corresponds
C to the Nth value in the table, so X=0 is the first value in the table.
C Beyond the table limits the values are just the end values.
C
C The function values and derivatives at midpoints is exact and given
C by: F = (Y(N+2)+Y(N+1))/2  ,  F' = (Y(N+2)-Y(N+1))/XSTEP
C
	REAL Y(NY)
C
	IF(NY .LT. 3) STOP 'BUG(calc_quad_spline): NY < 3'
C
C Determine which knot we will use and relative shift in X
C
	IX=FLOOR(X/XSTEP+0.5)
	DX=X/XSTEP+0.5-IX
C
C Evaluate spline for general and special cases
C
C General case, calculate coeffs and evaluate quadratic
	IF(IX.GT.0 .AND. IX.LT.NY-1) THEN
	  A0=(Y(IX)+Y(IX+1))/2.0
	  A1=Y(IX+1)-Y(IX)
	  A2=(Y(IX)+Y(IX+2))/2.0-Y(IX+1)
	  CALC_QUAD_SPLINE=(A2*DX+A1)*DX+A0
C If beyond table limits, fix to the limit values
	ELSE IF(IX .LT. 0) THEN
	  CALC_QUAD_SPLINE=Y(1)
	ELSE IF(IX .GT. NY-1) THEN
	  CALC_QUAD_SPLINE=Y(NY)
C If at table limits, use quadratic decay to the limit values
	ELSE IF(IX .EQ. 0) THEN
	  CALC_QUAD_SPLINE=Y(1)+0.5*(Y(2)-Y(1))*DX**2
	ELSE IF(IX .EQ. NY-1) THEN
	  CALC_QUAD_SPLINE=Y(NY)+0.5*(Y(NY-1)-Y(NY))*(1.0-DX)**2
C Huh, what?
	ELSE
	  STOP 'BUG(calc_quad_spline): Impossible IX'
	ENDIF
C
	RETURN
	END


	SUBROUTINE UP2LOW_CASE(OUT, IN)
C
	CHARACTER OUT*(*),IN*(*)
C
	CHARACTER CH*1
C
	IDIFF=ICHAR('a')-ICHAR('A')
	DO I=1,LEN(IN)
	  CH=IN(I:I)
	  IF(CH.GE.'A' .AND. CH.LE.'Z') CH=CHAR( ICHAR(CH)+IDIFF )
	  OUT(I:I)=CH
	ENDDO
C
	DO I=LEN(IN)+1,LEN(OUT)
	  OUT(I:I)=' '
	ENDDO
C
	RETURN
	END


	FUNCTION CUMUL_INV(P)
C
C Calculates the inverse cumulative distribution function
C
C Return the value (in sigma below the mean) where we expect a
C probability, P, of a normal distribution to be less than the value.
C
C ALGORITHM AS 111, APPL.STATIST., VOL.26, 118-121, 1977.
C
	DATA A0,A1,A2,A3/2.506628,-18.61500,41.39112,-25.44106/
	DATA B1,B2,B3,B4/-8.473511,23.08337,-21.06224,3.130829/
	DATA C0,C1,C2,C3/-2.787189,-2.297964,4.850141,2.321212/
	DATA D1,D2/3.543889,1.637068/
C
	IF(P.GT.0.08 .AND. P.LT.0.92) THEN		! 0.08 < P < 0.92
	  R=(P-0.5)**2
	  CUMUL_INV= (P-0.5)*(((A3*R + A2)*R + A1)*R + A0) / 
	1		(1.0 + (((B4*R + B3)*R + B2)*R + B1)*R)
	ELSE								! P < 0.08 OR P > 0.92
	  R=MAX(1E-8, MIN(P,1.0-P) )
	  R=SQRT(-LOG(R))
	  CUMUL_INV= (((C3*R + C2)*R + C1)*R + C0)/(1.0 + (D2*R + D1)*R)
	  IF(P .LT. 0.5) CUMUL_INV = -CUMUL_INV
	ENDIF
C
	RETURN
	END
