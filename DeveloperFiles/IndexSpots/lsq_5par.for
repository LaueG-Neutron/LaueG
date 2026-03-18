	SUBROUTINE ITER_LSQ(UB,DXY,RMS,NSUM, XSPOT,YSPOT,HKLS,CUTOFF,NPEAKS)
C
	REAL UB(3,3),DXY(2),XSPOT(NPEAKS),YSPOT(NPEAKS),HKLS(3,NPEAKS)
C
	PARAMETER (NHKLS_MAX=1000)
	COMMON /LSQ_COM/ HINV(5,5),JAC(2*NHKLS_MAX,5)
C
	LOGICAL LEXCLUDE(NHKLS_MAX)
	REAL RES_IN(2*NHKLS_MAX),RES_OUT(2*NHKLS_MAX)
	REAL HVEC(3),VEC5(5),DPARS(5),ROT(3,3)
C
C Calculate the initial residual (X/Y obs - calc), and
C those peaks to not refine LEXCLUDE() because outside CUTOFF
C
	NSUM=0
	DO I=1,NPEAKS
	  CALL VC3_MULT(HVEC, UB,HKLS(1,I))
	  CALL CALC_HVEC_TO_PIXEL(X0,Y0, HVEC)
	  RES_IN(2*I-1)=XSPOT(I)-X0-DXY(1)
	  RES_IN(2*I  )=YSPOT(I)-Y0-DXY(2)
	  DSQ=RES_IN(2*I-1)**2 + RES_IN(2*I)**2
	  LEXCLUDE(I)=(DSQ .GT. CUTOFF**2)
	  IF( LEXCLUDE(I) ) THEN
	    RES_IN(2*I-1)=0.0
	    RES_IN(2*I  )=0.0
	  ELSE
	    NSUM=NSUM+1
	  ENDIF
	ENDDO
C
	IF(NSUM .LT. 3) THEN
	  NSUM=0
	  RETURN
	ENDIF
C
C Do the maths to solve the LSQ
C
C VEC5 = JAC' * RES_IN
	DO I1=1,5
	  VEC5(I1)=0.0
	  DO I2=1,2*NPEAKS
	    VEC5(I1)=VEC5(I1) +JAC(I2,I1)*RES_IN(I2)
	  ENDDO
	ENDDO
C DPARS = HINV * (JAC' * RES_IN)
	DO I1=1,5
	  DPARS(I1)=0.0
	  DO I2=1,5
	    DPARS(I1)=DPARS(I1) +HINV(I1,I2)*VEC5(I2)
	  ENDDO
	ENDDO
C Approx. compensation of DPARS for the excluded peaks
	DO I=1,5
	  DPARS(I)=DPARS(I)*NPEAKS/NSUM
	ENDDO
C RES_OUT = JAC * (HINV * JAC' * RES_IN)
	DO I1=1,2*NPEAKS
	  RES_OUT(I1)=0.0
	  DO I2=1,5
	    RES_OUT(I1)=RES_OUT(I1) + JAC(I1,I2)*DPARS(I2)
	  ENDDO
	ENDDO
C
C Calculate the approximate rms from the IN & OUT residuals
C
	RMS=0.0
	DO I=1,NPEAKS
	  IF( .NOT.LEXCLUDE(I) ) RMS=RMS + (RES_IN(2*I-1)-RES_OUT(2*I-1))**2 +
	1									   (RES_IN(2*I)-RES_OUT(2*I))**2
	ENDDO
	RMS=SQRT(RMS/NSUM)
C
C Add LSQ results to DXY
C
	DXY(1)=DXY(1)+DPARS(1)
	DXY(2)=DXY(2)+DPARS(2)
C
C Apply LSQ rotation (in degrees) to UB
C NB: Not an exact rotation, we'll fix this later
C
      SIN_XROT=DPARS(3)/57.296
      SIN_YROT=DPARS(4)/57.296
      SIN_ZROT=DPARS(5)/57.296
	ROT(1,1)=1.0
	ROT(1,2)=			-SIN_ZROT
	ROT(1,3)=						-SIN_YROT
	ROT(2,1)=SIN_ZROT
	ROT(2,2)=			1.0
	ROT(2,3)=						-SIN_XROT
	ROT(3,1)=SIN_YROT
	ROT(3,2)=			SIN_XROT
	ROT(3,3)=						1.0
C
	CALL MX3_MULT(UB, ROT,UB)
C
	RETURN
	END



	SUBROUTINE SETUP_LSQ(XSPOT,YSPOT,NPEAKS)
C
	REAL XSPOT(NPEAKS),YSPOT(NPEAKS)
C
	PARAMETER (NHKLS_MAX=1000)
	COMMON /LSQ_COM/ HINV(5,5),JAC(2*NHKLS_MAX,5)
C
	REAL ROT(3,3),HVEC(3),HVEC2(3),HESS(5,5)
C
C JAC(IPAR,I) stores derivatives (dX/dPar or dY/dPar) for parameter IPAR,
C where odd/even I corresponds to dX/dY of peak number (I+1)/2.
C Parameters 1..5 are DXCEN,DYCEN,XROT,YROT,ZROT
C
C These rotations (in degrees) produce spot movements of ~1 pixel
C
	XROT=0.035
	YROT=0.044
	ZROT=0.085
C
C Calculate the Jacobian
C
	DO I=1,NPEAKS
C DXCEN
	  JAC(2*I-1,1)=1.0
	  JAC(2*I ,1)=0.0
C DYCEN
	  JAC(2*I-1,2)=0.0
	  JAC(2*I  ,2)=1.0
C Calculate the scattering vector HVEC
	  CALL CALC_PIXEL_TO_HVEC(HVEC, XSPOT(I),YSPOT(I))
C XROT
	  CALL CALC_XROT_MX(ROT, XROT)
	  CALL VC3_MULT(HVEC2, ROT,HVEC)
	  CALL CALC_HVEC_TO_PIXEL(X,Y, HVEC2)
	  JAC(2*I-1,3)=(X-XSPOT(I))/XROT
	  JAC(2*I  ,3)=(Y-YSPOT(I))/XROT
C YROT
	  CALL CALC_YROT_MX(ROT, YROT)
	  CALL VC3_MULT(HVEC2, ROT,HVEC)
	  CALL CALC_HVEC_TO_PIXEL(X,Y, HVEC2)
	  JAC(2*I-1,4)=(X-XSPOT(I))/YROT
	  JAC(2*I  ,4)=(Y-YSPOT(I))/YROT
C ZROT
	  CALL CALC_ZROT_MX(ROT, ZROT)
	  CALL VC3_MULT(HVEC2, ROT,HVEC)
	  CALL CALC_HVEC_TO_PIXEL(X,Y, HVEC2)
	  JAC(2*I-1,5)=(X-XSPOT(I))/ZROT
	  JAC(2*I  ,5)=(Y-YSPOT(I))/ZROT
C
	ENDDO
C
C Calculate the Hessian from the Jacobian
C
	DO I1=1,5
	  DO I2=1,5
	    HESS(I1,I2)=0.0
	    DO I3=1,2*NPEAKS
	      HESS(I1,I2)=HESS(I1,I2)+JAC(I3,I1)*JAC(I3,I2)
	    ENDDO
	  ENDDO
	ENDDO
C
C Invert the Hessian
C
	CALL MX5_INVERT(HINV, HESS,0.1)
C
	RETURN
	END



	SUBROUTINE MX5_INVERT(HINV, HESS,CUTOFF)
C
	REAL HINV(5,5),HESS(5,5)
C
	REAL E_VAL(5),E_VEC(5,5),MX5(5,5)
C
C DXCEN=1.0 DYCEN=1.0 XROT=0.035 YROT=0.044 ZROT=0.085
	REAL PARS_NORM(5)
	DATA PARS_NORM/1.0, 1.0, 0.035, 0.044, 0.085/
C
	DO I1=1,5
	  DO I2=1,5
	    HESS(I1,I2)=HESS(I1,I2)*PARS_NORM(I1)*PARS_NORM(I2)
	  ENDDO
	ENDDO
C
	CALL MX5_EIGEN(E_VAL,E_VEC, HESS)
C
	DO I1=1,5
	  DO I2=1,5
	    HINV(I1,I2)=0.0
	  ENDDO
	  IF(ABS(E_VAL(I1)) .LT. CUTOFF) E_VAL(I1)=CUTOFF
	  HINV(I1,I1)=1.0/E_VAL(I1)
	ENDDO
C
	DO I1=1,5
	  DO I2=1,5
	    MX5(I1,I2)=0.0
	    DO I3=1,5
	      MX5(I1,I2)=MX5(I1,I2)+HINV(I1,I3)*E_VEC(I2,I3)
	    ENDDO
	  ENDDO
	ENDDO
C
	DO I1=1,5
	  DO I2=1,5
	    HINV(I1,I2)=0.0
	    DO I3=1,5
	      HINV(I1,I2)=HINV(I1,I2)+E_VEC(I1,I3)*MX5(I3,I2)
	    ENDDO
	    HINV(I1,I2)=HINV(I1,I2)
	  ENDDO
	ENDDO
C
	DO I1=1,5
	  DO I2=1,5
	    HINV(I1,I2)=HINV(I1,I2)*PARS_NORM(I1)*PARS_NORM(I2)
	  ENDDO
	ENDDO
C
	RETURN
	END



	SUBROUTINE MX5_EIGEN(VALUES,VECTORS, MATRIX)
C
C Calculate eigen-values VALUES() and eigen-vectors VECTORS()
C for the 5 x 5 matrix MATRIX().
C
	REAL VALUES(5),VECTORS(5,5),MATRIX(5,5)
C
C Looks like we only need WORK() to be 5*N+1, but use 100 to be safe.
C
	REAL WORK(100),MATRIX2(5,5)
C
C Make a copy of MATRIX() so that EA06C doesn't corrupt it.
C
	DO I1=1,5
		DO I2=1,5
			MATRIX2(I1,I2)=MATRIX(I1,I2)
		ENDDO
	ENDDO
C
	CALL EA06C(MATRIX2,VALUES,VECTORS,5,5,5,WORK)
	RETURN
	END
