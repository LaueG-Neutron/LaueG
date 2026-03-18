C ***************** This file must be compiled with R8 as the default *************

C ===================== Routines in this file ==============================
C ---------------- Calculate the EFH ---------------------
C	SUBROUTINE CALC_FITTED_EFH(E4,F4,H4, X4,Y4)
C	SUBROUTINE CALC_SPOT_EFH(E,F,H, XY_SPOT, PIX_CEN,PIX_SIZE,DRUM_RAD, EIJ,ETA,ETA_ISO)
C	SUBROUTINE CALC_SHAPE_EFH(E,F,H, ALPHA_X,ALPHA_Y, EIJ)
C	SUBROUTINE CALC_MOSAIC_EFH(E,F,H, ALPHA_X,ALPHA_Y,TTH, ETA,ETA_ISO,DRUM_PIX)
C ---------------- Fit EFH and return fitted values  -------------
C	SUBROUTINE FIT_EFH(PLIB4,NLIB)
C----------------- Run the LSQ refinement ---------------
C	SUBROUTINE RUN_LSQ(SSQ,NPARS, NSPOTS,LPRINT)
C----------------- The function called by LSQ refinement ---------------
C	SUBROUTINE LSQ_FUNC(CALC,NOBS, PARS,NPARS)
C-------------------- Load data used for LSQ refinement ------------------
C	SUBROUTINE LOAD_LSQ_DATA(OBS,SIG,NSPOTS)
C----------------- Parameter specific routines ---------------
C	SUBROUTINE CLEAR_LSQ_PARS
C	SUBROUTINE REFINE_LSQ_SHAPE
C	SUBROUTINE REFINE_LSQ_ANISO
C	SUBROUTINE REFINE_LSQ_ISO
C	SUBROUTINE PACK_LSQ_PARS(PARS,NPARS)
C	SUBROUTINE UNPACK_LSQ_PARS(PARS,ESDS,LPRINT)
C	SUBROUTINE GET_PAR_NAME(NAME, IPAR)
C	SUBROUTINE GET_PAR_ESDS(ESDS, COV,NPARS)
C-------------------- NLSCON SPECIFIC ROUTINES --------------------
C	SUBROUTINE START_NLSCON(COVAR,ISIG, PAR,NPARS, IOBS,NOBS)
C	SUBROUTINE SET_NLSCON_OPTIONS(IOPT,IW,RW, IIW,IRW)
C	SUBROUTINE NLSCON_DUMMY_FUNC(N,M,MCON,X,DFX,IFAIL)
C	SUBROUTINE NLSCON_FUNC(NPARS, NOBS,NCON, PAR,ICALC, IFAIL)
C ================= R*8 versions of some MX3UTILS.FOR routines =============
C	SUBROUTINE VC3_SET_R8(VC, V1,V2,V3)
C	SUBROUTINE MX3_SET_R8(MX, V11,V12,V13, V21,V22,V23, V31,V32,V33)
C	SUBROUTINE MX3_XROT_R8(ROT, ANG)
C	SUBROUTINE MX3_YROT_R8(ROT, ANG)
C	SUBROUTINE VC3_UNIT_R8(VC)
C	SUBROUTINE MX3_MULT_R8(M_OUT, M1,M2)
C	SUBROUTINE MX3_TRANS_R8(A,B)
C	FUNCTION VC3_DOT_R8(A,B)
C ============ OLD EFH calculation routines ==================
C	SUBROUTINE CALC_SPOT_EFH_OLD(E,F,H, XY_SPOT, XY_CEN,DRUM_PIX, EIJ,ETA,ETA_ISO)
C =========================================================================


C ---------------- Calculate the EFH ---------------------

	SUBROUTINE CALC_FITTED_EFH(E4,F4,H4, X4,Y4)
C
	REAL*4 E4,F4,H4,X4,Y4
C
	REAL*4 PIX_CEN4,PIX_SIZE4,DRUM_RAD4
	COMMON /GEOM_COM/ PIX_CEN4(2),PIX_SIZE4(2),DRUM_RAD4
C
	LOGICAL LPARS
	REAL*8 EIJ8,ETA8,ETA_ISO8
	COMMON /FIT_EFH_COM/ EIJ8(6),ETA8(3),ETA_ISO8,LPARS(50)
C
	REAL*8 E8,F8,H8,XY_SPOT8(2),PIX_CEN8(2),PIX_SIZE8(2),DRUM_RAD8
C
C Convert R*4 arguments to local R*8 variables
C
	XY_SPOT8(1)=X4
	XY_SPOT8(2)=Y4
C
C Convert R*4 common variables to local R*8 variables
C
	PIX_CEN8(1)=PIX_CEN4(1)
	PIX_CEN8(2)=PIX_CEN4(2)
	PIX_SIZE8(1)=PIX_SIZE4(1)
	PIX_SIZE8(2)=PIX_SIZE4(2)
	DRUM_RAD8=DRUM_RAD4
C
C Calculate e,f,h using arguments of R*8
C
	CALL CALC_SPOT_EFH(E8,F8,H8, XY_SPOT8, PIX_CEN8,PIX_SIZE8,DRUM_RAD8, EIJ8,ETA8,ETA_ISO8)
C
C Convert output arguments from R*8 to R*4 arguments
C
	E4=E8
	F4=F8
	H4=H8
C
	RETURN
	END


	SUBROUTINE CALC_SPOT_EFH(E,F,H, XY_SPOT, PIX_CEN,PIX_SIZE,DRUM_RAD, EIJ,ETA,ETA_ISO)
C
C NB: E,F,H have the units of 1/(X pixel size)**2
C
	REAL XY_SPOT(2),PIX_CEN(2),PIX_SIZE(2),EIJ(6),ETA(3)
C
C Calculate basic angles
C
	ALPHA_X=      (XY_SPOT(1)-PIX_CEN(1))*PIX_SIZE(1)/DRUM_RAD
	ALPHA_Y=ATAN( (XY_SPOT(2)-PIX_CEN(2))*PIX_SIZE(2)/DRUM_RAD )
	TTH=ACOS( COS(ALPHA_Y) * COS(ALPHA_X) )
C
C Calculate ellipse parameters (E1,F1,H1) due to crystal shape
C
	CALL CALC_SHAPE_EFH(E1,F1,H1, ALPHA_X,ALPHA_Y, EIJ)
C
C Calculate ellipse parameters (E2,F2,H2) due to mosaicity
C
	DRUM_XPIX=DRUM_RAD/PIX_SIZE(1)
	CALL CALC_MOSAIC_EFH(E2,F2,H2, ALPHA_X,ALPHA_Y,TTH, ETA,ETA_ISO,DRUM_XPIX)
C
C Convolute the two ellipses
C
	E=( E1*E2*(F1+F2) - E1*H2**2 - E2*H1**2 ) / ( (E1+E2)*(F1+F2) - (H1+H2)**2 )
	F=( F1*F2*(E1+E2) - F1*H2**2 - F2*H1**2 ) / ( (E1+E2)*(F1+F2) - (H1+H2)**2 )
	H=( H1*(E2*F2-H2**2) + H2*(E1*F1-H1**2) ) / ( (E1+E2)*(F1+F2) - (H1+H2)**2 )
C
C Convert final ellipse to "cylindrical detector" coords
C
	F=F*COS(ALPHA_Y)**2
	H=H*COS(ALPHA_Y)
C
	RETURN
	END


	SUBROUTINE CALC_SHAPE_EFH(E,F,H, ALPHA_X,ALPHA_Y, EIJ)
C
	REAL EIJ(6)		! contains in order E11,E22,E33,E12,E13,E23
C
	REAL EMAT(3,3),UMAT(3,3),RX(3,3),RY(3,3)
C
C Load 3D ellipsoidal shape as a 3x3 matrix
C KLUDGE 1: Force Eii > 0
C
CC	CALL MX3_SET_R8(EMAT,  ABS(EIJ(1)),EIJ(4),EIJ(5),  EIJ(4),ABS(EIJ(2)),EIJ(6),  EIJ(5),EIJ(6),ABS(EIJ(3))  )
	CALL MX3_SET_R8(EMAT,  EIJ(1),EIJ(4),EIJ(5),  EIJ(4),EIJ(2),EIJ(6),  EIJ(5),EIJ(6),EIJ(3)  )
C
C Rotation matrix for lab to "spherical detector" coords
C
	CALL MX3_YROT_R8(RY,-ALPHA_X)
	CALL MX3_XROT_R8(RX,-ALPHA_Y)
	CALL MX3_MULT_R8(UMAT, RY,RX)
C
C Convert to spherical detector coords using EMAT=UMAT'*(EMAT*UMAT)
C
	CALL MX3_MULT_R8(EMAT, EMAT,UMAT)
	CALL MX3_TRANS_R8(UMAT, UMAT)
	CALL MX3_MULT_R8(EMAT, UMAT,EMAT)
C
C Calculate 2D ellipse parameters from ellipsoidal matrix
C
	E=EMAT(1,1)-EMAT(1,3)*EMAT(1,3)/EMAT(3,3)
	F=EMAT(2,2)-EMAT(2,3)*EMAT(2,3)/EMAT(3,3)
	H=EMAT(1,2)-EMAT(1,3)*EMAT(2,3)/EMAT(3,3)
C
C KLUDGE 2: Force ellipse to be non-negative
C
CC	E=ABS(E)
CC	F=ABS(F)
CC	IF(H**2 .GT. E*F) H=H*SQRT(E*F/H**2)
C
C Always add 0.5 to E & F, which is equivalent to a 0.1 mm radius crystal
C
	E=E+0.5
	F=F+0.5
C
	RETURN
	END


	SUBROUTINE CALC_MOSAIC_EFH(E,F,H, ALPHA_X,ALPHA_Y,TTH, ETA,ETA_ISO,DRUM_XPIX)
C
	REAL ETA(3)
C
	REAL SPERP(3),S10(3)
C
C Create unit vector along S0 x S1
C
	CALL VC3_SET_R8(SPERP, -SIN(ALPHA_Y), SIN(ALPHA_X)*COS(ALPHA_Y), 0.0 )
	CALL VC3_UNIT_R8(SPERP)
C
C Create unit vector along S1 - S0
C
	CALL VC3_SET_R8(S10, SIN(ALPHA_X)*COS(ALPHA_Y), SIN(ALPHA_Y), COS(ALPHA_X)*COS(ALPHA_Y) )
	S10(3)=S10(3)-1.0
	CALL VC3_UNIT_R8(S10)
C
C Get components of 3D mosaicity along SPERP & S10
C
	ETA_PERP=VC3_DOT_R8(ETA,SPERP)
	ETA_10  =VC3_DOT_R8(ETA,S10)
C
C Add isotropic contribution to components and convert to radians
C Calculate ellipse parameters for "spherical detector" coords
C
	A=2.0*DRUM_XPIX*(ABS(ETA_PERP)+ABS(ETA_ISO))/57.296             ! 57.296 as ETA are in degrees
	B=2.0*DRUM_XPIX*(SIN(TTH*0.5)*ABS(ETA_10)+ABS(ETA_ISO))/57.296
	A=MAX(0.01,A)
	B=MAX(0.01,B)
C
	SIN_ANGLE=COS(ALPHA_X)*SIN(ALPHA_Y)/SIN(TTH)
	COS_ANGLE=SIN(ALPHA_X)/SIN(TTH)
	E=(COS_ANGLE/A)**2 + (SIN_ANGLE/B)**2
	F=(SIN_ANGLE/A)**2 + (COS_ANGLE/B)**2
	H=(A**(-2) - B**(-2))*SIN_ANGLE*COS_ANGLE
C
	RETURN
	END


C ---------------- Fit EFH and return fitted values  -------------

	SUBROUTINE FIT_EFH(PLIB4,NLIB)
C
	REAL*4 PLIB4(1000,6)
C
	REAL*4 PIX_CEN4,PIX_SIZE4,DRUM_RAD4
	COMMON /GEOM_COM/ PIX_CEN4(2),PIX_SIZE4(2),DRUM_RAD4
C
	COMMON /FIT_EFH_DATA_COM/ EFH(3,1000),XYPIX(2,1000),PIX_CEN(2),PIX_SIZE(2),DRUM_RAD
C
	LOGICAL LPARS
	COMMON /FIT_EFH_COM/ EIJ(6),ETA(3),ETA_ISO,LPARS(50)
C
	PARAMETER (NDATA_MAX=3000)
	COMMON /LSQ_FUNC_COM/ CALC(NDATA_MAX),OBS(NDATA_MAX),SIG(NDATA_MAX)
C
C Copy R*4 variables in /GEOM_COM/ to local R*8 variables
C
	PIX_CEN(1)=PIX_CEN4(1)
	PIX_CEN(2)=PIX_CEN4(2)
	PIX_SIZE(1)=PIX_SIZE4(1)
	PIX_SIZE(2)=PIX_SIZE4(2)
	DRUM_RAD=DRUM_RAD4
C
C Copy R*4 argument PLIB4() to R*8 variables XYPIX() & EFH()
C
	DO I=1,NLIB
	  XYPIX(1,I)=PLIB4(I,1)
	  XYPIX(2,I)=PLIB4(I,2)
	  EFH(1,I)=PLIB4(I,3)
	  EFH(2,I)=PLIB4(I,4)
	  EFH(3,I)=PLIB4(I,5)
	ENDDO
C
C Try refining 7 different combinations of parameters and record the best one
C Add a bias to favour the models with less parameters
C On the 8th loop, switch to the best one so to load its parameters
C
	WRITE(30,'(/,A)') 'Fitting EFH values using least-squares'
	SSQ_BEST=1E9
	DO ITRY=1,8
	  I=ITRY
C On ITRY=8, reload the best model from ITRY=1-7
	  IF(ITRY .EQ. 8) THEN
	    I=IBEST
	    WRITE(30,'(A,I2,A)') 'Using Model',IBEST,', Parameters:'
	  ENDIF
C Setup the parameters to refine
	  CALL CLEAR_LSQ_PARS
	  IF(I.EQ.1 .OR. I.EQ.3 .OR. I.EQ.5 .OR. I.EQ.7) CALL REFINE_LSQ_ISO
	  IF(I.EQ.2 .OR. I.EQ.3 .OR. I.EQ.6 .OR. I.EQ.7) CALL REFINE_LSQ_ANISO
	  IF(I .GE. 4) CALL REFINE_LSQ_SHAPE
C If looping through models, output summary and record which is best
	  IF(ITRY .LT. 8) THEN
	    CALL RUN_LSQ(SSQ,NPARS, NLIB,.FALSE.)
	    WRITE(30,'(2X,A,I2,A,F6.1,A,I3,A)',IOSTAT=IDUMMY) 'Model',I,
	1					':   SSQ=',SSQ,' for',NPARS,' parameters'
C Add a bias to SSQ for less parameters
	    SSQ=SSQ*(1.0+NPARS/20.0)
	    IF(SSQ .LT. SSQ_BEST) THEN
	      IBEST=I
	      SSQ_BEST=SSQ
	    ENDIF
	  ELSE
C Refine model and output parameters plus esds
	    CALL RUN_LSQ(SSQ,NPARS, NLIB,.TRUE.)
	  ENDIF
C
	ENDDO
C
C Calculated the EFH and copy to the R*4 PLIB()
C
	DO I=1,NLIB
	  CALL CALC_SPOT_EFH(E,F,H, XYPIX(1,I), PIX_CEN,PIX_SIZE,DRUM_RAD, EIJ,ETA,ETA_ISO)
	  PLIB4(I,3)=E
	  PLIB4(I,4)=F
	  PLIB4(I,5)=H
	ENDDO
C
	RETURN
	END


C----------------- Run the LSQ refinement ---------------

	SUBROUTINE RUN_LSQ(SSQ,NPARS, NSPOTS,LPRINT)
C
	LOGICAL LPRINT
C
C Start the LSQ refinements of NSPOTS model ellipses.
C The number of refined parameters is returned in NPARS.
C In this case, all warnings (except software bugs) have been disabled
C and the sum-of-squared residuals are returned in SSQ. If a problem
C occurs in the refinements it is assumed it will be rejected for
C having a large SSQ value.
C
	PARAMETER (NDATA_MAX=3000)
C
	COMMON /LSQ_PARS_COM/ PARS(50),ESDS(50),COVAR(50,50)
C
	COMMON /LSQ_FUNC_COM/ CALC(NDATA_MAX),OBS(NDATA_MAX),SIG(NDATA_MAX)
C
C Sanity check on data size
C
	NOBS=3*NSPOTS
	IF(NOBS .GT. NDATA_MAX) CALL QUIT('BUG(run_lsq): NSPOTS > 1000')
C
C Load LSQ parameters into PARS, and check on the total number
C
	CALL PACK_LSQ_PARS(PARS,NPARS)
	IF(NPARS .GT. 50) CALL QUIT('BUG(run_lsq): NPARS > 50')
C
C Load the observations with esds
C
	CALL LOAD_LSQ_DATA(OBS,SIG,NSPOTS)
C
C Call routine to setup and start NLSCON refinement routines
C
	CALL START_NLSCON(COVAR,SIG, PARS,NPARS, OBS,NOBS)
C
C Calculate ICALC() with the final parameters.
C
	CALL LSQ_FUNC(CALC,NOBS, PARS,NPARS)
C
C Calculate the parameter esds and output any bad correlations
C
	CALL GET_PAR_ESDS(ESDS, COVAR,NPARS)
C
C Unpack parameters, and print values if LPRINT is true
C
	CALL UNPACK_LSQ_PARS(PARS,ESDS,LPRINT)
C
C Sum the parameter esds and set huge SSQ if sum > 1
C
	SUME=0.0
	DO I=1,NPARS
	  SUME=SUME+ABS(ESDS(I))
	ENDDO
C
	IF(SUME .GT. 1.0) THEN
	   SSQ=1E6
	   RETURN
	ENDIF
C
C If esds ore OK, calculate sum of residuals squared
C
	SSQ=0.0
	DO I=1,3*NSPOTS
	  SSQ=SSQ+(CALC(I)-OBS(I))**2
	ENDDO
C
	RETURN
	END


C----------------- The function called by LSQ refinement ---------------

	SUBROUTINE LSQ_FUNC(CALC,NOBS, PARS,NPARS)
C
	REAL CALC(NOBS),PARS(NPARS)
C
	COMMON /FIT_EFH_DATA_COM/ EFH(3,1000),XYPIX(2,1000),PIX_CEN(2),PIX_SIZE(2),DRUM_RAD
C
	LOGICAL LPARS
	COMMON /FIT_EFH_COM/ EIJ(6),ETA(3),ETA_ISO,LPARS(50)
C
C Unpack new PARS() values to /FIT_EFH_COM/
C NB: PARS() in the second argument is a dummy instead of ESDS()
C
	CALL UNPACK_LSQ_PARS(PARS,PARS,.FALSE.)
C
C Calculate ellipse shape and load E,F,H into CALC()
C
	DO I=1,NOBS/3
	  CALL CALC_SPOT_EFH(E,F,H, XYPIX(1,I), PIX_CEN,PIX_SIZE,DRUM_RAD, EIJ,ETA,ETA_ISO)
C Do not refine directly on E,F,H (must match LOAD_LSQ_DATA)
	  CALC(3*I-2)=1.0/SQRT(MAX(1E-5,E))
	  CALC(3*I-1)=1.0/SQRT(MAX(1E-5,F))
	  CALC(3*I  )=H/SQRT(MAX(1E-10,E*F))
CC	  CALC(3*I-2)=E
CC	  CALC(3*I-1)=F
CC	  CALC(3*I  )=H
	ENDDO
C
	RETURN
	END


C-------------------- Load data used for LSQ refinement ------------------

	SUBROUTINE LOAD_LSQ_DATA(OBS,SIG,NSPOTS)
C
	PARAMETER (NDATA_MAX=3000)
	REAL OBS(NDATA_MAX),SIG(NDATA_MAX)
C
	COMMON /FIT_EFH_DATA_COM/ EFH(3,1000),XYPIX(2,1000),PIX_CEN(2),PIX_SIZE(2),DRUM_RAD
C
C Load OBS() with observed E,F,H (or related values)
C
	DO I=1,NSPOTS
C Do not refine directly on E,F,H (must match LSQ_FUNC)
	  OBS(3*I-2)=1.0/SQRT(EFH(1,I))
	  OBS(3*I-1)=1.0/SQRT(EFH(2,I))
	  OBS(3*I  )=EFH(3,I)/SQRT(EFH(1,I)*EFH(2,I))
C	  OBS(3*I-2)=EFH(1,I)
C	  OBS(3*I-1)=EFH(2,I)
C	  OBS(3*I  )=EFH(3,I)
	ENDDO
C
C Load unit weights
C
	DO I=1,3*NSPOTS
	  SIG(I)=1.0
	ENDDO
C
	RETURN
	END


C----------------- Parameter specific routines ---------------

	SUBROUTINE CLEAR_LSQ_PARS
C
	LOGICAL LPARS
	COMMON /FIT_EFH_COM/ EIJ(6),ETA(3),ETA_ISO,LPARS(50)
C
C Default is 0.1mm radius crystal, and no anisotropic mosaicity
C
	DO I=1,3
	  EIJ(I)  =0.5
	  EIJ(I+3)=0.0
	  ETA(I)  =0.0
	ENDDO
C
C Default is 0.05 degree isotropic mosaicity
C
	ETA_ISO=0.05
C
	DO I=1,50
	  LPARS(I)=.FALSE.
	ENDDO
C
	RETURN
	END


	SUBROUTINE REFINE_LSQ_SHAPE
C
	LOGICAL LPARS
	COMMON /FIT_EFH_COM/ EIJ(6),ETA(3),ETA_ISO,LPARS(50)
C
	DO I=1,3
	  EIJ(I)  =0.1
	  EIJ(I+3)=0.0
	ENDDO
C
	DO I=1,6
	  LPARS(I)=.TRUE.
	ENDDO
C
	RETURN
	END


	SUBROUTINE REFINE_LSQ_ANISO
C
	LOGICAL LPARS
	COMMON /FIT_EFH_COM/ EIJ(6),ETA(3),ETA_ISO,LPARS(50)
C
	DO I=1,3
	  ETA(I)=0.05
	ENDDO
C
	DO I=7,9
	  LPARS(I)=.TRUE.
	ENDDO
C
	RETURN
	END


	SUBROUTINE REFINE_LSQ_ISO
C
	LOGICAL LPARS
	COMMON /FIT_EFH_COM/ EIJ(6),ETA(3),ETA_ISO,LPARS(50)
C
	ETA_ISO=0.01
	LPARS(10)=.TRUE.
C
	RETURN
	END


	SUBROUTINE PACK_LSQ_PARS(PARS,NPARS)
C
C Pack all refined LSQ parameters into PARS, also return
C the number of refined parameters in NPARS
C
	REAL PARS(*)
C
	LOGICAL LPARS
	COMMON /FIT_EFH_COM/ EIJ(6),ETA(3),ETA_ISO,LPARS(50)
C
	NPARS=0
C
	DO I=1,6
	  IF( LPARS(I) ) THEN
	    NPARS=NPARS+1
	    PARS(NPARS)=EIJ(I)
	  ENDIF
	ENDDO
C
	DO I=1,3
	  IF( LPARS(I+6) ) THEN
	    NPARS=NPARS+1
	    PARS(NPARS)=ETA(I)
	  ENDIF
	ENDDO
C
	IF( LPARS(10) ) THEN
	  NPARS=NPARS+1
	  PARS(NPARS)=ETA_ISO
	ENDIF
C
	RETURN
	END


	SUBROUTINE UNPACK_LSQ_PARS(PARS,ESDS,LPRINT)
C
C Unpack all refined LSQ parameters from PARS into /FIT_EFH_COM/
C If LPRINT is true, print out refined parameter values and esds.
C
	LOGICAL LPRINT
	REAL PARS(50),ESDS(50)
C
	LOGICAL LPARS
	COMMON /FIT_EFH_COM/ EIJ(6),ETA(3),ETA_ISO,LPARS(50)
C
C Unpack the refined parameters and output the refined result and esds
C
	NPARS=0
C
	DO I=1,6
	  IF( LPARS(I) ) THEN
	    NPARS=NPARS+1
	    EIJ(I)=PARS(NPARS)
	    IF( LPRINT ) WRITE(30,'(2X,A,I1,A,2F10.4)',IOSTAT=IDUM)
	1						'EIJ(',I,') =',PARS(NPARS),ESDS(NPARS)
	  ENDIF
	ENDDO
C
	DO I=1,3
	  IF( LPARS(I+6) ) THEN
	    NPARS=NPARS+1
	    ETA(I)=PARS(NPARS)
	    IF( LPRINT ) WRITE(30,'(2X,A,I1,A,2F10.4)',IOSTAT=IDUM)
	1						'ETA(',I,') =',PARS(NPARS),ESDS(NPARS)
	  ENDIF
	ENDDO
C
	IF( LPARS(10) ) THEN
	  NPARS=NPARS+1
	  ETA_ISO=PARS(NPARS)
	  IF( LPRINT ) WRITE(30,'(2X,A,2F10.4)',IOSTAT=IDUM)
	1						'ETA_ISO=',PARS(NPARS),ESDS(NPARS)
	ENDIF
C
	RETURN
	END


	SUBROUTINE GET_PAR_NAME(NAME, IPAR)
C
	CHARACTER NAME*(*)
C
	LOGICAL LPARS
	COMMON /FIT_EFH_COM/ EIJ(6),ETA(3),ETA_ISO,LPARS(50)
C
	NAME=REPEAT(' ',LEN(NAME))
	IF(LEN(NAME) .LT. 7) CALL QUIT('BUG(get_par_name): LEN(NAME) < 7')
C
	NPARS=0
C
	DO I=1,6
	  IF( LPARS(I) ) THEN
	    NPARS=NPARS+1
	    IF(IPAR .EQ. NPARS) WRITE(NAME,'(A,I3)') 'EIJ',I
	  ENDIF
	ENDDO
C
	DO I=1,3
	  IF( LPARS(6+I) ) THEN
	    NPARS=NPARS+1
	    IF(IPAR .EQ. NPARS) WRITE(NAME,'(A,I3)') 'ETA',I
	  ENDIF
	ENDDO
C
	IF( LPARS(10) ) THEN
	  NPARS=NPARS+1
	  IF(IPAR .EQ. NPARS) NAME='ETA_ISO'
	ENDIF
C
C Sanity check
C
	IF(NAME(1:1) .EQ. ' ') THEN
	  PRINT *,'IPAR=',IPAR
	  CALL QUIT('BUG(get_par_name): Invalid IPAR')
	ENDIF
C
	RETURN
	END


	SUBROUTINE GET_PAR_ESDS(ESDS, COV,NPARS)
C
	REAL ESDS(NPARS),COV(NPARS,NPARS)
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
	RETURN
	END


C-------------------- NLSCON SPECIFIC ROUTINES --------------------

	SUBROUTINE START_NLSCON(COVAR,ISIG, PAR,NPARS, IOBS,NOBS)
C
C NB: COVAR() is the returned covariance matrix
C NB: ISIG() is the uncertainty in IOBS() used in a weighted refinement
C
	REAL COVAR(NPARS*NPARS),ISIG(NOBS),PAR(NPARS),IOBS(NOBS)
C
C Size of work arrays, hopefully big enough!
C
	PARAMETER (IRW=170000)
	PARAMETER (IIW=3500)
	REAL RW(IRW)
	INTEGER IW(IIW)
C
	INTEGER IOPT(50)
C
C Two more arrays hopefully big enough, will check later.
C
	PARAMETER (NPARS_MAX=300)
	PARAMETER (NOBS_MAX=3000)
	REAL PAR_SCALE(NPARS_MAX)
C
	EXTERNAL NLSCON_FUNC,NLSCON_DUMMY_FUNC
C
C Complain and die if we haven't compiled in double-precision
C
	IF(1.0 .EQ. 1.0+1E-10) CALL QUIT('BUG: RECOMPILE USING DOUBLE PRECISION')
C
C Sanity check on array sizes.
C
	IF(NPARS .GT. NPARS_MAX) CALL QUIT('BUG(start_nlscon): NPARS > NPARS_MAX')
	IF(NOBS .GT. NOBS_MAX) CALL QUIT('BUG(start_nlscon): NOBS > NOBS_MAX')
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
      ITER=0
100	ITER=ITER+1
C Have editted NLSCON so RW(50...) returns esd's of parameters.
	CALL NLSCON(NPARS, NOBS,NOBS, NLSCON_FUNC,NLSCON_DUMMY_FUNC,
	1		PAR,PAR_SCALE, IOBS,ISIG, EPS, IOPT, IERR, IIW,IW,IRW,RW)
	IF (IERR.EQ.-1) GOTO 100
C
C Copy covariance matrix into COVAR().
C NB: The NLSCON documentation says RW(...) contains the
C     correlation matrix, this is incorrect.
C
	DO I=1,NPARS*NPARS
		COVAR(I)=RW(50+2*NPARS+I)
	ENDDO
C
C Only stop for software bugs (IERR>3)
C
	IF(IERR .EQ. 10) CALL QUIT('BUG(start_nlscon): IWK or RWK too small')
	IF(IERR .EQ. 20) CALL QUIT('BUG(start_nlscon): Invalid N,M,MFIT')
	IF(IERR .EQ. 21) CALL QUIT('BUG(start_nlscon): RTOL < 0')
	IF(IERR .EQ. 22) CALL QUIT('BUG(start_nlscon): XSCAL() < 0')
	IF(IERR .EQ. 30) CALL QUIT('BUG(start_nlscon): Invalid IOPT()')
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
      IW(31)=50
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
	CALL QUIT('BUG(nlscon_dummy_func): Function should never be called!')
      END


	SUBROUTINE NLSCON_FUNC(NPARS, NOBS,NCON, PAR,ICALC, IFAIL)
C
	REAL PAR(*),ICALC(*)
C
C A sanity check (there should be no constraints)
C
	IF(NCON .NE. 0) CALL QUIT('BUG(f): NCON .NE. 0')
C
	CALL LSQ_FUNC(ICALC,NOBS, PAR,NPARS)
C
	IFAIL=0
	RETURN
	END


C ================= R*8 versions of some mx3utils.for routines =============

	SUBROUTINE VC3_SET_R8(VC, V1,V2,V3)
C
	REAL VC(3)
C
	REAL TMP(3)
C
	TMP(1)=V1
	TMP(2)=V2
	TMP(3)=V3
C
	DO I=1,3
	  VC(I)=TMP(I)
	ENDDO
C
	RETURN
	END

	
	SUBROUTINE MX3_SET_R8(MX, V11,V12,V13, V21,V22,V23, V31,V32,V33)
C
	REAL MX(3,3)
C
	REAL TMP(3,3)
C
	TMP(1,1)=V11
	TMP(1,2)=	V12
	TMP(1,3)=		V13
	TMP(2,1)=V21
	TMP(2,2)=	V22
	TMP(2,3)=		V23
	TMP(3,1)=V31
	TMP(3,2)=	V32
	TMP(3,3)=		V33
C
	DO I1=1,3
	  DO I2=1,3
		MX(I1,I2)=TMP(I1,I2)
	  ENDDO
	ENDDO
C
	RETURN
	END

	
	SUBROUTINE MX3_XROT_R8(ROT, ANG)
C
C Create matrix ROT for a rotation of ANG (in degrees) around the X axis.
C                                    
      REAL ROT(3,3)
C
      COS_ANG=COSD(ANG)
      SIN_ANG=SIND(ANG)   
      ROT(1,1)=1.0
      ROT(1,2)=			0.0
      ROT(1,3)=						0.0
      ROT(2,1)=0.0
      ROT(2,2)=			COS_ANG
      ROT(2,3)=						-SIN_ANG
      ROT(3,1)=0.0
      ROT(3,2)=			SIN_ANG
      ROT(3,3)=						COS_ANG
C
      RETURN
      END


      SUBROUTINE MX3_YROT_R8(ROT, ANG)
C
C Create matrix ROT for a rotation of ANG (in degrees) around the Y axis.
C                                    
      REAL ROT(3,3)
C
      COS_ANG=COSD(ANG)
      SIN_ANG=SIND(ANG)   
      ROT(1,1)=COS_ANG
      ROT(1,2)=			0.0
      ROT(1,3)=						-SIN_ANG
      ROT(2,1)=0.0
      ROT(2,2)=			1.0
      ROT(2,3)=						0.0
      ROT(3,1)=SIN_ANG
      ROT(3,2)=			0.0
      ROT(3,3)=						COS_ANG
C
      RETURN
      END


	SUBROUTINE VC3_UNIT_R8(VC)
C
	REAL VC(3)
C
	SIZE=SQRT( VC(1)**2 + VC(2)**2 + VC(3)**2 )
	IF(SIZE .NE. 0) THEN
	  VC(1)=VC(1)/SIZE
	  VC(2)=VC(2)/SIZE
	  VC(3)=VC(3)/SIZE
	ENDIF
C
	RETURN
	END


	SUBROUTINE MX3_MULT_R8(M_OUT, M1,M2)
C
	REAL M_OUT(3,3),M1(3,3),M2(3,3)
C
	REAL SAVE(3,3)
C
	DO I=1,3
		DO J=1,3
			SAVE(I,J)=M1(I,1)*M2(1,J)+M1(I,2)*M2(2,J)+M1(I,3)*M2(3,J)
		ENDDO
	ENDDO
C     
	DO I=1,3
		DO J=1,3
			M_OUT(I,J)=SAVE(I,J)
		ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE MX3_TRANS_R8(A,B)
C
	REAL A(3,3),B(3,3)
C
	REAL SAVE(3,3)
C
	DO I1=1,3
	   DO I2=1,3
	      SAVE(I1,I2)=B(I1,I2)
	   ENDDO
	ENDDO
C
	DO I1=1,3
	   DO I2=1,3
	      A(I1,I2)=SAVE(I2,I1)
	   ENDDO
	ENDDO
C
	RETURN
	END


	FUNCTION VC3_DOT_R8(A,B)
C
	REAL A(3),B(3)
C
	VC3_DOT_R8=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
C
	RETURN
	END


C ============ Old EFH calculation routines ==================

	SUBROUTINE CALC_SPOT_EFH_OLD(E,F,H, XY_SPOT, XY_CEN,DRUM_PIX, EIJ,ETA,ETA_ISO)
C ===== DRUM_PIX is the radius in muliples of X pixel size
C
	REAL XY_SPOT(2),XY_CEN(2),EIJ(6),ETA(3)
C
	REAL EMAT(3,3),UMAT(3,3),TEMP(3,3),SPERP(3),S01(3)
C
C Setup diffraction geometry
C
C Calculate basic angles
	ALPHA_X=(XY_SPOT(1)-XY_CEN(1))/DRUM_PIX*57.296
	ALPHA_Y=ATAND( (XY_SPOT(2)-XY_CEN(2))/DRUM_PIX )
	TTH=ACOSD( COSD(ALPHA_Y) * COSD(ALPHA_X) )
	THETA=TTH/2
C Unit vector along S0 x S1
	SPERP(1)=-SIND(ALPHA_Y)/SIND(TTH)
	SPERP(2)=COSD(ALPHA_Y)*SIND(ALPHA_X)/SIND(TTH)
	SPERP(3)=0.0
C Unit vector along S1 - S0
	S01(1)=COSD(ALPHA_Y)*SIND(ALPHA_X)/(2*COSD(THETA))
	S01(2)=SIND(ALPHA_Y)/(2*COSD(THETA))
	S01(3)=1+COSD(ALPHA_Y)*COSD(ALPHA_X)/(2*COSD(THETA))
C Rotation matrix from lab to "spherical detector" coords
	UMAT(1,1)=COSD(ALPHA_X)
	UMAT(1,2)=				-SIND(ALPHA_Y)*SIND(ALPHA_X)
	UMAT(1,3)=								COSD(ALPHA_Y)*SIND(ALPHA_X)
	UMAT(2,1)=0.0
	UMAT(2,2)=				COSD(ALPHA_Y)
	UMAT(2,3)=								SIND(ALPHA_Y)
	UMAT(3,1)=-SIND(ALPHA_X)
	UMAT(3,2)=				-SIND(ALPHA_Y)*COSD(ALPHA_X)
	UMAT(3,3)=								COSD(ALPHA_Y)*COSD(ALPHA_X)
C
C Calculate "crystal shape" ellipse parameters (E1,F1,H1)
C
C Load 3D ellipsoidal matrix
C NB: EIJ contains E11,E22,E33,E12,E13,E23
C KLUDGE 1: Effectively force Eii	> 0
	EMAT(1,1)=ABS(EIJ(1))
	EMAT(1,2)=		EIJ(4)
	EMAT(1,3)=				EIJ(5)
	EMAT(2,1)=EIJ(4)
	EMAT(2,2)=		ABS(EIJ(2))
	EMAT(2,3)=				EIJ(6)
	EMAT(3,1)=EIJ(5)
	EMAT(3,2)=		EIJ(6)
	EMAT(3,3)=				ABS(EIJ(3))
C Convert to spherical detector coords, EMAT=UMAT'*(EMAT*UMAT)
	DO I=1,3
	  DO J=1,3
	    TEMP(I,J)=0.0
	    DO K=1,3
	      TEMP(I,J)=TEMP(I,J)+EMAT(I,K)*UMAT(K,J)
	    ENDDO
	  ENDDO
	ENDDO
C
	DO I=1,3
	  DO J=1,3
	    EMAT(I,J)=0.0
	    DO K=1,3
	      EMAT(I,J)=EMAT(I,J)+UMAT(K,I)*TEMP(K,J)
	    ENDDO
	  ENDDO
	ENDDO
C Calculate 2D ellipse parameters from ellipsoidal matrix
	E1=EMAT(1,1)-EMAT(1,3)*EMAT(1,3)/EMAT(3,3)
	F1=EMAT(2,2)-EMAT(2,3)*EMAT(2,3)/EMAT(3,3)
	H1=EMAT(1,2)-EMAT(1,3)*EMAT(2,3)/EMAT(3,3)
C KLUDGE 2: Force ellipse to be non-negative
	E1=ABS(E1)
	F1=ABS(F1)
	IF(H1**2 .GT. E1*F1) H1=H1*SQRT(E1*F1/H1**2)
C Always add 0.5 to E1 & F1, which is equivalent to a 0.1 mm radius crystal
	E1=E1+0.5
	F1=F1+0.5
C
C Calculate "mosaicity" ellipse parameters (E2,F2,H2)
C
C Get components of anisotropic mosaicity along SPERP & S01
	ETA_PERP=ETA(1)*SPERP(1) + ETA(2)*SPERP(2) + ETA(3)*SPERP(3)
	ETA_01  =ETA(1)*S01(1) + ETA(2)*S01(2) + ETA(3)*S01(3)
C Add isotropic contribution to components and convert to radians
C Calculate ellipse parameters for "spherical detector" coords
	A=2.0*DRUM_PIX*(ABS(ETA_PERP)+ABS(ETA_ISO))/57.296
	B=2.0*DRUM_PIX*(SIND(THETA)*ABS(ETA_01)+ABS(ETA_ISO))/57.296
	A=MAX(0.01,A)
	B=MAX(0.01,B)
	SIN_ANGLE=COSD(ALPHA_X)*SIND(ALPHA_Y)/SIND(TTH)
	COS_ANGLE=SIND(ALPHA_X)/SIND(TTH)
	E2=(COS_ANGLE/A)**2 + (SIN_ANGLE/B)**2
	F2=(SIN_ANGLE/A)**2 + (COS_ANGLE/B)**2
	H2=(A**(-2) - B**(-2))*SIN_ANGLE*COS_ANGLE
C
C Convolute the E1,F1,H1 and E2,F2,H2 ellipses
C
	E=( E1*E2*(F1+F2) - E1*H2**2 - E2*H1**2 ) / ( (E1+E2)*(F1+F2) - (H1+H2)**2 )
	F=( F1*F2*(E1+E2) - F1*H2**2 - F2*H1**2 ) / ( (E1+E2)*(F1+F2) - (H1+H2)**2 )
	H=( H1*(E2*F2-H2**2) + H2*(E1*F1-H1**2) ) / ( (E1+E2)*(F1+F2) - (H1+H2)**2 )
C
C Convert final ellipse to "cylindrical detector" coords
C
	F=F*COSD(ALPHA_Y)**2
	H=H*COSD(ALPHA_Y)
C
	RETURN
	END