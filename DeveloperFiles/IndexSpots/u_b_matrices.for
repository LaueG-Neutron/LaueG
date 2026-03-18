	SUBROUTINE ALIGN_VECTORS(U_ROT,H_PEAK,H_HKL,NPEAKS, ISTATUS)
C
C Align two sets of NPEAKS vectors (H_PEAK & H_HKL) by rotating
C the second set. LSQ is used to create a linearised rotation
C matrix. It is assumed that all vectors are unit vectors.
C A true rotation matrix corresponding to the rotation of the
C hkl vectors is returned in U_ROT().
C ISTATUS=1 on succes or =0 on failure (due to coplanar triplet)
C
	REAL U_ROT(3,3),H_PEAK(3,NPEAKS),H_HKL(3,NPEAKS)
C
	PARAMETER (NPEAKS_MAX=100)
	REAL JAC(3,3*NPEAKS_MAX),HESS(3,3),AINV(3,3),DY(3),ROT(3,3),PARS(3)
C
	IF(NPEAKS .GT. NPEAKS_MAX) CALL QUIT('BUG(align_vectors) NPEAKS > 100')
C
C Set the status flag for possible failure
C
	ISTATUS=0
C
C Can't work if less than 3 peaks
C
	IF(NPEAKS .LT. 3) RETURN
C
C Initialise U_ROT() as a unit matrix
C
	DO I1=1,3
	  DO I2=1,3
	    U_ROT(I1,I2)=0.0
	  ENDDO
	  U_ROT(I1,I1)=1.0
	ENDDO
C
C Create the Jacobian
C
	NITER=0
1	DO I=1,NPEAKS
	  JAC(1,3*I-2)=0.0
	  JAC(2,3*I-2)=				-H_HKL(3,I)
	  JAC(3,3*I-2)=							-H_HKL(2,I)
	  JAC(1,3*I-1)=H_HKL(3,I)
	  JAC(2,3*I-1)=				0.0
	  JAC(3,3*I-1)=							H_HKL(1,I)
	  JAC(1,3*I  )=-H_HKL(2,I)
	  JAC(2,3*I  )=				H_HKL(1,I)
	  JAC(3,3*I  )=							0.0
	ENDDO
C
C HESS = JAC * JAC^T
C
	DO I1=1,3
	  DO I2=1,3
	    HESS(I1,I2)=0.0
	    DO I3=1,3*NPEAKS
	      HESS(I1,I2)=HESS(I1,I2)+JAC(I1,I3)*JAC(I2,I3)
	    ENDDO
	  ENDDO
	ENDDO
C
C DY = JAC * (H_PEAK-H_HKL)
C
	DO I1=1,3
	  DY(I1)=0.0
	  DO I2=1,NPEAKS
	    DO I3=1,3
	      DY(I1)=DY(I1)+JAC(I1,I3+3*I2-3)*(H_PEAK(I3,I2)-H_HKL(I3,I2))
	    ENDDO
	  ENDDO
	ENDDO
C
C Solve for PARS the linear equation HESS * PARS = DY
C
	CALL MX3_INVERT(AINV,DET, HESS)
	IF(ABS(DET) .LT. 1E-4) RETURN
	CALL VC3_MULT(PARS,AINV,DY)
C
C Create the linearised rotation matrix
C
	ROT(1,1)=1.0
	ROT(1,2)=			-PARS(3)
	ROT(1,3)=						-PARS(2)
	ROT(2,1)=PARS(3)
	ROT(2,2)=			1.0
	ROT(2,3)=						PARS(1)
	ROT(3,1)=PARS(2)
	ROT(3,2)=			-PARS(1)
	ROT(3,3)=						1.0
C
C Rotate the H_HKL vectors (make sure they are unit vectors)
C
	DO I=1,NPEAKS
	  CALL VC3_MULT(H_HKL(1,I), ROT,H_HKL(1,I))
	  CALL VC3_UNIT(H_HKL(1,I))
	ENDDO
C
C Pre-multiply U_ROT to simulate the effect on the vectors
C
	CALL MX3_MULT(U_ROT, ROT,U_ROT)
C
C If routine hasn't converged then iterate up to 5 times
C
	DPAR=SQRT( PARS(1)**2 + PARS(2)**2 + PARS(3)**2 )
	NITER=NITER+1
	IF(DPAR.GT.1E-4 .AND. NITER.LT.5) GOTO 1
C
C Massage the U_ROT matrix to make it orthogonal, and retain
C the status flag value.
C
	CALL MAKE_U_ORTHOGONAL(U_ROT,ISTATUS)
C
	RETURN
	END



	SUBROUTINE CALC_U_FROM_TRIPLET(U, IHKLS,IPEAKS, ISTATUS)
C
C Calculate the U matrix from 3 vectors calculated by B * (hkl)
C and 3 observed "scattering vectors that are indexed in H_HKLS()
C and H_PEAKS() by the arrays IHKLS() & IPEAKS(). The scattering
C vectors are actually unit vectors, which is essential for this
C routine to work.
C ISTATUS=1 if success, =0 if failure (due to coplanar triplet)
C
	INTEGER IHKLS(3),IPEAKS(3)
	REAL U(3,3)
C
	PARAMETER (NHKLS_MAX=1000)
	COMMON /HKL_VECTORS_COM/ NHKLS,H_HKLS(3,NHKLS_MAX)
C
	PARAMETER (NPEAKS_MAX=10000)
	REAL H_PEAKS
	COMMON /PEAKS_COM/ H_PEAKS(3,NPEAKS_MAX)
C
	REAL MAT_PEAKS(3,3),MAT_HKLS(3,3),MAT_INV(3,3)
C
C Set the status flag to possible failure
C
	ISTATUS=0
C
C Copy the HKL & PEAK "scattering vectors" into matrices
C
	DO I1=1,3
	  DO I2=1,3
	    MAT_HKLS(I1,I2)=H_HKLS(I1,IHKLS(I2))
	    MAT_PEAKS(I1,I2)=H_PEAKS(I1,IPEAKS(I2))
	  ENDDO
	ENDDO
C
C Solve for U the linear equation MAT_PEAKS = U * MAT_HKLS.
C This works because were are actually using unit vectors.
C If impossible return ISTATUS=0 (which shouldn't happen).
C
	CALL MX3_INVERT(MAT_INV,DET, MAT_HKLS)
	IF(ABS(DET) .LT. 1E-3) RETURN
	CALL MX3_MULT(U, MAT_PEAKS,MAT_INV)
C
C Massage the U matrix to make it orthogonal, and retain
C the status flag value in case it is a coplanar triplet.
C
	CALL MAKE_U_ORTHOGONAL(U,ISTATUS)
	IF(ISTATUS.EQ.0) CALL MX3_PRINT('MAT_PEAKS',MAT_PEAKS)
C
	RETURN
	END



	SUBROUTINE MAKE_U_ORTHOGONAL(U,ISTATUS)
C
C Massage U matrix until it is a perfect rotation matrix.
C ISTATUS=1 if success, =0 if failed (due to coplanar vectors)
C
	REAL U(3,3)
C
	PARAMETER (TINY=1.0E-7)
	REAL TEMP(3,3)
C
C Set status flag for possible success
C
	ISTATUS=1
C
C Do a maximum of 10 iterations of the algorithm.
C
	DO ITER=1,10
C
C Invert U() and transpose and put into TEMP().
C For a pure rotation matrix this should equal U().
C If we can't invert the matrix, just give up.
C
	  CALL MX3_INVERT(TEMP,DET, U)
	  IF(ABS(DET) .LT. 1.E-3) THEN
	    ISTATUS=0
	    RETURN
	  ENDIF
	  CALL MX3_TRANS(TEMP,TEMP)
C
C Sum the differences squared between elements in U() and TEMP().
C
	  SUM=0.0
	  DO I1=1,3
	    DO I2=1,3
	      SUM=SUM+(U(I1,I2)-TEMP(I1,I2))**2
	    ENDDO
	  ENDDO
	  SUM=SQRT( SUM/9.0 )
C
C
C Average U() & TEMP() and load into U()
C
	  CALL MX3_ADD(TEMP,TEMP,U)
	  CALL MX3_SCALE(U,TEMP,0.5)
C
C Finish algorithm if SUM is close to the numeric precision.
C
	  IF(SUM .LT. 10.0*TINY) RETURN
C
	ENDDO
C
	RETURN
	END



	SUBROUTINE CALC_B_MATRIX(B, DIRECT)
C
C Calculate the Busing & Levy B matrix from the direct
C space cell dimensions
C
	REAL B(3,3),DIRECT(6)
C
	REAL RECIP(6)
C
C First calculate reciprocal lattice constants.
C
	CALL CONV_DIRECT_RECIP(RECIP,DIRECT)
C
C Busing & Levy's equation for B().
C
	B(1,1)=RECIP(1)
	B(1,2)=		RECIP(2)*COSD(RECIP(6))
	B(1,3)=				RECIP(3)*COSD(RECIP(5))
	B(2,1)=0.0
	B(2,2)=		RECIP(2)*SIND(RECIP(6))
	B(2,3)=				-RECIP(3)*SIND(RECIP(5))*COSD(DIRECT(4))
	B(3,1)=0.0
	B(3,2)=		0.0
	B(3,3)=				1.0/DIRECT(3)
C
	RETURN
	END



	SUBROUTINE CONV_DIRECT_RECIP(RECIP, DIRECT)
C
C Routine accepts a vector of the direct lattice constants
C in DIRECT() and returns the recip. lattice constants
C in RECIP().
C NB:  This routine can also convert from recip. to direct
C      lattice constants.
C
	REAL DIRECT(6),RECIP(6)
C
	A=DIRECT(1)
	B=DIRECT(2)
	C=DIRECT(3)
	ALPHA=DIRECT(4)
	BETA =DIRECT(5)
	GAMMA=DIRECT(6)
C
C Using formulae in Int. Tables. II, p106
C
	S=(ALPHA+BETA+GAMMA)/2.0
	TEMP=SIND(S)*SIND(S-ALPHA)*SIND(S-BETA)*SIND(S-GAMMA)
	V=2.0*A*B*C*SQRT( MAX(0.0,TEMP) )
C
	RECIP(1)=B*C*SIND(ALPHA)/V
	RECIP(2)=A*C*SIND(BETA )/V
	RECIP(3)=A*B*SIND(GAMMA)/V
C
	TEMP=(COSD(BETA)*COSD(GAMMA)-COSD(ALPHA))/(SIND(BETA)*SIND(GAMMA))
	RECIP(4)=ACOSD( MIN(1.0,TEMP) )
	TEMP=(COSD(ALPHA)*COSD(GAMMA)-COSD(BETA))/(SIND(ALPHA)*SIND(GAMMA))
	RECIP(5)=ACOSD( MIN(1.0,TEMP) )
	TEMP=(COSD(ALPHA)*COSD(BETA)-COSD(GAMMA))/(SIND(ALPHA)*SIND(BETA))
	RECIP(6)=ACOSD( MIN(1.0,TEMP) )
C
	RETURN
	END



	SUBROUTINE CALC_FROM_LAUEGEN_ANGLES(U, CELL,PHI_X,PHI_Y,PHI_Z)
C
C Given the CELL parameters and the corresponding PHI_* angles for Lauegen,
C calculate	our U matrix. We need the cell parameters as Lauegen uses
C the Buerger formulae for B instead of Busing & Levy. It only matters
C when some angles are not 90 degrees.
C
	REAL U(3,3),CELL(6)
C
	REAL B_BUSING(3,3),B_LAUEGEN(3,3),B_INV(3,3),U_LAUEGEN(3,3),TEMP(3,3)
C
	CALL CALC_B_MATRIX(B_BUSING, CELL)
C
	CALL CALC_B_MATRIX_LAUEGEN(B_LAUEGEN, CELL)
C
	CALL MX3_INVERT(B_INV,DET, B_BUSING)
	CALL MX3_MULT(TEMP, B_LAUEGEN,B_INV)
C
	CALL LAUEGEN_TO_U_MX(U_LAUEGEN, PHI_X,PHI_Y,PHI_Z)
	CALL MX3_MULT(U, U_LAUEGEN,TEMP)
C
	RETURN
	END



	SUBROUTINE CALC_LAUEGEN_ANGLES(PHI_X,PHI_Y,PHI_Z, U,CELL)
C
C Given our U matrix and the CELL parameters calculate the corresponding
C PHI_* angles for Lauegen. We need the cell parameters as Lauegen uses
C the Buerger formulae for B instead of Busing & Levy. It only matters
C when some angles are not 90 degrees.
C
	REAL U(3,3),CELL(6)
C
	REAL B_BUSING(3,3),B_LAUEGEN(3,3),B_INV(3,3),U_LAUEGEN(3,3),TEMP(3,3)
C
	CALL CALC_B_MATRIX(B_BUSING, CELL)
C
	CALL CALC_B_MATRIX_LAUEGEN(B_LAUEGEN, CELL)
C
	CALL MX3_INVERT(B_INV,DET, B_LAUEGEN)
	CALL MX3_MULT(TEMP, B_BUSING,B_INV)
C
	CALL MX3_MULT(U_LAUEGEN, U,TEMP)
	CALL U_MX_TO_LAUEGEN(PHI_X,PHI_Y,PHI_Z, U_LAUEGEN)
C
	RETURN
	END



	SUBROUTINE CALC_B_MATRIX_LAUEGEN(B, CELL)
C
	REAL B(3,3),CELL(6)
C
	REAL RCELL(6)
C
C First calculate reciprocal lattice constants.
C
	CALL CONV_DIRECT_RECIP(RCELL,CELL)
C
C====== Buerger's Xray Crystallography, p348
C
	COS_PSI=( COSD(RCELL(6)) - COSD(RCELL(5)) * COSD(RCELL(4)) )/SIND(RCELL(5))
	COS_RHO=SQRT( 1.0 - COSD(RCELL(6))**2 - COSD(RCELL(4))**2 - COSD(RCELL(5))**2
     *		+2.0*COSD(RCELL(6))*COSD(RCELL(4))*COSD(RCELL(5)) ) / SIND(RCELL(5))
C
	B(1,1)=RCELL(1)*SIND(RCELL(5))
	B(1,2)=RCELL(2)*COS_PSI
	B(1,3)=0.0
	B(2,1)=0.0
	B(2,2)=RCELL(2)*COS_RHO
	B(2,3)=0.0
	B(3,1)=RCELL(1)*COSD(RCELL(5))
	B(3,2)=RCELL(2)*COSD(RCELL(4))
	B(3,3)=RCELL(3)
C
	RETURN
	END



	SUBROUTINE LAUEGEN_TO_U_MX(U, PHI_X,PHI_Y,PHI_Z)
C
C U = COORDS * RZ(PHI_Z) * RY(-PHI_Y) * RX(PHI_X)
C
	REAL U(3,3)
C
	REAL COORDS(3,3)
C
	DATA COORDS /0.0,0.0,1.0, 1.0,0.0,0.0, 0.0,1.0,0.0/ 
C
	CALL CALC_RZYX(U, PHI_X,-PHI_Y,PHI_Z)
	CALL MX3_MULT(U, COORDS,U)
C
	RETURN
	END



	SUBROUTINE U_MX_TO_LAUEGEN(PHI_X,PHI_Y,PHI_Z, U)
C
C U = COORDS * RZ(PHI_Z) * RY(-PHI_Y) * RX(PHI_X)
C
	REAL U(3,3)
C
	REAL COORDS_INV(3,3),RZYX(3,3)
C
	DATA COORDS_INV /0.0,1.0,0.0, 0.0,0.0,1.0, 1.0,0.0,0.0/ 
C
	CALL MX3_MULT(RZYX, COORDS_INV,U)
	CALL CALC_RZYX_ANGLES(PHI_X,PHI_Y,PHI_Z, RZYX)
	PHI_Y=-PHI_Y
C
	RETURN
	END



	SUBROUTINE CALC_RXYZ_ANGLES(XROT,YROT,ZROT, RXYZ)
C
C Given a rotation matrix, RXYZ, made by RX * RY *RZ this
C routine calculates the X,Y,Z rotation angles
C
	REAL RXYZ(3,3)
C
	IF(RXYZ(1,3) .GT. 0.99999) THEN
	  XROT=0.0
	  ZROT=ATAN2D(RXYZ(2,1)+RXYZ(3,2),RXYZ(2,2)-RXYZ(3,1))
	ELSE IF(RXYZ(1,3) .LT. -0.99999) THEN
	  XROT=0.0
	  ZROT=ATAN2D(RXYZ(2,1)-RXYZ(3,2),RXYZ(2,2)+RXYZ(3,1))
	ELSE
	  XROT=ATAN2D(-RXYZ(2,3),RXYZ(3,3))
	  ZROT=ATAN2D(-RXYZ(1,2),RXYZ(1,1))
	ENDIF
	YROT=ASIND(-RXYZ(1,3))
	IF(SIND(XROT)*RXYZ(2,3) .GT. 0.0) YROT=180.0-YROT
C
	RETURN
	END



	SUBROUTINE CALC_RZYX_ANGLES(XROT,YROT,ZROT, RZYX)
C
C Given a rotation matrix, RZYX, made by RZ * RY *RX this
C routine calculates the X,Y,Z rotation angles
C
	REAL RZYX(3,3)
C
	IF(RZYX(3,1) .GT. 0.99999) THEN
	  XROT=0.0
	  ZROT=ATAN2D(-RZYX(1,2)-RZYX(2,3),RZYX(2,2)-RZYX(1,3))
	ELSE IF(RZYX(3,1) .LT. -0.99999) THEN
	  XROT=0.0
	  ZROT=ATAN2D(-RZYX(1,2)+RZYX(2,3),RZYX(2,2)+RZYX(1,3))
	ELSE
	  XROT=ATAN2D(RZYX(3,2),RZYX(3,3))
	  ZROT=ATAN2D(RZYX(2,1),RZYX(1,1))
	ENDIF
	YROT=ASIND(RZYX(3,1))
	IF(SIND(XROT)*RZYX(3,2) .LT. 0.0) YROT=180.0-YROT
C
	RETURN
	END



	SUBROUTINE CALC_RXYZ(RXYZ, XROT,YROT,ZROT)
C
	REAL RXYZ(3,3)
C
	REAL RX(3,3),RY(3,3),RZ(3,3)
C
	CALL CALC_XROT_MX(RX,XROT)
	CALL CALC_YROT_MX(RY,YROT)
	CALL CALC_ZROT_MX(RZ,ZROT)
C
	CALL MX3_MULT(RXYZ, RY,RZ)
	CALL MX3_MULT(RXYZ, RX,RXYZ)
C
	RETURN
	END



	SUBROUTINE CALC_RZYX(RZYX, XROT,YROT,ZROT)
C
	REAL RZYX(3,3)
C
	REAL RX(3,3),RY(3,3),RZ(3,3)
C
	CALL CALC_XROT_MX(RX,XROT)
	CALL CALC_YROT_MX(RY,YROT)
	CALL CALC_ZROT_MX(RZ,ZROT)
C
	CALL MX3_MULT(RZYX, RY,RX)
	CALL MX3_MULT(RZYX, RZ,RZYX)
C
	RETURN
	END



      SUBROUTINE CALC_XROT_MX(MX, ROT)
C
C Create matrix MX for a rotation of ROT around the X axis.
C                                    
      REAL MX(3,3)
C
      COS_ROT=COSD(ROT)
      SIN_ROT=SIND(ROT)   
      MX(1,1)=1.0
      MX(1,2)=			0.0
      MX(1,3)=						0.0
      MX(2,1)=0.0
      MX(2,2)=			COS_ROT
      MX(2,3)=						-SIN_ROT
      MX(3,1)=0.0
      MX(3,2)=			SIN_ROT
      MX(3,3)=						COS_ROT
C
      RETURN
      END



      SUBROUTINE CALC_YROT_MX(MX, ROT)
C
C Create matrix MX for a rotation of ROT around the Y axis.
C                                    
      REAL MX(3,3)
C
      COS_ROT=COSD(ROT)
      SIN_ROT=SIND(ROT)   
      MX(1,1)=COS_ROT
      MX(1,2)=			0.0
      MX(1,3)=						-SIN_ROT
      MX(2,1)=0.0
      MX(2,2)=			1.0
      MX(2,3)=						0.0
      MX(3,1)=SIN_ROT
      MX(3,2)=			0.0
      MX(3,3)=						COS_ROT
C
      RETURN
      END



      SUBROUTINE CALC_ZROT_MX(MX, ROT)
C
C Create matrix MX for a rotation of ROT around the Z axis.
C                                    
      REAL MX(3,3)
C
      COS_ROT=COSD(ROT)
      SIN_ROT=SIND(ROT)   
      MX(1,1)=COS_ROT
      MX(1,2)=			-SIN_ROT
      MX(1,3)=						0.0
      MX(2,1)=SIN_ROT
      MX(2,2)=			COS_ROT
      MX(2,3)=						0.0
      MX(3,1)=0.0
      MX(3,2)=			0.0
      MX(3,3)=						1.0
C
      RETURN
      END
