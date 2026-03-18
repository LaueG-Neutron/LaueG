	PROGRAM ORIENT_TWINS
C
	COMMON /U_B_INIT_COM/ U1_0(3,3),U2_0(3,3),B1_0(3,3),B2_0(3,3)
C
	REAL UB1(3,3),UB2(3,3),CELL1(6),CELL2(6)
	real mx1(3,3),mx2(3,3)
C
	LOGICAL LFITCELL
	CHARACTER FILE_NAME*132
C
C Output a simple banner
C
	PRINT '(1X,A)','Twin/pair spot orienting program for LaueG'//
	1			' (Ross Piltz, 2/10/2024)'
C
C Load parameters from SCILAB file
C
	CALL READ_PAR_FILE(UB1,UB2,MAX_HKL,WAV_MIN,DXY_MAX,NITER,LFITCELL,FILE_NAME)
C
C Load image from	TIF file and strip off the background
C
	CALL LOAD_STRIPPED_IMAGE(TRIM(FILE_NAME))
C
C Extract U, B & unit cell from initial UBs, and store initial U1 & U2 in /U_INIT_COM/
C
	CALL UNPACK_UB(U1_0,B1_0,CELL1, UB1)
	CALL UNPACK_UB(U2_0,B2_0,CELL2, UB2)
	call mx3_print('b1',b1_0)
	call mx3_print('b2',b2_0)
	call mx3_print('u1',u1_0)
	call mx3_print('u2',u2_0)
	call mx3_trans(mx1, u1_0)
	call mx3_mult(mx2, mx1,u2_0)
	call mx3_print('u1_T *u2',mx2)
C
C Output how UBs will transformed by iterations
C
	PRINT '(/,1X,A)','Iterating UBs to fit XY of spot pairs'
	IF( LFITCELL ) THEN
		PRINT '(1X,A)','Cell dimensions will be fitted as triclinic and independent'
		PRINT '(/,1X,A,2(/,1X,3F8.3,3F8.2,A))','Initial cells:',CELL1,' (1)',CELL2,' (2)'
	ELSE
		PRINT '(1X,A)','Cell dimensions will not be fitted'
	ENDIF
C
C Generate HKLs & calculated XYs for pairs of twin spots
C
	call mx3_print('ub1',ub1)
	call mx3_print('ub2',ub2)
	DO ITER=1,NITER
		PRINT '(/,1X,A,I2)','Iteration',ITER
C Run in INDEX mode for ITER = 1
	  IF(ITER .EQ. 1) THEN
			PRINT '(3X,A)','Iteration running in INDEX mode'
			CALL RUN_ITERATION(UB1,UB2, MAX_HKL,WAV_MIN,DXY_MAX, LFITCELL,.TRUE.)
		ELSE
			CALL RUN_ITERATION(UB1,UB2, MAX_HKL,WAV_MIN,DXY_MAX, LFITCELL,.FALSE.)
		ENDIF
	call mx3_print('ub1',ub1)
	call mx3_print('ub2',ub2)
	ENDDO
C
	CALL OUTPUT_RESULTS(UB1,UB2,LFITCELL)
C
CCCC	CALL DELETE_FILE('___laueg_orient_twins.in')
C
	PRINT *
	PRINT *,'SUCCESSFUL COMPLETION'
	END


	SUBROUTINE RUN_ITERATION(UB1,UB2, MAX_HKL,WAV_MIN,DXY_MAX, LFITCELL,LINDEX)
C
	LOGICAL LINDEX,LFITCELL
	REAL UB1(3,3),UB2(3,3)
C
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	REAL HKLLIST(16,1000),VEC(3),TMP(3,3)
	REAL U1_1(3,3),U1_2(3,3),U2_1(3,3),U2_2(3,3)
	REAL B1(3,3),B2(3,3),CELL1(6),CELL2(6)
C
C Generate HKLs & calculated XYs for pairs of twin spots
C
	CALL GEN_HKLLIST_PAIRS(HKLLIST,NHKLS, UB1,UB2,WAV_MIN,MAX_HKL,DXY_MAX)
	NHKLS0=NHKLS
	PRINT '(I7,A)',NHKLS,' spot pairs generated'
C
C Fit observed XY for spot pairs and store in HKLLIST(11..14,*)
C Amplitudes of spot pairs stored in HKLLIST(15..16,*)
C Also remove pairs from HKLLIST() that may be singletons
C
	CALL FIT_ALL_PAIRS(HKLLIST,NHKLS,LINDEX)
	PRINT '(I7,A)',NHKLS,' spot pairs fitted'
C
C Output some statistics on the fit, and prune outliers if LINDEX is false
C
	CALL CALC_STATS_PRUNE(HKLLIST,NHKLS, (.NOT.LINDEX) )
C
C Calc U matrices from UB matrices before changing
C
	CALL UNPACK_UB(U1_1,B1,CELL1, UB1)
	CALL UNPACK_UB(U2_1,B2,CELL2, UB2)
C
C Transform UB1 & UB2 matrices to fit observed XY of spot pairs
C
	IF( LINDEX ) THEN
		CALL TRANSFORM_UB_MX(UB1,UB2, HKLLIST,NHKLS, PHI,.FALSE.)
	ELSE
		CALL TRANSFORM_UB_MX(UB1,UB2, HKLLIST,NHKLS, PHI,LFITCELL)
	ENDIF
C
C Calc U matrices from UB matrices after changing
C
	CALL UNPACK_UB(U1_2,B1,CELL1, UB1)
	CALL UNPACK_UB(U2_2,B2,CELL2, UB2)
C
C Calculate angles that U matrices were rotated
C
	CALL MX3_DIVIDE(TMP, U1_2,U1_1)
	CALL CALC_MATRIX_ROTATION(VEC,ANGLE1, TMP)
	CALL MX3_DIVIDE(TMP, U2_2,U2_1)
	CALL CALC_MATRIX_ROTATION(VEC,ANGLE2, TMP)
	PRINT '(3X,A,2(F7.3,A))','UB1 & UB2 rotated by',ANGLE1,' &',ANGLE2,' degrees'
C
C If refined, output new cell dimensions
C
	IF( LFITCELL .AND. (.NOT.LINDEX) ) THEN
		PRINT '(3X,A,2(/,3X,3F8.3,3F8.2,A))','Refined cells:',CELL1,' (1)',CELL2,' (2)'
	ENDIF
C
	RETURN
      END


	SUBROUTINE CALC_STATS_PRUNE(HKLLIST,NHKLS, LPRUNE)
C
C Output some statistics on the fit, and prune outliers if LPRUNE is true
C
	LOGICAL LPRUNE
	REAL HKLLIST(16,1000)
C
	COMMON /RMS_COM/ NHKLS0,RMS1,RMS2
C
C Do 2 iterations (at most) of stats and pruning
C
	DO ITER=1,2
C Calculate and output statistics on r.m.s. and worst XY differences
		RMS1=0.0
		RMS2=0.0
		BAD1=0.0
		BAD2=0.0
		DO I=1,NHKLS
			DSQ1=(HKLLIST(11,I)-HKLLIST(7,I))**2 + (HKLLIST(12,I)-HKLLIST( 8,I))**2
			DSQ2=(HKLLIST(13,I)-HKLLIST(9,I))**2 + (HKLLIST(14,I)-HKLLIST(10,I))**2
			RMS1=RMS1+DSQ1
			RMS2=RMS2+DSQ2
			BAD1=MAX(BAD1,SQRT(DSQ1))
			BAD2=MAX(BAD2,SQRT(DSQ2))
		ENDDO
		RMS1=SQRT(RMS1/NHKLS)
		RMS2=SQRT(RMS2/NHKLS)
		PRINT '(7X,2(A,F4.1),A)','Twin 1: rms=',RMS1,', worst=',BAD1,' pixels'
		PRINT '(7X,2(A,F4.1),A)','Twin 2: rms=',RMS2,', worst=',BAD2,' pixels'
C Exit loop if not pruning outliers, or second iteration
		IF( .NOT.LPRUNE .OR. (ITER .EQ. 2) ) EXIT
C Remove any really bad outliers from HKLLIST
		NSKIP=0
		DO I=1,NHKLS
			DSQ1=(HKLLIST(11,I)-HKLLIST(7,I))**2 + (HKLLIST(12,I)-HKLLIST( 8,I))**2
			DSQ2=(HKLLIST(13,I)-HKLLIST(9,I))**2 + (HKLLIST(14,I)-HKLLIST(10,I))**2
			IF( (SQRT(DSQ1) .GT. 3.0*RMS1) .OR. (SQRT(DSQ2) .GT. 3.0*RMS2) ) THEN
				NSKIP=NSKIP+1
			ELSE
				DO I2=1,16
					HKLLIST(I2,I-NSKIP)=HKLLIST(I2,I)
				ENDDO
			ENDIF
		ENDDO
		PRINT '(I6,A)',NSKIP,' outliers removed'
C Exit loop if no outliers were removed
		IF(NSKIP .EQ. 0) EXIT
C Loop back to recalculate stats
		NHKLS=NHKLS-NSKIP
	ENDDO
C
C Save NHKLS (along with RMS1,2) in /RMS_COM/
C
	NHKLS0=NHKLS
C
	RETURN
	END


	SUBROUTINE OUTPUT_RESULTS(UB1,UB2,LFITCELL)
C
	LOGICAL LFITCELL
	REAL UB1(3,3),UB2(3,3)
C
	COMMON /U_B_INIT_COM/ U1_0(3,3),U2_0(3,3),B1_0(3,3),B2_0(3,3)
C
	COMMON /RMS_COM/ NHKLS0,RMS1,RMS2
C
	REAL U1(3,3),U2(3,3),B1(3,3),B2(3,3)
	REAL ROT1(3,3),ROT2(3,3),ROT21(3,3),TMP(3,3),VEC(3),CELL1(6),CELL2(6)
C
C Extract U, B & unit cell from final UBs
C
	CALL UNPACK_UB(U1,B1,CELL1, UB1)
	CALL UNPACK_UB(U2,B2,CELL2, UB2)
C
C Calculate rotation matrix of final Us compared to initial Us
C
	CALL MX3_DIVIDE(ROT1, U1,U1_0)
	CALL MX3_DIVIDE(ROT2, U2,U2_0)
C
C Output rotation axis in hkl for Twin 1, and angle of rotation for ROT1 & ROT2
C
	CALL CALC_MATRIX_ROTATION(VEC,ANGLE, ROT1)
	CALL MX3_INVERT(TMP,DET, UB1)
	CALL VC3_MULT(VEC, TMP,VEC)
	SIZE=MAX( ABS(VEC(1)), ABS(VEC(2)), ABS(VEC(3)) )
	PRINT '(/,1X,A,/,3X,F6.2,A,3F5.1)','UB1 total rotation:',ANGLE,
	1				' degrees around HKL(1) = ',VEC/SIZE*6.0
C
	CALL CALC_MATRIX_ROTATION(VEC,ANGLE, ROT2)
	CALL MX3_INVERT(TMP,DET, UB1)
	CALL VC3_MULT(VEC, TMP,VEC)
	SIZE=MAX( ABS(VEC(1)), ABS(VEC(2)), ABS(VEC(3)) )
	PRINT '(1X,A,/,3X,F6.2,A,3F5.1)','UB2 total rotation:',ANGLE,
	1				' degrees around HKL(1) = ',VEC/SIZE*6.0
C
C Calculate rotation matrix of final U2 compared to final U1
C
	CALL MX3_DIVIDE(ROT21, U2,U1)
C
C Get vector (in reciprocal space) and angle of rotation
C
	CALL CALC_MATRIX_ROTATION(VEC,ANGLE, ROT21)
	CALL MX3_INVERT(TMP,DET, UB1)
	CALL VC3_MULT(VEC, TMP,VEC)
	SIZE=MAX( ABS(VEC(1)), ABS(VEC(2)), ABS(VEC(3)) )
	PRINT '(1X,A,/,3X,F6.2,A,3F5.1)','Rotation of UB2 relative to UB1:',ANGLE,
	1				' degrees around HKL(1) = ',VEC/SIZE*6.0
C
C If refined, output final refined cell dimensions
C
	IF( LFITCELL ) PRINT '(/,1X,A,2(/,1X,3F8.3,3F8.2,A))',
	1									'Final cells:',CELL1,' (1)',CELL2,' (2)'
C
C Write new cells, UBs and some statistics to output file
C
	OPEN(UNIT=2,FILE='___laueg_orient_twins.out',STATUS='UNKNOWN')
C
	WRITE(2,'(6F12.5)') CELL1
	WRITE(2,'(6F12.5)') CELL2
C
	WRITE(2,'(3F12.7)') (UB1(1,K),K=1,3)
	WRITE(2,'(3F12.7)') (UB1(2,K),K=1,3)
	WRITE(2,'(3F12.7)') (UB1(3,K),K=1,3)
C
	WRITE(2,'(3F12.7)') (UB2(1,K),K=1,3)
	WRITE(2,'(3F12.7)') (UB2(2,K),K=1,3)
	WRITE(2,'(3F12.7)') (UB2(3,K),K=1,3)
C
	WRITE(2,'(I8,2F8.2)') NHKLS0,RMS1,RMS2
C
	CLOSE(UNIT=2)
C
	RETURN
	END
C===================================================================


C========= TRANSFORM UB MATRICES FROM CALC AND OBS SPOT XY =========
	SUBROUTINE TRANSFORM_UB_MX(UB1,UB2, HKLLIST,NHKLS, PHI,LFITCELL)
C
C Use observed and calculated spots XY to transform UB1 & UB2. If LFITCELL
C is false, the UBs are only rotated, no change to the cell dimensions.
C If LFITCELL is true, the UBs are rotated and the cell dimensions are
C fitted as if twins are independent triclinic lattices.
C
	LOGICAL LFITCELL
	REAL UB1(3,3),UB2(3,3),HKLLIST(16,1000)
C
	REAL ROT1(3,3),ROT2(3,3)
C
	COMMON /CELL_COM/ CELL(6),ILATT,ICEN
C
	REAL HVEC_1C(3,100),HVEC_2C(3,100),HVEC_1O(3,100),HVEC_2O(3,100)
	REAL TMP_1C_1C(3,3),TMP_1C_1O(3,3),TMP_1C_2O(3,3),PHI_ROT(3,3)
C
C Calculate unit scatt. vectors for calc & obs XY of all pairs
C
	DO I=1,NHKLS
		CALL PIX2HVEC(HVEC_1C(1,I), HKLLIST( 7,I) )
		CALL PIX2HVEC(HVEC_2C(1,I), HKLLIST( 9,I) )
		CALL PIX2HVEC(HVEC_1O(1,I), HKLLIST(11,I) )
		CALL PIX2HVEC(HVEC_2O(1,I), HKLLIST(13,I) )
	ENDDO
C
C Do SCILAB operations TMP_p_q = HVEC_p \ HVEC_q using
C A\B = inv(A'*A) * (A'*B), which is a LSQ solution to the problem
C
	DO I1=1,3
		DO I2=1,3
		  TMP_1C_1C(I1,I2)=0.0
		  TMP_1C_1O(I1,I2)=0.0
		  TMP_1C_2O(I1,I2)=0.0
			DO I3=1,NHKLS
	      TMP_1C_1C(I1,I2)=TMP_1C_1C(I1,I2)+HVEC_1C(I1,I3)*HVEC_1C(I2,I3)
				TMP_1C_1O(I1,I2)=TMP_1C_1O(I1,I2)+HVEC_1C(I1,I3)*HVEC_1O(I2,I3)
				TMP_1C_2O(I1,I2)=TMP_1C_2O(I1,I2)+HVEC_1C(I1,I3)*HVEC_2O(I2,I3)
			ENDDO
		ENDDO
	ENDDO
C
	CALL MX3_INVERT(TMP_1C_1C,DET, TMP_1C_1C)
	CALL MX3_MULT(TMP_1C_1O, TMP_1C_1C,TMP_1C_1O)
	CALL MX3_MULT(TMP_1C_2O, TMP_1C_1C,TMP_1C_2O)
C
C Correct for PHI rotation: ROT=PHI_ROT' * TMP' * PHI_ROT
C
	CALL CALC_YROT_MX(PHI_ROT, -PHI)
C
	CALL MX3_MULT(ROT1, TMP_1C_1O,PHI_ROT)
	CALL MX3_TRANS(ROT1, ROT1)
	CALL MX3_MULT(ROT1, ROT1,PHI_ROT)
C
	CALL MX3_MULT(ROT2, TMP_1C_2O,PHI_ROT)
	CALL MX3_TRANS(ROT2, ROT2)
	CALL MX3_MULT(ROT2, ROT2,PHI_ROT)
C
C If LFITCELL is false, iterate ROT1 & ROT2 into pure rotation matrices
C
	IF( .NOT.LFITCELL ) THEN
		DO ITER=1,3
			CALL CLEAN_ROTATION_MATRIX(ROT1, ROT1)
			CALL CLEAN_ROTATION_MATRIX(ROT2, ROT2)
		ENDDO
	ENDIF
C
C Calculate new UB1 & UB2
C NB: Both ROT1 & ROT2 are relative to the UB1 matrix
C
	CALL MX3_MULT(UB2, ROT2,UB1)
	CALL MX3_MULT(UB1, ROT1,UB1)
C
	RETURN
      END


	SUBROUTINE PIX2HVEC(HVEC, XYPIX )
C
	REAL HVEC(3),XYPIX(2)
C
	REAL PIXWAV(3)
C
	PIXWAV(1)=XYPIX(1)
	PIXWAV(2)=XYPIX(2)
	PIXWAV(3)=1.0
C
	CALL CALC_PIXWAVS_TO_HVECS(HVEC, PIXWAV,1)
	SIZE=SQRT( HVEC(1)**2 + HVEC(2)**2 + HVEC(3)**2 )
	HVEC(1)=HVEC(1)/SIZE
	HVEC(2)=HVEC(2)/SIZE
	HVEC(3)=HVEC(3)/SIZE
C
	RETURN
	END
C===================================================================


C============== SOME MATRIX MANIPULATION ROUTINES ==================
	SUBROUTINE UNPACK_UB(U,B,CELL, UB)
C
	REAL U(3,3),B(3,3),CELL(6),UB(3,3)
C
C Extract U, B & unit cell from UB
C
	CALL CALC_UB_TO_CELL(CELL, UB)
	CALL CALC_B_MATRIX(B, CELL)
	CALL MX3_DIVIDE(U, UB,B)
C
	RETURN
	END


	SUBROUTINE CLEAN_ROTATION_MATRIX(ROT2, ROT)
C
	REAL ROT2(3,3),ROT(3,3)
C
	REAL TMP(3,3)
C
	CALL MX3_INVERT(TMP,DET, ROT)
	CALL MX3_TRANS(TMP, TMP)
	CALL MX3_ADD(TMP, TMP,ROT)
	CALL MX3_SCALE(ROT2, TMP,0.5)
C
	RETURN
	END


	SUBROUTINE CALC_MATRIX_ROTATION(VEC,ANGLE, ROT)
C
	REAL VEC(3),ROT(3,3)
C
	REAL VPERP(3),VPERP2(3),VEC2(3)
C
C Get rotation axis vector
C
	VEC(1)=ROT(3,2)-ROT(2,3)
	VEC(2)=ROT(1,3)-ROT(3,1)
	VEC(3)=ROT(2,1)-ROT(1,2)
	CALL VC3_UNIT(VEC)
C
C Make largest indices positive
C
	IMAX=1
	IF(ABS(VEC(2)) .GT. ABS(VEC(IMAX))) IMAX=2
	IF(ABS(VEC(3)) .GT. ABS(VEC(IMAX))) IMAX=3
	IF(VEC(IMAX) .LT. 0.0) CALL VC3_SCALE(VEC, VEC,-1.0)
C
C Make perpendicular vector roughly along 001
C
	CALL VC3_INIT(VPERP, 0.0,0.0,-1.0)
	CALL VC3_CROSS(VPERP, VEC,VPERP)
	CALL VC3_CROSS(VPERP, VEC,VPERP)
	CALL VC3_UNIT(VPERP)
C
C Make second perpendicular vector
C
	CALL VC3_CROSS(VPERP2, VEC,VPERP)
C
C Calc rotation angle of VPERP in direction of VPERP2
C
	CALL VC3_MULT(VEC2, ROT,VPERP)
	ANGLE=ACOSD(VC3_DOT(VPERP2,VEC2))-90.0
C
	RETURN
	END
C===================================================================


C============== FIT PAIRS OF SPOTS TO INTENSITIES ==================
	SUBROUTINE FIT_ALL_PAIRS(HKLLIST,NHKLS,LINDEX)
C
	LOGICAL LINDEX
      REAL HKLLIST(16,1000)
C
	REAL CON(-30:30,-30:30)
	REAL AFAC(1000)
	REAL XY0(2),XY1(2),XY2(2),XYCEN(2)
C
C Loop through HKLLIST() fitting observed spot XY
C
	N=0
	DO I=1,NHKLS
C
C Load calculated XYs as starting position for Twin 1 spot
		XY1(1)=HKLLIST(7,I)
		XY1(2)=HKLLIST(8,I)
C
C Initial "index" mode: (only XY1 is valid)
		IF(LINDEX) THEN
C Calculate COG "center" of intensities around XY1
			CALL CALC_CEN_XY(XYCEN, XY1)
C Calculate XY2 assuming XYCEN = (XY1 + XY2)/2
			XY2(1)=2.0*XYCEN(1)-XY1(1)
			XY2(2)=2.0*XYCEN(2)-XY1(2)
C
C Normal mode: (both XY1 & XY2 are valid)
		ELSE
C Load calculated XYs as starting position for Twin 2 spot
			XY2(1)=HKLLIST(9,I)
			XY2(2)=HKLLIST(10,I)
C Move the XY1 & XY2 midpoint to match COG "center"
			XY0(1)=(XY1(1)+XY2(1))/2.0
			XY0(2)=(XY1(2)+XY2(2))/2.0
			CALL CALC_CEN_XY(XYCEN, XY0)
			XY1(1)=XY1(1) + XYCEN(1) - XY0(1)
			XY1(2)=XY1(2) + XYCEN(2) - XY0(2)
			XY2(1)=XY2(1) + XYCEN(1) - XY0(1)
			XY2(2)=XY2(2) + XYCEN(2) - XY0(2)
		ENDIF
C
C Round XY1 & XY2 midpoint to integers
		IX0=NINT(XYCEN(1))
		IY0=NINT(XYCEN(2))
C Convolute spot with ellipse for 61x61 pixels centred on IX0,IY0
		CALL CONVOL_SMOOTH(CON,IERR, IX0,IY0)
	  IF(IERR .NE. 0) CYCLE
C Fit XY of spots pairs, and amplitudes (height) of pair
	  CALL FIT_SPOT_PAIR(XY1,XY2,A1,A2,IERR, CON,IX0,IY0)
	  IF(IERR .NE. 0) CYCLE
C Add info to HKLLIST(), allowing for skipped  pairs
		N=N+1
		DO I2=1,10
			HKLLIST(I2,N)=HKLLIST(I2,I)
		ENDDO
		HKLLIST(11,N)=XY1(1)
		HKLLIST(12,N)=XY1(2)
		HKLLIST(13,N)=XY2(1)
		HKLLIST(14,N)=XY2(2)
	  HKLLIST(15,N)=A1
	  HKLLIST(16,N)=A2
C
	ENDDO
	NHKLS=N
C
C Remove pairs with outlier amplitude ratios
C
	CUTOFF=1.0
C Calculate log of amplitude ratios, and do sum of values
	SUM=0.0
	DO I=1,NHKLS
	  AFAC(I)=LOG(MAX(0.1,MIN(10.0, HKLLIST(15,I)/HKLLIST(16,I) )))
		SUM=SUM+AFAC(I)
	ENDDO
C Remove spots outside CUTOFF from HKLLIST()
	N=0
	DO I=1,NHKLS
		IF(ABS(AFAC(I)-SUM/NHKLS) .GE. CUTOFF)	CYCLE
		N=N+1
		DO I2=1,16
			HKLLIST(I2,N)=HKLLIST(I2,I)
		ENDDO
	ENDDO
	NHKLS=N
C
	RETURN
	END


	SUBROUTINE CALC_CEN_XY(XYCEN, XY)
C
C Move XY1 & XY2 so their average is the local intensity COG
C
	REAL XYCEN(2),XY(2)
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500,IUNIT=6)
	COMMON /RIMAGE_COM/ RIMAGE(MAX_NUMX,MAX_NUMY),NUMX,NUMY
C
	IX0=NINT(XY(1))
	IY0=NINT(XY(2))
C
C Calculate pair centre from moments around IX0,IY0
C
	S=0.0
	SX=0.0
	SY=0.0
	DO IX=MAX(1,IX0-30),MIN(NUMX,IX0+30)
	  DO IY=MAX(1,IY0-30),MIN(NUMY,IY0+30)
			RI=MAX(0.1, RIMAGE(IX,IY) )				! to prevent divide by zero
			S=S+RI
			SX=SX+IX*RI
			SY=SY+IY*RI
		ENDDO
	ENDDO
	XYCEN(1)=SX/S
	XYCEN(2)=SY/S
C
	RETURN
	END


	SUBROUTINE FIT_SPOT_PAIR(XY1,XY2,A1,A2,IERR, CON,IX0,IY0)
C
	REAL XY1(2),XY2(2),CON(-30:30,-30:30)
C
	REAL DXYLO(2),DXYHI(2),DXY1(2),DXY2(2)
C	
C Setup return values for failure
C
	A1=0.0
	A2=0.0
C
C Convolute spot with ellipse for 61x61 pixels centred on IX0,IY0
C	
C Get offsets of calculated spot positions relative to array centre
C
	X1_OFF=NINT(XY1(1))-IX0
	Y1_OFF=NINT(XY1(2))-IY0
	X2_OFF=NINT(XY2(1))-IX0
	Y2_OFF=NINT(XY2(2))-IY0
C
C Calculate the start and end of a line crossing XY1 & XY2
C
	DX12=X2_OFF-X1_OFF
	DY12=Y2_OFF-Y1_OFF
	SIZE=SQRT(DX12**2 + DY12**2)
	IF(SIZE .LT. 1.0) THEN
	  DX12=0.0
		DY12=10.0
	ELSEIF(SIZE .LT. 10.0) THEN
	  DX12=DX12*10.0/SIZE
	  DY12=DY12*10.0/SIZE
	ENDIF
	DXYLO(1)=(X1_OFF+X2_OFF)/2.0-DX12
	DXYLO(2)=(Y1_OFF+Y2_OFF)/2.0-DY12
	DXYHI(1)=(X1_OFF+X2_OFF)/2.0+DX12
	DXYHI(2)=(Y1_OFF+Y2_OFF)/2.0+DY12
C	
C Fit both spots fixed to a line from DXYLO to DXYHI
	CALL FIT_PAIR_LINE(DXY1,DXY2,A1,A2,IERR, CON,DXYLO,DXYHI)
	IF(IERR .NE. 0) RETURN
C	
C Fit both spots along the perpendicular to the separation
C
	CALL FIT_PAIR_PERP(DXY1,DXY2, CON,DXYLO,DXYHI)
C	
	XY1(1)=DXY1(1) + IX0
	XY1(2)=DXY1(2) + IY0
	XY2(1)=DXY2(1) + IX0
	XY2(2)=DXY2(2) + IY0
C
	RETURN
	END


	SUBROUTINE FIT_PAIR_LINE(DXY1,DXY2,A1,A2,IERR, CON,DXYLO,DXYHI)
C
	REAL DXY1(2),DXY2(2),CON(-30:30,-30:30),DXYLO(2),DXYHI(2)
C
	INTEGER IDXY(2,120)
	REAL VALS(120)
C
		IERR=-1
C
C Store in VALS() the CON() intensities along a path from XYLO to XYHI
C in steps of 1/2 pixel. The X,Y of the path are stored in IDXY().
C
	DIST=SQRT( (DXYLO(1)-DXYHI(1))**2 + (DXYLO(2)-DXYHI(2))**2 )
	DIST2=FLOAT(NINT(2.0*DIST))
	NVALS=0
	DO R=0.0,DIST2
	  IDX=NINT( DXYLO(1) + (DXYHI(1)-DXYLO(1))*R/DIST2 )
	  IDY=NINT( DXYLO(2) + (DXYHI(2)-DXYLO(2))*R/DIST2 )
	  IF(MAX(ABS(IDX),ABS(IDY)) .LE. 30) THEN
			NVALS=NVALS+1
			IF(NVALS .GT. 120) RETURN
	    VALS(NVALS)=CON(IDX,IDY)
			IDXY(1,NVALS)=IDX
			IDXY(2,NVALS)=IDY
		ENDIF
	ENDDO
C
C Find maxima either side of weighted mean
C
	SY=0.0
	SXY=0.0
	DO I=1,NVALS
	  SY=SY+VALS(I)
	  SXY=SXY+I*VALS(I)
	ENDDO
	IMID=NINT(SXY/SY)
	IF( (IMID .LE. 1) .OR. (IMID .GE. NVALS) ) RETURN
C
	ILO=1
	DO I=1,IMID
	  IF(VALS(I) .GT. VALS(ILO)) ILO=I
	ENDDO
C
	IHI=IMID+1
	DO I=IMID+1,NVALS
	  IF(VALS(I) .GT. VALS(IHI)) IHI=I
	ENDDO
	IMX1=ILO
	IMX2=IHI
C
C
	CALL FIT_2GAUSS(X1,X2,A1,A2,E1,E2, VALS,NVALS, IMX1,IMX2,IMID)
C
	IF( (MIN(X1,X2) .LT. 1.0) .OR. (MAX(X1,X2) .GT. FLOAT(NVALS)) ) RETURN
C
	DXY1(1)=IDXY(1,NINT(X1))
	DXY1(2)=IDXY(2,NINT(X1))
	DXY2(1)=IDXY(1,NINT(X2))
	DXY2(2)=IDXY(2,NINT(X2))
C
	IERR=0
	RETURN
	END


	SUBROUTINE FIT_2GAUSS(X1,X2,A1,A2,E1,E2, FOBS,NOBS, IMX1,IMX2,IMID)
C
C Fit two gaussians to FOBS(1..NOBS) intensities and return XY,
C amplitudes, and XY esd of pair in X1,X2,A1,A2,E1,E2
C
	REAL FOBS(NOBS)
C
	COMMON /LSQ_REF_COM/ NPARS
C
	COMMON /FIT_PEAKS_COM/ WAXAX(5)
C
	WAXAX(1)=8.0			! Fixed HWHM value
	WAXAX(2)=10.0
	WAXAX(3)=MIN(IMID-4,IMX1)
	WAXAX(4)=10.0
	WAXAX(5)=MAX(IMID+4,IMX2)
C
	CALL RUN_LSQ(X1,X2,E1,E2, FOBS,NOBS)
C
	A1=2.0**WAXAX(2)
	A2=2.0**WAXAX(4)
C
	RETURN
	END


	SUBROUTINE FIT_PAIR_PERP(DXY1,DXY2, CON,DXYLO,DXYHI)
C
C Calculate unit vector perpendicular to the line used in FIT_PAIR_LINE
C
	REAL DXY1(2),DXY2(2),CON(-30:30,-30:30),DXYLO(2),DXYHI(2)
C
	REAL DPERP(2),VAL1(-18:18),VAL2(-18:18),DXY(2)
C
	DPERP(1)=+(DXYHI(2)-DXYLO(2))
	DPERP(2)=-(DXYHI(1)-DXYLO(1))
	SIZE=SQRT( DPERP(1)**2 + DPERP(2)**2 )
	DPERP(1)=DPERP(1)/SIZE
	DPERP(2)=DPERP(2)/SIZE
C
C Store in val1 & val2[] the convolution values for spots 1 & 2 displaced along dperp by +/-9 pixels
C
	SUM1=0
	SUM2=0
	DO IOFF=-18,18
	  DXY(1)=DPERP(1)*(0.5*IOFF)
	  DXY(2)=DPERP(2)*(0.5*IOFF)
	  VAL1(IOFF)=CON( NINT(DXY1(1)+DXY(1)), NINT(DXY1(2)+DXY(2)) )
	  VAL2(IOFF)=CON( NINT(DXY2(1)+DXY(1)), NINT(DXY2(2)+DXY(2)) )
	  SUM1=SUM1+VAL1(IOFF)
	  SUM2=SUM2+VAL2(IOFF)
	ENDDO
C
C Calc centroids of val1 & val2[]
C
	CEN1=0.0
	CEN2=0.0
	DO IOFF=-18,18
		CEN1=CEN1+(0.5*IOFF)*VAL1(IOFF)/SUM1
		CEN2=CEN2+(0.5*IOFF)*VAL2(IOFF)/SUM2
	ENDDO
C
C Return XY positions of VAL1 & VAL2 centroids
C
	DXY1(1)=NINT( DXY1(1) + DPERP(1)*CEN1 )
	DXY1(2)=NINT( DXY1(2) + DPERP(2)*CEN1 )
	DXY2(1)=NINT( DXY2(1) + DPERP(1)*CEN2 )
	DXY2(2)=NINT( DXY2(2) + DPERP(2)*CEN2 )
C
	RETURN
	END
C===================================================================


C====================== CONVOLUTION STUFF ==========================
	SUBROUTINE CONVOL_SMOOTH(CON,IERR, IX0,IY0)
C
	REAL CON(-30:30,-30:30)
C
	PARAMETER (MAX_NUMX=8000,MAX_NUMY=2500,IUNIT=6)
	COMMON /RIMAGE_COM/ RIMAGE(MAX_NUMX,MAX_NUMY),NUMX,NUMY
C
C
	REAL Z(-30:30,-30:30),ZT(-30:30,-30:30),Z0(-30:30,-30:30),ELL(-30:30,-30:30)
C
	IERR=-1
C
C Return with error if within 30 pixels of image boundary
C
	IF(MIN(IX0,IY0,NUMX-IX0,NUMY-IY0) .LT. 30) RETURN
C
C Make copy of background stripped image around xy0()
C
	NHALF=30
	DO I1=-30,NHALF
		DO I2=-NHALF,NHALF
			Z(I1,I2)=RIMAGE(I1+IX0,I2+IY0)
			ZT(-I1,-I2)=Z(I1,I2)
		ENDDO
	ENDDO
C
C Make autocorrelation
C
	CALL CONVOL2D(Z0, Z,ZT,NHALF)
C
C Calc moments of auto-correlation cut at 80% peak-height
	RMIN0=+1E9
	RMAX0=-1E9
	DO I1=-NHALF,NHALF
		DO I2=-NHALF,NHALF
			RMIN0=MIN(RMIN0, Z0(I1,I2) )
			RMAX0=MAX(RMAX0, Z0(I1,I2) )
		ENDDO
	ENDDO
	CUTOFF=0.2*RMIN0+0.8*RMAX0
C
	S00=0.0
	S20=0.0
	S11=0.0
	S02=0.0
	DO I1=-NHALF,NHALF
		DO I2=-NHALF,NHALF
	    IF(Z0(I1,I2) .GT. CUTOFF ) THEN
	       S00=S00+1.0
	       S20=S20+I1**2
	       S11=S11+I1*I2
	       S02=S02+I2**2
	    ENDIF
	  ENDDO
	ENDDO
C
C If less than 15 pixels above 80%, return as failure
C
	IF(S00 .LT. 15.0) RETURN
	S20=S20/S00
	S11=S11/S00
	S02=S02/S00
C	
C Convert moments to ellipse parameters, and increase size by sqrt(2)
	DEN=4*(S20*S02 - S11**2)
	E= S20/DEN *0.5
	F= S02/DEN *0.5
	H=-S11/DEN *0.5
C Create ellipse using quadratic radial dependence
	SUM=0.0
	DO I1=-NHALF,NHALF
		DO I2=-NHALF,NHALF
	    Z0(I1,I2)=MAX(0.0, 1.0 - (E*I1**2 + F*I2**2 + H*I1*I2) )
	    SUM=SUM+Z0(I1,I2)
	  ENDDO
	ENDDO
	DO I1=-NHALF,NHALF
		DO I2=-NHALF,NHALF
	    ELL(I2,I1)=Z0(I1,I2)/SUM
	  ENDDO
	ENDDO
C	
C Convolute data with ellipse
	CALL CONVOL2D(CON, Z,ELL,NHALF)
C
	IERR=0
	RETURN
	END


	SUBROUTINE CONVOL2D(Z12, Z1,Z2,NHALF)
C
	REAL Z1(-30:30,-30:30),Z2(-30:30,-30:30),Z12(-30:30,-30:30)
C
	DO M1=-NHALF,NHALF
	  DO M2=-NHALF,NHALF
	    S=0.0
	    DO I1=MAX(-NHALF,-NHALF+M1),MIN(NHALF,NHALF+M1)
	      DO I2=MAX(-NHALF,-NHALF+M2),MIN(NHALF,NHALF+M2)
	        S=S + Z1(I1,I2) * Z2(M1-I1,M2-I2)
	      ENDDO
	    ENDDO
	    Z12(M1,M2)=S
	  ENDDO
	ENDDO
C
	RETURN
	END
C===================================================================


C=========== GENERATE HKLS FOR PAIRS OF ADJACENT SPOTS =============
	SUBROUTINE GEN_HKLLIST_PAIRS_OLD(HKLLIST,NHKLS, UB1,UB2,WAV_MIN,MAX_HKL,DXY_MAX)
C
      REAL HKLLIST(16,1000),UB1(3,3),UB2(3,3)
C
	INTEGER HKLS
	COMMON /HKLGEN_COM/ HKLS(3,100000),PIXWAVS(3,100000),IMULT(100000)
C
	REAL HKL2(3),PIXWAV2(3)
C
C Generate HKLs for Twin 1 using WAV_MIN and put results in /HKLGEN_COM/
C
	CALL GEN_HKLS(NHKLS, UB1,WAV_MIN,0.5*WAV_MIN)
C

C
	RETURN
	END


	SUBROUTINE GEN_HKLLIST_PAIRS_OLD(HKLLIST,NHKLS, UB1,UB2,WAV_MIN,MAX_HKL,DXY_MAX)
C
      REAL HKLLIST(16,1000),UB1(3,3),UB2(3,3)
C
	INTEGER HKLS
	COMMON /HKLGEN_COM/ HKLS(3,100000),PIXWAVS(3,100000),IMULT(100000)
C
	REAL HKL2(3),PIXWAV2(3)
C
C Generate HKLs for Twin 1 using WAV_MIN and put results in /HKLGEN_COM/
C
	CALL GEN_HKLS(NHKLS, UB1,WAV_MIN,0.5*WAV_MIN)
C
C Find nearby Twin 2 reflection, and load both reflections into HKLLIST()
C
	N=0
	DO I=1,NHKLS
C Prune Twin 1 HKL to MAX_HKL limit
	  IF(MAX( ABS(HKLS(1,I)), ABS(HKLS(2,I)), ABS(HKLS(3,I)) ) .GT. MAX_HKL) CYCLE
C Find best integer HKLs for Twin 2 (obeys MAX_HKL, but ignores WAV_MIN)
		CALL CALC_PIXWAVS_TO_HKLS(HKL2, UB2,PIXWAVS(1,I),1)
		SIZE=MAX( ABS(HKL2(1)), ABS(HKL2(2)), ABS(HKL2(3)) )
		SUM_MIN=1E9
		DO IMUL=1,MAX_HKL
			TMP1=HKL2(1)/SIZE*IMUL
			TMP2=HKL2(2)/SIZE*IMUL
			TMP3=HKL2(3)/SIZE*IMUL
			SUM=ABS(TMP1-NINT(TMP1)) + ABS(TMP2-NINT(TMP2)) + ABS(TMP3-NINT(TMP3))
			IF(SUM .LT. SUM_MIN) THEN
				SUM_MIN=SUM
				IBEST=IMUL
			ENDIF
		ENDDO
		HKL2(1)=NINT( HKL2(1)/SIZE*IBEST )
		HKL2(2)=NINT( HKL2(2)/SIZE*IBEST )
		HKL2(3)=NINT( HKL2(3)/SIZE*IBEST )
C Calculate X,Y for HKL2
		CALL CALC_HKLS_TO_PIXWAVS(PIXWAV2, UB2,HKL2,1)			
C Ignore if twin spots are more than DXY_MAX pixels apart
		DXY=SQRT( (PIXWAVS(1,I)-PIXWAV2(1))**2 + (PIXWAVS(2,I)-PIXWAV2(2))**2 )
		IF(DXY .GT. DXY_MAX) CYCLE
C Load HKLs & calculated XYs for twin pairs into HKLLIST()
		N=N+1
		DO I2=1,3
			HKLLIST(I2,N)=HKLS(I2,I)
			HKLLIST(I2+3,N)=HKL2(I2)
		ENDDO
		DO I2=1,2
			HKLLIST(I2+6,N)=PIXWAVS(I2,I)
			HKLLIST(I2+8,N)=PIXWAV2(I2)
		ENDDO
	ENDDO
	NHKLS=N
C
	RETURN
	END


	SUBROUTINE GEN_HKLS(NHKLS, UB,WAV_MIN,D_MIN)
C
C Copy to /HKLGEN_COM/ a list of unique hkls which strike the detector.
C Minimum d-spacing and wavelength is set by D_MIN and WAV_MIN.
C	INTEGER HKLS
C	COMMON /HKLGEN_COM/ HKLS(3,100000),PIXWAVS(3,100000),IMULT(100000)
C
	REAL UB(3,3)
C
	COMMON /CELL_COM/ CELL(6),ILATT,ICEN
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	CALL GEN_HEMISPHERE(NHKLS, UB,ICEN,D_MIN,1E6)
C	PRINT *,NHKLS,' hkls generated in a hemisphere'
	CALL PRUNE_DETECTED(NHKLS, UB,WAV_MIN,1E6, 1.0,FLOAT(NUMX),1.0,FLOAT(NUMY))
C	PRINT *,NHKLS,' hkls after pruning to detection limits'
	CALL PRUNE_MULTIPLES(NHKLS)
C	PRINT *,NHKLS,' hkls after pruning multiples'
C
	RETURN
	END
C===================================================================


C===================================================================
	SUBROUTINE DELETE_FILE(FILE_NAME)
C
C Delete the file if the "delete files" is missing or contains a TRUE
C
	CHARACTER FILE_NAME*(*)
C
	LOGICAL LDELETE
C
	OPEN(UNIT=17,FILE='___laueg_delete_files.in',STATUS='OLD',ERR=100)
	READ(17,*,ERR=100,END=100) LDELETE
	IF( .NOT.LDELETE ) RETURN
100	CLOSE(UNIT=17,IOSTAT=IDUMMY)
C
	OPEN(UNIT=17,FILE=FILE_NAME,STATUS='OLD',IOSTAT=IDUMMY)
	CLOSE(UNIT=17,STATUS='DELETE',IOSTAT=IDUMMY)
C
	RETURN
	END


	SUBROUTINE QUIT(TEXT)
C
C Workaround as SCILAB doesn't return STOP messages as they
C go to stderr not stdout
C
	CHARACTER*(*) TEXT
C
	PRINT *,TEXT
	CALL EXIT()
C
	END
