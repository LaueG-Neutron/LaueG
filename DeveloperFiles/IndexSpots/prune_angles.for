C ============ Routines specific to SOLN_PRUNE_ANGLES =============

	SUBROUTINE SOLN_PRUNE_ANGLES(DXY_SOLN, U_SOLN,NSOLN,B, NPEAKS)
C
C Loop through solutions in U_SOLN, searching for the best agreement
C for 20% of the peaks. The search includes X & Y offsets up to
C 12 pixels from the nominal center. The optimal XY offset and U
C matrix found for each search is returned in DXY_SOLN & U_SOLN.
C Any solution that cannot fit 20% of the peaks to an average of
C 0.25 degrees or better is removed.
C
	REAL B(3,3)
	PARAMETER (NSOLN_MAX=3000)
	REAL DXY_SOLN(2,NSOLN_MAX),U_SOLN(3,3,NSOLN_MAX)
C
	REAL H_PEAKS
	COMMON /PEAKS_COM/ H_PEAKS(3,10000)
C
	COMMON /HKL_VECTORS_COM/ NHKLS,H_HKLS(3,1000)
C
	COMMON /PIXEL_CEN_COM/ XCEN,YCEN
C
C
	REAL DANG_SOLN(NSOLN_MAX)
C
	REAL U(3,3),U_BEST(3,3)
C	
	FRAC=0.2
C
	DO I=1,NSOLN
	  DXY_SOLN(1,NSOLN)=0.0
	  DXY_SOLN(2,NSOLN)=0.0
	ENDDO
C
	DO I=1,NSOLN
	  DANG_SOLN(I)=1E6
	  RMULT=4.0
	  DO ITER=1,4
	    IF(ITER.GT.1 .AND. DANG_SOLN(I).GT.0.5) CYCLE
	    DX0=DXY_SOLN(1,I)
	    DY0=DXY_SOLN(2,I)
C
	    DO IDX=-2,2
	      DO IDY=-2,2
	        CALL MX3_COPY(U, U_SOLN(1,1,I))
	        DX2=DX0+RMULT*IDX
	        DY2=DY0+RMULT*IDY
	        CALL CHECK_ROT_XY_FRAC(DANG, U,B,FRAC, DX2,DY2, NPEAKS)
	        IF(DANG .LT. DANG_SOLN(I)) THEN
	          DANG_SOLN(I)=DANG
	          DXY_SOLN(1,I)=DX2
	          DXY_SOLN(2,I)=DY2
	          CALL MX3_COPY(U_BEST,U_SOLN(1,1,I))
	        ENDIF
	      ENDDO
	    ENDDO
C
	    CALL MX3_COPY(U_SOLN(1,1,I),U_BEST)
C
	    RMULT=RMULT*0.5
	  ENDDO
C
	ENDDO
C
C
C
	NSOLN2=0
	DO I=1,NSOLN
	  IF(DANG_SOLN(I) .LE. 0.25) THEN
	    NSOLN2=NSOLN2+1
	    CALL MX3_COPY(U_SOLN(1,1,NSOLN2), U_SOLN(1,1,I))
	    DXY_SOLN(1,NSOLN2)=DXY_SOLN(1,I)
	    DXY_SOLN(2,NSOLN2)=DXY_SOLN(2,I)
	    DANG_SOLN(NSOLN2)=DANG_SOLN(I)
	  ENDIF
	ENDDO
	NSOLN=NSOLN2
C
	PRINT '(1X,I5,A)',NSOLN,' solutions with best angles'
C
	RETURN
	END






	SUBROUTINE CHECK_ROT_XY_FRAC(DANG, U,B,FRAC, DX,DY, NPEAKS)
C
	REAL U(3,3),B(3,3)
C
	EXTERNAL DIFF_ANGLE_COMPARE
C
	PARAMETER (NTRIP_MAX=100000)
	COMMON /DIFF_ANGLE_COM/ DIFF_ANG(NTRIP_MAX)
C
	PARAMETER (NPEAKS_MAX=100)
	INTEGER ITAGS(NPEAKS_MAX)
C
	REAL H_PEAKS
	COMMON /PEAKS_COM/ H_PEAKS(3,10000)
C
	REAL HKLS(3,NPEAKS_MAX)
	REAL H_HKLS(3,NPEAKS_MAX),H_HKLS2(3,NPEAKS_MAX)
	REAL H_PEAKS1(3,NPEAKS_MAX),H_PEAKS2(3,NPEAKS_MAX)
C
	REAL U_ROT(3,3),UB(3,3),VEC(3)
C
C Calculate the UB matrix
C
	CALL MX3_MULT(UB, U,B)
C
C Calculate H_PEAKS() with X & Y shifted, store in H_PEAKS1
C
	DO I=1,NPEAKS
	  CALL CALC_HVEC_TO_PIXEL(X,Y, H_PEAKS(1,I))
	  CALL CALC_PIXEL_TO_HVEC(H_PEAKS1(1,I), X-DX,Y-DY)
	ENDDO
C
C Calculate the hkl list for the X & Y shifted peaks
C
	CALL CALC_BEST_HKLS(HKLS, UB,H_PEAKS1,NPEAKS)
C
C Calculate (unit) scattering vectors for HKLS
C
	DO I=1,NPEAKS
	  CALL VC3_MULT(H_HKLS(1,I), UB,HKLS(1,I))
	  CALL VC3_UNIT(H_HKLS(1,I))
	ENDDO
C
C Calculate angles between H_HKLS & H_PEAKS1 and put in DIFF_ANG()
C
	DO I=1,NPEAKS
	  CALL VC3_SUB(VEC,H_PEAKS1(1,I),H_HKLS(1,I))
	  DIFF_ANG(I)=ASIND( VC3_SIZE(VEC) )
	ENDDO
C
C Tagged sort on angle differences
C
	CALL TAG_SORT(DIFF_ANGLE_COMPARE,ITAGS,NPEAKS)
C
C Copy the best FRAC*NPEAKS peaks in new lists
C
	NPEAKS2=NINT(NPEAKS*FRAC)
	DO I=1,NPEAKS2
	  CALL VC3_COPY(H_HKLS2(1,I),H_HKLS(1,ITAGS(I)))
	  CALL VC3_COPY(H_PEAKS2(1,I),H_PEAKS1(1,ITAGS(I)))
	ENDDO
C
C Rotate H_HKLS2 vectors (and U) for the best match with H_PEAKS2
C If impossible due to coplanar triplet, return with failure
C
	CALL ALIGN_VECTORS(U_ROT,H_PEAKS2,H_HKLS2,NPEAKS2, ISTATUS)
	IF(ISTATUS .EQ. 0) RETURN
	CALL MX3_MULT(U, U_ROT,U)
C
C Calculate the average angle between H_PEAKS2 & H_HKLS2 vectors
C
	SUMSQ=0.0
	DO I=1,NPEAKS2
	  CALL VC3_SUB(VEC,H_PEAKS2(1,I),H_HKLS2(1,I))
	  DANG=ASIND( VC3_SIZE(VEC) )
	  SUMSQ=SUMSQ+DANG**2
	ENDDO
	DANG=SQRT(SUMSQ/NPEAKS2)
C
	RETURN
	END
