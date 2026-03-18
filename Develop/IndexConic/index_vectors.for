	SUBROUTINE INDEX_NODALS(V_NODAL1,V_NODAL2, HVECS,NVECS)
C
C Index the scattering vectors, HVECS, for spots on a conic
C between two nodal spots, V_NODAL1 & V_NODAL2.
C
	REAL V_NODAL1(3),V_NODAL2(3),HVECS(3,10000)
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS0
C
	REAL RH(100),RK(100)
C
C
	INTEGER IH(100),IK(100)
	REAL VEC1(3),PIXWAV(3)
C
C Calculate RH() & RK() where HVECS || RH * V_NODAL1 + RK * V_NODAL2
C
	DOT12=VC3_DOT(V_NODAL1,V_NODAL2)
	DO I=1,NVECS
	  DOT1=VC3_DOT(HVECS(1,I),V_NODAL1)
	  DOT2=VC3_DOT(HVECS(1,I),V_NODAL2)
	  RH(I)=DOT1 - DOT12*DOT2
	  RK(I)=DOT2 - DOT12*DOT1
	ENDDO
C
C
	ITEST=MIN(10, (NVECS-2)/2)
C
C Determine a likely ratio for the sizes of nodal vectors 1 & 2
C
	CALL GET_NODAL_RATIO(RATIO, RH,RK,NVECS,ITEST)
	CALL VC3_SCALE(V_NODAL2, V_NODAL2,RATIO)
C
	DO I=1,NVECS
	  RK(I)=RK(I)/RATIO
	ENDDO
C
C
	PRINT '(/,1X,A,/)','Initial indexing'
	PRINT '(1X,A,3F7.3,3X,3F7.3)','a*, b* =',V_NODAL1,V_NODAL2
	A_B=VC3_SIZE(V_NODAL1)/VC3_SIZE(V_NODAL2)
	GAMMA=ACOSD( VC3_DOT(V_NODAL1,V_NODAL2)/VC3_SIZE(V_NODAL1)/VC3_SIZE(V_NODAL2) )
	PRINT '(4X,A,F6.3,3X,A,F6.1)','|a*| / |b*| =',A_B,'angle =',GAMMA
C
C
C
	CALL CALC_INITIAL_HKS(IH,IK, RH,RK,NVECS)
C
	CALL PRUNE_BEST_RATIOS(IH(3),IK(3),RH(3),RK(3),NVECS,ITEST)
C
	PRINT *
C
	DO I=1,NVECS
C
	  IF(IH(I).NE.0 .OR. IK(I).NE.0) THEN
	    DO I2=1,3
	      VEC1(I2)=IH(I)*V_NODAL1(I2) + IK(I)*V_NODAL2(I2)
	    ENDDO
	    CALL VC3_UNIT(VEC1)
	    CALL CALC_HVECS_TO_PIXWAVS(PIXWAV, VEC1,1)
	    DIFF=SQRT( (SPOTS(1,I)-PIXWAV(1))**2 + (SPOTS(2,I)-PIXWAV(2))**2 )
	    PRINT '(I4,I3,3X,2F7.1,3X,3F7.1)',IH(I),IK(I),(SPOTS(K,I),K=1,2),PIXELS,DIFF
	  ENDIF
C
	ENDDO
C
	RETURN
	END



	SUBROUTINE PRUNE_BEST_RATIOS(IH,IK,RH,RK,NVECS,ITEST)
C
C Prune IH & IK(3..NVECS) to only include the ITEST elements
C with the best agreement between RH/RK = IH/IK. The pruning
C is done by marking IH & IK = 0 for the rejected spots.
C
	INTEGER IH(NVECS),IK(NVECS)
	REAL RH(NVECS),RK(NVECS)
C
	REAL RDIFF(1000),RDIFF2(1000)
C
C Calculate the fractional difference between each RH()/RK()
C and IH()/IK() and put into RDIFF() & RDIFF2().
C
	DO I=3,NVECS
	  IF(IH(I).EQ.0 .OR. IK(I).EQ.0) THEN
	    RDIFF(I)=1E6
	  ELSE
	    RDIFF(I)=ABS( RK(I)/RH(I)/IK(I)*IH(I) - 1.0 )
	    IF(ABS(IK(I)) .GT. 9) CALL QUIT('BUG: WEIRD')
	  ENDIF
	  RDIFF2(I)=RDIFF(I)
	ENDDO
C
C Do bubble sort of RDIFF2() and set CUTOFF to
C the ITEST-smallest value
C
	DO I1=3,NVECS
	  DO I2=I1+1,NVECS
	    IF(RDIFF2(I1) .GT. RDIFF2(I2)) THEN
	      RSAVE=RDIFF2(I1)
	      RDIFF2(I1)=RDIFF2(I2)
	      RDIFF2(I2)=RSAVE
	    ENDIF
	  ENDDO
	ENDDO
	CUTOFF=RDIFF2(ITEST)
C
C Change INUMER(),IDENOM(),RAT(), and NRAT to only include values
C where RDIFF is less than CUTOFF.
C
	DO I=3,NVECS
	  IF(ABS(RDIFF(I)) .GT. CUTOFF) THEN
	    IH(I)=0
	    IK(I)=0
	  ENDIF
	ENDDO
C
	RETURN
	END


	SUBROUTINE CALC_INITIAL_HKS(IH,IK, RH,RK,NVECS)
C
C Calculate H & K for spots from RH & RK
C
	INTEGER IH(NVECS),IK(NVECS)
	REAL RH(NVECS),RK(NVECS)
C
	INTEGER IHK(2)
C
C First 2 spots are the nodals
C
	IH(1)=1
	IK(1)=0
	IH(2)=0
	IK(2)=1
C
C Calculate remaining H & K from RH & RK
C
	IHK_MAX=9
	IHK_PROD_MAX=30
	DO I=3,NVECS
	  CALL GET_BEST_HK(IHK, RH(I),RK(I),IHK_MAX,IHK_PROD_MAX)
	  IH(I)=IHK(1)
	  IK(I)=IHK(2)
	ENDDO
C
	RETURN
	END


	SUBROUTINE GET_BEST_HK(IHK, RH,RK,IHK_MAX,IHK_PROD_MAX)
C
C Find the best integers up to +/- IHK_MAX that approximate the
C the fraction of RH & RK. The product of the integers must not
C exceed +/- IHK_PROD_MAX.
C
	INTEGER IHK(2)
C
	RHK_MAX=MAX(ABS(RH),ABS(RK))
	RH2=RH/RHK_MAX
	RK2=RK/RHK_MAX
C
	IBEST=0
	DBEST=1E6
	DO I=1,IHK_MAX
	  IH=NINT(RH2*I)
	  IK=NINT(RK2*I)
	  DIFF=ABS(IH*RK2-IK*RH2)*SQRT(FLOAT(I))
	  IF(DIFF .LT. DBEST) THEN
	    IF(ABS(IH*IK) .LE. IHK_PROD_MAX) THEN
	      IBEST=I
	      DBEST=DIFF
	    ENDIF
	  ENDIF
	ENDDO
C
	IHK(1)=NINT(RH2*IBEST)
	IHK(2)=NINT(RK2*IBEST)
C
	RETURN
	END


	SUBROUTINE GET_HK_SUM(ISUM, V_NODAL1,V_NODAL2, HVECS,NVECS)
C
C Calculate the sum of the magnitude of integer H & K for all "good" spots
C and add a penalty if the pixel position of a spot is poorly estimated
C
	REAL V_NODAL1(3),V_NODAL2(3),HVECS(3,10000)
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS0
C
	INTEGER IHK(2)
	REAL VEC1(3),PIXWAV(3)
C
C
	SPOT_MIN=5.0
	IHK_MAX=20
	IHK_PROD_MAX=50
	DIFF_MAX=10.0
C
C Precalculate some geometric factors
C
	DOT11=VC3_DOT(V_NODAL1,V_NODAL1)
	DOT12=VC3_DOT(V_NODAL1,V_NODAL2)
	DOT22=VC3_DOT(V_NODAL2,V_NODAL2)
C
C Loop over the spots
C
	ISUM=0
	DO I=1,NVECS
C Ignore weak spots
	  IF(SPOTS(3,I) .LT. SPOT_MIN) CYCLE
C Calculate the "real unscaled" values of H & K
	  DOT1H=VC3_DOT(V_NODAL1,HVECS(1,I))
	  DOT2H=VC3_DOT(V_NODAL2,HVECS(1,I))
	  RH=DOT1H*DOT22-DOT2H*DOT12
	  RK=DOT2H*DOT11-DOT1H*DOT12
C Calculate pixel position from RH & RK values
	  VEC1(1)=RH*V_NODAL1(1) + RK*V_NODAL2(1)
	  VEC1(2)=RH*V_NODAL1(2) + RK*V_NODAL2(2)
	  VEC1(3)=RH*V_NODAL1(3) + RK*V_NODAL2(3)
	  CALL CALC_HVECS_TO_PIXWAVS(PIXWAV, VEC1,1)
C If pixel positions are poorly estimated, ignore this spot
	  DIFF=SQRT( (SPOTS(1,I)-PIXWAV(1))**2 + (SPOTS(2,I)-PIXWAV(2))**2 )
	  IF(DIFF .GT. DIFF_MAX) CYCLE
C Calculate the best integers to approximate RH & RK
	  CALL GET_BEST_HK(IHK, RH,RK,IHK_MAX,IHK_PROD_MAX)
C Calculate pixel position from integer H & K values
	  VEC1(1)=IHK(1)*V_NODAL1(1) + IHK(2)*V_NODAL2(1)
	  VEC1(2)=IHK(1)*V_NODAL1(2) + IHK(2)*V_NODAL2(2)
	  VEC1(3)=IHK(1)*V_NODAL1(3) + IHK(2)*V_NODAL2(3)
	  CALL CALC_HVECS_TO_PIXWAVS(PIXWAV, VEC1,1)
C Add the absolute values of IHK to ISUM, unless the pixel positions
C were poorly estimated and then use a constant penalty of 100
	  DIFF=SQRT( (SPOTS(1,I)-PIXWAV(1))**2 + (SPOTS(2,I)-PIXWAV(2))**2 )
	  IF(DIFF .GT. DIFF_MAX) THEN
	    ISUM=ISUM+100
	  ELSE
	    ISUM=ISUM+ABS(IHK(1))+ABS(IHK(2))
	  ENDIF
C
	ENDDO
C
	RETURN
	END


	SUBROUTINE SCALE_NODALS_MIN_HK(V_NODAL1,V_NODAL2, HVECS,NVECS)
C
C Scale V_NODAL2 by a simple ratio that minimises the absolute sum of H & K
C
	REAL V_NODAL1(3),V_NODAL2(3),HVECS(3,10000)
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS0
C
	REAL V_TRY(3),SCALE(11)
	DATA SCALE /0.25,0.33333,0.5,0.66666,0.75,1.0,1.33333,1.5,2.0,3.0,4.0/
C
C Test scaling nodal2 by simple ratios
C
	ISUM_BEST=100000
	DO I=1,11
	  CALL VC3_SCALE(V_TRY, V_NODAL2, SCALE(I))
	  CALL GET_HK_SUM(ISUM, V_NODAL1,V_TRY, HVECS,NVECS)
	  IF(ISUM .LT. ISUM_BEST) THEN
	    ISUM_BEST=ISUM
	    IBEST=I
	  ENDIF
	ENDDO
C
C Output results
C
	CALL VC3_SCALE(V_NODAL2, V_NODAL2, SCALE(IBEST))
	PRINT '(/,1X,A,F7.3,A)','Scaling Nodal2 by',SCALE(IBEST),' to minimize HK values'
	PRINT '(/,1X,A,3F7.3,3X,3F7.3)','a*, b* =',V_NODAL1,V_NODAL2
	A_B=VC3_SIZE(V_NODAL1)/VC3_SIZE(V_NODAL2)
	GAMMA=ACOSD( VC3_DOT(V_NODAL1,V_NODAL2)/VC3_SIZE(V_NODAL1)/VC3_SIZE(V_NODAL2) )
	PRINT '(4X,A,F6.3,3X,A,F6.1,/)','|a*| / |b*| =',A_B,'angle =',GAMMA
C
	RETURN
	END



	SUBROUTINE CALC_INDEX_HK(IHK, V_NODAL1,V_NODAL2,HVECS,NVECS)
C
C Index the scattering vectors, HVECS(), for spots on the conic
C between the two nodal spots, V_NODAL1 & V_NODAL2.
C
	INTEGER IHK(2,1000)
	REAL V_NODAL1(3),V_NODAL2(3),HVECS(3,10000)
C
	IHK_MAX=20
	IHK_PROD_MAX=50
C
	DOT11=VC3_DOT(V_NODAL1,V_NODAL1)
	DOT12=VC3_DOT(V_NODAL1,V_NODAL2)
	DOT22=VC3_DOT(V_NODAL2,V_NODAL2)
C
	DO I=1,NVECS
C Calculate the "real unscaled" values of H & K
	  DOT1H=VC3_DOT(V_NODAL1,HVECS(1,I))
	  DOT2H=VC3_DOT(V_NODAL2,HVECS(1,I))
	  RH=DOT1H*DOT22-DOT2H*DOT12
	  RK=DOT2H*DOT11-DOT1H*DOT12
C Calculate the best integers to approximate RH & RK
	  CALL GET_BEST_HK(IHK(1,I), RH,RK,IHK_MAX,IHK_PROD_MAX)
	ENDDO
C
	RETURN
	END



	SUBROUTINE SCALE_NODALS_WAVS(WAVS,IHK, V_NODAL1,V_NODAL2,NVECS)
C
C Index the scattering vectors, HVECS, for spots on a conic
C between two nodal spots, V_NODAL1 & V_NODAL2.
C
	INTEGER IHK(2,1000)
	REAL V_NODAL1(3),V_NODAL2(3),WAVS(1000)
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS0
C
	REAL VEC(3),PIXWAV(3)
C
C Calculate the wavelengths for reasonable spots (else set to 1E6)
C
	DO I=1,NVECS
C
	  WAVS(I)=1E6
	  IF(SPOTS(3,I) .LE. 5.0) CYCLE
C
	  DO I2=1,3
	    VEC(I2)=IHK(1,I)*V_NODAL1(I2) + IHK(2,I)*V_NODAL2(I2)
	  ENDDO
	  CALL CALC_HVECS_TO_PIXWAVS(PIXWAV, VEC,1)
C
	  DIFF=SQRT( (SPOTS(1,I)-PIXWAV(1))**2 + (SPOTS(2,I)-PIXWAV(2))**2 )
	  IF(DIFF .GE. 10.0) CYCLE
C
	  WAVS(I)=-2.0*VEC(3)/VC3_SIZE(VEC)
	  PRINT '(2I4,3X,2F7.1,3X,3F7.1)',IHK(1,I),IHK(2,I),
	1			SPOTS(1,I),SPOTS(2,I),PIXWAV(1),PIXWAV(2),DIFF
C
	ENDDO
C
C Average the 10% smallest WAVS() values
C
	NOK=0
	DO I=1,NVECS
	  IF( WAVS(I) .LT. 1E5 ) NOK=NOK+1
	ENDDO
C
	NTRY=1+NOK/10
	SUM=0.0
	DO I1=1,NTRY
	  DO I2=I1+1,NVECS
	    IF( WAVS(I2) .LT. WAVS(I1)) THEN
	      WSAVE=WAVS(I1)
	      WAVS(I1)=WAVS(I2)
	      WAVS(I2)=WSAVE
	    ENDIF
	  ENDDO
	  SUM=SUM+WAVS(I1)
	ENDDO
	AVE=SUM/NTRY
C
C Assume shortest wavelengths are ~1 Angstrom
C
	SCALE=AVE/1.00
C
C Apply scale to all cell lengths
C
	CALL VC3_SCALE(V_NODAL1, V_NODAL1, SCALE)
	CALL VC3_SCALE(V_NODAL2, V_NODAL2, SCALE)
	PRINT '(/,1X,A,F7.3,A)','Scaling Nodals by',SCALE,' from minimum wavelengths'
	PRINT '(/,1X,A,3F7.3,3X,3F7.3)','a*, b* =',V_NODAL1,V_NODAL2
C
C Scale all wavelengths for the new nodal vectors
C
	DO I=1,NVECS
	  IF( WAVS(I) .LT. 1E5 ) WAVS(I)=WAVS(I)/SCALE
	ENDDO
C
C Output direct-space cell dimensions
C
	A_B=VC3_SIZE(V_NODAL1)/VC3_SIZE(V_NODAL2)
	A=1.0/VC3_SIZE(V_NODAL1)
	B=1.0/VC3_SIZE(V_NODAL2)
	GAMMA=ACOSD( VC3_DOT(V_NODAL1,V_NODAL2)/VC3_SIZE(V_NODAL1)/VC3_SIZE(V_NODAL2) )
	PRINT '(1X,A,2F6.2,F6.1,/)','a, b, angle =',A,B,GAMMA
C
	RETURN
	END



	SUBROUTINE GET_NODAL_RATIO(RATIO, RH,RK,NVECS,ITEST)
C
C Try dividing all ratios RAT(1..NRAT) by each element of RAT()
C to give RAT2(). Test how well RAT2() are approximated by simple
C fractions and return in SCALE the best RAT() value used as a
C divisor for RAT2(). ITEST is used by GET_DENOM_MERIT().
C
	REAL RH(NVECS),RK(NVECS)
C
	REAL RAT(1000),RAT2(1000)
C
C Create ratios of RK & RH ignoring the first two spots
C which are the nodals
C
	NRAT=NVECS-2
	DO I=1,NRAT
	  RAT(I)=RK(I+2)/RH(I+2)
	ENDDO
C
C The scaling the ratios by each value in RAT() and save
C the best results.
C
	RBEST=1E6
	DO I1=1,NRAT
C Make new ratios scaled so that element I1 is equal to 1
	  DO I2=1,NRAT
	    RAT2(I2)=RAT(I2)/RAT(I1)
	  ENDDO
C Test the merit of the scaled ratios and save the best one
	  CALL GET_DENOM_MERIT(RMERIT,RAT2,NRAT,ITEST)
	  IF(RMERIT .LT. RBEST) THEN
	    RBEST=RMERIT
	    IBEST=I1
	  ENDIF
C
	ENDDO
C
C Return the best scale factor
C
	RATIO=ABS(RAT(IBEST))
	RETURN
	END



	SUBROUTINE GET_DENOM_MERIT(RMERIT,RAT,NRAT,ITEST)
C
C Given an array of potential ratios, RAT(1..NRAT) return
C RMERIT which is how well the best ITEST ratios are fitted
C by fractions containing numbers up to +/- 9.
C
	REAL RAT(NRAT)
C
	REAL DBEST(1000)
C
C For each ratio, try to fit to a fraction using integers less
C than +/- NMULT_MAX. Put in DBEST() the difference between the
C integer and real value of the smaller number of the fraction.
C
	NMULT_MAX=9
	DO I=1,NRAT
	  RATIO=ABS(RAT(I))
	  IF(RATIO .GT. 1.0) RATIO=1.0/RATIO
	  DBEST(I)=1E6
	  DO IMULT=1,NMULT_MAX
	    DIFF=ABS(RATIO*IMULT-NINT(RATIO*IMULT))
	    DBEST(I)=MIN(DBEST(I),DIFF)
	  ENDDO
	ENDDO
C
C Do bubble sort of DBEST()
C
	DO I1=1,ITEST
	  DO I2=I1+1,NRAT
	    IF(DBEST(I1) .GT. DBEST(I2)) THEN
	      DSAVE=DBEST(I1)
	      DBEST(I1)=DBEST(I2)
	      DBEST(I2)=DSAVE
	    ENDIF
	  ENDDO
	ENDDO
C
C Return the ITEST-smallest value of DBEST
C
	RMERIT=DBEST(ITEST)
	RETURN
	END
