C ========== Order functions needed for TAG_SORT =================
C	FUNCTION ORDER_REF(IREF1,IREF2)
C	FUNCTION ORDER_SEQ(HKLM1,HKLM2)
C	FUNCTION ORDER_HKL(HKL1,HKL2)
C	FUNCTION SORT_WEAK_FUNC(I1,I2)
C	FUNCTION SORT_BIN_DATA_FUNC(I1,I2)
C ========== Base function to do tagged binary sort ==============
C	SUBROUTINE TAG_SORT(SORT_FUNC,TAG,NTOT)
C ================================================================


C ========== Order functions needed for TAG_SORT =================

	FUNCTION ORDER_REF(IREF1,IREF2)
C
C Returns +|0|- if IREF1 should appear before|identical|after IREF2.
C Compares in order:	sequence number; twin number, if applicable;
C									actual HKL; wavelength (reversed)
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
C First order spots by the sequences they belong to
C
	ORDER_REF=ORDER_SEQ(HKLMS(1,IREF1),HKLMS(1,IREF2))
	IF(ORDER_REF .NE. 0.0) RETURN
C
C If twinned, order versus twin number
C
	IF(IDATA_MODE .EQ. 2) THEN
		ORDER_REF=HKLMS(4,IREF2)-HKLMS(4,IREF1)
		IF(ORDER_REF .NE. 0.0) RETURN
	ENDIF
C
C Order versus HKL indices (ignoring M)
C
	ORDER_REF=ORDER_HKL(HKLMS(1,IREF1),HKLMS(1,IREF2))
	IF(ORDER_REF .NE. 0.0) RETURN
C
C Finally, order on wavelengths (reversed)
C
	ORDER_REF=WAVS(IREF2)-WAVS(IREF1)
C
	RETURN
	END


	FUNCTION ORDER_SEQ(HKLM1,HKLM2)
C
C Return +|0|- if HKLM1's sequence is before/equal/after HKLM2's sequence
C
	INTEGER HKLM1(4),HKLM2(4)
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	INTEGER HKL1(3),HKL2(3),EQV1(3),EQV2(3)
	REAL VEC1(3),VEC2(3)
C
C Copy hkl indices to HKL1 & 2
C
	IF(IDATA_MODE .NE. 3) THEN
		DO I=1,3
			HKL1(I)=HKLM1(I)
			HKL2(I)=HKLM2(I)
		ENDDO
	ELSE
C For modulated, use indices = 1000 * (HKL + modulation)
C NB: Works for modulation precision down to 0.001
		CALL CALC_MOD_QVEC(VEC1, HKLM1(4) )
		CALL CALC_MOD_QVEC(VEC2, HKLM2(4) )
		DO I=1,3
		  HKL1(I)=NINT( 1E3*( HKLM1(I) + VEC1(I) ) )
		  HKL2(I)=NINT( 1E3*( HKLM2(I) + VEC2(I) ) )
		ENDDO
	ENDIF
C
C Calculate unique equivalents of HKL1 & 2
C
	CALL CALC_EQUIV_HKL(EQV1,HKL1)
	CALL CALC_EQUIV_HKL(EQV2,HKL2)
C
C Simply order HKLs by indices
C
	ORDER_SEQ=ORDER_HKL(EQV1,EQV2)
C
	RETURN
	END


	FUNCTION ORDER_HKL(HKL1,HKL2)
C
C Returns +|0|- if HKL1 should appear before|identical|after HKL2
C Compares the HKLs given as arguments
C NB: Ignores the M value
C
	INTEGER HKL1(3),HKL2(3)
C
	ICOMP=HKL2(1)-HKL1(1)
	IF(ICOMP .EQ. 0) ICOMP=HKL2(2)-HKL1(2)
	IF(ICOMP .EQ. 0) ICOMP=HKL2(3)-HKL1(3)
C
	ORDER_HKL=ICOMP
C
	RETURN
	END


	FUNCTION SORT_WEAK_FUNC(I1,I2)
C
C A positive return value means that item I1 should appear before I2.
C Used by PRUNE_WEAK_SEQ.
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SORT_WEAK_COM/ CSEQ(NSEQ_MAX)
C
	SORT_WEAK_FUNC=CSEQ(I2)-CSEQ(I1)
C
	RETURN
      END


	FUNCTION SORT_BIN_DATA_FUNC(I1,I2)
C
C A positive return value means that item I1 should appear before I2.
C
	PARAMETER (NDATA_MAX=2000000)
	COMMON /SORT_BIN_COM/ BIN_DATA(NDATA_MAX)
C
	SORT_BIN_DATA_FUNC=BIN_DATA(I2)-BIN_DATA(I1)
C
	RETURN
	END


C ========== Base function to do tagged binary sort ==============
	SUBROUTINE TAG_SORT(SORT_FUNC,TAG,NTOT)
C
C Sort NTOT items and return TAG() which is the tag-index to
C the items sorted by the real function SORT_FUNC.
C SORT_FUNC(I1,I2) returns positive if item I1 should appear
C before item I2 in the tagged list, zero if items are identical,
C and negative if I1 should appear after I2.
C
	INTEGER TAG(NTOT)
C
	INTEGER TAG0
C
	IF(NTOT .EQ. -1) THEN
	  TAG(1)=1
	  RETURN
	ENDIF
C
	IF(NTOT .EQ. -2) THEN
	  IF(SORT_FUNC(1,2) .GE. 0.0) THEN
	    TAG(1)=1
	    TAG(2)=2
	  ELSE
	    TAG(1)=2
	    TAG(2)=1
	  ENDIF
	  RETURN
	ENDIF
C
C INITIALIZE TAGS
C
	DO I=1,NTOT
	  TAG(I)=I
	ENDDO
C
C TREE SORT (PART A)
C
	DO I1=2,NTOT/2
	  I2=NTOT/2-I1+2
	  TAG0=TAG(I2)
C
100	  I3=2*I2
	  IF(I3.LT. NTOT) THEN
	    IF(SORT_FUNC(TAG(I3+1),TAG(I3)) .LE. 0.0) I3=I3+1
	  ENDIF
	  IF(I3 .LE. NTOT) THEN
	    IF(SORT_FUNC(TAG(I3),TAG0) .LE. 0.0) THEN
	      TAG(I2)=TAG(I3)
	      I2=I3
	      GO TO 100
	    ENDIF
	  ENDIF
C
	  TAG(I2)=TAG0
	ENDDO
C
C TREE SORT (PART B)
C
	DO I1=2,NTOT
	  I2=NTOT-I1+2
	  I3=1
	  TAG0=TAG(I3)
C
200	  I4=2*I3
	  IF(I4 .LT. I2) THEN
	    IF(SORT_FUNC(TAG(I4+1),TAG(I4)) .LE. 0.0) I4=I4+1
	  ENDIF
	  IF(I4 .LE. I2) THEN
	    IF(SORT_FUNC(TAG(I4),TAG0) .LE. 0.0) THEN
	      TAG(I3)=TAG(I4)
	      I3=I4
	      GO TO 200
	    ENDIF
	  ENDIF
C
	  TAG(I3)=TAG0
	  TAG0=TAG(1)
	  TAG(1)=TAG(I2)
	  TAG(I2)=TAG0
	ENDDO
C
	RETURN
	END
