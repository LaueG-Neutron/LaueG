	SUBROUTINE TAG_SORT_REALS(ITAGS,REALS0,NREALS)
C
	INTEGER ITAGS(NREALS)
	REAL REALS0(NREALS)
C
	EXTERNAL REALS_SORT_COMPARE
C
	PARAMETER (NREALS_MAX=10000)
	COMMON /REALS_SORT_COM/ REALS(NREALS_MAX)
C
	IF(NREALS .GT. NREALS_MAX) CALL QUIT('BUG(tag_sort_reals): NREALS too big')
C
	DO I=1,NREALS
	  REALS(I)=REALS0(I)
	ENDDO
C
	CALL TAG_SORT(REALS_SORT_COMPARE,ITAGS,NREALS)
C
	RETURN
	END



	FUNCTION REALS_SORT_COMPARE(I1,I2)
C
C A positive return value means that item I1 should appear before I2.
C
	PARAMETER (NREALS_MAX=10000)
	COMMON /REALS_SORT_COM/ REALS(NREALS_MAX)
C
	REALS_SORT_COMPARE=REALS(I2)-REALS(I1)
C
	RETURN
	END



	FUNCTION ESDS_COMPARE(I1,I2)
C
C A positive return value means that item I1 should appear before I2.
C
	PARAMETER (NSPOTS_MAX=10000)
	COMMON /SPOTS_COM/ NSPOTS,XSPOT(NSPOTS_MAX),YSPOT(NSPOTS_MAX),ESDS(NSPOTS_MAX)
C
	ESDS_COMPARE=ESDS(I1)-ESDS(I2)
C
	RETURN
	END



	FUNCTION DIFF_ANGLE_COMPARE(I1,I2)
C
C Compares solutions for  triplets I1 & I2. The comparision is the
C averaged difference in angles between observed scattering vectors
C and those generated from HKLs for the given triplets.
C
	PARAMETER (NTRIP_MAX=100000)
	COMMON /DIFF_ANGLE_COM/ DIFF(NTRIP_MAX)
C
	DIFF_ANGLE_COMPARE=DIFF(I2)-DIFF(I1)
	RETURN
	END



	FUNCTION CHECK_ANGLE_COMPARE(I1,I2)
C
C Compares the difference in angles between observed scattering vectors
C and those generated from HKLs for a given triplet "solution".
C I1,I2 refers to the indices of the peaks we are "checking".
C
	PARAMETER (NCHECK_MAX=100)
	COMMON /CHECK_COMPARE_COM/ CHECK_ANG(NCHECK_MAX)
C
	CHECK_ANGLE_COMPARE=CHECK_ANG(I2)-CHECK_ANG(I1)
C
	RETURN
	END



	FUNCTION AZI_ROT_COMPARE(I1,I2)
C
C Compares the azimuthal rotation of two vectors
C
	PARAMETER (NPAIRS_MAX=1000)
	COMMON /AZI_ROT_COM/ AZI_ROT(NPAIRS_MAX)
C
	AZI_ROT_COMPARE=AZI_ROT(I2)-AZI_ROT(I1)
	RETURN
	END



	FUNCTION PAIR_ANGLE_COMPARE(I1,I2)
C
C Compare the angles between pairs of HKL generated vectors
C
	PARAMETER (NANGLES_MAX=150000)
	COMMON /ANGLES_COM/ NANGLES,ANGLES(NANGLES_MAX),IANGLES(2,NANGLES_MAX)
C
	PAIR_ANGLE_COMPARE=ANGLES(I2)-ANGLES(I1)
	RETURN
	END



	SUBROUTINE TAG_SORT(SORT_FUNC,TAG,NTOT)
C
C Sort NTOT items and return TAG() which is the tag-index to
C the items sorted by the real function SORT_FUNC().
C SORT_FUNC(I1,I2) returns positive if item I1 should appear
C before item I2 in the tagged list, zero if items are identical,
C and negative if I1 should appear after I2.
C
	INTEGER TAG(NTOT)
C
	INTEGER TAG0
C
C INITIALIZE TAGS
C
	DO I=1,NTOT
	  TAG(I)=I
	ENDDO
C
C TREE SORT (PART A)
C
	DO 240 I1=2,NTOT/2
	I2=NTOT/2-I1+2
	TAG0=TAG(I2)
200	I3=I2+I2
	IF(I3-NTOT) 220,230,240
220	IF(SORT_FUNC(TAG(I3+1),TAG(I3)) .GT. 0.0) GO TO 230
	I3=I3+1
230	IF(SORT_FUNC(TAG(I3),TAG0) .GT. 0.0) GO TO 240
	TAG(I2)=TAG(I3)
	I2=I3
	GO TO 200
240	TAG(I2)=TAG0
C
C TREE SORT (PART B)
C
	DO 340 I1=2,NTOT
	I2=NTOT-I1+2
	I3=1
	TAG0=TAG(I3)
300	I4=I3+I3
	IF(I4-I2) 310,320,330
310	IF(SORT_FUNC(TAG(I4+1),TAG(I4)) .GT. 0.0) GO TO 320
	I4=I4+1
320	IF(SORT_FUNC(TAG(I4),TAG0) .GT. 0.0) GO TO 330
	TAG(I3)=TAG(I4)
	I3=I4
	GO TO 300
330	TAG(I3)=TAG0
	TAG0=TAG(1)
	TAG(1)=TAG(I2)
	TAG(I2)=TAG0
340	CONTINUE
C
	RETURN
	END
