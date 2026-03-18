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
