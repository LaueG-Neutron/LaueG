C ==================
C	SUBROUTINE SORT_PFRAC_LIST()
C	FUNCTION SORT_PFRAC_FUNC(I1,I2)
C	SUBROUTINE TAG_SORT(SORT_FUNC,TAG,NTOT)
C ==================

	SUBROUTINE SORT_PFRAC_LIST()
C
	EXTERNAL SORT_PFRAC_FUNC
C
	COMMON /PFRAC_CHECK_COM/ PFRAC_LIST(3,20000),PFRAC_TARGET,NPFRAC
C
	INTEGER ITAGS(20000)
	REAL TEMP(20000)
C
C Sort the reflections using a tagged sort
C
	CALL TAG_SORT(SORT_PFRAC_FUNC,ITAGS,NPFRAC)
C
C Reorder the reflections according to ITAGS()
C
	DO I1=1,3
	  DO I2=1,NPFRAC
	    TEMP(I2)=PFRAC_LIST(I1,ITAGS(I2))
	  ENDDO
	  DO I2=1,NPFRAC
	    PFRAC_LIST(I1,I2)=TEMP(I2)
	  ENDDO
	ENDDO
C
	RETURN
	END


	FUNCTION SORT_PFRAC_FUNC(I1,I2)
C
	COMMON /PFRAC_CHECK_COM/ PFRAC_LIST(3,20000),PFRAC_TARGET,NPFRAC
C
	RAT1=(PFRAC_LIST(1,I1)-PFRAC_TARGET)/PFRAC_LIST(3,I1)
	RAT2=(PFRAC_LIST(1,I2)-PFRAC_TARGET)/PFRAC_LIST(3,I2)
	SORT_PFRAC_FUNC=RAT2-RAT1
C
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

