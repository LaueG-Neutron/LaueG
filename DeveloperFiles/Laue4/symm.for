C----- Symmetry related routines -------------------------------------------
c	SUBROUTINE ASK_CELL_TYPE
c	SUBROUTINE GET_CELL_TYPE(STRING)
c	SUBROUTINE CALC_EQUIV_HKL(EQV, HKL)
C---------------------------------------------------------------------------

	SUBROUTINE ASK_CELL_TYPE
C
	COMMON /CELL_TYPE_COM/ ITYPE,IFRIEDEL,ITYPE_ORIG,ICEN,CELL(6)
C
	DATA ICUB,ITETA,ITETB,ITETC,IORTH,IHEX,ITRIG,IRHOM,IMONA,
	1		IMONB,IMONC,IDENT/1,2,3,4,5,6,7,8,9,10,11,12/
C
C Output the merging options available for the symmetry type
C
100	NOPTS=4
	PRINT '(/,1X,A)','Merge reflection using:'
	PRINT '(1X,A)',' 1  Identical HKLs'
	PRINT '(1X,A)',' 2  Identical HKLs + Friedels'
	IF     (ITYPE_ORIG .EQ. ICUB) THEN
	  PRINT '(1X,A)',' 3  Cubic symmetry'
	  PRINT '(1X,A)',' 4  Cubic symmetry + Friedels'
	ELSE IF(ITYPE_ORIG.EQ.ITETA .OR.
	1		  ITYPE_ORIG.EQ.ITETB .OR. ITYPE_ORIG.EQ.ITETC) THEN
	  PRINT '(1X,A)',' 3  Tetragonal symmetry'
	  PRINT '(1X,A)',' 4  Tetragonal symmetry + Friedels'
	ELSE IF(ITYPE_ORIG .EQ. IORTH) THEN
	  PRINT '(1X,A)',' 3  Orthorhombic symmetry'
	  PRINT '(1X,A)',' 4  Orthorhombic symmetry + Friedels'
	ELSE IF(ITYPE_ORIG.EQ.IHEX .OR. ITYPE_ORIG.EQ.ITRIG) THEN
	  NOPTS=6
	  PRINT '(1X,A)',' 3  Trigonal symmetry'
	  PRINT '(1X,A)',' 4  Trigonal symmetry + Friedels'
	  PRINT '(1X,A)',' 5  Hexagonal symmetry'
	  PRINT '(1X,A)',' 6  Hexagonal symmetry + Friedels'
	ELSE IF(ITYPE_ORIG .EQ. IRHOM) THEN
	  PRINT '(1X,A)',' 3  Rhombohedral symmetry'
	  PRINT '(1X,A)',' 4  Rhombohedral symmetry + Friedels'
	ELSE IF(ITYPE_ORIG.EQ.IMONA .OR.
	1		  ITYPE_ORIG.EQ.IMONB .OR. ITYPE_ORIG.EQ.IMONC) THEN
	  PRINT '(1X,A)',' 3  Monoclinic symmetry'
	  PRINT '(1X,A)',' 4  Monoclinic symmetry + Friedels'
	ELSE IF(ITYPE_ORIG .EQ. IDENT) THEN
	  NOPTS=2
	ELSE
	  STOP 'BUG(ask_cell_type): Invalid ITYPE'
	ENDIF
C
C Ask for merging option, loop back to ask again if answer invalid
C
	PRINT '(1X,A,I1,A,$)','Selection option 1 - ',NOPTS,' : '
	READ(*,*,ERR=100) IOPT
	IF(IOPT.LT.1 .OR. IOPT.GT.NOPTS) GOTO 100
C
C Set switch for averaging Friedel pairs
C
	IFRIEDEL=0
	IF(IOPT .EQ. IOPT/2*2) IFRIEDEL=1
C
C Change symmetry if user requests it
C
	IOPT=(IOPT+1)/2
	IF(IOPT .EQ. 1) THEN
	  ITYPE=IDENT
	ELSE IF(NOPTS .EQ. 6) THEN
	  IF(IOPT .EQ. 2) ITYPE=ITRIG
	  IF(IOPT .EQ. 3) ITYPE=IHEX
	ELSE
	  ITYPE=ITYPE_ORIG
	ENDIF
C
	RETURN
	END



	SUBROUTINE GET_CELL_TYPE(STRING)
C
	CHARACTER*(*) STRING
C
	COMMON /CELL_TYPE_COM/ ITYPE,IFRIEDEL,ITYPE_ORIG,ICEN,CELL(6)
C
	DATA ICUB,ITETA,ITETB,ITETC,IORTH,IHEX,ITRIG,IRHOM,IMONA,
	1		IMONB,IMONC,IDENT/1,2,3,4,5,6,7,8,9,10,11,12/
C
C Fill STRING with blanks
C
	STRING=REPEAT(' ',LEN(STRING))
C
C Copy cell system into STRING
C
	IF     (ITYPE .EQ. ICUB) THEN
	  STRING='Cubic'
	ELSE IF(ITYPE .EQ. ITETA) THEN
	  STRING='Tetragonal(a axis)'
	ELSE IF(ITYPE .EQ. ITETB) THEN
	  STRING='Tetragonal(b axis)'
	ELSE IF(ITYPE .EQ. ITETC) THEN
	  STRING='Tetragonal(c axis)'
	ELSE IF(ITYPE .EQ. IORTH) THEN
	  STRING='Orthorhombic'
	ELSE IF(ITYPE .EQ. IHEX) THEN
	  STRING='Hexagonal'
	ELSE IF(ITYPE .EQ. ITRIG) THEN
	  STRING='Trigonal(hex)'
	ELSE IF(ITYPE .EQ. IRHOM) THEN
	  STRING='Trigonal(rhom)'
	ELSE IF(ITYPE .EQ. IMONA) THEN
	  STRING='Monoclinic(a axis)'
	ELSE IF(ITYPE .EQ. IMONB) THEN
	  STRING='Monoclinic(b axis)'
	ELSE IF(ITYPE .EQ. IMONC) THEN
	  STRING='Monoclinic(c axis)'
	ELSE IF(ITYPE .EQ. IDENT) THEN
	  STRING='Triclinic'
	ELSE
	  PRINT *,'ITYPE=',ITYPE
	  STOP 'BUG(get_cell_type): Invalid ITYPE'
	ENDIF
C
C Add "+ Friedels"
C
	IF(IFRIEDEL .EQ. 1) THEN
	  ILEN=LEN_TRIM(STRING)
	  STRING(ILEN+1:ILEN+11)=' + Friedels'
	ENDIF
C
	RETURN
	END


	SUBROUTINE CALC_EQUIV_HKL(EQV, HKL)
C NB: Ignores the M values
	INTEGER EQV(3),HKL(3)
C
	COMMON /CELL_TYPE_COM/ ITYPE,IFRIEDEL,ITYPE_ORIG,ICEN,CELL(6)
C
	DATA ICUB,ITETA,ITETB,ITETC,IORTH,IHEX,ITRIG,IRHOM,IMONA,
	1		IMONB,IMONC,IDENT/1,2,3,4,5,6,7,8,9,10,11,12/
C
C Simplifies the code if we copy the indices at the beginning
C
	EQV(1)=HKL(1)
	EQV(2)=HKL(2)
	EQV(3)=HKL(3)
C
C Cubic(P23): HKL,LHK,H-K-L
C
	IF(ITYPE .EQ. ICUB) THEN
C Rearrange H,K,L until (|H|>|K| & |H|>|L|) or |H|=>|K|=>|L|
	  DO IPERM=1,3
C Stop permuting if |H|>|K| & |H|>|L|
	    IF(ABS(EQV(1)) .GT. MAX(ABS(EQV(2)),ABS(EQV(3)))) EXIT
C Stop permuting if |H|=|K| and |K|>=|L|
	    IF(ABS(EQV(1)).EQ.ABS(EQV(2)) .AND. ABS(EQV(2)).GE.ABS(EQV(3))) EXIT
C Cyclic permutation of H,K,L
	    ITEMP=EQV(1)
	    EQV(1)=EQV(2)
	    EQV(2)=EQV(3)
	    EQV(3)=ITEMP
C Bug check in case EXIT wasn't activated
	    IF(IPERM .EQ. 3) STOP 'BUG(calc_equiv_hkl): Impossible P23'
	  ENDDO
C If H<0, invert H,K
	  IF(EQV(1) .LT. 0) THEN
	    EQV(1)=-EQV(1)
	    EQV(2)=-EQV(2)
	  ENDIF
C If K<0, invert K,L
	  IF(EQV(2) .LT. 0) THEN
	    EQV(2)=-EQV(2)
	    EQV(3)=-EQV(3)
	  ENDIF
C If K=0, we can invert L if negative
	  IF(EQV(2) .EQ. 0) EQV(3)=ABS(EQV(3))
C If Friedel, we can invert L if negative
	  IF(IFRIEDEL .EQ. 1) EQV(3)=ABS(EQV(3))
C
C Tetragonal(P4,c-axis): HKL,-H-KL,-KHL
C
	ELSE IF(ITYPE .EQ. ITETC) THEN	! tetragonal (unique c)
C If |H| < |K|, swap HK to -KH
	  IF(ABS(HKL(1)) .LT. ABS(HKL(2))) THEN
	    EQV(1)=-HKL(2)
	    EQV(2)=HKL(1)
	  ENDIF
C Invert HK if H<0
	  IF(EQV(1) .LT. 0) THEN
	    EQV(1)=-EQV(1)
	    EQV(2)=-EQV(2)
	  ENDIF
C If H = -K, swap HK to -KH (=HH) 
	  IF(EQV(1) .EQ. -EQV(2)) EQV(2)=EQV(1)
C If Friedel, we can invert L if negative
	  IF(IFRIEDEL .EQ. 1) EQV(3)=ABS(EQV(3))
C
C Tetragonal(P4,b-axis): HKL,-HK-L,-LKH
C
	ELSE IF(ITYPE .EQ. ITETB) THEN	! tetragonal (unique b)
C If |H| < |L|, swap HL to -LH
	  IF(ABS(HKL(1)) .LT. ABS(HKL(3))) THEN
	    EQV(1)=-HKL(3)
	    EQV(3)=HKL(1)
	  ENDIF
C Invert HL if H<0
	  IF(EQV(1) .LT. 0) THEN
	    EQV(1)=-EQV(1)
	    EQV(3)=-EQV(3)
	  ENDIF
C If H = -L, swap HL to -LH (=HH) 
	  IF(EQV(1) .EQ. -EQV(3)) EQV(3)=EQV(1)
C If Friedel, we can invert K if negative
	  IF(IFRIEDEL .EQ. 1) EQV(2)=ABS(EQV(2))
C
C Tetragonal(P4,a-axis): HKL,H-K-L,H-LK
C
	ELSE IF(ITYPE .EQ. ITETA) THEN	! tetragonal (unique a)
C If |K| < |L|, swap KL to -LK
	  IF(ABS(HKL(2)) .LT. ABS(HKL(3))) THEN
	    EQV(2)=-HKL(3)
	    EQV(3)=HKL(2)
	  ENDIF
C Invert KL if K<0
	  IF(EQV(2) .LT. 0) THEN
	    EQV(2)=-EQV(2)
	    EQV(3)=-EQV(3)
	  ENDIF
C If K = -L, swap KL to -LK (=KK) 
	  IF(EQV(2) .EQ. -EQV(3)) EQV(3)=EQV(2)
C If Friedel, we can invert H if negative
	  IF(IFRIEDEL .EQ. 1) EQV(1)=ABS(EQV(1))
C
C Orthorhombic(P222): HKL,-H-KL,-HK-L
C
	ELSE IF(ITYPE .EQ. IORTH) THEN
C If H<0, invert H,K
	  IF(EQV(1) .LT. 0) THEN
	    EQV(1)=-EQV(1)
	    EQV(2)=-EQV(2)
	  ENDIF
C If K<0, invert K,L
	  IF(EQV(2) .LT. 0) THEN
	    EQV(2)=-EQV(2)
	    EQV(3)=-EQV(3)
	  ENDIF
C If H=0 or K=0, we can invert L if negative
	  IF(EQV(1).EQ.0 .OR. EQV(2).EQ.0) EQV(3)=ABS(EQV(3))
C If Friedel, we can invert L if negative
	  IF(IFRIEDEL .EQ. 1) EQV(3)=ABS(EQV(3))
C
C Hexagonal(P6): HKL,-H-KL,KIL
C
	ELSE IF(ITYPE .EQ. IHEX) THEN
	  I=-(HKL(1)+HKL(2))	! hexagonal I index
C If H,K or I is zero, cyclically rotate the zero to K
	  IF(HKL(1) .EQ. 0) THEN
	    EQV(1)=I
	    EQV(2)=HKL(1)
	  ELSE IF(HKL(2) .EQ. 0) THEN
	    CONTINUE		! nothing to do
	  ELSE IF(I      .EQ. 0) THEN
	    EQV(1)=HKL(2)
	    EQV(2)=I
C No zeroes, so choose rotation with largest absolute I
	  ELSE IF(ABS(I) .GT. MAX( ABS(HKL(1)) , ABS(HKL(2)) ) ) THEN
	    CONTINUE		! nothing to do
	  ELSE IF(ABS(HKL(1)) .GT. MAX( ABS(I) , ABS(HKL(2)) ) ) THEN
	    EQV(1)=HKL(2)
	    EQV(2)=I
	  ELSE IF(ABS(HKL(2)) .GT. MAX( ABS(I) , ABS(HKL(1)) ) ) THEN
	    EQV(1)=I
	    EQV(2)=HKL(1)
	  ELSE
	    STOP 'BUG(calc_equiv_hkl): Impossible P6'
	  ENDIF
C Invert HK if H<0
	  IF(EQV(1) .LT. 0) THEN
	    EQV(1)=-EQV(1)
	    EQV(2)=-EQV(2)
	  ENDIF
C If Friedel, we can invert L if negative
	  IF(IFRIEDEL .EQ. 1) EQV(3)=ABS(EQV(3))
C
C Trigonal(P3): HKL,KIL
C
	ELSE IF(ITYPE .EQ. ITRIG) THEN
	  I=-(HKL(1)+HKL(2))	! hexagonal I index
C If H,K or I is zero, cyclically rotate the zero to K
	  IF(HKL(1) .EQ. 0) THEN
	    EQV(1)=I
	    EQV(2)=HKL(1)
	  ELSE IF(HKL(2) .EQ. 0) THEN
	    CONTINUE		! nothing to do
	  ELSE IF(I      .EQ. 0) THEN
	    EQV(1)=HKL(2)
	    EQV(2)=I
C No zeroes, so choose rotation with largest absolute I
	  ELSE IF(ABS(I) .GT. MAX(ABS(HKL(1)) , ABS(HKL(2)) ) ) THEN
	    CONTINUE		! nothing to do
	  ELSE IF(ABS(HKL(1)) .GT. MAX(ABS(I) , ABS(HKL(2)) ) ) THEN
	    EQV(1)=HKL(2)
	    EQV(2)=I
	  ELSE IF(ABS(HKL(2)) .GT. MAX( ABS(I) , ABS(HKL(1)) ) ) THEN
	    EQV(1)=I
	    EQV(2)=HKL(1)
	  ELSE
	    STOP 'BUG(calc_equiv_hkl): Impossible P3'
	  ENDIF
C If Friedel, invert HKL if H<0 or invert L if a 00L
	  IF(IFRIEDEL .EQ. 1) THEN
	    IF(EQV(1) .LT. 0) THEN
	      EQV(1)=-EQV(1)
	      EQV(2)=-EQV(2)
	      EQV(3)=-EQV(3)
	    ELSE IF(EQV(1) .EQ. 0) THEN
	      EQV(3)=ABS(EQV(3))
	    ENDIF
	  ENDIF
C
C Rhombohedral(R3): HKL,KLH
C
	ELSE IF(ITYPE .EQ. IRHOM) THEN
C Special case where H,K,L are all equal magnitude
	  IF(ABS(HKL(1)).EQ.ABS(HKL(2)) .AND. ABS(HKL(1)).EQ.ABS(HKL(3))) THEN
C Cyclically rotate h,k,l so we have +++, ++-, +--, or ---
	    IF(HKL(1).EQ.HKL(2) .AND. HKL(1).EQ.HKL(3)) THEN
	      IF(IFRIEDEL .EQ. 1) THEN	! +++ or --- case
	        EQV(1)=ABS(HKL(1))
	        EQV(2)=ABS(HKL(1))
	        EQV(3)=ABS(HKL(1))
	      ENDIF
	    ELSE IF(HKL(1)*HKL(2)*HKL(3).LT.0 .OR. IFRIEDEL.EQ.1) THEN
	      EQV(1)=+ABS(HKL(1))			! ++- case or Friedel
	      EQV(2)=+ABS(HKL(1))
	      EQV(3)=-ABS(HKL(1))
	    ELSE
	      EQV(1)=+ABS(HKL(1))			! +-- case
	      EQV(2)=-ABS(HKL(1))
	      EQV(3)=-ABS(HKL(1))
	    ENDIF
C General case:
	  ElSE
C Rearrange H,K,L until (|H|>|K| & |H|>|L|) or |H|=>|K|=>|L|
	    DO IPERM=1,3
	      IF(ABS(EQV(1)) .GT. MAX(ABS(EQV(2)),ABS(EQV(3)))) EXIT
	      IF(ABS(EQV(1)).EQ.ABS(EQV(2)) .AND. ABS(EQV(2)).GE.ABS(EQV(3))) EXIT
	      ITEMP=EQV(1)
	      EQV(1)=EQV(2)
	      EQV(2)=EQV(3)
	      EQV(3)=ITEMP
C Bug check in case EXIT wasn't activated
	      IF(IPERM .EQ. 3) STOP 'BUG(calc_equiv_hkl): Impossible R3'
	    ENDDO
C If Friedel, invert HKL if H<0
	    IF(IFRIEDEL.EQ.1 .AND. EQV(1).LT.0) THEN
	      EQV(1)=-EQV(1)
	      EQV(2)=-EQV(2)
	      EQV(3)=-EQV(3)
	    ENDIF
	  ENDIF
C
C Monoclinic(P2,a-axis): HKL,H-K-L
C
	ELSE IF(ITYPE .EQ. IMONA) THEN
	  IF(HKL(2) .LT. 0) THEN
	    EQV(2)=-HKL(2)
	    EQV(3)=-HKL(3)
	  ELSE IF(HKL(2) .EQ. 0) THEN
	    EQV(3)=ABS(HKL(3))
	  ENDIF
C If Friedel: HKL,-HKL
	  IF(IFRIEDEL .EQ. 1) EQV(1)=ABS(HKL(1))
C
C Monoclinic(P2,b-axis): HKL,-HK-L
C
	ELSE IF(ITYPE .EQ. IMONB) THEN
	  IF(HKL(1) .LT. 0) THEN
	    EQV(1)=-HKL(1)
	    EQV(3)=-HKL(3)
	  ELSE IF(HKL(1) .EQ. 0) THEN
	    EQV(3)=ABS(HKL(3))
	  ENDIF
C If Friedel: HKL,H-KL
	  IF(IFRIEDEL .EQ. 1) EQV(2)=ABS(HKL(2))
C
C Monoclinic(P2,c-axis): HKL,-H-KL
C
	ELSE IF(ITYPE .EQ. IMONC) THEN
	  IF(HKL(1) .LT. 0) THEN
	    EQV(1)=-HKL(1)
	    EQV(2)=-HKL(2)
	  ELSE IF(HKL(1) .EQ. 0) THEN
	    EQV(2)=ABS(HKL(2))
	  ENDIF
C If Friedel: HKL,HK-L
	  IF(IFRIEDEL .EQ. 1) EQV(3)=ABS(HKL(3))
C
C Identical hkls, triclinic (P1): HKL
C
	ELSE IF(ITYPE .EQ. IDENT) THEN
	  IF(IFRIEDEL .EQ. 1) THEN
C If Friedel: HKL,-H-K-L
	    IF(EQV(1) .LT. 0) THEN
	      EQV(1)=-EQV(1)
	      EQV(2)=-EQV(2)
	      EQV(3)=-EQV(3)
	    ELSE IF(EQV(1).EQ.0 .AND. EQV(2).LT.0) THEN
	      EQV(2)=-EQV(2)
	      EQV(3)=-EQV(3)
	    ELSE IF(EQV(1).EQ.0 .AND. EQV(2).EQ.0) THEN
	      EQV(3)=ABS(EQV(3))
	    ENDIF
	  ENDIF
C
C Unknown type:
C
	ELSE
	  STOP 'BUG(calc_equiv_hkl): ITYPE invalid'
	ENDIF
C
	RETURN
	END
