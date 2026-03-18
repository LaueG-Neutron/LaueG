C================= Routines for XYZ vectors and matrices ====================
C============= Elementary functions for (3) vectors
C	SUBROUTINE VC3_SET(VC, V1,V2,V3)
C	FUNCTION VC3_SIZE(VC)
C	SUBROUTINE VC3_UNIT(VC)
C	SUBROUTINE VC3_COPY(A,B)
C	SUBROUTINE VC3_ADD(A,B,C)
C	SUBROUTINE VC3_SUB(A,B,C)
C	SUBROUTINE VC3_SCALE(A,B,SCALE)
C	LOGICAL FUNCTION VC3_EQ(A,B)
C	FUNCTION VC3_DOT(A,B)
C	SUBROUTINE VC3_CROSS(A, B,C)
C	SUBROUTINE VC3_PRINT(NAME,VC)
C============= Rotate vector about X,Y,Z axis using angles in radians or degrees
C	SUBROUTINE VC3_XROTR(VC1, VC2,ANG)
C	SUBROUTINE VC3_YROTR(VC1, VC2,ANG)
C	SUBROUTINE VC3_ZROTR(VC1, VC2,ANG)
C	SUBROUTINE VC3_XROTD(VC1, VC2,ANG)
C	SUBROUTINE VC3_YROTD(VC1, VC2,ANG)
C	SUBROUTINE VC3_ZROTD(VC1, VC2,ANG)
C============= Elementary functions for (3,3) matrices
C	SUBROUTINE MX3_SET(MX, V11,V12,V13, V21,V22,V23, V31,V32,V33)
C	SUBROUTINE MX3_COPY(A,B)
C	SUBROUTINE MX3_TRANS(A,B)
C	SUBROUTINE MX3_ADD(A,B,C)
C	SUBROUTINE MX3_SUB(A,B,C)
C	SUBROUTINE MX3_SCALE(A,B,SCALE)
C	SUBROUTINE MX3_INVERT(AINV,DET0, A0)
C	LOGICAL FUNCTION MX3_EQ(A,B)
C	SUBROUTINE MX3_MULT(M_OUT, M1,M2)
C	SUBROUTINE MX3_PRINT(NAME,MX)
C============= Matrix for X,Y,Z rotations using angles in radians or degrees
C	SUBROUTINE MX3_XROTR(ROT, ANG)
C	SUBROUTINE MX3_YROTR(ROT, ANG)
C	SUBROUTINE MX3_ZROTR(ROT, ANG)
C	SUBROUTINE MX3_XROTD(ROT, ANG)
C	SUBROUTINE MX3_YROTD(ROT, ANG)
C	SUBROUTINE MX3_ZROTD(ROT, ANG)
C============= Matrix * Vector multiplication
C	SUBROUTINE VC3_MULT(V_OUT, M,V)
C============= Extract eigen-vectors from matrix
C	SUBROUTINE MX3_EIGEN(VALUES,VECTORS, MATRIX)
C============= Conversion between radians and degrees
C	FUNCTION ANG2DEG(ANG)
C	FUNCTION ANG2RAD(ANG)
C============= Library routines used by MX3_EIGEN()
C      SUBROUTINE EA06C(A,VALUE,VECTOR,M,IA,IV,W)
C      SUBROUTINE EA08C(A,B,VALUE,VEC,M,IV,W)
C      SUBROUTINE EA09C(A,B,VALUE,M,OFF)
C      SUBROUTINE MC04B(A,ALPHA,BETA,M,IA,Q)
C=======================================================================

C============= Elementary functions for (3) vectors
	SUBROUTINE VC3_SET(VC, V1,V2,V3)
C
	REAL VC(3)
C
	TMP1=V1
	TMP2=V2
	TMP3=V3
C
	VC(1)=TMP1
	VC(2)=TMP2
	VC(3)=TMP3
C
	RETURN
	END

	
	FUNCTION VC3_SIZE(VC)
C
	REAL VC(3)
C
	VC3_SIZE=SQRT( VC(1)**2 + VC(2)**2 + VC(3)**2 )
C
	RETURN
	END


	SUBROUTINE VC3_UNIT(VC)
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


	SUBROUTINE VC3_COPY(A,B)
C
	REAL A(3),B(3)
C
	DO I=1,3
	  A(I)=B(I)
	ENDDO
C
	RETURN
	END


	SUBROUTINE VC3_ADD(A,B,C)
C
	REAL A(3),B(3),C(3)
C
	DO I=1,3
	  A(I)=B(I)+C(I)
	ENDDO
C
	RETURN
	END


	SUBROUTINE VC3_SUB(A,B,C)
C
	REAL A(3),B(3),C(3)
C
	DO I=1,3
	  A(I)=B(I)-C(I)
	ENDDO
C
	RETURN
	END


	SUBROUTINE VC3_SCALE(A,B,SCALE)
C
	REAL A(3),B(3)
C
	DO I=1,3
	  A(I)=B(I)*SCALE
	ENDDO
C
	RETURN
	END


	LOGICAL FUNCTION VC3_EQ(A,B)
C
	REAL A(3),B(3)
C
	PARAMETER (TINY=1.0E-7)
C
	A_MAX=0.0
	DO I=1,3
	  A_MAX=MAX(A_MAX,ABS(A(I)))
	ENDDO
C
	VC3_EQ=.FALSE.
	DO I=1,3
	  DIFF=ABS(A(I)-B(I))
	  IF(DIFF .GT. 10.0*TINY*A_MAX) GOTO 100
	ENDDO
	VC3_EQ=.TRUE.
C
100	RETURN
	END

	
	FUNCTION VC3_DOT(A,B)
C
	REAL A(3),B(3)
C
	VC3_DOT=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
C
	RETURN
	END


	SUBROUTINE VC3_CROSS(A, B,C)
C
	REAL A(3),B(3),C(3)
C
	A_1=B(2)*C(3)-B(3)*C(2)
	A_2=B(3)*C(1)-B(1)*C(3)
	A_3=B(1)*C(2)-B(2)*C(1)
C
	A(1)=A_1
	A(2)=A_2
	A(3)=A_3
C
	RETURN
	END


	SUBROUTINE VC3_PRINT(NAME,VC)
C
	CHARACTER*(*) NAME
	REAL VC(3)
C
	PRINT '(1X,A,3F9.3)',NAME//'=',VC
C
	RETURN
	END

C============= Rotate vector about X,Y,Z axis using angles in radians or degrees
	SUBROUTINE VC3_XROTR(VC1, VC2,ANG)
C
C Rotate vector VC2 by ANG radians about the X axis
C Put the result in VC1
C
	REAL VC1(3),VC2(3)
C
	COS_ANG=COS(ANG)
      	SIN_ANG=SIN(ANG)   
	VC1(1)=VC2(1)
	TEMP  =COS_ANG*VC2(2) - SIN_ANG*VC2(3)
	VC1(3)=SIN_ANG*VC2(2) + COS_ANG*VC2(3)
	VC1(2)=TEMP
C
	RETURN
	END


	SUBROUTINE VC3_YROTR(VC1, VC2,ANG)
C
C Rotate vector VC2 by ANG radians about the Y axis
C Put the result in VC1
C
	REAL VC1(3),VC2(3)
C
	COS_ANG=COS(ANG)
	SIN_ANG=SIN(ANG)   
	VC1(2)=VC2(2)
	TEMP  =COS_ANG*VC2(1) - SIN_ANG*VC2(3)
	VC1(3)=SIN_ANG*VC2(1) + COS_ANG*VC2(3)
	VC1(1)=TEMP
C
	RETURN
	END


	SUBROUTINE VC3_ZROTR(VC1, VC2,ANG)
C
C Rotate vector VC2 by ANG radians about the Z axis
C Put the result in VC1
C
	REAL VC1(3),VC2(3)
C
	COS_ANG=COS(ANG)
	SIN_ANG=SIN(ANG)   
	VC1(3)=VC2(3)
	TEMP  =COS_ANG*VC2(1) - SIN_ANG*VC2(2)
	VC1(2)=SIN_ANG*VC2(1) + COS_ANG*VC2(2)
	VC1(1)=TEMP
C
	RETURN
	END


	SUBROUTINE VC3_XROTD(VC1, VC2,ANG)
	REAL VC1(3),VC2(3)
	CALL VC3_XROTR(VC1, VC2,ANG2RAD(ANG))
	RETURN
	END


	SUBROUTINE VC3_YROTD(VC1, VC2,ANG)
	REAL VC1(3),VC2(3)
	CALL VC3_YROTR(VC1, VC2,ANG2RAD(ANG))
	RETURN
	END


	SUBROUTINE VC3_ZROTD(VC1, VC2,ANG)
	REAL VC1(3),VC2(3)
	CALL VC3_ZROTR(VC1, VC2,ANG2RAD(ANG))
	RETURN
	END


C============= Elementary functions for (3,3) matrices
	SUBROUTINE MX3_SET(MX, V11,V12,V13, V21,V22,V23, V31,V32,V33)
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


		SUBROUTINE MX3_COPY(A,B)
C
	REAL A(3,3),B(3,3)
C
	DO I1=1,3
	  DO I2=1,3
	    A(I1,I2)=B(I1,I2)
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE MX3_TRANS(A,B)
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


	SUBROUTINE MX3_ADD(A,B,C)
C
	REAL A(3,3),B(3,3),C(3,3)
C
	DO I1=1,3
	  DO I2=1,3
	    A(I1,I2)=B(I1,I2)+C(I1,I2)
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE MX3_SUB(A,B,C)
C
	REAL A(3,3),B(3,3),C(3,3)
C
 	DO I1=1,3
	  DO I2=1,3
	    A(I1,I2)=B(I1,I2)-C(I1,I2)
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE MX3_SCALE(A,B,SCALE)
C
	REAL A(3,3),B(3,3)
C
	DO I1=1,3
	  DO I2=1,3
	    A(I1,I2)=B(I1,I2)*SCALE
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE MX3_INVERT(AINV,DET0, A0)
C				 
C Invert A0(3,3) and put in AINV(3,3).
C Uses classical adjoint method with double precision maths.
C
	REAL AINV(3,3),A0(3,3)
C
	REAL*8 ADJ(3,3),A(3,3),DET
C
C Copy the REAL*4 A0() into the REAL*8 A().
C
	DO I1=1,3
	  DO I2=1,3
	    A(I1,I2)=A0(I1,I2)
	  ENDDO
	ENDDO
C
C Explicitly calculate adjoint matrix.
C
	ADJ(1,1)= A(2,2)*A(3,3)-A(2,3)*A(3,2)
	ADJ(2,1)=-A(2,1)*A(3,3)+A(2,3)*A(3,1)
	ADJ(3,1)= A(2,1)*A(3,2)-A(2,2)*A(3,1)
C
	ADJ(1,2)=-A(1,2)*A(3,3)+A(1,3)*A(3,2)
	ADJ(2,2)= A(1,1)*A(3,3)-A(1,3)*A(3,1)
	ADJ(3,2)=-A(1,1)*A(3,2)+A(1,2)*A(3,1)
C
	ADJ(1,3)= A(1,2)*A(2,3)-A(1,3)*A(2,2)
	ADJ(2,3)=-A(1,1)*A(2,3)+A(1,3)*A(2,1)
	ADJ(3,3)= A(1,1)*A(2,2)-A(1,2)*A(2,1)
C
C Calculate determinant and give up on inversion if too small.
C
	DET=A(1,1)*ADJ(1,1)+A(1,2)*ADJ(2,1)+A(1,3)*ADJ(3,1)
	DET0=DET
	IF(ABS(DET) .LT. 1E-20) RETURN
C
C Trivially calculate AINV() from the adjoint and determinant.
C
	DO I1=1,3
	  DO I2=1,3
	    AINV(I1,I2)=ADJ(I1,I2)/DET
	  ENDDO
	ENDDO
C
	RETURN
	END


	LOGICAL FUNCTION MX3_EQ(A,B)
C
	REAL A(3,3),B(3,3)
C
	PARAMETER (TINY=1.0E-7)
C
	A_MAX=0.0
	DO I1=1,3
	  DO I2=1,3
	    A_MAX=MAX(A_MAX,ABS(A(I1,I2)))
	  ENDDO
	ENDDO
C
	MX3_EQ=.FALSE.
	DO I1=1,3
	  DO I2=1,3
	    DIFF=ABS(A(I1,I2)-B(I1,I2))
	    IF(DIFF .GT. 10.0*TINY*A_MAX) GOTO 100
	  ENDDO
	ENDDO
	MX3_EQ=.TRUE.
C
100	RETURN
	END


	SUBROUTINE MX3_MULT(M_OUT, M1,M2)
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


	SUBROUTINE MX3_PRINT(NAME,MX)
C
	CHARACTER*(*) NAME
	REAL MX(3,3)
C
	PRINT *,NAME//'='
	PRINT '(1X,3F9.5)',((MX(K1,K2),K2=1,3),K1=1,3)
C
	RETURN
	END

	
C============= Matrix for X,Y,Z rotations using angles in radians or degrees
	SUBROUTINE MX3_XROTR(ROT, ANG)
C
C Create rotation matrix ROT of ANG radians about the X axis.
C                                    
	REAL ROT(3,3)
C
	COS_ANG=COS(ANG)
	SIN_ANG=SIN(ANG)   
	ROT(1,1)=1.0
	ROT(1,2)=			0.0
	ROT(1,3)=						0.0
	ROT(2,1)=0.0
	ROT(2,2)=			COS_ANG
	ROT(2,3)=						-SIN_NG
	ROT(3,1)=0.0
	ROT(3,2)=			SIN_ANG
	ROT(3,3)=						COS_ANG
C
	RETURN
	END


	SUBROUTINE MX3_YROTR(ROT, ANG)
C
C Create rotation matrix ROT of ANG radians about the Y axis.
C                                    
	REAL ROT(3,3)
C
	COS_ANG=COS(ANG)
	SIN_ANG=SIN(ANG)   
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


	SUBROUTINE MX3_ZROTR(ROT, ANG)
C
C Create rotation matrix ROT of ANG radians about the Z axis.
C                                    
	REAL ROT(3,3)
C
	COS_ANG=COS(ANG)
	SIN_ANG=SIN(ANG)   
	ROT(1,1)=COS_ANG
	ROT(1,2)=			-SIN_ANG
	ROT(1,3)=						0.0
	ROT(2,1)=SIN_ANG
	ROT(2,2)=			COS_ANG
	ROT(2,3)=						0.0
	ROT(3,1)=0.0
	ROT(3,2)=			0.0
	ROT(3,3)=						1.0
C
	RETURN
	END


	SUBROUTINE MX3_XROTD(ROT, ANG)
	REAL ROT(3,3)
	CALL MX3_XROTR(ROT, ANG2RAD(ANG))
	RETURN
	END


	SUBROUTINE MX3_YROTD(ROT, ANG)
	REAL ROT(3,3)
	CALL MX3_YROTR(ROT, ANG2RAD(ANG))
	RETURN
	END


	SUBROUTINE MX3_ZROTD(ROT, ANG)
	REAL ROT(3,3)
	CALL MX3_ZROTR(ROT, ANG2RAD(ANG))
	RETURN
	END


C============= Matrix * Vector multiplication
	SUBROUTINE VC3_MULT(V_OUT, M,V)
C
	REAL V_OUT(3),M(3,3),V(3)
C
	REAL SAVE(3)
C
	DO I=1,3
	  SAVE(I)=M(I,1)*V(1)+M(I,2)*V(2)+M(I,3)*V(3)
	ENDDO
C     
	DO I=1,3
	  V_OUT(I)=SAVE(I)
	ENDDO
C
	RETURN
	END


C============= Extract eigen-vectors from matrix
	SUBROUTINE MX3_EIGEN(VALUES,VECTORS, MATRIX)
C
C Calculate eigen-values VALUES() and eigen-vectors VECTORS()
C for the 3 x 3 matrix MATRIX().
C
	REAL VALUES(3),VECTORS(3,3),MATRIX(3,3)
C
C Looks like we only need WORK() to be 3*N+1, but use 100 to be safe.
C
	REAL WORK(100),MATRIX2(3,3)
C
C Make a copy of MATRIX() so that EA06C doesn't corrupt it.
C
	DO I1=1,3
	  DO I2=1,3
	    MATRIX2(I1,I2)=MATRIX(I1,I2)
	  ENDDO
	ENDDO
C
	CALL EA06C(MATRIX2,VALUES,VECTORS,3,3,3,WORK)
	RETURN
	END


C============= Conversion between radians and degrees
	FUNCTION ANG2DEG(ANG)
	ANG2DEG=ANG*ATAND(1.0)/ATAN(1.0)
	RETURN
	END


	FUNCTION ANG2RAD(ANG)
	ANG2RAD=ANG*ATAN(1.0)/ATAND(1.0)
	RETURN
	END

	
C============= Library routines used by MX3_EIGEN()
C-----------------------------------------------------------------------
      SUBROUTINE EA06C(A,VALUE,VECTOR,M,IA,IV,W)
C-----------------------------------------------------------------------
C###### 18/05/70 LAST LIBRARY UPDATE
      REAL A(IA,M),VALUE(M),VECTOR(IV,M),W(1)
      M1=M+1
      W(1)=A(1,1)
      IF(M-2)30,10,20
   10 W(2)=A(2,2)
      W(4)=A(2,1)
      GOTO 30
   20 CALL MC04B(A,W,W(M1),M,IA,W(M+M1))
   30 CALL EA08C(W,W(M1),VALUE,VECTOR,M,IV,W(M+M1))
      IF(M.LE.2)RETURN
      DO 70 L=1,M
      DO 70 II=3,M
      I=M-II+1
      IF(W(M1+I))40,70,40
   40 PP=0.
      I1=I+1
      DO 50 K=I1,M
   50 PP=PP+A(I,K)*VECTOR(K,L)
      PP=PP/(A(I,I+1)*W(M1+I))
      DO 60 K=I1,M
   60 VECTOR(K,L)=VECTOR(K,L)+PP*A(I,K)
   70 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE EA08C(A,B,VALUE,VEC,M,IV,W)
C-----------------------------------------------------------------------
C###### 15/06/70 LAST LIBRARY UPDATE
      REAL A(M),B(M),VALUE(M),VEC(1),W(1)
      DATA A34/0./,EPS/1.E-6/
C
C---- THIS USES QR ITERATION TO FIND THE EIGENVALUES AND EIGENVECTORS
C---- OF THE SYMMETRIC TRIDIAGONAL MATRIX WHOSE DIAGONAL ELEMENTS ARE
C---- A(I),I=1,M AND OFF-DIAGONAL ELEMENTS ARE B(I),I=2,M.  THE ARRAY
C---- W IS USED FOR WORKSPACE AND MUST HAVE DIMENSION AT LEAST 2*M.
C---- WE TREAT VEC AS IF IT HAD DIMENSIONS (IV,M).
C
      CALL EA09C(A,B,W(M+1),M,W)
C
C---- SET VEC TO THE IDENTITY MATRIX.
C
      DO 20 I=1,M
      VALUE(I)=A(I)
      W(I)=B(I)
      K=(I-1)*IV+1
      L=K+M-1
      DO 10 J=K,L
   10 VEC(J)=0.
   20 VEC(K+I-1)=1.
      ITER=0
      IF(M.EQ.1)RETURN
      N2=M
C
C---- EACH QR ITERATION IS PERFORMED OF ROWS AND COLUMNS N1 TO N2
C
   30 DO 40 II=2,N2
      N1=2+N2-II
      IF(ABS(W(N1)).LE.(ABS(VALUE(N1-1))+ABS(VALUE(N1)))*EPS)GOTO 50
   40 CONTINUE
      N1=1
   50 IF(N2.NE.N1)GOTO 60
      N2=N2-1
      IF(N2-1)90,90,30
   60 ROOT=W(M+N2)
      ITER=ITER+1
      N2M1=N2-1
      A22=VALUE(N1)
      A12=A22-ROOT
      A23=W(N1+1)
      A13=A23
      DO 80 I=N1,N2M1
      A33=VALUE(I+1)
      IF(I.NE.N2M1)A34=W  (I+2)
      S=SIGN(SQRT(A12*A12+A13*A13),A12)
      SI=A13/S
      CO=A12/S
      JK=I*IV+1
      J1=JK-IV
      J2=J1+MIN0(M,I+ITER)-1
      DO 70 JI=J1,J2
      V1=VEC(JI)
      V2=VEC(JK)
      VEC(JI)=V1*CO+V2*SI
      VEC(JK)=V2*CO-V1*SI
   70 JK=JK+1
      IF(I.NE.N1)  W(I)=S
      A11=CO*A22+SI*A23
      A12=CO*A23+SI*A33
      A13=SI*A34
      A21=CO*A23-SI*A22
      A22=CO*A33-SI*A23
      A23=CO*A34
      VALUE(I)=A11*CO+A12*SI
      A12=-A11*SI+A12*CO
      W(I+1)=A12
   80 A22=A22*CO-A21*SI
      VALUE(N2)=A22
      GOTO 30
C
C---- RAYLEIGH QUOTIENT
C
   90 DO 110 J=1,M
      K=(J-1)*IV
      XX=VEC(K+1)**2
      XAX=XX*A(1)
      DO 100 I=2,M
      XX=XX+VEC(K+I)**2
  100 XAX=XAX+VEC(K+I)*(2.*B(I)*VEC(K+I-1)+A(I)*VEC(K+I))
  110 VALUE(J)=XAX/XX
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE EA09C(A,B,VALUE,M,OFF)
C-----------------------------------------------------------------------
C###### 15/06/70 LAST LIBRARY UPDATE
      REAL A(M),B(M),VALUE(M),OFF(M)
      DATA A34/0./,EPS/0.6E-7/
      VALUE(1)=A(1)
      IF(M.EQ.1)RETURN
      DO 10 I=2,M
      VALUE(I)=A(I)
   10 OFF(I)=B(I)
C
C---- EACH QR ITERATION IS PERFORMED OF ROWS AND COLUMNS N1 TO N2
C
      N2=M
   20 IF(N2.LE.1)RETURN
      DO 30 II=2,N2
      N1=2+N2-II
      IF(ABS(OFF(N1)).LE.(ABS(VALUE(N1-1))+ABS(VALUE(N1)))*EPS)GOTO 40
   30 CONTINUE
      N1=1
   40 IF(N2.NE.N1)GOTO 50
      N2=N2-1
      GOTO 20
C
C---- ROOT  IS THE EIGENVALUE OF THE BOTTOM 2*2 MATRIX THAT IS NEAREST
C---- TO THE LAST MATRIX ELEMENT AND IS USED TO ACCELERATE THE
C---- CONVERGENCE
C
   50 BB=(VALUE(N2)-VALUE(N2-1))*0.5
      CC=OFF(N2)*OFF(N2)
      SBB=1.
      IF(BB.LT.0.)SBB=-1.
      ROOT=VALUE(N2)+CC/(BB+SBB*SQRT(BB*BB+CC))
      N2M1=N2-1
   60 A22=VALUE(N1)
      A12=A22-ROOT
      A23=OFF(N1+1)
      A13=A23
      DO 70 I=N1,N2M1
      A33=VALUE(I+1)
      IF(I.NE.N2M1)A34=OFF(I+2)
      S=SQRT(A12*A12+A13*A13)
      SI=A13/S
      CO=A12/S
      IF(I.NE.N1)OFF(I)=S
      A11=CO*A22+SI*A23
      A12=CO*A23+SI*A33
      A13=SI*A34
      A21=CO*A23-SI*A22
      A22=CO*A33-SI*A23
      A23=CO*A34
      VALUE(I)=A11*CO+A12*SI
      A12=-A11*SI+A12*CO
      OFF(I+1)=A12
   70 A22=A22*CO-A21*SI
      VALUE(N2)=A22
      GOTO 20
      END
C-----------------------------------------------------------------------
      SUBROUTINE MC04B(A,ALPHA,BETA,M,IA,Q)
C-----------------------------------------------------------------------
C###### 27/03/72 LAST LIBRARY UPDATE
      DIMENSION A(IA,1),ALPHA(1),BETA(1),Q(1)
      ALPHA(1)=A(1,1)
      DO 20 J=2,M
      J1=J-1
      DO 10 I=1,J1
      A(I,J)=A(J,I)
   10 CONTINUE
      ALPHA(J)=A(J,J)
   20 CONTINUE
      M1=M-1
      M2=M-2
      DO 170 I=1,M2
      PP=0.0
      I1=I+1
      DO 30 J=I1,M
      PP=PP+A(I,J)**2
   30 CONTINUE
      PP1=SQRT(PP)
      IF(A(I,I+1))50,40,40
   40 BETA(I+1)=-PP1
      GOTO 60
   50 BETA(I+1)=PP1
   60 IF(PP)170,170,70
   70 H=PP-BETA(I+1)*A(I,I+1)
      A(I,I+1)=A(I,I+1)-BETA(I+1)
      DO 120 KI=I1,M
      QJ=0.0
      DO 80 KJ=I1,KI
      QJ=QJ+A(KJ,KI)*A(I,KJ)
   80 CONTINUE
      IF(KI-M)90,110,110
   90 I2=KI+1
      DO 100 KJ=I2,M
      QJ=QJ+A(KI,KJ)*A(I,KJ)
  100 CONTINUE
  110 Q(KI)=QJ/H
  120 CONTINUE
      BIGK=0.0
      DO 130 KJ=I1,M
      BIGK=BIGK+A(I,KJ)*Q(KJ)
  130 CONTINUE
      BIGK=BIGK/(2.0*H)
      DO 140 KJ=I1,M
      Q(KJ)=Q(KJ)-BIGK*A(I,KJ)
  140 CONTINUE
      DO 160 KI=I1,M
      DO 150 KJ=KI,M
      A(KI,KJ)=A(KI,KJ)-Q(KI)*A(I,KJ)-Q(KJ)*A(I,KI)
  150 CONTINUE
  160 CONTINUE
  170 CONTINUE
      DO 180 I=2,M
      H=ALPHA(I)
      ALPHA(I)=A(I,I)
      A(I,I)=H
  180 CONTINUE
      BETA(M)=A(M-1,M)
      RETURN
      END
