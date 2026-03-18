C ========================== Routines in this file ============================
C	SUBROUTINE APPLY_REJECT_AREAS()
C	SUBROUTINE EDGE_ELLIPSE(IOUT, XCEN,YCEN,E,F,H,RMULT)
C	SUBROUTINE EDGE_RADIUS(IOUT, XCEN,YCEN,RAD)
C	SUBROUTINE REJECT_ELLIPSE(IOUT, XCEN,YCEN,E,F,H,RMULT)
C	SUBROUTINE REJECT_RADIUS(IOUT, XCEN,YCEN,RAD)
C	SUBROUTINE REJECT_POINT(IOUT, IX,IY)
C =============================================================================

	SUBROUTINE APPLY_REJECT_AREAS()
C
C NB: This routine can be a hog, so put IY in outer loops
C
	COMMON /RIMAGE_COM/ RIMAGE(8000,2500)
	COMMON /EXCLUDE_COM/ NCIRC,ICIRC(3,999),NRECT,IRECT(4,999)
C
	LOGICAL LREJECT
	COMMON /REJECT_COM/ LREJECT(8000,2500)
C
C Output the exclusion areas to the log file
C
	WRITE(30,'(/,A)') 'All spots must lie completely within the image limits and'
	WRITE(30,'(A)') 'not include any pixels in the following exclusion areas:'
	WRITE(30,'(I3,A)') NCIRC,' circular exclusion areas (Xcen,Ycen,Radius)'
	WRITE(30,'(3(4X,3(A,I4),A))') ('(',ICIRC(1,K),',',ICIRC(2,K),',',ICIRC(3,K),')',K=1,NCIRC)
	WRITE(30,'(I3,A)') NRECT,' rectangular exclusion areas (Xmin,Xmax,Ymin,Ymax)'
	WRITE(30,'(2(4X,4(A,I4),A))') ('(',IRECT(1,K),',',IRECT(2,K),',',
	1									IRECT(3,K),',',IRECT(4,K),')',K=1,NRECT)
C
C Start with all pixels not rejected
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
	DO IY=1,NUMY
	  DO IX=1,NUMX
	    LREJECT(IX,IY)=.FALSE.
	  ENDDO
	ENDDO
C
C Zero all pixels that are on or within the circular exclusions areas
C
	DO I=1,NCIRC
	  DO IY=MAX(1,ICIRC(2,I)-ICIRC(3,I)),MIN(NUMY,ICIRC(2,I)+ICIRC(3,I))
	    DO IX=MAX(1,ICIRC(1,I)-ICIRC(3,I)),MIN(NUMX,ICIRC(1,I)+ICIRC(3,I))
	      IF((IX-ICIRC(1,I))**2 + (IY-ICIRC(2,I))**2 .LE. ICIRC(3,I)**2) THEN
	        RIMAGE(IX,IY)=0.0
	        LREJECT(IX,IY)=.TRUE.
	      ENDIF
	    ENDDO
	  ENDDO
	ENDDO
C
C Zero all pixels that are on or within the rectangular exclusions areas
C
	DO I=1,NRECT
	  DO IY=MAX(1,IRECT(3,I)),MIN(NUMY,IRECT(4,I))
	    DO IX=MAX(1,IRECT(1,I)),MIN(NUMX,IRECT(2,I))
	      RIMAGE(IX,IY)=0.0
	      LREJECT(IX,IY)=.TRUE.
	    ENDDO
	  ENDDO
	ENDDO
C
C Zero the pixels beyond NUMX x NUMY
C
	DO IY=1,2500
	  DO IX=NUMX+1,8000
	    RIMAGE(IX,IY)=0.0
	  ENDDO
	ENDDO
C
	DO IY=NUMY+1,2500
	  DO IX=1,8000
	    RIMAGE(IX,IY)=0.0
	  ENDDO
	ENDDO
C
	RETURN
	END


	SUBROUTINE EDGE_ELLIPSE(IOUT, XCEN,YCEN,E,F,H,RMULT)
C
C Check if any pixels within ellipse cross the image edge
C Return IOUT =-1 if ellipse is non-physical, =1 if any pixels
C cross the edge, =0 if all pixels are valid
C
C Reject ellipses with non-physical parameters
C
	IOUT=-1
	IF(E.LE.0.0 .OR. F.LE.0.0 .OR. E*F.LE.H**2) RETURN
C
C Calculate tangents of ellipse
C
	TY=SQRT(E/(E*F-H**2))*RMULT
	TX=SQRT(F/(E*F-H**2))*RMULT
C
C Return with IOUT=1 if edge is crossed
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
	IF(XCEN+TX.GT.NUMX-0.1 .OR. XCEN-TX.LT.1.1 .OR. 
	1   YCEN+TY.GT.NUMY-0.1 .OR. YCEN-TY.LT.1.1) IOUT=1
C
	IOUT=0
	RETURN
	END


	SUBROUTINE EDGE_RADIUS(IOUT, XCEN,YCEN,RAD)
C
C Return IOUT=1 if there an edge is within RAD of XCEN,YCEN
C Return IOUT=0 otherwise
C
	COMMON /EXCLUDE_COM/ NCIRC,ICIRC(3,999),NRECT,IRECT(4,999)
C
C If XCEN,YCEN within RAD pixels of an image edge, set flag to reject
C NB: Extra 0.1 margin added to ensure never outside array limits
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
	IF(	XCEN-RAD.LT.1.1 .OR. XCEN+RAD.GT.NUMX-0.1 .OR.
	1	YCEN-RAD.LT.1.1 .OR. YCEN+RAD.GT.NUMY-0.1     ) IOUT=1
C
	IOUT=0
	RETURN
	END


	SUBROUTINE REJECT_ELLIPSE(IOUT, XCEN,YCEN,E,F,H,RMULT)
C
C Check if any pixels within ellipse are no edges, holes, etc.
C Return IOUT =-1 if ellipse is non-physical, =1 if any pixels
C are rejected, =0 if all pixels are valid
C
C Reject ellipses with non-physical parameters
C
	IOUT=-1
	IF(E.LE.0.0 .OR. F.LE.0.0 .OR. E*F.LE.H**2) RETURN
C
C Calculate tangents of ellipse
C
	TY=SQRT(E/(E*F-H**2))
	TX=SQRT(F/(E*F-H**2))
C
C If encompassing circle not rejected, then must be OK
C
	RADIUS=MAX(TX,TY)*RMULT
	CALL REJECT_RADIUS(IOUT, XCEN,YCEN,RADIUS)
	IF(IOUT .EQ. 0) RETURN
C
C Potential reject, so check each pixel in ellipse
C
	DX=TX*RMULT
	NTOT=0
	DO IX=CEILING(XCEN-DX),FLOOR(XCEN+DX)
	  DY=SQRT(MAX(0.0, (1.0 - ((IX-XCEN)/DX)**2 )/F ))*RMULT
	  Y0=YCEN-(IX-XCEN)*H/F
	  DO IY=CEILING(Y0-DY),FLOOR(Y0+DY)
	    CALL REJECT_POINT(IOUT, IX,IY)
	    IF(IOUT .EQ. 1) RETURN
	    NTOT=NTOT+1
	  ENDDO
	ENDDO
C
C Also capture ellipses with very little area
C
	IOUT=0
	IF(NTOT .LT. 10) IOUT=-1
C
	RETURN
	END


	SUBROUTINE REJECT_RADIUS(IOUT, XCEN,YCEN,RAD)
C
C Return IOUT=1 if there are rejected pixels within RAD of XCEN,YCEN
C Return IOUT=0 otherwise
C
	COMMON /EXCLUDE_COM/ NCIRC,ICIRC(3,999),NRECT,IRECT(4,999)
C
C Set the flag to accept pixel
C
	IOUT=0
C
C If XCEN,YCEN within RAD pixels of an image edge, set flag to reject
C NB: Extra 0.1 margin added to ensure never outside array limits
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
	IF(	XCEN-RAD.LT.1.1 .OR. XCEN+RAD.GT.NUMX-0.1 .OR.
	1	YCEN-RAD.LT.1.1 .OR. YCEN+RAD.GT.NUMY-0.1     ) IOUT=1
C
C If XCEN,YCEN within DIST_REJECT pixels of a circular exclusion area, set flag to reject
C
	DO I=1,NCIRC
	  IF((XCEN-ICIRC(1,I))**2 + (YCEN-ICIRC(2,I))**2 .LE. (RAD+ICIRC(3,I))**2) IOUT=1
	ENDDO
C
C If a square of size 2*DIST_REJECT pixels centered on XCEN,YCEN is within a rectangular
C exclusion area, set flag to reject
C
	DO I=1,NRECT
	  IF(	XCEN.GE.IRECT(1,I)-RAD .AND. XCEN.LE.IRECT(2,I)+RAD .AND. 
	1	    YCEN.GE.IRECT(3,I)-RAD .AND. YCEN.LE.IRECT(4,I)+RAD       ) IOUT=1
	ENDDO
C
	RETURN
	END


	SUBROUTINE REJECT_POINT(IOUT, IX,IY)
C
C Return IOUT=1 if pixel XCEN,YCEN should be rejected, return IOUT=0 otherwise
C
	LOGICAL LREJECT
	COMMON /REJECT_COM/ LREJECT(8000,2500)
C
C Check IX,IY are valid
C
	CALL GET_IMAGE_NUMXY(NUMX,NUMY)
	IF(IX.LT.1 .OR. IX.GT.NUMX .OR. IY.LT.1 .OR. IY.GT.NUMY) THEN
	  IOUT=1
	  RETURN
	ENDIF
C
C Check if IX,IY is rejected
C
	IOUT=0
	IF( LREJECT(IX,IY) ) IOUT=1
C
	RETURN
	END
