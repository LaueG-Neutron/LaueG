C Commons used to define instrument geometry:
C   ITYPE: 1=KOALA1, 2=IMAGINE, 3=CYCLOPS ,4=KOALA2
C	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
C	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
C	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C	COMMON /CELL_COM/ CELL(6),ILATT,ICEN,CELL_PROD
C	COMMON /U_COM/ U_BASE(3,3),RXYZ(3)

C Sections of this file:
C	SUBROUTINE CALC_PIXELS_TO_HVECS(HVECS, PIXELS,NHKLS)
C	SUBROUTINE CALC_HVECS_TO_PIXELS(PIXELS, HVECS,NHKLS)
C	SUBROUTINE CALC_HKLS_TO_PIXWAVS(PIXWAVS, UB,HKLS,NHKLS)
C	SUBROUTINE CALC_HVECS_TO_PIXWAVS(PIXWAVS, HVECS,NVECS)
C	SUBROUTINE CALC_PIXWAVS_TO_HKLS(HKLS, UB,PIXWAVS,NHKLS)
C	SUBROUTINE CALC_PIXWAVS_TO_HVECS(HVECS, PIXWAVS,NVECS)
C================== GEOMETRY DEFINITIONS ============================
C====================== GEOMETRIC RELATIONSHIPS ======================
C===================== PARAMETER CORRELATIONS ======================
C====================================================================
C	SUBROUTINE CALC_UB_TO_CELL(CELL, UB)
C	SUBROUTINE CALC_UB_MATRIX(UB)
C	SUBROUTINE CALC_ZYXROT_MX(ROT, ANGX,ANGY,ANGZ)
C	SUBROUTINE CALC_XYZROT_MX(ROT, ANGX,ANGY,ANGZ)
C	SUBROUTINE CALC_XROT_MX(ROT, ANG)
C	SUBROUTINE CALC_YROT_MX(ROT, ANG)
C	SUBROUTINE CALC_ZROT_MX(ROT, ANG)
C	SUBROUTINE ALIGN_VECTORS(U_ROT,H_PEAK,H_HKL,NPEAKS)
C	SUBROUTINE CALC_U_FROM_VECTORS(U, H_HKLS,H_PEAKS, I_HKLS,I_PEAKS)
C	SUBROUTINE MAKE_U_ORTHOGONAL(U)
C	SUBROUTINE CALC_B_MATRIX(B, DIRECT)
C	SUBROUTINE CONV_DIRECT_RECIP(RECIP, DIRECT)


	SUBROUTINE CALC_PIXELS_TO_HVECS(HVECS, PIXELS,NHKLS)
C
C NB: Assumes the UB matrix has not been corrected for PHI 
C
	REAL HVECS(3,*),PIXELS(2,*)
C
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
	COMMON /GEOM_COM/ DRUM_RAD,PHI,NUMX,NUMY
C
	REAL ROT_PHI(3,3),INV_ROT(3,3)
	REAL XTAL_LAB(3),OFFSET(3)
C
	S0Y=SIND(BEAM_VERT)
	S0Z=COSD(BEAM_VERT)
C
	CALL CALC_YROT_MX(ROT_PHI, -PHI)
	CALL VC3_MULT(XTAL_LAB, ROT_PHI,XTAL_OFF)
	CALL VC3_SUB(OFFSET, DRUM_OFF,XTAL_LAB)
	CALL VC3_SCALE(OFFSET, OFFSET,1.0/DRUM_RAD)
C
	CALL CALC_YROT_MX(INV_ROT, +PHI)
	CALL MX3_INVERT(INV_ROT,DET, ROT_PHI)
C
	ANGLE0=( XTAL_LAB(1) - PIX_SIZE(1)*PIX_CEN(1) -
	1			PIX_SKEW*PIX_SIZE(2)*PIX_CEN(2) )/DRUM_RAD
	FAC1=PIX_SIZE(1)/DRUM_RAD
	FAC2=PIX_SKEW*PIX_SIZE(2)/DRUM_RAD
	FAC3=PIX_SIZE(2)/DRUM_RAD
C
	DO I=1,NHKLS
	  ANGLE=ANGLE0 + FAC1*PIXELS(1,I) + FAC2*PIXELS(2,I)
	  S1X=OFFSET(1) + SIN(ANGLE)
	  S1Y=OFFSET(2) + FAC3*(PIXELS(2,I)-PIX_CEN(2))
	  S1Z=OFFSET(3) + COS(ANGLE)
	  SIZE=SQRT( S1X*S1X + S1Y*S1Y + S1Z*S1Z )
	  S1X= S1X/SIZE
	  S1Y=(S1Y/SIZE -S0Y)
	  S1Z=(S1Z/SIZE -S0Z)
	  SIZE=SQRT( S1X*S1X + S1Y*S1Y + S1Z*S1Z )
	  HVECS(1,I)=(INV_ROT(1,1)*S1X + INV_ROT(1,2)*S1Y + INV_ROT(1,3)*S1Z)/SIZE
	  HVECS(2,I)=(INV_ROT(2,1)*S1X + INV_ROT(2,2)*S1Y + INV_ROT(2,3)*S1Z)/SIZE
	  HVECS(3,I)=(INV_ROT(3,1)*S1X + INV_ROT(3,2)*S1Y + INV_ROT(3,3)*S1Z)/SIZE
	ENDDO
C
	RETURN
	END


	SUBROUTINE CALC_HVECS_TO_PIXELS(PIXELS, HVECS,NHKLS)
C
C NB: Assumes the UB matrix has not been corrected for PHI 
C
	REAL HVECS(3,*),PIXELS(2,*)
C
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
	COMMON /GEOM_COM/ DRUM_RAD,PHI,NUMX,NUMY
C
	REAL ROT_PHI(3,3)
	REAL XTAL_LAB(3),OFFSET(3),H(3)
C
	S0Y=SIND(BEAM_VERT)
	S0Z=COSD(BEAM_VERT)
C
	CALL CALC_YROT_MX(ROT_PHI, -PHI)
	CALL VC3_MULT(XTAL_LAB, ROT_PHI,XTAL_OFF)
	CALL VC3_SUB(OFFSET, DRUM_OFF,XTAL_LAB)
C
	PIX_XCEN=PIX_CEN(1)-XTAL_LAB(1)/PIX_SIZE(1)
	SKEW=PIX_SKEW/PIX_SIZE(1)
	FAC1=DRUM_RAD/PIX_SIZE(1)
C
	DO I=1,NHKLS
	  H(1)=ROT_PHI(1,1)*HVECS(1,I)+ROT_PHI(1,2)*HVECS(2,I)+ROT_PHI(1,3)*HVECS(3,I)
	  H(2)=ROT_PHI(2,1)*HVECS(1,I)+ROT_PHI(2,2)*HVECS(2,I)+ROT_PHI(2,3)*HVECS(3,I)
	  H(3)=ROT_PHI(3,1)*HVECS(1,I)+ROT_PHI(3,2)*HVECS(2,I)+ROT_PHI(3,3)*HVECS(3,I)
	  SIZE2=(H(1)**2 + H(2)**2 + H(3)**2)
	  WAV=-2.0*(S0Y*H(2) + S0Z*H(3))/SIZE2
	  S1X=H(1)*WAV
	  S1Y=H(2)*WAV +S0Y
	  S1Z=H(3)*WAV +S0Z
	  OBLIQ=MAX(1E-4, 1.0-S1Y*S1Y )	! stops crashes for impossibly oblique spots
	  DIST=( DRUM_RAD*SQRT(OBLIQ) + OFFSET(1)*S1X + OFFSET(3)*S1Z )/OBLIQ
	  DX=DIST*S1X-OFFSET(1)
	  DY=DIST*S1Y-OFFSET(2)
	  DZ=DIST*S1Z-OFFSET(3)
	  PIXELS(1,I)=PIX_XCEN + FAC1*ATAN2(DX,DZ) -SKEW*DY
	  PIXELS(2,I)=PIX_CEN(2) + DY/PIX_SIZE(2)
	ENDDO
C
	RETURN
	END


	SUBROUTINE CALC_HKLS_TO_PIXWAVS(PIXWAVS, UB,HKLS,NHKLS)
C
	PARAMETER (NMAX_HKLS=300000)
C
	REAL HKLS(3,*),PIXWAVS(3,*),UB(3,3)
C
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	REAL UB2(3,3),ROT_PHI(3,3),HVECS(3,NMAX_HKLS)
C
C Create the UB adjusted for PHI
C
	CALL CALC_YROT_MX(ROT_PHI, -PHI)
	CALL MX3_MULT(UB2, ROT_PHI,UB)
C
C Calculate the scatterings vectors for the reflections HKLS(*,I)
C
	DO I=1,NHKLS
	  HVECS(1,I)=UB2(1,1)*HKLS(1,I) +UB2(1,2)*HKLS(2,I) +UB2(1,3)*HKLS(3,I)
	  HVECS(2,I)=UB2(2,1)*HKLS(1,I) +UB2(2,2)*HKLS(2,I) +UB2(2,3)*HKLS(3,I)
	  HVECS(3,I)=UB2(3,1)*HKLS(1,I) +UB2(3,2)*HKLS(2,I) +UB2(3,3)*HKLS(3,I)
	ENDDO
C
	CALL CALC_HVECS_TO_PIXWAVS(PIXWAVS, HVECS,NHKLS)
C
	RETURN
	END


	SUBROUTINE CALC_HVECS_TO_PIXWAVS(PIXWAVS, HVECS,NVECS)
C
	REAL PIXWAVS(3,*),HVECS(3,*)
C
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	REAL OFFSET(3)
C
	S0Y=SIND(BEAM_VERT)
	S0Z=COSD(BEAM_VERT)
C
	OFFSET(1)=COSD(-PHI)*XTAL_OFF(1) -SIND(-PHI)*XTAL_OFF(3) -DRUM_OFF(1)
	OFFSET(2)=XTAL_OFF(2)                                    -DRUM_OFF(2)
	OFFSET(3)=SIND(-PHI)*XTAL_OFF(1) +COSD(-PHI)*XTAL_OFF(3) -DRUM_OFF(3)
C An empirical correction that reduces correlations between XCEN & X offsets
	XTAL_OFFX=COSD(-PHI)*XTAL_OFF(1) -SIND(-PHI)*XTAL_OFF(3)
C
	DO I=1,NVECS
C Calculate wavelength for diffraction
	  SIZE2=HVECS(1,I)**2 + HVECS(2,I)**2 + HVECS(3,I)**2
	  WAV=-2.0*(S0Y*HVECS(2,I) + S0Z*HVECS(3,I))/MAX(1E-10,SIZE2)
C Calculate the (unit vector) scattering beam vector
	  S1X=HVECS(1,I)*WAV
	  S1Y=HVECS(2,I)*WAV+S0Y
	  S1Z=HVECS(3,I)*WAV+S0Z
C
	  IF(ITYPE .NE. 3) THEN
C=========== The calculation for cylindrical detector
C Calculate the distance from the sample to the spot on the detector
	    S1_RADIAL=SQRT(S1X**2+S1Z**2)
	    DIST=DRUM_RAD/S1_RADIAL -(S1X*OFFSET(1)+S1Z*OFFSET(3))/S1_RADIAL**2
C Calculate the position of the spot in laboratory coords
	    XLAB=DIST*S1X+OFFSET(1)
	    YLAB=DIST*S1Y+OFFSET(2)
	    ZLAB=DIST*S1Z+OFFSET(3)
C Convert from laboratory (3D) to plate (2D) coords
C NB: must use the radian verions of ATAN
	    DX=DRUM_RAD*ATAN2(XLAB,ZLAB)
	    DY=YLAB
C Adjust X for pixel skewness and XTAL_OFFX (see above)
	    DX=DX - DY*PIX_SKEW + XTAL_OFFX
C Convert to pixels and add the pixel centers of the plate
	    XPIX=DX/PIX_SIZE(1) +PIX_CEN(1)
	    YPIX=DY/PIX_SIZE(2) +PIX_CEN(2)
	  ELSE
C=========== The calculation for octagonal detector
C Calculate which detector-plate the spot is on
	    IPLATE=5+NINT( ATAN2D(S1X,S1Z)/45.0 )
	    IF(IPLATE .EQ. 9) IPLATE=1
	    YROT=45.0*(IPLATE-5)
C Rotate S1 from laboratory to detector-plate coords
	    TMP=COSD(YROT)*S1X -SIND(YROT)*S1Z
	    S1Z=SIND(YROT)*S1X +COSD(YROT)*S1Z
	    S1X=TMP
C Rotate OFFSET() from laboratory to detector-plate coords
	    OFFSET_X=COSD(YROT)*OFFSET(1) -SIND(YROT)*OFFSET(3)
	    OFFSET_Z=SIND(YROT)*OFFSET(1) +COSD(YROT)*OFFSET(3)
C Calculate the position of the spot in detector-plate coords
	    DIST=(DRUM_RAD-OFFSET_Z)/S1Z
	    DX=DIST*S1X+OFFSET_X
	    DY=DIST*S1Y+OFFSET(2)
C Adjust X for pixel skewness and XTAL_OFFX (see above)
	    DX=DX - DY*PIX_SKEW + XTAL_OFFX
C Convert to pixels and add the pixel centers of plate 5
	    XPIX=DX/PIX_SIZE(1) +PIX_CEN(1) +(IPLATE-5)*NUMX/8.0
	    YPIX=DY/PIX_SIZE(2) +PIX_CEN(2)
	  ENDIF
C Load x,y,wav into PIXWAVS()
	  PIXWAVS(1,I)=XPIX
	  PIXWAVS(2,I)=YPIX
	  PIXWAVS(3,I)=WAV
C
	ENDDO
C
	RETURN
	END


	SUBROUTINE CALC_PIXWAVS_TO_HKLS(HKLS, UB,PIXWAVS,NHKLS)
C
	REAL HKLS(3,*),PIXWAVS(3,*),UB(3,3)
C
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	REAL UB_ROT(3,3),ROT_PHI(3,3),UB_INV(3,3),HVECS(3,1000)
C
	CALL CALC_PIXWAVS_TO_HVECS(HVECS, PIXWAVS,NHKLS)
C
C Create the UB adjusted for PHI, then invert it
C
	CALL CALC_YROT_MX(ROT_PHI, -PHI)
	CALL MX3_MULT(UB_ROT, ROT_PHI,UB)
	CALL MX3_INVERT(UB_INV,DET, UB_ROT)
C
C Calculate the scatterings vectors for the reflections HKLS(*,I)
C
	DO I=1,NHKLS
	  HKLS(1,I)=UB_INV(1,1)*HVECS(1,I) +
	1		UB_INV(1,2)*HVECS(2,I) + UB_INV(1,3)*HVECS(3,I)
	  HKLS(2,I)=UB_INV(2,1)*HVECS(1,I) +
	1		UB_INV(2,2)*HVECS(2,I) + UB_INV(2,3)*HVECS(3,I)
	  HKLS(3,I)=UB_INV(3,1)*HVECS(1,I) +
	1		UB_INV(3,2)*HVECS(2,I) + UB_INV(3,3)*HVECS(3,I)
	ENDDO
C
	RETURN
	END


	SUBROUTINE CALC_PIXWAVS_TO_HVECS(HVECS, PIXWAVS,NVECS)
C
C
	REAL HVECS(3,*),PIXWAVS(3,*)
C
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	REAL OFFSET(3)
C
	S0Y=SIND(BEAM_VERT)
	S0Z=COSD(BEAM_VERT)
C
	OFFSET(1)=COSD(-PHI)*XTAL_OFF(1) -SIND(-PHI)*XTAL_OFF(3) -DRUM_OFF(1)
	OFFSET(2)=XTAL_OFF(2)                                    -DRUM_OFF(2)
	OFFSET(3)=SIND(-PHI)*XTAL_OFF(1) +COSD(-PHI)*XTAL_OFF(3) -DRUM_OFF(3)
C An empirical correction that reduces correlations between XCEN & X offsets
	XTAL_OFFX=COSD(-PHI)*XTAL_OFF(1) -SIND(-PHI)*XTAL_OFF(3)
C
	DO I=1,NVECS
	  IF(ITYPE .NE. 3) THEN
C=========== The calculation for cylindrical detector
	    DX=(PIXWAVS(1,I)-PIX_CEN(1))*PIX_SIZE(1)
	    DY=(PIXWAVS(2,I)-PIX_CEN(2))*PIX_SIZE(2)
	    ANGLE=( DX + DY*PIX_SKEW - XTAL_OFFX )/DRUM_RAD
	    S1X=DRUM_RAD*SIN(ANGLE) -OFFSET(1)
	    S1Y=DY                  -OFFSET(2)
	    S1Z=DRUM_RAD*COS(ANGLE) -OFFSET(3)
	  ELSE
C=========== The calculation for octagonal detector
C Calculate which detector-plate the spot is on
	    IPLATE=1+IFIX( PIXWAVS(1,I) / (NUMX/8.0) )
C Calculate spot position (in mm) relative to the center of plate 5
	    DX=(PIXWAVS(1,I) -PIX_CEN(1) -(IPLATE-5)*NUMX/8.0)*PIX_SIZE(1)
	    DY=(PIXWAVS(2,I) -PIX_CEN(2)		            )*PIX_SIZE(2)
C Adjust X for pixel skewness and XTAL_OFFX (see above)
	    DX=DX + DY*PIX_SKEW - XTAL_OFFX
C Rotate OFFSET() from laboratory to detector-plate coords
	    YROT=45.0*(IPLATE-5)
	    OFFSET_X=COSD(YROT)*OFFSET(1) -SIND(YROT)*OFFSET(3)
	    OFFSET_Z=SIND(YROT)*OFFSET(1) +COSD(YROT)*OFFSET(3)
C Calculate S1 the vector from the sample to the spot
	    S1X=DX-OFFSET_X
	    S1Y=DY-OFFSET(2)
	    S1Z=DRUM_RAD-OFFSET_Z
C Rotate S1 from the detector-plate to laboratory coords
	    TMP=COSD(-YROT)*S1X -SIND(-YROT)*S1Z
	    S1Z=SIND(-YROT)*S1X +COSD(-YROT)*S1Z
	    S1X=TMP
	  ENDIF
C Make S1 a unit vector and calculate HVEC from S0, S1, and wavelength
	  SIZE=SQRT(S1X**2 + S1Y**2 + S1Z**2)
	  HVECS(1,I)= S1X/SIZE      /PIXWAVS(3,I)
	  HVECS(2,I)=(S1Y/SIZE -S0Y)/PIXWAVS(3,I)
	  HVECS(3,I)=(S1Z/SIZE -S0Z)/PIXWAVS(3,I)
	ENDDO
C
	RETURN
	END


C================== GEOMETRY DEFINITIONS ============================
C
C We fix nominal XYZ axes as: left (looking downstream), up, downstream.
C Unspecified positions in mm, all angles anti-clockwise in degrees.
C
C	BEAM_ANGLE			Angle of incident beam from horizontal
C	XTAL_OFFSET(3)		Position offset of crystal center
C	DRUM_OFFSET(3)		Position offset of IP drum axis
C	DRUM_RADIUS			Radius of IP surface inside the drum
C	PIXEL_SIZE(2)			X & Y pixel size in microns
C	PIXEL_SKEW			X/Y pixel skewness factor
C	PIXEL_CEN(2)			Position of TTH=0 in pixels
C	PHI					Spindle rotation angle
C	CELL(6)				Sample cell parameters
C	UB_RXYZ(3)			Rotation angles to correct UB matrix
C	U_BASE(3,3)			Approximate U matrix fixed at program start
C
C	X_PIXEL,Y_PIXEL		Position on the image plate in terms of pixels
C
C====================== GEOMETRIC RELATIONSHIPS ======================
C
C Position of crystal center in lab coords
C We ignore an offset of the PHI axis as DRUM_OFFSET has the same effect
C
C	XTAL_LAB = ROT_Y(-PHI) * XTAL_OFFSET
C
C Rotation of UB matrix into lab coords
C
C	UB_LAB = ROT_Y(-PHI) * ROT_XYZ(UB_RXYZ) * U_BASE * B
C
C Calculating scattering vector for a given HKL
C
C	H_VEC = UB_LAB * HKL
C
C Calculate diffracted beam for a scattering vector
C
C	S1_VEC = S0_VEC + WAVELENGTH * H_VEC
C		for some value of WAVELENGTH such that |S1_VEC| = 1
C
C Connection between sample position, diffracted beam and pixel position
C
C	PIXEL_LAB = XTAL_LAB + DISTANCE * S1_VEC
C		for some value of DISTANCE
C
C Direction of incident beam
C
C	S0_VEC = [0,SIN(BEAM_ANGLE),COS(BEAM_ANGLE)]
C
C Position of a pixel on IP relative to TTH=0
C
C	PIXEL_IP(2) = PIXEL_SIZE(2) * ( Y_PIXEL - PIXEL_CEN(2) )
C	PIXEL_IP(1) = PIXEL_SIZE(1) * ( X_PIXEL - PIXEL_CEN(1) )
C					+ PIXEL_SKEW * PIXEL_IP(2)
C
C Position of a point on IP relative to TTH=0 in drum coords
C
C	PIXEL_DRUM = [0, PIXEL_IP(2),0]
C					+ DRUM_RADIUS * [SIN(ANGLE),0,COS(ANGLE)]
C		where ANGLE = 180/PI * PIXEL_IP(1) / DRUM_RADIUS
C
C Position of point in drum coords to lab coords
C
C	PIXEL_LAB = DRUM_OFFSET + PIXEL_DRUM
C
C The above geometry assumes the PIXEL_CEN is defined by a TTH=0
C diffracted beam (which can come from an offset sample) instead of
C the incident beam. The following corrects this effect.
C	X_PIXEL(corrected) = X_PIXEL(above eqns) + XTAL_LAB[1] / PIXEL_SIZE
C
C===================== PARAMETER CORRELATIONS ======================
C 100% correlations:
C	PIXEL_SIZE(1-2) & DRUM_RADIUS
C	PIXEL_CEN(2) & XTAL_OFFSET(2) & DRUM_OFFSET(2)
C
C 100% correlation if only one PHI used:
C	XTAL_OFFSET(1) & DRUM_OFFSET(1)
C	XTAL_OFFSET(3) & DRUM_OFFSET(3)
C
C Strong correlations when limited range in X
C	PIXEL_SIZE(1) & XTAL_OFFSET(3)
C	PIXEL_CEN(2) & BEAM_ANGLE
C	PIXEL_CEN(1) & XTAL_OFFSET(1)
C
C High correlation, but better fit if used
C	RXYZ_ANGLES(2),PIXEL_CEN(1),XTAL_OFFSET(1)
C
C====================================================================


	SUBROUTINE CALC_UB_TO_CELL(CELL, UB)
C
	REAL CELL(6),UB(3,3)
C
	REAL RECIP(6),VECH(3),VECK(3),VECL(3)
C
	CALL VC3_COPY(VECH, UB(1,1))
	CALL VC3_COPY(VECK, UB(1,2))
	CALL VC3_COPY(VECL, UB(1,3))
C
	RECIP(1)=SQRT( VC3_DOT(VECH,VECH) )
	RECIP(2)=SQRT( VC3_DOT(VECK,VECK) )
	RECIP(3)=SQRT( VC3_DOT(VECL,VECL) )
	RECIP(4)=ACOSD( VC3_DOT(VECK,VECL)/RECIP(2)/RECIP(3) )
	RECIP(5)=ACOSD( VC3_DOT(VECH,VECL)/RECIP(1)/RECIP(3) )
	RECIP(6)=ACOSD( VC3_DOT(VECH,VECK)/RECIP(1)/RECIP(2) )
C
	CALL CONV_DIRECT_RECIP(CELL, RECIP)
C
	RETURN
	END


	SUBROUTINE CALC_UB_MATRIX(UB)
C
C Rotation of UB matrix into lab coords
C
	REAL UB(3,3)
C
	COMMON /CELL_COM/ CELL(6),ILATT,ICEN,CELL_PROD
	COMMON /U_COM/ U_BASE(3,3),RXYZ(3)
C
	REAL RX(3,3),RY(3,3),RZ(3,3),B(3,3)
C
C	UB = RX * RY * RZ * U_BASE * B
C
	CALL CALC_B_MATRIX(B, CELL)
	CALL MX3_MULT(UB, U_BASE,B)
	CALL CALC_XROT_MX(RX, RXYZ(1))
	CALL CALC_YROT_MX(RY, RXYZ(2))
	CALL CALC_ZROT_MX(RZ, RXYZ(3))
	CALL MX3_MULT(UB, RZ,UB)
	CALL MX3_MULT(UB, RY,UB)
	CALL MX3_MULT(UB, RX,UB)
C
	RETURN
	END


	SUBROUTINE CALC_ZYXROT_MX(ROT, ANGX,ANGY,ANGZ)
C
	REAL ROT(3,3)
C
	REAL XROT(3,3),YROT(3,3),ZROT(3,3)
C
	CALL CALC_XROT_MX(XROT, ANGX)
	CALL CALC_YROT_MX(YROT, ANGY)
	CALL CALC_ZROT_MX(ZROT, ANGZ)
C
	CALL MX3_MULT(ROT, ZROT,YROT)
	CALL MX3_MULT(ROT, ROT, XROT)
C
	RETURN
	END

	
	SUBROUTINE CALC_XYZROT_MX(ROT, ANGX,ANGY,ANGZ)
C
	REAL ROT(3,3)
C
	REAL XROT(3,3),YROT(3,3),ZROT(3,3)
C
	CALL CALC_XROT_MX(XROT, ANGX)
	CALL CALC_YROT_MX(YROT, ANGY)
	CALL CALC_ZROT_MX(ZROT, ANGZ)
C
	CALL MX3_MULT(ROT, XROT,YROT)
	CALL MX3_MULT(ROT, ROT, ZROT)
C
	RETURN
	END


      SUBROUTINE CALC_XROT_MX(ROT, ANG)
C
C Create matrix ROT for a rotation of ANG around the X axis.
C                                    
      REAL ROT(3,3)
C
      COS_ANG=COSD(ANG)
      SIN_ANG=SIND(ANG)   
      ROT(1,1)=1.0
      ROT(1,2)=			0.0
      ROT(1,3)=						0.0
      ROT(2,1)=0.0
      ROT(2,2)=			COS_ANG
      ROT(2,3)=						-SIN_ANG
      ROT(3,1)=0.0
      ROT(3,2)=			SIN_ANG
      ROT(3,3)=						COS_ANG
C
      RETURN
      END


      SUBROUTINE CALC_YROT_MX(ROT, ANG)
C
C Create matrix ROT for a rotation of ANG around the Y axis.
C                                    
      REAL ROT(3,3)
C
      COS_ANG=COSD(ANG)
      SIN_ANG=SIND(ANG)   
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


      SUBROUTINE CALC_ZROT_MX(ROT, ANG)
C
C Create matrix MX for a rotation of ANG around the Z axis.
C                                    
      REAL ROT(3,3)
C
      COS_ANG=COSD(ANG)
      SIN_ANG=SIND(ANG)   
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


	SUBROUTINE ALIGN_VECTORS(U_ROT,H_PEAK,H_HKL,NPEAKS)
C
C Align two sets of NPEAKS vectors (H_PEAK & H_HKL) by rotating
C the second set. LSQ is used to create a linearised rotation
C matrix. It is assumed that all vectors are unit vectors.
C A true rotation matrix corresponding to the rotation of the
C hkl vectors is returned in U_ROT().
C
	REAL U_ROT(3,3),H_PEAK(3,NPEAKS),H_HKL(3,NPEAKS)
C
	PARAMETER (NPEAKS_MAX=100)
	REAL JAC(3,3*NPEAKS_MAX),HESS(3,3),AINV(3,3),DY(3),ROT(3,3),PARS(3)
C
	IF(NPEAKS .GT. NPEAKS_MAX) CALL QUIT('BUG(align_vectors) NPEAKS > 100')
C
C Initialise U_ROT() as a unit matrix
C
	DO I1=1,3
	  DO I2=1,3
	    U_ROT(I1,I2)=0.0
	  ENDDO
	  U_ROT(I1,I1)=1.0
	ENDDO
C
C Create the Jacobian
C
	NITER=0
1	DO I=1,NPEAKS
	  JAC(1,3*I-2)=0.0
	  JAC(2,3*I-2)=				-H_HKL(3,I)
	  JAC(3,3*I-2)=							-H_HKL(2,I)
	  JAC(1,3*I-1)=H_HKL(3,I)
	  JAC(2,3*I-1)=				0.0
	  JAC(3,3*I-1)=							H_HKL(1,I)
	  JAC(1,3*I  )=-H_HKL(2,I)
	  JAC(2,3*I  )=				H_HKL(1,I)
	  JAC(3,3*I  )=							0.0
	ENDDO
C
C HESS = JAC * JAC^T
C
	DO I1=1,3
	  DO I2=1,3
	    HESS(I1,I2)=0.0
	    DO I3=1,3*NPEAKS
	      HESS(I1,I2)=HESS(I1,I2)+JAC(I1,I3)*JAC(I2,I3)
	    ENDDO
	  ENDDO
	ENDDO
C
C DY = JAC * (H_PEAK-H_HKL)
C
	DO I1=1,3
	  DY(I1)=0.0
	  DO I2=1,NPEAKS
	    DO I3=1,3
	      DY(I1)=DY(I1)+JAC(I1,I3+3*I2-3)*(H_PEAK(I3,I2)-H_HKL(I3,I2))
	    ENDDO
	  ENDDO
	ENDDO
C
C Solve for PARS the linear equation HESS * PARS = DY
C If poorly conditioned, complain and give up
C
	CALL MX3_INVERT(AINV,DET, HESS)
	IF(ABS(DET) .LT. 1E-4) THEN
	  PRINT *,'WARNING: ALIGN_VECTORS unable to rotate vectors'
	  RETURN
	ENDIF
	CALL VC3_MULT(PARS,AINV,DY)
C
C Create the linearised rotation matrix
C
	ROT(1,1)=1.0
	ROT(1,2)=			-PARS(3)
	ROT(1,3)=						-PARS(2)
	ROT(2,1)=PARS(3)
	ROT(2,2)=			1.0
	ROT(2,3)=						PARS(1)
	ROT(3,1)=PARS(2)
	ROT(3,2)=			-PARS(1)
	ROT(3,3)=						1.0
C
C Rotate the H_HKL vectors (make sure they are unit vectors)
C
	DO I=1,NPEAKS
	  CALL VC3_MULT(H_HKL(1,I), ROT,H_HKL(1,I))
	  CALL VC3_UNIT(H_HKL(1,I))
	ENDDO
C
C Pre-multiply U_ROT to simulate the effect on the vectors
C
	CALL MX3_MULT(U_ROT, ROT,U_ROT)
C
C If routine hasn't converged then iterate up to 5 times
C
	DPAR=SQRT( PARS(1)**2 + PARS(2)**2 + PARS(3)**2 )
	NITER=NITER+1
	IF(DPAR.GT.1E-4 .AND. NITER.LT.5) GOTO 1
C
C Massage U_ROT to be a true rotation matrix
C
	CALL MAKE_U_ORTHOGONAL(U_ROT)
C
	RETURN
	END


	SUBROUTINE CALC_U_FROM_VECTORS(U, H_HKLS,H_PEAKS, I_HKLS,I_PEAKS)
C
C Calculate the U matrix from 3 vectors corresponding to B * (hkl) and
C 3 scattering vectors from the observed peaks.
C The arrays I_HKLS() and I_PEAKS() contain the indices for H_HKLS() & H_PEAKS().
C
	INTEGER I_HKLS(3),I_PEAKS(3)
C
	PARAMETER (NMAX_HKLS=300000)
	REAL U(3,3),H_PEAKS(3,NMAX_HKLS),H_HKLS(3,NMAX_HKLS)
C
	REAL MAT_PEAKS(3,3),MAT_HKLS(3,3),MAT_INV(3,3)
C
C Copy the HKL & PEAK vectors into matrices
C
	DO I1=1,3
	  DO I2=1,3
	    MAT_HKLS(I1,I2) =H_HKLS(I1,I_HKLS(I2))
	    MAT_PEAKS(I1,I2)=H_PEAKS(I1,I_PEAKS(I2))
	  ENDDO
	ENDDO
C
C Solve for U the linear equation MAT_PEAKS = U * MAT_HKLS
C
	CALL MX3_INVERT(MAT_INV,DET, MAT_HKLS)
	IF(ABS(DET) .LT. 1E-6) CALL QUIT('BUG(calc_u_from_vectors): |DET| < 1e-6')
	CALL MX3_MULT(U, MAT_PEAKS,MAT_INV)
C
C Massage the U matrix to make it orthogonal
C
	CALL MAKE_U_ORTHOGONAL(U)
C
	RETURN
	END


	SUBROUTINE MAKE_U_ORTHOGONAL(U)
C
C Massage U matrix until it is a perfect rotation matrix.
C
	REAL U(3,3)
C
	PARAMETER (TINY=1.0E-7)
	REAL TEMP(3,3)
C
C Do a maximum of 10 iterations of the algorithm.
C
	DO ITER=1,10
C
C Invert U() and transpose and put into TEMP().
C For a pure rotation matrix this should equal U().
C
	  CALL MX3_INVERT(TEMP,DET, U)
	  IF(ABS(DET) .LT. 1.E-20) CALL QUIT('BUG(make_u_orthogonal): DET too small')
	  CALL MX3_TRANS(TEMP,TEMP)
C
C Sum the differences squared between elements in U() and TEMP().
C
	  SUM=0.0
	  DO I1=1,3
	    DO I2=1,3
	      SUM=SUM+(U(I1,I2)-TEMP(I1,I2))**2
	    ENDDO
	  ENDDO
	  SUM=SQRT( SUM/9.0 )
C
C
C Average U() & TEMP() and load into U()
C
	  CALL MX3_ADD(TEMP,TEMP,U)
	  CALL MX3_SCALE(U,TEMP,0.5)
C
C Finish algorithm if SUM is close to the numeric precision.
C
	  IF(SUM .LT. 10.0*TINY) RETURN
C
	ENDDO
C
	RETURN
	END


	SUBROUTINE CALC_B_MATRIX(B, DIRECT)
C
	REAL B(3,3),DIRECT(6)
C
	REAL RECIP(6)
C
C First calculate reciprocal lattice constants.
C
	CALL CONV_DIRECT_RECIP(RECIP,DIRECT)
C
C Busing & Levy's equation for B().
C
	B(1,1)=RECIP(1)
	B(1,2)=		RECIP(2)*COSD(RECIP(6))
	B(1,3)=				RECIP(3)*COSD(RECIP(5))
	B(2,1)=0.0
	B(2,2)=		RECIP(2)*SIND(RECIP(6))
	B(2,3)=				-RECIP(3)*SIND(RECIP(5))*COSD(DIRECT(4))
	B(3,1)=0.0
	B(3,2)=		0.0
	B(3,3)=				1.0/DIRECT(3)
C
	RETURN
	END


	SUBROUTINE CONV_DIRECT_RECIP(RECIP, DIRECT)
C
C Routine accepts a vector of the direct lattice constants
C in DIRECT() and returns the recip. lattice constants
C in RECIP().
C NB:  This routine can also convert from recip. to direct
C      lattice constants.
C
	REAL DIRECT(6),RECIP(6)
C
	A=DIRECT(1)
	B=DIRECT(2)
	C=DIRECT(3)
	ALPHA=DIRECT(4)
	BETA =DIRECT(5)
	GAMMA=DIRECT(6)
C
C Using formulae in Int. Tables. II, p106
C
	S=(ALPHA+BETA+GAMMA)/2.0
	TEMP=SIND(S)*SIND(S-ALPHA)*SIND(S-BETA)*SIND(S-GAMMA)
	V=2.0*A*B*C*SQRT( MAX(0.0,TEMP) )
C
	RECIP(1)=B*C*SIND(ALPHA)/V
	RECIP(2)=A*C*SIND(BETA )/V
	RECIP(3)=A*B*SIND(GAMMA)/V
C
	TEMP=(COSD(BETA)*COSD(GAMMA)-COSD(ALPHA))/(SIND(BETA)*SIND(GAMMA))
	RECIP(4)=ACOSD( MIN(1.0,TEMP) )
	TEMP=(COSD(ALPHA)*COSD(GAMMA)-COSD(BETA))/(SIND(ALPHA)*SIND(GAMMA))
	RECIP(5)=ACOSD( MIN(1.0,TEMP) )
	TEMP=(COSD(ALPHA)*COSD(BETA)-COSD(GAMMA))/(SIND(ALPHA)*SIND(BETA))
	RECIP(6)=ACOSD( MIN(1.0,TEMP) )
C
	RETURN
	END