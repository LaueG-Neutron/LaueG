C ============= Same as in ORIENT_SPOTS, HKL_GEN, INDEX_CONICS  ==========

C ============= Convenience routine for HKLS to HVECS ====================
C	SUBROUTINE CALC_HKLS_TO_HVECS(HVECS, UB,HKLS,NHKLS)
C ============= Convenience routines for PIXWAVS to/from HKLS ============
C	SUBROUTINE CALC_HKLS_TO_PIXWAVS(PIXWAVS, UB,HKLS,NHKLS)
C	SUBROUTINE CALC_PIXWAVS_TO_HKLS(HKLS, UB,PIXWAVS,NHKLS)
C ============= Convenience routines for PIXWAVS to/from HVECS ===========
C	SUBROUTINE CALC_PIXWAVS_TO_HVECS(HVECS, PIXWAVS,NHKLS)
C	SUBROUTINE CALC_HVECS_TO_PIXWAVS(PIXWAVS, HVECS,NHKLS)
C ========== Base Routines for PIXWAVS >> SVECS >> HVECS >> HKLS =========
C	SUBROUTINE CALC_PIXEL_TO_SVEC(S1X,S1Y,S1Z, XPIX,YPIX,OFFSETS,XTAL_OFFX)
C	SUBROUTINE CALC_SVEC_TO_HVEC(HVX,HVY,HVZ, S1X,S1Y,S1Z,WAV)
C	SUBROUTINE CALC_HVEC_TO_HKL(RH,RK,RL, UB_ROT_INV,HVX,HVY,HVZ)
C ========== Base Routines for HKLS >> HVECS >> SVECS >> PIXWAVS =========
C	SUBROUTINE CALC_HKL_TO_HVEC(HVX,HVX,HVX, UB_ROT,RH,RK,RL)
C	SUBROUTINE CALC_HVEC_TO_SVEC(S1X,S1Y,S1Z,WAV,DSPACE, HVX,HVY,HVZ)
C	SUBROUTINE CALC_SVEC_TO_PIXEL(XPIX,YPIX, OFFSETS,XTAL_OFFX,S1X,S1Y,S1Z)
C ========== Correction to reduce LSQ correlations =======================
C	SUBROUTINE XTAL_OFFX_ENABLE()
C	SUBROUTINE XTAL_OFFX_DISABLE()
C=========== Lattice/UB Related Routines =================================
C	SUBROUTINE CALC_UB_TO_CELL(CELL, UB)
C	SUBROUTINE CALC_UB_MATRIX(UB)
C	SUBROUTINE ALIGN_VECTORS(U_ROT,H_PEAK,H_HKL,NPEAKS)
C	SUBROUTINE CALC_U_FROM_VECTORS(U, H_HKLS,H_PEAKS, I_HKLS,I_PEAKS)
C	SUBROUTINE MAKE_U_ORTHOGONAL(U)
C	SUBROUTINE CALC_B_MATRIX(B, DIRECT)
C	SUBROUTINE CONV_DIRECT_RECIP(RECIP, DIRECT)
C============ Matrix rotations ===========================================
C	SUBROUTINE CALC_ZYXROT_MX(ROT, ANGX,ANGY,ANGZ)
C	SUBROUTINE CALC_XYZROT_MX(ROT, ANGX,ANGY,ANGZ)
C	SUBROUTINE CALC_XROT_MX(ROT, ANG)
C	SUBROUTINE CALC_YROT_MX(ROT, ANG)
C	SUBROUTINE CALC_ZROT_MX(ROT, ANG)
C==================== Geometry Definitions ===============================
C==================== Geometric Relationships ============================
C==================== Parameter Correlations =============================
C==================== COMMON Parameters  =================================
C==================== PARAMETER Values ===================================
C=========================================================================

C ============= Convenience routine for HKLS to HVECS ====================

	SUBROUTINE CALC_HKLS_TO_HVECS(HVECS, UB,HKLS,NHKLS)
C
C Convert hkls to reciprocal lattice vectors (length = 1/d-spacing)
C
	REAL HVECS(3,*),UB(3,3),HKLS(3,*)
C
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	REAL ROT_PHI(3,3),UB_ROT(3,3)
C
C Rotate UB by -PHI
C
	CALL CALC_YROT_MX(ROT_PHI, -PHI)
	CALL MX3_MULT(UB_ROT, ROT_PHI,UB)
C
	DO I=1,NHKLS
	  CALL CALC_HKL_TO_HVEC(HVECS(1,I),HVECS(2,I),HVECS(3,I), UB_ROT,HKLS(1,I),HKLS(2,I),HKLS(3,I))
	ENDDO
C
	RETURN
	END


C ============= Convenience routines for PIXWAVS to/from HKLS ==============

	SUBROUTINE CALC_PIXWAVS_TO_HKLS(HKLS, UB,PIXWAVS,NHKLS)
C
C Convert pixel X,Y and wavelengths to diffracted beam vectors
C Beam vectors have lengths of 1/wavelength, as given by PIXWAVS(3,*)
C
	REAL HKLS(3,*),UB(3,3),PIXWAVS(3,*)
C
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	LOGICAL LXTAL_OFFX
	COMMON /LXTAL_OFFX_COM/ LXTAL_OFFX
C
	REAL OFFSETS(3),ROT_PHI(3,3),UB_ROT(3,3),UB_ROT_INV(3,3)
C
C Calculate sample position offset in LAB coords
C
	OFFSETS(1)=COSD(-PHI)*XTAL_OFF(1) -SIND(-PHI)*XTAL_OFF(3) -DRUM_OFF(1)
	OFFSETS(2)=XTAL_OFF(2)                                    -DRUM_OFF(2)
	OFFSETS(3)=SIND(-PHI)*XTAL_OFF(1) +COSD(-PHI)*XTAL_OFF(3) -DRUM_OFF(3)
C
C Optional empirical correction that reduces correlations between XCEN & X offset
C
	XTAL_OFFX=0.0
	IF( LXTAL_OFFX ) XTAL_OFFX=( COSD(-PHI)*XTAL_OFF(1) -SIND(-PHI)*XTAL_OFF(3) )
C
C Rotate UB by -PHI, then invert the matrix
C
	CALL CALC_YROT_MX(ROT_PHI, -PHI)
	CALL MX3_MULT(UB_ROT, ROT_PHI,UB)
	CALL MX3_INVERT(UB_ROT_INV,DET, UB_ROT)
C
	DO I=1,NHKLS
	  CALL CALC_PIXEL_TO_SVEC(S1X,S1Y,S1Z, PIXWAVS(1,I),PIXWAVS(2,I),OFFSETS,XTAL_OFFX)
	  CALL CALC_SVEC_TO_HVEC(HVX,HVY,HVZ, S1X,S1Y,S1Z,PIXWAVS(3,I))
	  CALL CALC_HVEC_TO_HKL(HKLS(1,I),HKLS(2,I),HKLS(3,I), UB_ROT_INV,HVX,HVY,HVZ)
	ENDDO
C
	RETURN
	END


	SUBROUTINE CALC_HKLS_TO_PIXWAVS(PIXWAVS, UB,HKLS,NHKLS)
C
C Convert pixel X,Y and wavelengths to diffracted beam vectors
C Beam vectors have lengths of 1/wavelength, as given by PIXWAVS(3,*)
C
	REAL HKLS(3,*),UB(3,3),PIXWAVS(3,*)
C
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	LOGICAL LXTAL_OFFX
	COMMON /LXTAL_OFFX_COM/ LXTAL_OFFX
C
	REAL OFFSETS(3),ROT_PHI(3,3),UB_ROT(3,3)
C
C Combine drum and xtal offsets for lab coords
C
	OFFSETS(1)=COSD(-PHI)*XTAL_OFF(1) -SIND(-PHI)*XTAL_OFF(3) -DRUM_OFF(1)
	OFFSETS(2)=XTAL_OFF(2)                                    -DRUM_OFF(2)
	OFFSETS(3)=SIND(-PHI)*XTAL_OFF(1) +COSD(-PHI)*XTAL_OFF(3) -DRUM_OFF(3)
C
C Optional empirical correction that reduces correlations between XCEN & X offset
C
	XTAL_OFFX=0.0
	IF( LXTAL_OFFX ) XTAL_OFFX=( COSD(-PHI)*XTAL_OFF(1) -SIND(-PHI)*XTAL_OFF(3) )
C
C Rotate UB by -PHI
C
	CALL CALC_YROT_MX(ROT_PHI, -PHI)
	CALL MX3_MULT(UB_ROT, ROT_PHI,UB)
C
	DO I=1,NHKLS
	  CALL CALC_HKL_TO_HVEC(HVX,HVY,HVZ, UB_ROT,HKLS(1,I),HKLS(2,I),HKLS(3,I))
	  CALL CALC_HVEC_TO_SVEC(S1X,S1Y,S1Z,PIXWAVS(3,I), HVX,HVY,HVZ)
	  CALL CALC_SVEC_TO_PIXEL(PIXWAVS(1,I),PIXWAVS(2,I), OFFSETS,XTAL_OFFX,S1X,S1Y,S1Z)
	ENDDO
C
	RETURN
	END


C ============= Convenience routines for PIXWAVS to/from HVECS ==============

	SUBROUTINE CALC_PIXWAVS_TO_HVECS(HVECS, PIXWAVS,NHKLS)
C
C Convert pixel X,Y and wavelengths to diffracted beam vectors
C Beam vectors have lengths of 1/wavelength, as given by PIXWAVS(3,*)
C
	REAL HVECS(3,*),PIXWAVS(3,*)
C
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	LOGICAL LXTAL_OFFX
	COMMON /LXTAL_OFFX_COM/ LXTAL_OFFX
C
	REAL OFFSETS(3)
C
C Calculate sample position offset in LAB coords
C
	OFFSETS(1)=COSD(-PHI)*XTAL_OFF(1) -SIND(-PHI)*XTAL_OFF(3) -DRUM_OFF(1)
	OFFSETS(2)=XTAL_OFF(2)                                    -DRUM_OFF(2)
	OFFSETS(3)=SIND(-PHI)*XTAL_OFF(1) +COSD(-PHI)*XTAL_OFF(3) -DRUM_OFF(3)
C
C Optional empirical correction that reduces correlations between XCEN & X offset
C
	XTAL_OFFX=0.0
	IF( LXTAL_OFFX ) XTAL_OFFX=( COSD(-PHI)*XTAL_OFF(1) -SIND(-PHI)*XTAL_OFF(3) )
C
	DO I=1,NHKLS
	  CALL CALC_PIXEL_TO_SVEC(S1X,S1Y,S1Z, PIXWAVS(1,I),PIXWAVS(2,I),OFFSETS,XTAL_OFFX)
	  CALL CALC_SVEC_TO_HVEC(HVECS(1,I),HVECS(2,I),HVECS(3,I), S1X,S1Y,S1Z,PIXWAVS(3,I))
	ENDDO
C
	RETURN
	END


	SUBROUTINE CALC_HVECS_TO_PIXWAVS(PIXWAVS, HVECS,NHKLS)
C
C Convert diffracted beam vectors (reciprocal lattice points) to pixel X,Y and wavelength
C
	REAL PIXWAVS(3,*),HVECS(3,*)
C
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	LOGICAL LXTAL_OFFX
	COMMON /LXTAL_OFFX_COM/ LXTAL_OFFX
C
	REAL OFFSETS(3)
C
C Combine drum and xtal offsets for lab coords
C
	OFFSETS(1)=COSD(-PHI)*XTAL_OFF(1) -SIND(-PHI)*XTAL_OFF(3) -DRUM_OFF(1)
	OFFSETS(2)=XTAL_OFF(2)                                    -DRUM_OFF(2)
	OFFSETS(3)=SIND(-PHI)*XTAL_OFF(1) +COSD(-PHI)*XTAL_OFF(3) -DRUM_OFF(3)
C
C Optional empirical correction that reduces correlations between XCEN & X offset
C
	XTAL_OFFX=0.0
	IF( LXTAL_OFFX ) XTAL_OFFX=COSD(-PHI)*XTAL_OFF(1) - SIND(-PHI)*XTAL_OFF(3)
C
	DO I=1,NHKLS
	  CALL CALC_HVEC_TO_SVEC(S1X,S1Y,S1Z,PIXWAVS(3,I), HVECS(1,I),HVECS(2,I),HVECS(3,I))
	  CALL CALC_SVEC_TO_PIXEL(PIXWAVS(1,I),PIXWAVS(2,I), OFFSETS,XTAL_OFFX,S1X,S1Y,S1Z)
	ENDDO
C
	RETURN
	END


C ========== Base Routines for PIXWAVS >> SVECS >> HVECS >> HKLS =========

	SUBROUTINE CALC_PIXEL_TO_SVEC(S1X,S1Y,S1Z, XPIX,YPIX,OFFSETS,XTAL_OFFX)
C
C Convert pixel X,Y and wavelengths to diffracted beam unit-vectors
C
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	REAL OFFSETS(3)
C
C Spot position in pixels relative to the pixel center
	DX=XPIX -PIX_CEN(1)
	DY=YPIX -PIX_CEN(2)
C If CYCLOPS, get the plate number and shift X pixel relative to plate 5
	IF(ITYPE .EQ. 3) THEN
	  IPLATE=1+IFIX( XPIX / (NUMX/8.0) )
	  DX=DX -(IPLATE-5)*NUMX/8.0
	ENDIF
C Convert spot position from pixels to mm across the detector
	DX=DX*PIX_SIZE(1)
	DY=DY*PIX_SIZE(2)
C Correct for skewness
	DX=DX + PIX_SKEW*DY
C Possibly adjust X for XTAL_OFFX used to reduce parameter correlations in LSQ
	DX=DX - XTAL_OFFX
C
	IF(ITYPE .NE. 3) THEN
C=========== Drum detectors
C Calculate S1, the vector from the sample to the spot
	  ANGLE=DX/DRUM_RAD
	  S1X=DRUM_RAD*SIN(ANGLE) -OFFSETS(1)
	  S1Y=DY                  -OFFSETS(2)
	  S1Z=DRUM_RAD*COS(ANGLE) -OFFSETS(3)
	ELSE
C=========== Octagonal CYCLOPS detector
C Rotate sample offset from laboratory to detector-plate coords
	  YROT=45.0*(IPLATE-5)
	  OFFSET_X=COSD(YROT)*OFFSETS(1) -SIND(YROT)*OFFSETS(3)
	  OFFSET_Z=SIND(YROT)*OFFSETS(1) +COSD(YROT)*OFFSETS(3)
C Calculate S1, the vector from the sample to the spot
	  S1X=DX-OFFSET_X
	  S1Y=DY-OFFSETS(2)
	  S1Z=DRUM_RAD-OFFSET_Z
C Rotate S1 from the detector-plate to laboratory coords
	  TMP=COSD(-YROT)*S1X -SIND(-YROT)*S1Z
	  S1Z=SIND(-YROT)*S1X +COSD(-YROT)*S1Z
	  S1X=TMP
	ENDIF
C
C Scale S1 to a unit vector
C
	SIZE=SQRT(S1X**2 + S1Y**2 + S1Z**2)
	S1X=S1X/SIZE
	S1Y=S1Y/SIZE
	S1Z=S1Z/SIZE
C
	RETURN
	END


	SUBROUTINE CALC_SVEC_TO_HVEC(HVX,HVY,HVZ, S1X,S1Y,S1Z,WAV)
C
C Convert diffracted beam unit-vector and wavelength to the scattering vector of length = 1/d-spacing
C
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
C
C Incident-beam vector (S0X is zero)
C
	S0Y=SIND(BEAM_VERT)
	S0Z=COSD(BEAM_VERT)
C
C Calculate scattering vectors, HVEC = (SVEC - S0)/WAV
C
	HVX=(S1X     )/WAV
	HVY=(S1Y -S0Y)/WAV
	HVZ=(S1Z -S0Z)/WAV
C
	RETURN
	END


	SUBROUTINE CALC_HVEC_TO_HKL(RH,RK,RL, UB_ROT_INV,HVX,HVY,HVZ)
C
C Convert scattering vectors to Miller indices
C
	REAL UB_ROT_INV(3,3)
C
C HKLS() = UB_ROT_INV() * HVECS()
C
	RH=UB_ROT_INV(1,1)*HVX + UB_ROT_INV(1,2)*HVY + UB_ROT_INV(1,3)*HVZ
	RK=UB_ROT_INV(2,1)*HVX + UB_ROT_INV(2,2)*HVY + UB_ROT_INV(2,3)*HVZ
	RL=UB_ROT_INV(3,1)*HVX + UB_ROT_INV(3,2)*HVY + UB_ROT_INV(3,3)*HVZ
C
	RETURN
	END


C ========== Base Routines for HKLS >> HVECS >> SVECS >> PIXWAVS =========

	SUBROUTINE CALC_HKL_TO_HVEC(HVX,HVY,HVZ, UB_ROT,RH,RK,RL)
C
	REAL UB_ROT(3,3)
C
C Calculate the scatterings vector for the reflection (RH,RK,RL)
C
	HVX=UB_ROT(1,1)*RH +UB_ROT(1,2)*RK +UB_ROT(1,3)*RL
	HVY=UB_ROT(2,1)*RH +UB_ROT(2,2)*RK +UB_ROT(2,3)*RL
	HVZ=UB_ROT(3,1)*RH +UB_ROT(3,2)*RK +UB_ROT(3,3)*RL
C
	RETURN
	END


	SUBROUTINE CALC_HVEC_TO_SVEC(S1X,S1Y,S1Z,WAV, HVX,HVY,HVZ)
C
C Calculates diffracted-beam unit-vectors from the scattering vectors (of length 1/d-spacing)
C Also returns the wavelength
C
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
C
C Incident-beam vector S0, as a unit-vector
C
C	S0X=0.0
	S0Y=SIND(BEAM_VERT)
	S0Z=COSD(BEAM_VERT)
C
C Calculate wavelength using Braggs' Law
C
	SIZE=SQRT( HVX**2 + HVY**2 + HVZ**2 )
	SIZE=MAX(1E-10, SIZE )
	SIN_TH=-( S0Y*HVY + S0Z*HVZ )/SIZE
	DSPACE=1/SIZE
	WAV=2.0*SIN_TH*DSPACE
C
C Combine incident-beam and scattering vectors to give the diffracted-beam unit-vector
C
	S1X=HVX*WAV
	S1Y=HVY*WAV + S0Y
	S1Z=HVZ*WAV + S0Z
C
	RETURN
	END


	SUBROUTINE CALC_SVEC_TO_PIXEL(XPIX,YPIX, OFFSETS,XTAL_OFFX,S1X,S1Y,S1Z)
C
C Calculate pixel positions from the diffracted-beam unit-vectors
C
	REAL OFFSETS(3)
C
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	IF(ITYPE .NE. 3) THEN
C=========== The calculation for cylindrical detector
C Distance from the sample to the spot on the detector
	  FAC1=S1X**2+S1Z**2
	  FAC2=OFFSETS(1)*S1X + OFFSETS(3)*S1Z
	  FAC3=OFFSETS(1)*S1Z - OFFSETS(3)*S1X
	  DIST =( SQRT( FAC1*DRUM_RAD**2 - FAC3**2 ) - FAC2 ) / FAC1
C Position of the spot relative to sample origin
	  X=DIST*S1X+OFFSETS(1)
	  Y=DIST*S1Y+OFFSETS(2)
	  Z=DIST*S1Z+OFFSETS(3)
C Convert X,Y,Z to spot position in mm across the IP cylinder
C NB: must use the radian verions of ATAN
	  DX=DRUM_RAD*ATAN2(X,Z)
	  DY=Y
	ELSE
C=========== The calculation for octagonal detector
C Plate the beam will hit
	  IPLATE=5+NINT( ATAN2D(S1X,S1Z)/45.0 )
	  IF(IPLATE .EQ. 9) IPLATE=1
C Rotate S1 & offset from laboratory to detector coords
	  YROT=45.0*(IPLATE-5)
	  TMP=COSD(YROT)*S1X -SIND(YROT)*S1Z
	  S1Z=SIND(YROT)*S1X +COSD(YROT)*S1Z
	  S1X=TMP
	  OFFSET_X=COSD(YROT)*OFFSETS(1) -SIND(YROT)*OFFSETS(3)
	  OFFSET_Z=SIND(YROT)*OFFSETS(1) +COSD(YROT)*OFFSETS(3)
C Position of the spot in detector-plate coords
	  DIST=(DRUM_RAD-OFFSET_Z)/S1Z
	  DX=DIST*S1X+OFFSET_X
	  DY=DIST*S1Y+OFFSETS(2)
	ENDIF
C
C Possibly, use empirical correction to reduce correlations between XCEN & XTAL_OFFSET
	DX=DX + XTAL_OFFX
C Convert to pixels and add the pixel centers of the plate
	XPIX=PIX_CEN(1) + (DX - DY*PIX_SKEW)/PIX_SIZE(1)
	YPIX=PIX_CEN(2) +  DY               /PIX_SIZE(2)
C For CYCLOPS, add the pixel centre of plate 5
	IF(ITYPE .EQ. 3) XPIX=XPIX +(IPLATE-5)*NUMX/8.0
C
	RETURN
	END


C ========== Correction to reduce LSQ correlations ===================

	SUBROUTINE XTAL_OFFX_ENABLE()
C
C Enable correction to reduce correlations between XCEN & XTAL_OFFSET
C
	LOGICAL LXTAL_OFFX
	COMMON /LXTAL_OFFX_COM/ LXTAL_OFFX
C
	LXTAL_OFFX=.TRUE.
C
	RETURN
	END


	SUBROUTINE XTAL_OFFX_DISABLE()
C
C Disable correction to reduce correlations between XCEN & XTAL_OFFSET
C
	LOGICAL LXTAL_OFFX
	COMMON /LXTAL_OFFX_COM/ LXTAL_OFFX
C
	LXTAL_OFFX=.FALSE.
C
	RETURN
	END


C =================== Lattice/UB Related =============================

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
	PARAMETER (NHKLS_MAX=300000)
	REAL U(3,3),H_PEAKS(3,NHKLS_MAX),H_HKLS(3,NHKLS_MAX)
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


C==================== Matrix rotations ===============================

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


C==================== Geometry Definitions ===========================
C We fix nominal XYZ axes as: left (looking downstream), up, downstream.
C Unspecified positions in mm, all angles anti-clockwise in degrees.
C
C	BEAM_ANGLE			Angle of incident beam from horizontal
C	XTAL_OFFSET(3)		Position offset of crystal center
C	DRUM_OFFSET(3)		Position offset of IP drum axis
C	DRUM_RADIUS			Radius of IP surface inside the drum
C	PIXEL_SIZE(2)		X & Y pixel size in microns
C	PIXEL_SKEW			X/Y pixel skewness factor
C	PIXEL_CEN(2)		Position of TTH=0 in pixels
C	PHI					Spindle rotation angle
C	CELL(6)				Sample cell parameters
C	UB_RXYZ(3)			Rotation angles to correct UB matrix
C	U_BASE(3,3)			Approximate U matrix fixed at program start
C
C	X_PIXEL,Y_PIXEL		Position on the image plate in terms of pixels
C====================================================================

C==================== Geometric Relationships ========================
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
C====================================================================

C==================== Parameter Correlations =========================
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
C====================================================================

C==================== COMMON Parameters  =============================
C To define instrument geometry
C	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C ! ITYPE: 1=KOALA1,2=IMAGINE,3=CYCLOPS,4=KOALA2
C	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
C	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
C To define Lattice or to make UB
C	COMMON /CELL_COM/ CELL(6),ILATT,ICEN,CELL_PROD
C	COMMON /U_COM/ U_BASE(3,3),RXYZ(3)
C=====================================================================

C==================== PARAMETER Values ===============================
C	PARAMETER (NHKLS_MAX=300000)
C=====================================================================
