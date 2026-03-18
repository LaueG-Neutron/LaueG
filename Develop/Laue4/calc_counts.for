C----- Calculate intensities and related values -------------------------------
c	SUBROUTINE CALC_ALL_COUNTS_SEQ
c	SUBROUTINE CALC_ALL_DCOUNTS_SEQ(DCOUNTS_SEQ)
c	FUNCTION CALC_DCOUNTS_SEQ(ISEQ)
c	FUNCTION CALC_COUNTS_SEQ(COUNTS_START, IFIRST,ILAST)
c	SUBROUTINE DERIV_COUNTS_SEQ(VALUE,DERIV, COUNTS_SEQ,  IREF)
C----- Twins ratio correction -------------------------------------------------
c	FUNCTION CALC_TWINS_FACTOR(ITWIN)
C----- Extinction correction --------------------------------------------------
c	FUNCTION CALC_EXTI_FACTOR(COUNTS,WAV,TTH)
C----- Absorption correction --------------------------------------------------
c	FUNCTION CALC_ABS_FACTOR(WAV,TTH)
C----- Wavelength distribution correction -------------------------------------
c	FUNCTION CALC_WAV_SCALE(WAV)
c	FUNCTION CALC_WAV_SPLINE(WAV)
c	FUNCTION CALC_WAV_ROUGH(WAV)
C----- Refined detector efficiency corrections --------------------------------
c	FUNCTION CALC_EFFIC_FACTOR(IREF)
c	FUNCTION CALC_EFFIC_REFINE(Y)
c	FUNCTION CALC_EFFIC_POLY(Y, POLY,NPOLY)
c	FUNCTION CALC_EFFIC_FIXED(WAV,Y)
c	FUNCTION CALC_IP_EFFIC(WAV,Y)
c	FUNCTION CALC_GD_CROSS(WAV)
C----- File scale-factor correction -------------------------------------------
c	FUNCTION CALC_FILE_SCALE(IFILE)
C------------------------------------------------------------------------------


C----- Calculate intensities and related values -------------------------------

	SUBROUTINE CALC_ALL_COUNTS_SEQ
C
C Returns the corrected counts values for all sequences 
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
C Call CALC_COUNTS_SEQ to iterate a value for COUNTS_SEQ()
C
	DO ISEQ=1,NSEQ
	  COUNTS_SEQ(ISEQ)=CALC_COUNTS_SEQ(COUNTS_SEQ(ISEQ), IFIRST(ISEQ),ILAST(ISEQ))
	ENDDO
C
	RETURN
	END


	SUBROUTINE CALC_ALL_DCOUNTS_SEQ(DCOUNTS_SEQ)
C
C Returns the esd of the corrected counts values for all sequences 
C
	PARAMETER (NSEQ_MAX=100000)
	REAL DCOUNTS_SEQ(NSEQ_MAX)
C
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	DO ISEQ=1,NSEQ
	  DCOUNTS_SEQ(ISEQ)=CALC_DCOUNTS_SEQ(ISEQ)
	ENDDO
C
	RETURN
	END


	FUNCTION CALC_DCOUNTS_SEQ(ISEQ)
C
C Returns the esd of the corrected counts value of the sequence ISEQ.
C
	PARAMETER (NDATA_MAX=2000000)
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	NSUM=0
	SDD=0.0
	DO I=IFIRST(ISEQ),ILAST(ISEQ)
C Completely ignore the spot if ISIG() < 0 (used for outliers)
	  IF(ISIG(I) .LT. 0.0) CYCLE
	  CALL DERIV_COUNTS_SEQ(ICALC(I),FDER, COUNTS_SEQ(ISEQ),I)
	  NSUM=NSUM+1
	  SDD=SDD+FDER**2/ISIG(I)**2
	ENDDO
C
C Return the esd of the least-squares refinement of COUNTS_SEQ
C
	IF(NSUM .GE. 1) THEN
	  CALC_DCOUNTS_SEQ=1.0/SQRT(MAX(1E-30,SDD))
	ELSE
C Special case for all spots rejected as outliers, return -1
	  CALC_DCOUNTS_SEQ=-1.0
	ENDIF
C
	RETURN
	END


	FUNCTION CALC_COUNTS_SEQ(COUNTS_START, IFIRST,ILAST)
C
C Returns the "corrected counts" value of a sequence defined by ILAST & IEND.
C The value is calculated iteratively by optimised individual intensities
C ICALC(IFIRST..ILAST) versus IOBS() using weights according to ISIG().
C COUNTS_START is the starting value for the iteration.
C
C Special care is needed for an extinction correction as the optimisation
C may be highly non-linear for the strongest intensities.
C
	PARAMETER (NDATA_MAX=2000000)
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	REAL FDER(NDATA_MAX)
C
C Do LSQ iterations on the COUNTS_SEQ value starting from COUNTS_START
C
	COUNTS_SEQ=COUNTS_START
	DO ITER=1,10
C Calculate sums used for LSQ
	  SDF=0.0
	  SDD=0.0
	  DO I=IFIRST,ILAST
C Completely ignore the spot if ISIG() < 0 (signals an outlier)
	    IF(ISIG(I) .GE. 0.0) THEN
	      CALL DERIV_COUNTS_SEQ(ICALC(I),FDER(I), COUNTS_SEQ,I)
	      SDF=SDF+FDER(I)*(IOBS(I)-ICALC(I))/ISIG(I)**2
	      SDD=SDD+FDER(I)**2/ISIG(I)**2
	    ENDIF
	  ENDDO
C Adjust COUNTS_SEQ using LSQ
	  DEV=SDF/MAX(1E-8,SDD)
	  COUNTS_SEQ=COUNTS_SEQ+DEV
C For no extinction the problem is linear and converges in 1 iteration,
C however, it may need two if the first step is large (i.e. the initial refinement).
C Jump out of loop if the adjustment is small (less small if no extinction).
	  IF( LEXTI_CORR ) THEN
	    IF(ABS(DEV) .LT. 1E-8) GOTO 100
	  ELSE
	    IF(ABS(DEV) .LT. 1E-1) GOTO 100
	  ENDIF
C
	ENDDO
C
C Apply linear correction to ICALC() instead of a full calculation
C
100	DO I=IFIRST,ILAST
	  ICALC(I)=ICALC(I)+FDER(I)*DEV
	ENDDO
C
	CALC_COUNTS_SEQ=COUNTS_SEQ
	RETURN
	END


	SUBROUTINE DERIV_COUNTS_SEQ(VALUE,DERIV, COUNTS_SEQ,  IREF)
C
C Calculate the uncorrected individual spot intensity, VALUE, for spot IREF.
C Uses as input COUNTS_SEQ, the corrected counts value of the sequence which
C IREF is part of.
C Also calculates DERIV, the derivative of VALUE w.r.t. COUNTS_SEQ.
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
C Load parameters for reflection IREF
C
	WAV=WAVS(IREF)
	TTH=TTHS(IREF)
	IFILE=IFILES(IREF)
	ITWIN=HKLMS(4,IREF)
C
C Calculate the extinction term where y = (1 + EX_FAC)**(-1/2) 
C
	EX_FAC=CALC_EXTI_FACTOR(COUNTS_SEQ,WAV,TTH)
C
C Calculate the value and its derivative w.r.t. COUNTS_SEQ
C
	FACTOR=WAV**4 /( SIND(TTH/2) )**2		! Lorentz factor for Laue
	1		* CALC_FILE_SCALE(IFILE)				! scale for each file
	2		* CALC_WAV_SCALE(WAV)						! wavelength correction
	3		* CALC_ABS_FACTOR(WAV,TTH)			! sample absorption
	4		* CALC_EFFIC_FACTOR(IREF)				! total efficiency factor
	5		/ SQRT( MAX(0.5,1.0+EX_FAC) )		! extinction
     6		* CALC_TWINS_FACTOR(ITWIN)			! twins ratio correction
C
      VALUE=FACTOR*COUNTS_SEQ
	DERIV=FACTOR*(1.0+0.5*EX_FAC)/MAX(0.5,1.0+EX_FAC)
C
	RETURN
	END


C----- Twins ratio correction -------------------------------------------------

	FUNCTION CALC_TWINS_FACTOR(ITWIN)
C
	LOGICAL LFILE,LWAV,LEFF,LABS,LEXTI,LTWIN_RATIO
	COMMON /LSQ_LPAR_COM/ LFILE,LWAV,LEFF,LABS,LEXTI,LTWIN_RATIO
C
      LOGICAL LTWIN1,LTWIN2,LTWIN11,LTWIN12,LRATIO
      COMMON /TWINS_CORR_COM/ ITWIN_OPT,TWIN_RATIO,LTWIN1,LTWIN2,LTWIN11,LTWIN12,LRATIO
C
	CALC_TWINS_FACTOR=1.0
cccccccccccccccccccccccc      IF( .NOT.LTWIN_RATIO ) RETURN
C
C Correct for Twin 2 contribution
C
	IF(ITWIN .EQ. 2) THEN
	  CALC_TWINS_FACTOR=TWIN_RATIO
	ELSEIF(ITWIN .EQ. 11) THEN
	  CALC_TWINS_FACTOR=1.0+TWIN_RATIO
	ELSEIF(ITWIN .EQ. 12) THEN
	  CALC_TWINS_FACTOR=1.0-TWIN_RATIO
	ENDIF
C
	RETURN
	END


C----- Extinction correction --------------------------------------------------

	FUNCTION CALC_EXTI_FACTOR(COUNTS,WAV,TTH)
C
C Calculate the extinction term EX_FAC where y = (1 + EX_FAC)**(-1/2) 
C NB: The B&C formulae are not really linear in EXTI so it takes
C	more iterations in CALC_COUNTS_SEQ for COUNTS_SEQ to
C	converge and are therefore much slower.
C
	COMMON /EXTI_CORR_COM/ IEXTI_OPT,EXTI
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
C If no extinction or negative intensity, return with a value of zero
C
	IF( .NOT.LEXTI_CORR .OR. EXTI.EQ.0.0 .OR. COUNTS.LE.0.0) THEN
	  CALC_EXTI_FACTOR=0.0
	  RETURN
	ENDIF
C
	SIN_TTH=SIND(TTH)
	COS_TTH=COSD(TTH)
	X=1.0E-6*EXTI*COUNTS*WAV**2
C
	IF(IEXTI_OPT .EQ. 1) THEN			! Zach Type 1, also SHELX
	  EX_FAC=2.0*X*WAV/SIN_TTH
	ELSE IF(IEXTI_OPT .EQ. 2) THEN	! Zach Type 2
	  EX_FAC=2.0*X/SIN_TTH
	ELSE IF(IEXTI_OPT .EQ. 3) THEN	! Zach Type 1 with B&C correction
	  EX_FAC=2.0*X
	ELSE IF(IEXTI_OPT .EQ. 4) THEN	! B&C I,G
	  A=0.58 + 0.48*COS_TTH + 0.24*COS_TTH**2
	  B=0.02 - 0.025*COS_TTH
	  EX_FAC=2.12*X + A*X*X/(1.0+B*X)
	ELSE IF(IEXTI_OPT .EQ. 5) THEN	! B&C I,L
	  A=0.025 + 0.285*COS_TTH
	  IF(COS_TTH .GE. 0.0) THEN
	    B=0.15 - 0.2*(0.75-COS_TTH)**2
	  ELSE
	    B=-0.45*COS_TTH
	  ENDIF
	  EX_FAC=2.0*X + A*X*X/(1.0+B*X)
	ELSE IF(IEXTI_OPT .EQ. 6) THEN	! B&C II
	  A=0.20 + 0.45*COS_TTH
	  B=0.22 - 0.12*(0.5-COS_TTH)**2
	  EX_FAC=2.0*X + A*X*X/(1.0+B*X)
	ELSE
	  STOP 'BUG(calc_exti_factor): Invalid IMODEL'
	ENDIF
C
	CALC_EXTI_FACTOR=EX_FAC
C
	RETURN
	END


C----- Absorption correction --------------------------------------------------

	FUNCTION CALC_ABS_FACTOR(WAV,TTH)
C
C Calculate the absorption factor for a spherical sample
C
C For hydrogenated organics we expect mu = (0.14 - 0.22) * ( wav + 1.07).
C For radius = 0.5 - 1 mm, mu * R = (0.075 - 0.24) + (0.07 - 0.22) * wav,
C which at wavelength of 0.85 - 1.7 we get mu * R = 0.13 - 0.61. 
C
	COMMON /ABS_CORR_COM/ IABS_OPT,UR0,UR_LIN
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
C Return 1.0 if the correction is turned off
C
	IF( .NOT.LABS_CORR ) THEN
	  CALC_ABS_FACTOR=1.0
	  RETURN
	ENDIF
C
C Calculate uR (mu * radius) as a function of WAV
C NB: UNPACK_ABS_PARS calculates UR0 & UR_LIN depending on IABS_OPT
C     IABS_OPT = 1	UR0=0, UR_LIN refined
C     IABS_OPT = 2	UR0=1.07*UR_LIN, UR_LIN refined
C     IABS_OPT = 3	UR0 and UR_LIN refined
C
	uR=UR0+UR_LIN*WAV
C
C Inverse of absorption factor for a sphere. Approximation to the 
C Int.Tables Vol II and good for uR < 2 & TTH=20-150.
C NB: A simpler approximation would be EXP(-uR*pi/2)
C
	ABS=1.0 - (0.205-0.059*uR-0.006*uR**2)*uR**2 * COSD(TTH)
CC Ignore the non-TTH dependent part of the correction
CC	ABS=ABS * EXP(-1.5146*uR) * (1.0+0.13*uR*EXP(uR))
	CALC_ABS_FACTOR=ABS
C
	RETURN
	END


C----- Wavelength distribution correction -------------------------------------

	FUNCTION CALC_WAV_SCALE(WAV)
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	COMMON /WAV_CORR_COM/ IWAV_OPT,NWAV,WAV_GAMMA,WAV_MIN,WAV_STEP,WAV_PAR(20)
C
C Return the default distribution if the full correction is turned off
C
	IF( .NOT.LWAV_CORR ) THEN
	  CALC_WAV_SCALE=CALC_WAV_ROUGH(WAV)
	  RETURN
	ENDIF
C
C IOPT = 1	Quadratic spline
C        2	Quadratic spline, then non-param correction
C
	IF(IWAV_OPT .EQ. 1) THEN
	  CORR=CALC_WAV_SPLINE(WAV)
	ELSE IF(IWAV_OPT .EQ. 2) THEN
	  CORR=CALC_WAV_SPLINE(WAV) * CALC_WAV_NONPAR(WAV)
     	ELSE
	  STOP' BUG(calc_wav_scale): Invalid IWAV_OPT'
	ENDIF
C
C Multiply by rough wavelength and the non-param correction
C
	CALC_WAV_SCALE=CORR*CALC_WAV_ROUGH(WAV)
C
C
	RETURN
	END


	FUNCTION CALC_WAV_SPLINE(WAV)
C
C Evaluates a wavelength distribution factor based on a spline curve
C
	COMMON /WAV_CORR_COM/ IWAV_OPT,NWAV,WAV_GAMMA,WAV_MIN,WAV_STEP,WAV_PAR(20)
C
	CALC_WAV_SPLINE=CALC_QUAD_SPLINE(WAV-WAV_MIN,WAV_STEP,WAV_PAR,NWAV)
C
	RETURN
	END


	FUNCTION CALC_WAV_ROUGH(WAV)
C
	COMMON /WAV_FILE_COM/ NWAV_FILE,WAV_FILE_MIN,WAV_FILE_STEP,WAV_FILE(10000)
C
	ROUGH=CALC_LIN_SPLINE(WAV-WAV_FILE_MIN, WAV_FILE_STEP,WAV_FILE,NWAV_FILE)
C Never let it get to zero as things start crashing!
	CALC_WAV_ROUGH=MAX(1E-6,ROUGH)
C
	RETURN
	END


C----- Refined detector efficiency corrections --------------------------------

	FUNCTION CALC_EFFIC_FACTOR(IREF)
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	Y=XY_PIX(2,IREF)
	CALC_EFFIC_FACTOR = CALC_EFFIC_FIXED(WAVS(IREF),Y) * CALC_EFFIC_REFINE(Y)
C
	RETURN
	END


	FUNCTION CALC_EFFIC_REFINE(Y)
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
	COMMON /EFF_CORR_COM/ IEFF_OPT,EFF_POLY(10)
C
C Return 1.0 if the correction is not being refined,
C otherwise, calculate the polynomial efficiency correction
C
	IF( .NOT.LEFF_CORR ) THEN
	  CALC_EFFIC_REFINE=1.0
	ELSE
	  CALC_EFFIC_REFINE = CALC_EFFIC_POLY(Y, EFF_POLY,IEFF_OPT)
	ENDIF
C
	RETURN
	END


	FUNCTION CALC_EFFIC_POLY(Y, POLY,NPOLY)
C
	REAL POLY(NPOLY)
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
C Calculate the polynomial variation in efficiency versus Y
C
C The Legendre polynomials are functions of DY which roughly span
C the range -1 to +1 with a constant density of spots
C
	CALL GET_NUMXY(NUMX,NUMY)
	YCEN=NUMY/2.0
C
	DY=Y/YCEN-1.0
	DY=3.0*DY/(2.0+ABS(DY))
	CORR=1.0
	DO I=1,NPOLY
	  CORR=CORR+POLY(I)*POLY_ORTHO(DY,I)
	ENDDO
	CALC_EFFIC_POLY=CORR
C
	RETURN
	END


	FUNCTION CALC_EFFIC_FIXED(WAV,Y)
C
	COMMON /EFF_CORR_COM/ IEFF_OPT,EFF_POLY(10)
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
C Apply obliqueness correction for IP efficiency (relative to the Y=EFF_YCEN efficiency)
C
	CALL GET_NUMXY(NUMX,NUMY)
	EFF_YCEN=0.5*NUMY
	CALC_EFFIC_FIXED=CALC_IP_EFFIC(WAV,Y) / CALC_IP_EFFIC(WAV,EFF_YCEN)
C
	RETURN
	END


	FUNCTION CALC_IP_EFFIC(WAV,Y)
C
	COMMON /EFF_CORR_COM/ IEFF_OPT,EFF_POLY(10)
C
	COMMON /MISC_PARS_COM/ NSEQ_MIN,SEQ_FRAC,X_LO,X_HI,Y_LO,Y_HI,WAV_LO,WAV_HI,
	1	CELL_MULT,EFF_WIDTH,EFF_THICK,EXTI_CORR_MAX,NBINS_NPAR,
	2	ALL_REL_ERR,ALL_SIG_MULT,WEAK_SIG_MULT,EXTI_ERR,WAV_ERR1,WAV_ERR2,
	3	HARM_CUTOFF_MULT,IVERBOSE
C
C Use the attenuation measurement value
	RMU=8.04E-5*CALC_GD_CROSS(WAV)
C ROBLIQ is the multiplier of the IP path length due to an oblique beam
	CALL GET_NUMXY(NUMX,NUMY)
	EFF_YCEN=0.5*NUMY
	ROBLIQ=SQRT( 1.0 + ( (Y-EFF_YCEN)/EFF_WIDTH )**2 )
	CALC_IP_EFFIC=1.0 - EXP(-RMU*EFF_THICK*ROBLIQ)
C
	RETURN
	END


	FUNCTION CALC_GD_CROSS(WAV)
C
C Calculate Gd total cross-section in barns for a wavelength in Angstroms.
C Fits the barn-book graph within 4% over 0.62 to 9.1 Angstroms.
C
	ALOG_E=4.4045-2.0*LOG(WAV)
	Y1=11.3982-0.374484*ALOG_E
	Y2=19.7400-2.37120*ALOG_E
	ARG=( Y1**(-24.7016) + Y2**(-24.5142) )**(-0.0415175)
	CALC_GD_CROSS=EXP(ARG)
C
	RETURN
	END


C----- File scale-factor correction -------------------------------------------

	FUNCTION CALC_FILE_SCALE(IFILE)
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_CORR_COM/ FSCALE(NFILES_MAX)
C
	LOGICAL LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
	COMMON /LCORR_COM/ LFILE_CORR,LWAV_CORR,LEFF_CORR,LEXTI_CORR,LHARM_CORR,LABS_CORR
C
C Return with 1.0 if the correction is turned off
C
	IF( .NOT.LFILE_CORR ) THEN
	  CALC_FILE_SCALE=1.0
	  RETURN
	ENDIF
C
C The file scale factor is simply the FSCALE() value
C
	CALC_FILE_SCALE=FSCALE(IFILE)
C
	RETURN
	END
