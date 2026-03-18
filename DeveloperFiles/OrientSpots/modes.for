C =================================================
C	SUBROUTINE RUN_AUTO_MODE(UB, LCEN_FIXED)
C	SUBROUTINE RUN_MANUAL_MODE(UB,IMODE)
C	SUBROUTINE RUN_HKL_MODE(UB)
C	SUBROUTINE SETUP_UB_MATRIX(UB)
C =================================================


	SUBROUTINE RUN_AUTO_MODE(UB, LCEN_FIXED)
C
	LOGICAL LCEN_FIXED
	REAL UB(3,3)
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
C
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
C
	COMMON /MODE_COM/ HKLGEN1,HKLGEN2,DIST_MATCH,IHKLGEN,IREFS(6)
C
C Output starting values for pixel center and UB matrix
C
	PRINT '(/,1X,A,2F8.2)','Initial pixel(xy) =',PIX_CEN
	PRINT '(1X,A,3(/,8X,3F10.6))','Initial UB matrix:',
	1				((UB(K2,K1),K1=1,3),K2=1,3)
C
C ================== Phase 1 ================
C
C Try to roughly optimise the pixel centers and the UB rotation
C using a small number of low index hkls and strong spots
C
	PRINT '(/,1X,A)','Initial Spot Match and Optimisation'
	CALL ORIENT_SPOTS_LEVEL1(RMS, UB, LCEN_FIXED)
C
C Setup the cell dimensions to use in LSQ refinement
C
	CALL INIT_CELL_REFINE()
C
C Recalculate UB matrix from the cell dimensions
C
	CALL SETUP_UB_MATRIX(UB)
C
C Do a LSQ refinement just on pixel center and UB rotation
C
	PRINT '(/,1X,A)','Initial Orientation Refinement'
	CALL SET_LSQ_REFINE(1,0,1,0,0,0)
	CALL RUN_LSQ()
C
C ================== Phase 2 ================
C
C Do a LSQ refinement of pixel center, rotations and cell for
C the low/medium index hkls and strong/medium spots
C
C Create a medium size match list (~50 matches)
C
	PRINT '(/,1X,A)','Second Orientation Refinement'
	CALL MAKE_PARTIAL_MATCH(RMS,NMATCH,NTRY)
	IF(NMATCH .LT. 5) GOTO 900
C
C Run the LSQ refining the cell, rotations, and pixel centers
C
	CALL SET_LSQ_REFINE(1,1,1,0,0,0)
	CALL RUN_LSQ()
C
C ================== Phase 3 ================
C
C Create a full list of matches, but restrict the parameters
C
	PRINT '(/,1X,A)','Third Orientation Refinement'
	CALL MAKE_FULL_MATCH()
	IF(NMATCH .LT. 10) GOTO 900
C
C Run the LSQ refining the cell, rotations, and pixel centers
C
	CALL SET_LSQ_REFINE(1,1,1,0,0,0)
	CALL RUN_LSQ()
C
C ================== Phase 4 ================
C
C Create a full list of matches, and try to refine all parameters
C
C Create a large list of matches
C
	PRINT '(/,1X,A)','Final Orientation Refinement'
	CALL MAKE_FULL_MATCH()
	IF(NMATCH .LT. 10) GOTO 900
C
C Run the LSQ refining different parameter models
C
	PRINT '(/,1X,A)','Testing Model 1 Refinement'
	CALL SET_LSQ_REFINE(1,1,1,0,0,0)
	CALL RUN_LSQ()
	RMS1=CALC_RMS_ERR(UB)
C
C Save parameters not refined in Model 1
C
	SAVE1=PIX_SIZE(1)
	SAVE2=PIX_SIZE(2)
	SAVE3=PIX_SKEW
	SAVE4=DRUM_OFF(1)
	SAVE5=DRUM_OFF(2)
	SAVE6=DRUM_OFF(3)
	SAVE7=BEAM_VERT
C
	PRINT '(/,1X,A)','Testing Model 2 Refinement'
	CALL SET_LSQ_REFINE(1,1,1,0,1,0)
	CALL RUN_LSQ()
	RMS2=CALC_RMS_ERR(UB)
C
	PRINT '(/,1X,A)','Testing Model 3 Refinement'
	CALL SET_LSQ_REFINE(1,1,1,1,1,1)
	CALL RUN_LSQ()
	RMS3=CALC_RMS_ERR(UB)
C
C Restore saved parameters
C
	PIX_SIZE(1)=SAVE1
	PIX_SIZE(2)=SAVE2
	PIX_SKEW=SAVE3
	DRUM_OFF(1)=SAVE4
	DRUM_OFF(2)=SAVE5
	DRUM_OFF(3)=SAVE6
	BEAM_VERT=SAVE7
C
C Decide which model is best, then make final refinement
C More complex models have to be at least 10% better
C Also, use fewer matches for very bad refinements
C
	IF(MIN(RMS1,RMS2,RMS3) .GT. 1.0) THEN
	  PRINT '(/,1X,A)','Using Model 1 Refinement with fewer spots'
	  NMATCH=MAX(10, NMATCH/2)
	  CALL SET_LSQ_REFINE(1,1,1,0,0,0)
	  CALL RUN_LSQ()
	ELSEIF(RMS1 .LE. MIN(RMS2*1.1,RMS3*1.2)) THEN
	  PRINT '(/,1X,A)','Repeating Model 1 Refinement'
	  CALL SET_LSQ_REFINE(1,1,1,0,0,0)
	  CALL RUN_LSQ()
	ELSEIF(RMS2 .LE. RMS3*1.1) THEN
	  PRINT '(/,1X,A)','Repeating Model 2 Refinement'
	  CALL SET_LSQ_REFINE(1,1,1,0,1,0)
	  CALL RUN_LSQ()
	ELSE
	  PRINT '(/,1X,A)','Repeating Model 3 Refinement'
	  CALL SET_LSQ_REFINE(1,1,1,1,1,1)
	  CALL RUN_LSQ()
	ENDIF
C
	RETURN
C
900	PRINT '(I4,A)',NMATCH,' matches insufficient for refinement'
	CALL QUIT('ERROR: Insufficient matches for refinement')
	END


	SUBROUTINE RUN_MANUAL_MODE(UB,IMODE)
C
C Run the program in its various manual modes
C IMODE=2	manual, full refinement mode
C       3	manual, hklgen and spot matching mode
C       4	manual, hklgen only mode
C
	REAL UB(3,3)
C
	COMMON /MODE_COM/ HKLGEN1,HKLGEN2,DIST_MATCH,IHKLGEN,IREFS(6)
C
	CHARACTER RESULTS*132
	COMMON /RESULTS_COM/ RESULTS
C
	REAL HKLS(3,10000),PIXWAV(3,10000)
C
C Always run hklgen and make a RESULTS string
C
	RESULTS=REPEAT(' ',LEN(RESULTS))
	IF(IHKLGEN .EQ. 1) THEN
	  CALL GEN_HKLS_SIMPLE(HKLS,PIXWAV,NHKLS, UB,HKLGEN1,HKLGEN2)
	  WRITE(RESULTS,'(I4,A)') NHKLS,' hkls generated using simple algorithm'
	ELSE IF(IHKLGEN .EQ. 2) THEN
	  CALL GEN_HKLS_FULL(HKLS,PIXWAV,NHKLS, UB,HKLGEN2,HKLGEN1)
	  WRITE(RESULTS,'(I5,A)') NHKLS,' hkls generated using full algorithm'
	ELSE
	  CALL QUIT('ERROR: Invalid IHKLGEN in input file')
	ENDIF
	PRINT *,TRIM(RESULTS)
C
C Match spots if IMODE=2,3
C
	IF(IMODE .LE. 3) THEN
C Make the matches
	  CALL FIND_MATCH(NMATCH,NOVER, HKLS,PIXWAV,NHKLS, 10000,DIST_MATCH)
C Output results to the console
	  RESULTS=REPEAT(' ',LEN(RESULTS))
	  IF(NOVER .EQ. 0) THEN
	    WRITE(RESULTS,'(I5,A)') NMATCH,' matches found'
	  ELSE
		  WRITE(RESULTS,'(I5,A,I4,A)') NMATCH,' matches found   (',
	1			  NOVER,' ambiguous matches were rejected)'
	  ENDIF
	  PRINT *,TRIM(RESULTS)
	ENDIF
C
C Run the LSQ if IMODE=2
C
	IF(IMODE .EQ. 2) THEN
	  PRINT *
C Load cell constraints and setup UB for refinement
	  CALL INIT_CELL_REFINE()
	  CALL SETUP_UB_MATRIX(UB)
C Turn on the requested refinements
	  CALL SET_LSQ_REFINE(IREFS(1),IREFS(2),IREFS(3),
	1			IREFS(4),IREFS(5),IREFS(6))
C Run the refinement
	  CALL RUN_LSQ()
	ENDIF
C
	RETURN
	END


	SUBROUTINE RUN_HKL_MODE(UB)
C
C Run the program in the HKL mode
C
	REAL UB(3,3)
C
C Load cell constraints and setup UB for refinement
C
	CALL INIT_CELL_REFINE()
	CALL SETUP_UB_MATRIX(UB)
C
C Turn on the refinement of cell and rotations, only
C
	CALL SET_LSQ_REFINE(0,1,1,0,0,0)
C
C Run the refinement
C
	CALL RUN_LSQ()
C
	RETURN
	END


	SUBROUTINE SETUP_UB_MATRIX(UB)
C
C Calculate the starting (base) U matrix from the UB and the
C cell dimensions. Also, set the U matrix rotations to zero.
C
	REAL UB(3,3)
C
	COMMON /CELL_COM/ CELL(6),ILATT,ICEN,CELL_PROD
	COMMON /U_COM/ U_BASE(3,3),RXYZ(3)
C
	REAL B(3,3),B_INV(3,3)
C
	CALL CALC_B_MATRIX(B, CELL)
	CALL MX3_INVERT(B_INV,DET, B)
	CALL MX3_MULT(U_BASE, UB,B_INV)
	CALL MAKE_U_ORTHOGONAL(U_BASE)
C
	CALL VC3_SET(RXYZ, 0.0,0.0,0.0)
C
	RETURN
	END
