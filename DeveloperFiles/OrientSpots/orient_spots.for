C ==================================================
C	PROGRAM ORIENT_SPOTS
C	SUBROUTINE READ_IN_FILE(UB,IMODE,NSPOTS2)
C	SUBROUTINE OUTPUT_RESULTS(IMODE)
C	SUBROUTINE DELETE_FILE(FILE_NAME)
C	SUBROUTINE QUIT(TEXT)
C ==================================================

	PROGRAM ORIENT_SPOTS
C
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
	COMMON /MODULATE_COM/ NHKLS_MOD,DMIN_MOD,IMOD_STRONG,RHKLS_MOD(3,100)
	COMMON /SPECIALS_COM/ SPECS(3,100),NSPECS
C
	LOGICAL LCEN_FIXED
	REAL UB(3,3)
C
C Output a simple banner
C
	PRINT '(1X,A)','Spot orienting program for LaueG (Ross Piltz, 5/9/2024)'
	PRINT *
C
C Turn on X-offset trick to reduce LSQ correlations
C
	CALL XTAL_OFFX_ENABLE()
C
C Turn off modulations and HKL special conditions
C
	NHKLS_MOD=0
	NSPECS=0
C
	CALL READ_IN_FILE(UB,IMODE,NSPOTS)
	PRINT '(1X,A,I5,A)','Read in spot list of',NSPOTS,' spots'
C
C Output instrument type
C
	IF(ITYPE .EQ. 1) THEN
	  PRINT '(1X,A)','Type 1: Drum IP, unstable XY, nonzero basecounts'
	ELSEIF(ITYPE .EQ. 2) THEN
	  PRINT '(1X,A)','Type 2: Drum IP, stable XY, zero basecounts'
	ELSEIF(ITYPE .EQ. 3) THEN
	  PRINT '(1X,A)','Type 3: Octagonal CCD'
	ELSEIF(ITYPE .EQ. 4) THEN
	  PRINT '(1X,A)','Type4 : Drum IP, stable XY, nonzero basecounts'
	ELSE
	  CALL QUIT('BUG(gen_hkls): Unknown Instrument Type')
	ENDIF
C
C Output instrument type and program mode, then start that mode
C IMODE=1	automatic mode
C       2	manual mode: hklgen, spot matching, and refinement
C       3	manual mode: hklgen and spot matching
C       4	manual mode: hklgen only
C       5	HKL mode
C
	IF(IMODE .EQ. 1) THEN
	  PRINT *,'Running in Automatic mode'
		LCEN_FIXED=.TRUE.
		IF(ITYPE .EQ. 1)	LCEN_FIXED=.FALSE.
	  CALL RUN_AUTO_MODE(UB, LCEN_FIXED)
	ELSEIF(IMODE .LE. 4) THEN
	  PRINT *,'Running in Manual mode'
	  CALL RUN_MANUAL_MODE(UB,IMODE)
	ELSEIF(IMODE .EQ. 5) THEN
	  PRINT *,'Running in HKL mode'
	  CALL RUN_HKL_MODE(UB)
	ENDIF
C
	CALL OUTPUT_RESULTS(IMODE)
C
	CALL DELETE_FILE('___laueg_orient_spots.in')
C
	PRINT *,'SUCCESSFUL COMPLETION'
	END


	SUBROUTINE READ_IN_FILE(UB,IMODE,NSPOTS2)
C
	REAL UB(3,3)
C
	COMMON /CELL_COM/ CELL(6),ILATT,ICEN,CELL_PROD
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
C
	COMMON /HKL_MODE_COM/ IMODE5,CELL0(6)
	COMMON /MODE_COM/ HKLGEN1,HKLGEN2,DIST_MATCH,IHKLGEN,IREFS(6)
C
	OPEN(UNIT=1,FILE='___laueg_orient_spots.in',STATUS='OLD',ERR=900)
C
C Read the program mode
C
	READ(1,*,ERR=901) IMODE
	IF(IMODE.LT.1 .OR. IMODE.GT.5) CALL QUIT('ERROR: Invalid IMODE in input file')
C
C Parameters common to all modes
C
C Input instrument type: 1=IP, 2=CCD
	READ(1,*,ERR=901) ITYPE
C Input crystallography info + UB
	READ(1,*,ERR=901) ILATT,ICEN
	READ(1,*,ERR=901) CELL
	READ(1,*,ERR=901) (UB(1,K),K=1,3)
	READ(1,*,ERR=901) (UB(2,K),K=1,3)
	READ(1,*,ERR=901) (UB(3,K),K=1,3)
C Input info probably used by all instruments
	READ(1,*,ERR=901) PHI
	READ(1,*,ERR=901) NUMX,NUMY
	READ(1,*,ERR=901) PIX_CEN
	READ(1,*,ERR=901) XTAL_OFF
	READ(1,*,ERR=901) DRUM_OFF
	READ(1,*,ERR=901) BEAM_VERT
C Input info that may vary for instruments
	READ(1,*,ERR=901) PIX_SIZE,PIX_SKEW
	READ(1,*,ERR=901) DRUM_RAD
C
C Read manual mode specific parameters and switches
C
	IF(IMODE.GE.2 .AND. IMODE.LE.4) THEN
	  READ(1,*,ERR=901) IHKLGEN,HKLGEN1,HKLGEN2
	  READ(1,*,ERR=901) DIST_MATCH
	  READ(1,*,ERR=901) IREFS
	ENDIF
C
C Read in the found spot list (and HKLs for MODE=5)
C
	READ(1,*,ERR=901) NSPOTS
	IF(IMODE .NE. 5) THEN
	  DO I=1,NSPOTS
	    READ(1,*,ERR=901) (SPOTS(K,I),K=1,3)
	  ENDDO
	ELSE
C For HKL mode, read XY and HKLs directly to the *_MATCH arrays
	  DO I=1,NSPOTS
	    READ(1,*,ERR=901) (XY_MATCH(K,I),K=1,2),(HKL_MATCH(K,I),K=1,3)
	  ENDDO
	  NMATCH=NSPOTS
	ENDIF
C
	CLOSE(UNIT=1)
	NSPOTS2=NSPOTS
C
C Calculate the product of cell lengths (which will be conserved)
C
	CELL_PROD=CELL(1)*CELL(2)*CELL(3)
C
C For orienting spots, rhom(hex) & hex are the same
C NB: R centering is handled by ICEN not ILATT
C
	IF(ILATT .EQ. 8) ILATT=7
C
C Save values for later use in HKL mode
C
	IMODE5=IMODE
	DO I=1,6
	  CELL0(I)=CELL(I)
	ENDDO
C
	RETURN
C
900	CALL QUIT('ERROR: Unable to open input file')
901	CALL QUIT('ERROR: Unable to read input file')
	END


	SUBROUTINE OUTPUT_RESULTS(IMODE)
C
	COMMON /CELL_COM/ CELL(6),ILATT,ICEN,CELL_PROD
	COMMON /U_COM/ U_BASE(3,3),RXYZ(3)
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS
	COMMON /MATCH_COM/ XY_MATCH(2,10000),HKL_MATCH(3,10000),NMATCH
C
	CHARACTER RESULTS*132
	COMMON /RESULTS_COM/ RESULTS
C
	REAL UB(3,3),UB2(3,3)
C
C The cell lengths are scaled by CELL_SCALE so the
C product of cell lengths are conserved (=CELL_PROD)
C
	CELL_SCALE=( CELL_PROD / ( CELL(1)*CELL(2)*CELL(3) ) )**(1.0/3.0)
	CELL(1)=CELL(1)*CELL_SCALE
	CELL(2)=CELL(2)*CELL_SCALE
	CELL(3)=CELL(3)*CELL_SCALE
C
C Calculate the final rms error
C
	CALL CALC_UB_MATRIX(UB)
	RMS1=CALC_RMS_ERR(UB)
C
	OPEN(UNIT=2,FILE='___laueg_orient_spots.out',STATUS='UNKNOWN')
C
C Write the cell & beam geometry info
C
	WRITE(2,'(3F10.5,3F10.4)') CELL
	CALL CALC_UB_MATRIX(UB)
	WRITE(2,'(3F12.7)') (UB(1,K),K=1,3)
	WRITE(2,'(3F12.7)') (UB(2,K),K=1,3)
	WRITE(2,'(3F12.7)') (UB(3,K),K=1,3)
	WRITE(2,'(2F10.3)') PIX_CEN
	WRITE(2,'(2F10.6)') PIX_SIZE
	WRITE(2,'(F10.6)') PIX_SKEW
	WRITE(2,'(2F10.6)') XTAL_OFF(1),XTAL_OFF(3)
	WRITE(2,'(F10.6)') BEAM_VERT
C
C If an auto, or a manual full refinement, make a new RESULTS string 
C
	IF(IMODE .LE. 2) THEN
C Reread parameters so we can calculate the "before" rms
	  CALL READ_IN_FILE(UB2,IDUM,NDUM)
	  CALL SETUP_UB_MATRIX(UB2)
	  CALL CALC_UB_MATRIX(UB2)
	  RMS2=CALC_RMS_ERR(UB2)
	  RESULTS=REPEAT(' ',LEN(RESULTS))
C In auto mode write just the number, for manual mode make it pretty
	  IF(IMODE .EQ. 1) THEN
	    WRITE(RESULTS,'(2I8,2F8.2)') NSPOTS,NMATCH,RMS2,RMS1
	  ELSE
	    WRITE(RESULTS,'(A,F6.3,A,3X,A,F6.3,A)') 'Final rms error =',
	1			RMS1,' mm','(initial orientation ',RMS2,' mm)'
	  ENDIF
	ENDIF
C
C Write the results string to the output file, then close the file
C
	WRITE(2,*) TRIM(RESULTS)
	CLOSE(UNIT=2)
C
C If an auto, or a manual full refinement, output a summary
C
	IF(IMODE .LE. 2) THEN
	  PRINT '(/,1X,A,F6.3,A,3X,A,F6.3,A)','Final rms error =',
	1		RMS1,' mm','(initial orientation ',RMS2,' mm)'
	  PRINT '(/,1X,A,3(/,8X,3F10.6),/)','Refined UB matrix:',
	1				((UB(K2,K1),K1=1,3),K2=1,3)
	ENDIF
C
	RETURN
	END



	SUBROUTINE DELETE_FILE(FILE_NAME)
C
C Delete the file if the "delete files" is missing or contains a TRUE
C
	CHARACTER FILE_NAME*(*)
C
	LOGICAL LDELETE
C
	OPEN(UNIT=17,FILE='___laueg_delete_files.in',STATUS='OLD',ERR=100)
	READ(17,*,ERR=100,END=100) LDELETE
	IF( .NOT.LDELETE ) RETURN
100	CLOSE(UNIT=17,IOSTAT=IDUMMY)
C
	OPEN(UNIT=17,FILE=FILE_NAME,STATUS='OLD',IOSTAT=IDUMMY)
	CLOSE(UNIT=17,STATUS='DELETE',IOSTAT=IDUMMY)
C
	RETURN
	END


	SUBROUTINE QUIT(TEXT)
C
C Workaround as SCILAB doesn't return STOP messages as they
C go to stderr not stdout
C
	CHARACTER*(*) TEXT
C
	PRINT *,TEXT
	CALL EXIT()
C
	END
