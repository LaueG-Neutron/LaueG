C	PROGRAM GEN_HKLS
C	SUBROUTINE MAIN()
C	SUBROUTINE READ_INPUT_FILE(UB,ICEN, FILE_NAME)
C ================= Misc. routines =============
C	SUBROUTINE DELETE_FILE(FILE_NAME)
C	SUBROUTINE QUIT(TEXT)


	PROGRAM GEN_HKLS
C
C Output a simple banner
C
	PRINT '(1X,A,/)','HKL generation program for LaueG (Ross Piltz, 4/9/2024)'
C
C Start it up
C
	CALL MAIN()
C
C LaueG needs this output to check the program didn't crash
C
	PRINT *
	PRINT *,'SUCCESSFUL COMPLETION'
	END


	SUBROUTINE MAIN()
C
C Generate a file of unique hkls which strike within the detector limits.
C The file contains hkl, xy, wav, d-space, multiples and satellite number.
C A satellite number of zero means the "normal" (or "main") reflection.
C	
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
C
	COMMON /LIMITS_COM/ D_MIN,D_MAX,WAV_MIN,WAV_MAX,X_MIN,X_MAX,Y_MIN,Y_MAX
C
	LOGICAL LMOD_WEAK
	COMMON /MODULATE_COM/ NHKLS_MOD,DMIN_MOD,LMOD_WEAK,RHKLS_MOD(3,100)
C
	REAL UB(3,3)
C
	PARAMETER (NHKLS_MAX=300000)
	INTEGER IHKLS(3,NHKLS_MAX),IMODS(NHKLS_MAX),IMULTS(NHKLS_MAX)
	REAL HKLS(3,NHKLS_MAX),DSPACES(NHKLS_MAX),PIXWAVS(3,NHKLS_MAX)
C
C Delete output file, if it still exists
C
	OPEN(UNIT=1,FILE='___laueg_gen_hkls.out',STATUS='OLD',IOSTAT=IDUMMY)
	CLOSE(1,STATUS='DELETE',IOSTAT=IDUMMY)
C
C Read parameters from the input file (mostly copied to COMMON)
C
	CALL READ_INPUT_FILE(UB,ICEN, '___laueg_gen_hkls.in')
C
C Increase D_MIN if it is impossible due to sin(TH) > 1
C
	D_MIN=MAX(D_MIN,WAV_MIN/2.0)
C
C Turn on the X-offset trick
C
	CALL XTAL_OFFX_ENABLE()
C
C Generate hkl list according to symmetry and within d-spacing limits
C NB: Only one of the Friedel pair is generated
C
	CALL GEN_HEMISPHERE(HKLS,NHKLS, UB,ICEN,D_MIN,D_MAX)
	PRINT '(/,I8,A)',NHKLS,' hkls generated within d-space limits'
C
C For modulations create satellite HKLs, otherwise zero modulation indices
C
	IF(NHKLS_MOD .GT. 0) THEN
	  CALL ADD_SAT_HKLS(HKLS,IMODS,NHKLS, UB,D_MIN,D_MAX)
	  PRINT '(I8,A)',NHKLS,' hkls after generating satellites'
	ELSE
	  DO I=1,NHKLS
	    IMODS(I)=0
	  ENDDO
	ENDIF
C
C Store main HKLs as integers in IHKLS(), add modulation fractions to HKLS()
C
	DO I=1,NHKLS
	  IHKLS(1,I)=NINT(HKLS(1,I))
	  IHKLS(2,I)=NINT(HKLS(2,I))
	  IHKLS(3,I)=NINT(HKLS(3,I))
	  IF(IMODS(I) .GT. 0) THEN
	    HKLS(1,I)=HKLS(1,I)+RHKLS_MOD(1,IMODS(I))
	    HKLS(2,I)=HKLS(2,I)+RHKLS_MOD(1,IMODS(I))
	    HKLS(3,I)=HKLS(3,I)+RHKLS_MOD(1,IMODS(I))
	  ENDIF
	ENDDO
C
C Prune to X, Y, and wavelength limits
C NB: Switches to Friedel (-h,-k,-l) if wavelength is negative
C
	CALL PRUNE_XY_WAV(HKLS,IHKLS,IMODS,NHKLS, UB, X_MIN,X_MAX,Y_MIN,Y_MAX, WAV_MIN,WAV_MAX)
	PRINT '(I8,A)',NHKLS,' hkls after pruning to X, Y, and wavelength limits'
C
C Prune/merge multiples of spots and return multiplicity in IMULTS()
C Satellite spots are treated according to the variable LMOD_WEAK
C
	CALL PRUNE_MULTIPLES(HKLS,IHKLS,IMODS,IMULTS,NHKLS, UB)
	PRINT '(I8,A)',NHKLS,' hkls after merging multiples'
C
C Calculate final pixel X,Y, wavelength
C
	CALL CALC_HKLS_TO_PIXWAVS(PIXWAVS, UB,HKLS,NHKLS)
C
C Calculate d-spaces from hkls
C NB: Using unrotated UB is OK for d-spaces
C
	DO I=1,NHKLS
	  CALL CALC_HKL_TO_HVEC(HVX,HVY,HVZ, UB,HKLS(1,I),HKLS(2,I),HKLS(3,I))
	  DSPACES(I)=1.0/SQRT(HVX**2 + HVY**2 + HVZ**2)
	ENDDO
C
C Write reflection list to output file
C
	OPEN(UNIT=2,FILE='___laueg_gen_hkls.out',STATUS='UNKNOWN',IOSTAT=IERR)
	IF(IERR .NE. 0) CALL QUIT('ERROR: Unable to open output file')
	DO I=1,NHKLS
	  WRITE(2,'(3I5,F8.3,2F8.1,F8.3,2I4,3F8.3)') (IHKLS(K,I),K=1,3),
	1			DSPACES(I),(PIXWAVS(K,I),K=1,3),IMULTS(I),IMODS(I),(HKLS(K,I),K=1,3)
	ENDDO
	CLOSE(2)
C
C Delete the input file and return
C
	CALL DELETE_FILE('___laueg_gen_hkls.in')
	RETURN
	END

	SUBROUTINE CALC_DSPACES(DSPACES, UB,HKLS,NHKLS)
C
	REAL DSPACES(*),UB(3,3),HKLS(3,*)
C
	RETURN
	END



	SUBROUTINE READ_INPUT_FILE(UB,ICEN, FILE_NAME)
C
	CHARACTER FILE_NAME*(*)
	REAL UB(3,3)
C
C Generate a file of unique hkls which strike within the detector limits.
C The file contains hkl, xy, wav, d-space, multiples and satellite number.
C A satellite number of zero means the "normal" reflection.
C
	COMMON /GEOM_COM/ ITYPE,NUMX,NUMY,DRUM_RAD,PHI
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
C	
	COMMON /LIMITS_COM/ D_MIN,D_MAX,WAV_MIN,WAV_MAX,X_MIN,X_MAX,Y_MIN,Y_MAX
C
	LOGICAL LMOD_WEAK
	COMMON /MODULATE_COM/ NHKLS_MOD,DMIN_MOD,LMOD_WEAK,RHKLS_MOD(3,100)
C
	COMMON /SPECIALS_COM/ SPECS(3,100),NSPECS
C
	CHARACTER*4 STR_CEN(7)
	DATA STR_CEN/'none','A   ','B   ','C   ','I   ','F   ','R   '/
C
C Read parameters from the input file, then delete the file
C
	OPEN(UNIT=1,FILE=FILE_NAME,STATUS='OLD',IOSTAT=IERR)
	IF(IERR .NE. 0) CALL QUIT('ERROR: Unable to open output file')
C
C Instrument type: 1=KOALA1, 2=IMAGINE, 3=CYCLOPS, 4=KOALA2
C

	READ(1,*) ITYPE
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
C Centering type: 1-6 for A,B,C,I,F,R centering, =0 for none
C
	READ(1,*) ICEN
C
C D-spacing and wavelength limits for all spots
C
	READ(1,*) WAV_MIN,WAV_MAX,D_MIN,D_MAX
C
C Various geometry parameters
C
	READ(1,*) (UB(1,K),K=1,3)
	READ(1,*) (UB(2,K),K=1,3)
	READ(1,*) (UB(3,K),K=1,3)
	READ(1,*) PHI
	READ(1,*) X_MIN,X_MAX,Y_MIN,Y_MAX
	READ(1,*) PIX_CEN
	READ(1,*) XTAL_OFF
	READ(1,*) DRUM_OFF
	READ(1,*) BEAM_VERT
	READ(1,*) PIX_SIZE,PIX_SKEW
	READ(1,*) DRUM_RAD
C
C Read list of special conditions for HKLs
	READ(1,*) NSPECS
	IF(NSPECS .GT. 100) CALL QUIT('ERROR: >100 special conditions')
	DO I=1,NSPECS
	  READ(1,'(3I1)') I1,I2,I3
	  SPECS(1,I)=I1
	  SPECS(2,I)=I2
	  SPECS(3,I)=I3
	ENDDO
C
C Read list of modulation satellites
	LMOD_WEAK=.FALSE.
	READ(1,*) NHKLS_MOD
	IF(NHKLS_MOD .GT. 100) CALL QUIT('ERROR: >100 modulation satellites')
	IF(NHKLS_MOD .GT. 0) THEN
	  READ(1,*) ((RHKLS_MOD(K1,K2),K1=1,3),K2=1,NHKLS_MOD)
	  READ(1,*) DMIN_MOD
	ENDIF
C
	CLOSE(1)
C
	PRINT '(1X,A)','Instrument limits:'
	PRINT '(11X,2(A,F6.2),/,11X,2(A,F6.2))','Wavelength =',
	1		WAV_MIN,' to',WAV_MAX,'d-space =',D_MIN,' to',D_MAX
	PRINT '(11X,2(A,I5))','X =',NINT(X_MIN),' to',NINT(X_MAX)
	PRINT '(11X,2(A,I5))','Y =',NINT(Y_MIN),' to',NINT(Y_MAX)
	PRINT '(/,1X,2A)','Centering = ',STR_CEN(ICEN+1)
C
	IF(NSPECS .EQ. 0) THEN
	  PRINT '(1X,A)','No HKL special conditions'
	ELSE
	  PRINT '(I4,1X,A)',NSPECS,'HKL special conditions (3-digit format):'
	  PRINT '(10X,20(3I1,A1))',((NINT(SPECS(K1,K2)),K1=1,3),' ',K2=1,NSPECS)
	ENDIF
C
	IF(NHKLS_MOD .EQ. 0) THEN
	  PRINT '(/,1X,A)','No modulation satellites'
	ELSE
	  PRINT '(/,I4,1X,A)',NHKLS_MOD,'modulation satellite HKLs in list'
	  PRINT '(1X,A,F6.2)','Satellites added to main spots with d-space >',DMIN_MOD
	ENDIF
C
	RETURN
	END


C ================= Misc. routines =============

	SUBROUTINE DELETE_FILE(FILE_NAME)
C
C Delete the file if the "delete files" is missing or contains a TRUE
C
	CHARACTER FILE_NAME*(*)
C
	LOGICAL LDELETE
C
	LDELETE=.TRUE.
	OPEN(UNIT=17,FILE='___laueg_delete_files.in',STATUS='OLD',ERR=100)
	LDELETE=.FALSE.
	READ(17,*,ERR=100,END=100) LDELETE
100	CLOSE(UNIT=17,IOSTAT=IDUMMY)
C
	IF( LDELETE ) THEN
	  OPEN(UNIT=17,FILE=FILE_NAME,STATUS='OLD',IOSTAT=IDUMMY)
	  CLOSE(UNIT=17,STATUS='DELETE',IOSTAT=IDUMMY)
	ENDIF
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
