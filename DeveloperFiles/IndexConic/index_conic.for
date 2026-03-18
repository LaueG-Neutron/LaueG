	PROGRAM INDEX_CONIC
C
	INTEGER IHK(2,1000)
	REAL HVECS(3,10000),V_NODAL1(3),V_NODAL2(3),V_CONIC(3),WAVS(1000),PIXWAV(3)
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS
C
C Output a simple banner
C
	PRINT '(1X,A)','Conic indexing program for LaueG (Ross Piltz, 9/9/2024)'
	PRINT *
C
C Delete the output file, if it exists
C
	OPEN(UNIT=1,FILE='___laueg_index_conic.out',STATUS='OLD',IOSTAT=IDUMMY)
	CLOSE(UNIT=1,STATUS='DELETE',IOSTAT=IDUMMY)
C
C Read the input file
C
	CALL READ_IN_FILE(V_CONIC,INODAL1,INODAL2)
C
C Prune SPOTS(1..3,1..NSPOTS) to only include those spots on
C the conic and move the nodals spots to spots 1 & 2
C
	CALL PRUNE_CONIC_SPOTS(V_CONIC,INODAL1,INODAL2)
	NCONIC=NSPOTS
	PRINT '(I5,1X,A)',NCONIC,'spots lie on the conic line'
C
C Calculate scattering unit-vectors for spots on the conic
C
	DO I=1,NCONIC
		PIXWAV(1)=SPOTS(1,I)
		PIXWAV(2)=SPOTS(2,I)
		PIXWAV(3)=1.0
	  CALL CALC_PIXWAVS_TO_HVECS(HVECS(1,I), PIXWAV,1)
	  CALL VC3_UNIT(HVECS(1,I))
	ENDDO
C
C Copy scattering unit-vectors for the 2 nodal spots
C
	DO I=1,3
	  V_NODAL1(I)=HVECS(I,1)
	  V_NODAL2(I)=HVECS(I,2)
	ENDDO
C
C Do initial indexing using nodal vectors
C
	CALL INDEX_NODALS(V_NODAL1,V_NODAL2, HVECS,NCONIC)
C
C Scale relative lengths of nodal vectors by simple fractions
C to obtain the set of smallest H & K indices
C
	CALL SCALE_NODALS_MIN_HK(V_NODAL1,V_NODAL2, HVECS,NCONIC)
C
C Calculate H & K indices, IHK, for nodal vectors
C
	CALL CALC_INDEX_HK(IHK, V_NODAL1,V_NODAL2, HVECS,NCONIC)
C
C Scale both nodal vectors to give reasonable minimum wavelengths
C for the spots on the conic. Also calculate the spot wavelengths
C using the scaled nodal vectors.
C
	CALL SCALE_NODALS_WAVS(WAVS,IHK, V_NODAL1,V_NODAL2, NCONIC)
C
C Write output file for LaueG to read in the results
C
	OPEN(UNIT=2,FILE='___laueg_index_conic.out',STATUS='UNKNOWN')
	WRITE(2,'(3F8.4)') V_NODAL1
	WRITE(2,'(3F8.4)') V_NODAL2
C
	DO I=1,NCONIC
	  IF(WAVS(I) .LT. 1E5) WRITE(2,'(2F7.1,2I4)') (SPOTS(K,I),K=1,2),(IHK(K,I),K=1,2)
	ENDDO
C
	CLOSE(UNIT=2)
C
C Delete the input file
C
	CALL DELETE_FILE('___laueg_index_conic.in')
C
	PRINT *,'SUCCESSFUL COMPLETION'
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



	SUBROUTINE PRUNE_CONIC_SPOTS(V_CONIC,INODAL1,INODAL2)
C
	REAL V_CONIC(3)
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS
C
	LOGICAL LOK(10000)
	REAL PIXWAV(3),HVEC(3),SPOTS2(3,1000)
C
C Normalise the conic vector
C
	CALL VC3_UNIT(V_CONIC)
C
C Set LOK() for all spots within CUTOFF pixels of the conic line
C
	CUTOFF=10.0
	DO I=1,NSPOTS
C
C Calculate scattering unit-vectors for all observed spots
C
		PIXWAV(1)=SPOTS(1,I)
		PIXWAV(2)=SPOTS(2,I)
		PIXWAV(3)=1.0
	  CALL CALC_PIXWAVS_TO_HVECS(HVEC, PIXWAV,1)
C
C Move scattering vectors onto the conic line and normalise them
C
	  CALL MAKE_UNIT_PERP(HVEC, HVEC,V_CONIC)
C
C Calculate pixel position from the "moved" scattering vectors
C
	  CALL CALC_HVECS_TO_PIXWAVS(PIXWAV, HVEC,1)
C
C Remove from PIXELS() any spots where the pixel position moved >10
C pixels. Use ISPOTS() to keep track of the original SPOT() index.
C
	  DX=SPOTS(1,I)-PIXWAV(1)
	  DY=SPOTS(2,I)-PIXWAV(2)
	  DSQ=DX**2 + DY**2
	  LOK(I)=(DSQ .LT. CUTOFF**2)
C
	ENDDO
C
C Copy INODAL1 & 2 from SPOTS() to SPOTS2()
C
	DO I=1,3
	  SPOTS2(I,1)=SPOTS(I,INODAL1)
	  SPOTS2(I,2)=SPOTS(I,INODAL2)
	ENDDO
C
C Copy remaining OK from SPOTS() to SPOTS2()
C
	LOK(INODAL1)=.FALSE.
	LOK(INODAL2)=.FALSE.
	NOK=2
	DO I=1,NSPOTS
	  IF( LOK(I) ) THEN
	    NOK=NOK+1
	    SPOTS2(1,NOK)=SPOTS(1,I)
	    SPOTS2(2,NOK)=SPOTS(2,I)
	  ENDIF
	ENDDO
	NSPOTS=NOK
C
C Copy SPOTS2() back to SPOTS()
C
	DO I=1,NSPOTS
	  SPOTS(1,I)=SPOTS2(1,I)
	  SPOTS(2,I)=SPOTS2(2,I)
	ENDDO
C
	RETURN
	END



	SUBROUTINE MAKE_UNIT_PERP(VOUT, VIN,VPERP)
C
	REAL VOUT(3),VIN(3),VPERP(3)
C
	DOT=VIN(1)*VPERP(1) +VIN(2)*VPERP(2) +VIN(3)*VPERP(3)
C
	VOUT(1)=VIN(1)-DOT*VPERP(1)
	VOUT(2)=VIN(2)-DOT*VPERP(2)
	VOUT(3)=VIN(3)-DOT*VPERP(3)
C
	RETURN
	END



	SUBROUTINE READ_IN_FILE(V_CONIC,INODAL1,INODAL2)
C
	REAL V_CONIC(3)
C
	COMMON /CELL_COM/ CELL(6),ILATT,ICEN,CELL_PROD
	COMMON /PIXEL_COM/ PIX_CEN(2),PIX_SIZE(2),PIX_SKEW
	COMMON /OFFSETS_COM/ XTAL_OFF(3),DRUM_OFF(3),BEAM_VERT
	COMMON /GEOM_COM/ DRUM_RAD
C
	COMMON /SPOTS_COM/ SPOTS(3,10000),NSPOTS
C
	OPEN(UNIT=1,FILE='___laueg_index_conic.in',STATUS='OLD',ERR=900)
C
C Read in instrument geometry parameters
C
	READ(1,*,ERR=901) PHI
	READ(1,*,ERR=901) NUMX,NUMY
	READ(1,*,ERR=901) PIX_CEN
	READ(1,*,ERR=901) PIX_SIZE
	READ(1,*,ERR=901) PIX_SKEW
	READ(1,*,ERR=901) XTAL_OFF
	READ(1,*,ERR=901) BEAM_VERT
	READ(1,*,ERR=901) DRUM_OFF
	READ(1,*,ERR=901) DRUM_RAD
C
C Read in conic vector and nodal indices
C
	READ(1,*,ERR=901) V_CONIC
	READ(1,*,ERR=901) INODAL1,INODAL2
C
C Read in the found spot list
C
	READ(1,*,ERR=901) NSPOTS
	DO I=1,NSPOTS
	  READ(1,*,ERR=901) (SPOTS(K,I),K=1,3)
	ENDDO
C
C Sanity check on the nodal indices
C
	IF( MIN(INODAL1,INODAL2).LT.1 .OR.
	1	MAX(INODAL1,INODAL2).GT.NSPOTS )
	2		CALL QUIT('ERROR: Nodal indices invalid')
C
C Brief output of spots read in, then return
C
	CLOSE(UNIT=1)
	PRINT '(1X,I5,1X,A)',NSPOTS,'spots read from input file'
	RETURN
C
900	CALL QUIT('ERROR: Unable to open input file')
901	CALL QUIT('ERROR: Unable to read input file')
	END



	SUBROUTINE DELETE_IN_FILE()
C
	OPEN(UNIT=1,FILE='___laueg_index_conic.in',STATUS='UNKNOWN')
	CLOSE(UNIT=1,STATUS='DELETE')
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
