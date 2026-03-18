	PROGRAM INDEX_SPOTS
C
	REAL CELL(6),B(3,3),UB(3,3),XYCEN(2)
C
C Output a simple banner
C
	PRINT '(1X,A)','Spot indexing program for LaueG (Ross Piltz, 9/9/2024)'
	PRINT *
C
C Delete the output file, if it exists
C
	OPEN(UNIT=1,FILE='___laueg_index_spots.out',STATUS='OLD',IOSTAT=IDUMMY)
	CLOSE(1,STATUS='DELETE',IOSTAT=IDUMMY)
C
C Read parameters & spots list (delete input file after reading)
C
	PRINT *,'Reading input file'
	CALL READ_INPUT_FILE(CELL,XCEN0,YCEN0,NPEAKS_MATCH,NPEAKS_SOLN)
C
C Sort the observed spots in order of decreasing "strength", except
C all "guide spots" are placed at the top of the list
C
	CALL SORT_OBS_SPOTS
C
C Calculate the reciprocal space vectors for all spots
C
	CALL CALC_PEAK_VECTORS(NPEAKS_SOLN)
C
C Calculate B matrix from the cell parameters
C
	CALL CALC_B_MATRIX(B, CELL)
C
C Calculate a list of scattering vectors H_HKLS for low order hkls
C
	PRINT '(/,1X,A)','Matching observed and HKL generated spots'
	CALL GENERATE_HKL_VECTORS(B)
C
C Make a table of the angles between all vectors in H_HKLS
C
	CALL MAKE_HKL_ANGLES_TABLE
C
C Find all triplets of spots that match a triplet of hkls.
C Use the first NPEAKS_MATCH observed spots of the sorted list.
C
	CALL FIND_TRIPLETS(NPEAKS_MATCH)
C
C Remove equivalent triplet solutions and find the best solution
C
	PRINT '(/,1X,A)','Testing for valid solutions'
	CALL FIND_SOLUTIONS(UB,XYCEN, B,NPEAKS_MATCH,NPEAKS_SOLN)
C
C The output file should have been deleted, so crash if it exists
C
	OPEN(UNIT=1,FILE='___laueg_index_spots.out',STATUS='UNKNOWN')
	WRITE(1,'(2F8.1)') XYCEN
	WRITE(1,'((3F10.5))') ((UB(K1,K2),K2=1,3),K1=1,3)
	CLOSE(UNIT=1)
C
C Output the change in pixel centers
C
	PRINT '(1X,A,4(I4,A))','Pixel center changed from (' ,NINT(XCEN0),
	1		',',NINT(YCEN0),') to (',NINT(XYCEN(1)),',',NINT(XYCEN(2)),')'
C
C Delete the input file
C
	CALL DELETE_FILE('___laueg_index_spots.in')
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



	SUBROUTINE READ_INPUT_FILE(CELL,XCEN2,YCEN2,NPEAKS_MATCH,NPEAKS_SOLN)
C
	REAL CELL(6)
C
	PARAMETER (NSPOTS_MAX=10000)
	COMMON /SPOTS_COM/ NSPOTS,XSPOT(NSPOTS_MAX),YSPOT(NSPOTS_MAX),ESDS(NSPOTS_MAX)
C
	PARAMETER (NPEAKS_MAX=100)
C
	COMMON /GUIDE_SPOTS_COM/ NGUIDES,IGUIDES(100)
C
	COMMON /INSTRUM_COM/ ITYPE,NUMXY(2),XY_CEN(2),XY_SIZE(2),DRUM_RAD
C
C IHSQ_MAX IS max H*H+K*K+L*L used to generate hkl list
C PAIR_ANG_MIN is the min angle between spots within a pair
C PAIR_DANG_MAX is the max +/- tolerance in angles between obs & calc pairs
C AZI_DANG_MAX is the max difference in azimuthal angle for a triplet group
C MIN_AZI_GROUP is the min number of pairs in a triplet group
C DOT_CROSS_MIN is the min dot-cross between the spots in a triplet
C IHKL_MAX1 is the max abs HKL for solutions in the initial tests
C IHKL_MAX2 is the max abs HKL for solutions in the final tests
C ISOLN_SELECT is the number of the highest merit solution to output
C NPEAKS_MATCH is the max number of peaks to use in vector matching
C NPEAKS_SOLN is the max number of peaks to use when testing solutions
C
	COMMON /HKLGEN_COM/ IHSQ_MAX
	COMMON /PAIR_MATCH_COM/ PAIR_ANG_MIN,PAIR_DANG_MAX
	COMMON /AZI_MATCH_COM/ AZI_DANG_MAX,MIN_AZI_GROUP
	COMMON /TRIP_MATCH_COM/ DOT_CROSS_MIN
	COMMON /HKL_BEST_COM/ IHKL_MAX1,IHKL_MAX2
	COMMON /SOLN_SELECT_COM/ ISOLN_SELECT
C
	OPEN(UNIT=1,FILE='___laueg_index_spots.in',STATUS='OLD',IOSTAT=IERR)
	IF(IERR .NE. 0) CALL QUIT('ERROR: Unable to open input file')
C
C Read cell dimensions
C
	READ(1,*) CELL
C
C Read the instrument type and its parameters
C
	READ(1,*) ITYPE
	READ(1,*) NUMXY
	READ(1,*) XY_CEN
	READ(1,*) XY_SIZE
	READ(1,*) DRUM_RAD
C
C Read the indexing specific parameters
C
	READ(1,*) IHSQ_MAX
	READ(1,*) PAIR_ANG_MIN,PAIR_DANG_MAX
	READ(1,*) AZI_DANG_MAX,MIN_AZI_GROUP
	READ(1,*) DOT_CROSS_MIN
	READ(1,*) IHKL_MAX1,IHKL_MAX2
	READ(1,*) ISOLN_SELECT
	READ(1,*) NPEAKS_MATCH,NPEAKS_SOLN
C
C Read in the observed spots
C
	READ(1,*) NSPOTS
	IF(NSPOTS .GT. NSPOTS_MAX) CALL QUIT('ERROR: Too many observed spots')
	DO I=1,NSPOTS
	  READ(1,*) XSPOT(I),YSPOT(I),ESDS(I)
	ENDDO
C
C Read in the guide spots which index the observed spots list
C
	READ(1,*) NGUIDES
	IF(NGUIDES .GT. 100) CALL QUIT('ERROR: Too many guide spots')
	DO I=1,NGUIDES
	  READ(1,*) IGUIDES(I)
	  IF(IGUIDES(I).LT.1 .OR. IGUIDES(I).GT.NSPOTS)
	1		CALL QUIT('ERROR: Invalid guide spot index')
	ENDDO
C
100	CLOSE(1)
C
C Enforce NSPOTS limit on NPEAK_* parameters
C
	NPEAKS_MATCH=MIN(NPEAKS_MATCH,NSPOTS)
	NPEAKS_SOLN=MIN(NPEAKS_SOLN,NSPOTS)
C
C Output a summary to the user
C
	IF(ITYPE .EQ. 3) THEN
		PRINT *,' Detector type: octagonal CCD'
	ELSE
		PRINT *,' Detector type: cylindrical IP'
	ENDIF
	PRINT '(1X,I4,A)',NSPOTS,' observed spots read from file'
	PRINT '(1X,I4,A)',NPEAKS_MATCH,' spots will be used for spot matching'
	PRINT '(1X,I4,A)',NPEAKS_SOLN,' spots will be used for testing solutions'
	IF(MAX(NPEAKS_MATCH,NPEAKS_SOLN) .GT. NPEAKS_MAX)
	1	CALL QUIT('ERROR: Too many matching and/or testing spots')
	IF(NGUIDES .GT. 0) PRINT '(1X,I4,A)',NGUIDES,
	1		' spots are designated as guide spots'
C
	XCEN2=XY_CEN(1)
	YCEN2=XY_CEN(2)
	RETURN
	END



	SUBROUTINE SORT_OBS_SPOTS
C
C Sort the observed spots lists (XSPOT,YSPOT,ESDS) in order of
C decreasing ESDS(), but with any guide spots placed before
C the normal spots.
C NB: ESDS() are invalid after the sort, but they aren't used again
C
	EXTERNAL ESDS_COMPARE
C
	PARAMETER (NSPOTS_MAX=10000)
	COMMON /SPOTS_COM/ NSPOTS,XSPOT(NSPOTS_MAX),YSPOT(NSPOTS_MAX),ESDS(NSPOTS_MAX)
C
	COMMON /GUIDE_SPOTS_COM/ NGUIDES,IGUIDES(100)
C
	
	INTEGER ITAGS(NSPOTS_MAX)
	REAL TEMP(3,NSPOTS_MAX)
C
C Change guide spot ESDS so they are the most important
C
      DO I=1,NGUIDES
        ESDS(IGUIDES(I))=ESDS(IGUIDES(I))+1E6
      ENDDO
C
C Do a tagged sort on the ESDS() values
C
	CALL TAG_SORT(ESDS_COMPARE,ITAGS,NSPOTS)
C
C Sort X,Y & esds in descending order of esds
C
	DO I=1,NSPOTS
	  TEMP(1,I)=XSPOT(ITAGS(I))
	  TEMP(2,I)=YSPOT(ITAGS(I))
	  TEMP(3,I)=ESDS(ITAGS(I))
	ENDDO
C
	DO I=1,NSPOTS
	  XSPOT(I)=TEMP(1,I)
	  YSPOT(I)=TEMP(2,I)
	  ESDS(I) =TEMP(3,I)
	ENDDO
C
C Change guide spot ESDS back to the original values
C
      DO I=1,NGUIDES
        ESDS(I)=ESDS(I)-1E6
      ENDDO
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
