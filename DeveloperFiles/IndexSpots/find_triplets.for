	SUBROUTINE FIND_TRIPLETS(NPEAKS)
C
	PARAMETER (NHKLS_MAX=1000)
	COMMON /HKL_VECTORS_COM/ NHKLS,H_HKLS(3,NHKLS_MAX)
C
	PARAMETER (NPAIRS_MAX=1000)
	INTEGER IPAIRS(2,NPAIRS_MAX)
C
	REAL H_PEAKS
	COMMON /PEAKS_COM/ H_PEAKS(3,10000)
C
	PARAMETER (NTRIP_MAX=100000)
	COMMON /TRIPLETS_COM/ NTRIP,IPEAKS_TRIP(3,NTRIP_MAX),
	1								IHKLS_TRIP(3,NTRIP_MAX)
C
C
C Search for pairs, then triplets, of matching peaks and hkls where
C IPEAK1 & IHKL1 is the first match. The triplet matches are stored
C in IPEAKS_TRIP & IHKLS_TRIP.
C
	NTRIP=0
	NPAIRS_TOTAL=0
C
C Zero the counter for number of rejected triplets due to guide spots
C
	NREJ=0
C
C Loop over possible observed peaks and hkl values for the first
C spots of a possible pairs of spots.
C
	DO IPEAK1=1,NPEAKS
	  DO IHKL1=1,NHKLS
C Find potential "second parts" of the pairs in IPAIRS(1..2,NPAIRS)
	    CALL FIND_PAIR_MATCHES(IPAIRS,NPAIRS,NPAIRS_MAX, IPEAK1,IHKL1, NPEAKS)
	    NPAIRS_TOTAL=NPAIRS_TOTAL+NPAIRS
	    CALL ADD_TRIPLET_MATCHES(IPAIRS,NPAIRS,NREJ, IPEAK1,IHKL1, H_PEAKS)
	  ENDDO
	ENDDO
C
	PRINT '(1X,I7,A)',NPAIRS_TOTAL,' matching pairs of spots'
	IF(NREJ .NE. 0) PRINT '(1X,I7,A)',NREJ,' triplets do not match any guide spots'
	PRINT '(1X,I7,A)',NTRIP,' matching triplets found'
C
	RETURN
	END



	SUBROUTINE ADD_TRIPLET_MATCHES(IPAIRS,NPAIRS,NREJ, IPEAK1,IHKL1, H_PEAKS)
C
C Assuming that the peak IPEAK1 corresponds to the hkl IHKL1, find two other
C peaks and hkls that correspond, add these three peaks and hkls to the
C list of triplet solutions in IPEAKS_TRIP() & IHKL_TRIP(). Each call to
C this routine may generate many (or zero) triplets.
C
	REAL H_PEAKS(3,10000)
C
	PARAMETER (NTRIP_MAX=100000)
	COMMON /TRIPLETS_COM/ NTRIP,IPEAKS_TRIP(3,NTRIP_MAX),
	1								IHKLS_TRIP(3,NTRIP_MAX)
C
	COMMON /HKL_VECTORS_COM/ NHKLS,H_HKLS(3,1000)
C
	PARAMETER (NPAIR_MAX=1000)
	INTEGER IPAIRS(2,NPAIR_MAX)
C
	EXTERNAL AZI_ROT_COMPARE
C
	COMMON /AZI_ROT_COM/ AZI_ROT(NPAIR_MAX)
C
	INTEGER ITAGS(NPAIR_MAX)
C
	CALL CALC_AZIMUTH_ANGLES(AZI_ROT, IPAIRS,NPAIRS,
	1							H_PEAKS,H_HKLS, IPEAK1,IHKL1)
C
C We will add 9 elements below, and ADD_AZIMUTH_MATCHES need 1 more
C
	IF(NPAIRS .GT. NPAIR_MAX+10) CALL QUIT('BUG(add_triplet_matches) NPAIRS too big')
C
C Do a tagged ascending sort on the AZI_ROT() values
C
	CALL TAG_SORT(AZI_ROT_COMPARE,ITAGS,NPAIRS)
C
C Fudge a +/- 180 wraparound by adding 9 elements to lists
C
	DO I=NPAIRS+1,NPAIRS+9
	  ITAGS(I)=I
	  AZI_ROT(I)=AZI_ROT(ITAGS(I-NPAIRS))+360.0
	  IPAIRS(1,I)=IPAIRS(1,ITAGS(I-NPAIRS))
	  IPAIRS(2,I)=IPAIRS(2,ITAGS(I-NPAIRS))
	ENDDO
	NPAIRS=NPAIRS+9
C
C Search for groups of AZI_ROT() values that are similar
C
	CALL ADD_AZIMUTH_MATCHES(IPEAKS_TRIP,IHKLS_TRIP,NTRIP,NREJ,
	1			AZI_ROT, IPAIRS,NPAIRS,ITAGS, H_HKLS, IPEAK1,IHKL1)
C
	RETURN
	END



	SUBROUTINE ADD_AZIMUTH_MATCHES(IPEAKS_TRIP,IHKLS_TRIP,NTRIP,NREJ,
	1				AZI_ROT, IPAIRS,NPAIRS,ITAGS, H_HKLS, IPEAK1,IHKL1)
C
C Search AZI_ROT() values for groups of at least NGROUP that differ by < DAZI degrees.
C Find the best two peak/hkl pairs to complement IPEAK1/IHKL1 and add to triplet list
C if the dot-product exceeds DOT_CUTOFF in magnitude.
C If guides spots are present, then the group of pairs must contain a guide spot
C NREJ is a counter for the number of rejected triplets due to guide spots
C
	PARAMETER (NTRIP_MAX=100000)
	INTEGER IPEAKS_TRIP(3,NTRIP_MAX),IHKLS_TRIP(3,NTRIP_MAX)
C
	COMMON /AZI_MATCH_COM/ AZI_DANG_MAX,MIN_AZI_GROUP
C
	COMMON /TRIP_MATCH_COM/ DOT_CROSS_MIN
C
	COMMON /GUIDE_SPOTS_COM/ NGUIDES,IGUIDES(100)
C
	REAL H_HKLS(3,1000)
C
	PARAMETER (NPAIR_MAX=1000)
	INTEGER IPAIRS(2,NPAIR_MAX),ITAGS(NPAIR_MAX)
	REAL AZI_ROT(NPAIR_MAX)
C
C Make sure we don't think the last value is part of a larger group
C
	ITAGS(NPAIRS+1)=NPAIRS+1
	AZI_ROT(NPAIRS+1)=AZI_ROT(NPAIRS+1)+2.0*AZI_DANG_MAX
C
C Search for AZI_ROT groups by looping through possible ends of the group
C
	DO IEND=2,NPAIRS
        AZI_END=AZI_ROT(ITAGS(IEND))
C Search backwards in AZI_ROT values to find possible starts of the group
	  ISTART=IEND
	  DO ITRY=IEND-1,1,-1
	    AZI_TRY=AZI_ROT(ITAGS(ITRY))
	    IF(AZI_TRY .LT. AZI_END-AZI_DANG_MAX) EXIT
	    ISTART=ITRY
	  ENDDO
C Reject the group if it is too small
	  IF(IEND-ISTART .LT. MIN_AZI_GROUP) CYCLE
C Reject the group if it is part of a larger group
		IF(AZI_ROT(ITAGS(IEND+1)) .LT.
	1					AZI_ROT(ITAGS(ISTART))+AZI_DANG_MAX) CYCLE
C If using guide spots, reject group if it doesn't contain one
	  IF(NGUIDES .GT. 0) THEN
	    N=0
	    DO I=ISTART,IEND
	      IF(IPAIRS(1,ITAGS(I)) .LE. NGUIDES) N=N+1
	    ENDDO
	    IF(N .EQ. 0) THEN
	    	NREJ=NREJ+1
	      CYCLE
	    ENDIF
	  ENDIF
C Find the best three pairs in the group
	  CALL FIND_BEST_TRIPLET(IPAIR2,IPAIR3,DOT_CROSS,
	1						IPAIRS,ITAGS,ISTART,IEND, H_HKLS,IHKL1)
C If the triplet has a good dot-cross product, add to the triplet list
	  IF(DOT_CROSS .GE. DOT_CROSS_MIN) THEN
	    NTRIP=NTRIP+1
	    IPEAKS_TRIP(1,NTRIP)=IPEAK1
	    IPEAKS_TRIP(2,NTRIP)=IPAIRS(1,IPAIR2)
	    IPEAKS_TRIP(3,NTRIP)=IPAIRS(1,IPAIR3)
	    IHKLS_TRIP(1,NTRIP)=IHKL1
	    IHKLS_TRIP(2,NTRIP)=IPAIRS(2,IPAIR2)
	    IHKLS_TRIP(3,NTRIP)=IPAIRS(2,IPAIR3)
	  ENDIF
C
	ENDDO
C
	RETURN
	END



	SUBROUTINE FIND_BEST_TRIPLET(IPAIR2,IPAIR3,DOT_MAX,
	1						IPAIRS,ITAGS,ISTART,IEND, H_HKLS,IHKL1)
C
C Choose the two matches of hkls/peaks in the list IPAIRS(ISTART..IEND) to
C make a triplet with IHKL1/IPEAK1. The algorithm trys to maximise the angles
C between the triplet vectors, but it is not necessarily the best solution.
C The routine returns indices to IPAIRS(), IPAIR2 & IPAIR3, for the two
C other pairs, and returns the dot-product of the three vectors in DOT_MAX.
C
	PARAMETER (NPAIR_MAX=1000)
	INTEGER IPAIRS(2,NPAIR_MAX),ITAGS(NPAIR_MAX)
C
	REAL H_HKLS(3,1000)
C
	REAL VEC(3)
C
C Find the hkl most nearly perpendicular to HKL1
C
	ITAG2=0
	DOT_MIN=1E6
	DO ITAG=ISTART,IEND
	  I=ITAGS(ITAG)
	  IHKL2=IPAIRS(2,I)
	  DOT=ABS( VC3_DOT(H_HKLS(1,IHKL2),H_HKLS(1,IHKL1)) )
	  IF(DOT .LT. DOT_MIN) THEN
	    DOT_MIN=DOT
	    ITAG2=ITAG
	  ENDIF
	ENDDO
	IPAIR2=ITAGS(ITAG2)
C
C Find the hkl with the largest dot-cross product with HKL1,HKL2
C
	IHKL2 =IPAIRS(2,IPAIR2)
	CALL VC3_CROSS(VEC, H_HKLS(1,IHKL1), H_HKLS(1,IHKL2) )
	ITAG3=0
	DOT_MAX=-1E6
	DO ITAG=ISTART,IEND
	  IF(ITAG .NE. ITAG2) THEN
	    I=ITAGS(ITAG)
	    IHKL3=IPAIRS(2,I)
	    DOT=ABS( VC3_DOT(H_HKLS(1,IHKL3),VEC) )
	    IF(DOT .GT. DOT_MAX) THEN
	      DOT_MAX=DOT
	      ITAG3=ITAG
	    ENDIF
	  ENDIF
	ENDDO
	IPAIR3=ITAGS(ITAG3)
C
	RETURN
	END



	SUBROUTINE CALC_AZIMUTH_ANGLES(AZI_ROT,	IPAIRS,NPAIRS,
	1								H_PEAKS,H_HKLS, IPEAK1,IHKL1)
C
C Calculate the azimuthal angle (rotation about IPEAK1 or IHKL1) for all
C the second peaks in the matched pairs in IPAIRS, then subtract the angle
C for the peaks from the angle for the hkls and store in AZI_ROT().
C The zeroes for the azimuthal angles are arbitrary. Matched pairs with
C similar values of AZI_ROT are possibly from the same solution.
C
	PARAMETER (NPAIRS_MAX=1000)
	INTEGER IPAIRS(2,NPAIRS_MAX)
	REAL AZI_ROT(NPAIRS_MAX)
C
	REAL H_PEAKS(3,10000),H_HKLS(3,1000)
C
	REAL PERP1_PEAKS(3),PERP2_PEAKS(3),PERP1_HKLS(3),PERP2_HKLS(3)
C
C Create arbitrary vectors perpendicular to the scattering vectors
C for IPEAK1 & IHKL1
C
	CALL MAKE_PERP_VECS(PERP1_PEAKS,PERP2_PEAKS, H_PEAKS(1,IPEAK1))
	CALL MAKE_PERP_VECS(PERP1_HKLS,PERP2_HKLS, H_HKLS(1,IHKL1))
C
C Calculate azimuthal rotation of each matched peak and hkl relative
C to PEAK1 and HKL1. The difference between the two rotations is
C stored in AZI_ROT().
C
	DO I=1,NPAIRS
	  IPEAK2=IPAIRS(1,I)
	  IHKL2=IPAIRS(2,I)
C
	  DOT1=VC3_DOT(PERP1_PEAKS,H_PEAKS(1,IPEAK2))
	  DOT2=VC3_DOT(PERP2_PEAKS,H_PEAKS(1,IPEAK2))
	  ROT_PEAK=ATAN2D(DOT1,DOT2)
C
	  DOT1=VC3_DOT(PERP1_HKLS,H_HKLS(1,IHKL2))
	  DOT2=VC3_DOT(PERP2_HKLS,H_HKLS(1,IHKL2))
	  ROT_HKL=ATAN2D(DOT1,DOT2)
C
	  AZI_ROT(I)=ROT_PEAK-ROT_HKL
	  IF(AZI_ROT(I) .GT.  180.0) AZI_ROT(I)=AZI_ROT(I)-360.0
	  IF(AZI_ROT(I) .LT. -180.0) AZI_ROT(I)=AZI_ROT(I)+360.0
	ENDDO
C
	RETURN
	END



	SUBROUTINE MAKE_PERP_VECS(PERP1,PERP2, VEC_IN)
C
	REAL PERP1(3),PERP2(3),VEC_IN(3)
C
	REAL VEC1(3),VEC2(3),VEC3(3)
C
C Make a cross-product of VEC_IN with (100) & (010)
C and use the result with the larger norm
C
	CALL VC3_SET(VEC1, 1.0,0.0,0.0)
	CALL VC3_CROSS(VEC2, VEC_IN,VEC1)
	CALL VC3_SET(VEC1, 0.0,1.0,0.0)
	CALL VC3_CROSS(VEC3, VEC_IN,VEC1)
	IF(VC3_DOT(VEC2,VEC2) .GT. VC3_DOT(VEC3,VEC3)) THEN
	  CALL VC3_COPY(PERP1,VEC2)
	ELSE
	  CALL VC3_COPY(PERP1,VEC3)
	ENDIF
C
C Make second perpendicular vector
C
	CALL VC3_CROSS(PERP2, VEC_IN,PERP1)
C
C Make the perpendicular vectors unit vectors
C
	CALL VC3_UNIT(PERP1)
	CALL VC3_UNIT(PERP2)
C
	RETURN
	END



	SUBROUTINE FIND_PAIR_MATCHES(IPAIRS,NPAIRS,NPAIRS_MAX,
	1										IPEAK1,IHKL1, NPEAKS)
C
C Loop through IPEAK2, matching HKL angles to the angle between peaks
C IPEAK1 & IPEAK2, within +/- DANGLE, and store any matches in IPAIRS()
C
	INTEGER IPAIRS(2,NPAIRS_MAX)
C
	REAL H_PEAKS
	COMMON /PEAKS_COM/ H_PEAKS(3,10000)
C
	PARAMETER (NANG_MAX=150000)
	COMMON /ANGLES_COM/ NANGLES,ANGLES(NANG_MAX),IANGLES(2,NANG_MAX)
C
	COMMON /PAIR_MATCH_COM/ PAIR_ANG_MIN,PAIR_DANG_MAX
C
C Loop through IPEAK2
C
	NPAIRS=0
	DOT_CUTOFF=COSD(PAIR_ANG_MIN)
	DO IPEAK2=1,NPEAKS
C Ignore this pair if the angle between IPEAK1 & IPEAK2 is too large or small
	  DOT=VC3_DOT(H_PEAKS(1,IPEAK1),H_PEAKS(1,IPEAK2))
	  IF(ABS(DOT) .GT. DOT_CUTOFF) CYCLE
C Calculate the angle between IPEAK1 & IPEAK2
	  ANGLE=ACOSD(DOT)
C Search the angles list for the first and last indices with an angle
C within ANGLE +/- PAIR_DANG_MAX. Return zeroes on failure.
	  CALL FIND_ANGLES(IFIRST,ILAST, ANGLE,PAIR_DANG_MAX)
	  IF(IFIRST .NE. 0) THEN
C Loop through the list finding any first hkl corresponding to IHKL1
	    DO I=IFIRST,ILAST
	      IF(IANGLES(1,I) .EQ. IHKL1) THEN
C Save the second hkl & peak in IPAIRS()
	        NPAIRS=NPAIRS+1
	        IF(NPAIRS .GT. NPAIRS_MAX)
	1				CALL QUIT('BUG(find_pair_matches) NPAIRS > NPAIRS_MAX')
	        IPAIRS(1,NPAIRS)=IPEAK2
	        IPAIRS(2,NPAIRS)=IANGLES(2,I)
	      ENDIF
	    ENDDO
	  ENDIF
C
	ENDDO
C
	RETURN
	END
