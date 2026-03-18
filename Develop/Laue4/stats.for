C----- Routines to handle indentifying outliers -------------------------------
c	SUBROUTINE PROCESS_OUTLIERS
c	SUBROUTINE MARK_OUTLIERS(LBAD_SEQ, CUTOFF)
c	SUBROUTINE GET_OUTLIERS(OUTLIERS, CUTOFF)
C ---------------- Output statistics about reflection redundancies --------------
c	SUBROUTINE OUTPUT_HKL_REDUNDANCY
C-------------------- Output intensity statistics versus several factors --------------
c	SUBROUTINE OUTPUT_INTENSITY_STATS
c	SUBROUTINE OUTPUT_RATIO_BINS(SPAR)
c	SUBROUTINE OUTPUT_RATIO_BINS_FILES
c	SUBROUTINE OUTPUT_RATIO_BINS_PLATES
c	SUBROUTINE CALC_BINS_STATS(RAT,ESD,GOF,RBINS,ISTATUS, NBINS,LREAL)
C-------------------- Output intensity statistics versus d-spacing --------------
c	SUBROUTINE OUTPUT_DSPACE_STATS
c	SUBROUTINE CALC_DSPACE_RFACTOR(R1,R2,GOOF,NUNIQ,
c	1				NSIG3,NSIG10,NSIG30, D_LO,D_HI, RSIG,DSPACES)
C------------------------------------------------------------------------------


C----- Routines to handle indentifying outliers -------------------------------

	SUBROUTINE PROCESS_OUTLIERS
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	LOGICAL LOPEN,LBAD_SEQ(NSEQ_MAX)
	CHARACTER FILE_NAME*255
	REAL OUTLIERS(NDATA_MAX)
C
C
	PRINT *
	WRITE(10,'(/,A,/)') '====== Analysis of intensity outliers ======'
C
C Set flag that we have not yet opened the "laue4_bad.lis" file
C
	LOPEN=.FALSE.
C
C Mark all sequences as not containing any outliers
C
	DO ISEQ=1,NSEQ
	  LBAD_SEQ(ISEQ)=.FALSE.
	ENDDO
C
C Print out information on potential outliers
C
	CALL GET_OUTLIERS(OUTLIERS,0.0)		! cutoff = 0.0, so all data
	OUT_MIN=+1E6
	OUT_MAX=-1E6
	DO I=1,ILAST(NSEQ)
	  OUT_MIN=MIN(OUT_MIN,OUTLIERS(I))
	  OUT_MAX=MAX(OUT_MAX,OUTLIERS(I))
	ENDDO
C
C Start loop with user input on outlier cutoff (or exit loop)
C
	PRINT '(1X,A,2(SP,F6.1,A))','Intensity outliers found ',OUT_MIN,
	1					' to',OUT_MAX,' esds from merged values'	
C
100	PRINT '(/,1X,A,$)','Cutoff (in esds) for rejecting outliers [RETURN = skip]: '
	READ(*,'(F10.0)',ERR=100) CUTOFF
C
C If a valid cutoff, identify and mark outliers, then loop back for more
C
	IF(CUTOFF .GT. 0.0) THEN
C If laue4_bad.lis is not open, do so on UNIT=3
	  IF( .NOT.LOPEN ) THEN
	    LOPEN=.TRUE.
	    OPEN(UNIT=3,STATUS='UNKNOWN',FILE='laue4_bad.lis')
	    IF(IDATA_MODE .EQ. 1) THEN
	      WRITE(3,'(A,//,A)') 'LaueG Rejects File Version 2, Option 1',
	1	  'Base-name                                  h   k   l   Comments ...'
	    ELSE
	      WRITE(3,'(A,//,A)') 'LaueG Rejects File Version 2, Option 2',
	1	  'Base-name                                  h   k   l   m   Comments ...'
	    ENDIF
        ENDIF
C Find new outliers, write them to UNIT=3 and mark ISIG() and LBAD_SEQ()
	  CALL MARK_OUTLIERS(LBAD_SEQ, CUTOFF)
C Loop back to get user input
	  GOTO 100
	ENDIF
C
C If no outliers written to laue4_bad.lis, log this and then exit
C
	PRINT *
	IF( .NOT.LOPEN ) THEN
	  WRITE(10,'(A)') 'No outliers identified'
        RETURN
      ENDIF
C
C Close "laue4_bad.lis", and tell user about the file
C
	CLOSE(UNIT=3)
	PRINT '(1X,A)','Individual outliers written to "laue4_bad.lis"'
C
C Open "laue4_sus.lis" and write the header
C
	OPEN(UNIT=3,STATUS='UNKNOWN',FILE='laue4_sus.lis')
	IF(IDATA_MODE .EQ. 1) THEN
	  WRITE(3,'(A,//,A)') 'LaueG Rejects File Version 2, Option 1',
	1	  'Base-name                                  h   k   l   Comments ...'
	ELSE
	  WRITE(3,'(A,//,A)') 'LaueG Rejects File Version 2, Option 2',
	1	  'Base-name                                  h   k   l   m   Comments ...'
	ENDIF
C
C Calculate all "outliers" with a cutoff of zero
C
	CALL GET_OUTLIERS(OUTLIERS,0.0)
C
C Output all sequences with an outlier to "laue4_sus.lis"
C
	NOUT_SEQ=0
	NOUT=0
	DO ISEQ=1,NSEQ
	  IF( LBAD_SEQ(ISEQ) ) THEN
	    NOUT_SEQ=NOUT_SEQ+1
	    DO I=IFIRST(ISEQ),ILAST(ISEQ)
	      NOUT=NOUT+1
	      CALL MAKE_FILE_NAME(FILE_NAME, IFILE_NUM(IFILES(I)),'')
	      RSIGMA=MAX(-99.9,MIN(99.9, OUTLIERS(I) ))
	      IX=NINT(XY_PIX(1,I))
	      IY=NINT(XY_PIX(2,I))
	      LEN_NAME=MAX(40,LEN_TRIM(FILE_NAME))
	      IF(IDATA_MODE .EQ. 1) THEN
	        WRITE(3,'(A,3I4,A,F5.2,A,2I5,A,F5.1,A)') FILE_NAME(1:LEN_NAME),
	1		  (HKLMS(K,I),K=1,3),'  ! wav=',WAVS(I),
	2		  '  x,y=(',IX,IY,')  I(hkl) - Imerge =',RSIGMA,' esds'
	      ELSE
	        WRITE(3,'(A,4I4,A,F5.2,A,2I5,A,F5.1,A)') FILE_NAME(1:LEN_NAME),
	1		  (HKLMS(K,I),K=1,4),'  ! wav=',WAVS(I),
	2		  '  x,y=(',IX,IY,')  I(hkl) - Imerge =',RSIGMA,' esds'
		  ENDIF
	    ENDDO
	  ENDIF
	ENDDO
	CLOSE(UNIT=3)
C
C If no outliers actually written to files, nothing to do, simply return
C
	IF(NOUT_SEQ .EQ. 0) RETURN
C
C Explain what we have written to output files
C
	PRINT '(1X,A,I5,A)','Total of',NOUT_SEQ,
	1	' suspect groups of equivalents with outliers'
	PRINT '(1X,A,I6,A,/)','All',NOUT,
	1	' reflections in these groups written to "laue4_sus.lis"'
C
C Calculate new merge statistics excluding any outliers
C
	PRINT '(1X,A)','Merge statistics with outliers excluded:'
	CALL PRINT_LSQ_MERGE
C
C Ask if outliers are to be removed from data
C
	PRINT *
	IVAL=1	
	CALL ASK_YES_NO(IVAL,'Exclude outliers from output data?')
C Count number of outliers, and reinstate outliers if answer is NO
	NREJ=0
	DO I=1,ILAST(NSEQ)
	  IF(ISIG(I) .LT. 0.0) THEN
	    NREJ=NREJ+1
C If NOT excluding outliers, negate ISIG() to remove outlier status
	    IF(IVAL .EQ. 0) ISIG(I)=-ISIG(I)
	  ENDIF
	ENDDO
C Final output about outliers
	IF(IVAL .EQ. 0) THEN
	  PRINT '(1X,I6,A)',NREJ,' outliers will be included in final output'
	ELSE
	  PRINT '(1X,I6,A)',NREJ,' outliers have been removed'
	ENDIF
C
C Recalculate COUNTS_SEQ() (ignores any spots with ISIG() < 0)
C
	CALL CALC_ALL_COUNTS_SEQ
C
C Output what has happened to the log file
C
	WRITE(10,'(/,I6,A)') NREJ,
	1		' outliers identified and written to "laue4_bad.lis"'
	WRITE(10,'(A,I5,A)') 'Outliers are part of',NOUT_SEQ,
	1						' suspect groups of equivalents'
	WRITE(10,'(A,I6,A)') 'All',NOUT,
	1		' reflections in these groups written to "laue4_sus.lis"'
C
	IF(IVAL .EQ. 0) THEN
	  WRITE(10,'(/,A)') 'Outliers have NOT been removed from final intensities'
	ELSE
	  WRITE(10,'(/,A)') 'Outliers have been removed from final intensities'
	  WRITE(10,'(/,A)') 'Merge statistics with outliers removed:'
	  CALL LOG_LSQ_MERGE
	ENDIF
C
	RETURN
	END



	SUBROUTINE MARK_OUTLIERS(LBAD_SEQ, CUTOFF)
C
C Find new group outliers with DIFF/ESD > CUTOFF, mark the group in LBAD_SEQ(),
C negate ISIG(), and add worst outliers to laue4_bad.lis.
C
	LOGICAL LBAD_SEQ(*)
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	CHARACTER FILE_NAME*255
	REAL OUTLIERS(NDATA_MAX)
C
C Find new outliers above CUTOFF and return in OUTLIERS() the
C DIFF/ESD of the worst reflection. Otherwise, return OUTLIERS() = 0
C
      CALL GET_OUTLIERS(OUTLIERS, CUTOFF)
C
C Zero counters for outlier statistics
C
	NOUT=0
	OUT_SUM=0.0
	OUT_SUM2=0.0
	OUT_MIN=+1E6
	OUT_MAX=-1E6
C
C Loop through all reflections finding non-zero OUTLIERS()
C
	DO ISEQ=1,NSEQ
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    IF(OUTLIERS(I) .EQ. 0.0) CYCLE
C Mark this sequence as containing an outlier
	    LBAD_SEQ(ISEQ)=.TRUE.
C Update statistics about outliers
		NOUT=NOUT+1
	    OUT_SUM=OUT_SUM+OUTLIERS(I)
	    OUT_SUM2=OUT_SUM2+OUTLIERS(I)**2
	    OUT_MIN=MIN(OUT_MIN,OUTLIERS(I))
	    OUT_MAX=MAX(OUT_MAX,OUTLIERS(I))
C Write the outlier to "laue4_bad.lis" on UNIT=3
	    CALL MAKE_FILE_NAME(FILE_NAME, IFILE_NUM(IFILES(I)),'')
	    RSIGMA=MAX(-99.9,MIN(99.9, OUTLIERS(I) ))
	    IX=NINT(XY_PIX(1,I))
	    IY=NINT(XY_PIX(2,I))
	    LEN_NAME=MAX(40,LEN_TRIM(FILE_NAME))
	    IF(IDATA_MODE .EQ. 1) THEN
	      WRITE(3,'(A,3I4,A,F5.2,A,2I5,A,F5.1,A)') FILE_NAME(1:LEN_NAME),
	1		  (HKLMS(K,I),K=1,3),'  ! wav=',WAVS(I),
	2		  '  x,y=(',IX,IY,')  I(hkl) - Imerge =',RSIGMA,' esds'
	    ELSE
	      WRITE(3,'(A,4I4,A,F5.2,A,2I5,A,F5.1,A)') FILE_NAME(1:LEN_NAME),
	1		  (HKLMS(K,I),K=1,4),'  ! wav=',WAVS(I),
	2		  '  x,y=(',IX,IY,')  I(hkl) - Imerge =',RSIGMA,' esds'
	    ENDIF
C
	  ENDDO
	ENDDO
C
C Output statistics on the outliers
C
	AVE=OUT_SUM/MAX(1,NOUT)
	STD=SQRT(MAX(0.0, OUT_SUM2/MAX(1,NOUT) - AVE**2) )
	PRINT '(I6,A)',NOUT,' groups with a worst outlier above cutoff'
	IF(NOUT .GT. 0) PRINT '(3X,4(A,F6.1,2X),A)','Ave=',AVE,'Spread=',STD,
	1				'Min=',OUT_MIN,'Max=',OUT_MAX,'(esds)'

C
      WRITE(10,'(A,I4,A,F4.1,A)') 'Found',NOUT,' groups with an outlier >',MIN(99.9,CUTOFF),' esds,'
C
C Negate ISIG() if OUTLIERS() has a non-zero value
C
	DO I=1,ILAST(NSEQ)
	  IF(OUTLIERS(I) .NE. 0.0) ISIG(I)=-ABS(ISIG(I))
	ENDDO
C
C Recalculate COUNTS_SEQ() (ignores any spots with ISIG() < 0)
C
	CALL CALC_ALL_COUNTS_SEQ
C
	RETURN
	END


	SUBROUTINE GET_OUTLIERS(OUTLIERS, CUTOFF)
C
C Return in OUTLIERS() the relative errors (Obs-Calc)/Sig for the worst
C "new" outlier in a sequence if the absolute relative error exceeds CUTOFF.
C Otherwise, OUTLIERS() is zero. In the special case where a sequence has
C only two reflections not already outliers, if one exceeds CUTOFF then
C both return their relative errors in OUTLIERS().
C
	REAL OUTLIERS(*)
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
C Loop through sequences:
C
	DO ISEQ=1,NSEQ
C
C Load OUTLIERS() with (Obs-Calc)/Sig, or 0 if already an outlier
C Also, count how many non-outliers there are
C
	  NOK=0
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    IF(ISIG(I) .LT. 0.0) THEN		! Already an outlier ?
	      OUTLIERS(I)=0.0
	    ELSE							! Not an outlier ?
	      OUTLIERS(I)=(IOBS(I)-ICALC(I))/ISIG(I)
	      NOK=NOK+1
	    ENDIF
	  ENDDO
C
C Identify the worst outlier
C
	  IWORST=IFIRST(ISEQ)
	  DO I=IFIRST(ISEQ)+1,ILAST(ISEQ)
	    IF(ABS(OUTLIERS(I)) .GT. ABS(OUTLIERS(IWORST))) IWORST=I
	  ENDDO
C
C Zero the OUTLIERS(), except the worst one
C
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    IF(I .NE. IWORST) OUTLIERS(I)=0.0
	  ENDDO
C
C If the worst outlier is below the cutoff, zero its OUTLIERS() value.
C If above the cutoff, and only two non-outliers, both have
C OUTLIERS() values equal but opposite.
C
	  IF(ABS(OUTLIERS(IWORST)) .LT. CUTOFF) THEN
	    OUTLIERS(IWORST)=0.0
	  ELSE IF(NOK .EQ. 2) THEN
	    DO I=IFIRST(ISEQ),ILAST(ISEQ)
	      IF(ISIG(I).GE.0.0 .AND. I.NE.IWORST) ISECOND=I
	    ENDDO
C Make the other non-outlier have an opposite OUTLIERS() value
	    OUTLIERS(ISECOND)=-OUTLIERS(IWORST)
	  ENDIF
C
C Loop back for the next sequence
C
	ENDDO
C
	RETURN
	END



C ---------------- Output statistics about reflection redundancies --------------

	SUBROUTINE OUTPUT_HKL_REDUNDANCY
C
      COMMON /MODES_COM/ IFILE_MODE,IDATA_MODE
C
      PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	PARAMETER (NDATA_MAX=2000000)
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
C
	INTEGER NBIN1(10),NBIN2(10)
C
	REAL RBIN(10)
	DATA RBIN/0.0,1.4,2.4,3.5,4.5,5.5,6.5,7.5,10.5,19.5/
C
C Zero the bins
C
	DO I=1,10
	  NBIN1(I)=0
	  NBIN2(I)=0
	ENDDO
C
C Loop through all sequences adding up data redundancies
C
	NTOTAL=0
	SUM_REFF=0.0
	DO ISEQ=1,NSEQ
C Sum ISIG**2 & ISIG**4 to calculate the effective redundancy
	  NSUM=0
	  SUM_SIG2=0.0
	  SUM_SIG4=0.0
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
C Ignore a peak if marked as an outlier by ISIG < 0
	    IF(ISIG(I) .GT. 0.0) THEN
	      NSUM=NSUM+1
	      SUM_SIG2=SUM_SIG2+ISIG(I)**2
	      SUM_SIG4=SUM_SIG4+ISIG(I)**4
	    ENDIF
	  ENDDO
C
	  IF(NSUM .GT. 0) THEN
C Calculate the effective redundancy for a weighted merge
	    REFF=SUM_SIG2**2 / SUM_SIG4
C Find the "numeric" and "effective" bin numbers to use
	    DO IBIN=1,10
	      IF(FLOAT(NSUM) .GT. RBIN(IBIN)) IBIN1=IBIN
	      IF(REFF        .GT. RBIN(IBIN)) IBIN2=IBIN
	    ENDDO
C Increment the "numeric" and "effective" bins

	    NBIN1(IBIN1)=NBIN1(IBIN1)+1
	    NBIN2(IBIN2)=NBIN2(IBIN2)+1
C Increment the sums used for the mean redundancy
	    NTOTAL=NTOTAL+NSUM
	    SUM_REFF=SUM_REFF+REFF
C
	  ENDIF
C
	ENDDO
C
C Output the average redundancy values
C
	RED=FLOAT(NTOTAL)/MAX(1,NSEQ)
	REFF=SUM_REFF/MAX(1,NSEQ)
	PRINT '(/,1X,A,2(F6.1,A),/)','Data redundancy:',RED,
	1							' (unweighted)',REFF,' (weighted)'
C
C Convert # in bins to % of overall number
C
	DO I=1,10
	  NBIN1(I)=NINT( 100.0*NBIN1(I)/NSEQ )
	  NBIN2(I)=NINT( 100.0*NBIN2(I)/NSEQ )
	ENDDO
C
C Write the binned redundancy statistics to the log file
C
	WRITE(10,'(/,A,//,A,2(/,A,7I4,3I5,A),/,2X,A,2(F6.1,A))')
	1	'====== Redundancy of data (equivalents and repeats) ======',
	2	'           1   2   3   4   5   6   7  8-10 11-20 >20',
	3	'Fraction',NBIN1,' %  (unweighted)',
	4	'Fraction',NBIN2,' %   (weighted)',
	5	'mean values',RED,' (unweighted)',REFF,' (weighted)'
	IF(IDATA_MODE .EQ. 3) WRITE(10,'(/,A)')
     1        'WARNING: Redundancies are incorrect for modulated structures'
C
	RETURN
	END



C-------------------- Output intensity statistics versus several factors --------------

	SUBROUTINE OUTPUT_INTENSITY_STATS
C
C Output to the log file data statistics versus various experimental parameters
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	REAL DCOUNTS_SEQ(NSEQ_MAX)
C
C CALC_BINS_STATS & SORT_BIN_DATA_FUNC use data from:
	COMMON /SORT_BIN_COM/ BIN_DATA(NDATA_MAX)
C
	WRITE(10,'(/,A)') '====== Comparison of observed and equivalent'//
	1														' intensities ======'
C
C
	DO I=1,NDATA
	  BIN_DATA(I)=WAVS(I)
	ENDDO
	CALL OUTPUT_RATIO_BINS('wavelength')
C
C
	DO I=1,NDATA
	  BIN_DATA(I)=TTHS(I)
	ENDDO
	CALL OUTPUT_RATIO_BINS('two-theta')
C
C
	DO I=1,NDATA
	  BIN_DATA(I)=XY_PIX(1,I)
	ENDDO
	CALL OUTPUT_RATIO_BINS('X-pixel')
C
C
	DO I=1,NDATA
	  BIN_DATA(I)=XY_PIX(2,I)
	ENDDO
	CALL OUTPUT_RATIO_BINS('Y-pixel')
C
C
	DO I=1,NDATA
	  BIN_DATA(I)=SIND(TTHS(I)/2.0)/WAVS(I)
	ENDDO
	CALL OUTPUT_RATIO_BINS('(sin(theta)/wav)')
C
C
	DO ISEQ=1,NSEQ
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    BIN_DATA(I)=COUNTS_SEQ(ISEQ)
	  ENDDO
	ENDDO
	CALL OUTPUT_RATIO_BINS('corrected-intensity')
C
C
	CALL CALC_ALL_DCOUNTS_SEQ(DCOUNTS_SEQ)
	DO ISEQ=1,NSEQ
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    BIN_DATA(I)=DCOUNTS_SEQ(ISEQ)
	  ENDDO
	ENDDO
	CALL OUTPUT_RATIO_BINS('corrected-intensity esds')
C
C
	DO I=1,NDATA
	  BIN_DATA(I)=ICALC(I)/ISIG(I)
	ENDDO
	CALL OUTPUT_RATIO_BINS('(intensity/esd)')
C
C
	DO ISEQ=1,NSEQ
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    EX_FAC=CALC_EXTI_FACTOR(COUNTS_SEQ(ISEQ),WAVS(I),TTHS(I))
	    BIN_DATA(I)=1.0 - 1.0/SQRT(MAX(1E-4,1.0+EX_FAC))
	  ENDDO
	ENDDO
	CALL OUTPUT_RATIO_BINS('extinction correction')
C
C
	DO I=1,NDATA
	    BIN_DATA(I)=1.0 - CALC_ABS_FACTOR(WAVS(I),TTHS(I))
	ENDDO
	CALL OUTPUT_RATIO_BINS('absorption correction')
C
C Special case of stats per image file
	CALL OUTPUT_RATIO_BINS_FILES
C
C Explicitly say what the statistics are and the bin size
C
	WRITE(10,'(/,4(2X,A,/))')
	1	'Ratio = Sum{ Iobs } / Sum{ Iequiv }',
	1	'(esd) = sqrt[ Sum{ esd(Iobs)^2 } ] / Sum{ Iequiv }',
	2	'GOOF  = sqrt[ Sum{ (Iobs-Iequiv)^2 / esd(Iobs)^2 } / N ]',
	3	'       where N is the degrees of freedom of the merge'
	WRITE(10,'(2X,A)')
	1	'Data is binned in groups of approximately equal size'
C
	RETURN
	END



	SUBROUTINE OUTPUT_RATIO_BINS(SPAR)
C
	CHARACTER SPAR*(*)
C
C CALC_BINS_STATS & SORT_BIN_DATA_FUNC use data from:
C	COMMON /SORT_BIN_COM/ BIN_DATA(NDATA_MAX)
C
	CHARACTER LINE_OUT*132
	REAL RAT(12),ESD(12),GOF(12),RBINS(12+1)
C
C Calculate stats for 12 bins of the array in /SORT_BIN_COM/
C Just return if failure signalled by ISTATUS=0
C
	NBINS=12
	CALL CALC_BINS_STATS(RAT,ESD,GOF,RBINS,ISTATUS, NBINS,.TRUE.)
	IF(ISTATUS .EQ. 0) RETURN
C
C Output relevant header line to log file
C
	ISCALE=1
	IF(RBINS(NBINS+1) .GT. 9999.0) ISCALE=10
	IF(RBINS(NBINS+1) .GT. 99999.0) ISCALE=100
	IF(ISCALE .EQ. 1) THEN
	  WRITE(10,'(/,2A)') 'Binned versus ',SPAR
	ELSE
	  WRITE(10,'(/,3A,I4)') 'Binned versus ',SPAR,' /',ISCALE
	ENDIF
C
C Output to log file the limits of the bins in terms of BIN_DATA
C Change sig-figs of output depending on maximum BIN_DATA value
C
	IF     (RBINS(NBINS) .LE. 9.9) THEN
	  WRITE(LINE_OUT,'(3X,21F5.2)',IOSTAT=IDUM) (RBINS(K)/ISCALE,K=1,NBINS)
	ELSE IF(RBINS(NBINS) .LE. 99.9) THEN
	  WRITE(LINE_OUT,'(3X,21F5.1)',IOSTAT=IDUM) (RBINS(K)/ISCALE,K=1,NBINS)
	ELSE
	  WRITE(LINE_OUT,'(3X,21I5)',IOSTAT=IDUM) (NINT(RBINS(K)/ISCALE),K=1,NBINS)
	ENDIF
C The final value is often much larger, so output it individually
	IF     (RBINS(NBINS+1) .LE. 9.9) THEN
	  WRITE(10,'(A,F5.2)',IOSTAT=IDUM) LINE_OUT(1:3+5*NBINS),RBINS(NBINS+1)
	ELSE IF(RBINS(NBINS+1) .LE. 99.9) THEN
	  WRITE(10,'(A,F5.1)',IOSTAT=IDUM) LINE_OUT(1:3+5*NBINS),RBINS(NBINS+1)
	ELSE
	  WRITE(10,'(A,I5)',IOSTAT=IDUM) LINE_OUT(1:3+5*NBINS),
	1					NINT(RBINS(NBINS+1)/ISCALE)
	ENDIF
C
C Output the relative average intensities and esds
C
	WRITE(10,'(A,20F5.2)') 'Ratio ',(MAX(-9.99,MIN(9.99,RAT(K))),K=1,NBINS)
	WRITE(10,'(A,20F5.2)') ' (esd)',(MAX(-9.99,MIN(9.99,ESD(K))),K=1,NBINS)
	WRITE(10,'(A,20F5.2)') 'GOOF  ',(MAX(-9.99,MIN(9.99,GOF(K))),K=1,NBINS)
C
C
	RETURN
	END



	SUBROUTINE OUTPUT_RATIO_BINS_FILES
C
C Special case of OUTPUT_RATIO_BINS for image file numbers
C
	EXTERNAL SORT_BIN_DATA_FUNC
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
C CALC_BINS_STATS & SORT_BIN_DATA_FUNC use data from:
	COMMON /SORT_BIN_COM/ BIN_DATA(NDATA_MAX)
C
	REAL RAT(NFILES_MAX),ESD(NFILES_MAX),GOF(NFILES_MAX),RBINS(NFILES_MAX+1)
C
C Create array to bin versus IFILES values
C
	DO I=1,ILAST(NSEQ)
	  BIN_DATA(I)=IFILES(I)
	ENDDO
C
C Create binned statistics
C
	CALL CALC_BINS_STATS(RAT,ESD,GOF,RBINS,ISTATUS, NFILES,.FALSE.)
C
C Output relevant header line to log file
C
	WRITE(10,'(/,2A)') 'Binned versus File Number'
C
	DO ILO=1,NFILES,12
	  IHI=MIN(ILO+11,NFILES)
	  IF(ILO .NE. 1) WRITE(10,*)
	  WRITE(10,'(5X,21I5)') (IFILE_NUM(K),K=ILO,IHI) 
	  WRITE(10,'(A,20F5.2)') 'Ratio ',(MAX(-9.99,MIN(9.99,RAT(K))),K=ILO,IHI)
	  WRITE(10,'(A,20F5.2)') ' (esd)',(MAX(-9.99,MIN(9.99,ESD(K))),K=ILO,IHI)
	  WRITE(10,'(A,20F5.2)') 'GOOF  ',(MAX(-9.99,MIN(9.99,GOF(K))),K=ILO,IHI)
	ENDDO
C
	RETURN
	END



	SUBROUTINE OUTPUT_RATIO_BINS_PLATES
C
C Special case of OUTPUT_RATIO_BINS for CYCLOPS plate/CCD
C
	EXTERNAL SORT_BIN_DATA_FUNC
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	PARAMETER (NFILES_MAX=400)
	COMMON /FILE_NUM_COM/ IFILE_NUM(NFILES_MAX),NFILES
C
C CALC_BINS_STATS & SORT_BIN_DATA_FUNC use data from:
	COMMON /SORT_BIN_COM/ BIN_DATA(NDATA_MAX)
C
	REAL RAT(NFILES_MAX),ESD(NFILES_MAX),GOF(NFILES_MAX),RBINS(NFILES_MAX+1)
C
C Create array to bin versus PLATE (1-8 is lower plate, 9-16 is upper plate)
C
	DO I=1,ILAST(NSEQ)
	  BIN_DATA(I)=1+IFIX(XY_PIX(1,I)/480.0)
	  IF(XY_PIX(2,I) .GT. 600.0) BIN_DATA(I)=BIN_DATA(I)+8
	ENDDO
C
C Create binned statistics
C
	NBINS=16
	CALL CALC_BINS_STATS(RAT,ESD,GOF,RBINS,ISTATUS, NBINS,.FALSE.)
C
C Output relevant header line to log file
C
	DO I=0,8,8
	  IF(I .EQ. 0) THEN
	    WRITE(10,'(/,2A)') 'Binned versus Plate Number (lower half)'
	  ELSE
	    WRITE(10,'(/,2A)') 'Binned versus Plate Number (upper half)'
	  ENDIF
	  WRITE(10,'(5X,21I5)') (K,K=1,8) 
	  WRITE(10,'(A,20F5.2)') 'Ratio ',(MAX(-9.99,MIN(9.99,RAT(K+I))),K=1,8)
	  WRITE(10,'(A,20F5.2)') ' (esd)',(MAX(-9.99,MIN(9.99,ESD(K+I))),K=1,8)
	  WRITE(10,'(A,20F5.2)') 'GOOF  ',(MAX(-9.99,MIN(9.99,GOF(K+I))),K=1,8)
	ENDDO
C
	RETURN
	END


	SUBROUTINE CALC_BINS_STATS(RAT,ESD,GOF,RBINS,ISTATUS, NBINS,LREAL)
C
C Calc the RATIO & esd, GOF and the bin limits RBINS(i) - RBINS(i+1) for
C NBINS bins of increasing value, the values are BIN_DATA() in /SORT_BIN_COM/.
C If LREAL is true,  BIN_DATA() are real values and bins are made of
C                    approx. equal number of members
C If LREAL is false, BIN_DATA() are integer values from 1 to NBINS
C		       which are the bins number to use
C
	LOGICAL LREAL
	REAL RAT(NBINS),ESD(NBINS),GOF(NBINS),RBINS(NBINS+1)
C
	EXTERNAL SORT_BIN_DATA_FUNC
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C      
      REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	COMMON /SORT_BIN_COM/ BIN_DATA(NDATA_MAX)
C
	INTEGER ITAGS(NDATA_MAX),ISEQS(NDATA_MAX)
C
C Do tagged ascending sort of BIN_DATA 
C
	CALL TAG_SORT(SORT_BIN_DATA_FUNC,ITAGS,NDATA)
C
C If no variation in BIN_DATA, return with ISTATUS=0 to signal failure
C for the case when absorption/extinction corections not used
C
	ISTATUS=1
	IF(ABS(BIN_DATA(ITAGS(1))-BIN_DATA(ITAGS(NDATA))) .LT. 1E-4) THEN
	  ISTATUS=0
	  RETURN
	ENDIF
C
C Create a list of the sequence numbers for each refln
C
	ILIST=0
	DO ISEQ=1,NSEQ
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    ILIST=ILIST+1
	    ISEQS(ILIST)=ISEQ
	  ENDDO
	ENDDO
C
C Prune ITAGS() to remove any reflns in a sequence of 1 (no equivs)
C
	NTAGS=0
	DO I=1,NDATA
	  ISEQ=ISEQS(ITAGS(I))
	  IF(IFIRST(ISEQ) .NE. ILAST(ISEQ)) THEN
	    NTAGS=NTAGS+1
	    ITAGS(NTAGS)=ITAGS(I)
	  ENDIF
	ENDDO
C
C Average ICALC, IOBS, ISIG in NBINS bins of ascending BIN_DATA
C For real data the bins are approx. equally spaced in terms on BIN_DATA,
C otherwise, (ie IFILE,IPLATE) the bins correspond to the BIN_DATA integer values
C
	IHI=0
	DO IBIN=1,NBINS
C Determine limits of data in the bin in terms of an index in ITAGS()
	  IF( LREAL ) THEN
	    ILO=NINT( FLOAT(IBIN-1)*NTAGS/NBINS ) +1
	    IHI=NINT( FLOAT(IBIN  )*NTAGS/NBINS )
	  ELSE
	    ILO=IHI+1
	    DO I=ILO,NTAGS
	      IF(NINT(BIN_DATA(ITAGS(I))) .GT. IBIN) EXIT
	      IHI=I
	    ENDDO
	  ENDIF
C Form the sums for averages (DOF is the degrees of freedom)
	  NSUM=0
	  SUMC=0.0
	  SUMO=0.0
	  SUMV=0.0
	  SUMG=0.0
	  DOF=0.0
	  DO I=ILO,IHI
C ITAG contains the actual index of the refln in ICALC(), etc.
	    ITAG=ITAGS(I)
	    NSUM=NSUM+1
	    SUMC=SUMC+ICALC(ITAG)
	    SUMO=SUMO+IOBS(ITAG)
	    SUMV=SUMV+ISIG(ITAG)**2
	    SUMG=SUMG + ( (IOBS(ITAG)-ICALC(ITAG))/ISIG(ITAG) )**2
C Calculating the degrees of freedom, DOF, is a bit tricky as we
C are not summing over complete sequences.
	    ISEQ=ISEQS(ITAG)
	    DOF=DOF +1.0 -1.0/(ILAST(ISEQ)-IFIRST(ISEQ)+1)
	  ENDDO
C Fudges to prevent divide by zero errors
	  NSUM=MAX(1,NSUM)
	  DOF=MAX(0.5,DOF)
	  SUMS=SQRT(SUMV)
	  IF(ABS(SUMC) .LT. 0.01*ABS(SUMS)) SUMC=0.01*SUMS
C Calculate the averages
	  RAT(IBIN)=ABS(SUMO)/MAX(1E-6,SUMC)
	  ESD(IBIN)=SUMS/MAX(1E-6,SUMC)
	  GOF(IBIN)=SQRT(SUMG/MAX(1E-6,DOF))
C Store the start value of BIN_DATA for this bin
	  IF( LREAL ) THEN
	    RBINS(IBIN)=BIN_DATA(ITAGS(ILO))
	  ELSE
	    RBINS(IBIN)=IBIN
	  ENDIF
	ENDDO
C Store the end value of BIN_DATA for the last bin
	IF( LREAL ) RBINS(NBINS+1)=BIN_DATA(ITAGS(NTAGS))
C
	RETURN
	END


C-------------------- Output intensity statistics versus d-spacing --------------

	SUBROUTINE OUTPUT_DSPACE_STATS
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	PARAMETER (NSEQ_MAX=100000)
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
	REAL DSPACES(NSEQ_MAX),RSIG(NSEQ_MAX)
C
	PARAMETER(IBIN_MAX=11)
	INTEGER NUNIQ(IBIN_MAX),NSIG3(IBIN_MAX)
	INTEGER NSIG10(IBIN_MAX),NSIG30(IBIN_MAX)
	REAL D_BIN(IBIN_MAX+1),R1(IBIN_MAX),R2(IBIN_MAX),GOOF(IBIN_MAX)
C
C Calculate the obs/esd for COUNTS_SEQ() and put in RSIG()
C
	CALL CALC_ALL_DCOUNTS_SEQ(RSIG)
	DO I=1,NSEQ
	  RSIG(I)=COUNTS_SEQ(I)/RSIG(I)
	ENDDO
C
C Calculate d-spacings for all sequences and the min & max values
C
	D_MIN=1E6
	D_MAX=0.0
	DO ISEQ=1,NSEQ
	  DSPACES(ISEQ)=WAVS(IFIRST(ISEQ)) / ( 2.0*SIND( 0.5*TTHS(IFIRST(ISEQ)) ) )
	  D_MIN=MIN(D_MIN,DSPACES(ISEQ))
	  D_MAX=MAX(D_MAX,DSPACES(ISEQ))
	ENDDO
C
C Calculate the d-min to use for each bin
C
	D_HI=D_MIN-1E-3
	DO I=1,IBIN_MAX
	  D_LO=D_HI
	  D_HI=D_MIN * 1.03**(I + 0.2*I**2)
	  D_HI=0.01*NINT(D_HI/0.01)
	  IF(I .GE. IBIN_MAX) D_HI=D_MAX+1E-3
	  D_HI=MIN(D_MAX+1E-3,D_HI)
	  D_BIN(I)=D_LO
	  IF(D_HI .GE. D_MAX-1E-3 ) EXIT
	ENDDO
	NBINS=MAX(IBIN_MAX,I)
	D_BIN(NBINS+1)=D_MAX+1E-3
C
C Calculate statistics for the various bins
C
	D_HI=D_MIN
	DO I=1,NBINS
	  CALL CALC_DSPACE_RFACTOR(R1(I),R2(I),GOOF(I),NUNIQ(I),
	1		NSIG3(I),NSIG10(I),NSIG30(I), D_BIN(I),D_BIN(I+1),
	2											RSIG,DSPACES)
	ENDDO
C
C Output results to the user
C
	PRINT '(1X,A,15F6.2)','d-space'  ,(D_BIN(K),K=NBINS+1,1,-1)
	PRINT '(1X,A,4X,14F6.1)','R1(%) ',(R1(K)   ,K=NBINS,1,-1)
	PRINT '(1X,A,4X,14F6.1)','R2(%) ',(R2(K)   ,K=NBINS,1,-1)
	PRINT '(1X,A,4X,14F6.1)','GOOF  ',(GOOF(K) ,K=NBINS,1,-1)
	PRINT '(1X,A,4X,14I6)'  ,'Nuniq ',(NUNIQ(K),K=NBINS,1,-1)
	PRINT '(1X,A,4X,14I6)'  ,'>3sig ',(NSIG3(K),K=NBINS,1,-1)
	PRINT '(1X,A,4X,14I6)'  ,'>10sig',(NSIG10(K),K=NBINS,1,-1)
	PRINT '(1X,A,4X,14I6)'  ,'>30sig',(NSIG30(K),K=NBINS,1,-1)
	PRINT *
C
	RETURN
	END


	SUBROUTINE CALC_DSPACE_RFACTOR(R1,R2,GOOF,NUNIQ,
	1				NSIG3,NSIG10,NSIG30, D_LO,D_HI, RSIG,DSPACES)
C
C Do the sums to calculate R_MERGE for intensities > RSIG * esd
C Ignore any reflections marked as outliers (ISIG < 0)
C
	PARAMETER (NSEQ_MAX=100000)
	REAL DSPACES(NSEQ_MAX),RSIG(NSEQ_MAX)
C
	PARAMETER (NDATA_MAX=2000000)
	INTEGER HKLMS(4,NDATA_MAX),IFILES(NDATA_MAX)
	REAL WAVS(NDATA_MAX),COUNTS(NDATA_MAX),DCOUNTS(NDATA_MAX)
	REAL TTHS(NDATA_MAX),XY_PIX(2,NDATA_MAX)
	COMMON /DATA_COM/ HKLMS,IFILES,WAVS,COUNTS,DCOUNTS,TTHS,XY_PIX,NDATA
C
	COMMON /SEQ_COM/ IFIRST(NSEQ_MAX),ILAST(NSEQ_MAX),COUNTS_SEQ(NSEQ_MAX),NSEQ
C
	REAL ICALC,IOBS,ISIG
	COMMON /LSQ_FUNC_COM/ ICALC(NDATA_MAX),IOBS(NDATA_MAX),ISIG(NDATA_MAX)
C
C Perform the usual sums to determine statistics.
C
	DOF=0.0
	NUNIQ=0
	NSUM=0
	SUM_I1=0.0
	SUM_D1=0.0
	SUM_S1=0.0
	SUM_I2=0.0
	SUM_D2=0.0
	NSIG3=0
	NSIG10=0
	NSIG30=0
	DO ISEQ=1,NSEQ
C Ignore sequence if incorrect d-spacing
	  IF( (DSPACES(ISEQ).LT.D_LO) .OR. (DSPACES(ISEQ).GE.D_HI) ) CYCLE
C Collect statistics for merged intensities
	  NUNIQ=NUNIQ+1
	  IF(RSIG(ISEQ) .GE. 3.0) NSIG3 =NSIG3+1
	  IF(RSIG(ISEQ) .GE.10.0) NSIG10=NSIG10+1
	  IF(RSIG(ISEQ) .GE.30.0) NSIG30=NSIG30+1
C Collect statistics for unmerged intensities
	  DO I=IFIRST(ISEQ),ILAST(ISEQ)
	    NSUM=NSUM+1
	    SUM_I1=SUM_I1+IOBS(I)
	    SUM_D1=SUM_D1+ABS(IOBS(I)-ICALC(I))
	    SUM_S1=SUM_S1+ABS(ISIG(I))
	    SUM_I2=SUM_I2+IOBS(I)**2/ISIG(I)**2
	    SUM_D2=SUM_D2+(IOBS(I)-ICALC(I))**2/ISIG(I)**2
	  ENDDO
C
	ENDDO
	DOF=NSUM-NUNIQ
C
C Ensure no divide by zero errors
C
	NSUM=MAX(1,NSUM)
	DOF=MAX(1.0,DOF)
	SUM_I1=MAX(1E-6,SUM_I1)
	SUM_I2=MAX(1E-6,SUM_I2)
C
C R1 is the unweighted R = SUM{|Ihkl -Imerge|} / SUM{Ihkl}
C R2 is the weighted-rms R = SQRT(  SUM{w * |Ihkl -Imerge|^2} / SUM{w * Ihkl^2}  )
C GOOF is the goodness-of-fit = SQRT(  SUM{w * |Ihkl -Imerge|^2} / DOF  )
C                   where w is the weight = 1.0 / variance(Ihkl)
C
	R1=        SUM_D1/SUM_I1
	R2=  SQRT( SUM_D2/SUM_I2 )
	GOOF=SQRT( SUM_D2/DOF )
C
C Convert R1 & R2 to percentages, and limit the range of all values
C
	R1    =MIN(99.9,100.0*R1)
	R2    =MIN(99.9,100.0*R2)
	GOOF  =MIN(9.9, GOOF )
C
	RETURN
	END
