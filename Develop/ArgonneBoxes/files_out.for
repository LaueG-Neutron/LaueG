C ================ Routines in this file =========================
C	SUBROUTINE OPEN_LIST_FILE(FILE_NAME)
C	SUBROUTINE OUTPUT_RESULTS(BASE_NAME,SCALE,NREFS)
C	SUBROUTINE OUTPUT_ELL_RECORD(IR, E,F,H,RAD_MULT)
C	SUBROUTINE OUTPUT_INT_FILES(RINTINT,RSIGMA,IOPTION,IOVER,NREFS, SCALE, FILE_NAME)
C	SUBROUTINE CHECK_SPOTS_PFRAC()
C------------ Debug routines ------------
C	SUBROUTINE WRITE_BACK_SMOOTH_IMAGES
C ================================================================

	SUBROUTINE OPEN_LIST_FILE(FILE_NAME)
C
	CHARACTER FILE_NAME*(*)
C
	CHARACTER VERSION_STRING*80
	COMMON /VERSION_STRING_COM/ VERSION_STRING
C
C Open the main log file and write the program name and version
C
	OPEN(UNIT=30,FILE=FILE_NAME,STATUS='UNKNOWN')
	WRITE(30,'(10X,A,5X,A,/)') 'ARGONNE_BOXES',TRIM(VERSION_STRING)
C
C Write an acknowledgement banner to the main log file
C
C Please do not remove the following output statement
	WRITE(30,'(A)')
     1  '========================================================================',
     2  ' Published use of this program should refer to:',
     3  ' "Integration of Single-Crystal Reflections using Area Multidetectors",',
     4  '    C. Wilkinson, H.W. Khamis, R.F.D. Stansfield and G.J. McIntyre,',
     5  '    J. Appl. Cryst. 21 (1988) 471-478.'
C
C Update the following as required
C
	WRITE(30,'(A)')
     1  ' Code by C.Wilkinson, modified by G.J.McIntyre and R.O.Piltz',
     2  ' This version maintained by R.O.Piltz (rop@ansto.gov.au)',
     3  '========================================================================'
C
	RETURN
	END


	SUBROUTINE OUTPUT_RESULTS(BASE_NAME,SCALE,NREFS)
C
	CHARACTER BASE_NAME*80
C
	COMMON /OPTIONS_COM/ IOPT_HKL,NHKLM
	COMMON /LIB_COM/ PLIB(1000,6),ILIB(1000),NLIB
C
	COMMON /SPOT2_COM/ IHKLM(4,20000),WAV(20000),TTH(20000),MULT(20000)
	COMMON /INTINT_COM/ RINTINT(20000),RSIGMA(20000),IOPTION(20000),IOVER(20000)
C
	COMMON /INT_COUNTERS_COM/ NSUCCESS,NFAIL,NFULL,NWEAK1,NWEAK2,NWEAK3,
	1                      NWEAK4, NPAIR,NEDGE,NNOBACK,NPIXEL,NNEIGH,NFAR,NELLIP
C
C Output a summary of results
C
	WRITE(30,'(A)') '**************************************'
	WRITE(30,'(/,A)') 'Integration Summary:'
	WRITE(30,'(I6,A)') NFULL,' full peak integrations'
	WRITE(30,'(I6,A)') NSUCCESS-NFULL,' sig(I)/I integrations'
	WRITE(30,'(I6,A)') NREFS-NSUCCESS,' failed integrations'
C
C Output total successful, and number of overlapped twins or satellites
C
	IF(IOPT_HKL .EQ. 1) THEN
	  WRITE(30,'(/,A,I6)') 'Successful spot integrations:',NSUCCESS
	ELSEIF(IOPT_HKL .EQ. 2) THEN
	  NOVER=0
	  DO IR=1,NREFS
	    IF(IOPTION(IR) .GT. 0) THEN
	      IF(IHKLM(4,IR) .GE. 10) NOVER=NOVER+1
	    ENDIF
	  ENDDO
	  WRITE(30,'(/,A,2(/,I6,A))') 'Successful spot integrations:',
	1		NSUCCESS-NOVER,' non-overlapping or separable twins',
	1		NOVER,' overlapping and inseparable twins'
	ELSE
	  NSAT=0
	  DO IR=1,NREFS
	    IF(IOPTION(IR) .GT. 0) THEN
	      IF(IHKLM(4,IR) .GE. 1) NSAT=NSAT+1
	    ENDIF
	  ENDDO
	  WRITE(30,'(/,A,2(/,I6,A))') 'Successful spot integrations:',
	1		NSUCCESS-NSAT,' normal reflections',
	1		NSAT,' satellite reflections'
	ENDIF
C
C Check, and output to list file, how well the estimated peak-fractions worked
C
	CALL CHECK_SPOTS_PFRAC()
C
C Output the *.int and *.flg files.
C All intensities and esd's are divided by SCALE.
C
	WRITE(30,'(//,A)') '================ Intensity Output ==============='
	CALL OUTPUT_INT_FILES(RINTINT,RSIGMA,IOPTION,IOVER,NREFS, SCALE,TRIM(BASE_NAME))
C
C Close log and *.ell files
C
	CLOSE(30)
	CLOSE(40)
C
C Write summary of results to console and the file "___laueg_argonne_boxes.out"
C
	NMISC=NNOBACK+NELLIP
	PRINT '(3(1X,A))','Model','Pass(Full,Weak,PkPk,CoPk,CoCo,Twin)',
	1			'Fail(Edge,Neigh,Far,Pixl,Misc)'
	PRINT '(1X,I4,I6,6I5,1X,6I5)',NLIB,
	1		NSUCCESS,NFULL,NWEAK1,NWEAK2,NWEAK3,NWEAK4,NPAIR,
	2		NFAIL,NEDGE,NNEIGH,NFAR,NPIXEL,NMISC
C
	OPEN(UNIT=1,FILE='___laueg_argonne_boxes.out',STATUS='UNKNOWN',IOSTAT=IDUMMY)
	WRITE(1,'(14I7)',IOSTAT=IDUMMY) NLIB,NSUCCESS,NFULL,NWEAK1,NWEAK2,
	1		NWEAK3,NWEAK4,NPAIR,NFAIL,NEDGE,NNEIGH,NFAR,NPIXEL,NMISC
	CLOSE(UNIT=1,IOSTAT=IDUMMY)
C
	RETURN
	END


	SUBROUTINE OUTPUT_ELL_RECORD(IR, E,F,H,RAD_MULT)
C
	REAL RAD_MULT(3)
C
	COMMON /OPTIONS_COM/ IOPT_HKL,NHKLM
	COMMON /LIB_COM/ PLIB(1000,6),ILIB(1000),NLIB
C
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /SPOT2_COM/ IHKLM(4,20000),WAV(20000),TTH(20000),MULT(20000)
	COMMON /INTINT_COM/ RINTINT(20000),RSIGMA(20000),IOPTION(20000),IOVER(20000)
C	
C Determine spot type
C  =1   full integration
C  =2   sigI/I integration, no overlap 
C  =3-5 sigI/I integration, (P-P,C-P,C-C) overlap
C  =6   sigI/I integration, overlapped twin
C  =9   failed integration
C  >10  same as status-10 but also a model spot
C
		IF(IOPTION(IR) .LT. 0) THEN
		  ITYPE=9
		ELSEIF(IOPTION(IR) .EQ. 3) THEN
		  ITYPE=1
		ELSEIF(IOVER(IR) .EQ. 5) THEN
		  ITYPE=6
		ELSE
		  ITYPE=IOVER(IR)+1
		ENDIF
C
		DO I=1,NLIB
		  IF(IR .EQ. ILIB(I)) ITYPE=ITYPE+10
		ENDDO
C
C Output ellipse info to the *.ell file
C NB: M for modulated spots is ignored
C
	  IF(IOPT_HKL .EQ. 1) THEN
	  WRITE(40,'(3I4,2F7.1,1P,3E11.3,0P,1X,3F4.1,I3)')
	1              (IHKLM(K1,IR),K1=1,3),X(IR),Y(IR),E,F,H,
	2              (RAD_MULT(K2)**2,K2=1,3),ITYPE
C
	ELSE
	    IM=IHKLM(4,IR)
C Convert M indices to 1,2 for separated twins and 11,12 for overlapped twins
	    IF(IOPT_HKL .EQ. 2) THEN
	      IM=IM-IM/10*10
	      IF(IHKLM(4,IR) .GT. 50) IM=IM+10
	    ENDIF
	  WRITE(40,'(3I4,2F7.1,1P,3E11.3,0P,1X,3F4.1,2I3)')
	1              (IHKLM(K1,IR),K1=1,3),X(IR),Y(IR),E,F,H,
	2              (RAD_MULT(K2)**2,K2=1,3),ITYPE,IM
	ENDIF
C
	RETURN
	END


	SUBROUTINE OUTPUT_INT_FILES(RINTINT,RSIGMA,IOPTION,IOVER,NREFS, SCALE, FILE_NAME)
C
	CHARACTER FILE_NAME*(*)
	REAL RINTINT(20000),RSIGMA(20000)
	INTEGER IOPTION(20000),IOVER(20000)
C
	COMMON /BGLOCAL_COM/ BKG(20000),BKG_ESD(20000),BKG_VAR(20000)
	COMMON /SPOT1_COM/ X(20000),Y(20000),RMAX(20000)
	COMMON /SPOT2_COM/ IHKLM(4,20000),WAV(20000),TTH(20000),MULT(20000)
	COMMON /OPTIONS_COM/ IOPT_HKL,NHKLM
C
C Output final info to the log file
C
	WRITE(30,'(/,2A)') 'Integration status flags written to ',FILE_NAME//'.flg'
	WRITE(30,'(2A)') 'Intensities for Laue4 written to ',FILE_NAME//'.int'
	WRITE(30,'(2X,A,F7.2)') 'Output intensities & esd''s divided by',SCALE
C
C Output the count values and the meaning of the M index
C
	IF(IOPT_HKL .NE. 1) WRITE(30,'(/,A)') 'Intensities are assigned H,K,L,M indices'
C
	IF(IOPT_HKL .EQ. 2) THEN
	  WRITE(30,'(A)') 'For non-overlapping or separable spots, '//
	1		'M=1 or 2 is the twin number.',
	2		'For overlapping and inseparable spots, M=11 is the total intensity',
	3		'of both twins, M=12 is the approximate intensity for twin 2.'
	ELSEIF(IOPT_HKL .EQ. 3) THEN
	  WRITE(30,'(2X,A)') 'M is the satellite number'
	ENDIF
C
C Open the intensity output file on unit=1
C	
	OPEN(UNIT=1,FILE=FILE_NAME//'.int',STATUS='UNKNOWN')
	IF(IOPT_HKL .EQ. 1) THEN
	  WRITE(1,'(A,/)') 'LaueG Intensities File Version 1, Option 1'
	  WRITE(1,'(A)') '   h   k   l   Wav  Mult Xpix   Ypix    TTH    Counts  esd'
	ELSEIF(IOPT_HKL .EQ. 2) THEN
	  WRITE(1,'(A,/)') 'LaueG Intensities File Version 1, Option 2'
	  WRITE(1,'(A)') '   h   k   l   m   Wav  Mult Xpix   Ypix    TTH    Counts  esd'
	ELSE
	  WRITE(1,'(A,/)') 'LaueG Intensities File Version 1, Option 3'
	  WRITE(1,'(A)') '   h   k   l   m   Wav  Mult Xpix   Ypix    TTH    Counts  esd'
	ENDIF
C	
C Open the integration status flag file on unit=2
C
	OPEN(UNIT=2,FILE=FILE_NAME//'.flg',STATUS='UNKNOWN')
	WRITE(2,'(3X,A,/)') 'ARGONNE_BOXES integration status and flags file'
	WRITE(2,'(A)')
     1	'Overlap number 1 = clear, no spot overlap',
     2	'               2 = peak overlap with neighbouring peak',
     3	'               3 = core overlap with neighbouring peak',
     4	'               4 = core overlap with neighbouring core',
     5	'               5 = overlapping twin integration'
	WRITE(2,*)
	WRITE(2,'(A)')
     1	'Integration option -1 = spot on edge or join',
     2	'                   -2 = pixel overload',
     3	'                   -3 = strong neighbour overlap',
     4	'                   -4 = invalid ellipse or background calc',
     5	'                   -5 = strong spot far from calc position',
     6	'                    1 = sig(i)/i due to weakness',
     7	'                    2 = sig(i)/i due to peak overlap',
     8	'                    3 = full peak integration from model'
	IF(IOPT_HKL .EQ. 1) THEN
	  WRITE(2,'(/,2A)') '  Ref num   h     k    l mult  x     y   ',
     1			'over  opt   wav        Int  sigma'
	ELSE
	  WRITE(2,'(/,2A)') '  Ref num   h     k    l mult  x     y   ',
     1			'over  opt   wav        Int  sigma   M'
	ENDIF
C
C Output to intensities and status/flags files
C
	DO IR=1,NREFS
C Convert M indices to 1,2 for separated twins and 11,12 for overlapped twins
	  IF(IOPT_HKL .EQ. 2) THEN
	    ITWIN=IHKLM(4,IR) - IHKLM(4,IR)/10*10
	    IF(IHKLM(4,IR) .GT. 50) ITWIN=ITWIN+10
	    IHKLM(4,IR)=ITWIN
	  ENDIF
C Scale the intensities & esd's by SCALE (except the failures)
	  IF(IOPTION(IR) .LE. 0) THEN
	    INTP=-9999
	    ISIGP=-9999
	  ELSE
	    INTP=NINT( RINTINT(IR)/SCALE )
	    ISIGP=NINT( RSIGMA(IR)/SCALE )
	  ENDIF
C Write a line to the *.int file
	  WRITE(1,'(<NHKLM>I4,F7.3,I3,F8.1,F7.1,F8.2,I7,I6)') (IHKLM(K,IR),K=1,NHKLM),
	1		WAV(IR),MULT(IR),X(IR),Y(IR),TTH(IR),INTP,ISIGP
C Write a line to the *.flg file (include M index if IOPT_HKL>1)
	  IF(IOPT_HKL .EQ. 1) THEN
	    WRITE(2,'(I5,4X,3I5,I3,2I6,2I5,F8.3,I10,I6)')  IR,
     1		(IHKLM(J,IR),J=1,3),MULT(IR),NINT(X(IR)),NINT(Y(IR)),
     2		IOVER(IR),IOPTION(IR),WAV(IR),INTP,ISIGP
	  ELSE
	    WRITE(2,'(I5,4X,3I5,I3,2I6,2I5,F8.3,I10,I6,I4)',IOSTAT=IDUMMY) IR,
     1		(IHKLM(J,IR),J=1,3),MULT(IR),NINT(X(IR)),NINT(Y(IR)),
     2		IOVER(IR),IOPTION(IR),WAV(IR),INTP,ISIGP,IHKLM(4,IR)
	  ENDIF
C
	ENDDO
C
C Close the output files
C
	CLOSE(UNIT=1)
	CLOSE(UNIT=2)
C
	RETURN
	END


	SUBROUTINE CHECK_SPOTS_PFRAC()
C
	COMMON /OPTIONS_COM/ IOPT_HKL,NHKLM
      COMMON /PFRAC_CHECK_COM/ PFRAC_LIST(3,20000),PFRAC_TARGET,NPFRAC
C
	WRITE(30,'(/,A)') 'Comparison of observed to estimated peak-fractions:'
	IF(IOPT_HKL .EQ. 2) WRITE(30,'(A)') '  WARNING: These statistics are unreliable in Twin mode'
	WRITE(30,'(I5,A)') NPFRAC,' suitable peak-fractions from peak integration'
C================================================================
C
C Remove outliers of RAT = integration_esd(Pfrac) / model_esd(Pfrac)
C
	NUM1=NPFRAC
	DO ITER1=1,20
C
C Calculate average and st-dev of RAT
	  AVE=0.0
	  SDEV=0.0
	  DO I=1,NUM1
	    RAT=PFRAC_LIST(2,I)/PFRAC_LIST(3,I)
	    AVE=AVE+RAT
	    SDEV=SDEV+RAT**2
	  ENDDO
	  AVE=AVE/MAX(1,NUM1)
	  SDEV=SQRT( SDEV/MAX(1,NUM1) - AVE**2 )
C
C Remove PFRAC_LIST() values with RAT outside AVE +/- 2*ST_DEV
	  NUM2=0
	  DO I=1,NUM1
	    RAT=PFRAC_LIST(2,I)/PFRAC_LIST(3,I)
		IF(RAT.GE.AVE-2*SDEV .AND. RAT.LE.AVE+2*SDEV) THEN
		  NUM2=NUM2+1
	      PFRAC_LIST(1,NUM2)=PFRAC_LIST(1,I)
	      PFRAC_LIST(2,NUM2)=PFRAC_LIST(2,I)
	      PFRAC_LIST(3,NUM2)=PFRAC_LIST(3,I)
	    ENDIF
	  ENDDO
C
C Stop if iterations complete or lost half the values
	  IF(NUM1.EQ.NUM2 .OR. NUM2.LT.NPFRAC/2) EXIT
	  NUM1=NUM2
	ENDDO
	NPFRAC=NUM1
C
	WRITE(30,'(I5,A,I3,A)') NPFRAC,' after',ITER1,' iterations removing esd outliers'
	WRITE(30,'(A,F5.2,A,F4.2,A)',IOSTAT=IDUMMY) '    esd ratio =',AVE,'(',SDEV,')'
C
C================================================================
C
C Remove outliers of RAT = (Pfrac - Target) / model_esd(Pfrac)
C
C Skip if less than 3 points
      IF(NPFRAC .GE. 3) THEN
C
C Sort PFRAC_LIST() in terms of RAT
        CALL SORT_PFRAC_LIST()
C
C Iteration removing outliers
	  NUM=NPFRAC
        DO ITER2=1,10
C Do linear fit (in indices) to RAT values
	    SUMX=0.0
	    SUMX2=0.0
	    SUMY=0.0
	    SUMXY=0.0
	    DO I=1,NUM
	      RAT=(PFRAC_LIST(1,I)-PFRAC_TARGET)/PFRAC_LIST(3,I)
	      SUMX=SUMX+I
	      SUMX2=SUMX2+I**2
	      SUMY=SUMY+RAT
	      SUMXY=SUMXY+I*RAT
          ENDDO
	    RAT_ZERO=(SUMY*SUMX2-SUMX*SUMXY)/(NUM*SUMX2-SUMX**2)
	    RAT_GRAD=(NUM*SUMXY-SUMX*SUMY)/(NUM*SUMX2-SUMX**2)
C Find low end above linear fit
	    DO I1=1,NUM/10
	      RAT=(PFRAC_LIST(1,I1)-PFRAC_TARGET)/PFRAC_LIST(3,I1)
	      IF(RAT .GT. RAT_ZERO+RAT_GRAD*I1) EXIT
	    ENDDO
C Find high end below linear fit
	    DO I2=NUM,NUM*9/10,-1
	      RAT=(PFRAC_LIST(1,I2)-PFRAC_TARGET)/PFRAC_LIST(3,I2)
	      IF(RAT .LT. RAT_ZERO+RAT_GRAD*I2) EXIT
          ENDDO
C Exit before pruning if result would be less than 3 or 67%
          IF(I2-I1+1 .LT. MAX(3,NPFRAC*2/3)) EXIT         
C Prune PFRAC_LIST() to remove leading/trailing outliers
	    DO I=I1,I2
	      PFRAC_LIST(1,I-I1+1)=PFRAC_LIST(1,I)
	      PFRAC_LIST(2,I-I1+1)=PFRAC_LIST(2,I)
	      PFRAC_LIST(3,I-I1+1)=PFRAC_LIST(3,I)
	    ENDDO
	    NUM=I2-I1+1
C
        ENDDO
	  NPFRAC=NUM
C
	  WRITE(30,'(I5,A,I3,A)') NPFRAC,' after',ITER2,' iterations removing (value-target)/esd outliers'
C
C From the last linear fit of RAT = RAT_ZERO+RAT_GRAD*I for I=I1 to I2
C
	  AVE=RAT_ZERO+RAT_GRAD*(I1+I2)/2
	  SDEV=0.42*RAT_GRAD*(I2-I1)
	  WRITE(30,'(A,F5.2,A,F4.2,A)',IOSTAT=IDUMMY) '    Average RAT2 =',AVE,'(',SDEV,')'
C
      ENDIF
C================================================================
C
C Calculate average and st-dev of peak-fraction
C
	AVE=0.0
	SDEV=0.0
	DO I=1,NPFRAC
	  AVE=AVE+PFRAC_LIST(1,I)
	  SDEV=SDEV+PFRAC_LIST(1,I)**2
	ENDDO
	NUM=MAX(1,NPFRAC)
	AVE=AVE/NUM
	SDEV=SQRT( SDEV/NUM - AVE**2 )
C
	WRITE(30,'(A,F5.2,A,F4.2,A)',IOSTAT=IDUMMY) '    Averaged peak-fraction =',AVE,'(',SDEV,')'
C
	RETURN
	END


C------------ Debug routines ------------

	SUBROUTINE WRITE_BACK_SMOOTH_IMAGES
C
	CHARACTER*2000 COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,ADC_INFO
	COMMON /TIFFINFO_COM/ COMMENT,SAMPLE,USER,DATETIME,HOST,BEAM,
	1			IOFF_DATA,NUMX,NUMY,IXYRES,ISTARTX,ISTARTY,ISPEED,
	2			EXPOSE_TIME,EXPOSE_PHI,TEMP_BEGIN,TEMP_END,TEMP_MIN,TEMP_MAX,
	3			IADC,ADC_ADD,ADC_DIV,ADC_INFO,IAPER_1,IAPER_2,IKAPPA
C
	COMMON /BACK_COM/ BACK(8000,2500)
	COMMON /SMOOTH_COM/ SMOOTH_IMAGE(8000,2500)
C
	INTEGER*2 IMAGE2_OUT(8000*2500)
C
	WRITE(30,'(/,A)') '*** DEBUG: Output background image to _BACK_IMAGE.TIF'
	WRITE(30,'(A)') '*** DEBUG: Output smoothed image to _SMOOTH_IMAGE.TIF'
C
C Convert background image to I*2 and output
C
	IOUT=0
	DO IY=1,NUMY
	  DO IX=1,NUMX
	    IOUT=IOUT+1
	    ISTRIP=MAX(0, NINT(BACK(IX,IY)) )
	    IF(ISTRIP .GE. 128*256) ISTRIP=ISTRIP-256*256
	    IMAGE2_OUT(IOUT)=ISTRIP
	  ENDDO
	ENDDO
C
	CALL WRITE_LAUE_TIFF('_BACK_IMAGE.TIF', IMAGE2_OUT,NUMX,NUMY,
	1			TRIM(HOST),TRIM(USER),TRIM(SAMPLE),
	2			DATETIME,TRIM(COMMENT),EXPOSE_PHI)
C
C Convert smoothed image to I*2 and output
C
	IOUT=0
	DO IY=1,NUMY
	  DO IX=1,NUMX
	    IOUT=IOUT+1
	    ISTRIP=MAX(0, NINT(SMOOTH_IMAGE(IX,IY)) )
	    IF(ISTRIP .GE. 128*256) ISTRIP=ISTRIP-256*256
	    IMAGE2_OUT(IOUT)=ISTRIP
	  ENDDO
	ENDDO
C
	CALL WRITE_LAUE_TIFF('_SMOOTH_IMAGE.TIF', IMAGE2_OUT,NUMX,NUMY,
	1			TRIM(HOST),TRIM(USER),TRIM(SAMPLE),
	2			DATETIME,TRIM(COMMENT),EXPOSE_PHI)
C
	RETURN
	END
