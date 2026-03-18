global ImageFigure ImageInfo DisplayData
global LatticeInfo DataDir ModulateInfo
global ArgboxSet BatchFiles

// ================== Process image files in batch mode =================
//function LoadBatchFileList()
//function BatchOrientSpots(bMemory)
//function BatchOrientTwins(bMemory)
//function BatchArgonneBoxes(bOldMode,bTwinMode,bModulateMode)
//function BatchLaue4(bTwin,bMod)
//function TestIntegDmin()
//function PlotWavelengthFiles()
//======================= Status for Laue4 =========================
//function [bTwin,bMod]=CheckLaue4Status(bRun)
//function [bOK,sOut]=CheckLaue4Files(bStatus)
//===================== Status for batch files =====================
//function sStatusTable=GetBatchStatus()
//function OutputBatchStatus(sStatusTable)
//function bStatus=ParseBatchStatus(sStatusTable)
//function sOut=BatchStatusString(bStatus)
//function [bTwins,bModulate,levels,ilatts,icens,vols]=GetIndexStatus(sFile)


// ================== Process image files in batch mode =================

function LoadBatchFileList()
global BatchFiles
// Ask for file names, sort and check them, then set BatchFiles
// and DataDir, and load ImageInfo from the first *.tif file
//
// NB: Assumes file names are not case-sensitive, i.e. using Windows
//
// Clear BatchFiles in case of failure or cancel
  BatchFiles=[]
//
// Ask for batch files
  sFiles=uigetfile(["*.tif"],DataDir,"Select image files to process",%t)
  if( isempty(sFiles) ) then return; end
//
// Create an example string to use on failures
   example=[" ";
            "Examples of valid file names:"; ...
            "   C:\W1\SiC1.tif C:\W1\SiC2.tif C:\W1\SiC3.tif"; ...
            "   D:\Q3\d1_001.tif D:\Q3\d1_002.tif D:\Q3\d1_003.tif"]
//
// Cut sFiles[] into path,non-numeric base,numeric,extension into parts[]
  parts=[]
  for file_name=sFiles
    [pname,bname,ename]=fileparts(file_name)
    ipos=regexp(bname,"/[0-9]+$/")
    if(ipos == []) then
      AbortBox(["Files must have a base name followed by a numeric";example])
    end
    bname=strsplit(bname,ipos-1)
    parts($+1,1:4)=[pname,bname(1),bname(2),ename]
  end
//
// Save path and non-numeric base names for the first file
  path_name=parts(1,1)
  base_name=parts(1,2)
//
// Abort if path & non-numeric base names are not unique, or if not a .tif file
  parts=convstr(parts,"l")
  if( or( parts(:,1) ~= parts(1,1) ) ) then
    AbortBox(["Files must have the same path";example])
  elseif( or( parts(:,2) ~= parts(1,2) ) ) then
    AbortBox(["Files must have the same base name";example])
  elseif( or( parts(:,4) ~= ".tif" ) ) then
    AbortBox(["File must have a .tif extension";example])
  end
//
// Abort if numeric parts are invalid
  if( ~and(isnum(parts(:,3))) ) then
    AbortBox(["File names must end in a numeric part";example])
  end
//
// If path is not the default, change to it and try to load settings file
// bDminValid=%t if d_min, etc. settings are valid
  bDminValid=(ImageInfo.type ~= 0)
  if(DataDir ~= pathconvert(path_name)) then
    ChangeDataDir(path_name,%t)
    WriteIniFile()
    bDminValid=LoadSettingsFile()
  end
//
// Sort numeric parts as integers since uigetfile() sorts as text
  [nums,keys]=unique(strtod(parts(:,3)))
  sNums=parts(keys,3)

// Combine non-numeric and numeric parts to BatchFiles[]
  BatchFiles=base_name+sNums
//
// Log the start of batch mode and output the file list
  WriteLogFileHeader("Start Batch Mode")
  WriteLogFile("Files:")
  out=strcat(BatchFiles,",")
  for i=1:70:length(out)
    WriteLogFile( "  "+part(out,i:i+69) )
  end
//
// Load ImageInfo (and partially OrientInfo) from the first batch file and
// instrument setup file. Also load the reduced image, but ignore it
  ReadImageHeader(BatchFiles(1),bDminValid)
//
endfunction


function BatchOrientSpots(bMemory)
global LatticeInfo OrientInfo ModulateInfo
//
// "Refine Orientation" for Batch Files
//
// This routine changes OrientInfo & LatticeInfo(1), while ignoring twins.
//
// If bMemory is true, use the index in memory (possi nly the index
// from the last image. Otherwise, always start from the saved
// index for that image, and skip any image without a saved file.
//
//
// Complain and give up if tif files are missing (How?)
  if( ~and(isfile(BatchFiles+".tif"))  ) then
    AbortBox("Image (*.tif) files are missing")
  end
//
// Write headers to console and log file
  LockMenus("Orienting spots")
  mprintf("\nSTART: Orient Spots\n")
// Make left-justified array of "Files" and base names
  ilen=max(length(BatchFiles))
  sFiles=strncpy(["File",BatchFiles']+blanks(80),ilen)
  mprintf("%s  %s\n",sFiles(1), ...
        "Found Matched  Start Final   Cell dimensions")
  WriteLogFile("Start Refine Orientation")
//
// Loop over batch files, refining each one in turn
  for base_name=BatchFiles'
// Output to console the next left-justified base name
    sFiles(1)=[]
    mprintf("%s",sFiles(1))
//
// Load index from saved file, if that is the option
    if( ~bMemory )
      [LatticeInfo,OrientInfo,ModulateInfo]=ReadIndexFile(base_name+".idx")
      if(LatticeInfo == []) then
        mprintf("ERROR: Index file missing, or invalid\n")
        continue
      end
    end
//
// Give up if we don't have a UB matrix for the refinement
    if(LatticeInfo(1).level < 2) then
      mprintf("ERROR: Starting orientation is missing\n")
      continue
    end
//
// Calculate det0 (used for initial cell volume)
    det0=det( LatticeInfo(1).ub )
//
// Read in image and find spots
    ReadImageHeader(base_name,%t)
    RunFindSpots(15)
    if( size(DisplayData.found_spots,1) < 20 ) then
      mprintf("ERROR: Less than 20 observed spots\n")
      continue
    end
//
// Refine orientation in auto-mode
    [LatticeInfo(1),OrientInfo,text]= ...
                 RunOrientSpots(LatticeInfo(1),OrientInfo,1,[])
//
// Conserve initial cell volume and set level to fully oriented
    rfac=(det0/det(LatticeInfo(1).ub))^(1/3)
    LatticeInfo(1).ub=LatticeInfo(1).ub*rfac
    LatticeInfo(1).cell(1:3)=LatticeInfo(1).cell(1:3)/rfac
    LatticeInfo(1).level=3
//
// Output summary to console and log file
    toks=tokens(text)
    mprintf("%7s%7s%8s%6s   %.2f %.2f %.2f %.1f %.1f %.1f\n", ...
                                     toks',LatticeInfo(1).cell)
    str=msprintf("  %s %s %.3f %.3f %.3f %.2f %.2f %.2f", ...
                               toks([2,4])',LatticeInfo(1).cell)
    WriteLogFile(str)
//
// Update index file
    SaveIndexFile(base_name+".idx",%f)
//
  end
//
// Output a "completed" message
  printf("\nFINISH: Orient Spots\n")
  WriteLogFile("End Refine Orientation")
  UnlockMenus()
//
endfunction


function BatchOrientTwins(bMemory)
global LatticeInfo OrientInfo ModulateInfo
//
// "Refine Split Spots" for Batch Files
//
// This routine changes the UBs of the twin lattices, while LatticeInfo(1)
// and OrientInfo are unaffected.
//
// If bMemory is true, use the index in memory and the index
// from the last image. Otherwise, always start from the saved
// index for that image, and skip any image without a saved file.
//
//
// Write headers to console and log file
  LockMenus("Orienting twins")
  if( bMemory ) then
    stmp="Orient Split Spots (initial)"
  else
    stmp="Orient Split Spots (repeat)"
  end
  mprintf("\nSTART: %s\n",stmp)
// Make left-justified array of "Files" and base names
  ilen=max(length(BatchFiles))
  sFiles=strncpy(["File",BatchFiles']+blanks(80),ilen)
  mprintf("%s  %s\n",sFiles(1), ...
        "Pairs   rms for 1 & 2    Rotation angle and axis")
  WriteLogFile("Start "+stmp)
//
// If using memory, get rotation of twin UBs relative to default
  if( bMemory ) then
    rot21=LatticeInfo(2).ub/LatticeInfo(1).ub
    rot31=LatticeInfo(3).ub/LatticeInfo(1).ub
    rot21=rot21/det(rot21)^(1/3)
    rot31=rot31/det(rot31)^(1/3)
  end
//
// Loop over batch files, refining each one in turn
  for base_name=BatchFiles'
// Output to console the next left-justified base name
    sFiles(1)=[]
    mprintf("%s",sFiles(1))
//
// Read index file to get OrientInfo, and LatticeInfo in Linfo[]
    [Linfo,OrientInfo,ModulateInfo]=ReadIndexFile(base_name+".idx")
    if(Linfo == []) then
      AbortBox("Index file missing, or invalid")
    end
// If using memory:
    if( bMemory ) then
// Load twin lattices from default read in
      Linfo(2:3)=Linfo(1)
// Calculate twin UBs by rotation from the default
      Linfo(2).ub=rot21*Linfo(1).ub
      Linfo(3).ub=rot31*Linfo(1).ub
    end

// Abort on errors
    if(size(Linfo,1) ~= 3) then
      AbortBox("Require 2 twins to refine pair separation")
    elseif(min(Linfo(2:3).level) < 2) then
      AbortBox("Starting orientation is missing")
    end
//
// Read image header to get PHI value, etc.
    ReadImageHeader(base_name,%f)
//
// Refine pair separation using hardwired value
    max_hkl=5, wav_min=1.8, dxy_max=20, niter=4, bFitCell=%f
    [LatticeInfo,sOut]=RunOrientTwins(Linfo,OrientInfo, base_name+".tif", ...
                                    max_hkl,wav_min,dxy_max,niter,bFitCell)
//
// Set level of twin lattices to fully oriented
    LatticeInfo(2:3).level=3
//
// Output summary to console and log file
    [ndum,nspots,rms1,rms2]=msscanf(sOut,"%d %f %f")
    rms1=rms1*abs(OrientInfo.pixsize(1))
    rms2=rms2*abs(OrientInfo.pixsize(1))
    [hkl1,hkl2,angle]=CalcSplitRotation(LatticeInfo(2).ub,LatticeInfo(3).ub)
    hkl1=6*hkl1/max(abs(hkl1))
    str=msprintf("%4d%8.2f%5.2f (mm)%10.2f%7.1f%5.1f%5.1f",nspots,rms1,rms2,angle,hkl1')
    mprintf("%s\n",str)
    WriteLogFile(str)
//
// Update index file
    SaveIndexFile(base_name+".idx",%f)
//
  end
//
// Output a "completed" message
  printf("\nFINISH: Orient Split Spots\n")
  WriteLogFile("End Orient Split Spots")
  UnlockMenus()
//
endfunction


function BatchArgonneBoxes(bOldMode,bTwinMode,bModulateMode)
global LatticeInfo OrientInfo ModulateInfo
//
  LockMenus("Integrating intensities")
//
// Output headers to log file
  WriteLogFile("Start Integrate Spots (argonne_boxes)")
// Output options used
  if( bOldMode ) then
    text="  Options: Old-Mode"
  else
    text="  Options: New-Mode"
  end
  if( bTwinMode ) then
    text=text+", Twin-mode"
  end
  if( bModulateMode ) then
    text=text+", Modulated-mode"
  end
  WriteLogFile(text)
// Output integration parameter names and values used by the mode
  if( bOldMode ) then
    str=["model_cutoff:"  string(ArgboxSet.model_cutoff)
         "strong_cutoff:" string(ArgboxSet.strong_cutoff)
         "cont_level:"    string(ArgboxSet.cont_level)
         "core_mult:"     string(ArgboxSet.core_mult)
         "peak_mult:"     string(ArgboxSet.peak_mult)
         "neigh_toler:"   string(ArgboxSet.neigh_toler)
         "min_fill:"      string(ArgboxSet.min_fill)
         "d_min:"         string(ArgboxSet.d_min)
         "wav_min:"       string(ArgboxSet.wav_min)
         "wav_max:"       string(ArgboxSet.wav_max)]
    str=strcat(str,"","c")
    WriteLogFile( "    "+strcat(str(1:5)," ") )
    WriteLogFile( "    "+strcat(str(6:10)," ") )
  else
    str=["pfrac_target:"    string(ArgboxSet.pfrac_target)
         "peak_mult:"       string(ArgboxSet.peak_mult)
         "pfrac_uncertain:" string(ArgboxSet.pfrac_uncertain)
         "neigh_toler:"     string(ArgboxSet.neigh_toler)
         "model_center:"    string(ArgboxSet.model_center)
         "model_outer:"     string(ArgboxSet.model_outer)
         "cont_level:"      string(ArgboxSet.cont_level)
         "min_fill:"        string(ArgboxSet.min_fill)
         "recenter:"        string(ArgboxSet.recenter)
         "d_min:"           string(ArgboxSet.d_min)
         "wav_min:"         string(ArgboxSet.wav_min)
         "wav_max:"         string(ArgboxSet.wav_max)]
    str=strcat(str,"","c")
    WriteLogFile( "    "+strcat(str(1:4)," ") )
    WriteLogFile( "    "+strcat(str(5:9)," ") )
    WriteLogFile( "    "+strcat(str(10:12)," ") )
  end
// Header for batch file output to log file
  WriteLogFile("  File Model Pass(Full,Weak,PkPk,CoPk,CoCo,Twin)" + ...
                            " Fail(Edge,Neigh,Far,Pixl,Misc)" )
//
// Make left-justified array of "Files" and base names
  ilen=max(length(BatchFiles))
  sFiles=strncpy(["File",BatchFiles']+blanks(80),ilen)
//
// Output headers to console
  printf("\nSTART: Integrate Spots (argonne_boxes)\n\n")
  mprintf("%s%s  %s  %s\n",sFiles(1),"Model", ...
              "Pass(Full,Weak,PkPk,CoPk,CoCo,Twin)", ...
              "Fail(Edge,Neigh,Far,Pixl,Misc)" )
//
// Loop over batch files
  for base_name=BatchFiles'
// Output to console the next left-justified base name
    sFiles(1)=[]
    mprintf("%s",sFiles(1))
//
// Load ImageInfo (and partially OrientInfo) from the image file and
// instrument setup file (also loads the reduced image which is ignored)
    ReadImageHeader(base_name,%t)
// Load the corresponding index file (without log file), abort on failure
    [LatticeInfo,OrientInfo,ModulateInfo]=ReadIndexFile(base_name+".idx")
    if(LatticeInfo == []) then
      BugBox("Index file is missing")
    end
// Run argonne_boxes (with twin_dist=30 for twin mode)
    stats=RunArgonneBoxes(base_name,bOldMode,bTwinMode,bModulateMode, 30)
// Output the integration summary to console and log file
    mprintf("%4d%7d%5d%5d%5d%5d%5d%5d%7d%5d%5d%5d%5d%5d\n",stats')
    text=msprintf(" %d %d %d %d %d %d %d %d %d %d %d %d %d %d",stats')
    WriteLogFile("  "+base_name+text)
  end
//
// Output completion to console and log file
  printf("\nFINISH: Integrate Spots\n")
  WriteLogFile("End Integrate Spots")
//
  UnlockMenus()
endfunction


function BatchLaue4(bTwin,bMod)
global LatticeInfo OrientInfo ModulateInfo
//
// If no wavelength file, load the nominal one for the instrument
  if( isempty(fileinfo("laue4_wav.dat")) ) then
    LoadNominalWaveFile()
  end
//
// Read first index file to get symmetry info, complain and die if missing
  [Linfo,Oinfo,Minfo]=ReadIndexFile(BatchFiles(1)+".idx")
  if(Linfo == []) then
    AbortBox("First index file  is missing")
  end
//
// Update globals
  LatticeInfo=Linfo
  OrientInfo=Oinfo
  ModulateInfo=Minfo
//
// Remove the outlier files
  mdelete("laue4_bad.lis")
  mdelete("laue4_sus.lis")
//
// Make an array of the numeric part of the batch file names
  nums=[]
  for bname=BatchFiles'
    [istart,iend,match]=regexp(bname,"/[0-9]*$/")
    nums($+1)=strtod(match)
  end
  no_num_name=strncpy(BatchFiles(1),istart-1)
//
// Run laue4 using the "non-numeric name" and the array of numeric parts
// and flags for twins or modulations
  RunLaue4(no_num_name,nums, bTwin,bMod)
//
endfunction


function TestIntegDmin()
//
// Complain and die if some *.int files are missing
  if( ~and(isfile(BatchFiles+".int"))  ) then
    AbortBox("Spot intensity (*.int) files are missing")
  end
//
// Read the *.int files and concatenate d,esds,bfails arrays
  d=[]
  esds=[]
  bfails=[]
  for base_name=BatchFiles'
    [hkl,imul,wav,tth,xy,counts]=ReadIntFile(base_name+".int")
    d=[d; (wav ./ (2*sin( tth/114.592 ))) ]
    esds=[esds; ( counts(:,1) ./ max(counts(:,2),1) ) ]
    bfails=[bfails; (counts(:,2) == -9999) ]
  end
//
// Create suitable bins for the range in d-spaces
  d_min=min(d);
  d_max=max(d);
  bins=1 ./ (  1/d_min + 1.001*(1/d_max-1/d_min)*[0:0.1:1]  )
//
// Count how many spots of each type in each bin
// > 9 esds
  ibins=find( ~bfails(:) & (esds(:) > 9) )
  [ind, nbins(1,:), info] = dsearch(d(ibins), bins)
// > 5 esds
  ibins=find( ~bfails(:) & (esds(:) > 5) )
  [ind, nbins(2,:), info] = dsearch(d(ibins), bins)
// > 3 esds
  ibins=find( ~bfails(:) & (esds(:) > 3) )
  [ind, nbins(3,:), info] = dsearch(d(ibins), bins)
// failures
  ifail=find( bfails(:) )
  [ind, nbins(4,:), info] = dsearch(d(ifail), bins)
// total number
  [ind, ntot, info] = dsearch(d, bins)
//
// Output the bin results to the console in a nice table
  mprintf("            Summed results for all batch files\n\n")
  mprintf("       d-spacing     Num   Fails  >3 esd  >5 esd  >9 esd\n")
  for i=1:10
    ntoti=max(1,ntot(i))
    mprintf("     %5.2f -%5.2f%7i%7i%%%7i%%%7i%%%7i%%\n", ...
             bins(i),bins(i+1),ntot(i), ...
             100*nbins(4,i)/ntoti,100*nbins(3,i)/ntoti, ...
             100*nbins(2,i)/ntoti,100*nbins(2,i)/ntoti      )
  end
//
endfunction



function PlotWavelengthFiles()
//
//// Plot "old", "original", and possibly "new", spectra versus wavelength
//
//// Read wav/old/orig arrays from laue4_wav.dat
// Try to open laue4_wav.dat, complain and die if missing
  [fd,err]=mopen("laue4_wav.dat")
  if(err ~= 0) then
    AbortBox("laue4_wav.dat is missing")
  end
// Skip the header, read lines into buffer, close the file
  mgetl(fd,3)
  buffer=mgetl(fd,-1)
  mclose(fd)
// Decode buffer as 2 or 3 columns (wav,old,orig)
  for i=1:size(buffer,1)
    [npars wav(i) old(i) orig(i)]=msscanf(buffer(i),"%f %f %f");
  end
// If orig is missing, copy old to orig
  if(npars == 2) then
    orig=old
  end
//
//// Read new array from laue4_wav.new
// Try to open new wavelength file
  [fd,err]=mopen("laue4_wav.new")
// Start with empty new in case the file is missing
  new=[]
// If new file does exist, check if compatible and read new
  if(err == 0) then
// Skip the header, read lines into buffer, close the file
    mgetl(fd,3)
    buffer=mgetl(fd,-1)
    mclose(fd)
// Decode buffer as (wav2,new)
    for i=1:size(buffer,1)
      [npars wav2(i) new(i)]=msscanf(buffer(i),"%f %f");
    end
// Complain and die if wavs & wav2 are different
    if( ~and(wav==wav2) ) then
      AbortBox("laue4_wav.dat and laue4_wav.new are incompatible")
    end
  end
//
//// Plot old, original, and possibly new, values versus wavelength
// Change the figure for the plots
  f=CreateEmptyFigure([600,400])
  f.figure_name="Normalisation Wavelength Distributions"
  f.infobar_visible="off"
  f.background=-2
//
// Choose suitable wavelength range to plot
  wav_min=ArgboxSet.wav_min
  wav_max=wav_min*2.0
  wav_min=wav_min*0.8
  wav_max=wav_max*1.2
//
// If new array is empty, plot old/orig data versus wav
  if(new == []) then
    dat=[old,orig]
    plot2d(wav,dat,[1,2],leg="old@original",rect=[wav_min,0,wav_max,max(dat)])
//
// Otherwise, plot new/old/orig data versus wav
  else
    dat=[new,old,orig]
    plot2d(wav,dat,[1,2,3],leg="new@old@original",rect=[wav_min,0,wav_max,max(dat)])
  end
//
// Customise new/old/original plots lines
  a=gca()
  a.children(1).children(:).thickness=2
  a.children(1).children(1).foreground=9    // "original"
  a.children(1).children(1).line_style=7
  a.children(1).children(2).foreground=13   // "old"
  a.children(1).children(2).line_style=8
  if(new ~= []) then
    a.children(1).children(3).foreground=21   // "new"
    a.children(1).children(3).line_style=1
  end
// Customise plots legend
  a.children(2).font_size=3
  a.children(2).font_style=8
  a.children(2).legend_location="in_upper_left"
//
endfunction


//======================= Status for Laue4 =========================

function [bTwin,bMod]=CheckLaue4Status(bRun)
//
// Write number of batch files to sOut1[]
  sOut1=msprintf("Files to normalise: %d",size(BatchFiles,1))
// Add info about laue4_wav.dat
  if( isfile("laue4_wav.dat") ) then
    sOut1($+1)="Wavelength spectra file (laue4_wav.dat) exists"
  end
// Add info about laue4_rej.dat
  [fd,ierr]=mopen("laue4_rej.dat","r")
  if(ierr == 0) then
    buff=mgetl(fd,-1)
    mclose(fd)
    nrej=size(buff,1)-3
    sOut1($+1)=msprintf("%d Reject Spots found in laue4_rej.dat",nrej)
  end
//
// Get the status of memory and batch file info
  sStatusTable=GetBatchStatus()
  bStatus=ParseBatchStatus(sStatusTable)
//
// Check files for any errors
  [bOK,sOut2]=CheckLaue4Files(bStatus)
//
// Display any errors and abort
  if( ~bOK ) then
    OutputBatchStatus(sStatusTable)
    AbortBox([sOut2;"";sOut1;"";"Additional information on console window"])
  end
//
// Copy for return if intensities are Twinned or Modulated
  bTwin=bStatus.bIntIntTwin
  bMod=bStatus.bIntIntMod
//
// If bRun=%t, just return
  if( bRun ) then
    return
  end
//
// Create status summary for InfoBox display
  sOut3=BatchStatusString(bStatus)
//
// 
  if( bStatus.bIntIntTwin ) then
    sOut4="Normalisation ready to run in Twinned mode"
  elseif( bStatus.bIntIntMod ) then
    sOut4="Normalisation ready to run in Modulated mode"
  elseif( bStatus.bIntAll ) then
    sOut4="Normalisation ready to run in HKL mode"
  else
    BugBox("Invalid I.I. for sOut4")
  end
//
// Print status table to console window
  OutputBatchStatus(sStatusTable)
//
// Output results to InfoBox
  InfoBox([sOut4;"";sOut1;sOut2;"";sOut3;""; ...
             "Additional information on console window"])
//
endfunction


function [bOK,sOut]=CheckLaue4Files(bStatus)
//
  bOK=%f
  if( ~bStatus.bIdxFull ) then
    sOut="All batch files must be fully oriented"
  elseif( bStatus.bIdxBad ) then
    sOut="Files have inconsistent lattices or cell volumes"
  elseif( ~bStatus.bIntAll ) then
    sOut="All batch files must be integrated"
  elseif( bStatus.bIntIntMix ) then
    sOut="Integrated intensities must be the same type"
// Special case for twinned
  elseif( bStatus.bIntIntTwin ) then
    sOut="Integrated intensity type: Twinned"
    if( ~bStatus.bTwinFull ) then
      sOut(2)="Some twinned lattices not fully oriented"
    elseif( bStatus.bTwinBad ) then
      sOut(2)="Index files have inconsistent twin lattices or cell volumes"
    elseif( bStatus.bTwinSwap ) then
      sOut(2)="Use ""Orient/Fix Twin Indices"" to fix swapped twins"
    else
      bOK=%t
    end
// Special case for modulated
  elseif( bStatus.bIntIntMod ) then
    sOut="Integrated intensity type: Modulated"
    if( ~bStatus.bModAll ) then
      sOut(2)="Some index files are not modulated"
    elseif( bStatus.bModBad ) then
      sOut(2)="Index files have inconsistent modulations"
    else
      bOK=%t
    end
// Special case for normal HKL
  else
    sOut="Integrated intensity type: HKL (normal)"
    bOK=%t
  end
//
// Prepend ERROR to last line of sOut
  if( ~bOK ) then
    sOut($)="ERROR: "+sOut($)
  end
//
  endfunction


//===================== Status for batch files =====================

function sStatusTable=GetBatchStatus()
//
// Create table of strings for status of in-memory and batch-files
//
  sLatt=["","?????","Tricl","Mono ","Ortho","Tetra", + ...
      "Cubic","Trg-R","Trg-H","Hexag"]
  sCen=["","P","A","B","C","I","F","R"]
  sLev=["None","None","Cell","Part","Full"]
  YesNo=["N","Y"]
//
// In-memory status, put in sStatusTable(1,:)
  [bTwins,bMod,levels,ilatts,icens,vols]=GetIndexStatus([])
  sIdx=[sLev(levels+2);sLatt(ilatts+2);sCen(icens+2);string(round(vols))]
  sIdx=matrix(sIdx,1,-1)
  sMod=YesNo( double(1+bMod) )
  sStatusTable=["[in memory]","","","","",sIdx,sMod,""]
//
// Batch-files status, add to sStatusTable
  for sBase=BatchFiles'
// File existences
    bTif=isfile(sBase+".tif")
    bIdx=isfile(sBase+".idx")
    bInt=isfile(sBase+".int")
    bEll=isfile(sBase+".ell")
    sExists=YesNo( 1+double([bTif,bIdx,bInt,bEll]) )
// Detailed index file info (works if file is missing)
    [dum,bMod,levels,ilatts,icens,vols]=GetIndexStatus(sBase+".idx")
    sIdx=[sLev(levels+2);sLatt(ilatts+2);sCen(icens+2);string(round(vols))]
    sIdx=matrix(sIdx,1,-1)
    sMod=YesNo( 1+double(bMod) )
// Read integrated intensity mode from *.int header
    sIntType="-1"
    [fd,ierr]=mopen(sBase+".int","r")
    if(ierr == 0) then
      head=mgetl(fd,1)
      if(part(head,1:$-1) == "LaueG Intensities File Version 1, Option ") then
        sIntType=part(head,$)
      end
      mclose(fd)
    end
// Put all together in a table of strings
    sStatusTable($+1,:)=[sBase,sExists,sIdx,sMod,sIntType]
  end
//
// Blank out modulations if default level = "None"
  i=find(sStatusTable(:,6)  == "None" )
  if(i ~= []) then
    sStatusTable(1,18)=""
  end
// Blank out lattices where level = "None"
  i=find( sStatusTable == "None" )
  if(i ~= []) then
    nrows=size(sStatusTable,1)
    sStatusTable([i,i+nrows,i+2*nrows,i+3*nrows])=""
  end
// Blank out any "-1"
  sStatusTable( find(sStatusTable == "-1") )=""
//
// Check if twins are swapped compared to first file
  nsize=size(sStatusTable,1)
  sum_rot=zeros(nsize,1)
  rot0=zeros(3,3)
  Linfo=LatticeInfo
  for i=1:nsize
    if( or(sStatusTable(i,10) == ["Part","Full"]) & ...
                    or(sStatusTable(i,14) == ["Part","Full"]) ) then
      if(i > 1) then
        [Linfo,Oinfo,Minfo]=ReadIndexFile(sStatusTable(i,1)+".idx");
      end
      rot=Linfo(2).ub/Linfo(3).ub;
      if(rot0(1,1) == 0) then
        rot0=rot;
      end
      sum_rot(i)=sum( (rot-eye(3,3)) .* (rot0-eye(3,3)) )
    end
  end
// Ensure first file is non-negative by negating sum_rot
  if(sum_rot(2) < 0) then
    sum_rot=-sum_rot
  end
// Set column 20 to N for positive, Y for negative, blank for zero
  stmp=["Y","","N"]
  sStatusTable(:,20)=stmp( sign(sum_rot)+2 )'
//
// Rearrange around last 3 columns
  sStatusTable(:,[18,19,20])=sStatusTable(:,[20,18,19])
//
endfunction


function OutputBatchStatus(sStatusTable)
//
// Output status table, with explanation, to console
//
// Centre justify table
  sTable=justify(sStatusTable,'c')
// Swap some columns
  sTable(:,2:9)=sTable(:,[6:9,2:5])
// Convert empty columns to a single blank
  sTable( find(sTable == "") )=" "
// Concat columns separated by "|"
  sTable='|'+strcat(sTable,'|','c')+'|'
//
// Create a header of same width as table
  sBar=part("=",ones(1:200))
  sHead=part(sBar,1:15)+" Memory and Batch File Status "+sBar
  sHead=part(sHead, 1:length(sTable(1)) )
//
// Output header line and table
  mprintf("%s\n",[sHead;sTable])
//
// Output explanation
  mprintf("     %s\n", ...
  ["Table columns are:"; ...
   "  (1)    File name"; ...
   " (2-5)   Orientation status, lattice type, centering, cell volume"; ...
   " (6-9)   Existence of *.tif,*.idx,*.int,*.ell files"; ...
   "(10-17)  Same as (2-5) for Twin 1 & 2 lattices"; ...
   " (18)    Twins are swapped?"; ...
   " (19)    Modulated structure?"; ...
   " (20)    Integrated intensities: 1 (normal), 2 (twin), 3 (modul.)"] ...
  )
//
// Output trailing line of "="
  mprintf("%s\n",part(sBar, 1:length(sTable(1)) ))

//
endfunction


function bStatus=ParseBatchStatus(sStatusTable)
//
// Make local status table without blanks and all lower case
  sTable=convstr(stripblanks(sStatusTable))
//
// Do in-memory status, then remove it from local status table
  bStatus.bIdxMem=or(sTable(1,6) == ["part","full"])
  bStatus.bTwinMem=and( (sTable(1,[10,14]) == "part") | ...
                      (sTable(1,[10,14]) == "full") )
  bStatus.bModMem=(sTable(1,18) == "y")
  sTable(1,:)=[]
//
// Count the number of batch files
  nBatch=size(BatchFiles,1)
  nIdx=sum(sTable(:,3) == "y")
  nInt=sum(sTable(:,4) == "y")
  nEll=sum(sTable(:,5) == "y")
//
// Check how many .int & .ell files
  bStatus.bIntAll=(nInt+nEll == 2*nBatch)
  bStatus.bIntNone=(nInt+nEll == 0)
//
// Count number of fully oriented "default" (non-twin) lattices
  nFull=sum(sTable(:,6) == "full")
  nPart=sum(sTable(:,6) == "part")
// Count consistent lattices, centering, cell volumes for default lattice
  nLatt=sum(sTable(:,7) == sTable(1,7))
  nCen =sum(sTable(:,8) == sTable(1,8))
  vols=strtod(sTable(:,9))
  nVol=sum(abs(vols-vols(1)) < 10)
//
// Check default lattice in index files
  bStatus.bIdxNone=(nIdx == 0)
  bStatus.bIdxFull=(nFull == nBatch)
  bStatus.bIdxPart=(nPart+nFull == nBatch)
  bStatus.bIdxBad=(nLatt+nCen+nVol ~= 3*nBatch)
//
// Repeat above for twin lattices
  bStatus.bTwinNone=or(sTable(:,[10,14]) == "")
  bStatus.bTwinFull=and(sTable(:,[10,14]) == "full")
  bStatus.bTwinPart=and( (sTable(:,[10,14]) == "part") | ...
                         (sTable(:,[10,14]) == "full")     )
  nLatt=sum(sTable(:,[11,15]) == sTable(1,7))
  nCen =sum(sTable(:,[12,16]) == sTable(1,8))
  vols1=strtod(sTable(:,13))
  vols2=strtod(sTable(:,17))
  nVol=sum(abs([vols1,vols2]-vols(1)) < 10)
  bStatus.bTwinBad=(nLatt+nCen+nVol ~= 6*nBatch)
//
// Are twins a mix of swapped and not swapped
  bStatus.bTwinSwap=or( sTable(:,18) ~= sTable(1,18) )
//
// Modulated structures
  nMod=sum(sTable(:,19) == "y")
  bStatus.bModNone=(nMod == 0)
  bStatus.bModAll=(nMod == nBatch)
// Check that modulations are consistent
  bStatus.bModBad=%f
  if( bStatus.bModAll ) then
    [Linfo,Oinfo,Minfo1]=ReadIndexFile(sTable(2,1)+".idx")
    bStatus.bModBad=%t
    for sBase=sTable(2:$,1)'
      [Linfo,Oinfo,Minfo2]=ReadIndexFile(sBase+".idx")
      if(norm(Minfo1.vecs-Minfo2.vecs)+norm(Minfo1.mults-Minfo2.mults) > 1e-3) then
        bStatus.bModBad=%f
      end
    end
  end
//
// Integrated intensity types
  types=unique( sTable(:,20) )
  bStatus.bIntIntMix= (size(types,2) > 1)
  bStatus.bIntIntTwin= (~bStatus.bIntIntMix) & (types(1) == "2")
  bStatus.bIntIntMod= (~bStatus.bIntIntMix) & (types(1) == "3")
//
endfunction


function sOut=BatchStatusString(bStatus)
//
// In-memory status
  if( bStatus.bIdxMem ) then
    sOut="Orientation solution found in memory"
  else
    sOut="No orientation solution in memory"
  end
  if( bStatus.bTwinMem ) then
    sOut($+1)="Twin lattice information found in memory"
  elseif( bStatus.bModMem ) then
    sOut($+1)="Modulation information found in memory"
  else
    sOut($+1)="No twin or modulation information found in memory"
  end
//
// Index files
// Normal HKLs
  sOut($+1)=""
  if( bStatus.bIdxPart ) then
    if( bStatus.bIdxFull ) then
      sOut($+1)="All index files fully oriented for normal HKLs"
    else
      sOut($+1)="All index files at least partially oriented"
    end
    if( bStatus.bIdxBad ) then
      sOut($+1)="Index files have inconsistent lattices or cell volumes"
    end
  elseif( bStatus.bIdxNone ) then
    sOut($+1)="No batch files oriented"
  else
    sOut($+1)="Some batch files not oriented"
  end
// Twins
  if( bStatus.bTwinPart ) then
    if( bStatus.bTwinFull ) then
      sOut($+1)="All index files fully oriented for twins"
    else
      sOut($+1)="All index files at least partially oriented for twins"
    end
    if( bStatus.bTwinBad ) then
      sOut($+1)="Index files have inconsistent twin lattices or cell volumes"
    elseif( bStatus.bTwinSwap ) then
      sOut($+1)="Some index files have swapped twins"
    end
  elseif( ~bStatus.bTwinNone ) then
    sOut($+1)="Some index files have oriented twins"
  end
// Modulations
  if( bStatus.bModAll ) then
    sOut($+1)="All index files have modulation information"
    if( bStatus.bModBad ) then
      sOut($+1)="Index files have inconsistent modulations"
    end
  elseif(~bStatus.bModNone ) then
    sOut($+1)="Some index files have modulation information"
  end
// No twins or modulation
  if( ~bStatus.bTwinPart & ~bStatus.bModAll ) then
    sOut($+1)="No twin or modulation information in index files"
  end
//
// Integration files and intensity types of *.int files
  sOut($+1)=""
  if( bStatus.bIntNone ) then
    sOut($+1)="No integration files exist"
  elseif( ~bStatus.bIntAll ) then
    sOut($+1)="Some integration files are missing"
  else
    if( bStatus.bIntIntMix ) then
      sOut($+1)="Integration files have mixed types"
    elseif( bStatus.bIntIntTwin ) then
      sOut($+1)="All integration files exist for Twinned mode"
    elseif( bStatus.bIntIntMod ) then
      sOut($+1)="All integration files exist for Modulated mode"
    else
      sOut($+1)="All integration files exist for normal HKL mode"
    end
  end
//
endfunction


function [bTwins,bModulate,levels,ilatts,icens,vols]=GetIndexStatus(sFile)
//
// Return values in case of failure
  bModulate=%f
  levels=[-1,-1,-1]
  ilatts=[-1,-1,-1]
  icens =[-1,-1,-1]
  vols  =[-1,-1,-1]
//
// Load Linfo & Minfo from memory, or from an index file
  if(sFile == []) then
    Linfo=LatticeInfo
    Minfo=ModulateInfo
  else
// Try to read index file, return bTwins=[] if failed
    [Linfo,Oinfo,Minfo]=ReadIndexFile(sFile)
    if(Linfo == []) then
      bTwins=[]
      return
    end
  end
//
// Check for modulation information
  bModulate= ~isempty(Minfo)
  if( bModulate ) then
    bModulate= ~isempty(Minfo.vecs)
  end
//
// If 3 lattices, set bTwins
  bTwins=( size(Linfo,1) == 3 )
//
// Copy Linfo[].level/ilatt/icen values to return arrays
  for i=1:size(Linfo,1)
    levels(i)=Linfo(i).level
    ilatts(i)=Linfo(i).ilatt
    icens (i)=Linfo(i).icen
    if(Linfo(i).ub ~= []) then
      vols(i)=1/abs(det( Linfo(i).ub ))
    end
  end
//
endfunction
