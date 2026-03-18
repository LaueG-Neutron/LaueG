global ImageFigure ImageInfo DisplayData
global LatticeInfo DataDir ModulateInfo
global ArgboxSet BatchFiles

// ============  Menu definitions and callbacks ============
//function SetupMainMenus()
//function MainStart_Menu(iMenu)
//function MainOrient_Menu(iMenu)
//function MainInteg_Menu(iMenu)
//function MainNorm_Menu(iMenu)
//function MainTools_Menu(iMenu)
//function MainHelp_Menu(iMenu)
// ==================== Helper routines  ================
//function DoBatchChecks(bFails,sFails)
//function LockMainMenus()
//function UnlockMainMenus()
//function UpdateMainMenus()
//function FixTwinsIndex()
//function RenameImageFiles()
//function angle=GetIndexFileRotation()
//function WriteBatchSpots()
//function CombineMergedFiles()


// ============  Menu definitions and callbacks ============

function SetupMainMenus()
// Remove the default Scilab popup menus and the toolbar
  RemoveMenuBar([])
// Remove any remaining LaueG pulldown menus
  delmenu("Start")
  delmenu("Orient")
  delmenu("Integrate")
  delmenu("Normalise")
  delmenu("Tools")
  delmenu("Help")
// Add popup menus
  AddConsoleMenu("Start", ...
           ["Single Image Mode","Start Batch Mode", ...
            "Abort, Unlock Menus","Exit LaueG"], ...
           "MainStart_Menu")
//
  AddConsoleMenu("Orient", ...
           ["Load Index File","Refine Orientation (initial)", ...
            "Refine Orientation (repeat)", "Change Lattice Centering", ...
            "Refine Twin Separation (initial)", ...
            "Refine Twin Separation (repeat)", ...
            "Fix Twins Index","Batch File Status"], ...
           "MainOrient_Menu")
//
  AddConsoleMenu("Integrate", ...
           ["Integrate Spots","Intensity vs d-spacing", ...
            "Display Ellipses","Integrate Spots (Twinned)", ...
            "Integrate Spots (Modulated)","Integrate Spots (OLD)", ...
            "Batch File Status"], ...
           "MainInteg_Menu")
//
  AddConsoleMenu("Normalise", ...
           ["Normalise Intensities", "Read Normalise Output", ...
            "Check Outlier Spots","Check Outliers + Equivs", ...
            "Edit Reject Spots","Delete Reject Spots", ...
            "Load Wavelength File","Plot Wavelength Files", ...
            "Batch File Status"], ...
           "MainNorm_Menu")
//
  AddConsoleMenu("Tools", ...
           ["Rename Image Files","Compare Orientations", ...
            "Convert CYCLOPS Files","Output Observed Spots", ...
            "Flicker Image Display","Merge *_mrg.hkl Files"], ...
           "MainTools_Menu")
//
// Add a "help" popup menu
  AddConsoleMenu("Help",["Help Manual","Software Versions","LaueG Updates", ...
                  "References","Support and Feedback"], ...
           "MainHelp_Menu")
//
// Turn off the batch-mode menus
  UpdateMainMenus()
//
endfunction


function MainStart_Menu(iMenu)
global BatchFiles
// Start single mode: disable batch mode and tools, display an image
  if(iMenu == 1) then
    BatchFiles=[]       // turns off batch mode
    UpdateMainMenus()
    CreateSingleImageWindow(%f)     // %f for new window, don't reuse
  end
// Start batch mode: load files, enable menus, delete image
  if(iMenu == 2) then
    if IsImageOn() then delete(ImageFigure); end
    LoadBatchFileList()
    UpdateMainMenus()
  end
// Abort a running operation and restore menus
  if(iMenu == 3) then
    AbortLaueG()
  end
// Exit LaueG & SCILAB
  if(iMenu == 4) then
    exit
  end
endfunction


function MainOrient_Menu(iMenu)
global LatticeInfo
//
// "Load Index File"
  if(iMenu == 1) then
    LoadSelectedIndexFile(%t)
    return
  end
//
// Demote any in-memory level=3
  for i=1:size(LatticeInfo,1)
    LatticeInfo(i).level=min(2,LatticeInfo(i).level)
  end
//
// Get the status of memory and batch file info
  sTable=GetBatchStatus()
  bStatus=ParseBatchStatus(sTable)
//
// "Refine Orientation (initial)"
  if(iMenu == 2) then
    DoBatchChecks(~bStatus.bIdxMem,"No orientation solution found in memory")
    BatchOrientSpots(%t)     // bMemory
  end
//
//"Refine Orientation (repeat)"
  if(iMenu == 3) then
    DoBatchChecks([~bStatus.bIdxPart,(bStatus.bIdxPart & bStatus.bIdxBad)], ...
        ["All batch files must be at least partially oriented", ...
         "Inconsistent lattices or cell volumes"])
    BatchOrientSpots(%f)
  end
//
// "Change Lattice Centering"
  if(iMenu == 4) then
    DoBatchChecks(~bStatus.bIdxPart, ...
        "All batch files must be at least partially oriented")
    BatchCenter_Popup(LatticeInfo(1).icen)
  end
//
// "Refine Twin Separation (initial)"
  if(iMenu == 5) then
    DoBatchChecks([~bStatus.bIdxFull,~bStatus.bTwinMem], ...
        ["All batch files must be fully oriented", ...
         "No twin lattice information found in memory"])
    BatchOrientTwins(%t)     // bMemory
  end
//
// "Refine Twin Separation (repeat)"
  if(iMenu == 6) then
    DoBatchChecks([~bStatus.bIdxPart,~bStatus.bTwinPart, ...
                      (bStatus.bIdxPart & bStatus.bIdxBad), ...
                      (bStatus.bTwinPart & bStatus.bTwinBad)], ...
        ["All batch files must be at least partially oriented",...
         "All batch files must have at least partially oriented twins",...
         "Inconsistent lattices or cell volumes",...
         "Inconsistent twin lattices or cell volumes"])
    BatchOrientTwins(%f)
  end
//
// "Fix Twins Index"
  if(iMenu == 7) then
    FixTwinsIndex()
  end
//
// "Batch File Status"
  if(iMenu == 8) then
    OutputBatchStatus(sTable)
    sOut=BatchStatusString(bStatus)
    InfoBox([sOut;"";"More information shown on console window"])
  end
//
endfunction


function MainInteg_Menu(iMenu)
//
// Get the status of memory and batch file info
  sTable=GetBatchStatus()
  bStatus=ParseBatchStatus(sTable)
//
// "Integrate Spots"
  if(iMenu == 1) then
    DoBatchChecks([~bStatus.bIdxFull,(bStatus.bIdxPart & bStatus.bIdxBad)], ...
        ["All batch files must be fully oriented", ...
         "Inconsistent lattices or cell volumes"])
    ArgonneBoxes_Popup(%f,%f,%f)    // bOldMode,bTwinMode,bModulateMode
  end
//
// "Intensity vs d-spacing"
  if(iMenu == 2) then
    DoBatchChecks([~bStatus.bIntAll,(bStatus.bIntAll & bStatus.bIntIntMix)], ...
        ["All batch files must be integrated", ...
         "Cannot mix integrated intensity types"])
    TestIntegDmin()
  end
//
// "Display Ellipses"
  if(iMenu == 3) then
    DoBatchChecks([~bStatus.bIntAll,(bStatus.bIntAll & bStatus.bIntIntMix)], ...
        ["All batch files must be integrated", ...
         "Cannot mix integrated intensity types"])
// Create the image window for the first batch file
    CreateIntegImageWindow(1)
  end
//
// "Integrate Spots (Twinned)"
  if(iMenu == 4) then
    DoBatchChecks([~bStatus.bIdxFull,~bStatus.bTwinFull, ...
                     (bStatus.bIdxPart & bStatus.bIdxBad), ...
                     (bStatus.bTwinFull & bStatus.bTwinBad),~bStatus.bModNone], ...
        ["All batch files must be fully oriented", ...
         "All batch files must have fully oriented twins", ...
         "Inconsistent lattices or cell volumes", ...
         "Inconsistent twin lattices or cell volumes", ...
         "Cannot mix twinned and modulated structures"])
    ArgonneBoxes_Popup(%f,%t,%f)
  end
//
// "Integrate Spots (Modulated)"
  if(iMenu == 5) then
    DoBatchChecks([~bStatus.bIdxFull,~bStatus.bModAll, ...
                       (bStatus.bIdxPart & bStatus.bIdxBad), ...
                       (bStatus.bModAll & bStatus.bModBad),~bStatus.bTwinNone], ...
        ["All batch files must be fully oriented", ...
         "All batch files must have modulations", ...
         "Inconsistent lattices or cell volumes", ...
         "Batch files have inconsistent modulations", ...
         "Cannot mix twinned and modulated structures"])
    ArgonneBoxes_Popup(%f,%f,%t)
  end
//
// "Integrate Spots (OLD)"
  if(iMenu == 6) then
    DoBatchChecks([~bStatus.bIdxFull,(bStatus.bIdxPart & bStatus.bIdxBad)], ...
        ["All batch files must be fully oriented", ...
         "Inconsistent lattices or cell volumes"])
    ArgonneBoxes_Popup(%t,%f,%f)
  end
//
// "Batch File Status"
  if(iMenu == 7) then
    OutputBatchStatus(sTable)
    sOut=BatchStatusString(bStatus)
    InfoBox([sOut;"";"More information shown on console window"])
  end
//
endfunction


function MainNorm_Menu(iMenu)
//
// "Normalise Intensities"
  if(iMenu == 1) then
    [bTwin,bMod]=CheckLaue4Status(%t)    // bRun
    BatchLaue4(bTwin,bMod)
  end
//
// "Read Normalise Output"
  if(iMenu == 2) then
    if( isempty(fileinfo("laue4.lis")) ) then
      WarnBox("Output file, laue4.lis, not found")
    else
      dos("start notepad laue4.lis")
    end
  end
//
// "Check Outlier Spots"
  if(iMenu == 3) then
    if( isempty(fileinfo("laue4_bad.lis")) ) then
      WarnBox("No intensity outlier spots from Laue4")
    else
      CreateRejectsImageWindow("laue4_bad.lis",%t)
    end
  end
//
// "Check Outliers + Equivs"
  if(iMenu == 4) then
    if( isempty(fileinfo("laue4_sus.lis")) ) then
      InfoBox("No intensity outlier spots from Laue4")
    else
      CreateRejectsImageWindow("laue4_sus.lis",%t)
    end
  end
//
// "Edit Reject Spots"
  if(iMenu == 5) then
    if( isempty(fileinfo("laue4_rej.dat")) ) then
      WarnBox("No reject spots, so far")
    else
      CreateRejectsImageWindow("laue4_rej.dat",%f)
    end
  end
//
// "Delete Reject Spots"
  if(iMenu == 6) then
    deletefile("laue4_rej.dat")
  end
//
// "Load Wavelength File"
  if(iMenu == 7) then
    WaveFileSelect_Popup()
  end
//
// "Plot Wavelength Files"
  if(iMenu == 8) then
    PlotWavelengthFiles()
  end
//
// "Batch File Status"
  if(iMenu == 9) then
    CheckLaue4Status(%f)
  end
endfunction


function MainTools_Menu(iMenu)
global BatchFiles
//
// Rename image files to base-name + numeric format
  if(iMenu == 1) then
    RenameImageFiles()
  end
//
// Output rotation angle between the UBs of two index files
  if(iMenu == 2) then
    angle=GetIndexFileRotation()
    str=msprintf("Angle between PHI axes = %.f degrees",angle)
    InfoBox(str)
  end
//
// Run the CYCLOPS converter program
  if(iMenu == 3) then
    temp=uigetdir(DataDir, "Select directory with CYCLOPS files")
    if( isdir(temp) ) then
        ChangeDataDir(temp,%t)
        RunCyclopsPrep()
    end
  end
//
// Run find spots on batch files and output spots to files
  if(iMenu == 4) then
    if IsImageOn() then
      WarnBox("Close displayed image, then restart tool")
    else
      WriteBatchSpots()
    end
  end
//
// Load list of flick files, and create flicker image display
  if(iMenu == 5) then
    BatchFiles=[]   // Turn off batch mode
    LoadFlickList()
    CreateFlickImageWindow()
  end
//
// Combine *_mrg.hkl files
  if(iMenu == 6) then
    CombineMergedFiles()
  end
//
endfunction


function MainHelp_Menu(iMenu)
//
// Help Manual
  if(iMenu == 1) then
    if( ~dos( "start """" """ +ProgDir+ "LaueG Manual.htm""" ) ) then
      ErrorBox("Unable to open the LaueG manual in a web browser")
    end
  end
//
// Software Versions
  if(iMenu == 2) then
//
// Get table of console *.exe files and version dates
    LockMenus("Collecting Version Dates")
    sTable=RunExeVersions()
    UnlockMenus()
// Add LaueGVersion and Scilab versions and display
    text=[LaueGVersion(); ...
          getversion()+" running on "+getos(); ""; ...
          msprintf("%-30.30s %s\n",sTable)]
    InfoTitleBox(text,"LaueG Software Versions")
  end
//
// LaueG Updates
  if(iMenu == 3) then
    text=["Updates available from:";""; ...
          "Koala Data Analysis page of the www.ansto.gov.au website";""; ...
          "Pressing the OK button will attempt to take you to that web page"]
    InfoTitleBox(text,"LaueG Updates")
    site="https://www.ansto.gov.au/research/user-office/"+ ...
         "instruments/neutron-scattering-instruments/koala/services"
    dos( "start """" """ +site+ """" )
  end
//
// References
  if(iMenu == 4) then
    text=["If you find LaueG useful in processing your Laue data, please", ...
          "reference the papers 1 & 2. If you wish to reference the analysis", ...
          "technique, please consider referencing papers 2 & 3.", ...
          "", "", ...
          "LaueG software for displaying and processing neutron Laue images", ...
          "Piltz, R. O. (2018). J. Appl. Cryst. 51, 963-965.", ...
          "https:/"+"/doi.org/10.1107/S1600576718005046", ...
          "", ...
          "Spot integration using ""argonne_boxes""", ...
          " C Wilkinson, H W Khamis, R F D Stansfield and", ...
          " G J McIntyre (1988), J. Appl. Cryst., 21, 471", ...
          "https:/"+"/doi.org/10.1107/S0021889888005400", ...
          "", ...
          "Accurate data processing for neutron Laue diffractometers", ...
          "J. Appl. Cryst. (2018). 51, 635-645.", ...
          "https:/"+"/doi.org/10.1107/S1600576718005058"]
    InfoTitleBox(text,"LaueG References")
  end
//
// Support and Feedback
  if(iMenu == 5) then
    InfoTitleBox("Please send any bugs or comments to rop@ansto.gov.au","LaueG Support")
  end
endfunction


// ==================== Helper routines  ================

function DoBatchChecks(bFails,sFails)
//
  sOut=""
  for i=1:size(bFails,2)
    if( bFails(i) ) then
      sOut($+1)=sFails(i)
    end
  end
  sOut(1)=[]
//
  if(sOut ~= []) then
    AbortBox([sOut;"";"See ""Batch File Status"" for more information"])
  end
//
endfunction


function LockMainMenus()
//
// Lock batch-mode menus
  ChangeMenus(get(0),"Orient","","off")
  ChangeMenus(get(0),"Integrate","","off")
  ChangeMenus(get(0),"Normalise","","off")
//
// Lock Tools menus
  ChangeMenus(get(0),"Tools","","off")
//
// Lock "change mode" items in Start menu
  ChangeMenus(get(0),"Start","Single Image Mode","off")
  ChangeMenus(get(0),"Start","Start Batch Mode","off")
//
endfunction


function UnlockMainMenus()
//
// Unlock "change mode" items in Start menu
  ChangeMenus(get(0),"Start","Single Image Mode","on")
  ChangeMenus(get(0),"Start","Start Batch Mode","on")
//
// Unlock Tools menus
  ChangeMenus(get(0),"Tools","","on")
//
// Unlock batch-mode menus if in batch mode
  UpdateMainMenus()
//
endfunction


function UpdateMainMenus()
//
    if( IsBatchMode() ) then
      ChangeMenus(get(0),"Orient","","on")
      ChangeMenus(get(0),"Integrate","","on")
      ChangeMenus(get(0),"Normalise","","on")
    else
      ChangeMenus(get(0),"Orient","","off")
      ChangeMenus(get(0),"Integrate","","off")
      ChangeMenus(get(0),"Normalise","","off")
    end
//
endfunction


function FixTwinsIndex()
global LatticeInfo OrientInfo ModulateInfo
//
// Get integrated intensity type, and check which twins are swapped
  sStatusTable=GetBatchStatus()
  if( or(sStatusTable(2:$,[10,14]) == "") ) then
    ErrorBox("All files must have oriented twins")
    return
  end
//
  bSwaps=( convstr(sStatusTable(2:$,18)) == "y" )
  iSwaps=find( double(bSwaps) ~= median(double(bSwaps)) )
  if(iSwaps == []) then
    InfoBox(["Nothing to do";"All files have consistent twin numbers"])
    return
  end
//
  InfoBox(["Swapping twin numbers for files:"; ...
               blanks(5)+BatchFiles(iSwaps)'+".idx"; ...
               "";"Affected integration files will be deleted"])
  WriteLogFile(["Swapping twins for:",BatchFiles(iSwaps)'])
//
  Lsave=LatticeInfo
  Osave=OrientInfo
  Msave=ModulateInfo
  for iFile=iSwaps
    sFile=BatchFiles(iFile)+".idx"
    if( ~LoadIndexFile(sFile,%f) ) then
      AbortBox("Unable to open "+sFile)
    end
    LatticeInfo([3,2],:)=LatticeInfo([2,3],:)
    SaveIndexFile(sFile,%f)
    mdelete(sFile+".int")
    mdelete(sFile+".ell")
  end
  LatticeInfo=Lsave
  OrientInfo=Osave
  ModulateInfo=Msave
//
//
endfunction


function RenameImageFiles()
//
// Get the file names, bail out if none selected
  get_files=uigetfile(["*.tif"],DataDir,"Select image files to rename",%t);
  if( isempty(get_files) ) then return; end
//
// Create popup to get base name, start number, and options
  RenameFiles_Popup(get_files)
//
endfunction


function angle=GetIndexFileRotation()
global LatticeInfo OrientInfo ModulateInfo
//
// Calculate angle between PHI axes from two orientations read from index files
//
// Read two index files and save UBs in UB(1:2).ub
  UB=[]
  for iFile=1:2
// Ask for, load and check UB from index file
    fname=uigetfile(["*.idx"],DataDir, "Select Index File "+string(iFile), %f)
    if( isempty(fname) ) then abort; end
// Read and check index file
    [Linfo,Oinfo,Minfo]=ReadIndexFile(fname)
    if(Linfo == []) then
      AbortBox("Index file is missing or corrupt")
    end
    if(Linfo(1).level < 2) then
      AbortBox("Index file is not oriented")
    end
// Copy the UB
    UB(iFile).ub=Linfo(1).ub
// Loop back for second index file
  end
//
// Calculate and check the rotation matrix from UB1 to UB2
  rot_mx=UB(2).ub / UB(1).ub
  if( abs( abs(det(rot_mx)) -1) > 0.10) then
    AbortBox("The two UB matrices are from different cells")
  end
  rot_mx=( rot_mx + inv(rot_mx') )/2.0
//
// For a HKL with a vector parallel to Y using UB 1, how
// far away in angle is the vector for HKL using UB 2.
  angle=acosd(max(-1,min(1, rot_mx(2,2) )))
//
endfunction


function WriteBatchSpots()
//
// Get the file names, bail out if none selected
  get_files=uigetfile(["*.tif"],DataDir,"Select image files to process",%t)
  if( isempty(get_files) ) then return; end
//
// Banner message to console
  mprintf("\nPixel X & Y, plus spot merit, written to *.xym files\n")
//
// Loop over selected files, finding spots and writing them to a file
  LockMenus("Finding Spots for Batch Files")
  for fname=get_files

// Load ImageInfo from the image file and find spots
    [pname,bname,ename]=fileparts(fname)
    ReadImageHeader(pname+bname,%t)
    RunFindSpots(15)
    nspots=size(DisplayData.found_spots,1)
//
// Write spots to *.xym file
    [fd,ierr]=mopen(bname+".xym","w")
    if(ierr ~= 0) then
      AbortBox("Unable to open output file")
    end
    mfprintf(fd,"%8.1f%8.1f%8.1f\n",DisplayData.found_spots)
    mclose(fd)
//
    mprintf("%5d spots output to file %s.xym\n",nspots,bname)
//
  end
  UnlockMenus()
//
endfunction


function CombineMergedFiles()
//
  sNames=[]
  tmp=uigetfile(["*_mrg.hkl"],DataDir, "Files from first Folder", %t)
  while (tmp(1) ~= "") then
    sNames=[sNames,tmp]
    tmp=uigetfile(["*_mrg.hkl"],DataDir, "Files from another Folder [press Cancel if none]", %t)
  end
  sNames=unique(sNames)'
//
  nFiles=size(sNames,1)
  if(nFiles == 0) then
    return
  elseif(nFiles == 1) then
    WarnBox("Need at least two files to merge together")
  end
//
// Check for modulation indices by reading first line of first file
  fd=mopen(sNames(1),'rt')
  lines=mgetl(fd,1)
  mclose(fd)
  idot=strindex(lines(1),'.')
  nMods=(idot(1)-19)/4
  idx=(12+4*nMods)
//
  files=[]
  hklms=[]
  ints=[]
  for iFile=1:nFiles
// Read complete file into buf[]
    fd=mopen(sNames(iFile),'rt')
    buf=mgetl(fd)
    mclose(fd)
// Extract hklm as strings, and counts/dcounts as doubles
// Normal hkl & twin cases:
    strs=part(buf,[1:idx])
    if(nMods == 0) then
      strs=strs+part(buf,29:32)
    end
    dbls=[ strtod(part(buf,idx+[1:12])), strtod(part(buf,2*idx+[1:16])) ]
// Append hklm & counts to arrays
    hklms=[hklms;strs]
    ints =[ints ;dbls]
// Append iFile to files[]
    files=[files; ones(buf)*iFile ]
  end
//
// Create list of unique indices (hkl + twin + modul)
  uniqs=unique(hklms)
  nUniqs=size(uniqs,1)
//
// Create tables of counts & weights for unique hklm & ifile
// Missing observations have weights=0
  counts=zeros(nUniqs,nFiles);
  weights=zeros(nUniqs,nFiles);
  for iuniq=1:nUniqs
    i=find(hklms == uniqs(iuniq));
    counts(iuniq,files(i))=ints(i,1)';
    weights(iuniq,files(i))=ints(i,2)'.^(-2);
  end
//
// Iteratively solve weighted merge starting with equal scale factors
  scales=ones(nFiles,1);
// Iterate solve merged intensites then scale factors
  for iter=1:100
    counts0=((weights.*counts) * scales) ./ (weights * (scales.^2));
    scales2=((weights.*counts)' * counts0) ./ (weights' * (counts0.^2));
    if(max(abs(scales2-scales)) < 1e-4) then break; end
// Add convergence accelerator of 1.7 (must be < 2)
    scales=scales + (scales2-scales)*1.7;
    mprintf("%3d %10.6f %10.6f %10.6f\n",iter,scales')
  end
//
// Calculate esds and add as second column in counts0 & scales
  counts0(:,2)=( weights * (scales.^2) ) .^ (-0.5)
  scales(:,2)=( (weights' * (counts0(:,1).^2)) ) .^ (-0.5)
//
// Perform reasonable factor-of-ten scaling of merged counts (as used in Laue4)
// Multiply counts0 & scales by 10 so that smallest dcounts0 > 0.3
  for iter=1:3
    if(min(counts0(:,2)) <= 0.3) then
      counts0=counts0*10.0
      scales=scales*10.0
    end
  end
// Multiply counts0 & /scales by 0.1 so that largest counts0 < 1e5
  for iter=1:3
    if(max(counts0(:,1)) >= 1e5) then
      counts0=counts0*0.1
      scales=scales*0.1
    end
  end
//
// Convert "divisor" scale factors to "multiplier" ones
  scales(:,1)=1 ./ scales(:,1)
  scales(:,2)=scales(:,2) .* scales(:,1).^2
//
// Output the scale factor used for each file
  mprintf("   Scale & esd         Input File\n")
  for i=1:nFiles
    mprintf("%8.4f(%6.4f)    %s\n",scales(i,:),sNames(i))
  end
//
// Open output file and save merged intensities
  fname=OutputFileBox("merged_mrg.hkl")
  if(fname == "") then return; end
  fd=mopen(fname,'wt')
  if(nMods == 0) then
    mfprintf(fd,"%s%8.1f%8.1f%s\n",part(uniqs,1:12),counts0,part(uniqs,13:16))
  else
    mfprintf(fd,"%s%8.1f%8.1f\n",uniqs,counts0)
  end
  mclose(fd)
//
endfunction
