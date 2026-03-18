global RawImage MainImage ImageInfo ImageFigure
global LatticeInfo OrientInfo ModulateInfo
global DisplayMode DisplayData DataDir

// ===============  Menu definitions and callbacks ===============
//function SetupImageMenus()
//function ImageFile_Menu(iMenu,iWin)
//function ImageDisplay_Menu(iMenu,iWin)
//function ImageIndex_Menu(iMenu,iWin)
//function ImageObsSpots_Menu(iMenu,iWin)
//function ImageCalcSpots_Menu(iMenu,iWin)
//function ImageIndex2_Menu(iMenu,iWin)
//function ImageConics_Menu(iMenu,iWin)
//function ImageTwins_Menu(iMenu,iWin)
//function ImageModul_Menu(iMenu,iWin)
// ============  Helper routines for refining orientations ==========
//function SingleImageOrient()
//function SingleImageOrientTwins()
// ==============  Helper routines for observed spots ============
//function ClearObsSpots()
//function FitAllObsSpots()
//function xy2=FitSpotsConvol(xy1)
// ============  Helper routines for fitting conics ==========
//function IndexNodals()
//function Vconic=FitConicSpots(xy)
//function FitPixcenConic()
//function [f, g, ind]=FitPixcenConic_Cost(x, ind)
//function ssq=FitPixcenConic_SSQ(dxy)
// ==================  Information message box routines ===============
//function ShowImageInfo()
//function ShowIndexInfo()
//function ShowInstrumInfo()


// ===============  Menu definitions and callbacks ===============

function SetupImageMenus()
// Manage title/toolbar/pulldowns for images window in single image mode
//
// Add the new menus
  AddImageMenu("File", ...
                ["Load Image File","Load Index File","Save Index File", ...
                 "Save Displayed Image","Save Obs. Spots","Load Obs. Spots", ...
                 "Show Image Info","Show Index Info","Show Instrument Info", ...
                 "Close Image"], ...
               "ImageFile_Menu")
//
// >>> This menu list must agree with SetupIntegMenus() <<<
  AddImageMenu("Display", ...
                ["Lighter Image","Darker Image","Brightness Levels", ...
                 "Invert Grey Levels","Strip Background","With Background", ...
                 "Smoothed Image","No Smoothing"], ...
                "ImageDisplay_Menu")
//
  AddImageMenu("Index", ...
                ["Input Unit Cell","Index Spots","Refine Orientation", ...
                "Transform Lattice","Reset Orientation","Refinement Options"], ...
               "ImageIndex_Menu")
//
  AddImageMenu("Obs. Spots", ...
                ["Show Obs. Spots","Hide Obs. Spots","Spot Search", ...
                 "Prune Spots","Add/Remove Spot (left click)", ...
                 "Recenter Spot (left click)","Recenter All Spots", ...
                 "Turn Off Left Click"], ...
                "ImageObsSpots_Menu")
//
  AddImageMenu("Calc. Spots", ...
            ["Show Calc. Spots", "Hide Calc. Spots", "Find HKL Spot", ...
             "Show Unmatched Spots","Hide Unmatched Spots"], ...
            "ImageCalcSpots_Menu")
//
  AddImageMenu("Index2", ...
            ["Index Spots","Refine Orient","Prune Unmatched Spots", ...
             "Index from 3 Spots","Manual Rotation + Shift", ...
             "Input Lauegen Angles"], ...
            "ImageIndex2_Menu")
//
  AddImageMenu("Conics", ...
            ["Mark/Unmark Conic Spot (left click)","Refine Pixel Center", ...
             "Mark Nodal Spots (left click)","Index from 3 Nodals", ...
             "Delete all Conic/Nodal Spots","Turn Off Left Click"], ...
            "ImageConics_Menu")
//
  AddImageMenu("Twins", ...
            ["Index Split Spots","Refine Split Spots","Default to Twin 1", ...
             "Default to Twin 2","Twin 1 to Default","Twin 2 to Default", ...
             "Display/Hide Twin 1 Spots","Display/Hide Twin 2 Spots", ...
             "Delete Twins"], ...
            "ImageTwins_Menu")
//
  AddImageMenu("Modulated", ...
            ["Input Modulation","Change Modulation", ...
             "Show Satellites","Hide Satellites", ...
             "Cursor Info: Default HKL","Cursor Info: Vector HKL", ...
             "Cursor Info: Satellite HKL","Delete Modulations"], ...
            "ImageModul_Menu")
endfunction


function ImageFile_Menu(iMenu,iWin)
//
// Select a file and load a new image
  if(iMenu == 1) then
    CreateSingleImageWindow(%t)     // %t for possibly reuse window
  end
//
// Read in a *.idx file
  if(iMenu == 2) then
    bOK=LoadSelectedIndexFile(%t)
    if( bOK ) then
// If a valid indexed solution, display calc. spots
      if(LatticeInfo(1).level > 1) then
        DrawGenHKLs(0)
      end
    end
  end
//
// Write orientation info to a *.idx file
  if(iMenu == 3) then
    if(LatticeInfo(1).level == 0) then
      AbortBox("No information to save")
    end
// Ask user for the file name
    fname=OutputFileBox(ImageInfo.basename+".idx")
    if(fname == "") then return; end
// Save index file, and write summary to logfile
    SaveIndexFile(fname,%t)
  end
//
// Write the image to a GIF file
  if(iMenu == 4) then
// Ask user for the output file name
    fname=OutputFileBox(ImageInfo.basename+".gif")
    if(fname == "") then return; end
// Output the file
    xs2gif(ImageFigure,fname)
  end
//
// Write observed spots to a *.xym file
  if(iMenu == 5) then
// Ask user for the output file name
    fname=OutputFileBox(ImageInfo.basename+".xym")
    if(fname == "") then return; end
// Output the file
    [fd,ierr]=mopen(fname,"w")
    if(ierr ~= 0) then AbortBox("Unable to open output file"); end
    mfprintf(fd,"%8.1f%8.1f%8.1f\n",DisplayData.found_spots)
    mclose(fd)
  end
//
// Load observed spots from a *.xym file
  if(iMenu == 6) then
// Ask for the file name, and open it
    fname=uigetfile(["*.xym"],DataDir, "Select spots file", %f)
    if( isempty(fname) ) then return; end
    [fd,ierr]=mopen(fname,"r")
    if(ierr ~= 0) then return; end
// Read spots from *.xym file
    [n,x,y,m]=mfscanf(1,fd,"%f%f%f")
    mclose(fd)
    if(modulo(n,3) ~= 0) then 
      AbortBox("Spot file "+fname+" is empty or corrupt")
    end
// Load and display spots if all OK
    DisplayData.found_spots=[x,y,m]
    DrawObsSpots(DisplayData.found_spots)
  end
//
// Show info about the image
  if(iMenu == 7) then
    ShowImageInfo()
  end
//
// Show index info
  if(iMenu == 8) then
    ShowIndexInfo()
  end
//
// Show instrument info
  if(iMenu == 9) then
    ShowInstrumInfo()
  end
//
// Close the image window
  if(iMenu == 10) then
    delete(ImageFigure)
  end
endfunction


function ImageDisplay_Menu(iMenu,iWin)
global DisplaySet FileImage StripImage RawImage
// This routine is used for both single and integ images
//
// Change max & min intensities for a brighter display
  if(iMenu == 1) then
    fac=1.0/1.5
    if(ImageFigure.color_map(1,1) > 0.5) then fac=1.5; end
    RedrawImageNow( MainImage.user_data*fac )
  end
//
// Change max & min intensities for a darker display
  if(iMenu == 2) then
    fac=1.0/1.5
    if(ImageFigure.color_map(1,1) > 0.5) then fac=1.5; end
    RedrawImageNow( MainImage.user_data/fac )
  end
//
// Change max & min intensities for display
  if(iMenu == 3) then
    Brightness_Popup()
  end
//
// Invert black & white (for Garry)
  if(iMenu == 4) then
    InvertColorMap()
    DisplaySet.invert_grey=~(DisplaySet.invert_grey)
    SaveSettingsFile()
  end
//
// Load stripped images, if currently unstripped
  if(iMenu == 5) then
    if( ~DisplaySet.strip_back ) then
// Load the stripped version of RawImage[], and create ReducedImage[]
      RunStripBack()
      RawImage=StripImage
      CalcReducedImage()
// Optionally, blur the main image (i.e. ReducedImage[])
      if( DisplaySet.blur_main ) then
        BlurMainImage()
      end
// Display image with new intensity limits
      levs=CalcImageHisto(0.02,0.0005)
      RedrawImageNow(levs)
// Save settings
      DisplaySet.strip_back=%t
      SaveSettingsFile()
    end
  end
//
// Load unstripped images, if currently stripped
  if(iMenu == 6) then
    if( DisplaySet.strip_back ) then
// Load FileImage[] from file, and create RawImage[] and ReducedImage[]
      ReadImageData()
      RawImage=FileImage
      CalcReducedImage()
// Optionally, blur the main image (i.e. ReducedImage[])
      if( DisplaySet.blur_main ) then
        BlurMainImage()
      end
// Display image with new intensity limits
      levs=CalcImageHisto(0.01,0.01)
      RedrawImageNow(levs)
// Save settings
      DisplaySet.strip_back=%f
      SaveSettingsFile()
    end
  end
//
// Blur the main image, if currently unblurred
  if(iMenu == 7) then
    if( ~DisplaySet.blur_main ) then
// Blur the main image (i.e. ReducedImage[])
      BlurMainImage()
// Display blurred main image using old intensity limits
      RedrawImageNow(MainImage.user_data)
// Save settings
      DisplaySet.blur_main=%t
      SaveSettingsFile()
    end
  end
//
// Restore main image without the blur
  if(iMenu == 8) then
    if( DisplaySet.blur_main ) then
// Make ReducedImage[] from RawImage[]
      CalcReducedImage()
// Display unblurred main image using old intensity limits
      RedrawImageNow(MainImage.user_data)
// Save settings
      DisplaySet.blur_main=%f
      SaveSettingsFile()
    end
  end
endfunction


function ImageIndex_Menu(iMenu,iWin)
global LatticeInfo OrientInfo
//
// Check if index solution is required, complain and die if none
  bNoIndex=(LatticeInfo(1).level < 2)
  if( bNoIndex & or(iMenu == [3,4]) ) then
    AbortBox("Require an index solution. Run Index Spots or load an index file.")
  end
//
// Check if cell and observed spots are requirements and are missing
  bNoSpots=(DisplayData.found_spots == [])
  bNoCell=(LatticeInfo(1).level < 1)
  bNeedSpots=( bNoSpots & or(iMenu == [2,3]) )
  bNeedCell=( bNoCell & (iMenu == 2) )
//
// Find spots, if needed
  if ( bNeedSpots ) then
    LockMenus("Finding spots from image")
    ClearObsSpots()
    RunFindSpots(15)
    DrawObsSpots(DisplayData.found_spots)
    UnlockMenus()
  end
//
// Enter cell parameters
  if( (iMenu == 1) | bNeedCell ) then
// Ask user to input the cell/type/centering, then set .level to 1
    InputCell_Popup( (iMenu == 2) )
    return
  end
//
// Index spots
// NB: Check InputCell_CB() before changing this code
  if(iMenu == 2) then
    LockMenus("Indexing spots")
// Run the indexer using default options
    [ub,pixcen]=RunIndexSpots([],DisplayData.found_spots,%f)
// Indexer didn"t abort, so save the new index values
    OrientInfo.pixcen=pixcen
    LatticeInfo(1).ub=CalcYRot(ImageInfo.phi)*ub
// Set level to indexed, and remove any twins
    LatticeInfo=LatticeInfo(1)
    LatticeInfo.level=2
// Draw calculated spots and unlock menus
    DrawGenHKLs(0)
    UnlockMenus()
  end
//
// Refine spot orientation
  if(iMenu == 3) then
    SingleImageOrient()
  end
//
// Transform the lattice and set cell type and centering
  if(iMenu == 4) then
    TransLattice_Popup()
  end
//
// Reset orientation parameters
  if(iMenu == 5) then
// Set DisplaySet with the instrument defaults for this image file
    LoadSetupOrient()
    for i=1:size(LatticeInfo,1)
      LatticeInfo(i).level=min(2,LatticeInfo(i).level)
    end
    InfoBox(["The following parameters were reset:"; ""; ...
             "pixel center, pixel size & skewness,"; ...
             "crystal offsets, beam vert. angle"])
    DrawGenHKLs(0)
  end
//
// Refinement Options
  if(iMenu == 6) then
    WarnBox("Index refinement option not yet implemented")
  end
endfunction


function ImageObsSpots_Menu(iMenu,iWin)
global DisplayMode DisplayData LatticeInfo OrientInfo
//
// Check for prerequisites, complain and die if any missing
  bNoSpots=(DisplayData.found_spots == [])
  if( bNoSpots & or(iMenu == [1:2,4,6:7]) ) then
    AbortBox("Require observed spots")
  end
//
// Create the popup to manage obs. spots
  if(iMenu == 1) then
    ObsSpots_Popup()
  end
//
// Hide observed spots
  if(iMenu == 2) then
    DrawObsSpots([])
  end
//
// Find spots
  if(iMenu == 3) then
    text="Minimum separation between spots? (5 to 50 pixels)"
    [radius,trail]=strtod(x_dialog(text,"30"))
    if( (radius ~= []) & (trail == "") ) then
      radius=max(5,min(50, round(radius) ))
      LockMenus("Searching for spots")
      ClearObsSpots()
      RunFindSpots(radius)
      DrawObsSpots(DisplayData.found_spots)
      UnlockMenus()
    end
  end
//
// Prune the found spots list to a user given number
  if(iMenu == 4) then
    LockMenus("")
    nspots=size(DisplayData.found_spots,1)
    text=string(nspots)+" observed spots in list, number to keep?"
    [nkeep,trail]=strtod(x_dialog(text,"0"))
    if( trail == "" ) then
      nkeep=max(0,min(nspots, round(nkeep) ))
      DisplayData.found_spots = DisplayData.found_spots(1:nkeep,:)
    end
    DrawObsSpots(DisplayData.found_spots)
    UnlockMenus()
  end
//
// Add/remove an observed spot using left-click
  if (iMenu == 5) then
    DisplayMode.mouse_mode=1
  end
//
// Recenter spots using left-click
  if (iMenu == 6) then
    DisplayMode.mouse_mode=2
  end
//
// Recenter all spots
  if (iMenu == 7) then
    FitAllObsSpots()
    DrawObsSpots(DisplayData.found_spots)
  end
//
// Turn off left click
  if (iMenu == 8) then
    DisplayMode.mouse_mode=0
  end
//
endfunction


function ImageCalcSpots_Menu(iMenu,iWin)
global DisplayMode DisplayData LatticeInfo OrientInfo
//
// Check for prerequisites, complain and die if any missing
  bNoSpots=(DisplayData.found_spots == [])
  bNoIndex=(LatticeInfo(1).level < 2)
  if( bNoIndex ) then
    AbortBox("Require index solution")
  elseif( bNoSpots & or(iMenu == [4:5]) ) then
    AbortBox("Require observed spots")
  end
//
// Create the popup to manage calc. spots
  if(iMenu == 1) then
    CalcSpots_Popup()
  end
//
// Hide calculated spots
  if(iMenu == 2) then
    DrawCalcSpots([])
  end
//
// Show position of an individual HKL
  if(iMenu == 3) then
    FindHKLSpot_Popup()
  end
//
// Draw the unmatched obs. spots (erase obs & calc spots)
  if(iMenu == 4) then
    DrawUnmatchedSpots(%t)
  end
//
// Erase unmatched spots
  if(iMenu == 5) then
    DrawUnmatchedSpots(%f)
  end
//
endfunction


function ImageIndex2_Menu(iMenu,iWin)
global LatticeInfo
//
// Check for prerequisites
  bNoSpots=(DisplayData.found_spots == [])
  bNoCell=(LatticeInfo(1).level < 1)
  bNoIndex=(LatticeInfo(1).level < 2)
// Add any missing prerequisites to "text"
  text=[]
  if( bNoCell & or(iMenu == [1,6]) ) then
    text=["Requires unit cell or index solution"; ""]
  elseif( bNoIndex & or(iMenu == [2:5]) ) then
    text=["Requires initial index solution"; ""]
  end
  if( bNoSpots &  or(iMenu == [1:4]) ) then
    text=[text;"Requires observed spots"; "Run ""Spot Search"" in ""Obs.Spots"""; ""]
  end
// Complain and die if any missing prerequisites
  if(text ~= []) then
    AbortBox(text)
  end
//
// Advanced version of index spots
  if (iMenu == 1) then
    AdvIndex_Popup()
  end
//
// Advanced version of refine orientation
  if (iMenu == 2) then
    RefineOrient_Popup()
  end
//
// Prune any obs. spots too far for calc. spots
  if (iMenu == 3) then
    AdvPrune_Popup()
  end
//
// Reindex using 3 spots by assigning X,Y,Wav to spots
  if (iMenu == 4) then
    HKLIndex3_Popup()
  end
//
// Manual rotation of the lattice, and shift of xy_cen
  if (iMenu == 5) then
// If no UB, create one from cell and an identity U matrix
    if( LatticeInfo(1).level == 1 ) then
      LatticeInfo(1).level=2
      LatticeInfo(1).ub=Cell2Bmatrix(LatticeInfo(1).cell)
      InfoBox("Creating UB matrix from cell dimensions")
    end
// Do the manual orientation
    ManualOrient_Popup()
  end
//
// Create UB from Lauegen angles
  if (iMenu == 6) then
    ConvertLauegen_Popup()
    DrawGenHKLs(0)
  end
//
endfunction


function ImageConics_Menu(iMenu,iWin)
global DisplayMode DisplayData LatticeInfo
//
// Check for prerequisites, complain and die if any missing
  bNoSpots=isempty(DisplayData.found_spots)
  b3Conics=(size(DisplayData.conic_spots,1) > 2)
  b3Nodals=(size(DisplayData.nodal_spots,1) > 2)
  if( bNoSpots & (iMenu ~= 5) ) then
    AbortBox("Require observed spots")
  elseif( ~b3Nodals & (iMenu == 4) ) then
    AbortBox("Need 3 nodal spots to index using conics")
  elseif( ~b3Conics & (iMenu == 2) ) then
    AbortBox("Need 3 conic spots to refine pixel center")
  end
//
// Mark/Unmark Conic Spot, and remove any nodal spots
  if (iMenu == 1) then
    DisplayMode.mouse_mode=3
    DisplayData.nodal_spots=[]
  end
//
// Refine Pixel Center
  if (iMenu == 2) then
    FitPixcenConic()
  end
//
// Mark Nodal Spot, and remove any conic spots
  if (iMenu == 3) then
    DisplayMode.mouse_mode=4
    DisplayData.conic_spots=[]
  end
//
// Index from 3 nodals, then return
  if (iMenu == 4) then
    IndexNodals()
    return
  end
//
// Delete all conic/nodal spots
  if (iMenu == 5) then
    DisplayData.conic_spots=[]
    DisplayData.nodal_spots=[]
  end
//
// Turn off left click, then return
  if (iMenu == 6) then
    DisplayMode.mouse_mode=0
    return
  end
//
// Redraw marks (unused by iMenu=4,6)
  DrawConicLines()
  DrawConicSpots()
  DrawNodalSpots()
//
endfunction


function ImageTwins_Menu(iMenu,iWin)
global LatticeInfo
//
// Check for prerequisites, complain and die if any missing
  levs=[LatticeInfo(1).level,0,0]
  execstr("levs(2)=LatticeInfo(2).level","errcatch")
  execstr("levs(3)=LatticeInfo(3).level","errcatch")
  if( or(iMenu == [1,3,4]) & (levs(1) < 2) ) then
    AbortBox("Require indexed default lattice")
  elseif( (iMenu == 2) & or(levs(2:3) < 2) ) then
    AbortBox("Require indexed twin 1 & 2 lattices")
  elseif( or(iMenu == [5,7]) & (levs(2) < 2) ) then
    AbortBox("Require indexed twin 1 lattice")
  elseif( or(iMenu == [6,8]) & (levs(3) < 2) ) then
    AbortBox("Require indexed twin 2 lattice")
  end
//
// "Index" Twin Orientations from default lattice
  if (iMenu == 1) then
    OrientTwins_Popup(%t)
  end
//
// Refine Twins Orientations from twins lattices
  if (iMenu == 2) then
    OrientTwins_Popup(%f)
  end
//
// Default lattice to Twin 1
  if (iMenu == 3) then
    LatticeInfo(2)=LatticeInfo(1)
  end
//
// Default lattice to Twin 2
  if (iMenu == 4) then
    LatticeInfo(3)=LatticeInfo(1)
  end
//
// Twin 1 to default lattice
  if (iMenu == 5) then
    LatticeInfo(1)=LatticeInfo(2)
  end
//
// Twin 2 to default lattice
  if(iMenu == 6) then
    LatticeInfo(1)=LatticeInfo(3)
  end
//
// Display, or hide, Twin 1 spots
  if (iMenu == 7) then
// If twin spots are currently drawn, erase them
    if( NumDrawnSpots(1) > 0 ) then 
      DrawTwin1Spots([])
// Otherwise, generate and draw twin spots
    else
      DrawGenHKLs(1)
    end
  end
//
// Display, or hide, Twin 2 spots
  if (iMenu == 8) then
// If twin spots are currently drawn, erase them
    if( NumDrawnSpots(2) > 0 ) then 
      DrawTwin2Spots([])
// Otherwise, generate and draw twin spots
    else
      DrawGenHKLs(2)
    end
  end
//
// Delete twins
  if (iMenu == 9) then
    LatticeInfo=LatticeInfo(1)
    DrawTwin1Spots([])
    DrawTwin2Spots([])
  end
//
endfunction


function ImageModul_Menu(iMenu,iWin)
global DisplayMode ModulateInfo
//
// Check for prerequisites, complain and die if any missing
  bNoModul= ~IsModulated()
  bNoIndex=(LatticeInfo(1).level < 2)
  if( bNoIndex & or(iMenu == [1:4,6,7]) ) then
    AbortBox("Requires indexed main lattice")
  elseif( bNoModul & or(iMenu == [3,4,7]) ) then
    AbortBox("Requires modulation info")
  end
//
// Input modulation for satellites
  if (iMenu == 1) then
    InputModul()
  end
//
// Change modulations vectors
  if (iMenu == 2) then
    ChangeModul_Popup()
  end
//
// Show popup to draw satellite spots
  if (iMenu == 3) then
    ShowModul_Popup()
  end
//
// Hide satellite spots
  if (iMenu == 4) then
    DisplayMode.modul_show=%f
    DrawGenHKLs(0)
  end
//
// Switch mouse display to normal integer HKL
  if (iMenu == 5) then
    DisplayMode.cursor_info=1
  end
//
// Switch mouse display to non-integer HKL vector
  if (iMenu == 6) then
    DisplayMode.cursor_info=2
  end
//
// Switch mouse display to satellite HKL vector
  if (iMenu == 7) then
    ModulMainSpot_Popup()
  end
//
// Remove all modulation stuff and redraw display
  if (iMenu == 8) then
    ModulateInfo.vecs=[]
    ModulateInfo.mults=[]
    DisplayMode.modul_show=%f
    DisplayMode.cursor_info=1
    DrawGenHKLs(0)
  end
//
endfunction


// ============  Helper routine for refining orientation ==========

function SingleImageOrient()
global LatticeInfo OrientInfo
//
// Refine spot orientation
//
// Run orient_spots in automatic mode (mode=1)
  nfound=size(DisplayData.found_spots,1)
  if(nfound < 50) then
    str=msprintf("Using only %d observed spots",nfound)
    WarnBox(str)
  end
  LockMenus("Refining orientation")
  WriteLogFile("Refining Orientation")
  [LatticeInfo(1),OrientInfo,text]=RunOrientSpots(LatticeInfo(1),OrientInfo,1,[])
//
// Set level to fully oriented
  LatticeInfo(1).level=3
//
// Generate and draw calculated spots
  DrawGenHKLs(0)
//
// Must unlock before displaying message box
  UnlockMenus()
//
// Output summary to message box and log file
  toks=tokens(text)
  text=msprintf("%s observed spots, %s matches, "+ ...
                "rms= %s (from %s) mm",toks([1,2,4,3])')
  InfoTitleBox(text,"Final Refinement Results")
  WriteLogFile( msprintf("  %s %s %.3f %.3f %.3f %.2f %.2f %.2f", ...
                               toks([2,4])',LatticeInfo(1).cell) )
//
endfunction


// ==============  Helper routines for observed spots ============

function ClearObsSpots()
global DisplayData
// Clear all observed/marked/conic spots from arrays and display
  DisplayData.found_spots=[]
  DrawObsSpots(DisplayData.found_spots)
  DisplayData.conic_spots=[]
  DisplayData.nodal_spots=[]
  DrawConicLines()
  DrawConicSpots()
  DrawNodalSpots()
endfunction


function FitAllObsSpots()
global DisplayData
// Inform user if they are trying without any observed spots
  if( isempty(DisplayData.found_spots) ) then
    WarnBox("No observed spots to recenter")
    return
  end
// Fit the x & y for all observed spots
  xy=FitSpotsConvol( DisplayData.found_spots(:,1:2) )
  DisplayData.found_spots=[xy DisplayData.found_spots(:,3)]
endfunction


function xy2=FitSpotsConvol(xy1)
// Search for the peak position xy2 within +/-10 pixels of xy1.
// Convolutes the image with a 2D gaussian of FWHM=5 pixels,
// then uses the maximum of the convolution as the peak position.
  width=3; isize=5; isize2=3*isize
//
// Create the gaussian matrix
  dsq=[-isize:isize].^2
  dsq=dsq'*ones(dsq)
  gauss=exp(-(dsq+dsq')/width^2)
//
  nspots=size(xy1,1)
  for ispot=1:nspots
// Get image limits +/-isize2 of current xy
    ixy=round(xy1(ispot,1:2))
    xs=ixy(1)+[-isize2:isize2]
    ys=ixy(2)+[-isize2:isize2]
// Extrapolate past image boundaries by repeating values
    xs=max(1,min(ImageInfo.numxy(1), xs ))
    ys=max(1,min(ImageInfo.numxy(2), ys ))
// Convolute the local image data with the gaussian
    image=conv2( RawImage(xs,ys), gauss, "valid" )
// Get maximum of the convolution and convert to pixel x,y
    [val,pos]=max(image)
    xy2(ispot,1:2)=ixy + pos - (1+size(image))/2
  end
endfunction


// ============  Helper routines for conics and nodals ==========

function IndexNodals()
global LatticeInfo
//
// Fit nodal 1 spots to find the conic vector, then run indexer
  xy_spots=DisplayData.found_spots(DisplayData.nodal_spots([1,2]),1:2)
  Vconic=FitConicSpots(xy_spots)
  Vnodals1=RunIndexConic(Vconic,DisplayData.nodal_spots([1,2]))
// Fit nodal 2 spots to find the conic vector, then run indexer
  xy_spots=DisplayData.found_spots(DisplayData.nodal_spots([1,3]),1:2)
  Vconic=FitConicSpots(xy_spots)
  Vnodals2=RunIndexConic(Vconic,DisplayData.nodal_spots([1,3]))
// Calculate and update UB and Cell from nodal solutions
  LatticeInfo=[]
  LatticeInfo.itype=1
  LatticeInfo.level=1
  LatticeInfo.ilatt=1
  LatticeInfo.icen=0
  LatticeInfo.special_pairs=[]
  scale=norm(Vnodals1(1,:))/norm(Vnodals2(1,:))
  LatticeInfo.ub=CalcYRot(ImageInfo.phi) * [Vnodals1;scale*Vnodals2(2,:)]'/1
  LatticeInfo.cell=UB2UnitCell(LatticeInfo.ub)
//
  cell_vol=1/abs(det(LatticeInfo.ub))
  text=msprintf("Cell = %.2f %.2f %.2f  %.1f %.1f %.1f\n       "+ ...
                    "(Volume = %.1f)",LatticeInfo.cell,cell_vol)
  InfoTitleBox(text,"Final Results from Nodal Indexing")
//
endfunction


function Vconic=FitConicSpots(xy)
// Create the "LSQ" matrix
  mvec=zeros(3,3)
  for i=1:size(xy,1)
    vec=Pix2Hvec(xy(i,:))
    mvec=mvec+vec'*vec

  end
// The non-trivial solution is the eigenvector
// corresponding to the smallest eigenvalue
  [evecs,evals] =spec(mvec)
  [val,ipos]=min(abs(diag(evals)))
  Vconic=real( evecs(:,ipos)' )
endfunction


function FitPixcenConic()
//
  xy0=OrientInfo.pixcen
  [fopt,xyopt]=optim(FitPixcenConic_Cost,xy0,"ar",30,10,1,1)
endfunction


function [f, g, ind]=FitPixcenConic_Cost(x, ind)
  f=FitPixcenConic_SSQ(x')
  g=numderivative( FitPixcenConic_SSQ , x )
endfunction


function ssq=FitPixcenConic_SSQ(dxy)
global OrientInfo
  OrientInfo.pixcen=dxy'
// Create an array of pixel x,y
//  xy=DisplayData.found_spots(ispots,1:2)
  xy=DisplayData.found_spots(DisplayData.conic_spots,1:2)
// Create the list of vectors and the "LSQ" matrix
  vec=Pix2Hvec(xy)
  mvec=vec' * vec
// The non-trivial solution is the eigenvector
// corresponding to the smallest eigenvalue
  [evecs,evals] =spec(mvec)
  [val,ipos]=min(abs(diag(evals)))
  Vconic=real( evecs(:,ipos)' )
//
  vec2=vec-vec*(Vconic' * Vconic)
  xy2=Hvec2Pix(vec2)
//
  if( xy2 == [] ) then
    ssq=1e6
  else
    ssq=sum( (xy-xy2).^2 )
  end
endfunction


// ============  Helper routine to input modulations ==========

function InputModul()
global ModulateInfo
//
// Input the modulation vectors
  buff=x_dialog('Enter modulation vectors',['0 0 0';'';''])
  [vecs,ierr]=evstr('['+buff+']')
  if(buff == []) then
    return
  elseif( (ierr ~= 0) | (size(vecs,2) ~= 3) ) then
    AbortBox("Invalid values")
  elseif(or( sum(abs(vecs),'c') == 0 )) then
    AbortBox("Vectors cannot be 0 0 0")
  end
//
// Input modulation multipliers
  nvecs=size(vecs,1)
  str1=strcat(repmat(" 0",1,nvecs-1))
  str2=["1"+str1;"-1"+str1;""]
  buff=x_dialog('Enter integer multiplier(s), '+string(nvecs)+' per line',str2)
  [mults,ierr]=evstr('['+buff+']')
  if(ierr ~= 0) then
    AbortBox("Invalid values")
  elseif(norm( mults-round(mults) ) > 0)
    AbortBox("Multipliers must be integers")
  elseif(size(mults,2) ~= nvecs)
    AbortBox("Not '+string(nvecs)+' values per line")
  end
//
// Input d_min limit
  buff=x_dialog('Minimum d-spacing for main spot',"2.0")
  [d_min,ierr]=evstr(buff)
  if(d_min == []) then
    AbortBox("Invalid value")
  end
//
// Everything seems valid, so load in ModulateInfo
  ModulateInfo.vecs=vecs
  ModulateInfo.mults=mults
  ModulateInfo.d_min=d_min
  ModulateInfo.d_max=99.9
//
// Output values to console
  str=msprintf('(%g,%g,%g)\n',vecs);
  mprintf("\nModulation vector(s): %s\n",strcat(str,'  '))
  mprintf("Satellites surrounding main spots:\n")
  mod_hkls=mults * vecs;
  str=strcat(repmat( "%3d", 1,size(vecs,1) ))
  mprintf("  "+str+"   (%6.3f,%6.3f,%6.3f)\n",[mults,mod_hkls])
  mprintf('Satellites not displayed for main spots with d-spacing < %g',d_min)
//
endfunction


// ==================  Information message box routines ===============

function ShowImageInfo()
// Build title and text to display
  text     ="User: "+ImageInfo.user
  text($+1)="Sample: "+ImageInfo.sample
  text($+1)="Comments: "+ImageInfo.comment
  text($+1)="Instrument: "+ImageInfo.instrum
  text($+1)="Experiment date: "+ImageInfo.date
  text($+1)=msprintf("Exposure time: %d sec",ImageInfo.expo)
  text($+1)=msprintf("Phi angle: %.3f",ImageInfo.phi)
  text($+1)=msprintf("Intensity base counts: %i",ImageInfo.basecounts)
// Output information in a popup box
  InfoTitleBox(text,ImageInfo.basename+" information")
endfunction


function ShowIndexInfo()
//
// Setup some string arrays
  latts=["Default","Twin 1","Twin 2"]
  levels=["No unit cell","Not indexed","Indexed solution","Fully refined"]
  types=["Triclinic","Monoclinic","Orthorhombic","Tetragonal", + ...
         "Cubic","Trigonal (rhom)","Trigonal (hex)","Hexagonal"];
  cens=["P","A","B","C","I","F","R"]
  spaces="          "
//
// Loop through the lattices
  text=[]
  for ilatt=1:size(LatticeInfo,1)
// Select the orientation to use
    L=LatticeInfo(ilatt)
// Add level of orientation solution
    if(size(LatticeInfo,1) > 1) then
      text($+1)="==== "+latts(ilatt)+" Lattice ===="
    end
    text($+1)=levels(L.level+1)
// Add cell dimensions and type, if possible
    if (L.level > 0) then
      text($+1)=types(L.ilatt)+", "+cens(L.icen+1)
      text($+1)=msprintf("%.3f %.3f %.3f %.2f %.2f %.2f",L.cell)+spaces
    end
// Add UB, if possible
    if (L.level > 1) then
      text($+1)="UB:"
      text($+1)=spaces+msprintf("%.5f %.5f %.5f",L.ub(1,:))
      text($+1)=spaces+msprintf("%.5f %.5f %.5f",L.ub(2,:))
      text($+1)=spaces+msprintf("%.5f %.5f %.5f",L.ub(3,:))
    end
// Add blanks between lattices
    text($+1)=""
  end
//
  text($+1)="==== "+"Orientation Parameters ===="
  text($+1)=msprintf("Pixel center:  %.2f %.2f",OrientInfo.pixcen)
  text($+1)=msprintf("Pixel size:  %.2f %.2f",1e3*OrientInfo.pixsize)
  text($+1)=msprintf("Pixel skewness:  %.5f",OrientInfo.pixskew)
  text($+1)=msprintf("Crystal offset:  %.2f %.2f %.2f",OrientInfo.xtaloff)
  text($+1)=msprintf("Beam vert. angle:  %.3f",OrientInfo.beamvert)
// Output information in a popup box
  InfoTitleBox(text,"Index & Orientation Information")
//
endfunction


function ShowInstrumInfo()
// Build the text array
  text     ="Instrument: "+ImageInfo.instrum
  text($+1)="Detector: "+ImageInfo.detector
  text($+1)=msprintf("Image size: %d x %d",ImageInfo.numxy)
  text($+1)=msprintf("Drum radius:  %.3f",OrientInfo.drumrad)
  text($+1)=msprintf("Drum offset:  %.3f %.3f %.3f",OrientInfo.drumoff)
  text($+1)=msprintf("Central hole (X,Y,R): (%.0f,%.0f,%.0f)",ImageInfo.circs(1,:))
//
// Add the list of circular exclusion areas
  if( isempty(ImageInfo.circs(2:$,:)) ) then
    text($+1)="No circular exclusion areas"
  else
    text($+1)=""
    text($+1)="Circular exclusion areas (X,Y,R):"
    ncircs=size(ImageInfo.circs,1)
    for i=2:6:ncircs    // NB: the first circle is the central given above
      text($+1)=msprintf("  (%.0f,%.0f,%.0f) ", ...
                        ImageInfo.circs(i:min(i+5,ncircs),:))
    end
  end
//
// Add the list of rectangular exclusion areas
  if( isempty(ImageInfo.rects) ) then
    text($+1)="No rectangular exclusion areas"
  else
    text($+1)=""
    text($+1)="Rectangular exclusion areas (X1,X2,Y1,Y2):"
    nrects=size(ImageInfo.rects,1)
    for i=1:4:nrects
      text($+1)=msprintf("  (%.0f,%.0f,%.0f,%.0f) ", ...
                        ImageInfo.rects(i:min(i+3,nrects),:))
    end
  end
// Output information in a popup box
  InfoTitleBox(text,"Instrument information")
endfunction
