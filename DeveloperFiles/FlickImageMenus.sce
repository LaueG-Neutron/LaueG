global ImageInfo ImageFigure
global ReducedImage MainImage FlickData

// ===============  Menu definitions and callbacks ===============
//function SetupFlickMenus()
//function FlickFile_Menu(iMenu,iWin)
//function FlickDisplay_Menu(iMenu,iWin)
//function FlickOrder_Menu(iMenu,iWin)
// ========================  Helper routines  ===================
//function LoadFlickFileList()
//function LoadFlickData(bStrip,bBlur)

// ===============  Menu definitions and callbacks ===============

function SetupFlickMenus()
// Manage title/toolbar/pulldowns for images window in flick image mode
//
// Add the new menus
  AddImageMenu("File", ...
                ["Load Image Files","Load Index File","Close Image"], ...
               "FlickFile_Menu")
//
  AddImageMenu("Display", ...
                ["Lighter Image","Darker Image","Brightness Levels", ...
                 "Invert Grey Levels","Strip Background","With Background", ...
                 "Smoothed Image","No Smoothing"], ...
                "FlickDisplay_Menu")
//
  AddImageMenu("Image Order", ...
                ["Sort by PHI","Sort by time","Sort by name"], ...
                "FlickOrder_Menu")
//
endfunction


function FlickFile_Menu(iMenu,iWin)
//
// Select image files and reload flicker images
  if(iMenu == 1) then
    LoadFlickList()
    CreateFlickImageWindow()
  end
//
// Read in a *.idx file
  if(iMenu == 2) then
    LoadSelectedIndexFile(%f)
  end
//
// Close the image window
  if(iMenu == 3) then
    delete(ImageFigure)
  end
endfunction


function FlickDisplay_Menu(iMenu,iWin)
global DisplaySet ReducedImage MainImage
//
// Change max & min intensities for a brighter/darker display
  if(iMenu <= 2) then
    fac=2/3
    if(iMenu == 2) then
      fac=1/fac
    end
    if(ImageFigure.color_map(1,1) > 0.5) then
      fac=1/fac
    end
// Calculate all flick pixmaps for new intensity levels
    levs=MainImage.user_data*fac
    LockMenus("Creating flicker images")
    CalcFlickPixmaps(levs)
    UnlockMenus()
// Draw image from pixmap of the selected image file
    MainImage.data=squeeze( FlickData.pixmaps(FlickData.iflick,:,:) )
// Store the brightness levels used for the pixmaps
    MainImage.user_data=levs
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
// Do menu items > 4 (which all change ReducedImage[])
  if(iMenu > 4) then
//
// Decide stripped and blurred options to use
    bStrip=DisplaySet.strip_back
    bBlur=DisplaySet.blur_main
// Change bStrip & bBlur to what we want
    if(iMenu == 5) then
      bStrip=%t
    elseif(iMenu == 6) then
      bStrip=%f
    elseif(iMenu == 7) then
      bBlur=%t
    elseif(iMenu == 8) then
      bBlur=%f
    end
// If "want" is different from "current":
    if( (bStrip ~= DisplaySet.strip_back) | ...
        (bBlur ~= DisplaySet.blur_main)         ) then
// Load reduced images and other info for flick files
// LoadFlickData() locks the menus and displays a progress bar
      LoadFlickData(bStrip,bBlur)
// Calculate intensity levels from selected file
      ReducedImage=squeeze( FlickData.images(FlickData.iflick,:,:) )
// NB: Use much lower intensities levels in flicker mode
      if( DisplaySet.strip_back ) then
        levs=CalcImageHisto(1e-5,0.02)
      else
        levs=CalcImageHisto(1e-3,0.01)
      end
// Calculate all flick pixmaps for new intensity levels
      LockMenus("Creating flicker images")
      CalcFlickPixmaps(levs)
// Draw image from pixmap of the selected image file
      MainImage.data=squeeze( FlickData.pixmaps(FlickData.iflick,:,:) )
// Store the brightness levels used for the pixmaps
      MainImage.user_data=levs
// Save new strip & blur settings
      DisplaySet.strip_back=bStrip
      DisplaySet.blur_main=bBlur
      SaveSettingsFile()
      UnlockMenus()
    end
//
  end
//
endfunction


function FlickOrder_Menu(iMenu,iWin)
global FlickData
//
  if(iMenu == 1) then
    [v,idx]=gsort(FlickData.phis,"g","i")
  elseif(iMenu == 2) then
    [v,idx]=gsort(FlickData.times,"g","i")
  else
    [v,idx]=gsort(FlickData.names,"g","i")
  end
//
  FlickData.iflick=find(idx == FlickData.iflick);
  FlickData.phis=FlickData.phis(idx);
  FlickData.times=FlickData.times(idx);
  FlickData.names=FlickData.names(idx);
  FlickData.images=FlickData.images(idx,:,:);
  FlickData.pixmaps=FlickData.pixmaps(idx,:,:);
//
endfunction


// ========================  Helper routines  ===================

function LoadFlickList()
global FlickData
//
// Ask user for the two file names, and return if an empty string
  fnames=uigetfile(["*.tif"],DataDir,"Select image files to flicker",%t)
  if( fnames == [] ) then return; end
//
// Run sanity checks on file names
  if(size(fnames,2) < 2) then
    AbortBox("You must select more than 1 file")
  elseif( or(convstr(fileext(fnames)) <> ".tif") ) then
    AbortBox("Images files must have "".tif"" extension")
  end
//
// If path changed, update path name and try to load settings file
  dname=dirname(fnames(1))
  if(DataDir ~= pathconvert(dname)) then
    ChangeDataDir(dname,%t)
    WriteIniFile()
    LoadSettingsFile()
  end
//
// Save base-names of flick files
  FlickData.names=basename(fnames)
//
endfunction


function LoadFlickData(bStrip,bBlur)
global FlickData ImageInfo
//
// Load reduced images, phi, time for all flick files
// Also load ImageInfo from first flick file
//
// Load ImageInfo from flick file 1
  ReadImageHeader(FlickData.names(1),%f)
//
// Load phi, times, and reduced images into FlickData
  FlickData.images=[]
  for i=1:size(FlickData.names,2)
   LockMenus("Reading flicker image "+string(i))
// Read file headers for phi and time values
    bname=FlickData.names(i)
    info=RunImageInfo(bname,%f)
// Load phi angles
    FlickData.phis(i)=info.phi
// Load time in "days" (not strictly, but good enough for sorting)
    tim=strtod( strsplit(info.date,[":","/"," "]) )
    tim=400*tim(6)+31*tim(5)+tim(4)+( tim(1)+tim(2)/60 )/24
    FlickData.times(i)=tim
// Complain and die if images have different sizes
    if( or(ImageInfo.numxy <> info.numxy) ) then
      AbortBox("The files are different sizes")
    end
//
// Load RawImage (but only use to make ReducedImage)
    ImageInfo.basename=bname
    if( bStrip ) then
      RunStripBack(%f)        // Don't reuse
    else
      ReadImageData()
    end
// Create ReducedImage[] from RawImage[]
    CalcReducedImage()
// Optionally, blur ReducedImage[]
    if( bBlur ) then
      BlurMainImage()
    end
// Save ReducedImage[] in FlickData
    FlickData.images(i,:,:)=ReducedImage;
  end
//
endfunction
