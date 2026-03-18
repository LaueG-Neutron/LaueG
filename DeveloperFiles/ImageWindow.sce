global ImageFigure ImageInfo LatticeInfo
global DisplaySet DataDir BatchFiles

//function CreateSingleImageWindow(bReuse)
//function CreateIntegImageWindow(iFile)
//function CreateRejectsImageWindow(file_name,bOutliers)
//function CreateFlickImageWindow()
//function ResetDisplay()


function CreateSingleImageWindow(bReuse)
global FileImage StripImage RawImage
global LatticeInfo ImageInfo ImageFigure MainAxes DisplayMode LastEvent
//
//// Get the image file name, and set/check image parameters
// Ask for the file name, and return if an empty string
  file_name=uigetfile(["*.tif"],DataDir, "Select image file", %f)
  if( isempty(file_name) ) then return; end
// If path changed, update path name and try to load settings file
  dname=dirname(file_name)
  if(DataDir ~= pathconvert(dname)) then
    ChangeDataDir(dname,%t)
    WriteIniFile()
    bDminValid=LoadSettingsFile()
  end
// Construct the base name (without path)
  base_name=basename(file_name)
// Log starting image mode and the file opened
  WriteLogFileHeader("Start Image Mode")
  WriteLogFile("File: "+base_name+".tif")
// Save the current ImageInfo
  info_last=ImageInfo
// bDminValid=%t if d_min, etc. settings are valid
// FOR SOME REASON, ImageInfo can lose the .type
// IF SO, set type to 0 and instrum to 0 (dummy instrument)
  if( ~isfield(ImageInfo, "type") )
    ImageInfo.type=0
    ImageInfo.instrum=0
  end
  bDminValid=(ImageInfo.type ~= 0)
// Load ImageInfo (and partially OrientInfo) from the image file and
// instrument setup file. Update and save the settings to a file.
  ReadImageHeader(base_name,bDminValid)
// Do not allow level=3 when using the orientation in memory
  for i=1:size(LatticeInfo,1)
    LatticeInfo(i).level=min(2,LatticeInfo(i).level)
  end
//
// Do sanity checks on reusing the existing window
  if( bReuse ) then
    if( IsImageOn() ) then
      bReuse=( DisplayMode.image_mode &                   ...
               and(ImageInfo.numxy == info_last.numxy) &  ...
              (ImageInfo.reduce == info_last.reduce)        )
// Delete incompatible existing window
      if( ~bReuse ) then delete(ImageFigure); end
    end
  end
//
// Reset display modes and set image size
  ResetDisplay()
  DisplayMode.image_mode=%t
  image_nxy=ImageInfo.numxy/ImageInfo.reduce
// Lock the existing menus and display a progress bar
  LockMenus("Loading image")
// If needed, create a new image window with pulldown menus
  if( ~bReuse ) then
    ImageFigure=CreateEmptyFigure(image_nxy)
    SetupImageMenus()
  end
// Update the title of the image window
  ImageFigure.figure_name=file_name+"     Right-click & mouse-wheel to zoom image"
// Lock any new menus, but keep the same progress bar
  LockMenus("")
//
//// Read and manipulate the image data
// Load FileImage[] and StripImage[] from image file and RunStripBack()
  ReadImageData()
  RunStripBack()
// Copy stripped or original image to RawImage[]
  if( DisplaySet.strip_back ) then
    RawImage=StripImage
  else
    RawImage=FileImage
  end
// Calculate ReducedImage[] from RawImage[]
  CalcReducedImage()
// Optionally, blur the main image (i.e. ReducedImage[])
  if( DisplaySet.blur_main ) then
    BlurMainImage()
  end
// Calculate intensity limits from the intensity histogram
  if( DisplaySet.strip_back ) then
    levs=CalcImageHisto(0.02,0.0005)
  else
    levs=CalcImageHisto(0.01,0.01)
  end
//
//// Draw the main image and setup marks
// Reuse or draw a new image
  if( bReuse ) then
    RedrawImage(levs)
  else
    DrawNewImage(image_nxy,levs)
    SetupImageMarks()
// Wait for image to draw (which sizes the window) then turn off resize
    sleep(50)
    ImageFigure.resize="off"
  end
//
//// Final bit of housekeeping
// Make sure data_bounds is set to the centre of pixels
  MainAxes.data_bounds=[1,1; ImageFigure.axes_size] -0.5
// Setup the event handler
  LastEvent.time=0
  LastEvent.xy=[0 0]
  ImageFigure.event_handler="ImageEventHandler"
  ImageFigure.event_handler_enable ="on"
// Unlock the menus and remove the progress bar
  UnlockMenus()
endfunction


function CreateIntegImageWindow(iFile)
global FileImage StripImage RawImage
global ImageFigure DisplayMode MainAxes LastEvent
global LatticeInfo OrientInfo ModulateInfo
//
//// Get the image file name, and set/check image parameters
// Get the file name from the list
  if( (iFile < 1) | (iFile > size(BatchFiles,1))) then
    BugBox("Invalid iFile")
  end
// Save the batch file number
  DisplayMode.ibatch=iFile
// Get the *.tif file name and its base
  base_name=BatchFiles(iFile)
  file_name=base_name+".tif"
// Load ImageInfo (and partially OrientInfo) from the image file,
// load the reduced image into ReducedImage
  ReadImageHeader(base_name,%t)
// Load orientation & lattice parameters from index file
  [LatticeInfo,OrientInfo,ModulateInfo]=ReadIndexFile(base_name+".idx")
  if(LatticeInfo == []) then
    AbortBox("Index file "+base_name+".idx is missing")
  end
// Complain and die if not indexed
  if(LatticeInfo(1).level < 2) then
    AbortBox("Image "+file_name+" is not indexed")
  end
//
//// Prepare a new image window (don't try to re-use the old one)
  if( IsImageOn() ) then
    delete(ImageFigure)
  end
// Reset display modes and set image size
  ResetDisplay()
  DisplayMode.integ_mode=%t
  image_nxy=ImageInfo.numxy/ImageInfo.reduce
// Lock the existing menus and display a progress bar
  LockMenus("Drawing image")
// Create a new image window with pulldown menus
  ImageFigure=CreateEmptyFigure(image_nxy)
  SetupIntegMenus()
// Update the title of the image window
  ImageFigure.figure_name=file_name+"     Right-click & mouse-wheel to zoom image"
// Lock any new menus, but keep the same progress bar
  LockMenus("")
//
//// Read and manipulate the image data
// Clear saved background stripped image
  StripImage=[]
// Load RawImage[] directly from the image file or use RunStripBack()
  if( DisplaySet.strip_back ) then
    RunStripBack(%f)
    RawImage=StripImage
  else
    ReadImageData()
    RawImage=FileImage
  end
// Calculate ReducedImage[] from RawImage[]
  CalcReducedImage()
// Optionally, blur the main image (i.e. ReducedImage[])
  if( DisplaySet.blur_main ) then
    BlurMainImage()
  end
// Calculate intensity limits from the intensity histogram
  if( DisplaySet.strip_back ) then
    levs=CalcImageHisto(0.02,0.0005)
  else
    levs=CalcImageHisto(0.01,0.01)
  end
//
//// Draw the main image and setup marks
  DrawNewImage(image_nxy,levs)
  SetupImageMarks()
//
// Draw (delayed) ellipse pixels on the MainImage.data() array
  scf(ImageFigure)
  drawlater
  DrawIntegEllipses()
  drawnow
//
// Wait for image to draw (which sizes the window) then turn off resize
  sleep(50)
  ImageFigure.resize="off"
//
//// Final bit of housekeeping
// Make sure data_bounds is set to the centre of pixels
  MainAxes.data_bounds=[1,1; ImageFigure.axes_size] -0.5
// Setup the event handler
  LastEvent.time=0
  LastEvent.xy=[0 0]
  ImageFigure.event_handler="ImageEventHandler"
  ImageFigure.event_handler_enable ="on"
// Unlock the menus and remove the progress bar
  UnlockMenus()
endfunction


function CreateRejectsImageWindow(file_name,bOutliers)
global ImageFigure DisplayMode RejectsInfo MainAxes LastEvent
//
//// Get the reject spots list
// Lock the existing menus and display a progress bar
  LockMenus("Drawing spots")
// Load RejectsInfo from reject file & *.ell files (with warnings)
  RejectsInfo=RunRejectsData(file_name,%t)
// Complain and die if no reject files
  if(RejectsInfo.itwin == []) then
    AbortLaueG()
  end
//
//// Prepare the new image window and remove any existing one
// NB: the rejects window is "reused" via UpdateRejectsImage()
  if IsImageOn() then delete(ImageFigure); end
// Reset display modes and set image size
  ResetDisplay()
  DisplayMode.rejects_mode=%t
  image_nxy=[1000,500]
// Create a image window with pulldown menus
  ImageFigure=CreateEmptyFigure(image_nxy)
  SetupRejectsMenus(bOutliers)
// Update the title of the image window
  nspots=size(RejectsInfo.base_name,1)
  ImageFigure.figure_name=DataDir+file_name+"  ("+string(nspots)+ ...
                 " spots)     Left-click to mark a spot"
// Lock any new menus, but keep the same progress bar
  LockMenus("")
//
//// Read and manipulate the image data
// Load the images starting from spot #1 
  LoadRejectImagesData(1)
// Calculate intensity levels from histogram of ReducedImage[]
  levs=CalcImageHisto(0.001,0.01)
//
//// Draw the intensities image and setup marks
  DrawNewImage(image_nxy,levs)
  SetupImageMarks()
// Draw (delayed) ellipse pixels on the MainImage.data() array
  scf(ImageFigure)
  drawlater
  DrawRejectEllipses()
  drawnow
//
// Wait for image to draw (which sizes the window) then turn off resize
  sleep(50)
  ImageFigure.resize="off"
//
//// Final bit of housekeeping
// Make sure data_bounds is set to the centre of pixels
  MainAxes.data_bounds=[1,1; ImageFigure.axes_size] -0.5
// Setup the event handler
  LastEvent.time=0
  LastEvent.xy=[0 0]
  ImageFigure.event_handler="RejectsEventHandler"
  ImageFigure.event_handler_enable ="on"
// Unlock the menus and remove the progress bar
  UnlockMenus()
endfunction


function CreateFlickImageWindow()
global FileImage StripImage RawImage
global ImageFigure DisplayMode FlickData ImageInfo ReducedImage MainImage MainAxes LastEvent
//
// Load ImageInfo from flick file 1
  ReadImageHeader(FlickData.names(1),%f)
// Check if we can reuse an existing window
  bReuse=%f
  if( IsImageOn() ) then
    bReuse=DisplayMode.flick_mode
// If an incompatible image window, delete it
    if( ~bReuse ) then delete(ImageFigure); end
  end
//
// Reset display modes and set image size
  ResetDisplay()
  DisplayMode.flick_mode=%t
  image_nxy=ImageInfo.numxy/ImageInfo.reduce
//
// Load phi, times, and reduced images into FlickData
  for i=1:size(FlickData.names,2)
    base_name=FlickData.names(i)
    LockMenus("Reading flicker image "+base_name)
// Read file headers for phi and time values
    info=RunImageInfo(base_name,%f)
// Complain and die if images have different sizes
    if( or(ImageInfo.numxy <> info.numxy) ) then
      AbortBox("The files have different sizes")
    end
// Save phi angle
    FlickData.phis(i)=info.phi
// Save time in "days" (not strictly, but good enough for sorting)
    tim=strtod( strsplit(info.date,[":","/"," "]) )
    FlickData.times(i)=400*tim(6)+31*tim(5)+tim(4)+( tim(1)+tim(2)/60 )/24
// Load RawImage[]
    ImageInfo.basename=base_name
    if( DisplaySet.strip_back ) then
      RunStripBack(%f)
      RawImage=StripImage
    else
      ReadImageData()
      RawImage=FileImage
    end
// Create ReducedImage[] from RawImage[]
    CalcReducedImage()
// Optionally, blur ReducedImage[]
    if( DisplaySet.blur_main ) then BlurMainImage(); end
// Create ReducedImage[] and save to FlickData
    FlickData.images(i,:,:)=ReducedImage;
  end
//
// Set to display flick file 1
  FlickData.iflick=1
//
// Calculate intensity levels from first flick file
  ReducedImage=squeeze( FlickData.images(1,:,:) );
// NB: Use much lower intensity levels in flicker mode
  if( DisplaySet.strip_back ) then
    levs=CalcImageHisto(1e-4,0.02)
  else
    levs=CalcImageHisto(1e-3,0.01)
  end
//
// Calculate all pixmaps and display the first one
  if( bReuse ) then
// Calculate all pixmaps and copy first one to MainImage
    LockMenus("Creating flicker images")
    CalcFlickPixmaps(levs)
    MainImage.data=squeeze( FlickData.pixmaps(1,:,:) )
  else
// Create window with pulldown menus
    ImageFigure=CreateEmptyFigure(image_nxy)
    SetupFlickMenus()
// Lock newly created menus
    LockMenus("Creating flicker images")
// Calc pixmaps and draw the first one
    CalcFlickPixmaps(levs)
    DrawNewImage(image_nxy,levs)  // levs[] ignored in flicker mode
// Delay for image to draw and resize the window, then turn off resize
    sleep(50)
    ImageFigure.resize="off"
  end
// Store the brightness levels used for the pixmaps
  MainImage.user_data=levs
//
//// Final bit of housekeeping
// Make sure data_bounds is set to the centre of pixels
  MainAxes.data_bounds=[1,1; ImageFigure.axes_size] -0.5
// Setup the event handler
  LastEvent.time=0
  LastEvent.xy=[0 0]
  ImageFigure.event_handler="FlickEventHandler"
  ImageFigure.event_handler_enable ="on"
// Unlock the menus and remove the progress bar
  UnlockMenus()
// Update the title of the image window
  ImageFigure.figure_name=DataDir+FlickData.names(1)+ ...
             "     Left/Right-click to move Down/Up image list"
endfunction


function ResetDisplay()
global DisplayMode DisplayData FlickData MainMarks ZoomAxes MainAxes
// Set all display modes to off/default
  DisplayMode.image_mode=%f
  DisplayMode.integ_mode=%f
  DisplayMode.rejects_mode=%f
  DisplayMode.flick_mode=%f
  DisplayMode.ell_integ=%f
  DisplayMode.ell_model=%f
  DisplayMode.ell_rejects=%f
  DisplayMode.modul_show=%f
  DisplayMode.mouse_mode=0
  DisplayMode.cursor_info=1
// Reset all display data
  DisplayData.found_spots=[]
  DisplayData.conic_spots=[]
  DisplayData.nodal_spots=[]
  DisplayData.ell_data=[]
  DisplayData.ell_amults=[]
  DisplayData.modul_main_hkl=[0 0 0]
// Reset the flick image data
  FlickData.iflick=0
  FlickData.phis=[]
  FlickData.times=[]
  FlickData.images=[]
  FlickData.pixmaps=[]
// If no valid image, we have finished
  if( ~IsImageOn() ) then return; end
// Erase any marks drawn on the main image
  for i=1:size(MainMarks,1)
    if( is_handle_valid(MainMarks(i)) ) then
      MainMarks(i).data=[]
      MainMarks(i).user_data=[]
    end
  end
// Remove mouse-wheel zoom from main image
  if( is_handle_valid(MainAxes) ) then
    MainAxes.zoom_box=[]
  end
// Remove mouse-wheel zoom from zoom image, and turn off zoom image
  if( IsZoomOn() )
    ZoomAxes.zoom_box=[]
    ZoomAxes.visible="off"
  end
//
endfunction
