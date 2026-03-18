// Routines that perform pixel level graphics
global RawImage ReducedImage ImageFigure MainImage
global DisplayMode ZoomAxes ImageInfo FlickData DisplaySet

// ============= Routines to draw main image =============
//function DrawNewImage(nxy,levs)
//function RedrawImageNow(levs)
//function RedrawImage(levs)
// ============= Routines to prepare image data =============
//function CalcReducedImage()
//function pixmap=CalcPixmap(image,levs)
//function levs=CalcImageHisto(frac_lo,frac_hi)
//function BlurMainImage()
//function CalcFlickPixmaps(levs)
// ============= Routines to draw zoom image =============
//function SetupZoomImage()
//function ShowZoomImage(xip,yip, xcur,ycur)
//function RedrawZoomImage()
// ============= Routines for Colour Map =============
//function [ColorMap]=MakeColorMap()
//function InvertColorMap()


// ============= Routines to draw main image =============

function DrawNewImage(nxy,levs)
global ImageFigure MainAxes MainImage
//
// Display the main image with no borders, etc.
  MainAxes=gca()
  MainAxes.data_bounds =  [1, 1; nxy]-0.5
  MainAxes.axes_visible = ["off" "off" "off"]
  MainAxes.margins=[0,0,0,0]
  ImageFigure.axes_size = nxy
  ImageFigure.color_map = MakeColorMap()
  if( DisplaySet.invert_grey ) then
    InvertColorMap()
  end
//
// Load pixmap for main image
  if( DisplayMode.flick_mode ) then
// If flick mode, copy pixmap for selected image
    pixmap=squeeze( FlickData.pixmaps(FlickData.iflick,:,:) )
  else
// If not flick mode, calculate pixmap from ReducedImage[]
    pixmap=CalcPixmap( ReducedImage ,levs)
  end
//
// Issue command to draw the pixmap, which creates MainImage
// NB: The pixmap of the image is copied to MainImage.data[]
  Matplot(pixmap,"020")
  MainImage=gce()
//
// If not flick mode, store the brightness levels used for the pixmap
  if( ~DisplayMode.flick_mode ) then
    MainImage.user_data=levs
  end
//
// If not rejects or flick mode, setup the zoom image (but invisible)
  if( ~or([DisplayMode.rejects_mode DisplayMode.flick_mode]) ) then
    SetupZoomImage()
  end
//
endfunction


function RedrawImageNow(levs)
// Immediately redraw existing image with new high & low intensities limits
//
// Postpone drawing on ImageFigure
  scf(ImageFigure)
  drawlater
//
// Redraw the existing image
  RedrawImage(levs)
  drawnow
//
endfunction


function RedrawImage(levs)
global MainImage
// Redraw existing image with new high & low intensities limits
//
// Reload main image
  MainImage.data=CalcPixmap( ReducedImage ,levs)
//
// Store the brightness levels used for the main image
  MainImage.user_data=levs
//
// Reload the main image marks
  RedrawImageMarks()
// Draw the ellipses on the MainImage.data array
  DrawCurrentEllipses()
// Redraw the zoom image, if it is visible
  if( IsZoomOn() ) then 
    RedrawZoomImage()
  end
//
endfunction


// ============= Routines to prepare image data =============

function CalcReducedImage()
global ReducedImage
// Load ReducedImage[] from RawImage[] averaged over ImageInfo.reduce^2 pixels
//
// Sanity check
  numxy=ImageInfo.numxy
  ireduce=ImageInfo.reduce
  if( modulo(numxy,ireduce) ~= [0,0] ) then
    BugBox("ImageInfo.reduce="+string(ireduce)+ ...
            " is not a divisor of the image dimensions.")
  end
//
// Copy 4D equivalent of RawImage to tmp
  numxy=numxy/ireduce
  tmp=matrix(RawImage,ireduce,numxy(1),ireduce,numxy(2));
// 
// Average tmp over index 1
  tmp1=tmp(1,:,:,:);
  for i=2:ireduce
   tmp1=tmp1+tmp(i,:,:,:);
  end
  tmp=tmp1/ireduce;
//
// Average tmp over index 3
  tmp1=tmp(:,:,1,:);
  for i=2:ireduce
   tmp1=tmp1+tmp(:,:,i,:);
  end
  tmp=tmp1/ireduce;
//
// Copy 2D equivalent of tmp to ReducedImage
  ReducedImage=squeeze(tmp);
//
endfunction


function pixmap=CalcPixmap(image,levs)
//
// Calculate pixmaps for image data and intensity levels
// Used for main and zoom images
//
// Flip Y axis and transpose flick image
  pixmap=flipdim(image,2)'
//
// Scale values to intensity levels
  pixmap=(pixmap-levs(1)) / (levs(2)-levs(1))
//
// sqrt() and clamp values to 1-256
  pixmap=1+255*min(1,sqrt(max(0, pixmap )))
//
endfunction


function levs=CalcImageHisto(frac_lo,frac_hi)
// Calculate pixel intensities corresponding to the frac_lo
// or frac_hi fraction of intensities in ReducedImage().
//
// Sort the intensities from a 5x5 grid of pixels, avoiding edges
  grid=ReducedImage(10:5:$-10,10:5:$-10)
  ixycen=int( ImageInfo.circs(1,1:2)/ImageInfo.reduce/5 -2 )
// Remove a square around the center hole by setting values to zero
  grid(ixycen(1)-5:ixycen(1)+5,ixycen(2)-5:ixycen(2)+5)=0
  isort=gsort(matrix(grid,-1,1))
// Remove any zero intensities from the sorted list
  nzero=find(isort == 0,1)
  if(nzero > 1) then
    isort=isort(1:nzero-1)
  end
// Get frac_lo lowest intensity and frac_hi highest intensity
  ilo=floor(size(isort,1)*(1-frac_lo))
  ihi=ceil(size(isort,1)*frac_hi)
  levs(1)=isort(ilo)
  levs(2)=isort(ihi)
endfunction


function BlurMainImage()
global ReducedImage
//
// Blur ReducedImage[] using star-like convolution
// Takes ~0.09 sec on my old i5-2500
//
  cnv=[1;2;1]*[1,2,1]
  cnv=cnv/sum(cnv)
  ReducedImage=conv2(ReducedImage,cnv,"same")
//
// Make edge values same as 2 from the edge
  ReducedImage([1,2,$-1,$],:)=ReducedImage([3,3,$-2,$-2],:)
  ReducedImage(:,[1,2,$-1,$])=ReducedImage(:,[3,3,$-2,$-2])
//
endfunction


function CalcFlickPixmaps(levs)
global FlickData MainImage
//
// Calculate pixmaps for all flick images for given intensity levels
//
// Calculate from FlickData.images and store in FlickData.pixmaps
  FlickData.pixmaps=[]
  for i=1:size(FlickData.names,2)
    image=squeeze( FlickData.images(i,:,:) )
    FlickData.pixmaps(i,:,:)=CalcPixmap(image,levs)
  end
//
endfunction


// ============= Routines to draw zoom image =============

function SetupZoomImage()
global ZoomAxes ZoomImage
// Setup but don't display the zoom image
//
// The zoom image is a 200 x 200 pixel image (plus overlays)
// drawn over part of the main image. The main image is displayed
// at 1:(ImageInfo.reduce) while the zoom image is 1:1
//
// Create an invisible zoom image within the current window
  ZoomAxes=newaxes()
  ZoomAxes.margins=[0,0,0,0]
// axes_bounds[] is in units of fractions of the main image
  dx=200.0/(ImageInfo.numxy(1)/ImageInfo.reduce)
  dy=200.0/(ImageInfo.numxy(2)/ImageInfo.reduce)
  ZoomAxes.axes_bounds=[0,0,dx,dy]
  ZoomAxes.user_data=[0 0]
// Create the (invisible) zoom image and store the handle
  Matplot(zeros(201,201),"022")
  ZoomImage=gce()
////////    ZoomImage.rect=[0,0,199,199]    ????? still needed !!!!!
// Setup marks for the zoom image
  SetupZoomMarks()
// Must delay these until after marks are setup
  ZoomAxes.auto_ticks = ["off" "off" "off"]
  ZoomAxes.axes_visible = ["off" "off" "off"]
  ZoomAxes.visible="off"
endfunction


function ShowZoomImage(xip,yip, xcur,ycur)
// Show the zoom image, which is in fact a 200 x 200 pixel image
// (plus overlays) drawn over part of the main image.
// xip,yip are the image-plate coords to zoom on, while xcur,ycur
// are the mouse cursor coords on the main image where we try to draw
// the zoom image. When xcur,ycur approach the image edges we have to
// compromise to keep the zoom image within the main image limits
global ZoomAxes
//
  [numx numy]=size(RawImage)
  [numx4 numy4]=size(ReducedImage)
// Store the IP coords for the lower left corner of the zoom image
  ix0=max(1,min(numx-200, xip-100 ))
  iy0=max(1,min(numy-200, yip-100 ))
  ZoomAxes.user_data=int([ix0 iy0]);
// Position the lower-left corner of the zoom image on the main image
// Positions are in fractions of the main image
  dx=ZoomAxes.axes_bounds(3)
  dy=ZoomAxes.axes_bounds(4)
  xmin=min(1-dx,max(0, (xcur-100)/numx4 ))
  ymin=min(1-dy,max(0, (ycur-100)/numy4 ))
// Declare the section of the the main image to draw on
  ZoomAxes.axes_bounds=[xmin,ymin,dx,dy]
// Draw the zoom image plus overlays
  RedrawZoomImage()
endfunction



function RedrawZoomImage()
global ZoomImage
// Draws the "zoom image" which is a 200 x 200 pixel section of
// the main image with a unscaled image-plate image and with
// overlays corresponding to the overlays on the main image
//
// Retrieve the IP pixel location of the center of the zoom image
  ix0=ZoomAxes.user_data(1)
  iy0=ZoomAxes.user_data(2)
// Create a 200x200 piece of RawImage
  image=RawImage(ix0+(1:200),iy0+(1:200))
// Draw zoom image with the same brightness/contrast as main image
  ZoomImage.data=CalcPixmap(image,MainImage.user_data)
// Draw ellipses in the zoom image as in the main image
  DrawZoomEllipses()
// Draw marks in the zoom image
  RedrawZoomMarks()
endfunction


// ============= Routines for Colour Map =============

function [ColorMap]=MakeColorMap()
// Create a grey colour map in 1-256 plus some extra colours
  ColorMap = graycolormap(256)

  ColorMap(257,1:3)=[1 .1 .3]    // red
  ColorMap(258,1:3)=[.1 .9 .1]   // green
  ColorMap(259,1:3)=[.4 .5 1]    // blue
  ColorMap(260,1:3)=[.8 .9 0]    // yellow
  ColorMap(261,1:3)=[.8 0 .8]    // magenta
  ColorMap(262,1:3)=[.1 .8 .8]   // cyan
  ColorMap(263,1:3)=[1 .6 .3]    // orange
  ColorMap(264,1:3)=[.5 .7 1]    // sky blue
  ColorMap(265,1:3)=[.7 .3 1]    // mauve
  ColorMap(266,1:3)=[1 .4 .7]    // pink
  ColorMap(267,1:3)=[.4 .8 .2]   // olive
endfunction


function InvertColorMap()
global ImageFigure
// Invert the grey part of the colour map
  ImageFigure.color_map(1:256,1:3)=1-ImageFigure.color_map(1:256,1:3)
endfunction
