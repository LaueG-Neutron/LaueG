global RawImage ZoomAxes ImageInfo LatticeInfo DisplayMode ImageFigure BatchFiles
global IntegRectSave LastEvent

// ============= Event Handlers for image/integ and rejects modes ===========
//function ImageEventHandler(win, x,y, ibut)
//function FlickEventHandler(win, x,y, ibut)
//function RejectsEventHandler(win, x,y, ibut)
// =========== Routines to handle mouse events ========
//function MouseMove(x,y,bZoom)
//function LeftMouseClick(x,y,bZoom)
//function LeftMousePress(x,y)
//function LeftMouseRelease(x,y)
//function RightMousePress(x2,y2, x,y)
// ==========  Specific mouse/spot interaction routines =========
//function AddRemoveObsSpot(x,y,bZoom)
//function FitObsSpot(x,y)
//function MarkUnmarkConicSpot(x,y)
//function MarkUnmarkNodalSpot(x,y)
//function [ispot,dist]=GetClosestSpot(x,y)
// =========== Routines to calc x,y & HKL from mouse position ========
//function [x1,y1, x2,y2]=CalcIPCoords(x,y)
//function sTemp=MouseHKLMessage(x,y,bZoom)
//function sTemp=MouseHKLString(x,y,radius)
//function [hkl_sat,wav,dxy]=CalcBestSatellite(xy)
//function hkl=CalcHKLGuess(xy,radius)


// ============= Event Handlers for image/integ and rejects modes ===========

function ImageEventHandler(win, x,y, ibut)
global MainAxes LastEvent
//
// Event handler for the main images window which contains a scaled
// main image of the image-plate image, and a 200x200 zoom image which
// appears/disappears on right-click. Both the main and zoom images
// have drawn overlays (marks and ellipses) corresponding to the same
// image-plate coordinates. The position of the mouse is converted to
// the image-plate coordinates depending on the mouse being within the
// main or zoom image.
//
// Ignore the "close window" event
   if( ibut == -1000 ) then return; end
// Do nothing if handle is invalid (happens occasionally)
   if( ~is_handle_valid(ImageFigure) ) then return; end
//
// If a simple mouse move, handle event pileup problem
  if(ibut == -1) then
// Get number of msec since last processed event
    dat=getdate()
    dtim=dat($)+dat($-1)*1000 - LastEvent.time
    if(dtim < 0) then
      dtim=dtim+60000
    end
// Only process "fast move" events every 100 msec
    if(dtim < 100) then
      if(norm( [x y] - LastEvent.xy ) >= 0.03*dtim) then
        return
      end
    end
  end
//
// Calculate IP coords x1,y1 from the cursor position, also
// x2,y2 which is the value if the zoom image was invisible.
// On failure, CalcIPCoords returns x1=[]
  [x1,y1, x2,y2]=CalcIPCoords(x,y)
  if(x1 == []) then return; end
//
// Do all mouse movement events (also includes button drags)
  bZoom=or([ x1 ~= x2 , y1 ~= y2 ])
  MouseMove(x1,y1,bZoom)
//
// Handle events for left mouse press, release or click
  if( ibut == 0 ) then
    LeftMousePress(x1,y1)
  elseif( ibut == -5 ) then
    LeftMouseRelease(x1,y1)
  elseif( or(ibut == 3) ) then
    LeftMouseClick(x1,y1,bZoom)
  end
//
// Toggle the zoom image on right mouse press/click
  if( or(ibut == [2,5]) ) then
    RightMousePress(x2,y2, x,y)
  end
//
// Record current time and mouse x,y
  dat=getdate()
  LastEvent.time=dat($)+dat($-1)*1000
  LastEvent.xy=[x y]
//
endfunction


function FlickEventHandler(win, x,y, ibut)
global ImageFigure FlickData ReducedImage MainImage ImageInfo LastEvent
//
// Event handler for the flick images window which contains a scaled
// image of one of two "flick" images files. The left/right click are
// used to switch the displayed images.
//
//
// Ignore the "close window" event
   if( ibut == -1000 ) then return; end
// Do nothing if handle is invalid (happens occasionally)
   if( ~is_handle_valid(ImageFigure) ) then return; end
//
// If a simple mouse move, handle event pileup problem
  if(ibut == -1) then
// Get number of msec since last processed event
    dat=getdate()
    dtim=dat($)+dat($-1)*1000 - LastEvent.time
    if(dtim < 0) then
      dtim=dtim+60000
    end
// Only process "fast move" events every 100 msec
    if(dtim < 100) then
      if(norm( [x y] - LastEvent.xy ) >= 0.03*dtim) then
        return
      end
    end
  end
//
// For left,right click set imove=-1,+1
  imove=0
  if( or(ibut == [0,3,10]) ) then
    imove=-1
  elseif( or(ibut == [2,5,12]) ) then
    imove=+1
  end
//
// Change displayed image depending on imove
  if( imove ~= 0) then
    i=FlickData.iflick+imove
    if( (i<1) | (i>size(FlickData.names,2)) ) then
      beep()
    else
      FlickData.iflick=i
      ReducedImage=squeeze( FlickData.images(i,:,:) )
      MainImage.data=squeeze( FlickData.pixmaps(i,:,:) )
      ImageInfo.phi=FlickData.phis(i)
// Change figure title to the displayed image
      ImageFigure.figure_name=DataDir+FlickData.names(i)+ ...
             "     Left/Right-click to move Down/Up image list"
    end
  end
//
// Calculate IP coords x1,y1 from the cursor position
// x2,y2 (the values for the zoom image) are ignored
  [x1,y1, x2,y2]=CalcIPCoords(x,y)
//
// Do mouse movement events
  MouseMove(x1,y1,%f)
//
// Record current time and mouse x,y
  dat=getdate()
  LastEvent.time=dat($)+dat($-1)*1000
  LastEvent.xy=[x y]
//
endfunction


function RejectsEventHandler(win, x,y, ibut)
global ImageFigure RejectsInfo LastEvent
// Special handler for image events when in "reject spots" mode
//
// Ignore the window closed event
   if( ibut == -1000 ) then return; end
// Do nothing if handle is invalid (happens occasionally)
   if( ~is_handle_valid(ImageFigure) ) then return; end
//
// If a simple mouse move, handle event pileup problem
  if(ibut == -1) then
// Get number of msec since last processed event
    dat=getdate()
    dtim=dat($)+dat($-1)*1000 - LastEvent.time
    if(dtim < 0) then
      dtim=dtim+60000
    end
// Only process "fast move" events every 100 msec
    if(dtim < 100) then
      if(norm( [x y] - LastEvent.xy ) >= 0.03*dtim) then
        return
      end
    end
  end
//

// Should display the hkl and file_name for the reject spot
  ix=floor(x/100)
  iy=floor(y/100)
  ispot=RejectsInfo.ifirst+ix+10*iy
  nspots=size(RejectsInfo.itwin,1)
  if(ispot <= nspots) then
    if(RejectsInfo.itwin(1) == 0) then
      ImageFigure.info_message=msprintf( "%s  (%4i%4i%4i)  %s", ...
          RejectsInfo.base_name(ispot),RejectsInfo.hkl(ispot,1:3) , ...
          RejectsInfo.comment(ispot)       )
    else
      ImageFigure.info_message=msprintf(  "%s  (%4i%4i%4i%4i)  %s", ...
          RejectsInfo.base_name(ispot),RejectsInfo.hkl(ispot,1:3) , ...
          RejectsInfo.itwin(ispot),RejectsInfo.comment(ispot)    )
    end
// Handle events for left mouse press, click or double-click
// Should "mark" the reject spot and draw a red cross in it
    if( (ibut == 0) | (ibut == 3) | (ibut == 10) ) then
      RejectsInfo.marked(ispot)=~RejectsInfo.marked(ispot)
      scf(ImageFigure)
      drawlater
      DrawRejectSpots()
      drawnow
    end
  else
    ImageFigure.info_message=""
  end
//
// Record current time and mouse x,y
  dat=getdate()
  LastEvent.time=dat($)+dat($-1)*1000
  LastEvent.xy=[x y]
//
endfunction


// =========== Routines to handle mouse events ========

function MouseMove(x,y,bZoom)
global ImageFigure IntegRectSave
//
// IF doing a "rubber" rectangle for integ-mode (NYI)
  if( DisplayMode.mouse_mode > 11 ) then
    IntegRectSave(3:4)=([x,y])
    DrawRectNow(IntegRectSave)
    ImageFigure.info_message==msprintf("xy=(%4i,%4i)",round(x),round(y))
    return
  end
//
// ELSE, write information message on image window
//
// Get intensity from ReducedImage or RawImage
try
  if( DisplayMode.flick_mode ) then
    ix=1 + round( (x-1)/ImageInfo.reduce )
    iy=1 + round( (y-1)/ImageInfo.reduce )
    inten=ReducedImage(ix,iy)
  else
    inten=RawImage(int(x),int(y))
  end
// Write x,y and intensity to string "sTemp"
  sTemp=msprintf( "xy=(%4i,%4i)  int=%5i  ",x,y,inten)
catch
// Write x,y to string "sTemp"
  sTemp=msprintf( "xy=(%4i,%4i)  ",x,y)
end
//
// Append to sTemp information on none/single/twin hkls
// KLUDGE: break up command into 2 lines due to TRY/CATCH bug
  sTemp2=MouseHKLMessage(x,y,bZoom)
  sTemp=sTemp+sTemp2
//
// Update information message on image window
  ImageFigure.info_message=sTemp
//
endfunction


function LeftMouseClick(x,y,bZoom)
global DisplayMode
// Handle events for left mouse click depending on mouse_mode
//
// mouse_mode=0     Nothing to do, so return
  if    ( DisplayMode.mouse_mode ==  0 ) then
    return
// mouse_mode=1     Add/remove an observed spot
  elseif( DisplayMode.mouse_mode ==  1 ) then
    AddRemoveObsSpot(x,y,bZoom)
// mouse_mode=2     Fit x,y for an observed spot
  elseif( DisplayMode.mouse_mode ==  2 ) then
    FitObsSpot(x,y)
// mouse_mode=3    Add a conic spot
  elseif( DisplayMode.mouse_mode == 3 ) then
    MarkUnmarkConicSpot(x,y)
// mouse_mode=4    Add a nodal spot
  elseif( DisplayMode.mouse_mode == 4 ) then
    MarkUnmarkNodalSpot(x,y)
// mouse_mode>30    Reset for next rubber rectangle
  elseif( DisplayMode.mouse_mode > 10 ) then
    DisplayMode.mouse_mode = 11
  end
//
// Update the observed, marked or conic spot display
  if( DisplayMode.mouse_mode <=  2 ) then
    DrawObsSpots(DisplayData.found_spots(:,1:2))
  elseif( DisplayMode.mouse_mode <= 4 ) then
    DrawConicLines()
    DrawConicSpots()
    DrawNodalSpots()
  end
endfunction


function LeftMousePress(x,y)
global DisplayMode IntegRectSave
//
// Set/reset to draw rubber rectangles
  if( DisplayMode.mouse_mode > 10 ) then
    DisplayMode.mouse_mode = 12
    IntegRectSave=[x,y,x,y]
    DrawRectNow(IntegRectSave)
  end
//
endfunction


function LeftMouseRelease(x,y)
global DisplayMode IntegRectSave
//
// Reset mouse mode for next rubber rectangle
  if( DisplayMode.mouse_mode == 12 ) then
    DisplayMode.mouse_mode = 11
    IntegRectSave(3:4)=([x,y])
    DrawRectNow(IntegRectSave)
    irect=round(IntegRectSave([1,3,2,4]))
    text=msprintf("spots with X = %i - %i, Y = %i - %i",irect)
    ibut=messagebox(msprintf("Reject %s ?",text), "Reject spots", ...
            "question",["This image only","All images","Cancel"],"modal")
    if(ibut == 1) then
      mprintf("Rejecting %s for current image\n",text)
    elseif(ibut == 2) then
      mprintf("Rejecting %s for all %i images\n", text,size(BatchFiles,1))
    end
//
    DrawRectNow([])
  end
//
endfunction


function RightMousePress(x2,y2, x,y)
global ZoomAxes
//
// If zoom image is off, draw it (and marks) over the main image
  if(ZoomAxes.visible == "off") then
    scf(ImageFigure)
    drawlater
    ShowZoomImage(x2,y2, x,y)
    ZoomAxes.visible="on"
    drawnow
  else
// Turn off the zoom image
    ZoomAxes.visible="off"
  end
//
endfunction


// ==========  Specific mouse/spot interactions =========

function AddRemoveObsSpot(x,y,bZoom)
global DisplayData
//
// Find index and distance to nearest observed spot
  [ispot,dist]=GetClosestSpot(x,y)
  if(ispot == 0) then
    dist=11
  end
//
// Double distance value (halve cutoff) if in zoom mode
  if( bZoom ) then
    dist=dist*2
  end
//
// If distance > 10 pixels add a new spot, else remove it
  if(dist > 10) then
//
// Add spot to start of found_spots and give it a large merit value
    DisplayData.found_spots=[x,y,1000.0; DisplayData.found_spots]
//
// Increase by 1 any conic/nodal indices
    if(DisplayData.conic_spots ~= []) then
      DisplayData.conic_spots=DisplayData.conic_spots+1
    end
    if(DisplayData.nodal_spots ~= []) then
      DisplayData.nodal_spots=DisplayData.nodal_spots+1
    end
//
  else
//
// Remove ispot from found_spots()
    DisplayData.found_spots=[ DisplayData.found_spots([1:ispot-1],:) ; ...
                              DisplayData.found_spots([ispot+1:$],:) ]
//
// Remove ispot from marked/conic/nodal lists
    ilist=find(DisplayData.conic_spots ~= ispot)
    DisplayData.conic_spots=DisplayData.conic_spots(ilist)
    ilist=find(DisplayData.nodal_spots ~= ispot)
    DisplayData.nodal_spots=DisplayData.nodal_spots(ilist)
//
// Decrease by 1 any conic/nodal indices after ispot
    ilist=find(DisplayData.conic_spots > ispot)
    if(ilist ~= []) then
      DisplayData.conic_spots(ilist)=DisplayData.conic_spots(ilist)-1
    end
    ilist=find(DisplayData.nodal_spots > ispot)
    if(ilist ~= []) then
      DisplayData.nodal_spots(ilist)=DisplayData.nodal_spots(ilist)-1
    end
//
// Redraw any conic or nodal lines
    DrawConicLines()
//
  end
//
endfunction


function FitObsSpot(x,y)
global DisplayData
//
// Find index of nearest observed spot, return if none
  ispot=GetClosestSpot(x,y)
  if(ispot == 0) then
    WarnBox("No observed spots to fit")
    return
  end
//
// Change the x & y values for the spot
  DisplayData.found_spots(ispot,1:2)=FitSpotsConvol([x,y])
//
endfunction


function MarkUnmarkConicSpot(x,y)
global DisplayData
//
// Find index of nearest observed spot, return if none
  ispot=GetClosestSpot(x,y)
  if(ispot == 0) then
    WarnBox("No observed spots to mark")
    return
  end
//
// If the spot is not in the list add it, else remove it
  ifind=find(DisplayData.conic_spots == ispot)
  if(ifind == []) then
    DisplayData.conic_spots($+1)=ispot
  else
    ifind=find(DisplayData.conic_spots ~= ispot)
    DisplayData.conic_spots=DisplayData.conic_spots(ifind)
  end
//
endfunction


function MarkUnmarkNodalSpot(x,y)
global DisplayData
//
// Find index of nearest observed spot, return if none
  ispot=GetClosestSpot(x,y)
  if(ispot == 0) then
    WarnBox("No observed spots to mark")
    return
  end
//
// If the spot is not in the list add it, else remove it
  ifind=find(DisplayData.nodal_spots == ispot)
  if(ifind == []) then
    DisplayData.nodal_spots($+1)=ispot
  else
    ifind=find(DisplayData.nodal_spots ~= ispot)
    DisplayData.nodal_spots=DisplayData.nodal_spots(ifind)
  end
//
// Remove spots if more than 3 in the list (FIFO behaviour)
  if( size(DisplayData.nodal_spots,1) > 3) then
    DisplayData.nodal_spots=DisplayData.nodal_spots($-2:$)
  end
// 
endfunction


function [ispot,dist]=GetClosestSpot(x,y)
//
// If no observed spots, return ispot=0 & dist=1e10
  if( isempty(DisplayData.found_spots) ) then
    ispot=0
    dist=11
    return
  end

// Get index and distance to closest observed spot
  dsq=(x-DisplayData.found_spots(:,1)).^2 + (y-DisplayData.found_spots(:,2)).^2
  [dist,ispot]=min(sqrt(dsq))
//
endfunction


// =========== Routines to calc x,y & HKL from mouse position ========

function [x1,y1, x2,y2]=CalcIPCoords(x,y)
//
// Return IP pixel x1,y1 from the mouse x,y coords including if
// the mouse is within the zoom image. Also returns x2,y2 which
// is calculated ignoring any zoom image.
// Return x1=[] to signal a failure
//
try
//
//
  numx=ImageInfo.numxy(1)
  numy=ImageInfo.numxy(2)
  numxc=numx/ImageInfo.reduce
  numyc=numy/ImageInfo.reduce
//
// Shift X,Y and flip Y so that (1,1) is at the lower-left corner
  x2=x+1
  y2=numyc-y
//
// Correct for mouse-wheel zooming of main image
  zb=MainAxes.zoom_box   
  if(zb ~= []) then
    x2=zb(1) + (zb(3)-zb(1))*(x2-1)/(numxc-1)
    y2=zb(4) + (zb(2)-zb(4))*(numyc-y2)/(numyc-1)
  end
//
// Scale X,Y from reduced image to full-sized IP image
  x2=x2*ImageInfo.reduce -1
  y2=y2*ImageInfo.reduce -1
//
// Copy x2,y2 to x1,y1 in case we are not within a visible zoom image
  x1=x2
  y1=y2
//
// If mouse within a visible zoom, calculate x2,y2 from the zoom image
  if( IsZoomOn() ) then
// Shift X,Y and flip Y so that (1,1) is the lower-left corner of the zoom image
    x3=x+1 - ZoomAxes.axes_bounds(1)*numxc
    y3=201-y + ZoomAxes.axes_bounds(2)*numyc
// Test if mouse is inside the zoom image boundaries
    if( max(abs( [x3 y3]-101 )) < 100.5 ) then
// Correct for mouse-wheel zooming of zoom image
      zb=ZoomAxes.zoom_box   
      if(zb ~= []) then
        x3=zb(1) + (zb(3)-zb(1))*(x3-1)/200
        y3=zb(4) + (zb(2)-zb(4))*(201-y3)/200
      end
// Add the zoom image offset (plus 1/2)
      x1=ZoomAxes.user_data(1)+x3+0.5
      y1=ZoomAxes.user_data(2)+y3+0.5
    end

  end
//
// Ensure results are within data bounds
  x1=max(1,min(numx, x1 ))
  y1=max(1,min(numy, y1 ))
  x2=max(1,min(numx, x2 ))
  y2=max(1,min(numy, y2 ))
//
//
// Return empty arrays if an error occurs
catch
  x1=[]; y1=[]; x2=[]; y2=[]
end
//
endfunction


function sTemp=MouseHKLMessage(x,y,bZoom)
global LatticeInfo
//
try
//
  sTemp=""
  bUB=%f
  bTwins=%f
//
// Check if we have a valid default UB
  if( LatticeInfo ~= [] ) then
    if( LatticeInfo(1).level > 1 ) then
      if( LatticeInfo(1).ub ~= [] ) then
        bUB=%t
      end
    end
//
// Check if we have both twin UBs
    if( size(LatticeInfo,1) >= 3 ) then
      if( min( LatticeInfo(2:3).level ) > 1) then
        if( and(LatticeInfo(2:3).ub ~= []) ) then
          bTwins=%t
        end
      end
    end
//
  end
//
// Set radius around (x,y) to find an integer HKL
  radius=30
  if( bZoom ) then
    radius=10
  end
//
// If valid twins, write twin information to sTemp
  if( bTwins ) then
// Calculate hkl messages for both twins
    Linfo1=LatticeInfo(1)
    LatticeInfo(1)=LatticeInfo(2)
    sTemp1=MouseHKLString(x,y,radius)
    LatticeInfo(1)=LatticeInfo(3)
    sTemp2=MouseHKLString(x,y,radius)
    LatticeInfo(1)=Linfo1
// Combine messages after stripping sTemp2 of labels
    sTemp2=strsubst(strsubst(strsubst(sTemp2,"hkl=",""),"wav=",""),"d-space=","")
    sTemp="Twin1: "+sTemp1+"  Twin2: "+sTemp2
//
// Else, if valid default UB, add default lattice information to sTemp
  elseif( bUB ) then
    sTemp=MouseHKLString(x,y,radius)
  end
//
end
//
endfunction


function sTemp=MouseHKLString(x,y,radius)
//
try
//
// Try to guess hkl
  hkl=CalcHKLGuess([x,y],radius)
// If routine fails, write warning to sTemp
  if( hkl == [] ) then
    sTemp=">> ERROR: Restart LaueG <<"
    return
  end
//
// Calculate wavelength and d-spacing from hkl
  if(norm(hkl) ~= 0) then
    pixwav=HKL2PixWav(hkl)
    dspace=HKL2DSpacing(hkl)
  else
    pixwav=[0,0,0]
    dspace=0
  end
//
// For cursor_info=1 and hkl are integers:
  bIntHKL=( max(abs(hkl-round(hkl))) < 0.001 )
  if( (DisplayMode.cursor_info == 1) & bIntHKL ) then
    sTemp=msprintf( "hkl=(%3i,%3i,%3i) wav=%5.2f d-space=%5.2f", ...
                                  hkl,pixwav(3),dspace )
// For cursor_info=3: (display satellites around main spot)
  elseif( DisplayMode.cursor_info == 3 ) then
    [hkl_sat,wav,dxy]=CalcBestSatellite([x,y])
    sTemp=msprintf("hkl=(%5.3f,%5.3f,%5.3f) wav=%5.2f,  Diff(x,y)=%d %d", ...
                                  hkl_sat,wav,dxy )
// Otherwise: (display as fractional hkl vector)
  else
    sTemp=msprintf( "hkl-vector=(%5.3f,%5.3f,%5.3f) wav=%5.2f d-space=%5.2f", ...
                                  hkl,pixwav(3),dspace )
  end
//
end
//
endfunction


function [hkl_sat,wav,dxy]=CalcBestSatellite(xy)
//
// Calculate hkl from xy[], and normalise to the largest index
  hkl=Pix2HKL(xy)
  hkl=hkl/( max(abs(hkl)) + 1e-9)
//
// Generate satellites hkls around main spot
  hkl_sat=ModulateInfo.mults * ModulateInfo.vecs
  hkl_sat=hkl_sat + ones(hkl_sat(:,1))*DisplayData.modul_main_hkl
// Add main spot to start of hkl_sat
  hkl_sat=[DisplayData.modul_main_hkl; hkl_sat]
// Normalise satellite hkls to largest h,k,l index
  sizes=max( abs(hkl_sat), "c") + 1e-9
  hkl_norm=hkl_sat ./ ( sizes * [1,1,1] )
// 
// Find (possibly multiple) best matches of hkl & hkl_norm
  dhkl=ones(hkl_norm(:,1))*hkl - hkl_norm
  dsq=sum( dhkl.^2, 'c')
  ifind=find( abs(dsq-min(dsq)) < 1e-6)
//
// Accept "first" best match (so to preference main spot)
  hkl_sat=hkl_sat(ifind(1),:)
//
// Calculate wavelength and obs. to calc. pixel difference
  pixwav=HKL2PixWav(hkl_sat)
  wav=pixwav(3)
  dxy=xy-pixwav(1:2)
//
endfunction


function hkl=CalcHKLGuess(xy,radius)
//
// Return the best guess HKL for the image-plate coords x,y
// It now enforces centering rules if an integer hkl is found
//
// Scale h-vector so largest index is +/-1
  hvec=PixWav2HKL([xy,1])
  hvec=hvec/max([1e-6,abs(hvec)])
//
// Simply return normalised vector if cursor_info > 1
  if( DisplayMode.cursor_info > 1 ) then
    hkl=hvec
    return
  end
//
// Try multiplying h-vector by 1 to 9, and round to integer hkl
  hkl_try=round(hvec'*[1:9])
// Calculate coords for integer hkls
  pixwav=HKL2PixWav(hkl_try')
// Find first calculated coords within "radius" of (x,y)
  ifind=find( (pixwav(:,1)-xy(1)).^2 + (pixwav(:,2)-xy(2)).^2 < radius^2 ,1)
//
// If no valid integer HKL found, try multipliers 10 to 29
  if(ifind == []) then
    hkl_try=round(hvec'*[10:29])
    pixwav=HKL2PixWav(hkl_try')
    ifind=find( (pixwav(:,1)-xy(1)).^2 + (pixwav(:,2)-xy(2)).^2 < radius^2 ,1)
  end
//
// If no valid integer HKL found, try multipliers 30 to 99
  if(ifind == []) then
    hkl_try=round(hvec'*[30:99])
    pixwav=HKL2PixWav(hkl_try')
    ifind=find( (pixwav(:,1)-xy(1)).^2 + (pixwav(:,2)-xy(2)).^2 < radius^2 ,1)
  end
//
// If no valid integer HKL found, give up and return with a unit vector
  if(ifind == []) then
    hkl=hvec/norm(hvec)
    return
  end
//
// Load best integer values into hkl
  hkl=hkl_try(1:3,(ifind))'
//
// If hkl not allowed by centering, scale it so it is
  bOK=AllowedHKL( hkl, LatticeInfo(1) )
  if( ~bOK ) then
// If R centering (icen=6), then scale by 3
    if(LatticeInfo(1).icen == 6) then
      hkl=hkl*3
    else
// Else scale by 2
      hkl=hkl*2
    end
  end
//
// If scaled hkl outside +/-99, then revert hkl to a unit vector
  if(max(abs(hkl)) > 99) then
    hkl=hvec/norm(hvec)
  end
//
endfunction
