global ReducedImage MainImage ImageInfo DisplayData
global DisplaySet DisplayMode ZoomAxes RejectsInfo

//function DrawCurrentEllipses()
//function DrawIntegEllipses()
//function DrawModelEllipses()
//function DrawRejectEllipses()
//function LoadEllipses(imults,itypes)
//function DrawMainEllipses()
//function DrawZoomEllipses()
//function ClearEllipses()


function DrawCurrentEllipses()
// Draw none/integ./model ellipses depending on DisplayMode.*
  if( DisplayMode.ell_integ ) then
    DrawIntegEllipses()
  elseif( DisplayMode.ell_model ) then
    DrawModelEllipses()
  elseif( DisplayMode.ell_rejects ) then
    DrawRejectEllipses()
  end
endfunction


function DrawIntegEllipses()
global DisplayMode
// Don't use ellipses in image-mode
  if( DisplayMode.image_mode ) then return; end
// Draw integration ellipses in colours depending on status
  itypes(1,:)=[1,259]    // blue = full integration
  itypes(2,:)=[2,262]    // cyan = sigI/I, no overlap
  itypes(3,:)=[3,258]    // green = sigI/I, P-P overlap
  itypes(4,:)=[4,260]    // yellow = sigI/I, C-P overlap
  itypes(5,:)=[5,263]    // orange = sigI/I, C-C overlap
  itypes(6,:)=[9,257]    // red = failed integration
  itypes(7:12,1)=itypes(1:6,1)+10  // add the model spots
  itypes(7:12,2)=itypes(1:6,2)
// Read the ellipse file and load the ellipses
  if( ~DisplaySet.ell_back ) then
    LoadEllipses([1,2],itypes)
  else
    LoadEllipses([1,2,3],itypes)
  end
// Draw the ellipses on the main image
  DrawMainEllipses()
// Redraw the zoom image with ellipses, if zoom is visible
  if( IsZoomOn() ) then 
    RedrawZoomImage()
  end
// Save which type of ellipse is being displayed
  DisplayMode.ell_integ=%t
  DisplayMode.ell_model=%f
  DisplayMode.ell_rejects=%f
endfunction


function DrawModelEllipses()
global DisplayMode
// Don't use ellipses in image-mode
  if( DisplayMode.image_mode ) then return; end
// Draw the model ellipses in green
  itypes(:,1)=[11:15,19]
  itypes(:,2)=258
// Read the ellipse file and load the ellipses
  if( ~DisplaySet.ell_back ) then
    LoadEllipses([1,2],itypes)
  else
    LoadEllipses([1,2,3],itypes)
  end
// Draw the ellipses on the main image
  DrawMainEllipses()
// Redraw the zoom image with ellipses, if zoom is visible
  if( IsZoomOn() ) then 
    RedrawZoomImage()
  end
// Save which type of ellipse is being displayed
  DisplayMode.ell_integ=%f
  DisplayMode.ell_model=%t
  DisplayMode.ell_rejects=%f
endfunction


function DrawRejectEllipses()
global MainImage DisplayMode
//
// Calculate first and last reject to display
  ifirst=RejectsInfo.ifirst
  nspots=size(RejectsInfo.itwin,1)
  ilast=min(ifirst+49,nspots)
//
// Get the ellipse parameters for displayed rejects
  if( ~DisplaySet.ell_back ) then
    amult=RejectsInfo.cont(ifirst:ilast,1:2)
  else
    amult=RejectsInfo.cont(ifirst:ilast,1:3)
  end
  cens=RejectsInfo.xycen(ifirst:ilast,1:2)
  pars=RejectsInfo.efh(ifirst:ilast,1:3)
//
// Calculate the X,Y limits for each reject image
  ixcen=round(cens(:,1))
  iycen=round(cens(:,2))
  ixy_lim=[ixcen-49,ixcen+49,iycen-49,iycen+49]
//
// Calculate pixels on MainImage that are drawn as ellipses
  ipoint=RunMakeRejectEllipses(amult,cens,pars,ixy_lim)
//
// Update the ellipses pixels with colour #258 (green)
  MainImage.data(ipoint)=258
// Save which type of ellipse is being displayed
  DisplayMode.ell_integ=%f
  DisplayMode.ell_model=%f
  DisplayMode.ell_rejects=%t
endfunction


function LoadEllipses(imults,itypes)
global ImageInfo DisplayData
// Load ellipses from the *.ell file for specified types
// of spots with each type having a specified colour.
//
// Spots with an "ellipse status" matching itypes(:,1) are
// given the colour of the corresponding itypes(:,2) value.
// Spots not matching any itypes(:,1) are ignored.
// imults[] can contain values of 1, 2 or 3 which correspond
// to core, peak and background ellipses being drawn. These
// values will be stored in DisplayData.ell_amults.
//
  file_name=ImageInfo.basename+".ell"
  [hkl,xy_efh,amult,imodel,itwin]=ReadEllFile(file_name)
// Store the requested ellipse area multipliers
  DisplayData.ell_amults=amult(1,imults)
// Store the required ellipses and colours
  DisplayData.ell_data=[]
  for ityp=itypes'
    itag=find(imodel == ityp(1))
    if(itag ~= []) then
      xy_efh2=xy_efh(itag,1:5)
      xy_efh2(:,6)=ityp(2)
      DisplayData.ell_data=[DisplayData.ell_data;xy_efh2]
    end
  end
endfunction


function DrawMainEllipses()
// Draw the ellipses on the main image directly into MainImage
global MainImage
// Don't use ellipses in image-mode
  if( DisplayMode.image_mode ) then return; end
// Nothing to do if no ellipse information
  if(DisplayData.ell_data == []) then return; end
// Reload main image
  MainImage.data=CalcPixmap(ReducedImage,MainImage.user_data)
// Get the ellipse area multipliers
  amult=DisplayData.ell_amults
// Load the X,Y limits of the main image
  ixy_lim=round([1,ImageInfo.numxy(1)/ImageInfo.reduce, ...
                  1,ImageInfo.numxy(2)/ImageInfo.reduce])
// Loop through each of the unique ellipse colours
  for icol=unique(DisplayData.ell_data(:,6))'
// Get the ellipse data for a particular colour
// Scale cen & pars for the 1/ImageInfo.reduce size of the main image
    itags=find(DisplayData.ell_data(:,6) == icol)
    cens=DisplayData.ell_data(itags,1:2)/ImageInfo.reduce
    pars=DisplayData.ell_data(itags,3:5)*ImageInfo.reduce^2
// Calculate and draw the ellipses for the amult() values
    ixy=RunMakeEllipses(cens,pars,amult,ixy_lim)
    MainImage.data( ixy_lim(4)*round(ixy(:,1))-round(ixy(:,2))+1 )=icol
  end
endfunction


function DrawZoomEllipses()
// Draw the ellipses on the zoom image directly into MainImage
global ZoomImage
// Don't use ellipses in image-mode
  if( DisplayMode.image_mode ) then return; end
// Nothing to do if no ellipse information
  if(DisplayData.ell_data == []) then return; end
// Calculate the X,Y limits for the zoom image
  ix0=ZoomAxes.user_data(1)
  iy0=ZoomAxes.user_data(2)
  ixy_lim=[ix0,ix0-1+200,iy0,iy0-1+200];
// Get the ellipse area multipliers
  amult=DisplayData.ell_amults
// Loop through each of the unique ellipse colours
  for icol=unique(DisplayData.ell_data(:,6))'
// Get the ellipse data for a particular colour
    itags=find(DisplayData.ell_data(:,6) == icol)
    cens=DisplayData.ell_data(itags,1:2)
    pars=DisplayData.ell_data(itags,3:5)
// Calculate and draw the ellipses for the amult() values
    ixy=RunMakeEllipses(cens,pars,amult,ixy_lim)
    if(ixy ~= []) then
      ixy(:,1)=round(ixy(:,1)-ix0+1)
      ixy(:,2)=round(ixy(:,2)-iy0+1)
      ZoomImage.data( 200*ixy(:,1)-ixy(:,2)+1 )=icol
    end
  end
endfunction


function ClearEllipses()
global DisplayData DisplayMode
// Clear any ellipses data
  DisplayData.ell_data=[]
// Turn off both types of ellipse being displayed
  DisplayMode.ell_integ=%f
  DisplayMode.ell_model=%f
// Redraw the images using the current display limits
  RedrawImageNow(MainImage.user_data)
endfunction
