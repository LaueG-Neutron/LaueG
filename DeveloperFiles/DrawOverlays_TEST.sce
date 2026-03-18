global MainAxes MainMarks ImageFigure DisplayData RejectsInfo
global LatticeInfo ImageInfo MainLines ZoomMarks

// ===== Generate and draw, or erase, calculated spots =====
//function DrawGenHKLs(iTwin)
// ============= Routine to draw conic lines =============
//function DrawConicLines()
//function xy=CalcConicLine(Vconic)
// ============= Routines to draw various spots as "marks" =============
//function DrawObsSpots(spots)
//function DrawCalcSpots(spots)
//function DrawTwin1Spots(spots)
//function DrawTwin2Spots(spots)
//function DrawCalcModSpots(spots)
//function DrawTwin1ModSpots(spots)
//function DrawTwin2ModSpots(spots)
//function DrawUnmatchedSpots(bDraw)
//function DrawConicSpots()
//function DrawNodalSpots()
//function DrawRejectSpots()
// ================= Routines to retrieve X,Y or number of drawn marks =============
//function ndraw=NumDrawnSpots(itype)
//function xy_pos=GetMarksXY(imark)
// ================= Routines store and retrieve the marks =============
//function DrawMarksNow(imark,xy_pos)
//function RedrawImageMarks()
//function RedrawZoomMarks()
// ================= Routines to setup marks and lines =============
//function SetupImageMarks()
//function SetupZoomMarks()

// Overlay marks:
// 1    Observed spots
// 2    Calculated spots (default)
// 3    Calculated spots (twin 1)
// 4    Calculated spots (twin 2)
// 5    Calc. satellite spots (default)
// 6    Calc. satellite spots (twin 1)
// 7    Calc. satellite spots (twin 2)
// 8    Unmatched spots
// 9    Conic spots
// 10   Nodal spots
// 11   Rejected spots


// ===== Generate and draw, or erase, calculated spots =====

function DrawGenHKLs(iTwin)
//
// Generate and draw HKLs for:
//  iTwin=0     default lattice
//  iTwin=1,2   Twins 1 or 2
//  iTwin=12    Twins 1 and 2
//
// Convert iTwin to imode=1,2,3,4
  imode=find(iTwin == [0,1,2,12])
  if(imode == []) then BugBox("Invalid iTwin"); end
//
// Ask for confirmation when a large number of spots
  nguess=GuessCalcSpotsNumber(LatticeInfo(min(3,imode)))
  if(nguess > 1e4) then
    text=[sprintf("Estimated number of spots is large (~%d,000)", ...
                               nguess/1000),"Do you wish to display them?"]
    button=messagebox(text,"Confirmation", "question", ["Yes","No"], "modal")
// Answer is no, so erase corresponding spots and return
    if(button == 2) then
      if(imode < 4) then
        DrawMarksNow(imode+1,[])
        DrawMarksNow(imode+4,[])
      else
        DrawMarksNow(3,[])
        DrawMarksNow(4,[])
        DrawMarksNow(6,[])
        DrawMarksNow(7,[])
      end
      UnlockMenus()
      return
    end
  end
//
  LockMenus("Generating calculated spots")
//
// Generate calculated spots using current display settings
  if(imode < 4) then
    hkllist=RunGenHKLs(LatticeInfo(imode), ...
                       DisplaySet.wav_min,DisplaySet.wav_max, ...
                       DisplaySet.d_min,DisplaySet.d_max, ...
                       DisplayMode.modul_show)
////////////////// ??????????????????? <<<<<<<<<<<<<<<<<<<<
// Test code for new peak shape algorithm
if( %f ) then
//
    global OrientInfo
//
    off0 = OrientInfo.xtaloff;
//    doff = [-1.000,-0.5,0.265];    // Si1_vert_phis_001
    doff = [-1.100,-1.2,0.1];      // Si1_obliq_30deg_phis_001

    OrientInfo.xtaloff=off0-doff;
    hkllist1=RunGenHKLs(LatticeInfo(1), DisplaySet.wav_min, DisplaySet.wav_max, ...
                   DisplaySet.d_min, DisplaySet.d_max, DisplayMode.modul_show);
    OrientInfo.xtaloff=off0+doff;
    hkllist2=RunGenHKLs(LatticeInfo(1), DisplaySet.wav_min, DisplaySet.wav_max, ...
                   DisplaySet.d_min, DisplaySet.d_max, DisplayMode.modul_show);
    OrientInfo.xtaloff=off0;
//
    DrawMarksNow(1,[hkllist1(:,5:6);hkllist2(:,5:6)])
else
    ifind=find( hkllist(:,9) == 0 )             // Main spots
    DrawMarksNow(imode+1,hkllist(ifind,5:6))
    ifind=find( hkllist(:,9) ~= 0 )             // Satellite spots
    DrawMarksNow(imode+4,hkllist(ifind,5:6))
end
////////////////// ??????????????????? <<<<<<<<<<<<<<<<<<<<
  else
    hkllist=RunGenHKLs(LatticeInfo(2), ...
                       DisplaySet.wav_min,DisplaySet.wav_max, ...
                       DisplaySet.d_min,DisplaySet.d_max, ...
                       DisplayMode.modul_show)
    ifind=find( hkllist(:,9) == 0 )
    DrawMarksNow(3,hkllist(ifind,5:6))
    ifind=find( hkllist(:,9) ~= 0 )
    DrawMarksNow(6,hkllist(ifind,5:6))
    hkllist=RunGenHKLs(LatticeInfo(3), ...
                       DisplaySet.wav_min,DisplaySet.wav_max, ...
                       DisplaySet.d_min,DisplaySet.d_max, ...
                       DisplayMode.modul_show)
    ifind=find( hkllist(:,9) == 0 )
    DrawMarksNow(4,hkllist(ifind,5:6))
    ifind=find( hkllist(:,9) ~= 0 )
    DrawMarksNow(7,hkllist(ifind,5:6))
  end
//
  UnlockMenus()
//
endfunction


// ============= Routines to draw conic lines =============

function DrawConicLines()
global MainLines
// Draws both the conic or nodal lines
//
// Ignore if no valid image
  if( ~is_handle_valid(ImageFigure) ) then
    return
  end
//
// Delay image draw
  scf(ImageFigure)
  drawlater
//
// Make sets of xy coords from nodal or conic spots
// 2 sets if 3 nodal spots, 1 set if 2 nodal spots, 1 set if >1 conic spots
  xy_spots1=[]
  xy_spots2=[]
  if(size(DisplayData.nodal_spots,1) > 1) then
    xy_spots1=DisplayData.found_spots(DisplayData.nodal_spots([1,2]),1:2)
    if(size(DisplayData.nodal_spots,1) > 2) then
      xy_spots2=DisplayData.found_spots(DisplayData.nodal_spots([1,3]),1:2)
    end
  elseif(size(DisplayData.conic_spots,1) > 1) then
    xy_spots1=DisplayData.found_spots(DisplayData.conic_spots,1:2)
  end
//
// Erase any existing lines
  MainLines(1).user_data=[]
  MainLines(1).data=[]
  MainLines(2).user_data=[]
  MainLines(2).data=[]
//
// Fit spots to find conic vector, create polyline coords, copy polyline to figure data
// Draw line 1 if xy_spots1 is not empty
  if(xy_spots1 ~= []) then
    Vconic=FitConicSpots(xy_spots1)
    xy_line=CalcConicLine(Vconic)
    MainLines(1).user_data=xy_line
    MainLines(1).data=xy_line/ImageInfo.reduce
  end
// Draw line 2 if xy_spots2 is not empty
  if(xy_spots2 ~= []) then
    Vconic=FitConicSpots(xy_spots2)
    xy_line=CalcConicLine(Vconic)
    MainLines(2).user_data=xy_line
    MainLines(2).data=xy_line/ImageInfo.reduce
  end
//
// Redraw the zoom image (if visible)
  if(IsZoomOn()) then
    RedrawZoomImage()
  end
//
  drawnow
//
endfunction


function xy=CalcConicLine(Vconic)
// Make 3 arbitrary vectors perpendicular to the conic
  v1=Vconic(1)*Vconic-[1,0,0]
  v2=Vconic(2)*Vconic-[0,1,0]
  v3=Vconic(3)*Vconic-[0,0,1]
// If v1 or v2 is tiny, use v3 instead
  if(norm(v1) < 0.1) then v1=v3; end
  if(norm(v2) < 0.1) then v2=v3; end
// Make v1 & v2 orthonormal
  v1=v1/norm(v1)
  v2=v2-(v1*v2')*v1
  v2=v2/norm(v2)
// Make up 180 points on the conic line
  v=cosd(1:180)'*v1+sind(1:180)'*v2
  xy=Hvec2Pix(v)
// Find the biggest step in x
  dx=xy(2:$,1)-xy(1:$-1,1)
  [val,ipos]=max(abs(dx))
// If the step > image width, move the step to start of the array
  if( val > ImageInfo.numxy(1)) then
    xy=[xy(ipos+1:$,:); xy(1:ipos,:)]
  else      // No big step, so make the line a closed loop
    xy($+1,:)=xy(1,:)
  end
//
endfunction


// ============= Routines to draw various spots as "marks" =============

function DrawObsSpots(spots)
// Draw all observed spots in the array spots[]
  DrawMarksNow(1,spots(:,1:2))
endfunction


function DrawCalcSpots(spots)
// Draw all calculated spots in the array spots[]
  DrawMarksNow(2,spots)
endfunction


function DrawTwin1Spots(spots)
// Draw all Twin 1 spots in the array spots[]
  DrawMarksNow(3,spots)
endfunction


function DrawTwin2Spots(spots)
// Draw all Twin 2 spots in the array spots[]
  DrawMarksNow(4,spots)
endfunction


function DrawCalcModSpots(spots)
// Draw satellite spots for default lattice
  DrawMarksNow(5,spots)
endfunction


function DrawTwin1ModSpots(spots)
// Draw satellite spots for Twin 1 lattice
  DrawMarksNow(6,spots)
endfunction


function DrawTwin2ModSpots(spots)
// Draw satellite spots for Twin 2 lattice
  DrawMarksNow(7,spots)
endfunction


function DrawUnmatchedSpots(bDraw)
// If bDraw=%f, erase the spots instead of draw them
  if( bDraw ) then
    DrawMarksNow(8,CalcUnmatchedSpots())
  else
    DrawMarksNow(8,[])
  end
endfunction


function DrawConicSpots()
// Draw the conic spots on the main (and zoom) image
  xy_pos=DisplayData.found_spots( DisplayData.conic_spots(:) , 1:2 )
  DrawMarksNow(9,xy_pos)
endfunction


function DrawNodalSpots()
// Draw the nodal spots on the main (and zoom) image
  xy_pos=DisplayData.found_spots( DisplayData.nodal_spots(:) , 1:2 )
  DrawMarksNow(10,xy_pos)
endfunction


function DrawRejectSpots()
//
// Get the marked spots that are being displayed
  nspots=size(RejectsInfo.itwin,1)
  ifirst=RejectsInfo.ifirst

  ilast=min(ifirst+49,nspots)
  imarks=find(RejectsInfo.marked(ifirst:ilast))
//
// If no marks, erase existing marks and return
  if( imarks == [] ) then
    DrawMarksNow(11,[])
    return
  end
//
// Load pixels positions of marks as if a normal image
  ix=round(modulo(imarks-0.1,10))
  iy=1+floor((imarks-1)/10)
  xy_pos=[-50+100*ix; 550-iy*100]
  xy_pos=ImageInfo.reduce*xy_pos'
//
// Draw the marks
  DrawMarksNow(11,xy_pos)
//
endfunction


// ================= Routines to retrieve X,Y or number of drawn marks =============

function ndraw=NumDrawnSpots(itype)
// Returns the number of drawn observed/calc/twin1/twin2 spots
// for itype = -1/0/1/2
  if( itype == -1 ) then
    ndraw=size(  MainMarks(1).user_data ,1)
  elseif( itype <= 2 ) then
// Include both normal and modulated spots
    ndraw=size(  MainMarks(itype+2).user_data ,1) + ...
          size(  MainMarks(itype+5).user_data ,1)
  else
    BugBox("Invalid itype = "+string(itype))
  end
endfunction


function xy_pos=GetMarksXY(imark)
// Returns the X,Y of marks of type "imark"
  if(MainMarks == []) then
    xy_pos=[]
  else
    xy_pos=MainMarks(imark).user_data
  end
//
endfunction


// ================= Routines store and draw the marks =============

function DrawMarksNow(imark,xy_pos)
global MainMarks
// Draw marks of type "imark" at x,y given by xy_pos
// Marks are stored in "MainMarks(imark).user_data"
//
// Ignore if no valid image
  if( ~is_handle_valid(ImageFigure) ) then
    return
  end
//
// Delay image draw triggered by change to marks data
  scf(ImageFigure)
  drawlater
//
// Update marks data
  if(xy_pos == []) then
    MainMarks(imark).user_data=[]
    MainMarks(imark).data=[]
  else
    MainMarks(imark).user_data=xy_pos
// NB: Add 1/4 to the drawn X,Y (for reduce=4)
    MainMarks(imark).data=(xy_pos+1)/ImageInfo.reduce
  end
//
// Redraw the zoom image (if visible)
  if(IsZoomOn()) then
    RedrawZoomImage()
  end
//
  drawnow
//
endfunction


function RedrawImageMarks()
global MainMarks
// Draw all main image marks stored in user_data arrays
  for i=1:size(MainMarks,1)
// NB: Add 1/4 pixel to the drawn X,Y
    if(MainMarks(i).data ~= []) then
      MainMarks(i).data=(MainMarks(i).user_data+1)/ImageInfo.reduce
    end
  end
endfunction


function RedrawZoomMarks()
global ZoomMarks
// Draw all zoom image marks stored in user_data arrays
  for i=1:size(ZoomMarks,1)
    if(MainMarks(i).user_data == []) then
      ZoomMarks(i).data=[]
    else
      ix0=ZoomAxes.user_data(1)
      iy0=ZoomAxes.user_data(2)
      ZoomMarks(i).data=[ MainMarks(i).user_data(:,1)-ix0 , ...
                          MainMarks(i).user_data(:,2)-iy0 ]
    end
  end
endfunction


// ================= Routines to setup marks and lines =============

function SetupImageMarks()
global MainMarks MainLines
// Setup marks to be displayed on the main image
//
// Set the current axes to that of the main image
  sca(MainAxes)
//
//Available colours:
// 257 red, 258 green, 259 blue, 260 yellow, 261 magenta, 262 cyan
// 263 orange, 264 sky blue, 265 mauve, 266 pink, 267 olive
//
// Different types of marks
marks=[...
         9   5    264  // blue circle
         5   6    263  // orange diamond
         5   6    260  // yellow diamond
         5   6    261  // mauve diamond
         1   6    263  // orange cross "+"
         1   6    260  // yellow cross "+"
         1   6    261  // mauve cross "+"
         9   5    261  // mauve circle
         2   7    264  // blue cross "x"
         9   11   264  // blue circle (large)
         2   50   257  // red cross "x" (very large)
      ]
// Draw polylines for each mark, and save results in MainMarks
  MainMarks=[]
  for i=1:size(marks,1)
    xpoly(0,0,"lines") // polyline with 1 point (0 not allowed)
    MainMarks(i)=gce()
    MainMarks(i).line_mode = "off"
    MainMarks(i).mark_mode = "on"
    MainMarks(i).mark_style = marks(i,1)
    MainMarks(i).mark_size = marks(i,2)
    MainMarks(i).mark_size_unit = "point"
    MainMarks(i).mark_foreground = marks(i,3)
    MainMarks(i).mark_background = 0
    MainMarks(i).data=[]
  end
//
// Draw polylines for each line, and save results in MainLines
// All 3 lines are dashed and blue
  MainLines=[]
  for i=1:3
    xpoly(0,0,"lines")
    MainLines(i)=gce()
    MainLines(i).line_mode = "on"
    MainLines(i).mark_mode = "off"
    MainLines(i).foreground = 259
    MainLines(i).line_style = 8
    MainLines(i).thickness=2
    MainLines(i).data=[]
  end
//
endfunction


function SetupZoomMarks()
global ZoomMarks
// Setup marks to be displayed on the zoom image
// NB: Assumes ZoomAxes is the current axes
//
// Setup arrays marks as in the main image
// NB: The marks are larger in the zoom image
marks=[...
         9   10   264  // blue circle
         5   9    263  // orange diamond
         5   9    260  // yellow diamond
         5   9    261  // mauve diamond
         1   9    263  // orange cross "+"
         1   9    260  // yellow cross "+"
         1   9    261  // mauve cross "+"
         9   8    261  // mauve circle
         2   12   264  // blue cross "x"
         9   26   264  // blue circle (large)
         2   50   257  // red cross "x" (very large)
      ]
//
// Draw polylines for each mark, and save results in ZoomMarks
  ZoomMarks=[]
  for i=1:size(marks,1)
    xpoly(0,0,"lines") // polyline with 1 point (0 not allowed)
    ZoomMarks(i)=gce()
    ZoomMarks(i).line_mode = "off"
    ZoomMarks(i).mark_mode = "on"
    ZoomMarks(i).mark_style = marks(i,1)
    ZoomMarks(i).mark_size = marks(i,2)
    ZoomMarks(i).mark_size_unit = "point"
    ZoomMarks(i).mark_foreground = marks(i,3)
    ZoomMarks(i).mark_background = 0
    ZoomMarks(i).data=[]
  end
endfunction
