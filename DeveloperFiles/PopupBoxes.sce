// Creation and callback routines for all popup boxes
global OrientInfo LatticeInfo ModulateInfo ImageInfo
global DisplayData DisplayMode DisplaySet ArgboxSet
global MainImage MainMarks BatchFiles FlickData ProgDir DataDir


// -------- Popups to input and modify lattices -------
//function TransLattice_Popup()
//function TransLattice_CB()
//function HKLIndex3_Popup()
//function HKLIndex3_CB()
//function InputCell_Popup(bIndex)
//function InputCell_VisCB()
//function InputCell_CB()
//function [bOK,ilatt,icen,cell_dims]=ReadCellPopupValues(f)
//function cell_dims=ConvertCellDims(cell_dims,ilatt)
//function BatchCenter_Popup(icen)
//function BatchCenter_CB()
//function ConvertLauegen_Popup()
//function ConvertLauegen_CB()
// ------------ Popups to process spots -----------
//function ManualOrient_Popup()
//function ManualOrient_CB()
//function AdvIndex_Popup()
//function AdvIndex_CB()
//function RefineOrient_Popup()
//function RefineOrient_CB()
//function f=OrientTwins_Popup(bIndex)
//function OrientTwins_CB()
//function ArgonneBoxes_Popup(bOldMode,bTwinMode,bModulateMode)
//function f=ArgBoxNew_Popup()
//function f=ArgBoxOld_Popup()
//function ArgonneBoxes_CB()
// ---------- Popups to change image & marks --------
//function Brightness_Popup()
//function Brightness_CB()
//function AdvPrune_Popup()
//function AdvPrune_CB()
//function FindHKLSpot_Popup()
//function FindHKLSpot_CB()
//function ObsSpots_Popup()
//function ObsSpots_CB()
//function CalcSpots_Popup()
//function CalcSpots_CB()
// ------- Popups to select files from a list -------
//function RenameFiles_Popup(get_files)
//function RenameFiles_CB()
//function IntegFileSelect_Popup(iFile)
//function IntegFileSelect_CB()
//function WaveFileSelect_Popup()
//function WaveFileSelect_CB()
// -------- Popups to manage modulation vectors  --------
//function ModulMainSpot_Popup()
//function ModulMainSpot_CB()
//function ShowModul_Popup()
//function ShowModul_CB()
//function ChangeModul_Popup()
//function ChangeModul_CB()
// ------- Convenience routines for popup boxes  -------
//function f=CreateEmptyPopup(siz,str)
//function CloseExistingPopup(name)
//function SetPopupEdits(f,strings)
//function [vals]=ReadPopupEditsNumber(f,nvals)
//function [strings]=ReadPopupEdits(f)
//function [vals]=ReadPopupChecks(f)
//function SetupHintUI()
//function HintUI_CB()
// ------- Routines for adding widgets to popup boxes  -------
//function h=AddTextUI(fig,pos,str,siz)
//function h=AddEditUI(fig,pos,str,siz)
//function h=AddListUI(fig,pos,str,siz,val)
//function h=AddCheckUI(fig,posxy,bVal)
//function h=AddButtonUI(fig,pos,str,siz)
//function AddUICallback(h,cb,tag)
//function AddHintUI(fig,posxy,text)


///////////////////////////////////////////////////////////
// -------- Popups to input and modify lattices -------- //
///////////////////////////////////////////////////////////

function f=TransLattice_Popup()
//
// Create a popup to transform the unit cell
//
// Create a popup window without the menus and toolbar
  f=CreateEmptyPopup([240 425],"Lattice Transformation")
//
// Output original cell dimensions and volume
  AddTextUI(f,[50 400 160 25],"Original Cell Dimensions",13)
  str=msprintf("%.3f %.3f %.3f %.2f %.2f %.2f",LatticeInfo(1).cell)
  AddTextUI(f,[20 385 200 20],str,11)
  str=msprintf("( Cell Volume = %.f )",1/abs(det(LatticeInfo(1).ub)))
  AddTextUI(f,[65 370 200 20],str,11)
//
// Create input for matrix and checkbox for inverse transform
  AddTextUI(f,[ 10 345 130 25],"Transformation Matrix:",13)
  AddEditUI(f,[150 340 80 25],"1 0 0",13)
  AddEditUI(f,[150 315 80 25],"0 1 0",13)
  AddEditUI(f,[150 290 80 25],"0 0 1",13)
  AddTextUI(f,[ 30 320 60 20],"Inverse",13)
  AddTextUI(f,[ 20 300 100 20],"HKL mode",13)
  AddCheckUI(f,[85 320],%f)
  AddCheckUI(f,[85 300],%f)
//
// Create first listbox for cell type
  str="Triclinic|Monoclinic|Orthorhombic|Tetragonal|" + ...
      "Cubic|Trigonal (rhom)|Trigonal (hex)|Hexagonal"
  AddListUI(f,[40 105 95 165],str,13,1)
  f.children(1).value=LatticeInfo(1).ilatt
//
// Create listbox for centering
  AddListUI(f,[160 95 35 190],"| P| A| B| C|  I| F| R",18,1)
  f.children(1).value=LatticeInfo(1).icen+1
//
//// Output the transformed cell dimensions and volume
  AddTextUI(f,[55 75 130 20],"Final Cell Dimensions",13)
  str=msprintf("%.3f %.3f %.3f %.2f %.2f %.2f",LatticeInfo(1).cell)
  AddTextUI(f,[20 60 240 20],str,11)
  str=msprintf("( Cell Volume = %.f )",1/abs(det(LatticeInfo(1).ub)))
  AddTextUI(f,[65 45 200 20],str,11)
//
// Pushbuttons for accept/test/cancel, and add callbacks
  h1=AddButtonUI(f,[ 10  7 70 30],"Accept",18)
  h2=AddButtonUI(f,[ 90  7 60 30],"Test",18)
  h3=AddButtonUI(f,[160  7 70 30],"Cancel",18)
  AddUICallback(h1,"TransLattice_CB","accept")
  AddUICallback(h2,"TransLattice_CB","test")
  AddUICallback(h3,"TransLattice_CB","cancel")
//
// Prevent window closing except by the pushbuttons
  f.closerequestfcn="InfoBox(""You must use the push buttons to exit this popup"")"
//
// Pass original LatticeInfo(1) to callback
  f.user_data=LatticeInfo(1)
//
endfunction


function TransLattice_CB()
global LatticeInfo
// Get graphics object, parent figure
  g=gcbo()
  f=g.parent
//
// Load the original lattice
  LatticeInfo(1)=f.user_data
//
// If "cancel", draw original spots and close popup
  if(g.tag == "cancel") then
    DrawGenHKLs(0)
    close(f)
    return
  end
//
// Read matrix, complain and return if invalid
  vals=ReadPopupEditsNumber(f,9)
  if(vals == []) then return; end
  mat=matrix(vals, 3,3)
// If singular matrix, complain and return
  if(abs(det(mat)) < 1e-6) then
    ErrorBox("Invalid transformation matrix")
    return
  end
//
// Get check-box values for "inverse" and "hkl mode"
  check=ReadPopupChecks(f)
//
// Invert matrix, if inverse selected
  if( check(1) ) then
    mat=1/mat
  end
//
// Invert and transpose matrix, if hkl mode selected
  if( check(2) ) then
    mat=1/(mat')
  end
//
// If "Accept" and we have a left hand lattice, suggest a fix
  if(g.tag == "accept") then
    if(det(LatticeInfo(1).ub / mat) < 0) then
      text=[ "Use an equivalent matrix to", ...
             "create a right handed system?"]
      if( QuestionBox(text) ) then
        mat=-mat;
      end
    end
  end
//
// Transform ub, calculate new cell and U matrix
  LatticeInfo(1).ub=LatticeInfo(1).ub / mat
  LatticeInfo(1).cell=UB2UnitCell(LatticeInfo(1).ub)
  u_mat=LatticeInfo(1).ub / Cell2Bmatrix(LatticeInfo(1).cell)
//
// Get the selected lattice & centering
  LatticeInfo(1).ilatt=f.children(8).value
  LatticeInfo(1).icen=f.children(7).value-1
//
// If the cell changes a lot against lattice type (>1 Angstrom, >10 degrees),
// assume the lattice type is incorrect and set it to triclinic
  dcell=LatticeInfo(1).cell - ConvertCellDims(LatticeInfo(1).cell,LatticeInfo(1).ilatt)
  itmp=0
  for i=6:-1:4
    if(abs(dcell(i)) > 10) then itmp=i; end
  end
  for i=3:-1:1
    if(abs(dcell(i)) > 1) then itmp=i; end
  end
  if( itmp > 0) then
    LatticeInfo(1).ilatt=1
    f.children(8).value=1
    tmp=["a","b","c","alpha","beta","gamma"]
    WarnBox("Using Triclinic due to inconsistent """+tmp(itmp)+""" value")
  end
//
// Adjust the cell for the lattice type, and recalculate UB
  LatticeInfo(1).cell=ConvertCellDims(LatticeInfo(1).cell,LatticeInfo(1).ilatt)
  LatticeInfo(1).ub=u_mat * Cell2Bmatrix(LatticeInfo(1).cell)
//
// Update the final cell dimensions and volume on popup box
  f.children(5).String=msprintf("%.3f %.3f %.3f %.2f %.2f %.2f",LatticeInfo(1).cell)
  f.children(4).String=msprintf("( Cell Volume = %.f )",1/abs(det(LatticeInfo(1).ub)))
//
// Generate and draw calculated spots
  DrawGenHKLs(0)
//
// If "Accept", close popup
  if(g.tag == "accept") then
    close(f)
    return
  end
//
endfunction


function HKLIndex3_Popup()
//
// Create a popup to transform the UB using 3 spots
//
// Create a popup window and strip off the menus and toolbar
  f=CreateEmptyPopup([360 300],"Change orientation using 3 Spots")
//
// Output original cell dimensions and volume
  AddTextUI(f,[140 270 120 25],"Original Cell",15)
  AddTextUI(f,[80 255 290 20],"? ? ? ? ? ?",11)
  AddTextUI(f,[120 240 120 20],"( Cell Volume = ? )",11)
  if(LatticeInfo(1).level > 1) then
    f.children(2).string=msprintf("%.3f %.3f %.3f %.2f %.2f %.2f",LatticeInfo(1).cell)
    f.children(1).string=msprintf("( Cell Volume = %.f )",1/abs(det(LatticeInfo(1).ub)))
  end
//
// Create inputs for HKL,XY,WAV values and outputs for calculated values
  AddTextUI(f,[ 30 200 40 25],"H,K,L",13)
  AddTextUI(f,[100 200 110 25],"X, Y, Wav. (calc.)",13)
  AddTextUI(f,[230 200 110 25],"X, Y, Wav. (new)",13)
  for i=1:3
    iy=175-(i-1)*30
    AddEditUI(f,[10 iy 70 25],"0 0 0",13)
    AddTextUI(f,[100 iy 100 25],"? ? ?",13)
    AddEditUI(f,[220 iy 130 25],"",13)
  end
//
// Output original cell dimensions and volume
  AddTextUI(f,[150 75 120 25],"New Cell",15)
  AddTextUI(f,[80 60 290 20],"? ? ? ? ? ?",11)
  AddTextUI(f,[120 45 120 20],"( Cell Volume = ? )",11)
  if(LatticeInfo(1).level > 1) then
    f.children(2).string=msprintf("%.3f %.3f %.3f %.2f %.2f %.2f",LatticeInfo(1).cell)
    f.children(1).string=msprintf("( Cell Volume = %.f )",1/abs(det(LatticeInfo(1).ub)))
  end
//
// Pushbuttons for accept/test/cancel, and add callbacks
  h1=AddButtonUI(f,[ 30  7 70 30],"Accept",18)
  h2=AddButtonUI(f,[150  7 60 30],"Test",18)
  h3=AddButtonUI(f,[260  7 70 30],"Cancel",18)
  AddUICallback(h1,"HKLIndex3_CB","accept")
  AddUICallback(h2,"HKLIndex3_CB","test")
  AddUICallback(h3,"HKLIndex3_CB","cancel")
//
// Prevent window closing except by the pushbuttons
  f.closerequestfcn="InfoBox(""You must use the push buttons to exit this popup"")"
//
// Pass original lattice to callback
  f.user_data=LatticeInfo(1)
//
endfunction

function HKLIndex3_CB()
global LatticeInfo
//
// Get graphics object, parent figure
  g=gcbo()
  f=g.parent
//
// Get original lattice
    Linfo=f.user_data
//
// If "cancel", revert to original lattice, redraw spots and close popup
  if(g.tag == "cancel") then
    LatticeInfo(1)=Linfo
    DrawGenHKLs(0)
    close(f)
    return
  end
//
// If we have a UB:
// Read HKLs and calc. X,Y,Wav (output but ignore any problems)
  if(LatticeInfo(1).level > 1 ) then
    strs=ReadPopupEdits(f)
    for i=1:3
      [hkl,bErr]=strvec2reals(strs(2*i-1),3)
      if( bErr ) then
        str="< syntax error"
        bOK=%f
      else
        pixwav=HKL2PixWav(hkl')
        if( (pixwav == []) | (pixwav(3) < 0) ) then
          str="< impossible hkl"
          bOK=%f
        else
          str=msprintf("%.f %.f %.2f",pixwav)
          bOK=%t
        end
      end
      f.children(17-3*i).string=str
    end
  end
//
// Try to read XYs and HKLs (complain and return if invalid)
  vals=ReadPopupEditsNumber(f,18)
  if(vals == []) then return; end
  vals=matrix(vals,6,3)'
//
// Calculate new UB and cell from the 3 spots
  LatticeInfo(1).ub = CalcYRot(ImageInfo.phi) * ...
                          PixWav2Hvec(vals(:,4:6))' * inv(vals(:,1:3)')
  LatticeInfo(1).cell=UB2UnitCell(LatticeInfo(1).ub)
//
// Change to triclinic with no centering, and set ilevel to not refined
  LatticeInfo(1).level=2
  LatticeInfo(1).ilatt=1
  LatticeInfo(1).icen=0
//
// Update the final cell dimensions and volume on popup box
  f.children(5).String=msprintf("%.3f %.3f %.3f %.2f %.2f %.2f",LatticeInfo(1).cell)
  f.children(4).String=msprintf("( Cell Volume = %.f )",1/abs(det(LatticeInfo(1).ub)))
//
// Generate and draw calculated spots
  DrawGenHKLs(0)
//
// If "Accept", close popup
  if(g.tag == "accept") then
    close(f)
    return
  end
//
// Must be "Test", so revert to original UB
  LatticeInfo(1)=Linfo
//
endfunction


function InputCell_Popup(bIndex)
//
// If bIndex is set, index cell after "accept" pressed
//
// Create a popup window and strip off the menus and toolbar
  f=CreateEmptyPopup([300 235],"Enter cell parameters:")
// Convert cell_dims to a string array (on error use defaults)
  values=["10.0";"10.0";"10.0";"90.0";"90.0";"90.0"]
// Create the 6 text & 6 edit boxes for cell dimensions
  names=["        a","        b","        c","  alpha","    beta","gamma"]
  for i=1:6
    iy=220-28*i
    AddTextUI(f,[120 iy 45 25],names(i),13)
    if(i <= 3) then
      str=msprintf("%.3f",LatticeInfo(1).cell(i))
    else
      str=msprintf("%.2f",LatticeInfo(1).cell(i))
    end
    AddEditUI(f,[170 iy-3 60 25],str,13)
  end
//
// Listbox for centering
  icen=LatticeInfo(1).icen+1
  AddListUI(f,[255 30 35 190],"| P| A| B| C|  I| F| R",18,icen)
//
// Pushbuttons for accept and cancel
  h1=AddButtonUI(f,[60 5 70 30],"Accept",18)
  h2=AddButtonUI(f,[160 5 70 30],"Cancel",18)
// Add callbacks and tags for the pushbuttons
  AddUICallback(h1,"InputCell_CB","accept")
  AddUICallback(h2,"InputCell_CB","cancel")
//
// *** This must be last created graphics object ***
// Listbox for lattice type
  str="Triclinic|Monoclinic|Orthorhombic|Tetragonal|" + ...
      "Cubic|Trigonal (rhom)|Trigonal (hex)|Hexagonal"
  h=AddListUI(f,[10 50 95 165],str,13,1)
  f.children(1).value=LatticeInfo(1).ilatt
  AddUICallback(h,"InputCell_VisCB","")
//
// Call callback to force correct cell dims for lattice type
  InputCell_VisCB()
//
// Pass bIndex to callback
  f.user_data=bIndex
//
endfunction

function InputCell_VisCB()
// Set popup visibilities for selected lattice type
// Get popup figure and value of lattice type ListUI
  f=gcf()
  ilatt=f.children(1).value
//
// Put in vis() which cell dims are on or off
  if(ilatt == 1) then
    vis=["on" "on"  "on"  "on"  "on"  "on"]
  elseif(ilatt == 2) then
    vis=["on" "on"  "on"  "off" "on"  "off"]
  elseif(ilatt == 3) then
    vis=["on" "on"  "on"  "off" "off" "off"]
  elseif(ilatt == 4) then
    vis=["on" "off" "on"  "off" "off" "off"]
  elseif(ilatt == 5) then
    vis=["on" "off" "off" "off" "off" "off"]
  elseif(ilatt == 6) then
    vis=["on" "off" "off" "on"  "off" "off"]
  else
    vis=["on" "off" "on"  "off" "off" "off"]
  end
// Turn on/off the corresponding text and edit widgets

  for i=1:6
    f.children(18-2*i).visible=vis(i)
    f.children(17-2*i).visible=vis(i)
  end
//
endfunction

function InputCell_CB()
global LatticeInfo
// Get graphics object and parent figure
  g=gcbo()
  f=g.parent
//
// Get bIndex value
  bIndex=f.user_data
//
// If "cancel", close popup
  if(g.tag == "cancel") then
    close(f)
    return
  end
//
// Read cell parameters, return if invalid
  [bOK,ilatt,icen,cell_dims]=ReadCellPopupValues(f)
  if (~bOK) then return; end
// Close popup
  close(f)
//
// Change the cell/type/centering in memory
  LatticeInfo(1).ilatt=ilatt
  LatticeInfo(1).icen=icen
  LatticeInfo(1).cell=ConvertCellDims(cell_dims,ilatt)
//
// Set level to "cell given" and remove any twins
  LatticeInfo(1).level=1
  LatticeInfo=LatticeInfo(1)

// Index the cell if bIndex is set, else erase the calc spots
  if( bIndex ) then
    ImageIndex_Menu(2,0)
  else
    DrawCalcSpots([])
  end
endfunction

function [bOK,ilatt,icen,cell_dims]=ReadCellPopupValues(f)
//
// Read cell dims, if invalid complain and return with bOK=%f
  cell_dims=ReadPopupEditsNumber(f,6)
  bOK=( cell_dims ~= [] )
  if( ~bOK ) then
    ErrorBox("Invalid number")
    return
  end
//
// Read the cell type and centering from the popup
  ilatt=f.children(1).value
  icen=f.children(4).value-1
endfunction


function cell_dims=ConvertCellDims(cell_dims,ilatt)
//
// Convert cell dims for the given lattice type
  if(ilatt == 2) then
    cell_dims([4,6])=[90,90]
  elseif(ilatt == 3) then
    cell_dims([4,5,6])=[90,90,90]
  elseif(ilatt == 4) then
    cell_dims([2,4:6])=[cell_dims(1),90,90,90]
  elseif(ilatt == 5) then
    cell_dims([2:6])=[cell_dims([1,1]),90,90,90]
  elseif(ilatt == 6) then
    cell_dims([2:3,5:6])=cell_dims([1,1,4,4])
  elseif(ilatt >= 7) then
    cell_dims([2,4:6])=[cell_dims(1),90,90,120]

  end

endfunction


function BatchCenter_Popup(icen)
//
// Create a popup window and strip off the menus and toolbar
  f=CreateEmptyPopup([180 235],"Select centering:")
//
// Second listbox for centering
  AddListUI(f,[70 40 35 190],"| P| A| B| C|  I| F| R",18,1)
  f.children(1).value=icen+1
//
// Pushbuttons for accept and cancel
  h1=AddButtonUI(f,[10 5 70 30],"OK",18)
  AddUICallback(h1,"BatchCenter_CB","ok")
  h2=AddButtonUI(f,[100 5 70 30],"Cancel",18)
  AddUICallback(h2,"BatchCenter_CB","cancel")
//
// Prevent window closing except by the pushbuttons
  f.closerequestfcn="InfoBox(""You must use the push buttons to exit this popup"")"
//
endfunction

function BatchCenter_CB()
global LatticeInfo OrientInfo ModulateInfo
// Get graphics object and parent figure
  g=gcbo()
  f=g.parent
//
// If "cancel", close popup
  if(g.tag == "cancel") then
    close(f)
    return
  end
//
// Read centering from popup, then close popup
  icen=f.children(3).value-1
  close(f)
//
// Write info to console and log file
  sCen=["P","A","B","C","I","F","R"]
  str=msprintf("Changing index files to %s",sCen(icen+1))
  mprintf("%s\n",str)
  WriteLogFile(str)
//
// Change centering in index files
  text=[]
  for fname=BatchFiles'
    [LatticeInfo,OrientInfo,ModulateInfo]=ReadIndexFile(fname+".idx")
    if(LatticeInfo == []) then
      text($+1)="Index file "+fname+".idx is missing"
    else
      LatticeInfo(:).icen=icen
      SaveIndexFile(fname+".idx",%f)
    end
  end
//
// Output if any files were missing
  if( text ~= [] ) then
    ErrorBox(text)
  end
//
endfunction


function ConvertLauegen_Popup()
// Create a UB using current cell and Lauegen angles from user
//
// Create a popup without the menus and toolbar
  f=CreateEmptyPopup([220 125],"UB from Lauegen angles")
//
// Text & edit boxes
  AddTextUI(f,[ 40 85  20 25],"X",14)
  AddTextUI(f,[100 85  20 25],"Y",14)
  AddTextUI(f,[160 85  20 25],"Z",14)
  AddEditUI(f,[ 20 60  50 25],"0",14)
  AddEditUI(f,[ 80 60  50 25],"0",14)
  AddEditUI(f,[140 60  50 25],"0",14)
//
// Pushbuttons to copntinue, or quit
  h1=AddButtonUI(f,[ 20 10 70 30],"OK",18)
  h2=AddButtonUI(f,[120 10 70 30],"Cancel",18)
  AddUICallback(h1,"ConvertLauegen_CB","ok")
  AddUICallback(h2,"ConvertLauegen_CB","cancel")
//
endfunction

function ConvertLauegen_CB()
global LatticeInfo
//

// Get graphics object, parent figure and user_data
  g=gcbo()
  f=g.parent
//
// If "ok", read Lauegen angles and calculate UB
  if(g.tag == "ok") then
    vals=ReadPopupEditsNumber(f,3)
    if(vals == []) then return; end
    LatticeInfo(1).ub=CalcUBfromLauegen(vals(1),vals(2),vals(3))
    LatticeInfo(1).level=2
    InfoBox("New UB created")
  end
//
  close(f)
//
endfunction


//////////////////////////////////////////////////////
// ------------ Popups to process spots ----------- //
//////////////////////////////////////////////////////

function ManualOrient_Popup()
//
// Create a popup window to manually rotate and shift lattice
  f=CreateEmptyPopup([250 280],"Manual Lattice Orientation")
//
  AddTextUI(f,[10 240 70 25],"X Y Pivot",13)
  AddEditUI(f,[85 240 60 25],"0",13)
  AddEditUI(f,[150 240 60 25],"0",13)
  f.children(1:2).enable="off"
  h=AddCheckUI(f,[220 240],%f)
  AddUICallback(h,"ManualOrient_CB","pivot")
//
  AddTextUI(f,[100 205 180 25],"Value",13)
  AddTextUI(f,[195 205 180 25],"Step",13)
//
  for iy=1:5
    if(iy < 4) then
      str=char(ascii("X")+iy-1)+" rotate"
      y=210-30*iy
      step=1
    else
      str=char(ascii("X")+iy-4)+" shift"
      y=200-30*iy
      step=10
    end
    AddTextUI(f,[5 y 50 25],str,13)
    AddEditUI(f,[90 y 50 25],"0",13)
    h1=AddButtonUI(f,[66 y+1 20 23],"<",15)
    h2=AddButtonUI(f,[145 y+1 20 23],">",15)
    AddEditUI(f,[190 y 40 25],string(step),13)
    AddUICallback(h1,"ManualOrient_CB",string(-iy))
    AddUICallback(h2,"ManualOrient_CB",string(+iy))
  end
//
// Pushbuttons for test and cancel
  h1=AddButtonUI(f,[10 10 70 30],"Accept",18)
  h2=AddButtonUI(f,[90 10 70 30],"Test",18)
  h3=AddButtonUI(f,[170 10 70 30],"Cancel",18)
  AddUICallback(h1,"ManualOrient_CB","accept")
  AddUICallback(h2,"ManualOrient_CB","test")
  AddUICallback(h3,"ManualOrient_CB","cancel")
//
// Prevent window closing except by the pushbuttons
  f.closerequestfcn="InfoBox(""You must use the push buttons to exit this popup"")"
//
// Copy original pixel center and UB to CB
  f.user_data=list(OrientInfo.pixcen,LatticeInfo(1).ub)
//
endfunction


function ManualOrient_CB()
global OrientInfo LatticeInfo
//
// Get graphics object and parent figure
    g=gcbo()
    f=g.parent
//
// If "cancel", restore orginal values and close popup
  if(g.tag == "cancel") then
    OrientInfo.pixcen=f.user_data(1)
    LatticeInfo(1).ub=f.user_data(2)
    DrawGenHKLs(0)
    close(f)
    return
  end
//
// Pivot check-box has changed
  if(g.tag == "pivot") then
    ixyz=[4:23]
    ipiv=[32:34]
    ixrot=28
    if( ReadPopupChecks(f) )
      f.children(ixyz).Enable="off"
      f.children(ipiv).enable="on"
      f.children(ixrot).String="Rotate"
    else
      f.children(ixyz).Enable="on"
      f.children(ipiv).enable="off"
      f.children(ixrot).String="X rotate"
    end
  end
//
// Get edit boxes strings
  buff=ReadPopupEdits(f)
// Check for weird SCILAB bug
  if(buff($) == "") then
    ErrorBox(["A SCILAB bug has occurred that I can''t fix!"; ...
              "The popup will close, then you can try again."])
    close(f)
    return
  end
//
// Decode vector values
  [fValue,bErr]=strvec2reals(buff(3:2:$),5)
  if( bErr ) then
    ErrorBox("Invalid orientation value")
    return
  end
// Decode vector steps
  [fStep,bErr]=strvec2reals(buff(4:2:$),5)
  if( bErr ) then
    ErrorBox("Invalid orientation step")
    return
  end
//
// If an arrow button, step the appropriate vector value
  itag=strtod(g.tag)
  if( ~isnan(itag) ) then
    ival=abs(itag)
    fValue(ival)=fValue(ival)+fStep(ival)*(itag/ival)
    f.children(32-5*ival).String=msprintf("%.1f",fValue(ival))
  end
//
// If pivot check-box is true, create rotation matrix for pivot
  if( ReadPopupChecks(f) )
// Read x & y values
    [xy_pivot,bErr]=strvec2reals(buff(1:2),2)
    if( bErr ) then
      ErrorBox("Invalid pivot value")
      return
    end
// Get hkl for pivot point
    hkl=PixWav2HKL([xy_pivot,1])
// Calculate h-vector (in PHI=0 coords)
    hvec=f.user_data(2) * hkl'
    hvec=hvec/norm(hvec)
// Calculate perpendiculars (in PHI=0 coords)
    hperp1=cross(hvec,[1;0;0])
    if(norm(hperp1) < 0.1) then hperp1=cross(hvec,[0;1;0]); end
    hperp1=hperp1/norm(hperp1)
    hperp2=cross(hvec,hperp1)
    hperp2=hperp2/norm(hperp2)
// Calculate U matrix from hvec to X axis
    U=[ hvec, hperp1, hperp2 ]
// Calculate rotation matrix around h-vector
    rot_mx = U * CalcXRot(fValue(1)) * U'
  else
// Else, create rotation matrix from X,Y,Z angles
    rot_mx=CalcZRot(fValue(1)) * CalcYRot(fValue(2)) * CalcXRot(fValue(3))
  end
//
// Change pixcen and UB
  dxy=fValue(4:5)
  OrientInfo.pixcen=f.user_data(1)+dxy
  LatticeInfo(1).ub=rot_mx*f.user_data(2)
//
// Display new calculated spots
  DrawGenHKLs(0)
//
// If "Accept", close popup
  if(g.tag == "accept") then
    close(f)
  end
//
endfunction


function AdvIndex_Popup()
// Create a popup to run indexer in advanced mode
//
// Create a popup window and strip off the menus and toolbar
  f=CreateEmptyPopup([490 265],"Advanced Spot Indexing");
// Text & edit boxes for 
  AddTextUI(f,[ 10 235 250 20],"Matching Spots to Pairs",14);
  AddTextUI(f,[ 30 215 160 20],"No of obs. spots to match",14);
  AddEditUI(f,[200 215  40 20],"20",14);
  AddTextUI(f,[ 30 195 160 20],"Max. h*h+k*k+l*l for calc.",14);
  AddEditUI(f,[200 195  40 20],"15",14);
  AddTextUI(f,[ 30 175 160 20],"Min. pair angle",14);
  AddEditUI(f,[200 175  40 20],"15",14);
  AddTextUI(f,[ 30 155 160 20],"Max. angle mismatch",14);

  AddEditUI(f,[200 155  40 20],"2.0",14);
// Text & edit boxes for calc. spots method 2
  AddTextUI(f,[ 10 120 200 20],"Matching Pairs to Triplets",14);
  AddTextUI(f,[ 30 100 160 20],"Max. angle mismatch",14);
  AddEditUI(f,[200 100  40 20],"2.0",14);
  AddTextUI(f,[ 30  80 160 20],"Min. number of matches",14);
  AddEditUI(f,[200  80  40 20],"4",14);
  AddTextUI(f,[ 30  60 160 20],"Min. triplet dot-cross",14);
  AddEditUI(f,[200  60  40 20],"0.25",14);
// Text & edit boxes for matching obs. & calc. spots
  AddTextUI(f, [260 235 200 20],"Initial Acceptance Tests",14);
  AddTextUI(f, [280 215 160 20],"No of obs. spots to match",14);
  AddEditUI(f, [450 215  30 20],"25",14);
  AddTextUI(f, [280 195 160 20],"Max. hkl for calc. spots",14);
  AddEditUI(f, [450 195  30 20],"5",14);
//
// Text & edit boxes for parameters refined
  AddTextUI(f, [260 160 200 20],"Final Test & Solution Selection",14);
  AddTextUI(f, [280 140 160 20],"Max. hkl for calc. spots",14);
  AddEditUI(f, [450 140  30 20],"5",14);
  AddTextUI(f, [280 120 160 20],"Solution number",14);
  AddEditUI(f, [450 120  30 20],"1",14);
//
// Pushbuttons to run indexer or cancel
  h1=AddButtonUI(f,[100 10 130 30],"Run Indexer",18);
  h2=AddButtonUI(f,[270 10 130 30],"Cancel",  18);
  AddUICallback(h1,"AdvIndex_CB","accept")
  AddUICallback(h2,"AdvIndex_CB","cancel")
//
endfunction


function AdvIndex_CB()
global LatticeInfo ImageInfo OrientInfo
// Get graphics object and parent figure
  g=gcbo()
  f=g.parent
// If "Cancel" was pressed, close the popup
  if(g.tag == "cancel") then
    close(f)
    return
  end
// Decode numeric parameters, return if invalid
  vals=ReadPopupEditsNumber(f,11)
  if(vals == []) then return; end
// Copy the values to pars
  pars.nmatch          =round(vals(1))
  pars.hsq_max         =round(vals(2))
  pars.pair_ang_min    =      vals(3)
  pars.pair_dang_max   =      vals(4)
  pars.azi_dang_max    =      vals(5)
  pars.min_azi_group   =round(vals(6))
  pars.dot_cross_min   =      vals(7)
  pars.ntest           =round(vals(8))
  pars.hkl_max1        =round(vals(9))
  pars.hkl_max2        =round(vals(10))
  pars.iselect         =round(vals(11))
// Everything looks OK, so start by removing popup box
  close(f)
// Run the indexer and display all program output
  LockMenus("Indexing spots")
  [ub,pixcen]=RunIndexSpots(pars,DisplayData.found_spots,%t)
// Indexer didn"t abort, so save values
  LatticeInfo(1).level=2
  OrientInfo.pixcen=pixcen
  LatticeInfo(1).ub=CalcYRot(ImageInfo.phi)*ub
// Display the calculated spots and unlock the menus
  DrawGenHKLs(0)
  UnlockMenus()
endfunction


function RefineOrient_Popup()
// Global saves parameters and modes for callback
global SaveAdvOrientPars
//
// Reset the saved parameters and modes
  SaveAdvOrientPars=[]
  SaveAdvOrientPars.ihklgen=0
  SaveAdvOrientPars.hklgen=[0,0]
  SaveAdvOrientPars.imatch=0
  SaveAdvOrientPars.nobs=0
  SaveAdvOrientPars.max_sep=0
  SaveAdvOrientPars.irefs=[0,0,0,0,0,0]
//
// Create a popup without the menus and toolbar
  f=CreateEmptyPopup([470 305],"Orientation Refinement");
//
// Text & edit boxes for obs. spot selection
  AddTextUI(f,[ 10 275 250 20],"Observed Spot Selection",14);
  AddTextUI(f,[ 30 255 150 20],"Max. number to use",14);
  nobs=size(DisplayData.found_spots,1)
  AddEditUI(f,[200 255  30 20],string(nobs),14);
// Text & edit boxes for calc. spots method 1
  AddTextUI(f,[ 10 225 200 20],"Calculate Spots Method 1",14);
  AddTextUI(f,[ 30 205 170 20],"Max. h*h + k*k +l*l",14);
  AddEditUI(f,[200 205  30 20],"6",14);
  AddTextUI(f,[ 30 185 150 20],"Min. wavelength",14);
  AddEditUI(f,[200 185  30 20],"2.0",14);
// Text & edit boxes for calc. spots method 2
  AddTextUI(f,[ 10 155 200 20],"Calculate Spots Method 2",14);
  AddTextUI(f,[ 30 135 170 20],"Min. d-spacing",14);
  AddEditUI(f,[200 135  30 20],string(DisplaySet.d_min),14);
  AddTextUI(f,[ 30 115 150 20],"Min. wavelength",14);
  AddEditUI(f,[200 115  30 20],string(DisplaySet.wav_min),14);
// Text & edit boxes for matching obs. & calc. spots
  AddTextUI(f, [260 275 200 20],"Matching Obs. & Calc. Spots",14);
  AddTextUI(f, [280 255 170 20],"Max. separation (mm)",14);
  AddEditUI(f, [420 255  30 20],"1.3",14);
// Text & edit boxes for parameters refined
  AddTextUI(f, [260 225 200 20],"Parameters Refined",14);
  AddTextUI(f, [280 205 170 20],"Pixel centers",14);
  AddCheckUI(f,[420 205],%t);
  AddTextUI(f, [280 185 150 20],"Cell dimensions",14);
  AddCheckUI(f,[420 185],%t);
  AddTextUI(f, [280 165 150 20],"Sample orientation",14);
  AddCheckUI(f,[420 165],%t);
  AddTextUI(f, [280 145 150 20],"Pixel size & shape",14);
  AddCheckUI(f,[420 145],%t);
  AddTextUI(f, [280 125 150 20],"Sample offset",14);
  AddCheckUI(f,[420 125],%t);
  AddTextUI(f, [280 105 150 20],"Beam direction",14);
  AddCheckUI(f,[420 105],%t);
//
// Pushbuttons to run the algorithms
  h1=AddButtonUI(f,[ 20  60  90 30],"Run Calc. 1",16);
  h2=AddButtonUI(f,[130  60  90 30],"Run Calc. 2",16);
  h3=AddButtonUI(f,[240  60  90 30],"Run Match"  ,16);
  h4=AddButtonUI(f,[350  60  90 30],"Run Refine" ,16);
  AddUICallback(h1,"RefineOrient_CB","calc1")
  AddUICallback(h2,"RefineOrient_CB","calc2")
  AddUICallback(h3,"RefineOrient_CB","match")
  AddUICallback(h4,"RefineOrient_CB","refine")
//
// Pushbuttons to accept or ignore the results
  h1=AddButtonUI(f,[ 80 10 130 30],"Accept & Exit",18);
  h2=AddButtonUI(f,[250 10 130 30],"Reject & Exit",  18);
  AddUICallback(h1,"RefineOrient_CB","accept")
  AddUICallback(h2,"RefineOrient_CB","cancel")
//
// Prevent window closing except by the pushbuttons
  f.closerequestfcn="InfoBox(""You must use the push buttons to exit this popup"")"
//
// Set nobs to use all obs. spots
//
// Set d_min & wav_min for Method 2 to display values
//
// Save the orientation data and calculated spots
  f.user_data=list(OrientInfo,LatticeInfo,MainMarks(2))
//
endfunction


function RefineOrient_CB()
global DisplayData OrientInfo LatticeInfo MainMarks
// Global saves the previous parameters and modes
global SaveAdvOrientPars
//
// Get the current graphic object and its parent figure
  g=gcbo();
  f=g.parent;
//
// If "accept", close popup
  if(g.tag=="accept" ) then
    close(f)
    return
  end
//
// If "cancel", restore orientation data and calculated spots, then close popup
  if( g.tag=="cancel" ) then
    OrientInfo=f.user_data(1)
    LatticeInfo=f.user_data(2)
    MainMarks(2)=f.user_data(3)
    close(f)
    return
  end
//
// Restore previous parameters and modes for this callback
  pars=SaveAdvOrientPars
//
// Read the popup and checkbox values, return if invalid
  vals=ReadPopupEditsNumber(f,6)
  if(vals == []) then return; end
  checks=ReadPopupChecks(f)
//
// Ensure nobs is within limits and update the popup box
  nobs=size(DisplayData.found_spots,1)
  nobs=max(1,min(nobs, vals(1)))
  f.children(33).string=string(nobs)
//
// Pressed "Calc 1", calculate hkls (method 1)
  imode=0
  if( g.tag=="calc1" ) then
    imode=4
    pars.ihklgen=1
    pars.imatch=0
    pars.hklgen=[vals(2),vals(3)]
// Pressed "Calc 2", calculate hkls (method 2)
  elseif( g.tag=="calc2" ) then
    imode=4
    pars.ihklgen=2
    pars.hklgen=[vals(4),vals(5)]
    pars.imatch=0
// Pressed "match", calculate hkl and spot match
  elseif( g.tag=="match" ) then
    if(pars.ihklgen == 0) then
      ErrorBox("Must generate hkls before matching")
      return
    end
    imode=3
    pars.imatch=1
    pars.nobs=nobs
    pars.max_sep=vals(6)/abs(OrientInfo.pixsize(1))
// Pressed "match", calculate hkl, spot match then refinement
  elseif( g.tag=="refine" ) then
// Check that we have matched hkls
    if(pars.ihklgen == 0) then
      ErrorBox("Must generate hkls before matching")
      return
    elseif(pars.imatch == 0) then
      ErrorBox("Must match spots before refining")
      return
    end
// Load params and modes for refinement
    imode=2
    pars.irefs=double(checks')
  else
// Sanity check
    BugBox("invalid tag="+g.tag)
  end
//
// Run orient_spots in manual mode
  [Linfo,Oinfo,text]=RunOrientSpots(LatticeInfo(1),OrientInfo,imode,pars)
// Output results, except in refine mode
  if(imode ~= 2) then
    InfoBox(text)
  end
//
// If refine mode, update orientation/lattice and draw spots
  if(imode == 2) then
    Linfo.level=3         // now fully oriented (we presume)
    LatticeInfo(1)=Linfo
    OrientInfo=Oinfo
    DrawGenHKLs(0)
  end
//
// Save current parameters and modes for next call of callback
  SaveAdvOrientPars=pars
//
endfunction


function f=OrientTwins_Popup(bIndex)
//
// Copy LatticeInfo to Linfo
  Linfo=LatticeInfo
//
// Create a popup to index or refine the split spot
// positions relative to each other
//
// Create a popup window without the menus and toolbar
  if( bIndex ) then
    f=CreateEmptyPopup([260 370-80],"Index Split Spots")
// If indexing, copy twins from default with slightly rotated UBs
    Linfo(2:3)=Linfo(1)
    rot=eye(3,3)+[0,-1,-2;1,0,-3;2,3,0]*3e-4
    rot=(rot+inv(rot'))/2
    Linfo(2).ub=rot*Linfo(2).ub
    Linfo(3).ub=inv(rot)*Linfo(2).ub
// Output initial cell parameters and the lattice rotation
    AddTextUI(f,[ 10 345-80 160 25],"Initial Cell Dimensions",13)
    str=msprintf("%.3f %.3f %.3f %.2f %.2f %.2f",Linfo(2).cell)
    AddTextUI(f,[ 15 332-80 210 20],str,11)
  else
    f=CreateEmptyPopup([260 370],"Refine Split Spots")
// Output initial cell parameters and the lattice rotation
    AddTextUI(f,[ 10 345 160 25],"Initial Cell Dimensions",13)
    str=msprintf("%.3f %.3f %.3f %.2f %.2f %.2f",Linfo(2).cell)
    AddTextUI(f,[ 15 332 210 20],str,11)
    AddTextUI(f,[230 332 50 20],"#1",11)
    str=msprintf("%.3f %.3f %.3f %.2f %.2f %.2f",Linfo(3).cell)
    AddTextUI(f,[ 15 315 210 20],str,11)
    AddTextUI(f,[230 315 50 20],"#2",11)
// Calculate and output rotations between lattices
    AddTextUI(f,[ 10 293 240 20],"Rotation between Lattices",13)
    [hkl1,hkl2,ang]=CalcSplitRotation(Linfo(2).ub,Linfo(3).ub)
    str=msprintf("%.2f degree about hkl = %.1f,%.1f,%.1f",ang,6*hkl1')
    AddTextUI(f,[ 30 277 200 20],str,11)
    AddTextUI(f,[230 277 50 20],"#1",11)
    str=msprintf("%.2f degree about hkl = %.1f,%.1f,%.1f",ang,6*hkl2')
    AddTextUI(f,[ 30 260 200 20],str,11)
    AddTextUI(f,[230 260 50 20],"#2",11)
  end
//
// Spot search parameters to use
  AddTextUI(f,[ 10 225 160 20],"Spot Search Parameters",13)
  AddTextUI(f,[ 30 205 70 20],"hkl max.",13)
  AddEditUI(f,[ 90 205 20 20],"5",13)
  AddTextUI(f,[150 205 70 20],"wav. min.",13)
  AddEditUI(f,[210 205 30 20],"1.8",13)
  AddTextUI(f,[ 30 185 130 20],"split max. (pixels)",13)
  AddEditUI(f,[145 185 25 20],"20",13)
//
// Are cell parameters fixed?
  AddTextUI(f,[ 10 160 130 20],"Fixed Cell Parameters",13)
  AddCheckUI(f,[150 155],%t)
//
// Output refined cell dimensions and lattice rotations as unknown
  AddTextUI(f,[10 120 160 25],"Final Cell Dimensions",13)
  AddTextUI(f,[15 107 210 20],"    unknown",11)
  AddTextUI(f,[230 107 50 20],"#1",11)
  AddTextUI(f,[15  90 210 20],"    unknown",11)
  AddTextUI(f,[230 90 50 20],"#2",11)
  AddTextUI(f,[10  70 240 20],"Rotation between Lattices",13)
  AddTextUI(f,[30  55 200 20],"unknown",11)
  AddTextUI(f,[230 55 50 20],"#1",11)
  AddTextUI(f,[30  38 200 20],"unknown",11)
  AddTextUI(f,[230 38 50 20],"#2",11)
//
// Pushbuttons for accept/test/cancel, and add callbacks
  h1=AddButtonUI(f,[ 10  7 70 30],"Accept",18)
  h2=AddButtonUI(f,[100  7 60 30],"Test",18)
  h3=AddButtonUI(f,[180  7 70 30],"Cancel",18)
  AddUICallback(h1,"OrientTwins_CB","accept")
// "test" may do indexing or refining
  if(bIndex) then
    AddUICallback(h2,"OrientTwins_CB","index")
  else
    AddUICallback(h2,"OrientTwins_CB","refine")
  end
  AddUICallback(h3,"OrientTwins_CB","cancel")
//
// Disable "Accept" until an index is performed
  f.children(3).Enable="off"
//
// Window closing done in callback
  f.closerequestfcn=msprintf("scf(%d); OrientTwins_CB()",f.figure_id)
//
// Store initial lattice and calculated spot marks for callback
  f.user_data=list(Linfo)
  for i=2:7
    f.user_data(i)=GetMarksXY(i)
  end
  f.user_data(8)=bIndex
//
// Remove calculated (and satellite) spots
  DrawCalcSpots([])
  DrawCalcModSpots([])
//
endfunction


function OrientTwins_CB()
global LatticeInfo
// Private globals to store last refinement results
global OrientTwins_Linfo OrientTwins_sOut
//
// Get graphics object
  g=gcbo()
//
// Get parent figure in a way that "closerequestfcn()" mimics
// a callback with "cancel"
  if(g.children == []) then
    f=g.parent
    tag=g.tag
  else
    f=g
    tag="cancel"
  end
//
// If "cancel" pressed, redraw calc. spots and close popup
  if(tag == "cancel") then
    for i=2:7
      DrawMarksNow(i,f.user_data(i))
    end
    close(f)
    return
  end
//
// If "accept", load refinement results, update log file, close popup
  if(tag == "accept") then
    LatticeInfo=OrientTwins_Linfo
    WriteLogFile("Manual Orient Twins for "+ImageInfo.basename+".tif")
    WriteLogFile( "   "+OrientTwins_sOut )
    WriteLogFile( "   "+f.children(12).String )
    WriteLogFile( "   "+f.children(10).String )
    WriteLogFile( "   "+f.children( 7).String )
    WriteLogFile( "   "+f.children( 5).String )
    close(f)
    return
  end
//
// Get initial lattice values and bIndex from user_data[]
  Linfo0=f.user_data(1)
  bIndex=f.user_data(8)
//
// Get checkbox value
  bFitCell= ~ReadPopupChecks(f)
//
// Get search parameter values, or exit if invalid
  vals=ReadPopupEditsNumber(f,3)
  if(vals == []) then return; end
  hkl_max=vals(1)
  wav_min=vals(2)
  dxy_max=vals(3)
//
// Orient spot pairs
  LockMenus("Orienting Split Spots")
  sFile=ImageInfo.basename+".tif"
  niter=4
  if( bIndex ) then niter=2; end
  [Linfo,sOut]=RunOrientTwins(Linfo0,OrientInfo, sFile, ...
                      hkl_max,wav_min,dxy_max,niter,bFitCell)
// Save refined values
  OrientTwins_Linfo=Linfo
  OrientTwins_sOut=sOut
//
// Force unit cell volumes to be conserved
  mult=( det(Linfo0(2).ub)/det(Linfo(2).ub) )^(1.0/3.0)
  Linfo(2).ub=Linfo(2).ub * mult
  Linfo(2).cell(1:3)=Linfo(2).cell(1:3) / mult
  mult=( det(Linfo0(3).ub)/det(Linfo(3).ub) )^(1.0/3.0)
  Linfo(3).ub=Linfo(3).ub * mult
  Linfo(3).cell(1:3)=Linfo(3).cell(1:3) / mult
//
// Output refined cell values
  f.children(12).String=msprintf("%.3f %.3f %.3f %.2f %.2f %.2f",Linfo(2).cell)
  f.children(10).String=msprintf("%.3f %.3f %.3f %.2f %.2f %.2f",Linfo(3).cell)
//
// Output lattice rotations for both lattices (in case hkls are equiv but not ident)
  [hkl1,hkl2,ang]=CalcSplitRotation(Linfo(2).ub,Linfo(3).ub)
  f.children(7).String=msprintf("%.2f degree about hkl = %.1f,%.1f,%.1f",ang,6*hkl1')
  f.children(5).String=msprintf("%.2f degree about hkl = %.1f,%.1f,%.1f",ang,6*hkl2')
//
// Set twins" levels to fully oriented
  Linfo(2:3).level=3
//
// If cell is refined, change twin lattices to triclinic
  if( bFitCell ) then
    Linfo(2:3).ilatt=1
  end
//
// Generate and draw calculated twin spots (also unlocks the menus)
  Lsave=LatticeInfo
  LatticeInfo=Linfo
  DrawGenHKLs(12)
  LatticeInfo=Lsave
//
// Show fit results in info popup
  [ndum,nspots,rms1,rms2]=msscanf(OrientTwins_sOut,"%d %f %f")
  rms1=rms1*abs(OrientInfo.pixsize(1))
  rms2=rms2*abs(OrientInfo.pixsize(1))
  text=msprintf("%d fitted spot pairs",nspots)
  text(2:3)=msprintf("      Twin %d  diff (rms) = %.2f mm\n",[1,rms1;2,rms2])
  InfoTitleBox(text,"Orientation Results")
// Turn on "accept" since "test" has been performed
  f.children(3).Enable="on"
//
endfunction


function ArgonneBoxes_Popup(bOldMode,bTwinMode,bModulateMode)
//
// Reset list of hints
  SetupHintUI()
//
// Create mode dependent part at the top of the popup
  if ( bOldMode ) then
    f=ArgBoxOld_Popup()
  else
    f=ArgBoxNew_Popup()
  end
//
// Add items common to both modes
  AddTextUI(f,[ 10 106 230 25],"Min. d-spacing for hkls",13)
  AddEditUI(f,[250 103  40 25],string(ArgboxSet.d_min),13)
  AddHintUI(f,[295 103],"Check values using ""Intensity vs d-spacing""")
  AddTextUI(f,[ 10  78 230 25],"Min. wavelength for hkls",13)
  AddEditUI(f,[250  75  40 25],string(ArgboxSet.wav_min),13)
  AddTextUI(f,[ 10  50 230 25],"Max. wavelength for hkls",13)
  AddEditUI(f,[250  47  40 25],string(ArgboxSet.wav_max),13)
//
// Pushbuttons for accept and cancel
  h1=AddButtonUI(f,[20 5 70 30]," Run ",18)
  h2=AddButtonUI(f,[115 5 80 30],"Defaults",18)
  h3=AddButtonUI(f,[220 5 70 30],"Cancel",18)
// Add callbacks and tags for the pushbuttons
  AddUICallback(h1,"ArgonneBoxes_CB","accept")
  AddUICallback(h2,"ArgonneBoxes_CB","defaults")
  AddUICallback(h3,"ArgonneBoxes_CB","cancel")
//
// Pass to callback the modes in use
  f.user_data=[bOldMode,bTwinMode,bModulateMode]
//
endfunction


function f=ArgBoxNew_Popup()
//
// Create the popup for the new-mode parameters
  f=CreateEmptyPopup([315 398],"Enter spot integration parameters:")
//
  AddTextUI(f,[ 10 365 230 25],"Target peak fraction",13)
  AddEditUI(f,[250 363  40 25],string(ArgboxSet.pfrac_target),13)
  AddHintUI(f,[295 363],"Decrease for weak spots, increase for poor spot shape")
  AddTextUI(f,[ 10 337 230 25],"Area multiplier for peak ellipse",13)
  AddEditUI(f,[250 335  40 25],string(ArgboxSet.peak_mult),13)
  AddHintUI(f,[295 335],"Maybe, increase if ""Target peak fraction"" is decreased")
  AddTextUI(f,[ 10 309 230 25],"Peak fraction uncertainty (0-100)",13)
  AddEditUI(f,[250 307  40 25],string(ArgboxSet.pfrac_uncertain),13)
  AddHintUI(f,[295 307],"Increase for poor spot shape")
  AddTextUI(f,[ 10 281 230 25],"Neighbour overlap tolerance (0-100)",13)
  AddEditUI(f,[250 279  40 25],string(ArgboxSet.neigh_toler),13)
  AddHintUI(f,[295 279],"Increase to reject fewer overlapped spots")
  AddTextUI(f,[ 10 253 230 25],"Model ellipses in middle zone",13)
  AddEditUI(f,[250 251  40 25],string(ArgboxSet.model_center),13)
  AddHintUI(f,[295 251],"Maybe, increase for good spots")
  AddTextUI(f,[ 10 225 230 25],"Model ellipses in outer zones",13)
  AddEditUI(f,[250 223  40 25],string(ArgboxSet.model_outer),13)
  AddHintUI(f,[295 223],"Maybe, increase for good spots")
  AddTextUI(f,[ 10 197 230 25],"Cutoff level of contour ellipse",13)
  AddEditUI(f,[250 195  40 25],string(ArgboxSet.cont_level),13)
  AddHintUI(f,[295 195],"Increase for weak spots or poor shape (Try this one first)")
  AddTextUI(f,[ 10 169 230 25],"Min. fill fraction for contour",13)
  AddEditUI(f,[250 167  40 25],string(ArgboxSet.min_fill),13)
  AddHintUI(f,[295 167],"Decrease for weak spots (Try this one second)")
//
  AddTextUI(f,[ 20 140 165 25],"Recenter spots",13)
  AddCheckUI(f,[120 140],double(ArgboxSet.recenter))
  AddHintUI(f,[143 140],"May help for weak spots with good shape")
//
endfunction


function f=ArgBoxOld_Popup()
//
// Create the popup for the old-mode parameters
  f=CreateEmptyPopup([315 345],"Enter spot integration parameters:")
//
  AddTextUI(f,[ 10 312 230 25],"Peak/Backg. cutoff for model peaks",13)
  AddEditUI(f,[250 309  40 25],string(ArgboxSet.model_cutoff),13)
  AddHintUI(f,[295 309],"Decrease for weak spots")
  AddTextUI(f,[ 10 284 230 25],"Peak/Backg. cutoff for strong peaks",13)
  AddEditUI(f,[250 281  40 25],string(ArgboxSet.strong_cutoff),13)
  AddHintUI(f,[295 281],"Maybe, decrease for weak spots")
  AddTextUI(f,[ 10 256 230 25],"Cutoff level for contour ellipse",13)
  AddEditUI(f,[250 253  40 25],string(ArgboxSet.cont_level),13)
  AddHintUI(f,[295 253],"Increase for weak spots (Try this one first)")
  AddTextUI(f,[ 10 228 230 25],"Area multiplier for core ellipse",13)
  AddEditUI(f,[250 225  40 25],string(ArgboxSet.core_mult),13)
  AddHintUI(f,[295 225],"Decrease for weak spots, increase for poor shape")
  AddTextUI(f,[ 10 200 230 25],"Area multiplier for peak ellipse",13)
  AddEditUI(f,[250 197  40 25],string(ArgboxSet.peak_mult),13)
  AddHintUI(f,[295 197],"Decrease for weak spots, increase for poor shape")
  AddTextUI(f,[ 10 172 230 25],"Neighbour overlap tolerance (0-100)",13)
  AddEditUI(f,[250 169  40 25],string(ArgboxSet.neigh_toler),13)
  AddHintUI(f,[295 169],"Increase to reject fewer overlapped spots")
  AddTextUI(f,[ 10 144 230 25],"Min. fill fraction for contour",13)
  AddEditUI(f,[250 141  40 25],string(ArgboxSet.min_fill),13)
  AddHintUI(f,[295 141],"Decrease for weak spots or poor shape (Try this one second)")
//
endfunction


function ArgonneBoxes_CB()
global ArgboxSet
// Get graphics object and parent figure
  g=gcbo()
  f=g.parent
// Get the user_data information
  bOldMode=f.user_data(1)
  bTwinMode=f.user_data(2)
  bModulateMode=f.user_data(3)
//
// If "cancel", close the popup
  if(g.tag == "cancel") then
    close(f)
    return
  end
//
// If "defaults", reset the edit widgets
  if(g.tag == "defaults") then
    str=ReadPopupEdits(f)
    if ( bOldMode ) then
      str(1:7)=["2","2","0.1","1","4","10","0.8"]
    else
      str(1:8)=["0.8","4","25","10","20","10","0.1","0.8"]
    end
    SetPopupEdits(f,str)
    return
  end
//
// Decode the numeric parameters for old- & new-modes
  if ( bOldMode ) then
// Read parameters, return if invalid
    pars=ReadPopupEditsNumber(f,10)
    if(pars == []) then return; end
// Save the values to ArgboxSet and the settings file
    ArgboxSet.model_cutoff =pars(1)
    ArgboxSet.strong_cutoff=pars(2)
    ArgboxSet.cont_level   =pars(3)
    ArgboxSet.core_mult    =pars(4)
    ArgboxSet.peak_mult    =pars(5)
    ArgboxSet.neigh_toler  =pars(6)
    ArgboxSet.min_fill     =pars(7)
    ArgboxSet.d_min        =pars(8)
    ArgboxSet.wav_min      =pars(9)
    ArgboxSet.wav_max      =pars(10)
  else
// For new-mode
// Read parameters, return if invalid
    pars=ReadPopupEditsNumber(f,11)
    if(pars == []) then return; end
// Save the values to ArgboxSet and the settings file
    ArgboxSet.pfrac_target   =pars(1)
    ArgboxSet.peak_mult      =pars(2)
    ArgboxSet.pfrac_uncertain=pars(3)
    ArgboxSet.neigh_toler    =pars(4)
    ArgboxSet.model_center   =pars(5)
    ArgboxSet.model_outer    =pars(6)
    ArgboxSet.cont_level     =pars(7)
    ArgboxSet.min_fill       =pars(8)
    ArgboxSet.d_min          =pars(9)
    ArgboxSet.wav_min        =pars(10)
    ArgboxSet.wav_max        =pars(11)
  end
//
// Read check box
  pars=ReadPopupChecks(f)
  if ( bOldMode ) then
    ArgboxSet.recenter=%f
  else
    ArgboxSet.recenter=pars(1)
  end
//
// Remove the popup
  close(f)
// Save the argbox settings
  SaveSettingsFile()
//
// Start argonne_boxes on batch files
  BatchArgonneBoxes(bOldMode,bTwinMode,bModulateMode)
//
endfunction


/////////////////////////////////////////////////////////
// --------- Popups to change images & marks --------- //
/////////////////////////////////////////////////////////

function Brightness_Popup()
//
// Create a popup without the menus and toolbar
  f=CreateEmptyPopup([270 125],"Intensity Brightness Levels")
// Get current levels
  levs=round(MainImage.user_data)
// Text & edit boxes for max. number of obs.
  AddTextUI(f,[ 60 85 100 25],"Max. Intensity",14)
  AddEditUI(f,[160 85  50 25],string(levs(2)),14)
  AddTextUI(f,[ 60 55 100 25],"Min. Intensity",14)
  AddEditUI(f,[160 55  50 25],string(levs(1)),14)
// Pushbuttons to accept, test and cancel setting
  h1=AddButtonUI(f,[ 20 10 70 30],"Accept",18)
  h2=AddButtonUI(f,[100 10 70 30],"Test",  18)
  h3=AddButtonUI(f,[180 10 70 30],"Cancel",18)
// Add callbacks and tags for the pushbuttons
  AddUICallback(h1,"Brightness_CB","accept")
  AddUICallback(h2,"Brightness_CB","test")
  AddUICallback(h3,"Brightness_CB","cancel")
//
// Save levs for Brightness_CB() using user_data
  f.user_data=levs
//
endfunction


function Brightness_CB()
global MainImage
//
// Get graphics object, parent figure and user_data
  g=gcbo()
  f=g.parent
  levs=f.user_data
//
// If not cancel, read levs from the popup
  if(g.tag ~= "cancel") then
    vals=ReadPopupEditsNumber(f,2)
    if(vals == []) then return; end
// Enforce lev_lo <= lev_hi
// NB: levs & vals have indices swapped
    levs(2)=vals(1)
    levs(1)=min(vals(2),vals(1)-1)
    f.children(4).string=string(levs(1))
  end
//
// If not "test", close popup



  if(g.tag ~= "test") then
    close(f)
  end
//
// Update image depending on mode
  if( DisplayMode.flick_mode ) then
// Calculate all flick pixmaps for new intensity levels
    LockMenus("Creating flicker images")
    CalcFlickPixmaps(levs)
    UnlockMenus()
// Draw image from pixmap of the selected image file
    MainImage.data=squeeze( FlickData.pixmaps(FlickData.iflick,:,:) )
// Store the brightness levels used for the pixmaps
    MainImage.user_data=levs
  else
    RedrawImageNow(levs)
  end
//
endfunction


function AdvPrune_Popup()
//
// Create a popup without the menus and toolbar
  f=CreateEmptyPopup([265 240],"Advanced Pruning of Obs. Spots")
//
  AddTextUI(f,[ 20 210 250 25],"Remove any observed spots more than",13)
  AddTextUI(f,[ 20 190 250 25],"""Spot Dist"" pixels from a calculated spot",13)
//
// Text & edit boxes for d-spacing & wavelength limits
  AddTextUI(f,[125 160 40 25],"min",14)
  AddTextUI(f,[195 160 40 25],"max",14)
  AddTextUI(f,[ 30 135 70 25],"d-spacing",14)
  AddEditUI(f,[110 135 50 25],string(DisplaySet.d_min),14)
  AddEditUI(f,[180 135 55 25],string(DisplaySet.d_max),14)
  AddTextUI(f,[ 30 105 70 25],"wavelength",14)
  AddEditUI(f,[110 105 50 25],string(DisplaySet.wav_min),14)
  AddEditUI(f,[180 105 55 25],string(DisplaySet.wav_max),14)
//
  AddTextUI(f,[ 30 70 110 25],"Spot Dist. (pixels)",14)
  AddEditUI(f,[160 70 50 25],"20",14)
//
// Pushbuttons to accept, test and cancel setting
  h1=AddButtonUI(f,[ 20 10 70 30],"Accept",18)
  h2=AddButtonUI(f,[100 10 70 30],"Test",  18)
  h3=AddButtonUI(f,[180 10 70 30],"Cancel",18)
//
// Add callbacks and tags for the pushbuttons
  AddUICallback(h1,"AdvPrune_CB","accept")
  AddUICallback(h2,"AdvPrune_CB","test")
  AddUICallback(h3,"AdvPrune_CB","cancel")
//
// Prevent window closing except by the pushbuttons
  f.closerequestfcn="InfoBox(""You must use the push buttons to exit this popup"")"
//
// Save X,Y for calc. and twins marks for the callback
  f.user_data=list(GetMarksXY(2),GetMarksXY(3),GetMarksXY(4))
//
// Remove any twins marks (too confusing to implement for twins)
  DrawTwin1Spots([])
  DrawTwin2Spots([])
//
endfunction


function AdvPrune_CB()
global DisplayData
//
// Get graphics object and parent figure
  g=gcbo()
  f=g.parent
//
// If "accept" or "test", calculate matching calc & obs spots
  if(g.tag ~= "cancel") then
//
// Read parameters from the popup window
    buff=ReadPopupEdits(f)
    [vals,ierr]=evstr(buff)
    if(size(vals,"*") ~= 5) then
      ErrorBox("Invalid values")
      return
    end
//
// Give up if too many calculated spots
    nguess=GuessCalcSpotsNumber2(LatticeInfo(1), vals(1),vals(3),vals(4))
    if(nguess > 10000) then
      ErrorBox("Too many calculated spots")
      return
    end
//
// Generate hkl list
    hkllist=RunGenHKLs(LatticeInfo(1),vals(3),vals(4),vals(1),vals(2),%f)
//
// Find obs and calc spots that match within vals(5) pixels
    match=MatchXY(DisplayData.found_spots(:,1:2),hkllist(:,5:6), vals(5))
//
  end
// 
// If "test", output summary to popup, draw spots, and return
  if(g.tag == "test") then
    nmatch=size(match,1)
    ncalc=size(hkllist,1)
    nobs=size(DisplayData.found_spots,1)
    f.children(5).String=msprintf("%d of %d observed spots within",nmatch,nobs)
    f.children(4).String=msprintf("%d pixels of %d calculated spots",vals(5),ncalc)
//
// Draw the matching obs spots and all calc spots
    DrawMarksNow(1,DisplayData.found_spots(match(:,1),1:2))
    DrawCalcSpots(hkllist(:,5:6))
//
    return
  end
//
//
// If "accept", reduce obs spots to those matched
  if(g.tag == "accept") then
    DisplayData.found_spots=DisplayData.found_spots(match(:,1),:)
  end
//
// Draw the observed spots
  DrawObsSpots(DisplayData.found_spots(:,1:2))
//
// Restore the original calc & twins spots
  DrawCalcSpots(f.user_data(1))
  DrawTwin1Spots(f.user_data(2))
  DrawTwin2Spots(f.user_data(3))
//
  close(f)
  return
//
endfunction


function FindHKLSpot_Popup()
//
// Create a popup without the menus and toolbar
  f=CreateEmptyPopup([250 150],"Find HKL Spot")
//
  AddTextUI(f,[10 110 30 25],"HKL:",13)
  AddEditUI(f,[50 110 100 25],"0 0 0",13)
//
  AddTextUI(f,[40 80 180 25],"",13)
  AddTextUI(f,[40 55 180 25],"",13)
//
// Pushbuttons for test and cancel
  h1=AddButtonUI(f,[30 10 70 30],"Find",18)
  h2=AddButtonUI(f,[150 10 70 30],"Close",18)
  AddUICallback(h1,"FindHKLSpot_CB","find")
  AddUICallback(h2,"FindHKLSpot_CB","close")
//
// Prevent window closing except by the pushbuttons
  f.closerequestfcn="InfoBox(""You must use the push buttons to exit this popup"")"
//
endfunction


function FindHKLSpot_CB()
//
// Get graphics object and parent figure
  g=gcbo()
  f=g.parent
//
// Clear the marked HKL
  DrawMarksNow(11,[])
//
// If "cancel", close popup
  if(g.tag == "close") then
    close(f)
    return
  end
//
// Read HKL, complain and return if invalid
  vals=ReadPopupEditsNumber(f,3)
  if(vals == []) then return; end
//
// Calculate X,Y of HKL, check if invalid, then output results
  pix_wav=HKL2PixWav(vals)
  f.children(3:4).string=""
  if(pix_wav == []) then
    f.children(3).string="Impossible Calculation"
  else
    f.children(4).string=msprintf("XY = %.f %.f",pix_wav(1:2))
    bOff=( or(pix_wav(1:2)<0) | or(pix_wav(1:2)>ImageInfo.numxy(1:2)) )
    if(pix_wav(3) < 0) then
      f.children(3).string=msprintf("For Friedel: %d %d %d",-vals(1:3))
      if( bOff ) then
        f.children(3).string=f.children(3).string+" & Off Screen"
      end
    elseif( bOff ) then
      f.children(3).string="Off Screen"
    end
    if( ~bOff ) then
      DrawMarksNow(11,pix_wav(1:2))
    end
  end
//
endfunction


function ObsSpots_Popup()
//
// Sanity check
  if( isempty(DisplayData.found_spots) ) then
//    WarnBox("No observed spots to display")
//    return
  end
//
// Create a popup without the menus and toolbar
  f=CreateEmptyPopup([250 125],"Observed Spots")
//
// Get number of observed spots, and number displayed
  nspots=size(DisplayData.found_spots,1)
  ndisp=size(MainMarks(1).user_data,1)
//
// Text & edit boxes for max. number of obs.
  str=string(nspots)+" observed spots were found"
  AddTextUI(f,[ 10 85 250 25],str,16)
//
// Text & edit boxes for max. number of obs.
  AddTextUI(f,[ 40 55 130 25],"Number to display",14)
  AddEditUI(f,[170 55  40 25],string(ndisp),14)
//
// Pushbuttons to accept, test and cancel setting
  h1=AddButtonUI(f,[ 10 10 70 30],"Accept",18)
  h2=AddButtonUI(f,[ 90 10 70 30],"Test",  18)
  h3=AddButtonUI(f,[170 10 70 30],"Cancel",18)
//
// Add callbacks and tags for the pushbuttons
  AddUICallback(h1,"ObsSpots_CB","accept")
  AddUICallback(h2,"ObsSpots_CB","test")
  AddUICallback(h3,"ObsSpots_CB","cancel")
//
// Prevent window closing except by the pushbuttons
  f.closerequestfcn="InfoBox(""You must use the push buttons to exit this popup"")"
//
// Save num spots displayed for ObsSpots_CB() using user_data
  f.user_data=ndisp
//
endfunction


function ObsSpots_CB()
// Get graphics object and parent figure
  g=gcbo()
  f=g.parent
// Get the user_data information
//
// If "cancel", revert to original ndisp
  if(g.tag == "cancel") then
    ndisp=f.user_data
  else
// Else, read ndisp, return if invalid
    ndisp=ReadPopupEditsNumber(f,1)
    if(ndisp == []) then return; end
  end
//
// Clean up number and redisplay on popup
  nspots=size(DisplayData.found_spots,1)
  ndisp=max(0,min(nspots, round(ndisp) ))
  f.children(4).string=string(ndisp)
//
// Display the spots
  DrawObsSpots(DisplayData.found_spots(1:ndisp,1:2))
//
// If not "test", close popup
  if(g.tag ~= "test") then
    close(f)
  end
endfunction


function CalcSpots_Popup()
//
// Create a popup without the menus and toolbar
  f=CreateEmptyPopup([270 145],"Calculated Spots")
//
// Text & edit boxes for d-spacing & wavelength limits
  AddTextUI(f,[125 120 40 25],"min",14)
  AddTextUI(f,[195 120 40 25],"max",14)
  AddTextUI(f,[ 30  95 70 25],"d-spacing",14)
  AddEditUI(f,[110  95 50 25],string(DisplaySet.d_min),14)
  AddEditUI(f,[180  95 55 25],string(DisplaySet.d_max),14)
  AddTextUI(f,[ 30  55 70 25],"wavelength",14)
  AddEditUI(f,[110  55 50 25],string(DisplaySet.wav_min),14)
  AddEditUI(f,[180  55 55 25],string(DisplaySet.wav_max),14)
//
// Pushbuttons to accept, test and cancel setting
  h1=AddButtonUI(f,[ 20 10 70 30],"Accept",18)
  h2=AddButtonUI(f,[100 10 70 30],"Test",  18)
  h3=AddButtonUI(f,[180 10 70 30],"Cancel",18)
//
// Add callbacks and tags for the pushbuttons
  AddUICallback(h1,"CalcSpots_CB","accept")
  AddUICallback(h2,"CalcSpots_CB","test")
  AddUICallback(h3,"CalcSpots_CB","cancel")
//
// Prevent window closing except by the pushbuttons
  f.closerequestfcn="InfoBox(""You must use the push buttons to exit this popup"")"
//
// Save the displayed calc. spots in user_data
  f.user_data=MainMarks(2).user_data
//
endfunction


function CalcSpots_CB()
global DisplaySet
// Get graphics object and parent figure
  g=gcbo()
  f=g.parent
//
// If "cancel", restore displayed spots and close popup
  if(g.tag == "cancel") then
    DrawCalcSpots(f.user_data)
    close(f)
    return
  end
//
// Read d_min & wav_min from the popup, return in invalid
  vals=ReadPopupEditsNumber(f,4)
  if(vals == []) then return; end
  d_min=vals(1)
  d_max=vals(2)
  wav_min=vals(3)
  wav_max=vals(4)
//
// Override unreasonable values
  d_min=max(0.3, d_min )
  wav_min=max(0.6, wav_min )
  f.children(8).string=string(d_min)
  f.children(5).string=string(wav_min)
  d_max=max(d_min+0.01, d_max )
  wav_max=max(wav_min+0.01, wav_max )
  f.children(7).string=string(d_max)
  f.children(4).string=string(wav_max)
//
// Save current values, and then update DisplaySet
  Dsave=DisplaySet
  DisplaySet.d_min=d_min
  DisplaySet.d_max=d_max
  DisplaySet.wav_min=wav_min


  DisplaySet.wav_max=wav_max
//
// Generate and draw calculated spots
  DrawGenHKLs(0)
//
// If accept, save setings and close popup
  if(g.tag == "accept") then
    SaveSettingsFile()
    close(f)
// ELSE, restore saved DisplaySet values
  else
    DisplaySet=Dsave
  end
endfunction


////////////////////////////////////////////////////////
// ------- Popups to select files from a list ------- //
////////////////////////////////////////////////////////

function RenameFiles_Popup(get_files)
//
// Create a popup without the menus and toolbar
  f=CreateEmptyPopup([240 285],"Rename Files")
//
// New base name and start index
  AddTextUI(f,[ 10 255 100 25],"Base file name",14)
  AddEditUI(f,[110 255 120 25],"new",14)
  AddTextUI(f,[ 10 225 90 25],"Start number",14)
  AddEditUI(f,[100 225  45 25],"1",14)
//
// How to order files
  AddTextUI(f,[ 10 200 90 25],"Sort files by:",14)
  str="File names|Phi increasing|Phi decreasing|Scan date incr.|Scan date decr."
  AddListUI(f,[80 90 100 110],str,13,1)
  f.children(1).value=1
//
// Do a rename or a copy
  AddTextUI(f,[ 60 60 120 25],"Keep original files",14)
  AddCheckUI(f,[180 60],0)
//
// Pushbuttons to accept, test and cancel setting
  h1=AddButtonUI(f,[ 20 10 60 30],"OK",18)
  h2=AddButtonUI(f,[140 10 80 30],"Cancel",18)
// Add callbacks and tags for the pushbuttons
  AddUICallback(h1,"RenameFiles_CB","rename")
  AddUICallback(h2,"RenameFiles_CB","cancel")
//
// Pass file list to CB
  f.user_data=get_files
//
endfunction


function RenameFiles_CB()
//
// Get graphics object, parent figure and user_data
  g=gcbo()
  f=g.parent
  get_files=f.user_data;
// If "cancel", close popup
  if(g.tag == "cancel") then
    close(f)
    return
  end
//
// Get new base name and starting number
  vals=ReadPopupEdits(f)
  if(vals(1) == "") then
    ErrorBox("Missing base file name")
    return
  end
// String 1 is the base name
  base_name=vals(1)
// String 2 is istart, a non-negative integer
  [iStart,bErr]=strvec2reals(vals(2),1)
  if( ~bErr ) then
    bErr= (iStart < 0) | (modulo(iStart,1) ~= 0)
  end
  if( bErr ) then
    ErrorBox("Invalid integer for start number")
    return
  end
//
// Get option to keep original files
  bKeep=ReadPopupChecks(f)
//
// Get the sort method
  iSort=f.children(5).value
//
// Finished with user input, so close popup
  close(f)
//
// Read phi and scan times from each file header (if needed)
  if(iSort > 1) then
    phis=[];
    times=[];
    for file_name=get_files
      [pname,bname,ename]=fileparts(file_name);
      Info=RunImageInfo(pname+bname,%f);
      phis=[phis,Info.phi];
// Total seconds as if 100 sec. in min., 100 min. in hour, ...
      isec=sum( strtod(tokens(Info.date,[":"," ","/"])) .* (10^[0;2;4;6;8;10]) )
      times=[times,isec]
    end
  end
//
// Get the sort index "k"
  if(iSort == 1) then
    [v,k]=gsort(get_files,'g','i')
  elseif(iSort == 2) then
    [v,k]=gsort(phis,'g','i')
  elseif(iSort == 3) then
    [v,k]=gsort(phis,'g','d')
  elseif(iSort == 4) then
    [v,k]=gsort(times,'g','i')
  elseif(iSort == 5) then
    [v,k]=gsort(times,'g','d')
  else
    BugBox("Impossible iSort")
  end
//
// Copy or rename files in the sorted order
  for file_name=get_files(k)
// Create new name and increment index number
    [pname,bname,ename]=fileparts(file_name)
    name_out=msprintf("%s%s_%03d%s",pname,base_name,iStart,ename)
    iStart=iStart+1
// Do actual copy, or rename, of file
    if( bKeep ) then
      [iOK,mess]=copyfile(file_name,name_out)
      if(iOK == 0) then AbortBox(["Unable to copy to file "+name_out;mess]); end
    else
      [iOK,mess]=movefile(file_name,name_out)
      if(iOK == 0) then AbortBox(["Unable to rename to file "+name_out;mess]); end
    end
  end
//
// Display results to prove we did something
  if( bKeep ) then
    InfoBox(msprintf("Copied %d files",size(get_files,2)))
  else
    InfoBox(msprintf("Renamed %d files",size(get_files,2)))
  end
//
endfunction


function IntegFileSelect_Popup(iFile)
// Create a popup to select a file from the batch files list,
// then display the integration ellipses in the callback
//
// Create a popup without the menus and toolbar
  f=CreateEmptyPopup([300 235],"Select an image file")
//
// Create listbox for image files
  AddListUI(f,[20,50,260,170],BatchFiles+".tif",12,1)
  f.children(1).value=iFile
//
// Pushbuttons to accept and cancel setting
  h1=AddButtonUI(f,[ 60 5 70 30],"Accept",18)
  h2=AddButtonUI(f,[160 5 70 30],"Cancel",18)
  AddUICallback(h1,"IntegFileSelect_CB","accept")
  AddUICallback(h2,"IntegFileSelect_CB","cancel")
//
endfunction


function IntegFileSelect_CB()
// Get the current graphic object and its parent figure
  g=gcbo();
  f=g.parent;
// If "accept", store selected file
  if(g.tag == "accept") then
    CreateIntegImageWindow(f.children(3).value)
  end
// Remove the popup
  close(f);
// If necessary, turn off Next or Previous menu items
  ManageIntegMenus()
endfunction


function WaveFileSelect_Popup()
//
// Read list of wavelength files
  CloseAllFiles()
  [fd,ierr]=mopen(ProgDir+"wavs.lis","r")
  if(ierr ~= 0) then
    AbortBox(["Unable to find ""wavs.lis"" file"; ...
            "Has LaueG been correctly installed?"])
  end
//
// Read the file list and comments from the wavs.lis
  names=[]
  comments=[]
  while (%t) then
    s=mfscanf(fd,"%s")
    if(s == []) then break; end
    names($+1)=ProgDir+s
    comments($+1)=stripblanks(mgetl(fd,1))
  end
  mclose(fd)
//
// Get the line in "comments" that matches the instrument
  ival=grep(convstr(comments),convstr(ImageInfo.instrum))
  if(ival == []) then
    ival=1
  else
    ival=ival(1)
  end
//
// If it exists, add the new wavelength file and set ival to it
  if( ~isempty(fileinfo("laue4_wav.new")) ) then
    names($+1)="laue4_wav.new"
    comments($+1)="New wavelength file created by Laue4"
    ival=size(names,1)
  end
//
// Create popup to let user select the file to use
  f=CreateEmptyPopup([300 235],"Select a wavelength file to use")
//
// Create listbox for image files
  h=AddListUI(f,[20,50,260,170],comments,12,1)
  f.children(1).value=ival
//
// Pushbuttons to accept and cancel setting
  h1=AddButtonUI(f,[ 60 5 70 30],"Accept",18)
  h2=AddButtonUI(f,[160 5 70 30],"Cancel",18)
// Add callbacks and tags for the pushbuttons
  AddUICallback(h1,"WaveFileSelect_CB","accept")
  AddUICallback(h2,"WaveFileSelect_CB","cancel")
//
// Copy Laue4 flag and file names for callback
  f.user_data=names
//
endfunction


function WaveFileSelect_CB()
//
// Get gcbo, its parent figure, and user_data
  g=gcbo();
  f=g.parent;
  names=f.user_data
//
// If "cancel", close popup
  if(g.tag == "cancel") then
    close(f)
    return
  end
//
// Get selected line before closing popup
  n=f.children(3).value
  close(f)
//
// Complain and die if file doesn"t exist
  if( isempty(fileinfo(names(n))) ) then
    AbortBox(["Unable to find wavelength file:";names(n)])
  end
//
// Copy the file into the current path
  copyfile(names(n),"laue4_wav.dat")
  printf( "Loaded wavelength file: %s\n",names(n) )
  WriteLogFile( "Loaded wavelength file: "+names(n) )
//
endfunction



////////////////////////////////////////////////////////////
// -------- Popups to manage modulation vectors  -------- //
////////////////////////////////////////////////////////////


function ModulMainSpot_Popup()
// Used to switch mouse display to satellite HKL vector
//
// Make a popup to set display of cursor position in terms of satellites
// around a central spot
//
// Create a popup window and strip off the menus and toolbar
  f=CreateEmptyPopup([260 120],"Satellite HKL Mode")
//
  AddTextUI(f,[10 70 140 25],"Central Spot hkl:",13)
  str=msprintf("%d %d %d",DisplayData.modul_main_hkl)
  AddEditUI(f,[150 70 50 25],str,13)
//
// Pushbuttons for test and cancel
  h1=AddButtonUI(f,[25 10 70 30],"Accept",18)
  h2=AddButtonUI(f,[120 10 70 30],"Cancel",18)
  AddUICallback(h1,"ModulMainSpot_CB","accept")
  AddUICallback(h2,"ModulMainSpot_CB","cancel")
//
endfunction


function ModulMainSpot_CB()
global DisplayMode DisplayData
//
// Get graphics object and parent figure
  g=gcbo()
  f=g.parent
//
// If "cancel", close popup
  if(g.tag == "cancel") then
    close(f)
    return
  end
//
// Read numbers, complain if invalid, non-integer, or zero
  vals=ReadPopupEditsNumber(f,3)
  if(vals == []) then return; end
  hkl=vals(1:3)
//
  if( max(abs(hkl-round(hkl))) > 0.001 ) then
    ErrorBox("Main spot HKL is non-integer")
    return
  end
  if(norm(hkl) == 0) then
    ErrorBox("Main spot HKL is zero")
    return
  end
  if(nmax == 0) then
    ErrorBox("Max satellite indices is zero")
    return
  end
//
// Set cursor to display in satellite mode
  DisplayMode.cursor_info=3
  DisplayData.modul_main_hkl=hkl
//
// Remove popup
  close(f)
endfunction


function ShowModul_Popup()
//
// Give warning if current d-spacing limits won"t display spots
  if(ModulateInfo.d_min >= ModulateInfo.d_max)
    WarnBox("Reset d-spacing min & max to display modulation spots")
  end
//
// Create a popup window and strip off the menus and toolbar
  f=CreateEmptyPopup([235 120],"Show Modulation Spots")
//
  AddTextUI(f,[10 90 180 25],"d-spacing limits (for display)",13)
  str=string(ModulateInfo.d_min)+" "+string(ModulateInfo.d_max)
  AddEditUI(f,[140 60 70 25],str,13)
//
// Pushbuttons for test and cancel
  h1=AddButtonUI(f,[5 10 70 30],"Accept",18)
  h2=AddButtonUI(f,[85 10 60 30],"Test",18)
  h3=AddButtonUI(f,[155 10 70 30],"Cancel",18)
  AddUICallback(h1,"ShowModul_CB","accept")
  AddUICallback(h2,"ShowModul_CB","test")
  AddUICallback(h3,"ShowModul_CB","cancel")
//
// Prevent window closing except by the pushbuttons
  f.closerequestfcn="InfoBox(""You must use the push buttons to exit this popup"")"
//

// Copy current ModulateInfo and display mode to user_data
  f.user_data=list(ModulateInfo,DisplayMode.modul_show)
//
endfunction


function ShowModul_CB()
global ModulateInfo DisplayMode
//
// Get graphics object and parent figure
  g=gcbo()
  f=g.parent
//
// If "cancel", restore original ModulateInfo and close popup
  if(g.tag == "cancel") then
    ModulateInfo=f.user_data(1)
    DisplayMode.modul_show=f.user_data(2)
    DrawGenHKLs(0)
    close(f)
    return
  end
//
// Read d-spacing limits for main spot
  buff=ReadPopupEdits(f)
  [dspace,bErr]=strvec2reals(buff,2)
  if( bErr ) then
    ErrorBox("Invalid d-spacing")
    return
  end
//
// Display satellite spots
  DisplayMode.modul_show=%t
  ModulateInfo.d_min=dspace(1)
  ModulateInfo.d_max=dspace(2)
  DrawGenHKLs(0)
//
// If "accept", save to settings file and close popup
  if(g.tag == "accept") then
    SaveSettingsFile()
    close(f)
  end
//
endfunction


function ChangeModul_Popup()
global DisplayMode
//
  if(ModulateInfo == []) then
    AbortBox("Cannot proceed. No modulation data to change.")
  end
//
// Create a popup window and strip off the menus and toolbar
  f=CreateEmptyPopup([320 230],"Change Modulation Vectors")
  AddTextUI(f,[75 200 180 25],"Value",13)
  AddTextUI(f,[235 200 180 25],"Step",13)
//
// Add vector values, step size, and <> buttons, and disable unused ones
  nvecs=size(ModulateInfo.vecs,1)
  for i=1:3
    y=210-35*i
    str1="? ? ?"
    str2="? ? ?"
    if(i <= nvecs) then
      str1=msprintf("%.3f %.3f %.3f",ModulateInfo.vecs(i,1:3))
      str2=msprintf("%.3f %.3f %.3f",ModulateInfo.vecs(i,1:3)/10)
    end
    AddEditUI(f,[30 y 120 25],str1,13)
    AddEditUI(f,[190 y 120 25],str2,13)
    h1=AddButtonUI(f,[6 y+1 20 23],"<",15)
    h2=AddButtonUI(f,[155 y+1 20 23],">",15)
    AddUICallback(h1,"ChangeModul_CB",string(-i))
    AddUICallback(h2,"ChangeModul_CB",string(i))
    if(i > nvecs) then
       f.children(1:4).enable="off"
    end
  end
//
  AddTextUI(f,[10 60 170 25],"d-spacing limits (for display)",13)
  str=string(ModulateInfo.d_min)+" "+string(ModulateInfo.d_max)
  AddEditUI(f,[200 60 70 25],str,13)
//
// Pushbuttons for test and cancel
  h1=AddButtonUI(f,[25 10 70 30],"Accept",18)
  h2=AddButtonUI(f,[120 10 70 30],"Test",18)
  h3=AddButtonUI(f,[215 10 70 30],"Cancel",18)
  AddUICallback(h1,"ChangeModul_CB","accept")
  AddUICallback(h2,"ChangeModul_CB","test")
  AddUICallback(h3,"ChangeModul_CB","cancel")
//
// Prevent window closing except by the pushbuttons
  f.closerequestfcn="InfoBox(""You must use the push buttons to exit this popup"")"
//
// Turn on display of modulations
  DisplayMode.modul_show=%t
//
// Copy ModulateInfo to user_data
  f.user_data=ModulateInfo
//
endfunction


function ChangeModul_CB()
global ModulateInfo
//
// Get graphics object and parent figure
    g=gcbo()
    f=g.parent
//
// If "cancel", restore orginal values and close popup
  if(g.tag == "cancel") then
    ModulateInfo=f.user_data
    DrawGenHKLs(0)
    close(f)
    return
  end
//
// Decode the edit boxes
  buff=ReadPopupEdits(f)
  nvecs=size(ModulateInfo.vecs,1)
//
// Check for weird SCILAB bug
  if(buff($) == "") then
    ErrorBox(["A SCILAB bug has occurred that I can""t fix!"; ...
              "The popup will close, then you can try again."])
    close(f)
    return
  end
// Vector values
  [fVec,bErr]=strvec2reals(buff(2*(1:nvecs)-1),3*nvecs)
  fVec=fVec'
  if( bErr ) then
    ErrorBox("Invalid vector value")
    return
  end
// Vector steps
  [fStep,bErr]=strvec2reals(buff(2*(1:nvecs)),3*nvecs)
  fStep=fStep'
  if( bErr ) then
    ErrorBox("Invalid vector step")
    return
  end
// d-spacing limits
  [dspace,bErr]=strvec2reals(buff($),2)
  if( bErr ) then
    ErrorBox("Invalid d-spacing")
    return
  end
//
// If an arrow button, step the appropriate vector value
  itag=strtod(g.tag)
  if( ~isnan(itag) ) then
    ivec=abs(itag)
    fVec(ivec,:)=fVec(ivec,:)+fStep(ivec,:)*(itag/ivec)
    f.children($+2-4*ivec).String=msprintf("%.3f %.3f %.3f",fVec(ivec,:))
  end
//
// Update values and display modulation spots
  ModulateInfo.vecs=fVec
  ModulateInfo.d_min=dspace(1)
  ModulateInfo.d_max=dspace(2)
  DrawGenHKLs(0)
//
// If "Accept", save to settings file and close popup
  if(g.tag == "accept") then
    SaveSettingsFile()
    close(f)
  end
//
endfunction


///////////////////////////////////////////////////////////
// ------- Convenience routines for popup boxes  ------- //
///////////////////////////////////////////////////////////

function f=CreateEmptyPopup(siz,str)
// Create a window without menus or toolbar and a canvas of "siz" pixels
//
// Close any existing figure with the same title
  CloseExistingPopup(str)
//
// Set default figure to about the right size
  df=gdf()
  df.axes_size=siz
  df.figure_size=siz+[8,71]
//
// Create figure
  f=createWindow()
  f.figure_name=str
  f.axes_size=siz
  f.auto_resize="off"
//
//
// Set figure to accommodate an axis of "siz" pixels
  sleep(50)
  f.auto_resize="on"
  f.axes_size=siz
//
endfunction


function CloseExistingPopup(name)
// Close figure with id=0 to 10 and matching title
  for i=0:10
    f=get_figure_handle(i)
    if(f == []) then
      continue
    end
    if(f.figure_name == name) then
      close(f)
    end
  end
endfunction


function SetPopupEdits(f,strings)
//
// Find all "edit" widgets
  iedit=find(f.children.style=="edit")
//
  if(size(iedit,"*") ~= size(strings,"*")) then
    BugBox("Wrong size for strings[]")
  end
//
// Copy string vector in reverse order
  for i=1:size(iedit,"*")
    f.children(iedit(i)).string=strings($-i+1)
  end
//
endfunction


function SetPopupChecks(f,bools)
//
// Find all "checkbox" widgets
  iedit=find(f.children.style=="checkbox")
//
  if(size(iedit,"*") ~= size(bools,"*")) then
    BugBox("Wrong size for bools[]")
  end
//
// Copy 0/1 values in reverse order
  for i=1:size(iedit,"*")
    f.children(iedit(i)).Value=double(bools($-i+1))
  end
//
endfunction


function [vals]=ReadPopupEditsNumber(f,nvals)
//
// Concatenate all "edit" strings together
  str=strcat(ReadPopupEdits(f)," ")
//
// Convert to a vector of numbers
  [vals,tails]=strtod(tokens(str))
//
// If any problem decoding numbers, complain and return []
  if( (size(vals,1) ~= nvals) | or(isnan(vals)) | (strcat(tails) ~= "") ) then
    ErrorBox("Invalid/missing number(s)")
    vals=[]
    return
  end
//
  vals=vals'
//
endfunction


function [strings]=ReadPopupEdits(f)
//
// Read the string value for all "edit" widgets
  i=find(f.children.style=="edit")
  if(i == []) then
    strings=[]
  else
    strings=f.children(i).string
  end
// Reverse the order of the string vector
  strings=strings($:-1:1)  
//
endfunction


function [vals]=ReadPopupChecks(f)
//
// Read the values for all "checkbox" widgets
  i=find(f.children.style=="checkbox")
  if(i == []) then
    vals=[]
  else
    vals=(f.children(i).value == 1)
  end
// Reverse the order of the values
// NB: SCILAB ($:-1:1) fails for unitary boolean
  if(size(vals,"*") > 1) then
    vals=vals($:-1:1)
  end
//
endfunction


function SetupHintUI()
global HintUI
  HintUI=[]
  HintUI.text=[]
endfunction


function HintUI_CB()
global HintUI
//
// Get graphics object, number of tag, and text of hint
  g=gcbo()
  num=strtod(g.tag)
  str=HintUI.text(num)
// Convert string with \n to array of strings
  strs=strsplit(str,"/\\n/")
// Display hint in modal message box
  messagebox(strs, "modal", "info", "OK")
//
endfunction



/////////////////////////////////////////////////////////////////
// ------- Routines for adding widgets to popup boxes  ------- //
/////////////////////////////////////////////////////////////////

function h=AddTextUI(fig,pos,str,siz)
// Add a "text" UI to a figure
  if(str == []) then
    str=" "
  end
//
  h=uicontrol(fig,"style","text","position",pos,"string",str, ...
                 "FontSize",siz,"FontUnits","points", ...
                 "FontName","helvetica","HorizontalAlignment","left", ...
                 "BackgroundColor",[0.8 0.8 0.8])
endfunction


function h=AddEditUI(fig,pos,str,siz)
// Add an "edit" UI to a figure
  if(str == []) then
    str=" "
  end
//
  h=uicontrol(fig,"style","edit","position",pos,"string",str, ...
                 "FontSize",siz,"FontUnits","points", ...
                 "FontName","helvetica","HorizontalAlignment","left", ...
                 "Relief","flat", ...
                 "BackgroundColor",[1 1 1])
endfunction


function h=AddListUI(fig,pos,str,siz,val)
// Add a "listbox" UI to a figure
  if(str == []) then
    str=" "
  end
//
  h=uicontrol(fig,"style","listbox","position",pos,"string",str, ...
                 "value",val,"relief","sunken", ...
                 "FontSize",siz,"FontUnits","points", ...
                 "FontName","helvetica","HorizontalAlignment","left", ...
                 "BackgroundColor",[1 1 1])
endfunction


function h=AddCheckUI(fig,posxy,bVal)
// Add a "pushbutton" UI to a figure
  val=double(bVal)
  pos=[posxy(1),posxy(2)+3,20,20]
  h=uicontrol(fig,"style","checkbox","position",pos,"value",val,...
                 "BackgroundColor",[0.8 0.8 0.8])
endfunction


function h=AddButtonUI(fig,pos,str,siz)
// Add a "pushbutton" UI to a figure
  if(str == []) then
    str=" "
  end
//
  h=uicontrol(fig,"style","pushbutton","position",pos,"string",str, ...
                 "FontSize",siz,"FontUnits","points", ...
                 "FontName","helvetica","HorizontalAlignment","center", ...
                 "Relief","raised", ...
                 "BackgroundColor",[0.8 0.8 0.8])
//
endfunction


function AddUICallback(h,cb,tag)
// Add a callback and tag to the UI with handle "h"
  if(cb ~= "") then
    set(h,"callback",cb)
  end
  if(tag ~= "") then
    set(h,"tag",tag)
  end
endfunction


function AddHintUI(fig,posxy,text)
global HintUI
//
//
// Create small "?" button
  pos=[posxy(1),posxy(2)+5,15,15]
  h=uicontrol(fig,"style","pushbutton","position",pos,"string","?", ...
                 "FontSize",9,"FontUnits","points", ...
                 "FontName","helvetica","HorizontalAlignment","center", ...
                 "Relief","raised", ...
                 "BackgroundColor",[0.8 0.8 0.8])
//
// Add text to array of hints, and register tag with callback
  num=size(HintUI.text,1)+1
  HintUI.text(num)=text
  tag=string(num)
  set(h,"callback","HintUI_CB")
  set(h,"tag",tag)
//
endfunction
