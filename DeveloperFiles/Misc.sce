global LatticeInfo ModulateInfo BatchFiles bDevelopMode ZoomAxes
global ImageFigure ProgBarId

// ========= Create figure without menus =========
//function f=CreateEmptyFigure(siz)
//function RemoveMenuBar(f)
// ========= Manage menus + progress bar =========
//function AddConsoleMenu(label,items,callback)
//function AddImageMenu(label,items,callback)
//function LockMenus(str)
//function UnlockMenus()
//function ChangeMenus(fig,sMenu,sItem,sOn)
//function RemoveProgBar()
//function AbortLaueG()
// ============= Message box routines ===============
//function ModalMessageBox(sText,sTitle,sType)
//function InfoBox(text)
//function InfoTitleBox(text,name)
//function WarnBox(text)
//function ErrorBox(text)
//function AbortBox(text)
//function BugBox(text)
//function result=QuestionBox(text)
//function file_name_out=OutputFileBox(file_name_in)
//function text=AddHtmlSpaces(text_in)
// ============= Status and information ===============
//function result=IsImageOn()
//function result=IsZoomOn()
//function result=IsBatchMode()
//function result=IsModulated()
//function str=LaueGVersion()
// ============= Decode strings ===============
//function [vals,bErr]=strvec2reals(str,num)
//function [vals,bErr]=strmat2reals(strvec)

// ========= Create figure without menus =========

function f=CreateEmptyFigure(siz)
// Create window to contain an image of size "siz" without menus or toolbar
//
// Create figure
// NB: Cannot use dockable="off" due to unsetmenu() bug
  f=figure("axes_size",siz,"resize","off","auto_resize","on", ...
                 "menubar","none","toolbar","none")
//
// After a delay, adjust figure size to get correct axes size
  sleep(50)
  f.figure_size=f.figure_size + (siz-f.axes_size)
  f.axes_size=siz
//
endfunction


function RemoveMenuBar(f)
// Remove the default menus and the toolbar from figure "f"
// The SCILAB "console" window is selected if f=[]
//
  if( f == [] ) then
    delmenu("File")
    delmenu("Edit")
    delmenu("Preferences")
    delmenu("Control")
    delmenu("Applications")
    delmenu("?")
    toolbar(-1,"off")
  else
// Loop through children removing any menu items
    delmenu(f.figure_id,"File")
    delmenu(f.figure_id,"Tools")
    delmenu(f.figure_id,"Edit")
    delmenu(f.figure_id,"?")
// Remove the toolbar
    toolbar(f.figure_id,"off")
  end
endfunction


// ========= Manage menus + progress bar =========

function AddConsoleMenu(label,items,callback)
//
// Add menu for console window
  addmenu(label,items,list(2,callback))
//
// Make text a nice blue
  fig=get(0)
  fig.children(1).Foregroundcolor=[0,0.3,1]
  fig.children(1).children(:).Foregroundcolor=[0,0.3,1]
//
endfunction


function AddImageMenu(label,items,callback)
//
// Add menu for ImageFigure
  fig=ImageFigure
  addmenu(fig.figure_id,label,items,list(2,callback))
//
// Make text a nice blue
  fig.children(1).Foregroundcolor=[0,0.3,1]
  fig.children(1).children(:).Foregroundcolor=[0,0.3,1]
//
endfunction


function LockMenus(str)
global ProgBarId
//
// If not an empty string, create a new progress bar
  if( str ~= "" ) then
    RemoveProgBar()
    ProgBarId=progressionbar(str)
  end
//
// Lock some console menus and items
  LockMainMenus()
//
// Lock top level menus of image window, if it exists
  if( IsImageOn() ) then
    ChangeMenus(ImageFigure,"*","","off")
  end
//
//
endfunction


function UnlockMenus()
//
// Close progress bar and unlock menus
  RemoveProgBar()
//
// Unlock console menus according to batch mode
  UnlockMainMenus()
//
// Unlock top level menus of image window, if it exists
  if( IsImageOn() ) then
    ChangeMenus(ImageFigure,"*","","on")
  end
//
endfunction


function ChangeMenus(fig,sMenu,sItem,sOn)
//
// Enables/disables menus and changes label colour
//
// DEBUG
//  disp(fig)
//  mprintf("   sMenu,sItem,sOn = %s,%s,%s\n",sMenu,sItem,sOn)
//
// Search the top menu items, expand "*" search to a string array
  i=find(fig.children.type == "uimenu")
  if(sMenu ~= "*") then
    i=i( find(fig.children(i).label == sMenu) )
    if(i == []) then BugBox( msprintf("1 M,I=%s,%s",sMenu,sItem) ); end
  end
  sCommand=msprintf("fig.children(%d)\n",i')
//
// Search the sub-menu items, implement "*" search using "(:)"
  if(sItem == "*") then
    sCommand=sCommand+".children(:)"
  else
    if(sItem ~= "") then
      i=find( evstr(sCommand+".children(:).label") == sItem )
      if(i == []) then BugBox( msprintf("2 M,I=%s,%s",sMenu,sItem) ); end
      sCommand=sCommand+msprintf(".children(%d)",i)
    end
  end
//
// DEBUG
//  for sLabel=sCommand'+".label"
//    mprintf("sLabels: %s = %s\n",sLabel,strcat(evstr(sLabel),","))
//  end
//
  if(convstr(sOn) == "on") then
    execstr(sCommand' +".Foregroundcolor=[0,0.3,1]")
    execstr(sCommand' +".Enable=%t")
  else
    execstr(sCommand' +".Foregroundcolor=[0.5,0.5,0.5]")
    execstr(sCommand' +".Enable=%f")
  end
//
endfunction


function RemoveProgBar()
//
// Close any existing progress bar
  if( type(ProgBarId) == 9 ) then
    if( is_handle_valid(ProgBarId) ) then
      execstr("close(ProgBarId)","errcatch")    // sometimes fails
    end
  end
//
endfunction


function AbortLaueG()
  UnlockMenus()
  abort
endfunction


// ============= Message box routines ===============

function ModalMessageBox(sText,sTitle,sType)
// Important to first remove the progress bar as I have
// hacked it to always be "on top"
  RemoveProgBar()
//
  messagebox(AddHtmlSpaces(sText), sTitle, sType, "modal")
endfunction


function InfoBox(text)
  ModalMessageBox(text,"LaueG Information","info")
endfunction


function InfoTitleBox(text,name)
  ModalMessageBox(text, name,"info")
endfunction


function WarnBox(text)
  ModalMessageBox(text, "LaueG Warning","warning")
// Save warning message if in batch mode
  if( IsBatchMode() )
    WriteLogFile("LaueG Warning: "+text)
  end
endfunction


function ErrorBox(text)
  ModalMessageBox(text, "LaueG Error","error")
// Save error message if in batch mode
  if( IsBatchMode() )
    WriteLogFile("LaueG Error: "+text)
  end
endfunction


function AbortBox(text)
// Remove progress bar in case it hides the message box
  RemoveProgBar()
// Write modal message
  ModalMessageBox([matrix(text,-1,1);"";blanks(15)+"Aborting Operation"], ...
                                                    "LaueG Error","error")
// Save abort message if in batch mode
  if( IsBatchMode() )
    WriteLogFile("LaueG Abort: "+text)
  end
// Do full LaueG abort
  AbortLaueG()
endfunction


function BugBox(text)
// Remove progress bar in case it hides the message box
  RemoveProgBar()
// Create string array for message with traceback
  [lnum,mac]=where()
  if(size(lnum,1) < 2) then
    lnum=[0,0]
    mac=["","""unknown routine"""]
  end
  calls=msprintf("Bug detected at line %i of %s",lnum(2),mac(2))
// Write modal message
  messagebox([calls,"",AddHtmlSpaces(text)], "LaueG BUG", "error", "modal")
// Write bug messages to the log file
  WriteLogFile("")
  WriteLogFileHeader("LaueG "+calls)
  WriteLogFile("     "+text)
  WriteLogFile("")
// Do full LaueG abort
  AbortLaueG()
endfunction


function result=QuestionBox(text)
  button=messagebox(AddHtmlSpaces(text), "LaueG Question", "question", ["Yes" "No"], "modal")
  result=(button == 1)
endfunction


function file_name_out=OutputFileBox(file_name_in)
//
// Ask if use suggested name, other name, or cancel
  button=messagebox("File name: "+file_name_in, ...
      "Output file name", "question", ["Yes" "Rename" "No"], "modal")
//
// Default for "no" and "cancel"
  file_name_out=""
//
// If "yes", make output name from DataDir path and input name
  if( button == 1 ) then
    file_name_out=DataDir+file_name_in
// If "rename", get name using the file-name-popup
  elseif( button == 2 ) then
    [fpath,fname,fext]=fileparts(file_name_in)
    file_name_out=uiputfile(fext,DataDir,"Name of output file")
// Convert a [] to a "" for the output name
    if(file_name_out == []) then
      file_name_out=""
    end

// If missing, add the extension from the input name
    if(file_name_out ~= "") then
      [fpath,fname,fext2]=fileparts(file_name_out)
      if(fext2 == "") then
        file_name_out=file_name_out+fext
      end
    end
  end
endfunction


function text=AddHtmlSpaces(text_in)
// Replaces double spaces with " &nbsp;" as message box uses HTML
 text=strsubst(text_in,"  "," &nbsp;")
endfunction


// ============= Status and information ===============

function result=IsImageOn()
  if( type(ImageFigure) == 9 ) then
    result= is_handle_valid(ImageFigure)
  else
    result=%f;
  end
endfunction

function result=IsZoomOn()
   result=%f;
  if( type(ZoomAxes) == 9 ) then
    if( is_handle_valid(ZoomAxes) ) then
      result=(ZoomAxes.visible == "on")
    end
  end
endfunction

function result=IsBatchMode()
  result=(size(BatchFiles,1) > 0)
endfunction

function result=IsModulated()
  if( ModulateInfo == [] )
    result=%f
  else
    result=(ModulateInfo.vecs ~= [])
  end
endfunction

function str=LaueGVersion()
// Used to output LaueG release version information
  if( bDevelopMode ) then
    str="LaueG (development version)"
  else
    str="LaueG (build: "+LaueGReleaseDate()+")"
  end
endfunction


// ============= Decode strings ===============

function [vals,bErr]=strvec2reals(str,num)
// Convert a vector of strings to a vector of reals
// Return bErr=%t if syntax error or not "num" reals found
// If num<1, accept any number of reals
//
  [vals,ierr]=evstr(strcat(str," "))
  bErr=( (ierr ~= 0) | (num>0 & (size(vals,2) ~= num)) )
//
endfunction


function [vals,bErr]=strmat2reals(strvec)
// Convert a vector of strings to a matrix of reals
// Return bErr=%t if syntax error or varying number of columns
// If strvec == "", return with [], but no error
//
  [vals,ierr]=evstr( matrix(strvec,-1,1) )
  bErr=(ierr ~= 0)
//
endfunction

