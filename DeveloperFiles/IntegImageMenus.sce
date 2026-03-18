global ImageFigure BatchFiles DisplayMode DataDir

//function SetupIntegMenus()
//function IntegFile_Menu(iMenu,iWin)
//function IntegEll_Menu(iMenu,iWin)
//function ManageIntegMenus()


// ========================  Menu definitions  ===================

function SetupIntegMenus()
// Manage title/toolbar/pulldowns for images window to
// display the image integration ellipses in batch mode
//
// Create the file menu options
  fig_id=ImageFigure.figure_id
//
// Create the Next/Previous/Select image menu options
  AddImageMenu("File", ...
            ["Next Image","Previous Image","Select Image", ...
             "Save Displayed Image", ...
             "Show File Info.","Show Index Info.", ...
             "Show Instrument Info.", ...
             "Close Image"], ...
            "IntegFile_Menu")
//
// Create image display menu options
// >>> This menu list must agree with SetupImageMenus() <<<
  AddImageMenu("Display", ...
            ["Lighter Image","Darker Image","Brightness Levels", ...
             "Invert Grey Levels","Strip Background","Reload Background", ...
             "Blur Main Image","Remove Blur"], ...
            "ImageDisplay_Menu")        // NB: uses ImageDisplay_Menu
//
// Create types of ellipses menu options
  AddImageMenu("Ellipses", ...
            ["Integration Ellipses","Model Ellipses","No Ellipses", ...
             "Toggle Background Ellipses"], ...
            "IntegEll_Menu")
//
// Possibly turn off menu items
  ManageIntegMenus()
//
endfunction


// ========================  Menu callbacks  ===================



function IntegFile_Menu(iMenu,iWin)
//
// Display the "Next" batch file
  if (iMenu == 1) then
    CreateIntegImageWindow(DisplayMode.ibatch+1)
  end
//
// Display the "Previous" batch file
  if (iMenu == 2) then
    CreateIntegImageWindow(DisplayMode.ibatch-1)
  end
//
// "Select" a specific batch file from the list
  if (iMenu == 3) then
    IntegFileSelect_Popup(DisplayMode.ibatch);
  end
//
// Output the image to a GIF file
  if(iMenu == 4) then
    name=uiputfile(["*.gif"],DataDir,"Name of GIF file")
    if( isempty(name)) then return; end
    [fpath,fname,fext]=fileparts(name)
    if(fext == "") then
      name=name+".gif"
    end
    xs2gif(ImageFigure,name)
  end
//
// Show info about the image
  if(iMenu == 5) then
    ShowImageInfo()
  end
//
// Show index info
  if(iMenu == 6) then
    ShowLatticeInfo()
  end
//
// Show instrument info
  if(iMenu == 7) then
    ShowInstrumInfo()
  end
//
// Close the image window
  if(iMenu == 8) then
    delete(ImageFigure)
  end
//
// If necessary, turn off Next or Previous menu items
  ManageIntegMenus()
endfunction


function IntegEll_Menu(iMenu,iWin)
global DisplaySet
// Draw integration ellipses
  if (iMenu == 1) then
// drawlater/drawnow only works on the current figure!
    scf(ImageFigure)
    drawlater
    DrawIntegEllipses()
    drawnow
  end
// Draw model ellipses
  if (iMenu == 2) then
// drawlater/drawnow only works on the current figure!
    scf(ImageFigure)
    drawlater
    DrawModelEllipses()
    drawnow
  end
// Clear all ellipses
  if (iMenu == 3) then
    ClearEllipses()
  end
// Toggle background ellipses, and redraw ellipses
  if (iMenu == 4) then
    DisplaySet.ell_back=~DisplaySet.ell_back
    SaveSettingsFile()
// drawlater/drawnow only works on the current figure!
    scf(ImageFigure)
    drawlater
    DrawCurrentEllipses()
    drawnow
  end
endfunction


function ManageIntegMenus()
// Turn off Previous image if at the first batch file
  if(DisplayMode.ibatch == 1) then
    ChangeMenus(ImageFigure,"File","Previous Image","off")
  else
    ChangeMenus(ImageFigure,"File","Previous Image","on")
  end
//
// Turn off Next image if at the last batch file
  if(DisplayMode.ibatch == size(BatchFiles,1)) then
    ChangeMenus(ImageFigure,"File","Next Image","off")
  else
    ChangeMenus(ImageFigure,"File","Next Image","on")
  end
//
// Turn off Select if only one batch file
  if(size(BatchFiles,1) == 1) then
    ChangeMenus(ImageFigure,"File","Select Image","off")
  end
endfunction
