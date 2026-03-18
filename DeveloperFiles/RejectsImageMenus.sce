global ImageFigure MainImage RejectsInfo

// =================  Menu definitions and callbacks  =============
//function SetupRejectsMenus(bOutliers)
//function RejectsFile_Menu(iMenu,iWin)
//function RejectsDisplay_Menu(iMenu,iWin)
//function RejectsRejects_Menu(iMenu,iWin)
//function RejectsOutliers_Menu(iMenu,iWin)
// ========================  Helper routines  ===================
//function LoadRejectImagesData(ifirst)
//function UpdateRejectsImage(sMenu,ifirst)
//function CreateRejectFile(itags)
//function AppendRejectFile(itags)



// =================  Menu definitions and callbacks  =============

function SetupRejectsMenus(bOutliers)
// Manage title/toolbar/pulldowns for the "reject spots" image window
// bOutliers=TRUE for "Check Outlier Spots", FALSE for "Edit Reject Spots"
//
// Total number of spots
// Add the menus common to "rejects" and "outliers"
  AddImageMenu("File", ...
                ["Save Displayed Image", "Close Image"], ...
               "RejectsFile_Menu")
  AddImageMenu("Display", ...
                ["Invert Image","Brighter Image","Darker Image", ...
                 "Brightness Levels","Toggle Background Ellipses"], ...
                "RejectsDisplay_Menu")
// Add the menus specific to "outliers" or "rejects"
  if( bOutliers ) then
    AddImageMenu("Outliers", ...
                  ["Next 50","Previous 50","Add Marked to Rejects", ...
                   "Add Unmarked to Rejects"], ...
                  "RejectsOutliers_Menu")
  else
    AddImageMenu("Rejects", ...
                  ["Next 50","Previous 50","Remove Marked from Rejects", ...
                   "Remove Unmarked from Rejects"], ...
                  "RejectsRejects_Menu")
  end
//
// Turn off unneeded navigation menu items
  if( bOutliers ) then
    str="Outliers"
  else
    str="Rejects"
  end
  if(size(RejectsInfo.itwin,1) <= 50) then
    ChangeMenus(ImageFigure,str,"Next 50","off")
  end
  ChangeMenus(ImageFigure,str,"Previous 50","off")
endfunction


function RejectsFile_Menu(iMenu,iWin)
// Save the image to a GIF file
  if(iMenu == 1) then
    name=uiputfile(["*.gif"],DataDir,"Name of GIF file")
    if( isempty(name)) then return; end
    [fpath,fname,fext]=fileparts(name)
    if(fext == "") then
      name=name+".gif"
    end
    xs2gif(ImageFigure,name)
  end
// Close image window
  if(iMenu == 2) then
    delete(ImageFigure);
  end
endfunction


function RejectsDisplay_Menu(iMenu,iWin)
global DisplaySet
//
// Invert black & white (for Garry)
  if(iMenu == 1) then
    InvertColorMap();
    DisplaySet.invert_grey=~(DisplaySet.invert_grey)
    SaveSettingsFile()
  end
//
// Change max & min intensities for a brighter display
  if(iMenu == 2) then
    fac=1.0/1.5
    if(ImageFigure.color_map(1,1) > 0.5) then fac=1.5; end
    RedrawImageNow( MainImage.user_data*fac )
  end
//
// Change max & min intensities for a darker display
  if(iMenu == 3) then
    fac=1.0/1.5
    if(ImageFigure.color_map(1,1) > 0.5) then fac=1.5; end
    RedrawImageNow( MainImage.user_data/fac )
  end
//
// Change max & min intensities for display
  if(iMenu == 4) then
    Brightness_Popup()
  end
//
// Toggle background ellipses
  if(iMenu == 5) then
    DisplaySet.ell_back=~DisplaySet.ell_back
    SaveSettingsFile()
// drawlater/drawnow only works on the current figure!
    scf(ImageFigure)
    drawlater
    DrawCurrentEllipses()
    drawnow
  end
endfunction


function RejectsRejects_Menu(iMenu,iWin)
//
// Next 50 spots
  if(iMenu == 1) then
    UpdateRejectsImage("Rejects",RejectsInfo.ifirst+50)
  end
// Previous 50 spots
  if(iMenu == 2) then
    UpdateRejectsImage("Rejects",RejectsInfo.ifirst-50)
  end
// Remove Marked spots from Rejects file
  if(iMenu == 3) then
    itags=find(~RejectsInfo.marked);
    CreateRejectFile(itags)
  end
// Remove Unmarked spots from Rejects file
  if(iMenu == 4) then
    itags=find(RejectsInfo.marked);
    CreateRejectFile(itags)
  end
endfunction


function RejectsOutliers_Menu(iMenu,iWin)
//
// Next 50 spots
  if(iMenu == 1) then
    UpdateRejectsImage("Outliers",RejectsInfo.ifirst+50)
  end
// Previous 50 spots
  if(iMenu == 2) then
    UpdateRejectsImage("Outliers",RejectsInfo.ifirst-50)
  end
// Add Marked spots to Rejects file
  if(iMenu == 3) then
    itags=find(RejectsInfo.marked);
    AppendRejectFile(itags)
  end
// Add Unmarked spots to Rejects file
  if(iMenu == 4) then
    itags=find(~RejectsInfo.marked);
    AppendRejectFile(itags)
  end
endfunction



// ========================  Helper routines  ===================

function LoadRejectImagesData(ifirst)
global RejectsInfo
// Load ReducedImage with 100x100 images starting from reject "ifirst"
//
// Store the number of the first image displayed
  RejectsInfo.ifirst=ifirst
// Create the rejects image and copy to ReducedImage
  nspots=size(RejectsInfo.itwin,1)
  ilast=min(ifirst+49,nspots)
  RunRejectsImage(RejectsInfo.base_name(ifirst:ilast), ...
                             RejectsInfo.xycen(ifirst:ilast,1:2))
endfunction


function UpdateRejectsImage(sMenu,ifirst)
//
// Load up to 50 rejects starting from ifirst
  LoadRejectImagesData(ifirst)
//
// Set/unset Next & Previous menu items
  nspots=size(RejectsInfo.itwin,1)
  if(ifirst+49 >= nspots) then
    ChangeMenus(ImageFigure,sMenu,"Next 50","off")
  else
    ChangeMenus(ImageFigure,sMenu,"Next 50","on")
  end
  if(ifirst == 1) then
    ChangeMenus(ImageFigure,sMenu,"Previous 50","off")
  else
    ChangeMenus(ImageFigure,sMenu,"Previous 50","on")
  end
//
// Draw the new image, ellipses, and marks
  levs=CalcImageHisto(0.001,0.01)
// drawlater/drawnow only works on the current figure!
  scf(ImageFigure)
  drawlater
  RedrawImage(levs)
  DrawRejectEllipses()
  DrawRejectSpots()
  drawnow
//
endfunction


function CreateRejectFile(itags)
//
// Create a new file with a 3 line header
  CloseAllFiles()
  [fd,ierr]=mopen("laue4_rej.dat","w");
  if(ierr ~= 0) then
    AbortBox("Unable to open laue4_rej.dat")
  end
//
// Write a 3 line header depending on M index
  if(RejectsInfo.itwin(1) == 0) then
    mfprintf(fd,"LaueG Rejects File Version 2, Option 1\n" + ...
                "\n" + ...
                "Base-name                                  " + ...
                "h   k   l   Comments ...\n");
  else
    mfprintf(fd,"LaueG Rejects File Version 2, Option 2\n" + ...
                "\n" + ...
                "Base-name                                  " + ...
                "h   k   l   m   Comments ...\n");
  end
//
// Write out the tagged spots of RejectsSpotsInfo
  for i=itags
    if(RejectsInfo.itwin(1) == 0) then
      mfprintf(fd,"%-40s%4i%4i%4i  %s\n", ...
                  RejectsInfo.base_name(i),RejectsInfo.hkl(i,:), ...
                  RejectsInfo.comment(i)                        )
    else
      mfprintf(fd,"%-40s%4i%4i%4i%4i  %s\n", ...
                  RejectsInfo.base_name(i),RejectsInfo.hkl(i,:), ...
                  RejectsInfo.itwin(i),RejectsInfo.comment(i)   )
    end
  end
  mclose(fd)
//
// Give a summary to the user
  nspots=size(itags,2)
  InfoBox(msprintf("%i spots written to laue4_rej.dat",nspots))
endfunction


function AppendRejectFile(itags)
//
// If no tagged rejects, nothing to do
  if(itags == []) then
    return
  end
//
// Have to create the rejects file if it doesn't exist
  file_name="laue4_rej.dat";
  if( ~isfile(file_name) ) then
    CreateRejectFile(itags);
    return
  end
//
// Make the equivalent of RejectsInfo from the rejects file
  RInfo=RunRejectsData(file_name,%f)
//
// If rejects file is empty, just write out the new rejects
  if(RInfo.itwin == []) then
    CreateRejectFile(itags);
    return
  end
//
// If the two RejectsInfo's are incompatible, ask what to do
  bTwin=(RejectsInfo.itwin(1) ~= 0)
  bTwin2=(RInfo.itwin(1) ~= 0)
  if(bTwin ~= bTwin2) then
    text=["Old and new rejected spots are incompatible in terms";
          "of twins. Continue and delete old rejected spots?"]
    if( QuestionBox(text) ) then
      CreateRejectFile(itags);
    end
    return
  end
//
// Append tagged rejects to the rejects file
  CloseAllFiles()
  [fd,ierr]=mopen(file_name,"a");
  if(ierr ~= 0) then
    AbortBox("Unable to open laue4_rej.dat")
  end
//
// Write out the tagged spots of RejectsInfo
  for i=itags
    if(RejectsInfo.itwin(1) == 0) then
      mfprintf(fd,"%-40s%4i%4i%4i  %s\n", ...
                  RejectsInfo.base_name(i),RejectsInfo.hkl(i,:), ...
                  RejectsInfo.comment(i)                        )
    else
      mfprintf(fd,"%-40s%4i%4i%4i%4i  %s\n", ...
                  RejectsInfo.base_name(i),RejectsInfo.hkl(i,:), ...
                  RejectsInfo.itwin(i),RejectsInfo.comment(i)   )
    end
  end
  mclose(fd)
//
// Give a summary to the user
  nspots=size(itags,2)
  InfoBox(msprintf("%i spots added to laue4_rej.dat",nspots))
endfunction
