///////////////////////mode(-1)
//
// LaueG startup file
//
//
// Crash out if we are not running v5.5.2
if(getversion() ~= "scilab-5.5.2") then
  mprintf("\n\n\n")
  mprintf("    ##############################################\n")
  mprintf("    #                                            #\n")
  mprintf("    #  This LaueG requires SCILAB version 5.5.2  #\n")
  mprintf("    #                                            #\n")
  mprintf("    ##############################################\n\n")
  mprintf("Download and install SCILAB 5.5.2 (64-bit) from\n")
  mprintf("https://www.scilab.org/download/5.5.2 or the Koala page at ANSTO.\n")
  mprintf("Unfortunately, version 6 of SCILAB has many bugs and is unreliable.\n\n")
  x=input("Press <Enter> to exit LaueG")
  quit
end
//
// Need the largest stacksize
stacksize("max")
//
// Output start of wait message
mprintf("\n====== Initialising LaueG, please wait ======\n")
//
// Check if I am running the code in development mode
global bDevelopMode
bDevelopMode=%f
if( isfile("C:\DevApps\LaueG\MyVersion.sce") ) then
  if( convstr(PWD) == convstr("C:\DevApps\LaueG") ) then
    bDevelopMode=%t
  end
end
//
// Set ProgDir
if( bDevelopMode ) then
  ProgDir="C:\DevApps\LaueG\"
else
  ProgDir=pathconvert("SCI\laueg\")
  if( ~isfile(ProgDir+"Release.sce") ) then
    ProgDir=pathconvert("SCI\..\..\laueg\")
    if( ~isfile(ProgDir+"Release.sce") ) then
      disp("ERROR: Cannot find Release.sce")
      abort
    end
    ProgDir=pathconvert(cd(ProgDir))
  end
end
//
// Load functions from ProgDir
if( bDevelopMode ) then
// Names of all *.sce files, including this one (which must be first)
  global LaueGFiles
  LaueGFiles=["LaueG.sce", ...
              "MainMenus.sce","SingleImageMenus.sce", ...
              "IntegImageMenus.sce","RejectsImageMenus.sce", ...
              "FlickImageMenus.sce","PopupBoxes.sce", ...
              "DrawImage.sce","DrawOverlays.sce","DrawEllipses.sce", ...
              "ImageWindow.sce","EventHandler.sce","BatchMode.sce", ...
              "CalcSpots.sce","RunConsole.sce","Files.sce", ...
              "Init.sce","Misc.sce","MyVersion.sce"]
// Execute each *.sce file from the program directory
  funcprot(0)    // warnings off when functions changed
  mprintf("  Loading Script Files\n ")
  [str,n,line,func]=lasterror(%t)
  for name=LaueGFiles(2:$)
    mprintf("\r%s              ",name)
    exec(ProgDir+name,-1)
  end
  funcprot(1)    // warnings on when functions changed
  mprintf("\r  Finished Loading                   \n")
else
// Execute Release.sce in the program directory
  mprintf("  Loading Script File\n")
  exec(ProgDir+"Release.sce",-1)
end
//
// Setup the pulldown menus on the main (console) window
mprintf("  Setup Menus\n")
SetupMainMenus()
//
// Add development menu for me
if( bDevelopMode ) then
  SetupDevelMenu()
end
//
// Set the directory for ini & wavelength files
mprintf("  Setup Defaults\n")
cd("SCIHOME\..\")
mkdir("laueg")      // harmless if it already exists
IniDir=pathconvert(cd("laueg"))
//
// Load the laueg.ini file
ReadIniFile()
//
// Move to the data directory
ChangeDataDir(DataDir,%f)
//
// Initialise main structures
LoadSetupsFile()      // read instrument setups file
LoadDummySetup()      // set instrument to the dummy setup
DummyLatticeInfo()    // set lattice, simple cubic, no twins
DummyDisplaySet()     // set defaults for image display
DummyArgboxSet()      // set defaults for argonnes_boxes
ResetDisplay()        // clear all display modes
//
// Try to load settings file
LoadSettingsFile()
//
// Output end of wait message and some useful information
mprintf("====== Initialisation complete ======\n\n")
mprintf("%s running on %s (%s)\n",LaueGVersion(),convstr(getversion(),"u"),getos())
if(InstrumDefault > 1) then
  mprintf("Default instrument setup: %s\n", ...
            InstrumSetups(InstrumDefault).name)
end
mprintf("Working directory: %s\n",DataDir)
//
// Turn off warnings for users, turn on for me
if( bDevelopMode ) then
  warning("on")
else
  warning("off")
end
//
mode(1)
