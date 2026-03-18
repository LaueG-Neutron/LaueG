; Create an Install program for LaueG

#define MyDateTimeString GetDateTimeString('dd-mmm-yyyy', '-', '');

[Setup]
AppVersion={#MyDateTimeString}
AppName=LaueG
OutputDir=C:\DevApps\LaueG\_Release\Installs\
OutputBaseFilename=LaueG_setup_{#MyDateTimeString}
SourceDir=C:\DevApps\LaueG\_Release\WindowsBuild\
DisableDirPage=yes
UsePreviousAppDir=no
DefaultDirName="C:\Program Files\"
InfoAfterFile=setup_finished.txt
SetupLogging=yes
ArchitecturesAllowed= x64 ia64
PrivilegesRequired=none
Uninstallable=no
WindowVisible=yes


[Files]
; Copy the executables to the "prog" folder
Source: "argonne_boxes.exe"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "find_spots.exe"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "gen_hkls.exe"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "image_info.exe"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "index_conic.exe"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "index_spots.exe"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "laue4_exe.exe"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "laueg_cyclops.exe"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "make_ellipses.exe"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "orient_spots.exe"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "rejects_data.exe"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "rejects_image.exe"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "strip_back.exe"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "run_pause.bat"; DestDir: "{code:GetScilabHome}\laueg\";
; Copy the SCILAB script files to the "prog" folder
Source: "LaueG.sce"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "Release.sce"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "LaueG_setup.sce"; DestDir: "{code:GetScilabHome}\laueg\";
; Copy the wavelength distribution files to the "prog" folder,
; also create the desktop shortcut
Source: "wavs.lis"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "koala_wav.dat"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "imagine_wav.dat"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "vivaldi_wav.dat"; DestDir: "{code:GetScilabHome}\laueg\"; AfterInstall: MakeShortcut();
; Copy the instruments setups file to "prog" folder
Source: "setups.dat"; DestDir: "{code:GetScilabHome}\laueg\";
; Copy the icon file to the "prog" folder
Source: "laueg_icon.ico"; DestDir: "{code:GetScilabHome}\laueg\";
; Copy the LaueG Manual htm file and folder to the "prog" folder
Source: "LaueG Manual.htm"; DestDir: "{code:GetScilabHome}\laueg\";
Source: "LaueG Manual_files\*"; DestDir: "{code:GetScilabHome}\laueg\LaueG Manual_files\";
; Overwrite jar files for SCILAB v5.5.2
Source: "Scilab-5.5.2\org.scilab.modules.renderer.jar"; DestDir: "{code:GetScilabHome}\modules\renderer\jar\"; Check: IsScilab5_5_2;
Source: "Scilab-5.5.2\org.scilab.modules.gui.jar"; DestDir: "{code:GetScilabHome}\modules\gui\jar\"; Check: IsScilab5_5_2;


[Run]
  Filename: "{code:GetScilabHome}\bin\WScilex.exe"; Parameters: "-nb -f ""{code:GetScilabHome}\laueg\\LaueG_setup.sce"""; StatusMsg: "Please wait while SCILAB starts up"


[Code]
var
  ScilabExe: String;
  ScilabVersion: String;
  ScilabHome: String;
  Page: TWizardPage;
  FolderTreeView: TFolderTreeView;
  NewLine: String;

// ============ InitializeSetup() and subroutines ==========

function FindScilabRegistry(): String;
var
  temp: String;
  temp2: String;
begin
// Fall back values for Scilab folder
  temp := 'scilab-5.5.2';
  temp2 := 'C:\Program Files\scilab-5.5.2';
// Get folder for 64 bit Scilab
// Following gives error if not running 64 bit, so first check with IsWin64
  if RegKeyExists(HKLM64, 'SOFTWARE\Scilab') then
  begin
    if RegQueryStringValue(HKLM64, 'SOFTWARE\Scilab', 'LASTINSTALL', temp) then
      RegQueryStringValue(HKLM64, 'SOFTWARE\Scilab\'+temp, 'SCIPATH', temp2);
// Return name (and folder) of main Scilab executable
    result := temp2+'\bin\wscilex.exe';
  end
end;


function FindScilabAssociation(): String;
var
  temp: String;
  temp2: String;
begin
// Fall back values for Scilab icon
  temp := 'Scilab5.sce';
  temp2 := 'C:\Program Files\scilab-5.5.2\bin\wscilex.exe,7';
// Get file association of *.sce files
  RegQueryStringValue(HKCR, '.sce', '', temp);
// Get default icon for Scilab
  RegQueryStringValue(HKCR, temp, 'DefaultIcon', temp2);
// Return name (and folder) of main Scilab executable
  result := ChangeFileExt(temp2,'.exe');
end;


function InitializeSetup(): Boolean;
begin
//
// Convenient place to initialise NewLine
  NewLine := Chr(13)+Chr(10);
//
// If logged in as admin, try to find the Scilab executable
// using the registry or file associations.
  ScilabExe := '';
  if ( IsAdminLoggedOn() ) then
    begin
    ScilabExe := FindScilabRegistry();
    if ( not FileExists(ScilabExe) ) then
      ScilabExe := FindScilabAssociation();
    end;
//
// If executable missing, try guessing "local application data" folders
  if ( not FileExists(ScilabExe) ) then
    ScilabExe := GetEnv('USERPROFILE')+'\AppData\Local\scilab-5.5.2\bin\wscilex.exe';

  result := true;
end;


// ============ InitializeWizard() and subroutines ==========

procedure FolderSelected(Sender: TObject);
begin
  ScilabExe := TFolderTreeView(Sender).Directory + '\bin\wscilex.exe';
end;


function NextButtonClick(PageID: Integer): Boolean;
begin
  result := true;
// Give Welcome popup depending on ScilabExe being found
  if ( PageID = wpWelcome ) then
  begin
    if ( not FileExists(ScilabExe) ) then
      MsgBox('Cannot find SCILAB installation!' +NewLine+NewLine+
             'If SCILAB was not installed "as administrator", you ' +
             'will have to manually find a folder named "scilab-5.5.2", ' +
             'or similar.' +Newline+Newline+ 
             'If SCILAB was installed "as administrator", LaueG must ' +
             'also be installed "as administrator".',
              mbInformation, MB_OK);
  end
// Pressed the 'Next' button, so check if ScilabExe exists.
// Ask to try again, or extract the Scilab path and version number.
  else if ( PageID = Page.ID ) then
    begin
    if ( FileExists(ScilabExe) ) then
    begin
      ScilabHome    := ExpandFileName(ExtractFilePath(ScilabExe)+'\..\');
      ScilabVersion := ExtractFileName(RemoveBackslash(ScilabHome));
      ScilabVersion := Copy(ScilabVersion,8,5);
      if ( ScilabVersion = '5.5.2' ) then
        MsgBox('Make sure SCILAB is not running', mbInformation, MB_OK)
      else
        begin
        MsgBox('LaueG only supports SCILAB version 5.5.2', mbInformation, MB_OK);
        result := false;
        end;
      end
    else
      begin
      MsgBox('SCILAB path is invalid, try again or cancel installation', mbError, MB_OK);
      result := false;
      end;
    end
end;


procedure InitializeWizard();
begin
  Page := CreateCustomPage(wpWelcome, 'Select the SCILAB home folder', 'TFolderTreeView');
  FolderTreeView := TFolderTreeView.Create(Page);
  FolderTreeView.Width := Page.SurfaceWidth;
  FolderTreeView.Height := Page.SurfaceHeight;
  FolderTreeView.Parent := Page.Surface;
  if ( FileExists(ScilabExe) ) then
    FolderTreeView.Directory := ExpandFileName(ExtractFilePath(ScilabExe)+'\..\')
  else
    FolderTreeView.Directory := ExpandConstant('{localappdata}');

  FolderTreeView.OnChange := @FolderSelected;

end;
 

// ========== Routines called in the [Files] section ============

function GetScilabHome(DummyParam: string): String;
begin
  result := ScilabHome;
end;

function IsScilab5_5_2(): Boolean;
// Return true if version 5.5.2
begin
  result := ( ScilabVersion =  '5.5.2' );
end;

procedure MakeShortcut();
var
  Filename: String;
begin
// Set the filename for the shortcut
  Filename := ExpandConstant('{userdesktop}\LaueG.lnk');
// Create the shortcut
  CreateShellLink(
    Filename,
    'LaueG',
    ScilabHome + '\bin\WScilex.exe',
    '-nb -f "' + ScilabHome + 'laueg\laueg.sce"',
    ScilabHome + 'laueg\',
    ScilabHome + '\laueg\laueg_icon.ico',
    0,
    SW_SHOWNORMAL);
// Tell user about the shortcut
  MsgBox('Shortcut to run LaueG created on the desktop of user "' +
          ExpandConstant('{username}')+'"',
          mbInformation,  MB_OK)
end;
