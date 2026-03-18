mode(-1)


  //////////////////////////////////////////////////
  //                                              //
  //   This script is not for manual execution.   //
  //   You should close this window and follow    //
  //   the installation instructions.             //
  //                                              //
  //                                              //
  //////////////////////////////////////////////////



// Set default language to English
setlanguage('en_US')
if( getos() == 'Windows' ) then
  warning('off')
  setdefaultlanguage('en_US')
  warning('on')
end

// Move to the directory with the ini file
cd('SCIHOME\..\')
mkdir('laueg')
cd('laueg')

// If no laueg.ini file, create one, complain and die if it fails
if( ~isfile('laueg.ini') ) then
  [fd,ierr]=mopen('laueg.ini','w')
  if(ierr ~= 0) then
    messagebox(['Cannot create ""laueg.ini"" file in ""'+pwd()+'""', ...
              '', 'See readme file'], ...
          "LaueG Install Error", "error", "modal")
    quit
  end
  DataDir=pathconvert(cd('~'))
  mfprintf(fd,'%s\n',DataDir)
  mfprintf(fd,'%d\n',0)
  mclose(fd)
end

// Give manual instructions to remove side-panels
messagebox(['LaueG Windows Layout','', ...
            'Remove the  ""File Browser"", ""Variable Browser"" and', ...
            '""Command History"" side panels from the SCILAB window.', ...
            'This will leave a single ""Scilab Console"" panel.'], ...
                     "LaueG Install", "information", "continue", "modal")

// Give manual instruction to quit
mprintf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
mprintf("   ######################################################\n")
mprintf("   #                                                    #\n")
mprintf("   #   Exit SCILAB once side panels have been removed   #\n")
mprintf("   #                                                    #\n")
mprintf("   ######################################################\n\n\n")
