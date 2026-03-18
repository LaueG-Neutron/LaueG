global ImageInfo LatticeInfo OrientInfo ModulateInfo DataDir IniDir
global ArgboxSet DisplaySet InstrumDefault

// ========= Read Tif, Ell and Int files ==============
//function ReadImageHeader(base_name,bDminValid)
//function ReadImageData()
//function [hkl,xy_efh,amult,itype,itwin]=ReadEllFile(file_name)
//function [hkl,imul,wav,tth,xy,counts]=ReadIntFile(file_name)
// ======= Load nominal wavelength file for instrument ========
//function LoadNominalWaveFile()
// ========= Loading/updating INI and Log files ==========
//function ReadIniFile()
//function WriteIniFile()
//function WriteLogFile(text)
//function WriteLogFileHeader(text)
// ========= Saving/loading parameter files ==========
//function LoadSetupsFile()
//function [Linfo,Oinfo,Minfo]=ReadIndexFile(sFile)
//function bOK=LoadSelectedIndexFile(bLog)
//function bOK=LoadIndexFile(sFile,bLog)
//function [levels,bTwins,bModulate]=CheckIndexFile(sFile)
//function SaveIndexFile(sFile,bLog)
//function bOK=LoadSettingsFile()
//function SaveSettingsFile()
// ========= ASCII method for saving/loading parameters to files ==========
//function [lStructs,bOK]=LoadParamFile(sFile,sHeader,sStructs)
//function bOK=SaveParamFile(sFile,sHeader,sStructs)
// ========= Other file related routines ==========
//function ChangeDataDir(dname,bOutput)
//function CloseAllFiles()
//function DeleteTempFile(fname)


// ========= Read Tif, Ell and Int files ==============

function ReadImageHeader(base_name,bDminValid)
// Read image header and update info/data structures from the instrument
// setups. Update the settings and save to the file.
// NB: No longer loads the reduced image into ReducedImage.
// 
// bDminValid=%t if d_min, etc. settings are valid
global ImageInfo
//
// Needed in case ImageInfo is corrupted
  if(or( [~isfield(ImageInfo,"setup"), ~isfield(ImageInfo,"type")])) then
    ImageInfo.setup=0
    ImageInfo.type=0
  end
//
// Store current setup and type
  last_setup=ImageInfo.setup
  last_type=ImageInfo.type
//
// Read the image file base_name. Returns file info
// It no longer loads a reduced image into ReducedImage
  Info=RunImageInfo(base_name,%t)
//
// Load ImageInfo from Info and the corresponding instrument setup.
// Also set LatticeInfo.pixcen for "unstable" VIVALDI & KOALA.
  LoadInstrumSetup(Info)
//
// Add the basename of the image file to ImageInfo
  ImageInfo.basename=base_name
//
// If previous setup was "dummy", try loading DisplaySet & ArgboxSet from
// the settings file, else use the instrument defaults
  if(last_type == 0) then
    if( ~LoadSettingsFile() ) then
       LoadSetupSettings()
       SaveSettingsFile()
    end
// If bDminValid=%t, d_min, etc. settings are valid, so don't change them
// If bDminValid=%f, load default d_min, etc. settings if setup has changed
  elseif( ~bDminValid ) then
    if((last_setup ~= ImageInfo.setup) ) then
      LoadSetupSettings()
      SaveSettingsFile()
    end
  end
//
endfunction


function ReadImageData()
global FileImage
// Close all open files in case any were left open
  CloseAllFiles()
// Open the TIF file, complain and abort if unable to
  file_name=ImageInfo.basename+".tif"
  numx=ImageInfo.numxy(1)
  numy=ImageInfo.numxy(2)
  [fd,ierr]=mopen(file_name,"rb");
  if(ierr ~= 0) then
    AbortBox("Unable to open "+file_name)
  end
// Skip the offset bytes before the image data
  mseek(ImageInfo.offset)
// Read the image as unsigned 16bit integers and convert to a matrix of doubles
  FileImage=double( matrix(mgeti(numx*numy,"us",fd),numx,numy) )
// Close image file
  mclose(fd);
//
endfunction


function [hkl,xy_efh,amult,itype,itwin]=ReadEllFile(file_name)
// Read information from the *.ell file
//
// Open and read the header line
  CloseAllFiles()
  [fd,err]=mopen(file_name,"r")
  if(err ~= 0) then
    AbortBox("Unable to open "+file_name)
  end
// Read header lines, complain and abort if invalid
  header=mgetl(fd,1)
  [nread,iver,iopt]=msscanf(header,"LaueG Ellipse File Version %i, Option %i")
  if(nread ~= 2) then
    AbortBox(file_name+" does not have a valid header line")
  elseif( (iver ~= 1) | (iopt < 1) | (iopt > 2) ) then
    AbortBox(file_name+" has an unknown version or option")
  end
// Skip 2 comment lines
  mgetl(fd,2)
// Read remaining lines into buffer() and close the file
  buffer=mgetl(fd,-1)
  mclose(fd)
//
// Decode buffer() depending on iopt value
  if(iopt == 1) then
    [nread,ih,ik,il,x,y,e,f,h,c1,c2,c3,itype]= ...
      msscanf(-1,buffer,"%i %i %i %f %f %f %f %f %f %f %f %i")
      itwin=zeros(itype)
  else
    [nread,ih,ik,il,x,y,e,f,h,c1,c2,c3,itype,itwin]= ...
      msscanf(-1,buffer,"%i %i %i %f %f %f %f %f %f %f %f %i %i")
  end
// Complain and die if not all lines were decoded
  if(size(buffer,1) ~= size(itwin,1)) then
    AbortBox("Invalid data line in "+file_name)
  end
// Pack the data into matrices for convenience
  hkl=[ih,ik,il]
  xy_efh=[x,y,e,f,h]
  amult=[c1,c2,c3]
endfunction


function [hkl,imul,wav,tth,xy,counts]=ReadIntFile(file_name)
// Read information from the *.int file
//
// Open and read the header line
  CloseAllFiles()
  [fd,err]=mopen(file_name,"r")
  if(err ~= 0) then
    AbortBox("Unable to open "+file_name)
  end
// Read header lines, abort on errors
  header=mgetl(fd,1);
  [nread,iver,iopt]=msscanf(header,"LaueG Intensities File Version %i, Option %i")
  if(nread ~= 2) then
    AbortBox(file_name+" does not have a valid header line")
  end
// Check version and option, abort if invalid
  if(iver ~= 1) then
    AbortBox(file_name+" has an unknown version="+string(iver))
  elseif( (iopt < 1) & (iopt > 3) ) then
    AbortBox(file_name+" is version 1 with unknown option="+string(iopt))
  end
// Skip 2 comment lines
  mgetl(fd,2)
// Read remaining lines into buffer() and close the file
  buffer=mgetl(fd,-1)
  mclose(fd)
//
// Decode buffer() (ignore any errors)
  if(iopt > 1) then
    [nread,ih,ik,il,im,wav,imul,x,y,tth,iCounts,iSig]= ...
      msscanf(-1,buffer,"%4i%4i%4i%4i%7f%3i%8f%7f%8f%7i%6i");
  else
    [nread,ih,ik,il,wav,imul,x,y,tth,iCounts,iSig]= ...
      msscanf(-1,buffer,"%4i%4i%4i%7f%3i%8f%7f%8f%7i%6i");
  end
// Pack the data into matrices for convenience
  if(iopt > 1) then
    hkl=[ih,ik,il,im]
  else
    hkl=[ih,ik,il]
  end
  xy=[x,y]
  counts=[iCounts,iSig]
endfunction


// ======= Load nominal wavelength file for instrument ========

function LoadNominalWaveFile()
//
// Open wavs.lis
  CloseAllFiles()
  [fd,ierr]=mopen(ProgDir+"wavs.lis","r")
  if(ierr ~= 0) then
    AbortBox(["Unable to find ""wavs.lis"" file"; ...
            "Has LaueG been correctly installed?"])
  end
//
// Read list of wavelength files and comments
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
  end
//
// Copy the wavelength file and log to the console
  copyfile(names(ival(1)),"laue4_wav.dat")
  printf( "Loaded wavelength file: %s\n",names(ival(1)) )
  WriteLogFile( "Loaded wavelength file: "+names(ival(1)) )
//
endfunction


// ========= Loading/updating INI and Log files ==========

function ReadIniFile()
global DataDir InstrumDefault
// Read DataDir from the laueg.ini file
//
// Try to open the laueg.ini file
  [fd,ierr]=mopen(IniDir+"laueg.ini","r")
  if(ierr == 0) then
// Read the directory location from laueg.ini
    DataDir=mgetl(fd,1)
// Read the default instrument setup
// If missing, use zero.
    [n,InstrumDefault]=mfscanf(fd,"%d")
    if(n ~= 1) then
      InstrumDefault=0
    end
    mclose(fd)
  else
// If no ini file, try to make a new one
// Can't use WriteIniFile() here as it calls CloseAllFiles()
// No ini file, so give warning and set DataDir to C:\ or ~
    mprintf("Unable to find ""laueg.ini"" file, assuming new installation\n")
    DataDir="~"
    if(getos() == "Windows") then 
        DataDir="C:\"
    end
// Make sure directory exists (its safe)
    mkdir(IniDir)
// Try to create a new ini file, abort if it fails
    [fd,ierr]=mopen(IniDir+"laueg.ini","w")
    if(ierr ~= 0) then
      AbortBox("Cannot create a new ""laueg.ini"" file")
    end
// Write DataDir to ini file
    mfprintf(fd,"%s\n",DataDir)
// Set the default setup and write to the ini file
    InstrumDefault=0
    mfprintf(fd,"%d\n",InstrumDefault)
    mclose(fd)
  end
//
endfunction


function WriteIniFile()
// Write DataDir to the ini file
  CloseAllFiles()
// Try to open the ini file and write DataDir to it
  ini_name=IniDir+"laueg.ini"
  [fd,ierr]=mopen(ini_name,"w")
  if(ierr ~= 0) then
    AbortBox("Cannot open the ""laueg.ini"" file")
  end
  mfprintf(fd,"%s\n",DataDir)
  mfprintf(fd,"%d\n",InstrumDefault)
  mclose(fd)
endfunction


function WriteLogFile(text)
// Open the log file, complain if we can't
  [fd,ierr]=mopen("laueg.log","a")
  if(ierr ~= 0) then
    WarnBox("Unable to open ""laueg.log""")
    return
  end
// Output the line of text to the log file
  mfprintf(fd,"%s\n",text)
  mclose(fd)
endfunction


function WriteLogFileHeader(text)
// Create string with the current time
  t=clock()
  time=msprintf(" %.2i:%.2i:%.2i",t(4:6))
// Write the header text and date/time to the log file
  WriteLogFile("==== "+text+" ("+date()+time+") ====")
endfunction


// ========= Saving/loading parameter files ==========

function LoadSetupsFile()
global InstrumSetups
//
// Try to load named structures from sFile
  sFile=ProgDir+"setups.dat"
  sStructs="InstrumSetups"
  [lStructs,bOK]=LoadParamFile(sFile,"LaueG Instrument Setups File",sStructs)
//
// Complain and abort on failure
  if( ~bOK ) then
    AbortBox("Invalid instrument setups file")
  end
//
// Finally, copy setups to the global variable
  InstrumSetups=lStructs(1)
//
endfunction


function [Linfo,Oinfo,Minfo]=ReadIndexFile(sFile)
//
// Try to load named structures from sFile
  sStructs=["LatticeInfo","OrientInfo","ModulateInfo"]
  [lStructs,bOK]=LoadParamFile(sFile,"LaueG Index File",sStructs)
//
// If unable to load file, return with empty structures
  if( ~bOK )
    Linfo=[]
    Oinfo=[]
    Minfo=[]
    return
  end
//
// Copy structures to return values
  Linfo=lStructs(1)
  Oinfo=lStructs(2)
  Minfo=lStructs(3)
//
// For early index files without ModulateInfo or special_pairs
  if( Minfo == [] )
    Minfo.vecs=[]
    Minfo.mults=[]
    Minfo.d_min=1
    Minfo.d_max=10
  end
  if( ~isfield(Linfo,"special_pairs") )
    Linfo.special_pairs=[]
  end
//
endfunction


function bOK=LoadSelectedIndexFile(bLog)
global LatticeInfo
//
// Ask for the file name, and return if an empty string
// DO NOT change DataDir to this new path
  fname=uigetfile(["*.idx"],DataDir, "Select index file", %f)
  bOK= ~isempty(fname)
  if( ~bOK ) then return; end
//
// Load globals from the index file
  if( ~LoadIndexFile(fname,bLog) ) then
    AbortBox("Index file "+fname+" is missing or corrupt")
  end
//
// Demote any level=3 in batch mode, or if idx and tif base names differ
  idx_bname=pathconvert(dirname(fname))+basename(fname)
  tif_bname=pathconvert(DataDir)+ImageInfo.basename
  if( (BatchFiles ~= []) | (idx_bname ~= tif_bname) ) then
    for i=1:size(LatticeInfo,1)
      LatticeInfo(i).level=min(2,LatticeInfo(i).level)
    end
  end
//
endfunction


function bOK=LoadIndexFile(sFile,bLog)
global LatticeInfo OrientInfo ModulateInfo
//
// Try to read index file, return with bOK=%f if failed
  [Linfo,Oinfo,Minfo]=ReadIndexFile(sFile)
  if(Linfo == []) then
    bOK=%f
    return
  end
//
// Copy values to globals and update return status
// Retain any twin information in LatticeInfo[]
  LatticeInfo(1)=Linfo
  OrientInfo=Oinfo
  ModulateInfo=Minfo
  bOK=%t
//
// Update log file, if bLog is set
  if( bLog ) then
    WriteLogFile("Loading Index File: "+sFile)
  end
//
endfunction


function [levels,bTwins,bModulate]=CheckIndexFile(sFile)
//
// Check if the index file has twins or is modulated, and
// return the index levels for LatticeInfo(1:3).
//
// Try to read index file, abort if failed
  [Linfo,Oinfo,Minfo]=ReadIndexFile(sFile)
  if(Linfo == []) then
    AbortBox("Cannot read index file = "+sFile)
  end
//
// Copy Linfo[].level values to levels[]
  levels=[]
  for i=1:size(Linfo,1)
    levels(i)=Linfo(i).level
  end
//
// Check if any modulation information
  if( Minfo == [] ) then
    bModulate=%f
  else
    bModulate= ~isempty(Minfo.vecs)
  end
//
// If more than one lattice, set bTwins
  bTwins=( size(levels,1) > 1 )
//
endfunction


function SaveIndexFile(sFile,bLog)
//
// Try to load named structures from sFile
  sStructs=["LatticeInfo","OrientInfo","ModulateInfo"]
  bOK=SaveParamFile(sFile,"LaueG Index File",sStructs)
//
// Output info about success or failure
  if( bOK )
    if( bLog ) then
      WriteLogFile("Saving Index File: "+sFile)
    end
  else
    AbortBox("Error writing to index file: "+sFile)
  end
//
endfunction


function bOK=LoadSettingsFile()
global ArgboxSet DisplaySet
//
// Try to load named structures from sFile
  sFile=DataDir+"laueg.set"
  sStructs=["ArgboxSet","DisplaySet"]
  [lStructs,bOK]=LoadParamFile(sFile,"LaueG Settings File",sStructs)
//
// If unable to load file, simply return
  if( ~bOK )
    return
  end
//
// Update global structures from the output list
  ArgboxSet=lStructs(1)
  DisplaySet=lStructs(2)
//
// Update log file
  WriteLogFile("Loading Settings File: "+sFile)
//
endfunction


function SaveSettingsFile()
//
// Try to load named structures from sFile
  sFile=DataDir+"laueg.set"
  sStructs=["ArgboxSet","DisplaySet"]
  bOK=SaveParamFile(sFile,"LaueG Settings File",sStructs)
//
// Output info about success or failure
  if( bOK )
    WriteLogFile("Saving Settings File: "+sFile)
  else
    ErrorBox("Error writing to settings file: "+sFile)
  end
//
endfunction



// ========= ASCII method for saving/loading structures to files ==========

// MERDE! MERDE! MERDE!
//
// I made a bad move ditching my old ASCII method of saving Settings and
// and Index information in files. SCILAB keeps changing the format for
// their SAVE files, so I am reverting to a (new) ASCII method. This
// also allows the instrument setups file to be in the same format.


function [lStructs,bOK]=LoadParamFile(sFile,sHeader,sStructs)
// lStructs is a list of structures named in vector sStructs[]
//
// Set bOK to false
  bOK=%f
// Remove twins from local copy of LatticeInfo (in case it is used)
  LatticeInfo=LatticeInfo(1)
// Save structures named in sStructs to output list
  execstr( "lStructs=list("+strcat(sStructs,",")+")" )
//
// Read parameter file to determine which type it is
  [fd,err]=mopen(sFile,"r")
  if(err ~= 0) then
    return
  end
// Check header for my LaueG format or SCILAB's HDF format
  str=mgetl(fd,1)
  bLaueG=(str == "// "+sHeader)
  bHDF=(str == "‰HDF")
// If my format, read the remainder of the file into str[]
  if( bLaueG ) then
    str=mgetl(fd,-1)
  end
  mclose(fd)
//
// "Load" parameter file depending on its type
// If my format, evaluate the string
  if( bLaueG ) then
    if ( execstr(str,"errcatch") ~= 0 ) then 
      AbortBox("Cannot decode Parameter File = "+sFile)
    end
// If a HDF file, use load() to read the file
  elseif( bHDF ) then
    if ( execstr( "load(""" +sFile+ """)", "errcatch") ~= 0 ) then 
      AbortBox("Cannot load Parameter File = "+sFile)
    end
// Else, give up
  else
    AbortBox("Unknown format for Parameter File = "+sFile)
 end
//
// Ensure all fields in saved structures also exist in loaded ones
  for i=1:size(sStructs,"*")
// Find field names in the old structure that are not in the new
    fold=fieldnames(lStructs(i))
    fnew=fieldnames(evstr(sStructs(i)))
    fdiff=setdiff(fold,fnew)
// Copy any "difference" fields from old to new structure
    if(fdiff ~= []) then
      execstr(sStructs(i)+"."+fdiff+"=lStructs(i)."+fdiff)
    end
  end
//
// Copy local structures to make new output list, and set bOK to true
  execstr( "lStructs=list("+strcat(sStructs,",")+")" )
  bOK=%t
//
endfunction


function bOK=SaveParamFile(sFile,sHeader,sStructs)
// Save structures named in sStructs to an ascii save file, sFile,
// with a header line of sHeader.
//
// Convert structures in sStructs to a vector of strings
//
  sLines=[]
//
// Loop through structures named in sStructs
  for sName=sStructs
    Struct=evstr(sName)
//
// Loop through elements of the structure
    for i1=1:size(Struct,1)
      for i2=1:size(Struct,2)
//
// Loop through field names
        sFields=fieldnames(Struct)
        for sField=sFields'
//
// Add string expression for field to end of vector of strings
          sPar=msprintf("%s(%d,%d).%s",sName,i1,i2,sField)
          Par=evstr(sPar)
          sLines($+1)=sci2exp(Par,sPar)
//
        end
      end
    end
  end
//
// Output results to sFile
//
// Open output file
  [fd,err]=mopen(sFile,"w")
  if(err ~= 0) then
    bOK=%f
    return
  end
// Write header line, plus blank line, and vector of strings
  mputl("// "+sHeader,fd)
  mputl("//",fd)
  bOK=mputl(sLines,fd)
// Close file
  mclose(fd)
//
endfunction


// ========= Other file related routines ==========

function ChangeDataDir(dname,bOutput)
global DataDir
//
// Bomb out if trying to use a UNC path name in Windows
  if (getos() == "Windows") then
    if( grep(dname,"\\") ~= [] ) then
      AbortBox(["Unable to use UNC path names such as "+dname;""; ...
      "Copy files to a local drive, or map the drive to a letter (i.e. Z:)"])
    end
  end
//
// Try to change to the new directory
  if( ~chdir(dname) ) then
    ErrorBox("Unable to change to the requested directory")
// If it fails change to the root directory
    if (getos() == "Windows") then
      dname="C:\"
    else
      dname="/"
    end
    if( ~chdir(dname) ) then
      return
    end
  end
// Update DataDir
  DataDir=pathconvert(dname)
// Revert to deleting LaueG temporary files (in this directory)
  mdelete("___laueg_delete_files.in")
// If requested, write directory name to the console
  if( bOutput ) then
    mprintf("Working directory: %s\n",DataDir)
  end
endfunction


function CloseAllFiles()
// Closing the 8 lowest fd's should do it
  sWarn=warning("query")
  warning("off")
  for i=[1:4,7:10]
    mclose(i)
  end
  warning(sWarn)
endfunction


function DeleteTempFile(fname)
// Delete file if "___laueg_delete_files.in" doesn't exist
  if( isempty(fileinfo("___laueg_delete_files.in")) ) then
    mdelete(fname)
  end
endfunction
