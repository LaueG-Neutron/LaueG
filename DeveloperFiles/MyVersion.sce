//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
//                                                                     //
//  Given x,y on the perimeter of an ellipse:                          //
//    Calculate r0=sqrt(x.^2 + y.^2), ang0=atand(y,x)                  //
//    Fourier decompose r0.^-2 to A + B*cosd(2*ang0) + C*sind(2*ang0)  //
//    Then e=A+B, f=A-B, h=C for equation e*x^2 + f*y^2 + 2*h*x*y = 1  //
//                                                                     //
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//


/////////=============================================
function junko()
/////////=============================================

ub1=[ ...
      0.100306   0.000095   0.023973
     -0.067687   0.130657  -0.004215
     -0.051128  -0.013525   0.064028 ]
ub2=[ ...
      0.010006   0.109973  -0.023909
      0.130149  -0.068417   0.004317
      0.014418  -0.023353  -0.063945 ]


ub1_ub2=ub1\ub2

ub1_ub2=ub1_ub2*sign(det(ub1_ub2))/abs(det(ub1_ub2))^(1/3)


ub1_ub2=round(ub1_ub2*12)/12
ub1-(ub2*ub1_ub2)


[H,K,L]=ndgrid(0:3,0:3,0:3)
hkls=[matrix(H,-1),matrix(K,-1),matrix(L,-1)]
hkls(1,:)=[]



ub1*hkls
ub1_ub2*hkls'


//
//The SCILAB functions
//    [Hvec]=PixWav2Hvec(pixwav)
//    [pixwav]=Hvec2PixWav(Hvec)
//readily convert between
//    pixwav=[ 239.     323.    1.3]
//    pixwav=[1239.    3082.    1.8]
//and
//    Hvec=[-0.4812083  -0.4926137  -1.1152420]
//    Hvec=[-0.1603395   0.5173230  -0.4389141]
//using these parameters
    ImageInfo.type= 1
    ImageInfo.phi= -89.989998
    OrientInfo.pixcen= [1986.932,984.26001]
    OrientInfo.pixsize= [0.200174,0.199388]
    OrientInfo.pixskew= -0.000434
    OrientInfo.xtaloff= [-0.227274,0,-0.415665]
    OrientInfo.drumoff= [-0.11,0,0.14]
    OrientInfo.beamvert= 0.169727
    OrientInfo.drumrad= 159.16
//LatticeInfo.cell= [12.838697,12.838697,12.838697,51.1721,51.1721,51.1721]
//LatticeInfo.ub=[...
//    0.0813219  - 0.0332875    0.0213319  
//    0.0319983  - 0.1007455    0.014934   
//    0.0640528    0.0219589  - 0.1051753  ]
//
/////////////////////////////

Hvec=[-0.4812083  -0.4926137  -1.1152420]
pixwav=Hvec2PixWav(Hvec)

OrientInfo.drumoff=OrientInfo.drumoff+[1,0,0];
shiftx=Hvec2PixWav(Hvec) - pixwav
OrientInfo.drumoff=OrientInfo.drumoff-[1,0,0];

OrientInfo.drumoff=OrientInfo.drumoff+[0,1,0];
shifty=Hvec2PixWav(Hvec) - pixwav
OrientInfo.drumoff=OrientInfo.drumoff-[0,1,0];

OrientInfo.drumoff=OrientInfo.drumoff+[0,0,1];
shiftz=Hvec2PixWav(Hvec) - pixwav
OrientInfo.drumoff=OrientInfo.drumoff-[0,0,1];

shiftx(1:2)- OffsetShift(pixwav(1:2),[1,0,0])
shifty(1:2)- OffsetShift(pixwav(1:2),[0,1,0])
shiftz(1:2)- OffsetShift(pixwav(1:2),[0,0,1])

//
//
// Ellipse formula: e*x^2 + f*y^2 + 2*h*x*y = 1
efh=[1,2,0.9]

// Create a list of X,Y within the ellipse given by efh[]
dxy=0.001;
g=efh(1)*efh(2) - efh(3)^2;
ylo = -sqrt( efh(1)/g );
yhi = +sqrt( efh(1)/g );
xys=[0,0];
for y=ylo:dxy:yhi
  xlo = ( -efh(3)*y - sqrt(efh(1) - g*y^2) )/efh(1);
  xhi = ( -efh(3)*y + sqrt(efh(1) - g*y^2) )/efh(1);
//  xys2=[xlo:dxy:xhi]';
//  xys2(:,2)=y;
//  xys=[xys;xys2];
  xys=[xys;[xlo,y];[xhi,y]];
end
xys(1,:)=[];

efh2=xy_2_efh(xys)
////========================


dxyz=0.01

mx=[1,0,0; 0,2,0; 0,0,3]
Rot=CalcXRot(12.8) * CalcYRot(42.0) * CalcZRot(67.3)
mx = Rot' * mx * Rot

xyzs=Make3DVolume(mx,dxyz);
xyzs=Make3DSurface(mx,dxyz);

ssq=sum( xyzs .* (mx * xyzs')' ,'c');
tabul(int(ssq*1e3)/1e3)



/////////=============================================
endfunction
/////////=============================================



function xyzs=Make3DVolume(mx,dxyz)
//
  M11=mx(1,1)
  M22=mx(2,2)
  M33=mx(3,3)
  M12=mx(1,2)
  M13=mx(1,3)
  M23=mx(2,3)
//
  xyzs=[0,0,0]
//
  A3 = M11*M22*M33 -M11*M23^2 -M22*M13^2 -M33*M12^2 +2*M23*M12*M13
  C3 = M11*M22 - M12^2
  dz=sqrt(max(0, C3/A3 ))
  nz=floor(dz/dxyz)
  for z=[-1:1/nz:+1]*dz
//
    xys=[0,0]
//
    A2=M11*M22 - M12^2
    B2=(M11*M23 - M12*M13)*z
    C2=M11 - (M11*M33 - M13^2)*z^2
    dy=sqrt(max(0, (B2/A2)^2 + (C2/A2) ))
    ny=floor(dy/dxyz)
    y_tmp=[0]
    if(ny > 0) then y_tmp= [-1:1/ny:+1]*dy; end
    for y= -(B2/A2) + y_tmp
//
      A1=M11
      B1=M12*y + M13*z
      C1=1 - M22*y^2 - M33*z^2 - 2*M23*y*z
	  dx=sqrt(max(0, (B1/A1)^2 + (C1/A1) ))
      nx=floor(dx/dxyz)
      x_tmp=[0]
      if(nx > 0) then x_tmp= [-1:1/nx:+1]*dx; end
      xs = -(B1/A1) + x_tmp';
      xs(:,2)=y;
	  xys=[xys; xs];
//
    end
    xys(1,:)=[]
//
    xys(:,3)=z
    xyzs=[xyzs; xys]
//
  end
  xyzs(1,:)=[];
//
  return
endfunction


function xyzs=Make3DSurface(mx,dxyz)
//
  M11=mx(1,1)
  M22=mx(2,2)
  M33=mx(3,3)
  M12=mx(1,2)
  M13=mx(1,3)
  M23=mx(2,3)
//
  xyzs=[0,0,0]
//
  A3 = M11*M22*M33 -M11*M23^2 -M22*M13^2 -M33*M12^2 +2*M23*M12*M13
  C3 = M11*M22 - M12^2
  dz=sqrt(max(0, C3/A3 ))
  nz=floor(dz/dxyz)
  for z=[-1:1/nz:+1]*dz
//
    xys=[0,0]
//
    A2=M11*M22 - M12^2
    B2=(M11*M23 - M12*M13)*z
    C2=M11 - (M11*M33 - M13^2)*z^2
    dy=sqrt(max(0, (B2/A2)^2 + (C2/A2) ))
    ny=floor(dy/dxyz)
    y_tmp=[0]
    if(ny > 0) then y_tmp= [-1:1/ny:+1]*dy; end
    for y= -(B2/A2) + y_tmp
//
      A1=M11
      B1=M12*y + M13*z
      C1=1 - M22*y^2 - M33*z^2 - 2*M23*y*z
	  dx=sqrt(max(0, (B1/A1)^2 + (C1/A1) ))
      xs = -(B1/A1) - dx
      if(dx > 0) then xs(2) = -(B1/A1) + dx; end
//
      xs(:,2)=y
	  xys=[xys; xs]
//
    end
    xys(1,:)=[]
//
    xys(:,3)=z
    xyzs=[xyzs; xys]
//
  end
  xyzs(1,:)=[];
//
  return
endfunction


function efh=xy_2_efh(xys)
//
// Convert x,y to polar coords
  ang=(180/%pi)*atan(xys(:,2),xys(:,1));
  rsq=xys(:,1).^2 + xys(:,2).^2;
//
// Bin polar coords in groups of increasing angle
  nxy=size(xys,1);
  nang=ceil( sqrt(nxy) );
  angs=linspace(-180,180,nang);
  ibin=dsearch(ang,angs);
//
// From each bin find perimeter by choosing the point with the largest radius
  perim=[]
  for ib=unique(ibin)'
    i=find(ibin == ib);
    [v,k]=max(rsq(i));
    perim=[perim; ang(i(k)) , rsq(i(k))];
  end
//
// Convert perim[] to ang[] & rinvsq[]
  ang=perim(:,1)
  rinvsq=perim(:,2) .^ (-1.0)
//
// Integrate values for approx. Fourier decomp. to an ellipse
// Inaccurate due to uneven ang[] spacing
  IR=mean( rinvsq );
  IRCos2=mean( cosd(2*ang) .* rinvsq ) / mean( cosd(2*ang).^2 );
  IRSin2=mean( sind(2*ang) .* rinvsq ) / mean( sind(2*ang).^2 );
//
// Correct Fourier components by scaling sin & cos parts
  calc0=IRCos2 .* cosd(2*ang) + IRSin2 .* sind(2*ang);
  IR=mean(rinvsq - calc0)
  rfac=sum(abs( rinvsq - IR)) / sum(abs( calc0 ))
  IRCos2=IRCos2*rfac
  IRSin2=IRSin2*rfac
//
//  calc=IR + IRCos2 .* cosd(2*ang) + IRSin2 .* sind(2*ang);
//  clf; plot(ang,[rinvsq,calc])
//
// Convert to e,f,h values
  efh(1) = 1/(IR+IRCos2)^2
  efh(2) = 1/(IR-IRCos2)^2
  efh(3) = 1/(IR+IRSin2)^2 - (efh(1)+efh(2))/2
//
  return
endfunction


function [dxy_pix]=OffsetShift(pix,offset)
// Shift in pixel position of a spot due to sample offset (in mm)
// Sample offset is in laboratory coords, i.e. after PHI rotation 
//
// Spot position in mm relative to the pixel center
  xy_mm=( pix - OrientInfo.pixcen ) .* OrientInfo.pixsize
// Angle that x pixel goes around the drum
  angle=xy_mm(1)/OrientInfo.drumrad
// Shift in mm of x (= angle_shift * OrientInfo.drumrad) of spot due to offset
  dx_mm=offset(3)*sin(angle) - offset(1)*cos(angle)
// Shift in height of spot
  dy_mm=-offset(2) + ( offset(1)*sin(angle) + offset(3)*cos(angle) )*xy_mm(2)/OrientInfo.drumrad
// Difference in terms of pixel size
  dxy_pix=[dx_mm , dy_mm] ./ OrientInfo.pixsize
//
endfunction


function [Dpix]=OffsetShift_OLD(pix,offset)
// Shift in position for a spot at pix due to sample offset (in mm)
// Sample offset is in laboratory coords, i.e. after PHI rotation 
//
// Spot position in mm relative to the pixel center
  dxy=( pix - OrientInfo.pixcen ) .* OrientInfo.pixsize
// Diffracted beam vector from the offset sample to the pixel position
  angle=dxy(1)/OrientInfo.drumrad
  S1=OrientInfo.drumrad*[sin(angle),0,cos(angle)] + [0,dxy(2),0]
  S1=S1-offset
// Recalculate spot position for same beam direction but no sample offset
  xyz=OrientInfo.drumrad * S1 / norm( S1([1,3]) )
// Convert to spot position in mm around the IP cylinder
  angle=atan(xyz(1),xyz(3))
  dxy2=[OrientInfo.drumrad*angle , xyz(2)]
// Return difference in terms of pixel size
  Dpix=(dxy2 -dxy) ./ OrientInfo.pixsize
//
endfunction



// LaueG globals
global ArgboxSet BatchFiles DataDir DisplayData DisplayMode
global DisplaySet FlickData ImageFigure ImageInfo IniDir
global InstrumDefault InstrumSetups IntegRectSave LastEvent
global LatticeInfo MainAxes MainImage MainLines MainMarks
global ModulateInfo OrientInfo ProgBarId ProgDir RawImage
global ReducedImage StripImage RejectsInfo ZoomAxes
global bDevelopMode
// My globals
global LaueGFiles


// ==================  Routines I only want in my version  ===============
//function SetupDevelMenu()
//function MyStuff_Menu(iMenu)
//function SaveSourceFiles()
//function sCode=CreateReleaseCode()
//function SetupBuildFolder()
//
// ======================= Experimental Functions =======================
//function UBOptim() & more appear at the end of this file
//
// ===================  Program Notes to Myself  ================
//
//// The main variable structures are:
// ImageInfo      Information from the TIFF header and values that
//                are set when the image is first opened, i.e. edges.
// OrientInfo     Information needed to orient the sample and are
//                common to all lattices.
// LatticeInfo()  An array where each element describes the information
//                needed to orient individual lattices.
// LatticeInfo.level = 0 no cell info, 1 no UB, 2 rough UB, 3 refined UB
// ModulateInfo   Information on modulations to the lattice
// RejectsInfo    Information specific to the reject spots
// ArgboxSet      User set parameters for argonne_boxes
// DisplaySet     User set parameters to display images
// DisplayMode    Options and modes used to display images
// DisplayData    Arrays of data used to display images, i.e.
//                found_spots, marked_spots, ell_mults, ell_data
// FlickData      Information, double() & 8-bit images for the flick files
// InstrumSetups  Vector of instruments setups read from "setups.dat"
// BatchFiles     Vector of base file names for batch processing
// LaueGFiles     Vector of source file names (unused in release version)
//
//// SCILAB structures created for displaying the main and zoom images:
//   ImageFigure   MainAxes     MainImage     MainRects    MainArcs
//   MainLines     MainMarks    ZoomAxes      ZoomMarks
//
//// Large matices containing intensity data:
// FileImage       double(), read directly from image file
// StripImage      double(), as RawImage but stripped of background
// RawImage        double(), created from FileImage or StripImage for display
// ReducedImage    double(), RawImage averaged over ImageInfo.reduce^2 pixels
// NB: FlickData.images are ReducedImage for flick images
// NB: FlickData.pixmaps are MainImage.data for flick images
// NB: StripImage[] is used for orienting twin & split spots
//
//// Other global variables:
// IniDir DataDir ProgDir ProgBarId
//
//// Under development:
// IntegRectSave   Used to reject spots in a rectangular area
// InstrumDefault  Default setup for ambiguous cases


// ==================  Routines I only want in my version  ===============

function SetupDevelMenu()
//
// Add the MyStuff menu
  delmenu("MyStuff")
  AddConsoleMenu("MyStuff", ...
                  ["Edit Code","Save Source","Create Release Code", ...
                   "Generate Build Folder","Edit Graphics","Keep Temp Files"], ...
                 "DevelMenu")
//
endfunction


function DevelMenu(iMenu)
//
// Start Scinotes with *.sce files loaded
  if(iMenu == 1) then
    editor(ProgDir+LaueGFiles)
  end
//
// Save source files to an archive
  if(iMenu == 2) then
    SaveSourceFiles()
  end
//
// Compile functions
// Save compiled (stripped down) functions to Release.sce in ProgDir.
// NB: Tried writing the compiled functions but SCILAB failed.
  if(iMenu == 3) then
    CreateReleaseCode()
  end
//
// Copy files to Build folder ready to run InnoSetup
  if(iMenu == 4) then
    SetupBuildFolder()
  end
//
// Turn on graphics editor
  if(iMenu == 5) then
    if IsImageOn() then
      ged(8,0)
    end
  end
//
// Stop temporary files being deleted in current directory
// Create "___laueg_delete_files.in" in DataDir
  if(iMenu == 6) then
    fd=mopen("___laueg_delete_files.in","wt")
    mfprintf(fd,".false.\n")
    mclose(fd)
  end
//
endfunction


function SaveSourceFiles()
//
// Copies all source files and other (mainly text) files needed to
// reconstruct LaueG in case of a disaster. Files and folders are
// copied to D:\Work\LaueG\_Release\Source\ and a dated 7zip file
// is then added to D:\Work\LaueG_Backup\
//
  RelDir=ProgDir+"_Release\"
  SrcDir=RelDir+"Source\"
  ArcDir="C:\DevApps\LaueG_Backup\"
//
// Remove old Source folder, if it exists
  if( isdir(SrcDir) ) then
    if( ~removedir( SrcDir ) ) then
      mprintf("Unable to delete old source folder %s\n",SrcDir)
      return
    end
  end
//
// Create a new Source folder
  sleep(100)        // Apparently needed by Windows
  if( ~createdir( part(SrcDir,1:$-1) ) ) then
    mprintf("Unable to create new source folder %s\n",SrcDir)
    return
  end
  mprintf("\nCreated new source code folder %s\n",SrcDir)
//
// Copy source files in ProgDir to the Source folder
  sTypes=["*.sce","*.bat","*.dat","*.lis"];
  for sType=sTypes;
    sFiles=ls(ProgDir+sType);
    for fnam=sFiles';
      if( fnam == ProgDir+"Release.sce" ) then
        continue
      end
      [iOK,message]=copyfile(fnam,SrcDir);
      if(iOK ~= 1) then
        mprintf("Error copying files from %s to %s\n",ProgDir,SrcDir)
        return
      end
    end
    mprintf("Copied %d %s files from %s\n",size(sFiles,1),sType,ProgDir)
  end
//
// Copy subfolders with *.for and *.vfproj files to Source subfolders
  S=dir(ProgDir)
  sDirs=S.name( find(S.isdir) )
//
  nDir=0
  nFor=0
  nProj=0
  for dnam=sDirs'
    dnam=dnam+"\"
// Always ignore _Release folder
    if( dnam == "_Release\" ) then
      continue
    end
// Create subfolder corresponding to ProgDir in SrcDir
    if( ~createdir( part(SrcDir+dnam,1:$-1) ) ) then
      mprintf("Unable to create folder %s\n",SrcDir+dnam)
      return
    end
    nDir=nDir+1
// Copy *.for,*.vfproj files from subfolder in ProgDir to SrcDir
    sFiles=[ls(ProgDir+dnam+"*.for"); ls(ProgDir+dnam+"*.vfproj")]
    for fnam=sFiles'
      iOK=copyfile(fnam,SrcDir+dnam);
      if(iOK ~= 1) then
        mprintf("Error copying files from %s to %s\n",ProgDir,SrcDir+dnam)
        return
      end

      if( grep(fnam,".vfproj") ) then
        nProj=nProj+1
      else
        nFor=nFor+1
      end
    end
  end
  mprintf("Copied %d *.for, %d *.vproj files from %d %s subfolders\n", ...
                                  nFor,nProj,nDir,ProgDir)
//
// Copy all files and subfolders from _Release\BuildFiles\
  iOK=copyfile(RelDir+"BuildFiles\",SrcDir+"BuildFiles\")
  if(iOK ~= 1) then
    mprintf("Error copying files from %s to %s\n", ...
                  RelDir+"BuildFiles\",SrcDir+"BuildFiles\")
    return
  end
  mprintf("Copied all files and subfolders from %s\n",RelDir+"BuildFiles\")
//
// Copy _Notes and InnoSetup files to BuildFiles\
  iOK1=copyfile(RelDir+"_Notes.txt",SrcDir+"BuildFiles\")
  iOK2=copyfile(RelDir+"WindowsInstall.iss",SrcDir+"BuildFiles\")
  if( (iOK1 ~= 1) | (iOK2 ~= 1)) then
    mprintf("Error copying _Notes.txt or WindowInstall.iss to %s\n", ...
                                  SrcDir+"BuildFiles\")
    return
  end
  mprintf("Copied _Notes.txt and WindowInstall.iss to %s\n", ...
                                  SrcDir+"BuildFiles\")
//
// Try to create .7z archive of SrcDir
  sRun="""C:\Program Files\7-Zip\7z.exe"" a"
  sDest=ArcDir+"Source\Source_"+strsubst(date(),"-","")+".7z"
  sLog="/DevApps/LaueG/_Release/Source/zip.log"
  sDos=strcat([sRun,sDest,SrcDir]," ") + " > "+sLog
  mdelete(sLog)
  dos(sDos)
  fd=mopen(sLog,"r")
  sOut=mgetl(fd)
  mclose(fd)
  mdelete(sLog)

  bFail=isempty(grep(sOut,"Everything is Ok"))
  if( bFail ) then
    mprintf("Error creating 7zip archive of Source files\n")
  else
    mprintf("Source files archived to %s\n\n",sDest)
  end
//
endfunction


function CreateReleaseCode()
//
// Read source files and extract all functions and globals.
// Write stripped code lines to Release.sce in current folder.
//
// Sanity check on source files to ignore
  if(LaueGFiles(1) ~= "LaueG.sce") then
    AbortBox("LaueG.sce is not first file in LaueGFiles")
  elseif(LaueGFiles($) ~= "MyVersion.sce")
    AbortBox("MyVersion.sce is not last file in LaueGFiles")
  end
//
// Copy code from the source files to sCode[] (not LaueG,MyPopups,MyVersion)
  mprintf("\nSTART copying SCILAB code\n")
  sCode=[]
  for fname=LaueGFiles(2:$-1)
    [fd,ierr]=mopen(ProgDir+fname,"r");
    txt=mgetl(fd);
    sCode=[sCode;txt];
    mclose(fd);
  end
  mprintf("%d lines of code read from %d source files\n", ...
                         size(sCode,1),size(LaueGFiles,"*")-2)
//
// Blank out all comments
// Find lines containing //
  for iline=grep(sCode,"//")
// Split line at //, ", or '
    [res1,res2]=strsplit(sCode(iline),["//","""","''"]);
// Find which splits are due to //
    icomm=find(res2 == "//");
// If first // split is within a string constant:
    if( modulo(icomm(1),2) == 0 ) then
// If there is a second //, we complain, else do nothing
      if( size(icomm,"*") > 1 ) then
        disp("Possible comment exception: "+sCode(iline))
      end
// Else, we only keep characters before the first //
    else
      sCode(iline)=strsubst(sCode(iline),"/\/\/.*/","","r");
    end
  end
//
// Remove leading and trailing blanks and tabs
  sCode=stripblanks(sCode,%t);
//
// Replace "..." continuation with a blank line and a joined line
  iCont=find(part(strrev(sCode),1:3) == "...");
  for i=iCont
    sCode(i+1)=part(sCode(i),1:($-3))+sCode(i+1);
    sCode(i)="";
  end
//
// Remove blank lines
  i=find(sCode == "");
  sCode(i)=[];
  mprintf("%d lines after removing blanks, comments, continuations\n",size(sCode,1))
//
// Find start and end lines for each function
  iStart=grep(sCode,"/^function /","r");
  iEnd=grep(sCode,"endfunction");
  if(size(iStart,2) ~= size(iEnd,2)) then
    disp("Mismatch between number of ""function"" and ""endfunction""")
    abort
  end
  nFuncs=size(iStart,2)
  mprintf("%d functions identified\n",nFuncs)
//
// Keep lines within functions in sCode[], copy the rest to sCode2[]
  sCode2=sCode;
  sCode="";
  for i=1:nFuncs
    for iline=iStart(i):iEnd(i)
      sCode($+1)=sCode2(iline);
      sCode2(iline)="";
    end
  end
  sCode(1)=[];
  i=find(sCode2 == "");
  sCode2(i)=[];
//
// Find all globals in sCode2[] and save unique names in sGlobals[]
  i=grep(sCode2,"/^global /","r");
  sGlobals=strsubst(sCode2(i),"/^global/","","r");
  sGlobals=strcat(" "+sGlobals);
  sGlobals=unique(tokens(sGlobals));
  mprintf("%d unique global variables identified\n",size(sGlobals,1))
// Remove the global lines from sCode2, and complain if anything is left
  sCode2(i)=[];
  if( sCode2 ~= [] ) then
    AbortBox("Code found not in a function or a global")
  end
//
// Extract function names
  iStart=grep(sCode,"/^function /","r");
  sFuncs=strsubst(sCode(iStart),"/^function /","","r");
  sFuncs=strsubst(sFuncs,"/\(.*/","","r");
  sFuncs=strsubst(sFuncs,"/.*=/","","r");
// Check that all sFuncs[] refer to compiled functions
  for sFunc=sFuncs'
    if(typeof(evstr(sFunc)) ~= "function") then
      AbortBox("Unloaded function found: "+sFunc)
    end
  end
//
// Create global declarations and add to start of sCode[]
  sLines="global";
  for stmp=sGlobals'
    sLine=sLines($)+" "+stmp;
    if(length(sLine) < 70) then
      sLines($)=sLine;
    else
      sLines($+1)="global "+stmp;
    end
  end
  sCode=[sLines;sCode];
//
// Add function hardwired with release date to end of sCode[]
  sCode=[sCode;
       "function str=LaueGReleaseDate()";
       "str=""" + date() + """";
       "endfunction"];
//
// Write sCode[] to Release.sce in ProgDir
  [fd,ierr]=mopen(ProgDir+"Release.sce","w")
  if(ierr ~= 0) then
    AbortBox("Unable to open "+ProgDir+"Release.sce")
  end
  mprintf("Writing %d lines of code to %sRelease.sce\n",size(sCode,1),ProgDir)
  mfprintf(fd,"%s\n",sCode)
  mclose(fd)
//
  mprintf("FINISHED copying SCILAB code\n\n")
//
endfunction


function SetupBuildFolder()
//
//
  RelDir=ProgDir+"_Release\"
  BldDir=RelDir+"WindowsBuild\"
//
// Check that wavs.lis and *.dat files are consistent
  sFiles=["wavs.lis"; findfiles(RelDir+"BuildFiles","*.dat")]
  for fnam=sFiles'
// Read the file from RelDir
    [fd,ierr]=mopen(RelDir+"BuildFiles\"+fnam,"r")
    if(ierr ~= 0) then
      mprintf("Cannot find %s in %s\n",fnam,RelDir+"BuildFiles\")
      return
    end
    buff1=mgetl(fd,-1)
    mclose(fd)
// Read the file from ProgDir
    [fd,ierr]=mopen(ProgDir+fnam,"r")
    if(ierr ~= 0) then
      mprintf("Cannot find %s in %s\n",fnam,ProgDir)
      return
    end
    buff2=mgetl(fd,-1)
    mclose(fd)
//
    if( ~and(buff1 == buff2) ) then
      mprintf("%s is different in %s and %s\n",fnam, ...
                              ProgDir,RelDir+"BuildFiles\")
      return
    end
  end
//
// Remove old WindowsBuild folder
  if( isdir(BldDir) ) then
    if( ~removedir( BldDir ) ) then
      mprintf("Unable to delete old build folder %s\n",BldDir)
      return
    end
  end
//
// Create a new WindowsBuild folder
  sleep(100)  // apparently needed by Windows
  if( ~createdir( part(BldDir,1:$-1) ) ) then
    mprintf("Unable to create new build folder %s\n\n",BldDir)
    return
  end
  mprintf("Created new build folder %s\n",BldDir)
//
// Create list of executables used by RunConsole.sce
// Read RunConsole.sce into buff[]
  [fd,ierr]=mopen(ProgDir+"\RunConsole.sce")
  buff=mgetl(fd,-1)
  mclose(fd)
// Extract names of executables that are "run"
  buff=stripblanks(buff,%t)
  i=grep(buff,"/^RunConsole.*\(/","r")
  buff=strsubst(buff(i),"/^RunConsole.*\(/","","r")
  buff=stripblanks(buff,%t)
  buff=strsubst(buff,"/^[''""]/","","r")
  buff=strsubst(buff,"/[''""].*/","","r")
  buff=unique(buff)
// Create list of executable file names (with .exe)
  sFiles=buff+".exe"
  mprintf("%d *.exe files referenced in RunConsole.sce\n",size(sFiles,1))
//
// Copy each exe file to the build folder
  for fnam=sFiles'
    iOK=copyfile(ProgDir+fnam, BldDir);
    if(iOK ~= 1) then
      mprintf("Error copying %s to build folder\n",fnam)
      return
    end
  end
  mprintf("Copied %d *.exe files from %s\n",size(sFiles,1),ProgDir)
//
// Copy the Scilab code files
  iOK1=copyfile(ProgDir+"LaueG.sce", BldDir)
  iOK2=copyfile(ProgDir+"Release.sce", BldDir)
  if( iOK1+iOK2 ~= 2 ) then
    mprintf("Error copying LaueG.sce or Release.sce\n")
    return
  end
  mprintf("Copied LaueG.sce and Release.sce from %s\n",ProgDir)
//
// Copy all files from _Release\BuildFiles\
  iOK=copyfile(RelDir+"BuildFiles\",BldDir)
  if(iOK ~= 1) then
    mprintf("Error copying files from %s\n",RelDir+"BuildFiles\")
    return
  end
  mprintf("Copied all files from %s\n",RelDir+"BuildFiles\")
//
// Remove the copied Windows subfolders (don't care if it fails)
   removedir( BldDir+"\Windows" )
//
// Copy files from _Release\BuildFiles\Windows
  iOK=copyfile(RelDir+"BuildFiles\Windows\",BldDir)
  if(iOK ~= 1) then
    mprintf("Error copying files from %s\n",RelDir+"BuildFiles\")
    return
  end
  mprintf("Copied all files from %s\n",RelDir+"BuildFiles\Windows\")
//
  mprintf("\nTo finish Windows build, COMPILE WindowsInstall.iss in %s\n",RelDir)
  mprintf("The installation file will be created in %s\n\n",RelDir+"Installs\")
//
endfunction


// ======================= Experimental Functions =======================

function UBOptim()
global LatticeInfo
//
// The UB from modals often matches only a portion of the spots.
// This routine tries matrix transforms on the UB consisting of
// simple transformations such as doubling/halving cell lengths
// and shearing axes relative to the ones in the UB.
// It has worked at least once to convert from a partial match
// to an excellent match after a few iterations.
// The transformations include a rescaling that ensures the unit cell
// volume is preserved.
//
  mprintf("Processing: 0%%")
//
// Set number of observed spots, and wav_min & d_min for calculated spots
  nobs=100
  wav_min=1.4
  d_min=3.0
//
// Create vector of possible transformation matrices
  imat=0
  mats=[]
  for Rx=[1,2,0.5]
    for Ry=[1,2,0.5]
      for Rz=[1,2,0.5]
        mat_diag=[Rx,0,0;0,Ry,0;0,0,Rz];
        imat=imat+1; mats(:,:,imat)=mat_diag;
        for iOff1=1:3
          for iOff2=1:3
            if(iOff1 ~= iOff2) then

              for Roff=[-0.5,0.5]
                imat=imat+1; mats(:,:,imat)=mat_diag;
                mats(iOff1,iOff2,imat)=Roff
              end
            end
          end
        end
      end
    end
  end
//
// Normalise matrices to determinant = 1
  nmats=size(mats,3)
  mats2=[]
  for i=1:nmats
    mats2(:,:,i)=mats(:,:,i) * det(mats(:,:,i))^(-1/3)
  end
//
// Start with default lattice, triclinic and no centering
  Linfo=LatticeInfo(1)
  LInfo.ilatt=1
  LInfo.icen=0
  LInfo.special_pairs=[]
//
// Count number of matched spots and summed merits for matches using new UBs
  nmatch=[]
  merit=[]
  for i=1:nmats
    if( modulo(i,round(nmats/10.9)) == 0 ) then
      mprintf(" %i%%",100*i/nmats)
    end
    LInfo.ub = LatticeInfo(1).ub / mats2(:,:,i)
    hkllist=RunGenHKLs(Linfo,wav_min,1e4, d_min,1e4, %f)
    match=MatchXY(DisplayData.found_spots(1:nobs,1:2),hkllist(:,5:6), 10)
    nmatch(i)=size(match,1)

    merit(i)=sum(DisplayData.found_spots(match(:,1),3))
  end
  mprintf(" 100%%\n\n")
//
// Print out information on 9 solutions with maximum merit
  [vals,idx]=gsort(merit);
  mprintf("Top 9 solutions:\n#  Merit1   Merit2         Cell Dimensions\n")
  for i0=1:9
    i=idx(i0)
    ucell=UB2UnitCell( LatticeInfo(1).ub / mats2(:,:,i) )
    mprintf("%1i%8i%7i  %6.2f%6.2f%6.2f %6.1f%6.1f%6.1f\n",i0,merit(i),nmatch(i),ucell)
  end
//
// Print out information on the untransformed matrix
  ucell=UB2UnitCell( LatticeInfo(1).ub );
  mprintf("*%8i%7i  %6.2f%6.2f%6.2f %6.1f%6.1f%6.1f\n",merit(1),nmatch(1),ucell)
  mprintf("     * denotes the original UB matrix\n")
//
// Ask which solution to load
  sVal=input("Choose solution to use (or <RET> for none):","string")
  iVal=0
  if( isnum(sVal) ) then
    iVal=strtod(sVal)
  end
//
// Exit if using the original UB
  if( (iVal < 1) | (iVal > 9) )
    mprintf("Reloading original UB\n")
    return
  end
//
// Update LatticeInfo and output results
  LatticeInfo.ub=LatticeInfo(1).ub / mats2(:,:,idx(iVal))
  LatticeInfo.cell=UB2UnitCell(LatticeInfo(1).ub)
  LatticeInfo.ilatt=1
  LatticeInfo.icen=0
  LatticeInfo.special_pairs=[]
  mprintf("Loading solution %d\n",iVal)
  mprintf("Transformation Matrix:\n")
  mprintf("    %9.5f%9.5f%9.5f\n",mats2(:,:,idx(iVal)))
  mprintf("Cell Dimensions:\n")
  mprintf("    %6.2f%6.2f%6.2f %6.1f%6.1f%6.1f\n",LatticeInfo.cell)
//
endfunction


//========================================================//
//########################################################//
//========================================================//

function jacob_merge

fd=mopen("CCT4_1_4K_fast_mrg.hkl","rt")
buff1=mgetl(fd)
mclose(fd)
data1=msscanf(-1,buff1,"%4d%4d%4d%8f%8f\n")
hkl1=data1(:,1:3)
int1=data1(:,4:5)

fd=mopen("CCT4_1_4K_slow_mrg.hkl","rt")
buff2=mgetl(fd)
mclose(fd)
data2=msscanf(-1,buff2,"%4d%4d%4d%8f%8f\n")
hkl2=data2(:,1:3)
int2=data2(:,4:5)

[nb,loc]=members(hkl1,hkl2,"rows")
i1=find(loc>0)
i2=loc(i1)
[hkl1(i1,:)==hkl2(i2,:)]

for nord=1:3
  sum1=sum( ( int1(i1,1) ./ int2(i2,2) ).^nord )
  sum2=sum( ( int2(i2,1) ./ int2(i2,2) ).^nord )
  rats(nord)=(sum2/sum1).^(1/nord)
end

if(rats(1) < 0) then
  for nord=1:3
    sum1=sum( ( int1(i1,1) ./ int1(i1,2) ).^nord )
    sum2=sum( ( int2(i2,1) ./ int1(i1,2) ).^nord )
    rats(nord)=(sum2/sum1).^(1/nord)
  end
end


if(rats(1) > 1) then
  int2=int2/rats(3)
else
  int1=int1*rats(3)
end



sum12=int1(i1,1) .* int1(i1,2).^(-2) + int2(i2,1) .* int2(i2,2).^(-2)
var12=( int1(i1,2).^(-2) + int2(i2,2).^(-2) ) .^(-1)
data12=[hkl1(i1,:), sum12 .* var12 , var12 .^ (1/2)]

data1(i1,:)=[]
data2(i2,:)=[]
data12=[data12; data1; data2]

fd=mopen("CCT4_1_4K_both_mrg.hkl","wt")
mfprintf(fd,"%4d%4d%4d%8.1f%8.1f\n",data12)
mclose(fd)


endfunction





function run_siggi_run()
global ModulateInfo DisplayMode DisplayData

CreateSingleImageWindow()
LockMenus("Searching for spots")
ClearObsSpots()
RunFindSpots(15)
DrawObsSpots(DisplayData.found_spots(:,1:2))
SingleImageOrient()
SingleImageOrient()
ModulateInfo.d_min=0.5
ModulateInfo.d_max=10
DisplayMode.modul_show=%t
DrawGenHKLs(0)
  

nobs=400
// Copy found spots to xy_obs
xy_obs=DisplayData.found_spots(1:nobs,1:2);
// Read calc. spots file into hkllist
fd=mopen("___laueg_gen_hkls.out","r");
buff=mgetl(fd);
mclose(fd);
hkllist=evstr(buff);
// Find obs to calc matches within 10 pixels
match=MatchXY(xy_obs,hkllist(:,5:6), 10);
// Find unique matches
ihist=histc(nobs,match(:,1),%f)
ifind=find(ihist(match(:,1)) == 1);
// Prune obs & calc lists to the unique matches
xy_obs=xy_obs(match(ifind,1),:);
hkllist=hkllist(match(ifind,2),:);
// Calculate (fractional) hkl for observed spots
pixwav=[xy_obs hkllist(:,7)]
hkl_obs=PixWav2HKL(pixwav)
// Split fractional hkl into integer hkl and modulation offset
hkl_main=hkllist(:,1:3)
hkl_mod=hkllist(:,$-2:$)-hkl_main
// Prune obs & calc lists to the satellites only
ifind=find(abs(hkl_mod(:,2)) > 0.001)
xy_obs=xy_obs(ifind,:)
hkl_obs=hkl_obs(ifind,:)
hkllist=hkllist(ifind,:)
hkl_main=hkl_main(ifind,:)
hkl_mod=hkl_mod(ifind,:)
// Display the pruned observed spots
DisplayData.found_spots=xy_obs
DisplayData.found_spots(:,3)=1000
DrawObsSpots(DisplayData.found_spots(:,1:2))
// Create normalised fractional obs. hkls
size_obs=sqrt(sum(hkl_obs.^2,2))
norm_obs=hkl_obs ./ (size_obs*[1,1,1])


// Find rmin a minimum of rmult within 0.001
rmults=[0.90:0.001:1.10]
vals=[];
for rmult=rmults
  hkl_try=hkl_main + rmult*hkl_mod;
  size_try=sqrt(sum(hkl_try.^2,2));
  norm_try=hkl_try ./ (size_try*[1,1,1]);
  diffs=norm_obs-norm_try;
  vals= [vals , sum( diffs(:,2).^2 ) ];
end
[val,idx]=min(vals); rmin=rmults(idx);
// Create values +/-0.0015 around rmin
rmults=[-.0015:0.0001:0.0015]+rmin
vals=[];
for rmult=rmults
  hkl_try=hkl_main + rmult*hkl_mod;
  size_try=sqrt(sum(hkl_try.^2,2));
  norm_try=hkl_try ./ (size_try*[1,1,1]);
  diffs=norm_obs-norm_try;
  vals= [vals , sum( diffs(:,2).^2 ) ];
end
[val,idx]=min(vals); rmin1=rmults(idx);

// Prune spots with diffs > 0.005
rmult=rmin1
hkl_try=hkl_main + rmult*hkl_mod;
size_try=sqrt(sum(hkl_try.^2,2));
norm_try=hkl_try ./ (size_try*[1,1,1]);
diffs=norm_obs-norm_try;
sdiffs=sum(abs(diffs(:,:)),2)
f=figure();
plot(sdiffs);
ifind=find(sdiffs < 0.005);

// Create values +/-0.0015 around rmin
rmults=[-.0015:0.0001:0.0015]+rmin
vals=[];
for rmult=rmults
  hkl_try=hkl_main + rmult*hkl_mod;
  size_try=sqrt(sum(hkl_try.^2,2));
  norm_try=hkl_try ./ (size_try*[1,1,1]);
  diffs=norm_obs-norm_try;
  vals= [vals , sum( diffs(ifind,2).^2 ) ];
end
[val,idx]=min(vals); rmin2=rmults(idx);

// Prune spots with diffs > 0.003
rmult=rmin2
  hkl_try=hkl_main + rmult*hkl_mod;
  size_try=sqrt(sum(hkl_try.^2,2));
  norm_try=hkl_try ./ (size_try*[1,1,1]);
  diffs=norm_obs-norm_try;
sdiffs=sum(abs(diffs(:,:)),2)
ifind=find(sdiffs < 0.003);

// Create values +/-0.0015 around rmin
rmults=[-.0015:0.0001:0.0015]+rmin
vals=[];
for rmult=rmults
  hkl_try=hkl_main + rmult*hkl_mod;
  size_try=sqrt(sum(hkl_try.^2,2));
  norm_try=hkl_try ./ (size_try*[1,1,1]);
  diffs=norm_obs-norm_try;
  vals= [vals , sum( diffs(ifind,2).^2 ) ];
end
[val,idx]=min(vals); rmin3=rmults(idx);

mprintf("nu = %.5f %.5f %.5f,  nsatt=%d\n", ...
   [rmin1 rmin2 rmin3]*max(hkl_mod),size(ifind,2))


endfunction




function junko_twins()

hkl_mod=ModulateInfo.idxs*ModulateInfo.vecs
hkl_200=hkl_mod
hkl_200(:,2)=hkl_mod(:,2)-2
hkl_400=hkl_mod
hkl_400(:,2)=hkl_mod(:,2)-4
hkl_600=hkl_mod
hkl_600(:,2)=hkl_mod(:,2)-6
hkl_200=hkl_200 ./ (-hkl_200(:,2) * [1,1,1])
hkl_400=hkl_400 ./ (-hkl_400(:,2) * [1,1,1])
hkl_600=hkl_600 ./ (-hkl_600(:,2) * [1,1,1])




hkllist=RunGenHKLs(LatticeInfo(1),1,99,3,99,%f)
ifind1=find(abs(hkllist(:,1)) < 3)
ifind2=find(abs(hkllist(:,2)) < 3)
ifind3=find(abs(hkllist(:,3)) < 3)
ifind=intersect(ifind1,ifind2)
ifind=intersect(ifind,ifind3)
hkllist=hkllist(ifind,:)

hvec=HKL2Hvec(hkllist(:,1:3))
sizes=sqrt(sum(hvec(:,:).^2,2))
hvec(:,1)=hvec(:,1) ./ sizes
hvec(:,2)=hvec(:,2) ./ sizes
hvec(:,3)=hvec(:,3) ./ sizes

a=hvec * hvec'
ang=round(acosd(a*(1-1e-6)))
int32(ang)

[i10a,i10b]=find(ang == 118)

ang

aaa=[]
for i=1:size(i10a,2)/2
  aaa(i)=a(i10a(i),i10b(i))
end
[i10a',i10b',aaa]



[i10a,i10b]=find(abs(a) < 0.1)
gsort(i10a)



[i10a;i10b]


ifind=find(i10a ~= i10b)
[i10a(ifind);i10b(ifind)]






hkl=[]
for i=1:size(DisplayData.found_spots,1)
  hkl(i,:)=CalcHKLGuess(DisplayData.found_spots(i,1:2),20)
end
hvec=HKL2Hvec(hkl)


hvec=Pix2Hvec(DisplayData.found_spots(:,1:2))
sizes=sqrt(sum(hvec(:,:).^2,2))
hvec(:,1)=hvec(:,1) ./ sizes
hvec(:,2)=hvec(:,2) ./ sizes
hvec(:,3)=hvec(:,3) ./ sizes


a=hvec * hvec'

round(acosd(a*(1-1e-6)))

[i10a,i10b]=find(abs(a) < 0.03)



///////////////////////////////////////


ub1=LatticeInfo(2).ub;
ub2=LatticeInfo(3).ub;
ub2=ub2 /(det(ub2)/det(ub1))^(1/3);
R0=ub1\ub2;
R=round(R0*12)/12;
mprintf("Transformation Matrix from Twin 1 to Twin 2:\n")
mprintf("%8.4f%8.4f%8.4f\n",inv(R))


h1=[1;2;-2]
h2=inv(R)*h1
v1=ub1*h1; e1=v1/norm(v1);
v2=ub2*h2; e2=v2/norm(v2);
RR=eye(3,3)+(e1-e2)*e2';
LatticeInfo(3).ub=RR*ub2

endfunction
///////////////////////////////////////////////////////////////////////////


function junko4


  rot_ub=CalcYRot(-165) * Linfo.ub
  rot_ub*pars.hkls'


// Calculate new UB using orient_spots




[LatticeInfo,OrientInfo,text]=RunOrientSpots(Linfo,Oinfo,5,pars)
LatticeInfo.cell=UB2UnitCell(LatticeInfo.ub)


hvecs=Linfo.ub * pars.hkls'



hvecs=HKL2Hvec(pars.hkls)'
pvecs=Pix2Hvec(pars.xys)'
scales=mean(hvecs ./ pvecs,"r")
spvecs=([1 1 1]' * scales) .* pvecs


Hmat=pars.hkls' * pars.hkls

ub2=inv(Hmat)*(spvecs*hvecs')


ub2=(spvecs*hvecs') * inv(Hmat)


(ub2*pars.hkls') ./ pvecs


endfunction




function junko3

ub0=LatticeInfo(1).ub

hmax=3
kmax=3
lmax=3

hkls=[]
for h=0:hmax
  for k=-kmax:kmax
    for l=-lmax:lmax
//
      igcd=gcd([h,k,l]);
      if(igcd == 1 ) then
        hkls=[hkls;[h,k,l]];
      end
//
    end
  end
end
hkls=[hkls,hkls*ub0']
dots=hkls(:,4:6) * hkls(:,4:6)'
d_dots=sqrt(diag(dots))
dots=dots ./ (d_dots*d_dots')

[i1,i2]=find(abs(dots) < cosd(89.5))
for i=1:size(i1,2)
  if(i2(i) > i1(i)) then
    a=1/d_dots(i1(i));
    b=1/d_dots(i2(i));
    alpha=acosd(dots(i1(i),i2(i)));
    mprintf("%6.2f%6.2f%6.1f%3i%3i%3i   %3i%3i%3i\n",a,b,alpha,hkls(i1(i),1:3),hkls(i2(i),1:3))
  end
end


u_dots=unique(-matrix(dots-triu(dots),-1,1))

ifind=find(dots > cosd(1.6))



acosd(u_dots($:-1:$-5))'









endfunction



function junko2




a=[1 2; 3 5; 4 2]
b=[7 9; 6 1; 8 1]




hkl1=CalcHKLGuess(DisplayData.found_spots(1,1:2),10)
hkl2=CalcHKLGuess(DisplayData.found_spots(2,1:2),10)
hkl3=CalcHKLGuess(DisplayData.found_spots(3,1:2),10)
hkl4=CalcHKLGuess(DisplayData.found_spots(4,1:2),10)

vec1=Pix2Hvec(DisplayData.found_spots(1,1:2))
vec2=Pix2Hvec(DisplayData.found_spots(2,1:2))
vec3=Pix2Hvec(DisplayData.found_spots(3,1:2))
vec4=Pix2Hvec(DisplayData.found_spots(4,1:2)+[0,0])
obs_vec=[vec1' vec2' vec3' vec4']

hvec1=HKL2Hvec(hkl1); hvec1=hvec1/norm(hvec1)
hvec2=HKL2Hvec(hkl2); hvec2=hvec2/norm(hvec2)
hvec3=HKL2Hvec(hkl3); hvec3=hvec3/norm(hvec3)
hvec4=HKL2Hvec(hkl4); hvec4=hvec4/norm(hvec4)
calc_vec=[hvec1' hvec2' hvec3' hvec4']

M=obs_vec/calc_vec

hvec1=CalcYRot(-phi) * ub * hkl1; hvec1=hvec1/norm(hvec1)




endfunction

function junko()
// Ask for new image to draw, do spot search and draw calculated
// spots according to current display settings
CreateSingleImageWindow()
RunFindSpots(15)
DrawObsSpots(DisplayData.found_spots(:,1:2))
DrawGenHKLs(0)
SpotsObs=DisplayData.found_spots;
SpotsCalc=MainMarks(2).user_data;


// Remove any observed or calculated spots not fitting to "radius" pixels
radius=5
match=MatchXY(SpotsCalc,SpotsObs(:,1:2), radius);
[val,idx]=gsort(-match(:,2));
isize=size(match,1)/2
match=match(idx(1:isize),:);
SpotsCalc=SpotsCalc(match(:,1),:);
SpotsCalc(:,3:5)=0;
for i=1:isize
  SpotsCalc(i,3:5)=CalcHKLGuess(SpotsCalc(i,1),SpotsCalc(i,2),radius);
end
SpotsObs=SpotsObs(match(:,2),:);
DisplayData.found_spots=SpotsObs;
DrawObsSpots(DisplayData.found_spots(:,1:2))
DrawCalcSpots(SpotsCalc(:,1:2))

////////////////

// Ask for new image to draw, do spot search and draw calculated spots

// from previous image
CreateSingleImageWindow()
RunFindSpots(15)
DrawObsSpots(DisplayData.found_spots(:,1:2))
DrawCalcSpots(SpotsCalc(:,1:2))

// Prune observed spots to be within "radius" pixels of a calculated,
// and of merit greater than the "cutoff" time the maximum of
// observed spots associated with the same calculated spots.
radius=70
cutoff=0.1
SpotsObs=DisplayData.found_spots;
match=MatchXY(SpotsCalc(:,1:2),SpotsObs(:,1:2), radius);
match(:,3)=SpotsObs(match(:,2),3);
match2=[];
for icalc=unique(match(:,1))'
  iobs=find(match(:,1) == icalc);
  ifind=find( match(iobs,3) > max(match(iobs,3))*cutoff );
  match2=[match2 ; match(iobs(ifind),:)];
end
SpotsObs=SpotsObs(match2(:,2),:);
DisplayData.found_spots=SpotsObs;
DrawObsSpots(DisplayData.found_spots(:,1:2))

// For each calculated spot: remove all observed spots within "radius"
// pixels if the spots are not within "rad_match" pixels of being
// symmetrical about the mean x,y point. Calculated spots x,y are
// changed to the mean observed x,y point. Removing a group of observed
// spots also removes the associated calculated spot.
rad_match=5
SpotsObs=DisplayData.found_spots;
match=MatchXY(SpotsCalc(:,1:2),SpotsObs(:,1:2), radius);
obs2=[];
calc2=[];

for icalc=unique(match(:,1))'
  iobs=find(match(:,1) == icalc);
  obs=SpotsObs(match(iobs,2),:);
  xyave=mean(obs(:,1:2),"r");
  reflect=[2*xyave(1)-obs(:,1) , 2*xyave(2)-obs(:,2) ];
  match3=MatchXY(obs(:,1:2),reflect,rad_match);
  if( size(match3,1) == size(obs,1) ) then
    obs2=[ obs2 ; obs ];
    calc2=[ calc2 ; xyave , SpotsCalc(icalc,3:5) ];
  end
end
SpotsObs=obs2;
DisplayData.found_spots=SpotsObs;
DrawObsSpots(DisplayData.found_spots(:,1:2))
SpotsCalc=calc2;
DrawCalcSpots(SpotsCalc(:,1:2))



// Update current UB with triclinic cell from SpotsCalc which
// contains the x,y averaged over the associated observed spots
//
// For pixwav, x,y from observed and wav from hkl
pixwav=HKL2PixWav(SpotsCalc(:,3:5));
pixwav(:,1:2)=SpotsCalc(:,1:2);
Hvec_obs=PixWav2Hvec(pixwav);
// Calculate new UB
rot_ub_tran=SpotsCalc(:,3:5) \ Hvec_obs;
ub=CalcYRot(ImageInfo.phi) * rot_ub_tran';
ucell=UB2UnitCell(ub)
// Modify cell and ub for orthorhombic
ucell(4:6)=90
b=Cell2Bmatrix(ucell)
u=ub/b
u=(u+1/u')/2
u=(u+1/u')/2
u=(u+1/u')/2
ub=u*b
// Load UB & cell, and set centering to P
LatticeInfo(1).ub=ub;
LatticeInfo(1).cell=ucell;
LatticeInfo(1).icen=0;


match=MatchXY(SpotsCalc(:,1:2),SpotsObs(:,1:2), radius);
HvecObs=SpotsObs(:,1:2);
for icalc=unique(match(:,1))'
  pixwav=HKL2PixWav( SpotsCalc(icalc,3:5) ); wav=pixwav(3);
  iobs=find(match(:,1) == icalc);
  obs=SpotsObs(match(iobs,2),:); obs(:,3)=wav;
  Hvec=PixWav2Hvec(obs);
  HvecObs(match(iobs,2),3)=SpotsCalc(icalc,3);
  HvecObs(match(iobs,2),4)=SpotsCalc(icalc,4);
  HvecObs(match(iobs,2),5)=SpotsCalc(icalc,5);
  HvecObs(match(iobs,2),6)=Hvec(:,1)-mean(Hvec(:,1));
  HvecObs(match(iobs,2),7)=Hvec(:,2)-mean(Hvec(:,2));
  HvecObs(match(iobs,2),8)=Hvec(:,3)-mean(Hvec(:,3));
end
HvecObs



ub0=LatticeInfo(1).ub;
for i=1:size(HvecObs,1)
  hkl=HvecObs(i,3:5);
  Hvec=(ub0*hkl')';
  HvecObs(i,9:11)=Hvec;
end


for icalc=1:3
  imatch=MatchXY(SpotsCalc(icalc,1:2),SpotsObs(:,1:2), radius);
  nspots(icalc)=size(imatch,1)
end

ival=0;
vals=[];
ubs=[];

for i1=1:nspots(1)
  for i2=1:nspots(2)
    for i3=1:nspots(3)
      idx=[i1,i2,i3];
      ival=ival+1;
      [ndiff,ub]=test_idx(idx);
      vals(:,ival)=[i1,i2,i3,ndiff]';
      ubs(1:3,1:3,ival)=ub;
    end
  end
end
vals

[idum,idx2]=gsort(-vals(4,:));
vals(:,idx2(1:10))
ubs(:,:,idx2(1:10))


LatticeInfo(1).ub=ubs(:,:,idx2(5))


endfunction


function [n2,ub]=test_idx(idx)
//
Hvec=[]
hkl=[]
for icalc=1:3
  imatch=MatchXY(SpotsCalc(icalc,1:2),HvecObs(:,1:2), radius)
  iobs=imatch(idx(icalc),2)
// Calculate Hvec for calculated HKL + observed twin offsets in HKL
  hkl(icalc,:)=HvecObs(iobs,3:5)
  Hvec(icalc,:)=HvecObs(iobs,9:11)+HvecObs(iobs,6:8)
end
ub=Hvec' / hkl'
//
for icalc=4:6
  imatch=MatchXY(SpotsCalc(icalc,1:2),HvecObs(:,1:2), radius)
  if(imatch == []) then continue; end
  iobs=imatch(:,2)
// Calculate Hvec for calculated HKL + observed twin offsets in HKL
  hkl($+1,:)=HvecObs(iobs(1),3:5)
  nbest=1e6
  for i=iobs'
    htry=HvecObs(i,9:11)+HvecObs(i,6:8)
    n=norm( inv(ub) * htry(1,:)' - hkl(icalc,:)' )
    if(n < nbest) then
      nbest=n
      ibest=i
    end
  end
  Hvec($+1,:)=HvecObs(i,9:11)+HvecObs(i,6:8)
end
ub=Hvec' / hkl'
//
//
for icalc=13  //7:size(SpotsCalc,1)
  imatch=MatchXY(SpotsCalc(icalc,1:2),HvecObs(:,1:2), radius)
  if(imatch == []) then continue; end
  iobs=imatch(:,2)
// Calculate Hvec for calculated HKL + observed twin offsets in HKL
  hkl($+1,:)=HvecObs(iobs(1),3:5)
  nbest=1e6
  for i=iobs'
    htry=HvecObs(i,9:11)+HvecObs(i,6:8)
    n=norm( inv(ub) * htry(1,:)' - hkl(icalc,:)' )
    if(n < nbest) then
      nbest=n
      ibest=i
    end
  end
  Hvec($+1,:)=HvecObs(i,9:11)+HvecObs(i,6:8)
end
ub=Hvec' / hkl'
n2=norm( inv(ub) * Hvec' - hkl' )
//
endfunction




//vec=calc_obs(match(:,1),1:2)-SpotsObs(match(:,2),1:2);
//match(:,3)=sqrt(sum(vec.^2,"c"));
//match(:,4)=SpotsObs(match(:,2),3);






//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_212121()

ang=0.58
hkl_piv=[0,1,-1]

para=ub*hkl_piv'
para=para/norm(para)
perp1=[0;para(3);-para(2)]
perp2=[para(3);0;-para(1)]
if(norm(perp1) < norm(perp2)) then
  perp1=perp2
end
perp2=cross(para,perp1)
para=para/norm(para)
perp1=perp1/norm(perp1)
perp2=perp2/norm(perp2)
Urot=[para,perp1,perp2]

rot=Urot * [1,0,0;0,cosd(ang),sind(ang);0,-sind(ang),cosd(ang)] * Urot'


vec12=Pix2Hvec(hkllist(5:6,4:5))'
Hvec2HKL( cross(vec12(:,1),vec12(:,2))' )*180

endfunction
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// Separate spots list into unique hklXY1[] and pairs of spots hklXY1XY2[]
// Calculate rotation matrix where hvec2 ~ hvec1 * rot
function [hklXY1,hklXY1XY2,rot]=IndexSplitSpots(spots,maxhkl,minwav,maxdev)
//
// Estimate hkl & wavs for observed spots
  hkllist=CalcObsHKLs(spots,maxhkl,minwav);
//
// Remove spots more than 3*maxdev pixels from calc. XY
  pix=HKL2Pix(hkllist(:,1:3));
  dev=sqrt(sum( (pix - hkllist(:,4:5)).^2, "c"));
  i=find(dev > 3*maxdev);
  hkllist(i,:)=[];
//
// Create HKL & XY lists for unique and split pairs of spots
  isingle=[];
  isplit=[];
  for i1=1:size(hkllist,1)
    hkl1=hkllist(i1,1:3);
    i2=find( (hkl1(1) == hkllist(:,1)) & (hkl1(2) == hkllist(:,2)) & ...
                                            (hkl1(3) == hkllist(:,3)) );
    if(i2(1) == i1) then
      if(size(i2,2) == 1) then
        isingle($+1)=i1;
      else
        isplit($+1,1:2)=i2(1:2);   // only use first two matches
      end
    end
  end
  hklXY1=hkllist(isingle,1:5);
  hklXY1XY2=[hkllist(isplit(:,1),1:5),hkllist(isplit(:,2),4:5)];
//
// Swap 1 & 2 in HKL2PixWav[] if spot #2 closer to calc. XY
  pixwav=HKL2PixWav(hklXY1XY2(:,1:3));
  tmp=hklXY1XY2(:,4:7) - pixwav(:,[1,2,1,2]);
  dsq1=tmp(:,1).^2 + tmp(:,2).^2;
  dsq2=tmp(:,3).^2 + tmp(:,4).^2;
  i=find(dsq2 < dsq1);
  hklXY1XY2(i,4:7)=hklXY1XY2(i,[6:7,4:5]);
//
//
// Progessively prune hklXY1XY2[] using decreasing mismatch cutoffs:
  phi_rot=CalcYRot(-ImageInfo.phi);
  for cutoff=[3.0,2.0,1.5,1.2,1.0]*maxdev
// Create unit scatt. vectors for XY of spot pairs
    hvec1=Pix2Hvec( hklXY1XY2(:,4:5) );
    hvec2=Pix2Hvec( hklXY1XY2(:,6:7) );
// Create UB rotation matrix to match the #2 spots
    rot=phi_rot' * (hvec1\hvec2)' * phi_rot;
// Massage rot[] to a pure rotation matrix
/////////////////////    rot=(rot+inv(rot'))/2;
/////////////////////    rot=(rot+inv(rot'))/2;
// Prune spots so XY of spots #2 calculated using rot[]
// are within cutoff pixels of their observed XY
    pix2=HKL2Pix_UB(hklXY1XY2(:,1:3),rot*LatticeInfo.ub)
    tmp=pix2 - hklXY1XY2(:,6:7);
    tmp=sqrt( tmp(:,1).^2 + tmp(:,2).^2 );
    i=find(tmp < cutoff);
    hklXY1XY2=hklXY1XY2(i,:);
  end
//
endfunction


function [xy]=HKL2Pix_UB(hkls,ub)
global LatticeInfo
  ub_save=LatticeInfo.ub;
//
// Calc XY from hlk using given UB
  LatticeInfo.ub=ub;
  xy=HKL2Pix(hkls);
//
  LatticeInfo.ub=ub_save;
endfunction


// Calculate hkls and wavelength for spot XY for current orientation
// Spot with large hkls or small wavelengths are removed
function hkllist=CalcObsHKLs(spots,maxhkl,minwav)
//
// Convert x,y of spots into hkls using wav. = 1.0
  pixwav=spots;
  pixwav(:,3)=1.0;
  hkl=PixWav2HKL(pixwav);
//
// Find smallest indice of hkl not "maxhkl" smaller than the largest
  hkl_sort=gsort(abs(hkl),"c");
  hkl_max=hkl_sort(:,1);
  hkl_mid=hkl_sort(:,2);
  hkl_min=hkl_sort(:,3);
  i=find(hkl_min*maxhkl < hkl_max); hkl_min(i)=hkl_mid(i);
  i=find(hkl_min*maxhkl < hkl_max); hkl_min(i)=hkl_max(i);
//
// Scale hkls to make smallest indice +/- 1
  hkl=hkl ./ (hkl_min*[1,1,1]);
//
// Scale hkl to get best integer values
  hklr=[]
  for r=1:(maxhkl+1)
    hklr(:,:,r)=hkl * r;
  end
  [v,k]=min( squeeze( sum( (hklr-round(hklr)).^2 ,2) ) ,"c");
  hkl=round( hkl .* (k*[1,1,1]) );
//
// Recalculate wavs for integer hkls
  pixwav(:,3)=HKL2Wav(hkl);
//
// Package hkls, xy_obs, wav_calc into one array
  hkllist=[hkl,pixwav]
//
// Remove hkls larger than maxhkl
  i=find(max(abs(hkllist(:,1:3)),"c") > maxhkl);
  hkllist(i,:)=[];
//
// Remove wavelengths less than minwav
  i=find(hkllist(:,6) < minwav);
  hkllist(i,:)=[];
//
endfunction


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

function try_try_try3()

////////////////////

nXY=20;
spots=[];
for tmp = DisplayData.found_spots'
  x0=tmp(1); y0=tmp(2);
// Copy local pixels to Z0[]
  dX=[-nXY:nXY]; dY=[-nXY:nXY];
  Z0=double(RawImage(x0+dX,y0+dY));
// Calculate autocorrelation to Z[]
  Z1=flipdim(flipdim(Z0,1),2);
  Z1=conv2(Z0,Z1,"same");
// Cutoff peak below 10% of maximum
  Z1=max(0, Z1-0.1*max(Z1) );
// Fit ellipse to autocorrelation
  Z1=Z1/sum(Z1);
  x1=ones(dX')*dX; y1=(ones(dY')*dY)'; x2=x1.^2; y2=y1.^2; x1y1=x1 .*y1;
  M20=sum(Z1 .*x2); M02=sum(Z1 .*y2); M11=sum(Z1 .*x1y1);
  den=M20*M02-M11^2; e=M20/den/4; f=M02/den/4; h=-M11/den/2;

// Calculate "radius" of ellipse minor axis
  r_min=1.0/sqrt( e+f + sqrt( (e-f)^2 + h^2 ) );
// Create 71% sized ellipse, with sharp boundary and normalised sum
  Z2=max(0, 1 - 2*(e*x2 + f*y2 + h*x1y1) );
  Z2=Z2/sum(Z2);
// Convolute ellipse with data
  Z3=conv2(Z0,Z2,"same");
// Get pixel x,y of convolution maximum
  [v,k]=max(Z3);
// Save pixel position, ellipse centre and shape
  spots($+1,:)=[x0,y0,dX(k(1)),dY(k(2)),e,f,h,r_min];
//
end

// Remove spots more than 3 pixels from calc position
i=find(sqrt(sum(spots(:,3:4).^2,"c")) < 3)
spots=spots(i,:);

// Convert spots[] from x,y to pixel distance & minor-radius
dist=sqrt( (spots(:,1)-2000).^2 + (spots(:,2)-1000).^2 );
spots=[dist,spots(:,$)];

// Sort spots[] into increasing distance
[v,k]=gsort(-spots(:,1));
spots=spots(k,:);

// Remove any ridiculous values
i=find( (spots(:,2) > 1) & (spots(:,2) < 20) );
spots=spots(i,:);

// Linear fit versus distance, then refit to those below the line
[a,b,dev]=reglin(spots(:,1)',spots(:,2)')
i=find(spots(:,2) < (a*spots(:,1) + b));
[a,b,dev]=reglin(spots(i,1)',spots(i,2)')
figure(); plot(spots(i,1),spots(i,2))


endfunction


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


function [pairs,fails,rot1,rot2]=FitAllSplitXYs_TEST(spots,rot1,rot2,nXY)
//
// Iterate optimizing rotation matrices
  for iter=1:3
    mprintf("Iteration%2d    ",iter);
//
// Create list of fitted pair positions in pairs[]
    pairs=[];
    fails=[];
    for ispot=1:size(spots,1)
      hkl=spots(ispot,1:3)
// Calc scatt. vectors by rotating nominal vector
      hvec0=HKL2Hvec(hkl);
      hvec1=hvec0 * rot1;
      hvec2=hvec0 * rot2;
// Calc pixel positions from vectors
      xy1=Hvec2Pix(hvec1');
      xy2=Hvec2Pix(hvec2');
// Try to find local maxima for split pair
      [xy1,xy2]=FitSplitXY(xy1,xy2,nXY);
// If fit failed, save hkl and initial xy1 in fails[]
      if( xy1 == [] ) then
        fails($+1,:)=spots(ispot,1:5)
// If successful, save hkl and positions of maxima in pairs[]
      else
        pairs($+1,:)=[hkl,xy1,xy2];
      end
    end
// Fix v5 empty array bug
    if( and(pairs(1,:) == 0) ) then
      pairs(1,:)=[];
      fails(1,:)=[];
    end
//
// Recalculate rotation matrices from pairs[]
    mprintf("Spot pairs fitted:%3d    ",size(pairs,1))
// Calculate nominal scatt. vectors
    hvec0=HKL2Hvec(pairs(:,1:3));
// Create scatt. vectors from spot positions and dummy wavelength = 1
    pixwav1=pairs(:,4:5); pixwav1(:,3)=1;
    pixwav2=pairs(:,6:7); pixwav2(:,3)=1;
    hvec1=PixWav2Hvec(pixwav1);
    hvec2=PixWav2Hvec(pixwav2);
// Normalise vector lengths (i.e. ignore wavelengths)
    hvec0=hvec0 ./ ( sqrt(sum(hvec0.^2,"c")) * [1,1,1] );

    hvec1=hvec1 ./ ( sqrt(sum(hvec1.^2,"c")) * [1,1,1] );
    hvec2=hvec2 ./ ( sqrt(sum(hvec2.^2,"c")) * [1,1,1] );
// Create rotation matrices to relate directions of vectors
    rot1=hvec0 \ hvec1;
    rot2=hvec0 \ hvec2;
// Scale matrices to det=1
    rot1=rot1 * det(rot1)^(-1.0/3.0)
    rot2=rot2 * det(rot2)^(-1.0/3.0)
//
// Loop back for next iteration
  end
//
endfunction


function [xy1_ret,xy2_ret]=FitSplitXY(xy1,xy2,nXY)
//
// Setup null return values to signal fitting failure
  xy1_ret=[]
  xy2_ret=[]
//
// Setup for local data map of size 2*nXY+1
  dXY=[-nXY:nXY]
  [dY,dX]=ndgrid(dXY,dXY)
//
// Create local data map around mean position
  xy0=round( (xy1+xy2)/2 );
  Z0=double(RawImage(xy0(1)+dXY,xy0(2)+dXY));
// Convolute map with circle of radius 4 to 8
  dist=norm(xy0 - [2000,1000]);
  rad=round(3.5+2.1*dist/2000.0)
  [dYc,dXc]=ndgrid([1-rad+1:rad-1],[1-rad:rad-1])
  Zc=max(0, 1 - (dXc.^2 + dYc.^2)/rad^2 );
  Z0=conv2(Z0,Zc,"same");
//
// Round the pair of x,y displacements from map centre
  dxy1=round(xy1)-xy0;
  dxy2=round(xy2)-xy0;
//
// Calc vectors perpendicular and parallel to line from dxy1 to dxy2
  para=dxy2-dxy1; para=para/norm(para);
  perp=[para(2),-para(1)];
//
// Search for maximum -4:4 pixels perpendicular from dxy1
  for i=-4:4
    ixy=dxy1+perp*i+nXY+1;
    Zperp(i+5)=Z0(ixy(1),ixy(2));
  end
  [v,k]=max(Zperp);
  dxy1=dxy1 +perp*(k-5);
// If maximum at either end, return as failed
  if( (k == 1) | (k == 9) ) then return; end
//
// Search for maximum -4:4 pixels perpendicular from dxy2
  for i=-4:4
    ixy=dxy2+perp*i+nXY+1;
    Zperp(i+5)=Z0(ixy(1),ixy(2));
  end
  [v,k]=max(Zperp);
  dxy2=dxy2 + perp*(k-5);
// If maximum at either end, return as failed
  if( (k == 1) | (k == 9) ) then return; end
//
// Search for maximum -4:4 pixels parallel from dxy1
  for i=-4:4
    ixy=dxy1+para*i+nXY+1;
    Zpara(i+5)=Z0(ixy(1),ixy(2));
  end
  [v,k]=max(Zpara);
  dxy1=dxy1 + para*(k-5);
// If maximum at either end, return as failed
  if( (k == 1) | (k == 9) ) then return; end
//
// Search for maximum -4:4 pixels parallel from dxy2
  for i=-4:4
    ixy=dxy2+para*i+nXY+1;
    Zpara(i+5)=Z0(ixy(1),ixy(2));
  end
  [v,k]=max(Zpara);
  dxy2=dxy2 + para*(k-5);
// If maximum at either end, return as failed
  if( (k == 1) | (k == 9) ) then return; end
//
// Return with fitted positions
  xy1_ret=xy0+round(dxy1)
  xy2_ret=xy0+round(dxy2)
//
endfunction


//////////////////////////


function hkldata=read_int_file(fname)
//
// Open and read file as text, ignoring 3 line header
  fd=mopen(fname,"r");
  buff=mgetl(fd);
  mclose(fd)
  buff(1:3)=[];
//
  hkldata=msscanf(-1,buff,"%4d%4d%4d%7f%3d%8f%8f%8f%7f%6f");
// Only keep  esd > 0
  i=find(hkldata(:,10) > 0);
  hkldata=hkldata(i,:);
// Only keep intensities more than 10 esd
  i=find(hkldata(:,9) > 10*hkldata(:,10));
  hkldata=hkldata(i,:);
// Reduce data to hkl,wav,I,esd
  hkldata=hkldata(:,[1:4,9:10]);
//
endfunction


function try_try_try()
//
hkldata1=read_int_file("120K#1_fast2_1.int");
hkldata2=read_int_file("120K#1_fast2_2.int");
hkldata3=read_int_file("120K#1_fast2_3.int");
hkldata4=read_int_file("120K#1_fast2_4.int");
hkldata5=read_int_file("120K#1_fast2_5.int");
hkldata6=read_int_file("120K#1_fast2_6.int");
hkldata7=read_int_file("120K#1_fast2_7.int");
hkldata8=read_int_file("120K#1_fast2_8.int");
//
// Find hkls common to all data
hkl1=hkldata1(:,1:3) * [1e8,1e4,1]';
hkl2=hkldata2(:,1:3) * [1e8,1e4,1]';
hkl3=hkldata3(:,1:3) * [1e8,1e4,1]';
hkl4=hkldata4(:,1:3) * [1e8,1e4,1]';
hkl5=hkldata5(:,1:3) * [1e8,1e4,1]';
hkl6=hkldata6(:,1:3) * [1e8,1e4,1]';
hkl7=hkldata7(:,1:3) * [1e8,1e4,1]';
hkl8=hkldata8(:,1:3) * [1e8,1e4,1]';
tab=tabul([hkl1;hkl2;hkl3;hkl4;hkl5;hkl6;hkl7;hkl8]);
i=find( tab(:,2) == max(tab(:,2)) );
hkl0=tab(i,1);
//
// Prune data to only the common hkls
[n,k]=members(hkl0,hkl1); hkldata1=hkldata1(k,:);
[n,k]=members(hkl0,hkl2); hkldata2=hkldata2(k,:);
[n,k]=members(hkl0,hkl3); hkldata3=hkldata3(k,:);
[n,k]=members(hkl0,hkl4); hkldata4=hkldata4(k,:);
[n,k]=members(hkl0,hkl5); hkldata5=hkldata5(k,:);
[n,k]=members(hkl0,hkl6); hkldata6=hkldata6(k,:);
[n,k]=members(hkl0,hkl7); hkldata7=hkldata7(k,:);
[n,k]=members(hkl0,hkl8); hkldata8=hkldata8(k,:);
//
// Make list of hkls, wavs and ints
hkldata0=[hkldata1,hkldata2,hkldata3,hkldata4,hkldata5,hkldata6,hkldata7,hkldata8];
hkldata0=hkldata0(:,[1:5,10:11,16:17,22:23,28:29,34:35,40:41,46:47]);
hkldata0=hkldata0(:,[1:3,4:2:19,5:2:19]);
//
// Sort wavs and ints to make wavs increasing
for i=1:size(hkldata0,1)
  if(hkldata0(i,4) > hkldata0(i,11)) then
    hkldata0(i,4:11)=hkldata0(i,11:-1:4);
    hkldata0(i,12:19)=hkldata0(i,19:-1:12);
  end
end
//
// Remove if first int is within 20% of max int
i=find(hkldata0(:,12) < 0.8*max(hkldata0(:,12:19),"c"));
hkldata0=hkldata0(i,:);//
// Remove if last int is within 20% of max int
i=find(hkldata0(:,19) < 0.8*max(hkldata0(:,12:19),"c"));
hkldata0=hkldata0(i,:);
//
// Get wavelength of max int
wav_max=[];
for i=1:size(hkldata0,1)
  [v,k]=max(hkldata0(i,12:19));
  wav_max(i)=hkldata0(i,k+3);
end


clf();
for i=1:size(hkldata0,1)
  plot(hkldata0(i,4:11)/median(wav_max),hkldata0(i,12:19))
end


i=find(wav_max > 1.5*median(wav_max));
[hkldata0(i,1:3),wav_max(i),wav_max(i)/median(wav_max)]
wav_scale=ones(hkldata0(:,1))
wav_scale(i)=[2,3,2,2]'


clf();
for i=1:size(hkldata0,1)
  plot(hkldata0(i,4:11)/median(wav_max) ./ wav_scale(i),hkldata0(i,12:19))
end



i=find(wav_max < 1.5*median(wav_max));
[hkldata0(i,1:3),wav_max(i)/median(wav_max)]


imin=find(wav_max/median(wav_max) ./ wav_scale < 0.7)


//
// Normalise ints to maximum per hkl, and scale wavelength for maximum to be 1.3 Ang
int_max=[]; wav_max=[];
for i=1:size(hkldata0,1)
  [v,k]=max(hkldata0(i,8:11));
  int_max(i)=hkldata0(i,k+7);
  wav_max(i)=hkldata0(i,k+3);
end
hkldata0(:,4:7) = 1.3 * hkldata0(:,4:7) ./ (wav_max * [1,1,1,1]);
hkldata0(:,8:11) = hkldata0(:,8:11) ./ (int_max * [1,1,1,1]);
// Save maximum int and its wavelength
hkldata0(:,12:13)=[wav_max,int_max];
//
// Only keep if max ints not at ends
i=find( or( hkldata0(:,9:10) == 1 ,"c") );
hkldata0=hkldata0(i,:);
// Only keep if int(2)<int(1) and int(3)>int(4)
i=find( (hkldata0(:,8) < hkldata0(:,9)) & (hkldata0(:,10) > hkldata0(:,11)) );
hkldata0=hkldata0(i,:);
//
// Only keep if max int > 2000
i=find( hkldata0(:,13) > 5000 );
hkldata0=hkldata0(i,:);
//
// Only keep if max int wav < 2.6
i=find( hkldata0(:,12) < 2.6 );
hkldata0=hkldata0(i,:);
//
// Only keep if wav. spread >0.2 and < 1
dwav= abs( hkldata0(:,4) - hkldata0(:,7) );
i=find( (dwav > 0.2) & (dwav < 1) );
hkldata0=hkldata0(i,:);
//
// Only keep if wav. spread >0.2 and < 1
dwav= abs( hkldata0(:,4) - hkldata0(:,7) );
i=find( (dwav > 0.2) & (dwav < 1) );
hkldata0=hkldata0(i,:);





for i=1:size(hkldata0,1)
  plot(hkldata0(i,4:7),hkldata0(i,8:11))
end





mean(hkldata0(:,8:11),"r")


endfunction

