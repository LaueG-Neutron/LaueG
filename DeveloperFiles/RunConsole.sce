global ImageInfo LatticeInfo ProgDir OrientInfo
global ArgboxSet DisplaySet ModulateInfo

// ============ Main data processing =============
//function RunFindSpots(radius)
//function [UB,pixcen]=RunIndexSpots(pars,spots,bShowOutput)
//function [Linfo,Oinfo,text]=RunOrientSpots(Linfo,Oinfo,imode,pars)
//function [Linfo,sOut]=RunOrientTwins(Linfo,Oinfo, sFile, ...
//                           max_hkl,wav_min,dxy_max,niter,bFitCell)
//function stats=RunArgonneBoxes(base_name,bOldMode,bTwinMode,bModulateMode,twin_dist)
//function RunLaue4(base_name,nums, bTwin,bMod)
// ============ Displayed image =============
//function RunRejectsImage(base_names,xy_cen)
//function Info=RunImageInfo(base_name,bFull)
//function RunStripBack()
// ============ Ellipses and rejects =============
//function [ixy]=RunMakeEllipses(xy_cen,ell_pars,ellim,ixy_lim)
//function ipoint=RunMakeRejectEllipses(amult,cens,pars,ixy_lim)
//function RInfo=RunRejectsData(file_name,bVerbose)
// ============ Miscellaneous ==============
//function hkllist=RunGenHKLs(Linfo,wav_min,wav_max,d_min,d_max,bModulate)
//function v_nodals=RunIndexConic(Vconic,Nodals)
//function RunCyclopsPrep()
//function sTable=RunExeVersions()
// =========== Low level ==========
//function RunConsoleAbort(Command,SuccessComment,FailComment)
//function bOK=RunConsole(Command, SuccessComment,FailComment)
//function [sDos,sOut,sErr]=RunDOSCommand(Command)
//function RunConsolePopup(Command)


// ============ Main data processing =============

function RunFindSpots(radius)
global DisplayData
// Create input file for find_spots
  CloseAllFiles()
  [fd,ierr]=mopen("___laueg_find_spots.in","w")
  if(ierr ~= 0) then
    AbortBox("Unable to create input file")
  end
// Output base name for file
  mfprintf(fd,"%s\n",ImageInfo.basename)
// Output image properties
  mfprintf(fd,"%i %i\n",ImageInfo.numxy)
  mfprintf(fd,"%i\n",ImageInfo.basecounts)
  mfprintf(fd,"%i\n",radius)
// The pixel center and drum radius in mm
  mfprintf(fd,"%f %f %f\n",OrientInfo.pixcen,OrientInfo.drumrad)
// ipixcen=0 for VIVALDI/KOALA as pixcen not stable
  ipixcen=1
  if(ImageInfo.type == 1) then ipixcen=0; end
  mfprintf(fd,"%f %f %i\n",OrientInfo.pixsize,ipixcen)
// Output edges of image and exclusion regions
  mfprintf(fd,"%i %i %i %i\n",ImageInfo.edges)
  ncirc=size(ImageInfo.circs,1)
  mfprintf(fd,"%i\n",ncirc)
  for i=1:ncirc
    mfprintf(fd,"%i %i %i\n",ImageInfo.circs(i,:))
  end
  nrect=size(ImageInfo.rects,1)
  mfprintf(fd,"%i\n",nrect)
  for i=1:nrect
    mfprintf(fd,"%i %i %i %i\n",ImageInfo.rects(i,:))
  end
  mfprintf(fd,"%i\n",ImageInfo.type)
  mclose(fd)
// Run find_spots and abort if it crashes
  RunConsoleAbort("find_spots", ...
        "","Unable to find spots for image")
// Read the output file from find_spots
  [fd,ierr]=mopen("___laueg_find_spots.out","r")
  if(ierr ~= 0) then
    AbortBox("Unable to find output file")
  end
  spots2=mfscanf(-1,fd,"%f %f %f\n")
  mclose(fd)
// Delete the output file
  DeleteTempFile("___laueg_find_spots.out")
// Sort in terms of esds and find those with esds>1.5
  [temp,tags]=gsort(spots2(:,3),"r","d")
  itag=find((spots2(tags(:),3)<1.5),1)
  if isempty(itag) then itag=size(spots2,1); end
// Copy observed spots to DisplayData
  DisplayData.found_spots=spots2(tags(1:itag),:)
endfunction


function [UB,pixcen]=RunIndexSpots(pars,spots,bShowOutput)
// Complain and die if no spots to index
  if(spots == []) then
    AbortBox("Observed spots list is missing")
  end
// Create the index_spots input file
  CloseAllFiles()
  [fd,ierr]=mopen("___laueg_index_spots.in","w")
  if(ierr ~= 0) then
    AbortBox("Unable to create input file")
  end
//
// Unit cell and some geometry parameters
  mfprintf(fd,"%.4f %.4f %.4f %.3f %.3f %.3f\n",LatticeInfo(1).cell)
// Write the instrument parameters
  mfprintf(fd,"%i\n",ImageInfo.type)
  mfprintf(fd,"%i %i\n",ImageInfo.numxy)
  mfprintf(fd,"%f %f\n",OrientInfo.pixcen)
  mfprintf(fd,"%f %f\n",OrientInfo.pixsize)
  mfprintf(fd,"%f\n",OrientInfo.drumrad)
//
// Run the indexer using default parameters
  if(pars == []) then
    pars.hsq_max=15
    pars.pair_ang_min=15.0;  pars.pair_dang_max=1.5
    pars.azi_dang_max=2.0;   pars.min_azi_group=4
    pars.dot_cross_min=0.25
    pars.hkl_max1=5;         pars.hkl_max2=5
    pars.iselect=1
    pars.nmatch=20;          pars.ntest=25
  end
//
// Indexing control parameters
  mfprintf(fd,"%i\n",pars.hsq_max)
  mfprintf(fd,"%.2f %.3f\n",pars.pair_ang_min,pars.pair_dang_max)
  mfprintf(fd,"%.3f %i\n",pars.azi_dang_max,pars.min_azi_group)
  mfprintf(fd,"%.3f\n",pars.dot_cross_min)
  mfprintf(fd,"%i %i\n",pars.hkl_max1,pars.hkl_max2)
  mfprintf(fd,"%i\n",pars.iselect)
  mfprintf(fd,"%i %i\n",pars.nmatch,pars.ntest)
//
// Observed spots list
  mfprintf(fd,"%i\n",size(spots,1))
  if( spots ~= []) then
    mfprintf(fd,"%.1f %.1f %.1f\n",spots)
  end
//
// Don"t use any "guide spots"
  mfprintf(fd,"0\n")
//
  mclose(fd)
//
// Run index_spots, display a progress bar
// Show output and abort if indexing not successful
// If requested, always show the program output
  if( bShowOutput ) then
    RunConsoleAbort("index_spots", ...
          "Indexing program results","Indexing program failed")
  else
    RunConsoleAbort("index_spots","","Indexing program failed")
  end
// Read the UB and pixel center from the index_spots output file
  [fd,ierr]=mopen("___laueg_index_spots.out","r")
  if(ierr ~= 0) then
    AbortBox("Unable to find output file")
  end
  pixcen=mfscanf(fd,"%f %f\n")
  UB=mfscanf(3,fd,"%f %f %f\n")
  mclose(fd);
// Delete the output file
  DeleteTempFile("___laueg_index_spots.out")
endfunction


function [Linfo,Oinfo,text]=RunOrientSpots(Linfo,Oinfo,imode,pars)
//
// Run orient_spots in auto-mode (imode=1), manual modes (2-4),
// HKL index mode (5), or manual mode 2 with no output (6).
//
// Manual mode: (uses pars for parameter values)
//  Calculate hkls: imode=4
//    For Calc 1 method: pars.ihklgen=1; pars.hklgen=(HSQ_MAX,WAV_MIN)
//    For Calc 2 method: pars.ihklgen=2; pars.hklgen=(WAV_MIN,D_MIN)
//  Match hkls: imode=3 (after imode=4 to calc hkls)
//    imode=3; pars.nobs; pars.max_sep=(separation in mm)
//  Refine: imode=2 (after imode=4 & 3 to calc matched hkls)
//    pars.irefs=1 to refine (pix_cen,cell,rot,pix_siz,xtal_off,beam_vert)
//
// Create the input file
  CloseAllFiles()
  [fd,ierr]=mopen("___laueg_orient_spots.in","w")
  if(ierr ~= 0) then
    AbortBox("Unable to create input file")
  end
// Write mode to run orient_spots and instrument type
  if(imode == 6 ) then
    mfprintf(fd,"%i\n",2)
  else
    mfprintf(fd,"%i\n",imode)
  end
  mfprintf(fd,"%d\n",ImageInfo.type)
// Write the xtal information + UB
// If imode=5, always treat as triclinic P1
  if(imode == 5) then
    mfprintf(fd,"1 0\n")
  else
    mfprintf(fd,"%i %i\n",Linfo.ilatt,Linfo.icen)
  end
  mfprintf(fd,"%f %f %f %f %f %f\n",Linfo.cell)
  mfprintf(fd,"%f %f %f\n",Linfo.ub(1,1:3))
  mfprintf(fd,"%f %f %f\n",Linfo.ub(2,1:3))
  mfprintf(fd,"%f %f %f\n",Linfo.ub(3,1:3))
// Write instrument info common to instrument types
  mfprintf(fd,"%f\n",ImageInfo.phi)
  mfprintf(fd,"%i %i\n",ImageInfo.numxy)
  mfprintf(fd,"%f %f\n",Oinfo.pixcen)
  mfprintf(fd,"%f %f %f\n",Oinfo.xtaloff)
  mfprintf(fd,"%f %f %f\n",Oinfo.drumoff)
  mfprintf(fd,"%f\n",Oinfo.beamvert)
// Write instrument info specific to instrument types
  mfprintf(fd,"%f %f %f\n",Oinfo.pixsize,Oinfo.pixskew)
  mfprintf(fd,"%f\n",Oinfo.drumrad)
// Get number of observed spots to fit
  nspots=size(DisplayData.found_spots,1)
// If imode=2,3,4,6: write manual parameters and reset nspots
  if( or(imode == [2,3,4,6]) ) then
    mfprintf(fd,"%i %f %f\n",pars.ihklgen,pars.hklgen)
    mfprintf(fd,"%f\n",pars.max_sep)
    mfprintf(fd,"%i %i %i %i %i %i\n",matrix(pars.irefs,1,-1))
    nspots=pars.nobs
  end
// If imode=5, always uses 4 spots
  if(imode == 5) then
    nspots=4
  end
// Write the observed spots list
  mfprintf(fd,"%i\n",nspots)
  for i=1:nspots
    if(imode == 5) then
      mfprintf(fd,"%.2f %.2f %.3f %.3f %.3f\n",pars.xys(i,:),pars.hkls(i,:))
    else
      mfprintf(fd,"%.2f %.2f %.3f\n",DisplayData.found_spots(i,:))
    end
  end
// Finished writing input file
  mclose(fd)
//
// Run orient_spots, abort on failure depending on imode
  if(imode == 1) then
    bOK=RunConsole("orient_spots", "","Orientation program failed")
    if( ~bOK ) then
      DeleteTempFile("___laueg_orient_spots.out")
      AbortLaueG()
    end
  elseif(imode == 2) then
    RunConsole("orient_spots", ...
       "Orientation program completed","Orientation program failed")
  else
    RunConsole("orient_spots", "","Operation failed")
  end
//
// Open the output file
  [fd,ierr]=mopen("___laueg_orient_spots.out","r")
  if(ierr ~= 0) then
    AbortBox("Unable to find output file")
  end
//
// Overwite record items refined by orient_spots
  Linfo.cell=mfscanf(fd,"%f %f %f %f %f %f\n")
  Linfo.ub=mfscanf(3,fd,"%f %f %f\n")
  Oinfo.pixcen=mfscanf(fd,"%f %f\n")
  Oinfo.pixsize=mfscanf(fd,"%f %f\n")
  Oinfo.pixskew=mfscanf(fd,"%f\n")
  Oinfo.xtaloff([1,3])=mfscanf(fd,"%f %f\n")
  Oinfo.beamvert=mfscanf(fd,"%f\n")
//
// Read the results string
  text=mgetl(fd,1)
//
// Close and delete the file
  mclose(fd)
  DeleteTempFile("___laueg_orient_spots.out")
//
endfunction


function [Linfo,sOut]=RunOrientTwins(Linfo,Oinfo, sFile, ...
                           max_hkl,wav_min,dxy_max,niter,bFitCell)
//
// Create the input file
  CloseAllFiles()
  [fd,ierr]=mopen("___laueg_orient_twins.in","w")
  if(ierr ~= 0) then
    AbortBox("Unable to create input file")
  end
// File name of data file
  mfprintf(fd,"%s\n",sFile)
// Parameters for hkl generation
  mfprintf(fd,"%d %f %f %d\n",max_hkl,wav_min,dxy_max,niter)
// Fit all cell parameters, or keep the same
  mfprintf(fd,"%d\n",double(bFitCell))
// Instrument type
  mfprintf(fd,"%d\n",ImageInfo.type)
// Lattice info (assumes same for both twins)
  mfprintf(fd,"%i %i\n",Linfo(2).ilatt,Linfo(2).icen)
  mfprintf(fd,"%f %f %f %f %f %f\n",Linfo(2).cell)
// UB for Twin 1
  mfprintf(fd,"%f %f %f\n",Linfo(2).ub(1,1:3))
  mfprintf(fd,"%f %f %f\n",Linfo(2).ub(2,1:3))
  mfprintf(fd,"%f %f %f\n",Linfo(2).ub(3,1:3))
// UB for Twin 2
  mfprintf(fd,"%f %f %f\n",Linfo(3).ub(1,1:3))
  mfprintf(fd,"%f %f %f\n",Linfo(3).ub(2,1:3))
  mfprintf(fd,"%f %f %f\n",Linfo(3).ub(3,1:3))
// Instrument info common to instrument types
  mfprintf(fd,"%f\n",ImageInfo.phi)
  mfprintf(fd,"%i %i\n",ImageInfo.numxy)
  mfprintf(fd,"%f %f\n",Oinfo.pixcen)
  mfprintf(fd,"%f %f %f\n",Oinfo.xtaloff)
  mfprintf(fd,"%f %f %f\n",Oinfo.drumoff)
  mfprintf(fd,"%f\n",Oinfo.beamvert)
// Instrument info specific to instrument types
  mfprintf(fd,"%f %f %f\n",Oinfo.pixsize,Oinfo.pixskew)
  mfprintf(fd,"%f\n",Oinfo.drumrad)
// Finished writing input file
  mclose(fd)
//
// Run orient_twins and abort on failure
  bOK=RunConsole("orient_twins","","Orientation program failed")
  if( ~bOK ) then
    DeleteTempFile("___laueg_orient_twins.out")
    AbortLaueG()
  end
//
// Open and read the output file
  [fd,ierr]=mopen("___laueg_orient_twins.out","r")
  if(ierr ~= 0) then
    AbortBox("Unable to find output file")
  end
// Overwite cells & UBs
  Linfo(2).cell=mfscanf(fd,"%f %f %f %f %f %f\n")
  Linfo(3).cell=mfscanf(fd,"%f %f %f %f %f %f\n")
  Linfo(2).ub=mfscanf(3,fd,"%f %f %f\n")
  Linfo(3).ub=mfscanf(3,fd,"%f %f %f\n")
//
// Read the results string
  sOut=mgetl(fd,1)
//
// Close and delete the file
  mclose(fd)
  DeleteTempFile("___laueg_orient_twins.out")
//
endfunction


function stats=RunArgonneBoxes(base_name,bOldMode,bTwinMode,bModulateMode,twin_dist)
// Run argonne_boxes in normal, twinned, or modulated modes
//
// Generate the hkls to integrate
  if( bTwinMode ) then
// Generate hkls for twins 1 & 2, and the overlapped spots
    iopt_hkl=2
// Calculate twin hkls which contain (hkl, d-space, xy, wav, mult, imod)
// hkllist contains unpaired twin 1 spots, unpaired 2, interleaved paired spots
    hkllist=TwinGenHKLs(ArgboxSet.d_min,ArgboxSet.wav_min,twin_dist)
  else
// Generate hkllist which contains (hkl, d-space, xy, wav, mult, imod)
    hkllist=RunGenHKLs(LatticeInfo(1),ArgboxSet.wav_min,ArgboxSet.wav_max, ...
                                       ArgboxSet.d_min,1e4,bModulateMode)
    if( bModulateMode ) then
// Sort the list on hklm
      iopt_hkl=3
      [val tag]=gsort(hkllist(:,[1:3,9]),"lr","i")
      hkllist=hkllist(tag,:)
    else
// Sort the list on hkl
      iopt_hkl=1
      [val tag]=gsort(hkllist(:,1:3),"lr","i")
      hkllist=hkllist(tag,:)
    end
  end
//
// Open parameter file
// NB: Have to do this after RunGenHKLs() which closes all open files
  [fd,ierr]=mopen("___laueg_argonne_boxes.in","w")
  if(ierr ~= 0) then
    AbortBox("Unable to create input file")
  end
// Write header with hkl option (1=HKL, 2=HKLM) and lines used for parameters and hkls
  ncircs=size(ImageInfo.circs,1)
  nrects=size(ImageInfo.rects,1)
  nlines=12+ncircs+nrects   // nrects+4  rectangles are output
// Add 1 lines for pixsize & ipixcen
  if( ~bOldMode ) then nlines=nlines+1; end
  nhkls=size(hkllist,1)
  mfprintf(fd,"%i %i %i\n",iopt_hkl,nlines,nhkls)
//
// Output base name of image file
  mfprintf(fd,"%s\n",base_name)
// Output algorithm dependant parameters
  if( bOldMode ) then
    mfprintf(fd,"1\n")          // unused (I think?)
    mfprintf(fd,"%f %f %f\n",ArgboxSet.model_cutoff, ...
                ArgboxSet.strong_cutoff,ArgboxSet.core_mult)
  else
    mfprintf(fd,"%d %d %f %f %d\n", ...
             ArgboxSet.model_center,ArgboxSet.model_outer, ...
             ArgboxSet.pfrac_target,ArgboxSet.pfrac_uncertain/100, ...
             double(ArgboxSet.recenter) )
  end
//
// Output option independent parameters
  mfprintf(fd,"%f %f\n",ArgboxSet.cont_level,ArgboxSet.peak_mult)
  mfprintf(fd,"%f %f\n",ArgboxSet.neigh_toler,ArgboxSet.min_fill)
//
// The pixel center and drum radius in X pixels
  mfprintf(fd,"%f %f %f\n",OrientInfo.pixcen, ...
             OrientInfo.drumrad/abs(OrientInfo.pixsize(1)))
// If using the new mode, output X,Y pixel size and if PIXCEN is stable
  if( ~bOldMode ) then
// Pixcen for VIVALDI/KOALA is not stable
    ipixcen=1
    if(ImageInfo.type == 1) then ipixcen=0; end
    mfprintf(fd,"%f %f %i\n",OrientInfo.pixsize,ipixcen)
  end
//
// Output 3 crosstalk factors for instrument
  mfprintf(fd,"%f %f %f\n",ImageInfo.crosstalks)
//
// Output circular exclusion regions
  mfprintf(fd,"%i\n",ncircs)
  for i=1:ncircs
    mfprintf(fd,"%i %i %i\n",ImageInfo.circs(i,:))
  end
//
// Output number of rectanglar exclusion regions (plus image edges)
  mfprintf(fd,"%i\n",nrects+4)
// Output edges as rectangles (left,right,top bottom)
  mfprintf(fd,"%i %i %i %i\n",1,ImageInfo.edges(1),1,ImageInfo.numxy(2))
  mfprintf(fd,"%i %i %i %i\n",ImageInfo.edges(2), ...
                           ImageInfo.numxy(1),1,ImageInfo.numxy(2))
  mfprintf(fd,"%i %i %i %i\n",1,ImageInfo.numxy(1),1,ImageInfo.edges(3))
  mfprintf(fd,"%i %i %i %i\n",1,ImageInfo.numxy(1), ...
                           ImageInfo.edges(4),ImageInfo.numxy(2))
// Output list of exclusion rectangles
  for i=1:nrects
    mfprintf(fd,"%i %i %i %i\n",ImageInfo.rects(i,:))
  end
//
// Output reflection lines with hkl(m),x-pix,y-pix,wav,mult,tth
  for i=1:nhkls
    sin_th=0.5*hkllist(i,7)/hkllist(i,4)
    if(sin_th > 1.01) then
      mprintf("WARNING: sin(theta) > 1  for hkl %i %i %i\n",hkllist(i,1:3))
    end
    tth=2*asind( min(1.0,sin_th) )
    if(iopt_hkl == 1) then
      mfprintf(fd,"%i %i %i %.2f %.2f %.4f %i %.2f\n", ...
                                   hkllist(i,[1:3,5:8]),tth)
    else
      mfprintf(fd,"%i %i %i %i %.2f %.2f %.4f %i %.2f\n", ...
                                    hkllist(i,[1:3,9,5:8]),tth)
    end
  end
//
  mclose(fd)
// Run argonne_boxes, complain and abort if it crashes

  if( bOldMode ) then
    RunConsoleAbort("argonne_boxes_OLD", ...
       "","Aborting Integrate Spots due to an error from argonne_boxes")


  else
    RunConsoleAbort("argonne_boxes", ...
       "","Aborting Integrate Spots due to an error from argonne_boxes")
  end
//
// Read the output file from argonne_boxes
  [fd,ierr]=mopen("___laueg_argonne_boxes.out","r")
  if(ierr ~= 0) then
    AbortBox("Unable to find output file")
  end
  stats=mfscanf(14,fd,"%7d")
  mclose(fd)
//
// Delete the output file
  DeleteTempFile("___laueg_argonne_boxes.out")
//
endfunction


function RunLaue4(base_name,nums, bTwin,bMod)
//
// Create laue4 input file
  CloseAllFiles()
  [fd,ierr]=mopen("___laueg_laue4.in","w")
  if(ierr ~= 0) then
    AbortBox("Unable to create input file")
  end
// Write version number of file
  mfprintf(fd,"%i\n",6)
// Write file names to process
  nfiles=size(nums,1)
  mfprintf(fd,"%i\n",nfiles)           // number of file names
  mfprintf(fd,"%s\n",base_name);       // base name of all files
  mfprintf(fd,"%i\n",nums)             // the list of file numbers
// Write image size, instrument name, and experiment date
  mfprintf(fd,"%i %i\n",ImageInfo.numxy)
  mfprintf(fd,"%s\n",ImageInfo.instrum)
  temp=tokens(ImageInfo.date)
  mfprintf(fd,"%s\n",temp($)+"")
// Write correction on/off switches
  mfprintf(fd,"%s\n","T")       // file scale flag
  mfprintf(fd,"%s\n","T")       // wav. distrib. flag
  mfprintf(fd,"%s\n","T")       // efficiency flag
  mfprintf(fd,"%s\n","F")       // extinction flag
  mfprintf(fd,"%s\n","F")       // harmonic overlap flag
  mfprintf(fd,"%s\n","F")       // absorb flag
// Convert from laueg to laue4 numbering of cell types
  convert=[12,10,5,4,1,8,7,6]
  icryst=convert(LatticeInfo(1).ilatt)
// Write cell type, 1 (merge Friedels), and cell type again
  mfprintf(fd,"%i %i %i\n",icryst,1,icryst)
// Write centering and cell dimensions (for *.cif & *.inf files)
  mfprintf(fd,"%i %.3f %.3f %.3f %.2f %.2f %.2f\n", ...
                    LatticeInfo(1).icen+1,LatticeInfo(1).cell)
// Write options and parameters for various corrections
  mfprintf(fd,"2 10 300 700\n")     // 10 point spline & non-parametric
  mfprintf(fd,"2\n")                // quadratic efficiency
  mfprintf(fd,"1 0.1\n")            // extinction type & est. max. corr.
  mfprintf(fd,"1 0 0.1\n")          // absorb
// Write minimum number and relative intensity for refined sequences
  mfprintf(fd,"2 0.75\n")
// Write X & Y and wavelength limits
  mfprintf(fd,"%i %i  %i %i\n",ImageInfo.edges)
  mfprintf(fd,"%.4f %.4f\n",ArgboxSet.wav_min,ArgboxSet.wav_min*2.0)
// Write multiplier for cell dimensions
  mfprintf(fd,"1.00\n")
// Write all_rel_err,weak_sig_mult parameters
  if(ImageInfo.instrum == "KOALA") then
    mfprintf(fd,"1.0 0.005 1.0\n")
  elseif(ImageInfo.instrum == "IMAGINE") then
    mfprintf(fd,"0.7 0.02 1.0\n")
  else
    mfprintf(fd,"1.3 0.03 1.0\n")
  end
// For blue plates:  ip_thick=0.69
  ip_thick=1.0
  ip_width=OrientInfo.drumrad/abs(OrientInfo.pixsize(2))
// Write ip width and thickness
  mfprintf(fd,"%.1f %.3f\n",ip_width,ip_thick)
// Write wav_err1,wav_err2,exti_err
  mfprintf(fd,"0.2 0.2 0.25 0\n")
// Write iverbose,harm_cutoff_mult
  mfprintf(fd,"1 1.0\n")
//
// Add flags for twins or modulation, "0 0" if neither
  mfprintf(fd,"%d %d\n",double([bTwin,bMod]))
// Add modulation info
  if( bMod ) then
    [nidx,nvec]=size(ModulateInfo.mults)
    mfprintf(fd,"%d %d\n",nvec,nidx)
    sout=strcat(string(ModulateInfo.vecs)," ","c")
    mfprintf(fd,"%s\n",sout)
    sout=strcat(string(ModulateInfo.mults)'," ","c")
    mfprintf(fd,"%s\n",sout)
  end
// Close laue4 input file
  mclose(fd)
//
// Run laue4 in a popup console, and ignore any errors
  RunConsolePopup("laue4_exe")
endfunction



// ============ Displayed image =============

function RunRejectsImage(base_names,xy_cen)
global ReducedImage
// Run rejects_image to load ReducedImage[]
//
// Create the input file
   CloseAllFiles()
  [fd,ierr]=mopen("___laueg_rejects_image.in","w");
  if(ierr ~= 0) then
    AbortBox("Unable to create input file")
  end
  nfiles=size(xy_cen,1)
  mfprintf(fd,"%d\n",nfiles)
  for i=1:nfiles
    mfprintf(fd,"%s\n",base_names(i)+".tif")
    mfprintf(fd,"%.1f %.1f\n",xy_cen(i,1:2))
  end

  mclose(fd);
//
// Run rejects_image and abort if it crashes
  RunConsoleAbort("rejects_image","","Unable to make reject spots image")
//
  DeleteTempFile("___laueg_rejects_image.in")
//
// Read the rejects image file into ReducedImage[]
// Open the rejects image file, complain and die on an error
  [fd,ierr]=mopen("___laueg_rejects_image.img","rb")
  if(ierr ~= 0) then
    AbortBox("Unable to find rejects image file")
  end
// Read the *.img file into ReducedImage[] as a matrix of doubles
  ReducedImage=double(matrix(mgeti(1000*500,"us",fd),1000,500))
// Close and delete the *.img file
  mclose(fd)
  DeleteTempFile("___laueg_rejects_image.img")
endfunction


function Info=RunImageInfo(base_name,bFull)
//
// Run image_info and return the values in Info, which is a subset of ImageInfo.
// If bFull == %t, use image data to find central hole and base intensity.
//
// Create the input file
   CloseAllFiles()
  [fd,ierr]=mopen("___laueg_image_info.in","w")
  if(ierr ~= 0) then
    AbortBox("Unable to create input file")
  end
  mfprintf(fd,"%s\n",base_name)
  if( bFull ) then
    mfprintf(fd,"2\n")
  else
    mfprintf(fd,"1\n")
  end
  mclose(fd)
//
// Run image_info and abort if it crashes
  RunConsoleAbort("image_info", ...
        "","Unable to read information header for image")
//
// Read the output text file into Info.*
// Open the output text file, complain and die on an error
  [fd,ierr]=mopen("___laueg_image_info.out","r")
  if(ierr ~= 0) then
    AbortBox("Unable to find output text file")
  end
// Copy the output text file results to a clean ImageInfo
  Info=[]
  Info.basename=base_name
  Info.comment=mgetl(fd,1)
  Info.sample=mgetl(fd,1)
  Info.user=mgetl(fd,1)
  Info.host=mgetl(fd,1)
  Info.date=mgetl(fd,1)
  Info.phi=mfscanf(fd,"%f\n")
  Info.expo=mfscanf(fd,"%i\n")
  Info.numxy=mfscanf(fd,"%i %i\n")
  Info.offset=mfscanf(fd,"%i\n")
  Info.startxy=mfscanf(fd,"%i %i\n")    // ignored
  Info.resol=mfscanf(fd,"%i\n")         // ignored
  Info.hole=mfscanf(fd,"%i %i %i\n")
  Info.basecounts=mfscanf(fd,"%i\n")
  Info.reduce=mfscanf(fd,"%i\n")
// Additional parameters for KOALA2
  tmp=[0,0,0,0]
  if(Info.host == 'KOALA2') then
    tmp=mfscanf(fd,"%i %i %i %i\n")
    Info.adc=tmp(1)
    Info.aper=tmp(2:3)
    Info.kappa=tmp(4)
  end
// Close and delete the output text file
  mclose(fd)
  DeleteTempFile("___laueg_image_info.out")
//
// If strings are empty, return [] to indicate an invalid file
  if(Info.user+Info.host+Info.date == "") then
    Info=[]
    return
  end
//
endfunction


function RunStripBack()
global FileImage StripImage RawImage
//
// Write input file for strip_back.exe
  CloseAllFiles()
  [fd,ierr]=mopen("___laueg_strip_back.in","w")
  if(ierr ~= 0) then
    AbortBox("Unable to create input file")
  end
// Output base filename and image properties
  mfprintf(fd,"%s\n",ImageInfo.basename)
  mfprintf(fd,"%i %i\n",ImageInfo.numxy)
  mfprintf(fd,"%i\n",ImageInfo.basecounts)
  mfprintf(fd,"%i %i\n",OrientInfo.pixcen)
  mfprintf(fd,"%f %f\n",OrientInfo.pixsize)
  mfprintf(fd,"%f\n",OrientInfo.drumrad)
// Pixcen is unstable for VIVALDI/KOALA
  ipixcen=1
  if(ImageInfo.type == 1) then ipixcen=0; end
  mfprintf(fd,"%i\n",ipixcen)
// Output image edges and exclusion regions
  mfprintf(fd,"%i %i %i %i\n",ImageInfo.edges)
  ncirc=size(ImageInfo.circs,1)
  mfprintf(fd,"%i\n",ncirc)
  for i=1:ncirc
    mfprintf(fd,"%i %i %i\n",ImageInfo.circs(i,:))
  end
  nrect=size(ImageInfo.rects,1)
  mfprintf(fd,"%i\n",nrect)
  for i=1:nrect
    mfprintf(fd,"%i %i %i %i\n",ImageInfo.rects(i,:))
  end
  mclose(fd)
//
// Run strip_back and abort if it crashes
  RunConsoleAbort("strip_back", ...
        "","Unable to strip background from image")
//

// Load RawImage from the background stripped raw image file
  [fd,ierr]=mopen("___laueg_strip_back.raw","rb");
  if(ierr ~= 0) then
    AbortBox("Cannot find stripped image file")
  end
// Read the *.raw file into RawImage[] as a matrix of double
  numx=ImageInfo.numxy(1)
  numy=ImageInfo.numxy(2)
  RawImage=matrix(mget(numx*numy,"f",fd),numx,numy);
// Close and (possibly) delete the *.raw file
  mclose(fd);
  DeleteTempFile("___laueg_strip_back.raw")
//
// Copy RawImage[] to StripImage[]
  StripImage=RawImage;
//
endfunction



// ============ Ellipses and rejects =============

function [ixy]=RunMakeEllipses(xy_cen,ell_pars,ellim,ixy_lim)
// Run make_ellipses to find pixels that lie on the ellipse edge
// Create the input file
  CloseAllFiles()
  [fd,ierr]=mopen("___laueg_make_ellipses.in","w")
  if(ierr ~= 0) then
    AbortBox("Unable to create input file")
  end
  nspots=size(xy_cen,1)
  for ell=ellim
    mfprintf(fd,"%i %.2f %i %i %i %i\n",nspots,ell,ixy_lim)
    mfprintf(fd,"%.1f %.1f %.5f %.5f %.5f\n",[xy_cen ell_pars])
  end
  mclose(fd)
// Run make_ellipses
// Show output and abort if not successful
  RunConsoleAbort("make_ellipses", ...
              "","Program failed, unable to draw ellipses")
// Open the output file
  [fd,ierr]=mopen("___laueg_make_ellipses.out","r")
  if(ierr ~= 0) then
    AbortBox("Unable to find output file")
  end
// Read all pixel X,Y values and the ellipse index
  ixy=mfscanf(-1,fd,"%i %i\n");
  mclose(fd)
// Remove the index from ixy
  ipos=find(ixy(:,1) == 99999, 1)
  ixy=ixy(1:ipos-1,1:2)
// Delete the output file
  DeleteTempFile("___laueg_make_ellipses.out")
endfunction


function ipoint=RunMakeRejectEllipses(amult,cens,pars,ixy_lim)
//
// Create the input file
  CloseAllFiles()
  [fd,ierr]=mopen("___laueg_make_ellipses.in","w")
  if(ierr ~= 0) then
    AbortBox("Unable to create input file")
  end
// Write the ellipse parameters
  [nspots nmults]=size(amult)
  for i=1:nspots
    for ell=amult(i,:)
      mfprintf(fd,"%i %.2f %i %i %i %i\n",1,ell,ixy_lim(i,:))
      mfprintf(fd,"%.1f %.1f %.5f %.5f %.5f\n",cens(i,:),pars(i,:))
    end
  end
  mclose(fd)
//
// Run make_ellipses, show output and abort if failed
  RunConsoleAbort("make_ellipses", ...
             "","Program failed, unable to draw ellipses")
//
// Read the output file
  [fd,ierr]=mopen("___laueg_make_ellipses.out","r")
  if(ierr ~= 0) then
    AbortBox("Unable to find output file")
  end
// Read all X,Y values and the ellipse index
  ixy=mfscanf(-1,fd,"%i %i\n")
  mclose(fd)
// Extract the index from the list
  index=ixy( find(ixy(:,1) == 99999) , 2)
//
// Add offsets due to tiling to the pixel positions
  ifirst=1
  for i_tile=1:nspots
// Get the end ixy for this ellipse
    ilast=index(i_tile*nmults)
// Calculate the X,Y offsets of the "tiled" image
    ix_tile=round(modulo(i_tile-0.1,10))
    iy_tile=1+floor((i_tile-1)/10)
    idx=ix_tile*100-100 +1
    idy=500-iy_tile*100 +1
// Add the image offsets to ixy
    ixy(ifirst:ilast,1)=ixy(ifirst:ilast,1)+idx-ixy_lim(i_tile,1)
    ixy(ifirst:ilast,2)=ixy(ifirst:ilast,2)+idy-ixy_lim(i_tile,3)
// Update the start for the next ellipse
    ifirst=ilast+1
  end
//
// Calculate the indices in the main image of the pixel positions
  ipoint=500*ixy(1:ilast,1)-ixy(1:ilast,2)
//
// Delete the output file
  DeleteTempFile("___laueg_make_ellipses.out")
endfunction


function RInfo=RunRejectsData(file_name,bVerbose)
//
// Initialise returned structure
  RInfo=struct("itwin",[],"base_name",[],"hkl",[], ...
     "xycen",[],"efh",[],"cont",[],"comment",[],"marked",[])
//
// Create the input file
  CloseAllFiles()
  [fd,ierr]=mopen("___laueg_rejects_data.in","w")
  if(ierr ~= 0) then
    AbortBox("Unable to create input file")
  end
  mfprintf(fd,"%s\n",file_name)
  mclose(fd)
//
// Run rejects_data
  RunConsoleAbort("rejects_data", ...
         "","Program failed, unable to load reject spots")
//
// Read the output file and save the data for spots in the
// input file that match a record in the corresponding *.ell file
  [fd,ierr]=mopen("___laueg_rejects_data.out","r")
  if(ierr ~= 0) then
    AbortBox("Unable to find output file")
  end
//
  nread=0
  nrej=0
  while (%t) then
// Read most parameters
    base=stripblanks(mgetstr(40,fd))
    [nscan,h,k,l,itwin,x,y,e,f,h2,c1,c2,c3,itype]= ...
          mfscanf(fd,"%i %i %i %i %f %f %f %f %f %f %f %f %i\n")
    if(nscan == -1) then
      break
    end
    if(nscan ~= 13) then
      AbortBox("Error reading reject spots data")
    end
    nread=nread+1
// Skip this spot if unable to find in *.ell file
    if(itype < 0) then
      continue
    end
// Save the parameters
    nrej=nrej+1
    RInfo.itwin(nrej)=itwin
    RInfo.base_name(nrej)=base
    RInfo.hkl(nrej,:)=[h,k,l]
    RInfo.xycen(nrej,:)=[x,y]
    RInfo.efh(nrej,:)=[e,f,h2]/1000
    RInfo.cont(nrej,:)=[c1,c2,c3]
// Read comment that trails the parameters
    RInfo.comment(nrej)=mgetl(fd,1)
// "Unmark" all spots
    RInfo.marked(nrej)=%f
  end
  mclose(fd)
//
// Delete the output file
  DeleteTempFile("___laueg_rejects_data.out")
//
// In verbose mode, output a warning if required
  if( bVerbose ) then
    if ( nread == 0 ) then
      WarnBox("No spots to display")
    elseif ( nrej == 0 ) then
      WarnBox(msprintf("Read %i spots, found no matches",nread))
    elseif ( nrej ~= nread ) then
      WarnBox(msprintf("Read %i spots, found %i matches",nread,nrej))
    end
  end
//
endfunction


// ============ Miscellaneous ==============

function hkllist=RunGenHKLs(Linfo,wav_min,wav_max,d_min,d_max,bModulate)
//
// Create the gen_hkls input file and returns hkllist, an array
// containing [hkl, d-space(4), xy, wav, mult, imod, Rhkl] where
// Rhkl are real-valued Miller indices used for satellite spots.
// "Linfo" is LatticeInfo(1) or some related structure
//
// Sanity check on unindexed lattice
  if( Linfo.level < 2 ) then
    hkllist=[]
    return
  end
//
// Write output file for Fortran program
  CloseAllFiles()
  [fd,ierr]=mopen("___laueg_gen_hkls.in","w")
  if(ierr ~= 0) then
    AbortBox("Unable to create input file")
  end
// Write the instrument type
  mfprintf(fd,"%i\n",ImageInfo.type)
// Write the centering, wav_min,d_min & UB
  mfprintf(fd,"%i\n",Linfo.icen)
  mfprintf(fd,"%f %f %f %f\n",wav_min,wav_max,d_min,d_max)
  mfprintf(fd,"%f %f %f\n",Linfo.ub(1,1:3))
  mfprintf(fd,"%f %f %f\n",Linfo.ub(2,1:3))
  mfprintf(fd,"%f %f %f\n",Linfo.ub(3,1:3))
// Write phi, image edges, pixel center
  mfprintf(fd,"%f\n",ImageInfo.phi)
  mfprintf(fd,"%i %i %i %i\n",ImageInfo.edges)
  mfprintf(fd,"%f %f\n",OrientInfo.pixcen)
// Write xtal & drum offsets, beam vertical angle
  mfprintf(fd,"%f %f %f\n",OrientInfo.xtaloff)
  mfprintf(fd,"%f %f %f\n",OrientInfo.drumoff)
  mfprintf(fd,"%f\n",OrientInfo.beamvert)
// Write pixel size & skew, and drum-radius
  mfprintf(fd,"%f %f %f\n",OrientInfo.pixsize,OrientInfo.pixskew)
  mfprintf(fd,"%f\n",OrientInfo.drumrad)
// Write number of, and list of, HKL special-condition vector-pairs
  n=size(Linfo.special_pairs,1)
  mfprintf(fd,"%i\n",n)
  for i=1:n
    mfprintf(fd,"%i %i %i   %i %i %i\n",Linfo.special_pairs(i,:))
  end
//
// Write modulation information
  if( bModulate ) then
// Write number of, and list of, satellite hkls
    mod_hkls=ModulateInfo.mults * ModulateInfo.vecs
    n=size(mod_hkls,1)
    mfprintf(fd,"%i\n",n)
    for i=1:n
      mfprintf(fd,"%f %f %f\n",mod_hkls(i,:))
    end
// Write d-spacing limit of main reflection for calculating satellites
    mfprintf(fd,"%f\n",ModulateInfo.d_min)
  else
// Write an empty satellite list to turn off modulation
    mfprintf(fd,"0\n")
  end
// Close the input file
  mclose(fd)
//
// Run gen_hkls, complain and abort if it crashes
// Show output and abort if not successful
  RunConsoleAbort("gen_hkls", ...
              "","Program failed, unable to generate HKL list")
//
// Read the output file
  CloseAllFiles()
  [fd,ierr]=mopen("___laueg_gen_hkls.out","r")
  if(ierr ~= 0) then
    AbortBox("Unable to find output file")
  end
  buff=mgetl(fd,-1)
  mclose(fd)
// Decode buff to hkllist
  [hkllist,bErr]=strmat2reals(buff)
  if(size(hkllist,1) == 0) then hkllist=[]; end
// Delete the output file
  DeleteTempFile("___laueg_gen_hkls.out")
endfunction


function v_nodals=RunIndexConic(Vconic,Nodals)
// Complain and die if no spots to index
  if(DisplayData.found_spots == []) then
    AbortBox("Observed spots list is missing")
  end
// Create the index_conic input file
  CloseAllFiles()
  [fd,ierr]=mopen("___laueg_index_conic.in","w")
  if(ierr ~= 0) then
    AbortBox("Unable to create input file")
  end
// Write the geometry information for the image
  mfprintf(fd,"%f\n",ImageInfo.phi)
  mfprintf(fd,"%i %i\n",ImageInfo.numxy)
  mfprintf(fd,"%f %f\n",OrientInfo.pixcen)
  mfprintf(fd,"%f %f\n",OrientInfo.pixsize)
  mfprintf(fd,"%f\n",OrientInfo.pixskew)
  mfprintf(fd,"%f %f %f\n",OrientInfo.xtaloff)
  mfprintf(fd,"%f\n",OrientInfo.beamvert)
  mfprintf(fd,"%f %f %f\n",OrientInfo.drumoff)
  mfprintf(fd,"%f\n",OrientInfo.drumrad)
// Write the conic vector
  mfprintf(fd,"%.4f %.4f %.4f\n",Vconic)
// Write the indices of the two nodal spots
  mfprintf(fd,"%i %i\n",Nodals(1),Nodals(2))
// Write the observed spots list
  mfprintf(fd,"%i\n",size(DisplayData.found_spots,1))
  mfprintf(fd,"%.1f %.1f %.1f\n",DisplayData.found_spots)
// Close the input file
  mclose(fd)
//
// Run index_conic, display a progress bar
// Show output and abort if indexing not successful
  RunConsoleAbort("index_conic", ...
          "Conic Indexing Results","Indexing program failed")
//
// Open and read the index_conic output file
  [fd,ierr]=mopen("___laueg_index_conic.out","r")

  if(ierr ~= 0) then
    AbortBox("Unable to find output file")
  end
// Read the two nodal vectors
  v_nodals(1,1:3)=mfscanf(1,fd,"%f %f %f\n")
  v_nodals(2,1:3)=mfscanf(1,fd,"%f %f %f\n")
// Read, but ignore, the calculated pixel positions
  xy_conic=mfscanf(-1,fd,"%f %f\n")
  mclose(fd)
//
// Delete the output file
  DeleteTempFile("___laueg_index_conic.out")
endfunction


function RunCyclopsPrep()
// Run laueG_cyclops in a popup console, and ignore any errors
  RunConsolePopup("laueg_cyclops")
endfunction


function sTable=RunExeVersions()
// Create table of *.exe file names and version dates
//
// Create a list of *.exe file names used for RunConsole, etc.
  str="findstr /i RunConsole """+ProgDir+"*.sce"""
  [sOut,bOK]=dos(str)
  sOut=strsubst(sOut,"''","""")
  sOut=strsubst(sOut,"RunConsole(""","XYZZY")
  sOut=strsubst(sOut,"RunConsoleAbort(""","XYZZY")
  sOut=strsubst(sOut,"RunConsolePopup(""","XYZZY")
  sOut=sOut( grep(sOut,"XYZZY") )
  sOut=strsubst(sOut,"/.*XYZZY/","","r")
  sOut=strsubst(sOut,"/"".*/","","r")
  progs=unique(convstr(sOut))
  if(progs(1) == "") then progs(1)=[]; end
//
// Delete all LaueG input files so executables crash immediately
  dos("del ___laueg_*.in")
//
// Run executables and collect version dates in sTable[]
  sTable=["",""]
  for prog=progs'
    [sDos,sOut,sErr]=RunDOSCommand(prog)
    sTable($+1,:)=[prog,"("+strsubst(sOut(1),"/.* /","","r")]
  end
  sTable(1,:)=[]
//
endfunction


// =========== Low level ==========

function RunConsoleAbort(Command,SuccessComment,FailComment)
// Run an executable and abort on a failure
//
  bOK=RunConsole(Command,SuccessComment,FailComment)
  if( ~bOK ) then
    AbortLaueG()
  end
endfunction


function bOK=RunConsole(Command, SuccessComment,FailComment)
global ImageFigure
// Run "Command" from ProgDir, and if the comment strings are
// not empty output some information at the end
//
// Complain and die if using an UNC path name
  if( grep(DataDir,"\\") ~= [] ) then
    AbortBox(["Unable to use UNC path names such as "+DataDir;""; ...
    "Copy files to a local drive, or map the drive to a letter (i.e. Z:)"])
  end
//
// Run Command in a hidden dos-box and write summary to log file
  WriteLogFile("RunConsole: "+Command)
  [sDos,sOut,sErr]=RunDOSCommand(Command)
  bOK= ~isempty(grep(sOut,"SUCCESSFUL COMPLETION"))
  if( ~bOK ) then WriteLogFile("RunConsole: FAILED"); end
//
// Ignore sDos as it keeps giving: "The system cannot find the path specified."
//
// Shouldn't happen, but check anyway
  if( bOK & (sErr ~= []) ) then
    BugBox("bOK=True and sErr="+sErr)
  end
//
// If success, output success message (if it exists)
  if( bOK ) then
    if(SuccessComment ~= "") then
      stmp=SuccessComment
      if(sOut ~= []) then
        stmp=[stmp;"";"Program output:";sOut]
      end
      InfoBox(stmp)
    end
  else
//
// If failure, output the failure message (if it exists)
// and then any program output or error messages
    if( ~bOK & (FailComment ~= "") ) then
      stmp=FailComment
      if(sOut ~= []) then
        stmp=[stmp;"";"Program output:";sOut]
      end
      if(sErr ~= []) then
        stmp=[stmp;"";"System error for command:";"   "+[ProgDir+Command;sErr]]
      end
      ErrorBox(stmp)
    end
  end
//
endfunction


function [sDos,sOut,sErr]=RunDOSCommand(Command)
//
  outfile="___laueg_dos.out"
  errfile="___laueg_dos.err"
  Q='""'
// NB: PRINT & WRITE(6) goes to "1", WRITE(0) goes to "2"
  dos_line=Q+Q+Q+ProgDir+Command+Q+" 1>"+outfile+" 2>"+errfile+Q+Q
  [sDos,bOK]=dos(dos_line)
  sOut=mgetl(outfile)
  sErr=mgetl(errfile)
  mdelete(outfile)
  mdelete(errfile)
//
endfunction


function RunConsolePopup(Command)
// Run "Command" in a popup console in ProgDir, pause, then exit
//
  WriteLogFile("RunConsolePopup: "+Command)
//
  Q='""'
  dos_line="start /wait cmd /c "+Q+ProgDir+"run_pause.bat"+Q+" "+Command
  dos(dos_line)
//
endfunction
