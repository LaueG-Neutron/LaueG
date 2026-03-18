global InstrumSetups InstrumDefault ImageInfo OrientInfo ProgDir

// ====== Load Instrument Setup given clues in Info ======
//function LoadInstrumSetup(Info)
// ====== Load main structures with dummy values  ======
//function LoadDummySetup()
//function DummyLatticeInfo()
//function DummyDisplaySet()
//function DummyArgboxSet()
// ======== Load info/settings with instrument defaults =======
//function LoadSetupOrient()
//function LoadSetupSettings()
// ====== Routines to help choose Instrument Setup ======
//function iSetup=ChooseInstrumSetup(hostname,numxy,sdate)
//function Match=MatchSetups(hostname,numxy,sdate)
//function n=ChooseInstrumDefault()


// ====== Load Instrument Setup given clues in Info ======

function LoadInstrumSetup(Info)
global ImageInfo
//
// Info is a version of ImageInfo which only uses the
// fields .host, .numxy, and .date (plus .hole if "type" = 1).
//
// Choose the instrument setup to use
  iSetup=ChooseInstrumSetup(Info.host,Info.numxy,Info.date)
  setup=InstrumSetups(iSetup)
//
// If dummy instrument, or setup has changed
  if(setup.type == 0) then
    bChanged=%t
  else
    bChanged=( ImageInfo.setup ~= setup.setup )
  end
//
// ==== Create a new ImageInfo ====
// Copy fields from Info to ImageInfo
  ImageInfo=Info
//
// Add fields to ImageInfo from the instrument setup
  ImageInfo.type=setup.type
// ImageInfo.type:
//  0   dummy instrument
//  1   KOALA/VIVALDI, drum IP, unstable XY, non-zero basecounts
//  2   IMAGINE, drum IP, stable XY, zero basecounts
//  3   CYCLOPS, octagonal CCD
//  4   KOALA2, drum IP, stable XY, nonzero basecounts
//
  ImageInfo.instrum=setup.instrum
  ImageInfo.setup=setup.setup
  ImageInfo.detector=setup.detector
  ImageInfo.edges=setup.edges
  ImageInfo.circs=[setup.hole;setup.circs]
  ImageInfo.rects=setup.rects
  ImageInfo.crosstalks=setup.crosstalks

//
// For KOALA/VIVALDI apply a shift in pixel x,y
  if(ImageInfo.type == 1) then
    dx=Info.hole(1)-setup.hole(1)
    dy=Info.hole(2)-setup.hole(2)
    ImageInfo.edges=ImageInfo.edges
    ImageInfo.circs(:,1)=ImageInfo.circs(:,1)+dx
    ImageInfo.circs(:,2)=ImageInfo.circs(:,2)+dy
    if( ImageInfo.rects ~= [] ) then
      ImageInfo.rects(:,1)=ImageInfo.rects(:,1)+dx
      ImageInfo.rects(:,2)=ImageInfo.rects(:,2)+dx
      ImageInfo.rects(:,3)=ImageInfo.rects(:,3)+dy
      ImageInfo.rects(:,4)=ImageInfo.rects(:,4)+dy
    end
  end
//
// For IMAGINE/CYCLOPS basecounts is zeroed
  if( or(ImageInfo.type == [2,3]) )
    ImageInfo.basecounts=0
  end
//
// For unstable VIVALDI/KOALA always reset pixcen
///  if(ImageInfo.type == 1)
///    OrientInfo.pixcen=setup.pixcen
///  end
//
// Reload OrientInfo if setup has changed, and update
// log file if not the dummy setup
  if( bChanged ) then
    LoadSetupOrient()
    if(ImageInfo.type ~= 0) then
      WriteLogFile("Loading Instrument Setup: "+ImageInfo.setup)
    end
  end
//
endfunction


// ====== Load main structures with dummy values  ======

function LoadDummySetup()
// Load "dummy" setup from instrument setups
  Info=[]
  Info.host= "laueg_dummy"
  Info.numxy=[4000 2000]
  Info.date="00:00:00 1/1/2010"
  LoadInstrumSetup(Info)
endfunction


function DummyLatticeInfo()
global LatticeInfo ModulateInfo
//
// Load the boring default cell
  LatticeInfo=[]
  LatticeInfo.itype=1
  LatticeInfo.level=0
  LatticeInfo.ilatt=1
  LatticeInfo.icen=0
  LatticeInfo.cell=[10,10,10,90,90,90]
  LatticeInfo.ub=[]
  LatticeInfo.special_pairs=[]
// Load no modulation
  ModulateInfo.vecs=[]
  ModulateInfo.mults=[]
  ModulateInfo.d_min=1
  ModulateInfo.d_max=10
endfunction


function DummyDisplaySet()
global DisplaySet
// Image & spots display params
  DisplaySet.blur_main=%f
  DisplaySet.invert_grey=%f
  DisplaySet.strip_back=%f
  DisplaySet.d_min=0.7
  DisplaySet.d_max=99
  DisplaySet.wav_min=0.85
  DisplaySet.wav_max=99
  DisplaySet.ell_back=%f
endfunction


function DummyArgboxSet()
global ArgboxSet
// Argonne_boxes params
  ArgboxSet.strong_cutoff=2
  ArgboxSet.model_cutoff=2
  ArgboxSet.cont_level=0.1
  ArgboxSet.core_mult=1
  ArgboxSet.peak_mult=4
  ArgboxSet.neigh_toler=10
  ArgboxSet.min_fill=0.8
  ArgboxSet.d_min=0.8
  ArgboxSet.wav_min=0.85
  ArgboxSet.wav_max=99
  ArgboxSet.model_center=20
  ArgboxSet.model_outer=10
  ArgboxSet.pfrac_target=0.8
  ArgboxSet.pfrac_uncertain=25
  ArgboxSet.recenter=%f
endfunction


// ======== Load info/settings with instrument defaults =======

function LoadSetupOrient()
global OrientInfo
//
// Find current setup in the instrument list
  iSetup=0
  for i=1:size(InstrumSetups,1)
    if(InstrumSetups(i).setup==ImageInfo.setup) then
      iSetup=i
    end
  end
//
// Reload OrientInfo values from the instrument defaults
  setup=InstrumSetups(iSetup)
  OrientInfo.pixcen=setup.pixcen
  OrientInfo.pixsize=setup.pixsize
  OrientInfo.pixskew=setup.pixskew
  OrientInfo.xtaloff=setup.xtaloff
  OrientInfo.drumoff=setup.drumoff
  OrientInfo.beamvert=setup.beamvert
  OrientInfo.drumrad=setup.drumrad
//
endfunction


function LoadSetupSettings()
global DisplaySet ArgboxSet
//
// Find current setup in the instrument list
  iSetup=0
  for i=1:size(InstrumSetups,1)
    if(InstrumSetups(i).setup==ImageInfo.setup) then
      iSetup=i
    end
  end
//
// Reload the settings contained in the instrument defaults
  setup=InstrumSetups(iSetup)
  DisplaySet.d_min=setup.display_dspaces(1)
  DisplaySet.d_max=setup.display_dspaces(2)
  DisplaySet.wav_min=setup.display_wavs(1)
  DisplaySet.wav_max=setup.display_wavs(2)
  ArgboxSet.d_min=setup.argbox_dwavs(1)
  ArgboxSet.wav_min=setup.argbox_dwavs(2)
  ArgboxSet.wav_max=setup.argbox_dwavs(3)
//
endfunction


// ====== Routines to help choose Instrument Setup ======

function iSetup=ChooseInstrumSetup(hostname,numxy,sdate)
//
// Select (hopefully automatically) a setup compatible with the current image
//
// Special case if the "dummy" instrument
  if( hostname == InstrumSetups.host(1) ) then
    iSetup=1
    return
  end
//
// Get merit of matches between image and known instruments
  Match=MatchSetups(hostname,numxy,sdate)
//
// If there is an exact match, load it and return
  iMatch=find(Match == 1)
  if( iMatch ~= [] ) then
    iSetup=iMatch(1)
    return
  end
//
// Makes a vector of strings of compatible setup names
  iMatch=find(Match > 0)
  names=[]
  for i=iMatch
    names($+1)=InstrumSetups(i).setup;
  end
//
// Give up if no clue
  if(names == []) then
    AbortBox("Image file from unknown instrument type")
  end
//
// Let the user choose the setup to use (or not)
  iSetup = x_choose(names,["Select instrument setup to use:"
                      "(To avoid this message use"
                      """Set Default"" in the ""Start"" menu)"])
  if(iSetup == 0) then
    AbortBox("No instrument setup selected")
  end
//
endfunction


function Match=MatchSetups(hostname,numxy,sdate)
//
// Create numeric equivalent of image date
// The image date (unlike setups.dates) is preceded by a time
  idate=msscanf(sdate,"%d:%d:%d %d/%d/%d")
  idate0=idate(4)+100*idate(5)+10000*idate(6)
//
// For each setup, score 1 for matching: instrum, numxy, dates
  Match=[]
  for i=1:size(InstrumSetups,1)
  Match(i)=0
// Does numxy match? If so, score at least 0.5
  if( and(InstrumSetups.numxy(i) == numxy) ) then
    Match(i)=0.5
// Does instrument match?
    if( InstrumSetups.host(i) == hostname ) then
// Calculate numeric equivalent of start and end dates
      idate=msscanf(InstrumSetups.dates(i)(1),"%d/%d/%d")
      idate1=idate(1)+100*idate(2)+10000*idate(3)
      idate=msscanf(InstrumSetups.dates(i)(2),"%d/%d/%d")
      idate2=idate(1)+100*idate(2)+10000*idate(3)
// Are dates compatible?
      if( (idate0 >= idate1) & (idate0 <= idate2) ) then
        Match(i)=1.0
      end
    end
  end
// Never return a match to the dummy setup
  Match(1)=0
end
//
endfunction


function n=ChooseInstrumDefault()
global InstrumDefault
//
// Make a vector of strings from the setup names
  names=[]
  for i=2:size(InstrumSetups,1)
    names(i-1)=InstrumSetups(i).setup
  end
//
// Let the user choose the setup to use, or not
  n = x_choose(names,["Select instrument setup to always use:"
                      "(Do not use this procedure unless"
                      "LaueG has suggested doing so!)"])
//
// Set (or clear) default setup
  if( (n == 0) & (InstrumDefault > 0) ) then
    InfoBox("The default instrument setup has been cleared")
    InstrumDefault=0
  else
    InstrumDefault=n+1
  end
//
endfunction
