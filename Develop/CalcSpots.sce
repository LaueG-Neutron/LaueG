global ImageInfo LatticeInfo OrientInfo DisplayData DisplayMode

// ===== Convert between pixels, vectors, and HKLs =====
//function [Hvec]=PixWav2Hvec(pixwav)
//function [pixwav]=Hvec2PixWav(Hvec)
//function [Hvec]=Pix2Hvec(xy_pix)
//function [xy_pix]=Hvec2Pix(Hvec)
//function [hkl]=PixWav2HKL(pixwav)
//function [pixwav]=HKL2PixWav(hkls)
//function [hkls]=Pix2HKL(pix)
//function [pix]=HKL2Pix(hkls)
//function [wav]=HKL2Wav(hkls)
//function [Hvec]=HKL2Hvec(hkls)
//function [hkls]=Hvec2HKL(Hvec)
// ===== Convert between d-spacing, cell & UB =====
//function dspace=HKL2DSpacing(hkl)
//function B=Cell2Bmatrix(d)
//function B=Cell2BmatrixBuerger(d)
//function UB=CalcUBfromLauegen(phi_x,phi_y,phi_z)
//function ucell=UB2UnitCell(ub)
//function r=Direct2Recip(d)
// ===== Symmetry/Centering routine =====
//function bOK=AllowedHKL(hkl, Linfo)
// ===== Miscellaneous routines =====
//function [hkl1,hkl2,angle]=CalcSplitRotation(ub1,ub2)
//function xy_unmatch=CalcUnmatchedSpots()
//function nguess=GuessCalcSpotsNumber(linfo)
//function nguess=GuessCalcSpotsNumber2(linfo,d_min,wav_min,wav_max)
//function match=MatchXY(xy1,xy2, dxy)
//function mx=CalcXRot(rot)     // rot is in degrees
//function mx=CalcYRot(rot)     // rot is in degrees
//function mx=CalcZRot(rot)     // rot is in degrees
// ===== Generate HKLs for twins =====
//function hkls=TwinGenHKLs(d_min,wav_min,twin_dist)
// ===== Unused single vector versions without CYCLOPS code =====
//function [Hvec]=PixWav2Hvec_ALT(pixwav)
//function [pixwav]=Hvec2PixWav_ALT(Hvec)


//=====================================================================
//
//Instrument axis: Y Phi rotation axis (upwards)
//                 Z Incident beam direction (horizontal plane projection)
//                 X Cross-product Y by Z (left if facing upstream)
//OrientInfo.beamvert     Angle of incident beam from horizontal
//OrientInfo.xtaloff      Offset of sample to instrument origin at Phi=0
//OrientInfo.drumoff      Offset of IP drum to instrument origin
//OrientInfo.drumrad      Radius of IP surface inside the drum
//OrientInfo.pixsize      X & Y pixel sizes in microns
//OrientInfo.pixskew      X/Y pixel skewness factor
//OrientInfo.pixcen       Pixel position where two-theta=0
//LatticeInfo.cell        Sample lattice cell dimensions
//LatticeInfo.ub          Sample lattice ub matrix
//ImageInfo.phi           Spindle rotation angle (CW from above)
//
//=====================================================================


// ===== Convert between pixels, vectors, and HKLs =====

function [Hvec]=PixWav2Hvec(pixwav)
// From pixel X,Y and wavelength to reciprocal space 
//
try
//
// Return null if input if null
  if(pixwav == []) then Hvec=[]; return; end
//
// Non-vectorized stuff
  xtallab=CalcYRot(-ImageInfo.phi) * OrientInfo.xtaloff'
  offset=xtallab - OrientInfo.drumoff'
  S0=[0; sind(OrientInfo.beamvert); cosd(OrientInfo.beamvert)]
// Spot position in pixels relative to the pixel center
  dx=pixwav(:,1) - OrientInfo.pixcen(1)
  dy=pixwav(:,2) - OrientInfo.pixcen(2)
// If CYCLOPS, get the plate number and subtract the offset from x
  if( ImageInfo.instrum == "CYCLOPS") then
    iPlate=1+int( pixwav(:,1)/480 )
    dx=dx -(iPlate-5)*480
  end
// Convert spot position from pixels to mm
  dx=dx * OrientInfo.pixsize(1)
  dy=dy * OrientInfo.pixsize(2)
// Apply skewness correction
  dx=dx +OrientInfo.pixskew*dy
// Apply empirical correction to x
  dx=dx -xtallab(1)
// Diffracted beam vector from the sample to the pixel position
  if( ImageInfo.instrum == "CYCLOPS") then
// Rotate offset from laboratory to detector coords
    yrot=45*(iPlate-5)
    offset=CalcYRot(yrot)*offset
// Calculate the vector S1 from the sample to the spot
    S1=[ dx - offset(1) , dy - offset(2) ]
    S1(:,3)=OrientInfo.drumrad - offset(3)
// Rotate S1 from detector to laboratory coords (NOT VECTORISED)
    S1=CalcYRot(-yrot)*S1
  else
// Calculate the vector S1 from the sample to the spot
    angle=dx/OrientInfo.drumrad
    S1x=OrientInfo.drumrad*sin(angle) - offset(1)
    S1y=dy                            - offset(2)
    S1z=OrientInfo.drumrad*cos(angle) - offset(3)
  end
// Calculate diffraction vector from wavelength
  S1_sizes=sqrt(S1x.^2 + S1y.^2 + S1z.^2)
  Hvec(:,1)=( S1x ./ S1_sizes - S0(1) ) ./ pixwav(:,3)
  Hvec(:,2)=( S1y ./ S1_sizes - S0(2) ) ./ pixwav(:,3)
  Hvec(:,3)=( S1z ./ S1_sizes - S0(3) ) ./ pixwav(:,3)
//
catch
  Hvec=[]
end
//
endfunction


function [pixwav]=Hvec2PixWav(Hvec)
// From reciprocal space scattering vector to pixel X,Y and wavelength
//
try
//
// Return null if input is null
  if(Hvec == []) then pixwav=[]; return; end
//
// Non-vectorized stuff
  xtallab=CalcYRot(-ImageInfo.phi) * OrientInfo.xtaloff'
  offset=xtallab-OrientInfo.drumoff'
  S0=[0; sind(OrientInfo.beamvert); cosd(OrientInfo.beamvert)]
// Calculate wavelengths and the diffracted beam unit-vectors S1
  wav=-2*(Hvec*S0) ./ sum(Hvec.^2,'c')
  S1x=S0(1) + Hvec(:,1) .* wav
  S1y=S0(2) + Hvec(:,2) .* wav
  S1z=S0(3) + Hvec(:,3) .* wav
// Calculate (x,y) the spot position, in mm relative to the detector center
  if( ImageInfo.instrum == "CYCLOPS") then
// Calculate the plate the beam will hit
    iplate=5+round( atand(S1x,S1z) / 45 )
    if(iplate == 9) then iplate=1; end
// Rotate S1 & offset from laboratory to detector coords
    yrot=45*(iplate-5)
    S1=CalcYRot(yrot)*[S1x;S1y;S1z]
    offset=CalcYRot(yrot)*offset
// Calculate the position of the spot on the plate relative to the center
    x=S1(1,:)./S1(3,:)*(OrientInfo.drumrad-offset(3)) + offset(1)
    y=S1(2,:)./S1(3,:)*(OrientInfo.drumrad-offset(3)) + offset(2)
  else
// Intersection of S1 with IP surface relative to sample origin
    fac1=S1x.^2+S1z.^2
    fac2=offset(1)*S1x + offset(3)*S1z
    fac3=offset(1)*S1z - offset(3)*S1x
    dist =( sqrt( fac1 .* OrientInfo.drumrad.^2 - fac3.^2 ) - fac2 ) ./ fac1
    x=offset(1) + dist .* S1x
    y=offset(2) + dist .* S1y
    z=offset(3) + dist .* S1z
// Convert (x,z) to spot position in mm around the IP cylinder
    x=OrientInfo.drumrad*atan(x,z)
  end
//
// Apply empirical correction to x
  x = x - xtallab(1)
// Add skewness correction to x
  x = x - OrientInfo.pixskew*y
// Convert to pixel positions and add wavelength
  pixwav(:,1)=OrientInfo.pixcen(1) + x / OrientInfo.pixsize(1)
  pixwav(:,2)=OrientInfo.pixcen(2) + y / OrientInfo.pixsize(2)
  pixwav(:,3)=wav
// If CYCLOPS, add the plate offset to X
  if( ImageInfo.instrum == "CYCLOPS") then
    pixwav(:,1)=pixwav(:,1) +(iplate-5)*(ImageInfo.numxy(1)/8)
  end
//
catch
  pixwav=[]
end
//
endfunction


function [Hvec]=Pix2Hvec(xy_pix)
// Return reciprocal space (unit) scattering vector for pixel X,Y
//
// Return null if input if null
  if(xy_pix == []) then Hvec=[]; return; end
//
  pixwav=xy_pix
  pixwav(:,3)=1
  Hvec=PixWav2Hvec(pixwav)
  Hnorm=sqrt(Hvec(:,1).^2 + Hvec(:,2).^2 + Hvec(:,3).^2)
  Hvec=Hvec./[Hnorm,Hnorm,Hnorm]
endfunction


function [xy_pix]=Hvec2Pix(Hvec)
//
// Return null if input if null
  if(Hvec == []) then xy_pix=[]; return; end
//
  pixwav=Hvec2PixWav(Hvec)
  xy_pix(:,1)=pixwav(:,1)
  xy_pix(:,2)=pixwav(:,2)
endfunction


function [hkl]=PixWav2HKL(pixwav)
// Return hkl for pixel X,Y and wavelength
//
// Return null if input if null
  if(pixwav == []) then hkl=[]; return; end
//
  rot_phi=CalcYRot(-ImageInfo.phi)
  trans_mx=inv( rot_phi * LatticeInfo(1).ub )'
  Hvec=PixWav2Hvec(pixwav)
  hkl=Hvec*trans_mx
endfunction


function [pixwav]=HKL2PixWav(hkls)
//
// Return null if input if null
  if(hkls == []) then pixwav=[]; return; end
//
// Return [] if no UB matrix
  if( isempty(LatticeInfo(1).ub) ) then
    pixwav=[]
    return
  end
// Rotate the UB by Phi angle
  rot_ub=CalcYRot(-ImageInfo.phi) * LatticeInfo(1).ub
// Now all the vectorized stuff
  Hvec=rot_ub*hkls'
  pixwav=Hvec2PixWav(Hvec')
endfunction


function [hkls]=Pix2HKL(pix)
// This works with pix=[]
  pixwav=[pix,ones(pix(:,1))]
  hkls=PixWav2HKL(pixwav)
endfunction


function [pix]=HKL2Pix(hkls)
// This works with hkls=[]
  pixwav=HKL2PixWav(hkls)
  pix=pixwav(:,1:2)
endfunction


function [wav]=HKL2Wav(hkls)
// This works with hkls=[]
  pixwav=HKL2PixWav(hkls)
  wav=pixwav(:,3)
endfunction


function [Hvec]=HKL2Hvec(hkls)
//
// Return null if input if null
  if(hkls == []) then Hvec=[]; return; end
//
// Return [] if no UB matrix
  if( isempty(LatticeInfo(1).ub) ) then
    Hvec=[]
    return
  end
// Rotate the UB by Phi angle
  rot_ub=CalcYRot(-ImageInfo.phi) * LatticeInfo(1).ub
// Now all the vectorized stuff
  Hvec=(rot_ub*hkls')'
endfunction


function [hkls]=Hvec2HKL(Hvec)
//
// Return null if input if null
  if(Hvec == []) then hkls=[]; return; end
//
// Return [] if no UB matrix
  if( isempty(LatticeInfo(1).ub) ) then
    hkls=[]
    return
  end
// Rotate the UB by Phi angle
  rot_ub=CalcYRot(-ImageInfo.phi) * LatticeInfo(1).ub
// Now all the vectorized stuff
  hkls=( inv(rot_ub) * Hvec' )'
endfunction


// ===== Convert between d-spacing, cell & UB =====

function dspace=HKL2DSpacing(hkl)
  hvec=LatticeInfo(1).ub*hkl';
  dspace=1.0/sqrt( hvec'*hvec );
endfunction


function B=Cell2Bmatrix(d)
// convert from direct to reciprocal cell
  r=Direct2Recip(d)
// Busing & Levy's equation for B().
  B=[r(1)  r(2)*cosd(r(6))  r(3)*cosd(r(5))
     0     r(2)*sind(r(6))  r(3)*sind(r(5))*cosd(d(4))
     0     0                1/d(3)                    ]
endfunction


function B=Cell2BmatrixBuerger(d)
// convert from direct to reciprocal cell
  r=Direct2Recip(d)
// Buerger's Xray Crystallography, p348
  cos_psi=( cosd(r(6)) - cosd(r(5)) * cosd(r(4)) )/sind(r(5))
  cos_rho=sqrt( 1 - cosd(r(6))^2 - cosd(r(4))^2 - cosd(r(5))^2 ...
                + 2*cosd(r(6))*cosd(r(4))*cosd(r(5)) ) / sind(r(5))
//
  B=[  r(1)*sind(r(5))   r(2)*cos_psi      0
       0                 r(2)*cos_rho      0
       r(1)*cosd(r(5))   r(2)*cosd(d(4))   r(3)  ]
endfunction


function UB=CalcUBfromLauegen(phi_x,phi_y,phi_z)
  B=Cell2Bmatrix(LatticeInfo(1).cell)
  U = CalcZRot(phi_z) * CalcYRot(-phi_y) * CalcXRot(phi_x)
  UB = [1 0 0; 0 0 1; 0 -1 0] * U * B
endfunction


function ucell=UB2UnitCell(ub)
  abc=1/(ub' * ub)
  a=sqrt(abc(1,1))
  b=sqrt(abc(2,2))
  c=sqrt(abc(3,3))
  alp=acosd( abc(2,3)/(b*c) )
  bet=acosd( abc(1,3)/(a*c) )
  gam=acosd( abc(1,2)/(a*b) )
  ucell=[a,b,c,alp,bet,gam]
endfunction


function r=Direct2Recip(d)
// Formulae from Int. Tables. II, p106
  s=(d(4)+d(5)+d(6))/2.0
  temp=sind(s)*sind(s-d(4))*sind(s-d(5))*sind(s-d(6))
  vol=2*d(1)*d(2)*d(3)*sqrt(max(0, temp ))
  r(1:3)=[d(2)*d(3)*sind(d(4)) d(1)*d(3)*sind(d(5)) d(1)*d(2)*sind(d(6))]/vol
  temp=[(cosd(d(5))*cosd(d(6))-cosd(d(4)))/(sind(d(5))*sind(d(6)))
        (cosd(d(4))*cosd(d(6))-cosd(d(5)))/(sind(d(4))*sind(d(6)))
        (cosd(d(4))*cosd(d(5))-cosd(d(6)))/(sind(d(4))*sind(d(5)))]
  r(4:6)=acosd(min(1, temp' ))
endfunction


// ===== Symmetry/Centering routine =====

function bOK=AllowedHKL(hkl, Linfo)
// LaueG only uses centering info, not spacegroups
  bOK=%t
  if    (Linfo.icen == 1) then
// A centering: (K+L=2N)
    bOK=(modulo(hkl(2)+hkl(3),2) == 0)
  elseif(Linfo.icen == 2) then
// B centering: (H+L=2N)
    bOK=(modulo(hkl(1)+hkl(3),2) == 0)
  elseif(Linfo.icen == 3) then
// C centering: (H+K=2N)
    bOK=(modulo(hkl(1)+hkl(2),2) == 0)
  elseif(Linfo.icen == 4) then
// I centering: (H+K+L=2N)
    bOK=(modulo(hkl(1)+hkl(2)+hkl(3),2) == 0)
  elseif(Linfo.icen == 5) then
// F centering: (H+K=2N), (K+L=2N)
    bOK= (modulo(hkl(1)+hkl(2),2) == 0) & (modulo(hkl(2)+hkl(3),2) == 0)
  elseif(Linfo.icen == 6) then
// R centering: (-H+K+L=3N)
    bOK=(modulo(-hkl(1)+hkl(2)+hkl(3),3) == 0)
  end
//
endfunction


// ===== Miscellaneous routines =====

function [hkl1,hkl2,angle]=CalcSplitRotation(ub1,ub2)
//
// Calculates rotation axis and angle for U matrices extracted from UBs.
// Rotation axis returned as a HKL in both lattices.
  cell1=UB2UnitCell(ub1)
  b1=Cell2Bmatrix(cell1)
  u1=ub1/b1
  cell2=UB2UnitCell(ub2)
  b2=Cell2Bmatrix(cell2)
  u2=ub2/b2
//
// Massage rotation in recip space to a true rotation
  rot=u2/u1
  for iter=1:4
    rot=( rot + inv(rot') )/2
  end
//
// Get vector of rotation axis
  vec(1)=rot(3,2)-rot(2,3)
  vec(2)=rot(1,3)-rot(3,1)
  vec(3)=rot(2,1)-rot(1,2)
  vec=vec/norm(vec)
//
// Make largest indices positive
  [v,k]=max(abs(vec))
  vec=vec*sign(vec(k))
//
// Make a perpendicular vector roughly along 001
  vperp=cross(vec, cross(vec,[0.0;0.0;-1.0]))
  vperp=vperp/norm(vperp)
//
// Make second perpendicular vector
  vperp2=cross(vec,vperp)
//
// Calc rotation of vperp in direction of vperp2
  vec2=rot*vperp
  angle=acosd(vperp2'*vec2) -90.0
//
// Ensure angle > 0 by inverting axis
  if(angle < 0) then
    angle=-angle
    vec=-vec
  end
//
// Convert rotation axis to HKL for pair 1 & 2
  hkl1=inv(ub1)*vec
  hkl2=inv(ub2)*vec
//
// Scale HKLs to make largest index = +/-1
  hkl1=hkl1/max(abs(hkl1))
  hkl2=hkl2/max(abs(hkl2))
//
endfunction


function xy_unmatch=CalcUnmatchedSpots()
// Generate calculated spots
  hkllist=RunGenHKLs(LatticeInfo(1), ...
                     DisplaySet.wav_min,DisplaySet.wav_max, ...
                     DisplaySet.d_min,DisplaySet.d_max, ...
                     DisplayMode.modul_show)
// Match obs & calc spots within 5 pixels apart
  xy_calc=hkllist(:,5:6)
  xy_obs=DisplayData.found_spots(:,1:2)
  match=MatchXY(xy_obs,xy_calc,5.0)
// Return x,y of any unmatched "obs" spots


  nobs=size(xy_obs,1)
  unmatch=setdiff([1:nobs],match(:,1)')
  xy_unmatch=xy_obs(unmatch,1:2)
endfunction


function nguess=GuessCalcSpotsNumber(linfo)
// Guess how many spots will be generated for LatticeInfo = linfo
//
  d_min=DisplaySet.d_min
  wav_min=DisplaySet.wav_min
  wav_max=DisplaySet.wav_max
  nguess=GuessCalcSpotsNumber2(linfo,d_min,wav_min,wav_max)
//
endfunction


function nguess=GuessCalcSpotsNumber2(linfo,d_min,wav_min,wav_max)
// Guess how many spots will be generated for LatticeInfo = linfo
  Vcell=1/abs( det(linfo.ub) )
  div=[1,2,2,2,2,4,3]
//
  nguess=0.658*Vcell*(  sqrt( (d_min/wav_min) / (0.06*wav_min^6+d_min^6) ) ...
                      - sqrt( (d_min/wav_max) / (0.06*wav_max^6+d_min^6) )  )
  nguess=nguess/div(linfo.icen+1)
//
endfunction


function match=MatchXY(xy1,xy2, dxy)
// Return match() where coords xy1(match(i,1),1:2) match
// xy2(match(i,2),1:2) within a distance of dxy pixels
//
  num1=size(xy1,1)
  num2=size(xy2,1)
//
// No matches if either list is empty
  if(min(num1,num2) == 0) then
    match=[]
    return
  end
//
// Do tagged sorts on the x values of xy1 & xy2
  [val tag1]=gsort(xy1(:,1),"r","i")
  [val tag2]=gsort(xy2(:,1),"r","i")
//
// Loop through each sorted element of xy1()
  match=[]
  i2lo=1
  for i1=1:num1
    x1=xy1(tag1(i1),1)
    y1=xy1(tag1(i1),2)
// Increment i2lo until x2(i2lo) >= x1(i1)-dxy
    for i2=i2lo:num2
      if(xy2(tag2(i2),1) >= x1-dxy) then
        break
      end
    end
    i2lo=i2;
// Find sorted xy2() with |x2-x1| < dxy, then if
// xy2() within dxy pixels of xy1() add to match()
    for i2=i2lo:num2
      x2=xy2(tag2(i2),1)
      if(x2 > x1+dxy) then
        break
      end
      y2=xy2(tag2(i2),2)
      if((x2-x1).^2 + (y2-y1).^2 < dxy.^2) then
        match($+1,1:2)=[tag1(i1) tag2(i2)]
      end
    end
//
  end
//
endfunction


function mx=CalcXRot(rot)     // rot is in degrees
      mx=[  1.0    0.0        0.0
            0.0  cosd(rot) -sind(rot)
            0.0  sind(rot)  cosd(rot)]
endfunction


function mx=CalcYRot(rot)     // rot is in degrees
      mx=[cosd(rot) 0.0 -sind(rot)
            0.0     1.0    0.0
          sind(rot) 0.0  cosd(rot)]
endfunction


function mx=CalcZRot(rot)     // rot is in degrees
      mx=[cosd(rot) -sind(rot) 0.0
          sind(rot)  cosd(rot) 0.0
            0.0        0.0     1.0]
endfunction


// ===== Generate HKLs for twins (used by argonne_boxes) =====

function hkls=TwinGenHKLs(d_min,wav_min,twin_dist)
//
global LatticeInfo
//
// Generate hkl lists of (hkl, d-space, xy, wav, mult, imod)
// List is ordered in unpaired spots of twin 1, then twin 2,
// then interleaved pairs of twinned spots from twins 1 & 2
// imod =1,2 for unpaired twins 1,2
// imod =11,12 for paired twins 1,2
//
// Calculate matrix for HKL equivalence of overlapping spots
  rat12=LatticeInfo(2).ub \ LatticeInfo(3).ub
  rat12=round(rat12*12)/12
//
// Generate hkl lists for Twin 1 & 2
  hkls1=RunGenHKLs(LatticeInfo(2),wav_min,999.99,d_min,999.99,%f);
  hkls2=RunGenHKLs(LatticeInfo(3),wav_min,999.99,d_min,999.99,%f);
//
// Copy twin-related pairs of spots to hkls11[] & hkls12[]
  [val,idx1,idx2]=intersect(hkls1(:,1:3),round( hkls2(:,1:3) * rat12' ),"r")
  hkls11=hkls1(idx1,:)
  hkls12=hkls2(idx2,:)
// Remove paired spots to get unpaired ones
  hkls1(idx1,:)=[]
  hkls2(idx2,:)=[]
//
// Find any paired twins further apart than twin_dist pixels
  dist=sqrt(sum( (hkls11(:,5:6)-hkls12(:,5:6)).^2, "c"))
  imove=find(dist>twin_dist)
// Move far spots from hkls11 & hkls12 to hkls1 & hkls2
  hkls1=[hkls1;hkls11(imove,:)]
  hkls2=[hkls2;hkls12(imove,:)]
  hkls11(imove,:)=[]
  hkls12(imove,:)=[]
//
// Find unpaired spots closer than twin_dist pixels
  idxs=[0,0]
  for i1=1:size(hkls1,1)
    dsq=(hkls1(i1,5)-hkls2(:,5)).^2 + (hkls1(i1,6)-hkls2(:,6)).^2
    [val,k]=min(dsq)
    if(val < twin_dist^2) then
      idxs($+1,:)=[i1,k]
    end
  end
  idxs(1,:)=[]
// Move close spots from hkls1 & hkls2 to hkls11 & hkls12
  if(idxs ~= []) then
    hkls11=[hkls11;hkls1(idxs(:,1),:)]
    hkls12=[hkls12;hkls2(idxs(:,2),:)]
    hkls1(idxs(:,1),:)=[]
    hkls2(idxs(:,2),:)=[]
  end
//
// Add imod values to hkl lists
  if(hkls1 ~= []) then hkls1(:,9)=1; end
  if(hkls2 ~= []) then hkls2(:,9)=2; end
  if(hkls11 ~= []) then hkls11(:,9)=11; end
  if(hkls12 ~= []) then hkls12(:,9)=12; end
//
// Concatenate unpaired spots 1 & 2 with interleaved pairs of spots
  hkls=[hkls1; hkls2; matrix( [hkls11,hkls12]', 12,-1)' ]
//
endfunction


// ===== Unused single vector versions without CYCLOPS code =====

function [Hvec]=PixWav2Hvec_ALT(pixwav)
// From pixel X,Y and wavelength to reciprocal space
// Only works for single pixwav[] vector, no CYCLOPS code
//
// Spot position in pixels relative to the pixel center
  dxy=pixwav(1:2) - OrientInfo.pixcen
// Wavelength of spot
  wav=pixwav(3)
// Convert spot position from pixels to mm, with skewness correction
  dxy=dxy .* OrientInfo.pixsize
  dxy(1)=dxy(1) +OrientInfo.pixskew*dxy(2)
// Vector from instrument origin to spot position
  dxy_rad=dxy/OrientInfo.drumrad
  S1=OrientInfo.drumrad * [sin(dxy_rad(1)), dxy_rad(2), cos(dxy_rad(1))]
// Total sample offset due to xtaloff & drumoff
  offset=CalcYRot(-ImageInfo.phi) * OrientInfo.xtaloff' - OrientInfo.drumoff'
// Diffracted beam unit-vector from sample to spot position
  S1=S1 - offset'
  S1=S1/norm(S1)
// Calculate diffraction vector from incident beam unit-vector and wavelength
  S0=[0, sind(OrientInfo.beamvert), cosd(OrientInfo.beamvert)]
  Hvec=( S1 - S0 ) / wav
//
endfunction


function [pixwav]=Hvec2PixWav_ALT(Hvec)
// From reciprocal space scattering vector to pixel X,Y and wavelength
// Only works for single Hvec[] vector, no CYCLOPS code
//
// Incident beam unit-vector S0
  S0=[0, sind(OrientInfo.beamvert), cosd(OrientInfo.beamvert)]
// Calculate wavelength and the diffracted beam unit-vector S1
  wav=-2*(Hvec*S0') / (Hvec*Hvec')      // Braggs' Law for vectors
  S1=S0 + Hvec * wav                    // Must be unit-vector
// Total sample offset due to xtaloff & drumoff
  offset=CalcYRot(-ImageInfo.phi) * OrientInfo.xtaloff' - OrientInfo.drumoff'
// Distance from sample to intersection of S1 with IP surface
  fac1=S1(1)^2 + S1(3)^2
  fac2=offset(1)*S1(1) + offset(3)*S1(3)
  fac3=offset(1)*S1(3) - offset(3)*S1(1)
  dist =( sqrt( fac1 * OrientInfo.drumrad^2 - fac3^2 ) - fac2 ) / fac1
// Spot position relative to instrument origin
  xyz=dist * S1 + offset'
// Convert xyz to spot position in mm around the IP cylinder
  xpix=OrientInfo.drumrad*atan(xyz(1),xyz(3))
  ypix=xyz(2)
// Convert from mm to pixel positions
  xpix=xpix - OrientInfo.pixskew*ypix
  pixwav=OrientInfo.pixcen + [xpix,ypix] ./ OrientInfo.pixsize
// Add wavelength
  pixwav(3)=wav
//
endfunction
