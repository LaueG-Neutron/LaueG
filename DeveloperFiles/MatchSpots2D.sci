nspots=30000
ix_width=30
iy_width=30
numx=4000
numy=2000

xy=[rand(nspots,1)*numx,rand(nspots,1)*numy];

ipairs=MatchSpots2D(xy,ix_width,iy_width);
size(ipairs)

//
ibrute=[];
for i1=1:nspots
  for i2=i1+1:nspots
    if( (abs(xy(i1,1)-xy(i2,1)) < ix_width) & (abs(xy(i1,2)-xy(i2,2)) < iy_width) ) then
      ibrute($+1,:)=[i1,i2];
    end
  end
end
ibrute(1,:)=[];
size(ibrute)


//
ibrute=[];
for i1=1:nspots
  i2=[i1+1:nspots];
  i=find( (abs(xy(i1,1)-xy(i2,1)) < ix_width) & (abs(xy(i1,2)-xy(i2,2)) < iy_width) );
  i12=[i1*ones(i'),i2(i)'];
  ibrute=[ibrute;i12];
end
ibrute(1,:)=[];
size(ibrute)





[xy(23,:);xy(92,:)]
[xy(26,:);xy(91,:)]
[xy(47,:);xy(100,:)]
[xy(86,:);xy(87,:)]

xy2=[xy(17,:);xy(83,:)]


ixy=round( [xy2(:,1)/ix_width + 0.0 , xy2(:,2)/iy_width + 0.5] )







function ipairs=MatchSpots2D(xy,ix_width,iy_width)
//
// Match spots using double sized X bins
  ipairs1=MatchSpots2D_half(xy,2*ix_width,iy_width,0.5)
// Repeat with X bins offset by 1/2 a bin
  ipairs2=MatchSpots2D_half(xy,2*ix_width,iy_width,0.0)
//
// Combine and remove any repeats
  ipairs=unique([ipairs1;ipairs2],'r')
//
// Remove any pairs with X move than 1 width apart
  i=find(abs(xy(ipairs(:,1),1)-xy(ipairs(:,2),1)) > ix_width)
  ipairs(i,:)=[]
//
endfunction


function ipairs=MatchSpots2D_half(xy,ix_width,iy_width,xoffset)
//
// Allocate spots to X & Y bins
  ixy=round( [xy(:,1)/ix_width + xoffset , xy(:,2)/iy_width + 0.5] )
//
// Get maximum indices for X & Y bins
  num_ix=max(ixy(:,1))
  num_iy=max(ixy(:,2))
//
// Do tagged sort in X bins, then in Y bins within each X bin
  isort=1000*ixy(:,1) + ixy(:,2)
  [v,k]=gsort(isort,'r','i')
//
// Create list of last indices to last X bins in sorted order
// List displaced by 1 to include index to non-existing bin "0"
  ix_last=zeros(num_ix+1,1)
  for ix=1:nspots
    ix_last( ixy(k(ix),1) +1 )=ix;
  end
//
// Fill in indices for empty X bins
  for i=1:num_ix
    if(ix_last(i+1) == 0) then
      ix_last(i+1)=ix_last(i)
    end
  end
//
// Loop over X bins, finding start and end indices for each X bin
  ipairs=[0,0];
  for ix=1:num_ix
    ix_start=ix_last(ix)+1;
    ix_end=ix_last(ix+1);
//
// Loop over spot #1 for all spots in this X bin
    for ix1=ix_start:ix_end
      iix1=k(ix1);
// Get X & Y for spot #1
      x1=xy(iix1,1);
      y1=xy(iix1,2);
// Loop over spot #2 for all spots later than spot #1
      for ix2=ix1+1:ix_end
        iix2=k(ix2);
// If spot #2 is later than the next Y bin to spot #1, terminate loop
        if(ixy(iix2,2) > ixy(iix1,2) + 1) then break; end
// Get X & Y for spot #2
        x2=xy(iix2,1);
        y2=xy(iix2,2);
// Ignore spot #2 is more than 1 width in X or Y from spot #1
        dx=(x2-x1);
        dy=(y2-y1);
        if( (abs(dx) > ix_width) | (dy > iy_width) ) then continue; end
// Output spot pair and distance between them
        ipairs($+1,:)=[iix1,iix2]
      end
//
    end
//
  end
  ipairs(1,:)=[]
//
endfunction

