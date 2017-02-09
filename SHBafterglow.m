function [appmag,offset]=SHBafterglow(d,time,offset)

  % from Fig 2, Kann et al 'The Afterglows of Swift -era Gamma-Ray
  % Bursts. II. Short/Hard (Type I) vs. Long/Soft (Type II) Optical 
  % Afterglows', arXiv: 0804.1959

%   error(mynargchk(2,3,nargin,'struct'));

  if nargin==2
    offset=[];
  end

  if isempty(offset)
    offset=8*rand;
  end
  
  if (offset<0 || offset>8)
    warning('SHBafterglow: magnitude offset outside physical band');
  end

  days=time/86400;
  
  refdist=6634.3; %% Mpc at z=1

  refmag=Inf(size(days));
  idx=find(days>1e-4); %% observation before this time

  % Correction to obtain the magnitude in the cosmological frame
  % Cosmological frame time= oberver frame*(1+z)
  daysOF=days*2;

  refmag(idx)=23+offset+8/3*log10(daysOF(idx));  %% app magnitude for source at ref distance

  appmag=refmag+5*log10(d/refdist);

end

  
