%% Copyright (C) 2010 Eric Chassande-Mottin, CNRS (France)
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, see
%% <http://www.gnu.org/licenses/>.

%% general routine to read ASCII data from a file and store into a cell struct array

function out=readcatalog(file,fields,delimiter,deblank)

  error(mynargchk(2,4,nargin,'struct'));
  
  if (nargin<3)
    delimiter=' ';
  endif

  if (nargin<4)
    deblank=false;
  endif

  %% determine local tmp dir
  tmpdir=getenv('TMPDIR');
  if isempty(tmpdir)
    tmpdir='/tmp';
  endif

  if exist(file,'file')~=2
    error('readcatalog: file %s not found',file);
  endif

  descr=dir(file);
  if descr.bytes == 0
    warning('readcatalog: file %s is empty',file);
    out=cell;
    return
  endif

  %% remove comments and blank lines
  tmpfile=sprintf('%s/tmp-%s.txt',tmpdir,randtag(9,'num'));
  if deblank
    deblank_cmd='| sed \'s/^ *//;s/ *$//;s/ \\{1,\\}/ /g\'';
  else
    deblank_cmd='';
  endif
  command=sprintf('grep -v \'^%\' %s %s > %s',file,deblank_cmd,tmpfile);
  system(command);

  try
    array=csv2cell(tmpfile,delimiter);
  catch
    [msg]=lasterr;
    printf('reading %s:\n',tmpfile);
    error(msg);
  end_try_catch

  if size(array,2)~=length(fields)
    error('file %s has %d columns while %d fields are requested',file,size(array,2),length(fields));
  endif
  out=cell2struct(array,fields,2);

  system(sprintf('rm -f %s',tmpfile));

  end
