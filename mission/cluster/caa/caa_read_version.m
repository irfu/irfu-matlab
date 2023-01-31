function [v,vdate] = caa_read_version(sp)
%CAA_READ_VERSION  read caa file version
%
% [v,vdate] = caa_read_version([sp])
%

% Copyright 2005 Yuri Khotyaintsev

if nargin<1, sp=pwd; end

old_pwd = pwd;
cd(sp);
if exist('./.version','file')
  [vers, d] = textread('./.version','%s %s',-1);
  if nargout==0
    for j=1:length(vers), disp([vers{j} ' : ' d{j}]), end
  elseif nargout==1
    v = str2num(vers{end});
  elseif nargout==2
    v = str2num(vers{end});
    vdate = d{end};
  else
    error('wrong number of input arguments')
  end
else
  if nargout==0
    disp('cannot find .version')
  elseif nargout==1
    v = [];
  elseif nargout==2
    v = [];
    vdate = '';
  else
    error('wrong number of input arguments')
  end
end
cd(old_pwd)
