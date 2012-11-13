function c_update
% LOCAL.C_UPDATE update file information in CAA directory
% 
% See also:
%	LOCAL.C_READ_AUX
%

curDir=pwd;
dirCaa='~/data/CAA';
cd(dirCaa);
clear aux
tmp=dir([dirCaa '/CL_SP_AUX']);
isub = [tmp(:).isdir]; %# returns logical vector
tmp(isub)=[];
f=vertcat(tmp.name);
fn=f;fn(:,end+1)='=';fn=fn';
%tt=textscan(fn(:),'%*11s%4d%2d%2d_%2d%2d%2d_%4d%2d%2d_%2d%2d%2d%*s','delimiter','=');
tt=textscan(fn(:),'%*11s%4f%2f%2f_%2f%2f%2f_%4f%2f%2f_%2f%2f%2f%*s','delimiter','=');

aux.filename=[repmat('CL_SP_AUX/',size(f,1),1) f];
aux.tstart=irf_time([tt{1} tt{2} tt{3} tt{4} tt{5} tt{6}],'vector2epoch');
aux.tend=irf_time([tt{7} tt{8} tt{9} tt{10} tt{11} tt{12}],'vector2epoch');

save('caa','aux','-append'); 
cd(curDir);

% $Id$
% $Revision$  $Date$

