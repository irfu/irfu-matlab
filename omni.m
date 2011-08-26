function f = omni( tint, parameter )
%OMNI download omni data 
% 
% DEVELOPMENT VERSION !!! 
%
% OMNI(tint,parameter) download parameters for specified time interval
% 
% parameters - string, paramters separated by comma. 
%               'B'     - magnetic field magnitude [nT]
%               'Bx'    - Bx GSE (the same as 'BxGSE')
%               'By'
%               'Bz'
%               'ByGSM' - By GSM 
%               'BzGSM' - Bz GSM
%               'n'     - proton density
%               'v'     - bulk speed
%               'dst'   - DST index
%               'f10.7' - F10.7 flux
%               'al'    - DST index
%   Detailed explanation goes here
%
% Examples:
%   ff= omni(tint,'b,bx,bygsm')
%   ff= omni(tint,'f10.7')

httpreq='http://omniweb.gsfc.nasa.gov/cgi/nx1.cgi?activity=retrieve&spacecraft=omni2&';
start_date=irf_time(tint(1),'yyyymmdd');
end_date=irf_time(tint(2),'yyyymmdd');

i=strfind(parameter,',');
iend=[i-1 length(parameter)];
istart=[1 i+1];
vars='';number_var=0;
for jj=1:length(istart)
    variable=parameter(istart(jj):iend(jj));
    switch lower(variable)
        case 'f10.7', var_number=50;
        otherwise, var_number=0;
    end
    if var_number>0, 
        vars=[vars '&vars=' num2str(var_number)]; 
        number_var=number_var+1;
    end
end

url=[httpreq 'start_date=' start_date '&end_date=' end_date vars];
disp(['url:' url]);
c=urlread(url);

cstart=strfind(c,'YEAR');
cend=strfind(c,'</pre>')-1;
fmt='%f %f %f';
for jj=1:number_var, fmt=[fmt ' %f']; end
cc=textscan(c(cstart:cend),fmt,'headerlines',1);

xx=double([cc{1} repmat(cc{1}.*0+1,1,2) repmat(cc{1}.*0,1,3)]);
f(:,1)=irf_time(xx)+(cc{2}-1)*3600*24+cc{3}*3600;
for jj=1:number_var, 
    f(:,jj+1)=cc{jj+3};
end
f(f==9999.9)=NaN;


