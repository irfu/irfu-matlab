function [y]=irf_gse2gsm(x,flag_gse2gsm)
%function [y]=irf_gse2gsm(x,flag);
%
% x - [t,Xx,Xy,Xz]  column vector, time in isdat epoch
% y - [t,Yx,Yy,Yz]  column vector
% flag - if flag not given do GSE->GSM, if flag=-1 do GSM->GSE any other flag values do GSE->GSM
%
% Uses: IRF.GEOCENTRIC_COORDINATE_TRANSFORMATION
%

%function [y,tilt]=irf_gse2gsm(x,flag); NOT IMPLEMENTED, NEEDED?

if nargin==0, help irf_gse2gsm;return;end
conv='gse>gsm';
if nargin==2 && flag_gse2gsm==-1 
	conv='gsm>gse';
end
if isempty(x) 
	irf_log('fcal','empty input variable'); 
	y=x; 
	return; 
end % if empty input, empty output

y=irf.geocentric_coordinate_transformation(x,conv);
