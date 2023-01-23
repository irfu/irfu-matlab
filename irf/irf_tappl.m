function [z]=irf_tappl(x,s)
%IRF_TAPPL   Apply expression to data
%
% [z]=irf_tappl(x,s)
% x is time vector, first column time  other data
% s is what expression to apply to the data part,
%
% Example:
%    y=irf_tappl(x,'*2/1e3');
%
Units=irf_units;

z=x;
if isa(x,'TSeries')
  z.data=eval(['x.data',s]);
else
  z(:,2:end)=eval(['x(:,2:end)',s]);
end