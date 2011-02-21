function ht=irf_pl_number_subplots(h,position,varargin)
% ht=irf_pl_number_subplots(h,position,text_property,text_value,...)
% 
% number subplots by putting in each subplot 'a', 'b', 'c', etc.
% INPUT 
%   h           handles to subplots
%   position    position of legend
% any text property/value pairs can be passed 

if nargin==1, position=[0.01, 0.97];end

h=h(:);ht=h;
abc='a':'z';
for j=1:length(h)
  hleg=irf_legend(h(j),abc(j),position,varargin{:});
  ht(j)=hleg;
end
