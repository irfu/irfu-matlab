function ht=irf_pl_number_subplots(h,position,varargin)
%IRF_PL_NUMBER_SUBPLOTS alpha numbering of subplots
%  IRF_PL_NUMBER_SUBPLOTS adds a),b) to current figure
%
%  IRF_PL_NUMBER_SUBPLOTS(AX) number specified axis 
%
%  IRF_PL_NUMBER_SUBPLOTS(AX,position) put numbering at specified X,Y position
%  X,Y units are normalized, values between 0 and 1.
%
%  IRF_PL_NUMBER_SUBPLOTS(AX,position,'Property1',PropertyValue1,...)
%  sets the text property values of the numbering 
% 
%  See also IRF_LEGEND
% 
%  Example:
%     irf_pl_number_subplots(h,[0.02,0.97],'fontsize',14);

if nargin<=1, position=[0.01, 0.97];end
if nargin==0, h=irf_plot_get_subplot_handles; end

h=h(:);ht=h;
abc='a':'z';
for j=1:length(h)
  hleg=irf_legend(h(j),[abc(j) ')'],position,varargin{:});
  ht(j)=hleg;
end

if nargout==0, clear ht;end
