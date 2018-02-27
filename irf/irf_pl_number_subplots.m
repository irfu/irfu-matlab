function ht=irf_pl_number_subplots(h,position,varargin)
%IRF_PL_NUMBER_SUBPLOTS alpha numbering of subplots
%  IRF_PL_NUMBER_SUBPLOTS adds a),b) labeling to current figure
%
%  IRF_PL_NUMBER_SUBPLOTS(AX) number specified axis
%
%  IRF_PL_NUMBER_SUBPLOTS(AX,position) put numbering at specified X,Y position
%  X,Y units are normalized, values between 0 and 1.
%
%  IRF_PL_NUMBER_SUBPLOTS(AX,position,'Property1',PropertyValue1,...)
%  sets the text property values of the numbering
%
%  IRF_PL_NUMBER_SUBPLOTS(AX,position,'num','(?)',...) specify template for
%  labels. ? changed to a-z. In this example labels would be (a),(b),..
%
%  IRF_PL_NUMBER_SUBPLOTS(AX,position,'firstletter','b',...) specify which
%  letter is the first
%
%  See also IRF_LEGEND
%
%  Example:
%     irf_pl_number_subplots(h,[0.02,0.97],'fontsize',14);

if nargin<=1, position=[0.01, 0.97];end
if nargin==0, h=irf_plot_get_subplot_handles; end
%% Defaults
num='?)'; % template for numbering, question mark changed to a-z
firstLetter = 'a';
%% Check input options
inp_flag=ones(size(varargin));
for jj=1:numel(varargin)
    switch(lower(varargin{jj}))
        case 'num'
            num = varargin{jj+1};
            inp_flag(jj:jj+1)=0;
        case 'firstletter'
            firstLetter = varargin{jj+1};
            inp_flag(jj:jj+1)=0;
    end
end
args=varargin(inp_flag==1);

%% add numbering
h=h(:);ht=h;
abc=firstLetter:'z';
for j=1:length(h)
    hleg=irf_legend(h(j),irf_ssub(num,abc(j)),position,args{:});
    ht(j)=hleg;
end

if nargout==0, clear ht;end
