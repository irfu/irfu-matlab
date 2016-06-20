function payload_mass_vs_number_of_sc
%% Defaults
payloadLim = [1.5 500]; % kg

nScMax   = 12;
xLim = [0.8 nScMax/0.8];
yLim = payloadLim; 

Cluster = struct('nSc',4,'mPayload',71,'name','Cluster');
CrossScaleI = struct('nSc',10,'mPayload',40,'name','Cross-Scale I');
CrossScaleII= struct('nSc',7,'mPayload',30,'name','Cross-Scale II');
EIDOSCOPE   = struct('nSc',1, 'mPayload',40,'name','EIDOSCOPE');
MMS = struct('nSc',4,'mPayload',160,'name','MMS');
THEMIS      = struct('nSc',5, 'mPayload',24,'name','THEMIS');
TorESA      = struct('nSc',1, 'mPayload',26,'name','Tor');

M4BigMum      = struct('nSc',1,  'mPayload',150,'name','THOR','nameLong','Big Mother');
M4MumSonI     = struct('nSc',1.3,'mPayload',125,'name','MS','nameLong','Mother Son');
M4MumSonII    = struct('nSc',1.3,'mPayload',19, 'name','MS','nameLong','Mother Son');
M4MumDaughterI= struct('nSc',1.5,'mPayload',105,'name','MD','nameLong','Mother Daughters');
M4MumDaughterII= struct('nSc',1.5,'mPayload',12,'name','MD','nameLong','Mother Daughters');
M4MotherMother= struct('nSc',2,  'mPayload',85, 'name','MM','nameLong','Mother Mother');
%M4MotherFather= struct('nSc',2,  'mPayload',79, 'name','MF','nameLong','Mother Father');
%M4Sisters     = struct('nSc',3,  'mPayload',50, 'name','S', 'nameLong','Sisters');

Lclass = struct('mLim',300*[1 1.4],'name','L'); % payload mass for single s/c
Mclass = struct('mLim',140*[1 1.4],'name','M'); % payload mass for single s/c
Sclass = struct('mLim', 10*[1 1.4],'name','S'); % payload mass for single s/c

ScProps = {'MarkerSize',30,'MarkerEdgeColor','b'};
ScTextProps = {'FontWeight','demi','FontSize',18,'color',[0.4 0.4 0.4]};
M4Props = {'MarkerSize',30,'MarkerEdgeColor','r'};
M4TextProps = {'FontWeight','demi','FontSize',18};
ClassLine = {'Color',[0.7 0.7 0.7],'linewidth',6,'linestyle','-'};
%% Initialize figure
set(0,'defaultLineLineWidth', 1.5);
fn=figure(61);
clf reset;
clear h;
set(fn,'color','white'); % white background for figures (default is grey)
set(gcf,'PaperUnits','centimeters')
xSize = 14; ySize = 2/3*xSize;
xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xLeft xSize sLeft ySize yTop
% additional good options
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');

h=axes('position',[0.1 0.15 0.6 0.7]); % [x y dx dy]

% y axis
set(h,'ylim',payloadLim,...
	'yscale','log',...
	'ytick',[3 6 12 25 50 100 200 ],...
	'yminortick','off',...
	'ygrid','on','yminorgrid','off');
ylabel(h,'payload mass [kg]');
% x-axis
set(h,'xlim',xLim,...
	'xscale','log',...
	'xtick',[1 2 3 4 5 7 10 12]);
xlabel(h,'Number of spacecraft');
hold(h,'on');

% plot class lines
plot_class(Sclass);
plot_class(Mclass);
plot_class(Lclass);

% plot spacecrat
plot_spacecraft(Cluster)
plot_spacecraft(CrossScaleI)
plot_spacecraft(CrossScaleII)
plot_spacecraft(EIDOSCOPE)
plot_spacecraft(MMS)
plot_spacecraft(THEMIS,'verticalalignment','top')
plot_spacecraft(TorESA)

plot_m4(M4BigMum)
%legend_m4(M4BigMum)

irf_legend(h,[datestr(now,31) '  ' mfilename],[0.99 1.01],...
	'fontsize',6,'interpreter','none','color',[0.9 0.9 0.9])
%% Functions nested
	function plot_spacecraft(Sc,varargin)
		plot(h,Sc.nSc,Sc.mPayload,'.',ScProps{:});
		text(Sc.nSc*1.05,Sc.mPayload*1.1,...
			Sc.name,'parent',h,ScTextProps{:},varargin{:});
	end
	function plot_m4(Sc,varargin)
		plot(h,Sc.nSc,Sc.mPayload,'.',M4Props{:});
		text(Sc.nSc*1.05,Sc.mPayload*1.1,...
			Sc.name,'parent',h,M4TextProps{:},varargin{:});
	end
	function legend_m4(varargin)
		textStruct = cell(1,nargin);
		for j=1:nargin
			Sc = varargin{j};
			textStruct{j} = [Sc.name ' - ' Sc.nameLong];
		end
		text(xLim(2),yLim(1),textStruct,...
			'verticalalignment','bottom');
	end
	function plot_class(Class)
		m = mean(Class.mLim);
		line([1 xLim(2)],[m m/xLim(2)],ClassLine{:});
		line([1 1],Class.mLim,'parent',h,ClassLine{:});
		text(0.95,mean(Class.mLim),Class.name,...
			'horizontalalignment','right','fontsize',20,'fontweight','bold');
	end
end

