function hout = irf_pl_matrix(x,y,F,dx,dy,Flab,xlab,ylab)
%IRF_PL_MATRIX plot matrix type data, for example energy/pitchangle distribution
%
% [h] = irf_pl_matrix(specrec)
% [h] = irf_pl_matrix(x,y,F,[dx],[dy],[Flabel],[xlabel],[ylabel])
%
% Input:
%    specrec - structure including spectra
%              specrec.x  - x vector
%              specrec.y - y vector
%              specrec.F - matrix to plot (size(x)xsize(y))
%              specrec.dx - vector of dt interval for every t point, can be single number if the same for all (can be omitted)
%              specrec.dy - vector of dF interval for every frequency f point, can be single number if the same for all (can be omitted)
%              specrec.Flabel - label of matrix to be plotted (can be omitted)
%              specrec.xlabel - label of xaxis to be plotted (can be omitted)
%              specrec.ylabel - label of yaxis to be plotted (can be omitted)
%
%
% See also CAA_SPECTROGRAM
%

narginchk(1,6)

if nargin==1
  specrec = x; x = [];
  if ~isfield(specrec,'dx'), specrec.dx=[];end
  if ~isfield(specrec,'dy'), specrec.dy=[];end
elseif nargin==2 % not defined
  help irf_pl_matrix;return
elseif nargin>=3 % caa_spectrogram(x,y,F,..)
  if nargin<4, dx=x(2)-x(1);     end % dx not defined
  if nargin<5, dy=y(2)-y(1);     end % dy not defined
  if nargin<6, Flab='';     end % not defined
  if nargin<7, xlab='';     end % not defined
  if nargin<8, ylab='';     end % not defined
  specrec.x = double(x);
  specrec.y = double(y);
  specrec.F = double(F);
  specrec.dx = double(dx);
  specrec.dy = double(dy);
  specrec.Flabel = Flab;
  specrec.xlabel = xlab;
  specrec.ylabel = ylab;
end



Fsize=size(F);
xx=[specrec.x(1)-specrec.dx(1)/2; specrec.x(:)+specrec.dx(:)/2];
xx=repmat(xx,1,Fsize(2)+1);
yy=[specrec.y(1)-specrec.dy(1)/2  specrec.y(:)'+specrec.dy(:)'/2];
yy=repmat(yy,Fsize(1)+1,1);
FF=zeros(size(xx));FF(1:end-1,1:end-1)=specrec.F;

FF(FF==0)=NaN; % zero counts put to NaN
pcolor(xx,yy,log10(FF))
%load caa/cmap.mat
%colormap(cmap)
shading('flat')
%	colorbar('vert')
%	set(gca,'TickDir','out','YScale','log')
set(gca,'TickDir','out')
xlabel(specrec.xlabel)
ylabel(specrec.ylabel)
