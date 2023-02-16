function varargout = ds2nfu(varargin)
% DS2NFU  Convert data space units into normalized figure units.
%
% [Xf, Yf] = DS2NFU(X, Y) converts X,Y coordinates from
% data space to normalized figure units, using the current axes.  This is
% useful as input for ANNOTATION.
%
% POSf = DS2NFU(POS) converts 4-element position vector, POS from
% data space to normalized figure units, using the current axes.  The
% position vector has the form [Xo Yo Width Height], as defined here:
%
%      web(['jar:file:D:/Applications/MATLAB/R2006a/help/techdoc/' ...
%           'help.jar!/creating_plots/axes_pr4.html'], '-helpbrowser')
%
% [Xf, Yf] = DS2NFU(HAX, X, Y) converts X,Y coordinates from
% data space to normalized figure units, on specified axes HAX.
%
% POSf = DS2NFU(HAX, POS) converts 4-element position vector, POS from
% data space to normalized figure units, using the current axes.
%
% Ex.
%       % Create some data
% 		t = 0:.1:4*pi;
% 		s = sin(t);
%
%       % Add an annotation requiring (x,y) coordinate vectors
% 		plot(t,s);ylim([-1.2 1.2])
% 		xa = [1.6 2]*pi;
% 		ya = [0 0];
% 		[xaf,yaf] = ds2nfu(xa,ya);
% 		annotation('arrow',xaf,yaf)
%
%       % Add an annotation requiring a position vector
% 		pose = [4*pi/2 .9 pi .2];
% 		posef = ds2nfu(pose);
% 		annotation('ellipse',posef)
%
%       % Add annotations on a figure with multiple axes
% 		figure;
% 		hAx1 = subplot(211);
% 		plot(t,s);ylim([-1.2 1.2])
% 		hAx2 = subplot(212);
% 		plot(t,-s);ylim([-1.2 1.2])
% 		[xaf,yaf] = ds2nfu(hAx1,xa,ya);
% 		annotation('arrow',xaf,yaf)
% 		pose = [4*pi/2 -1.1 pi .2];
% 		posef = ds2nfu(hAx2,pose);
% 		annotation('ellipse',posef)

% Michelle Hirsch
% mhirsch@mathworks.com
% Copyright 2006-2014 The MathWorks, Inc

% Source: Matlab Central File Exchange
% https://se.mathworks.com/matlabcentral/fileexchange/10656-data-space-to-figure-units-conversion

% SPDX-License-Identifier: BSD-3-Clause
% LICENSE
% Copyright (c) 2006, The MathWorks, Inc.
% Copyright (c) 1997, Christophe COUVREUR
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
%    * Neither the name of the The MathWorks, Inc. nor the names
%      of its contributors may be used to endorse or promote products derived
%      from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% END OF LICENSE

%% Process inputs
narginchk(1, 3)

% Determine if axes handle is specified
if length(varargin{1})== 1 && ishandle(varargin{1}) && strcmp(get(varargin{1},'type'),'axes')
  hAx = varargin{1};
  varargin = varargin(2:end);
else
  hAx = gca;
end

errmsg = ['Invalid input.  Coordinates must be specified as 1 four-element \n' ...
  'position vector or 2 equal length (x,y) vectors.'];

% Proceed with remaining inputs
if length(varargin)==1	% Must be 4 elt POS vector
  pos = varargin{1};
  if length(pos) ~=4
    error(errmsg);
  end
else
  [x,y] = deal(varargin{:});
  if length(x) ~= length(y)
    error(errmsg)
  end
end


%% Get limits
axun = get(hAx,'Units');
set(hAx,'Units','normalized');
axpos = get(hAx,'Position');
axlim = axis(hAx);
axwidth = diff(axlim(1:2));
axheight = diff(axlim(3:4));


%% Transform data
if exist('x','var')
  varargout{1} = (x-axlim(1))*axpos(3)/axwidth + axpos(1);
  varargout{2} = (y-axlim(3))*axpos(4)/axheight + axpos(2);
else
  pos(1) = (pos(1)-axlim(1))/axwidth*axpos(3) + axpos(1);
  pos(2) = (pos(2)-axlim(3))/axheight*axpos(4) + axpos(2);
  pos(3) = pos(3)*axpos(3)/axwidth;
  pos(4) = pos(4)*axpos(4)/axheight;
  varargout{1} = pos;
end


%% Restore axes units
set(hAx,'Units',axun)

