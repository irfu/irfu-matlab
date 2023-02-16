function errorbar_tick(h,w,xtype)
%ERRORBAR_TICK Adjust the width of errorbars
%   ERRORBAR_TICK(H) adjust the width of error bars with handle H.
%      Error bars width is given as a ratio of X axis length (1/80).
%   ERRORBAR_TICK(H,W) adjust the width of error bars with handle H.
%      The input W is given as a ratio of X axis length (1/W). The result
%      is independent of the x-axis units. A ratio between 20 and 80 is usually fine.
%   ERRORBAR_TICK(H,W,'UNITS') adjust the width of error bars with handle H.
%      The input W is given in the units of the current x-axis.
%
%   See also ERRORBAR
%

% Author: Arnaud Laurent
% Creation : Jan 29th 2009
% MATLAB version: R2007a
%
% Notes: This function was created from a post on the french forum :
% http://www.developpez.net/forums/f148/environnements-developpement/matlab/
% Author : Jerome Briot (Dut)
%   http://www.mathworks.com/matlabcentral/newsreader/author/94805
%   http://www.developpez.net/forums/u125006/dut/
% It was further modified by Arnaud Laurent and Jerome Briot.

% SPDX-License-Identifier: BSD-2-Clause
% LICENSE
% Copyright (c) 2009, Arnaud Laurent
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% * Redistributions of source code must retain the above copyright
%   notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
%   notice, this list of conditions and the following disclaimer in
%   the documentation and/or other materials provided with the distribution
%
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%   POSSIBILITY OF SUCH DAMAGE.
% END OF LICENSE

% Check numbers of arguments
narginchk(1,3)

% Check for the use of V6 flag ( even if it is depreciated ;) )
flagtype = get(h,'type');

% Check number of arguments and provide missing values
if nargin==1
  w = 80;
end

if nargin<3
  xtype = 'ratio';
end

% Calculate width of error bars
if ~strcmpi(xtype,'units')
  dx = diff(get(gca,'XLim'));	% Retrieve x limits from current axis
  w = dx/w;                   % Errorbar width
end

% Plot error bars
if strcmpi(flagtype,'hggroup') % ERRORBAR(...)
  
  hh=get(h,'children');		% Retrieve info from errorbar plot
  x = get(hh(2),'xdata');		% Get xdata from errorbar plot
  
  x(4:9:end) = x(1:9:end)-w/2;	% Change xdata with respect to ratio
  x(7:9:end) = x(1:9:end)-w/2;
  x(5:9:end) = x(1:9:end)+w/2;
  x(8:9:end) = x(1:9:end)+w/2;
  
  set(hh(2),'xdata',x(:))	% Change error bars on the figure
  
else  % ERRORBAR('V6',...)
  
  x = get(h(1),'xdata');		% Get xdata from errorbar plot
  
  x(4:9:end) = x(1:9:end)-w/2;	% Change xdata with respect to the chosen ratio
  x(7:9:end) = x(1:9:end)-w/2;
  x(5:9:end) = x(1:9:end)+w/2;
  x(8:9:end) = x(1:9:end)+w/2;
  
  set(h(1),'xdata',x(:))	% Change error bars on the figure
  
end
