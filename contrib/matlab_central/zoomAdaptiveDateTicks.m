function zoomAdaptiveDateTicks(varargin)
% ZOOMADAPTIVEDATETICKS - Make date ticks adapt to zooming
%
% zoomAdaptiveDateTicks('on')
% Turns on the automatic adaptation of date ticks
% to user zooming for the current figure window
%
% zoomAdaptiveDateTicks('off')
% Turns off the automatic adaptation of date ticks
% to user zooming for the current figure window
%
% zoomAdaptiveDateTicks('demo')
% Opens a demo figure window to play with

% Source: Matlab Central File Exchange
% https://www.mathworks.com/matlabcentral/fileexchange/15342-zoom-adaptive-date-ticks

% SPDX-License-Identifier: BSD-3-Clause
% LICENSE
% Copyright (c) 2007, The MathWorks, Inc.
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
% * Neither the name of the The MathWorks, Inc. nor the names
%   of its contributors may be used to endorse or promote products derived
%   from this software without specific prior written permission.
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

if (nargin>0)
  switch varargin{1}
    case 'demo'
      % Create demo values
      dates = floor(now) - linspace(1169,0,15000)';
      values= randn(15000,1);
      % Show data with date ticks
      figure
      plot(dates,values)
      datetick('x')
      zoomAdaptiveDateTicks('on')
    case 'on'
      % Define a post zoom callback
      set(zoom(gcf),'ActionPostCallback', @adaptiveDateTicks);
    case 'off'
      % Delete the post zoom callback
      set(zoom(gcf),'ActionPostCallback', '');
    otherwise
      figure(gcf)
  end
end


function adaptiveDateTicks(figureHandle,eventObjectHandle)
% Resetting x axis to automatic tick mark generation
set(eventObjectHandle.Axes,'XTickMode','auto')
% using automaticallly generate date ticks
datetick(eventObjectHandle.Axes,'x','keeplimits')
