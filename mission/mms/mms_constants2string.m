function s = mms_constants2string(listName)
%MMS_CONSTANTS2STRING  Convert cell array to string
%
%   s = mms_constants2string(listName)
%
% Example:
%   s = mms_constants2string('SDCProcs')

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

global MMS_CONST, if isempty(MMS_CONST), MMS_CONST = mms_constants(); end

if ~isfield(MMS_CONST,listName)
  errStr = ['No such MMS_CONST : ' listName];
  irf.log('critical',errStr), error(errStr)
end
l = MMS_CONST.(listName);
if ~iscell(l)
  errStr = ['MMS_CONST.' listName ' is not a cell array'];
  irf.log('critical',errStr), error(errStr)
end
s=['''' l{1} '''']; for iS=2:numel(l), s=sprintf('%s, ''%s''',s,l{iS}); end