function obj=tt2000(inp)
%TT2000  Creates new EpochTT2000 object
%
%  OBJ = TT2000(SECS)
%
%  Creates new EpochTT2000 object OBJ from SECS time array in units of sec.
%  Zero corresponds to epoch TT2000=0.

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if ~isa(inp,'double'), 
  error('irf:tt2000:badInputs','input of type double expected')
end
inp = inp(:);
obj = EpochTT2000(int64(inp*1e9));