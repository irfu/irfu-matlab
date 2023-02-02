function [data_level,caa_vs,cl_id,DATA_VERSION,sp,st,dt]=caa_cefname2specs(s)
%CAA_CEFNAME2SPECS  extract details from a CEF file name
%
% [data_level,caa_vs,cl_id,DATA_VERSION,sp,st,dt]=caa_cefname2specs(cef_fname)
%
% See also CAA_EXPORT
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

cl_id = str2double(s(2));
data_level = str2double(s(12));

s = s(14:end);
if s(2)=='_', caa_vs = s(1); s = s(4:end);
elseif s(3)=='_', caa_vs = s(1:2); s = s(5:end);
else, caa_vs = s(1:3); s = s(6:end);
end

YY = s(1:4);
MM = s(5:6);
DD = s(7:8);
st_hh = s(10:11);
st_mm = s(12:13);
st_ss = s(14:15);
et_hh = s(17:18);
et_mm = s(19:20);
et_ss = s(21:22);
DATA_VERSION = s(25:26);

h3 = num2str( fix( str2double(st_hh)/3 )*3 );
if length(h3) == 1, h3 = [ '0' h3 ]; end

sp = [ '/data/caa/l1/' YY '/' YY MM DD '_' h3 '00/C' num2str(cl_id) '/' ...
  YY MM DD '_' st_hh st_mm ];

st = iso2epoch([ YY '-' MM '-' DD 'T' st_hh ':' st_mm ':' st_ss 'Z' ]);
et = iso2epoch([ YY '-' MM '-' DD 'T' et_hh ':' et_mm ':' et_ss 'Z' ]);
if str2double(et_hh)==0 && str2double(et_mm)==0 && str2double(et_ss)==0
  et = et + 86400;
end
dt = et - st;
