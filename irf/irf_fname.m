function out=irf_fname(tint,fmt)
%IRF_FNAME construct a filename string from time
%
% fname = irf_fname(tint, [FORMAT])
% Input:
% st - ISDAT epoch
%
% Output:
% fname - string formatted accordint to FORMAT:
%	0: YYYYMMDD_hhmm (default)
%	1: YYMMDDhhmmss
%	2: YYYYMMDD_hhmmss_hhmmss (CAA)
%	3: YYYYMMDD    (new CAA daily format)
%	4: YYMMDD_hhmmss_mmm
%	5: YYYYMMDD_hhmmss_YYYYMMDD_hhmmss (CAA/MAARBLE)
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(1,2)

if nargin < 2, fmt = 0; end
if isa(tint,'GenericTimeArray')
  tint = [tint.start.epochUnix tint.stop.epochUnix];
end

if fmt==2 && length(tint)<=1, error('ST must have two elements for FORMAT=2'), end

d = fromepoch(tint(1));
s{1} = num2str(d(1));

for k=2:6
  s{k} = num2str(fix(d(k)));
  if d(k)<10, s{k} = ['0' s{k}]; end
end

switch fmt
  case 0
    out = [s{1} s{2} s{3} '_' s{4} s{5}];
  case 1
    out = [s{1}(3:4) s{2} s{3} s{4} s{5} s{6}(1:2)];
  case 2
    d = fromepoch(tint(end));
    for k=4:6
      se{k-3} = num2str(fix(d(k)));
      if d(k)<10, se{k-3} = ['0' se{k-3}]; end
    end
    out = [s{1} s{2} s{3} '_' s{4} s{5} s{6} '_' se{1} se{2} se{3}];
  case 3
    out = [s{1} s{2} s{3}];
  case 4
    s{6}=num2str(d(6),'%06.3f');
    out = [s{1}(3:4) s{2} s{3} '_' s{4} s{5} s{6}(1:2) '_' s{6}(4:6)];
  case 5
    d = fromepoch(tint(end));
    se{1} = num2str(d(1));
    for k=2:6
      se{k} = num2str(fix(d(k)));
      if d(k)<10, se{k} = ['0' se{k}]; end
    end
    out = [s{1} s{2} s{3} '_' s{4} s{5} s{6} '_'...
      se{1} se{2} se{3} '_' se{4} se{5} se{6}];
  otherwise
    error('unknown format')
end
