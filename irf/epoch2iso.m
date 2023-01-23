function out = epoch2iso(t,fmt)
%EPOCH2ISO   Convert ISDAT epoch to ISO time string
%
% ISDAT epoch is the number of seconds since 1-Jan-1970 and
% ISO time string has a format yyyy-mm-ddTHH:MM:ss.wwwwwwZ
% as described in the CEF data file syntax
% http://www.space-plasma.qmul.ac.uk/csds/welcome.html
%
% t_iso_s = epoch2iso(t_epoch,[fmt])
%   fmt - 0: long (default), 1: short
%
%   See also EPOCH2DATE, ISO2EPOCH, TOEPOCH, FROMEPOCH
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin<2, fmt = 0; end

switch fmt % set rounding precision for different formats
  case 0
    dt_res = 5e-7;
  case 1
    dt_res = 5e-4;
end

if length(t)<5 || any(diff(t)<0) || (t(end)-t(1)) > length(t)*100
  d = fromepoch(t);
  switch fmt
    case 0
      out = num2str(d,'%04d-%02d-%02dT%02d:%02d:%09.6fZ');
    case 1
      out = num2str(d,'%04d-%02d-%02dT%02d:%02d:%06.3fZ');
  end
  
  ii = find(out(:,18)=='6'); % in case there has been rounding leading to 60.000 seconds
  if any(ii)
    out(ii,:) = epoch2iso(t(ii)+dt_res,fmt);
  end
else
  % This approach is faster for data with many samples per minute, as we run
  % from epoch only once per minute
  % This approach is only applicable to monotonically increasing time
  out = char(zeros(length(t),27-fmt*3));
  out(:,[5 8]) = '-';
  out(:,11) = 'T';
  out(:,[14 17]) = ':';
  out(:,end) = 'Z';
  
  tss = fromepoch(t(1));
  tee = fromepoch(t(end));
  
  mins = toepoch([tss(1:5) 0]):60:toepoch([tee(1:5) 0]);
  d = fromepoch(mins);
  
  s1 = {'', '', '', '',''}; sl=[0 4; 5 2; 8 2; 11 2; 14 2];
  j_start = 0;
  
  for j=1:5
    if d(1,j)==d(end,j)
      ss = num2str(d(1,j),'%02d');
      for jj=1:sl(j,2), out(:,sl(j,1)+jj) = ss(jj); end
    else
      j_start = j;
      for jj=j:5, s1(jj) = {num2str(d(:,jj),'%02d')}; end
      break
    end
  end
  
  for j=1:length(mins)
    if j==length(mins), ii = find(t>=mins(j));
    else, ii = find(t>=mins(j) & t<mins(j+1));
    end
    if isempty(ii), continue, end
    if j_start
      for kk=j_start:5
        for jj=1:sl(kk,2), out(ii,sl(kk,1)+jj) = s1{kk}(j,jj); end
      end
    end
    switch fmt
      case 0
        out(ii,18:end-1) = num2str(t(ii)-mins(j),'%09.6f');
      case 1
        out(ii,18:end-1) = num2str(t(ii)-mins(j),'%06.3f');
    end
  end
  ii = find(out(:,18)=='6'); % in case there has been rounding leading to 60.000 seconds
  if any(ii)
    out(ii,:) = epoch2iso(t(ii)+dt_res,fmt);
  end
end

