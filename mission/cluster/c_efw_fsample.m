function sf = c_efw_fsample(data,mode,cl_id)
%C_EFW_FSAMPLE  guess sampling frequency
%
% fs = c_efw_fsample(data,[mode, cl_id])
% returns sampling frequency in Hz
%
% data - first column is interpreted as time. Other columns are ignored.
% mode - restrict to particular type of data :
%        'any' - default 0.25/5/25/450/4500/9000 Hz
%        'lx'  - LX data 5 Hz
%        'hx'  - HX data 25/450 Hz
%        'ib'  - internal burst data 450/4500/9000/18000 Hz
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

data = data(:,1);

if length(data)<=1
  error('cannot compute sampling frequency from less than two points')
end
if nargin <3, fdwp = 900.005;
else
  switch cl_id
    case 1
      fdwp = 900.020; % Hz
    case 2
      fdwp = 900.013;
    case 3
      fdwp = 900.005;
    case 4
      fdwp = 900.006;
    otherwise
      error('Incorrect CL_ID')
  end
end

if nargin <2, mode='any'; end
if ~(strcmp(mode,'any') || strcmp(mode,'lx') || strcmp(mode,'hx') || strcmp(mode,'ib'))
  error('bad value for MODE')
end

sf = guess_fsample((length(data) - 1)/(data(end) - data(1)),mode,fdwp);
if sf, return
else
  sf = guess_fsample(1/(data(2) - data(1)),mode,fdwp);
  if sf, return
  elseif length(data)>2
    sf = guess_fsample(1/(data(3) - data(2)),mode,fdwp);
    if sf, return
    else
      sf = guess_fsample(1/(data(end) - data(end-1)),mode,fdwp);
      if sf, return, end
    end
  end
end

if ~sf, irf_log('proc','cannot guess sampling frequency'), end

%% help function
function sf = guess_fsample(f,mode,fdwp)
K_PLUS = 1.1;
K_MINUS = .9;

if f<K_PLUS*5 && f>K_MINUS*5 && (strcmp(mode,'any') || strcmp(mode,'lx'))
  sf = 5;     % LX
elseif f<K_PLUS*fdwp/36 && f>K_MINUS*fdwp/36 && (strcmp(mode,'any') || strcmp(mode,'hx'))
  sf = fdwp/36;    % NM
elseif f<K_PLUS*fdwp/2 && f>K_MINUS*fdwp/2 && (strcmp(mode,'any') || strcmp(mode,'hx') || strcmp(mode,'ib'))
  sf = fdwp/2;   % BM1/IB
elseif f<K_PLUS*2250 && f>K_MINUS*2250 && (strcmp(mode,'any') || strcmp(mode,'ib'))
  sf = 2250;  % IB
elseif f<K_PLUS*4500 && f>K_MINUS*4500 && (strcmp(mode,'any') || strcmp(mode,'ib'))
  sf = 4500;  % IB
elseif f<K_PLUS*9000 && f>K_MINUS*9000 && (strcmp(mode,'any') || strcmp(mode,'ib'))
  sf = 9000;  % IB
elseif f<K_PLUS*18000 && f>K_MINUS*18000 && (strcmp(mode,'any') || strcmp(mode,'ib'))
  sf = 18000;  % IB
elseif f<K_PLUS*.25 && f>K_MINUS*.25 && strcmp(mode,'any'), sf = .25;   % SPIN
else, sf = 0;
end
