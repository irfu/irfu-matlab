function [Ddsi,Damp] = c_efw_dsi_off(t,cl_id,Ps)
%C_EFW_DSI_OFF  get EFW offsets
%
% [Ddsi,Damp] = c_efw_dsi_off(t,[cl_id,Ps])
%
% [Ddsi,Damp] = c_efw_dsi_off(t,cl_id,'magnetosphere')
%    Get the magnetospheric offsets
%
% Ddsi is complex: Dx = real(Ddsi), Dy = imag(Ddsi)
%
% See also CAA_COROF_DSI

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(1,3)

SC_POT_LIM = -8;  % Above this we apply SW/SH correction, below - MS
TAV = 300; % Averaging window for SC potential
Damp = 1.1*ones(1,4);

% i is for Ey +i send curve down

% Table of SW/SH offsets
% t>=toepoch([2018 11 19 21 00 0]), Ddsi = [ x x x x ];
% elseif t>=toepoch([2018 06 23 06 00 0]), Ddsi = [ x x x x ]; % Force MS offset for months when don't enter SW.

if t>=toepoch([2022 01 01 00 00 0]), Ddsi = [ 00 2.43  00 0.22 ];
elseif t>=toepoch([2021 01 01 00 00 0]), Ddsi = [ 00 2.25  00 0.20 ];
elseif t>=toepoch([2020 01 01 00 00 0]), Ddsi = [ 00 2.09  00 0.09 ];
elseif t>=toepoch([2019 01 01 00 00 0]), Ddsi = [ -0.15 1.93  00 -0.12 ];
elseif t>=toepoch([2017 09 01 00 00 0]), Ddsi = [ -0.18 1.84  00 -0.09 ];
elseif t>=toepoch([2017 03 24 12 00 0]), Ddsi = [ -0.14 1.73  00 -0.13 ];
elseif t>=toepoch([2017 03 23 12 00 0]), Ddsi = [ -1.14 1.73  00 -0.13 ];
elseif t>=toepoch([2017 01 01 00 00 0]), Ddsi = [ -0.14 1.73  00 -0.13 ];
elseif t>=toepoch([2016 01 01 00 00 0]), Ddsi = [ -0.36 1.95  00 -1.15 ];
elseif t>=toepoch([2015 01 01 00 00 0]), Ddsi = [ 0.34  2.89  00  0.31 ];
elseif t>=toepoch([2014 10 16 00 00 0]), Ddsi = [ -0.11  1.67  1.22  0.37 ];
elseif t>=toepoch([2014 10 16 00 00 0]), Ddsi = [ -0.11  3.67  1.22  0.37 ];
elseif t>=toepoch([2014 10 11 00 00 0]), Ddsi = [ -0.11  2.67  1.22  0.37 ];
elseif t>=toepoch([2014 08 01 00 00 0]), Ddsi = [ -0.11  3.67  1.22  0.37 ];
elseif t>=toepoch([2014 01 01 00 00 0]), Ddsi = [ -0.11  2.67  1.22  0.37 ];
elseif t>=toepoch([2013 01 01 00 00 0]), Ddsi = [ -0.13  2.19  0.97  0.43 ];
elseif t>=toepoch([2012 03 03 00 00 0]), Ddsi = [ -0.28  1.28  0.7  0.11 ];
elseif t>=toepoch([2012 01 06 12 00 0]), Ddsi = [ -0.28  1.28  1.0  0.11 ];
elseif t>=toepoch([2011 11 01 00 00 0]), Ddsi = [ -0.28   0.1  1.1  0.11 ];
elseif t>=toepoch([2011 05 09 15 00 0]), Ddsi = [ 0.03  0.65+0.25i  2.04  0.51 ];
elseif t>=toepoch([2011 05 02 21 00 0]), Ddsi = [ 0.03  3.0+0.65+0.25i  2.04  0.51 ];
elseif t>=toepoch([2011 04 30 06 40 0]), Ddsi = [ 0.03  4.0+0.65+0.25i  2.04  0.51 ];
elseif t>=toepoch([2011 03 01 16 50 0]), Ddsi = [ 0.03  0.65+0.25i  2.04  0.51 ];
elseif t>=toepoch([2011 02 28 17 00 0]), Ddsi = [ 0.03  0.65+0.25i  2.04+0.6  0.51 ];
elseif t>=toepoch([2010 10 30 00 00 0]), Ddsi = [ 0.03  0.65+0.25i  2.04  0.51 ];
elseif t>=toepoch([2010 07 12 00 00 0]), Ddsi = [ 0.4   1.6  1.22 0.82 ]; % Force MS offset for months when don't enter SW.
elseif t>=toepoch([2010 07 08 18 00 0]), Ddsi = [ -0.27 0.8  1.66 0.25 ]; % These orbits need the variable offset
elseif t>=toepoch([2010 07 01 00 00 0]), Ddsi = [ 0.4  1.58  1.21 0.82 ]; % Force MS offsets.
elseif t>=toepoch([2010 06 01 00 00 0]), Ddsi = [ -0.27 0.8  1.66 0.25 ];
elseif t>=toepoch([2009 11 13 00 00 0]), Ddsi = [ -0.34 0.67 1.50-0.25i 0.17 ]; % Back to variable offsets, Adjust Ey on C3
elseif t>=toepoch([2009 07 01 00 00 0]), Ddsi = [ 0.46  1.31 1.23  0.64 ]; % Force MS offset for months when don't enter SW.
elseif t>=toepoch([2008 12 01 00 00 0]), Ddsi = [-0.18  0.40 1.33 -0.07 ]; % Back to variable offsets
elseif t>=toepoch([2008 07 01 00 00 0]), Ddsi = [ 0.59  1.32 1.38  0.69 ]; % Force MS offset for months when don't enter SW.
elseif t>=toepoch([2008 01 01 00 00 0]), Ddsi = [-0.18  0.40 1.33 -0.07 ];
elseif t>=toepoch([2007 11 01 01 01 0]), Ddsi = [ 0.20  0.76 1.77  0.28 ];
elseif t>=toepoch([2007 08 05 01 01 0]), Ddsi = [ 0.72  1.46 1.53  0.87 ]; % Force MS offset for months when don't enter SW.
elseif t>=toepoch([2007 08 01 00 00 0]), Ddsi = [0.72 1.46+1.4 1.53 0.87]; % Force MS offset for months when don't enter SW.
elseif t>=toepoch([2007 07 21 18 45 0]), Ddsi = [-0.08     .46+1.4 1.65 .13 ]; % problem with guard settings on C2
elseif t>=toepoch([2007 02 01 00 00 0]), Ddsi = [-0.08     .46 1.65 .13 ]; % very approximate due to high-speed solar wind streams
elseif t>=toepoch([2007 01 01 00 00 0]), Ddsi = [-0.08     .46 1.95 .13 ]; % very approximate due to high-speed solar wind streams
elseif t>=toepoch([2006 10 01 00 00 0]), Ddsi = [ .26      .79 2.0  .59 ];
elseif t>=toepoch([2006 10 20 00 00 0]), Ddsi = [ .26      .79 1.7  .59 ];
elseif t>=toepoch([2006 07 20 01 01 0]), Ddsi = [ 0.84  1.59 1.69  1.02 ]; % Force MS offset for months when don't enter SW.
elseif t>=toepoch([2006 07 01 00 00 0]), Ddsi = [ .26      .79 1.7  .59 ];
elseif t>=toepoch([2006 02 01 00 00 0]), Ddsi = [ .46     1.04 2.0  .59 ];
elseif t>=toepoch([2006 01 01 00 00 0]), Ddsi = [ .46     1.13 2.0  .59 ];
elseif t>=toepoch([2005 07 01 00 00 0]), Ddsi = [ .31      .60 .52  .64 ]; % Big jump to 2006-01-01.
elseif t>=toepoch([2005 03 01 00 00 0]), Ddsi = [ .35+0.2i .78 .51  .62 ];
elseif t>=toepoch([2004 11 01 00 00 0]), Ddsi = [ .35      .78 .51  .62 ];
elseif t>=toepoch([2004 07 01 00 00 0]), Ddsi = [ .35      .78 .41  .42 ];
elseif t>=toepoch([2004 05 01 00 00 0]), Ddsi = [ .65      .95 .41  .42 ]; % manually checked
elseif t>=toepoch([2004 01 01 00 00 0]), Ddsi = [ .23      .72 .50  .42 ];
elseif t>=toepoch([2003 07 02 23 30 0]), Ddsi = [ .15      .53 .47  .71 ];
elseif t>=toepoch([2003 07 01 12 40 0]), Ddsi = [1.42      .53 .47  .71 ]; % HXONLY on C1 in the magnetosphere
elseif t>=toepoch([2002 12 02 00 00 0]), Ddsi = [ .15      .53 .47  .71 ];
elseif t>=toepoch([2002 05 02 00 00 0]), Ddsi = [ .33      .69 .73  .33 ];
elseif t>=toepoch([2002 01 01 00 00 0]), Ddsi = [ .33      .69 .73  .92 ];
elseif t>=toepoch([2001 07 01 00 00 0]), Ddsi = [ .47      .82 .89 1.04 ];
elseif t>=toepoch([2001 06 01 00 00 0]), Ddsi = [ .31      .60 .52  .64 ];
elseif t>=toepoch([2001 05 25 00 00 0]), Ddsi = [ .31     1.35 .52 1.55 ];
elseif t>=toepoch([2001 04 25 00 00 0]), Ddsi = [ .31      .60 .52  .64 ];
elseif t>=toepoch([2001 03 01 00 00 0]), Ddsi = [ .69     1.36 .68  .34 ];
elseif t>=toepoch([2001 02 02 15 00 0]), Ddsi = [ .55      .77 .44  .1  ];
elseif t>=toepoch([2001 02 02 00 00 0]), Ddsi = [ .48      .77 .44 1.11 ]; % Special puck/guard ?
elseif t>=toepoch([2001 02 01 00 00 0]), Ddsi = [ .55      .8  .4   .1  ];
else
  Ddsi = [ .55 .8 .4 .1 ];
end

if nargin == 1, return, end

DdsiSW = Ddsi;
Ddsi = Ddsi(cl_id);
Damp = Damp(cl_id);

if nargin == 2 || isempty(Ps), return, end

flagAlwaysMagnetosphere = 0;
if isnumeric(Ps)
  ndata = ceil((Ps(end,1) - Ps(1,1))/TAV);
  ta = Ps(1,1) + (1:ndata)*TAV - TAV/2; ta = ta';
  Psr = irf_resamp( Ps( ~isnan(Ps(:,2)) ,:), ta, 'window',TAV);
  if isempty(Psr), return, end

  ii = find(Psr(:,2) < SC_POT_LIM);
  if isempty(ii), return, end
elseif ischar(Ps)
  if strcmpi(Ps, 'magnetosphere')
    flagAlwaysMagnetosphere = 1;
  else
    error('Unrecognazed value of region')
  end
end

% Table of MS offsets
if t>=toepoch([2022 01 01 00 0 0]), Ddsi = [ 00 3.65 00 0.57 ];
elseif t>=toepoch([2021 01 01 00 0 0]), Ddsi = [ 00 2.99 00 0.35 ];
elseif t>=toepoch([2020 01 01 00 0 0]), Ddsi = [ 00 2.59 00 0.29 ];
elseif t>=toepoch([2019 01 01 00 0 0]), Ddsi = [ 00 2.44 00 -0.03 ];
elseif t>=toepoch([2018 01 01 00 0 0]), Ddsi = [ 0.05 2.44 00 0.17 ];
elseif t>=toepoch([2017 01 01 00 0 0]), Ddsi = [ 0.26 2.44 00 0.25 ];
elseif t>=toepoch([2016 01 01 00 0 0]), Ddsi = [ -0.17 2.45 00 -0.51 ];
elseif t>=toepoch([2015 01 01 00 0 0]), Ddsi = [ -0.02 2.89 00 -0.44 ]; % curves fixed
elseif t>=toepoch([2014 10 16 00 0 0]), Ddsi = [ 0.18  3.7  0.71  0 ]; % Only C3 has a good curve
elseif t>=toepoch([2014 10 11 00 0 0]), Ddsi = [ 0.18  2.7  0.71  0 ]; % Only C3 has a good curve
elseif t>=toepoch([2014 08 01 00 0 0]), Ddsi = [ 0.18  3.7  0.71  0 ]; % Only C3 has a good curve
elseif t>=toepoch([2014 07 01 00 0 0]), Ddsi = [ 0.18  2.7  0.71  0 ]; % Only C3 has a good curve
elseif t>=toepoch([2014 01 01 00 0 0]), Ddsi = [ 0.18  2.7  0.71  0.4 ]; % Only C3 has a good curve
elseif t>=toepoch([2013 01 01 00 0 0]), Ddsi = [ 0.15  2.7  0.94  0.2 ]; % Only C3 has a good curve
elseif t>=toepoch([2012 05 01 00 0 0]), Ddsi = [ 0.1  2.4  0.50  0.43 ];
elseif t>=toepoch([2012 01 01 00 0 0]), Ddsi = [ 0.1  2.4  0.47  0.45 ];
elseif t>=toepoch([2011 11 01 00 0 0]), Ddsi = [ 0.49  0.78  1.18  0.84 ];
elseif t>=toepoch([2011 06 01 00 0 0]), Ddsi = [ 0.49  0.78-2.78  1.18  0.84 ];
elseif t>=toepoch([2011 01 01 00 0 0]), Ddsi = [ 0.49  0.78  1.18  0.84 ]; % C2 strange curve limited data
elseif t>=toepoch([2010 10 30 00 0 0]), Ddsi = [ 0.03  0.78  2.04  0.51 ]; % force to SW offsets
elseif t>=toepoch([2010 01 01 00 0 0]), Ddsi = [ 0.4   1.58 1.21  0.82 ];
elseif t>=toepoch([2009 01 01 00 0 0]), Ddsi = [ 0.46  1.31 1.23  0.64 ];
elseif t>=toepoch([2008 01 01 00 0 0]), Ddsi = [ 0.59  1.32 1.38  0.69 ];
elseif t>=toepoch([2007 11 01 00 0 0]), Ddsi = [ 0.72  1.46 1.69  0.87 ];
elseif t>=toepoch([2007 08 05 01 1 0]), Ddsi = [ 0.72  1.46 1.53  0.87 ];
elseif t>=toepoch([2007 07 21 18 45 0]), Ddsi = [ 0.72  1.46+1.4 1.53 0.87 ]; % problem with guard settings on C2
elseif t>=toepoch([2007 01 01 00 0 0]), Ddsi = [ 0.72  1.46 1.53  0.87 ];
elseif t>=toepoch([2006 01 01 00 0 0]), Ddsi = [ 0.84  1.59 1.69  1.02 ];
elseif t>=toepoch([2005 01 01 00 0 0]), Ddsi = [ 1.26  2.34 1.81  1.37 ];
elseif t>=toepoch([2004 01 01 00 0 0]), Ddsi = [ 1.35  2.06 1.45  1.15 ];
elseif t>=toepoch([2003 01 01 00 0 0]), Ddsi = [ 1.42  2.18 1.64  1.43 ];
elseif t>=toepoch([2002 12 02 00 0 0]), Ddsi = [ 1.42  1.98 1.64  2.00 ];
elseif t>=toepoch([2002 05 02 00 0 0]), Ddsi = [ 1.33  1.98 1.66  1.30 ];
elseif t>=toepoch([2002 01 01 00 0 0]), Ddsi = [ 1.33  1.98 1.66  2.00 ];
elseif t>=toepoch([2001 06 01 00 0 0]), Ddsi = [ 1.26  1.74 1.54  1.06 ];
elseif t>=toepoch([2001 01 01 00 0 0]), Ddsi = [ 1.21  1.92 1.25  1.02 ];
else
  Ddsi = DdsiSW;
end

% SC pot is all the time below SC_POT_LIM
if flagAlwaysMagnetosphere || ~any(Psr(:,2) >= SC_POT_LIM)
  Ddsi = Ddsi(cl_id); return
end

DdsiMS = Ddsi;

Dd = Psr; clear Psr
Dd(:,1) = Dd(:,1) - TAV/2; % offset is set at the start of the interval
Dd(:,2) = DdsiSW(cl_id);
Dd(ii,2) = DdsiMS(cl_id);

% Remove repeating points
d = [1; diff(Dd(:,2))]; d(d~=0) = 1;
Ddsi = Dd(d==1,:);

