function kext = onera_desp_lib_kext(kext)
%***************************************************************************************************
% Copyright 2006, T.P. O'Brien
%
% This file is part of IRBEM-LIB.
%
%    IRBEM-LIB is free software: you can redistribute it and/or modify
%    it under the terms of the GNU Lesser General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    IRBEM-LIB is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public License
%    along with IRBEM-LIB.  If not, see <http://www.gnu.org/licenses/>.
%
% Last modified by A. C. Kellerman - added TS07d functionality - March 2016
%***************************************************************************************************
%
% function kext = onera_desp_lib_kext(kext)
% converts string representation of kext to number
% leaves numeric input unchanged
% kext: long integer to select external magnetic field 
% 
% 
% '',IGRF        : 0   = no external field 
% MF             : 1   = Mead & Fairfield [1975] (uses 0?Kp?9 - Valid for rGEO?17. Re) 
% T87, T87short  : 2   = Tsyganengo short [1987] (uses 0?Kp?9 - Valid for rGEO?30. Re) 
% T87long        : 3   = Tsyganengo long [1987] (uses 0?Kp?9- Valid for rGEO?70. Re) 
% T89, T89c      : 4   = Tsyganenko [1989c] (uses 0?Kp?9- Valid for rGEO?70. Re) 
% OPQ            : 5   = Olson & Pfitzer quiet [1977] (default - Valid for rGEO?15. Re) 
% OPD            : 6   = Olson & Pfitzer dynamic [1988] (uses 5.?dens?50., 300.?velo?500., -100.?Dst?20. - Valid for rGEO?60. Re) 
% T96            : 7   = Tsyganenko [1996] (uses -100.?Dst (nT)?20., 0.5?Pdyn (nPa)?10., |ByIMF| (nT)?10., |BzIMF| (nT)?10. - Valid for rGEO?40. Re) 
% OM97           : 8   = Ostapenko & Maltsev [1997] (uses dst,Pdyn,BzIMF, Kp) 
% T01            : 9   = Tsyganenko [2001] (uses -50.?Dst (nT)?20.,0.5?Pdyn (nPa)?5., |ByIMF| (nT)?5., |BzIMF| (nT)?5., 0.?G1?10., 0.?G2?10. - Valid for xGSM?-15. Re) 
% T01S           : 10  = Tsyganenko [2001] storm  (uses Dst, Pdyn, ByIMF, BzIMF, G2, G3 - there is no upper or lower limit for those inputs - Valid for xGSM?-15. Re) 
% T04            : 11  = Tsyganenko [2004] storm  (uses Dst, Pdyn, ByIMF, BzIMF, W1, W2, W3, W4, W5, W6 - there is no upper or lower limit for those inputs - Valid for xGSM?-15. Re) 
% A00,Paraboloid : 12 =Alexeev [2000] - also known as Paraboloid model - Submitted to ISO  (uses Dens, velo, Dst, BzIMF, AL) 
% TS07d         : 13 - Tsyganenko and Sitnov [2007] - data-based empirical field model (uses Pdyn);
% 
% 
% Notes: 
% 
% 
% when the magnetic field model inputs are outside the allowed range bad data values are returned. 
% When solar wind inputs are required they must be taken in the vicinity of the day side magnetopause and not at L1.
% 
if isempty(kext) % do this first because isnumeric([])==1
    kext = 0;
end
if isnumeric(kext)
    return
end
switch(upper(kext))
    case {'','IGRF'}, kext = 0;
    case {'MF'}, kext = 1;
    case {'T87','T87SHORT'}, kext = 2;
    case {'T87LONG'}, kext = 3;
    case {'T89','T89C'}, kext = 4;
    case {'OPQ'}, kext = 5;
    case {'OPD'}, kext = 6;
    case {'T96'}, kext = 7;
    case {'OM97'}, kext = 8;
    case {'T01'}, kext = 9;
    case {'T01S'}, kext = 10;
    case {'T04'}, kext = 11;
    case {'A00','PARABOLOID'}, kext = 12;
    case {'TS07D'}, kext = 13;
    otherwise
        error('Unknown value for kext "%s" in %s',kext,mfilename);
end
