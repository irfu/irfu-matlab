function [coef1,coef2,coef3,coef4]=c_efw_calib(isdat_epoch);
%function [c1,c2,c3,c4]=c_efw_calib(isdat_epoch);
% get calibration coefficients for time isdat_epoch (can also be in [yyyy mm dd hh mm ss] format)
% c1..c4 - calibration coefficients [[A_12 E_offs_12_s E_offs_12_xy];[A_34 E_offs_34_s E_offs_34_xy]]
% for calibration coefficient definition see C_DESPIN

if size(isdat_epoch,2) == 6, isdat_epoch=toepoch(isdat_epoch);end % in case tiem given as tiem vector
if isdat_epoch > toepoch([2000 01 01 00 00 00]),
info='EFW calibration coefficients from 2000-01-01';
coef1=[[1,-0.82, -0.25+i*0];[1,-0.67, 0+i*0]]; % sc/1
coef2=[[1,-0.82, 1+i*0];[1,-0.67, 1+i*0]]; % sc/2
coef3=[[1,-0.82, 1+i*0];[1,-0.67, 1+i*0]]; % sc/3
coef4=[[1,-0.82, 1+i*0];[1,-0.67, 1+i*0]]; % sc/4
end

%31 Dec 2000, mp crossing, L-mode
if isdat_epoch > toepoch([2000 12 31 00 00 00]),
info='EFW calibration coefficients from 2000-12-31';
coef1=[[1,0.2, 1.3+0i];[1,-0.2, 1.3+0i]]; % sc/1
coef2=[[1,-0.2, 1.4+.15-0.1i];[1,-0.3, 1.4+0i]]; % sc/2
coef3=[[1,0, 1+.1-0.15i];[1,-.2, 1+0i]]; % sc/3
coef4=[[1,0.9, 0.9+1.1+0.2i];[1,-0.7, 0.9+i*0]]; % sc/4
end

% 2001 01  14 auroral zone, M-mode
if isdat_epoch > toepoch([2001 01 14 03 00 00]),
info='EFW calibration coefficients from 2001-01-14 03:00';
coef1=[[1,0.2, 1-.25+0i];[1,-.6, 1+0i]]; % sc/1
coef2=[[1,-0.15, 1-0.2-0.1i];[1,0.95, 1+0i]]; % sc/2
coef3=[[1,0.35, 1-0.4-0.1i];[1,0.15, 1+0i]]; % sc/3
coef4=[[1,0.35, 1+.35+0i];[1,-1.05, 1+0i]]; % sc/4
end

% 14 Jan 2001, mp crossing, L-mode
if isdat_epoch > toepoch([2001 01 14 08 00 00]),
info='EFW calibration coefficients from 2001-01-14 08:00';
coef1=[[1,0.2, 1.3-0.3-0.1i];[1,-0.3, 1.3+0i]]; % sc/1
coef2=[[1,-0.2, 1.4-0.2-0.2i];[1,-0.5, 1.4+0i]]; % sc/2
coef3=[[1,-0.1, 1.3-0.5-0.15i];[1,0, 1.3+0i]]; % sc/3
coef4=[[1,0.8, 0.9+0.4+i*0];[1,-0.75, 0.9+i*0]]; % sc/4
end

% 26 Jan 2001, mp crossing, L-mode
if isdat_epoch > toepoch([2001 01 26 00 00 00]),
info='EFW calibration coefficients from 2001-01-26 00:00';
coef1=[[1,0.2, 1-0.1-0.1i];[1,-0.3, 1+0i]]; % sc/1
coef2=[[1,-0.2, 1.2-0-0.2i];[1,-0.5, 1.2+0i]]; % sc/2
coef3=[[1,-0.1, 1-0.1-0.05i];[1,0, 1+0i]]; % sc/3
coef4=[[1,0.6, 0.9+0.8+i*0];[1,-0.5, 0.9+i*0]]; % sc/4
end

% 26 Jan 2001, mp crossing, L-mode
if isdat_epoch > toepoch([2001 01 26 10 00 00]),
info='EFW calibration coefficients from 2001-01-26 10:00';
coef1=[[1,0, 1-0.2-0.1i];[1,-0.1, 1+0i]]; % sc/1
coef2=[[1,-0.2, 1.2-0.1-0.2i];[1,-0.5, 1.2+0i]]; % sc/2
coef3=[[1,-.2, 1-0.1-0.05i];[1,-0.1, 1+0i]]; % sc/3
coef4=[[1,0.6, 0.9+0.8+i*0];[1,-0.5, 0.9+i*0]]; % sc/4
end

% 14 Feb 2001, az crossing S, L-mode
if isdat_epoch > toepoch([2001 02 14 00 00 00]),
info='EFW calibration coefficients from 2001-02-14 00:00';
coef1=[[1,0.1, 1.3-0.1+0.i];[1,-0.2, 1.3+0i]]; % sc/1
coef2=[[1,-0.1, 1.4-0.1-0.1i];[1,-0.4, 1.4+0.3+0i]]; % sc/2
coef3=[[1,-0.1, 1.3-0.1i];[1,-0.1, 1.3+0i]]; % sc/3
coef4=[[1,0.6, 0.9+1+0.0i];[1,-0.7, 0.9+0i]]; % sc/4
end

% 14 Feb 2001, az crossing B, L-mode
if isdat_epoch > toepoch([2001 02 14 02 00 00]),
info='EFW calibration coefficients from 2001-02-14 02:00';
coef1=[[1,0.1, 1.3+0.1+0.i];[1,-0.2, 1.3+0i]]; % sc/1
coef2=[[1,-0.1, 1.4-0.1-0.1i];[1,-0.4, 1.4+0.3+0i]]; % sc/2
coef3=[[1,-0.1, 1.3-0.1+0.05i];[1,-0.1, 1.3+0i]]; % sc/3
coef4=[[1,0.6, 0.9+0.8+0.2i];[1,-0.7, 0.9+0i]]; % sc/4
end

% 14 Feb 2001, mp crossing, L-mode
if isdat_epoch > toepoch([2001 02 14 05 00 00]),
info='EFW calibration coefficients from 2001-02-14 05:00';
coef1=[[1,0.2, 1.3-0.3-0.1i];[1,-0.3, 1.3+0i]]; % sc/1
coef2=[[1,-0.2, 1.4-0.1-0.2i];[1,-0.3, 1.4+0i]]; % sc/2
coef3=[[1,-0.1, 1.3-0.1+0.05i];[1,-0.2, 1.3+0i]]; % sc/3
coef4=[[1,0.8, 0.9+1.+0.1i];[1,-0.75, 0.9+0i]]; % sc/4
end

% 18 Feb 2001, az crossing S, L-mode
if isdat_epoch > toepoch([2001 02 18 00 00 00]),
info='EFW calibration coefficients from 2001-02-18 00:00';
coef1=[[1,0.05, 1.-0.1+0.i];[1,-0.1, 1+0i]]; % sc/1
coef2=[[1,-0.1, 1.4-0.1-0.1i];[1,-0.4, 1.4+0.3+0i]]; % sc/2
coef3=[[1,-0.1, 1.-0.15i];[1,-0.1, 1.+0i]]; % sc/3
coef4=[[1,0.5, 1+0.9+0i];[1,-0.6, 1+0i]]; % sc/4
end


% 23 Feb 2001, az crossing, L-mode
if isdat_epoch > toepoch([2001 02 23 00 00 00]),
info='EFW calibration coefficients from 2001-02-23 00:00';
coef1=[[1,0., 1.3-0.1-0.1i];[1,-0.2, 1.3+0i]]; % sc/1
coef2=[[1,-0.1, 1.4-0.1-0.2i];[1,-0.3, 1.4+0.3+0i]]; % sc/2
coef3=[[1,-0.1, 1.3-0.15i];[1,-0.1, 1.3+0i]]; % sc/3
coef4=[[1,0.5, 0.9+1.+0.1i];[1,-0.65, 0.9+0i]]; % sc/4
end

% 26 Feb 2001, az crossing, L-mode
if isdat_epoch > toepoch([2001 02 26 02 00 00]),
info='EFW calibration coefficients from 2001-02-26 00:00';
coef1=[[1,0., 1.3-0.1-0.1i];[1,-0.2, 1.3+0i]]; % sc/1
coef2=[[1,-0.1, 1.4-0.1-0.2i];[1,-0.3, 1.4+0.3+0i]]; % sc/2
coef3=[[1,-0.1, 1.3-0.1+0.05i];[1,-0.1, 1.3+0i]]; % sc/3
coef4=[[1,0.5, 0.9+1.+0.1i];[1,-0.65, 0.9+0i]]; % sc/4
end

% 3 Mar 2001, mp crossing, M-mode (sc2 25Hz sampling) CALIBR NOT DONE
if isdat_epoch > toepoch([2001 03 03 00 00 00]),
info='EFW calibration coefficients from 2001-03-03 00:00';
coef1=[[1,0.3, 1.-0.2+0i];[1,-0.5, 1.+0i]]; % sc/1  DONE
coef2=[[1,-0.1, 1.4-0.1-0.2i];[1,-0.3, 1.4+0.3+0i]]; % sc/2 M-mode 25 sampling
coef3=[[1,0.5, 1.-0.1+0.05i];[1,0, 1+0i]]; % sc/3  DONE
coef4=[[1,0.4, 1+1.+0.i];[1,-1.05, 1+0i]]; % sc/4 DONE
end

% 12 Mar 2001, az crossing, L-mode
if isdat_epoch > toepoch([2001 03 12 00 00 00]),
info='EFW calibration coefficients from 2001-03-12 00:00';
coef1=[[1,0., 1.3-0.1-0.1i];[1,-0.2, 1.3+0i]]; % sc/1
coef2=[[1,-0.1, 1.4-0.1-0.2i];[1,-0.3, 1.4+0.3+0i]]; % sc/2
coef3=[[1,-0.1, 1.3-0.1+0.05i];[1,-0.1, 1.3+0i]]; % sc/3
coef4=[[1,0.5, 0.9+1.+0.1i];[1,-0.65, 0.9+0i]]; % sc/4
end

% 19 Mar 2001, az crossing, L-mode
if isdat_epoch > toepoch([2001 03 19 00 00 00]),
info='EFW calibration coefficients from 2001-03-19 00:00';
coef1=[[1,0., 1.3-0.3-0.1i];[1,-0.2, 1.3+0i]]; % sc/1
coef2=[[1,-0.1, 1.4-0.7-0.1i];[1,-0.3, 1.4+0.3+0i]]; % sc/2
coef3=[[1,-0.1, 1.3-0.1+0.05i];[1,-0.1, 1.3+0i]]; % sc/3
coef4=[[1,1.6, 4.9+1.+0.1i];[1,-0.65, 4.9+0i]]; % sc/4
end

% 28 Mar 2001, az S crossing, L-mode
if isdat_epoch > toepoch([2001 03 28 18 00 00]),
info='EFW calibration coefficients from 2001-03-28 00:00';
coef1=[[1,0., 1.3-0.1-0.1i];[1,-0.2, 1.3+0i]]; % sc/1
coef2=[[1,-0.1, 1.4+0.1-0.2i];[1,-0.4, 1.4+0.3+0i]]; % sc/2
coef3=[[1,-0.1, 1.3-0.1+0.05i];[1,-0.1, 1.3+0i]]; % sc/3
coef4=[[1,0.5, 0.9+1.+0.1i];[1,-0.65, 0.9+0i]]; % sc/4
end

% 28 Mar 2001, az N crossing, L-mode
if isdat_epoch > toepoch([2001 03 28 22 00 00]),
info='EFW calibration coefficients from 2001-03-28 00:00';
coef1=[[1,0., 1.3-0.1-0.1i];[1,-0.2, 1.3+0i]]; % sc/1
coef2=[[1,-0.1, 1.4+0.2-0.2i];[1,-0.5, 1.4+0.3+0i]]; % sc/2
coef3=[[1,-0.1, 1.3-0.1+0.05i];[1,-0.1, 1.3+0i]]; % sc/3
coef4=[[1,0.5, 0.9+1.+0.1i];[1,-0.65, 0.9+0i]]; % sc/4
end

% 28 Apr 2001, az S crossing, L-mode
if isdat_epoch > toepoch([2001 04 28 00 00 00]),
info='EFW calibration coefficients from 2001-04-28 00:00';
coef1=[[1,0., 1.3-0.1-0.1i];[1,-0.2, 1.3+0i]]; % sc/1
coef2=[[1,-0.1, 1.4-0.1i];[1,-0.5, 1.4+0.3+0i]]; % sc/2
coef3=[[1,0.1, 1.3-0.5+0.05i];[1,-0.2, 1.3+0i]]; % sc/3
coef4=[[1,0.7, 0.9+.6+0.2i];[1,-0.8, 0.9+0i]]; % sc/4
end

% 28 Apr 2001, az N crossing, L-mode
if isdat_epoch > toepoch([2001 04 28 21 00 00]),
info='EFW calibration coefficients from 2001-04-28 21:00';
coef1=[[1,-0.1, 1.3-0.2-0.1i];[1,-0.3, 1.3+0i]]; % sc/1
coef2=[[1,-0.3, 1.4-0.1i];[1,-0.7, 1.4+0.3+0i]]; % sc/2
coef3=[[1,-0.1, 1.3-0.6-0.05i];[1,-0.1, 1.3+0i]]; % sc/3
coef4=[[1,0.5, 0.9+.6+0.1i];[1,-0.6, 0.9+0i]]; % sc/4
end

% 11 Jun 2001, mp crossing, L-mode
if isdat_epoch > toepoch([2001 06 11 00 00 00]),
info='EFW calibration coefficients from 2001-06-11 00:00';
coef1=[[1,-0.1, 1.-0.1i];[1,-0.2, 1.+0i]]; % sc/1
coef2=[[1,-0.2, .9-0.1i];[1,-0.4, .9+0.3+0i]]; % sc/2
coef3=[[1,0, 1.-0.4-0.05i];[1,0, 1.+0i]]; % sc/3
coef4=[[1,0.5, 0.9+.6+0i];[1,-0.6, 0.9+0i]]; % sc/4
end


% 27 Jun 2001, plasma sheet, L-mode
if isdat_epoch > toepoch([2001 06 27 00 00 00]),
info='EFW calibration coefficients from 2001-06-27 00:00';
coef1=[[1,-0.1, 1.-0.1i];[1,-0.2, 1.+0i]]; % sc/1
coef2=[[1,-0.3, .9-0.1i];[1,-0.4, .9+0.3+0i]]; % sc/2
coef3=[[1,0.2, 1.-0.4-0.05i];[1,-0.15, 1.+0i]]; % sc/3
coef4=[[1,0.6, 0.9+.6+0i];[1,-0.7, 0.9+0i]]; % sc/4
end

% 30 Jun 2001, magnetopause, L-mode
if isdat_epoch > toepoch([2001 06 27 00 00 00]),
info='EFW calibration coefficients from 2001-06-30 18:00 (magnetopause)';
coef1=[[1,-0.3, 1.-0.1i];[1,-0.25, 1.+0i]]; % sc/1
coef2=[[1,-0.4, .9-0.1i];[1,-0.55, .9+0.2+0i]]; % sc/2
coef3=[[1,-0.2, 1.-0.4-0.05i];[1,0.05, 1.+0i]]; % sc/3
coef4=[[1,0.4, 0.9+.6+0i];[1,-0.6, 0.9+0i]]; % sc/4
end

% 7 Jul 2001, magnetopause, L-mode
if isdat_epoch > toepoch([2001 07 05 00 00 00]),
info='EFW calibration coefficients from 2001-07-05 00:00';
coef1=[[1,-0.1, 1.-0.1i];[1,-0.2, 1.+0i]]; % sc/1
coef2=[[1,-0.3, .9-0.1i];[1,-0.4, .9+0.3+0i]]; % sc/2
coef3=[[1,0., 1.-0.4-0.05i];[1,0., 1.+0i]]; % sc/3
coef4=[[1,0.5, 0.9+.6+0i];[1,-0.5, 0.9+0i]]; % sc/4
end


% 8 Aug 2001, neutral sheet crossing, M,M,M,L-mode
if isdat_epoch > toepoch([2001 08 05 00 00 00]),
info='EFW calibration coefficients from 2001-08-05';
coef1=[[1,0.1, 1.3-0.05i];[1,-0.7, 1.3+0i]]; % sc/1
coef2=[[1,-0.2, 1.4-0.1i];[1,1.4, 1.4+0.3+0i]]; % sc/2
coef3=[[1,0.7, 1.3-0.4-0.05i];[1,0.2, 1.3+0i]]; % sc/3
coef4=[[1,0.35, 0.9+.8+0i];[1,-0.7, 0.9+0i]]; % sc/4
end

% 7 Sep 2001, neutral sheet crossing, M,M,M,M-mode
if isdat_epoch > toepoch([2001 09 07 19 00 00]),
info='EFW calibration coefficients from 2001-09-07 1900 UT';
coef1=[[1,.2, 1.3-0.05i];[1,-0.7, 1.3+0i]]; % sc/1
coef2=[[1,-0.3, 1.4-0.1i];[1,1.2, 1.4+0.3+0i]]; % sc/2
coef3=[[1,.5, 1.3-0.4-0.05i];[1,0.2, 1.3+0i]]; % sc/3
coef4=[[1,0.5, 0.9+.8+0i];[1,-1.0, 0.9+0i]]; % sc/4
end

% 7 Sep 2001, neutral sheet crossing, L,M,L,L-mode
if isdat_epoch > toepoch([2001 09 07 21 00 00]),
info='EFW calibration coefficients from 2001-09-07 2100 UT';
coef1=[[1,-.4, 1.3-0.05i];[1,-0.2, 1.3+0i]]; % sc/1
coef2=[[1,-0.2, 1.4-0.1i];[1,1.4, 1.4+0.3+0i]]; % sc/2
coef3=[[1,-.2, 1.3-0.4-0.05i];[1,-0.2, 1.3+0i]]; % sc/3
coef4=[[1,0.5, 0.9+.8+0i];[1,-0.6, 0.9+0i]]; % sc/4
end

% 1 Oct 2001, neutral sheet crossing, L,M,L,L-mode
if isdat_epoch > toepoch([2001 10 01 00 00 00]),
info='EFW calibration coefficients from 2001-10-01';
coef1=[[1,-.4, 1.3-0.05i];[1,-0.3, 1.3+0i]]; % sc/1
coef2=[[1,-0.3, 1.4-0.1i];[1,1.3, 1.4+0.3+0i]]; % sc/2
coef3=[[1,0, 1.3-0.4-0.05i];[1,-0.4, 1.3+0i]]; % sc/3
coef4=[[1,0.5, 0.9+.8+0i];[1,-0.6, 0.9+0i]]; % sc/4
end

% 11 Oct 2001, neutral sheet crossing, L,M,L,L-mode
if isdat_epoch > toepoch([2001 10 11 16 00 00]),
info='EFW calibration coefficients from 2001-10-11';
coef1=[[1,-.4, 1.3-0.05i];[1,-0.3, 1.3+0i]]; % sc/1
coef2=[[1,-0.3, 1.4-0.1i];[1,1.3, 1.4+0.3+0i]]; % sc/2
coef3=[[1.2,0, 0.5-0.4-0.05i];[1.3,-0.2, 0.5+0i]]; % sc/3
coef4=[[1,0.5, 0.9+.8+0i];[1,-0.6, 0.9+0i]]; % sc/4
end

% 20 Feb 2002, MP, L,M,L,L-mode
if isdat_epoch > toepoch([2002 02 20 00 00 00]),
info='EFW calibration coefficients from 2002-02-20';
coef1=[[1,0, 1.3-0.05i];[1,0.3, 1.3+0i]]; % sc/1
coef2=[[1,-0.3, 1.4-0.1i];[1,1.2, 1.4+0.3+0i]]; % sc/2
coef3=[[1,0.2, 1.3-0.4-0.05i];[1,-0.2, 1.3+0i]]; % sc/3
coef4=[[1,0.2, 0.9+.8+0i];[1,-0.5, 0.9+0i]]; % sc/4
end

% 2 Mar 2002, MP-Xing, L,M,L,L-mode
if isdat_epoch > toepoch([2002 3 2 00 00 00]),
info='EFW calibration coefficients from 2002-03-02';
coef1=[[1,0, 1.3-0.05i];[1,-0.1, 1.3+0i]]; % sc/1
coef2=[[1,-0.15, 1.4-0.1i];[1,1.1, 1.4+0.3+0i]]; % sc/2
coef3=[[1,0.5, 1.3-0.4-0.05i];[1,0.35, 1.3+0i]]; % sc/3
coef4=[[1,0.3, 0.9+.8+0i];[1,-1, 0.9+0i]]; % sc/4
end

% 4 Mar 2002, MP-Xing, L,M,L,L-mode
if isdat_epoch > toepoch([2002 3 4 00 00 00]),
info='EFW calibration coefficients from 2002-03-04';
coef1=[[1,0, 1.3-0.05i];[1,0.4, 1.3+0i]]; % sc/1
coef2=[[1,-0.3, 1.2-0.1i];[1,1.2, 1.2+0.3+0i]]; % sc/2
coef3=[[1,0, 1.3-0.4-0.05i];[1,-.3, 1.3+0i]]; % sc/3
coef4=[[1,0.4, 0.9+.8+0i];[1,-.4, 0.9+0i]]; % sc/4
end

% 6 Mar 2002, magnetopause in M-mode
if isdat_epoch > toepoch([2002 03 6 00 00 00]),
info='EFW calibration coefficients from 2002-03-06';
coef1=[[1,0.1, 1.3-0.05i];[1,-0.7, 1.3+0i]]; % sc/1
coef2=[[1,-0.2, 1.4-0.1i];[1,1.1, 1.4+0.3+0i]]; % sc/2
coef3=[[1,0.2, 1.3-0.4-0.05i];[1,-0.4, 1.3+0i]]; % sc/3
coef4=[[1,0.15, 0.9+.8+0i];[1,-0.6, 0.9+0i]]; % sc/4
end


% 30 Mar 2002, magnetopause in M-mode
if isdat_epoch > toepoch([2002 03 30 00 00 00]),
info='EFW calibration coefficients from 2002-03-30';
coef1=[[1,0.1, 1.3-0.05i];[1,-0.7, 1.3+0i]]; % sc/1
coef2=[[1,-0.2, 1.4-0.1i];[1,1.1, 1.4+0.3+0i]]; % sc/2
coef3=[[1,0.5, 1.3-0.4-0.05i];[1,0.2, 1.3+0i]]; % sc/3
coef4=[[1,0.35, 0.9+.8+0i];[1,-0.9, 0.9+0i]]; % sc/4
end

% 10 Jan 2003, magnetopause in L-mode
if isdat_epoch > toepoch([2003 01 10 00 00 00]),
info='EFW calibration coefficients from 2003-01-10';
coef1=[[1,0., 1.3-0.05i];[1,0., 1.3+0i]]; % sc/1
coef2=[[1,0., 1.4-0.1i];[1,1.1, 1.4+0.3+0i]]; % sc/2
coef3=[[1,0., 1.3-0.4-0.05i];[1,-0.4, 1.3+0i]]; % sc/3
coef4=[[1,0.35, 0.9+.8+0i];[1,-0.7, 0.9+0i]]; % sc/4
end

disp(info);
