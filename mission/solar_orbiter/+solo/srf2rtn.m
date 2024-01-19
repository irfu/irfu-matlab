function TSout = srf2rtn(TSin,direction)
%SOLO.SRF2RTN  Transfrom between SRF and RTN
%
%  TS_RTN = solo.srf2rtn(TS_SRF)
%  TS_SRF = solo.srf2rtn(TS_RTN, -1)
%
% Transform coordinates from the Spacecraft Reference Frame(SRF) to a
% Sun Radial-Tangential-Normal (RTN) frame (or from RTN -> SRF).

%
% Example
% First get an example TSeries of mag data
% file='/data/solo/soar/mag/L2/mag-srf-burst/2020/07/solo_L2_mag-srf-burst_20200731_V02.cdf';
% dObj=dataobj(file);
% t_srf=get_variable(dObj,'EPOCH'); b_srf=get_variable(dObj,'B_SRF');
% TSin_srf=irf.ts_vec_xyz(t_srf.data(1:50),b_srf.data(1:50,:));

%  TSout=srf2rtn(TSin,[direction]); %direction=1 for srf->rtn, direction=-1 for rtn->srf

if nargin<2, direction = 1; end

solo.db_get_metakernel('flown');

% Compute et (SPICE ephemeries time, make use of input in TT2000)
et=irf_time(TSin.time,'EpochTT>tt')';

% The rotation matrix M (note that this is time dependent so it should be a
% 3x3xlength(time) matrix)
if direction==1
  M=cspice_pxform('SOLO_SRF','SOLO_SUN_RTN',et);
elseif direction==-1
  M=cspice_pxform('SOLO_SUN_RTN','SOLO_SRF',et);
else
  error('DIRECTION MUST 1 or -1')
end

out = zeros(size(TSin.data));
for idx=1:3
  out (:,idx) = ...
    squeeze(M(idx,1,:)).*TSin.data(:,1) + ...
    squeeze(M(idx,2,:)).*TSin.data(:,2) + ...
    squeeze(M(idx,3,:)).*TSin.data(:,3);
end

TSout=irf.ts_vec_xyz(TSin.time,out);
TSout.units = TSin.units;
TSout.siConversion = TSin.siConversion;

end


