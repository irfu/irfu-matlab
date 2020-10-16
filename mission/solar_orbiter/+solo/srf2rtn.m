function TSout = srf2rtn(TSin,direction)
%A function to transform coordinates from the Spacecraft Reference Frame
%(SRF) to a Sun Radial-Tangential-Normal (RTN) frame (or from RTN -> SRF).

% Example
% First get an example TSeries of mag data
% file='/data/solo/soar/mag/L2/mag-srf-burst/2020/07/solo_L2_mag-srf-burst_20200731_V02.cdf';
% dObj=dataobj(file);
% t_srf=get_variable(dObj,'EPOCH'); b_srf=get_variable(dObj,'B_SRF');
% TSin_srf=irf.ts_vec_xyz(t_srf.data(1:50),b_srf.data(1:50,:));  

%  TSout=srf2rtn(TSin,[direction]); %direction=1 for srf->rtn, direction=-1 for rtn->srf

if nargin<2, direction = 1; end

% Find the most recent metakernel file and load it
sharedPath = '/share/SPICE/'; % SPICE kernels for different missions are found in this folder on IRFU servers
dirs = dir([sharedPath,'Solar-Orbiter/kernels/mk/*pred-mk_v*.tm']);
if size(dirs, 1) > 1
    % Multiple kernels could be found if executing this script at the same time as syncing new kernel files
    error('Found multiple metakernels, please check your folder.');
end
kernelFile = [dirs.folder, filesep, dirs.name];
cspice_furnsh(kernelFile);

% Compute et (SPICE ephemeries time, make use of input in TT2000)
et=irf_time(TSin.time,'EpochTT>tt')';

% The rotation matrix M (note that this is time dependent so it should be a
% 3x3xlength(time) matrix)
if direction==1
    M=cspice_pxform('SOLO_SRF','SOLO_SUN_RTN',et);
elseif direction==-1
    M=cspice_pxform('SOLO_SUN_RTN','SOLO_SRF',et);
end

%This for loop can probably be optimized...
%for j=1:length(TSin.time)
%    out(j,:)=M(:,:,j)*TSin.data(j,1:3)';
%end

%..to this perhaps, but for long time series the above might actually be better
out_tmp=pagetranspose(pagemtimes(M,'none',TSin.data,'transpose'));
out=out_tmp(:,:,1); %Someone better skilled with matrix multiplication might be able to do something more efficient here

TSout=irf.ts_vec_xyz(TSin.time,out);

end


