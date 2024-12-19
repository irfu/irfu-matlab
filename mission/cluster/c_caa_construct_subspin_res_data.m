function res = c_caa_construct_subspin_res_data(variable_name)
%C_CAA_CONSTRUCT_SUBSPIN_RES_DATA construct subspin resolution for PEACE,
%RAPID, CIS data
%
% res = c_caa_construct_subspin_res_data(variable_name)
%
% res.tt   - time axis
% res.en   - energy vector
% res.phi   - azimuth angle vector
% res.theta   - pitch angle vector
% res.data - full data matrix
% res.omni - omni directional energy spectra (summed over pitch angles)
% res.pitch_angle - data summed over energies
% res.phiphi - azimuthal/pitch angle matrices
% res.enlabel - energy label
% res.dataunits - dataunits (for colorbar)
%
% Supported dataset types:
%	PEACE - PADLAR,PADMAR,PITCH_3DR,PITCH_3DX,3DXP
%	RAPID - L3DD,E3DD
%	  CIS - CODIF_HS,CODIF_LS,HIA_HS,HIA_LS
%
% Example:
%   res=c_caa_construct_subspin_res_data('Data__C4_CP_PEA_PITCH_3DRH_PSD')
%
%

if any([strfind(variable_name,'PADLAR') strfind(variable_name,'PADMAR') strfind(variable_name,'PADHAR')])  % PEACE variables PADMARL and PADMARH
  %% PADMARL and PADMARH need to be treated specially in order to recover
  % the subspin timing
  % PEACE sector angles are given in SR2 reference frame!
  if any(strfind(variable_name,'PADLAR'))
    AR_MODE = 1;
  elseif any(strfind(variable_name,'PADMAR'))
    AR_MODE = 2;
  elseif any(strfind(variable_name,'PADHAR'))
    AR_MODE = 4;
  else
    error('c_caa_construct_subspin_res_data: should not be here')
  end
  %XXX : possibly here we should use different values for LAR/MAR/HAR
  theta=15:30:180; % pitch angles to which rebin
  [variable,dataobject,varmat,dataunits]=c_caa_var_get(variable_name); % check that it is loaded in memory
  if isempty(variable)
    error('Cannot load the requested CAA variable')
  end

  phivar=c_caa_var_get(variable.DEPEND_1);phivar=fillval_to_nan(phivar);phi=phivar.data;%nan_phi=isnan(phi);phi(nan_phi)=[];
  phi_dminus = getv(dataobject,phivar.DELTA_MINUS); phi_dminus = phi_dminus.data;
  phi_dplus = getv(dataobject,phivar.DELTA_PLUS); phi_dplus = phi_dplus.data;

  polarvar=c_caa_var_get(variable.DEPEND_2);polar=polarvar.data(1,:);
  envar=c_caa_var_get(variable.DEPEND_3);envar=fillval_to_nan(envar);en=envar.data(1,:);nan_en=isnan(en);en(nan_en)=[];
  enunits=getfield(getv(dataobject,variable.DEPEND_3),'UNITS');
  enlabel=enunits;
  t=varmat.t(:); ndata = length(t);

  dd=regexp(variable_name, '__', 'split');
  ic=str2double(dd{2}(2));
  %Phi Angle in the SR2 co-ordinate system at the center time of the spin
  phi_spin_center = c_caa_var_get(['Angle_SR2phi__' dd{2}]);
  phi_spin_center = phi_spin_center.data(1);

  spin_period = median(diff(t));
  if spin_period > 4.3 || spin_period < 3.7
    spin_period = 4;
    irf_log('proc',sprintf('Using spin period of %.2f sec',spin_period))
  end

  % Sub-spin time
  phi_0 = double(phi - phi_spin_center);
  phi_0(phi_0<-180) = phi_0(phi_0<-180) + 360;
  phi_0(phi_0>180) = phi_0(phi_0>180) - 360;
  tt = repmat(t,1,AR_MODE) + spin_period * phi_0/360; tt = reshape(tt',ndata*AR_MODE,1);
  tt_dminus = spin_period * double(phi_dminus)/360; tt_dminus = reshape(tt_dminus',ndata*AR_MODE,1);
  tt_dplus = spin_period * double(phi_dplus)/360; tt_dplus = reshape(tt_dplus',ndata*AR_MODE,1);

  variable.data(:,:,:,nan_en)=[]; % remove NaN energy data
  variable=fillval_to_nan(variable); % FILLVALs put to NaN
  newdata = zeros(ndata*AR_MODE,size(variable.data,3),size(variable.data,4));
  pitchangle = newdata;
  ii = 1:AR_MODE:ndata*AR_MODE;
  newdata(ii,:,:) = squeeze(variable.data(:,1,:,:));
  if AR_MODE>1
    newdata(ii+1,:,:) = squeeze(variable.data(:,2,:,:));
  end
  if AR_MODE==4
    newdata(ii+2,:,:) = squeeze(variable.data(:,3,:,:));
    newdata(ii+3,:,:) = squeeze(variable.data(:,4,:,:));
  end

  B = c_caa_var_get(irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic),'mat');
  B_SR2=irf_resamp(c_coord_trans('GSE','SR2',B,'cl_id',ic),tt,'nearest'); % sample to t
  b_SR2=irf_norm(B_SR2);

  phi = reshape(phi',ndata*AR_MODE,1);
  cosphi=cosd(phi);sinphi=sind(phi);
  cospolar=cosd(polar);sinpolar=sind(polar);
  for jpolar=1:length(polar)
    % n-vector of given sector in ISR2 ref rframe
    nsector = [sinpolar(jpolar).*cosphi sinpolar(jpolar).*sinphi ones(ndata*AR_MODE,1)*cospolar(jpolar)];
    nparticle = -nsector;
    pitchsector = acosd((dot(nparticle,b_SR2(:,2:4),2)));
    pitchangle(:,jpolar,:)=repmat(pitchsector,[1 1 length(en)]);
  end

  data=ftheta(newdata,pitchangle,theta);

  [en,ix]=sort(en); % sort energy in ascending order
  data=data(:,:,ix); % sort data accordingly

  %[tt,ix]=sort(tt); % sort time in ascending order
  %data=data(ix,:,:); % sort data accordingly


  ind_data = data; ind_data(~isnan(data)) = 1; ind_data(isnan(data)) = 0;
  data_with_nan=data;
  data(isnan(data)) = 0;
  data_omni=reshape(sum(data,2)./sum(ind_data,2),size(data,1),size(data,3));
  data_angle=reshape(sum(data,3)./sum(ind_data,3),size(data,1),size(data,2));

  res.tt = tt;                    % time axis
  res.tt_deltaplus = tt_dplus;
  res.tt_deltaminus = tt_dminus;
  res.en = en;                    % energy levels
  res.phi = phi;                  % azimuthal angles
  res.theta = theta;              % pitch angles
  res.data=data_with_nan;
  res.dtsampling=[];
  res.omni=data_omni;
  res.pitch_angle=data_angle;
  res.phiphi=[];
  res.enlabel=enlabel;
  res.dataunits=dataunits;
  return
elseif any([strfind(variable_name,'PITCH_3DR') strfind(variable_name,'PITCH_3DX')])  % PEACE variable
  [variable,dataobject,peace,dataunits]=c_caa_var_get(variable_name); % check that it is loaded in memory
  if isempty(variable)
    error('Cannot load the requested CAA variable')
  end
  enunits=getfield(getv(dataobject,variable.DEPEND_3),'UNITS');
  enlabel=getfield(getv(dataobject,variable.DEPEND_3),'LABLAXIS');
  enlabel=[enlabel ' [' enunits ']'];
  phi=peace.dep_x{1}.data(1,:);nan_phi=isnan(phi);phi(nan_phi)=[];
  theta=peace.dep_x{2}.data';
  en=peace.dep_x{3}.data(1,:);nan_en=isnan(en);en(nan_en)=[];
  peace.data(:,:,:,nan_en)=[]; % remove NaN energy data
  peace.data(:,nan_phi,:,:)=[]; % remove NaN energy data
  dataraw=peace.data;
  t=peace.t(:);
elseif any([strfind(variable_name,'3DXPL') strfind(variable_name,'3DXPH')])
  %% PEACE_3DXPH does not have pitch angle matrix data, therefore rebinning
  % necessary
  % PEACE sector angles are given in SR2 reference frame!
  theta=15:30:180; % pitch angles to which rebin
  [variable,dataobject,varmat,dataunits]=c_caa_var_get(variable_name); % check that it is loaded in memory
  if isempty(variable)
    error('Cannot load the requested CAA variable')
  end
  phivar=c_caa_var_get(variable.DEPEND_1);phivar=fillval_to_nan(phivar);phi=phivar.data(1,:);nan_phi=isnan(phi);phi(nan_phi)=[];
  polarvar=c_caa_var_get(variable.DEPEND_2);polar=polarvar.data(1,:);
  envar=c_caa_var_get(variable.DEPEND_3);envar=fillval_to_nan(envar);en=envar.data(1,:);nan_en=isnan(en);en(nan_en)=[];
  enunits=getfield(getv(dataobject,variable.DEPEND_3),'UNITS');
  enlabel=enunits;
  variable.data(:,:,:,nan_en)=[]; % remove NaN energy data
  variable.data(:,nan_phi,:,:)=[]; % remove NaN energy data
  % variable.data indices are in order 1) time 2) azimuth 3) polar 4) energy
  variable.data=permute(variable.data,[1 3 2 4]); % permute in order [time, polar, azimuth, energery]
  pitchangle=variable.data;
  t=varmat.t(:);
  [tt]=subspintime(dataobject,phi);
  dd=regexp(variable_name, '__', 'split');
  ic=str2double(dd{2}(2));
  c_eval('caa_load C?_CP_FGM_FULL;',ic);
  c_eval('B=getmat(C?_CP_FGM_FULL,''B_vec_xyz_gse__C?_CP_FGM_FULL'');',ic);
  B_ISR2=irf_resamp(c_coord_trans('GSE','SR2',B,'cl_id',ic),tt,'nearest'); % sample to t
  b_ISR2=irf_norm(B_ISR2);
  cosphi=cos(phi/180*pi);sinphi=sin(phi/180*pi);
  cospolar=cos(polar/180*pi);sinpolar=sin(polar/180*pi);
  for jphi=1:length(phi)
    for jpolar=1:length(polar)
      % nn vector of given sector in ISR2 ref rframe
      nsector=[sinpolar(jpolar).*cosphi(jphi) sinpolar(jpolar).*sinphi(jphi) cospolar(jpolar)];
      nparticle=-nsector;
      pitchsector=acos(dot(b_ISR2(jphi:length(phi):end,2:4),repmat(nparticle,[length(t) 1]),2))*180/pi;
      pitchangle(:,jpolar,jphi,:)=repmat(pitchsector,[1 1 1 length(en)]);
    end
  end
  variable=fillval_to_nan(variable); % FILLVALs put to NaN
  dataraw=ftheta(variable.data,pitchangle,theta);
  dataraw=permute(dataraw,[1 3 2 4]); % permute in order [time, azimuth, pitch, energy]
elseif any([strfind(variable_name,'RAP_L3DD') strfind(variable_name,'RAP_E3DD')]) % RAPID variable
  %% RAPID does not have pitch angle matrix data, therefore rebinning
  % necessary
  [variable,dataobject,rapid,dataunits]=c_caa_var_get(variable_name); % check that it is loaded in memory
  if isempty(variable)
    error('Cannot load the requested CAA variable')
  end
  enunits=getfield(getv(dataobject,variable.DEPEND_1),'UNITS');
  enlabel=getfield(getv(dataobject,variable.DEPEND_1),'LABLAXIS');
  enlabel=[enlabel ' [' enunits ']'];
  phi=rapid.dep_x{2}.data(:);
  if size(phi,1) > 1 && size(phi,2) > 1
    phi = phi(1,:);
  end
  theta=10:20:180; % pitch angles
  en=sqrt(rapid.dep_x{1}.data(1,:).*...
    (rapid.dep_x{1}.data(1,:)+rapid.dep_x{1}.DELTA_PLUS(1,:))); % DELTA_MINUS = 0
  nan_en=isnan(en);
  en(nan_en)=[];
  rapid.data=permute(rapid.data,[1 4 3 2]); % permute in the order time, polar angle, azimuth angle, energy
  rapid.data(:,:,:,nan_en)=[]; % remove NaN energy data
  % read pitch angle information
  variable_pitch_name=['Electron_Pitch_' variable_name(regexp(variable_name,'_C?_')+(1:4)) 'CP_RAP_EPITCH'];
  variable_pitch=c_caa_var_get(variable_pitch_name,'mat');
  rapid_pitch=rapid.data;
  [~,ivar,irap]=intersect(variable_pitch.t,rapid.t);
  variable_pitch.t=variable_pitch.t(ivar);
  variable_pitch.data=variable_pitch.data(ivar,:,:);
  rapid_pitch = rapid_pitch(irap,:,:,:);
  for j=1:size(rapid_pitch,4), rapid_pitch(:,:,:,j)=permute(variable_pitch.data,[1 3 2]);end
  dataraw=ftheta(rapid.data,rapid_pitch,theta);
  dataraw=permute(dataraw,[1 3 2 4]); % permute in order [time, azimuth, pitch, energery]
  t=rapid.t(:);
elseif any([strfind(variable_name,'CODIF_HS') strfind(variable_name,'CODIF_LS') strfind(variable_name,'HIA_HS') strfind(variable_name,'HIA_LS')]) % CIS variable
  %% CODIF_HS does not have pitch angle matrix data, therefore rebinning
  % necessary
  % CIS sector angels are given in ISR2 reference frame and show particle arrival direction!
  theta=11.25:22.5:180; % pitch angles to which rebin
  [variable,dataobject,varmat,dataunits]=c_caa_var_get(variable_name); % check that it is loaded in memory
  if isempty(variable)
    error('Cannot load the requested CAA variable')
  end
  phivar=c_caa_var_get(variable.DEPEND_2);phi=phivar.data';
  polarvar=c_caa_var_get(variable.DEPEND_1);polar=polarvar.data';
  envar=c_caa_var_get(variable.DEPEND_3);en=envar.data';
  enunits=getfield(getv(dataobject,variable.DEPEND_3),'UNITS');
  enlabel=enunits;
  % variable.data indices are in order 1) time 2) polar 3) azimuthal 4) energy
  pitchangle=variable.data;
  t=varmat.data(:,1);
  [tt]=subspintime(dataobject,phi);

  dd=regexp(variable_name, '__', 'split');
  ic=str2double(dd{2}(2));
  c_eval('caa_load C?_CP_FGM_FULL;',ic);
  c_eval('B=getmat(C?_CP_FGM_FULL,''B_vec_xyz_gse__C?_CP_FGM_FULL'');',ic);
  B_ISR2=irf_resamp(c_coord_trans('GSE','ISR2',B,'cl_id',ic),tt,'nearest'); % sample to t
  b_ISR2=irf_norm(B_ISR2);
  cosphi=cos(phi/180*pi);sinphi=sin(phi/180*pi);
  cospolar=cos(polar/180*pi);sinpolar=sin(polar/180*pi);
  for jphi=1:length(phi)
    for jpolar=1:length(polar)
      % nn vector of given sector in SR2 ref rframe
      nparticle=[cospolar(jpolar).*cosphi(jphi) cospolar(jpolar).*sinphi(jphi) sinpolar(jpolar)];
      pitchsector=acos(dot(b_ISR2(jphi:length(phi):end,2:4),repmat(nparticle,[length(t) 1]),2))*180/pi;
      pitchangle(:,jpolar,jphi,:)=repmat(pitchsector,[1 1 1 length(en)]);
    end
  end
  variable=fillval_to_nan(variable); % FILLVALs put to NaN
  dataraw=ftheta(variable.data,pitchangle,theta);
  dataraw=permute(dataraw,[1 3 2 4]); % permute in order [time, azimuth, pitch, energery]
else
  error('This dataset is not know to c_caa_construct_subspin_res_data(). Please add it :-)')
end

% assume dataraw in order [time, azimuth, pitch, energery]
dataraw(dataraw == variable.FILLVAL)=NaN; % FILLVALs put to NaN
[en,ix]=sort(en); % sort energy in ascending order
dataraw=dataraw(:,:,:,ix); % sort data accordingly
dataraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
data=reshape(dataraw,size(dataraw,1)*size(dataraw,2),size(dataraw,3),size(dataraw,4));
tt=repmat(t(:),1,length(phi));
phiphi=tt;


% Fix permutations of matrix to get high time resolution
[tt,dtsampling]=subspintime(dataobject,phi);

for j=length(phi):-1:1
  phiphi(:,j)=phi(j);
end

ind_data = data; ind_data(~isnan(data)) = 1; ind_data(isnan(data)) = 0;
data_with_nan=data;
data(isnan(data)) = 0;
data_omni=reshape(sum(data,2)./sum(ind_data,2),size(data,1),size(data,3));
data_angle=reshape(sum(data,3)./sum(ind_data,3),size(data,1),size(data,2));
phiphi=reshape(phiphi',numel(phiphi),1);

res.tt = tt;                    % time axis
res.tt_deltaplus = [];
res.tt_deltaminus = [];
res.en = en;                    % energy levels
res.phi = phi;                  % azimuthal angles
res.theta = theta;              % pitch angles
res.data=data_with_nan;
res.dtsampling=dtsampling;
res.omni=data_omni;
res.pitch_angle=data_angle;
res.phiphi=phiphi;
res.enlabel=enlabel;
res.dataunits=dataunits;

end

%----------------------------------------------------------
function ftheta=ftheta(fpol,thetapol,theta)
% rebinning to theta vector
% assuming fpol second dimension is polar angle
% rebinnning fpol to ftheta
% assuming thetapol exactly the same size as fpol
ftheta_dim=size(fpol);
ftheta_dim(2)=length(theta);
ftheta=zeros(ftheta_dim);
thetahalfstep=(theta(2)-theta(1))/2;
thetamin=theta-thetahalfstep;
thetamax=theta+thetahalfstep;
for j=1:length(theta)
  ind=(thetapol>thetamin(j)) & (thetapol<thetamax(j));
  fpoltemp=fpol.*ind;
  ind(isnan(fpoltemp))=0;
  fpoltemp(isnan(fpoltemp))=0;
  switch length(ftheta_dim)
    case 3
      ftheta(:,j,:)=sum(fpoltemp,2)./sum(ind,2);
    case 4
      ftheta(:,j,:,:)=sum(fpoltemp,2)./sum(ind,2);
  end
end

end

%----------------------------------------------------------
function [tt,dtsampling]=subspintime(dataobject,phi)
% construct subspin time vector
% phi are azimuthal angles (spin period is divided in the number of azimuth
% angles)
timevar=getv(dataobject,dataobject.VariableAttributes.DEPEND_0{1,2});
tt=timevar.data(:);
tt=repmat(tt,1,length(phi));

if isfield(timevar,'DELTA_PLUS') && isfield(timevar,'DELTA_MINUS')
  if ischar(timevar.DELTA_PLUS)
    deltaplus= getv(dataobject,timevar.DELTA_PLUS);
    dtplus=deltaplus.data(1,:);
    if dtplus>5 % temporary solution for CIS problems
      if (dtplus/2>3.5) && (dtplus/2 < 4.5), dtplus=dtplus/2;
      elseif (dtplus/3>3.5) && (dtplus/3 < 4.5), dtplus=dtplus/3;
      elseif (dtplus/4>3.5) && (dtplus/4 < 4.5), dtplus=dtplus/4;
      end
    end
  elseif isnumeric(timevar.DELTA_PLUS)
    dtplus=timevar.DELTA_PLUS;
  end
  if ischar(timevar.DELTA_MINUS)
    deltaminus= getv(dataobject,timevar.DELTA_MINUS);
    dtminus=deltaplus.data(1,:);
  elseif isnumeric(timevar.DELTA_MINUS)
    dtminus=timevar.DELTA_MINUS;
  end
else
  dtplus=2;
  dtminus=2;
end
spin_period=double(dtplus+dtminus);
dtsampling=spin_period/length(phi);
for j=length(phi):-1:1
  tt(:,j)=tt(:,1)+double(-dtminus+(j-0.5)*dtsampling);
end
tt=reshape(tt',numel(tt),1);
end

%----------------------------------------------------------
function out=fillval_to_nan(in)
% FILLVALs put to NaN
out=in;
out.data(in.data == in.FILLVAL)=NaN;
end