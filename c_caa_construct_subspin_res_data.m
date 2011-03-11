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
% res.pitch_angle - data summer over energies
% res.phiphi=phiphi - azimuthal/pitch angle matrices
% res.enlabel - energy label
% res.dataunits - dataunits (for colorbar)
%
% Example:
%   res=c_caa_construct_subspin_res_data('Data__C4_CP_PEA_PITCH_3DRH_PSD')

if strfind(variable_name,'PITCH_3DR'), % PEACE variable
  [variable,dataobject]=c_caa_var_get(variable_name); % check that it is loaded in memory
  peace=getmat(dataobject,variable_name);
  dataunits=getunits(dataobject,variable_name);
  enunits=getfield(getv(dataobject,variable.DEPEND_3),'UNITS');
  enlabel=getfield(getv(dataobject,variable.DEPEND_3),'LABLAXIS');
  enlabel=[enlabel ' [' enunits ']'];
  phi=peace.dep_x{1}.data(1,:);
  theta=peace.dep_x{2}.data(1,:);
  en=peace.dep_x{3}.data(1,:);nan_en=isnan(en);en(nan_en)=[];
  dataraw=peace.data;
  t=peace.t(:);
elseif strfind(variable_name,'RAP_L3DD'), % RAPID variable
  % RAPID does not have pitch angle matrix data, therefore rebinning
  % necessary
  [variable,dataobject]=c_caa_var_get(variable_name); % check that it is loaded in memory
  rapid=getmat(dataobject,variable_name);
  dataunits=getunits(dataobject,variable_name);
  enunits=getfield(getv(dataobject,variable.DEPEND_1),'UNITS');
  enlabel=getfield(getv(dataobject,variable.DEPEND_1),'LABLAXIS');
  enlabel=[enlabel ' [' enunits ']'];
  phi=rapid.dep_x{2}.data(1,:);
  theta=10:20:180; % pitch angles
  en=rapid.dep_x{1}.data(1,:);nan_en=isnan(en);en(nan_en)=[];
  rapid.data=permute(rapid.data,[1 4 3 2]); % permute in the order time, polar angle, azimuth angle, energy
  % read pitch angle information
  variable_pitch_name=['Electron_Pitch_' variable_name(regexp(variable_name,'_C?_')+(1:4)) 'CP_RAP_EPITCH'];
  [variable_pitch,dataobject_pitch]=c_caa_var_get(variable_pitch_name);
  rapid_pitch=rapid.data;
  for j=1:size(rapid_pitch,4), rapid_pitch(:,:,:,j)=permute(variable_pitch.data,[1 3 2]);end
  dataraw=ftheta(rapid.data,rapid_pitch,theta);
  dataraw=permute(dataraw,[1 3 2 4]); % permute in order [time, azimuth, pitch, energery]
  t=rapid.t(:);
end

dataraw(dataraw == variable.FILLVAL)=NaN;
dataraw(:,:,:,nan_en)=[]; % remove NaN energy data
[en,ix]=sort(en); % sort energy in ascending order
dataraw=dataraw(:,:,:,ix); % sort data accordingly
dataraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, pitch angle, energy
data=reshape(dataraw,size(dataraw,1)*size(dataraw,2),size(dataraw,3),size(dataraw,4));
tt=repmat(t(:),1,length(phi));
phiphi=tt;


% Fix permutations of matrix to get high time resolution

timevar=getv(dataobject,dataobject.VariableAttributes.DEPEND_0{1,2});
if isfield(timevar,'DELTA_PLUS') && isfield(timevar,'DELTA_MINUS')
  if ischar(timevar.DELTA_PLUS)
    deltaplus= getv(dataobject,timevar.DELTA_PLUS);
    deltaminus= getv(dataobject,timevar.DELTA_MINUS);
    dtplus=deltaplus.data(1,:);
    dtminus=deltaminus.data(1,:);
  elseif isnumeric(timevar.DELTA_PLUS)
    dtplus=getv(dataobject,timevar.DELTA_PLUS);
    dtminus=getv(dataobject,timevar.DELTA_MINUS);
  end
else
  dtplus=2;
  dtminus=2;
end
spin_period=dtplus+dtminus;
dtsampling=spin_period/length(phi);
for j=length(phi):-1:1,
  tt(:,j)=tt(:,1)+double(-dtminus+(j-0.5)*dtsampling);
  phiphi(:,j)=phi(j);
end
tt=reshape(tt',numel(tt),1);
ind_data = data; ind_data(~isnan(data)) = 1; ind_data(isnan(data)) = 0;
data_with_nan=data;
data(isnan(data)) = 0;
data_omni=reshape(sum(data,2)./sum(ind_data,2),size(data,1),size(data,3));
data_angle=reshape(sum(data,3)./sum(ind_data,3),size(data,1),size(data,2));
phiphi=reshape(phiphi',numel(phiphi),1);

res.tt = tt;
res.en = en;
res.phi = phi;
res.theta = theta;
res.data=data_with_nan;
res.dtsampling=dtsampling;
res.omni=data_omni;
res.pitch_angle=data_angle;
res.phiphi=phiphi;
res.enlabel=enlabel;
res.dataunits=dataunits;

end

function ftheta=ftheta(fpol,thetapol,theta)
% rebinning to theta vector
% assuming fpol second dimension is polar angle
% rebinnning fpol to ftheta
ftheta_dim=size(fpol);
ftheta_dim(2)=length(theta);
ftheta=zeros(ftheta_dim);
thetahalfstep=(theta(2)-theta(1))/2;
thetamin=theta-thetahalfstep;
thetamax=theta+thetahalfstep;
for j=1:length(theta),
  ind=(thetapol>thetamin(j)).*(thetapol<thetamax(j));
  switch length(ftheta_dim)
    case 4
      fpoltemp=fpol.*ind;
      fpoltemp(isnan(fpoltemp))=0;
      ftheta(:,j,:,:)=sum(fpoltemp,2)./sum(ind,2);
  end
end

end