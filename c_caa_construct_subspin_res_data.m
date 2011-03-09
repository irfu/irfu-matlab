function res = c_caa_construct_subspin_res_data(variable_name)
%C_CAA_CONSTRUCT_SUBSPIN_RES_DATA construct subspin resolution for PEACE,
%RAPID, CIS data
%
% res = c_caa_construct_subspin_res_data(variable_name)
%
% res.data - full data matrix 
% res.omni - omni directional energy spectra (summed over pitch angles)
% res.pitch_angle - data summer over energies
% res.phiphi=phiphi - azimuthal/pitch angle matrices
% res.enlabel - energy label
% res.dataunits - dataunits (for colorbar)
%
% Example:
%   res=c_caa_construct_subspin_res_data('Data__C4_CP_PEA_PITCH_3DRH_PSD')


[variable,dataobject]=c_caa_var_get(variable_name); % check that it is loaded in memory
peace=getmat(dataobject,variable_name);
dataunits=getunits(dataobject,variable_name);
enunits=getfield(getv(dataobject,variable.DEPEND_3),'UNITS');
enlabel=getfield(getv(dataobject,variable.DEPEND_3),'LABLAXIS');
enlabel=[enlabel ' [' enunits ']'];
phi=peace.dep_x{1}.data(1,:);
theta=peace.dep_x{2}.data(1,:);
en=peace.dep_x{3}.data(1,:);nan_en=isnan(en);en(nan_en)=[];
data=peace.data;
data(:,:,:,nan_en)=[]; % remove NaN energy data
en=sort(en); % sort energy in ascending order
data=sort(data,4,'ascend'); % sort data accordingly
data=permute(data,[2 1 3 4]); % permute the order azimuts, pitch angle, energy
data=reshape(data,size(data,1)*size(data,2),size(data,3),size(data,4));
tt=repmat(peace.t(:),1,length(phi));
phiphi=tt;
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
data(isnan(data)) = 0;
data_omni=reshape(sum(data,2)./sum(ind_data,2),size(data,1),size(data,3));
data_angle=reshape(sum(data,3)./sum(ind_data,3),size(data,1),size(data,2));
phiphi=reshape(phiphi',numel(phiphi),1);

res.data=data;
res.omni=data_omni;
res.pitch_angle=data_angle;
res.phiphi=phiphi;
res.enlabel=enlabel;
res.dataunits=dataunits;

end

