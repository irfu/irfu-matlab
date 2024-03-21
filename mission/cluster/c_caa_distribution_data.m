function specrec = c_caa_distribution_data(product,detector)
% C_CAA_DISTRIBUTION_DATA   Prepares data to be plotted with
%                           c_caa_distribution_function.
%   Creates a data structure containing necessary data to plot pitch angle
%   distribution functions for various PEACE, CIS and RAPID products.
%
%   data_structure = C_CAA_DISTRIBUTION_DATA(data_product,detector);
%       data_product - for example: 'C3_CP_PEA_3DXPH_PSD'
%       detector - 'LEEA' or 'HEEA', needs to be specified when using for
%           example data product 'C3_CP_PEA_PITCH_FULL_DEFlux' which
%           contains both LEEA and HEEA data separately. Default is [].
%
%    data_structure =
%
%         product: 'C3_CP_PEA_3DXPH_PSD'
%        detector: []
%               t: [69344x1 double]
%              dt: 0.0646
%        en_label: 'Energy [eV]'
%           en_cs: [26x1 single]
%          en_pol: [27x1 single]
%         f_label: 'Pitch angle'
%            f_cs: [15 45 75 105 135 165]
%           f_pol: [0 30 60 90 120 150 180]
%         p_label: 'Log PSD [s^3/km^6]'
%               p: [69344x6x26 double]
%            p_bg: [69344x6x26 double]
%
%   For cross-section plots, the directly given energy levels and pitch
%   angles are used. For polar plots, the energy ranges and pitch angles
%   ranges are used. The resulting data matrix, data_structure.p, has
%   dimensions (N_time,N_pitchangle,N_energylevels).
%
%   Examples:
%       data_structure = C_CAA_DISTRIBUTION_DATA('C3_CP_PEA_3DXPH_PSD');
%       h=c_caa_plot_distribution_function('tint',tint,'polar',data_structure);
%
% See also c_caa_plot_distribution_function.m

distr=product(max(strfind(product,'_'))+1:end);
specrec.product=product;
specrec.detector=[];

if any(strfind(product,'PITCH')) % Something fishy with background level
  disp('Something fishy with background level...')
  if any(strfind(product,'FULL')) && ~exist('detector','var') % contains both detectors
    % Ask for detector is if not given
    detector = input('Which detector? LEEA (1) or HEEA (2)?');
    switch detector
      case 1; detector_str='_LEEA';
      case 2; detector_str='_HEEA';
    end
    specrec.detector=detector_str(2:end);
  elseif any(strfind(product,'FULL')) && exist('detector','var')
    detector_str=['_',detector];
    specrec.detector=detector;
  else
    detector_str='';
  end

  % Load pitch angle data
  [caaData,dataobject,data,data_units]=c_caa_var_get(['Data',detector_str,'__',product]);
  [caabg,~,bg,bg_units]=c_caa_var_get(['BackgroundLevel',detector_str,'__',product]);

  % Pitch angles
  theta=data.dep_x{2}.data(1,:);
  theta_plus=data.dep_x{2}.df.plus;
  theta_minus=data.dep_x{2}.df.minus;
  %nan_theta=isnan(theta);theta(nan_theta)=[];

  % Energy levels
  en=data.dep_x{3}.data(1,:);
  en_plus=data.dep_x{3}.df.plus;
  en_minus=data.dep_x{3}.df.minus;
  nan_en=isnan(en);en(nan_en)=[];en_plus(nan_en)=[];en_minus(nan_en)=[];
  en_pol=[en+en_plus en(end)-en_minus(end)];
  [en, en_cs_order]=sort(en);
  [en_pol, en_pol_order]=sort(en_pol);

  % Sub spin angles
  phi=data.dep_x{1}.data(1,:);
  nan_phi=isnan(phi);phi(nan_phi)=[];

  % Construct subspin resolution data
  dataraw=data.data; dataraw(:,:,:,nan_en)=[];dataraw(:,nan_phi,:,:)=[];
  dataraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
  data=reshape(dataraw,size(dataraw,1)*size(dataraw,2),size(dataraw,3),size(dataraw,4));
  [tt,dtsampling]=subspintime(dataobject,phi); % define subspin time from angle phi.
  data(data==0)=NaN;

  % Same for background level
  bgraw=bg.data; bgraw(:,:,:,nan_en)=[];bgraw(:,nan_phi,:,:)=[];
  bgraw=permute(dataraw,[2 1 3 4]); % permute the order azimuts, time, pitch angle, energy
  bg=reshape(bgraw,size(bgraw,1)*size(bgraw,2),size(bgraw,3),size(bgraw,4));

  % Time
  specrec.t=tt;
  specrec.dt=dtsampling/2;

  % Energy
  specrec.en_label='Energy [eV]';
  specrec.en_cs=en;
  specrec.en_pol=en_pol;

  % Pitch angle
  specrec.f_label='Pitch angle';
  specrec.f_cs=theta;
  specrec.f_pol=[theta-theta_minus theta(end)+theta_plus(end)];

  % Data
  specrec.p_label=['Log ' distr ' [' data_units ']']; % data
  specrec.p=double(data(:,:,en_cs_order)); % average over time
  specrec.p_bg=double(bg(:,:,en_cs_order)); % average over time
end
if any(strfind(product,'PEA_3DXP'))
  % Load pitch angle data
  res=c_caa_construct_subspin_res_data(['Data__', product]);
  bg=c_caa_construct_subspin_res_data(['BackgroundLevel__',product]);

  SEDL=c_caa_var_get(['Sweep_Energy_DeltaLower__', product],'mat');
  SEDU=c_caa_var_get(['Sweep_Energy_DeltaUpper__', product],'mat');
  en_nan=isnan(SEDL(1,:));SEDL(:,en_nan)=[]; [SEDL, SEDL_order]=sort(SEDL(1,2:end));
  en_nan=isnan(SEDU(1,:));SEDU(:,en_nan)=[]; [SEDU, SEDU_order]=sort(SEDU(1,2:end));

  % Store all data in specrec
  specrec.t=res.tt;
  specrec.dt=res.dtsampling/2;

  specrec.en_label='Energy [eV]'; % energy
  specrec.en_cs=res.en(:);
  specrec.en_pol=[res.en-SEDL res.en(end)+SEDU(end)];

  specrec.f_label='Pitch angle'; % pitch angle
  specrec.f_cs=res.theta;
  dtheta=diff(res.theta(1:2))/2;
  specrec.f_pol=[res.theta-dtheta res.theta(end)+dtheta];

  specrec.p_label=['Log ' distr ' [' res.dataunits ']'];
  specrec.p=res.data;
  specrec.p(specrec.p==0)=NaN;
  specrec.p_bg=bg.data;
end
if any([strfind(product,'CODIF_HS'),...
    strfind(product,'CODIF_LS'),...
    strfind(product,'HIA_LS'),...
    strfind(product,'HIA_HS'),...
    strfind(product,'HIA') ])

  % Load pitch angle data
  res=c_caa_construct_subspin_res_data(['3d_ions__', product]);
  caaSEDU=c_caa_var_get(['delta_plus_energy_table__', product]);
  caaSEDL=c_caa_var_get(['delta_minus_energy_table__', product]);
  SEDL1=flipdim(caaSEDL.data(1,:),2); en_nan=isnan(SEDL1);SEDL1(en_nan)=[];
  SEDU1=flipdim(caaSEDU.data(1,:),2); en_nan=isnan(SEDU1);SEDU1(en_nan)=[];

  % Store all data in specrec
  % Time
  specrec.t=res.tt;
  specrec.dt=res.dtsampling/2;

  % Energy
  specrec.en_label='Energy [eV]';
  specrec.en_cs=res.en;
  specrec.en_pol=[res.en-SEDL1 res.en(end)+SEDU1(end)];

  % Pitch angle
  specrec.f_label='Pitch angle';
  specrec.f_cs=res.theta;
  dtheta=diff(res.theta(1:2))/2;
  specrec.f_pol=[res.theta-dtheta res.theta(end)+dtheta];

  % Data
  specrec.p_label=['Log ' distr ' [' res.dataunits ']'];
  specrec.p=res.data;
  specrec.p(specrec.p==0)=NaN;
  specrec.p_bg=NaN(size(specrec.p)); % No bg data?
end
if any(strfind(product,'RAP')) && any(strfind(product,'PAD'))
  [caaData,dataobject,Data,Data_units]=c_caa_var_get(['PAD_Electron_Dif_flux__',product]);

  % Time
  specrec.t=Data.t;
  specrec.dt=(Data.t(2)-Data.t(1))/2;

  % Energy
  en=Data.dep_x{1}.data(1,:);
  en_plus=Data.dep_x{1}.df.plus;
  en_minus=Data.dep_x{1}.df.minus;
  nan_en=isnan(en); en(nan_en)=[]; en_plus(nan_en)=[]; en_minus(nan_en)=[];

  specrec.en_label='Energy [keV]';
  specrec.en_cs=tocolumn(en);
  specrec.en_pol=tocolumn([en(1)-en_minus en+en_plus(1)]);

  % Pitch angles
  specrec.f_label='Pitch angle';
  specrec.f_cs=Data.dep_x{2}.data(1,:);
  specrec.f_pol=[Data.dep_x{2}.data(1,:)-Data.dep_x{2}.df.minus,...
    Data.dep_x{2}.data(1,end)+Data.dep_x{2}.df.plus(end)];

  % Data
  specrec.p_label=['Log ' distr ' [' Data_units ']'];
  dataraw=Data.data; dataraw(:,nan_en,:)=[];
  specrec.p=permute(double(dataraw(:,:,:)),[1 3 2]); % order: time, pitch angle, energy
  specrec.p(specrec.p==0)=NaN;
  specrec.p_bg=specrec.p;
  specrec.p_bg(:)=NaN; % No background data
end


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
    spin_period=dtplus+dtminus;
    dtsampling=spin_period/length(phi);
    for j=length(phi):-1:1
      tt(:,j)=tt(:,1)+double(-dtminus+(j-0.5)*dtsampling);
    end
    tt=reshape(tt',numel(tt),1);
  end
end