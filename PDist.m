classdef PDist < TSeries
  % Particle distributions, subclass of TSeries
  % TODO: 
  % e65: collect data into 64 energy levels instead of alternating 32 
  %
  %
  properties (Access = protected)
    type_
    species_
    depend_
    ancillary_
  end
  
  properties (Dependent = true)
    type
    species
    depend
    ancillary
  end
  
  properties (SetAccess = immutable,Dependent = true)
%     tensorOrder
%     tensorBasis
  end
  
  properties (Constant = true, Hidden = true)
%     MAX_TENSOR_ORDER = 2;
%     BASIS = {'xyz','rtp','rlp','rpz','xy','rp'};
%     BASIS_NAMES = {...
%       'Cartesian','Spherical,colatitude', 'Spherical,latitude','Cylindrical',...
%       'Cartesian 2D','Polar 2D'};
  end
  
  properties (SetAccess = protected)
%     representation
  end
  
  properties
%     name = '';
%     units = '';
%     siConversion = '';
%     userData = [];
  end
  
  methods
    function obj = PDist(t,data,varargin) % constructor
      if nargin<2, error('2 inputs required'), end
      
      obj@TSeries(t,data,'to',0);           
      
      args = varargin;     
      if isa(args{1},'char'); obj.type_ = args{1}; args(1) = [];
      else error('3rd input must specify distribution type')
      end
            
      % collect required data, depend        
      switch obj.type_
        case {'skymap'} % construct skymap distribution                
          obj.depend{1} = args{1}; args(1) = []; obj.representation{2} = {'energy'};
          obj.depend{2} = args{1}; args(1) = []; obj.representation{3} = {'phi'};
          obj.depend{3} = args{1}; args(1) = []; obj.representation{4} = {'theta'};             
        case {'pitchangle'} % construct pitchangle distribution
          obj.depend{1} = args{1}; args(1) = []; obj.representation{2} = {'energy'};
          obj.depend{2} = args{1}; args(1) = []; obj.representation{3} = {'pitchangle'};                            
        case {'omni'} % construct omni directional distribution
          obj.depend{1} = args{1}; args(1) = []; obj.representation{2} = {'energy'};
        otherwise 
          warning('Unknown distribution type')
      end
      
      % collect additional data into ancillary
      while ~isempty(args)
        x = args{1}; args(1) = [];
        switch lower(x)
          case {'energy0'}
            obj.ancillary.energy0 = args{1}; args(1) = [];
          case {'energy1'}
            obj.ancillary.energy1 = args{1}; args(1) = [];
          case {'esteptable'}
            obj.ancillary.esteptable = args{1}; args(1) = [];  
        end     
      end
    end    
    % set
    function obj = set.species(obj,value)
      obj.species_ = value;
    end
    function obj = set.type(obj,value)
      obj.type_ = value;
    end
    function obj = set.depend(obj,value)
      obj.depend_ = value;
    end
    function obj = set.ancillary(obj,value)
      obj.ancillary_ = value;
    end    
    % get
    function value = get.species(obj)
      value = obj.species_;
    end
    function value = get.type(obj)
      value = obj.type_;
    end
    function value = get.depend(obj)
      value = obj.depend_;
    end
    function value = get.ancillary(obj)
      value = obj.ancillary_;
    end
    
    function obj = tlim(obj,tint)
      %TLIM  Returns data within specified time interval
      %
      % Ts1 = TLIM(Ts,Tint)
      %
      % See also: IRF_TLIM
      
      % This needs to be modified from TSeries.m to include tlim on depend
      % variables too.
      [idx,obj.t_] = obj.time.tlim(tint);
      sizeData = size(obj.data_);
      nd = ndims(obj.data_);
      if nd>6, error('we cannot support more than 5 dimensions'), end % we cannot support more than 5 dimensions      
      switch nd
        case 2, obj.data_ = obj.data_(idx,:);
        case 3, obj.data_ = obj.data_(idx,:,:,:);
        case 4, obj.data_ = obj.data_(idx,:,:,:,:);
        case 5, obj.data_ = obj.data_(idx,:,:,:,:,:);
        case 6, obj.data_ = obj.data_(idx,:,:,:,:,:,:);
        otherwise, error('should no be here')
      end      
      % on depend data      
      nDepend = numel(obj.depend);
      for ii = 1:nDepend
        sizeDepend =  size(obj.depend{ii});
        if sizeDepend(1) == 1, % same dependence for all times
          obj.depend_{ii} = obj.depend{ii};
        elseif sizeDepend(1) == sizeData(1);                    
          obj.depend_{ii} = obj.depend_{ii}(idx,:);
        else
          error('Depend has wrong dimensions.')
        end
      end
      % on ancillary data
      nameFields = fieldnames(obj.ancillary);
      nFields = numel(nameFields);
      for iField = 1:nFields
        eval(['sizeField = size(obj.ancillary.' nameFields{iField} ');'])
        if sizeField(1) == sizeData(1);
          eval(['obj.ancillary.' nameFields{iField} ' = obj.ancillary.' nameFields{iField} '(idx,:);'])
        end
      end
    end
    
    function PD = palim(obj,palim)
      if strcmp(obj.type,'pitchangle'); error('PDist type must be pitchangle.'); end      
      pitchangles = obj.depend{2};
      
      elevels0 = intersect(find(pitchangles(1,:)>palim(1)),find(pitchangles(1,:)<palim(2)));
      
      if numel(elevels0) ~= numel(elevels1)
        warning('Energy levels differ for different times. Including the largest interval.')
        elevels = unique([elevels0,elevels1]);
      end
      disp(['Effective eint = [' num2str(min(min(energy(1:2,elevels))),'%g') ' ' num2str(max(max(energy(1:2,elevels))),'%g') ']'])
      
      tmpEnergy = energy(:,elevels);
      tmpData = obj.data(:,elevels,:,:);
      
      PD = obj;
      PD.data_ = tmpData;
      PD.depend{1} = tmpEnergy; 
    end
    function PD = elim(obj,eint)  
      energy = obj.depend{1};
      % Picks out energies in an interval, or the closest energy (to be implemented!)
      if numel(eint) == 2
        elevels0 = intersect(find(energy(1,:)>eint(1)),find(energy(1,:)<eint(2)));
        elevels1 = intersect(find(energy(2,:)>eint(1)),find(energy(2,:)<eint(2)));      
        if numel(elevels0) ~= numel(elevels1)
          warning('Energy levels differ for different times. Including the largest interval.')
          elevels = unique([elevels0,elevels1]);
        else
          elevels = elevels0;
        end
        disp(['Effective eint = [' num2str(min(min(energy(1:2,elevels))),'%g') ' ' num2str(max(max(energy(1:2,elevels))),'%g') ']'])
      else
        ediff0 = abs(energy(1,:)-eint);
        ediff1 = abs(energy(2,:)-eint);
        if min(ediff0)<min(ediff1); ediff = ediff0;
        else ediff = ediff1; end        
        elevels = find(ediff==min(ediff));
        disp(['Effective energies alternate in time between ' num2str(energy(1,elevels),'%g') ' and ' num2str(energy(2,elevels),'%g') ''])
      end
      
      tmpEnergy = energy(:,elevels);
      tmpData = obj.data(:,elevels,:,:);
      
      PD = obj;
      PD.data_ = tmpData;
      PD.depend{1} = tmpEnergy;      
    end
    function PD = omni(obj,varargin)
      % Makes omnidirectional distribution
      % Conserves the units
      
      if ~strcmp(obj.type_,'skymap'); error('PDist must be a skymap'); end
      units = irf_units;
      
      diste = obj;
      if isempty(strfind(obj.units,'km')) % Unit conversion from s^3/cm^6 to s^3/km^6    
        diste.data = diste.data*1e30; 
      end
      energyspec = obj.depend{1};
      if size(energyspec,1) == 1
        energyspec = repmat(energyspec,diste.length,1);
      end

      % define angles
      thetae = obj.depend{3};
      dangle = pi/16;
      lengthphi = 32;

      z2 = ones(lengthphi,1)*sind(thetae);
      solida = dangle*dangle*z2;      
      allsolide = zeros(size(diste.data));
     
      for ii = 1:length(diste.time);
          for jj=1:size(energyspec,2);
              allsolide(ii,jj,:,:) = solida;
          end
      end
      
      distes = diste.data.*allsolide;

      % Electron analysis - OMNI
      for ii = 1:length(diste.time);
          disttemp = squeeze(distes(ii,:,:,:));
          PSDomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3))/(mean(mean(solida)));
      end
      if isempty(obj.species)
        species = varargin{1};
      else 
        species = obj.species;
      end
      switch obj.species
        case {'e','electron','electrons'}
          efluxomni = PSDomni.*energyspec.^2;
          efluxomni = efluxomni/1e6/(5.486e-4)^2/0.53707; %convert to normal units
        case {'i','p','ion','ions'}
          efluxomni = PSDomni.*energyspec.^2;
          efluxomni = efluxomni/1e6/0.53707; %convert to normal units
        otherwise          
      end
      
      PD = obj;
      PD.type = 'omni';
      PD.data_ = efluxomni;
      PD.depend = {obj.depend{1};};
      PD.representation = {obj.representation{1},'energy'};
      PD.units = 'keV/(cm^2 s sr keV)';
      PD.name = 'Differential energy flux';
    end
    function spec = specrec(obj,varargin)      
      if isempty(varargin); spectype = 'energy'; else spectype = varargin{1}; end % set default
      switch spectype
        case 'energy'
          spec.t = obj.time.epochUnix;
          spec.p = double(obj.data);
          spec.p_label = {'dEF',obj.units};
          spec.f = single(obj.depend{1});
          spec.f_label = {['E_ ' obj.species(1) ' (eV)']};
        case {'pitchangle','pa'}
          spec.t = obj.time.epochUnix;
          spec.p = double(squeeze(nanmean(obj.data,2))); % nanmean over energies
          spec.p_label = {'dEF',obj.units};
          spec.f = single(obj.depend{2});
          spec.f_label = {'\theta (deg.)'};
        otherwise % energy is default          
          spec.t = obj.time.epochUnix;
          spec.p = double(obj.data);
          spec.p_label = {'dEF',obj.units};
          spec.f = single(obj.depend{1});
          spec.f_label = {'E (eV)'};
      end
    end
    function PD = deflux(obj)
      % Changes units to differential energy flux
      units = irf_units;
      switch obj.species
        case {'e','electrons','electron'}
          mm = units.me/units.mp;          
        case {'i','p','ions','ion'}
          mm = 1;
        otherwise
          error('Units not supported.')
      end  
      switch obj.units
        case {'s^3/cm^6'}
          tmpData = obj.data*1e30/1e6/mm^2/0.53707;
        case {'s^3/m^6'}
          tmpData = obj.data*1e18/1e6/mm^2/0.53707;
        case {'s^3/km^6'}
          tmpData = obj.data/1e6/mm^2/0.53707;
        otherwise
          error('Units not supported.')
      end      
      energy = obj.depend{1};
      sizeData = size(tmpData);
      reshapedData = reshape(tmpData,sizeData(1),sizeData(2),prod(sizeData(3:end)));
      matEnergy = repmat(energy,1,1,prod(sizeData(3:end)));
      reshapedData = reshapedData.*matEnergy.^2;
      tmpData = reshape(reshapedData,sizeData);
      
%       for ii = 1:length(obj.time),
%         energytemp = energy(ii,:)'*ones(1,length(thetae));
%         paddiste.data(ii,:,:) = squeeze(paddiste.data(ii,:,:)).*energytemp.^2;
%       end

      PD = obj;
      PD.data_ = tmpData;
      PD.units = 'keV/(cm^2 s sr keV)';
    end
    function out = peflux(obj)
      % Changes units to differential particle flux
      out = obj;
    end
    function out = psd(obj)
      % Changes units to phase space density (s^3/km^6)
      out = obj;
    end        
    function PD = pitchangles(obj,obj1,obj2)
      %PITCHANGLES Calculate pitchangle distribution
      % Distribution.pitchangles(pitchangles,B,[nangles])
      % Input: 
      %     B - TSeries of B in dmpa coordinates
      %     nangles - Number of pitch angles
      %   See also MMS.GET_PITCHANGLEDIST         
      if isempty(obj2),
          nangles = 12;
      else 
          nangles = obj2; 
      end       
      [PD,~,~,~] = mms.get_pitchangledist(obj,obj1,'angles',nangles); % - For v1.0.0 or higher data      
    end  
    function PD = e64(obj)
      % E64 collect data into 64 energy levels per time
      %   
      %   see also MMS.PSD_REBIN
      
      % MMS.PSD_REBIN should be made compatible with pitch angle
      % distrbutions!!!
      
      [pdistr,phir,energyr] = mms.psd_rebin(obj,TSeries(obj.time,obj.depend{2}),obj.ancillary.energy0,obj.ancillary.energy1,TSeries(obj.time,obj.ancillary.energyStepTable));
      PD = obj.clone(pdistr.time,pdistr.data);      
      PD.depend{1} = energyr;
      PD.depend{2} = phir.data;  
      %if isfield(PD.ancillary,'energy1'); PD.ancillary = rmfield(PD.ancillary,'energy1'); end
      if isfield(PD.ancillary,'energy0'); PD.ancillary = setfield(PD.ancillary,'energy0',PD.depend{1}); end
      if isfield(PD.ancillary,'energyStepTable'); PD.ancillary = setfield(PD.ancillary,'energyStepTable',zeros(PD.length,1)); end
%       energy = obj.depend{1};
%       [newEnergy,energyOrder] = sort([energy(1,:) energy(2,:)]);
%             
%       sizeData = size(obj.data);
%       sizeNewData = sizeData; sizeNewData(1) = fix(sizeNewData(1)/2); sizeNewData(2) = 64;
%       tmpData = nan(sizeNewData);
%       tmpData(:,1:32,:,:) = obj.data(1:2:end-1,:,:,:);
%       tmpData(:,33:64,:,:) = obj.data(2:2:end,:,:,:);
%       tmpData = tmpData(:,energyOrder,:,:);
%       newEnergy = nan(sizeNewData(1),64);
%       
%                   
% 
%       % Define new times
%       deltat = median(diff(obj.time.epochUnix))/2;
%       newTimes = pdist.time(1:2:end-1)+deltat;      
%       phir = nan(length(newtimes),32);
%       newelnum = 1;
% 
%       phis = circshift(phi.data,1,2);
%       phis(:,1) = phis(:,1)-360;
% 
%       for ii=1:2:sizeData(1)-1;
%           if phi.data(ii,1) > phi.data(ii+1,1), 
%               phir(newelnum,:) = (phi.data(ii,:)+phis(ii+1,:))/2;
%               pdisttemp = circshift(squeeze(pdist.data(ii+1,:,:,:)),1,2);       
%           else
%               phir(newelnum,:) = (phi.data(ii,:)+phi.data(ii+1,:))/2;
%           end
%           newelnum = newelnum+1;
%       end
% 
%       phir = TSeries(newtimes,phir);
%       pdistr = TSeries(newtimes,pdistr);
%       toc;

    end
    function m = mass(obj)
      % Get mass of species
      units = irf_units;
      switch obj.species
        case {'e','electrons','electron'}
          m = units.me;
        case {'i','p','ions','ion'}
          m = units.mp;
        otherwise
          error('Species not supported.')
      end 
    end
    function e = energy(obj)
      % Get energy of object
      %indE = find(strcmp(obj.representation,'energy'))
      e = obj.depend{1};
    end
    function moms = moments(obj,varargin)
      % MOMENTS compute moments from the FPI particle phase-space densities 
      %
      % For brst mode data
      % particlemoments = PDist.moments(phi,theta,stepTable,energy0,energy1,SCpot,particle,option,option_value)
      %
      % For fast mode data
      % particlemoments = PDist.moments(phi,theta,energy,SCpot,particle,'fast',option,option_value)
      %
      % Input:
      %   pdist - TSeries of the full particle distribution of electrons or ions
      %   (must be in s^3/cm^6) (burst and fast)
      %   phi - TSeries of all phi angles of distribution for burst data. 1D array or
      %   structure for fast data.
      %   theta - 1D array or structure of theta angles (burst and fast)
      %   stepTable - TSeries of stepping table between energies (burst)
      %   energy0 - 1D array or structure of energy table 0 (burst)
      %   energy1 - 1D array or structure of energy table 1 (burst)
      %   energy - 1D array or structure of energy table (fast)
      %   SCpot - TSeries of spacecraft potential (burst and fast). 
      %   (Make sure sign is correct, should be typically positive)
      %   particle - indicate particle type: 'electron' or 'ion'
      %
      %   See Example_MMS_EDRsignatures for example of loading the necessary data 
      %   and running the function.
      %
      % Optional Inputs:
      %   'energyrange' - set energy range in eV to integrate over [E_min E_max].
      %   energy range is applied to energy0 and the same elements are used for energy1 to 
      %   ensure that the same number of points are integrated over. 
      %   'noscpot' - set to 1 to set spacecraft potential to zero. Calculates moments without
      %   correcting for spacecraft potential. 
      %   'enchannels' - set energy channels to integrate over [min max]; min and max
      %   between must be between 1 and 32.
      %   'partialmoms' - use a binary array (or TSeries) (pmomsarr) to select which psd points are used
      %   in the moments calculation. pmomsarr must be a binary array (1s and 0s, 1s correspond to points used).
      %   Array (or data of TSeries) must be the same size as pdist.data. For
      %   examples see Example_MMS_partialmoments.
      %
      % Output: 
      %   psd_moments - structure containing the particle moments: density, bulk
      %   velocity, pressure, temperature, and particle heat flux (n_psd, V_psd, P_psd, T_psd, and H_psd,
      %   respectively) as TSeries'. For temperature and
      %   pressure tensors the order of the columns is XX, XY, XZ, YY, YZ, ZZ.
      %
      % See also MMS.PSD_MOMENTS
      %
      % Notes: 
      % Regarding the spacecraft potential, the best estimate of is -1.2*(probe
      % to spacecraft voltage)+MMSoffset. Note that in most plasmas the spacecraft
      % potential is positive. E.g.
      % ic = 1,2,3, or 4;
      % c_eval('do = dataobj(''data/mms?_edp_brst_l2_scpot_20151202011414_v1.0.0.cdf'');',ic);
      % c_eval('SCpot = mms.variable2ts(get_variable(tmpDataObj,''mms?_edp_psp''));',ic);
      % offset1 = 1.3; offset2 = 1.5; offset3 = 1.2; offset4 = 0.0; %For v1 data
      % c_eval('SCpot.data = -SCpot.data*1.2+offset?;',ic);
      % Apply correction for input. Correction is not applied in this script. 
      % This correction is applied to v2 spacecraft potential so use 
      % c_eval('SCpot = mms.variable2ts(get_variable(tmpDataObj,''mms?_edp_scpot_fast_l2''));',ic);
      %
      % Currently the heat flux vector does not match with the FPI ion moments. Currently
      % using Eq. (6.8) of Analysis Methods for Multi-Spacecraft Data. This needs
      % to be investigated further. 
 
    end
  end
  methods (Static)
    function newUnits = changeunits(from,to)
      
    end
  end
end