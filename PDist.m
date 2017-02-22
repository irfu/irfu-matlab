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
      else, error('3rd input must specify distribution type')
      end
            
      % collect required data, depend        
      switch obj.type_
        case {'skymap'} % construct skymap distribution                
          obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'energy'};
          obj.depend{2} = args{1}; args(1) = []; obj.representation{2} = {'phi'};
          obj.depend{3} = args{1}; args(1) = []; obj.representation{3} = {'theta'};             
        case {'pitchangle'} % construct pitchangle distribution
          obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'energy'};
          obj.depend{2} = args{1}; args(1) = []; obj.representation{2} = {'pitchangle'};                            
        case {'omni'} % construct omni directional distribution
          obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'energy'};
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
    
    function [varargout] = subsref(obj,idx)
    %SUBSREF handle indexing
    switch idx(1).type
      % Use the built-in subsref for dot notation
      case '.'
        [varargout{1:nargout}] = builtin('subsref',obj,idx);
      case '()'
        tmpEpoch = builtin('subsref',obj.time,idx(1));        
        obj.t_ = tmpEpoch;
        idxTmp = repmat({':'}, ndims(obj.data), 1);
        idxTmp(1) = idx(1).subs;
        sizeData = size(obj.data_);
        obj.data_ = obj.data_(idxTmp{:});
        % on depend data      
        
        nDepend = numel(obj.depend);
        for ii = 1:nDepend
          sizeDepend =  size(obj.depend{ii});
          if sizeDepend(1) == 1 % same dependence for all times
            obj.depend_{ii} = obj.depend{ii};
          elseif sizeDepend(1) == sizeData(1)
            obj.depend_{ii} = obj.depend_{ii}(idxTmp{:},:);
          else
            error('Depend has wrong dimensions.')
          end
        end
        if numel(idx) > 1
          obj = builtin('subsref',obj,idx(2:end));
        end
        [varargout{1:nargout}] = obj;
      case '{}'
        error('irf:TSeries:subsref',...
          'Not a supported subscripted reference')
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
        if sizeDepend(1) == 1 % same dependence for all times
          obj.depend_{ii} = obj.depend{ii};
        elseif sizeDepend(1) == sizeData(1)
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
        if sizeField(1) == sizeData(1)
          eval(['obj.ancillary.' nameFields{iField} ' = obj.ancillary.' nameFields{iField} '(idx,:);'])
        end
      end
    end
    
    function PD = palim(obj,palim,varargin)
      % PDIST.PALIM Picks out given pitchangles
      %   distribution type must be 'pitchangle'
      %   PADist = PADist.palim(palims,[arg])
      %     palims - pitchangles, is one angle is given, the closest one is
      %              chosen. If two are equally close, the average is taken,
      %              unless the additional argument 'noav' is given
      %   
      %   PADist.palim([0 90])
      %   PADist.palim(90)
      %   PADist.palim(90,'noav')
      
      if ~strcmp(obj.type,'pitchangle'); error('PDist type must be pitchangle.'); end      
      pitchangles = obj.depend{2};
      doAverage = 0;
        
      if numel(palim) == 1        
        indPA = find(abs(pitchangles-palim) == min(abs(pitchangles-palim)));
        if nargin>2 && ischar(varargin{1}) && strcmpi(varargin{1},'noav')
          doAverage = 0;
        else 
          doAverage = 1;
        end                
      else
        indPA = intersect(find(pitchangles(1,:)>palim(1)),find(pitchangles(1,:)<palim(2)));
      end                  
      
      if doAverage
        tmpPA = mean(pitchangles(indPA));
        tmpData = irf.nanmean(obj.data(:,:,indPA),3);
      else
        tmpPA = pitchangles(indPA);
        tmpData = obj.data(:,:,indPA);
      end      
      
      PD = obj;
      PD.data_ = tmpData;
      PD.depend{2} = tmpPA; 
    end
    function PD = elim(obj,eint)  
      energy = obj.depend{1};
      
      % Picks out energies in an interval, or the closest energy (to be implemented!)
      if numel(eint) == 2
       if or(isempty(obj.ancillary), or(~isfield(obj.ancillary, 'energy0'), ~isfield(obj.ancillary, 'energy1')))
            energytmp0 = energy(1,:);
            energytmp1 = energy(2,:);
            if energytmp0(1) > energytmp1(1)
                tmp = energytmp0;
                energytmp0 = energytmp1;
                energytmp1 = tmp;
            end
            elevels0 = intersect(find(energytmp0>eint(1)),find(energytmp0<eint(2)));
            elevels1 = intersect(find(energytmp1>eint(1)),find(energytmp1<eint(2)));            
       else
            elevels0 = intersect(find(obj.ancillary.energy0>eint(1)),find(obj.ancillary.energy0<eint(2)));
            elevels1 = intersect(find(obj.ancillary.energy1>eint(1)),find(obj.ancillary.energy1<eint(2)));        
       end
       if numel(elevels0) ~= numel(elevels1)
          warning('Energy levels differ for different times. Including the largest interval.')
          elevels = unique([elevels0,elevels1]);
        else
          elevels = elevels0;
        end         
        disp(['Effective eint = [' num2str(min(min(energy(:,elevels))),'%g') ' ' num2str(max(max(energy(:,elevels))),'%g') ']'])
      else
        ediff0 = abs(energy(1,:)-eint);
        ediff1 = abs(energy(2,:)-eint);
        if min(ediff0)<min(ediff1); ediff = ediff0;
        else, ediff = ediff1; end
        elevels = find(ediff==min(ediff));
        disp(['Effective energies alternate in time between ' num2str(energy(1,elevels),'%g') ' and ' num2str(energy(2,elevels),'%g') ''])
      end      
      tmpEnergy = energy(:,elevels);
      tmpData = obj.data(:,elevels,:,:);      
      
      PD = obj;
      PD.data_ = tmpData;
      PD.depend{1} = tmpEnergy;
      if or(isempty(PD.ancillary), or(~isfield(PD.ancillary, 'energy0'), ~isfield(PD.ancillary, 'energy1')))    
          PD.ancillary.energy0 = energytmp0(elevels);
          PD.ancillary.energy1 = energytmp1(elevels);      
      else
          PD.ancillary.energy0 = PD.ancillary.energy0(elevels);
          PD.ancillary.energy1 = PD.ancillary.energy1(elevels);
      end
    end
    function PD = omni(obj)
      % Makes omnidirectional distribution, conserving units.
      
      if ~strcmp(obj.type_,'skymap'); error('PDist must be a skymap.'); end      
      
      dist = obj;
      % define angles
      energysize = size(obj.depend{1});
      theta = obj.depend{3};
      dangle = pi/16;
      lengthphi = 32;

      z2 = ones(lengthphi,1)*sind(theta);
      solida = dangle*dangle*z2;      
      allsolida = repmat(solida,1,1,length(dist.time), energysize(2));
      allsolida = squeeze(permute(allsolida,[3 4 1 2]));
      dists = dist.data.*allsolida;
      omni = squeeze(irf.nanmean(irf.nanmean(dists,3),4))/(mean(mean(solida)));
      
      PD = obj;
      PD.type = 'omni';
      PD.data_ = omni;
      PD.depend = {obj.depend{1}};
      PD.representation = {obj.representation{1}};
      PD.units = obj.units;
      PD.name = 'omni';
    end
    function spec = specrec(obj,varargin)      
      if isempty(varargin); spectype = 'energy'; else, spectype = varargin{1}; end % set default
      
      switch obj.units
        case {'s^3/km^6','s^3/cm^6','s^3/m^6'}
          spec.p_label = {'PSD',obj.units};
        case {'keV/(cm^2 s sr keV)'}
          spec.p_label = {'DEF',obj.units};
        case {'1/(cm^2 s sr keV)'}
          spec.p_label = {'PEF',obj.units};  
        otherwise
          spec.p_label = {obj.units};
      end
      switch spectype
        case 'energy'
          spec.t = obj.time.epochUnix;
          spec.p = double(obj.data);          
          spec.f = single(obj.depend{1});
          spec.f_label = {['E_' obj.species(1) ' (eV)']};
        case {'pitchangle','pa'}
          spec.t = obj.time.epochUnix;
          spec.p = double(squeeze(nanmean(obj.data,2))); % nanmean over energies
          %spec.p_label = {'dEF',obj.units};
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
    function PD = deflux(obj,flagdir)
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
      
      if nargin<2 || flagdir ~= -1
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
      elseif flagdir == -1 && strcmp(obj.units,'keV/(cm^2 s sr keV)')
        irf.log('warning','Converting DEFlux to PSD in SI units');
        tmpData = obj.data/1e12*mm^2*0.53707;
      end    
      energy = obj.depend{1};
      sizeData = size(tmpData);
      reshapedData = reshape(tmpData,sizeData(1),sizeData(2),prod(sizeData(3:end)));
      if size(energy,1) == 1
        matEnergy = repmat(energy,obj.length,1,prod(sizeData(3:end)));
      elseif size(energy,1) == obj.length
        matEnergy = repmat(energy,1,1,prod(sizeData(3:end)));
      end
       
      if nargin<2 || flagdir ~= -1
        reshapedData = reshapedData.*matEnergy.^2;
        tmpData = reshape(reshapedData,sizeData);
        PD = obj;
        PD.data_ = tmpData;
        PD.units = 'keV/(cm^2 s sr keV)';
      elseif flagdir == -1 && strcmp(obj.units,'keV/(cm^2 s sr keV)')
        reshapedData = reshapedData./(matEnergy.^2);
        tmpData = reshape(reshapedData,sizeData);
        PD = obj;
        PD.data_ = tmpData;
        PD.units = 's^3/m^6';  
      else 
      	irf.log('warning','No change to PDist');
      	PD = obj;
      end
    end
    function PD = dpflux(obj,flagdir)
      % Changes units to differential particle flux
      units = irf_units;
      switch obj.species
        case {'e','electrons','electron'}
          mm = units.me/units.mp;          
        case {'i','p','ions','ion'}
          mm = 1;
        otherwise
          error('Units not supported.')
      end 
      
      if nargin<2 || flagdir ~= -1
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
      elseif flagdir == -1 && strcmp(obj.units,'1/(cm^2 s sr keV)')
        irf.log('warning','Converting DPFlux to PSD');
        tmpData = obj.data/1e12*mm^2*0.53707;
      end   
      
      energy = obj.depend{1};
      sizeData = size(tmpData);
      reshapedData = reshape(tmpData,sizeData(1),sizeData(2),prod(sizeData(3:end)));
      if size(energy,1) == 1
        matEnergy = repmat(energy,obj.length,1,prod(sizeData(3:end)));
      elseif size(energy,1) == obj.length
        matEnergy = repmat(energy,1,1,prod(sizeData(3:end)));
      end
      
      if nargin<2 || flagdir ~= -1
        reshapedData = reshapedData.*matEnergy;
        tmpData = reshape(reshapedData,sizeData);
        PD = obj;
        PD.data_ = tmpData;
        PD.units = '1/(cm^2 s sr keV)';  
      elseif flagdir == -1 && strcmp(obj.units,'1/(cm^2 s sr keV)')
        reshapedData = reshapedData./matEnergy;
        tmpData = reshape(reshapedData,sizeData);
        PD = obj;
        PD.data_ = tmpData;
        PD.units = 's^3/m^6';  
      else 
        irf.log('warning','No change to PDist');
        PD = obj;
      end
    end
    function PD = convertto(obj,newunits)
      % Changes units of Pdist. 
      % Accepted inputs 's^3/cm^6', 's^3/km^6', 's^3/m^6', 'keV/(cm^2 s sr keV)',
      % and '1/(cm^2 s sr keV)'
        
      PD = obj;
      % Convert to SI units
      switch obj.units
        case {'s^3/cm^6'}
          PD.data_ = obj.data*1e12;
        case {'s^3/km^6'}
          PD.data_ = obj.data*1e-18;
        case {'s^3/m^6'}
          %PD = PD;
        case {'keV/(cm^2 s sr keV)'}
          PD = obj.deflux(-1);
        case {'1/(cm^2 s sr keV)'}
          PD = obj.dpflux(-1);
        otherwise
          error('Unknown units.')
      end
      PD.units = 's^3/m^6';
      PD.siConversion = 1;
      % Convert to new units
      switch newunits
        case {'s^3/cm^6'}
        	PD.data_ = PD.data*1e-12;
          PD.units = 's^3/cm^6';
          PD.siConversion = 1e12;
        case {'s^3/km^6'}
          PD.data_ = PD.data*1e18;
          PD.units = 's^3/km^6';
          PD.siConversion = 1e-18;
        case {'s^3/m^6'}
          %PD = PD;
        case {'keV/(cm^2 s sr keV)'}
          PD = PD.deflux;
        case {'1/(cm^2 s sr keV)'}
          PD = PD.dpflux;
        otherwise
          error('Units not supported.');
      end
    end          
    function PD = pitchangles(obj,obj1,obj2)
      %PITCHANGLES Calculate pitchangle distribution
      % PitchangleDistribution = Distribution.pitchangles(B,[nangles])
      % PitchangleDistribution = pitchangles(Distribution,B,[nangles])
      % Input: 
      %     B - TSeries of B in dmpa coordinates
      %     nangles - Number of pitch angles or edges of pitchangle bins
      %               default number of pitchangles is 12
      %   See also MMS.GET_PITCHANGLEDIST     
      
      if nargin<3 || isempty(obj2)
        nangles = 12;
      else 
        nangles = obj2; 
      end       
      [PD,~,~,~] = mms.get_pitchangledist(obj,obj1,'angles',nangles); % - For v1.0.0 or higher data      
    end  
    function PD = e64(obj)
      % E64 recompile data into 64 energy channels. Time resolution is
      % halved. Only applies to skymap.
      %   
      %   see also MMS.PSD_REBIN
      
      if ~strcmp(obj.type_,'skymap'); error('PDist must be a skymap.'); end 
      
      [pdistr,phir,energyr] = mms.psd_rebin(obj,TSeries(obj.time,obj.depend{2}),obj.ancillary.energy0,obj.ancillary.energy1,TSeries(obj.time,obj.ancillary.esteptable));
      PD = obj.clone(pdistr.time,pdistr.data);      
      PD.depend{1} = energyr;
      PD.depend{2} = phir.data;  
      
      if isfield(PD.ancillary,'energy0')
        PD.ancillary.energy0 = PD.depend{1};
        PD.ancillary.energy1 = PD.depend{1};
      end
      if isfield(PD.ancillary,'esteptable'); PD.ancillary.esteptable = zeros(PD.length,1); end
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