classdef PDist < TSeries
  % Particle distributions, subclass of TSeries
  
  properties (Access = protected)
    type_
    depend_
    ancillary_
  end
  
  properties (Dependent = true)
    type
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
    function PD = omni(obj,species)
      % Makes omnidirectional distribution
      % Conserves the units
      
      if ~strcmp(obj.type_,'skymap'); error('PDist must be a skymap'); end
      units = irf_units;
      
      diste = obj;
      if isempty(strfind(obj.units,'km')) % Unit conversion from s^3/m^6 to s^3/km^6    
        diste.data = diste.data*1e30; 
      end
      energyspec = obj.depend{1};

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


      switch species
        case {'e','electron','electrons'}
          efluxomni = PSDomni.*energyspec.^2;
          efluxomni = efluxomni; %convert to normal units
        case {'i','p','ion','ions'}
          efluxomni = PSDomni.*energyspec.^2;
          efluxomni = efluxomni/1e6/0.53707; %convert to normal units
      end
      PD = obj;
      PD.type = 'omni';
      PD.data_ = efluxomni;
      PD.depend = {energyspec};
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
          spec.f_label = {'E (eV)'};
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
      switch obj.units
        case {'s^3/m^6'}
          tmpData = obj.data*1e30/1e6/(5.486e-4)^2/0.53707;
        case {'s^3/km^6'}
          tmpData = obj.data/1e6/(5.486e-4)^2/0.53707;
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
      % Distribution.pitchangles(pitchangles,B,[pitchangles])
      %   See also MMS.GET_PITCHANGLEDIST      
      
      %[paddist,theta,energy,tint] = mms.get_pitchangledist(obj,obj.depend{2},obj.depend{3},obj.ancillary.steptable,obj.ancillary.energy0,obj.ancillary.energy1,obj1); % - For v1.0.0 or higher data
      %Dist = PDist(paddist.time,paddist.data,energy,theta);
      
      if ~strcmp(obj.type_,'skymap'); error('PDist must be a skymap'); end
      if isempty(obj2); anglevec = [15 30 45 60 75 90 105 120 135 150 165 180];
      else anglevec = obj2; end
      numechannels = size(obj.data,2);
      lengthphi = size(obj.data,3);
      lengththeta = size(obj.data,4);
        
      pitcha = anglevec-(anglevec(2)-anglevec(1))*0.5;
      phi = obj.depend{2};
      theta = obj.depend{3};


      pdist = obj;
      B = obj1;
      B = B.resample(pdist);
      Bvec = B/B.abs;
      Bvecx = repmat(Bvec.data(:,1),1,numechannels,lengthphi,lengththeta);
      Bvecy = repmat(Bvec.data(:,2),1,numechannels,lengthphi,lengththeta);
      Bvecz = repmat(Bvec.data(:,3),1,numechannels,lengthphi,lengththeta);
      
      x = zeros(length(pdist.time),lengthphi,lengththeta);
      y = zeros(length(pdist.time),lengthphi,lengththeta);
      z = zeros(length(pdist.time),lengthphi,lengththeta);

      for ii = 1:length(pdist.time);
          x(ii,:,:) = -cosd(phi(ii,:)')*sind(theta);
          y(ii,:,:) = -sind(phi(ii,:)')*sind(theta);
          z(ii,:,:) = -ones(lengthphi,1)*cosd(theta);
      end
      
      xt = repmat(x,1,1,1,32);
      xt = squeeze(permute(xt,[1 4 2 3]));
      yt = repmat(y,1,1,1,32);
      yt = squeeze(permute(yt,[1 4 2 3]));
      zt = repmat(z,1,1,1,32);
      zt = squeeze(permute(zt,[1 4 2 3]));

      thetab = acosd(xt.*Bvecx+yt.*Bvecy+zt.*Bvecz);

      c_eval('dist? = pdist.data;',anglevec);
      dist15(thetab > 15) = NaN;
      for jj = 2:(length(anglevec)-1)
        c_eval('dist?(thetab < (?-15)) = NaN;',anglevec(jj));
        c_eval('dist?(thetab > ?) = NaN;',anglevec(jj));
      end 
      dist180(thetab < 165) = NaN;
      c_eval('dist? =  squeeze(irf.nanmean(irf.nanmean(dist?,4),3));',anglevec);

      paddistarr = cat(3,dist15,dist30,dist45,dist60,dist75,dist90,dist105,dist120,dist135,dist150,dist165,dist180);      
      theta = pitcha;
      
      PD = PDist(pdist.time,paddistarr,'pitchangle',obj.depend{1},theta);
      PD.units = obj.units;
    end  
  end
  
end