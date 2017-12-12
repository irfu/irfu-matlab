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
        if isfield(obj.ancillary,'esteptable') && size(obj.ancillary.esteptable,1) == sizeData(1)
          obj.ancillary.esteptable = obj.ancillary.esteptable(idxTmp{1},:);
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
      if ~isempty(obj.ancillary)
        nameFields = fieldnames(obj.ancillary);
        nFields = numel(nameFields);
        for iField = 1:nFields
          eval(['sizeField = size(obj.ancillary.' nameFields{iField} ');'])
          if sizeField(1) == sizeData(1)
            eval(['obj.ancillary.' nameFields{iField} ' = obj.ancillary.' nameFields{iField} '(idx,:);'])
          end
        end
      end
    end

    function [x,y,z] = xyz(obj,varargin)
      % PDIST.XYZ Get xyz coordinates of each detector bin. DSL
      % coordinates. PLEASE REPORT ERRORS.
      %
      %   [x,y,z] = PDIST.xyz(options);
      %    x, y, z - ntx32x16 matrices
      %    options:
      %     'ts' - return x, y, z as TSeries
      %     xyz - transform x,y,z to new xyz = 3x3:          [x,y,z] = PDIST.xyz(xyz);
      %     x,y,z - transform x,y,z to new x,y,z = 1x3 each: [x,y,z] = PDIST.xyz(x,y,z);
      %     'plot' - plots grid, color coded to polar angle 
      %     'squeeze' - squeezes output data [1 32 16] -> [32 16] if PDist      
      %                 only has one time index for example
      
      doReturnTSeries = 0;
      doSqueeze = 0;
      doRotation = 0;
      have_options = 0;
      
      nargs = numel(varargin);      
      if nargs > 0, have_options = 1; args = varargin(:); end
      
      while have_options
        l = 1;
        if isnumeric(args{l})
          if size(args{l}) == [3 3];
            newx = args{l}(1,:);
            newy = args{l}(2,:);
            newz = args{l}(3,:);
            args = args(l+1:end);  
            doRotation = 1;
          elseif numel(args{l}) == 3 && numel(args{l+1}) && numel(args{l+2})
            newx = args{l};
            newy = args{l+1};
            newz = args{l+2};
            args = args(l+3:end);
            doRotation = 1;
          end
        end
        if isempty(args), break, end
        switch(lower(args{1}))
          case 'ts'
            doReturnTSeries = 1;
            args = args(l+1:end);
          case 'squeeze'
            doSqueeze = 1;
            args = args(l+1:end);
          otherwise
            irf.log('warning',sprintf('Input ''%s'' not recognized.',args{1}))
            args = args(l+1:end);
        end        
        if isempty(args), break, end
      end

      nLength = obj.length;
      phi = obj.phi;
      if(size(phi,1)==nLength)
        % Burst data or single data point??
        staticPhi = false;
        azimuthal = phi*pi/180;
      else
        % Srvy data
        staticPhi = true;
        azimuthal = repmat(phi*pi/180, nLength, 1);
      end
      
      theta = obj.theta;
      polar = repmat(theta*pi/180, nLength, 1);
      
      x = NaN(nLength, size(azimuthal,2), size(polar,2));
      y = NaN(nLength, size(azimuthal,2), size(polar,2));
      z = NaN(nLength, size(azimuthal,2), size(polar,2));
      
      for ii = 1:nLength
        [POL, AZ] = meshgrid(polar(ii,:), azimuthal(ii,:));
        X = -sin(POL).*cos(AZ); % '-' because the data shows which direction the particles were coming from
        Y = -sin(POL).*sin(AZ);
        Z = -cos(POL);
        if doRotation % Transform into different coordinate system, (which could be time dependend)
          xX = reshape(X,size(X,1)*size(X,2),1);
          yY = reshape(Y,size(Y,1)*size(Y,2),1);
          zZ = reshape(Z,size(Z,1)*size(Z,2),1);
          newTmpX = [xX yY zZ]*newx';
          newTmpY = [xX yY zZ]*newy';
          newTmpZ = [xX yY zZ]*newz';
          X = reshape(newTmpX,size(X,1),size(X,2));
          Y = reshape(newTmpY,size(X,1),size(X,2));
          Z = reshape(newTmpZ,size(X,1),size(X,2));
        end
        x(ii,:,:) = X;
        y(ii,:,:) = Y;
        z(ii,:,:) = Z;
        if ii==1 && ~doRotation && staticPhi
          % Static
          x = permute(repmat(X, 1, 1, nLength), [3, 1, 2]);
          y = permute(repmat(Y, 1, 1, nLength), [3, 1, 2]);
          z = permute(repmat(Z, 1, 1, nLength), [3, 1, 2]);
          break % Only run through the first "ii"
        end
      end
      %x = permute(x,[1 3 2]);
      %y = permute(y,[1 3 2]);
      %z = permute(z,[1 3 2]);
      
      if doSqueeze
        x = squeeze(x);
        y = squeeze(y);
        z = squeeze(z);
      end
      if doReturnTSeries
        x = irf.ts_scalar(obj.time,x);
        y = irf.ts_scalar(obj.time,y);
        z = irf.ts_scalar(obj.time,z);
      end
    end

    function [vx,vy,vz] = v(obj,varargin)
      % PDIST.V Get velocity corresponding to each detector bin. DSL
      % coordinates. PLEASE REPORT ERRORS.
      %
      %   [vx,vy,vz] = PDIST.v(options);
      %    vx, vy, vz - ntx32x32x16 matrices - km/s
      %    options:
      %     'ts' - return x, y, z as TSeries
      %     xyz - transform x,y,z to new xyz = 3x3:          [x,y,z] = PDIST.xyz(xyz);
      %     x,y,z - transform x,y,z to new x,y,z = 1x3 each: [x,y,z] = PDIST.xyz(x,y,z);
      %     'plot' - plots grid, color coded to polar angle 
      %     'squeeze' - squeezes output data [1 32 32 16] -> [32 32 16] 
      %                 if PDist only has one time index for example
      %
      %   Example:
      %     f = ePDist(100).convertto('s^3/km^6'); % single time PDist
      %     f.data(f.data < 2e3) = NaN; % remove low values
      %     [vx,vy,vz] = f.v('squeeze');
      %     dotsize = 50;
      %     scatter3(vx(:)*1e-3,vy(:)*1e-3,vz(:)*1e-3,f.data(:)*0+dotsize,log10(f.data(:)),'filled'); 
      %     axis equal; colorbar;
      %     vlim = [-5 5]; clim = [3 5];
      %     set(gca,'clim',clim,'xlim',vlim,'ylim',vlim,'zlim',vlim)
      
      doReturnTSeries = 0;
      doSqueeze = 0;
      doRotation = 0;
      have_options = 0;
      
      nargs = numel(varargin);
      if nargs > 0, have_options = 1; args = varargin(:); end
      
      while have_options
        l = 1;
        if isnumeric(args{l})
          if size(args{l}) == [3 3];
            newx = args{l}(1,:);
            newy = args{l}(2,:);
            newz = args{l}(3,:);
            args = args(l+1:end);  
            doRotation = 1;
          elseif numel(args{l}) == 3 && numel(args{l+1}) && numel(args{l+2})
            newx = args{l};
            newy = args{l+1};
            newz = args{l+2};
            args = args(l+3:end);
            doRotation = 1;
          end
        end
        if isempty(args), break, end
        switch(lower(args{1}))
          case 'ts'
            doReturnTSeries = 1;
            args = args(l+1:end);
          case 'squeeze'
            doSqueeze = 1;
            args = args(l+1:end);
          otherwise
            irf.log('warning',sprintf('Input ''%s'' not recognized.',args{1}))
            args = args(l+1:end);
        end
        if isempty(args), break, end
      end

      nLength = obj.length;
      phi = obj.phi;
      if(size(phi,1)==nLength)
        % Burst data or single data point??
        staticPhi = false;
        azimuthal = phi*pi/180;
      else
        % Srvy data with constant Phi for all datapoints
        staticPhi = true;
        azimuthal = repmat(phi*pi/180, nLength, 1);
      end

      theta = obj.theta;
      polar = repmat(theta*pi/180, nLength, 1);

      energy = obj.energy;
      units = irf_units;
      % FIXME: Should we not make sure units of data is as expected? (km/s)
      % FIXME: Should we not use the obj.mass? (ie depending on species of
      % object)
      velocity = sqrt(energy*units.eV*2/units.me)/1000; % km/s
      if(size(unique(velocity),2)==1), staticEnergy=true; else, staticEnergy=false; end

      vx = NaN(size(obj.data));
      vy = NaN(size(obj.data));
      vz = NaN(size(obj.data));

      for ii = 1:nLength
        [VEL,AZ,POL] = meshgrid(velocity(ii,:),azimuthal(ii,:),polar(ii,:));
        %[AZ,VEL,POL] = meshgrid(azimuthal(ii,:),velocity(ii,:),polar(ii,:));
        
        VX = -VEL.*sin(POL).*cos(AZ); % '-' because the data shows which direction the particles were coming from
        VY = -VEL.*sin(POL).*sin(AZ);
        VZ = -VEL.*cos(POL);
        
        % meshgrid permutes the 1st and 2nd indices, 
        % see for example [I1,I2] = meshgrid(1:3,1:2); size(I1), size(I2)
        % the following permutes them back
        % (one can also leave this out and do the following above:
        % [AZ,VEL,POL] = meshgrid(azimuthal(ii,:),velocity(ii,:),polar(ii,:));      
        VX = permute(VX,[2 1 3]);
        VY = permute(VY,[2 1 3]);
        VZ = permute(VZ,[2 1 3]);
        
        if doRotation % Transform into different coordinate system
          VxX = reshape(VX,numel(VX),1);
          VyY = reshape(VY,numel(VX),1);
          VzZ = reshape(VZ,numel(VX),1);

          newTmpX = [VxX VyY VzZ]*newx';
          newTmpY = [VxX VyY VzZ]*newy';
          newTmpZ = [VxX VyY VzZ]*newz';

          VX = reshape(newTmpX,size(VX));
          VY = reshape(newTmpY,size(VY));
          VZ = reshape(newTmpZ,size(VZ));
        end
        vx(ii,:,:,:) = VX;
        vy(ii,:,:,:) = VY;
        vz(ii,:,:,:) = VZ;

        if(ii==1 && ~doRotation && staticEnergy && staticPhi)
          % No rotation (which could be time dependent) and no changing
          % energy/velocity as well as static Phi angles for all
          % datapoints. Use simple static values.
          vx = permute(repmat(VX, 1, 1, 1, nLength), [4, 1, 2, 3]);
          vy = permute(repmat(VY, 1, 1, 1, nLength), [4, 1, 2, 3]);
          vz = permute(repmat(VZ, 1, 1, 1, nLength), [4, 1, 2, 3]);
          break % exit for loop
        end
      end

      if 0 % Diagnostics
        step = 2; %#ok<UNRCH>
        subplot(1,3,1)
        scatter3(VX(1:step:end),VY(1:step:end),VZ(1:step:end),VZ(1:step:end)*0+10,VEL(1:step:end)); axis equal
        subplot(1,3,2)
        scatter3(VX(1:step:end),VY(1:step:end),VZ(1:step:end),VZ(1:step:end)*0+10,AZ(1:step:end)); axis equal
        subplot(1,3,3)
        scatter3(VX(1:step:end),VY(1:step:end),VZ(1:step:end),VZ(1:step:end)*0+10,POL(1:step:end)); axis equal
      end
      if doSqueeze
        vx = squeeze(vx);
        vy = squeeze(vy);
        vz = squeeze(vz);
      end
      if doReturnTSeries
        vx = irf.ts_scalar(obj.time,vx);
        vy = irf.ts_scalar(obj.time,vy);
        vz = irf.ts_scalar(obj.time,vz);
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
      [pitchangles, indPitch] = obj.pitchangle;
      doAverage = 0;
        
      if numel(palim) == 1        
        indPA = find(abs(pitchangles-palim) == min(abs(pitchangles-palim)));
        if nargin>2 && ischar(varargin{1}) && strcmpi(varargin{1},'noav')
          doAverage = 0;
        else 
          doAverage = 1;
        end
      else
        % verify palim
        if(palim(1)>=palim(2)), error('Unexpected palim, should be [minLimit maxLimit] was: [%d %d]', palim(1), palim(2)); end
        indPA = find(bitand(pitchangles(1,:)>=palim(1), pitchangles(1,:)<=palim(2)));
      end
      
      if doAverage
        tmpPA = mean(pitchangles(indPA));
        switch indPitch
          case 2
            tmpData = irf.nanmean(obj.data(:,indPA,:), indPitch);
          case 3
            tmpData = irf.nanmean(obj.data(:,:,indPA), indPitch);
          otherwise
            error('Not yet implemented');
        end
      else
        tmpPA = pitchangles(indPA);
        switch indPitch
          case 2
            tmpData = obj.data(:,indPA,:);
          case 3
            tmpData = obj.data(:,:,indPA);
          otherwise
            error('Not yet implemented');
        end
      end
      
      PD = obj;
      PD.data_ = tmpData;
      PD.depend{indPitch-1} = tmpPA;
    end

    function PD = elim(obj,eint)  
      energy = obj.energy;
      
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
      [energy, indEnergy] = obj.energy;
      energysize = size(energy);
      theta = obj.theta;
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
      PD.depend = {obj.energy};
      PD.representation = {obj.representation{indEnergy-1}};
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
          spec.f = single(obj.energy);
          spec.f_label = {['E_' obj.species(1) ' (eV)']};
        case {'pitchangle','pa'}
          spec.t = obj.time.epochUnix;
          spec.p = double(squeeze(nanmean(obj.data,2))); % nanmean over energies
          %spec.p_label = {'dEF',obj.units};
          spec.f = single(obj.pitchangle);
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
      energy = obj.energy;
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
      
      energy = obj.energy;
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

    function PD = pitchangles(obj,obj1,obj2) %,method
      %PITCHANGLES Calculate pitchangle distribution
      % PitchangleDistribution = Distribution.pitchangles(B,[nangles])
      % PitchangleDistribution = pitchangles(Distribution,B,[nangles])
      % Input: 
      %     B - TSeries of B in dmpa coordinates
      %     nangles - Number of pitch angles or edges of pitchangle bins
      %               default number of pitchangles is 12
      %   See also MMS.GET_PITCHANGLEDIST     
      
      if ~strcmp(obj.type_,'skymap'); error('PDist must be a skymap.'); end 
      
      if nargin<3 || isempty(obj2)
        nangles = 12;
      else 
        nangles = obj2; 
      end       
%       if method % try new method to try to get away the stripes         
%         data_size = size(obj.data);
%         B = obj1.resample(obj.time);
%         [VX,VY,VZ] = obj.v('squeeze');        
%         vx = squeeze(VX(:,1,:,:));
%         vy = squeeze(VY(:,1,:,:));
%         vz = squeeze(VZ(:,1,:,:));
%         vabs = sqrt(vx.^2 + vy.^2 + vz.^2);
%         vxnorm = vx./vabs;
%         vynorm = vy./vabs;
%         vznorm = vz./vabs;
%         Bnorm = irf_norm(B.data);
%         Bxnorm = squeeze(repmat(Bnorm(:,1),1,1,data_size(3),data_size(4)));
%         Bynorm = squeeze(repmat(Bnorm(:,2),1,1,data_size(3),data_size(4)));
%         Bznorm = squeeze(repmat(Bnorm(:,3),1,1,data_size(3),data_size(4)));
%         
%         % pitch angle for each bin (dimension only includes one energy level)
%         pitchangle = acosd(vxnorm.*Bxnorm + vynorm.*Bynorm + vznorm.*Bznorm);
%         pitchangles = nan(data_size);
%         for iE = 1:data_size(2)          
%           pitchangles(:,iE,:,:) = pitchangle;
%         end
%         % sum up f and sort them into the right pitch angle bin
%         pitchangle_edges = linspace(0,180,nangles+1);
% %         %[count,edges,mid,loc] = histcn(pitchangle,pitchangle_edges,pitchangle_edges,pitchangle_edges);        
% %         [count,edges,mid,loc] = histcn(pitchangles(:),pitchangle_edges);
% %         locs = reshape(loc,data_size);
% %         [loc_t,loc_E,loc_az,loc_pol] = ind2sub(data_size,loc);
%         % use irf.nanmean to sum up f for each new bin
%         new_data = nan(obj.length,size(obj.data,2),nangles);
%             
%         for it = 1:data_size(1)
%           for iE = 1:data_size(2)          
%             pitchangles_ = pitchangles(it,iE,:,:);
%             [count,edges,mid,loc] = histcn(pitchangles_(:),pitchangle_edges); 
%             locs = reshape(loc,data_size(3:4));
%             for ipa = 1:nangles         
%               locs_ipa = find(loc == ipa);
%               new_data(it,iE,ipa) = irf.nanmean(obj.data(it,iE,locs_ipa));    
%             end
%           end
%         end
%         
% %         for ipa = 1:nangles      
% %           locs_ = find(loc == ipa);
% %           new_data(:,:,ipa) = irf.nanmean(obj.data(loc==ipa));          
% %         end
%         PD = obj.clone(obj.time,new_data);                
%         PD.depend = {PD.depend{1},repmat(mid{1},obj.length,1)};        
%       else
        [PD,~,~,~] = mms.get_pitchangledist(obj,obj1,'angles',nangles); % - For v1.0.0 or higher data      
%       end
    end

    function PD = einterp(obj,varargin)
      % PDIST.EINTERP Interpolates f to 64 energy channels. 
      %   OBS: ONLY FOR COSMETICS, it makes pitchangle spectrograms 
      %   smoother. Use with caution and always compare to unmodified 
      %   spectrograms and PDist.e64.
      %
      %   PD = PDIST.EINTERP(method);
      %   PD = PDIST.EINTERP;  
      %   method - interpolation method, see interp1, if left empty, default
      %            is 'pchip' which preserves the shape better than 'linear'
      %            and therefore makes pitchangle spectrograms smoother
      %   
      %   Example:
      %     h = irf_plot(3);
      %     tind = 850:950; % memory consuming on long time intervals, and
      %                     % only meaningful to do on shorter times where
      %                     % one can see the jump between energylevels
      %     irf_spectrogram(h(1),ePDist1(tind).pitchangles(gseB1,15).specrec('pa'),'log');
      %     irf_spectrogram(h(2),ePDist1(tind).e64.pitchangles(gseB1,15).specrec('pa'),'log');
      %     irf_spectrogram(h(3),ePDist1(tind).einterp('pchip').pitchangles(gseB1,15).specrec('pa'),'log');
      
      if ~strcmp(obj.type_,'skymap'); error('PDist must be a skymap.'); end 
      if isempty(varargin); method = 'pchip'; else, method = varargin{1}; end

      nt = obj.length;
      old_energies = obj.energy;
      unique_energies = unique(old_energies,'rows');
      new_energy = sort(torow(unique_energies(:)));
      new_energies = repmat(new_energy,nt,1);
      old_data = obj.data;
      new_data = nan(size(old_data,1),numel(new_energy),size(old_data,3),size(old_data,4));
      for it = 1:nt
        for iaz = 1:size(new_data,3)
          for ipol = 1:size(new_data,4)
            new_data(it,:,iaz,ipol) = interp1(old_energies(it,:),old_data(it,:,iaz,ipol),new_energies(it,:),method);
          end
        end
      end
      new_data(new_data<0) = 0; % pchip sometimes give negative values, set these to zero
      PD = obj.clone(obj.time,new_data);
      PD.depend{1} = new_energies;
      PD.ancillary.energy = PD.depend{1};
      PD.ancillary.energy0 = new_energy;
      PD.ancillary.energy1 = new_energy;
    end

    function PD = e64(obj)
      % E64 recompile data into 64 energy channels. Time resolution is
      % halved. Only applies to skymap.
      %   
      %   see also MMS.PSD_REBIN
      
      if ~strcmp(obj.type_,'skymap'); error('PDist must be a skymap.'); end 
      if size(obj.energy,2) == 64; irf.log('warning','PDist already has 64 energy levels.'); end 
      
      if ~any([isfield(obj.ancillary,'energy0') isfield(obj.ancillary,'energy1') isfield(obj.ancillary,'esteptable')]) % construct energy0, energy1, and esteptable 
        esteptable = zeros(obj.length,1);
        [energies,~,esteptable] = unique(obj.energy,'rows'); % consider using legacy
        energy0 = obj.depend{1}(1,:);
        energy1 = obj.depend{1}(2,:);
      end
      
      [pdistr,phir,energyr] = mms.psd_rebin(obj, ...
        TSeries(obj.time, obj.phi), obj.ancillary.energy0, ...
        obj.ancillary.energy1, TSeries(obj.time, obj.ancillary.esteptable));
      PD = obj.clone(pdistr.time, pdistr.data);
      PD.depend{1} = repmat(energyr,PD.length,1);
      PD.ancillary.energy = PD.depend{1}; 
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

    function [en, indEn] = energy(obj)
      % Get energy (Depend_i) of object and its column number (i) in data
      indEn = find(cellfun(@(x) strcmp(x,'energy'),obj.representation));
      if isempty(indEn), error('Did not find "energy" among the dependencies.'); end
      en = obj.depend{indEn};
      if(nargout>1), indEn=indEn+1; end % If requested return its column number in the data (+1 since first column is time).
    end
    function [ph, indPh] = phi(obj)
      % Get phi (Depend_i) of object and its column number (i) in data
      indPh = find(cellfun(@(x) strcmp(x,'phi'),obj.representation));
      if isempty(indPh), error('Did not find "phi" among the dependencies.'); end
      ph = obj.depend{indPh};
      if(nargout>1), indPh=indPh+1; end % If requested return its column number in the data (+1 since first column is time).
    end
    function [th, indTh] = theta(obj)
      % Get phi (Depend_i) of object and its column number (i) in data
      indTh = find(cellfun(@(x) strcmp(x,'theta'),obj.representation));
      if isempty(indTh), error('Did not find "theta" among the dependencies.'); end
      th = obj.depend{indTh};
      if(nargout>1), indTh=indTh+1; end % If requested return its column number in the data (+1 since first column is time).
    end
    function [pitchangle, indPitch] = pitchangle(obj)
      % Get pitchangles (Depend_i) of object and its column number (i) in data
      indPitch = find(cellfun(@(x) strcmp(x,'pitchangle'), obj.representation));
      if isempty(indPitch), error('Did not find "pitchangle" among the dependencies.'); end
      pitchangle = obj.depend{indPitch};
      if(nargout>1), indPitch=indPitch+1; end % If requested return its column number in the data (+1 since first column is time).
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
