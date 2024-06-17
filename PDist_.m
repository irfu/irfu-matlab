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
        case {'line (reduced)'} % % construct 1D distribution, through integration over the other 2 dimensions
          obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'velocity'};          
        case {'plane (reduced)'} % construct 2D distribution, either through integration or by taking a slice
          obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'velocity1'};
          obj.depend{2} = args{1}; args(1) = []; obj.representation{2} = {'velocity2'};      
        case {'plane (slice)'} % construct 2D distribution, either through integration or by taking a slice
          obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'velocity1'};
          obj.depend{2} = args{1}; args(1) = []; obj.representation{2} = {'velocity2'};        
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
    
    function varargout = subsref(obj,idx)
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

        % pick out correct indices for ancillary data time tables, nb. this
        % assumes anything with 'number of rows' = PDist.length is a timetable
        ancillary_fieldnames = fieldnames(obj.ancillary);
        new_ancillary_data = obj.ancillary;
        for iField = 1:numel(ancillary_fieldnames)
          field_data = getfield(obj.ancillary,ancillary_fieldnames{iField});          
          if isnumeric(field_data) && size(field_data,1) == sizeData(1) % has the same number of rows as the PDist has time indices, assume each row corresponds to the same time index
            new_ancillary_data = setfield(new_ancillary_data,ancillary_fieldnames{iField},field_data(idxTmp{1},:,:,:,:,:,:)); % repeated :,:,:,:,:,:, used to support multidimensional data
          end
        end
        obj.ancillary = new_ancillary_data;
                
        %varargout{1} = obj(idx);
        %dbstack
        %idx.subs
        %nargout
        if numel(idx) > 1  
          nargout_str = [];
          for inout = 1:nargout
            c_eval('nargout_str = [nargout_str ''tmp_vout?,''];',inout)
          end
          nargout_str = ['[' nargout_str(1:end-1) ']'];
          varargout_str = ['{' nargout_str(2:end-1) '}'];
          %disp(sprintf('%s = builtin(''subsref'',obj,idx(2:end));',nargout_str))
          eval(sprintf('%s = builtin(''subsref'',obj,idx(2:end));',nargout_str))
          %disp(sprintf('varargout = %s;',varargout_str))
          eval(sprintf('varargout = %s;',varargout_str))
%           if nargout == 1
%             tmp_vout = builtin('subsref',obj,idx(2:end));
%             obj = {tmp_vout};
%           elseif nargout == 2
%             [tmp_vout1,tmp_vout2] = builtin('subsref',obj,idx(2:end));
%             obj = {tmp_vout1,tmp_vout2};
%           elseif nargout == 3
%             [tmp_vout1,tmp_vout2,tmp_vout3] = builtin('subsref',obj,idx(2:end));
%             obj = {tmp_vout1,tmp_vout2,tmp_vout3};
%           end
        else
          [varargout{1:nargout}] = obj;
        end                
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
          if all(size(args{l}) == [3 3])
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

      phi = TSeries(obj.time,obj.depend{1,2});
      azimuthal = phi.data*pi/180;      
      
      theta = obj.depend{1,3};
      polar = repmat(theta*pi/180,obj.length,1);            
      
      x = nan(obj.length,size(azimuthal,2),size(polar,2));
      y = nan(obj.length,size(azimuthal,2),size(polar,2));
      z = nan(obj.length,size(azimuthal,2),size(polar,2));
      
      
      for ii = 1:length(obj.time)
        [POL,AZ] = meshgrid(polar(ii,:),azimuthal(ii,:));
        X = -sin(POL).*cos(AZ); % '-' because the data shows which direction the particles were coming from
        Y = -sin(POL).*sin(AZ);
        Z = -cos(POL);

                
        if doRotation % Transform into different coordinate system
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
          if all(size(args{l}) == [3 3])
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

      phi = TSeries(obj.time,obj.depend{1,2});
      azimuthal = phi.data*pi/180;      
      
      theta = obj.depend{1,3};
      polar = repmat(theta*pi/180,obj.length,1);      
      
      energy = obj.depend{1};
      units = irf_units;
      velocity = sqrt(energy*units.eV*2/units.me)/1000; % km/s
      
      vx = NaN*obj.data;
      vy = NaN*obj.data;
      vz = NaN*obj.data;
            
      for ii = 1:length(obj.time)
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
    function PD = reduce(obj,dim,x,varargin) 
      %PDIST.REDUCE Reduces (integrates) 3D distribution to 1D (line).      
      %   Example:
      %     f1D = iPDist1.reduce('1D',dmpaB1,'vint',[0 10000]);
      %     irf_spectrogram(irf_panel('f1D'),f1D.specrec('velocity_1D'));
      %
      %   Options:
      %     'vint'   - set limits on the from-line velocity to get cut-like
      %                distribution
      %     'nMC'    - number of Monte Carlo iterations used for integration,
      %                for default number see IRF_INT_SPH_DIST
      %     'weight' - how the number of MC iterations per bin is weighted, can be
      %                'none' (default), 'lin' or 'log'
      %     'vg'     - array with center values for the projection velocity
      %                grid in [km/s], determined by instrument if omitted
      %     'scpot'  - sets all values below scpot to zero. If a
      %                significant photoelectron population remains, try
      %                applying some margin on scPot: scPot -> scPot*1.5
      %                To be implemented: correct energy table: energy -> energy-scpot
      %
      %   Reduction to 2D (plane) is to be implemented.
      %
      % This is a shell function for irf_int_sph_dist
      % See also: MMS.PLOT_INT_DISTRIBUTION, IRF_INT_SPH_DIST
      % MMS.PLOT_INT_PROJECTION
      
      %% Input
      [ax,args,nargs] = axescheck(varargin{:});
      irf.log('warning','Please verify that you think the projection is done properly!');
      if isempty(obj); irf.log('warning','Empty input.'); return; else, dist = obj; end
      
      % Check to what dimension the distribution is to be reduced
      if any(strcmp(dim,{'1D','2D'}))
        dim = str2num(dim(1)); % input dim can either be '1D' or '2D'
      else
        error('First input must be a string deciding projection type, either ''1D'' or ''2D''.')
      end
      
      if dim == 1 % 1D: projection to line
        if isa(x,'TSeries')
          xphat_mat = x.resample(obj).norm.data; 
        elseif isnumeric(x) && numel(size(x) == 3)
          xphat_mat = repmat(x,dist.length,1);
        elseif isnumeric(x) && all(numel(size(x) == [dist.length 3]))
          xphat_mat = x;        
        end
      elseif dim == 2 % 2D: projection to plane
        %phig = linsp;
        if isa(x,'TSeries') && isa(varargin{1},'TSeries')
          y = varargin{1}; varargin = varargin(2:end); % assume other coordinate for perpendicular plane is given after and in same format
          xphat_mat = x.resample(obj).norm.data;          
          yphat_mat = y.resample(obj).norm.data;
        elseif isnumeric(x) && numel(size(x) == 3)
          y = varargin{1}; varargin = varargin(2:end); % assume other coordinate for perpendicular plane is given after and in same format
          xphat_mat = repmat(x,dist.length,1);
          yphat_mat = repmat(y,dist.length,1);
        elseif isnumeric(x) && all(numel(size(x) == [dist.length 3]))
          y = varargin{1}; varargin = varargin(2:end); % assume other coordinate for perpendicular plane is given after and in same format
          xphat_mat = x;
          yphat_mat = y;
        else
          error('Can''t recognize second vector for the projection plane, ''y'': PDist.reduce(''2D'',x,y,...)')
        end          
        zphat_mat = cross(xphat_mat,yphat_mat); % it's x and z that are used as input to irf_int_sph_dist
        
        nargs = nargs - 1;
        args = args(2:end);
        
        % Set default projection grid, can be overriden by given input 'phig'
        nAzg = 32; 
        dPhig = 2*pi/nAzg;        
        phig = linspace(0,2*pi-dPhig,nAzg)+dPhig/2; % centers
      end
      % make input distribution to SI units, s^3/m^6
      dist = dist.convertto('s^3/m^6');
            
      
      %% Check for input flags
      % Default options and values
      doTint = 0;      
      nMC = 100; % number of Monte Carlo iterations
      vint = [-Inf,Inf];
      aint = [-180,180]; % azimuthal interval
      vgInput = 0;
      weight = 'none';      
      correct4scpot = 0;
      isDes = 1;
      if strcmp(dist.species,'electrons'); isDes = 1; else, isDes = 0; end
      
      ancillary_data = {};
      
      have_options = nargs > 1;
      while have_options
        switch(lower(args{1}))
          case {'t','tint','time'} % time
            l = 2;
            tint = args{2};
            doTint = 1;
          case 'nmc' % number of Monte Carlo iterations
            l = 2;
            nMC = args{2};
            ancillary_data{end+1} = 'nMC';
            ancillary_data{end+1} = nMC;
          case 'vint' % limit on transverse velocity (like a cylinder) [km/s]
            l = 2;
            vint = args{2};
          case 'aint'
            l = 2;
            aint = args{2};
          case 'phig'
            l = 2;
            phig = args{2};
          case 'vg' % define velocity grid
            l = 2;
            vgInput = 1;
            vg = args{2}*1e3;
          case 'weight' % how data is weighted
            l = 2;
            weight = args{2};
            ancillary_data{end+1} = 'weight';
            ancillary_data{end+1} = weight;    
          case 'scpot'
            l = 2;
            scpot = args{2};
            ancillary_data{end+1} = 'scpot';
            ancillary_data{end+1} = scpot; 
            correct4scpot = 1;
        end
        args = args((l+1):end);
        if isempty(args), break, end
      end
      
      % set vint ancillary data      
      ancillary_data{end+1} = 'vint';
      ancillary_data{end+1} = vint;
      ancillary_data{end+1} = 'vint_units';
      ancillary_data{end+1} = 'km/s';
      
      %% Get angles and velocities for spherical instrument grid, set projection
      %  grid and perform projection
      emat = double(dist.depend{1});
      if correct4scpot
        scpot = scpot.tlim(dist.time).resample(dist.time);
        scpot_mat = repmat(scpot.data, size(emat(1,:)));
        [it_below_scpot,ie_below_scpot] = find(emat < scpot_mat);
        %dist.data(it_below_scpot,ie_below_scpot,:,:) = 0;
        %dist.data(:,1:6,:,:) = 0;
        %emat = emat - scpot_mat;                
        
        % must also remove all tabler energies below zero
        
        %ind_below_scpot = find(emat<scpot_mat);
      end
      
      u = irf_units;

      if isDes == 1; M = u.me; else; M = u.mp; end

      if doTint % get time indicies
        if length(tint) == 1 % single time
          it = interp1(dist.time.epochUnix,1:length(dist.time),tint.epochUnix,'nearest');
        else % time interval
          it1 = interp1(dist.time.epochUnix,1:length(dist.time),tint(1).epochUnix,'nearest');
          it2 = interp1(dist.time.epochUnix,1:length(dist.time),tint(2).epochUnix,'nearest');
          it = it1:it2;
        end
      else % use entire PDist
        it = 1:dist.length;
      end
    
      nt = length(it);
      if ~nt % nt = 0
        error('Empty time array. Please verify the time(s) given.')
      end
      % try to make initialization and scPot correction outside time-loop
      
      % loop to get projection
      disp('Integrating distribution')
      fprintf('it = %4.0f/%4.0f\n',0,nt) % display progress
      for i = 1:nt
        if mod(i,1) == 0, fprintf([repmat('\b', 1, 10) '%4.0f/%4.0f\n'],i,nt); end % display progress
        if dim == 1
          xphat = xphat_mat(i,:);
        elseif dim == 2
          xphat = xphat_mat(i,:); % corresponding to phi = 0 in 'phig'
          zphat = zphat_mat(i,:); % normal to the projection plane
        end
        %fprintf('%g %g %g',xphat)
        
        % 3d data matrix for time index it
        F3d = double(squeeze(double(dist.data(it(i),:,:,:)))); % s^3/m^6
        energy = emat(it(i),:);
        
        if correct4scpot
          if 0
            dist.data(it(i),emat(it(i),:) < 0,:,:) = 0; %#ok<UNRCH>
          else
            ie_below_scpot = find(abs(emat(it(i),:)-scpot_mat(it(i),:))); % energy channel below 
            ie_below_scpot = find(abs(emat(it(i),:)-scpot_mat(it(i),:)) == min(abs(emat(it(i),:)-scpot_mat(it(i),:)))); % closest energy channel
            remove_extra_ind = 0; % for margin, remove extra energy channels
            F3d(it(i),1:(max(ie_below_scpot) + remove_extra_ind),:,:) = 0; 
            %disp(sprintf('%8.1g ',energy))
            energy = energy-scpot_mat(it(i),:);
            %disp(sprintf('%8.1g ',energy))
            energy(energy<0) = 0;
            %disp(sprintf('%8.1g ',energy))
          end
        end
            
        v = sqrt(2*energy*u.e/M); % m/s       

        if 0%length(v) ~= 32 % shopuld be made possible for general number, e.g. 64 (dist.e64)
            error('something went wrong') %#ok<UNRCH>
        end

        % azimuthal angle
        phi = double(dist.depend{2}(it(i),:)); % in degrees
        %phi = phi+180;
        %phi(phi>360) = phi(phi>360)-360;
        phi = phi-180;
        phi = phi*pi/180; % in radians

        if length(phi) ~= 32
            error('something went wrong')
        end

        % elevation angle
        th = double(dist.depend{3}); % polar angle in degrees
        th = th-90; % elevation angle in degrees
        th = th*pi/180; % in radi ans

        if length(th) ~= 16
            error('something went wrong')
        end

        % Set projection grid after the first distribution function
        % bin centers
        if ~vgInput
          if dim == 1
            vg = [-fliplr(v),v];
          elseif dim == 2
            vg = v;
          end
        end
        if i == 1            
            % initiate projected f
            if dim == 1
              Fg = zeros(length(it),length(vg));
              vel = zeros(length(it),2);
            elseif dim == 2
              Fg = zeros(length(it),length(vg),length(vg));
              vel = zeros(length(it),2);
            end
            dens = zeros(length(it),1);
        end
        % perform projection
        if dim == 1 
          % v, phi, th corresponds to the bins of F3d
          tmpst = irf_int_sph_dist(F3d,v,phi,th,vg,'x',xphat,'nMC',nMC,'vzint',vint*1e3,'aint',aint,'weight',weight);
          all_vg(i,:) = vg;
        elseif dim == 2
          %tmpst = irf_int_sph_dist_mod(F3d,v,phi,th,vg,'x',xphat,'z',zphat,'phig',phig,'nMC',nMC,'vzint',vint*1e3,'weight',weight);
          tmpst = irf_int_sph_dist(F3d,v,phi,th,vg,'x',xphat,'z',zphat,'phig',phig,'nMC',nMC,'vzint',vint*1e3,'weight',weight);
          all_vx(i,:,:) = tmpst.vx;
          all_vy(i,:,:) = tmpst.vy;
          all_vx_edges(i,:,:) = tmpst.vx_edges;
          all_vy_edges(i,:,:) = tmpst.vy_edges;
        end
        
        if dim == 1
          Fg(i,:,:) = tmpst.F;
        elseif dim == 2
          Fg(i,:,:) = tmpst.F_using_edges;
        end
         dens(i) = tmpst.dens;
         vel(i,:) = tmpst.vel; % dimension of projection, 1D if projection onto line, 2D if projection onto plane
         
      end
      % vg is m/s, transform to km/s
      if dim == 1
        PD = PDist(dist.time(it),Fg,'line (reduced)',all_vg*1e-3);      
      elseif dim == 2
        Fg_tmp = Fg(:,:,:);
        all_vx_tmp = permute(all_vx(:,:,:),[1 2 3])*1e-3;
        all_vy_tmp = permute(all_vy(:,:,:),[1 2 3])*1e-3;
        all_vx_edges_tmp = permute(all_vx_edges(:,:,:),[1 2 3])*1e-3;
        all_vy_edges_tmp = permute(all_vy_edges(:,:,:),[1 2 3])*1e-3;
        PD = PDist(dist.time(it),Fg_tmp,'plane (reduced)',all_vx_tmp,all_vy_tmp);
        PD.ancillary.vx_edges = all_vx_edges_tmp;
        PD.ancillary.vy_edges = all_vy_edges_tmp;
      end
      PD.species = dist.species;
      PD.userData = dist.userData;
     
      if dim == 1      
        PD.units = 's/m^4';
        PD.ancillary.projection_direction = xphat_mat(it,:);
        PD.ancillary
      elseif dim == 2
        PD.units = 's^2/m^5';
        PD.ancillary.projection_dir_1 = xphat_mat(it,:);
        PD.ancillary.projection_dir_2 = yphat_mat(it,:);
        PD.ancillary.projection_axis = zphat_mat(it,:);
      end
      
      while ~isempty(ancillary_data)
        PD.ancillary.(ancillary_data{1}) = ancillary_data{2};
        ancillary_data(1:2) = [];
      end
        
      
      % Must add xphat to ancillary data!
      
    end
    function varargout = surf(varargin)
      % Surface/pcolor plot of PDist of type 'plane (reduced)'.
      %   PDist.surf(...)
      %   h_surf = PDist.surf(...)
      %   [h_surf,h_axis] = PDist.surf(...)
      %   [h_surf,h_axis,h_all] = PDist.surf(...)
      %
      %   Output:
      %      h_all - structure with all the used handles, useful for 
      %              changing LineWidth of contour or background color of
      %              'printinfo' text object
      %      h_all = 
      %         struct with fields:
      % 
      %               Axes: [1x1 Axes]
      %            Surface: [1x1 Surface]
      %            Contour: [1x1 Contour]
      %           Colorbar: [1x1 ColorBar]
      %           Infotext: [1x1 Text]
      %    
      %   Input:
      %     'tint'/tint  - if tint.length = 1, chooses closest distribution
      %                    if tint.length > 1, takes the average of 
      %                       everything within this interval, so be
      %                       cautious if a varying projection plane/axis
      %                       is used            
      %     'nolog10' - does not take log10 of data, default is to take log10
      %     'contour'/contour_levels - contour levels drawn in black, if
      %             log is used, use the power, as shown in the plot, i.e.
      %             contour_levels = [-8 -5] for contour at log10(f) = [-8 -5]      
      %     '10^3 km/s' - plots x and y axis velocities in 10^3 km/s, this
      %             is default for electrons
      %     'km/s' - plots x and y axis velocities in km/s, this is
      %             default for ions
      %     'printinfo' - prints info - number of distributions used, time
      %             or time interval, vint used for distribution 
      %             integration (see PDist.reduce), v1/v2 directions
      %     'flim'/flim - sets values outside of this range to NaN -
      %             default is [0 Inf]
      %
      %     Example:
      %       tint = irf.tint('2017-07-06T13:53:50.00Z',25);
      %       tint_plot = irf.tint('2017-07-06T13:54:00.00Z',5);
      %       Eint = [000 40000];
      %       vint = [-Inf Inf];
      %       iDist = iPDist1.tlim(tint).elim(eint);
      %       iLine = dmpaB1.resample(iDist).norm;
      %       iPlane1 = iLine.cross(irf.ts_vec_xyz(iLine.time,repmat([1 0 0],iLine.length,1)));
      %       iPlane2 = iLine.cross(iPlane1);      
      %       if2D = iDist.reduce('2D',iPlane1,iPlane2,'vint',vint); % reduced distribution perp to B
      %       h1 = subplot(3,1,1);
      %       if2D.surf(h1,'printinfo','tint',tint_plot(1))
      %       h2 = subplot(3,1,2);
      %       if2D.surf(h2,'printinfo','tint',tint_plot(2))
      %       h3 = subplot(3,1,3);
      %       if2D.surf(h3,'printinfo','tint',tint_plot)
      %       
      %   See also mms.plot_int_projection, PDist.reduce
      
      %[ax,args,nargs] = axescheck(varargin{:});
      % axescheck only checks if the first argument is an axis handle, but
      % PDist.surf can be called as both PDist.surf(ax) and surf(ax,PDist),
      % where the first option has the axis as the second argument while
      % the second option has the axis as the first argument. Instead we
      % search untilwe find the first axis handle.
      have_axes = 0;
      orig_args = varargin;
      iArg = 0;
      while ~have_axes
        iArg = iArg + 1; 
        if ((isscalar(orig_args{iArg}) && isgraphics(orig_args{iArg},'axes')) ...
          || isa(orig_args{iArg},'matlab.graphics.axis.AbstractAxes') || isa(orig_args{iArg},'matlab.ui.control.UIAxes'))
          ax = orig_args{iArg};
          axes(ax); % set current axis to ax, needed for text(), which does not accept ax as input
          have_axes = 1;
          args = orig_args;
          args(iArg) = [];
          nargs = numel(args);
        end
      end
            
      if isempty(ax); ax = gca; end
      all_handles.Axes = ax;
            
      % Make sure first non axes-handle input is PDist of the right type.
      if isa(args{1},'PDist') && any(strcmp(args{1}.type,{'plane (reduced)','plane (slice)'}))
        dist_orig = args{1};
      else
        error('First input that is not an axes handle must be a PDist of type ''plane (reduced)'' or ''plane (slice)'', see PDist.reduce.')
      end
      args = args(2:end);
      nargs = nargs - 1;
      
      dist = dist_orig;
            
      % default plotting parameters
      doLog10 = 1;
      doColorbar = 1;
      doAxisLabels = 1;
      doPrintInfo = 0;
      doContour = 0;
      doCircles = 0;      
      doFLim = 1; flim = [0 Inf];
      
      if strcmp(dist.species,'electrons') 
        v_scale = 1e-3;
        v_label_units = '10^3 km/s';
      elseif strcmp(dist.species,'ions') 
        v_scale = 1;
        v_label_units = 'km/s';
      else
        error(sprintf('Species %s not supported',dist.species))
      end
      
      tId = 1:dist.length; % if tint is not given as input (check below) default is to include all the time indices of input distribution
      
      % check for input, try to keep it at a minimum, so that the
      % functionality is similar to Matlabs plot function, all the details
      % can then be fixed outside the function using ax.XLim, ax.YLim, 
      % ax.CLim, etc...
      if nargs > 0; have_options = 1; else have_options = 0; end
      while have_options
        l = 1;
        switch(lower(args{1}))   
          case 'tint'
            l = 2;
            notint = 0;
            tint = args{2};
            if tint.length == 1 % find closest time
              [~,tId] = min(abs(dist.time-tint));                         
            else % take everything within time interval
              [tId,~] = dist.time.tlim(tint);                            
            end
          case 'vectors'
            l = 2;
            vectors = args{2};
            have_vectors = 1;              
          case 'scpot'
            l = 2;
            scpot = args{2};
            if isa(scpot,'TSeries')
              includescpot = 1;
              irf.log('notice','Spacecraft potential passed.')
            else
                includescpot = 0;
                irf.log('notice','scpot not recognized. Not using it.')
            end
          case 'nolog10'
            l = 1;
            doLog10 = 0;   
          case 'contour'
            l = 2;
            contour_levels = args{2};  
            doContour = 1;
%           case 'circles_origin'
%             l = 2;
%             v_levels = args{2};
%             doCirclesOrigin = 1;
%           case 'circles_drifting'
%             l = 2;
%             v_levels = args{2};
%             doCirclesDrifting = 1;
          case 'circles' % same as circles drifting
            l = 2;
            v_levels = args{2};
            doCircles = 1;  
          case '10^3 km/s'
            l = 1;
            v_scale = 1e-3;
            v_label_units = '10^3 km/s';
          case 'km/s'
            l = 1;
            v_scale = 1;
            v_label_units = 'km/s';
          case 'printinfo'
            doPrintInfo = 1;
          case 'flim'
            l = 2;
            doFLim = 1;
            flim = args{2};
        end
        args = args(l+1:end);  
        if isempty(args), break, end    
      end
      
      % due to Matlab functionality, we must explicitly call the overloaded
      % subsref (defined within this subclass), otherwise it will call the 
      % builtin function
      subs.type = '()';
      subs.subs = {tId};
      dist = dist_orig.subsref(subs);
      if (length(dist.time)<1); irf.log('warning','No data for given time interval.'); return; end
      
      % prepare data to be plotted
      plot_data = squeeze(mean(dist.data,1)); % average data over time indices
      if doFLim % put values outside given interval to NaN, default is [0 Inf]
        plot_data(plot_data<=flim(1)) = NaN;
        plot_data(plot_data>flim(2)) = NaN;
      end
      if doLog10 % take log10 of data
        plot_data = log10(plot_data);      
      end
      
      % main surface plot
      % NOTE, PCOLOR and SURF uses flipped dimensions of (x,y) and (z), but PDist.reduce does not, there we need to flip the dim of the data
      plot_x_edges = squeeze(irf.nanmean(dist.ancillary.vx_edges,1))*v_scale; % v_scale, default 1e-3 for electrons to put axes in 10^3 km/s
      plot_y_edges = squeeze(irf.nanmean(dist.ancillary.vy_edges,1))*v_scale;
      plot_z_edges = plot_x_edges*0;                  
      ax_surface = surf(ax,plot_x_edges,plot_y_edges,plot_z_edges,plot_data'); 
      all_handles.Surface = ax_surface;      
      view(ax,[0 0 1])
      ax.Box = 'on';
      
      if doContour
        hold(ax,'on')
        if numel(contour_levels) == 1
          contour_levels = [contour_levels contour_levels];
        end
        plot_x = squeeze(irf.nanmean(dist.depend{1},1))*v_scale; % v_scale, default 1e-3 for electrons to put axes in 10^3 km/s
        plot_y = squeeze(irf.nanmean(dist.depend{2},1))*v_scale;        
        [~,h_contour] = contour(ax,plot_x,plot_y,plot_data',contour_levels,'k');
        h_contour.LineWidth = 1.5;
        hold(ax,'off')
        all_handles.Contour = h_contour;
      end            
      if doColorbar
        hcb = colorbar('peer',ax);
        if doLog10
          hcb.YLabel.String = sprintf('log_{10} f (%s)',dist.units);
        else
          hcb.YLabel.String = sprintf('f (%s)',dist.units);
        end
        all_handles.Colorbar = hcb;
      end
      if doCircles
        hold(ax,'on')
        nAngles = 100; 
        angles = linspace(0,2*pi,nAngles);
        vx_drift = v_levels(:,1);
        vy_drift = v_levels(:,2);
        v_circle = v_levels(:,3);
        if numel(v_circle) == 1
          vx_levels = v_circle*cos(angles);
          vy_levels = v_circle*sin(angles);
          vx_drifts = vx_drift;
          vy_drifts = vy_drift;
        else
          vx_levels = tocolumn(v_circle)*sin(angles);
          vy_levels = tocolumn(v_circle)*cos(angles);
          vx_drifts = repmat(vx_drift',nAngles',1);
          vy_drifts = repmat(vy_drift',nAngles',1);
        end
        h_levels = plot(ax,vx_drifts+vx_levels',vy_drifts+vy_levels','k','LineWidth',1.5);                          
        hold(ax,'off')
        all_handles.Circles = h_levels;
      end
      if doAxisLabels
        ax.XLabel.String = sprintf('v (%s)',v_label_units);
        ax.YLabel.String = sprintf('v (%s)',v_label_units);
      end
      if doPrintInfo
        s1 = sprintf('v int (out-of-plane) = [%g %g] %s',dist.ancillary.vint(1),dist.ancillary.vint(1),dist.ancillary.vint_units);
        if numel(tId) == 1
          s2 = sprintf('time = %s',dist.time.utc);
        else
          s2 = sprintf('tint = %s - %s',dist.time(1).utc,dist.time(end).utc);
        end
        s3 = sprintf('nDist = %g',numel(tId));
        
        % check projection direction to see if they are varying or not
        n_proj_dirs_1 = size(unique(dist.ancillary.projection_dir_1,'rows'),1);
        n_proj_dirs_2 = size(unique(dist.ancillary.projection_dir_2,'rows'),1);
        if n_proj_dirs_1 == 1
          s4 = sprintf('v_{1,dir} = [%.2f %.2f %.2f]',dist.ancillary.projection_dir_1);
        else
          s4 = 'v_{1,dir} = varying';
        end
        if n_proj_dirs_2 == 1
          s5 = sprintf('v_{2,dir} = [%.2f %.2f %.2f]',dist.ancillary.projection_dir_2);
        else
          s5 = 'v_{2,dir} = varying';
        end
        
        h_text = text(ax.XLim(1),ax.YLim(2),sprintf('%s\n%s\n%s\n%s\n%s',s3,s2,s1,s4,s5));
        h_text.VerticalAlignment = 'top';
        h_text.HorizontalAlignment = 'left';
        h_text.BackgroundColor = [1 1 1];
        all_handles.Infotext = h_text;
      end
      
      % output
      if nargout == 0
        varargout = {};
      elseif nargout == 1 % return surface handle, this is how pcolor does it       
        varargout = {ax_surface};
      elseif nargout == 2 % return surface handle and axis handle
        varargout = {ax_surface,ax};
      elseif nargout == 3 % return all handles in a structure
        varargout = {ax_surface,ax,all_handles};
      end
    end
    function varargout = plot_plane(varargin)
      varargout_tmp = surf(varargin{:});
      varargout{1:nargout} = varargout_tmp;
      %varargout = varargout_tmp;
    end
    function varargout = plot_pad_polar(varargin)
      have_axes = 0;
      orig_args = varargin;
      iArg = 0;
      if nargin >1
        while ~have_axes
          iArg = iArg + 1; 
          if ((isscalar(orig_args{iArg}) && isgraphics(orig_args{iArg},'axes')) ...
            || isa(orig_args{iArg},'matlab.graphics.axis.AbstractAxes') || isa(orig_args{iArg},'matlab.ui.control.UIAxes'))
            ax = orig_args{iArg};
            axes(ax); % set current axis to ax, needed for text(), which does not accept ax as input
            have_axes = 1;
            args = orig_args;
            args(iArg) = [];
            nargs = numel(args);
          end
        end
      end
            
      if isempty(ax); ax = gca; end
      all_handles.Axes = ax;
            
      % Make sure first non axes-handle input is PDist of the right type.
      if isa(args{1},'PDist') && any(strcmp(args{1}.type,{'pitchangle'}))
        dist_orig = args{1};
      else
        error('First input that is not an axes handle must be a PDist of type ''pitchangle'', see PDist.pitchangles.')
      end
      args = args(2:end);
      nargs = nargs - 1;
      
      dist = dist_orig;
            
      % default plotting parameters
      doMirrorData = 1;
      doAxesV = 0; % default is to do energy
      doLog10 = 1;
      doLogAxes = 1;
      doColorbar = 1;
      doAxisLabels = 1;
      doPrintInfo = 0;
      doContour = 0;
      doCircles = 0;      
      doFLim = 1; flim = [0 Inf];
      
      if strcmp(dist.species,'electrons') 
        v_scale = 1e-3;
        v_label_units = '10^3 km/s';
      elseif strcmp(dist.species,'ions') 
        v_scale = 1;
        v_label_units = 'km/s';
      else
        error(sprintf('Species %s not supported',dist.species))
      end
      
      tId = 1:dist.length; % if tint is not given as input (check below) default is to include all the time indices of input distribution
      
      % check for input, try to keep it at a minimum, so that the
      % functionality is similar to Matlabs plot function, all the details
      % can then be fixed outside the function using ax.XLim, ax.YLim, 
      % ax.CLim, etc...
      if nargs > 0; have_options = 1; else have_options = 0; end
      while have_options
        l = 1;
        switch(lower(args{1}))   
          case 'tint'
            l = 2;
            tint = args{2};
            if tint.length == 1 % find closest time
              [~,tId] = min(abs(dist.time-tint));                         
            else % take everything within time interval
              [tId,~] = dist.time.tlim(tint);                            
            end         
          case 'scpot'
            l = 2;
            scpot = args{2};
            if isa(scpot,'TSeries')
              includescpot = 1;
              irf.log('notice','Spacecraft potential passed.')
            else
                includescpot = 0;
                irf.log('notice','scpot not recognized. Not using it.')
            end
          case 'nolog10'
            l = 1;
            doLog10 = 0;   
          case 'contour'
            l = 2;
            contour_levels = args{2};  
            doContour = 1;
%           case 'circles_origin'
%             l = 2;
%             v_levels = args{2};
%             doCirclesOrigin = 1;
%           case 'circles_drifting'
%             l = 2;
%             v_levels = args{2};
%             doCirclesDrifting = 1;
          case '10^3 km/s'
            l = 1;
            v_scale = 1e-3;
            v_label_units = '10^3 km/s';
            doAxesV = 1;
          case 'km/s'
            l = 1;
            v_scale = 1;
            v_label_units = 'km/s';
            doAxesV = 1;
          case 'printinfo'
            doPrintInfo = 1;
          case 'flim'
            l = 2;
            doFLim = 1;
            flim = args{2};
        end
        args = args(l+1:end);  
        if isempty(args), break, end    
      end
      
      % select time indices
      % due to Matlab functionality, we must explicitly call the overloaded
      % subsref (defined within this subclass), otherwise it will call the 
      % builtin function
      subs.type = '()';
      subs.subs = {tId};
      dist = dist_orig.subsref(subs);
      %dist = dist_orig(tId);
      if (length(dist.time)<1); irf.log('warning','No data for given time interval.'); return; end
      
      % prepare data to be plotted
      data = squeeze(mean(dist.data,1)); % average data over time indices
      data(data==0) = NaN; % put zero values to NaN
      if doFLim % put values outside given interval to NaN, default is [0 Inf]
        plot_data(data<=flim(1)) = NaN;
        plot_data(data>flim(2)) = NaN;
      end
      if doLog10 % take log10 of data
        data = log10(data);      
      end
      
      
      % main surface plot
      % NOTE, PCOLOR and SURF uses flipped dimensions of (x,y) and (z), but PDist.reduce does not, there we need to flip the dim of the data
      rho_edges = [dist.depend{1}-dist.ancillary.delta_energy_minus dist.depend{1}(:,end)+dist.ancillary.delta_energy_plus(:,end)];
      theta = dist.depend{2}; dtheta = theta(2)-theta(1); theta_edges = [theta(1)-dtheta/2 theta+dtheta/2];
      if doLogAxes, rho_edges = log10(rho_edges); end
      theta_edges = theta_edges + 90; % rotate data      
      [RHO,THETA] = meshgrid(rho_edges,theta_edges);            
      X = RHO.*cosd(THETA);
      Y = RHO.*sind(THETA);
      if doMirrorData % mirror data
        plot_X = [X; -flipdim(X(1:end-1,:),1)];
        plot_Y = [Y; flipdim(Y(1:end-1,:),1)];
        plot_data = [data flipdim(data,2)];
      end
      ax_surface = surf(ax,plot_X,plot_Y,plot_X*0,plot_data');
      
%       plot_x_edges = squeeze(irf.nanmean(dist.ancillary.vx_edges,1))*v_scale; % v_scale, default 1e-3 for electrons to put axes in 10^3 km/s
%       plot_y_edges = squeeze(irf.nanmean(dist.ancillary.vy_edges,1))*v_scale;
%       plot_z_edges = plot_x_edges*0;                  
%       ax_surface = surf(ax,plot_x_edges,plot_y_edges,plot_z_edges,plot_data'); 
      all_handles.Surface = ax_surface;      
      view(ax,[0 0 1])
      ax.Box = 'on';
      

      if doContour
        hold(ax,'on')
        if numel(contour_levels) == 1
          contour_levels = [contour_levels contour_levels];
        end
        plot_x = squeeze(irf.nanmean(dist.depend{1},1))*v_scale; % v_scale, default 1e-3 for electrons to put axes in 10^3 km/s
        plot_y = squeeze(irf.nanmean(dist.depend{2},1))*v_scale;        
        [~,h_contour] = contour(ax,plot_x,plot_y,plot_data',contour_levels,'k');
        h_contour.LineWidth = 1.5;
        hold(ax,'off')
        all_handles.Contour = h_contour;
      end            
      if doColorbar
        hcb = colorbar('peer',ax);
        if doLog10
          hcb.YLabel.String = sprintf('log_{10} f (%s)',dist.units);
        else
          hcb.YLabel.String = sprintf('f (%s)',dist.units);
        end
        all_handles.Colorbar = hcb;
      end
      if doCircles
        hold(ax,'on')
        nAngles = 100; 
        angles = linspace(0,2*pi,nAngles);
        vx_drift = v_levels(:,1);
        vy_drift = v_levels(:,2);
        v_circle = v_levels(:,3);
        if numel(v_circle) == 1
          vx_levels = v_circle*cos(angles);
          vy_levels = v_circle*sin(angles);
          vx_drifts = vx_drift;
          vy_drifts = vy_drift;
        else
          vx_levels = tocolumn(v_circle)*sin(angles);
          vy_levels = tocolumn(v_circle)*cos(angles);
          vx_drifts = repmat(vx_drift',nAngles',1);
          vy_drifts = repmat(vy_drift',nAngles',1);
        end
        h_levels = plot(ax,vx_drifts+vx_levels',vy_drifts+vy_levels','k','LineWidth',1.5);                          
        hold(ax,'off')
        all_handles.Circles = h_levels;
      end
      if doAxisLabels
        ax.XLabel.String = sprintf('v (%s)',v_label_units);
        ax.YLabel.String = sprintf('v (%s)',v_label_units);
      end
      if doPrintInfo
        s1 = '';%sprintf('v int (out-of-plane) = [%g %g] %s',dist.ancillary.vint(1),dist.ancillary.vint(1),dist.ancillary.vint_units);
        if numel(tId) == 1
          s2 = sprintf('time = %s',dist.time.utc);
        else
          s2 = sprintf('tint = %s - %s',dist.time(1).utc,dist.time(end).utc);
        end
        s3 = sprintf('nDist = %g',numel(tId));
        
        % check projection direction to see if they are varying or not
        %n_proj_dirs_1 = size(unique(dist.ancillary.projection_dir_1,'rows'),1);
        %n_proj_dirs_2 = size(unique(dist.ancillary.projection_dir_2,'rows'),1);
%         if n_proj_dirs_1 == 1
%           s4 = sprintf('v_{1,dir} = [%.2f %.2f %.2f]',dist.ancillary.projection_dir_1);
%         else
%           s4 = 'v_{1,dir} = varying';
%         end
%         if n_proj_dirs_2 == 1
%           s5 = sprintf('v_{2,dir} = [%.2f %.2f %.2f]',dist.ancillary.projection_dir_2);
%         else
%           s5 = 'v_{2,dir} = varying';
%         end
        
        h_text = text(ax.XLim(1),ax.YLim(2),sprintf('%s\n%s\n%s\n%s\n%s',s3,s2,'','',''));
        h_text.VerticalAlignment = 'top';
        h_text.HorizontalAlignment = 'left';
        h_text.BackgroundColor = [1 1 1];
        all_handles.Infotext = h_text;
      end
      
      %dbstack        
      %nargout
        
      % output
%       if nargout == 0
%         varargout = {};
%       elseif nargout == 1 % return surface handle, this is how pcolor does it       
%         varargout = {ax_surface};
%       elseif nargout == 2 % return surface handle and axis handle
%         varargout = {ax_surface,ax};
%       elseif nargout == 3 % return all handles in a structure
        varargout = {ax_surface,ax,all_handles};
%       end
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
      unique_etables = unique(obj.depend{1},'rows','stable');
      netables = size(unique_etables,1); % netables = 2 for older dta and 1 for newer data
        
      % find new elevels
      if numel(eint) == 2 % energy interval        
        elevels = [];
        for ietable = 1:netables % loop over 1 or 2 and saves all the unique indices, i.e. max range
          tmp_elevels = intersect(find(unique_etables(ietable,:)>eint(1)),find(unique_etables(ietable,:)<eint(2)));
          elevels = unique([elevels tmp_elevels]);
        end
        disp(['Effective eint = [' num2str(min(min(energy(:,elevels))),'%g') ' ' num2str(max(max(energy(:,elevels))),'%g') ']'])      
      else % pick closest energy level
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
      
      % update ancillary data
      if isempty(PD.ancillary) % if ancillary data is empty, this is just for backwards compatibility
        if netables == 1          
          PD.ancillary.esteptable = zeros(size(energy,1),0);
          PD.ancillary.energy0 = unique_etables(1,elevels);
          PD.ancillary.energy1 = unique_etables(1,elevels);
        elseif netables == 2
          [esteptable,~] = ismember(energy,unique_etables(2,:),'rows'); % where row is not unique_etables(2,:), '0' is returned
          PD.ancillary.esteptable = esteptable;
          PD.ancillary.energy0 = unique_etables(1,elevels);
          PD.ancillary.energy1 = unique_etables(2,elevels);
        end
      else % check what fields ancillary have, and update them
        if isfield(PD.ancillary, 'energy0'), PD.ancillary.energy0 = PD.ancillary.energy0(elevels);
          else PD.ancillary.energy0 = energy(1,elevels); end
        if isfield(PD.ancillary, 'energy1'), PD.ancillary.energy1 = PD.ancillary.energy1(elevels);
          else PD.ancillary.energy1 = energy(2,elevels); end
        if ~isfield(PD.ancillary, 'esteptable')
          [esteptable,~] = ismember(energy,PD.ancillary.energy1,'rows');
          PD.ancillary.esteptable = esteptable;
        end          
        if isfield(PD.ancillary,'energy'), PD.ancillary.energy = PD.ancillary.energy(:,elevels); end                
        if isfield(PD.ancillary,'delta_energy_minus'), PD.ancillary.delta_energy_minus = PD.ancillary.delta_energy_minus(:,elevels); end                
        if isfield(PD.ancillary,'delta_energy_plus'), PD.ancillary.delta_energy_plus = PD.ancillary.delta_energy_plus(:,elevels); end                
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
      % PDIST.SPECREC Prepares structure to be used with irf_spectrogram or irf_plot
      %   sr = PDIST.SPECREC(spectype)
      %     spectype - 'energy' - default for PDist.type 'omni'
      %                   Ex: PDist.omni.SPECREC
      %                       PDist.pitchangles(dmpaB).SPECREC('energy',pitchangles_indices_to_average_over)
      %                             (pitchangles_indices_to_average_over default is all the indices)                             
      %                'pitchangle'/'pa' - default for PDist.type 'pitchangle'
      %                   Ex: PDist.pitchangles(dmpaB).SPECREC
      %                'velocity'/'1D_velocity'/'velocity_1D' - default for PDist.type 'line (reduced)' 
      %                   Ex: PDist.reduce('1D',[1 0 0]).SPECREC
      %                       PDist.reduce('1D',[1 0 0]).SPECREC('velocity_1D','10^3 km/s')
      %                'v_f1D*v' - multiplies f by v^2
      %                'v_f1D*v^2' - multiplies f by v
      
      % supported spectrogram types
      set_default_spectype = 0;
      supported_spectypes = {'energy','pitchangle','pa','velocity','1D_velocity','velocity_1D','v_f1D*v','v_f1D*v^2'};      
      
      if isempty(varargin) % no spectype given
        set_default_spectype = 1;
      elseif ~isempty(varargin) && ~any(strcmp(supported_spectypes,varargin{1})) % given spectype not supported
        set_default_spectype = 1;
      end      

      if set_default_spectype
        switch obj.type
          case 'omni'
            spectype = 'energy';
            irf.log('warning',sprintf('Spectype not given, default spectype for distribution type ''%s'' is ''%s''.',obj.type, spectype));
          case 'line (reduced)'
            spectype = '1D_velocity';
            irf.log('warning',sprintf('Spectype not given, default spectype for distribution type ''%s'' is ''%s''.',obj.type, spectype));
          case 'pitchangle'
             if numel(obj.depend{2}) == 1
               spectype = 'energy';
             else
               spectype = 'pitchangle';
             end
%            spectype = 'pitchangle';
            irf.log('warning',sprintf('Spectype not given, default spectype for distribution type ''%s'' is ''%s''.',obj.type, spectype));
          otherwise
            irf.log('warning',sprintf('Distribution type ''%s'' not supported. Using spectype ''energy'' for whatever backwards compatibility there might be. This option will be removed in a future version.',obj.type));
            spectype = 'energy';
        end  
      else
         spectype = varargin{1};  
         varargin = varargin(2:end); % remove from varargin
      end
      
      switch obj.units % set to p_label according to units of PDist
        case {'s^3/km^6','s^3/cm^6','s^3/m^6','s^2/km^5','s^2/cm^5','s^2/m^5','s^1/km^4','s^1/cm^4','s^1/m^4','s/km^4','s/cm^4','s/m^4'}
          spec.p_label = {'PSD',obj.units};
        case {'keV/(cm^2 s sr keV)'}
          spec.p_label = {'DEF',obj.units};
        case {'1/(cm^2 s sr keV)'}
          spec.p_label = {'PEF',obj.units};  
        otherwise
          spec.p_label = {obj.units};
      end
      
      switch spectype % make specrec structure
        case 'energy'
          switch obj.type
            case 'pitchangle'
              if ~isempty(varargin) % assume next argument is the pitchangle level/levels we want to average over
                iPA = varargin{2};              
              else
                iPA = 1:size(obj.depend{2},2);
              end
              irf.log('warning',['Averaging over pitch angles [' sprintf(' %g',obj.depend{2}(iPA)) ']']);
              spec.p = squeeze(double(nanmean(obj.data(:,:,iPA),3)));
            case 'omni'
              spec.p = double(obj.data);  
            otherwise
              error('Spectype ''%s'' not yet implemented for distribution type ''%s''.',spectype,obj.type);
          end
          spec.t = obj.time.epochUnix;          
          spec.f = single(obj.depend{1});
          spec.f_label = {['E_' obj.species(1) ' (eV)']};
        case {'pitchangle','pa'}
          if ~isempty(varargin) % assume next argument is the pitchangle level/levels we want to average over
            iE = varargin{1};              
          else
            iE = 1:size(obj.depend{1},2);
          end
          spec.t = obj.time.epochUnix;
          irf.log('warning',['Averaging over energy levels [' sprintf(' %g',obj.depend{1}(1,iE)) ']']);
          irf.log('warning',['Averaging over energy levels [' sprintf(' %g',obj.depend{1}(2,iE)) ']']);
          spec.p = double(squeeze(nanmean(obj.data(:,iE,:),2))); % nanmean over energies
          %spec.p_label = {'dEF',obj.units};
          spec.f = single(obj.depend{2});
          spec.f_label = {'\theta (deg.)'};
        case {'velocity','1D_velocity','velocity_1D'}
          if ~strcmp(obj.type_,'line (reduced)'); error('PDist must be projected unto a vector: type: ''line (reduced)'', see PDist.reduce.'); end      
          % check for additional argument given
          if ~isempty(varargin) && strcmp(varargin{1},'10^3 km/s') % make y (v) units in 10^3 km/s (because they often go up 10^4)
            spec.f = obj.depend{1}*1e-3;
            spec.f_label = {'v (10^3 km/s)'};
          else
            spec.f = obj.depend{1};
            spec.f_label = {'v (km/s)'};
          end
          spec.t = obj.time.epochUnix;
          spec.p = double(squeeze(obj.data));  
        case {'v_f1D*v'}
          if ~strcmp(obj.type_,'line (reduced)'); error('PDist must be projected unto a vector: type: ''line (reduced)'', see PDist.reduce.'); end      
          % check for additional argument given
          if ~isempty(varargin) && strcmp(varargin{1},'10^3 km/s') % make y (v) units in 10^3 km/s (because they often go up 10^4)
            spec.f = obj.depend{1}*1e-3;
            spec.f_label = {'v (10^3 km/s)'};            
          else
            spec.f = obj.depend{1};
            spec.f_label = {'v (km/s)'};
          end
          spec.t = obj.time.epochUnix;
         % vscale = 
          spec.p = double(squeeze(obj.data)).*obj.depend{1};
          spec.p_label{2} = [spec.p_label{2} '*km/s'];
          spec.p_label{1} = [spec.p_label{1} '*v'];
        case {'v_f1D*v^2'}
          if ~strcmp(obj.type_,'line (reduced)'); error('PDist must be projected unto a vector: type: ''line (reduced)'', see PDist.reduce.'); end      
          % check for additional argument given
          if ~isempty(varargin) && strcmp(varargin{1},'10^3 km/s') % make y (v) units in 10^3 km/s (because they often go up 10^4)
            spec.f = obj.depend{1}*1e-3;
            spec.f_label = {'v (10^3 km/s)'};            
          else
            spec.f = obj.depend{1};
            spec.f_label = {'v (km/s)'};
          end
          spec.t = obj.time.epochUnix;
          spec.p = double(squeeze(obj.data)).*obj.depend{1}.^2;
          spec.p_label{2} = [spec.p_label{2} '*(km/s)^2'];
          spec.p_label{1} = [spec.p_label{1} '*v^2'];
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
      old_energies = obj.depend{1};
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
      if size(obj.depend{1},2) == 64; irf_log(proc,'PDist already has 64 energy levels.'); end 
      
      if ~any([isfield(obj.ancillary,'energy0') isfield(obj.ancillary,'energy1') isfield(obj.ancillary,'esteptable')]) % construct energy0, energy1, and esteptable 
        esteptable = zeros(obj.length,1);
        [energies,~,esteptable] = unique(obj.depend{1},'rows'); % consider using legacy
        energy0 = obj.depend{1}(1,:);
        energy1 = obj.depend{1}(2,:);
      end
      
      [pdistr,phir,energyr] = mms.psd_rebin(obj,TSeries(obj.time,obj.depend{2}),obj.ancillary.energy0,obj.ancillary.energy1,TSeries(obj.time,obj.ancillary.esteptable));
      PD = obj.clone(pdistr.time,pdistr.data);      
      PD.depend{1} = repmat(energyr,PD.length,1);
      PD.ancillary.energy = PD.depend{1}; 
      PD.depend{2} = phir.data;  
      
      if isfield(PD.ancillary,'energy0')
        PD.ancillary.energy0 = PD.depend{1}(1,:);
        PD.ancillary.energy1 = PD.depend{1}(1,:);
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