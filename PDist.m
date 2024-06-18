classdef PDist < TSeries

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
      % PDIST Create PDIST object.
      %   Constructor method (PDIST.PDIST) for class PDist.
      %   Load and work with particle distributions from satellite missions
      %   such as MMS. For a complete set of methods type >> methods PDIST.
      %
      %   PD = PDIST(time,data,type,depend_var1,...,depend_varN);
      %     N is the dimension of the input data
      %     nt is the number of time steps
      %     time - time in EpochTT format
      %     data - matrix of data in format [nt, sz1, ..., szN]
      %     type - type of distribution: 'moms-tens0', 'moms-tens1',
      %            'moms-tens2', 'skymap', 'pitchangle', 'omni',
      %            'line (reduced)' (same as '1Dcart'), 'plane (reduced)',
      %            'plane (slice)', 'box' (same as '3Dcart')
      %     depend_var - dependent variables (should be as many depend_var
      %                  as dimensions of the data set excluding time), e.g.:
      %                  velocity (km/s),
      %                  energy (eV),
      %                  instrument azimuthal angle (deg),
      %                  instrument polar angle (deg),
      %                  pitchangle (deg)
      %
      %   Note: PDist objects are typically constructed using some
      %   construction function like mms.make_pdist, mms.get_data for
      %   skymaps, or different methods of PDIST, like PDIST.reduce or
      %   PDIST.omni
      %
      % Example:
      %   % Note: mms.db_list_files and mms.get_data requires you have
      %   % database initiated:
      %   % If MMS data are located in directory:
      %   % /path/to/your/data/mms1/...
      %   % /path/to/your/data/mms2/...
      %   % etc..
      %   % >> mms.db_init('local_file_db','/path/to/your/data');
      %
      %   % Create skymap distributions from MMS data using mms.get_data
      %   tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
      %   mms_id = 1;
      %   iPDist1 = mms.get_data('PDi_fpi_brst_l2',tint,mms_id);
      %
      %   % Create skymap distributions from MMS data using mms.make_pdist
      %   tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
      %   list_files = mms.db_list_files('mms1_fpi_brst_l2_dis-dist',tint);
      %   ifile = 1; % in case there are several files, chose one
      %   filepath = [list_files(ifile).path '/' list_files(ifile).name];
      %   % It is of course possible to also directly type the path to your file.
      %   [iPDist1,iPDistErr1] = mms.make_pdist(filepath);
      %
      %   % Create phony 3D distribution on cartesian grid
      %   m = 9.1094e-31; % electron mass
      %   n = 1*1e6; % 1/m3
      %   vd = 1000; % m/s
      %   T = 1000; % eV
      %   vt = @(T) sqrt(2*units.kB*T/m); % m/s, thermal speed from temperature
      %   % Maxwellian distribution
      %   f = @(vx,vy,vz,T,n,vd) n./((pi)^(3/2).*vt(T).^3).*exp(-(vx-vd).^2./vt(T).^2-(vy).^2./vt(T).^2-(vz).^2./vt(T).^2);
      %   % Set up velocity grid
      %   nvx = 50; nvy = 50; nvz = 50;
      %   vx = linspace(-50,50,nvx); % km/s, depend_var1
      %   vy = linspace(-50,50,nvy); % km/s, depend_var2
      %   vz = linspace(-50,50,nvz); % km/s, depend_var3
      %   t = 0:1; % s, depend_var0
      %   [~,VX,VY,VZ] = ndgrid(t,vx*1e-3,vy*1e-3,vz*1e-3);
      %   F = f(VX,VY,VZ,T,n,vd);
      %   time = EpochTT(t);
      %   % Create PDist object
      %   PD = PDIST(time,F,'3Dcart',vx,vy,vz);
      %
      % See also: TSeries, EpochTT, irf.ts_skymap, mms.get_data,
      % mms.make_pdist, PDIST.reduce, PDIST.omni, PDIST.pitchangles,
      % mms.db_init

      if nargin<2, error('2 inputs required'), end

      obj@TSeries(t,data,'to',0);

      args = varargin;
      if isa(args{1},'char'); obj.type_ = args{1}; args(1) = [];
      else, error('3rd input must specify distribution type')
      end

      % collect required data, depend
      switch obj.type_
        case {'moms-tens0'} % eg. density or scalar temperature partial moments
          obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'energy'};
        case {'moms-tens1'} % eg. velocitypartial moments
        case {'moms-tens2'} % eg. pressure or temperature partial moments
        case {'skymap'} % construct skymap distribution
          obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'energy'};
          obj.depend{2} = args{1}; args(1) = []; obj.representation{2} = {'phi'};
          obj.depend{3} = args{1}; args(1) = []; obj.representation{3} = {'theta'};
        case {'pitchangle'} % construct pitchangle distribution
          obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'energy'};
          obj.depend{2} = args{1}; args(1) = []; obj.representation{2} = {'pitchangle'};
        case {'omni'} % construct omni directional distribution
          obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'energy'};
        case {'line (reduced)','1Dcart'} % % construct 1D distribution, through integration over the other 2 dimensions
          obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'velocity'};
        case {'plane (reduced)','2Dcart'} % construct 2D distribution, either through integration or by taking a slice
          obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'velocity1'};
          obj.depend{2} = args{1}; args(1) = []; obj.representation{2} = {'velocity2'};
        case {'plane (slice)'} % construct 2D distribution, either through integration or by taking a slice
          obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'velocity1'};
          obj.depend{2} = args{1}; args(1) = []; obj.representation{2} = {'velocity2'};
        case {'box','3Dcart'}
          obj.depend{1} = args{1}; args(1) = []; obj.representation{1} = {'velocity1'};
          obj.depend{2} = args{1}; args(1) = []; obj.representation{2} = {'velocity2'};
          obj.depend{3} = args{1}; args(1) = []; obj.representation{3} = {'velocity3'};
        otherwise
          warning('Unknown distribution type')
      end

      % Enforce energy to be a timeseries (not much difference in data size anyways)
      obj = obj.enforce_depend_timeseries('energy');

      % Should check dimension of depends, and switch if they are wrong,
      % time should always be first index, and it can be 1 or obj.nt
      % This has only been partly implemented here...
      % Should be moved to a private method to make code more easily readable.
      size_data = size(obj.data);
      for idep = 1:numel(obj.depend)
        size_dep = size(obj.depend{idep});
        if not(size_dep(1) == 1)
          if size_dep(2) == 1
            obj.depend{idep} = obj.depend{idep}';
          end
        end
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
          if not(isempty(obj.ancillary))
            ancillary_fieldnames = fieldnames(obj.ancillary);
            new_ancillary_data = obj.ancillary;
            for iField = 1:numel(ancillary_fieldnames)
              field_data = getfield(obj.ancillary,ancillary_fieldnames{iField});
              if isnumeric(field_data) && size(field_data,1) == sizeData(1) % has the same number of rows as the PDist has time indices, assume each row corresponds to the same time index
                new_ancillary_data = setfield(new_ancillary_data,ancillary_fieldnames{iField},field_data(idxTmp{1},:,:,:,:,:,:)); % repeated :,:,:,:,:,:, used to support multidimensional data
              end
            end
            obj.ancillary = new_ancillary_data;
          end
          if numel(idx) > 1
            nargout_str = [];
            if nargout == 0 % dont give varargout
              obj = builtin('subsref',obj,idx(2:end));
            else
              for inout = 1:nargout % create [out1,out2,...outN] to get the correct number or nargout for rest of subsrefs (idx)
                c_eval('nargout_str = [nargout_str ''tmp_vout?,''];',inout)
              end
              nargout_str = ['[' nargout_str(1:end-1) ']'];
              varargout_str = ['{' nargout_str(2:end-1) '}'];
              eval(sprintf('%s = builtin(''subsref'',obj,idx(2:end));',nargout_str)) % disp(sprintf('%s = builtin(''subsref'',obj,idx(2:end));',nargout_str))
              eval(sprintf('varargout = %s;',varargout_str)); % varargout = {out1,out2,...outN}; % disp(sprintf('varargout = %s;',varargout_str))
            end
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
          obj.depend_{ii} = reshape(obj.depend_{ii}(idx,:),[numel(idx) sizeDepend(2:end)]);
        else
          error('Depend has wrong dimensions.')
        end
      end
      % on ancillary data
      if not(isempty(obj.ancillary))
        nameFields = fieldnames(obj.ancillary);
        nFields = numel(nameFields);
        for iField = 1:nFields
          eval(['sizeField = size(obj.ancillary.' nameFields{iField} ');'])
          if sizeField(1) == sizeData(1)
            eval(['obj.ancillary.' nameFields{iField} ' = reshape(obj.ancillary.' nameFields{iField} '(idx,:),[numel(idx) sizeField(2:end)]);'])
          end
        end
      end
    end
    function obj = resample_depend_ancillary(obj,NewTime,varargin)
      TsTmp = obj;
      tData = double(TsTmp.time.ttns - TsTmp.time.start.ttns)/10^9;
      dataTmp = double(TsTmp.data);
      newTimeTmp = double(NewTime.ttns - TsTmp.time.start.ttns)/10^9;

      %         % reshape data so it can be directly inserted into irf_resamp
      %         origDataSize = size(dataTmp);
      %         dataTmpReshaped = squeeze(reshape(dataTmp,[origDataSize(1) prod(origDataSize(2:end))]));
      %         newDataTmpReshaped = irf_resamp([tData dataTmpReshaped], newTimeTmp, varargin{:}); % resample
      %         newDataReshaped = squeeze(newDataTmpReshaped(:,2:end)); % take away time column
      %         newData = reshape(newDataReshaped,[length(newTimeTmp) origDataSize(2:end)]); % shape back to original dimensions

      % depend data
      sizeData = size(obj.data);
      nDepend = numel(obj.depend);
      for ii = 1:nDepend
        sizeDepend =  size(obj.depend{ii});
        if sizeDepend(1) == 1 % same dependence for all times
          obj.depend_{ii} = obj.depend{ii};
        elseif sizeDepend(1) == TsTmp.length
          dataTmp = obj.depend{ii};
          origDataSize = size(dataTmp);
          dataTmpReshaped = squeeze(reshape(dataTmp,[origDataSize(1) prod(origDataSize(2:end))]));
          newDataTmpReshaped = irf_resamp([tData dataTmpReshaped], newTimeTmp, varargin{:}); % resample
          newDataReshaped = squeeze(newDataTmpReshaped(:,2:end)); % take away time column
          newData = reshape(newDataReshaped,[length(newTimeTmp) origDataSize(2:end)]); % shape back to original dimensions

          obj.depend_{ii} = newData;
        else
          error('Depend has wrong dimensions.')
        end
      end

      % ancillary data
      if not(isempty(obj.ancillary))
        nameFields = fieldnames(obj.ancillary);
        nFields = numel(nameFields);
        for iField = 1:nFields
          eval(['sizeField = size(obj.ancillary.' nameFields{iField} ');'])
          if sizeField(1) == TsTmp.length
            old_ancillary = eval(['obj.ancillary.' nameFields{iField}]);

            % temporary fix for upsampling any non single or double data (esteptable!)
            if not(any([isa(old_ancillary,'double'),isa(old_ancillary,'single')]))
              old_ancillary = single(old_ancillary);
            end

            new_ancillary = irf_resamp([tData old_ancillary], newTimeTmp, varargin{:});
            eval(['obj.ancillary.' nameFields{iField} ' = new_ancillary(:,2:end);'])
          end
        end
      end
    end
    function obj = mtimes(obj,value)
      obj.data = obj.data*value;
    end
    function obj = times(obj,value)
      obj.data = obj.data.*value;
    end
    function obj = mdivide(obj,value)
      obj.data = obj.data/value;
    end
    function obj = divide(obj,obj2)
      obj.data = obj.data./obj2.data;
    end
    function [x,y,z] = xyz(obj,varargin)
      % PDIST.XYZ Get xyz coordinates of each detector bin.
      % PLEASE REPORT ERRORS.
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
      %     'scpot' - correct velocity for spacecraft potential
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

      units = irf_units;
      doScpot = 0;
      doReturnTSeries = 0;
      doSqueeze = 0;
      doRotation = 0;
      have_options = 0;
      doEdges = 0;

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
          case 'scpot'
            l = 2;
            doScpot = 1;
            scpot = args{2};
            args = args(l+1:end);
          case 'edges'
            doEdges = 1;
            args = args(l+1:end);
          otherwise
            irf.log('warning',sprintf('Input ''%s'' not recognized.',args{1}))
            args = args(l+1:end);
        end
        if isempty(args), break, end
      end

      switch obj.type
        case 'skymap'
          phi = TSeries(obj.time,obj.depend{1,2});
          azimuthal = phi.data*pi/180;

          theta = obj.depend{1,3};
          polar = repmat(theta*pi/180,obj.length,1);
          energy = obj.depend{1};

          if doEdges % get V at edges of bins
            azimuthal = [azimuthal(:,1) - 0.5*(azimuthal(:,2)-azimuthal(:,1)) azimuthal(:,1:end-1) + 0.5*diff(azimuthal,1,2) azimuthal(:,end) + 0.5*(azimuthal(:,end)-azimuthal(:,end-1))];
            polar = [polar(:,1) - 0.5*(polar(:,2)-polar(:,1)) polar(:,1:end-1) + 0.5*diff(polar,1,2) polar(:,end) + 0.5*(polar(:,end)-polar(:,end-1))];
            de_minus = obj.ancillary.delta_energy_minus;
            de_plus = obj.ancillary.delta_energy_plus;
            energy = [energy(:,1) - de_minus(:,1) energy + de_plus];
          end

          if doScpot
            if isa(scpot,'TSeries')
              scpot = scpot.resample(obj).data;
            elseif isa(scpot,'numeric') && numel(scpot) == 1
              scpot = repmat(scpot,obj.length,1);
            elseif isa(scpot,'numeric') && numel(scpot) == obj.length
              scpot = scpot;
            end
          else
            scpot = zeros(obj.length,1);
          end

          vx = NaN*obj.data;
          vy = NaN*obj.data;
          vz = NaN*obj.data;

          if doEdges
            sizedata = size(obj.data);
            sizedata(2:end) = sizedata(2:end) + 1;
            vx = nan(sizedata);
            vy = nan(sizedata);
            vz = nan(sizedata);
          end

          sp = obj.species;
          switch sp
            case 'ions'
              m = units.mp;
              mult = 1;
            case 'electrons'
              m = units.me;
              mult = -1;
          end

          for ii = 1:length(obj.time)
            % Adjust for spacecraft potential
            energy_tmp = energy(ii,:) + mult*scpot(ii);
            energy_tmp(energy_tmp<0) = 0; % if scpot is not used, energy_tmp = energy and nothing is changed

            % Energy -> speed
            velocity = sqrt((energy_tmp)*units.eV*2/m)/1000; % km/s

            % ndgrid of spherical coordinates
            [VEL,AZ,POL] = ndgrid(velocity,azimuthal(ii,:),polar(ii,:));

            % From spherical to cartesian coordinates
            % '-' because the data shows which direction the particles were coming from
            VX = -VEL.*sin(POL).*cos(AZ);
            VY = -VEL.*sin(POL).*sin(AZ);
            VZ = -VEL.*cos(POL);

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
        case 'pitch'
      end
    end
    function PD = d3v(obj,varargin)
      % Calculate phase space volume of FPI bins.
      %
      % Get partial density by doing: dn = pdist*pdist.d3v;
      %
      %   Options:
      %     'scpot',scpot - Corrects for spacecraft potential. For better
      %                     accordance with FPI, multiply scpot with 1.2,
      %                     see mms.psd_moments.
      %     'mat' - returns matrix (nt x nE x nAz x nPol) with phase space
      %             volume

      % Default return is f_fpi*d3v, i.e. PDist multiplied with volume
      % corresponding to each bin, giving the units of density.

      units = irf_units;
      doScpot = 0;
      doReturnMat = 0;
      nargs = numel(varargin);
      have_options = 0;
      if nargs > 0, have_options = 1; args = varargin(:); end
      while have_options
        l = 0;
        switch(lower(args{1}))
          case 'scpot'
            scpot = varargin{2};
            doScpot = 1;
            l = 2;
            args = args(l+1:end);
          case 'mat'
            doReturnMat = 1;
            l = 1;
            args = args(l+1:end);
          otherwise
            l = 1;
            irf.log('warning',sprintf('Input ''%s'' not recognized.',args{1}))
            args = args(l+1:end);
        end
        if isempty(args), break, end
      end

      switch obj.units % check units and if they are supported
        case 's^3/cm^6' % m^3/s^3 = m^3/s^3 * cm^3/cm^3 = cm^3/s^3 * m^3/cm^3 = cm^3/s^3 * (10^-2)^3
          d3v_scale = 1/10^(-2*3);
          new_units = 'cm^3/s^3';
        case 's^3/m^6' % m^3/s^3 = m^3/s^3 * m^3/m^3 = m^3/s^3 * m^3/m^3 = m^3/s^3 * (10^0)^3
          d3v_scale = 1/10^0;
          new_units = 'm^3/s^3';
        case 's^3/km^6' % m^3/s^3 = m^3/s^3 * km^3/km^3 = km^3/s^3 * m^3/km^3 = km^3/s^3 * (10^3)^3
          d3v_scale = 1/10^(3*3);
          new_units = 'km^3/s^3';
        otherwise
          error(sprintf('PDist.d3v not supported for %s',obj.units))
      end

      % Calculate velocity volume of FPI bin
      % int(sin(th)dth) -> x = -cos(th), dx = sin(th)dth -> int(dx) -> x = [-cos(th2) + cos(th1)] = [cos(th1) - cos(th1)]
      bin_edge_polar = [obj.depend{3} - 0.5*mean(diff(obj.depend{3})) obj.depend{3}(end) + 0.5*mean(diff(obj.depend{3}))];
      d_polar = cosd(bin_edge_polar(1:(end-1))) - cosd(bin_edge_polar(2:end));
      d_polar_mat = zeros(size(obj.data));
      c_eval('d_polar_mat(:,:,:,?) = d_polar(?);',1:16)

      % int(dphi) -> phi
      bin_azim = obj.depend{2}(1,2) - obj.depend{2}(1,1);
      d_azim = bin_azim*pi/180;

      % int(v^2dv) -> v^3/3
      if doScpot
        E_minus = (obj.depend{1} - obj.ancillary.delta_energy_minus) - repmat(scpot.data,1,size(obj.depend{1},2));
        E_plus = (obj.depend{1} + obj.ancillary.delta_energy_plus)   - repmat(scpot.data,1,size(obj.depend{1},2));
        E_minus(E_minus<0)= 0;
        E_plus(E_plus<0)= 0;
      else
        E_minus = (obj.depend{1} - obj.ancillary.delta_energy_minus);
        E_plus = (obj.depend{1} + obj.ancillary.delta_energy_plus);
      end
      v_minus = sqrt(2*units.e*E_minus/units.me); % m/s
      v_plus = sqrt(2*units.e*E_plus/units.me); % m/s
      d_vel = (v_plus.^3 - v_minus.^3)/3; % (m/s)^3
      d_vel_mat = repmat(d_vel,1,1,32,16);

      d3v = d_vel_mat.*d_azim.*d_polar_mat; % (m/s)^3

      if doReturnMat
        PD = d3v*d3v_scale;
      else
        PD = obj;
        PD.data = d3v*d3v_scale;
        PD.units = new_units;
        PD.name = 'd3v';
        PD.siConversion = num2str(str2num(PD.siConversion)/d3v_scale,'%e');
      end
    end
    function PD = solidangle(obj)
      % Solid angle of bins, can for example be used when working with
      % pitchangles, or fluxes (where units is flux/sr)
      %
      % The change in solid angle is only due to the changes in polar (or
      % pitch) angle, you therefore get all the unique values as follows:
      %
      %   squeeze(ePDist.solidangle.data(1,1,1,:))
      %   squeeze(ePDist(1).pitchangles(dmpaB1,15).solidangle.data(1,1,:))
      %
      % Total solid angle is 4*pi
      %
      %   sum(ePDist.solidangle.data(1,1,:))
      %   sum(ePDist(1).pitchangles(dmpaB1,15).solidangle.data(1,1,:))

      if strcmp(obj.type,'pitchangle')
        if isfield(obj.ancillary,'pitchangle_edges') && not(isempty(obj.ancillary.pitchangle_edges))
          bin_edge_polar = obj.ancillary.pitchangle_edges;
        else
          bin_edge_polar = [obj.depend{2} - 0.5*mean(diff(obj.depend{2})) obj.depend{2}(end) + 0.5*mean(diff(obj.depend{2}))];
        end
        d_polar = cosd(bin_edge_polar(1:(end-1))) - cosd(bin_edge_polar(2:end));
        d_polar_mat = zeros(size(obj.data));
        c_eval('d_polar_mat(:,:,?) = d_polar(?);',1:numel(obj.depend{2}))

        % int(dphi) -> phi
        d_azim = 2*pi; % all around

        sr_mat = d_polar_mat*d_azim;
      elseif strcmp(obj.type,'skymap')
        bin_edge_polar = [obj.depend{3} - 0.5*mean(diff(obj.depend{3})) obj.depend{3}(end) + 0.5*mean(diff(obj.depend{3}))];
        d_polar = cosd(bin_edge_polar(1:(end-1))) - cosd(bin_edge_polar(2:end));
        d_polar_mat = zeros(size(obj.data));
        c_eval('d_polar_mat(:,:,:,?) = d_polar(?);',1:16)

        % int(dphi) -> phi
        bin_azim = obj.depend{2}(1,2) - obj.depend{2}(1,1);
        d_azim = bin_azim*pi/180;

        sr_mat = d_azim.*d_polar_mat;
      else
        error(sprintf('PDist.type = %s not supported.',PDist.type))
      end

      PD = obj;
      PD.data = sr_mat;
      PD.units = 'sr';
      if isfield(PD.ancillary,'meanorsum'), PD.ancillary = rmfield(PD.ancillary,'meanorsum'); end
    end
    function PD = flux(obj,varargin)
      % Flux/sr [cm-2 s-1 sr-1] for skymaps and pitch angle distributions.
      %  j = int(fv d3v) = int(fv v^2dv sin(th)dth dphi)
      %    ~> (fv^4/4)*solidangle
      %
      %  Reduced distributions to be added.
      %
      %  To get flux in units [cm-2 s-1], multiply with solid angle:
      %   ePDist.flux.*ePDist.solidangle
      %
      %  FPI flux in EDI energy range
      %    dv_FPI_485 = 1760; % km/s
      %    dv_EDI_500 = 660; % km/s
      %    ePitch1 = ePDist1.pitchangles(dmpaB1,[168.5 180]); % antiparallel flux
      %    irf_plot(ePitch1.elim(500).flux*dv_EDI_500/dv_FPI_485)

      doScpot = 0;
      doPerSr = 1;
      doDiff = 0;

      nargs = numel(varargin);
      have_options = 0;
      if nargs > 0, have_options = 1; args = varargin(:); end

      while have_options
        l = 0;
        switch(lower(args{1}))
          case 'scpot'
            scpot = varargin{2};
            doScpot = 1;
            l = 2;
          case 'sr'
            doPerSr = varargin{2};
            l = 2;
          case 'diff'
            doDiff = 1;
            l = 1;
          otherwise
            l = 1;
            irf.log('warning',sprintf('Input ''%s'' not recognized.',args{1}))
        end
        args = args(l+1:end);
        if isempty(args), break, end
      end

      units = irf_units;

      if doDiff
        PD = obj;
        E_mat = repmat(PD.depend{1},1,1,size(PD.depend{2},2)); % eV
        E_mat_SI = E_mat*units.e;
        if 0
          %PD.data = PD.data*2.*E_mat_SI/PD.mass/PD.mass;
          PD.data = PD.data*2.*E_mat*units.e/PD.mass/PD.mass;
        else
          v_mat = sqrt(2*units.e*E_mat/obj.mass); % m/s
          v_mat = v_mat*1e2; % cm/s
          PD.data = PD.data.*v_mat.^2/PD.mass*1e-3;
        end
        PD.units = '1/(cm^2 s sr eV)';
        return
      end
      % int(v^3dv) -> v^4/4
      if doScpot
        E_minus = (obj.depend{1} - obj.ancillary.delta_energy_minus) - repmat(scpot.data,1,size(obj.depend{1},2));
        E_plus = (obj.depend{1} + obj.ancillary.delta_energy_plus)   - repmat(scpot.data,1,size(obj.depend{1},2));
        E_minus(E_minus<0)= 0;
        E_plus(E_plus<0)= 0;
      else
        E_minus = (obj.depend{1} - obj.ancillary.delta_energy_minus);
        E_plus = (obj.depend{1} + obj.ancillary.delta_energy_plus);
      end
      v_minus = sqrt(2*units.e*E_minus/units.me); % m/s
      v_plus = sqrt(2*units.e*E_plus/units.me); % m/s
      d_vel = (v_plus.^4 - v_minus.^4)/4; % (m/s)^3

      if strcmp(obj.type,'skymap')
        d_vel_mat = repmat(d_vel,1,1,size(obj.depend{2},2),size(obj.depend{3},2));
      elseif strcmp(obj.type,'pitchangle')
        d_vel_mat = repmat(d_vel,1,1,numel(obj.depend{2}));
      end

      if doPerSr
        vd3v = d_vel_mat;
        str_sr = '/sr';
      else
        solidangle = obj.solidangle;
        vd3v = d_vel_mat.*solidangle;
        str_sr = '';
      end

      old_units = obj.units;
      switch obj.units
        case 's^3/cm^6' % m^4/s^4 = m^4/s^4 * cm^4/cm^4 = cm^4/s^4 * m^4/cm^4 = cm^4/s^4 * (10^-2)^4
          d3v_scale = 1/10^(-2*4);
          new_units = sprintf('1/cm^2s%s',str_sr);
        case 's^3/m^6' % m^4/s^4 = m^4/s^4 * m^4/m^4 = m^4/s^4 * m^4/m^4 = m^4/s^4 * (10^0)^4
          d3v_scale = 1/10^0;
          new_units = sprintf('1/m^2s%s',str_sr);
        case 's^3/km^6' % m^4/s^4 = m^4/s^4 * km^4/km^4 = km^4/s^4 * m^4/km^4 = km^4/s^4 * (10^3)^4
          d3v_scale = 1/10^(3*4);
          new_units = sprintf('1/km^2s%s',str_sr);
      end

      PD = obj;
      PD.data = PD.data.*vd3v*d3v_scale;
      PD.units = new_units;
      PD.siConversion = num2str(str2num(PD.siConversion)/d3v_scale,'%e');
    end
    function PD = flux_red(obj,varargin)
      % Flux/sr [cm-2 s-1 sr-1], int(v^3dv) -> v^4/4, for skymaps and pitch angle distributions.
      %  Reduced distributions to be added.
      %
      %  To get flux in units [cm-2 s-1], multiply with solid angle:
      %   ePDist.flux.*ePDist.solidangle
      %
      %  FPI flux in EDI energy range
      %    dv_FPI_485 = 1760; % km/s
      %    dv_EDI_500 = 660; % km/s
      %    ePitch1 = ePDist1.pitchangles(dmpaB1,[168.5 180]); % antiparallel flux
      %    irf_plot(ePitch1.elim(500).flux*dv_EDI_500/dv_FPI_485)

      doScpot = 0;
      doPerSr = 1;

      nargs = numel(varargin);
      have_options = 0;
      if nargs > 0, have_options = 1; args = varargin(:); end

      while have_options
        l = 0;
        switch(lower(args{1}))
          otherwise
            l = 1;
            irf.log('warning',sprintf('Input ''%s'' not recognized.',args{1}))
            args = args(l+1:end);
        end
        if isempty(args), break, end
      end

      units = irf_units;

      v_minus = obj.ancillary.v_edges(1:end-1); % m/s
      v_plus = obj.ancillary.v_edges(2:end); % m/s
      d_vel = abs(v_plus.^2 - v_minus.^2)/2; % (m/s)^3
      d_vel_mat = repmat(d_vel,obj.length,1);

      str_sr = '';

      old_units = obj.units;
      switch obj.units
        case 's^3/cm^6' % m^4/s^4 = m^4/s^4 * cm^4/cm^4 = cm^4/s^4 * m^4/cm^4 = cm^4/s^4 * (10^-2)^4
          d3v_scale = 1/10^(-2*4);
          new_units = sprintf('1/cm^2s%s',str_sr);
        case 's^3/m^6' % m^4/s^4 = m^4/s^4 * m^4/m^4 = m^4/s^4 * m^4/m^4 = m^4/s^4 * (10^0)^4
          d3v_scale = 1/10^0;
          new_units = sprintf('1/m^2s%s',str_sr);
        case 's^3/km^6' % m^4/s^4 = m^4/s^4 * km^4/km^4 = km^4/s^4 * m^4/km^4 = km^4/s^4 * (10^3)^4
          d3v_scale = 1/10^(3*4);
          new_units = sprintf('1/km^2s%s',str_sr);
      end

      PD = obj;
      PD.data = PD.data.*d_vel_mat*1;
      PD.units = 's-1m-2';

      PD.data = PD.data*1e-4;
      PD.units = 's-1cm-2';
      PD.siConversion = '>1e4';%num2str(str2num(PD.siConversion)/d3v_scale,'%e');
    end
    function PD = reduce(obj,dim,x,varargin)
      %PDIST.REDUCE Reduces (integrates) 3D distribution to 1D (line).
      %   Example (1D):
      %     f1D = iPDist1.reduce('1D',dmpaB1,'vint',[0 10000]);
      %     irf_spectrogram(irf_panel('f1D'),f1D.specrec('velocity_1D'));
      %
      %   Example (2D):
      %     f2D = iPDist1.reduce('2D',[1 0 0],[0 1 0]);
      %     f2D(100).plot_plane
      %     [h_surf,h_axis,h_all] = f2D(100).plot_plane;
      %
      %   See more example uses in Example_MMS_reduced_ion_dist,
      %   Example_MMS_reduced_ele_dist, and Example_MMS_reduced_ele_dist_2D
      %
      %   Options:
      %     'nMC'    - number of Monte Carlo iterations used for integration,
      %                for default number see IRF_INT_SPH_DIST
      %     'base'   - set the base for the projection to cartesian 'cart'
      %                (default) or polar 'pol' (only valid for 2D planes)
      %     'vg'     - array with center values for the projection velocity
      %                grid in [km/s], determined by instrument if omitted
      %     'vg_edges' - array with edge values for the projection velocity
      %                grid in [km/s]
      %     'phig'   - array with center values for the projection
      %                azimuthal angle in [rad]
      %     'vint'   - set limits on the out-of-plane velocity to get
      %                cut-like distribution in 2D or a cylindrical shell
      %                in 1D in [km/s]
      %     'aint'   - angular limit in out-of-plane direction to make
      %                projection cut-like in 2D (only valid for 2D planes)
      %     'scpot'  - sets all values below scpot to zero and changes the
      %                energy correspondingly (only valid for electrons)
      %     'lowerelim' - sets all values below lowerelim to zero, does not
      %                change the energy. Can be single value, vector or
      %                Tseries, for example 2*scpot
      %     'weight' - how the number of MC iterations per bin is weighted,
      %                can be 'none' (default), 'lin' or 'log'
      %
      %
      %   The output is a PDist object with the reduced distribution where
      %   'data' is the integrated phase space density and 'depend'
      %   contains one (line) or two (plane) vectors of the velocity
      %   centers. The units of the velocity is [km/s].
      %
      % The integration itself is performed in irf_int_sph_dist.m
      %
      % See also: IRF_INT_SPH_DIST, PDIST.PLOT_PLANE, PDIST.SPECREC,
      % IRF_SPECTROGRAM

      %% Input
      [~,args,nargs] = axescheck(varargin{:});
      irf.log('warning','Please verify that you think the projection is done properly!');
      if isempty(obj); irf.log('warning','Empty input.'); return; else, dist = obj; end

      % Check to what dimension the distribution is to be reduced
      if any(strcmp(dim,{'1D','2D'}))
        dim = str2double(dim(1)); % input dim can either be '1D' or '2D'
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

        xphat_amplitude = sqrt(sum(xphat_mat.^2,2));
        if abs(mean(xphat_amplitude)-1) < 1e-2 && std(xphat_amplitude) > 1e-2 % make sure x are unit vectors,
          xphat_mat = xphat_mat./repmat(xphat_amplitude,1,3);
          irf.log('warning','|<x/|x|>-1| > 1e-2 or std(x/|x|) > 1e-2: x is recalculated as x = x/|x|.');
        end
      elseif dim == 2 % 2D: projection to plane
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

        % it's x and z that are used as input to irf_int_sph_dist
        % x and y are given, but might not be orthogonal
        % first make x and y unit vectors
        xphat_amplitude = sqrt(sum(xphat_mat.^2,2));
        yphat_amplitude = sqrt(sum(yphat_mat.^2,2));
        % These ifs are not really necessary, but could be there if one
        % wants to add some output saying that they were not put in
        % (inputted) as unit vectors. The definition of unit vectors is not
        % quite clear, due to tiny roundoff(?) errors
        if abs(mean(xphat_amplitude)-1) < 1e-2 && std(xphat_amplitude) > 1e-2 % make sure x are unit vectors,
          xphat_mat = xphat_mat./repmat(xphat_amplitude,1,3);
          irf.log('warning','|<x/|x|>-1| > 1e-2 or std(x/|x|) > 1e-2: x is recalculated as x = x/|x|.');
        end
        if abs(mean(yphat_amplitude)-1) < 1e-2 && std(yphat_amplitude) > 1e-2 % make sure y are unit vectors
          yphat_mat = yphat_mat./repmat(yphat_amplitude,1,3);
          irf.log('warning','|<y/|y|>-1| > 1e-2 or std(y/|y|) > 1e-2: y is recalculated as y = y/|y|.');
        end
        % make z orthogonal to x and y
        zphat_mat = cross(xphat_mat,yphat_mat,2);
        zphat_amplitude = sqrt(sum(zphat_mat.^2,2));
        zphat_mat = zphat_mat./repmat(zphat_amplitude,1,3);
        % make y orthogonal to z and x
        yphat_mat = cross(zphat_mat,xphat_mat,2);
        % check amplitude again, incase x and y were not orthogonal
        yphat_amplitude = sqrt(sum(yphat_mat.^2,2));
        if abs(mean(yphat_amplitude)-1) < 1e-2 && std(yphat_amplitude) > 1e-2  % make sure y are unit vectors
          yphat_mat = yphat_mat./repmat(yphat_amplitude,1,3);
          irf.log('warning','x and y were not orthogonal, y is recalculated as y = cross(cross(x,y),x)');
        end


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
      doLowerElim = 0;
      flag_dphi = 0;
      flag_dtheta = 0;
      nMC = 100; % number of Monte Carlo iterations
      vint = [-Inf,Inf];
      aint = [-180,180]; % azimuthal intherval
      vgInput = 0;
      vgInputEdges = 0;
      weight = 'none';
      correct4scpot = 0;
      base = 'cart'; % coordinate base, cart or pol

      if strcmp(dist.species,'electrons'); isDes = 1; else, isDes = 0; end

      ancillary_data = {};

      have_options = nargs > 1;
      while have_options
        switch(lower(args{1}))
          case {'t','tint','time'} % time (undocumented, can be removed?)
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
          case 'vg_edges'
            l = 2;
            vgInputEdges = 1;
            vg_edges = args{2}*1e3; % m/s
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
          case 'lowerelim'
            l = 2;
            lowerelim = args{2};
            ancillary_data{end+1} = 'lowerelim';
            ancillary_data{end+1} = lowerelim;
            doLowerElim = 1;
            if isnumeric(lowerelim) && numel(lowerelim) == 1
              lowerelim = repmat(lowerelim,dist.length,1);
            elseif isnumeric(lowerelim) && numel(lowerelim) == dist.length
              lowerlim = lowerelim;
            elseif isa(lowerelim,'TSeries')
              lowerelim = lowerelim.resample(dist).data;
            else
              error(sprintf('Can not recognize input for flag ''%s'' ',args{1}))
            end
          case 'base' %
            l = 2;
            base = args{2};
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
      units = irf_units;
      emat = double(dist.depend{1});
      if doLowerElim
        lowerelim_mat = repmat(lowerelim, size(emat(1,:)));
      end
      if correct4scpot
        %scpot = scpot.tlim(dist.time).resample(dist.time);
        scpot = scpot.resample(dist.time);
        scpot_mat = repmat(scpot.data, size(emat(1,:)));
      end
      if isDes == 1; M = units.me; else; M = units.mp; end
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

      if not(any([vgInput,vgInputEdges])) % prepare a single grid outside the time-loop
        %        emax = dist.ancillary.energy(1,end)+dist.ancillary.delta_energy_plus(1,end);
        emax = dist.depend{1}(1,end)+dist.ancillary.delta_energy_plus(1,end);
        vmax = units.c*sqrt(1-(emax*units.e/(M*units.c^2)+1).^(-2));
        nv = 100;
        vgcart_noinput = linspace(-vmax,vmax,nv);
        irf.log('warning',sprintf('No velocity grid specified, using a default vg = linspace(-vmax,vmax,%g), with vmax = %g km/s.',nv,vmax*1e-3));
      end
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

        if doLowerElim
          remove_extra_ind = 0; % for margin, remove extra energy channels
          ie_below_elim = find(abs(emat(it(i),:)-lowerelim_mat(it(i),:)) == min(abs(emat(it(i),:)-lowerelim_mat(it(i),:)))); % closest energy channel
          F3d(1:(max(ie_below_elim) + remove_extra_ind),:,:) = 0;
        end
        if correct4scpot
          if isfield(dist.ancillary,'delta_energy_minus') % remove all that satisfies E-Eminus<Vsc
            ie_below_scpot = find(emat(it(i),:)-dist.ancillary.delta_energy_minus(it(i),:)-scpot_mat(it(i),1)<0,1,'last');
            if 0 % disp energy channel that is removed, interferes with it = ... display
              disp(sprintf('Spacecraft potential = %g, Energy channel removed [E-Eminus,E,E+Eplus] = [%g,%g,%g]',...
                scpot_mat(it(i),1),...
                emat(it(i),ie_below_scpot)-dist.ancillary.delta_energy_minus(it(i),ie_below_scpot),...
                emat(it(i),ie_below_scpot),...
                emat(it(i),ie_below_scpot)+dist.ancillary.delta_energy_plus(it(i),ie_below_scpot)))
            end
          else
            ie_below_scpot = find(abs(emat(it(i),:)-scpot_mat(it(i),:)) == min(abs(emat(it(i),:)-scpot_mat(it(i),:)))); % closest energy channel
          end
          remove_extra_ind = 0; % for margin, remove extra energy channels
          F3d(1:(max(ie_below_scpot) + remove_extra_ind),:,:) = 0;
          %disp(sprintf('%8.1g ',energy))
          energy = energy-scpot_mat(it(i),:);
          %disp(sprintf('%8.1g ',energy))
          energy(energy<0) = 0;
          %disp(sprintf('%8.1g ',energy))
        end

        v = units.c*sqrt(1-(energy*units.e/(M*units.c^2)+1).^(-2)); % m/s

        % azimuthal angle
        if size(dist.depend{2},1)>1
          phi = double(dist.depend{2}(it(i),:)); % in degrees
        else % fast mode
          phi = double(dist.depend{2}); % in degrees
        end
        %phi = phi+180;
        %phi(phi>360) = phi(phi>360)-360;
        phi = phi-180;
        phi = phi*pi/180; % in radians

        % elevation angle
        th = double(dist.depend{3}); % polar angle in degrees
        th = th-90; % elevation angle in degrees
        th = th*pi/180; % in radians

        if isfield(dist.ancillary,'delta_phi_minus') && isfield(dist.ancillary,'delta_phi_plus')
          deltaphi = (dist.ancillary.delta_phi_plus+dist.ancillary.delta_phi_minus)*pi/180;
          if size(deltaphi,1) > size(deltaphi,2)
            deltaphi = deltaphi';
          end
          flag_dphi = 1;
        end
        if isfield(dist.ancillary,'delta_theta_minus') && isfield(dist.ancillary,'delta_theta_plus')
          deltatheta = (dist.ancillary.delta_theta_plus+dist.ancillary.delta_theta_minus)*pi/180;
          if size(deltatheta,1) > size(deltatheta,2)
            deltatheta = deltatheta';
          end
          flag_dtheta = 1;
        end

        % Set projection grid after the first distribution function
        % bin centers
        if vgInputEdges % redefine vg (which is vg_center)
          vg = vg_edges(1:end-1) + 0.5*diff(vg_edges);
        elseif vgInput
          vg = vg;
        else % define from instrument velocity bins
          if strcmp(base,'cart')
            vg = vgcart_noinput; % maybe just bypass this and go directly through input vg_edges?
          else
            if dim == 1
              vg = [-fliplr(v),v];
            elseif dim == 2
              vg = v;
            end
          end
        end

        % initiate projected f
        if i == 1
          if dim == 1
            Fg = zeros(length(it),length(vg));
            vel = zeros(length(it),1);
          elseif dim == 2 && strcmpi(base,'pol')
            Fg = zeros(length(it),length(phig),length(vg));
            vel = zeros(length(it),2);
          elseif dim == 2 && strcmpi(base,'cart')
            Fg = zeros(length(it),length(vg),length(vg));
            vel = zeros(length(it),2);
          end
          dens = zeros(length(it),1);
        end
        % perform projection
        if dim == 1 % 1D plane
          % v, phi, th corresponds to the bins of F3d
          if vgInputEdges
            if flag_dphi && flag_dtheta
              tmpst = irf_int_sph_dist(F3d,v,phi,th,vg,'x',xphat,'nMC',nMC,'vzint',vint*1e3,'aint',aint,'weight',weight,'vg_edges',vg_edges,'dphi',deltaphi,'dth',deltatheta);
            else
              tmpst = irf_int_sph_dist(F3d,v,phi,th,vg,'x',xphat,'nMC',nMC,'vzint',vint*1e3,'aint',aint,'weight',weight,'vg_edges',vg_edges);
            end
          else
            if flag_dphi && flag_dtheta
              tmpst = irf_int_sph_dist(F3d,v,phi,th,vg,'x',xphat,'nMC',nMC,'vzint',vint*1e3,'aint',aint,'weight',weight,'dphi',deltaphi,'dth',deltatheta);
            else
              tmpst = irf_int_sph_dist(F3d,v,phi,th,vg,'x',xphat,'nMC',nMC,'vzint',vint*1e3,'aint',aint,'weight',weight);
            end
          end
          all_vg(i,:) = tmpst.v; % normally vg, but if vg_edges is used, vg is overriden
          all_vg_edges(1,:) = tmpst.v_edges;
        elseif dim == 2
          %tmpst = irf_int_sph_dist_mod(F3d,v,phi,th,vg,'x',xphat,'z',zphat,'phig',phig,'nMC',nMC,'vzint',vint*1e3,'weight',weight);
          % is 'vg_edges' implemented for 2d?
          if flag_dphi && flag_dtheta
            tmpst = irf_int_sph_dist(F3d,v,phi,th,vg,'x',xphat,'z',zphat,'phig',phig,'nMC',nMC,'vzint',vint*1e3,'weight',weight,'base',base,'dphi',deltaphi,'dth',deltatheta);
          else
            tmpst = irf_int_sph_dist(F3d,v,phi,th,vg,'x',xphat,'z',zphat,'phig',phig,'nMC',nMC,'vzint',vint*1e3,'weight',weight,'base',base);
          end
          all_vx(i,:,:) = tmpst.vx;
          all_vy(i,:,:) = tmpst.vy;
          all_vx_edges(i,:,:) = tmpst.vx_edges;
          all_vy_edges(i,:,:) = tmpst.vy_edges;
        end

        % fix for special cases
        % dimension of projection, 1D if projection onto line, 2D if projection onto plane
        if dim == 1 || strcmpi(base,'cart')
          Fg(i,:,:) = tmpst.F;
        elseif dim == 2
          Fg(i,:,:) = tmpst.F_using_edges;
        end
        % set moments from reduced distribution (for debug)
        dens(i) = tmpst.dens;
        vel(i,:) = tmpst.vel;

      end

      % Construct PDist objects with reduced distribution
      % vg is m/s, transform to km/s
      if dim == 1
        PD = PDist(dist.time(it),Fg,'line (reduced)',all_vg*1e-3);
        PD.ancillary.v_edges = all_vg_edges;
      elseif dim == 2 && strcmpi(base,'pol')
        Fg_tmp = Fg(:,:,:);
        all_vx_tmp = permute(all_vx(:,:,1:end-1),[1 2 3])*1e-3;
        all_vy_tmp = permute(all_vy(:,:,1:end-1),[1 2 3])*1e-3;
        all_vx_edges_tmp = permute(all_vx_edges(:,:,:),[1 2 3])*1e-3;
        all_vy_edges_tmp = permute(all_vy_edges(:,:,:),[1 2 3])*1e-3;
        PD = PDist(dist.time(it),Fg_tmp,'plane (reduced)',all_vx_tmp,all_vy_tmp);
        PD.ancillary.vx_edges = all_vx_edges_tmp;
        PD.ancillary.vy_edges = all_vy_edges_tmp;
        PD.ancillary.base = 'pol';
      elseif dim == 2 && strcmpi(base,'cart')
        PD = PDist(dist.time(it),Fg,'plane (reduced)',squeeze(all_vx)*1e-3,squeeze(all_vx)*1e-3);
        PD.ancillary.vx_edges = all_vx_edges*1e-3;
        PD.ancillary.vy_edges = all_vx_edges*1e-3;
        PD.ancillary.base = 'cart';
      end
      PD.species = dist.species;
      PD.userData = dist.userData;
      PD.ancillary.v_units = 'km/s';
      PD.ancillary.energy = obj.depend{1};
      if isfield(obj.ancillary,'delta_energy_minus')
        PD.ancillary.delta_energy_minus = obj.ancillary.delta_energy_minus;
        PD.ancillary.delta_energy_plus = obj.ancillary.delta_energy_plus;
      end

      % set units and projection directions
      if dim == 1
        PD.units = 's/m^4';
        PD.ancillary.projection_direction = xphat_mat(it,:);
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

      if doLowerElim
        PD.ancillary.lowerelim = lowerelim_mat;
      end

    end % end of reduce function
    function PD = rebin(obj,base,grid,orient,varargin)
      % PDIST.REBIN Rebins energies of distribution function.
      %   Usage:
      %     PD = REBIN(dist,base,grid,orient);
      %       base - 'sph', 'cart' (cart quite slow, only typically use
      %                             single time step, say for example you
      %                             want to create a grid of test particles
      %                             based on observed distribution)
      %       orient - cartesian unit vectors [x,y,z] in DSL coordinates
      %       grid - 'sph', 'cart'
      %           if base is 'sph', grid should contain {energy,azimuthal_angle,polar_angle}
      %           if any is empty, it is kept as it is,for example, one
      %           can, only {energy,[],[]} implemented
      %           if base is 'cart' or 'cart_v', grid should be {vx,vy,vz}
      %
      %     Rebin to correspond to EDI energy interval.
      %     ePDist1_rebin_500 = ePDist1.rebin(''sph'',{[475 525],[],[]});',1);
      % See also IRF_INT_SPH_DIST

      %     if base is 'sph', grid should contain {energy,azimuthal_angle,polar_angle}
      %       if any is empty, it is kept as it is,for example, one can
      %       choose to only rebin in energies
      %     if base is 'cart' or 'cart_v', grid should be {vx,vy,vz}
      %     if base is 'cart_E', grid should be {Ex,Ey,Ez}
      %       v/E defines the edges of the bins
      %


      % Default values
      units = irf_units;
      doScpot = 0;
      nMC = 200;
      nt = obj.length;
      its = 1:nt;
      PD = [];

      % Check input
      nargs = numel(varargin);
      have_options = 0;
      if nargs > 0, have_options = 1; args = varargin(:); end
      while have_options
        l = 0;
        switch(lower(args{1}))
          case 'scpot'
            l = 2;
            scpot = args{2};
            doScpot = 1;
          otherwise
            l = 1;
            irf.log('warning',sprintf('Input ''%s'' not recognized.',args{1}))
        end
        args = args(l+1:end);
        if isempty(args), break, end
      end

      % Return output in same units as input
      % Input velocities are in km/s, use this v_scaling to transform all
      % velocities to the length scale units of the input f
      switch obj.units % check units and if they are supported
        case 's^3/cm^6'
          v_scale = 1e3/1e-2; % km/cm
        case 's^3/m^6'
          v_scale = 1e0/1e-2; % m/cm
        case 's^3/km^6'
          v_scale = 1e-2/1e-2; % cm/cm
        otherwise
          error(sprintf('PDist.d3v not supported for %s',obj.units))
      end

      % Start binning
      switch base
        case 'sph'
          old_az_num = size(obj.depend{2},2);
          old_pol_num = size(obj.depend{3},2);

          old_energy_minus = obj.depend{1} - obj.ancillary.delta_energy_minus;
          old_energy_plus = obj.depend{1} + obj.ancillary.delta_energy_plus;
          old_energy_num = size(old_energy_minus,2);

          old_v_minus = sqrt(2*units.e*old_energy_minus/units.me); % m/s
          old_v_plus = sqrt(2*units.e*old_energy_plus/units.me); % m/s
          old_v2dv = (old_v_plus.^3 - old_v_minus.^3)/3;

          old_data = obj.data;
          %old_dn = obj.d3v.data;
          old_d3v = obj.d3v.data; % volume of bin

          if not(isempty(grid{1}))
            new_energy_minus = grid{1}(1:end-1);
            new_energy_plus = grid{1}(2:end);
            new_energy_edges = unique([new_energy_minus,new_energy_plus]);

            new_v_minus = sqrt(2*units.e*new_energy_minus/units.me); % m/s
            new_v_plus = sqrt(2*units.e*new_energy_plus/units.me); % m/s
            new_v2dv = (new_v_plus.^3 - new_v_minus.^3)/3;
          end
          new_energy_num = numel(new_energy_minus);
          new_data = zeros(nt,new_energy_num,size(obj.depend{2},2),size(obj.depend{3},2));

          % loop through time
          nskip_erange = 0;
          nskip_fzero = 0;
          for it = its
            % loop through old instrument bins
            for ie = 1:old_energy_num
              if or(old_energy_minus(it,ie) > max(new_energy_plus),old_energy_plus(it,ie) < min(new_energy_minus))
                nskip_erange = nskip_erange + 1;
                %disp(sprintf('skipping: new_energy_channel max range = [%.0f, %.0f], old energy channel = [%.0f, %.0f]',min(new_energy_minus),max(new_energy_plus),old_energy_minus(it,ie),old_energy_plus(it,ie)))
                continue
              end
              for ipol = 1:old_pol_num
                for iaz = 1:old_az_num %
                  f_per_MC = old_data(it,ie,iaz,ipol)*old_v2dv(it,ie)/nMC;
                  %dn_per_MC = old_dn(it,ie,iaz,ipol)/nMC;
                  if f_per_MC == 0
                    nskip_fzero = nskip_fzero + 1;
                    continue
                  end
                  %try
                  energy_MC = rand(nMC,1)*(old_energy_plus(it,ie)-old_energy_minus(it,ie))+old_energy_minus(it,ie);

                  %f_per_MC = old_dn(it,ie,iaz,ipol)/nMC;
                  % divide particles into new bins
                  [N,EDGES] = histcounts(energy_MC,new_energy_edges);
                  new_data(it,:,iaz,ipol) = new_data(it,:,iaz,ipol) + N*f_per_MC./new_v2dv;
                  %catch
                  %  1;
                  %end
                end
              end
            end
          end
          disp(sprintf('nskip_erange = %g, nskip_fzero = %g',nskip_erange,nskip_fzero))
          PD = obj;
          PD.depend{1} = repmat(((new_energy_plus(:,:)+new_energy_minus(:,:))/2),nt,1);
          PD.ancillary.delta_energy_minus = abs(PD.depend{1}-new_energy_minus(:,:));
          PD.ancillary.delta_energy_plus = abs(PD.depend{1}-new_energy_plus(:,:));
          PD.ancillary.energy0 = PD.depend{1}(1,:);
          PD.ancillary.energy1 = PD.depend{1}(1,:);
          %new_dn = PD.d3v('mat');
          PD.data_ = new_data;%./new_dn;
        case 'cart'
          units = irf_units;

          % Get input
          new_vx_unit = orient(1,:);
          new_vy_unit = orient(2,:);
          new_vz_unit = orient(3,:);
          new_vbins_edges = grid;
          new_vx_edges = new_vbins_edges{1};
          new_vy_edges = new_vbins_edges{2};
          new_vz_edges = new_vbins_edges{3};
          nvx_new = numel(new_vx_edges);
          nvy_new = numel(new_vy_edges);
          nvz_new = numel(new_vz_edges);
          new_dvx = diff(new_vx_edges);
          new_dvy = diff(new_vy_edges);
          new_dvz = diff(new_vz_edges);
          new_vx_centers = new_vx_edges(1:end-1) + 0.5*new_dvx;
          new_vy_centers = new_vy_edges(1:end-1) + 0.5*new_dvy;
          new_vz_centers = new_vz_edges(1:end-1) + 0.5*new_dvz;
          [NDVX,NDVY,NDVZ] = ndgrid(new_dvx,new_dvy,new_dvz);
          new_vol = NDVX.*NDVY.*NDVZ*v_scale.^3; % phase space volume should be in same length units as inout f

          % Initialize new matrix for f
          dn_new = zeros(obj.length, nvx_new-1, nvy_new-1, nvz_new-1);
          f_new = zeros(obj.length, nvx_new-1, nvy_new-1, nvz_new-1);
          sizedata = [nvx_new-1, nvy_new-1, nvz_new-1];

          % Data from PDist in spherical coordinate system
          [old_vx,old_vy,old_vz] = obj.v; % v_xyz_DSL of original grid
          old_f = obj.data; % phase space density of each cell
          old_vol = obj.d3v.data; % Phase space volume of each cell
          old_dn = old_f.*old_vol;

          % Edges of energy bins, same for each time step
          energy_minus = obj.depend{1}(1,:) - obj.ancillary.delta_energy_minus;
          energy_plus = obj.depend{1}(1,:) + obj.ancillary.delta_energy_plus;
          denergy = energy_plus - energy_minus;
          energy_edges = [energy_minus energy_plus(end)];

          % Edges of polar angle bins, same for each time step
          polar_center = obj.depend{3};
          dpolar = polar_center(2) - polar_center(1);
          polar_edges = [polar_center(1:end-1)-0.5*dpolar polar_center(end)+0.5*dpolar];

          for it = 1:obj.length
            % Edges of azimuthal angle bins, changes for each time step
            azim_center = obj.depend{2};
            dazim = azim_center(2) - azim_center(1);
            azim_edges = [azim_center(1:end-1)-0.5*dazim azim_center(end)+0.5*dazim];

            % Create N particles within each bin that each recieve 1/N of
            % the phase space density. These are then rotated into the new
            % coordinate system and binned in the new grid.
            N = nMC;
            for iEnergy    = 1:size(old_f,2)
              for iAzim    = 1:size(old_f,3)
                for iPolar = 1:size(old_f,4)
                  % Skip bin if space spade density is zero
                  if old_dn(it,iEnergy,iAzim,iPolar) == 0
                    continue;
                  end

                  %                 if doScpot
                  %                     energy_minus = energy_minus ;
                  %                     energy_plus = energy_plus - scpot.data(it);
                  %                     energy_edges = energy_edges - scpot.data(it);
                  %                   end

                  tmp_energy = energy_edges(iEnergy) + denergy(iEnergy)*rand(N,1); % eV
                  if doScpot, tmp_energy = tmp_energy - scpot.data(it); end
                  tmp_azim   = azim_edges(iAzim)     + dazim*rand(N,1);   % deg
                  tmp_polar  = polar_edges(iPolar)   + dpolar*rand(N,1);  % deg
                  tmp_v = sqrt(tmp_energy*units.eV*2/units.me)/1000; % km/s

                  if doScpot % check if energy is negative, then skip
                    if iEnergy>6;
                      1;
                    end
                    tmp_azim(tmp_energy<0) = [];
                    tmp_polar(tmp_energy<0) = [];
                    tmp_v(tmp_energy<0) = [];
                    if isempty(tmp_v) continue; end
                  end

                  old_vx = -tmp_v.*sind(tmp_polar).*cosd(tmp_azim); % '-' because the data shows which direction the particles were coming from
                  old_vy = -tmp_v.*sind(tmp_polar).*sind(tmp_azim);
                  old_vz = -tmp_v.*cosd(tmp_polar);

                  % Rotate into new coordinate system
                  new_vx = old_vx*new_vx_unit(1) + old_vy*new_vx_unit(2) + old_vz*new_vx_unit(3);
                  new_vy = old_vx*new_vy_unit(1) + old_vy*new_vy_unit(2) + old_vz*new_vy_unit(3);
                  new_vz = old_vx*new_vz_unit(1) + old_vy*new_vz_unit(2) + old_vz*new_vz_unit(3);

                  % Assign particle density to each macro particle
                  tmp_dn = old_dn(it,iEnergy,iAzim,iPolar)/N;

                  % Bin into new grid
                  iVxg = discretize(new_vx,new_vx_edges);
                  iVyg = discretize(new_vy,new_vy_edges);
                  iVzg = discretize(new_vz,new_vz_edges);

                  % sizedata is like size(g_new) but excludes first index (time)
                  loc = sub2ind(sizedata,iVxg,iVyg,iVzg);
                  loc(isnan(loc)) = []; % values that fall outside of box becomes nan, remove these
                  hasdata = all(loc>0, 2);

                  %sum_dn_over_vol = accumarray(loc(hasdata,:), tmp_dn./new_vol(loc),[numel(f_new) 1]);
                  sum_dn = accumarray(loc(hasdata,:), tmp_dn,[numel(f_new) 1]);
                  dn_new(it,:,:,:) = dn_new(it,:,:,:) + reshape(sum_dn,[1 sizedata]);
                end
              end
            end

            if 0 % plot results.
              %%
              hca = subplot(3,1,1);
              surf(hca,new_vx_edges,new_vy_edges,zeros(nvx_new,nvy_new),log10(squeeze(sum(f_new,3))));
              view([0,0,1])
              hca = subplot(3,1,2);
              surf(hca,new_vx_edges,new_vz_edges,zeros(nvx_new,nvz_new),log10(squeeze(sum(f_new,2))));
              view([0,0,1])
              hca = subplot(3,1,3);
              surf(hca,new_vy_edges,new_vz_edges,zeros(nvy_new,nvz_new),log10(squeeze(sum(f_new,1))));
              view([0,0,1])
              1;
            end
          end
          % divide accumulated density (e.g. cm^{-3}) for each bin by the
          % velocity space volume of that bin (e.g. cm^{3}/s^{3})
          f_new = dn_new./reshape(repmat(new_vol,obj.length,1,1,1),[obj.length,nvx_new-1,nvz_new-1,nvz_new-1]);

          PD = PDist(obj.time,f_new,'3Dcart',new_vx_centers,new_vy_centers,new_vz_centers);
          PD.ancillary.vx_edges = new_vx_edges*1e-3; % km/s
          PD.ancillary.vy_edges = new_vy_edges*1e-3; % km/s
          PD.ancillary.vz_edges = new_vz_edges*1e-3; % km/s
          PD.ancillary.base = 'cart';
          PD.units = obj.units;
      end
    end
    function varargout = shift(pdist,v_nf,nMC,orient,sc,varargin)
      % PDIST.SHIFT  Rebin the distribution function to shifted reference
      %              frame and a rotated coordinate system.
      %
      %   PD = pdist.shift(v_nf,nMC,orient,sc,varargin)
      %
      %%%This function rebins the distribution function to shifted reference
      %%%frame and a rotated coordinate system.
      %%% Input:
      %%% pdist: 3D skymap distribution function
      %%% v_nf transformation velocity in km/s.
      %%% orient: 3x3 rotation matrix between the old and the new cooridnate system.
      %%% When rotating to a Field alligned coordinate system, if you want the
      %%% theta to be the pitch angle, the orient input has to be of the form
      %%% [t1x t1y t1z; t2x t2y t2z; bx by bz] where t1x, t1y, t1z, t2x, t2y,
      %%% t2z, bx, by, bz are the dsl components of the t1 and t2
      %%% vectors (perpendicular to the magnetic field) and b vector (parallel to
      %%% the magnetic field.
      %%%
      %%% The flag returndiff if set to 1 returns the descrete phase space volume
      %%% elements d3v. It is 0 by default.
      %%% The flag new_grid allows the user to input their own new grid. After the
      %%% flag the user have to give in order: edges of new energy bins (in eV), edges of
      %%% new polar bins (in degrees), and edges of new azimuthal bins (in degrees).

      Var = varargin;

      returndiff = 0;%%by default don't return differentials
      flag_newgrid = 0;%%by default take the same old grid but change first and last energy bins to capture full distribution in new frame.
      while ~isempty(Var)
        flag = Var{1};
        switch flag
          case 'returndiff'
            returndiff =Var{2};
            Var(1:2) =[];
          case 'new_grid'
            flag_newgrid = 1;
            Eedgesn = Var{2};
            thedgesn = Var{3}; % CN: illogical order because th is depend{3}
            phedgesn = Var{4}; % CN: illogical order because ph is depend{2}
            Var(1:4) = [];
          otherwise
            error(['undefined input flag: ' flag])

        end
      end


      u = irf_units;
      pdist = pdist.convertto('s^3/m^6');%put in SI units
      sp = pdist.species;
      switch sp
        case 'ions'
          m = u.mp;
        case 'electrons'
          m = u.me;
      end
      dat = squeeze(pdist.data);
      E_old = pdist.depend{1};
      dEo = diff(E_old);
      Eedgeso = [ E_old(1:end-1)-dEo/2 E_old(end)-dEo(end)/2 E_old(end)+dEo(end)/2];
      Eedgeso(Eedgeso<0) = 0;
      %%Define new coordinate system

      e1 = orient(1,:);
      e2 = orient(2,:);
      e3 = orient(3,:);


      %%use same old polar and azimuthal grid points for new grid
      th = pdist.depend{3};
      %fix for the case when theta and phi points are not in ascending order
      if ~issorted(th);[th,ithsort] = sort(th);dat = dat(:,:,ithsort);end
      dth = diff(th);
      thedges = [th(1:end-1)-dth/2 th(end)-dth(end)/2 th(end)+dth(end)/2];

      ph = pdist.depend{2};
      dph = diff(ph);
      phedges = [ph(1:end-1) - dph/2 ph(end)-dph(end)/2 ph(end)+dph(end)/2];
      %Incorporates SolO where phi is from -180 to 180, make it from 0 to 360
      if ~issorted(ph);[ph,iphsort] = sort(ph);dat = dat(:,iphsort,:);phedges = sort((phedges));end

      %%velocity grid before shifting reference frame
      vedges_old = sqrt(2*u.e*Eedgeso./m);
      dv_old = diff(vedges_old);%% velocity increment in the old bin

      %%calculate old d3v and fd3v
      dv3o = diff(vedges_old.^3)/3;%%dv = int v^2dv = (v^3(2) - v^3(1))/3
      dcosth = -diff(cosd(thedges));%%int sin(theta) dth = -( cos(theta2) - cos(theta1))
      dp = diff(phedges)*pi/180;%%int dphi = phi2-phi1
      [dV3o,dPHo,dcosTHo] = ndgrid(dv3o,dp,dcosth);
      d3vo = dV3o.*dcosTHo.*dPHo;
      fd3vo = dat.*d3vo;
      %% shift reference frame and rotate to new coordinate system
      [vx,vy,vz] = pdist.v;
      %%%the function PDist.v gives vx vy vz of the electrons. But when rebinning
      %%%I should use vx vy and vz of the instrument bins. Note the instrument
      %%%bins sees the direction  where an electron is coming from so an electron
      %%%with v = (vx,vy,vz) will be observed in the instrument bin -(vx,vy,vz).
      vxn = -(squeeze(vx) - v_nf(1));
      vyn = -(squeeze(vy) - v_nf(2));
      vzn = -(squeeze(vz) - v_nf(3));

      vxn = vxn*1000;vyn = vyn*1000;vzn = vzn*1000;%put in m/s
      %%rotate to new coordinate system

      vxnr =  vxn*e1(1) + vyn*e1(2) + vzn*e1(3);
      vynr =  vxn*e2(1) + vyn*e2(2) + vzn*e2(3);
      vznr =  vxn*e3(1) + vyn*e3(2) + vzn*e3(3);

      %%get the spherical coordinates of grid points in the shifted reference
      %%frame
      [~,~,v_o] = cart2sph(vxnr,vynr,vznr);


      %% generate old velocity grid for shifted reference frame and roatated coordinate system
      %%%determine the location of the grid edges in spherical coordinates in
      %%%the new frame.

      [VO,PHO,THO] = ndgrid(vedges_old,phedges,thedges);

      %% generate new velocity grid
      if ~flag_newgrid
        %%%the new velocity grid will have the same discretization as the old one
        %%%except the start and end point will be determined from the minimum and
        %%%maximum values of the speed in the shifted reference frame.

        vmin = min(min(min(v_o)));
        vmin = min([vedges_old(1),vmin-dv_old(1)/2]);
        % if norm(v_nf)>vedges_old(1)/1000; vmin =0;end
        if vmin<0; vmin =0;end

        vmax = max(max(max(v_o)));
        vmax = max([vedges_old(end),vmax+dv_old(end)/2]);
        Vedges = [vmin vedges_old(2:end-1) vmax];

      else
        Vedges = sqrt(2*u.e*Eedgesn./m);
        thedges = thedgesn;
        phedges = phedgesn;

      end
      if strcmpi(sc,'solo')
        Emin = 0.5*m*(norm(v_nf)*1000)^2/u.e;
        if Eedgeso(1)-Emin<0
          Vedges(1) = 0;
        end
        thedges = linspace(0,180,45);
        phedges = linspace(0,360,67);
      end

      %get new velocity and energies at the center of each new grid box
      %         dV = diff(Vedges);
      %         Vn = Vedges(1:end-1)+dV/2;
      %         En = 0.5*m*Vn.^2/u.e;

      Eedgesn = 0.5*m*Vedges.^2/u.e;
      En = Eedgesn(1:end-1) + diff(Eedgesn)/2;
      %         En(2:end-2) = E_old(2:end-2);


      dthn = diff(thedges);
      thn = thedges(1:end-1)+dthn/2;


      dphn = diff(phedges);
      phn = phedges(1:end-1)+dphn/2;

      %%new d3v
      dv3 = diff(Vedges.^3)/3;
      dcosth = -diff(cosd(thedges));
      dp = diff(phedges)*pi/180;
      [dV3,dPH,dcosTH] = ndgrid(dv3,dp,dcosth);
      d3vn = dV3.*dcosTH.*dPH;

      %% rebin to new grid in shifted reference frame
      l1o = size(v_o,1);l2o = size(v_o,2);l3o = size(v_o,3);
      l1 = size(d3vn,1);l2 = size(d3vn,2);l3 = size(d3vn,3);
      fd3vn = zeros(l1,l2,l3);
      for i = 1:l1o

        for j = 1:l2o

          for k = 1:l3o


            fd3vperm = fd3vo(i,j,k)/nMC;%% fd3vperm = fd3v per monte carlo point
            if fd3vperm == 0
              continue;
            end




            %%get bin edges in the non-rotated, non-transformed reference
            %%frame (from now on called F1)

            Vm_o = VO(i,j,k);
            Vp_o = VO(i+1,j,k);

            Phm_o = PHO(i,j,k);
            Php_o = PHO(i,j+1,k);

            Thm_o = THO(i,j,k);
            Thp_o = THO(i,j,k+1);

            %%generate Monte Carlo points in the bin in F1

            vo_MC = rand(nMC,1)*(Vp_o - Vm_o)+Vm_o;
            tho_MC = rand(nMC,1)*(Thp_o - Thm_o)+Thm_o;
            pho_MC = rand(nMC,1)*(Php_o - Phm_o)+Phm_o;

            %%Put each MC vector in cartesian coordinate, transform to new
            %%reference frame and rotated coordinate system (from now on
            %%called F2) then get the spherical coordinates in that new
            %%frame

            %%%as before the - sign is to make the velocities of the electrons instead
            %%%of the instrument. After shifting to the new reference frame I change it
            %%%back to that of the isntrument.

            VXE = -vo_MC.*sind(tho_MC).*cosd(pho_MC);
            VYE = -vo_MC.*sind(tho_MC).*sind(pho_MC);
            VZE = -vo_MC.*cosd(tho_MC);

            VXEnf = -(VXE - v_nf(1)*1000);
            VYEnf = -(VYE - v_nf(2)*1000);
            VZEnf = -(VZE - v_nf(3)*1000);

            VXEnfr = VXEnf*e1(1) + VYEnf*e1(2) + VZEnf*e1(3);
            VYEnfr = VXEnf*e2(1) + VYEnf*e2(2) + VZEnf*e2(3);
            VZEnfr = VXEnf*e3(1) + VYEnf*e3(2) + VZEnf*e3(3);

            [ph_MC,th_MC,v_MC] = cart2sph(VXEnfr,VYEnfr,VZEnfr);

            %cat2sph returns elevation angle, make it polar angle in the range
            %(0,180) instead of -pi/2 to pi/2
            th_MC = (pi/2-th_MC)*180/pi;

            %put in degrees and 0 to 360 range (instead of -pi to pi)
            ph_MC = ph_MC*180/pi;
            ph_MC = wrapTo360(ph_MC);

            if phedges(1)<0 && max(phedges)~=360
              ph_MC(ph_MC>max(phedges)) = ph_MC(ph_MC>max(phedges)) - 360;
            end





            %%Discretize MC into new grid
            iv = discretize(v_MC,Vedges);
            ith_n = discretize(th_MC,thedges);
            iph_n = discretize(ph_MC,phedges);

            loc = sub2ind([l1 l2 l3],iv,iph_n,ith_n);
            loc(isnan(loc)) = []; % values that fall outside of box becomes nan, remove these
            hasdata = all(loc>0, 2);

            %%sum all fd3vperm for each MC in each new bin from the old bin
            sum = accumarray(loc(hasdata,:),fd3vperm,[numel(fd3vn) 1]);
            fd3vn = fd3vn(:,:,:) + reshape(sum,[ l1 l2 l3]);


          end

        end

      end
      %%
      % fd3vn(fd3vn == 0 ) = nan; % CN: nan not compatible with other PDist methods
      fn = fd3vn./d3vn;
      Fn = zeros([1,size(fn)]);
      Fn(1,:,:,:) = fn;
      PDistn = PDist(pdist.time,Fn,'skymap',En,phn,thn);

      PDistn.ancillary.V_edges = Vedges/1000;%km/s
      PDistn.ancillary.phi_edges = phedges;
      PDistn.ancillary.theta_edges = thedges;
      PDistn.ancillary.base = 'sph';
      PDistn.units = pdist.units;
      PDistn.species = pdist.species;
      PDistn.ancillary.energy0 = En;
      PDistn.ancillary.energy1 = En;
      PDistn.ancillary.esteptable = 0;
      PDistn.ancillary.energy = En;
      PDistn.ancillary.delta_energy_plus = -En + Eedgesn(2:end);
      PDistn.ancillary.delta_energy_minus = En - Eedgesn(1:end-1);

      if returndiff

        differentials.dV3 = dV3;
        differentials.dcosTh = dcosTH;
        differentials.dPH = dPH;
        differentials.d3v = d3vn;
        varargout{1} = PDistn;
        varargout{2} = differentials;
      else
        varargout{1} = PDistn;
      end

    end
    function PD = smooth(obj,step)
      % PDIST.SMOOTH Running average
      % PD = smooth(obj,step)
      % method: It just applies obj.resample on a step-basis

      PD = obj;
      %dists = cell(step,1);
      for istep = 1:step
        dist_tmp = obj.resample(obj.time(istep:step:end));
        PD.data(istep:step:end,:) = dist_tmp.data(:,:);
      end
    end
    function [ax,args,nargs] = axescheck_pdist(varargin)
      %[ax,args,nargs] = axescheck_pdist(varargin{:});
      % MATLAB's axescheck only checks if the first argument is an axis handle, but
      % PDist.surf can be called as both PDist.surf(ax) and surf(ax,PDist),
      % where the first option has the axis as the second argument while
      % the second option has the axis as the first argument. Instead we
      % search until we find the first axis handle.
      have_axes = 0;
      ax = [];
      orig_args = varargin;
      args = orig_args;
      nargs = numel(args);
      iArg = 0;
      while ~have_axes
        iArg = iArg + 1;
        if iArg > nargs, break; end
        if ((isscalar(orig_args{iArg}) && isgraphics(orig_args{iArg},'axes')) ...
            || isa(orig_args{iArg},'matlab.graphics.axis.AbstractAxes') || isa(orig_args{iArg},'matlab.ui.control.UIAxes'))
          ax = orig_args{iArg};
          axes(ax); % set current axis to ax, needed for text(), which does not accept ax as input
          have_axes = 1;
          args(iArg) = []; % remove axes from args
          nargs = nargs - 1;
        end
      end
    end
    function varargout = plot_isosurface(varargin)
      % PDIST.PLOT_ISOSURFACE Plots isosurfaces of skymap.
      %   Might require some newer Matlab version and specific toolbox.
      %
      % Example:
      %   pdist = ePDist2(100);
      %   h(1) = subplot(2,2,1);
      %   h(2) = subplot(2,2,2);
      %   h(3) = subplot(2,2,3);
      %   h(4) = subplot(2,2,4);
      %   hs1 = pdist.plot_isosurface(h(1));
      %   hs2 = pdist.plot_isosurface(h(2),'val',(1:3)*1e-27);
      %   hs3 = pdist.plot_isosurface(h(3),'val',(1:3)*1e-27,'smooth',3);
      %   hs4 = pdist.plot_isosurface(h(4),'val',(1:3)*1e-27,'smooth',3,'fill');
      %   linkprop(h,{'View'})

      % Check for axes
      [ax,args,nargs] = axescheck_pdist(varargin{:});
      if isempty(ax); ax = gca; end
      all_handles.Axes = ax;

      % Make sure first non axes-handle input is PDist of the right type.
      if isa(args{1},'PDist') && any(strcmp(args{1}.type,{'skymap'}))
        dist_orig = args{1};
      else
        error('First input that is not an axes handle must be a PDist of type ''skymap''.')
      end
      args = args(2:end);
      nargs = nargs - 1;

      pdist = dist_orig;

      tId = 1:pdist.length; % if tint is not given as input (check below) default is to include all the time indices of input distribution


      % Default values
      doScPot = 0;
      doSmooth = 0;
      doPrintInfo = 0;
      doFillGap = 0;
      doRotate = 0;
      T = [1 0 0; 0 1 0; 0 0 1]; % no rotation

      % Default formatting
      faceAlpha = 0.5;
      iso_values = [];
      colors = get(ax,'colororder');


      % check for input, try to keep it at a minimum, so that the
      % functionality is similar to Matlabs plot function, all the details
      % can then be fixed outside the function using ax.XLim, ax.YLim,
      % ax.CLim, etc... and colorbar perhaps?
      if nargs > 0; have_options = 1; else have_options = 0; end
      while have_options
        l = 1;
        switch(lower(args{1}))
          case {'rotate'}
            doRotate = 1;
            T = args{2};
            l = 2;
          case {'facealpha'}
            faceAlpha = args{2};
            l = 2;
          case {'fill'}
            doFillGap = 1;
            l = 1;
          case {'val','iso_val','values'}
            l = 2;
            iso_values = args{2};
          case {'prctile','percentile'}
            l = 2;
            iso_values_percentile = args{2};
          case {'tint','time','t'} % not implemented
            l = 2;
            notint = 0;
            tint = args{2};
            if tint.length == 1 % find closest time
              [~,tId] = min(abs(pdist.time-tint));
            else % take everything within time interval
              [tId,~] = pdist.time.tlim(tint);
            end
          case 'vectors'
            l = 2;
            vectors = args{2};
            have_vectors = 1;
          case 'scpot'
            l = 2;
            scpot = args{2};
            % if isa(scpot,'TSeries')
            %   includescpot = 1;
            %   irf.log('notice','Spacecraft potential passed.')
            % else
            %   includescpot = 0;
            %   irf.log('notice','scpot not recognized. Not using it.')
            % end
            doScPot = 1;
          case '10^3 km/s' % not implemented
            l = 1;
            v_scale = 1e-3;
            v_label_units = '10^3 km/s';
          case 'km/s' % not implemented
            l = 1;
            v_scale = 1;
            v_label_units = 'km/s';
          case 'printinfo' % not implemented
            doPrintInfo = 1;
          case 'smooth'
            l = 2;
            doSmooth = 1;
            nSmooth = args{2};
          case {'colorbar','docolorbar'} % not implemented
            l = 2;
            doColorbar = args{2};
        end
        args = args(l+1:end);
        if isempty(args), break, end
      end


      % due to Matlab functionality, we must explicitly call the overloaded
      % subsref (defined within this subclass), otherwise it will call the
      % builtin function
      subs.type = '()';
      subs.subs = {tId};
      pdist = dist_orig.subsref(subs);
      if (length(pdist.time)<1); irf.log('warning','No data for given time interval.'); return; end

      F = pdist.data;
      F = mean(F,1); % average over time, if multiple times
      F = squeeze(F);

      if 0 % interpolate to edges of bins, NOT WORKING
        phi = TSeries(pdist.time,pdist.depend{1,2});
        azimuthal = phi.data*pi/180;

        theta = pdist.depend{1,3};
        polar = repmat(theta*pi/180,pdist.length,1);
        energy = pdist.depend{1};
        doEdges = 1;
        if doEdges % get V at edges of bins
          azimuthal_e = [azimuthal(1) - 0.5*(azimuthal(2)-azimuthal(1)) azimuthal(1:end-1) + 0.5*diff(azimuthal) azimuthal(end) + 0.5*(azimuthal(end)-azimuthal(end-1))];
          polar_e = [polar(1) - 0.5*(polar(2)-polar(1)) polar(1:end-1) + 0.5*diff(polar) polar(end) + 0.5*(polar(end)-polar(end-1))];
          de_minus = pdist.ancillary.delta_energy_minus;
          de_plus = pdist.ancillary.delta_energy_plus;
          energy_e = [energy(:,1) - de_minus(:,1) energy + de_plus];
        end
        [EN,AZ,POL] = meshgrid(energy,azimuthal,polar);
        [ENe,AZe,POLe] = meshgrid(energy_e,azimuthal_e,polar_e);
        Fe = interp3(EN,AZ,POL,F,ENe,AZe,POLe);
        Fe = permute(Fe,[2 1 3]); % from meshgrid to ndgrid
      end

      if doScPot
        [VX,VY,VZ] = pdist.v('scpot',scP);
        %[VXe,VYe,VZe] = pdist.v('scpot',scP,'edges');
      else
        [VX,VY,VZ] = pdist.v;
        %[VXe,VYe,VZe] = pdist.v('edges');
      end
      % Assume all energy levels are the same
      VX = squeeze(mean(VX,1));
      VY = squeeze(mean(VY,1));
      VZ = squeeze(mean(VZ,1));
      %VXe = squeeze(mean(VXe,1));
      %VYe = squeeze(mean(VYe,1));
      %VZe = squeeze(mean(VZe,1));


      %Fg = griddedInterpolant(VX,VY,VZ,F);
      %Fq = Fg(VXe,VYe,VZe);
%      Fq = interp3(VX,VY,VZ,F,VXe,VYe,VZe);

      %F = Fe;
      %VX = VXe;
      %VY = VYe;
      %VZ = VZe;


      if doFillGap % make circular in azimuth, to fill in the gap
        icat_az = 1:(nSmooth-1)+3;
        F = cat(2,F,F(:,icat_az,:));
        VX = cat(2,VX,VX(:,icat_az,:));
        VY = cat(2,VY,VY(:,icat_az,:));
        VZ = cat(2,VZ,VZ(:,icat_az,:));

        if 0
          icat_pol = 1;
          F = cat(3,F,F(:,:,end)*0+mean(F(:,:,end),3));
          VX = cat(3,VX,VX(:,:,end)*0+mean(VX(:,:,end),3));
          VY = cat(3,VY,VY(:,:,end)*0+mean(VY(:,:,end),3));
          VZ = cat(3,VZ,VZ(:,:,end)*0+mean(VZ(:,:,end),3));
        end
      end

      if doSmooth
        F = smooth3(F,'box',nSmooth);
      end

      if doFillGap
        F = F(:,2:end-icat_az,:);
        VX = VX(:,2:end-icat_az,:);
        VY = VY(:,2:end-icat_az,:);
        VZ = VZ(:,2:end-icat_az,:);
      end

      if isempty(iso_values)
        iso_values = prctile(F(:),[50 70 90]);
      end

      if doRotate
        VxX = reshape(VX,numel(VX),1);
        VyY = reshape(VY,numel(VX),1);
        VzZ = reshape(VZ,numel(VX),1);

        newTmpX = [VxX VyY VzZ]*T(1,:)';
        newTmpY = [VxX VyY VzZ]*T(2,:)';
        newTmpZ = [VxX VyY VzZ]*T(3,:)';

        VX = reshape(newTmpX,size(VX));
        VY = reshape(newTmpY,size(VY));
        VZ = reshape(newTmpZ,size(VZ));
        all_handles.Rotation = T;
      end

      nSurf = numel(iso_values);

      if numel(faceAlpha) == 1
        faceAlpha = repmat(faceAlpha,size(iso_values));
      end
      %holdon = 0;

      if strcmp(ax.NextPlot,'replace')
        delete(ax.Children)
      end

      hps = gobjects(0);
      for isurf = 1:nSurf
        %if not(holdon); hold(ax,'on'); holdon = 1; end
        Flev = iso_values(isurf);
        s = isosurface(VX,VY,VZ,F,Flev);
        %s = isocaps(VX,VY,VZ,F,Flev);
        %p = patch(ax,Faces=s.faces,Vertices=s.vertices);
        p = patch(ax,'Faces',s.faces,'Vertices',s.vertices);
        hps(isurf) = p;

        % Default formatting
        p.EdgeColor = 'none';
        p.FaceColor = colors(mod(isurf-1,size(colors,1))+1,:).^0.5; % cycle colors
        p.FaceAlpha = faceAlpha(isurf);
      end
      %hold(ax,'off') % not with patch objects

      % might be leftover from previous plottings, only add if no lights
      hlight = findobj(ax,'Type','Light');
      if isempty(hlight)
        hlight = camlight(ax);
      end

      lighting(ax,'gouraud')

      view(ax,[2 2 1])
      %view(ax,[0 0 1])

      all_handles.Patch = hps;
      all_handles.Light = hlight;

      ax.XLabel.String = 'v_{x} (km/s)';
      ax.YLabel.String = 'v_{y} (km/s)';
      ax.ZLabel.String = 'v_{z} (km/s)';
      axis(ax,'square')
      axis(ax,'equal')

      if 1 % Print title with time
        if pdist.length == 1
          ax.Title.String = pdist.time.utc('yyyy-mm-ddTHH:MM:SS:mmmZ');
        else
          t1 = pdist.time.start.utc('yyyy-mm-ddTHH:MM:SS:mmmZ');
          t2 = pdist.time.stop.utc('yyyy-mm-ddTHH:MM:SS:mmmZ');
          tcomp = t1 == t2;
          if sum(tcomp(1:11))==11
            ax.Title.String = {sprintf('%s - %s',t1,t2(12:end))};
          else
            ax.Title.String = {sprintf('%s - %s',t1,t2)};
          end
        end
      end
      if 1 % Print iso values as legend
        legs = arrayfun(@(x)sprintf('%g',x),iso_values,'UniformOutput',false);
        hleg = legend(hps,legs,'location','northeast','box','off');
        hleg.Title.String = sprintf('f_%s (%s)',pdist.species(1),pdist.units);
        all_handles.Legend = hleg;
      end


      varargout{1} = all_handles;
    end
    function varargout = plot_plane(varargin)
      % PDIST.PLOT_PLANE Surface/pcolor plot of PDist of type 'plane (reduced)'.
      %   PDist.PLOT_PLANE(...)
      %   h_surf = PDist.PLOT_PLANE(...)
      %   [h_surf,h_axis] = PDist.PLOT_PLANE(...)
      %   [h_surf,h_axis,h_all] = PDist.PLOT_PLANE(...)
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
      %     'log10'/value - if value = 0, does not take log10 of data,
      %             default is to take log10, i.e. value = 1,
      %     'contour'/contour_levels - contour levels drawn in black, if
      %             option 'log' is passed, the function will also do
      %             log10(contour_levels).
      %             If one value is passed, this is the number of contours
      %             that will be drawn at levels decided by Matlab.
      %             If contour_levels is empty [], a default of 10 levels
      %             will be drawn.
      %     contourf/contour_levels - overrides the pcolor plot with filled
      %             contour levels, using Matlabs contourf function
      %     '10^3 km/s' - plots x and y axis velocities in 10^3 km/s, this
      %             is default for electrons
      %     'km/s' - plots x and y axis velocities in km/s, this is
      %             default for ions
      %     'printinfo' - prints info - number of distributions used, time
      %             or time interval, vint used for distribution
      %             integration (see PDist.reduce), v1/v2 directions
      %     'flim'/flim - sets values outside of this range to NaN -
      %             default is [0 Inf]
      %     'colorbar'/value - value is 0 for no colorbar or 1 for
      %             colorbar, default is to plot colorbar
      %     'off-diag-pres-cont'/v1/v2 - plots contributions to pressure
      %             tensor components, v1 and v2 are TSeries with the bulk
      %             speeds from moments, bulk speed from distribution not
      %             yet implemented
      %
      %     Example:
      %       tint = irf.tint('2017-07-06T13:53:50.00Z',25);
      %       tint_plot = irf.tint('2017-07-06T13:54:00.00Z',5);
      %       eint = [000 40000];
      %       vint = [-Inf Inf];
      %       iDist = iPDist1.tlim(tint).elim(eint);
      %       iLine = dmpaB1.resample(iDist).norm;
      %       iPlane1 = iLine.cross(irf.ts_vec_xyz(iLine.time,repmat([1 0 0],iLine.length,1)));
      %       iPlane2 = iLine.cross(iPlane1);
      %       if2D = iDist.reduce('2D',iPlane1,iPlane2,'vint',vint); % reduced distribution perp to B
      %       h1 = subplot(3,1,1);
      %       if2D.PLOT_PLANE(h1,'printinfo','tint',tint_plot(1))
      %       h2 = subplot(3,1,2);
      %       if2D.PLOT_PLANE(h2,'printinfo','tint',tint_plot(2))
      %       h3 = subplot(3,1,3);
      %       if2D.PLOT_PLANE(h3,'printinfo','tint',tint_plot)
      %
      %   See also mms.plot_int_projection, PDist.reduce

      % Check for axes
      [ax,args,nargs] = axescheck_pdist(varargin{:});
      if isempty(ax); ax = gca; end
      all_handles.Axes = ax;

      % Make sure first non axes-handle input is PDist of the right type.
      if isa(args{1},'PDist') && any(strcmp(args{1}.type,{'plane (reduced)','plane (slice)','2Dcart'}))
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
      doContourFill = 0;
      doCircles = 0;
      doFLim = 1; flim = [0 Inf];
      doSmooth = 0;
      doP12 = 0;
      doB = 0;
      doStress = 0;
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
      % ax.CLim, etc... and colorbar perhaps?
      if nargs > 0; have_options = 1; else have_options = 0; end
      while have_options
        l = 1;
        switch(lower(args{1}))
          case 'stress'
            l = 1;
            doStress = 1;
          case 'off-diag-pres-cont'
            l = 3;
            doP12 = 1;
            v1 = args{2};
            v2 = args{3};
            if isempty(v1)
              doCalculateVbulk = 1;
            else
              doCalculateVbulk = 0;
            end
          case {'tint','time','t'}
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
          case 'nolog10' % backwards compatibility
            l = 1;
            doLog10 = 0;
          case 'log10'
            l = 2;
            doLog10 = args{2};
          case 'contour'
            l = 2;
            contour_levels = args{2};
            doContour = 1;
          case 'contourf'
            l = 2;
            contour_levels = args{2};
            doContour = 1;
            doContourFill = 1;
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
          case 'smooth'
            l = 2;
            doSmooth = 1;
            nSmooth = args{2};
          case {'colorbar','docolorbar'}
            l = 2;
            doColorbar = args{2};
          case 'b'
            l = 2;
            doB = 1;
            B = args{2};
        end
        args = args(l+1:end);
        if isempty(args), break, end
      end
      if doP12, doLog10 = 0; end
      if doStress, doLog10 = 0; end

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
      if doSmooth
        plot_data = smooth2(plot_data,nSmooth);
      end
      if doP12
        [V1,V2] = ndgrid(dist.depend{1}(1,:),dist.depend{2}(1,:)); % km/s
        if doCalculateVbulk % incase moments are bad, this might be a better firs/temporary option
          % NOT IMPLEMENTED
          vx = dist.units;
        else
          v1_bulk = mean(v1.resample(dist).data,1);
          v2_bulk = mean(v2.resample(dist).data,1);
        end

        V1V2 = (V1-v1_bulk).*(V2-v2_bulk)*1e3*1e3; % km/s -> m/s
        % units_scaling
        new_units = 'arb. units';
        vp12_scale = 1;
        switch dist.units
          case 's^2/m^5'
            vp12_scale = 1e0; % depend: m -> m
            new_units = 'm^{-3}'; % [f] m/s*m/s = s2m-5*m*m*s-1*s-1 = m-3
          otherwise
            disp('unsupported/unimplemented units, output is in arb. units')
        end
        plot_data = vp12_scale*plot_data.*V1V2; % e.g. 1*s2*m-5*ms-1*ms-1 = m-3
        dv1 = diff(squeeze(irf.nanmean(dist.ancillary.vx_edges,1)))*1e3; % km/s -> m/s
        dv2 = diff(squeeze(irf.nanmean(dist.ancillary.vy_edges,1)))*1e3; % km/s -> m/s

        integrated_p12 = dist.mass*nansum(nansum(plot_data.*(dv1*dv2'))); % kg*m-3*m/s*m/s = kg*m*s-2 = Pa
        % pascal is kg*m*s-2, so go to nPa
        integrated_p12 = integrated_p12*1e9; % Pa -> nPa
      end
      if doStress
        [V1,V2] = ndgrid(dist.depend{1}(1,:),dist.depend{2}(1,:)); % km/s
        V1V2 = V1.*V2*1e3*1e3; % km/s -> m/s
        % units_scaling
        new_units = 'arb. units';
        vp12_scale = 1;
        switch dist.units
          case 's^2/m^5'
            vp12_scale = 1e0; % depend: m -> m
            new_units = 'm^{-3}'; % [f] m/s*m/s = s2m-5*m*m*s-1*s-1 = m-3
          otherwise
            disp('unsupported/unimplemented units, output is in arb. units')
        end
        plot_data = vp12_scale*plot_data.*V1V2; % e.g. 1*s2*m-5*ms-1*ms-1 = m-3
        dv1 = diff(squeeze(irf.nanmean(dist.ancillary.vx_edges,1)))*1e3; % km/s -> m/s
        dv2 = diff(squeeze(irf.nanmean(dist.ancillary.vy_edges,1)))*1e3; % km/s -> m/s

        integrated_p12 = dist.mass*nansum(nansum(plot_data.*(dv1*dv2'))); % kg*m-3*m/s*m/s = kg*m*s-2 = Pa
        % pascal is kg*m*s-2, so go to nPa
        integrated_p12 = integrated_p12*1e9; % Pa -> nPa
      end
      % main surface plot
      % NOTE, PCOLOR and SURF uses flipped dimensions of (x,y) and (z), but PDist.reduce does not, there we need to flip the dim of the data
      plot_x_edges = squeeze(irf.nanmean(dist.ancillary.vx_edges,1))*v_scale; % v_scale, default 1e-3 for electrons to put axes in 10^3 km/s
      plot_y_edges = squeeze(irf.nanmean(dist.ancillary.vy_edges,1))*v_scale;

      if strcmpi(dist.ancillary.base,'pol')
        plot_z_edges = plot_x_edges*0;
      elseif strcmpi(dist.ancillary.base,'cart')
        plot_z_edges = zeros(length(plot_x_edges),length(plot_y_edges));
      end
      ax_surface = surf(ax,plot_x_edges,plot_y_edges,plot_z_edges,plot_data');
      all_handles.Surface = ax_surface;
      view(ax,[0 0 1])
      ax.Box = 'on';
      shading(ax,'flat');

      if doB

        if isa(B,'TSeries')
          dt = obj.time(2) - obj.time(1);
          B1 = B{1}.tlim(dist.time([1 end]) + 0.5*dt*[-1 1]);
          B2 = B{2}.tlim(dist.time([1 end]) + 0.5*dt*[-1 1]);
          B1 = mean(B1.data,1);
          B2 = mean(B2.data,1);
          b1 = B1/norm(B1);
          b2 = B2/norm(B2);
        else % single vector value
          b = B/norm(B);
        end

        hold(hca,'on')
        Bnorm = B;
        hold(hca,'off')
      end
      if doP12 % add info about integrated value
        if integrated_p12 >= 0
          sumf_color = 'r';
        else
          sumf_color = 'b';
        end
        %irf_legend(ax,sprintf('m*int f vv dv2 = %.6f nPa',integrated_p12),[0.02 0.02],'k')
        %irf_legend(ax,sprintf('p = %.6f nPa',integrated_p12),[0.02 0.02],'color',sumf_color)
        irf_legend(ax,sprintf('p = %.3f pPa',integrated_p12*1e3),[0.98 0.98],'color',sumf_color,'fontsize',12)
      end
      %if doP12 % add info about integrated value
        %irf_legend(ax,sprintf('m*int f vv dv2 = %.6f nPa',integrated_p12),[0.02 0.02],'k')
      %end

      if doContour
        hold(ax,'on')
        if isempty(contour_levels), contour_levels = 10;
        elseif numel(contour_levels) > 1 && doLog10, contour_levels = log10(contour_levels);
        end
        if strcmp(dist.ancillary.base,'pol')
          if 1
            plot_data_tmp = zeros(size(plot_data,1)+1,size(plot_data,2)+1);
            plot_data_tmp(2:end,2:end) = plot_data;
            plot_data_tmp(2:end,1) = plot_data(:,end);
            plot_data_tmp(1,2:end) = plot_data(end,:);
            plot_data_tmp(1,1) = plot_data(end,end);
            plot_data = plot_data_tmp;
          else
            plot_data(:,end+1) = plot_data(:,1);
            plot_data(end+1,:) = plot_data(1,:);
          end
          plot_x = plot_x_edges;
          plot_y = plot_y_edges;
        else
          plot_x = squeeze(irf.nanmean(dist.depend{1},1))*v_scale; % v_scale, default 1e-3 for electrons to put axes in 10^3 km/s
          plot_y = squeeze(irf.nanmean(dist.depend{2},1))*v_scale;
        end
        if doContourFill
          [~,h_contour] = contourf(ax,plot_x,plot_y,plot_data',contour_levels,'k');
        else
          [~,h_contour] = contour(ax,plot_x,plot_y,plot_data',contour_levels,'k');
        end
        h_contour.LineWidth = 1.5;
        hold(ax,'off')
        all_handles.Contour = h_contour;
      end
      if doColorbar
        hcb = colorbar('peer',ax);
        data_units = dist.units;
        if doP12
          hcb.YLabel.String = sprintf('f(v_1,v_2)(v_1-<v_1>)(v_2-<v_2>) (%s)',new_units);
        else
          if doLog10
            hcb.YLabel.String = sprintf('log_{10} f (%s)', dist.units);
          else
            hcb.YLabel.String = sprintf('f (%s)', dist.units);
          end
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
    function varargout = plot_pad_polar(varargin)
      % PDIST.PLOT_PAD_POLAR polar pitchangle plot

      % Check for axes
      [ax,args,nargs] = axescheck_pdist(varargin{:});
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

      units = irf_units;

      % default plotting parameters
      doMirrorData = 0;
      doAxesV = 1; v_scale = 1; v_label_units = 'km/s'; % default is to do energy
      doLog10 = 1;
      doLogAxes = 0;
      doColorbar = 1;
      doAxisLabels = 1;
      doPrintInfo = 0;
      doContour = 0;
      doCircles = 0;
      doScpot = 0;
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
            doScpot = 1;
            if isa(scpot,'TSeries')
              scpot = scpot.resample(dist).data;
              irf.log('notice','scpot was TSeries.')
            elseif isnumeric(scpot) && numel(scpot) == 1
              irf.log('notice','scpot was scalar.')
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
          case {'km/s','v'}
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
      data = squeeze(irf.nanmean(dist.data,1)); % average data over time indices
      data = reshape(data,[size(dist.depend{1},2) size(dist.depend{2},2)]);
      data(data==0) = NaN; % put zero values to NaN
      if doFLim % put values outside given interval to NaN, default is [0 Inf]
        data(data<=flim(1)) = NaN;
        data(data>flim(2)) = NaN;
      end
      if doLog10 % take log10 of data
        data = log10(data);
      end

      % main surface plot
      % NOTE, PCOLOR and SURF uses flipped dimensions of (x,y) and (z), but PDist.reduce does not, there we need to flip the dim of the data
      rho_edges = [dist.depend{1}-dist.ancillary.delta_energy_minus dist.depend{1}(:,end)+dist.ancillary.delta_energy_plus(:,end)];
      if isfield(dist.ancillary,'delta_pitchangle_minus') && not(isempty(dist.ancillary.delta_pitchangle_minus))
        ntheta = numel(dist.depend{2});
        theta_minus = dist.depend{2} - dist.ancillary.delta_pitchangle_minus;
        theta_plus = dist.depend{2} + dist.ancillary.delta_pitchangle_plus;
        theta_edges = sort(unique([theta_minus;theta_plus]));
        % Look for any edges that are not common, this should be a gap;
        gaps = setdiff(theta_minus(2:end),theta_plus(1:end-1));
        ngaps = numel(gaps);
        new_itheta = 1;
        for itheta = 2:ntheta
          new_itheta = new_itheta + 1;
          if theta_minus(itheta) == theta_plus(itheta-1) % ok, the edges concide

          else % edges does not concide, need to pad with NaN
            theta_edges = [theta_edges(1:new_itheta); NaN; theta_edges(new_itheta+1:end)];
            data = [data(:,1:new_itheta-1) nan(size(data,1),2) data(:,new_itheta:end)];
            new_itheta = new_itheta + 2;
          end
        end
        1;

      elseif isfield(dist.ancillary,'pitchangle_edges') && not(isempty(dist.ancillary.pitchangle_edges))
        theta_edges = dist.ancillary.pitchangle_edges;
        if not(size(theta_edges,2)-1 == size(dist.depend{2},2)) % there are gaps in the pitchangle, for example for EDI flux
          % find gaps and pad with nans, only adapted for equally wide
          % pitch angle bins, and only one gap
          diff_theta_edges = diff(theta_edges);
          unique_diff_theta_edges = sort(unique(diff_theta_edges));
          % assume the smallest on is the proper one
          ind_pad = find(diff_theta_edges==unique_diff_theta_edges(end));
          data = [data(:,1:ind_pad-1) nan(size(data,1),2) data(:,ind_pad:end)];
          theta_edges = [theta_edges(1:ind_pad) NaN theta_edges(ind_pad+1:end)]; % also pad grid, to avoid empty boxes
          %           ngaps = size(theta_edges,2) - 1 - size(dist.depend{2},2);
          %           ngaps_remaining = ngaps;
          %           %while ngaps_remaining
          %           for iedge = 1:size(theta_edges,2)-1
          %             theta_minus_correct = dist.depend{2}(iedge);
          %             theta_plus = theta_edges(iedge+1);
          %             theta_minus = theta_edges(iedge);
          %             theta_plus = theta_edges(iedge+1);
          %           end

        end
      else
        theta = dist.depend{2}; dtheta = theta(2)-theta(1);
        theta_edges = [theta(1)-dtheta/2 theta+dtheta/2];
      end
      if doScpot
        if isscalar(scpot)
          rho_edges = rho_edges - scpot;
        else
          rho_edges = rho_edges - repmat(scpot(tId),1,size(rho_edges,2));
        end
        %rho_edges = rho_edges - repmat(scpot,1,size(rho_edges,2));
        rho_edges(rho_edges<0) = NaN;
      end
      if doAxesV
        rho_edges = sqrt(2*units.e*rho_edges/units.me)*1e-3*v_scale;
        stringLabel = sprintf('v (%s)',v_label_units);
      else
        stringLabel = sprintf('E (%s)','eV');
      end
      rho_edges = irf.nanmean(rho_edges,1); % average over times, do after removing scpot
      data(isnan(rho_edges),:) = NaN; % rho_edges<scpot was put to NaN above
      if doLogAxes, rho_edges = log10(rho_edges); stringLabel = sprintf('log_{10}(%s)%s',stringLabel(1),stringLabel(2:end)); end
      theta_edges = theta_edges + 90; % rotate data
      [RHO,THETA] = meshgrid(rho_edges,theta_edges);
      X = RHO.*cosd(THETA);
      Y = RHO.*sind(THETA);
      if doMirrorData % mirror data
        plot_X = [X; -flipdim(X(1:end-1,:),1)];
        plot_Y = [Y; flipdim(Y(1:end-1,:),1)];
        plot_data = [data flipdim(data,2)];
      else
        plot_X = -X; % the rotation by 90 puts it in negative v's, so option to put a minus sign here to have positive vperp
        plot_Y = Y;
        plot_data = data;
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
        ax.XLabel.String = stringLabel;
        ax.YLabel.String = stringLabel;
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
      % Ancillary data, problematic for pitchangle_edges since we dont
      % immediately know where gaps can be, use instead
      % pitchangle_delta_minus/plus
      if isfield(PD.ancillary,'pitchangle_delta_minus')
        PD.ancillary.pitchangle_delta_minus = PD.ancillary.pitchangle_delta_minus(:,indPA);
      end
      if isfield(PD.ancillary,'pitchangle_delta_plus')
        PD.ancillary.pitchangle_delta_plus = PD.ancillary.pitchangle_delta_plus(:,indPA);
      end
      % changed to delta_pitchangle_minus/plus to follow fpi way: delta_energy_minus/plus
      % keep above for now for backwards compatability
      if isfield(PD.ancillary,'delta_pitchangle_minus')
        PD.ancillary.delta_pitchangle_minus = PD.ancillary.delta_pitchangle_minus(:,indPA);
      end
      if isfield(PD.ancillary,'delta_pitchangle_plus')
        PD.ancillary.delta_pitchangle_plus = PD.ancillary.delta_pitchangle_plus(:,indPA);
      end

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
        %         if ~isfield(PD.ancillary, 'esteptable')
        %           [esteptable,~] = ismember(energy,PD.ancillary.energy1,'rows');
        %           PD.ancillary.esteptable = esteptable;
        %         end
        if isfield(PD.ancillary,'energy'), PD.ancillary.energy = PD.ancillary.energy(:,elevels); end
        if isfield(PD.ancillary,'delta_energy_minus'), PD.ancillary.delta_energy_minus = PD.ancillary.delta_energy_minus(:,elevels); end
        if isfield(PD.ancillary,'delta_energy_plus'), PD.ancillary.delta_energy_plus = PD.ancillary.delta_energy_plus(:,elevels); end
      end
    end
    function PD = omni(obj,V0,varargin)
      % Makes omnidirectional distribution, conserving units.
      %
      % distOmni = dist.OMNI Returns the average omnidirectional distribution
      %   function in the spacecraft frame. Units are the same as in the
      %   output.
      %
      % distOmni = dist.OMNI(V) Returns the average omnidirectional
      %   distribution in an arbitrary frame defined by the velocity V.
      %   V is the velocity of the new frame in the spacecraft frame.
      %   V is either a vector or a TSeries object with unit of [km/s].
      %
      % distOmni = dist.OMNI(V,propertyFlag,propertyValue) defines the
      %   property of the resulting omni distribution.
      %
      %   Properties:
      %     'E' -  New energy grid for the output. Array with units [eV]
      %
      % Example:
      %   % Get the ion omni distribution function in the ion rest frame
      %   idistOmniRestFrame = iPDist.omni(Vi); % Vi is the ion velocity
      %
      %
      % Note: When there is no full sky coverage, calling PDist.omni
      % without a velocity does not calculate a proper spherical mean of
      % the distribution and therefore does not conserve i.e. density. This
      % happens for SolO distributions. Instead call, PDist.omni([0,0,0])
      % to get the proper averaged distribution in the spacecraft frame.
      %

      if ~strcmp(obj.type_,'skymap'); error('PDist must be a skymap.'); end

      if nargin == 1 % Normal, classical usage
        dist = obj;
        % define angles
        energysize = size(obj.depend{1});
        phi = obj.depend{2};
        theta = obj.depend{3};
        lengthphi = size(phi,2);
        if isfield(obj.ancillary,'delta_theta_plus') && isfield(obj.ancillary,'delta_theta_minus')
          dangletheta = obj.ancillary.delta_theta_plus + obj.ancillary.delta_theta_minus;
        else
          dangletheta = median(diff(obj.depend{3}));
        end

        if isfield(obj.ancillary,'delta_phi_plus') && isfield(obj.ancillary,'delta_phi_minus')
          danglephi = obj.ancillary.delta_phi_plus + obj.ancillary.delta_phi_minus;
        else
          danglephi = median(diff(obj.depend{2}(1,:)));
        end

        dangletheta = dangletheta*pi/180;
        danglephi = danglephi*pi/180;

        z2 = ones(lengthphi,1)*sind(theta);
        solida = (danglephi*dangletheta').*z2;
        allsolida = repmat(solida,1,1,length(dist.time), energysize(2));
        allsolida = (permute(allsolida,[3 4 1 2]));
        dists = (dist.data).*allsolida;
        omni = squeeze(irf.nanmean(irf.nanmean(dists,3),4))/(mean(mean(solida)));

        PD = obj;
        PD.type = 'omni';
        PD.data_ = omni;
        PD.depend = {obj.depend{1}};
        PD.representation = {obj.representation{1}};
        PD.units = obj.units;
        PD.name = 'omni';
      else % frame transformation
        % call the protected class method
        PD = obj.get_omni_dist_framedep(V0,varargin{:});
      end



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
        case {'1/(cm^2 s sr eV)'}
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
            case 'moms-tens0'
              spec.p = double(obj.data);
            otherwise
              error('Spectype ''%s'' not yet implemented for distribution type ''%s''.',spectype,obj.type);
          end
          spec.t = obj.time.epochUnix;
          spec.f = single(obj.depend{1});
          if not(isempty(obj.species))
            spec.f_label = {['E_' obj.species(1) ' (eV)']};
          end
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
          if ~any(strcmp(obj.type_,{'line (reduced)','1Dcart'})); error('PDist must be projected unto a vector: type: ''line (reduced)'', see PDist.reduce.'); end
          %if ~strcmp(obj.type_,'line (reduced)'); error('PDist must be projected unto a vector: type: ''line (reduced)'', see PDist.reduce.'); end
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
        case {'i','p','ions','ion','hplus'}
          mm = 1;
        case {'oplus'}
          mm = 16;
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
        case {'i','p','ions','ion','hplus'}
          mm = 1;
        case {'oplus'}
          mm = 16;
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
      elseif flagdir == -1
        switch obj.units
          case '1/(cm^2 s sr eV)'
            irf.log('warning','Converting DPFlux to PSD');
            tmpData = obj.data/1e12*mm^2*0.53707;
          case '1/(cm^2 s sr keV)'
            irf.log('warning','Converting DPFlux to PSD'); %% OBS, not correct units!!! Same above
            tmpData = 1e-3*obj.data/1e12*mm^2*0.53707;
          otherwise
            error(sprintf('Units: %s not supported',obj.units))
        end
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
        PD.units = '1/(cm^2 s sr eV)';
      elseif flagdir == -1
          reshapedData = reshapedData./matEnergy;
          tmpData = reshape(reshapedData,sizeData);
          PD = obj;
          PD.data_ = tmpData;
        switch obj.units
          case '1/(cm^2 s sr eV)'
            PD.units = 's^3/m^6';
          case '1/(cm^2 s sr keV)' %% OBS, not correct units!!! Same above
            PD.units = 's^3/m^6';
        end
      else
        irf.log('warning','No change to PDist');
        PD = obj;
      end
    end
    function PD = convertto(obj,newunits)
      % Changes units of Pdist.
      % Accepted inputs 's^3/cm^6', 's^3/km^6', 's^3/m^6', 'keV/(cm^2 s sr keV)',
      % and '1/(cm^2 s sr eV)'

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
        case {'1/(cm^2 s sr eV)'}
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
        case {'1/(cm^2 s sr eV)'}
          PD = PD.dpflux;
        otherwise
          error('Units not supported.');
      end
    end
    function PD = pitchangles(obj,obj1,obj2,varargin) %,method
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
        pitchangle_edges = 0:(180/nangles):180;
      elseif isnumeric(obj2) % angles or number of angles
        nangles = obj2;
        if numel(nangles) > 1
          pitchangle_edges = nangles;
        else % if nothing is passed, they are equidistanced
          pitchangle_edges = 0:(180/nangles):180;
        end
      else % obj2 is part of varargin to be passed on to mms.get_pitchangles
        varargin = {obj2,varargin{:}};
      end

      PD = mms.get_pitchangledist(obj,obj1,'angles',nangles,varargin{:}); % - For v1.0.0 or higher data

      % if the pitch angle bins are not equally spaced, we pass this for
      % plotting purposes, can be empty
      PD.ancillary.pitchangle_edges = pitchangle_edges;
    end
    function PD = pitchangle_diffusion(obj,Dth)
      % Calculates df/dt for pitch angle distribution with some assumptions
      % Under development. Only implemented for single distribution.
      %   dfdt = 1/sin(th)*d/dth(D*sin(th)*df/dth)
      %   th ('pitch angle') should be with respect to the phase speed of the
      %   wave, but here we assume it to be close to zero (or equivalently,
      %   very far away from the resonant energy, so th is stricly the
      %   pitch angle
      %
      %   Example:
      %     pitch_dfdt = pitch.pitchangle_diffusion(Dth);
      %     pitch_dfdt.plot_pad_polar(hca,'10^3 km/s','flim',[-inf inf],'nolog10')
      %
      th_edges = obj.ancillary.pitchangle_edges;
      dth = diff(th_edges);
      th = th_edges(1:end-1) + 0.5*dth;
      dt = mean(dth); % assumes they are equispaced
      TH = repmat(th,[numel(obj.depend{1}(1,:)),1]);

      f = squeeze(obj.data);
      dfdth = f*0;
      dfdth(:,2:end-1) = (f(:,3:end)-f(:,1:end-2))/(2*dt);
      dfdth(:,1) = (f(:,2)-f(:,1))/(dt);
      dfdth(:,end) = (f(:,end)-f(:,end-1))/(dt);
      %[~,dfdth_] = gradient(f',E,TH);
      Dsinthdfdth = Dth*sind(TH).*dfdth;
      %imagesc(dfdth')

      %[~,dfdt__] = gradient(Dsinthdfdth,E',TH');

      dfdt__ = dfdth*0;
      dfdt__(:,2:end-1) = (Dsinthdfdth(:,3:end)-Dsinthdfdth(:,1:end-2))/(2*dt);
      dfdt__(:,1) = diff(Dsinthdfdth(:,1:2),1,2)/(dt);
      dfdt__(:,end) = diff(Dsinthdfdth(:,end-1:end),1,2)/(dt);
      dfdt = dfdt__./sind(TH);

      PD = obj;
      PD.data = reshape(dfdt,[1 size(dfdt)]);
      PD.units = '...';
    end
    function PD = squeeze(obj)
      PD = obj;
      PD.data = squeeze(PD.data);
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
      if size(obj.depend{1},2) == 64; irf_log('proc','PDist already has 64 energy levels.'); end

      %if ~any([isfield(obj.ancillary,'energy0') isfield(obj.ancillary,'energy1') isfield(obj.ancillary,'esteptable')]) % construct energy0, energy1, and esteptable
      %  esteptable = zeros(obj.length,1);
      %  [energies,~,esteptable] = unique(obj.depend{1},'rows'); % consider using legacy
      %  energy0 = obj.depend{1}(1,:);
      %  energy1 = obj.depend{1}(2,:);
      %end
      if isequal(obj.depend{1}(1,1),obj.depend{1}(2,1))
        irf_log('proc','PDist only has one set of energies, returning original PDist.')
        PD = obj;
        return
      end

      [pdistr,phir,energyr] = mms.psd_rebin(obj,TSeries(obj.time,obj.depend{2}),obj.ancillary.energy0,obj.ancillary.energy1,TSeries(obj.time,obj.ancillary.esteptable));
      PD = obj.clone(pdistr.time,pdistr.data);
      PD.depend{1} = repmat(energyr,PD.length,1);
      PD.ancillary.energy = PD.depend{1};
      PD.depend{2} = phir.data;

      PD.ancillary.dt_minus = obj.ancillary.dt_minus*2;
      PD.ancillary.dt_plus = obj.ancillary.dt_plus*2;

      % update delta_energy
      if isfield(PD.ancillary,'delta_energy_minus') && isfield(PD.ancillary,'delta_energy_plus')
        delta_energy = diff(energyr);
        log_energy = log10(energyr);
        log10_energy = diff(log_energy);
        log10_energy_plus  = log_energy + 0.5*[log10_energy log10_energy(end)];
        log10_energy_minus = log_energy - 0.5*[log10_energy(1) log10_energy];
        energy_plus = 10.^log10_energy_plus;
        energy_minus = 10.^log10_energy_minus;
        delta_energy_plus = energy_plus - energyr;
        delta_energy_minus = abs(energy_minus - energyr);
        delta_energy_plus(end) = max(PD.ancillary.delta_energy_minus(:,end));
        delta_energy_minus(1) = min(PD.ancillary.delta_energy_minus(:,1));
        PD.ancillary.delta_energy_plus = delta_energy_plus;
        PD.ancillary.delta_energy_minus = delta_energy_minus;
      end
      if isfield(PD.ancillary,'energy0')
        PD.ancillary.energy0 = PD.depend{1}(1,:);
        PD.ancillary.energy1 = PD.depend{1}(1,:);
      end
      if isfield(PD.ancillary,'esteptable'); PD.ancillary.esteptable = zeros(PD.length,1); end
    end
    function PD = interp(obj,varargin)
      % PDIST.INTERP Interpolates from spherical to spherical or cartesian
      % grid.
      %
      % To add: Interpolation between all possible grids.

      % Default values
      % Check input
      nargs = numel(varargin);
      have_options = 0;
      if nargs > 0, have_options = 1; args = varargin(:); end
      while have_options
        l = 0;
        switch(lower(args{1}))
          case 'cart'
            scpot = varargin{2};
            doCart = 1;
            l = 3;
            new_vxyz = args{2};
            new_vx_unit = args{2}(1,:);
            new_vy_unit = args{2}(2,:);
            new_vz_unit = args{2}(3,:);
            new_vbins_edges = args{3}; % edges of new bin
            args = args(l+1:end);
          otherwise
            l = 1;
            irf.log('warning',sprintf('Input ''%s'' not recognized.',args{1}))
            args = args(l+1:end);
        end
        if isempty(args), break, end
      end

      % Get v_xyz_DSL of original grid
      [old_vx,old_vy,old_vz] = obj.v;
      old_f = obj.data;


      % Rotate old coordinates into new coordinates
      old_vx_in_new_vxyz = old_vx*new_vx_unit(1) + old_vy*new_vx_unit(2) + old_vz*new_vx_unit(3);
      old_vy_in_new_vxyz = old_vx*new_vy_unit(1) + old_vy*new_vy_unit(2) + old_vz*new_vy_unit(3);
      old_vz_in_new_vxyz = old_vx*new_vz_unit(1) + old_vy*new_vz_unit(2) + old_vz*new_vz_unit(3);


      % Set up new grid
      % Add one outer bin
      [new_vx,new_vy,new_vz] = meshgrid(new_vbins_edges{1},new_vbins_edges{2},new_vbins_edges{3});
      new_f = nan(size(old_f));


      if 0
        %%
        scatter3(old_vx(:),old_vy(:),old_vz(:),old_vx(:)*0+1,old_vx(:))
        figure; scatter3(old_vx_in_new_vxyz(:),old_vy_in_new_vxyz(:),old_vz_in_new_vxyz(:),old_vx(:)*0+1,old_vx(:))

      end

      % Interpolate data from old to new grid
      sizeData = size(old_f);
      %%
      for itime = 1:length(obj.time)
        X = squeeze(old_vy_in_new_vxyz(itime,:,:,:));
        Y = squeeze(old_vy_in_new_vxyz(itime,:,:,:));
        Z = squeeze(old_vz_in_new_vxyz(itime,:,:,:));
        V = squeeze(old_f(itime,:,:,:));
        Xq = new_vx(:,:,:);
        Yq = new_vy(:,:,:);
        Zq = new_vz(:,:,:);
        Vq = interp3(X,Y,Z,V,Xq,Yq,Zq);
        new_f(itime,:,:,:) = Vq;
      end

    end
    function m = mass(obj)
      % Get mass of species
      units = irf_units;
      switch obj.species
        case {'e','electrons','electron'}
          m = units.me;
        case {'i','p','ions','ion','hplus'}
          m = units.mp;
        case {'oplus'}
          m = 16*units.mp;
        otherwise
          error('Species not supported.')
      end
    end
    function particles = macroparticles(obj,varargin)
      % PDIST.MACROPARTICLES Creates array of macro particles based on
      %   PDIST skymap.
      %
      %   particles = PD.macroparticles(inp1,arg1,...);
      %   particles = macroparticles(PD,inp1,arg1,...);
      %
      %   Input:
      %      PD - PDist object. Currently only 'skymap' is implemented.
      %
      %   Options:
      %     'skipzero',0 or 1 - skip bins that have zero phase space
      %             density, default is 1.
      %     'scpot',scpot - spacecraft potential in TSeries format, adjust energy
      %     'ntot',ntot - Total number of macro particles, scalar integer.
      %             If less than number of bins with non-zero phase space
      %             density, it is adjusted so that each cell with non-zero
      %             phase space density gets one particle.
      %     'nbin,nbin - Number of macro particles per bin.
      %         Of ntot and nbin, default is nbin with one particles per
      %         cell. If both are given as input, last one applies.
      %     'positioning' - 'random' or 'center', position of particles
      %             within each cell, default is 'randomize'
      %
      %   Output:   1×obj.length struct array with fields:
      %             iDep1, iDep2, ..., iDepN, dn, vx, vy, vz
      %             where N is the numer of dependent variables of PD.
      %             N = 3 for skymap (energy, azimuthal angle, polar angle)
      %             iDep give the index of the corresponding bin of the
      %             macro particle
      %             dn - particle density of macro particle,
      %                  dn = f(iDep1,iDep2,iDep3)*vol(iDep1,iDep2,iDep3)/Nmacro(iDep1,iDep2,iDep3)
      %                  sum(particle(1).dn) should give the density for
      %                  this time step, however, the absence of proper
      %                  calibration often makes them differ from FPI
      %                  moments.
      %
      %   Example:
      %     particles = PD(1).macroparticles('ntot',5e3,'skipzero',1,'scpot',scpot);
      %     scatter3(particles.vx,particles.vy,particles.vz,5,particles.dn)
      %

      units = irf_units;

      % Default values
      Ntot = 16*32*32; % one per cell if input is skymap
      doNtot = 0;
      Nbin = 1; % one per cell
      doNbin = 1;
      doSkipZero = 1;
      doScpot = 0;
      method = 'random'; % other option is center (of the bin)

      % Check input
      nargs = numel(varargin);
      have_options = 0;
      if nargs > 0, have_options = 1; args = varargin(:); end
      while have_options
        l = 0;
        switch(lower(args{1}))
          case 'ntot'
            Ntot = args{2};
            doNtot = 1;
            doNbin = 0;
            l = 2;
            args = args(l+1:end);
          case 'nbin'
            Nbin = args{2};
            doNbin = 1;
            doNtot = 0;
            l = 2;
            args = args(l+1:end);
          case 'skipzero'
            doSkipZero = args{2};
            l = 2;
            args = args(l+1:end);
          case 'scpot'
            doScpot = 1;
            l = 2;
            scpot = args{2};
            scpot = scpot.resample(obj); % make sure they have the same timeline
            args = args(l+1:end);
          case 'positioning'
            method = args{2};
            l = 2;
            args = args(l+1:end);
          otherwise
            l = 1;
            irf.log('warning',sprintf('Input ''%s'' not recognized.',args{1}))
            args = args(l+1:end);
        end
        if isempty(args), break, end
      end

      % Data from PDist in spherical coordinate system
      sizedata = obj.datasize;

      switch method
        case 'center' % add particles to the center of the bin
          % v_xyz_DSL of original grid
          if doScpot
            [vx,vy,vz] = obj.v('scpot',scpot);
          else
            [vx,vy,vz] = obj.v;
          end

          % Phase space density of each cell, can be different units depending
          % on original PDist
          f = obj.data;
          % Phase space volume of each cell, same base length and time units as f
          if doScpot
            vol = obj.d3v('scpot',scpot).data;
          else
            vol = obj.d3v.data;
          end
          % Partial density of each cell, same length unit as f
          dn = f.*vol;

          % Partial density for each macroparticle
          dn_part = dn/Nbin;

          Nbin_mat = Nbin*ones(sizedata);

          % Do not make any particles if the phase space density is zero
          if doSkipZero
            Nbin_mat(dn_part==0) = 0;
          end

          % Repeat values the number of times Nbin_mat specifies
          for it = 1:obj.length
            dn_all = repelem(dn_part(it,:),Nbin_mat(it,:));
            vx_all = repelem(vx(it,:),Nbin_mat(it,:));
            vy_all = repelem(vy(it,:),Nbin_mat(it,:));
            vz_all = repelem(vz(it,:),Nbin_mat(it,:));
            f_all = repelem(f(it,:),Nbin_mat(it,:));

            [iDep1,iDep2,iDep3] = ndgrid(1:sizedata(2),1:sizedata(3),1:sizedata(4));
            iDep1_all = repelem(iDep1(:),Nbin_mat(it,:));
            iDep2_all = repelem(iDep2(:),Nbin_mat(it,:));
            iDep3_all = repelem(iDep3(:),Nbin_mat(it,:));

            % Collect results into structure
            p(it).iDep1 = iDep1_all;
            p(it).iDep2 = iDep2_all;
            p(it).iDep3 = iDep3_all;
            p(it).dn = tocolumn(dn_all);
            p(it).vx = tocolumn(vx_all);
            p(it).vy = tocolumn(vy_all);
            p(it).vz = tocolumn(vz_all);
          end
        case 'random' % initialize random positions within each bin
          % First calculate how many particles should g in each bin
          % Phase space density of each cell, can be different units depending
          % on original PDist
          f = obj.data;
          % Phase space volume of each cell, same base length and time units as f
          if doScpot
            vol = obj.d3v('scpot',scpot).data;
          else
            vol = obj.d3v.data;
          end

          % Partial density of each cell, same length unit as f
          dn = f.*vol;

          if doNtot
            % Total density of each timestep, something wrong here, or just badly
            % messed up by background noise and particle contamination?
            n_tot = sum(dn(:,:),2);

            % Fraction of density in each separate bin
            n_frac = dn./repmat(n_tot,[1 sizedata(2:end)]);

            % Target number of particles per bin
            Ntmp = n_frac*Ntot; % sum(Ntmp(:)) = Ntot
            % Round up
            Ntmp_roundup = ceil(Ntmp);
            if not(doSkipZero)
              Ntmp_roundup(Ntmp_roundup==0) = 1;
            end
          elseif doNbin
            Ntmp_roundup = Nbin*ones(sizedata);
          end

          % Partial density for each macroparticle
          dn_part = dn./Ntmp_roundup;

          % n_frac = 0 divided by Ntmp_roundup = 0 gives NaN
          dn_part(isnan(dn_part)) = 0;

          % Edges of energy bins, same for each time step
          energy_minus = obj.depend{1}(1,:) - obj.ancillary.delta_energy_minus;
          energy_plus = obj.depend{1}(1,:) + obj.ancillary.delta_energy_plus;
          if doScpot
            scpotmat = repmat(scpot.data,[1 sizedata(2:end)]);
            energy_minus = energy_minus - scpotmat;
            energy_plus = energy_plus - scpotmat;
          end
          denergy = energy_plus - energy_minus;
          energy_edges = [energy_minus(:,:) energy_plus(:,end)];

          % Edges of polar angle bins, same for each time step
          polar_center = obj.depend{3};
          dpolar = polar_center(2) - polar_center(1);
          polar_minus = polar_center - 0.5*dpolar;

          for it = 1:obj.length
            % Initialize arrays for particles
            iDep1_all = zeros(sum(Ntmp_roundup(it,:),2),1);
            iDep2_all = zeros(sum(Ntmp_roundup(it,:),2),1);
            iDep3_all = zeros(sum(Ntmp_roundup(it,:),2),1);
            vx_all = zeros(sum(Ntmp_roundup(it,:),2),1);
            vy_all = zeros(sum(Ntmp_roundup(it,:),2),1);
            vz_all = zeros(sum(Ntmp_roundup(it,:),2),1);
            dn_all = zeros(sum(Ntmp_roundup(it,:),2),1);

            i_part_count = 1;

            % Edges of azimuthal angle bins, changes for each time step
            azim_center = obj.depend{2}(it,:);
            dazim = azim_center(2) - azim_center(1);
            azim_minus = azim_center-0.5*dazim;

            % Create N particles within each bin that each recieve 1/N of
            % the phase space density. These are then rotated into the new
            % coordinate system and binned in the new grid.
            for iEnergy    = 1:sizedata(2)
              % If scpot is inside bin, skip entire bin
              %if energy_minus(it,iEnergy) < 0, continue, end
              for iAzim    = 1:sizedata(3)
                for iPolar = 1:sizedata(4)
                  N_bin = Ntmp_roundup(it,iEnergy,iAzim,iPolar);
                  % Possiblity (defined by doSkipZero) to skip bin if space space density is zero
                  if doSkipZero && dn_part(it,iEnergy,iAzim,iPolar) == 0, continue; end

                  % Assign N_bin randomized positions within energy and
                  % angle ranges
                  if iAzim == 32 % debug
                    1;
                  end
                  tmp_energy = energy_minus(it,iEnergy) + denergy(it,iEnergy)*rand(N_bin,1); % eV
                  tmp_azim   = azim_minus(iAzim)        + dazim*rand(N_bin,1);   % deg
                  tmp_polar  = polar_minus(iPolar)      + dpolar*rand(N_bin,1);  % deg
                  tmp_v = sqrt(tmp_energy*units.eV*2/units.me)/1000; % km/s

                  if doScpot % check if energy is negative, then skip
                    if iEnergy>8;  % debug
                      1;
                    end
                    tmp_azim(tmp_energy<0) = [];
                    tmp_polar(tmp_energy<0) = [];
                    tmp_v(tmp_energy<0) = [];
                    if isempty(tmp_v) continue; end
                    N_bin = numel(tmp_v);
                  end

                  % Transform into cartesian velocity components

                  tmp_vx = -tmp_v.*sind(tmp_polar).*cosd(tmp_azim); % '-' because the data shows which direction the particles were coming from
                  tmp_vy = -tmp_v.*sind(tmp_polar).*sind(tmp_azim);
                  tmp_vz = -tmp_v.*cosd(tmp_polar);
                  if any(tmp_vy>0) % debug
                    1;
                  end
                  if iPolar == 16 % debug
                    1;
                  end

                  %if tmp_polar

                  % Rotate into new coordinate system
                  %new_vx = old_vx*new_vx_unit(1) + old_vy*new_vx_unit(2) + old_vz*new_vx_unit(3);
                  %new_vy = old_vx*new_vy_unit(1) + old_vy*new_vy_unit(2) + old_vz*new_vy_unit(3);
                  %new_vz = old_vx*new_vz_unit(1) + old_vy*new_vz_unit(2) + old_vz*new_vz_unit(3);

                  % Assign particle density to each macro particle
                  tmp_dn = repelem(dn_part(it,iEnergy,iAzim,iPolar),N_bin);
                  tmp_iDep1 = repelem(iEnergy,N_bin);
                  tmp_iDep2 = repelem(iAzim,N_bin);
                  tmp_iDep3 = repelem(iPolar,N_bin);

                  % Collect particles in array
                  iDep1_all(i_part_count + (0:N_bin-1),:) = tmp_iDep1;
                  iDep2_all(i_part_count + (0:N_bin-1),:) = tmp_iDep2;
                  iDep3_all(i_part_count + (0:N_bin-1),:) = tmp_iDep3;
                  vx_all(i_part_count + (0:N_bin-1),:) = tmp_vx;
                  vy_all(i_part_count + (0:N_bin-1),:) = tmp_vy;
                  vz_all(i_part_count + (0:N_bin-1),:) = tmp_vz;
                  dn_all(i_part_count + (0:N_bin-1),:) = tmp_dn;

                  % Increase counter
                  i_part_count = i_part_count + N_bin;
                end % end polar angle loop
              end % end azimuthal angle loop
            end % end energy loop
            p(it).iDep1 = iDep1_all(1:i_part_count-1);
            p(it).iDep2 = iDep2_all(1:i_part_count-1);
            p(it).iDep3 = iDep3_all(1:i_part_count-1);
            p(it).vx = vx_all(1:i_part_count-1);
            p(it).vy = vy_all(1:i_part_count-1);
            p(it).vz = vz_all(1:i_part_count-1);
            p(it).dn = dn_all(1:i_part_count-1);
          end % end time loop
      end % end switch method
      particles = p;
    end
    %     function e = energy(obj)
    %       % Get energy of object when not knowing its index
    %       isEnergy = cellfun(@(s) strcmp(s,'energy'),obj.representation);
    %       iEnergy = find(isEnergy);
    %       if isempty(iEnergy) % no energy dependence, return
    %         e = [];
    %       else
    %         e = obj.depend{iEnergy};
    %       end
    %     end
    function moms = moments(obj,varargin)
      % PRELIMINARY VERSION
      % Currently does not include spacecraft potential and has no
      % pressure tensor or heat flux.
      %
      % MOMENTS get particle moments from PDist object
      %
      %   moms = PDIST.MOMENTS returns a structure containing TSeries
      %   objects of the moments recalculated from the PDist object.
      %       Moments:
      %           n   -   number density [cm^-3]
      %           V   -   bulk velocity [km/s]
      %           T   -   Temperature tensor [eV]
      %
      %   See also: MMS.PSD_MOMENTS
      %
      %   TODO:   -Implement heat flux and maybe pressure tensor.
      %           -Spacecraft potential should really NOT be an input to
      %           this function but should be called in a separate class
      %           method (maybe it already exists?)
      %

      % make sure it's PSD and in SI units
      dist = obj.convertto('s^3/m^6');

      % units
      u = irf_units;
      % particle mass
      if strcmp(dist.species,'electrons'); isDes = 1; else, isDes = 0; end
      if isDes; M = u.me; else; M = u.mp; end

      % get intstrument values (azimuthal angle is set in loop)
      % elevation angle
      th = double(dist.depend{3}); % polar angle in degrees
      th = th-90; % elevation angle in degrees
      th = th*pi/180; % in radians
      dth = median(diff(th)); % scalar
      % velocity
      if ~isempty(dist.ancillary)
        esteptable = dist.ancillary.esteptable;
      else % if there is no ancillary data, it's assumed to be fast mode
        % set steptable to all zeros
        esteptable = zeros(1,length(dist));
      end
      idEstep0First = find(esteptable==0,1);
      emat = double(dist.depend{1}); % [eV]
      e0 = emat(idEstep0First,:); e0 = e0(1,:); % [eV]
      e1 = emat(idEstep0First+1,:); e1 = e1(1,:); % [eV]
      % velocity (size = [2,nE]) per steptable
      % this works also for no energy table switching since e0 == e1
      v = sqrt((2*[e0;e1])*u.e/M); % [m/s]

      % velocity diffs from delta energy
      % energy diffs minus/plus
      if ~isempty(dist.ancillary)
        dEm = dist.ancillary.delta_energy_minus; % [eV]
        dEp = dist.ancillary.delta_energy_plus; % [eV]
      else % if there is no ancillary data, it's assumed to be fast mode
        % display warning
        irf.log('w','No information on delta_energy_minus/plus in PDist, guessing values')
        % guess energy diffs
        %dEm = repmat([e0(2)-e0(1),diff(e0)],length(dist),1)/2;
        dEm = ([diff(e0),e0(end)-e0(end-1)]+[e0(2)-e0(1),diff(e0)])/2;
        dEm = repmat(dEm,length(dist),1)/2;
        dEp = dEm;
      end
      % per esteptable
      dEm0 = dEm(idEstep0First,:);
      dEp0 = dEp(idEstep0First,:);
      dEm1 = dEm(idEstep0First+1,:);
      dEp1 = dEp(idEstep0First+1,:);
      % vel diffs
      v0lower = sqrt(2*(e0-dEm0)*u.e/M); % [m/s]
      v0upper = sqrt(2*(e0+dEp0)*u.e/M);
      v1lower = sqrt(2*(e1-dEm1)*u.e/M);
      v1upper = sqrt(2*(e1+dEp1)*u.e/M);
      % same structure as for v
      dv = [v0upper-v0lower; v1upper-v1lower]; % [m/s]

      % Number of instrument bins
      nEle = length(th);
      nV = length(v);

      % initialize arrays
      N = zeros(1,obj.length);
      NV = zeros(3,obj.length);
      Pressure = zeros(3,3,obj.length);

      % loop'n through time
      for it = 1:obj.length
        % 3d data matrix for time index it, [E,phi,th]
        F3d = double(squeeze(dist.data(it,:,:,:)));

        % azimuthal angle
        if size(dist.depend{2},1)>1 % brst mode
          phi = double(dist.depend{2}(it,:)); % in degrees
        else % fast mode
          phi = double(dist.depend{2});
        end
        phi = phi-180; % travel/arrival correction
        phi = phi*pi/180; % in radians
        dphi = median(diff(phi)); % scalar

        % Number of instrument bins
        nAz = length(phi);


        % 3D matrices for instrumental bin centers
        TH = repmat(th,nV,1,nAz);                       % [v,th,phi]
        TH = permute(TH,[1,3,2]);                       % [v,phi,th]
        PHI = repmat(phi,nV,1,nEle);                    % [v,phi,th]
        VEL = repmat(v(esteptable(it)+1,:),nAz,1,nEle); % [phi,v,th]
        VEL = permute(VEL,[2,1,3]);                     % [v,phi,th]
        DV = repmat(dv(esteptable(it)+1,:),nAz,1,nEle); % [phi,v,th]
        DV = permute(DV,[2,1,3]);                       % [v,phi,th]


        [VX,VY,VZ] = sph2cart(PHI,TH,VEL);

        % density
        N(it) = sum(sum(sum(F3d.*VEL.^2.*DV.*cos(TH)*dphi*dth)));

        % should improve preformance by finding indices with non-zero
        % psd?
        % idf = find(F3d);
        % [idfV,idfPhi,idfTh] = ind2sub(size(F3d),idf);

        % mega loop (should skip empty bins)
        for iv = 1:nV
          for iphi = 1:nAz
            for ith = 1:nEle
              % Ignore bin if value of F is zero to save computations
              if F3d(iv,iphi,ith) == 0
                continue;
              end
              % velocity
              NV(:,it) = NV(:,it)+...
                [VX(iv,iphi,ith);VY(iv,iphi,ith);VZ(iv,iphi,ith)]*...
                (F3d(iv,iphi,ith)*VEL(iv,iphi,ith)^2*DV(iv,iphi,ith)*...
                cos(TH(iv,iphi,ith))*dphi*dth);

            end
          end
        end

        Vtemp = NV(:,it)/N(it);

        % mega loop #2 (should skip empty bins) to get pressure tensor
        % separate loop because it requires velocity moments
        for iv = 1:nV
          for iphi = 1:nAz
            for ith = 1:nEle
              % Ignore bin if value of F is zero to save computations
              if F3d(iv,iphi,ith) == 0
                continue;
              end
              % pressure (6.9)
              Pressure(:,:,it) = Pressure(:,:,it)+...
                M*(([VX(iv,iphi,ith);VY(iv,iphi,ith);VZ(iv,iphi,ith)]-Vtemp)*...
                ([VX(iv,iphi,ith);VY(iv,iphi,ith);VZ(iv,iphi,ith)]-Vtemp)')*...
                (F3d(iv,iphi,ith)*VEL(iv,iphi,ith)^2*DV(iv,iphi,ith)*...
                cos(TH(iv,iphi,ith))*dphi*dth);
            end
          end
        end
      end

      % set output structure
      moms = [];

      % density
      moms.n = irf.ts_scalar(obj.time,N*1e-6);
      moms.n.name = [dist.name,'_moms_density'];
      moms.n.units = 'cm^-3';
      moms.n.siConversion = '1e6>m^-3';

      % velocity
      moms.V = irf.ts_vec_xyz(obj.time,(NV./N)'*1e-3);
      moms.V.name = [dist.name,'_moms_velocity'];
      moms.V.units = 'km/s';
      moms.V.siConversion = '1.0e3>m s^-1';
      % tensorOrder, representation, etc are read-only, how to add?

      % temperature
      moms.T = irf.ts_tensor_xyz(obj.time,permute(Pressure,[3,1,2])./(u.e*N'));
      moms.T.name = [dist.name,'_moms_temperature'];
      moms.T.units = 'eV';
      moms.T.siConversion = '11604.50520>K';
      % tensorOrder, representation, etc are read-only, how to add?
    end


  end
  % Plotting and other functions
  methods (Access = protected)
    function PD = enforce_depend_timeseries(obj,depend)
      % Find if 'depend' is a depend, and if yes, enforce it to be a
      % timeseries
      PD = obj;
      isSep = cellfun(@(s) strcmp(s,depend),obj.representation);
      iDep = find(isSep);
      if isempty(iDep) % no energy dependence, return
        return;
      end
      current_depend = obj.depend{iDep};
      dimDependData = size(obj.data,iDep+1); % +1 because dim 1 is time
      if size(current_depend) == [obj.length,dimDependData] % depend{iEnergy} is timeseries
        return
      elseif size(current_depend) == [1,dimDependData] % depend{iEnergy} is 1 x nt
        PD.depend{iDep} = repmat(current_depend,[obj.length,1]);
      elseif size(current_depend) == [1,dimDependData] % depend{iEnergy} is nt x 1
        PD.depend{iDep} = repmat(current_depend',[obj.length,1]);
      end
    end

    % help function for PDist.omni
    function fmean = get_omni_dist_framedep(obj,V0,varargin)
      % GET_OMNI_DIST_FRAMEDEP Get omnidirectional distribution in an arbitrary frame
      %
      %   fOmni = GET_OMNI_DIST_FRAMEDEP(obj,V) returns the omnidirectional
      %   distribution of the PDist object dist in the frame moving with velocity
      %   V. V can either be a 1x3 vector for a frame moving with constant
      %   velocity or a TSeries object for a changing frame. The unit of the
      %   input distribution is conserved in the output. The unit of V is [km/s].
      %
      %   In mathematical terms, this function calculates the spherical mean of
      %   the distribution function https://en.wikipedia.org/wiki/Spherical_mean
      %   centered on the velocity in the input. It accomplishes this with a
      %   Monte-Carlo method where a number of points (5e3 per channel) evenly
      %   distributed on a sphere centered on V and with radius corresponding to
      %   each energy value v = sqrt(2*E/m). The omnidirectional
      %   distribution is then the average value of the distribution the
      %   Monte-Carlo points for the given energy value.
      %
      % Tested but somewhat experimental.

      % TODO:
      %   - More inputs: energy grid, number of MC points, mass, ...
      %   - Add more info to the output PDist object?
      %   - Allow for varying dPhi, dTh
      %   - Remove flags for input, better to force all things to be in
      %   input

      % handle input
      args = varargin;
      nargs = length(args);

      inpE = 0; % default use obj energy table
      inpPhi = 0; % default assume equally-sized bins in angle
      inpTh = 0;
      have_options = nargs > 1;
      while have_options
        switch(lower(args{1}))
          case 'e' % Energy
            Ef = args{2};
            inpE = 1;
            args = args(3:end);
            % more cases here in the future
          otherwise
            irf.log('w','Unknown input')

        end
        if isempty(args), break, end
      end

      % handle angles
      % check if obj has dphi, dth
      if isfield(obj.ancillary,'delta_phi_plus')
        % fine to have both plus and minus in if statement because only
        % having one makes no sense
        dPhi_plus = obj.ancillary.delta_phi_plus*pi/180; % in radians
        dPhi_minus = obj.ancillary.delta_phi_minus*pi/180;
        inpPhi = 1;
      end
      if isfield(obj.ancillary,'delta_theta_plus')
        dTh_plus = obj.ancillary.delta_theta_plus*pi/180; % in radians
        dTh_minus = obj.ancillary.delta_theta_minus*pi/180;
        inpTh = 1;
      end


      % units and constants
      u = irf_units;

      % number of Monte-Carlo points
      nMC = 5e3; % per energy channel

      % handle inputs
      if nargin == 1
        V0 = [0,0,0];
      end

      % check species (user should be able to input M)
      switch obj.species
        case 'ions'
          M = u.mp;
        case 'electrons'
          M = u.me;
      end

      nt = length(obj);

      % Handle input velocity
      % make sure V0 is an array with size PDist.length
      if isa(V0,'TSeries')
        %if size(V0.data,1) == size(PDist.data,1)
        if V0.length == obj.length
          irf.log('w','assuming same time sampling')
          V0v = V0.data;
        else
          irf.log('w','resampling velocity')
          V0v = V0.resample(obj).data;
        end
      elseif isnumeric(V0)
        irf.log('w','using constant velocity')
        if size(V0,1)==3
          V0v = repmat(V0,1,n);
        else
          V0v = repmat(V0,nt,1); % hopefully right
        end
      end

      % for good measure
      V0v = double(V0v);

      % pre treat distribution object
      % start with the normal command so it has the same structure
      % fmean = obj.omni;

      % assume they are all the same
      E = obj.depend{1}; % [eV]
      % delta energies [eV]
      if isfield(obj.ancillary,'delta_energy_plus')
        irf.log('w','PDist ancillary information present')
        dEp = obj.ancillary.delta_energy_plus;
        dEm = obj.ancillary.delta_energy_minus;
      else
        irf.log('w','PDist ancillary information NOT present')
        dE = diff(E')';
        dEp = [dE,dE(:,end)]/2;
        dEm = [dE(:,1),dE]/2;
      end
      % energy edges
      Ee = [E-dEm,E(:,end)+dEp(:,end)]; % [eV]
      ve = sqrt(2*Ee*u.e/M)*1e-3; % [km/s]

      % elevation angle
      th = double(obj.depend{3}); % polar angle in degrees
      th = th-90; % elevation angle in degrees
      th = th*pi/180; % in radians
      % SolO PDists don't follow MMS convention and th angles are
      % decreasing, this causes errors so they must be flipped
      if ~issorted(th)
        % assume it's decreasing, if it jumps around the pdist object should feel bad
        th = fliplr(th);
        flipTheta = 1; % flag to flip distribution function
      else
        flipTheta = 0;
      end
      if inpTh
        % ensure row vector
        the = [th-dTh_minus(:)',th(end)+dTh_plus(end)]; % [radians]
      else
        dth = median(diff(th)); % assume same size
        the = [th-dth/2,th(end)+dth/2]; % [radians]
      end

      % If not set, energy table is the same as first in input distribution
      if ~ inpE
        Ef = E(1,:); % [eV]
      end
      nEf = length(Ef);
      % velocity of grid in [km/s]
      vf = sqrt(2*Ef*u.e/M)*1e-3; % [km/s]

      % pre allocate data matrix
      % nE = size(fmean.data,2);
      fmeanData = zeros(nt,nEf);

      % create mc particles
      % acceptence-rejection method (not so slow apparently)
      % this creates elevation angle points distributed like a cosine
      thp = zeros(1,nMC);
      count = 1;
      while count <= nMC
        rn1 = rand(1); % [0,1] flat
        rn1 = (rn1-0.5)*pi; % [-pi/2,pi/2] flat

        rn2 = rand(1);
        if rn2<abs(cos(rn1))
          % accept
          thp(count) = rn1;
          count = count+1;
        end % else reject and try again
      end

      phip = (rand(1,nMC)-0.5)*2*pi; % [-pi,pi] flat

      % velocities in spacecraft frame
      vxp = zeros(nEf,nMC);
      vyp = zeros(nEf,nMC);
      vzp = zeros(nEf,nMC);

      % use the same angle values for all energy values
      for ii = 1:nEf
        [vxp(ii,:),vyp(ii,:),vzp(ii,:)] = sph2cart(phip,thp,vf(ii));
      end

      % get spherical mean of distribution
      % loop over time steps
      fprintf('it = %4.0f/%4.0f\n',0,nt) % display progress
      for it = 1:nt
        if mod(it,1) == 0, fprintf([repmat('\b', 1, 10) '%4.0f/%4.0f\n'],it,nt); end % display progress

        % ---------------- instrument part ----------------
        % psd
        F3d = double(squeeze(double(obj.data(it,:,:,:)))); % whatever units

        % remember to flip the distribution if th was the wrong order
        if flipTheta
          F3d = flip(F3d,3); % hopefully not slow
        end

        % azimuthal angle
        if size(obj.depend{2},1) > 1 % burst mode mms
          phi = double(obj.depend{2}(it,:)); % in degrees
        else % fast mode mms or SolO
          phi = double(obj.depend{2});
        end
        phi = phi-180;
        phi = phi*pi/180; % in radians (can be outside [-pi,pi])

        % check phi angles (hopefully not slow)
        if ~issorted(phi)
          % assume it's decreasing, if it jumps around the pdist object should feel bad
          phi = fliplr(phi);
          F3d = flip(F3d,2); % hope and pray
        end

        % special case for SolO where phi at this stage can be less than
        % -pi because azimuth is defined between [-180,180] degrees unlike
        % for MMS where it is defined in [0,360] degrees.
        % If this happens put a flag to fix the angles in the -x,+y (second)
        % quadrant
        if numel(find(phi<-pi))>0
          if ~numel(find(phi>pi/2))==0 % can't have more than 2 pi coverage
            error('something is wrong with the azimuth angles')
          else
            fix2ndQuadrant = 1;
          end
        else
          fix2ndQuadrant = 0;
        end

        % edges of phi bins
        if inpPhi
          phie = [phi-dPhi_minus(:)',phi(end)+dPhi_plus(end)];
        else
          dphi = median(diff(phi));
          phie = [phi-dphi/2,phi(end)+dphi/2];
        end

        % ---------------- grid part ----------------
        % velocities of mc points in desired frame
        vxp0 = vxp+V0v(it,1);
        vyp0 = vyp+V0v(it,2);
        vzp0 = vzp+V0v(it,3);

        % back to the spherical frame
        [phip0,thp0,vp0] = cart2sph(vxp0,vyp0,vzp0);

        % phip0 is now distributed in [-pi,pi]
        % if angles in obj are messed up, redefine phip0 to [-3*pi/2,pi/2]
        % (SO ELEGANT!!!)
        if fix2ndQuadrant
          phip0(phip0>pi/2) = phip0(phip0>pi/2)-2*pi;
        end


        % get good indices
        idPhip = discretize(phip0,phie);
        idThp = discretize(thp0,the);
        idVp = discretize(vp0,ve(it,:));

        % ---------------- calculating mean psd ----------------
        for jj = 1:nEf
          % get the instrument bin indices of all MC points
          idMC = sub2ind(size(F3d),idVp(jj,:),idPhip(jj,:),idThp(jj,:));
          % sometimes there is no instrument bin corresponding to the MC
          % point, those indices become NaNs but should count as zero psd
          idMC = idMC(~isnan(idMC));
          fmeanData(it,jj) = sum(F3d(idMC))/nMC;
        end
      end

      % construct object
      fmean = PDist(obj.time,fmeanData,'omni',Ef);


      fmean.ancillary.V0 = V0v; % add velocity to ancillary

      % set useful things
      fmean.representation = obj.representation(1);
      fmean.units = obj.units;
      fmean.name = 'omni';
      fmean.units = obj.units;


    end
  end
  methods (Static)
    function newUnits = changeunits(from,to)

    end
  end
end

