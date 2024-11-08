classdef TSeries
  %TSeries Generic time dependent variable
  %   Time varibale has 2 fields: T - time [GenericTimeArray] and DATA
  %
  %   TS = TSeries(T,DATA,[ARGS])
  %
  %   ARGS:
  %      'to','TensorOrder' - 0 (scalar-default), 1 (vector), 2 (tensor)
  %      'tb','TensorBasis' -
  %                           xyz (Cartesian)
  %                           rtp (Spherical,colatitude)
  %                           rlp (Spherical,latitude)
  %                           rpz (Cylindrical)
  %                           xy (Cartesian 2D)
  %                           tp (Polar 2D)
  %      'repres' - Representation for each of the dimensions of the DATA
  %      as a cell array containing description of tensor dimensionality
  %      or {} if the DATA dimension correspond to variables associated
  %      with other independent variables of the data product, such and a
  %      particle energy or a spectral frequency channel.
  %
  %      For more info on representation see Cluster Metadata Dictionary
  %           https://caa.esac.esa.int/caa/doc_format_meta.xml
  %
  %  ARGS MACROs:
  %      'vec_xyz','vec_rtp','vec_rlp','vec_rpz' - 3D vectors
  %      'vec_xy','vec_rp'                       - 2D vectors
  %
  %  Example:
  %
  %  epoch = EpochUnix(iso2epoch('2002-03-04T09:30:00Z')+(0:3));
  %  data4x2 = [1 2; 3 4; 5 6; 7 8];
  %  data4x3 = [1 2 3; 4 5 6; 7 8 9; 10 11 12];
  %  data4x3x3 = rand(4,3,3);
  %
  %  % TensorOrder=1, 2 components of 2D vector
  %  TsVec2DXY = TSeries(epoch,data4x2,'TensorOrder',1,'TensorBasis','xy',...
  %          'repres',{'x','y'})
  %
  %  % Or using equvivalent macro
  %  TsVec2DXY = TSeries(epoch,data4x2,'vec_xy')
  %
  %  % TensorOrder=1, 3 components of 3D vector
  %  TsVec3DRTP = TSeries(epoch,data4x3,'vec_rtp')
  %
  %  % TensorOrder=1, 2 components (x,z) of incomplete 3D vector
  %  TsVec3DXZ = TSeries(epoch,data4x2,'TensorOrder',1,'TensorBasis','xyz',...
  %          'repres',{'x','z'})
  %
  %  % TensorOrder=2, 3x3 components of tensor (3D, default basis=xyz)
  %  TsTen3D = TSeries(epoch,data4x3x3,'TensorOrder',2,...
  %          'repres',{'x','y','z'},'repres',{'x','y','z'})
  %
  %  % TensorOrder=1, 3 energies x 3 components of vector (3D, spherical)
  %  TsPflux = TSeries(epoch,data4x3x3,'TensorOrder',1,'TensorBasis','rtp',...
  %          'repres',{},'repres',{'r','t','p'})


  % ----------------------------------------------------------------------------
  % SPDX-License-Identifier: Beerware
  % "THE BEER-WARE LICENSE" (Revision 42):
  % <yuri@irfu.se> wrote this file.  As long as you retain this notice you
  % can do whatever you want with this stuff. If we meet some day, and you think
  % this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
  % ----------------------------------------------------------------------------

  properties (Access=protected)
    coordinateSystem_ = '';
    data_
    t_ % GenericTimeArray
    fullDim_
    tensorOrder_ = 0;
    tensorBasis_ = '';
  end

  properties (Dependent = true)
    data
    time
    coordinateSystem
  end

  properties (SetAccess = immutable,Dependent = true)
    tensorOrder
    tensorBasis
  end

  properties (Constant = true, Hidden = true)
    MAX_TENSOR_ORDER = 2;
    BASIS = {'xyz','rtp','rlp','rpz','xy','rp'};
    BASIS_NAMES = {...
      'Cartesian','Spherical,colatitude', 'Spherical,latitude','Cylindrical',...
      'Cartesian 2D','Polar 2D'};
  end

  properties (SetAccess = protected)
    representation
  end

  properties
    name = '';
    units = '';
    siConversion = '';
    userData = [];
  end

  methods
    function obj = TSeries(t,data,varargin)
      if nargin == 0, obj.data_ = []; obj.t_ = []; return, end
      if nargin == 1 && isempty(t), obj.data_ = []; obj.t_ = []; return, end
      if nargin == 2 && or(isempty(t),isempty(data)), obj.data_ = []; obj.t_ = []; return, end

      if nargin<2, error('2 inputs required'), end
      if ~isa(t,'GenericTimeArray')
        error('irf:TSeries:TSeries:badInputs',...
          'T must be of TSeries type or derived from it')
      end
      if ~isnumeric(data)
        error('irf:TSeries:TSeries:badInputs',...
          'DATA must be numeric')
      end
      if size(data,1) ~= t.length()
        error('irf:TSeries:TSeries:badInputs',...
          'T and DATA must have the same number of records')
      end
      obj.data_ = data; obj.t_ = t;
      obj.fullDim_ = cell(ndims(data)-1,1); % time should not be included -> -1
      obj.representation = cell(ndims(data)-1,1); iDim = 1; % time should not be included -> -1
      %obj.representation = cell(ndims(data),1); iDim = 1;
      % Why is the dimension of represenation set already here? When
      % picking out a component of a Tensor: T.xy, the data becomes 2D
      % (time x T.xy), but the representation should have two 'dimensions':
      % 'x' and 'y'. This can be solved by resetting representation when
      % specifying tensor order, since tensor order hase to be specified
      % before settin representation.

      % remove since time should not be included
      %if ~isempty(obj.t_)
      %  obj.representation{iDim} = obj.t_([]); iDim = iDim + 1;
      %end

      args = varargin;
      flagTensorOrderSet = false; flagTensorBasisSet = false;
      while ~isempty(args)
        x = args{1}; args(1) = [];
        switch lower(x)
          case {'tensor_xyz'}
            if ndims(obj.data_)>3 %#ok<ISMAT>
              error('irf:TSeries:TSeries:badInputs',...
                'DATA has more than 3 dimensions (needed for Tensor with order 2)')
            elseif size(obj.data_,2)~=3
              error('irf:TSeries:TSeries:badInputs',...
                'size(DATA,2) must be 3 for xyz tensor')
            elseif size(obj.data_,3)~=3
              error('irf:TSeries:TSeries:badInputs',...
                'size(DATA,3) must be 3 for xyz tensor')
            end
            obj.tensorOrder_ = 2; flagTensorOrderSet = true;
            [~,iB] = intersect(obj.BASIS,x(8:10));
            obj.tensorBasis_ = iB; flagTensorBasisSet = true;
            obj.representation{1} = {x(8), x(9), x(10)}; % 1st dimension, (rows for 2D)
            obj.fullDim_{1} = true;
            obj.representation{2} = {x(8), x(9), x(10)}; % 2nd dimension, (cols for 2D)
            obj.fullDim_{2} = true;
          case {'vec_xyz','vec_rtp','vec_rlp','vec_rpz'}
            if ndims(obj.data_)>2 %#ok<ISMAT>
              error('irf:TSeries:TSeries:badInputs',...
                'DATA has more than 2 dimensions (needed for 3D vec)')
            elseif size(obj.data_,2)~=3
              error('irf:TSeries:TSeries:badInputs',...
                'size(DATA,2) must be 3 for 3D vec')
            end
            obj.tensorOrder_ = 1; flagTensorOrderSet = true;
            [~,iB] = intersect(obj.BASIS,x(5:7));
            obj.tensorBasis_ = iB; flagTensorBasisSet = true;
            obj.representation{1} = {x(5), x(6), x(7)};
            obj.fullDim_{1} = true;
          case {'vec_xy','vec_rp'}
            if ndims(obj.data_)>2 %#ok<ISMAT>
              error('irf:TSeries:TSeries:badInputs',...
                'DATA has more than 2 dimentions (needed for 2D vec)')
            elseif size(obj.data_,2)~=2
              error('irf:TSeries:TSeries:badInputs',...
                'size(DATA,2) must be 2 for 2D vec')
            end
            obj.tensorOrder_ = 1; flagTensorOrderSet = true;
            [~,iB] = intersect(obj.BASIS,x(5:6));
            obj.tensorBasis_ = iB; flagTensorBasisSet = true;
            obj.representation{1} = {x(5), x(6)}; obj.fullDim_{1} = true;
          case {'to','tensororder'}
            if flagTensorOrderSet
              error('irf:TSeries:TSeries:badInputs',...
                'tensorOrder already set')
            end
            if ~isempty(args), y = args{1}; args(1) = [];
            else
              error('irf:TSeries:TSeries:badInputs',...
                'tensorOrder requires a second argument')
            end
            if ~isempty(intersect(y,(0:1:obj.MAX_TENSOR_ORDER)))
              obj.tensorOrder_ = y; flagTensorOrderSet = true;
              if y>0
                % Check if we have any data dimension with 1..3 elements
                found = false;
                for i=2:ndims(obj.data_)
                  if size(obj.data_,i)<=3, found = true; break, end
                end
                if ~found
                  error('irf:TSeries:TSeries:badInputs',...
                    'cannot construct tensor:all data dimensions have size<=3')
                end
                obj.tensorBasis_ = 1; % set default basis to XYZ
              end
            else
              error('irf:TSeries:TSeries:badInputs',...
                'tensorOrder must be 0<=tensorOrder<=%d',...
                obj.MAX_TENSOR_ORDER)
            end
            obj.representation = cell(obj.tensorOrder,1);
            iDim = 1;
            %if ~isempty(obj.t_)
            %  obj.representation{iDim} = obj.t_; iDim = iDim + 1;
            %end
          case {'tb','basis','tensorbasis'}
            if flagTensorBasisSet
              error('irf:TSeries:TSeries:badInputs',...
                'tensorBasis already set')
            end
            if obj.tensorOrder_==0
              error('irf:TSeries:TSeries:badInputs',...
                'cannot set tensorBasis for a scalar (tensorOrder=0)')
            end
            if ~isempty(args), y = args{1}; args(1) = [];
            else
              error('irf:TSeries:TSeries:badInputs',...
                'tensorBasis requires a second argument')
            end
            [~,iB] = intersect(obj.BASIS,y);
            if ~isempty(iB)
              obj.tensorBasis_ = iB; flagTensorBasisSet = true;
            else
              error('irf:TSeries:TSeries:badInputs',...
                'tensorBasis value not recognized')
            end
          case {'rep','repres','representation'}
            if isempty(obj.tensorOrder_)
              error('irf:TSeries:TSeries:badInputs',...
                'Must specify TensorOrder first')
            end
            if iDim>(obj.tensorOrder) %ndims(data), % remove '+1' after obj.tensorOrder, since time should no longer be included
              error('irf:TSeries:TSeries:badInputs',...
                'Representation already set for all DATA dimensions')
            end
            if ~isempty(args), y = args{1}; args(1) = [];
            else
              error('irf:TSeries:TSeries:badInputs',...
                'Representation requires a second argument')
            end
            if isempty(y) && iDim<ndims(data)-1, iDim = iDim + 1;
            else
              [ok,msg] = validate_representation(y);
              if ~isempty(ok)
                obj.fullDim_{iDim}=ok; obj.representation{iDim} = y;
                iDim = iDim + 1;
              else
                error('irf:TSeries:TSeries:badInputs',msg)
              end
            end
          case 'depend'
          otherwise
            error('irf:TSeries:TSeries:badInputs',...
              'Unknown parameter')
        end
      end

      % XXX: TODO Validate if we have non-empty Representation set for
      % number of dimensions corresponding to tensor order

      function [ok,msg] = validate_representation(x)
        ok = []; msg = '';
        sDim = size(obj.data,iDim+1); tb = obj.BASIS{obj.tensorBasis_};
        if sDim>length(tb)
          msg = sprintf(...
            'Dimension %d size %d>%d (Basis=%s) cannot have Representation. Use Depend instead',...
            iDim,sDim,length(tb),tb);
          return
        end
        if ~iscell(x)
          msg = 'Representation requires a cell input'; return
        end
        if sDim>0 && sDim~=length(x)
          msg = sprintf('Representation requires a cell(%d) input',sDim);
          return
        end
        if ~all(cellfun(@(x) ischar(x) && length(x)==1,x))
          msg = sprintf(...
            'Representation requires char(%d) cells elements',...
            obj.tensorOrder_);
          return
        end
        s = sprintf('%s',char(x)');
        if length(s)~=length(unique(s))
          msg = sprintf(...
            'Representation (%s) contains repeating attributes',s);
          return
        end
        s = unique(s);
        c = length((intersect(tb,s)));
        if c == 0, msg = sprintf('Unrecognized representation');
        elseif length(s) == c, ok = true; % complete dim
        elseif length(s) < c, ok = false; % incomplete dim
        elseif length(s) > c
          d = setdiff(s,tb); s1 = sprintf('%s',char(d)');
          msg = sprintf(...
            'Unrecognized entries (%s) for representation ''%s''',s1,tb);
        else, error('should not be here')
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
          obj.data_ = obj.data_(idxTmp{:});
          if numel(idx) > 1
            obj = builtin('subsref',obj,idx(2:end));
          end
          [varargout{1:nargout}] = obj;
        case '{}'
          error('irf:TSeries:subsref',...
            'Not a supported subscripted reference')
      end
    end

    function value = get.coordinateSystem(obj)
      value = obj.coordinateSystem_;
    end

    function value = get.data(obj)
      value = obj.data_;
    end

    function value = get.time(obj)
      value = obj.t_;
    end

    function value = basis(obj)
      value = obj.BASIS{obj.tensorBasis_};
    end

    function value = get.tensorBasis(obj)
      value = [obj.BASIS{obj.tensorBasis_}...
        ' (' obj.BASIS_NAMES{obj.tensorBasis_} ')'];
    end

    function value = get.tensorOrder(obj)
      value = obj.tensorOrder_;
    end

    %Components
    function y = x(obj)
      %access X component
      [y,ok] = getComponent(obj,'x'); if ~ok, error('cannot get X'), end
    end
    function y = y(obj)
      %access Y component
      [y,ok] = getComponent(obj,'y'); if ~ok, error('cannot get Y'), end
    end
    function y = z(obj)
      %access Z component
      [y,ok] = getComponent(obj,'z'); if ~ok, error('cannot get Z'), end
    end
    function y = r(obj)
      %access R component
      [y,ok] = getComponent(obj,'r'); if ~ok, error('cannot get R'), end
    end
    function y = theta(obj)
      %access T(theta) component
      [y,ok] = getComponent(obj,'t'); if ~ok, error('cannot get THETA'), end
    end
    function y = phi(obj)
      %access P(phi) component
      [y,ok] = getComponent(obj,'p'); if ~ok, error('cannot get PHI'), end
    end
    function y = lambda(obj)
      %access L(lambda) component
      [y,ok] = getComponent(obj,'l'); if ~ok, error('cannot get LAMBDA'), end
    end
    function y = xx(obj)
      %access X component
      [y,ok] = getComponent(obj,'xx'); if ~ok, error('cannot get XX'), end
    end
    function y = xy(obj)
      %access X component
      [y,ok] = getComponent(obj,'xy'); if ~ok, error('cannot get XY'), end
    end
    function y = xz(obj)
      %access X component
      [y,ok] = getComponent(obj,'xz'); if ~ok, error('cannot get XZ'), end
    end
    function y = yx(obj)
      %access X component
      [y,ok] = getComponent(obj,'yx'); if ~ok, error('cannot get YX'), end
    end
    function y = yy(obj)
      %access X component
      [y,ok] = getComponent(obj,'yy'); if ~ok, error('cannot get YY'), end
    end
    function y = yz(obj)
      %access X component
      [y,ok] = getComponent(obj,'yz'); if ~ok, error('cannot get YZ'), end
    end
    function y = zx(obj)
      %access X component
      [y,ok] = getComponent(obj,'zx'); if ~ok, error('cannot get ZX'), end
    end
    function y = zy(obj)
      %access X component
      [y,ok] = getComponent(obj,'zy'); if ~ok, error('cannot get ZX'), end
    end
    function y = zz(obj)
      %access X component
      [y,ok] = getComponent(obj,'zz'); if ~ok, error('cannot get ZZ'), end
    end

    function obj = set.coordinateSystem(obj,value)
      if obj.tensorOrder_ < 1
        error('irf:TSeries:setcoordinateSystem:badInputs',...
          'coordinateSystem can only be set for a tensor')
      end
      if ~ischar(value)
        error('irf:TSeries:setcoordinateSystem:badInputs',...
          'expecting string input')
      end
      obj.coordinateSystem_ = value;
    end

    function obj = set.data(obj,value)
      if all(size(value) == size(obj.data_)), obj.data_ = value;
      else
        error('irf:TSeries:setdata:badInputs',...
          'size of DATA cannot be changed')
      end
    end

    function obj = set.time(obj,value)
      if ~isa(value,'GenericTimeArray')
        error('irf:TSeries:sett:badInputs',...
          'T must be of GenericTimeArray type or derived from it')
      end
      if value.length() ~= obj.length()
        error('irf:TSeries:sett:badInputs',...
          'T must have the same number of records')
      end
      obj.t_ = value;
    end

    function res = isempty(obj)
      res = isempty(obj.t_);
    end

    function Ts = clone(obj,t,data)
      % CLONE TSeries.CLONE(t,data)
      %  Clones a TSeries, but allows to add new time and data.
      if size(data,1) ~= t.length()
        error('irf:TSeries:TSeries:badInputs',...
          'T and DATA must have the same number of records')
      end
      Ts = obj;
      Ts.data_ = data;
      Ts.t_ = t;
    end

    function Ts = norm(obj)
      % NORM Make into unit vector.
      % Ts = Ts/Ts.abs;
      %
      % Example:
      %   B = B.resample(E);
      %   irf_plot({B,...                       % B xyz
      %             B.norm,...                  % B xyz norm
      %             E,...                       % E xyz
      %             E.dot(B.norm),...           % E par
      %             E.dot(B.norm)*B.norm,...    % E par xyz
      %             E - E.dot(B.norm)*B.norm})  % E perp xyz
      %

      if obj.tensorOrder~=1
        error('Only tensororder = 1 implemented');
      end
      Ts = obj/obj.abs;
    end

    function Ts = sqrt(obj)
      %SQRT Square root
      if obj.tensorOrder~=0, error('Square root requires tensorOrder = 0'); end
      obj.data_ = sqrt(obj.data); Ts = obj;
      Ts.tensorBasis_ = ''; Ts.representation{2} = [];
      if ~isempty(obj.name), Ts.name = sprintf('sqrt(%s)',obj.name); end
    end

    function Ts = cumsum(obj,option)
      %CUMSUM Cumsum
      % TS = TSeries.cumsum(option)
      %   option = 't' - cumsum over time
      %            n - integer - cumsum over data dimension n

      if obj.tensorOrder~=0, error('Cumsum requires tensorOrder = 0'); end
      if not(exist('option'))
        option = 1;
      end
      switch option
        case 't'
          obj.data_ = cumsum(obj.data,1); Ts = obj;
        otherwise % dimension is numeric value
          %if option> ndims(obj.data), error('cumsum dimension exceeds data dimension'); end
          if option>0
            obj.data_ = cumsum(obj.data,abs(option)); Ts = obj;
          elseif option<0
            obj.data_ = cumsum(obj.data(:,end:-1:1),abs(option)); Ts = obj;
          end
      end

      Ts.tensorBasis_ = ''; Ts.representation{2} = [];
      if ~isempty(obj.name), Ts.name = sprintf('cumsum(%s)',obj.name); end
    end

    function Ts = abs(obj)
      %ABS Magnitude
      if obj.tensorOrder==0
        Tmpdata = abs(obj.data);
      elseif obj.tensorOrder==1
        switch obj.basis
          case {'xy','xyz'}, Tmpdata = sqrt( sum(abs(obj.data).^2, 2) );
          case {'rtp','rlp','rp'}, Tmpdata = abs(obj.r.data);
          case 'rpz' % cylindrical
            Tmpdata = sqrt(abs(obj.r.data).^2 + abs(obj.z.data).^2);
          otherwise
            error('Unknown representation'); % should not be here
        end
        obj.representation{2} = [];
      else
        error('Not yet implemented');
      end
      obj.data_ = Tmpdata; Ts = obj;
      Ts.tensorOrder_=0; Ts.tensorBasis_ = ''; %Ts.representation{2} = []; % moved this to tensorOrder = 1
      if ~isempty(obj.name), Ts.name = sprintf('|%s|',obj.name); end
    end

    function Ts = abs2(obj)
      %ABS2 Magnitude squared
      if obj.tensorOrder~=1, error('Not yet implemented'); end
      switch obj.basis
        case {'xy','xyz'}, Tmpdata = sum(abs(obj.data).^2, 2);
        case {'rtp','rlp','rp'}, Tmpdata = abs(obj.r.data).^2;
        case 'rpz' % cylindrical
          Tmpdata = abs(obj.r.data).^2 + abs(obj.z.data).^2;
        otherwise
          error('Unknown representation'); % should not be here
      end
      obj.data_ = Tmpdata; Ts = obj;
      Ts.tensorOrder_=0; Ts.tensorBasis_ = ''; Ts.representation{2} = [];
      if ~isempty(obj.name), Ts.name = sprintf('|%s|^2',obj.name); end
      if ~isempty(obj.units), Ts.units = sprintf('(%s)^2',obj.units); end
    end

    function Ts = cross(obj,obj1)
      %CROSS cross product
      if ~isa(obj,'TSeries') ||  ~isa(obj1,'TSeries')
        error('Both imputs must be TSeries');
      elseif obj.tensorOrder~=1 || obj1.tensorOrder~=1
        error('Only scalars and vectors are supported');
      elseif obj.time~=obj1.time
        error('Input TS objects have different timelines, use resample()')
      end
      Ts = obj;
      vector_product()
      update_name_units()

      function vector_product()
        switch obj.basis
          case {'xy','rp'}
            Ts = obj.transform('xy'); d1 = Ts.data;
            d2 = obj1.transform('xy').data;
            Ts.data_ = d1(:,1).*d2(:,2) - d1(:,2).*d2(:,1);
          case {'xyz','rtp','rlp','rpz'}
            Ts = obj.transform('xyz'); d1 = Ts.data;
            d2 = obj1.transform('xyz').data;
            Ts.data_=[	d1(:,2).*d2(:,3)-d1(:,3).*d2(:,2), ...
              - (d1(:,1).*d2(:,3)-d1(:,3).*d2(:,1)), ...
              d1(:,1).*d2(:,2)-d1(:,2).*d2(:,1)];
          otherwise
            error('Unknown representation'); % should not be here
        end
      end

      function update_name_units()
        if ~isempty(obj.name) || ~isempty(obj1.name)
          if isempty(obj.name), s = 'untitled';
          else, s = obj.name;
          end
          if isempty(obj1.name), s1 = 'untitled';
          else, s1 = obj1.name;
          end
          Ts.name = sprintf('cross(%s,%s)',s,s1);
        end
        if ~isempty(obj.units) || ~isempty(obj1.units)
          if isempty(obj.units), Ts.units = obj1.units;
          elseif isempty(obj1.units), Ts.units = obj.units;
          else, Ts.units = [obj.units ' ' obj1.units];
          end
        end
      end
    end

    function l = length(obj)
      % LENGTH - number of data points in time series
      %
      % LENGTH(TS) - number of data points in TS
      if isempty(obj.t_), l = 0;
      else, l = obj.t_.length();
      end
    end

    function l = datasize(obj,option_time)
      % DATASIZE - size of data
      %
      % TS.DATASIZE - size(TS.data):
      %                    1x1 scalar -> [nTimes 1]
      %                    1x3 vector -> [nTimes 3]
      %                    3x3 tensor -> [nTimes 3 3]
      % TS.DATASIZE('dataonly') - remove leading index from time
      %                    1x1 scalar -> [1 1]
      %                    1x3 vector -> [1 3]
      %                    3x3 tensor -> [3 3]
      % TS.DATASIZE('time=1') - include leading 1
      %                    1x1 scalar -> [1 1 1]
      %                    1x3 vector -> [1 1 3]
      %                    3x3 tensor -> [1 3 3]

      if isempty(obj.t_)
        l = 0;
      elseif nargin == 1
        l = size(obj.data);
      elseif nargin>1 % option given
        switch option_time
          case 'dataonly'
            l = size(obj(1).data);
            if obj.tensorOrder == 2
              l(1) = [];
            elseif obj.tensorOrder == 1
              l(1) = 1;
            elseif length(l)>2
              l(1) = [];
            elseif length(l) == 2
              l(1) = 1;
            else
              error('Option not covered.')
            end
          case 'time=1'
            l = obj.datasize('dataonly');
            l = [1 l];
          otherwise
            error('Do not understand input.')
        end
      end
    end

    function Ts = plus(obj,obj1)
      % PLUS Addition of TS with TS/scalars.
      %
      %   TS.PLUS(TS)
      %   TS + TS
      %
      %   TS.PLUS(constant)
      %   TS + constant
      %   - Constant can be also an object of the same size
      %     as each data sample.
      %
      %   Examples:
      %     TS3 = TS1 + TS2;
      %     TS2 = TS1 + 0.1;
      %     TS2 = TS1 + [1 2 4];  % if TS1,TS2 are vector time series
      [ST,~] = dbstack; % see if plus() was called from within minus()
      if numel(ST)>1 && strcmp(ST(2).name,'TSeries.minus')
        operationStr = 'Minus'; operationSymbol = '-';
      else, operationStr = 'Plus';  operationSymbol = '+';
      end

      if isnumeric(obj1) % then obj is TSeries
        Ts = obj;
        if isempty(obj) || isempty(obj1)
          sz = size(Ts.data_); sz(1) = 0; Ts.data_ = zeros(sz);
          Ts.t_ = Ts.t_([]);
        else
          if numel(obj1) == 1
            Ts.data_ = Ts.data_ + obj1;
          else
            sizeInp = size(obj1);
            sizeObj = size(Ts.data);
            if isequal(sizeInp,sizeObj)
              Ts.data_ = Ts.data_ + obj1;
            elseif numel(sizeInp) == numel(sizeObj) ...   % same dimensions
                && sizeInp(1) == 1  ...                     % one row in obj1
                && all(sizeInp(2:end) == sizeObj(2:end))    % otherwise same number of elements
              Ts.data_ = Ts.data_ + repmat(obj1,[sizeObj(1) ones(1,numel(sizeInp)-1)]);
            else
              error([operationStr ' not defined']);
            end
          end
        end
      elseif isnumeric(obj) % then obj1 is a TSeries
        Ts = obj1;
        if isempty(obj1) || isempty(obj)
          sz = size(Ts.data_); sz(1) = 0; Ts.data_ = zeros(sz);
          Ts.t_ = Ts.t_([]);
        else
          if numel(obj) == 1
            Ts.data_ = Ts.data_ + obj;
          else
            sizeInp = size(obj);
            sizeObj = size(Ts.data);
            if isequal(sizeInp,sizeObj)
              Ts.data_ = Ts.data_ + obj;
            elseif numel(sizeInp) == numel(sizeObj) ...   % same dimensions
                && sizeInp(1) == 1  ...                     % one row in obj1
                && all(sizeInp(2:end) == sizeObj(2:end))    % otherwise same number of elements
              Ts.data_ = Ts.data_ + repmat(obj,[sizeObj(1) ones(1,numel(sizeInp)-1)]);
            else
              error([operationStr ' not defined']);
            end
          end
        end
      elseif isa(obj,'TSeries') && isa(obj1,'TSeries')
        if obj.tensorOrder~=obj1.tensorOrder
          error('Inputs must have the same tensor order');
          %elseif obj.tensorBasis_ ~= obj1.tensorBasis_
          %  error('Inputs must have the same tensor basis, use transform()');
        elseif ~strcmpi(obj.units,obj1.units)
          warning('Inputs do not have the same units')
        end
        Ts = obj;
        if isempty(obj) || isempty(obj1)
          sz = size(Ts.data_); sz(1) = 0; Ts.data_ = zeros(sz);
          Ts.t_ = Ts.t_([]);
        else
          if obj.time~=obj1.time
            error('Input TS objects have different timelines, use resample()')
          end
          Ts.data_ = obj.data + obj1.data;
        end
        update_name()
      else
        error([operationStr ' not defined']);
      end
      function update_name()
        if ~isempty(obj.name) || ~isempty(obj1.name)
          if isempty(obj.name), s = 'untitled';
          else, s = obj.name;
          end
          if ~isa(obj1,'TSeries') || isempty(obj1.name), s1 = 'untitled';
          else, s1 = obj1.name;
          end
          Ts.name = sprintf(['(%s' operationSymbol '%s)'],s,s1);
        end
        Ts.userData = [];
      end
    end

    function Ts = uplus(obj)
      % UPLUS The positive of the TS.
      %
      %   TS.UPLUS(TS)
      %   +TS
      Ts = obj;
    end

    function Ts = minus(obj,obj1)
      % MINUS Subtraction of TS with TS/scalars.
      %
      %   TS.MINUS(TS)
      %   TS - TS
      %
      %   TS.MINUS(constant)
      %   TS - constant
      %   - Constant can be also an object of the same size
      %     as each data sample.
      %
      %   Examples:
      %     TS3 = TS1 - TS2;
      %     TS2 = TS1 - 0.1;
      %     TS2 = TS1 - [1 2 4];  % if TS1,TS2 are vector time series

      Ts = plus(obj,(-1)*obj1);
    end

    function Ts = uminus(obj)
      % UMINUS The negative of the TS.
      %
      %   TS.UMINUS(TS)
      %   -TS
      Ts = -1*obj;
    end

    function Ts = dot(obj,obj1)
      %DOT  Scalar product
      %
      %  TsOut = DOT(TsA,TsB)
      %  TsOut = DOT(TsA,VECTOR), VECTOR is a constant vector
      if isa(obj1,'TSeries') && ~isa(obj,'TSeries')
        obj1Tmp = obj1; obj1 = obj; obj = obj1Tmp;
      end

      if isa(obj1,'TSeries')
        if obj.tensorOrder~=1 || obj1.tensorOrder~=1
          error('Only vectors are supported');
        elseif obj.time~=obj1.time
          error('Input TS objects have different timelines, use resample()')
        end
        Ts = obj;
        scalar_product()
      else % constant vector, not TS
        if ~isnumeric(obj1) || ~isvector(obj1)
          error('Second argument must be a vector or TS')
        end
        Ts = obj;
        switch  obj.basis
          case 'xy'
            if length(obj1) ~= 2
              error('Second argument must be a 2D vector')
            end
            Ts.data_ = obj.data(:,1).*obj1(1) + obj.data(:,2)*obj1(2);
          case 'xyz'
            if length(obj1) ~= 3
              error('Second argument must be a 3D vector')
            end
            Ts.data_ = obj.data(:,1).*obj1(1) + obj.data(:,2)*obj1(2) + ...
              + obj.data(:,3)*obj1(3);
          otherwise
            error('Only Cartesian basis supported, use transform()'); % should not be here
        end
      end
      update_name_units()
      Ts.tensorOrder_=0; Ts.tensorBasis_ = ''; Ts.representation{2} = [];

      function scalar_product()
        switch obj.basis
          case {'xy','rp'}
            d1 = obj.transform('xy').data;
            d2 = obj1.transform('xy').data;
            Ts.data_ = d1(:,1).*d2(:,1) + d1(:,2).*d2(:,2);
          case {'xyz','rtp','rlp','rpz'}
            d1 = obj.transform('xyz').data;
            d2 = obj1.transform('xyz').data;
            Ts.data_ = d1(:,1).*d2(:,1) + d1(:,2).*d2(:,2) + ...
              d1(:,3).*d2(:,3);
          otherwise
            error('Unknown representation'); % should not be here
        end
      end

      function update_name_units()
        if ~isempty(obj.name) || (isa(obj1,'TSeries') && ~isempty(obj1.name))
          if isempty(obj.name), s = 'untitled';
          else, s = obj.name;
          end
          if ~isa(obj1,'TSeries') || isempty(obj1.name), s1 = 'untitled';
          else, s1 = obj1.name;
          end
          Ts.name = sprintf('dot(%s,%s)',s,s1);
        end
        if isempty(obj.units), s = '';
        else, s = obj.units;
        end
        if ~isa(obj1,'TSeries') || isempty(obj1.units), s1 = '';
        else, s1 = obj1.units;
        end
        Ts.units = sprintf('%s %s',s,s1);
        Ts.userData = [];
      end
    end

    function Ts = times(obj,obj1)
      %TIMES element-wise multiplication '.*'.

      if ~isa(obj1,'TSeries') || ~isa(obj,'TSeries')
        if (isnumeric(obj1) || isnumeric(obj))
          Ts = obj*obj1; return
        else
          error('TSeries inputs are expected')
        end
      elseif obj.tensorOrder~=obj1.tensorOrder
        error('Inputs must have the same tensor order');
      elseif obj.tensorBasis_ ~= obj1.tensorBasis_
        error('Inputs must have the same tensor basis, use transform()');
      elseif obj.time~=obj1.time
        error('Input TS objects have different timelines, use resample()')
      end

      Ts = obj;
      Ts.data_ = obj.data.*obj1.data;
      update_name_units()

      function update_name_units()
        if ~isempty(obj.name) || ~isempty(obj1.name)
          if isempty(obj.name), s = 'untitled';
          else, s = obj.name;
          end
          if isempty(obj1.name), s1 = 'untitled';
          else, s1 = obj1.name;
          end
          Ts.name = sprintf('%s * %s',s,s1);
        end
        if isempty(obj.units), s = '';
        else, s = obj.units;
        end
        if isempty(obj1.units), s1 = '';
        else, s1 = obj1.units;
        end
        Ts.units = sprintf('%s %s',s,s1);
        Ts.userData = [];
      end
    end

    function Ts = mtimes(obj1,obj2)
      %MTIMES  Matrix and scalar multiplication '*'.

      if any([isempty(obj1),isempty(obj2)])
        Ts = TSeries([]);
        return;
      end

      % Check dimensions of input
      if isa(obj1,'TSeries')
        sizeData1 = obj1.datasize('dataonly');
        nTimes1 = obj1.length;
        data1 = obj1.data;
        tmpData1 = permute(obj1.data,[2:ndims(obj1.data) 1]); % move time to last index
        tmpData1 = reshape(tmpData1,[sizeData1 nTimes1]);
        to1 = obj1.tensorOrder;
        repres1 = obj1.representation;
        name1 = obj1.name;
        units1 = obj1.units;
      else
        sizeData1 = size(obj1);
        data1 = obj1;%size(reshape(obj1,[1 sizeObj1]));
        tmpData1 = data1;
        if isscalar(data1)
          to1 = 0;
        elseif ndims(data1)>2 % assume it's a scalar matrix
          to1 = 1;
        elseif any(sizeData1==1)
          to1 = 1;
        elseif all(size(data1)>1)
          to1 = 2;
        elseif any(size(data1)>1)
          to1 = 1;
        end
        if to1 == 0
          repres1 = {};
        elseif to1 == 1
          repres1 = {'x','y','z'};
        elseif to1 == 2
          repres1 = {{'x','y','z'};{'x','y','z'}};
        end
        name1 = '';
        units1 = '';
      end
      if isa(obj2,'TSeries')
        sizeData2 = obj2.datasize('dataonly');
        nTimes2 = obj2.length;
        data2 = obj2.data;
        tmpData2 = permute(obj2.data,[2:ndims(obj2.data) 1]); % move time to last index
        tmpData2 = reshape(tmpData2,[sizeData2 nTimes2]);
        to2 = obj2.tensorOrder;
        repres2 = obj2.representation;
        name2 = obj2.name;
        units2 = obj2.units;
      else
        sizeData2 = size(obj2);
        data2 = obj2;%reshape(obj2,[1 sizeObj2]);
        tmpData2 = data2;
        if isscalar(data2)
          to2 = 0;
        elseif ndims(data2)>2 % assume it's a scalar matrix
          to2 = 1;
        elseif any(sizeData2==1)
          to2 = 1;
        elseif all(size(data2)>1)
          to2 = 2;
        elseif any(size(data2)>1)
          to2 = 1;
        end
        if to2 == 0
          repres2 = {};
        elseif to2 == 1
          repres2 = {'x','y','z'};
        elseif to2 == 2
          repres2 = {{'x','y','z'};{'x','y','z'}};
        end
        name2 = '';
        units2 = '';
      end

      % Multiplication by a scalar, supported combinations (if dimensions allow):
      % T*s s*T s*S S*s s*V V*s S*V V*S S*T T*S
      % (upper case = TSeries, lower case = numeric)
      if to1 == 0 || to2 == 0
        if validate_dimensions('scalar') && validate_times
          [newTensorOrder,newRepresentation] = determine_new_tensororder_and_representation;
          newTime = new_time;
          newData = multiply('scalar');
          Ts = construct_new_ts;
          update_name_units;
        end
      elseif any([to1 to2])
        if validate_dimensions('noscalar') && validate_times
          [newTensorOrder,newRepresentation] = determine_new_tensororder_and_representation;
          newTime = new_time;
          newData = squeeze(multiply('tensor'));
          Ts = construct_new_ts;
          update_name_units;
        end
      else, error('Only scalar multiplication supported.')
      end

      % Nested functions
      function isOk = validate_dimensions(multiplication_type)
        switch multiplication_type
          case 'scalar'
            if numel(sizeData1) == numel(sizeData2) && all(sizeData1 == sizeData2) && all(sizeData2 ~= [1 1]) && all(sizeData1 ~= [1 1])
              warning('Interpreting both inputs as scalars.')
              to1 = 0;
              to2 = 0;
              repres1 = {};
              repres2 = {};
              isOk = 1;
            elseif (numel(sizeData1) == 2 && all( sizeData1 == [1 1])) || (numel(sizeData2) == 2 && all(sizeData2 == [1 1]))
              isOk = 1;
            else
              error('For scalar multiplication, data sizes must agree or either data set must have size [1 1]. Consider using ''.*''.')
            end
          otherwise
            if sizeData1(2) ~= sizeData2(1)
              error('For vector/matric multiplication, inner matrix dimensions must agree.')
            else
              isOk = 1;
            end
        end
      end
      function isOk = validate_times
        if isa(obj1,'TSeries') && isa(obj2,'TSeries') && any(ne(obj1.time,obj2.time))
          error('Input TS objects have different timelines, use resample().')
        else
          isOk = 1;
        end
      end
      function newTime = new_time
        if isa(obj1,'TSeries'), newTime = obj1.time;
        else, newTime = obj2.time;
        end
      end
      function newData = multiply(multiplication_type)
        switch multiplication_type
          case 'scalar' % scalar multiplication
            if isa(obj1,'TSeries') && ~isa(obj2,'TSeries')
              if isscalar(obj2)
                data2 = repmat(data2,obj1.datasize);
              else % is vector
                data2 = repmat(data2,nTimes1,1);
              end
              newData = data1.*data2;
            elseif ~isa(obj1,'TSeries') && isa(obj2,'TSeries')
              if isscalar(obj1)
                newData = data1*data2;
              else % is vector
                data1 = repmat(data1,nTimes2,1);
                newData = data1.*data2;
              end
            elseif obj1.datasize('dataonly') == obj2.datasize('dataonly')
              newData = obj1.data.*obj2.data;
            elseif all(obj1.datasize('dataonly') == [1 1]) && any(obj2.datasize('dataonly') ~= [1 1])
              data1 = repmat(data1, obj2.datasize('dataonly') );
              newData = data1.*data2;
            elseif any(obj1.datasize('dataonly') ~= [1 1]) && all(obj2.datasize('dataonly') == [1 1])
              data2 = repmat(data2, obj1.datasize('dataonly'));
              newData = data1.*data2;
            else
              error('Not supported.')
            end
          otherwise % tensor multiplication
            %error('Tensor multiplication not yet implemented.')
            %reshapedData1 = permute(data1,[]);
            %newData = zeros(sizeData1(1),sizeData2(2),newTime.length);
            newTmpData = zeros(sizeData1(1),sizeData2(2),newTime.length);
            for i1 = 1:sizeData1(1)
              for i2 = 1:sizeData2(2)
                for indSum = 1:sizeData1(2)
                  %                  newData(i1,i2,:) = newData(i1,i2,:) + data1(i1,indSum,:).*data2(indSum,i2,:);
                  newTmpData(i1,i2,:) = newTmpData(i1,i2,:) + tmpData1(i1,indSum,:).*tmpData2(indSum,i2,:);
                end
              end
            end
            sizeNewData = size(newTmpData);
            newTmpData = permute(newTmpData,[length(sizeNewData) 1:length(sizeNewData)-1]); % move time back to first index
            %tmpData2 = reshape(tmpData2,[sizeData2 nTimes2]);
            newData = newTmpData;
        end
      end
      function [newTensorOrder,newRepresentation] = determine_new_tensororder_and_representation
        switch to1
          case 0 % obj1 is a scalar
            newTensorOrder = to2;
            newRepresentation = repres2;
          case 1 % obj1 is a vector
            if to2 == 0
              newTensorOrder = to1;
              newRepresentation = repres1;
            elseif to2 == 1
              if sizeData1(1)>1 % obj1 is 3x1 vector
                newTensorOrder = 2;
                newRepresentation = {repres1;repres2};
              else % 1x3
                newTensorOrder = 0;
                newRepresentation = {};
              end
            elseif to2 == 2
              newTensorOrder = 1; % rowvector*matrix, anything else is not supported
              newRepresentation = repres1(1);
            else
              error('Can not determine new tensor order or representation.')
            end
          case 2 % obj1 is a matrix
            if to2 == 0
              newTensorOrder = to1;
              newRepresentation = repres1;
            elseif to2 == 1
              if sizeData1(1)>1 % obj2 is row vector, (1x3)*(3x3)
                newTensorOrder = 1;
                newRepresentation = repres1(2);
              else % column vector
                error('Can not determine new tensor order.')
              end
            elseif to2 == 2 % matrix*matrix
              newTensorOrder = 2;
              newRepresentation = repres1;
            else
              error('Can not determine new tensor order.')
            end
          otherwise
            error('Can not determine new tensor order.')
        end
      end
      function newTs = construct_new_ts
        switch newTensorOrder
          case 0
            newTs = TSeries(newTime,newData,'to',newTensorOrder);
          case 1
            newTs = TSeries(newTime,newData,'to',newTensorOrder,'repres',newRepresentation{1});
          case 2
            newTs = TSeries(newTime,newData,'to',newTensorOrder,'repres',newRepresentation{1},'repres',newRepresentation{2});
          otherwise
            error('Can not construct new TSeries.')
        end
      end
      function update_name_units()
        if all([isempty(name1) isempty(name2)])
          Ts.name = '';
          Ts.units = '';
        elseif any([isempty(name1) isempty(name2)])
          Ts.name = sprintf('%s%s',name1,name2);
          Ts.units = sprintf('%s%s',units1,units2);
        else
          Ts.name = sprintf('%s*%s',name1,name2);
          Ts.units = sprintf('%s*%s',units1,units2);
        end
        if isempty(Ts.units); Ts.units = ''; end
        if isempty(Ts.name); Ts.name = ''; end
      end
    end

    function Ts = power(obj,obj1)
      %MPOWER Elementwise power (.^) with scalar
      if ~isa(obj,'TSeries')
        error('First input must be a TSeries.')
      elseif ~isnumeric(obj1)
        error('Second argument must be numeric.')
      elseif ~isscalar(obj1)
        error('Second argument must be numeric scalar.')
      end

      Ts = obj;
      Ts.data_ = Ts.data.^obj1;
      update_name_units()

      function update_name_units()
        if ~isempty(obj.name)
          if isempty(obj.name), s = 'untitled';
          else, s = obj.name;
          end
          Ts.name = sprintf('(%s).^%g',s,obj1);
        end
        if isempty(obj.units), s = '';
        else, s = obj.units;
        end
        Ts.units = sprintf('(%s)*%g',s,obj1);
        Ts.userData = [];
      end
    end

    function Ts = mrdivide(obj,obj1)
      % Divide

      % Division with complete vectors or matrices is not supported.
      % Division between partial tensors of the same order is supported.
      % The new tensor order is then 0.

      if isa(obj1,'TSeries') && isnumeric(obj) && isscalar(obj) % e.g. 2/ne
        if obj1.tensorOrder > 0 && any(obj1.datasize('dataonly')~= [1 1])
          error('Can only divide by scalar TSeries or partial tensor with data size [1 1].')
        end
        Ts = obj1;
        Ts.data_ = repmat(obj,size(obj1.data))./obj1.data;
        Ts.name = sprintf('%s/(%s)','unknown',Ts.name);
        return
      elseif isa(obj,'TSeries') && isnumeric(obj1) && isscalar(obj1) % e.g. ne/2
        Ts = obj;
        Ts.data_ = obj.data/obj1;
        Ts.name = sprintf('%s/(%s)',Ts.name,'unknown');
        return
      elseif isa(obj,'TSeries') && isa(obj1,'TSeries') % both are TSeries
        if any(obj1.datasize('dataonly')~=[1 1])
          error('Denominator must be of size [1 1], e.g. either single scalar or single partial tensor.');
        elseif (obj.tensorOrder ~= obj1.tensorOrder) && obj1.tensorOrder ~= 0
          error('Tensor orders not compatible.');
        elseif (obj.tensorOrder == obj1.tensorOrder) ...
            && obj.tensorOrder>0 ...
            && any(obj.datasize('dataonly')~=[1 1]) ...
            && any(obj1.datasize('dataonly')~=[1 1])
          error('Division with partial tensors require datasizes [1 1].');
        end
        if obj.time~=obj1.time
          warning('tseries:resampling','resamplig TSeries')
          obj1 = obj1.resample(obj.time);
        end

        if (obj.tensorOrder == obj1.tensorOrder) ... % Partial tensors, new tensor order = 0
            && obj.tensorOrder>0 ...
            && all(obj.datasize('dataonly')==[1 1]) ...
            && all(obj1.datasize('dataonly')==[1 1])
          %sizeData = size(obj.data);
          newData  = obj.data./obj1.data;
          Ts = irf.ts_scalar(obj.time,newData);
        else % Retain tensor order of numerator
          Ts = obj;
          sizeData = size(obj.data);
          newData  = obj.data./repmat(obj1.data,[1 sizeData(2:end)]);
          Ts.data_ = newData;
        end
        update_name_units()

        return
      end

      function update_name_units()
        if ~isempty(obj.name) || ~isempty(obj1.name)
          if isempty(obj.name), s = 'untitled';
          else, s = obj.name;
          end
          if ~isa(obj1,'TSeries') || isempty(obj1.name), s1 = 'untitled';
          else, s1 = obj1.name;
          end
          Ts.name = sprintf('%s/(%s)',s,s1);
        end
        if isempty(obj.units), s = '';
        else, s = obj.units;
        end
        if ~isa(obj1,'TSeries') || isempty(obj1.units), s1 = '';
        else, s1 = obj1.units;
        end
        Ts.units = sprintf('%s/(%s)',s,s1);
        Ts.userData = [];
      end
    end

    function Ts = trace(obj)
      % TRACE TSeries.trace
      %   Takes the trace of tensor.
      if obj.tensorOrder ~= 2
        error('Trace only applicable to order 2 tensors')
      end

      newData = 0;
      for ii = 1:size(obj.data,2)
        newData = newData + squeeze(obj.data(:,ii,ii));
      end
      obj.data_ = newData; Ts = obj;
      Ts.tensorOrder_=0; Ts.tensorBasis_ = ''; Ts.representation{2} = [];
      if ~isempty(obj.name), Ts.name = sprintf('trace(%s)',obj.name); end
    end

    function Ts = combine(obj,obj1)
      % Combine two time series, with different times but same data type &
      % representation into a single timeseries sorted by unique timestamps
      % (EpochTT values). Note: It will keep userData, name, units,
      % coordinateSystem and siConversion from "obj".
      if ~isa(obj1,'TSeries') || ~isa(obj,'TSeries')
        error('TSeries inputs are expected')
        %elseif isa(obj1,'TSeries') && isempty(obj1)
        %  Ts =
      elseif obj.tensorOrder~=obj1.tensorOrder
        error('Inputs must have the same tensor order');
      elseif obj.tensorBasis_ ~= obj1.tensorBasis_
        error('Inputs must have the same tensor basis, use transform()');
      end

      timeClass = class(obj.time);
      obj1TimeTmp = convert_epoch(obj1.time,timeClass);
      epochTmp = [obj.time.epoch; obj1TimeTmp.epoch];
      [srt_time, srt] = sort(epochTmp);
      [epochTmp, usrt] = unique(srt_time);
      NewTime = feval(timeClass,epochTmp);

      dataTmp = [obj.data; obj1.data];
      nd = ndims(obj.data);
      switch nd
        case 2, dataTmp = dataTmp(srt(usrt), :);
        case 3, dataTmp = dataTmp(srt(usrt), :, :);
        case 4, dataTmp = dataTmp(srt(usrt), :, :, :);
        case 5, dataTmp = dataTmp(srt(usrt), :, :, :, :);
        case 6, dataTmp = dataTmp(srt(usrt), :, :, :, :, :);
        otherwise
          errStr = 'Cannot handle more than 6 dimensions.';
          irf.log('critical', errStr);
          error(errStr);
      end
      if obj.tensorOrder==0
        Ts = TSeries(NewTime, dataTmp, 'tensorOrder', obj.tensorOrder);
      elseif obj.tensorOrder==1
        Ts = TSeries(NewTime, dataTmp, 'tensorOrder', obj.tensorOrder, ...
          'tensorBasis', obj.BASIS{obj.tensorBasis_} , 'repres', obj.representation{1}); % Combined TSeries
      elseif obj.tensorOrder==2
        Ts = TSeries(NewTime, dataTmp, 'tensorOrder', obj.tensorOrder, ...
          'tensorBasis', obj.BASIS{obj.tensorBasis_} , 'repres', obj.representation{1} , 'repres', obj.representation{2}); % Combined TSeries
      else
        Ts = TSeries(NewTime, dataTmp, 'tensorOrder', obj.tensorOrder, ...
          'tensorBasis', obj.BASIS{obj.tensorBasis_}); % Combined TSeries
      end
      % Perhaps fix a better combination of metadata, for now keep "obj".
      Ts.name = obj.name;
      Ts.units = obj.units;
      Ts.siConversion = obj.siConversion;
      Ts.userData = obj.userData;
      Ts.fullDim_ = obj.fullDim_;
      if(~isempty(obj.coordinateSystem))
        Ts.coordinateSystem = obj.coordinateSystem;
      end
    end

    function Ts = resample(obj,NewTime,varargin)
      % RESAMPLE  Resample TSeries to a new timeline
      %
      % TsOut = RESAMPLE(Ts,NewTime, [ARGS])
      %     irf_resamp METHOD: 'nearest','linear','spline','pchip','cubic','v5cubic'
      % NewTime should be GeneralTimeArray (e.g. EpochTT.)
      % Resampled data type is double.
      % ARGS are given as input to irf_resamp()
      %
      % TsOut = RESAMPLE(Ts,Ts2, [ARGS])
      %
      % Resample Ts to timeline of Ts2
      %
      % See also: IRF_RESAMP

      %if isempty(obj), error('Cannot resample empty TSeries'), end
      if isempty(obj)
        irf.log('warning','Cannot resample empty TSeries.');
        Ts = TSeries([]);
        return;
      end
      if ~isa(NewTime,'TSeries') && ~isa(NewTime,'GenericTimeArray')
        error('NewTime must be of TSeries or GenericTimeArray type or derived from it')
      end
      if isa(NewTime,'TSeries'), NewTime = NewTime.time; end
      if NewTime == obj.time, Ts = obj; return, end

      if obj.tensorOrder==0
        resample_(obj);
      elseif obj.tensorOrder==1
        % For non-cartesian bases, in order to do a proper inte/extrapolarion
        % we first transform into cartesian basis, resample, and then
        % transform back to the original basis
        basis = obj.BASIS{obj.tensorBasis_};
        switch basis
          case {'xy','xyz'}, resample_(obj); return
          case {'rtp','rlp','rpz'}, resample_(obj.transform('xyz'));
          case 'rp', resample_(obj.transform('xy'));
          otherwise
            error('Unknown representation'); % should not be here
        end
        Ts = Ts.transform(basis);
      elseif obj.tensorOrder==2
        resample_(obj)
      else
        error('Not yet implemented');
      end

      function resample_(TsTmp)
        tData = double(TsTmp.time.ttns - TsTmp.time.start.ttns)/10^9;
        dataTmp = double(TsTmp.data);
        newTimeTmp = double(NewTime.ttns - TsTmp.time.start.ttns)/10^9;

        % reshape data so it can be directly inserted into irf_resamp
        origDataSize = size(dataTmp);
        dataTmpReshaped = squeeze(reshape(dataTmp,[origDataSize(1) prod(origDataSize(2:end))]));
        newDataTmpReshaped = irf_resamp([tData dataTmpReshaped], newTimeTmp, varargin{:}); % resample
        newDataReshaped = squeeze(newDataTmpReshaped(:,2:end)); % take away time column
        newData = reshape(newDataReshaped,[length(newTimeTmp) origDataSize(2:end)]); % shape back to original dimensions

        Ts = TsTmp;

        if isa(Ts,'PDist') % update ancillary and depend before time is updated
          Ts = Ts.resample_depend_ancillary(NewTime,varargin{:});
        end

        Ts.t_ = NewTime; Ts.data_ = newData;
      end
    end %RESAMPLE

    function Ts = filt(obj,varargin)
      % FILT  Filter TSeries
      %
      % TsOut = filt(Ts,fmin,fmax,[Fs],[order])
      % TsOut = Ts.filt(fmin,fmax,[Fs],[order])
      %
      % See also: IRF_FILT

      if isempty(obj), error('Cannot filter empty TSeries'), end
      if numel(varargin)<2, error('Input missing'), end

      if obj.tensorOrder==0
        filt_(obj);
      elseif obj.tensorOrder==1
        % For non-cartesian bases, in order to do a proper inter-/extrapolation
        % we first transform into cartesian basis, resample, and then
        % transform back to the original basis
        basis = obj.BASIS{obj.tensorBasis_};
        switch basis
          case {'xy','xyz'}, filt_(obj); return
          case {'rtp','rlp','rpz'}, filt_(obj.transform('xyz'));
          case 'rp', filt_(obj.transform('xy'));
          otherwise
            error('Unknown representation'); % should not be here
        end
        Ts = Ts.transform(basis);
      elseif obj.tensorOrder==2
        filt_(obj)
      else
        error('Not yet implemented');
      end

      function filt_(TsTmp)
        Ts = irf_filt(TsTmp,varargin{:}); % resample
      end
    end

    function Ts = smooth(obj,varargin)
      % SMOOTH Smooth TSeries.
      %
      %   Calls matlab function SMOOTH on each data column.
      %
      %   TsOut = smooth(Ts); % uses default smoothing span of 5
      %   TsOut = smooth(Ts,inp)
      %   TsOut = Ts.smooth;
      %
      % See also: TSERIES.FILT, SMOOTH

      if isempty(obj), error('Cannot smooth empty TSeries'), end

      data = obj.data;
      size_data = size(data);
      new_data = data*0;
      for idata = 1:prod(size_data(2:end))
        new_data(:,idata) = smooth(data(:,idata),varargin{:});
      end
      Ts = obj.clone(obj.time,new_data);
    end
    function obj = tlim(obj,tint, mode)
      %TLIM  Returns data within specified time interval
      %
      % Ts1 = TLIM(Ts,Tint, [MODE])
      %
      % Where MODE can be:
      %        'and', 0 (default)
      %        'xor', 1
      %
      % Ts1 is part of the Ts that is within interval
      % LIM.START <= X(:,1) < LIM.STOP for "AND" mode
      %
      % Ts1 is part of Ts outside the interval for "XOR" mode:
      % X(:,1) < LIM.START & X(:,1) > LIM.STOP
      %
      % See also: IRF_TLIM
      if nargin<3, mode = 0; end
      [idx,obj.t_] = obj.time.tlim(tint, mode);
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
    end

    function Ts = transform(obj, newBasis)
      % TRANSFORM  transform TSeries to a new basis
      %
      % Ts1 = TRANSFORM(Ts,newBasis)
      %
      % newBasis is string:
      %         'rlp' - Cartesian XYZ to spherical latitude
      %         'xyz' - Spherical latitude to cartesian XYZ
      %         'rpz' - Cartesian XYZ to cylindrical
      %         'xyz' - Cylidrical to cartesian XYZ
      %         'rtp' - Cartesian XYZ to spherical colatitude
      %         'xyz' - Spherical colatitude to cartesian XYZ
      %         'rlp' - Spherical colatitude to spherical latitude
      %         'rtp' - Spherical latitude to colatitude
      if ~ischar(newBasis) || min(size('adc'))>1
        error('newBasis must be a string')
      end
      oldBasis = obj.BASIS{obj.tensorBasis_};
      if strcmpi(newBasis,oldBasis), Ts = obj; return; end
      Ts = obj.changeBasis([oldBasis '>' newBasis]);
    end
  end

  methods (Access=protected)
    function [res, ok] = getComponent(obj,comp)
      res = []; ok = false;
      nd = ndims(obj.data_)-1; % first is time
      if nd>6, error('we cannot support more than 5 dimensions'), end % we cannot support more than 5 dimensions
      teno = obj.tensorOrder_;
      if length(comp)~=teno, return, end
      basis = obj.BASIS{obj.tensorBasis_};
      if ~any(basis==comp(1)), return
      elseif length(comp)==2 && ~any(char(obj.BASIS(obj.tensorBasis_))==comp(2))
        return
      end
      switch obj.tensorOrder_
        case 0, error('should no be here')
        case 1
          iDim = find((~cellfun(@isempty,obj.fullDim_)));
          iComp = find_iComp(comp);
          if isempty(iComp), return, end
          idx = cell(nd,1);
          for i = 1:nd
            idx{i} = 1:size(obj.data_,i);
            if i==iDim, idx{i} = iComp; end
          end
          % XXX: This is ugly, what is a better way of doing this?
          switch nd
            case 1, dataNew = obj.data_(:,idx{1});
            case 2, dataNew = obj.data_(:,idx{1},idx{2});
            case 3, dataNew = obj.data_(:,idx{1},idx{2},idx{3});
            case 4, dataNew = obj.data_(:,idx{1},idx{2},idx{3},idx{4});
            case 5
              dataNew = obj.data_(idx{1},idx{2},idx{3},idx{4},idx{5},idx{6});
            otherwise, error('should no be here')
          end
          args = {obj.t_,dataNew,'TensorOrder',teno,'TensorBasis',basis};
          for i=1:nd
            if i==iDim, args = [args {'repres',{comp}}]; %#ok<AGROW>
            else, args = [args {'repres',{}}]; %#ok<AGROW>
            end
          end
          res = TSeries(args{:}); res.name = sprintf('%s_%s',obj.name,comp);
          ok = true;
        case 2
          iDim = find((~cellfun(@isempty,obj.fullDim_)));
          iComp1 = find_iComp(comp(1));
          iComp2 = find_iComp(comp(2));
          iComp = [iComp1 iComp2];
          if any(isempty([iComp1 iComp2])), return, end
          idx = cell(nd,1);
          indDim = 1;
          for i = 1:nd
            idx{i} = 1:size(obj.data_,i);
            if i==iDim(indDim), idx{i} = iComp(indDim); indDim = indDim+1; end
          end
          % XXX: This is ugly, what is a better way of doing this?
          switch nd
            case 1, dataNew = obj.data_(:,idx{1});
            case 2, dataNew = obj.data_(:,idx{1},idx{2});
            case 3, dataNew = obj.data_(:,idx{1},idx{2},idx{3});
            case 4, dataNew = obj.data_(:,idx{1},idx{2},idx{3},idx{4});
            case 5
              dataNew = obj.data_(idx{1},idx{2},idx{3},idx{4},idx{5},idx{6});
            otherwise, error('should no be here')
          end
          args = {obj.t_,dataNew,'TensorOrder',teno,'TensorBasis',basis};
          indDim = 1;
          for i=1:nd
            if i==iDim(indDim),  args = [args {'repres',{comp(indDim)}}]; indDim = indDim+1;%#ok<AGROW>
            else, args = [args {'repres',{}}]; %#ok<AGROW>
            end
          end
          res = TSeries(args{:}); res.name = sprintf('%s_%s',obj.name,comp);
          ok = true;
        otherwise, error('should no be here')
      end
      function res = find_iComp(c)
        rep = obj.representation{iDim}; lRep = numel(rep);
        if rep{1}==c, res = 1;
        elseif lRep>1 && rep{2}==c, res = 2;
        elseif lRep>2 && rep{3}==c, res = 3;
        else, res = [];
        end
      end
      res.units = obj.units; res.siConversion = obj.siConversion;
    end

    function Ts = changeBasis(obj, flag)
      % Tranform from one coordinate system to another and return new
      % TimeSeries.
      % flag: = 'xyz>rlp' - Cartesian XYZ to spherical latitude
      %         'rlp>xyz' - Spherical latitude to cartesian XYZ
      %         'xyz>rpz' - Cartesian XYZ to cylindrical
      %         'rpz>xyz' - Cylidrical to cartesian XYZ
      %         'xyz>rtp' - Cartesian XYZ to spherical colatitude
      %         'rtp>xyz' - Spherical colatitude to cartesian XYZ
      %         'rtp>rlp' - Spherical colatitude to spherical latitude
      %         'rlp>rtp' - Spherical latitude to colatitude
      switch lower(flag)
        case 'xyz>rlp'
          [phi, lambda, r] = cart2sph(obj.x.data, obj.y.data, obj.z.data);
          Ts = TSeries(obj.time, [r, lambda, phi], 'vec_rlp');
        case 'rlp>xyz'
          [x, y, z] = sph2cart(obj.phi.data, obj.lambda.data, obj.r.data);
          Ts = TSeries(obj.time, [x, y, z], 'vec_xyz');
        case 'xyz>rpz'
          [phi, r, z] = cart2pol(obj.x.data, obj.y.data, obj.z.data);
          Ts = TSeries(obj.time, [r, phi, z], 'vec_rpz');
        case 'rpz>xyz'
          [x, y, z] = pol2cart(obj.phi.data, obj.r.data, obj.z.data);
          Ts = TSeries(obj.time, [x, y, z], 'vec_xyz');
        case 'xyz>rtp'
          [phi, lambda, r] = cart2sph(obj.x.data, obj.y.data, obj.z.data);
          theta = pi/2 - lambda;
          Ts = TSeries(obj.time, [r, theta, phi], 'vec_rtp');
        case 'rtp>xyz'
          lambda = pi/2 - obj.theta.data;
          [x, y, z] = sph2cart(obj.phi.data, lambda, obj.r.data);
          Ts = TSeries(obj.time, [x, y, z], 'vec_xyz');
        case 'rtp>rlp'
          lambda = pi/2 - obj.theta.data;
          Ts = TSeries(obj.time, [obj.r.data,lambda,obj.phi.data],'vec_rlp');
        case 'rlp>rtp'
          theta = pi/2 - obj.lambda.data;
          Ts = TSeries(obj.time, [obj.r.data,theta,obj.phi.data],'vec_rtp');
        case 'xy>rp'
          [phi, r] = cart2pol(obj.x.data, obj.y.data);
          Ts = TSeries(obj.time, [r, phi], 'vec_rp');
        case 'rp>xy'
          [x, y] = pol2cart(obj.phi.data, obj.r.data);
          Ts = TSeries(obj.time, [x, y], 'vec_xy');
        otherwise
          errStr='Invalid transformation'; error(errStr);
      end
    end
  end

end
