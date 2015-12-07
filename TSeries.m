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
  %      For more infor on representation see Cluster Metadata Dictionary
  %           http://caa.estec.esa.int/caa/doc_format_meta.xml
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
      if  nargin == 0, obj.data_ = []; obj.t_ = []; return, end
      
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
      obj.fullDim_ = cell(ndims(data),1);
      obj.representation = cell(ndims(data),1); iDim = 1;
      if ~isempty(obj.t_)
        obj.representation{iDim} = obj.t_; iDim = iDim + 1;
      end
      
      args = varargin;
      flagTensorOrderSet = false; flagTensorBasisSet = false;
      while ~isempty(args)
        x = args{1}; args(1) = [];
        switch lower(x)
          case {'vec_xyz','vec_rtp','vec_rlp','vec_rpz'}
            if ndims(obj.data_)>2, %#ok<ISMAT>
              error('irf:TSeries:TSeries:badInputs',...
                'DATA has more than 2 dimentions (needed for 3D vec)')
            elseif size(obj.data_,2)~=3
              error('irf:TSeries:TSeries:badInputs',...
                'size(DATA,2) must be 3 for 3D vec')
            end
            obj.tensorOrder_ = 1; flagTensorOrderSet = true;
            [~,iB] = intersect(obj.BASIS,x(5:7));
            obj.tensorBasis_ = iB; flagTensorBasisSet = true;
            obj.representation{2} = {x(5), x(6), x(7)};
            obj.fullDim_{2} = true;
          case {'vec_xy','vec_rp'}
            if ndims(obj.data_)>2, %#ok<ISMAT>
              error('irf:TSeries:TSeries:badInputs',...
                'DATA has more than 2 dimentions (needed for 2D vec)')
            elseif size(obj.data_,2)~=2
              error('irf:TSeries:TSeries:badInputs',...
                'size(DATA,2) must be 2 for 2D vec')
            end
            obj.tensorOrder_ = 1; flagTensorOrderSet = true;
            [~,iB] = intersect(obj.BASIS,x(5:6));
            obj.tensorBasis_ = iB; flagTensorBasisSet = true;
            obj.representation{2} = {x(5), x(6)}; obj.fullDim_{2} = true;
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
              if y>0,
                % Check if we have any data dimension with 1..3 elements
                found = false;
                for i=2:ndims(obj.data_)
                  if size(obj.data_,i)<=3, found = true; break, end
                end
                if ~found,
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
                'tensorBais value not recognized')
            end
          case {'rep','repres','representation'}
            if isempty(obj.tensorOrder_),
              error('irf:TSeries:TSeries:badInputs',...
                'Must specify TensorOrder first')
            end
            if iDim>ndims(data),
              error('irf:TSeries:TSeries:badInputs',...
                'Representation already set for all DATA dimensions')
            end
            if ~isempty(args), y = args{1}; args(1) = [];
            else
              error('irf:TSeries:TSeries:badInputs',...
                'Representation requires a second argument')
            end
            if isempty(y) && iDim<ndims(data), iDim = iDim + 1;
            else
              [ok,msg] = validate_representation(y);
              if ~isempty(ok),
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
        sDim = size(obj.data,iDim); tb = obj.BASIS{obj.tensorBasis_};
        if sDim>length(tb),
          msg = sprintf(...
            'Dimension %d size %d>%d (Basis=%s) cannot have Representation. Use Depend instead',...
            iDim,sDim,length(tb),tb);
          return
        end
        if ~iscell(x)
          msg = 'Representation requires a cell input'; return
        end
        if sDim~=length(x)
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
        else error('should not be here')
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
          if numel(idx) > 1,
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
      y = getComponent(obj,'x'); if isempty(y), error('cannot get X'), end
    end
    function y = y(obj)
      %access Y component
      y = getComponent(obj,'y'); if isempty(y), error('cannot get Y'), end
    end
    function y = z(obj)
      %access Z component
      y = getComponent(obj,'z'); if isempty(y), error('cannot get Z'), end
    end
    function y = r(obj)
      %access R component
      y = getComponent(obj,'r'); if isempty(y), error('cannot get R'), end
    end
    function y = theta(obj)
      %access T(theta) component
      y = getComponent(obj,'t'); if isempty(y), error('cannot get THETA'), end
    end
    function y = phi(obj)
      %access P(phi) component
      y = getComponent(obj,'p'); if isempty(y), error('cannot get PHI'), end
    end
    function y = lambda(obj)
      %access L(lambda) component
      y = getComponent(obj,'l'); if isempty(y), error('cannot get LAMBDA'), end
    end
    
    function obj = set.coordinateSystem(obj,value)
      if obj.tensorOrder_ < 1, 
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
    
    function Ts = abs(obj)
      %ABS Magnitude
      if obj.tensorOrder~=1, error('Not yet implemented'); end
      switch obj.basis
        case {'xy','xyz'}, Tmpdata = sqrt( sum(abs(obj.data).^2, 2) );
        case {'rtp','rlp','rp'}, Tmpdata = abs(obj.r.data);
        case 'rpz' % cylindrical
          Tmpdata = sqrt(abs(obj.r.data).^2 + abs(obj.z.data).^2);
        otherwise
          error('Unknown representation'); % should not be here
      end
      obj.data_ = Tmpdata; Ts = obj;
      Ts.tensorOrder_=0; Ts.tensorBasis_ = ''; Ts.representation{2} = [];
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
          else s = obj.name;
          end
          if isempty(obj1.name), s1 = 'untitled';
          else s1 = obj1.name;
          end
          Ts.name = sprintf('cross(%s,%s)',s,s1);
        end
        if ~isempty(obj.units) || ~isempty(obj1.units)
          if isempty(obj.units), Ts.units = obj1.units;
          elseif isempty(obj1.units), Ts.units = obj.units;
          else Ts.units = [obj.units ' ' obj1.units];
          end
        end
      end
    end
    
    function l = length(obj)
			% LENGTH - number of data points in time series
			%
			% LENGTH(TS) - number of data points in TS
      if isempty(obj.t_), l = 0;
      else l = obj.t_.length();
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
      [ST,I] = dbstack; % see if plus() was called from within minus()
      if numel(ST)>1 && strcmp(ST(2).name,'TSeries.minus')
            operationStr = 'Minus'; operationSymbol = '-';
      else, operationStr = 'Plus';  operationSymbol = '+';
      end
          
      if isnumeric(obj) && isa(obj1,'TSeries')
        Ts = plus(obj1,obj);        
      end
      if isnumeric(obj1)
        Ts = obj;
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
      elseif isa(obj,'TSeries') && isa(obj1,'TSeries')
        if obj.time~=obj1.time
          error('Input TS objects have different timelines, use resample()')
        elseif ~strcmpi(obj.units,obj1.units)
          error('Input TS objects have different units')
        end
        Ts = obj;
        Ts.data_ = obj.data + obj1.data;
        update_name()
      else        
        error([operationStr ' not defined']);
      end
      function update_name()
        if ~isempty(obj.name) || ~isempty(obj1.name)
          if isempty(obj.name), s = 'untitled';
          else s = obj.name;
          end
          if ~isa(obj1,'TSeries') || isempty(obj1.name), s1 = 'untitled';
          else s1 = obj1.name;
          end
          Ts.name = sprintf(['(%s' operationSymbol '%s)'],s,s1);                              
        end
        Ts.userData = [];
      end
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
        if ~isempty(obj.name) || ~isempty(obj1.name)
          if isempty(obj.name), s = 'untitled';
          else s = obj.name;
          end
          if ~isa(obj1,'TSeries') || isempty(obj1.name), s1 = 'untitled';
          else s1 = obj1.name;
          end
          Ts.name = sprintf('dot(%s,%s)',s,s1);
        end
        if isempty(obj.units), s = 'unknown';
        else s = obj.units;
        end
        if ~isa(obj1,'TSeries') || isempty(obj1.units), s1 = 'unknown';
          else s1 = obj1.units;
        end
        Ts.units = sprintf('%s %s',s,s1);
        Ts.userData = [];
      end
    end
    
    function Ts = times(obj,obj1)
      %TIMES  by element multiplication
      
      if ~isa(obj1,'TSeries') || ~isa(obj,'TSeries')
        error('TSeries inputs are expected')
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
          else s = obj.name;
          end
          if isempty(obj1.name), s1 = 'untitled';
          else s1 = obj1.name;
          end
          Ts.name = sprintf('%s * %s',s,s1);
        end
        if isempty(obj.units), s = 'unknown';
        else s = obj.units;
        end
        if isempty(obj1.units), s1 = 'unknown';
          else s1 = obj1.units;
        end
        Ts.units = sprintf('%s %s',s,s1);
        Ts.userData = [];
      end
    end
    
    function Ts = mtimes(obj,obj1)
      %MTIMES  multiply TS by scalar
      if isa(obj1,'TSeries') && ~isa(obj,'TSeries')
        obj1Tmp = obj1; obj1 = obj; obj = obj1Tmp;
      end
        
      % Multiplication by a scalar
      if isa(obj1,'TSeries')
        if obj.tensorOrder>0 && obj1.tensorOrder>0
          error('One of TS inputs must be a scalar');
        elseif obj.tensorOrder<obj1.tensorOrder
          obj1Tmp = obj1; obj1 = obj; obj = obj1Tmp;
        elseif obj.time~=obj1.time
          error('Input TS objects have different timelines, use resample()')
        end
          
        Ts = obj;
        switch obj.tensorOrder
          case 0
            sz = size(obj.data,2); sz1 = size(obj1.data,2); 
            if sz>1 && sz1>1 && sz~=sz1
              error('Data dimensions do not match')
            elseif sz==sz1
              Ts.data_ = obj.data.*obj1.data;
            else
              if sz==1, d = obj.data; d1 = obj1.data;
              else d1 = obj.data; d = obj1.data;
              end
              Ts.data_ = d.*repmat(d1,1,size(d,2));
            end
          otherwise
            if size(obj1.data,2)>1,
              error('TS1.data has more than one column')
            end
            Ts.data_ = obj.data.*repmat(obj1.data,1,size(obj.data,2));
        end
        update_name_units()
        return
      end
      if ~isnumeric(obj1)
        error('second argument must be numeric or TSeries')
      elseif ~isscalar(obj1)
        error('only scalars are supported')
      end
      Ts = obj; Ts.data = Ts.data*obj1;
      
      function update_name_units()
        if ~isempty(obj.name) || ~isempty(obj1.name)
          if isempty(obj.name), s = 'untitled';
          else s = obj.name;
          end
          if ~isa(obj1,'TSeries') || isempty(obj1.name), s1 = 'untitled';
          else s1 = obj1.name;
          end
          Ts.name = sprintf('%s * %s',s,s1);
        end
        if isempty(obj.units), s = 'unknown';
        else s = obj.units;
        end
        if ~isa(obj1,'TSeries') || isempty(obj1.units), s1 = 'unknown';
          else s1 = obj1.units;
        end
        Ts.units = sprintf('%s %s',s,s1);
        Ts.userData = [];
      end
    end
    
    function Ts = mrdivide(obj,obj1)
      if isa(obj1,'TSeries') && ~isa(obj,'TSeries')
        obj1Tmp = obj1; obj1 = obj; obj = obj1Tmp;
      end
        
      % Division
      if isa(obj1,'TSeries')
        if obj.tensorOrder>1 || obj1.tensorOrder>1
          error('Only scalars and vectors are supported'); 
        end
        if obj1.tensorOrder>obj.tensorOrder
          Ts = mtimes(obj1,obj); return;
        end
        if obj.time~=obj1.time
          warning('tseries:resampling','resamplig TSeries')
          obj1 = obj1.resample(obj.time);
        end
        Ts = obj;
        switch obj.tensorOrder
          case 0, Ts.data_ = obj.data./obj1.data; 
          case 1
            switch obj1.tensorOrder
              case 0
                Ts.data_ = obj.data./repmat(obj1.data,1,size(obj.data,2));
              otherwise
                error('Not supported')
            end
          otherwise
            error('Not supported')
        end
        update_name_units()
        return
      end
      if ~isnumeric(obj1)
        error('second argument must be numeric or TSeries')
      elseif ~isscalar(obj1)
        error('only scalars are supported')
      end
      Ts = obj; Ts.data = Ts.data/obj1;
      
      function update_name_units()
        if ~isempty(obj.name) || ~isempty(obj1.name)
          if isempty(obj.name), s = 'untitled';
          else s = obj.name;
          end
          if ~isa(obj1,'TSeries') || isempty(obj1.name), s1 = 'untitled';
          else s1 = obj1.name;
          end
          Ts.name = sprintf('%s/(%s)',s,s1);
        end
        if isempty(obj.units), s = 'unknown';
        else s = obj.units;
        end
        if ~isa(obj1,'TSeries') || isempty(obj1.units), s1 = 'unknown';
          else s1 = obj1.units;
        end
        Ts.units = sprintf('%s/(%s)',s,s1);
        Ts.userData = [];
      end
    end
    
    function Ts = combine(obj,obj1)
      % Combine two time series, with different times but same data type &
      % representation into a single timeseries sorted by unique timestamps
      % (EpochTT values). Note: It will keep userData, name, units,
      % coordinateSystem and siConversion from "obj".
      if ~isa(obj1,'TSeries') || ~isa(obj,'TSeries')
        error('TSeries inputs are expected')
      elseif obj.tensorOrder~=obj1.tensorOrder
        error('Inputs must have the same tensor order');
      elseif obj.tensorBasis_ ~= obj1.tensorBasis_
        error('Inputs must have the same tensor basis, use transform()');
      end
      Tmptime = [obj.time.epoch; obj1.time.epoch];
      Tmpdata = [obj.data; obj1.data];
      [srt_time, srt] = sort(Tmptime);
      [Tmptime, usrt] = unique(srt_time);
      Tmpdata = Tmpdata(srt(usrt),:);
      if(obj.tensorOrder==0)
        Ts = TSeries(EpochTT(Tmptime), Tmpdata, ...
          'tensorOrder', obj.tensorOrder);
      else
        Ts = TSeries(EpochTT(Tmptime), Tmpdata, ...
          'tensorOrder', obj.tensorOrder, ...
          'tensorBasis', obj.BASIS{obj.tensorBasis_}); % Combined TSeries
      end
      % Perhaps fix a better combination of metadata, for now keep "obj".
      Ts.name = obj.name;
      Ts.units = obj.units;
      Ts.siConversion = obj.siConversion;
      Ts.userData = obj.userData;
      if(~isempty(obj.coordinateSystem))
        Ts.coordinateSystem = obj.coordinateSystem;
      end
    end

    function Ts = resample(obj,NewTime,varargin)
      % RESAMPLE  Resample TSeries to a new timeline
      %
      % TsOut = RESAMPLE(Ts,NewTime, [ARGS])
      % 
			% NewTime should be GeneralTimeArray (e.g. EpochTT.)
			% Resampled data type is double. 
      % ARGS are given as input to irf_resamp()
      %
      % TsOut = RESAMPLE(Ts,Ts2, [ARGS])
      %
      % Resample Ts to timeline of Ts2
      %
      % See also: IRF_RESAMP
      if ~isa(NewTime,'TSeries') && ~isa(NewTime,'GenericTimeArray')
        error('NewTime must be of TSeries or GenericTimeArray type or derived from it')
      end
      if isa(NewTime,'TSeries'), NewTime = NewTime.time; end
      if NewTime == obj.time, Ts = obj; return, end 
            
      if obj.tensorOrder==0, 
          resample_(obj);
      elseif obj.tensorOrder==1,       
          % For non-cartesian bases, in order to do a proper inte/extrapolarion
          % we first transform into cartesian basis, resample, and then
          % transform back to teh original basis
          basis = obj.BASIS{obj.tensorBasis_};
          switch basis
            case {'xy','xyz'}, resample_(obj); return 
            case {'rtp','rlp','rpz'}, resample_(obj.transform('xyz'));
            case 'rp', resample_(obj.transform('xy'));
            otherwise
              error('Unknown representation'); % should not be here
          end
          Ts = Ts.transform(basis);
      else
          error('Not yet implemented'); 
      end
           
      function resample_(TsTmp)
        tData = double(TsTmp.time.ttns - TsTmp.time.start.ttns)/10^9;
        dataTmp = double(TsTmp.data);
        newTimeTmp = double(NewTime.ttns - TsTmp.time.start.ttns)/10^9;
        newData = irf_resamp([tData dataTmp], newTimeTmp, varargin{:});
        Ts = TsTmp; Ts.t_ = NewTime; Ts.data_ = newData(:,2:end);
      end
    end %RESAMPLE
    
    function obj = tlim(obj,tint)
      %TLIM  Returns data within specified time interval
      %
      % Ts1 = TLIM(Ts,Tint)
      %
      % See also: IRF_TLIM
      [idx,obj.t_] = obj.time.tlim(tint);
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
    function res = getComponent(obj,comp)
      res = [];
      nd = ndims(obj.data_);
      if nd>6, error('we cannot support more than 5 dimensions'), end % we cannot support more than 5 dimensions
      teno = obj.tensorOrder_;
      if length(comp)~=teno, return, end
      basis = obj.BASIS{obj.tensorBasis_};
      if ~any(basis==comp(1)), return
      elseif length(comp)==2 && ~any(obj.BASIS(obj.tensorBasis_)==comp(2))
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
            case 2, dataNew = obj.data_(idx{1},idx{2});
            case 3, dataNew = obj.data_(idx{1},idx{2},idx{3});
            case 4, dataNew = obj.data_(idx{1},idx{2},idx{3},idx{4});
            case 5, dataNew = obj.data_(idx{1},idx{2},idx{3},idx{4},idx{5});
            case 6
              dataNew = obj.data_(idx{1},idx{2},idx{3},idx{4},idx{5},idx{6});
            otherwise, error('should no be here')
          end
          args = {obj.t_,dataNew,'TensorOrder',teno,'TensorBasis',basis};
          for i=2:nd
            if i==iDim, args = [args {'repres',{comp}}]; %#ok<AGROW>
            else args = [args {'repres',{}}]; %#ok<AGROW>
            end
          end
          res = TSeries(args{:}); res.name = sprintf('%s_%s',obj.name,comp);
        case 2
          error('not implemented')
        otherwise, error('should no be here')
      end
      function res = find_iComp(c)
        rep = obj.representation{iDim}; lRep = length(rep);
        if rep{1}==c, res = 1;
        elseif lRep>1 && rep{2}==c, res = 2;
        elseif lRep>2 && rep{3}==c, res = 3;
        else res = [];
        end
      end
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

