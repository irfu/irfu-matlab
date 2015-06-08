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
  %      'vec_xy','vec_rt'                       - 2D vectors
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
    data_
    t_ % GenericTimeArray
    fullDim_
    tensorOrder_ = 0;
    tensorBasis_ = '';
  end
  
  properties (Dependent = true)
    data
    time
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
    userData = [];
  end
  
  methods
    function obj = TSeries(t,data,varargin)
      if  nargin == 0, obj.data_ = []; obj.t_ = []; return, end
      
      if nargin<2, error('2 inputs required'), end
      if ~isa(t,'GenericTimeArray')
        error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
          'T must be of GenericTimeArray type or derived from it')
      end
      if ~isnumeric(data)
        error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
          'DATA must be numeric')
      end
      if size(data,1) ~= t.length()
        error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
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
              error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
                'DATA has more than 2 dimentions (needed for 3D vec)')
            elseif size(obj.data_,2)~=3
              error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
                'size(DATA,2) must be 3 for 3D vec')
            end
            obj.tensorOrder_ = 1; flagTensorOrderSet = true;
            [~,iB] = intersect(obj.BASIS,x(5:7));
            obj.tensorBasis_ = iB; flagTensorBasisSet = true;
            obj.representation{2} = {x(5), x(6), x(7)};
            obj.fullDim_{2} = true;
          case {'vec_xy','vec_rt'}
            if ndims(obj.data_)>2, %#ok<ISMAT>
              error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
                'DATA has more than 2 dimentions (needed for 2D vec)')
            elseif size(obj.data_,2)~=2
              error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
                'size(DATA,2) must be 2 for 2D vec')
            end
            obj.tensorOrder_ = 1; flagTensorOrderSet = true;
            [~,iB] = intersect(obj.BASIS,x(5:6));
            obj.tensorBasis_ = iB; flagTensorBasisSet = true;
            obj.representation{2} = {x(5), x(6)}; obj.fullDim_{2} = true;
          case {'to','tensororder'}
            if flagTensorOrderSet
              error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
                'tensorOrder already set')
            end
            if ~isempty(args), y = args{1}; args(1) = [];
            else
              error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
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
                  error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
                    'cannot construct tensor:all data dimensions have size<=3')
                end
                obj.tensorBasis_ = 1; % set default basis to XYZ
              end
            else
              error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
                'tensorOrder must be 0<=tensorOrder<=%d',...
                obj.MAX_TENSOR_ORDER)
            end
            
          case {'tb','basis','tensorbasis'}
            if flagTensorBasisSet
              error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
                'tensorBasis already set')
            end
            if obj.tensorOrder_==0
              error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
                'cannot set tensorBasis for a scalar (tensorOrder=0)')
            end
            if ~isempty(args), y = args{1}; args(1) = [];
            else
              error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
                'tensorBasis requires a second argument')
            end
            [~,iB] = intersect(obj.BASIS,y);
            if ~isempty(iB)
              obj.tensorBasis_ = iB; flagTensorBasisSet = true;
            else
              error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
                'tensorBais value not recognized')
            end
          case {'rep','repres','representation'}
            if isempty(obj.tensorOrder_),
              error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
                'Must specify TensorOrder first')
            end
            if iDim>ndims(data),
              error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
                'Representation already set for all DATA dimensions')
            end
            if ~isempty(args), y = args{1}; args(1) = [];
            else
              error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
                'Representation requires a second argument')
            end
            if isempty(y) && iDim<ndims(data), iDim = iDim + 1;
            else
              [ok,msg] = validate_representation(y);
              if ~isempty(ok),
                obj.fullDim_{iDim}=ok; obj.representation{iDim} = y;
                iDim = iDim + 1;
              else
                error('irf:GenericTimeArray:GenericTimeArray:badInputs',msg)
              end
            end
          case 'depend'
          otherwise
            error('irf:GenericTimeArray:GenericTimeArray:badInputs',...
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
    
    function value = get.data(obj)
      value = obj.data_;
    end
    
    function value = get.time(obj)
      value = obj.t_;
    end
   
    function value = get.tensorBasis(obj)
      value = [obj.BASIS{obj.tensorBasis_}...
        ' (' obj.BASIS_NAMES{obj.tensorBasis_} ')'];
    end
    
    function value = get.tensorOrder(obj)
      value = obj.tensorOrder_;
    end
    
    function y = tranform(obj, flag)
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
      y = [];
      switch lower(flag)
        case 'xyz>rlp'
          [p, l, r] = cart2sph(obj.x.data, obj.y.data, obj.z.data);
          y = TSeries(obj.time, [r, l, p], 'vec_rlp');
        case 'rlp>xyz'
          [x, y, z] = sph2cart(obj.p.data, obj.l.data, obj.r.data);
          y = TSeries(obj.time, [x, y, z], 'vec_xyz');
        case 'xyz>rpz'
          [p, r, z] = cart2pol(obj.x.data, obj.y.data, obj.z.data);
          y = TSeries(obj.time, [r, p, z], 'vec_rpz');
        case 'rpz>xyz'
          [x, y, z] = pol2cart(obj.p.data, obj.r.data, obj.z.data);
          y = TSeries(obj.time, [x, y, z], 'vec_xyz');
        case 'xyz>rtp'
          [p, l, r] = cart2sph(obj.x.data, obj.y.data, obj.z.data);
          t = pi/2 - l;
          y = TSeries(obj.time, [r, t, p], 'vec_rtp');
        case 'rtp>xyz'
          l = pi/2 - obj.t.data;
          [x, y, z] = sph2cart(obj.p.data, l, obj.r.data);
          y = TSeries(obj.time, [x, y, z], 'vec_xyz');
        case 'rtp>rlp'
          l = pi/2 - obj.t.data;
          y = TSeries(obj.time, [obj.r.data, l, obj.p.data], 'vec_rlp');
        case 'rlp>rtp'
          t = pi/2 - obj.l.data;
          y = TSeries(obj.time, [obj.r.data, t, obj.p.data], 'vec_rtp');
        otherwise
          errStr='Incorrect usage or conversion not yet implemented!';
          error(errStr);
      end
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
    
    function obj = set.data(obj,value)
      if all(size(value) == size(obj.data_)), obj.data_ = value;
      else
        error('irf:GenericTimeArray:setdata:badInputs',...
          'size of DATA cannot be changed')
      end
    end
    
    function obj = set.time(obj,value)
      if ~isa(t,'GenericTimeArray')
        error('irf:GenericTimeArray:sett:badInputs',...
          'T must be of GenericTimeArray type or derived from it')
      end
      if value.length() ~= obj.length()
        error('irf:GenericTimeArray:sett:badInputs',...
          'T must have the same number of records')
      end
      obj.t_ = value;
    end
    
    function res = isempty(obj)
      res = isempty(obj.t_);
    end
    
    function l = length(obj)
      if isempty(obj.t_), l = 0;
      else l = obj.t_.length();
      end
    end
    function obj = plus(obj,inp)
      if isnumeric(inp) 
        if numel(inp) == 1
          obj.data_ = obj.data_ + inp;
        else
          sizeInp = size(inp);
          sizeObj = size(obj.data);
          if numel(sizeInp) == numel(sizeObj) && ...
              sizeInp(1) == 1 && ...
              all(sizeInp(2:end) == sizeObj(2:end))
            obj.data_ = obj.data_ + repmat(inp,[sizeObj(1) ones(1,numel(sizeInp)-1)]);
          end
        end
      else
        error('Plus not defined');
      end
    end
    function obj = minus(obj,inp)
      if isnumeric(inp) && ...
          ((numel(inp) == 1) || (size(inp,2) == size(obj.data_,2)))
        obj.data_ = obj.data_ - inp;
      else
        error('Plus not defined');
      end
    end
    function obj = mtimes(obj,inp)
      if isnumeric(inp) && numel(inp) == 1,
        obj.data_ = obj.data_ * inp;
      else
        error('mtimes not defined');
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
          error('irf:GenericTimeArray:subsref',...
            'Not a supported subscripted reference')
      end
    end

    function obj = tlim(obj,tint)
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
  end
  
  methods (Access=protected)
    function res = getComponent(obj,comp)
      res = [];
      nd = ndims(obj.data_);
      if nd>6, error('we cannot support more than 5 dimensions'), end % we cannot support more than 5 dimensions
      teno = obj.tensorOrder_;
      if length(comp)~=teno, return, end
      basis = obj.BASIS{teno};
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
          res = TSeries(args{:});
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
  end
  
end

