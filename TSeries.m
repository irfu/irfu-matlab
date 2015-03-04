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
  %      Fpr more infor on representation see Cluster Metadata Dictionary
  %           http://caa.estec.esa.int/caa/doc_format_meta.xml
  %
  %  Example:
  %
  %  epoch = EpochUnix(iso2epoch('2002-03-04T09:30:00Z')+(0:3));
  %  data4x2 = [1 2; 3 4; 5 6; 7 8];
  %  data4x3x3 = rand(4,3,3);
  %
  %  % TensorOrder=1, 2 components of imcomplete vector (2D)
  %  TsVec = TSeries(epoch,data4x2,'TensorOrder',1,'TensorBasis','xy',...
  %          'repres',{'x','y'})
  %  
  %  % TensorOrder=2, 3x3 components of tensor (3D, default basis=xyz)
  %  TsTen = TSeries(epoch,data4x3x3,'TensorOrder',2,...
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
    rep_
    tensorOrder_ = 0;
    tensorBasis_ = '';
  end
  
  properties (Dependent = true)
    data
    t
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
      obj.rep_ = cell(ndims(data),1);
      obj.representation = cell(ndims(data),1); iDim = 1;
      if ~isempty(obj.t_)
        obj.representation{iDim} = obj.t_; iDim = iDim + 1;
      end
      
      args = varargin;
      flagTensorOrderSet = false; flagTensorBasisSet = false;
      while ~isempty(args)
        x = args{1}; args(1) = [];
        switch lower(x)
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
                obj.rep_{iDim}=ok; obj.representation{iDim} = y;
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
        ok = ''; msg = '';
        sDim = size(obj.data,iDim); tb = obj.BASIS{obj.tensorBasis_};
        if sDim>length(tb),
          msg = sprintf(...
            'Dimension %d size %d>%d (Basis=%s) cannot have Representation. Use Depend instead',...
            iDim,sDim,length(tb),tb);
          return
        end
        if ~iscell(x), 
          msg = 'Representation requires a cell input'; return
        end
        if sDim~=length(x), 
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
        elseif length(s) == c, ok = 1; % complete dim
        elseif length(s) < c, ok = -1; % incomplete dim
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
    
    function value = get.t(obj)
      value = obj.t_;
    end
    
    function value = get.tensorBasis(obj)
      value = [obj.BASIS{obj.tensorBasis_}...
        ' (' obj.BASIS_NAMES{obj.tensorBasis_} ')'];
    end
    
    function value = get.tensorOrder(obj)
      switch obj.tensorOrder_
        case 0, value = '0 (Scalar)';
        case 1, value = '1 (Vector)';
        case 2, value = '2 (Tensor)';
      end
    end
    
    function obj = set.data(obj,value)
      if all(size(value) == size(obj.data_)), obj.data_ = value;
      else
        error('irf:GenericTimeArray:setdata:badInputs',...
          'size of DATA cannot be changed')
      end
    end
    
    function obj = set.t(obj,value)
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
  end
  
end

