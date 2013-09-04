classdef Data
  %DATA Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    data = [];
    name = '';
  end
  properties(Dependent)
    dimension
  end
  properties(Access=private)
      privateDimension = [];
   end
  
  methods
    
    function da = Data(inp,dim,name)
      if isempty(inp)
        da.data = [];
      else
        da.data = inp;
      end
      if nargin < 2 || isempty(dim)
        if isempty(da.data), da.privateDimension = [];
        elseif isscalar(da.data), da.privateDimension = 0;
        elseif isrow(da.data), da.privateDimension = 1;
        else  da.privateDimension = ndims(da.data);
        end
      else
        da.dimension = dim;
      end
      if nargin < 3 || isempty(name)
        symbols = ['a':'z' 'A':'Z' '0':'9'];
        MAX_ST_LENGTH = 50;
        stLength = randi(MAX_ST_LENGTH);
        nums = randi(numel(symbols),[1 stLength]);
        da.name = symbols (nums);
      else
        da.name = name;
      end
    end % costructor
    
    function da = set.dimension(da,dim)
      % Dimension of data
      % 0 - scalar
      % 1 - vector
      % 2 - thensor
      
      if dim==0 && isscalar(da.data) ||...
          dim == 1 && isrow(da.data) || ...
          dim == ndims(da.data)
        da.privateDimension = dim;
      else
        error('irf:Data:setDimension:badInputs',...
          'Specified dimension does not match the data')
      end
    end
    
    function res = get.dimension(da)
      res = da.privateDimension;
    end
    
    function da = set.name(da,name)
      if ~ischar(name)
        error('irf:Data:setName:badInputs',...
          'Name must be a string')
      end
      da.name = name;
    end
    
  end %methods
  
end

