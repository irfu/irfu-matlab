classdef Data
  %DATA Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    data
    nElements
    sizes
    format
    
  end
  
  methods
    
    function da = Data(inp)
      if isempty(inp)
        da.data = [];
      end
    end
    
    function n = dimension(da)
      % Dimension of data
      % 0 - scalar
      % 1 - vector
      % 2 - thensor
      if isempty(da.sizes), n = 0;
      else n = length(da.sizes);
      end
    end
  end %methods
  
end

