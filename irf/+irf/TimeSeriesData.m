classdef TimeSeriesData < irf.TimeArray
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties(Dependent)
    data
  end
  properties(Access=private)
    privateData
  end
  
  methods
    function tsd = TimeSeriesData(timeArray,data)
      tsd = tsd@irf.TimeArray(timeArray);
      tsd.data = data;  
    end
    
    function tsd = set.data(tsd,data)
      if length(tsd) == size(data,1)
        tsd.privateData = data;
      else
        error('irf:TimeSeriesData:setData:badInputs',...
            'Number of points in data(%d) does not match the time vector(%d)',...
          size(data,1),length(tsd))
      end
    end
  end % methods
  
end

