classdef TimeSeriesData
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties(Dependent)
    t
    data
  end
  properties(Access=private)
    privateTimeArray
    privateData
  end
  
  methods
    function tsd = TimeSeriesData(timeArray,data)
      if ~isa(timeArray,'irf.TimeArray')
        error('irf:TimeSeriesData:setData:badInputs',...
            'timeArray must be of class irf.TimeArray')
      end
      if size(data,1) ~= length(timeArray)
        error('irf:TimeSeriesData:setData:badInputs',...
          'timeArray and data must have the same number of records')
      end
      tsd.privateTimeArray = timeArray;
      tsd.privateData = data;
    end
    
    function data = get.data(tsd)
      data = tsd.privateData;
    end
    
    function timeArray = get.t(tsd)
      timeArray = tsd.privateTimeArray;
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
    
    function tsd = set.t(tsd,timeArray)
      if length(tsd) == length(timeArray)
        tsd.privateTimeArray = timeArray;
      else
        error('irf:TimeSeriesData:setT:badInputs',...
            'Number of points in data(%d) does not match the time vector(%d)',...
          size(data,1),length(tsd))
      end
    end
    
    function r = length(tsd)
      r = length(tsd.privateTimeArray);
    end
  end % methods
  
end

