classdef TimeVar
  %TimeVar Generic time dependent variable
  %   Time varibale has 2 fields: T - time [GenericTimeArray] and DATA 
  %
  %   a = TimeVar(t,data)
  %
  %  Example:
  %  a = TimeVar(EpochUnix(iso2epoch('2002-03-04T09:30:00Z')+(0:3)),(1:4)')
  
  
% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

  properties (Access=protected)
    data_
    t_ % GenericTimeArray
  end
  
  properties (Dependent = true)
    data
    t
  end
  
  methods
    function obj = TimeVar(t,data,varargin)
      if  nargin > 0
        if nargin<2, error('2 inputs required'), end
        if ~isa(t,'GenericTimeArray')
          error('T must be of GenericTimeArray type or derived from it')
        end
        if size(data,1) ~= t.length()
          error('T and DATA must have the same number of records')
        end
        obj.data_ = data; obj.t_ = t;
      else obj.data_ = []; obj.t_ = [];
      end
    end
    
    function value = get.data(obj)
      value = obj.data_;
    end
    
    function value = get.t(obj)
      value = obj.t_;
    end
    
    function obj = set.data(obj,value)
      if all(size(value) == size(obj.data_)), obj.data_ = value;
      else
        error('size of DATA cannot be changed')
      end
    end
    
    function obj = set.t(obj,value)
      if ~isa(t,'GenericTimeArray')
        error('T must be of GenericTimeArray type or derived from it')
      end
      if value.length() ~= obj.length()
        error('T must have the same number of records')
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

