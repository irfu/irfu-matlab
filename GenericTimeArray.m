classdef (Abstract) GenericTimeArray
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  methods 
    toUtc(obj)
  end
  
  methods (Static)
    function [ output_args ] = validate_iso_time_str( input_args )
      %verify_iso_time_str Summary of this function goes here
      %   Detailed explanation goes here
      
      output_args = true;
    end
  end
end
