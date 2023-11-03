function [variableMatName] = variable_mat_name(variableCdfName)
%VARIABLE_MAT_NAME get matlab variable name from cdf variable name
%   Matlab variable name cannot have '-','.', cannot start with number etc.

variableMatName=variableCdfName;

% Replace minuses with underscores
variableMatName(strfind(variableMatName,'-')) = '_';

% Remove training dots
while (variableMatName(end) == '.')
  variableMatName(end) = [];
end

% Take care of '...'
d3 = strfind(variableMatName,'...');
if d3
  variableMatName( d3 + (1:2) ) = [];
end

% Replace dots with underscores
variableMatName(strfind(variableMatName,'.')) = '_';

% Add "x" if the varible name starts with a number
if ~isletter(variableMatName(1))
  variableMatName=['x' variableMatName];
end

% Take care of names longer than 63 symbols (Matlab limit)
if length(variableMatName)>63
  %Try typical changes (based on FEEPS variable names)
  variableMatName = strrep(variableMatName,'electron','elec');
  variableMatName = strrep(variableMatName,'energy','ener');
  variableMatName = strrep(variableMatName,'sensorid','sensid');
  variableMatName = strrep(variableMatName,'contamination','contam');
  %continue with truncate if above fails
  if length(variableMatName)>63
    variableMatName = variableMatName(1:63);
  end
end

% if changes made write out log
if ~strcmp(variableMatName,variableCdfName)
  irf.log('notice',['orig var : ' variableCdfName]);
  irf.log('notice',[' new var : ' variableMatName]);
end

end

