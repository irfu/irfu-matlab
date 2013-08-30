%
% Returns default configuration parameters for the Cef toolbox
%
%
function params=cefget
%
% Corresponding values in globals.c 
% CheckSyntax == CAAStrict
% CheckCAAfileNameSyntax == CAAfilename
% TreatAsBinartData == binaryData
% IgnoreWarnings 
% CheckDataFieldStructure == simpleData
% CheckValueFields == allData
% ExportAsGzipped == gzip 
% 
params = struct('CheckSyntax','off',...
                'CheckCAAfileNameSyntax','on',...
                'TreatAsBinaryData', 'on', ...
                'IgnoreWarnings', {{'ceflib:Warning:CefSyntaxError', 'ceflib:Info','ceflib:Info:GeneralError'}}, ...
                'CheckDataFieldStructure', 'off',...
                'CheckValueFields', 'on', ... 
                'ExportAsGzipped', 'off');


end 