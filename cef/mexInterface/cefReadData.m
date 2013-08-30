%
% Syntax: 
% data=cefReadData(file,[1 inf],[],[],header_efw, header_glob);
% First argument must be a string path to a CEF file
% Second argument must be a scalar or vector
% Third argument must be a string or cell containing product names, like
% 'B_mag__C4_CP_FGM_SPIN' or empty []
% Fourth argument must be a string or cell array containing the file
% timespan or empty []
% Fifth argument must be a string path to a CEF header directory or empty
% More arguments to different directories for CEF header...
%
% Ex: Read time tag from one record only 
%  data=cefReadData(file2,[100 100],'time_tags',[],header_fgm, header_glob);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% This is a dummy file you should have gotten a mex file 
disp('Mex file not found please check your path');
