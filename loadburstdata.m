function OUT=loadburstdata(filename,cols)

%LOADBURSTDATA Load Cluster EFW internal burst data from file
%
% DATA=LOADBURSTDATA('FILENAME',COLS) loads Cluster EFW internal
% burst data from file with name 'FILENAME'. COLS is the number of
% columns in the output data matrix.
%  
% By Anders Tjulin, last update 7/11-2002

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',
...
mfilename,'c_efw_burt_load')

if nargin<2
  cols=4;
end

% Get the data

fid = fopen(filename,'r','l');
data = fread(fid,'int16');

% Remove header

data = data(87:end);

% Reshape the data matrix

rest = length(data)-floor(length(data)/cols)*cols;
data = [data;zeros(cols-rest,1)];
OUT = reshape(data,cols,length(data)/cols)';

fclose(fid);
