function OUT=c_efw_burst_st(filename)
%C_EFW_BURST_ST Check the starting time of the burst according to ISDAT
%
% OUT=c_efw_burst_st(filename)
%

% By Anders Tjulin 13/8-2002

if nargin<1
  error "You need at least a filename"
end

% Get the data

fid = fopen(filename,'r','s');
data = fread(fid,10,'uint32');

time1=data(5);time2=data(6);
OUT=fromepoch(time1);
OUT(end)=OUT(end)+time2/1e9;

fclose(fid);
