function varargout = onera_desp_lib_multi_Lstar_hmin(what,varargin)
% write/read binary files for use with multi_Lstar_hmin.mpix
% [Lm,Lstar,Bmin,Bmir,I,hmin,hmin_lon] = onera_desp_lib_multi_Lstar_hmin('read',filename)
% onera_desp_lib_multi_Lstar_hmin('write',filename,kext,options,sysaxes,dates,x1,x2,x3,alpha,maginput,R0)
% default alpha = 90, default maginput = [] = zeros(25,1), default R0=1

varargout = cell(1,nargout);
switch(lower(what))
    case {'write','save'}
        write_inputs(varargin{:});
    case {'read','load'}
        [varargout{:}] = read_outputs(varargin{:});
    otherwise
        error('Unknown action (what) = <%s>',what);
end

function write_inputs(filename,kext,options,sysaxes,dates,x1,x2,x3,alpha,maginput,R0)
if nargin < 9
    alpha = 90;
end

if nargin < 10
    maginput = [];
end

if nargin < 11
    R0 = 1;
end

if isempty(maginput)
    maginput = nan(25,1);
end

if numel(maginput)==25
    maginput = reshape(maginput,1,25); % row vector
end

if size(maginput,2) ~= 25
    error('maginput is wrong size. dim 2 should be 25 elements');
end

ntimes = numel(dates);

if numel(alpha)==1
    alpha = repmat(alpha,ntimes,1);
end

magic = 'multi_Lstar_hmin';
header = char(zeros(100,1)); % 100-byte header
header(1:length(magic)) = magic; % magic string

% control for maginput varying
% use first byte in header after magic string
% other header bytes will be used for similar things
maginput_varies = numel(maginput)>25;
if maginput_varies
    header(length(magic)+1) = char(1); % maginput varies with time
end

kext = onera_desp_lib_kext(kext);
options = onera_desp_lib_options(options);
sysaxes = onera_desp_lib_sysaxes(sysaxes);
[iyear,idoy,UT] = onera_desp_lib_matlabd2yds(dates);
fid = fopen(filename,'wb');
fwrite(fid,header,'char');
fwrite(fid,ntimes,'int32');
fwrite(fid,kext,'int32');
fwrite(fid,options,'int32');
fwrite(fid,sysaxes,'int32');
fwrite(fid,R0,'double');
fwrite(fid,maginput','double');
fwrite(fid,iyear,'int32');
fwrite(fid,idoy,'int32');
fwrite(fid,UT,'double');
fwrite(fid,x1,'double');
fwrite(fid,x2,'double');
fwrite(fid,x3,'double');
fwrite(fid,alpha,'double');
fclose(fid);

function [Lm,Lstar,Bmin,Bmir,I,hmin,hmin_lon] = read_outputs(filename)
fid = fopen(filename,'rb');
ntimes = fread(fid,1,'int32');
Lm = fread(fid,ntimes,'double');
Lstar = fread(fid,ntimes,'double');
Bmin = fread(fid,ntimes,'double');
Bmir = fread(fid,ntimes,'double');
I = fread(fid,ntimes,'double');
hmin = fread(fid,ntimes,'double');
hmin_lon = fread(fid,ntimes,'double');
fclose(fid);

% flag is -1.0000e+031

Lm(Lm < -1e+30) = nan;
Lstar(Lstar < -1e+30) = nan;
Bmin(Bmin < -1e+30) = nan;
Bmir(Bmir < -1e+30) = nan;
I(I < -1e+30) = nan;
hmin(hmin < -1e+30) = nan;
hmin_lon(hmin_lon < -1e+30) = nan;

