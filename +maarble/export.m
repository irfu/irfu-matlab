function export( cl_id, freqRange, timeVector,...
    BB_xxyyzz_fac,EESum_xxyy_isr2,EE_xxyyzz_fac,...
    Poynting_rThetaPhi_fac,k_thphSVD_fac,polSVD_fac,ellipticity)
%EXPORT  Export data to CEF
%
% export( cl_id, freqRange, timeVector,BB_xxyyzz_fac,...
%          EESum_xxyy_isr2,EE_xxyyzz_fac,Poynting_rThetaPhi_fac,...
%          k_thphSVD_fac,polSVD_fac,ellipticity)
%
% Will create a gzipped CEF file in the current directory
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

% This must be changed when we do any major changes to our processing software
DATASET_VERSION = '0';

% We do not track versions here, CAA will do this for us
DATA_VERSION = '00';

%% Prepare the data array
nFreq = size(BB_xxyyzz_fac,2);
nData = size(BB_xxyyzz_fac,1);

% Replace NaN with FILLVAL (specified in the CEF header)
FILLVAL            = -999;
FILLVAL_EXP        = -1.00E+31;

BB_xxyyzz_fac(isnan(BB_xxyyzz_fac)) = FILLVAL_EXP;
k_thphSVD_fac(isnan(k_thphSVD_fac)) = FILLVAL;
ellipticity(isnan(ellipticity)) = FILLVAL;
polSVD_fac(isnan(polSVD_fac)) = FILLVAL;
Poynting_rThetaPhi_fac(isnan(Poynting_rThetaPhi_fac(:,:,1))) = FILLVAL_EXP;
Poynting_rThetaPhi_fac(isnan(Poynting_rThetaPhi_fac)) = FILLVAL;
EESum_xxyy_isr2(isnan(EESum_xxyy_isr2)) = FILLVAL_EXP;
EE_xxyyzz_fac(isnan(EE_xxyyzz_fac)) = FILLVAL_EXP;

% Define formats for output
formatExp = '%9.2e,'; % Amplitudes
formatAng = '%4.0f,'; % Angles - integer values
formatDeg = '%6.1f,'; % Degree of ... -1..1 or 0..1

% Reformat B matrix
BB_2D = zeros(nData,nFreq*3);
for comp=1:3,BB_2D(:,((1:nFreq)-1)*3+comp) = BB_xxyyzz_fac(:,:,comp); end

dataToExport = {...
    {formatExp, BB_2D},...                         % BB_xxyyzz_fac
    {formatAng, k_thphSVD_fac(:,:,1)},...          % THSVD_fac
    {formatAng, k_thphSVD_fac(:,:,2)},...          % PHSVD_fac
    {formatDeg, ellipticity},...                   % ELLSVD
    {formatDeg, polSVD_fac},...                    % POLSVD
    {formatExp, Poynting_rThetaPhi_fac(:,:,1)},... % AMPV
    {formatAng, Poynting_rThetaPhi_fac(:,:,2)},... % THPV
    {formatAng, Poynting_rThetaPhi_fac(:,:,3)},... % PHPV
    {formatExp, EESum_xxyy_isr2}                   % ESUM
    };

% For Pc3-5 we also add E spectrum in FAC
if strcmpi(freqRange,'pc35')
    % Reformat E matrix
    EE_2D = zeros(nData,nFreq*3);
    for comp=1:3,EE_2D(:,((1:nFreq)-1)*3+comp) = EE_xxyyzz_fac(:,:,comp); end
    dataToExport = [dataToExport {formatExp, EE_2D}]; % EE_xxyyzz_fac
end
    
% Time array goes first
out_CharArray = epoch2iso(timeVector,1); % Short ISO format: 2007-01-03T16:00:00.000Z
out_CharArray(:,end+1)=',';

for i=1:length(dataToExport)
    tmp_CharArray = sprintf(dataToExport{i}{1},dataToExport{i}{2});
    tmp_CharArray = reshape(tmp_CharArray,length(tmp_CharArray)/nData,nData)';
    out_CharArray = [out_CharArray tmp_CharArray]; clear tmp_CharArray %#ok<AGROW>
end

% Add END_OF_RECORD markers in the end of each line
out_CharArray(:,end)='$';
out_CharArray(:,end+1)=sprintf('\n');
out_CharArray = out_CharArray';
out_CharArray=out_CharArray(:)';

%% Prepare the file header
switch lower(freqRange)
    case 'pc12'
        DT2 = 0.5;
        datasetID = 'MAARBLE_ULF_PC12';
    case 'pc35';
        DT2 = 30;
        datasetID = 'MAARBLE_ULF_PC35';
    otherwise
        error('freqRange must be ''pc12'' or ''pc35''')
end
fileName = sprintf('C%d_CP_AUX_%s', cl_id, datasetID);
header = [...
    sprintf('!-------------------- CEF ASCII FILE --------------------|\n')...
    sprintf('! created on %s\n', datestr(now))...
    sprintf('!--------------------------------------------------------|\n')...
    sprintf('FILE_NAME = "%s.cef"\n',fileName)...
    sprintf('FILE_FORMAT_VERSION = "CEF-2.0"\n')...
    sprintf('END_OF_RECORD_MARKER = "$"\n')...
    sprintf('include = "C%d_CH_AUX_%s.ceh"\n', cl_id, datasetID)...
    sprintf(pmeta('FILE_TYPE','cef'))...
    sprintf(pmeta('DATASET_VERSION',DATASET_VERSION))...
    sprintf(pmeta('LOGICAL_FILE_ID',fileName))...
    sprintf(pmeta('VERSION_NUMBER',DATA_VERSION))...
    sprintf('START_META     =   FILE_TIME_SPAN\n')...
	sprintf('   VALUE_TYPE  =   ISO_TIME_RANGE\n')...
    sprintf('   ENTRY       =   %s/%s\n', ...
        epoch2iso(timeVector(1)-DT2,1),epoch2iso(timeVector(end)+DT2,1))...
    sprintf('END_META       =   FILE_TIME_SPAN\n')...
    sprintf('START_META     =   GENERATION_DATE\n')...
    sprintf('   VALUE_TYPE  =   ISO_TIME\n')...
    sprintf('   ENTRY       =   %s\n', epoch2iso(date2epoch(now()),1))...
    sprintf('END_META       =   GENERATION_DATE\n')...
    sprintf('!\n')...
    sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')...
    sprintf('!                       Data                          !\n')...
    sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')...
    sprintf('DATA_UNTIL = "END_OF_DATA"\n')];

%% Write the file
if 0
    f = fopen([fileName '.cef'],'w');
    fwrite(f,header);
    fwrite(f,out_CharArray);
    fwrite(f,sprintf('END_OF_DATA\n'));
    fclose(f);
    return
else
    % Write directly GZIPed file
    fileOutStream = java.io.FileOutputStream(java.io.File([fileName '.cef.gz']));
    gzipOutStream = java.util.zip.GZIPOutputStream( fileOutStream );
    gzipOutStream.write(java.lang.String(header).getBytes());
    gzipOutStream.write(java.lang.String(out_CharArray).getBytes());
    gzipOutStream.write(java.lang.String(sprintf('END_OF_DATA\n')).getBytes());
    gzipOutStream.close;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obuf = pmeta(metaID,metaValue)
% Print META
if isnumeric(metaValue), q = ''; metaValue = num2str(metaValue); else q = '"'; end
obuf = [...
    'START_META     =   ' metaID '\n'...
    '   ENTRY       =   ' q metaValue q '\n'...
    'END_META       =   ' metaID '\n'];
end