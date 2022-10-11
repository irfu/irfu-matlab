function export(ebsp,tint,cl_id,freqRange)
%EXPORT  Export data to CEF
%
% export( ebsp, tint, cl_id, freqRange)
%
% Export EBSP into CEF
%
% export( rotMatrix, tint, cl_id)
%
% export FAC rotation matrix into CEF
%
% Will create a gzipped CEF file in the current directory

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------
%
% This software was developed as part of the MAARBLE (Monitoring,
% Analyzing and Assessing Radiation Belt Energization and Loss)
% collaborative research project which has received funding from the
% European Community's Seventh Framework Programme (FP7-SPACE-2011-1)
% under grant agreement n. 284520.

% This must be changed when we do any major changes to our processing software
DATASET_VERSION = '0';

% We do not track versions here, CAA will do this for us
DATA_VERSION = '00';

% Define formats for output
FORMAT_EXP = '%9.2e,'; % Amplitudes
FORMAT_ANG = '%6.0f,'; % Angles - integer values
FORMAT_DEG = '%6.1f,'; % Degree of ... -1..1 or 0..1
FORMAT_EXP_ANG_ANG = [FORMAT_EXP FORMAT_ANG FORMAT_ANG];
FORMAT_ROTMATR = '%11.6f,' ; % Rotation matrix, < 1
FORMAT_R = '%10.2f,' ; % distance in km < 20 RE
% Replace NaN with FILLVAL (specified in the CEF header)
FILLVAL            = -999;
FILLVAL_EXP        = -1.00E+31;

dataToExport = []; nData = [];
if ~isstruct(ebsp)
  error('expecting sturcture output of irf_ebsp or irf_convert_fac')
end

if any( ebsp.t<tint(1) | ebsp.t>tint(end) )
  irf.log('critical','data outside TINT')
  error('data outside TINT')
end

flagHasE = true;
if isnumeric(cl_id) % Cluster
  flagCluster = true;
else
  flagCluster = false;
  if cl_id(1)=='g' % GOES
    flagHasE = false;
  end
end

if isfield(ebsp,'rotMatrix')
  export_rotMatrix
elseif isfield(ebsp,'flagFac')
  export_ebsp
else
  error('expecting sturcture output of irf_ebsp or irf_convert_fac')
end

%% Write out the data
% Time array goes first
out_CharArray = epoch2iso(ebsp.t,1); % Short ISO format: 2007-01-03T16:00:00.000Z
out_CharArray(:,end+1)=',';

% Write out data by columns and then combine into a common char matrix
for i=1:length(dataToExport)
  tmp_CharArray = sprintf(dataToExport{i}{1},dataToExport{i}{2}');
  tmp_CharArray = reshape(tmp_CharArray,length(tmp_CharArray)/nData,nData)';
  out_CharArray = [out_CharArray tmp_CharArray]; clear tmp_CharArray %#ok<AGROW>
end

% Add END_OF_RECORD markers in the end of each line
out_CharArray(:,end)='$';
out_CharArray(:,end+1)=sprintf('\n'); %#ok<SPRINTFN>
out_CharArray = out_CharArray';
out_CharArray=out_CharArray(:)';

%% Prepare the file header
if flagCluster, filePrefix = sprintf('C%d', cl_id);
else
  filePrefix = 'CC';
  datasetID = [upper(cl_id) '_' datasetID];
end
fileName = [filePrefix '_CP_AUX_MAARBLE_' datasetID...
  '__' irf_fname(tint,5) '_V' DATA_VERSION];

header = [...
  sprintf('!-------------------- CEF ASCII FILE --------------------|\n')...
  sprintf('! created on %s\n', char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss")))...
  sprintf('!--------------------------------------------------------|\n')...
  sprintf('FILE_NAME = "%s.cef"\n',fileName)...
  sprintf('FILE_FORMAT_VERSION = "CEF-2.0"\n')...
  sprintf('END_OF_RECORD_MARKER = "$"\n')...
  sprintf('include = "%s_CH_AUX_MAARBLE_%s.ceh"\n', filePrefix, datasetID)...
  sprintf(pmeta('FILE_TYPE','cef'))...
  sprintf(pmeta('DATASET_VERSION',DATASET_VERSION))...
  sprintf(pmeta('LOGICAL_FILE_ID',fileName))...
  sprintf(pmeta('VERSION_NUMBER',DATA_VERSION))...
  sprintf('START_META     =   FILE_TIME_SPAN\n')...
  sprintf('   VALUE_TYPE  =   ISO_TIME_RANGE\n')...
  sprintf('   ENTRY       =   %s/%s\n', ...
  epoch2iso(tint(1),1),epoch2iso(tint(2),1))...
  sprintf('END_META       =   FILE_TIME_SPAN\n')...
  sprintf('START_META     =   GENERATION_DATE\n')...
  sprintf('   VALUE_TYPE  =   ISO_TIME\n')...
  sprintf('   ENTRY       =   %s\n', char(datetime("now","Format","uuuu-MM-dd'T'HH:mm:ss.SSS'Z'")))...
  sprintf('END_META       =   GENERATION_DATE\n')...
  sprintf('!\n')...
  sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')...
  sprintf('!                       Data                          !\n')...
  sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')...
  sprintf('DATA_UNTIL = "END_OF_DATA"\n')];

%% Write the file
if 0
  % Write to plain CEF
  f = fopen([fileName '.cef'],'w'); %#ok<UNRCH>
  fwrite(f,header);
  fwrite(f,out_CharArray);
  fwrite(f,sprintf('END_OF_DATA\n'));
  fclose(f);
  return
else
  % Write directly GZIPed CEF file
  fileOutStream = java.io.FileOutputStream(java.io.File([fileName '.cef.gz']));
  gzipOutStream = java.util.zip.GZIPOutputStream( fileOutStream );
  gzipOutStream.write(java.lang.String(header).getBytes());
  gzipOutStream.write(java.lang.String(out_CharArray).getBytes());
  gzipOutStream.write(java.lang.String(sprintf('END_OF_DATA\n')).getBytes());
  gzipOutStream.close;
end

  function export_rotMatrix
    nData = length(ebsp.t);
    rm2d = zeros(nData,9);
    for iRow=1:3
      for iCol=1:3
        rm2d(:,(iRow-1)*3+iCol) = ebsp.rotMatrix(:,iRow,iCol);
      end
    end
    rm2d(isnan(rm2d)) = FILLVAL;
    dataToExport = { {FORMAT_ROTMATR, rm2d} };
    datasetID = 'ULF_FACMATR';
    if ~flagCluster
      r = ebsp.r; r(isnan(r)) = FILLVAL_EXP;
      dataToExport = [dataToExport {{FORMAT_R, r}}];
    end
  end % export_rotMatrix()

  function export_ebsp
    %% Check the input
    if ~ebsp.flagFac, error('EBSP must be in FAC'), end
    switch lower(freqRange)
      case 'pc12'
        %DT2 = 0.5; % time resolution
        datasetID = 'ULF_PC12';
        if ischar(cl_id) && cl_id(1)=='g', numberOfFreq = 12; % GOES
        else, numberOfFreq = 21; % Cluster and THEMIS
        end
      case 'pc35'
        %DT2 = 30; % time resolution
        datasetID = 'ULF_PC35';
        numberOfFreq = 21;
      otherwise
        error('freqRange must be ''pc12'' or ''pc35''')
    end
    nFreq = length(ebsp.f); nData = length(ebsp.t);
    if nFreq~=numberOfFreq
      error('number of frequencies in ebsp.f must be %d (not %d!)',...
        numberOfFreq,nFreq)
    end
    
    %% Prepare data array
    % B0
    if isempty(ebsp.fullB), magB = ebsp.B0; else, magB = ebsp.fullB; end
    magB = irf_abs(magB); magB = magB(:,[1 5]); magB = irf_resamp(magB,ebsp.t);
    magB = magB(:,2);
    
    ebsp.k_tp(isnan(ebsp.k_tp)) = FILLVAL;
    ebsp.ellipticity(isnan(ebsp.ellipticity)) = FILLVAL;
    ebsp.planarity(isnan(ebsp.planarity)) = FILLVAL;
    ebsp.dop(isnan(ebsp.dop)) = FILLVAL;
    ebsp.dop2d(isnan(ebsp.dop2d)) = FILLVAL;
    ebsp.pf_rtp(isnan(ebsp.pf_rtp)) = FILLVAL;
    ebsp.planarity(isnan(ebsp.planarity)) = FILLVAL;
    magB(isnan(magB)) = FILLVAL_EXP;
    ebsp.bb_xxyyzzss(isnan(ebsp.bb_xxyyzzss)) = FILLVAL_EXP;
    if flagHasE
      if isempty(ebsp.ee_ss)
        ebsp.ee_ss = ones(nData,nFreq)*FILLVAL_EXP;
      else
        ebsp.ee_ss(isnan(ebsp.ee_ss)) = FILLVAL_EXP;
      end
    end
    
    % fliplr to make frequencies ascending
    ebsp.ellipticity = fliplr(ebsp.ellipticity);
    ebsp.planarity = fliplr(ebsp.planarity);
    ebsp.dop = fliplr(ebsp.dop);
    ebsp.dop2d = fliplr(ebsp.dop2d);
    if flagHasE
      ebsp.ee_ss = fliplr(ebsp.ee_ss);
    end
    
    % Reformat matrices/vectors and fliplr to make frequencies ascending
    BB_2D = zeros(nData,nFreq*3);
    for comp=1:3
      BB_2D(:,((1:nFreq)-1)*3+comp) = fliplr(ebsp.bb_xxyyzzss(:,:,comp));
    end
    K = zeros(nData,nFreq*2);
    for comp=1:2
      K(:,((1:nFreq)-1)*2+comp) = fliplr(ebsp.k_tp(:,:,comp));
    end
    if flagHasE
      PV = ones(nData,nFreq*3)*FILLVAL;
      if ~isempty(ebsp.pf_rtp)
        for comp=1:3
          PV(:,((1:nFreq)-1)*3+comp) = fliplr(ebsp.pf_rtp(:,:,comp));
        end
      end
    end
    
    % NOTE: This list must be consistent with the CEF header file
    if ~flagHasE % Only magnetic field is available
      dataToExport = {...
        {FORMAT_EXP, BB_2D},...              % BB_xxyyzz_fac
        {FORMAT_ANG, K},...                  % KSVD_fac
        {FORMAT_DEG, ebsp.ellipticity},...   % ELLSVD
        {FORMAT_DEG, ebsp.planarity},...     % PLANSVD
        {FORMAT_DEG, ebsp.dop},...           % DOP
        {FORMAT_DEG, ebsp.dop2d},...         % POLSVD
        {FORMAT_EXP, magB}                   % BMAG
        };
    else % Both E & B
      dataToExport = {...
        {FORMAT_EXP, BB_2D},...              % BB_xxyyzz_fac
        {FORMAT_ANG, K},...                  % KSVD_fac
        {FORMAT_DEG, ebsp.ellipticity},...   % ELLSVD
        {FORMAT_DEG, ebsp.planarity},...     % PLANSVD
        {FORMAT_DEG, ebsp.dop},...           % DOP
        {FORMAT_DEG, ebsp.dop2d},...         % POLSVD
        {FORMAT_EXP_ANG_ANG, PV},...         % PV
        {FORMAT_EXP, ebsp.ee_ss},...         % ESUM
        {FORMAT_EXP, magB}                   % BMAG
        };
      
      % For Pc3-5 we also add E spectrum in FAC
      if strcmpi(freqRange,'pc35')
        ebsp.ee_xxyyzzss(isnan(ebsp.ee_xxyyzzss)) = FILLVAL_EXP;
        % Reformat E matrix and fliplr to make frequencies ascending
        EE_2D = ones(nData,nFreq*3)*FILLVAL;
        if ~isempty(ebsp.ee_xxyyzzss)
          for comp=1:3
            EE_2D(:,((1:nFreq)-1)*3+comp) = fliplr(ebsp.ee_xxyyzzss(:,:,comp));
          end
        end
        dataToExport = [dataToExport {{FORMAT_EXP, EE_2D}}]; % EE_xxyyzz_fac
      end
    end
  end % export_ebsp()

end %export

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obuf = pmeta(metaID,metaValue)
% Print META
if isnumeric(metaValue), q = ''; metaValue = num2str(metaValue); else, q = '"'; end
obuf = [...
  'START_META     =   ' metaID '\n'...
  '   ENTRY       =   ' q metaValue q '\n'...
  'END_META       =   ' metaID '\n'];
end