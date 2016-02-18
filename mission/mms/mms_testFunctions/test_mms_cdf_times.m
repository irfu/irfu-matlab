function test_mms_cdf_times( fileWithFiles )
% Input argument is file with files to be tested. (full filenames on separate lines).
% ie the output of for instance "ls /data/mms/mms[1-4]/edp/*/ql/dce2d/2015/08/mms[1-4]_edp_*_ql_dce2d_20150803*_v0.1.0.cdf"
% to test all QL data for the 3rd of August.
narginchk(1,1);
if(~exist(fileWithFiles,'file')), error('FILE NOT FOUND'); end

fileID = fopen(fileWithFiles,'r');
files = textscan(fileID,'%s');
fclose(fileID);

for ii=1:length(files{1,1})
  eachCDF(files{1,1}{ii})
end

end

function eachCDF(file)
  if(~exist(file,'file')), error('FILE NOT FOUND'); end
  [~, name, ext] = fileparts(file);
  if(~strcmpi(ext,'.cdf')), error('File Not CDF'); end
  scId = name(4);
  disp(['Verify file: ', name, ext]);
  dataObj = dataobj(file);

  if( isfield(dataObj.data,['mms',scId,'_edp_dce_epoch']) )
    time = dataObj.data.(['mms',scId,'_edp_dce_epoch']).data;
  elseif(isfield(dataObj.data,['mms',scId,'_edp_scpot_epoch']) )
    time = dataObj.data.(['mms',scId,'_edp_scpot_epoch']).data;
  else
    error('Unknown file, what is time?');
  end

  monotone_time(time)
  disp_start_stop(time)
end


% Check time is increasing for each separate cdf file
function monotone_time(time)
  if(any(diff(time)<=0))
    error('Non increasing time in cdf file');
  else
    disp('Time (of cdf file) is monoton increasing, all(diff(time)>0)');
  end
end

% Start and stop times of the cdf file
function disp_start_stop(time)
  tint = irf.tint([time(1), time(end)]);
  % Display the times
  tint %#ok<NOPRT>
end
