function filePath=get_file(fileUrlLink,varargin)
% IRF.GET_FILE if the file is not already on the path download
% the link to temporary location and return path
%
% [filePath] = get_file(fileUrlLink);
% [filePath] = get_file(fileUrlLink,'appName','key');
%		check if files is not there based on 'appName' and 'key'.
%		save the location of the file so that it can be reused at
%		the next matlab session. Uses DATASTORE.
%
%	Example:
%		linkUrlFile = 'https://www.space.irfu.se/cluster/matlab/indexCaaMeta.mat';
%		fileIndexCaa = irf.get_file(linkUrlFile,'caa','indexFile');
%

ii = strfind(fileUrlLink,'/');
fileName = fileUrlLink(ii(end)+1:end);

if nargin==1
  useDatastore = false;
elseif nargin==3
  useDatastore = true;
  appName = varargin{1};
  keyName = varargin{2};
end
if useDatastore
  if ~exist(fileName,'file') % file is not on matlab path
    temp = datastore(appName,keyName); % returned saved file path
    if isempty(temp) % no path saved from earlier sessions
      [filePath,status] = get_index_file();
      if ~status, return; end
    else
      if ~exist(temp,'file') % check if there is no file at the specified location
        datastore('@delete',appName,keyName); % remove the keyname because points to nonexisting file
        [filePath,status] = get_index_file();
        if ~status, return; end
      else
        filePath = temp; % returned file path saved at earlier sessions
      end
    end
  end
else

end

if nargout == 0
  clear filePath;
end

  function [filePath,status] = get_index_file
    disp(['You do not have file ' fileName '!']);
    disp('The file is located at:');
    disp(['<a href="' fileUrlLink '">' fileUrlLink '</a>']);
    disp('You can download it to somewhere on your matlab path,');
    disp('or I can download it for you to some temporary path.');
    disp('WARNING!!! Files from temporary path may be removed by system after some time.');
    reply = irf_ask('Shall I download to temporary path? y/n [%]>','y','y');
    if strcmpi(reply,'y')
      fileDir = tempname;
      mkdir(fileDir);
      filePath = [fileDir filesep fileName];
      if verLessThan('matlab', '8.4')
        [f,status]=urlwrite(fileUrlLink,filePath); %#ok<URLWR> websave introduced in R2014b
      else
        try
          f = websave(filePath, fileUrlLink);
          status = true;
        catch
          status = false;
        end
      end
      if status
        disp('Success!');
        if useDatastore
          datastore(appName,keyName,f);
        end
        status = true;
        disp(['File downloaded to: ' filePath]);
        disp('I will remember it, but it will be removed by system after some time.');
        disp('If you want to keep it, move the file to somewhere on your matlab path.');
      else
        disp('Did not succeed. Maybe problems with internet connections.');
        disp('Try again later or complain...');
        status = false;
      end
    else
      disp('OK! Download to matlab path and rerun the program.')
      filePath = [];
      status = false;
    end
  end
end
