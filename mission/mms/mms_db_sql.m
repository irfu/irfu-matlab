classdef mms_db_sql < handle
	%MMS_DB_SQL handle MMS database of files and variables in SQLITE
	%   IN DEVELOPMENT, DO NOT USE!!!!
	
	properties (Access=protected)
		connection = [];
		statement
		databaseDirectory
	end
	
	properties
		databaseFile
	end
	
	properties (Dependent = true)
		isConnected
	end
	
	methods
		function obj = mms_db_sql(fileName)
          javaPath=[irf('path'), filesep, 'contrib', filesep, 'java', filesep];
          listDir=dir([javaPath, 'sqlite-jdbc*.jar']);
          if(isempty(listDir)), error('Missing sqlite-jdbc.jar'); end
          javaaddpath([javaPath, listDir(end).name]);
			if nargin == 1,
				[dirPath,file,ext] = fileparts(fileName);
				if isempty(dirPath) || strcmp(dirPath,'.')
					dirPath = pwd;
				end
				obj.databaseDirectory = dirPath;
				obj.databaseFile = [dirPath filesep file ext];
				obj.open_database;
			end
        end
		
        function open_database(obj)
          % Create tables if they do not exist
          sql = [...
            'PRAGMA foreign_keys = ON;', ...
            'CREATE TABLE IF NOT EXISTS "FileList" (', ...
              '"idFile" INTEGER, ', ...
              '"directory",', ...
              '"dataset",', ...
              '"date",', ...
              '"version",', ...
              '"fileNameFullPath" TEXT UNIQUE, ', ...
              'PRIMARY KEY(idFile) );', ...
            'CREATE TABLE IF NOT EXISTS "FileListToImport" (', ...
              '"directory",', ...
              '"dataset",', ...
              '"date",', ...
              '"version",', ...
              '"fileNameFullPath" TEXT UNIQUE);', ...
            'CREATE TABLE IF NOT EXISTS "Datasets" (', ...
              '"idDataset" INTEGER NOT NULL,', ...
              '"dataset" TEXT,', ...
              '"varNames" TEXT,', ...
              'PRIMARY KEY(idDataset) );', ...
            'CREATE TABLE IF NOT EXISTS "VarNames" (', ...
              '"idVar" INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,', ...
              '"varName" TEXT,', ...
              '"idDataset" INTEGER, ', ...
              'FOREIGN KEY(idDataset) REFERENCES Datasets(idDataset) ', ...
              'ON UPDATE CASCADE ON DELETE CASCADE );', ...
            'CREATE TABLE IF NOT EXISTS "VarIndex" (', ...
              '"idFile" INTEGER,', ...
              '"idDataset" INTEGER,', ...
              '"startTT" INTEGER,', ...
              '"endTT" INTEGER, ', ...
              'FOREIGN KEY(idFile) REFERENCES FileList(idFile) ', ...
              'ON UPDATE CASCADE ON DELETE CASCADE, ', ...
              'FOREIGN KEY(idDataset) REFERENCES Datasets(idDataset) ',...
              'ON UPDATE CASCADE ON DELETE CASCADE);', ...
              ];
            obj.sqlUpdate(sql);
        end
		
		function connect(obj)
			% CONNECT open connection to database
			if ~obj.isConnected
				d = org.sqlite.JDBC;
				p = java.util.Properties();
				obj.connection = d.createConnection(...
					['jdbc:sqlite:' obj.databaseFile],p);
				obj.statement = obj.connection.createStatement();
			end
		end
		
		function close(obj)
			% CLOSE close connection to database
			obj.connection.close % close connection
			obj.connection = [];
		end
		
		function value = get.isConnected(obj)
			if isempty(obj.connection)
				value = false;
			else
				value = true;
			end
		end
		
		function import_a_file_from_list(obj)
          % import one new file from FileListToImport
          sql = ['SELECT * FROM FileListToImport ',...
            'WHERE fileNameFullPath NOT IN ', ...
            '(SELECT fileNameFullPath FROM FileList) ', ...
            'ORDER BY rowid DESC limit 1'];
          rs = obj.sqlQuery(sql);
          while rs.next
            fileToImport = char(rs.getString('fileNameFullPath'));
            status = obj.import_file(fileToImport);
            if(status==0)
              irf.log('warning',['Failed to import file: ',fileToImport]);
            end
          end
          if(exist('fileToImport','var'))
            % remove file from FileListToImport
            sql = ['DELETE from FileListToImport ',...
              'WHERE fileNameFullPath="' fileToImport '"'];
            obj.sqlUpdate(sql);
          else
            % No previously files not previously known were found, simply
            % clear FileListToImport
            irf.log('notice','No new files (not already in database) was found.');
            sql = 'DELETE FROM FileListToImport';
            obj.sqlUpdate(sql);
          end
        end

        function clear_deleted_files(obj)
          % Clear up database by removing all entries of files which have
          % been deleted.
          irf.log('warning', 'Re-running add_all_files_to_import_list to check which files are present on system.');
          % Clear any old entries in FileListToImport.
          obj.sqlUpdate('DELETE FROM FileListToImport');
          obj.add_all_files_to_import_list;
          % Verify at least some files are found. (ie avoid deleting all entries).
          sql = 'SELECT * FROM FileListToImport ORDER BY rowid LIMIT 1';
          rs = obj.sqlQuery(sql);
          if(~rs.next), irf.log('warning','Nothing found. Aborting.'); return; end
          sql = ['DELETE FROM FileList WHERE idFile IN ',...
            '(SELECT idFile FROM FileList WHERE fileNameFullPath NOT IN ',...
            '(SELECT fileNameFullPath FROM FileListToImport))'];
          obj.sqlUpdate(sql);
        end

        function import_files_from_list(obj)
          % import all files from FileListToImport
          someFilesDidNotImport = false;
% This could possibly speed up adding each result to "fileList" by
% pre-allocating it (when it is really huge). But this should be combined
% using some cleaver UNION so it is only one SQL query first of all..
%           sql = ['SELECT count(*) as records FROM FileListToImport ', ...
%             'WHERE fileNameFullPath NOT IN ', ...
%             '(SELECT fileNameFullPath FROM FileList) ',...
%             'ORDER BY rowid DESC'];
%           rs = obj.sqlQuery(sql);
%           if(rs.next)
%             count = rs.getInt(1);
%             fileList = cell(count,1);
%           end

          sql = ['SELECT * FROM FileListToImport ', ...
            'WHERE fileNameFullPath NOT IN ', ...
            '(SELECT fileNameFullPath FROM FileList) ',...
            'ORDER BY rowid DESC'];
          rs = obj.sqlQuery(sql);
          ii = 1;
          while rs.next
            fileToImport = char(rs.getString('fileNameFullPath'));
            fileList{ii,1} = fileToImport; %#ok<AGROW>
            ii = ii + 1;
          end
          if(exist('fileList','var'))
            % We got at least one new file to be imported.
            for ii=1:length(fileList)
              irf.log('warning',['Importing ' fileList{ii,1}]);
              status = obj.import_file(fileList{ii,1});
              if(status==0)
                irf.log('warning',['******** Did not succeed to import:',fileList{ii,1}]);
                someFilesDidNotImport = true;
              end
            end
          else
            irf.log('warning','No new files (not already in database) was found.');
            someFilesDidNotImport = true;
          end
          if someFilesDidNotImport
            irf.log('critical','Some files did not import!');
            % irf.log('warning','Partially cleaning up FileListToImport.');
            % sql = ['DELETE FROM FileListToImport ',...
            %  'WHERE fileNameFullPath IN ',...
            %  '(SELECT fileNameFullPath FROM FileList)'];
            % obj.sqlUpdate(sql);
          else
            irf.log('warning','All files imported from FileListToImport.');
          end
          % Clean up FileListToImport
          irf.log('warning','Cleaning up FileListToImport.');
          sql = 'DELETE FROM FileListToImport';
          obj.sqlUpdate(sql);
        end

        function add_all_files_to_import_list(obj)
			% works only in unix or mac
            % Note: This function adds all files found on system to import
            % list, when importing only new (not previously imported files
            % will be processed). This is to allows for deleting files from
            % DB which has been deleted upstream and removed from system.
			% Verify sqlite3 is installed
			[status, ~] = system('command -v sqlite3 >/dev/null 2>&1 || { exit 100; }');
			if(status==100), error('It appears Sqlite3 is not installed/found on your system.'); end
			system(['cd ' obj.databaseDirectory ...
				'; find ./mms[1-4]* -name mms*cdf -type f |  ' ...
				'perl -pe ''s/(\.\/mms.*)(mms[1-4]?_[\w-]*)_(20\d\d\d\d\d\d\d*)_(v[\d\.]*)(.cdf)\n/$1,$2,$3,$4,$_/'' > delme.txt;'...
				'echo -e ".mod csv\n.import delme.txt FileListToImport\n" | sqlite3 ' obj.databaseFile ';'...
				'rm ./delme.txt'...
				]);
        end

        function status = import_file(obj,fileToImport)
			% import a single file
			% return status = 1 if file is imported sucessully, or it exists or
			% file with newer version exists. Otherwise return status = 0.
			%
			irf.log('notice',['File to import: ' fileToImport]);
			status = 0; % Assume it did not succeed yet.
			removeOlderVersionFile = false;
			% check if files with same data and date exist
			FileInfo = mms_db_sql.get_file_info(fileToImport);
			sql = ['select idFile,version from FileList ' ...
				'where directory = "' FileInfo.directory '"'...
				' and dataset = "' FileInfo.dataset '"'...
				' and date = "' FileInfo.date '"'];
			rs = obj.sqlQuery(sql);
			if rs.next % file with same dataset and date exists, compare versions
				existingVersion = char(rs.getString('version'));
				if is_version_larger(FileInfo.version(2:end),existingVersion(2:end))
					irf.log('notice',['File with older version ' existingVersion ' exsists!']);
					removeOlderVersionFile = true;
					existingFileID = char(rs.getString('idFile'));
				else
					irf.log('notice',['Not importing version ' FileInfo.version ...
						' because file with newer version ' existingVersion ' exsists!']);
					return;
				end
			end
			
			% add fileName to FileList and get idFile
			sql = ['insert into FileList (directory,dataset,date,version,fileNameFullPath)'...
				' SELECT * from FileListToImport where fileNameFullPath = "' fileToImport '"'];
			obj.sqlUpdate(sql);
			sql = ['select idFile,dataset from FileList where fileNameFullPath = "' fileToImport '"'];
			rs = obj.sqlQuery(sql);
			while rs.next % file already exist
				idFile = char(rs.getString('idFile'));
				dataset = char(rs.getString('dataset'));
			end
			
			% add to VarIndex list
			out = obj.get_science_variables(fileToImport);
			if isempty(out) % reading cdf file did not succeed
				status = 0;
				return;
			end
			for iDataset = 1:numel(out)
				if isempty(out(iDataset).startTT) || out(iDataset).startTT<1000, break;end % energy channels are put as DEPEND_0 for FPI
				varNames = out(iDataset).varNames;
				% add dataset to Datasets if needed
				idDataset=obj.add_var_names(dataset,varNames);
                % Sqlite does not have NaN but uses NULL values..
                if(isnan(out(iDataset).startTT) || isnan(out(iDataset).endTT))
                  irf.log('debug',['.. insert into VarIndex: idDataset=' dataset ...
					' : "NULL"/"NULL" ']);
				  sql = ['insert into VarIndex (idFile,idDataset,startTT,endTT) '...
					'values (' idFile ',"' idDataset '", "NULL", "NULL")'];
                else
                  irf.log('debug',['.. insert into VarIndex: idDataset=' dataset ...
					' : ' irf_time([out(iDataset).startTT out(iDataset).endTT],'tint>utc')]);
				  sql = ['insert into VarIndex (idFile,idDataset,startTT,endTT) '...
					'values (' idFile ',"' idDataset '",' ...
					num2str(out(iDataset).startTT) ',' num2str(out(iDataset).endTT) ')'];
                end
				obj.sqlUpdate(sql);
			end
			% If we have reached this point then insert went well
            status = 1;
			if removeOlderVersionFile
				irf.log('notice',['Deleting information of the file with older version ' existingVersion '!']);
				obj.sqlUpdate(['delete from FileList where idFile = "' existingFileID '"']);
			end
        end
		
        function idDataset = add_var_names(obj,dataset,varNames)
          % ADD_VAR_NAMES add variable names to VarNames table
          if ischar(varNames), varNames = {varNames}; end
          % Check if dataset with varNames exists
          varNamesCharArray = char(sort(varNames))';
          varNamesCharArray(:,end+1) = ' ';
          varNamesString = reshape(varNamesCharArray, 1, []);
          rs = obj.sqlQuery(['select idDataset,dataset from Datasets where varNames="' varNamesString '";']);
          flagAddDatasetToDb = true;
          while rs.next
            datasetDb = char(rs.getString('dataset'));
            idDataset = char(rs.getString('idDataset'));
            if strcmp(datasetDb,dataset)
              flagAddDatasetToDb = false;
              break;
            end
            irf.log('warning',['! ' datasetDb ' and ' dataset ' include the same variables!']);
          end

          % add dataset and variables if needed and get idDataset
          if flagAddDatasetToDb
            % add dataset
            obj.sqlUpdate(['insert into Datasets(dataset,varNames) values("' dataset '","' varNamesString '");']);
            rs = obj.sqlQuery(['select idDataset from Datasets where dataset="' dataset '";']);
            while rs.next
              idDataset = char(rs.getString('idDataset'));
            end
            for iVar=1:length(varNames)
              if(iVar==1)
                sqlValues = ['("', varNames{iVar}, '",',idDataset,')'];
              else
                sqlValues = [sqlValues, ', ("', varNames{iVar}, ...
                  '",',idDataset,')']; %#ok<AGROW>
              end
            end
            sqlValues = [sqlValues, ';'];
            sql = ['INSERT OR IGNORE INTO VarNames(varName,idDataset) VALUES', sqlValues];
            obj.sqlUpdate(sql);
          end
        end

        function [idDatasetList,DatasetList] = find_datasets_with_varname(obj,varName)
          % find Datasets with varName
          idDatasetList = {}; iDataset = 1;
          if(nargout==2), DatasetList = {}; end
          sql = ['SELECT idDataset,dataset FROM VarNames LEFT JOIN Datasets ', ...
            'USING (idDataset) WHERE varName = "' varName '"'];
          rs=obj.sqlQuery(sql);
          while rs.next
            idDatasetList{iDataset} = char(rs.getString('idDataset')); %#ok<AGROW>
            if(nargout==2)
              DatasetList{iDataset} =  char(rs.getString('dataset')); %#ok<AGROW>
            end
            iDataset = iDataset + 1;
          end
          if(isempty(idDatasetList))
            irf.log('warning',['There is no variable with name ' varName '.']);
          end
        end

        function idDatasetList = find_dataset_id(obj,dataset)
			% find Datasets with name "dataset"
			idDatasetList = {};iDataset = 1;
			sql = ['select idDataset from Datasets where dataset = "' dataset '"'];
			rs=obj.sqlQuery(sql);
			if ~rs.next
				irf.log('warning',['There is no dataset with name ' dataset '.']);
				return;
			else
				while true
					idDatasetList{iDataset} = char(rs.getString('idDataset'));  %#ok<AGROW>
					iDataset = iDataset +1 ;
					if ~rs.next, break; end
				end
			end
		end
		
		function tintArray = index_var(obj,varName)
            sql = ['SELECT startTT,endTT FROM VarIndex LEFT JOIN VarNames ', ...
              'USING (idDataset) WHERE VarNames.varName = "' varName '" ',...
              'ORDER BY startTT ASC'];
			rs= obj.sqlQuery(sql);
			tintArray = zeros(0,2,'int64');
			while rs.next
				tint = sscanf([char(rs.getString('startTT')) ' ' ...
					char(rs.getString('endTT'))],'%ld %ld');
				tintArray(end+1,:)= tint;       %#ok<AGROW>
			end
			if nargout == 0, % print time intervals
				nTint = size(tintArray,1);
				for ii = 1:min(nTint,5),       disp(irf_time(tintArray(ii,:),'tint>utc')); end
				if nTint>10,                   disp('...'); end
				for ii = max(5,nTint-5):nTint, disp(irf_time(tintArray(ii,:),'tint>utc')); end
				clear tintArray;
			end
		end
		function res = file_has_var(obj,fileName,varName)
			% find files
			if ischar(varName), varName={varName};end
			res=true; % default is true
			for iVarname = 1:length(varName)
				varString = varName{iVarname};
				sql = ['select v.idVar from VarNames AS v where v.varName = "' varString '" and '...
					'v.idDataset IN (select vind.idDataset from VarIndex AS vind where vind.idFile='...
					'(select fl.idFile from FileList AS fl where fl.fileNameFullPath = "' fileName '"))'];
				rs=obj.sqlQuery(sql);
				if ~rs.next
					res = false;
					return;
				end
			end
		end
		
		function fileNames = find_files(obj,varargin)
			% FIND FILES search files given variable name and/or dataset and/or time interval
			%
			% FIND_FILES(obj,'varName', variableName,..)
			% FIND_FILES(obj,'dataset', datasetName,..)
			% FIND__FILES(obj,'tint,   , tint,..)
			%  time interval can be UTC string or GenericTimeArray
%			startTT = int64(0);
%			endTT = intmax('int64');
			fileNames = cell(0);
			searchTint = false;
			searchVariable = false;
			searchDataset = false;
			args = varargin;
			while numel(args)>=2
				switch lower(args{1})
					case {'tint'}
						timeInterval = args{2};
						[startTT,endTT] = mms_db_sql.start_stop_in_ttns(timeInterval);
						searchTint = true;
					case {'varname','variable'}
						varName = args{2};
						if ~ischar(varName),
							irf.log('critical','varName should be text string');
							return;
						end
						searchVariable = true;
					case {'dataset','fileprefix'}
						dataset = args{2};
						if ~ischar(dataset),
							irf.log('critical','varName should be text string');
							return;
						end
						searchDataset = true;
					otherwise
						irf.log('critical','unrecognized input');
						return;
				end
				args(1:2)=[];
			end
			
			if searchTint
				sqlTime = [' and startTT <= ' num2str(endTT) ...
					         ' and   endTT >= ' num2str(startTT) ];
			else
				sqlTime = '';
			end
			
			% find Datasets with varName
			if ~searchVariable && ~searchDataset
				sqlDataset = 'idDataset ';
			else
				if searchVariable && ~searchDataset
					idDatasetList = find_datasets_with_varname(obj,varName);
				elseif ~searchVariable && searchDataset
					idDatasetList = find_dataset_id(obj,dataset);
				elseif searchVariable && searchDataset ;
					idDatasetList = intersect(find_datasets_with_varname(obj,varName),find_dataset_id(obj,dataset));
				end
				sqlDataset = ['idDataset IN ("', strjoin(idDatasetList, '","'), '")'];
			end
			% find files
			idFileArray = []; iFile = 1;
			sql = ['select idFile from VarIndex where ' sqlDataset sqlTime ...
				' order by startTT asc'];
			rs=obj.sqlQuery(sql);
			while rs.next
				idFileArray(iFile) = str2double(rs.getString('idFile')); %#ok<AGROW>
				irf.log('debug',['idFile = ' num2str(idFileArray(iFile))]);
				iFile = iFile + 1;
			end
			
			% get filenames
			fileNames = cell(numel(idFileArray),1);
			for iFile = 1:numel(idFileArray)
				sql = ['select fileNameFullPath from FileList where idFile = ' ...
					num2str(idFileArray(iFile))];
				rs = obj.sqlQuery(sql);
				if rs.next
					fileNames{iFile} = char(rs.getString('fileNameFullPath'));
				end
			end
		end
				
		function var(obj,varargin)
			% VAR seach variable names
			%  VAR('par1','par2',..)
			%  finds variable name which contains a text string par1 and par2 and ...
			sql = 'select * from VarNames ';
			if nargin > 1,
				sql = [sql 'where '];
				for iVar = 1:nargin-1,
					if ischar(varargin{iVar})
						if iVar > 1,
							sql = [sql ' and ']; %#ok<AGROW>
						end
						sql = [sql ' varName like "%' strrep(varargin{iVar},'*','%') '%"']; %#ok<AGROW>
					end
				end
			end
			objName = inputname(1);
			res=obj.sqlQuery(sql);
			while res.next
%				disp([char(res.getString('varName'))]);
				varName = char(res.getString('varName'));
				linkTxt=mms_db_sql.matlab_link(varName,[objName '.var_attributes(''' varName ''')']);
				disp(linkTxt{1});
			end
		end
		
		function var_attributes(obj,varName)
			fileNameList=obj.find_files('varname',varName);
			fileName=fileNameList{1};
			d=spdfcdfinfo(fileName,'varstruct',true);
			ivar=ismember(d.Variables.Name,varName);
			varProp=d.Variables;
			fieldN=fieldnames(varProp);
			for ii = 1: numel(fieldN)
				varProp.(fieldN{ii})(~ivar)=[];
			end
			disp '------ VARIABLE PROPERTIES -------';
			disp(varProp);
			varAttr=d.VariableAttributes;
			fieldN=fieldnames(varAttr);
			for ii = 1: numel(fieldN)
				i=strcmp(varName,varAttr.(fieldN{ii})(:,1));
				if isempty(i) || ~any(i)
					varAttr=rmfield(varAttr,fieldN{ii});
				else
					varAttr.(fieldN{ii})=varAttr.(fieldN{ii})(i,2);
				end
			end
			disp '------ VARIABLE ATTRIBUTES -------';
			disp(varAttr);
		end
		
		function out = sqlQuery(obj,sql)
			obj.connect;
			irf.log('debug',['sqlite> ' sql]);
			out=obj.statement.executeQuery(sql);
		end
		
		function out = sqlUpdate(obj,sql)
			obj.connect;
			irf.log('debug',['sqlite> ' sql]);
			out=obj.statement.executeUpdate(sql);
		end
	end
	methods (Static)
		function out = get_science_variables(cdfFileName)
			% OUT returns structure array with fields
			%     epochVarName - name of the epoch variable, string
			%     varNames - cell array of strings with variables
			%     startTT  - int64 (TT2000 values) with corresponding start time
			%     endTT    - int64 (TT2000 values) with corresponding end time
			narginchk(1,1); % One file name
            if(~exist(cdfFileName, 'file'))
              errStr = ['File not found: ', cdfFileName];
              irf.log('warning', errStr); warning(errStr);
              out = []; return;
            end
            irf.log('debug',['Reading: ' cdfFileName]);
            try
              inf = spdfcdfinfo(cdfFileName, 'VARSTRUCT', true);
            catch ME
              errStr = ['Cannot get file information from: ', cdfFileName];
              irf.log('warning', errStr); irf.log('warning', ME.message);
              warning(errStr); out = []; return;
            end
            dep = inf.VariableAttributes.DEPEND_0;
            [tVarNames, ~, IC] = unique(dep(:,2));
            indGoodTVarName = 1:numel(tVarNames);
            iBadTVarName = ~ismember(tVarNames, inf.Variables.Name); % time variable should be a variable in the cdf file
            if any(iBadTVarName)
              irf.log('notice',['!! bad time variable names in DEPEND_O: ' tVarNames{iBadTVarName}]);
              indGoodTVarName(iBadTVarName) = [];
              tVarNames{iBadTVarName} = [];
            end
            isBadTime = cellfun(@(x) ~any(strcmpi({'tt2000', 'epoch', 'epoch16'}, ...
              inf.Variables.DataType(strcmp(inf.Variables.Name, x)))), ...
              tVarNames(indGoodTVarName));
            if any(isBadTime)
              irf.log('notice',['! not accepted time format for time DEPEND_O variable: ' ...
                tVarNames{indGoodTVarName(isBadTime)}]);
              indGoodTVarName(isBadTime) = [];
            end
            try
              epoch = spdfcdfread(cdfFileName, 'DataOnly', true, ...
                'Variable', tVarNames(indGoodTVarName), ...
                'KeepEpochAsIs', true);
            catch ME
              errStr = ['Cannot read Epochs from file: ', cdfFileName];
              irf.log('warning', errStr); irf.log('warning', ME.message);
              warning(errStr); out = []; return;
            end
            if isinteger(epoch), epoch = {epoch};end % only one time variable
            for iT = numel(indGoodTVarName):-1:1
              out(iT).epochVarName = tVarNames{iT};
              out(iT).varNames = dep(IC == iT,1);
              startTT = min(epoch{indGoodTVarName(iT)});
              endTT   = max(epoch{indGoodTVarName(iT)});
              if any(startTT) && any(endTT) && ...
                  (min(startTT,endTT) < int64(479390467184000000)...% '2015-03-12T00:00:00.000000000'
                  ||  max(startTT,endTT) > int64(1262260868184000000)),%'2040-01-01T00:00:00.000000000'
                out(iT).startTT = NaN;
                out(iT).endTT   = NaN;
                irf.log('notice',['!!! In file:' cdfFileName]);
                irf.log('notice',['!!! ' out(iT).epochVarName ' has startTT = ' num2str(startTT) ' and endTT = ' num2str(endTT)]);
              else
                out(iT).startTT = startTT;
                out(iT).endTT   = endTT;
              end
            end
        end

		function [startTT,endTT] = start_stop_in_ttns(timeInterval)
			if ischar(timeInterval)
				timeInterval = irf.tint(timeInterval);
			end
			if isa(timeInterval,'GenericTimeArray')
				startTT = timeInterval.start.ttns;
				endTT = timeInterval.stop.ttns;
			elseif isinteger(timeInterval)
				startTT = timeInterval(1);
				endTT = timeInterval(2);
			else
				irf.log('warning','Time interval in unknown format');
				return;
			end
		end
		function fileInfo = get_file_info(fileName)
			% GET_FILE_INFO get values of directory, dataset, date, version
			% fileInfo is structure with fields
			% 'directory','dataset','date',version'
			fileInfo = regexp(fileName,['(?<directory>\.\/mms.*)'...
				'(?<dataset>mms[1-4]?_[\w-]*)_(?<date>20\d\d\d\d\d\d\d*)'...
				'_(?<version>v[\d\.]*)(.cdf)'],'names');
		end
		function outStr=matlab_link(linkText,linkCommandText)
			% MATLAB_LINK returns string with link to matlab text to execute
			% 
			% MATLAB_LINK(linkText,linkCommandText)
			%  linkText can be cell array
			%  In linkCommandText the '?' is substituted with linkText

			if ischar(linkText)
				linkText = {linkText};
			end
			if iscell(linkText)
				outStr=cell(size(linkText));
				for ii = 1:numel(linkText)
					linkStr = linkText{ii};
					linkCommandStr = strrep(linkCommandText,'?',linkStr);
					outStr{ii} = ['<a href="matlab: ' linkCommandStr '">' linkStr '</a>' ];
				end
			end
		end
	end
	
end
