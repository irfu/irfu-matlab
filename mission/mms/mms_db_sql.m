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
		function obj = mdb(fileName)
			if nargin == 1,
				[dir,file,ext] = fileparts(fileName);
				if isempty(dir)
					dir = pwd;
				end
				obj.databaseDirectory = dir;
				obj.databaseFile = [dir filesep file ext];
				if ~exist(obj.databaseFile,'file')
					obj.create_new_database;
				end
			end
		end
		
		function create_new_database(obj)
			sql = ['CREATE TABLE "FileList" ('...
				'"fileId"	INTEGER, "fileNameFullPath"	TEXT,'...
				'	PRIMARY KEY(fileId) )'];
			obj.sqlUpdate(sql);
			sql = 'CREATE TABLE "FileListToImport" ("fileNameFullPath"	TEXT UNIQUE)';
			obj.sqlUpdate(sql);
			sql = ['CREATE TABLE "VarIndex" ('...
				'"fileId"	INTEGER,"varName"	TEXT,"startTT"	INTEGER,"endTT"	INTEGER)'];
			obj.sqlUpdate(sql);
			sql = ['CREATE TABLE "VarNames" ('...
				'"varId"	INTEGER,"varName"	TEXT UNIQUE,PRIMARY KEY(varId))'];
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
		
		function import_files_from_list(obj)
			% import all files from FileListToImport
			while 1
				sql = 'select * from FileListToImport';
				rs = obj.sqlQuery(sql);
				if rs.next
					fileToImport = char(rs.getString('fileNameFullPath'));
					irf.log('debug',['Importing ' fileToImport]);
					obj.import_file(fileToImport);
					% remove file from FileListToImport
					sql = ['delete from FileListToImport where fileNameFullPath="' fileToImport '"'];
					obj.sqlUpdate(sql);
				else
					irf.log('warning','All files imported from FileListToImport');
					break;
				end
			end
		end
		
		function add_all_files_to_import_list(obj)
			% works only in unix or mac
			% Verify sqlite3 is installed
			[status, ~] = system('command -v sqlite3 >/dev/null 2>&1 || { exit 100; }');
		        if(status==100), error('It appears Sqlite3 is not installed/found on your system.'); end
			system([' find ' obj.databaseDirectory ' -name *cdf -type f | ' ...
				'sqlite3 ' obj.databaseFile ' ".import /dev/stdin FileListToImport"']);
		end
		
		function import_a_file_from_list(obj)
			% import one new file from FileListToImport
			sql = 'select * from FileListToImport LIMIT 1';
			rs = obj.sqlQuery(sql);
			while rs.next
				fileToImport = char(rs.getString('fileNameFullPath'));
				obj.import_file(fileToImport);
			end
			% remove file from FileListToImport
			sql = ['delete from FileListToImport where fileNameFullPath="' fileToImport '"'];
			obj.sqlUpdate(sql);
		end
		
		function import_file(obj,fileToImport)
			% import a file
			irf.log('notice',['File to import: ' fileToImport]);
			% check if file is not in db
			sql = ['select * from FileList ' ...
				'where fileNameFullPath = "' fileToImport '"'];
			rs=obj.sqlQuery(sql);
			while rs.next % file already exist
				irf.log('warning','File exists!');
				obj.close;
				return;
			end
			% add file to FileList
			sql = ['insert into FileList (fileNameFullPath) values ("' fileToImport '")'];
			obj.sqlUpdate(sql);
			sql = ['select fileId from FileList where  fileNameFullPath = "' fileToImport '"'];
			rs=obj.sqlQuery(sql);
			while rs.next % file already exist
				fileId = char(rs.getString('fileId'));
			end
			%
			[varNames,startTT,endTT] = obj.get_science_variables(fileToImport);
			obj.add_var_names(varNames);
			for iVar = 1:numel(varNames)
				irf.log('debug',['.. insert into VarIndex: ' varNames{iVar}...
					' : ' irf_time(EpochTT([startTT(iVar) endTT(iVar)]),'tint>utc')]);
				sql = ['insert into VarIndex (fileId,varName,startTT,endTT) '...
					'values (' fileId ',"' varNames{iVar} '",' ...
					num2str(startTT(iVar)) ',' num2str(endTT(iVar)) ')'];
				obj.sqlUpdate(sql);
			end
		end
		
		function add_var_names(obj,varNames)
			if ischar(varNames), varNames = {varNames};end
			for iVar = 1:length(varNames)
				sql = ['insert or ignore into VarNames(varName) values("' varNames{iVar} '")'];
				obj.sqlUpdate(sql);
			end
		end
		
		function fileNames = search_files(obj,varName,startTT,endTT)
			if nargin < 3
				startTT = int64(0);
			end
			if nargin < 4
				endTT = intmax('int64');
			end
			if ~ischar(varName),
				irf.log('critical','varName should be text string');
				return;
			elseif ~isinteger(startTT) && ~isinteger(endTT),
				irf.log('critical','startTT and endTT should be integers');
				return;
			end
			sql = ['select * from VarNames where varName = "' varName '"'];
			rs=obj.sqlQuery(sql);
			if ~rs.next
				irf.log('warning',['There is no variable with name ' varName '.']);
				return;
			else
				sql = ['select fileId from VarIndex where varName = "' varName '"'];
				sql = [sql ' and startTT <= ' num2str(endTT)   ];
				sql = [sql ' and   endTT >= ' num2str(startTT) ];
				rs=obj.sqlQuery(sql);
				fileIdArray = []; iFile = 1;
				while rs.next
					fileIdArray(iFile) = str2num(rs.getString('fileId')); %#ok<AGROW>
					irf.log('debug',['fileId = ' num2str(fileIdArray(iFile))]);
					iFile = iFile + 1;
				end
				fileNames = cell(numel(fileIdArray),1);
				for iFile = 1:numel(fileIdArray)
					sql = ['select fileNameFullPath from FileList where fileId = ' ...
						num2str(fileIdArray(iFile))];
					rs = obj.sqlQuery(sql);
					if rs.next
						fileNames{iFile} = char(rs.getString('fileNameFullPath'));
					end
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
						sql = [sql ' varName like "%' varargin{iVar} '%"']; %#ok<AGROW>
					end
				end
			end
			res=obj.sqlQuery(sql);
			while res.next
				disp(char(res.getString('varName')));
			end
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
      function [varNames,startTT,endTT] = get_science_variables(cdfFileName)
        % Get science variable names and duration from "cdfFileName".
        % varNames - cell array of strings with variables
        % startTT  - int64 (TT2000 values) with corresponding start time
        % endTT    - int64 (TT2000 values) with corresponding end time
        narginchk(1,1); % One file name
        if(~exist(cdfFileName, 'file'))
          errStr = ['File not found: ', cdfFileName];
          irf.log('critical', errStr); error(errStr);
        end
        varNames=[]; startTT=[]; endTT=[];
        irf.log('debug',['Reading: ' cdfFileName]);
        try
          inf = spdfcdfinfo(cdfFileName);
        catch ME
          errStr = ['Cannot get file information from: ', cdfFileName];
          irf.log('warning', errStr); irf.log('warning', ME.message);
          warning(errStr); return;
        end
        dep = inf.VariableAttributes.DEPEND_0;
        varNames = dep(:,1);
        nVar = numel(varNames);
        [tVarNames,~,IC] = unique(dep(:,2));
        startTT = ones(nVar,1,'int64');
        endTT = startTT;
        % read time variables
        try
          epoch = spdfcdfread(cdfFileName, 'Variable', tVarNames, ...
            'KeepEpochAsIs', true);
        catch ME
          errStr = ['Cannot read Epochs from file: ', cdfFileName];
          irf.log('warning', errStr); irf.log('warning', ME.message);
          warning(errStr); return;
        end
        if isinteger(epoch), epoch = {epoch}; end % only one time variable
        % Some files have incorrect Epcoh or FillVal (for instance some
        % hk101 sunpulses) so only trust Epochs inside of interval
        % 2015-03-12T00:00:00.000000000 to 2040-01-01T00:00:00.000000000
        % (which is well within MMS mission life but should discard strange
        % epoch such as 1706-... and other invalid values).
        validEpoch = irf.tint('2015-03-12T00:00:00.000000000', ...
          '2040-01-01T00:00:00.000000000');
        for iT = 1:numel(tVarNames)
          if isempty(epoch{iT}), epoch{iT} = NaN; end
          epoch{iT}(epoch{iT}<validEpoch.start.ttns) = NaN; % Note, int64(NaN) = 0.
          epoch{iT}(epoch{iT}>validEpoch.stop.ttns) = NaN;
          startTT(IC==iT) = min(epoch{iT});
          endTT(IC==iT)   = max(epoch{iT});
        end
      end % get_science_variables

    end

end