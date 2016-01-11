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
							sql = [sql ' and '];
						end
						sql = [sql ' varName like "%' varargin{iVar} '%"'];
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
			%     varNames - cell array of strings
			%     startTT  - int64
			%     endTT    - int64
			irf.log('debug',['Reading: ' cdfFileName]);
			inf = spdfcdfinfo(cdfFileName);
			dep = inf.VariableAttributes.DEPEND_0;
			varNames = dep(:,1);
			nVar = numel(varNames);
			[tVarNames,~,IC] = unique(dep(:,2));
			startTT = ones(nVar,1,'int64');
			endTT = startTT;
			% read time variables
			epoch = spdfcdfread(cdfFileName, ...
				'Variable', tVarNames, 'KeepEpochAsIs', true);
			if isinteger(epoch), epoch = {epoch};end % only one time variable
			for iT = 1:numel(tVarNames)
				if isempty(epoch{iT}), epoch{iT}=NaN;end
				startTT(IC==iT) = min(epoch{iT});
				endTT(IC==iT)   = max(epoch{iT});
			end
		end
	end
	
end
