classdef TimeTable
  % irf.TimeTable class to handle time tables
  %
  % TT=irf.TimeTable	- generate empty time table
  % TT=irf.TimeTable(tint) - generate time table from matrix tint.
  %			1st column start time and 2nd column end times in isdat epoch
  %     or nanoseconds of TT2000 if tint is int64
  % TT=irf.TimeTable('demo') - generate short example time table
  % TT=irf.TimeTable('https://..') - load ascii time table from internet link
  % TT=irf.TimeTable(filename) - load ascii time table from file
  % TT=irf.TimeTable(ascii_text) - load time table from cell string array ascii_text
  %
  % Methods: (see also help irf.TimeTable.(method))
  %   TT  = add(TT,..)       add time interval to table, see help irf.TimeTable.add
  %   out = ascii(TT)	       time table in ascii format
  %   TT	= common(T1,T2)    return common elements
  %   TT  = intersect(T1,T2) intersection of two time tables
  %   N   = numel(TT)        number of time intervals
  %   TT  = overlap(T1,T2)   returns elements in T1 that overlap with T2
  %   T2  = remove(T1,index) remove elements
  %   T2  = select(T1,index) select elements
  %   TT  = setdiff(T1,T2)   returns elements of T1 that are not in T2
  %   T2  = sort(T1)         sort according to start times
  %   T2  = unique(T1)       return unique time table (sorted nonoverlapping intervals)
  %    export_ascii(TT,filename)   export time table to ascii file
  %
  % Examples:
  %	TT = irf.TimeTable % initialize empty time table
  %	tint=irf_time([2001 1 1 0 0 0]) + [0 3600];
  %	TT=add(TT,tint,{'Burst events','all tail crossings 2001'},'this is best event');
  %	TT=add(TT,tint+7200,{},'this is another best');
  %	ascii(TT)
  %
  %   %read Cluster Master Science Plan
  %   msp=irf.TimeTable('http://jsoc1.bnsc.rl.ac.uk/msp/full_msp_ascii.lst');

  % The following properties can be set only by class methods
  properties (SetAccess = public)
    Header
    TimeInterval
    Comment
    Description
    UserData
  end
  methods
    function TT		= TimeTable(varargin) % construct time table
      TT.Header={};		%	header for file, cell string array
      TT.TimeInterval=[];	%	matrix with 2 columns, start/end
      TT.Comment={};		%	short comment for each time interval
      TT.Description={};	%	description lines for each interval
      TT.UserData=struct([]);	%	user can put different info in this structure
      if numel(varargin)==1
        source = varargin{1};
      elseif numel(varargin)==0
        return;
      else
        irf.log('critical','max 1 argument supported');
        return
      end
      if ischar(source) && size(source,1)==1
        if strcmpi(source,'demo') % return demo time table
          TT = irf.TimeTable; % initialize empty time table
          tint=irf_time([2001 1 1 0 0 0]) + [0 3600];
          TT.Header = {'This is header','Where general file information is given'};
          TT=add(TT,tint,{'These are several lines','describing the time interval'},'this is comment');
          TT=add(TT,tint+7200,{},'this is comment, there is no description for this time interval');
          TT=add(TT,tint+600,{},'this is overlapping interval');
        elseif regexp(source,'^https?://')==1 % start with http:// or https://
          tempTT = urlread(source); %#ok<URLRD> webread introduced in R2014b
          TT	 = irf.TimeTable(tempTT);
        elseif exist(source,'file') % check if file
          fid = fopen(source);
          fileContents=fscanf(fid,'%c',inf);
          fclose(fid);
          TT=irf.TimeTable(fileContents);
        else % assume source is ascii text of time table
          ttAscii=textscan(source,'%s','delimiter',sprintf('\n')); %#ok<SPRINTFN>
          TT=irf.TimeTable(ttAscii{1});
          return;
        end
      elseif iscellstr(source) % source is ascii text to be converted
        % TT   = from_ascii(ascii)
        % convert ascii to time table
        % function out=from_ascii(ttascii) % convert ascii tt to matlab format
        TT=irf.TimeTable;

        isHeader					= 1;
        iHeaderLine					= 0;
        nTimeInterval				= 0;
        header				    	= cell(size(source,1));
        description			    	= cell(size(source,1));
        iDescriptionLine			= 0;
        for iTtAscii=1:numel(source)
          str=source{iTtAscii};
          timeInterval=regexp(str,'^\s*(?<start>[\d-]*T[\d:\.]*Z?)\s*(?<end>[\d-]*T[\d:\.]*Z?)\s?(?<comment>.*)','names');
          if ~isempty(timeInterval) && isfield(timeInterval,'start') && isfield(timeInterval,'end')
            isHeader=0;
            nTimeInterval=nTimeInterval+1;
            TT.TimeInterval(nTimeInterval,:) = [irf_time(timeInterval.start,'utc>epoch') irf_time(timeInterval.end,'utc>epoch')];
            if isfield(timeInterval,'comment')
              TT.Comment{nTimeInterval}	= timeInterval.comment;
            else
              TT.Comment{nTimeInterval}	= [];
            end
            if nTimeInterval>1
              TT.Description{nTimeInterval-1}	= description(1:iDescriptionLine);
            end
            iDescriptionLine				= 0;  % reset counter of description lines for the next time interval
          else
            if isHeader
              iHeaderLine = iHeaderLine + 1;
              header{iHeaderLine} = str;
            else
              iDescriptionLine = iDescriptionLine + 1;
              description{iDescriptionLine} = str;
            end
          end
        end
        TT.Description{nTimeInterval}=description(1:iDescriptionLine);
        TT.Header = header(1:iHeaderLine);
      elseif isnumeric(source) && size(source,2)==2 % time interval
        TT.TimeInterval	= source;
        TT.Header		= {};
        TT.Comment		= cell(size(source,1),1);
        TT.Description	= cell(size(source,1),1);
      else
        irf.log('critical','Unknown argument.');
        return;
      end
    end
    function N		= numel(TT)
      N = numel(TT.TimeInterval)/2;
    end
    function TT		= add(TT,tint,description,comment)
      % TT = add(TT,tint,description,comment)
      % add element to time table
      % tint - time interval, vector with start and end time in epoch
      %			or time interval in iso format
      % description - cell string array describing time interval
      % comment - character vector commenting time interval (in ascii format
      %				located on the same line)
      if nargin<4, comment='';end
      if nargin<3, description={};end
      if nargin<2, irf.log('critical','not enough arguments');return;end
      nElement=size(TT.TimeInterval,1)+1;
      if isnumeric(tint)
        if numel(tint)==0
          irf.log('critical','Time interval not specified. Returning.');
          return;
        elseif numel(tint)==1 % only time instant give, end the same as start time
          tint(2)=tint(1);
        end
      elseif ischar(tint) % assume iso format
        tint=irf_time(tint,'utc>tint');
      end
      TT.TimeInterval(nElement,:)=tint;
      if ~iscellstr(description) % description should be cellstr
        if ischar(description)
          description={description};
        else
        end
      end
      TT.Description{nElement}=description;
      if ischar(comment)
        TT.Comment{nElement}=comment;
      else
        TT.Comment{nElement}=[];
      end
    end
    function out	= ascii(tt) % tt in ascii format
      % Return time table in ascii format (cell string array)
      header=tt.Header;
      tint=tt.TimeInterval;
      description=tt.Description;
      comment=tt.Comment;
      nElements=numel(tt);
      if isempty(comment), comment=cell(nElements,1);end
      if isempty(description), description=cell(nElements,1);end
      nHeaderLines=numel(header);
      nDescriptionLines = sum(cellfun(@(x) numel(x),description));
      out=cell(nHeaderLines+nElements+nDescriptionLines,1);
      out(1:nHeaderLines)=header;
      currentLine=nHeaderLines+1;
      fmt='%s %s %s';
      if isa(tint,'int64')
        tstartIso = irf_time(tint(:,1),'ttns>utc');
        tendIso = irf_time(tint(:,2),'ttns>utc');
      else
        tstartIso = irf_time(tint(:,1),'utc');
        tendIso = irf_time(tint(:,2),'utc');
      end
      for iTT=1:numel(tt)
        out{currentLine}=sprintf(fmt,tstartIso(iTT,:),...
          tendIso(iTT,:),comment{iTT});
        currentLine=currentLine+1;
        for j=1:numel(description{iTT})
          out{currentLine}=description{iTT}{j};
          currentLine=currentLine+1;
        end
      end
      if nargout == 0 % interactive calling from command line
        out=char(out);
        if size(out,1)>nHeaderLines+10
          disp(out(1:nHeaderLines+3,:));
          disp(['.... ' num2str(size(out,1)-nHeaderLines-5) ' lines omitted.....']);
          disp(out(end-2:end,:));
        else
          disp(out);
        end
        clear out;
      end
    end
    function ok		= export_ascii(TT,filename) % export to file
      % ok = export_ascii(TT,filename)
      %	export time table TT to file
      out=ascii(TT);
      if ischar(filename) % read from file
        fid = fopen(filename,'w');
        if fid == -1
          irf.log('critical',['Cannot open file:' filename]);
          ok=0;
          return;
        end
        for j=1:numel(out)
          fprintf(fid,'%s\n',out{j});
        end
        fclose(fid);
        ok=1;
        if nargout==0, clear ok, end
      end
    end
    function TTout  = remove(TTin,index) % remove entries
      % TTout  = remove(TTin,index) remove entries
      if isempty(index)
        TTout=TTin;
        return;
      end
      TTout=[];
      if ~isnumeric(index)
        irf.log('critical','Index not number');
        return;
      end
      if max(index(:)) > numel(TTin) || min(index(:)) < 1
        irf.log('critical','Index out of range');
        return;
      end
      TTout				           = TTin;
      TTout.TimeInterval(index(:),:) = [];
      if ~isempty(TTout.Comment)
        TTout.Comment(index(:))     = [];
      end
      if ~isempty(TTout.Description)
        TTout.Description(index(:))    = [];
      end
      TTout.Header                   = TTin.Header;
      if ~isempty(TTout.UserData)
        TTout.UserData(index)          = [];
      end
    end
    function TTout = select(TTin,index) % return index
      if islogical(index)
        if numel(index) == size(TTin.TimeInterval,1)
          index = find(index);
        else
          error('select(TTin,index): index wrong size');
        end
      end
      if ~isnumeric(index)
        irf.log('critical','Index not number');
        return;
      elseif isempty(index)
        TTout=irf.TimeTable;
        irf.log('critical','index empty, returning empty timetable');
        return;
      end
      if max(index(:)) > numel(TTin) || min(index(:)) < 1
        irf.log('critical','Index out of range');
        return;
      end
      TTout				= irf.TimeTable(TTin.TimeInterval(index(:),:));
      TTout.Header		= TTin.Header;
      if ~isempty(TTin.Comment),		TTout.Comment		= TTin.Comment(index(:));		end
      if ~isempty(TTin.Description),	TTout.Description   = TTin.Description(index(:));	end
      if ~isempty(TTin.UserData),	TTout.UserData		= TTin.UserData(index);			end
    end
    function TTout	= sort(TTin) % sort according to start time
      TTout=TTin;
      [~,inew]=sort(TTin.TimeInterval(:,1));
      TTout.TimeInterval = TTin.TimeInterval(inew,:);
      TTout.Comment = TTin.Comment(inew);
      TTout.Description = TTin.Description;
      if isfield(TTout,'UserData')
        TTout.UserData = TTin.UserData(inew);
      end
    end
    function TTout	= unique(TTin) % time table sorted and unique
      TTout=TTin;
      tt=TTin.TimeInterval;
      com=TTin.Comment;
      desc=TTin.Description;
      [~,inew]=sort(tt(:,1));
      ttnew=tt(inew,:);
      iJoinedLines=zeros(numel(inew),1);
      for j=numel(inew):-1:2
        if ttnew(j,1)<ttnew(j-1,2)
          iJoinedLines(j)=1;
          ttnew(j-1,2)=max(ttnew(j-1,2),ttnew(j,2));
          ttnew(j,:)=[];
        end
      end
      TTout.TimeInterval = ttnew;
      if ~isempty(com)
        com = com(inew);
        TTout.Comment = com(~iJoinedLines);
      end
      if ~isempty(desc)
        desc = desc(inew);
        TTout.Description = desc(~iJoinedLines);
      end
    end
    function TT		= intersect(TT1,TT2) % intersect of 2 time tables
      TT1 = unique(TT1);
      TT2 = unique(TT2);
      t1=TT1.TimeInterval;
      t2=TT2.TimeInterval;
      tout=zeros(max(size(t1,1),size(t2,1))+1,2);
      iInterval=0;
      i2=1;
      endOfTT2 = false; % flag to know when end of time series 2
      for j= 1:size(t1,1)
        while (t1(j,1)>t2(i2,2))
          i2=i2+1;
          endOfTT2 = i2>size(t2,1); % check end of 2nd time series
          if endOfTT2, break; end % if end of 2nd time series stop
        end
        if endOfTT2, break; end % if end of 2nd time series stop
        if t1(j,2) < t2(i2,1)
          continue;
        end
        while t2(i2,1)<t1(j,2)
          iInterval=iInterval+1;
          tout(iInterval,:)=[max(t1(j,1),t2(i2,1)) min(t1(j,2),t2(i2,2))];
          if t1(j,2) > t2(i2,2) % take next i2
            i2=i2+1;
            endOfTT2 = i2>size(t2,1); % check end of 2nd time series
            if endOfTT2, break; end % if end of 2nd time series stop
          else
            break;
          end
        end
        if endOfTT2, break; end % if end of 2nd time series stop
      end
      if iInterval == 0
        irf.log('warning','There were no intersecting intervals, returning empty time table.');
      end
      tout(iInterval+1:end,:)=[];
      TT = irf.TimeTable;
      TT.Header = [cellfun(@(x) {['# ' x]}, TT1.Header);...
        'INTERSECT';cellfun(@(x) {['# ' x]}, TT2.Header)];
      TT.TimeInterval=tout;
      TT.Comment=cell(1,size(tout,1));
      TT.Description=cell(1,size(tout,1));
    end
    function [TT,i]	= setdiff(TT1,TT2) % returns TT1 elements that are not in TT2
      % TT=setdiff(TT1,TT2) returns TT1 elements that are not in TT2
      % [TT,ii]=setdiff(TT1,TT2) returns also indices of TT1 that are different
      t1=TT1.TimeInterval;
      t2=TT2.TimeInterval;
      [~,ia1,~]=intersect(t1(:,1),t2(:,1));
      [~,ia2,]=intersect(t1(:,2),t2(:,2));
      ia=intersect(ia1,ia2); % ia intervals are common in both series
      i=setdiff(1:size(t1,1),ia);
      TT=select(TT1,i);
    end
    function [TT,i1]	= overlap(TT1,TT2) % returns elements in TT1 which overlap with elements in TT2
      % [TT,i] = overlap(TT1,TT2) % returns elements in TT1 which overlap with elements in TT2
      % i is indices of overlapping intervals of TT1
      t1=TT1.TimeInterval;
      t2=TT2.TimeInterval;
      if isempty(t1) || isempty(t2)
        TT=irf.TimeTable;
        i1=[];
        return;
      end
      i1 = false(1,size(t1,1));
      i1Open = i1;
      indt1 = (1:size(t1,1))';
      t1List = [t1(:,1) 0*indt1+1 0*indt1 indt1;...
        t1(:,2) 0*indt1-1 0*indt1 indt1];
      indt2 = (1:size(t2,1))';
      t2List = [t2(:,1) 0*indt2 0*indt2+1 indt2;...
        t2(:,2) 0*indt2 0*indt2-1 indt2];
      tList = sortrows([t1List;t2List]);
      tList(:,5) = cumsum(tList(:,2)); % number of open TT1 intervals
      tList(:,6) = cumsum(tList(:,3)); % number of open TT2 intervals
      for i = 1:size(tList,1)
        if tList(i,2) == 0 % TT2 interval
          if tList(i,3) == -1 && tList(i,6) == 0 % all TT2 intervals are closed
            if tList(i,5) > 0 % there are open TT1 intervals
              i1(i1Open) = true; % all open TT1 intervals are overlapping
            end
          end
        else %% TT1 interval
          ind=tList(i,4);
          if tList(i,2) == 1 % interval start
            i1Open(ind)=true;
          else % interval end
            i1Open(ind)=false;
            if tList(i,6) > 0, i1(ind)=true;end
          end
        end
      end
      i1=find(i1);
      TT=select(TT1,i1);
    end
    function [TT,ii1,ii2]	= common(TT1,TT2) % returns TT1 elements that are commmon in TT2			% TT=common(TT1,TT2)
      % [TT,ii1,ii2]=common(TT1,TT2) returns also indexes of those elements
      % The UserData information is taken from TT1
      t1=TT1.TimeInterval;
      t2=TT2.TimeInterval;
      [~,ia1,ib1]=intersect(t1(:,1),t2(:,1)); % common start times
      [~,ia2,ib2]=intersect(t1(:,2),t2(:,2)); % common end times
      ii1=intersect(ia1,ia2); % TT1 intervals with start/end times common
      ii2=intersect(ib1,ib2); % TT2 intervals with start/end times common
      TT=select(TT1,ii1);
    end
  end % methods
end % classdef
