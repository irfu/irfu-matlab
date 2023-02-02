function out=tt(varargin)
% IRF.TT work with time tables
%   time table is a class defined by irf.TimeTable.
%	See methods of irf.TimeTable class on how to generate and export time tables
%
%   IRF.TT(tt) display time table on screen
%
%   IRF.TT('list_IRF') list time tables in IRF repository
%
%   tt=IRF.TT('read_IRF',tt_id) read time table tt_id from IRF
%
%   IRF.TT(tt,'write_IRF',tt_id) write time table tt_id to IRF
%
%   [tt]=IRF.TT('read_AMDA',tt_id) read time table tt_id from AMDA
%		you can see AMDA shared time tables http://cdpp-amda.cesr.fr/DDHTML/SHARED/ttrepository.html
%
%	See also:
%		irf.TimeTable
%
%   Example:
%       TT=irf.TimeTable('demo');
%       irf.tt(TT); % display time table, same as ascii(TT)
%		TT=irf.tt('read_IRF','Example');
%		irf.tt(TT,'write_IRF','demo'); % write time table TT to brain.irfu.se with new id 'new_id'
%



%  =========== BELOW NOT IMPLEMENTED =========
%
%   hashid=IRF.TT(tt,'save_AMDA',tt_id,'user',username,'passw',password)
%     save time table tt_id. Return hashid, if saving did not succeed return NaN.
%
%   out=IRF.TT(tt,'save_AMDA',hashid) replace time table with hashid with a new one
%    return NaN if operation did not succeed.
%
%   out=IRF.TT(tt,'add_AMDA',hashid) add time intervals to the table hashid
%    return NaN if operation did not succeed.
%


%% Check inputs

if nargin==0
  help irf.tt;
  return;
end

if isa(varargin{1},'irf.TimeTable') % first argument is TT in matlab format (structure)
  tt=varargin{1};
  if nargin==1 && nargout == 0 % IRF_TT(tt) display time table on screen
    action='display';
  elseif nargin>1 && ischar(varargin{2}) % second argument is action
    action=varargin{2};
    if strcmp(action,'write_IRF') % [out]=IRF_TT(tt,'write_IRF',tt_id)
      if nargin==2
        disp('ERROR! dont understand irf_tt arguments. See help!');
        return
      else
        tt_id=varargin{3};
        if ~ischar(varargin{3}) % filename should be character array
          disp('ERROR! dont understand irf_tt arguments. See help!');
          return
        end
      end
    end
  end
elseif ischar(varargin{1}) % first argument is action
  action=varargin{1};
  if nargin>1 && ischar(varargin{2}) && strcmpi(varargin{1},'import_ascii')
    filename=varargin{2};
  end
  if nargin>1 && ischar(varargin{2}) && ...
      (strcmpi(varargin{1},'read_irf') ...
      || strcmpi(varargin{1},'read_www') ...
      || strcmpi(varargin{1},'read_amda'))
    tt_id=varargin{2};
  end
end

%% Act
switch lower(action)
  case 'display'
    ascii(tt);
  case 'read_irf'
    out=irf.TimeTable(['https://www.space.irfu.se/TT/' tt_id]);
  case 'read_amda'
    httpLink=['http://cdpp-amda.cesr.fr/DDHTML/SHARED/' tt_id '.txt'];
    out=irf.TimeTable(httpLink);
  case 'write_irf'
    remoteFile=['hq.irfu.se:/usr/home/www/space/TT/' tt_id];
    tempFile=tempname;
    export_ascii(tt,tempFile);
    eval(['!scp ' tempFile ' ' remoteFile]);
    delete(tempFile);
  case 'list_irf'
    s = urlread('https://www.space.irfu.se/TT/'); %#ok<URLRD> webread introduced in R2014b
    A = regexp(s,'a href="(?<entry>\w+)"','tokens');
    A = [A{:}]';
    disp(A)
    %		eval('!ssh hq.irfu.se ls /usr/home/www/space/TT');
end

if nargout==0
  clear out;
end

end

