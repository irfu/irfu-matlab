function res = db_index(flagIndexON,flagSave)
%MMS.DB_INDEX  get/set DB caching settings
%
% res = MMS.DB_INDEX([flagIndex],[flagSave])
%
% Examples:
%     mms.db_index()     - display caching status
%     mms.db_index('on') - turn caching ON
%     mms.db_index('off') - turn caching OFF
%     mms.db_index('on','save') - turn caching ON permanently
%
% See also: MMS.DB_INDEX_DISP


narginchk(0,2)

global MMS_DB;
if isempty(MMS_DB) ||  isempty(MMS_DB.databases)
	mms.db_init();
	if  isempty(MMS_DB.databases)
		strTxt = 'No MMS database initialized. See help mms.db_init.';
		irf.log('critical',strTxt); error(strTxt);
	end
end

if nargin == 0 % display status
	res = MMS_DB.index.enabled;
	return;
end

if nargin>0
  switch lower(flagIndexON)
    case 'on', val = true; 
    case 'off', val = false; 
    otherwise
      warning('unrecognized value for flagIndex')
      val = [];
  end
  if ~isempty(val)
    MMS_DB.index.enabled = val;
    if nargin<2, flagSave = false;
    else
      if strcmpi(flagSave,'save'), flagSave = true;
      else
        warning('unrecognized value for flagSave')
        flagSave = false;
      end
    end
    if flagSave
      irf.log('notice','saving indexing status')
      datastore('mms_db','db_index_enabled',val),
    end
	end
end

if MMS_DB.index.enabled, disp('use of the DB indexing enabled')
else, disp('use of DB indexing disabled')
end

if nargout>0, res = MMS_DB.index.enabled; end
