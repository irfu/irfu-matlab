function res = db_index(flagIndexON,flagSave)
%SOLO.DB_INDEX  get/set DB caching settings
%
% res = SOLO.DB_INDEX([flagIndex],[flagSave])
%
% Examples:
%     solo.db_index()     - display caching status
%     solo.db_index('on') - turn caching ON
%     solo.db_index('off') - turn caching OFF
%     solo.db_index('on','save') - turn caching ON permanently
%
% See also: SOLO.DB_INDEX_DISP


narginchk(0,2)

global SOLO_DB;
if isempty(SOLO_DB) ||  isempty(SOLO_DB.databases)
  solo.db_init();
  if  isempty(SOLO_DB.databases)
    strTxt = 'No SOLO database initialized. See help solo.db_init.';
    irf.log('critical',strTxt); error(strTxt);
  end
end

if nargin == 0 % display status
  res = SOLO_DB.index.enabled;
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
    SOLO_DB.index.enabled = val;
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
      datastore('solo_db','db_index_enabled',val),
    end
  end
end

if SOLO_DB.index.enabled, disp('use of the DB indexing enabled')
else, disp('use of DB indexing disabled')
end

if nargout>0, res = SOLO_DB.index.enabled; end
