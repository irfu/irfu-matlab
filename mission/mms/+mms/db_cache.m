function res = db_cache(flagCacheON,flagSave)
%MMS.DB_CACHE  get/set DB caching settings
%
% res = MMS.DB_CACHE([flagCache],[flagSave])
%
% Examples:
%     mms.db_cache()     - display caching status
%     mms.db_cache('on') - turn caching ON
%     mms.db_cache('off') - turn caching OFF
%     mms.db_cache('on','save') - turn caching ON permanently
%
% See also: MMS.DB_CACHE_DISP


narginchk(0,2)

global MMS_DB; if isempty(MMS_DB), mms.db_init(), end

if nargin>0
  switch lower(flagCacheON)
    case 'on', val = true; 
    case 'off', val = false; 
    otherwise
      warning('unrecognized value for flagCache')
      val = [];
  end
  if ~isempty(val)
    MMS_DB.cache.enabled = val;
    if nargin<2, flagSave = false;
    else
      if strcmpi(flagSave,'save'), flagSave = true;
      else
        warning('unrecognized value for flagSave')
        flagSave = false;
      end
    end
    if flagSave
      disp('saving caching status')
      datastore('mms_db','db_cache_enabled',val),
    end
  end
end

if MMS_DB.cache.enabled, disp('DB caching enabled')
else, disp('DB caching disabled')
end

if nargout>0, res = MMS_DB.cache.enabled; end
