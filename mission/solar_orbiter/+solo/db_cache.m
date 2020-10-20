function res = db_cache(flagCacheON,flagSave)
%SOLO.DB_CACHE  get/set DB caching settings
%
% res = SOLO.DB_CACHE([flagCache],[flagSave])
%
% Examples:
%     solo.db_cache()     - display caching status
%     solo.db_cache('on') - turn caching ON
%     solo.db_cache('off') - turn caching OFF
%     solo.db_cache('on','save') - turn caching ON permanently
%
% See also: SOLO.DB_CACHE_DISP


narginchk(0,2)

global SOLO_DB; if isempty(SOLO_DB), solo.db_init(), end

if nargin>0
  switch lower(flagCacheON)
    case 'on', val = true;
    case 'off', val = false;
    otherwise
      warning('unrecognized value for flagCache')
      val = [];
  end
  if ~isempty(val)
    SOLO_DB.cache.enabled = val;
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
      datastore('solo_db','db_cache_enabled',val),
    end
  end
end

if SOLO_DB.cache.enabled, disp('DB caching enabled')
else, disp('DB caching disabled')
end

if nargout>0, res = SOLO_DB.cache.enabled; end
