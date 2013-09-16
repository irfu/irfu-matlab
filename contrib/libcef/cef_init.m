%%
%	Initialisation of CEFLIB library
%	--------------------------------
%
if not (libisloaded ('libcef'))
	loadlibrary (['libcef_' lower(computer)], @libcef_mfile,'alias','libcef')
end


cef_read	= @(file)	calllib ('libcef', 'cef_read', file);

cef_close	= @()		calllib	('libcef', 'cef_close');

cef_verbosity	= @(level)	calllib ('libcef', 'cef_verbosity', level);

cef_metanames	= @()		calllib ('libcef', 'cef_metanames');

cef_meta	= @(meta)	calllib ('libcef', 'cef_meta', meta);

cef_gattributes	= @()		calllib ('libcef', 'cef_gattributes');

cef_vattributes	= @(var)	calllib ('libcef', 'cef_vattributes', var);

cef_gattr	= @(key)	calllib ('libcef', 'cef_gattr', key);

cef_vattr	= @(var,key)	calllib ('libcef', 'cef_vattr', var, key);

cef_varnames	= @()		calllib ('libcef', 'cef_varnames');

cef_var		= @(var) 	calllib ('libcef', 'cef_var', var);

cef_depends	= @(var)	calllib ('libcef', 'cef_depends', var);

cef_date 	= @(d)		d / 86400000.0 + 715146.0;

milli_to_isotime = @(var,n)	calllib ('libcef', 'milli_to_isotime', var, n);
