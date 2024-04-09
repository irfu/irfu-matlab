function ts = variable2ts(v)
%SOLO.VARIABLE2TS  Converts variable structure to TSeries
%
% ts = variable2ts(v)

if ~isstruct(v) || ~isfield(v,'data') || ~isfield(v,'dim') || ...
    ~isfield(v,'nrec') || ~isfield(v,'name')|| ~isfield(v,'DEPEND_0')
  errS = 'Expecting VAR structure as returned by solo.db_get_variable()';
  irf.log('critical', errS), error(errS)
end

% User data
ud = [];
ud.GlobalAttributes = v.GlobalAttributes;
ud.CATDESC          = v.CATDESC;
if isfield(v,'DISPLAY_TYPE'), ud.DISPLAY_TYPE     = v.DISPLAY_TYPE; end
ud.FIELDNAM         = v.FIELDNAM;
if isfield(v,'VALIDMIN'), ud.VALIDMIN = v.VALIDMIN; end
if isfield(v,'VALIDMAX'), ud.VALIDMAX = v.VALIDMAX; end
if isfield(v,'LABLAXIS'), ud.LABLAXIS = v.LABLAXIS; end
if isfield(v,'LABL_PTR_1')
  ud.LABL_PTR_1 = v.LABL_PTR_1;
elseif isfield(v,'LABL_PTR_2')
  ud.LABL_PTR_2 = v.LABL_PTR_2;
elseif isfield(v,'LABL_PTR_3')
  ud.LABL_PTR_3 = v.LABL_PTR_3;
end

data = v.data; siConversion = '';
if     v.dim(1) == 3 && v.dim(2) == 1, varType = 'vec_xyz';
elseif v.dim(1) == 3 && v.dim(2) == 3 && isfield(v,'TENSOR_ORDER') && v.TENSOR_ORDER == 2, varType = 'tensor_xyz';
elseif v.dim(1) == 2 && v.dim(2) == 1, varType = 'vec_xy';
elseif v.dim(1) == 1 && v.dim(2) == 1, varType = 'scalar';
else
  varType = 'scalar';
end

if isempty(varType)
  errS = 'Unrecognized VAR: cannot convert to TSeries';
  irf.log('critical', errS), error(errS)
end

if isfield(v,'FILLVAL'), data(data==v.FILLVAL) = NaN; end
ts = feval(['irf.ts_' varType],v.DEPEND_0.data,data);
ts.name = v.name;
if isfield(v,'UNITS')
  ts.units = v.UNITS;
else
  ts.units = 'unitless';
end
if ~isempty(siConversion)
  ts.siConversion = siConversion;
elseif isfield(v,'SI_CONVERSION')
  ts.siConversion = v.SI_CONVERSION;
end
ts.userData = ud;
end
