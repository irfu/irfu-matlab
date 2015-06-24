function ts = variable2ts(v)
%MMS.VARIABLE2TS  Converst variable structure to TSeries
%
% ts = variable2ts(v)

if ~isstruct(v) || ~isfield(v,'data') || ~isfield(v,'dim') || ...
    ~isfield(v,'nrec') || ~isfield(v,'name')|| ~isfield(v,'DEPEND_0')
  errS = 'Expecting VAR structure as returned by mms.db_get_variable()';
  irf.log('critical', errS), error(errS)
end

% User data
ud = [];
ud.GlobalAttributes = v.GlobalAttributes;
ud.CATDESC          = v.CATDESC;
ud.DISPLAY_TYPE     = v.DISPLAY_TYPE;
ud.FIELDNAM         = v.FIELDNAM;
if isfield(v,'SI_CONVERSION')
ud.SI_CONVERSION    = v.SI_CONVERSION;
end
ud.VALIDMIN         = v.VALIDMIN;
ud.VALIDMAX         = v.VALIDMAX;
if isfield(v,'LABLAXIS'), ud.LABLAXIS = v.LABLAXIS; end
if isfield(v,'LABL_PTR_1'), ud.LABL_PTR_1 = v.LABL_PTR_1;
elseif isfield(v,'LABL_PTR_2'), ud.LABL_PTR_2 = v.LABL_PTR_2;
elseif isfield(v,'LABL_PTR_3'), ud.LABL_PTR_3 = v.LABL_PTR_3;
end

varType = ''; data = v.data; 
if v.dim(1)==3 && v.dim(2)==1, varType = 'vec_xyz';
elseif v.dim(1)==2 && v.dim(2)==1, varType = 'vec_xy';
elseif v.dim(1)==1 && v.dim(2)==1, varType = 'scalar';
else
  % Special quirks for different instruments
  if regexp(v.name,'^mms[1-4]_[d,a]fg_srvy(_gsm)?_dmpa$') %AFG/DFGb1
    data = data(:,1:3); % strip Btot
    varType = 'vec_xyz';
    ud.LABL_PTR_1.data = ud.LABL_PTR_1.data(1:3,:);
    ud.LABL_PTR_1.dim(1) = 3;
    ud.VALIDMIN(4) = []; ud.VALIDMAX(4) = [];
    ud.SI_CONVERSION = '1.0e-9>T';
  end
end

if isempty(varType)
  errS = 'Unrecognized VAR: cannot converto to TSeries';
  irf.log('critical', errS), error(errS)
end

if isfield(v,'FILLVAL'), data(data==v.FILLVAL) = NaN; end
ts = feval(['irf.ts_' varType],v.DEPEND_0.data,data);
ts.name = v.name;
ts.units = v.UNITS;
ts.userData = ud;
