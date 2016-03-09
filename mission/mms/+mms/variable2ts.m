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
if isfield(v,'DISPLAY_TYPE'), ud.DISPLAY_TYPE     = v.DISPLAY_TYPE; end
ud.FIELDNAM         = v.FIELDNAM;
ud.VALIDMIN         = v.VALIDMIN;
ud.VALIDMAX         = v.VALIDMAX;
if isfield(v,'LABLAXIS'), ud.LABLAXIS = v.LABLAXIS; end
if isfield(v,'LABL_PTR_1'), ud.LABL_PTR_1 = v.LABL_PTR_1;
elseif isfield(v,'LABL_PTR_2'), ud.LABL_PTR_2 = v.LABL_PTR_2;
elseif isfield(v,'LABL_PTR_3'), ud.LABL_PTR_3 = v.LABL_PTR_3;
end

data = v.data; siConversion = '';
if v.dim(1)==3 && v.dim(2)==1, varType = 'vec_xyz';
elseif v.dim(1)==2 && v.dim(2)==1, varType = 'vec_xy';
elseif v.dim(1)==1 && v.dim(2)==1, varType = 'scalar';
else
  % Special quirks for different instruments
  if ~isempty(regexp(v.name,'^mms[1-4]_[d,a]fg_(srvy|brst)(_gsm)?_dmpa$', 'once'))|| ...
      ~isempty(regexp(v.name,'^mms[1-4]_[d,a]fg_(srvy|brst)_l2pre_(gse|gsm|dmpa|bcs)$', 'once')) || ...
      ~isempty(regexp(v.name,'^mms[1-4]_fgm_b_(gse|gsm|dmpa|bcs)_(srvy|brst)_l2$', 'once'))%AFG/DFGb1
    data = data(:,1:3); % strip Btot
    varType = 'vec_xyz';
    ud.LABL_PTR_1.data = ud.LABL_PTR_1.data(1:3,:);
    ud.LABL_PTR_1.dim(1) = 3;
    ud.VALIDMIN(4) = []; ud.VALIDMAX(4) = [];
    siConversion = '1.0e-9>T';
  elseif ~isempty(regexp(v.name,'^mms[1-4]_?(ql_)pos_gs(e|m)$', 'once'))
    data = data(:,1:3); % strip Rtot
    varType = 'vec_xyz';
    ud.LABL_PTR_1.data = ud.LABL_PTR_1.data(1:3,:);
    ud.LABL_PTR_1.dim(1) = 3;
    siConversion = '1.0e3>m';
  else varType = 'scalar';
  end
end

% Shift times to center of deltat- and deltat+ for l2 particle
% distributions and moments
if ~isempty(regexp(v.name,'^mms[1-4]_d[ei]s_','once'))
	if isfield(v.DEPEND_0,'DELTA_MINUS_VAR') && isfield(v.DEPEND_0,'DELTA_PLUS_VAR'),
        if isfield(v.DEPEND_0.DELTA_MINUS_VAR,'data') && isfield(v.DEPEND_0.DELTA_PLUS_VAR,'data'),
            irf.log('warning','Times shifted to center of dt-+. dt-+ are recalculated');
            toffset = (int64(v.DEPEND_0.DELTA_PLUS_VAR.data)-int64(v.DEPEND_0.DELTA_MINUS_VAR.data))*1e6/2;
            tdiff = (int64(v.DEPEND_0.DELTA_PLUS_VAR.data)+int64(v.DEPEND_0.DELTA_MINUS_VAR.data))*1e6/2;
            v.DEPEND_0.DELTA_MINUS_VAR.data = tdiff;
            v.DEPEND_0.DELTA_PLUS_VAR.data = tdiff;
            v.DEPEND_0.data = v.DEPEND_0.data+toffset;
        end
    end
end

if isempty(varType)
  errS = 'Unrecognized VAR: cannot convert to TSeries';
  irf.log('critical', errS), error(errS)
end

if isfield(v,'FILLVAL'), data(data==v.FILLVAL) = NaN; end
ts = feval(['irf.ts_' varType],v.DEPEND_0.data,data);
ts.name = v.name;
if isfield(v,'UNITS'), ts.units = v.UNITS;
else ts.units = 'unitless';
end
if ~isempty(siConversion), ts.siConversion = siConversion;
elseif isfield(v,'SI_CONVERSION'), ts.siConversion    = v.SI_CONVERSION;
end
ts.userData = ud;
