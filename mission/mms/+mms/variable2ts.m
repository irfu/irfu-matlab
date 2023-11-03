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
if isfield(v,'VALIDMIN'), ud.VALIDMIN = v.VALIDMIN; end
if isfield(v,'VALIDMAX'), ud.VALIDMAX = v.VALIDMAX; end
if isfield(v,'LABLAXIS'), ud.LABLAXIS = v.LABLAXIS; end
if isfield(v,'LABL_PTR_1'), ud.LABL_PTR_1 = v.LABL_PTR_1;
elseif isfield(v,'LABL_PTR_2'), ud.LABL_PTR_2 = v.LABL_PTR_2;
elseif isfield(v,'LABL_PTR_3'), ud.LABL_PTR_3 = v.LABL_PTR_3;
end

data = v.data; siConversion = '';
if v.dim(1)==3 && v.dim(2)==1, varType = 'vec_xyz';
elseif v.dim(1) == 3 && v.dim(2) == 3 && isfield(v,'TENSOR_ORDER') && v.TENSOR_ORDER == 2, varType = 'tensor_xyz';
elseif v.dim(1)==2 && v.dim(2)==1, varType = 'vec_xy';
elseif v.dim(1)==1 && v.dim(2)==1, varType = 'scalar';
else
  % Special quirks for different instruments
  if (~isempty(regexp(v.name,'^mms[1-4]_[d,a]fg_(srvy|brst)(_gsm)?_dmpa$', 'once'))|| ...
      ~isempty(regexp(v.name,'^mms[1-4]_[d,a]fg_(srvy|brst)_l2pre_(gse|gsm|dmpa|bcs)$', 'once')) || ...
      ~isempty(regexp(v.name,'^mms[1-4]_fgm_b_(gse|gsm|dmpa|bcs)_(srvy|brst)_l2$', 'once')) || ...
      ~isempty(regexp(v.name,'^mms[1-4]_[d,a]fg_b_(gse|gsm|dmpa|bcs)_(srvy|brst)_l2pre$', 'once')) )%AFG/DFGb1
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
  else, varType = 'scalar';
  end
end

% Shift times to center of deltat- and deltat+ for l2 particle
% distributions and moments
if ~isempty(regexp(v.name,'^mms[1-4]_d[ei]s_','once')) || ~isempty(regexp(v.name,'^mms[1-4]_hpca_','once'))
  if isfield(v.DEPEND_0,'DELTA_MINUS_VAR') && isfield(v.DEPEND_0,'DELTA_PLUS_VAR')
    if isfield(v.DEPEND_0.DELTA_MINUS_VAR,'data') && isfield(v.DEPEND_0.DELTA_PLUS_VAR,'data')
      irf.log('warning','Times shifted to center of dt-+. dt-+ are recalculated');
      flag_MINUS = 1e3;       flag_PLUS = 1e3;
      if isfield(v.DEPEND_0.DELTA_MINUS_VAR, 'UNITS') && isfield(v.DEPEND_0.DELTA_PLUS_VAR, 'UNITS')
        if strcmp(v.DEPEND_0.DELTA_MINUS_VAR.UNITS, 's')
          flag_MINUS = 1e3;           % s --> ms
        elseif strcmp(v.DEPEND_0.DELTA_MINUS_VAR.UNITS, 'ms')
          flag_MINUS = 1;
        else
          irf.log('warning','Epoch_minus_var units are not clear, assume s');
          flag_MINUS = 1e3;
        end
        if strcmp(v.DEPEND_0.DELTA_PLUS_VAR.UNITS, 's')
          flag_PLUS = 1e3;           % s --> ms
        elseif strcmp(v.DEPEND_0.DELTA_PLUS_VAR.UNITS, 'ms')
          flag_PLUS = 1;
        else
          irf.log('warning','Epoch_plus_var units are not clear, assume s');
          flag_PLUS = 1e3;
        end
      else
        irf.log('warning','Epoch_plus_var/Epoch_minus_var units are not clear, assume s');
      end
      % toffset = (int64(v.DEPEND_0.DELTA_PLUS_VAR.data)-int64(v.DEPEND_0.DELTA_MINUS_VAR.data))*1e6/2;
      % tdiff = (int64(v.DEPEND_0.DELTA_PLUS_VAR.data)+int64(v.DEPEND_0.DELTA_MINUS_VAR.data))*1e6/2;
      % v.DEPEND_0.DELTA_PLUS_VAR.data*flag_PLUS in ms
      toffset = int64((v.DEPEND_0.DELTA_PLUS_VAR.data*flag_PLUS-v.DEPEND_0.DELTA_MINUS_VAR.data*flag_MINUS)*1e6/2);
      tdiff = int64((v.DEPEND_0.DELTA_PLUS_VAR.data*flag_PLUS+v.DEPEND_0.DELTA_MINUS_VAR.data*flag_MINUS)*1e6/2);
      tdiff_data = median(diff(v.DEPEND_0.data)) / 2;                   % ns
      if ~(tdiff_data == mean(tdiff))
        str1 = num2str(mean(v.DEPEND_0.DELTA_PLUS_VAR.data)*flag_PLUS);
        str2 = num2str(mean(v.DEPEND_0.DELTA_MINUS_VAR.data)*flag_MINUS);
        str0 = num2str(double(tdiff_data * 2) * 1e-6);
        str_warning = ['Epoch_plus_var (' str1 ') and Epoch_minus_var (' str2 ') donot match data sampling time (' str0 '), assume tdiff_data; '];
        irf.log('warning', str_warning);
        tdiff = tdiff_data;
        toffset = tdiff_data;
      end
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
else, ts.units = 'unitless';
end
if ~isempty(siConversion), ts.siConversion = siConversion;
elseif isfield(v,'SI_CONVERSION'), ts.siConversion    = v.SI_CONVERSION;
end
ts.userData = ud;
