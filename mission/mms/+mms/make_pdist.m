function PD = make_pdist(file)
% Construct PDist skymap from file name + file path
%   ePDist = mms.make_pdist(FilePath)
%
%   db_info = datastore('mms_db');  
%   c_eval('ePDist? = mms.make_pdist([db_info.local_file_db_root ''/mms?/fpi/brst/l2/des-dist/2015/10/16/mms?_fpi_brst_l2_des-dist_20151016103254_v2.1.0.cdf''])',ic)
 
filePathAndName = file;
if strfind(filePathAndName,'des')
  vS = 'des';
  species = 'electrons';
elseif strfind(filePathAndName,'dis')
  vS = 'dis';
  species = 'ions';
else
    PD = NaN;
    disp('File not recognized as particle distribution.'); 
    return;
end
fileParts = strsplit(filePathAndName,'/');
mmsId = fileParts{end}(4);

%disp('Loading electron distribution...')
tmpDataObj = dataobj(filePathAndName);
tmpDist = get_variable(tmpDataObj,['mms' mmsId '_' vS '_dist_brst']);
Dist = tmpDist.data;

% Collect User data
ud = [];
ud.GlobalAttributes = tmpDist.GlobalAttributes;
ud.CATDESC          = tmpDist.CATDESC;
if isfield(tmpDist,'DISPLAY_TYPE'), ud.DISPLAY_TYPE     = tmpDist.DISPLAY_TYPE; end
ud.FIELDNAM         = tmpDist.FIELDNAM;
ud.VALIDMIN         = tmpDist.VALIDMIN;
ud.VALIDMAX         = tmpDist.VALIDMAX;
if isfield(tmpDist,'LABLAXIS'), ud.LABLAXIS = tmpDist.LABLAXIS; end
if isfield(tmpDist,'LABL_PTR_1'), ud.LABL_PTR_1 = tmpDist.LABL_PTR_1;
elseif isfield(tmpDist,'LABL_PTR_2'), ud.LABL_PTR_2 = tmpDist.LABL_PTR_2;
elseif isfield(tmpDist,'LABL_PTR_3'), ud.LABL_PTR_3 = tmpDist.LABL_PTR_3;
end


% Shift times to center of deltat- and deltat+ for l2 particle
% distributions and moments
if ~isempty(regexp(tmpDist.name,'^mms[1-4]_d[ei]s_','once'))
	if isfield(tmpDist.DEPEND_0,'DELTA_MINUS_VAR') && isfield(tmpDist.DEPEND_0,'DELTA_PLUS_VAR'),
        if isfield(tmpDist.DEPEND_0.DELTA_MINUS_VAR,'data') && isfield(tmpDist.DEPEND_0.DELTA_PLUS_VAR,'data'),
            irf.log('critical','Times shifted to center of dt-+. dt-+ are recalculated');
            toffset = (int64(tmpDist.DEPEND_0.DELTA_PLUS_VAR.data)-int64(tmpDist.DEPEND_0.DELTA_MINUS_VAR.data))*1e6/2;
            tdiff = (int64(tmpDist.DEPEND_0.DELTA_PLUS_VAR.data)+int64(tmpDist.DEPEND_0.DELTA_MINUS_VAR.data))*1e6;
            tmpDist.DEPEND_0.DELTA_MINUS_VAR.data = tdiff/2;
            tmpDist.DEPEND_0.DELTA_PLUS_VAR.data = tdiff/2;
            tmpDist.DEPEND_0.data = tmpDist.DEPEND_0.data+toffset;
        end
    end
end

time = tmpDist.DEPEND_0.data;
if ~isa(time,'GenericTimeArray'), time = EpochTT(time); end
dt_minus = tmpDist.DEPEND_0.DELTA_MINUS_VAR.data;
dt_plus = tmpDist.DEPEND_0.DELTA_MINUS_VAR.data;
energy0 = get_variable(tmpDataObj,['mms' mmsId '_' vS '_energy0_brst']);
energy1 = get_variable(tmpDataObj,['mms' mmsId '_' vS '_energy1_brst']);
energy0 = energy0.data; energy1 = energy1.data; 
phi = get_variable(tmpDataObj,['mms' mmsId '_' vS '_phi_brst']);
theta = get_variable(tmpDataObj,['mms' mmsId '_' vS '_theta_brst']);
stepTable = get_variable(tmpDataObj,['mms' mmsId '_' vS '_steptable_parity_brst']);
phi = phi.data; theta = theta.data; stepTable = stepTable.data;

% Make energytable from energy0, energy1 and energysteptable
energy = repmat(torow(energy0),numel(stepTable),1);
energy(stepTable==1,:) = repmat(energy1,sum(stepTable),1);

% Construct PDist
PD = PDist(time,Dist,'skymap',energy,phi,theta);
PD.userData = ud; 
PD.name = tmpDist.name; 
%PD.units = tmpDist.units;
PD.units = 's^3/cm^6'; 
PD.siConversion = '1e12';
PD.species = species;
PD.ancillary.dt_minus = dt_minus; 
PD.ancillary.dt_plus = dt_plus; 
PD.ancillary.energy0 = energy0; 
PD.ancillary.energy1 = energy1;
PD.ancillary.energyStepTable = stepTable;

end