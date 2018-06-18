function out = thor_get_omni(tint,vars)  
% irf_get_data_omni can only download one years 1-min data at a time, so
% this one downloads different years separately and then combines them into
% a TSeries object.
%   
%   omni_data = thor_get_omni(tint,vars);

years = tint.start:60*60*24*366:(tint.stop+60*60*24*366);
yearsUTC = years.utc;
yearsVec = str2num(yearsUTC(:,1:4));

startYear = yearsVec(1);
subYears = yearsVec-yearsVec(1);


tsub = 1;
c_eval('tintUTC{tsub} = ''?-01-01T00:00:00/?-12-31T23:59:00''; tsub = tsub+1;',startYear+[subYears(1):subYears(end)]);

omni_orig = [];
tic;
for iy = 1:numel(tintUTC)  
  tint = irf.tint(tintUTC{iy});
  tmp_omni = irf_get_data_omni(tint,vars,'omni_min');
  omni_orig = [omni_orig; tmp_omni];
  toc
end

out = omni_orig;  
%t0 = irf_time(omni_orig(:,1),'epoch>epochtt');