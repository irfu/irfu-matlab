cd ('/Users/yuri/Dropbox (IRFU)/projects-all/MAARBLE/EMIC_stats/tt')
band = 'H';
for cl_id=1:4
  ttFile = sprintf('C%d_MAARBLE_%s_band_wave_events',cl_id,band);
  tt = irf.TimeTable(ttFile);
  Out = maarble.wave_psd_ell_thetak(tt); %#ok<NASGU>
  dbName = sprintf('dbEMIC_%s_C%d',band,cl_id);
  eval([dbName '=Out;'])
  clear Out
end