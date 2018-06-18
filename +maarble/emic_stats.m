if ismac
cd ('/Users/yuri/Dropbox (IRFU)/projects-all/MAARBLE/EMIC_stats/tt')
else, cd ('/home/yuri/EMIC_stats/tt')
end

dbFileName = ['dbEMIC_' irf_fname(irf_time(now,'date>epoch'))];
bands = {'H','He','O'};
for iBand = 1:length(bands)
  band = bands{iBand};
  % THEMIS
  th = 'ABCDE';
  for th_id=1:5
    th_s = th(th_id);
    ttFile = sprintf('CC_TH%s_MAARBLE_%s_band_wave_events',th_s,band);
    tt = irf.TimeTable(ttFile);
    Out = maarble.wave_psd_ell_thetak(tt);
    dbName = sprintf('dbEMIC_%s_TH%s',band,th_s);
    eval([dbName '=Out;'])
    clear Out
    if exist([dbFileName '.mat'],'file'), save(dbFileName,dbName,'-append')
    else, save(dbFileName,dbName)
    end
  end
  % Cluster
  for cl_id=1:4
    ttFile = sprintf('C%d_MAARBLE_%s_band_wave_events',cl_id,band);
    tt = irf.TimeTable(ttFile);
    Out = maarble.wave_psd_ell_thetak(tt);
    dbName = sprintf('dbEMIC_%s_C%d',band,cl_id);
    eval([dbName '=Out;'])
    clear Out
    if exist([dbFileName '.mat'],'file'), save(dbFileName,dbName,'-append')
    else, save(dbFileName,dbName)
    end
  end
end