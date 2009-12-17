function caa_comp_sun_angle(cl_id)

efwa = getData(ClusterProc,cl_id,'efwa','nosave');
efwa = efwa{2};
atwo = c_load('Atwo?',cl_id,'var');
atwo = c_phase(efwa(:,1),atwo);
irf_plot({efwa,atwo},'comp','LineStyle',{'.','r.'})
