function dphi=caa_comp_sun_angle(cl_id)
dphi=0;
efwa = getData(ClusterProc,cl_id,'efwa','nosave');
if(length(efwa)<1), return, end
efwa = efwa{2};
atwo = c_load('Atwo?',cl_id,'var');
if(length(atwo)<10), return, end
atwo = c_phase(efwa(:,1),atwo);
%irf_plot({efwa,atwo},'comp')
if(length(atwo)<10), return, end
[ii1,ii2]=irf_find_comm_idx(efwa,atwo);
dphi=efwa(ii1,2)-atwo(ii2,2);
dphi=mod(dphi+900,360)-180;
dphi=[efwa(ii1,1) dphi];


 

