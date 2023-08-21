%% Check Huishan list for similar events (nothing good found)
%workDirectory='/Users/andris/data/cluster/huishan_DFs';
%cd(workDirectory);

function batch_run(TT,varargin)

% A function to download data and make plots of MAARBLE parameters from
% all four Cluster spacecraft, given a time table, TT.
%
% TT is a time table
% varagin to be either 'pc12' or 'pc35'
% NOTE: if you put something different a
% default frequency range of .02 to 5 will be used
%

workDirectory='/Users/meghanmella/Documents/MATLAB/Maarble';
cd(workDirectory);
tt=TT;
args = varargin;
freq_range = lower(args{1});
%tt=irf_tt('read_IRF','huishans_df_list');
nEvents=numel(tt); % number of events
j=1;
%for cl_id=1:4,
for cl_id=1
  cl_id_str=int2str(cl_id);
  workDirectory2=['/Users/meghanmella/Documents/MATLAB/Maarble/C' cl_id_str];
  emic=ascii(tt);
  emic.dirName=cell(size(tt));
  for iEvent=1:nEvents
    cd(workDirectory2);
    tint=[tt.TimeInterval(iEvent) tt.TimeInterval(iEvent+numel(tt))];
    disp('=================================');
    disp(['  ' num2str(iEvent) '. event']);
    disp(irf_time(tint,'tint>utc'));
    disp('=================================');

    dirName=[sprintf('maarble_emic_',iEvent) irf_time(tint(1),'epoch>utc')];
    emic.dirName{iEvent}=dirName;
    if(~exist(dirName,'dir'))
      disp(['creating directory - ' dirName]);
      mkdir(dirName);
    end
    cd(emic.dirName{iEvent});
    dataDownloaded=sprintf('dataDownloaded');
    if(~exist(dataDownloaded,'dir'))
      %caa_download(tint,['C' cl_id_str '_CP_FGM_FULL']);
      caa_download(tint,['C' cl_id_str '_CP_FGM_FULL_ISR2']);
      caa_download(tint,['C' cl_id_str '_CP_EFW_L2_E']);
      mkdir('dataDownloaded');
    end
    %B=c_caa_var_get(['B_vec_xyz_gse__C' cl_id_str '_CP_FGM_FULL'],'mat');
    B=c_caa_var_get(['B_vec_xyz_isr2__C' cl_id_str '_CP_FGM_FULL_ISR2'],'mat');
    %xyz=c_caa_var_get(['sc_pos_xyz_gse__C' cl_id_str '_CP_FGM_FULL'],'mat');
    xyz=c_caa_var_get(['sc_pos_xyz_isr2__C' cl_id_str '_CP_FGM_FULL_ISR2'],'mat');
    E=c_caa_var_get(['E_Vec_xy_ISR2__C' cl_id_str '_CP_EFW_L2_E'],'mat');
    %figure(255+iEvent+cl_id*4), clf
    %figure(255+iEvent+cl_id*6), clf
    hf=figure(386+j), clf
    set(hf,'Position',[0 0 680 800]);

    %[timeVector,frequencyVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,...
    %    EE_xxyyzz_FAC,Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,...
    %    k_thphSVD_fac,polarization,ellipticity]...
    %    =maarble.prepare_ULF_data(E,B,xyz,1,'freq',[.01,1]);
    %    [timeVector,frequencyVector,BVector,BB_xxyyzz_fac]...
    %        =maarble.prepare_ULF_data(E,B,xyz,1,'freq',[.02,5]);
    maarble.prepare_ULF_data(E,B,xyz,cl_id,freq_range);
    %caa_download(tint,'C?_CP_FGM_5VPS');
    %caa_download(tint,'C?_CP_RAP_ESPCT6');
    j=j+1;
  end
  cd(workDirectory2);
  save _maarble_emic emic;
end
cd(workDirectory);

end
