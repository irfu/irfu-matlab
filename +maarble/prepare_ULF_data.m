function [timeVector,frequencyVector,BVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,...
    Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,k_thphSVD_fac,polSVD_fac,ellipticity]=prepare_ULF_data(e,b,xyz,varargin)
% This routine prepares the ULF data to be analyzed for the Maarble
% project. The routine then calculates the various wave parameters using
% the irf_ebsp routine, and then plots them using irf_pl_ebsp.


% e in ISR2; b in gse
% get e from E_Vec_xy_ISR2__C?_CP_EFW_L2_E
% get b from B_vec_xyz_gse__C?_CP_FGM_FULL
% position vector xyz from sc_pos_xyz_gse__C?_CP_FGM_FULL 
% converts b and xyz to ISR2
% first entry in varargin to be cl_id
% second and third entries to be the frequency range in the form: 'freq',[.01,1]

% Example:
%   maarble.prepare_ULF_data(E,B,xyz,4,'freq',[.01,1]);


  %% Check input 
  [ax,args,nargs] = axescheck(varargin{:});
  %x=args{1};
  if isempty(args), % nothing to plot, first input parameter empty
    display('first entry of varargin must contain cluster ID');
    return;
  end
  
  cl_id=varargin{1};
  args=args(2:end);

  %% Convert b and xyz to ISR2
  b = c_coord_trans('GSE','ISR2',b,'cl_id',cl_id);
  xyz = c_coord_trans('GSE','ISR2',xyz,'cl_id',cl_id);

  original_args=args;
  var_desc{1} = '';
  flag_subplot = 0;
  have_options = 1;

  if nargs > 1, have_options = 1; end
  
  % TODO: this needs to be determined from data ?
  %sampl=22.5;
  
  %% Default values that can be overridden by options
% dt = 0;
% flag_yy = 0;
% scaleyy = 1;
% plot_type = '';
% marker = '-';
% flag_plot_all_data=1;
flag_colorbar=1;
colorbar_scale=1;

   %% Set some important parameters
%   freq_int=[.01 5];
%   freq_number=25;
%   Morlet_width=5.36;
% 
%   if isnumeric(args{end})
%       freq_int=args{end};
%       args=args(1:end-2);
%   elseif ismember({'freq'},args)
%       disp('frequency interval values missing. using default')
%       args=args(1:end-1);
%   end
%   
%   amin=log10(0.5*sampl/freq_int(2));amax=log10(0.5*sampl/freq_int(1));anumber=freq_number;
%   a=logspace(amin,amax,anumber);
%   w0=sampl/2; % The maximum frequency
%   newfreq=w0./a;

  %% get background magnetic field
  sampl_low=5;
  t_low=b(1,1):1/sampl_low:b(end,1); t_low=t_low';
  b_low=irf_resamp(b,t_low); %low sample frequency to avoid filter issues
  bf=irf_filt(b_low,1/600,0,[],5);
%  bf=irf_filt(b,1/300,0,[],5);
  b0=b_low;
  b0(:,2:4)=b_low(:,2:4)-bf(:,2:4);
  B=b0;
%  b=bf;  %use B with background field removed
% %% resample to 25 Hz
  sampl_b=1/(b(2,1)-b(1,1));
  sampl=25;
  t=b(1,1):1/sampl:b(end,1); t=t'; 
  B=irf_resamp(B,t);
  e=irf_resamp(e,t); b=irf_resamp(b,t); disp('resampling to 25 Hz');
  b(:,2:4)=b(:,2:4)-B(:,2:4);
  
  %% get third component of e
  
  %e(:,4)=-(b(:,2).*e(:,2)+b(:,3).*e(:,3))./b(:,4);
  for j=1,size(b(:,1)),
%      if b(j,4) > .00375
      %if b(j,4) > mean(b(:,2).*b(:,2)+b(:,3).*b(:,3)+b(:,4).*b(:,4))-10^-2.5
      %if b(j,4) > mean(b(:,2).*b(:,2)+b(:,3).*b(:,3)+b(:,4).*b(:,4))-10^1
      if b(j,4) > 0.5*(b(j,2).*b(j,2)+b(j,3).*b(j,3))^.5
          e(j,4)=-(b(j,2).*e(j,2)+b(j,3).*e(j,3))./b(j,4);
      else
          e(j,4)=NaN;
          %e(j,4)=0;
      end 
  end

  %% Calculations
  if nargout==4,
    [timeVector,frequencyVector,BVector,BB_xxyyzz_fac]=...
        irf_ebsp(e,b,B,xyz,varargin(1),varargin{2});
    h=irf_pl_ebsp(timeVector,frequencyVector,BVector,BB_xxyyzz_fac);
  elseif nargout==8,
    [timeVector,frequencyVector,BVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,...
        Poynting_xyz_FAC]=irf_ebsp(e,b,B,xyz,varargin(1),varargin{2});
    h=irf_pl_ebsp(timeVector,frequencyVector,BVector,BB_xxyyzz_fac,...
        EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,Poynting_xyz_FAC,Poynting_rThetaPhi_FAC);
  else
    [timeVector,frequencyVector,BVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,...
        Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,k_thphSVD_fac,polSVD_fac,ellipticity]=...
        irf_ebsp(e,b,B,xyz,varargin(2),varargin{3});
    %plot_params={{timeVector},{frequencyVector},{BB_xxyyzz_fac},{EESum_xxyyzz_ISR2},{Poynting_xyz_FAC},{polarization},{ellipticity}};
    %plot_params={timeVector,frequencyVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,Poynting_xyz_FAC,polarization,ellipticity};
    h=irf_pl_ebsp(timeVector,frequencyVector,BVector,BB_xxyyzz_fac,...
        EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,k_thphSVD_fac,polSVD_fac,ellipticity);
  end
%   plot_params=varargout{1};
%   for i=2,nargout,
%       plot_params=strcat(plot_params,varargout{i});
%   end
  
%  included_params=ismember({'e' 'b' 's' 'eb' 'pol' 'ellip'},lower(args));
%   if any(included_params(1:4) == 1) && all(included_params(5:6) == 0),
%     [e1,b1,powerEx_plot,powerEy_plot,powerEz_plot,power2E_plot,powerBx_plot,powerBy_plot,powerBz_plot,power2B_plot,Spar_plot,EtoB_plot,powerBx_SM_plot,powerBy_SM_plot,powerBz_SM_plot,power2B_SM_plot]=irf_ebs(e1,b1,xyz,freq_int);
% %    [e1,power2E_plot,Spar_plot,power2B_SM_plot]=irf_ebs(e1,b1,xyz);
%     disp('irf_ebs ... calculations');
%   elseif all(included_params(1:4) == 0) && any(included_params(5:6) == 1),
%     [e1,b1,Lp,Ls5]=irf_pol(e1,b1,xyz,freq_int);
% %    [e1,Lp,Ls5]=irf_pol(e1,b1,xyz);
%     disp('irf_pol ... calculations');
%   elseif all(included_params([1 3:4]) == 0) && any(included_params([2 5:6]) == 1),
%     [e1,b1,Lp,Ls5,powerBx_SM_plot,powerBy_SM_plot,powerBz_SM_plot,power2B_SM_plot]=irf_bpol(e1,b1,xyz,freq_int);
% %    [e1,Lp,Ls5,power2B_SM_plot]=irf_bpol(e1,b1,xyz);
%     disp('irf_bpol ... calculations');
%   else    
%     %for everything
%     [e1,b1,Lp,Ls5,powerEx_plot,powerEy_plot,powerEz_plot,power2E_plot,powerBx_plot,powerBy_plot,powerBz_plot,power2B_plot,Spar_plot,EtoB_plot,powerBx_SM_plot,powerBy_SM_plot,powerBz_SM_plot,power2B_SM_plot]=irf_ebsp(e1,b1,xyz,freq_int);
% %    [e1,Lp,Ls5,power2E_plot,Spar_plot,power2B_SM_plot]=irf_ebsp(e1,b1,xyz);
%     disp('irf_ebsp ... calculations');

%% Plot
  %[timeVector,frequencyVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,Poynting_xyz_FAC,polarization,ellipticity]=irf_pl_ebsp(e,newfreq,Morlet_width,plot_params);
  %h=irf_pl_ebsp(plot_params);
  
  end
  
  
