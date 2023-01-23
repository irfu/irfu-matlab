function [t_start_save,t_stop_save,fdm_save]=caa_efw_mode_tab(fdm,var_list_s)
%CAA_EFW_MODE_TAB  read EFW FDM and extract time intervals of constant FDM
%
% [t_start_save,t_stop_save,fdm_save] = caa_efw_mode_tab(fdm,var_list_s)
%
% fdm - FRM{cl_id} produced by ClusterDB/getData
% var_list_s is in form 'px|r|ss'
% where:
% px	Burst playback
% bbb	Burst internal state
% r		Whisper pulses present
% w		Sweep in progress
% eeee	E/D mode
% ss	Sampling mode
% tm	Tape mode
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

t_fdm = fdm(:,1);
fdm_is = fdm(:,2:end)';
var_list = tokenize(var_list_s,'|');
t_start_save = [];
t_stop_save = [];
fdm_save = {};

if ~isempty(t_fdm)
  % this info is needed if next interval has no data
  prev_ok = 1;
  t_end_last_interv = t_fdm(end);
  
  % convert to bits
  s_pbbbxrwc = dec2bin(fdm_is(1,:),8);
  s_ssqiiiii = dec2bin(fdm_is(2,:),8);
  s_mmmmeeee = dec2bin(fdm_is(4,:),8);
  fdm_tm=fdm_is(5,:)';
  
  clear fdm_is
  
  fdm_px = str2num( s_pbbbxrwc(:,1) );  %#ok<ST2NM>
  fdm_bbb = str2num( s_pbbbxrwc(:,2) );  %#ok<ST2NM>
  fdm_bbb(:,2) = str2num( s_pbbbxrwc(:,3) );  %#ok<ST2NM>
  fdm_bbb(:,3) = str2num( s_pbbbxrwc(:,4) );  %#ok<ST2NM>
  fdm_px(:,2) = str2num( s_pbbbxrwc(:,5) ); %#ok<ST2NM>
  fdm_r = str2num( s_pbbbxrwc(:,6) ); %#ok<ST2NM>
  fdm_w = str2num( s_pbbbxrwc(:,7) ); %#ok<ST2NM>
  fdm_s = str2num( s_ssqiiiii(:,1) ); %#ok<ST2NM>
  fdm_s(:,2) = str2num( s_ssqiiiii(:,2) ); %#ok<ST2NM>
  fdm_e = str2num( s_mmmmeeee(:,5) ); %#ok<ST2NM>
  fdm_e(:,2) = str2num( s_mmmmeeee(:,6) ); %#ok<ST2NM>
  fdm_e(:,3) = str2num( s_mmmmeeee(:,7) ); %#ok<ST2NM>
  fdm_e(:,4) = str2num( s_mmmmeeee(:,8) ); %#ok<ST2NM>
  
  clear s_pbbbxrwc s_ssqiiiii s_mmmmeeee
  
  t_fdm_last = t_fdm(1);
  for j=1:length(var_list)
    if strcmp(var_list(j),'px')
      fdm_last_px = fdm_px(1,:); %#ok<NASGU>
    elseif strcmp(var_list(j),'bbb')
      fdm_last_bbb = fdm_bbb(1,:); %#ok<NASGU>
    elseif strcmp(var_list(j),'r')
      fdm_last_r = fdm_r(1,:); %#ok<NASGU>
    elseif strcmp(var_list(j),'w')
      fdm_last_w = fdm_w(1,:); %#ok<NASGU>
    elseif strcmp(var_list(j),'s')
      fdm_last_s = fdm_s(1,:); %#ok<NASGU>
    elseif strcmp(var_list(j),'e')
      fdm_last_e = fdm_e(1,:); %#ok<NASGU>
    elseif strcmp(var_list(j),'tm')
      fdm_last_tm = fdm_tm(1); %#ok<NASGU>
    else
      error('unknown variable')
    end
  end
  
  % If we have more than one point
  if length(fdm(:,1))>1
    empty = 1;
    for j=1:length(var_list)
      eval(['ii{j} = irf_find_diff(fdm_' var_list{j} ',fdm_last_' var_list{j} ');']);
      if empty && ~isempty(ii{j}), empty = 0; end  %#ok<USENS>
    end
    
    % no changes (should be a typical situation)
    if ~empty
      ii_all = [];
      for j=1:length(var_list)
        if ~isempty(ii{j})
          n = length(ii{j});
          ii_all(end+1:end+n) = ii{j};
        end
      end
      ii_all = irf_rm_double_idx(ii_all);
      n_diff = length(ii_all);
      
      for jj = 1:n_diff
        l = length(t_start_save);
        
        % save a good interval
        t_start_save(l+1) = t_fdm_last;
        % if the change happened exactly at the boundary between intervals
        if ii_all(jj)-1 >= 1, t_stop_save(l+1) = t_fdm(ii_all(jj)-1);
        else, t_stop_save(l+1) = t_end_last_interv;
        end
        
        for j=1:length(var_list)
          eval(['if length(fdm_last_' var_list{j} ')>1,fdm_save_' var_list{j} '(l+1,:)=fdm_last_' var_list{j} '; else, fdm_save_' var_list{j} '(l+1)=fdm_last_' var_list{j} '; end ']);
        end
        
        t_fdm_last = t_fdm(ii_all(jj));
        for j=1:length(var_list)
          eval(['fdm_last_' var_list{j} '= fdm_' var_list{j} '(ii_all(jj),:);']);
        end
        
      end
    end
  end
end


if prev_ok %last interval was also good
  if (t_fdm_last <= t_fdm(end))
    l = length(t_start_save);
    
    % save the good interval
    t_start_save(l+1) = t_fdm_last;
    t_stop_save(l+1) = t_fdm(end);
    
    %fdm_save{l+1} = fdm_last;
    for j=1:length(var_list)
      eval(['if length(fdm_last_' var_list{j} ')>1,fdm_save_' var_list{j} '(l+1,:)=fdm_last_' var_list{j} '; else, fdm_save_' var_list{j} '(l+1)=fdm_last_' var_list{j} '; end ']);
    end
  end
end

if length(var_list)>1
  for j=1:length(var_list)
    eval(['fdm_save{j}=fdm_save_' var_list{j} ';']);
  end
else
  if exist(['fdm_save_' var_list{1}],'var')
    eval(['fdm_save=fdm_save_' var_list{1} ';']);
  else, fdm_save = [];
  end
end
