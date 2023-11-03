function caa_apply_manproblems(st,dt,cli,extraprocessing,redo_sp_dir)
%CAA_APPLY_MANPROBLEMS: apply any manproblems during specified interval
%
% caa_apply_manproblems(st,dt,[cli],[extraprocessing],[redo_sp_dir])
% Input:
%  st,dt: start and duration (or end time) for the intervals to be
%         scanned, either as ISO strings or epochs.
%  cli: cluster satellite number (default=1:4)'
%  extraprocessing: extra processing steps to pass to ClusterProc/getData
%         after manproblems (e.g. to run 'manproblems|p|ps', pass 'p|ps')
%  redo_sp_dir: If this is not empty, summary plots are re-done for the
%         affected intervals and put in the specified directory.
%
narginchk(2,5)
if nargin < 5,redo_sp_dir=[];end
if (nargin < 4) || isempty(extraprocessing),varsproc='manproblems'; else, varsproc=['manproblems|' extraprocessing]; end
if nargin < 3, cli=1:4; end
if isa(st,'char'),st=iso2epoch(st); end
if isa(dt,'char'),dt=iso2epoch(dt); end
if dt>iso2epoch('2000-01-01T00:00:00Z'), dt=dt-st; end
et=st+dt;
old_pwd=pwd;
touched=[];

for cl_id=cli
  d=[c_ctl('get', 5, 'data_path') '/caa-control'];
  f_name = [d '/manual_problems_c' num2str(cl_id) '.dat'];
  if ~exist(f_name,'file')
    irf_log('load',['file ' f_name ' not found']);
    cd(old_pwd), return
  end
  fid = fopen(f_name);
  C = textscan(fid, '%s %n %1[+-] %s','commentStyle', '%');
  fclose(fid);

  for i=1:length(C{1})
    st_mp=iso2epoch(C{1}{i});
    dt_mp=C{2}(i);
    if (st_mp<et && st_mp>st)
      dirs=caa_get_subdirs(epoch2iso(st_mp), dt_mp, cl_id);
      for j=1:length(dirs)
        cd(dirs{j});
        touched=[touched '|' dirs{j}]; %#ok<AGROW>
        c_get_batch(0,0,'sc_list',cl_id,'sp',pwd,'varsproc',varsproc,'check_caa_sh_interval',1,'nosrc')
      end
    end
  end
end

if ~isempty(redo_sp_dir) && ~isempty(touched)
  touched=tokenize(touched,'|');
  for i=1:length(touched), touched{i}=touched{i}(1:strfind(touched{i},'/C')-1); end
  touched=unique(touched);
  cd(redo_sp_dir);
  for i=1:length(touched),caa_pl_summary_l1(-1,-1,touched{i},'savepdf'); end
end

cd(old_pwd);
end