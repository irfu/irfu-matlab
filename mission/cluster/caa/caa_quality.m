function caa_quality(dflt)
%CAA_QUALITY  set quality factor for the CAA products
%
% caa_quality([q_factor])
%   q_factor 1..4. (default is 3)
%   1 - known problems, use at your own risk
%   2 - survey data, not for publication
%   3 - good for publication, subject to PI approval
%   4 - excellent data which has received special treatment
%

if nargin<1, dflt=3; end

mmm = ['h - help';...
  'q - quit';...
  'd - disp';...
  's - save'];
qfac_s  = ['1 - known problems         ';...
  '2 - survey data, not f/publ';...
  '3 - good for publication   ';...
  '4 - excellent data         '];
if exist('./mInfo.mat','file')
  load -mat mInfo caa_q
end
if ~exist('caa_q','var')
  irf_log('proc',['using default Q=' num2str(dflt)])
  caa_q.p = ones(4,3)*dflt;
  caa_q.e = ones(4,3)*dflt;
end
disp_q(caa_q)

disp('commands:')
disp(mmm)
disp(['current dir is ' pwd])
q='0';flag_save=1;
while(q ~= 'q') % ====== MAIN LOOP =========
  
  q = input('input command or data level 1..3/all(a)>','s');
  if isempty(q), q = '0'; end
  if strcmp(q,'q'), return,
  elseif q == '0' || strcmp(q,'h'), disp('commands:'), disp(mmm)
  elseif strcmp(q,'d'), disp_q(caa_q)
  elseif strcmp(q,'s')
    if exist('./mInfo.mat','file'), save mInfo caa_q -append
    else, save mInfo caa_q
    end
    irf_log('save','caa_q -> mInfo.mat')
  elseif strcmp(q,'1') || strcmp(q,'2') || strcmp(q,'3')
    disp(sprintf('modifying level %s',q))
    lev = str2num(q);
    ok = 0;
    while ~ok
      sc_list = input('cl_id 1..4, all(a)>','s');
      if strcmp(sc_list,'a') || strcmp(sc_list,'all') || isempty(sc_list)
        sc_list = 1:4;
        ok = 1;
      else
        try
          sc_list = str2num(sc_list);
          ii = find(sc_list>4 | sc_list<1);
          if isempty(ii), ok = 1;
          else,disp('wrong input')
          end
        end
      end
    end
    ok = 0;
    while ~ok
      vars = input('e/p/all(a)>','s');
      if strcmp(vars,'a') || strcmp(vars,'all'), vars = [1 2]; ok = 1;
      elseif strcmp(vars,'p'), vars = 1; ok = 1;
      elseif strcmp(vars,'e'), vars = 2; ok = 1;
      else,disp('wrong input')
      end
    end
    disp(qfac_s)
    ok = 0;
    while ~ok
      qfac = input('Q-factor 1..4>','s');
      try
        qfac = str2num(qfac);
        if qfac>4 || qfac<1, disp('wrong input')
        else, ok = 1;
        end
      end
    end
    for v=length(vars)
      if vars(v)==1, caa_q.p(sc_list,lev) = qfac;
      elseif vars(v)==2, caa_q.e(sc_list,lev) = qfac;
      end
    end
    disp_q(caa_q)
  end
end

function disp_q(caa_q)
disp('sc\lev  1(p/e)  2(p/e)  3(p/e)')
for cl_id=1:4
  disp(sprintf('c%d      %d/%d      %d/%d      %d/%d',...
    cl_id,caa_q.p(cl_id,1),caa_q.e(cl_id,1),caa_q.p(cl_id,2),caa_q.e(cl_id,2),...
    caa_q.p(cl_id,3),caa_q.e(cl_id,3)))
end
