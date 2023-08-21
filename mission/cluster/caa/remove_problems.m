% script REMOVE_PROBLEMS
%
% Remove problems from data
%
% Input: signal,probe,cl_id,problems
% Output: res

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

res = signal;
if probe>10
  allProbes = [12 34 32 42];
  switch probe
    case 12, p_list = [1,2];
    case 32, p_list = [3,2];
    case 34, p_list = [3,4];
    case 42, p_list = [4,2];
    otherwise, error('Unknown probe')
  end
elseif probe>0 && probe <=4, p_list = probe;
else, error('Unknown probe')
end

param = tokenize(problems,'|');
for i=1:length(param)

  switch lower(param{i})
    case 'reset'
      % Remove bad bias around EFW reset
      [ok,bbias,msg] = c_load('BADBIASRESET?',cl_id);
      if ok
        if ~isempty(bbias)
          irf_log('proc','blanking bad bias due to EFW reset')
          res = caa_rm_blankt(res,bbias);
        end
      else, irf_log('load',msg)
      end
      clear ok bbias msg

    case 'bbias'
      % Remove bad bias from bias current indication
      for kk = p_list
        [ok,bbias,msg] = c_load(irf_ssub('BADBIAS?p!',cl_id,kk));
        if ok
          if ~isempty(bbias)
            irf_log('proc',['blanking bad bias on P' num2str(kk)])
            res = caa_rm_blankt(res,bbias);
          end
        else, irf_log('load',msg)
        end
        clear ok bbias msg
      end

    case 'hbiassa'
      % Remove saturation due to too high bias current
      [ok,wake,msg] = c_load(irf_ssub('HBIASSA?p!',cl_id,probe));
      if ok
        if ~isempty(wake)
          irf_log('proc','blanking HB saturation')
          res = caa_rm_blankt(res,wake);
          clear wake
        end
      else, irf_log('load',msg)
      end
      clear ok wake msg
      % Remove saturation due to too high bias current from nonstandard operations intervals
      [ok,nsops,msg] = c_load('NSOPS?',cl_id);
      if ok
        if ~isempty(nsops)
          idx = nsops(:,3)==caa_str2errid('high_bias');
          if any(idx)
            irf_log('proc','blanking HB saturation (NS_OPS)')
            res = caa_rm_blankt(res,nsops(idx,1:2));
          end
        end
      else, irf_log('load',msg)
      end
      clear ok nsops msg

    case 'saasa'
      % Remove probe saturation due to SAA
      if probe>10
        [ok,saasa,msg] = c_load('SAASADI?',cl_id);
        if ~isempty(saasa)
          saasa = saasa(:,2*find(probe==allProbes) -[1 0]);
          irf_log('proc',['blanking SAA saturation on P' num2str(probe)])
          res = caa_rm_blankt(res,saasa);
        end
      else
        [ok,saasa,msg] = c_load('SAASASE?',cl_id);
        if ~isempty(saasa)
          saasa = saasa(:,2*probe-[1 0]);
          irf_log('proc',['blanking SAA saturation on P' num2str(probe)])
          res = caa_rm_blankt(res,saasa);
        end
      end

    case 'probesa'
      % Remove probe saturation
      for kk = p_list
        [ok,sa,msg] = c_load(irf_ssub('PROBESA?p!',cl_id,kk));
        if ok
          if ~isempty(sa)
            irf_log('proc',['blanking saturated P' num2str(kk)])
            res = caa_rm_blankt(res,sa);
          end
        else, irf_log('load',msg)
        end
        clear ok sa msg
      end

    case 'probeld'
      % Remove probe saturation due to low density
      for kk = p_list
        [ok,sa,msg] = c_load(irf_ssub('PROBELD?p!',cl_id,kk));
        if ok
          if ~isempty(sa)
            irf_log('proc',...
              ['blanking low density saturation on P' num2str(kk)])
            res = caa_rm_blankt(res,sa);
          end
        else, irf_log('load',msg)
        end
        clear ok sa msg
      end

    case 'whip'
      % Remove whisper pulses
      [ok,whip,msg] = c_load('WHIP?',cl_id);
      if ok
        if ~isempty(whip)
          irf_log('proc','blanking Whisper pulses')
          res = caa_rm_blankt(res,whip);
        end
      else, irf_log('load',msg)
      end
      clear ok whip msg

    case 'sweep'
      % Remove sweeps
      [ok,sweep,msg] = c_load('SWEEP?',cl_id);
      if ok
        if ~isempty(sweep)
          irf_log('proc','blanking sweeps')
          res = caa_rm_blankt(res,sweep);
          clear sweep
        end
      else, irf_log('load',msg)
      end
      clear ok sweep msg

    case 'spike'
      % Remove spike
      [ok,spike,msg] = c_load('SPIKE?',cl_id);
      if ok
        if ~isempty(spike)
          irf_log('proc','blanking spike')
          res = caa_rm_blankt(res,spike);
          clear spike
        end
      else, irf_log('load',msg)
      end
      clear ok spike msg

    case 'bdump'
      % Remove burst dumps
      [ok,bdump,msg] = c_load('BDUMP?',cl_id);
      if ok
        if ~isempty(bdump)
          irf_log('proc','blanking burst dumps')
          res = caa_rm_blankt(res,bdump);
          clear bdump
        end
      else, irf_log('load',msg)
      end
      clear ok bdump msg

    case 'nsops'
      % Remove nonstandard operations intervals
      [ok,nsops,msg] = c_load('NSOPS?',cl_id);
      if ok
        if ~isempty(nsops)
          if ~exist('nsops_errlist','var'),nsops_errlist=[];end
          for j=1:length(nsops(:,3))
            opcode=nsops(j,3);
            % If nsops_errlist is present, only match on those
            % opcodes in the list
            if (opcode<10 && opcode>0) || ((opcode>=11 && opcode<=14) && any(opcode-10==p_list)) || (~isempty(nsops_errlist) && any(opcode==nsops_errlist))
              irf_log('proc',['blanking nsops interval. opcode:' num2str(opcode) ' probe:' num2str(probe)]);
              res = caa_rm_blankt(res,nsops(j,:));
            end
          end
        end
      else, irf_log('load',msg)
      end
      clear ok nsops msg opcode

    case 'wake'
      % Remove wakes
      [ok,wake,msg] = c_load(irf_ssub('PSWAKE?p!',cl_id,probe));
      if ok
        if ~isempty(wake)
          irf_log('proc','blanking plasmaspheric wakes')
          res = caa_rm_blankt(res,wake);
          clear wake
        end
      else, irf_log('load',msg)
      end
      clear ok wake msg

      [ok,wake,msg] = c_load(irf_ssub('LOWAKE?p!',cl_id,probe));
      if ok
        if ~isempty(wake)
          irf_log('proc','blanking lobe wakes')
          res = caa_rm_blankt(res,wake);
          clear wake
        end
      else, irf_log('load',msg)
      end
      clear ok wake msg

    otherwise
      error('Unknown parameter')
  end
end
