function c_ri_run_events_into_pictures(s_t,e_t,p_MP,p_Bp,p_E,p_Out, period)
%
%Input:
%
%Output:
%
%Descrition of the function:
%
%Using:
%
%Work method:
%
%Error:
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
%writing the ls-command to a file, which can be opened by matlab
file_prefix = 'E_';
start_path = pwd;

cd(p_E);

ls_out = sprintf('%sls_E_files.txt',p_Out);
unix_command = sprintf('ls %s*.mat >%s' ,file_prefix, ls_out );
unix(unix_command);

fp = fopen(ls_out, 'r');
% continue until end of file
while feof(fp) == 0
  file_name = fgetl(fp);
  f_length = length(file_name);

  if f_length > 6
    %excluding bad filenames
    if strcmp(file_name(1:2),file_prefix) && strcmp(file_name(f_length-3:f_length), '.mat')

      if c_ri_timestr_within_intervall_E(file_name,s_t,e_t) == 1
        c_ri_events_into_pictures(p_E,file_name,period,p_Out,p_Bp,p_MP);
      end

    end

  end
end
cd(start_path);
