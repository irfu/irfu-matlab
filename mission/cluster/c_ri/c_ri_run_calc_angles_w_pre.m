function c_ri_run_calc_angles_w_pre(s_t, e_t, path_input, path_output)
%
% c_ri_run_calc_angles_w_pre(s_t, e_t, path_input, path_output)
%
%Input:
% s_t, e_t - ex:[2002 03 02 3 0 0]
% path_input, path_output - end with / ex: /share/ang/. Can use default
%
%Output:
% Saves to file: /path_output/Ap_020302_F030100_T045500_b.mat
%
%Descrition of the function:
% Finds all files within the time intervall, loads the files, claculates
% the angles and saves the result to file.
%
%Using:
% c_ri_timestr_within_intervall
% c_ri_calc_angle_w_preprocess
%
%Work method:
%
%Error:
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
if nargin == 3
  path_input = [pwd '/'];
  path_output = [pwd '/'];
end

%writing the ls-command to a file, which can be opened by matlab
file_prefix = 'Bp_';
start_path = pwd;

cd(path_input);

ls_out = sprintf('%sls_B_data.txt',path_output);
unix_command = sprintf('ls %s*.mat >%s' ,file_prefix, ls_out );
unix(unix_command);

fp = fopen(ls_out, 'r');
% continue until end of file
while feof(fp) == 0
  file_name = fgetl(fp);
  f_length = length(file_name);

  %excluding bad filenames
  if strcmp(file_name(1:3),file_prefix) && strcmp(file_name(f_length-3:f_length), '.mat')

    if c_ri_timestr_within_intervall(file_name,s_t,e_t) == 1

      c_ri_calc_angle_w_preprocess(path_input, file_name, path_output);

    end

  end

end

cd(start_path);
