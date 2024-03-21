function c_ri_load_B_for_preprocess(s_t, e_t, path_input, path_output)
%
% load_B_for_preprocess(s_t, e_t, path_input, path_output)
%
%Input:
% s_t -start time ex [2002 03 02 1 2 3]
% e_t -end time ex [2002 03 02 1 2 3]
% path_input -path to input cataloge (can use default)
% path_output -path to output cataloge (can use default)
%
%Output:
% preprocessed data saved to cataloge "path_output"
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
if nargin == 2
  path_input = [pwd '/'];
  path_output = [pwd '/'];
end

%writing the ls-command to a file, which can be opened by matlab
file_prefix = 'Bm_';
start_path = pwd;

cd(path_input);

ls_out = sprintf('%sls_B_data.txt',path_output);
unix_command = sprintf('ls %s*.* >%s' ,file_prefix, ls_out );
unix(unix_command);

fp = fopen(ls_out, 'r');

% continue until end of file
while feof(fp) == 0
  pack
  file_name = fgetl(fp);
  f_length = length(file_name);

  %excluding bad filenames
  if strcmp(file_name(1:3),file_prefix) && strcmp(file_name(f_length-3:f_length), '.mat')
    if c_ri_timestr_within_intervall(file_name,s_t,e_t) == 1

      disp(file_name)
      p_and_f = sprintf('%s%s',path_input,file_name);
      load(p_and_f)
      mode = file_name(27);

      [B1, B2, B3, B4] = preprocess_B(B1,B2,B3,B4,mode);

      n_p_and_f = sprintf('%sBp%s',path_output, file_name(3:f_length-4));

      disp(['saving to: ' n_p_and_f])
      save(n_p_and_f,'B1', 'B2' ,'B3', 'B4')

    end
  end

end

cd(start_path);
