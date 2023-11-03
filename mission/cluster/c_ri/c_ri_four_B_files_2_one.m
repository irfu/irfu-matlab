function c_ri_four_B_files_2_one(path_input, path_output , s_t, e_t)
%
%c_ri_four_B_files_2_one(path_input, path_output, s_t, e_t)
%
%Input:
% path_input -the path to the catalog where the files are. (include / in the end)
%            ex "Ba_020605_F050939_T191504_n.04"
% path_output -the path to where the outputfiles should be (include / in the end)
% s_t -start time ex [2002 03 02 0 0 0]
% e_t -end time ex [2002 03 03 0 0 0]
%
%Output:
% saves the files to path_output with the same filename exept for
% the extension. (*.mat instead of *.01 *.02 *.03 *.04)
%
%Descrition of the function:
% Compess 4 files, each file is the Ba-file form a cluster satellite, to
% one file containing
%
%Using:
% find_str
% load_file
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
ls_out = sprintf('%sBa_list.txt',path_output);
start_dir = pwd;
cd(path_input)

unix_command = sprintf('ls > %s' , ls_out )
unix(unix_command);

fp = fopen(ls_out, 'r');

% continue until end of file
while feof(fp) == 0

  %find a filename with the end .01 to start with satellite 1
  %and the starttime is not the end time and that the file is
  %within the intervall.
  correct_l = -1;
  temp_l = 0;
  while  feof(fp) == 0 && correct_l == -1
    temp_l = fgetl(fp);
    i = length(temp_l);
    %if the file is witin the time intervall

    %checks the filename
    if i == 30 && temp_l(1:3) == 'Ba_'
      if c_ri_timestr_within_intervall(temp_l,s_t,e_t) == 1
        %finds the file for cluster 1
        if strcmp(temp_l(i-2:i),'.01') == 1
          % if start and end are not the same
          if strcmp(temp_l(11:16),temp_l(19:24)) == 0
            correct_l = 1;
          end
        end
      end
    end
  end % end of while, correct_l ~= 1 or eof

  %if a single correct line is found
  if correct_l == 1
    b_file_01 = temp_l;
    b_file_02 = find_str(b_file_01,fp,'2');
    b_file_03 = find_str(b_file_01,fp,'3');
    b_file_04 = find_str(b_file_01,fp,'4');

    mode = b_file_01(length(b_file_01)-3);

    B1 = load_file(path_input, b_file_01); %#ok<NASGU>
    B2 = load_file(path_input, b_file_02); %#ok<NASGU>
    B3 = load_file(path_input, b_file_03); %#ok<NASGU>
    B4 = load_file(path_input, b_file_04); %#ok<NASGU>

    filename = b_file_01(3:length(b_file_01)-3);
    p_and_f = sprintf('%sBm%s',path_output,filename);
    save(p_and_f, 'B1', 'B2', 'B3', 'B4')

    disp(['saved to file: ' p_and_f]);

  end

end % eof

cd(start_dir);
