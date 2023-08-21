function c_ri_run_get_B(s_t,e_t,path_input,path_output)
%
% c_ri_run_get_B(s_t,e_t,path_input,path_output)
%
%Input:
% s_t,e_t -start and endtime in ex: [2002 03 02 0 0 0].
% path_input -path to MP-files
%           (MP-file ex: "MP_20020101_00:00:00_to_20020201_00:00:00")
% path_output -where the Ba.......01 file will be saved
%
%Output:
% saves to output
%
%Descrition of the function:
% Downloads the B-files for the time intervall given in input.
%
%Using:
% c_ri_get_many_B
%
%Work method:
% Makes a file where the result of the ls command is saved. This file is
% loaded and for every line the file is loaded. If the loaded "passing_MP"
% has a time intervall within the time intervall, given in input, the the section
% is loaded.
%
%Error:
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
if nargin == 0
  path_input = [pwd '/'];
  path_output = [pwd '/'];
  s_t = 0;
  e_t = 0;
end

if nargin == 2
  path_input = [pwd '/'];
  path_output = [pwd '/'];
end

s_t_e = toepoch(s_t);
e_t_e = toepoch(e_t);

%writing the ls-command to a file, which can be opened by matlab
ls_out = sprintf('%sls_list.txt',path_output);
unix_command = sprintf('ls %sMP_*.* >%s' , path_input, ls_out );
unix(unix_command);

fp = fopen(ls_out, 'r');

% continue until end of file
while feof(fp) == 0
  f_line = fgetl(fp);
  passing_MP = 0;

  if f_line ~= -1
    %assigns to value to passing_MP
    load(f_line);
  end

  if passing_MP ~= 0

    if s_t == 0 && e_t == 0
      tmp = passing_MP;

    else
      tmp1 = passing_MP(:,1);
      tmp2 = passing_MP(:,2) - passing_MP(:,1);
      tmp = [tmp1 tmp2];
      [s_row,e_row] = find_row(s_t_e, e_t_e, tmp,1);

      if s_row == -1 || e_row == -1
        p_mp = -1;
      else
        p_mp = passing_MP(s_row:e_row,:);
      end

    end

    if p_mp(1,1) ~= -1
      c_ri_get_many_B(p_mp,path_output);
    end

  end

end

fclose(fp);
