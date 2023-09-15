function c_ri_get_many_B(passing_MP,path_output)
%function c_ri_get_many_B(passing_MP,path_output)
%
%Input:
% passing_MP -matrix where every row contains [entering MP | leaving MP]
%             entering and leaving must be in the same day.
%
%Output:
% Saved to file:
% All B-files
%
%Descrition of the function:
% Load all available B-data for every M
%
%Using:
% datestring
% create timetable
% get_B
%
%Work method:
% Creates two timetables, containing the downloadperiods for data
% in burstmode and normalmode. Then the data is downloaded and
% saved to file.
%
%Error:
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
if nargin == 1
  path_output = [pwd '/'];
end

[nr_cross, col] = size(passing_MP);

for mp_cross = 1:nr_cross
  pack
  from = passing_MP(mp_cross,1);
  to = passing_MP(mp_cross,2);
  disp([ 'mp crossing' datestring(fromepoch(from)) ' ' datestring(fromepoch(to))]);

  for cl_nr = 1:4
    b_table = 0;
    n_table = 0;

    [b_table, n_table] = create_timetable(from,to,cl_nr);

    if b_table ~= 0
      [b_s, col] = size(b_table);
      for k = 1:b_s
        c_ri_get_B( b_table(k,1), b_table(k,2), cl_nr, 'b',path_output);
      end
    end

    if n_table ~= 0
      [n_s, col] = size(n_table);
      for i = 1:2:n_s
        c_ri_get_B( n_table(i,1), n_table(i,2), cl_nr, 'n',path_output);
      end
    end
  end

end
