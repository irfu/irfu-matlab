function name2 = find_str(name1, fp, f_end)
%
%Input:
% name1 -stringname of models file,ends with .01
% fp -filepointer
% f_end - 2,3,4
%
%Output:
% name2 -the filename
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

n1 = length(name1);
t_name = sprintf('%s%s',name1(1:n1-1),f_end);
name2 = 0;

cont = 1;
while feof(fp) == 0 && cont == 1
  temp_l = fgetl(fp);
  t = length(temp_l);

  if n1 == t
    if strcmp(t_name,temp_l) == 1
      name2 = t_name;
      cont = 0;
    end
  end

  if strcmp(temp_l(t), '1') == 1
    cont = 0;
  end

end
