function OUT = c_ri_datestring_file(date)

%C_RI_DATESTRING_FILE Makes a string containing the date
%
% C_RI_DATESTRING_FILE(date) make a string of the date.
% The date is of form [Year Month Day Hour Min Sec].
%
% By Anders Tjulin
%
% The change is in the output. Less date separators 
%
%modified by Robert Isaksson in the summer of -03 

  if size(date) ~=[1,6]
    error('Wrong format of date')
  end
  
  year=num2str(date(1));
  
  if date(2)<10
    temp=num2str(date(2));
    month=['0' temp];
  else
    month=num2str(date(2));
  end
  
  if date(3)<10
    temp=num2str(date(3));
    day=['0' temp];
  else
    day=num2str(date(3));
  end
  
  if date(4)<10
    temp=num2str(date(4));
    hour=['0' temp];
  else
    hour=num2str(date(4));
  end
  
  if date(5)<10
    temp=num2str(date(5));
    minute=['0' temp];
  else
    minute=num2str(date(5));
  end
  
  if date(6)<10
    temp=num2str(round(date(6)));
    second=['0' temp];
  else
    second=num2str(round(date(6)));
  end
  
  OUT=[year '' month '' day 't' hour '' minute '' second];
