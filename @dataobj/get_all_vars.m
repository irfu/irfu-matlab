function out = get_all_vars(dobj,varargin)

if nargin>1 && strcmp(varargin,'ts')
  doGetTs = 1;
else
  doGetTs = 0;
end

if isempty(dobj.data)
	isDataobjEmpty = true;
	disp('dataobj is empty!');
	disp(' ');
else
	isDataobjEmpty = false;
end

nvars = size(dobj.vars,1);
vars = cell(1,nvars);
if nvars>0
	for v=1:nvars
		if isDataobjEmpty
			disp(dobj.vars{v,1});
    else
      tmp_var = get_variable(dobj,dobj.vars{v,1});
      if doGetTs
        vars{1,v} = mms.variable2ts(tmp_var);
      else
        vars{1,v} = tmp_var;
      end
		end
	end
else
	disp('No variables in data object');
end
out = vars;