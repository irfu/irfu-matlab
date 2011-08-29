function res = getmat(dobj,var_s)
%GETMAT(dobj, var_s)  get a variable in the Cluster AV format
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

data = getv(dobj,var_s);
fillv = getfillval(dobj,var_s);
if ~ischar(fillv),
    data.data(data.data==fillv) = NaN;
else
    disp('fill value is character: discarding')
end

dep = getdep(dobj,var_s);
dim = length(data.variance(3:end));

if dim <=1
    plot_data = double(data.data)';
    
    if isfield(dep,'DEPEND_O'), res = [dep.DEPEND_O plot_data'];
    else res = plot_data;
    end
else
    for d = 1:size(dep.DEPEND_X,1)
        dep_x{d} = getv(dobj,dep.DEPEND_X{d,1});
        dep_x{d}.s = dep.DEPEND_X{d,1};
        dep_x{d}.fillv = getfillval(dobj,dep_x{d}.s);
        if isnumeric(dep_x{d}.fillv), % only implemented for numeric data 
            dep_x{d}.data(dep_x{d}.data==dep_x{d}.fillv) = NaN;
        end
        dep_x{d}.units = getunits(dobj,dep_x{d}.s);
        dep_x{d}.lab = getlablaxis(dobj,dep_x{d}.s);
    end
    if isnumeric(data.FILLVAL), % put fillvalues to NaN
        data.data(data.data==data.FILLVAL) = NaN;
    end
    res = struct('t',dep.DEPEND_O,'dep_x',[],'data',data.data);
    res.dep_x = dep_x;
end
