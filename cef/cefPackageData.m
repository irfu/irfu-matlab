%
% Builds a correct data structure for cefWriteData from raw data and a structure of meta variables 
% 
% By Josef Hook, jhook@rssd.esa.int
%   Ex:
%   ut=cefPackageData(pv, metadata_variables)
%
function ret=cefPackageData(data, mvars)


fields=fieldnames(mvars);
[row,col]=size(data);
if(row<col && not(row==col)) 
    data=data'; 
end

    

ret_m={};
n=1;
for m=1:size(fields,1)
    ds=getfield(getfield(mvars, fields{m}), 'SIZES');
   ret_m=[ret_m, {num2cell(data(:,n:n+ds-1)')}];
    n=n+ds;
end

ret=cell2struct(ret_m, fields,2);





