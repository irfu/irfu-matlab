function c=cefVar2cdfVar(v)

vnames=fieldnames(v);

a = v.(vnames{1});
b = v.(vnames{2});

fa=fieldnames(a);
fb=fieldnames(b);
[ins,rfa,rfb]=intersect(fa,fb);
%comvars=union(fa,fb);

comvars = fa;
for i=1:length(fb);
   if(isempty(strmatch(fb{i}, ins)));
comvars =       {comvars{:}, fb{i}};
   end        
end
comvars = comvars';

c=cell2struct(cell(1,length(comvars)), comvars,2);

%
% populate b values
%
l=1
for i=1:length(fb)
    
        if(iscell(b.(fb{i})))
          d{l,1} = [[vnames{2},'_',fb{i}], b.(fb{i})]
          c.(fb{i}) =  {vnames{2}, [vnames{2},'_',fb{i}]}
          l=l+1
        else        
          c.(fb{i}) =  {vnames{2}, b.(fb{i})};
        end
end

%
% populate a values
%
for i=1:length(fa)
        
        if(iscell(a.(fa{i})))
          d{l,1} = [[vnames{2},'_',fb{i}], a.(fa{i})]
          c.(fa{i}) =  {vnames{1}, [vnames{1},'_',fa{i}]}
          l=l+1
        else        
          c.(fa{i}) =  {vnames{1}, a.(fa{i})};
        end    
 end

%
% Populate intersecting values
%
for i=1:length(ins)
   tmp =  c.(ins{i});
   tmp2 = {vnames{2}, b.(ins{i})};
   c.(ins{i}) = { tmp{:}; tmp2{:}};
end


%
% Repackage c and d
%
c = {c, d}
