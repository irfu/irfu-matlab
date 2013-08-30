%
% R=cefChunk(D,T,S,E) Extracts chunk of data D given start time S and end time E and times T.
% The data must be in the format given by the cefRead... functions
%
% Ex: 
% ut = cefChunk(e,t,'2002-06-12T07:29:54.0Z','2002-06-12T07:29:58.0Z')
function R=cefChunk(D,T,S,E)

%
% Assume time given in cef string format 2002-02-02T00:00:00Z
%
if(isempty(S) && isempty(E))
    R=D;
else
    

if(isempty(S))
    S=cefTimeToMjs(T(1));
end
if(isempty(E))
    E=cefTimeToMjs(T(end));
end
    
if(isstr(S)) 
    S=cefTimeToMjs(S);
end

if(isstr(E))
    E=cefTimeToMjs(E);
end

if(iscell(T))
    T=cefTimeToMjs(T);
end

a=find(T>=S);
b=find(T(a)<=E);

if(size(D,1)>size(D,2))
    Da=D(a,:);
    R=Da(b,:);
else
    Da=D(:,a);
    R=Da(:,b);    
end


end