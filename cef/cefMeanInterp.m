%
% R=cefMeanInterp(A,B,HB,C,HC) performs a mean interpolation of matrix A given the time
% axis B and C. We also need the halfintervall HB for B and HC for C.
%
%
%
%
% See also ceflib
% Matlab version of cefMeanInterp
% 
function ret=cefMeanInterp(a,b,c)


if(not(isnumeric(a))) 
 error('First argument, has to be a numeric matrix')
end

if(ischar(b) || iscell(b)) 
    b = cefTimeToMjs(b);
elseif(not(isnumeric(b)))
        error('Second argument, time line must be given either as a CEF formated date string or as a numeric valued vector')
end

if(ischar(c) || iscell(c)) 
    c = cefTimeToMjs(c);
elseif(not(isnumeric(c)))
        error('Third argument, time line must be given either as a CEF formated date string or as a numeric valued vector')
end

[nb,mb]=size(b);
[nc,mc]=size(c);
if(nb ~= 1 && mb ~= 1)
    error('Second argument has to be a vector.')
end

if(nc ~= 1 && mc ~= 1)
    error('Third argument has to be a vector.')
end


[n,m]=size(a);
do_transp = 0;

if(n<m) 
    do_transp=1;
    a=a';
end


ret=[]
epsilon=1;

k=1
iold=0
for i=1:length(b)
    
    if(abs(b(i) - c(k))<epsilon) 
        
        ret=[ret; mean(a((iold+1):i,:))]  
            k=k+1;
            iold=i
    end
    if(k==length(c))
       ret=[ret; mean(a((i+1):end,:))]
       break 
    end
end



if(do_transp)
    ret=ret'; 
end