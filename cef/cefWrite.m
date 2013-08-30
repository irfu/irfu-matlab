%
% cefWrite(F,D,V)
%
%   Write data D and Variable metadata V to file F. Returns 0 if success.
%
% cefWrite(F,D,V,H)
%
% cefWrite(F,D,V,H,INC)
%
% cefWrite(F,D,V,H,INC,OPTS)
%
% Note: Options tag not used in current version
%
% See also cefRead
function cefWrite(F,D,V, varargin)

if(not(ischar(F))) 
    error('First argument not a file')
end

if(not(isstruct(V)))
    error('Variable metadata must be a struct of structures')
end
fn=fieldnames(V);
    
if(not(isstruct(V.(fn{1}))))
   error('Variable metadata must be a struct of structures')
end

%
% Data D is given in a raw format we then need to repackage it.
% We also need the variable V.
%
if(not(isstruct(D)))
disp('Data not in structure format repackaging')    
data=cefPackageData(D,V);
    
else 
    
data = D;

end

%
% Load options
%
if(length(varargin)>=3)
 OPTS = varargin{3};
else
 OPTS = cefGet;   
end


%
% If we have specified a Header write header data.
%
if(length(varargin)>=1)

    H=varargin{1};
    if(not(isstruct(H)))
       error('Metadata header must be a struct') 
    end
    

    %
    %  Write include tags to header? 
    %
    if(length(varargin)>=2)
   % disp('writing meta header with include tags')
    if(not(iscell(varargin{2})))
        error('Include tags must be given in a cell structure')
    end
        
    cefWriteMetaData(H,V,varargin{2},F,'w');
    else
  %          disp('writing meta header')
        cefWriteMetaData(H,V,[],F,'w');
    end

end
%disp('writing data')
cefWriteData(data,V,F);

