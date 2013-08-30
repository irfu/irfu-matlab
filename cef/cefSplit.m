%
% Function for splitting strings.
% Ex: 
%
%  [a,b]=cefSplit('C4_CP_EFW_L1_P12__20010706_060000_064341_V01.cef','00','g')
%  [a,b]=cefSplit('C4_CP_EFW_L1_P12__20010706_060000_064341_V01.cef','__')
% 
% Third argument may either be set to 'g' which stands for global or it may
% be undefined
%
function [left, right]=cefSplit(str, splt, varargin)


if(not(ischar(str)))
        error('First argument not a string')
end

if(not(ischar(splt))) 
        error('Second argument not a charachter or string')
end

if(isempty(varargin))
    varargin(1) = {'l'};
end

if(length(splt)==1)
   
    
    if(strcmp(varargin(1),'g')) 
    [left, right]=strtok(str,splt);
        right=right(2:end);
        while(any(right))
           [dleft, right]=strtok(right, splt);
           left = strvcat(left, dleft);
           right=right(2:end);
        end
    else 
        [left, right]=strtok(str,splt);
        right=right(2:end);
    end
else
    
    if(strcmp(varargin(1), 'g'))
        
    pos=findstr(splt, str);
    left=str(1:(pos-1));
    right=str((pos+length(splt)):end);
  
    if(length(pos)>1)
    
    for k=2:length(pos)
       left=strvcat(left,str(pos(k-1):(pos(k)-1)));
 %      right=right((pos+length(splt)):end)
    end
    left=strvcat(left, str((pos(k)+length(splt)):end));
  right=[];  
    end
    
    else
    pos=findstr(splt, str);
    left=str(1:(pos-1));
    right=str((pos+length(splt)):end);
    end
    
end





