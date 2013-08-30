%
%
% Function for plotting the entire data structure with subplots
%
%
% Ex: cefPlot(cefReadData(file2,[1 inf],[],[],header_fgm, header_glob),1) 
%     cefPlot(cefReadData(file2,[1 inf],'B_vec_xyz_gse__C4_CP_FGM_SPIN',[],header_fgm, header_glob),0)
%
%
% Second argument should be set to 1 if the input structure has time vector 
% as its first element. 
% Note: This behaviour is now deprecated you only need to specify onte
% argument like  cefPlot(data)
%
function cefPlot(data,varargin)

     if(nargin == 0)
     error('input is not defined')
     end
     if(isempty(data))
     error('input is empty')
     end
     if(isstruct(data) == 0)
     error('input is not a structure')
     end

     fields=fieldnames(data);

     figure,
     warning off
     for k=1:length(fields) 

     subplot(length(fields),1,k), 
     if(k==1) 
         try
           plot((cefTimeToMjs(getfield(data,cell2mat(fields(k)) )) )')
         catch
           plot(transpose(cell2mat(getfield(data, cell2mat(fields(k)) ) ) ))         
         end
     else
       plot(transpose(cell2mat(getfield(data, cell2mat(fields(k)) ) ) ))
     end
     
     h=title(cell2mat(fields(k)));
     set(h,'Interpreter','none');
     end
       warning on


end
