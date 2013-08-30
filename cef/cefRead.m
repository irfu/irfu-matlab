%
% D = cefRead(F)
%
%   Read a CEF file.
%
% [D,V] = cefRead(F) 
%
%   Read a CEF file F and return data D and variable information V.
%
% [D,V,H] = cefRread(F)
%
%   Read a CEF file F and return data D, variable information V and header
%   information H.
%
% [D,V,H] = cefRead(F, [N M])
%
%   Read only records between N<M for example [1 inf] would read all records.
%
% [D,V,H] = cefRead(F, [N M, DT])
% 
%   With sample intervall DT.
%   Note: This has not been implemented yet
%
% [D,V,H] = cefRead(F, {S, E, DT})
%
%   Another examples are {S, inf, DT}, {1, E, DT}
%   S start time tag, E end time tag.
%
% [D,V,H] = cefRead(F, {S, E}, PRODUCT)
% 
%   Select a specifig product,
%   Ex: PRODUCT = 'P12__C4_CP_EFW_L1_P12'
%   It is also possible to match multiple products by with the same
%   beginning. PRODUCT='B' would return 
%
%   B_vec_xyz_gse__C4_CP_FGM_SPIN: {3x1 cell}
%   B_mag__C4_CP_FGM_SPIN: {[600.5240]}
%
%   For a CAA/Cluster FGM file. Note that in current version 
%   this only works for the data field.
%
% [D,V,H] = cefRead(F, {S, E}, PRODUCT, PATH, OPTIONS)
%
%   PATH to header files can be a cell of multiple pathes like 
%   PATH = {'/tmp', '/tmp/efw/headers'} 
%
% [D,V,H] = cefRead(F, {S, E}, PRODUCT, PATH, OPTIONS)
%
%   OPTIONS give by cefGet
%
%   A couple of examples are given below:
%   To extract one single time_tag record from a efw file
%   which need headers. Also do not ignore any warnings.
%
%       opts=cefGet
%       opts.Ignorewarnings = {}
%       d=cefRead(file,[2 2],'time',{header_efw, header_global}, opts )
%
%       d =
%
%           time_tags__C4_CP_EFW_L1_P12: {'2001-07-06T06:00:00.062856Z'}
%
%
% See also cefWrite, cefGet
%

%
% By Josef Hook, joh@kth.se, jhook@rssd.esa.int, 2006/2007
% All rights reserved
function [D,varargout]=cefRead(F, varargin)

error(nargoutchk(0, 3, nargout))


%
% Default variables
%
NOUT = max(nargout, 1) - 1;
D=[];
H = {};
OPTS = cefGet;
I=[1 inf];
PROD=[];


%
% Do we have a specific product that we want to select
%
if(length(varargin)>=2)
    PROD=varargin{2};
end

%
% 3 varargin arguments 
%
if(length(varargin)>=3)             
    H = varargin{3};
%    disp('Setting header tag')
end % end of if(length(varargin)>=3)       


%
% Do we have options 
%
if(length(varargin)>=4)        
%
% 4 varargin arguments 
%
%disp('Setting options tag')
if(not(isempty(varargin{4})))
    OPTS = varargin{4};
end

end % end of if(length(varargin)>=4)       

%
% Set warning arguments
%
if(not(isempty(OPTS.IgnoreWarnings)))
    if(iscell(OPTS.IgnoreWarnings))
        for k=1:length(OPTS.IgnoreWarnings)
            warning('off', OPTS.IgnoreWarnings{k})
        end
    else 
        warning('off', OPTS.IgnoreWarnings)
    end
%    disp('setting warnings')
end




%
% Collect I param for cefReadData, cefReadMetaData
%
if(length(varargin)>=1)    
%
% 1 varargin argument 
%
        if(not(isempty(varargin{1}))) 
            
        %Check type of varargin{1}, we accept numeric and cell
        if(isnumeric(varargin{1})) 
            %disp('Read numerical records')
             if(min(size(varargin{1}))>1)
                 error('Second argument must be a vector with length 2 or 3')
             end
        
             switch(max(size(varargin{1}))) 
                    case 2
                    % We have an [n,m] vector 
                    %disp('We have an [n,m] vector')
                    I = varargin{1};   
                   case 3
                    % We have an [n,m,dt] vector 
                    %disp('We have an [n,m,dt] vector')
                    tmp = varargin{1};
                    I = tmp(1:2);
                    warning('We have not implemented the usage of dt yet')
                 otherwise
                 error('Second argument must be a vector with length 2 or 3')
             end
        
        elseif(iscell(varargin{1}))
              %disp('Cell records')
              %
              % Check if we have a string or a numeric, we allow syntax like
              % {'20020101T00:00:00Z', inf}
              %
               switch(max(size(varargin{1}))) 
                    case 2
                    %    
                    % We have an [n,m] vector 
                    %
                    %disp('We have an [n,m] vector')
                    vals = varargin{1};
%                    [num2str(ischar(vals{1})), num2str(ischar(vals{2}))]
                    %
                    % Switch over individuall cell values
                    %
                    switch [num2str(ischar(vals{1})), num2str(ischar(vals{2}))]
                        case ['1', '1']
                            start_end_times = [cefTimeToMjs(vals{1}), cefTimeToMjs(vals{2})]
                            
                            disp('Ok two strings try and convert these to times')
                            
                            %
                            % Search for index number corresponding to
                            % start time tag and end time tag. 
                            %
                            firstRecord = cefReadData(F,[1 1],'time_tags',[], H)
                            
                            disp('This feature is not implemented yet')    
                        case ['0', '1']
                            if(isnumeric(vals{1}))
                                I = [vals{1}, I(2)]; 
                                warning('Parsing of date Il not yet implemented')
                                disp('Ok first arg a value second arg a char')
                            else
                               error('Could not parse second argument')
                            end
                        case ['1', '0']
                            if(isnumeric(vals{2}))
                                I = [I(1), vals{2}]; 
                                warning('Parsing of date Il not yet implemented')
                                disp('Ok second arg a value first arg a char')
                            else
                               error('Could not parse second argument')
                            end
                            
                        otherwise
                            disp('Something is wrong in your syntax for the second argument Note: we dont accept {1,1} use [1, 1] instead')
                    end
                    
                   case 3
                   %
                   % We have an [n,m,dt] vector 
                   %
                   %disp('We have an [n,m,dt] vector')
                    vals = varargin{1};
%                    [num2str(ischar(vals{1})), num2str(ischar(vals{2}))]
                    if(not(isnumeric(vals{3})))
                        error('dt must be a integer')
                    end
                    dt = vals{3};    
                    %
                    % Switch over individuall cell values
                    %
                    switch [num2str(ischar(vals{1})), num2str(ischar(vals{2}))]
                        case ['1', '1']
                            disp('Ok two strings try and convert these to times')
                            
                            
                            disp('This feature is not implemented yet')    
                        case ['0', '1']
                            if(isnumeric(vals{1}))
                                I = [vals{1}, I(2)]; 
                                warning('Parsing of date Il not yet implemented')
                                disp('Ok first arg a value second arg a char')
                            else
                               error('Could not parse second argument')
                            end
                        case ['1', '0']
                            if(isnumeric(vals{2}))
                                I = [I(1), vals{2}]; 
                                warning('Parsing of date Il not yet implemented')
                                disp('Ok second arg a value first arg a char')
                            else
                               error('Could not parse second argument')
                            end
                            
                        otherwise
                            disp('Something is wrong in your syntax for the second argument Note: we dont accept {1,1} use [1, 1] instead')
                    end

                   
                   
                   otherwise
                     error('Second argument must be a vector with length 2 or 3')
               end
             
        else
            error('Could not parse second argument, must be a vector or a cell')
        end
        
        end % end of if(not(isempty(varargin{1})))
end % end of if(length(varargin)>=1)       

if(NOUT == 1) 
    
[header, var]=cefReadMetaData(F,H{:});
varargout(1) = {var};

end

if(NOUT == 2)
    
[header, var]=cefReadMetaData(F,H{:});
varargout(1) = {var};
varargout(2) = {header};    

end
D=cefReadData(F,I,PROD,[],H{:});



