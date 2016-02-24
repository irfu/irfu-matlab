function CBH = cbhandle(varargin)
%CBHANDLE   Gets the handle of current colorbar or its peer axes.
%
%   SYNTAX:
%     CBH = cbhandle;
%     CBH = cbhandle(H,...);
%     CBH = cbhandle(...,'force');
%     CBH = cbhandle(...,'unhide');
%     PAH = cbhandle(...,'peer');
%
%   INPUTS:
%     'force'  - Forces to create colorbars when not found.
%     'unhide' - Looks even in hidden handles.
%     'peer'   - Looks for its peer axes handle instead.
%     H        - Colorbars, axes, figures or uipanels handles to look for
%                colorbars or its peer handles.
%                DEFAULT: gca (current axes)
%
%   OUTPUTS:
%     CBH - Colorbar or its peer axes handle(s).
%
%   DESCRIPTION:
%     By default, colorbars are hidden objects. This function searches for
%     them basically by its 'axes' type and 'Colorbar' tag.
%
%   SEE ALSO:
%     COLORBAR
%     and
%     CBUNITS, CBLABEL, CBFREEZE by Carlos Vargas
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   cbhandle.m
%   VERSION: 2.1 (Jul 03, 2014) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>)
%   MATLAB:  8.2.0.701 (R2013b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Jun 08, 2009)
%   1.1      Fixed bug with colorbar handle input. (Aug 20, 2009)
%   2.0      Rewritten code. New 'create', 'force', 'unhide' and 'peer'
%            optional inputs. Changed application data to 'cbfreeze'. (Jun
%            05, 2014)
%   2.1      Fixed small bug with input reading. (Jul 03, 2014)

%   DISCLAIMER:
%   cbhandle.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2009-2014 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Parameters:
appName = 'cbfreeze';

% Sets default:
H      = get(get(0,'CurrentFigure'),'CurrentAxes');
FORCE  = false;
UNHIDE = false;
PEER   = false;

% Checks inputs/outputs:
assert(nargin<=5,'CAVARGAS:cbhandle:IncorrectInputsNumber',...
    'At most 5 inputs are allowed.')
assert(nargout<=1,'CAVARGAS:cbhandle:IncorrectOutputsNumber',...
    'Only 1 output is allowed.')

% Gets H: Version 2.1
if ~isempty(varargin) && ~isempty(varargin{1}) && ...
        all(reshape(ishandle(varargin{1}),[],1))
    H = varargin{1};
    varargin(1) = [];
end

% Gets UNHIDE:
while ~isempty(varargin)
    if isempty(varargin) || isempty(varargin{1}) || ~ischar(varargin{1})...
            || (numel(varargin{1})~=size(varargin{1},2))
        varargin(1) = [];
        continue
    end
    
    switch lower(varargin{1})
        case {'force','forc','for','fo','f'}
            FORCE  = true;
        case {'unhide','unhid','unhi','unh','un','u'}
            UNHIDE = true;
        case {'peer','pee','pe','p'}
            PEER   = true;
    end
    varargin(1) = [];
end

% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

% Show hidden handles:
if UNHIDE
    UNHIDE = strcmp(get(0,'ShowHiddenHandles'),'off');
    set(0,'ShowHiddenHandles','on')
end

% Forces colormaps
if isempty(H) && FORCE
    H = gca;
end
H = H(:);
nH = length(H);

% Checks H type:
newH = [];
for cH = 1:nH
    switch get(H(cH),'type')
        
        case {'figure','uipanel'}
            % Changes parents to its children axes
            haxes = findobj(H(cH), '-depth',1, 'Type','Axes', ...
                '-not', 'Tag','legend');
            if isempty(haxes) && FORCE
                haxes = axes('Parent',H(cH));
            end
            newH = [newH; haxes(:)];
            
        case 'axes'
            % Continues
            newH = [newH; H(cH)];
    end
    
end
H  = newH;
nH = length(H);

% Looks for CBH on axes:
CBH = NaN(size(H));
for cH = 1:nH
    
    % If its peer axes then one of the following is not empty:
    hin  = double(getappdata(H(cH),'LegendColorbarInnerList'));
    hout = double(getappdata(H(cH),'LegendColorbarOuterList'));
    
    if ~isempty([hin hout]) && any(ishandle([hin hout]))
        % Peer axes:

        if ~PEER
            if ishandle(hin)
                CBH(cH) = hin;
            else
                CBH(cH) = hout;
            end
        else
            CBH(cH) = H(cH);
        end
        
    else
        % Not peer axes:
        
        if isappdata(H(cH),appName)
            % CBFREEZE axes:
            
            appdata = getappdata(H(cH),appName);
            if ~PEER
                CBH(cH) = double(appdata.cbHandle);
            else
                CBH(cH) = double(appdata.peerHandle);
            end
            
        elseif strcmp(get(H(cH),'Tag'),'Colorbar')
            % Colorbar axes:
            
            if ~PEER
                
                % Saves its handle:
                CBH(cH) = H(cH);
                
            else
                
                % Looks for its peer axes:
                peer = findobj(ancestor(H(cH),{'figure','uipanel'}), ...
                    '-depth',1, 'Type','Axes', ...
                    '-not', 'Tag','Colorbar', '-not', 'Tag','legend');
                for l = 1:length(peer)
                    hin  = double(getappdata(peer(l), ...
                        'LegendColorbarInnerList'));
                    hout = double(getappdata(peer(l), ...
                        'LegendColorbarOuterList'));
                    if any(H(cH)==[hin hout])
                        CBH(cH) = peer(l);
                        break
                    end
                end
                
            end
            
        else
            % Just some normal axes:
            
            if FORCE
                temp = colorbar('Peer',H(cH));
                if ~PEER
                    CBH(cH) = temp;
                else
                    CBH(cH) = H(cH);
                end
            end
        end
            
    end
    
end

% Hidden:
if UNHIDE
    set(0,'ShowHiddenHandles','off')
end

% Clears output:
CBH(~ishandle(CBH)) = [];


end

% [EOF] CBHANDLE.M by Carlos A. Vargas A.