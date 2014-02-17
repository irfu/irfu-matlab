function bitmask = mms_bitmasking( sourceCDFobj1, sourceCDFobj2 )
% MMS_BITMASKING compares time intervall of the two dataobj and returns status code bitmask indicating overlap or not.
%	bitmask=MMS_BITMASKING( sourceCDFobj1, sourceCDFobj2 ) return a bitmask array of same length as sourceCDFobj1's 
%	number of datapoints. First column of bitmask corresponds to timestamps in sourceCDFobj1 and the second column 
%	indicate if sourceCDFobj1's time stamp is within the range of sourceCDFobj2s time.
%	If the times do not match bitmask value MMS_CONST.Bitmask.OnlyDCE will be retured.
%
%	MMS_BITMASKING( sourceCDFobj1 ) return a bitmask indicating only a single source for all timestamps and no overlap.
%
%       MMS_BITMASKING( sourceCDFobj1, sourceCDFobj1 ) return a bitmask indication full overlap. (i.e. times are identical).
%
%	Example:
%		bitmask = mms_bitmasking(sourceCDFobj1);
%		bitmask = mms_bitmasking(sourceCDFobj1, sourceCDFobj2);
%
% 	See also DATAOBJ, MMS_INIT.

global MMS_CONST;

narginchk(1,2);

% Get times for first priority, (tt2000 times are at varTTsrc1.data(1:end))
varTTsrc1 = getv(sourceCDFobj1, sourceCDFobj1.vars{1,1});

bitmask = zeros(varTTsrc1.nrec,2);
% Bitmask should cover the entire time of src1
bitmask(:,1) = varTTsrc1.data;

if(nargin==1)
    % If only one source file mark all times with having only one source.
    bitmask(:,2) = MMS_CONST.Bitmask.OnlyDCE;

elseif(nargin == 2)
    % Get times for source 2.
    varTTsrc2 = getv(sourceCDFobj2, sourceCDFobj2.vars{1,1});
    
    % Check if second data source begins after first and if so set bitmask
    % to only DCE for all of these data packet times.
    j = 1;
    while((varTTsrc1.data(j) < varTTsrc2.data(1)) && (j <= varTTsrc1.nrec))
        % Our src1 (given priority begin before src2)
        bitmask(j,2) = MMS_CONST.Bitmask.OnlyDCE;
        j = j + 1;
    end

    % Check if second data source ends before first and if so set bitmask
    % to only DCE for all of these data packet times.
    if(varTTsrc1.data(end)>varTTsrc2.data(end))
        j = varTTsrc1.nrec;
        while(varTTsrc1.data(j)>varTTsrc2.data(end))
            bitmask(j,2) = MMS_CONST.Bitmask.OnlyDCE;
            j = j - 1;
        end
    else
        % Second source ends at same time or later than first source. Use
        % times for first source until it ends.
    end
    
end

