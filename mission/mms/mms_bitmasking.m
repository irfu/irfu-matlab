function bitmask = mms_bitmasking( sourceCDFobj1, sourceCDFobj2 )
% Bitmasking 
% For now, bitmasking returns bitmask = MMS_CONST.Bitmask.OnlyDCE for any
% times when the sourceCDFobj1 is not covered by sourceCDFobj2.
% 
% Date of latest change: 2014/01/14

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

