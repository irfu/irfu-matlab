function bitmask = mms_sdc_sdp_bitmasking( dceTimes, dcvTimes )
% MMS_SDC_SDP_BITMASKING compares time intervall of the two time series and 
% returns status code bitmask indicating overlap or not.
%   First column of bitmask corresponds to timestamps of dceTimes and
%   the second column indicate if it was overlapping or not. If the times 
%   do not match bitmask value MMS_CONST.Bitmask.OnlyDCE will be retured
%   for those times.
%
%	MMS_SDC_SDP_BITMASKING( dceTimes ) return a bitmask indicating 
%   only a single source for all timestamps and no overlap.
%
%   MMS_SDC_SDP_BITMASKING( dceTimes, dceTimes ) return a bitmask
%   indication full overlap, i.e. "0". (i.e. times are identical).
%
%	Example:
%		bitmask = mms_sdc_sdp_bitmasking(dceTimes);
%		bitmask = mms_sdc_sdp_bitmasking(dceTimes, dcvTimes);
%
% 	See also MMS_SDC_SDP_DATAMANAGER.

global MMS_CONST;

narginchk(1,2);


bitmask = zeros(length(dceTimes),2);
% Bitmask should cover the entire time interval.
bitmask(:,1) = dceTimes;

if(nargin==1)
    % If only one source file mark all times with having only one source.
    irf.log('debug','Only one input to bitmasking. All of bitmask is OnlyDCE');
    bitmask(:,2) = MMS_CONST.Bitmask.OnlyDCE;

elseif(nargin == 2)
    
    irf.log('debug','Two inputs to bitmasking. Compare both time series and set OnlyDCE on relevant times.');
    % Check if second data source begins after first and if so set bitmask
    % to only DCE for all of these data packet times.
    j = 1;
    while((dceTimes(j) < dcvTimes(1)) && (j <= length(dceTimes)))
        % Our src1 (given priority begin before src2)
        bitmask(j,2) = MMS_CONST.Bitmask.OnlyDCE;
        j = j + 1;
    end

    % Check if second data source ends before first and if so set bitmask
    % to only DCE for all of these data packet times.
    if( dceTimes(end) > dcvTimes(end) )
        j = length(dceTimes);
        while(dceTimes(j)>dcvTimes(end))
            bitmask(j,2) = MMS_CONST.Bitmask.OnlyDCE;
            j = j - 1;
        end
    else
        % Second source ends at same time or later than first source. Use
        % times for first source until it ends.
    end
    
end

