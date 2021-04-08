%
% Derive official filename for one VHT dataset file.
%
% NOTE: Technically uses hardcoded DATASET_ID.
%
%
% ARGUMENTS
% =========
%
%
%
% RETURN VALUES
% =============
% filename
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2021-03-30.
%
function filename = derive_dataset_filename(yearMonth, verNbr)
    % PROPOSAL: Automatic test code.
    %
    % Cf. "solo_L1_rpw-bia-current-cdag_20200401-20200430_V04.cdf"
    
    assert((length(yearMonth) == 2) && isnumeric(yearMonth))
    assert(isscalar(verNbr) && isnumeric(verNbr))
    
    % Date vectors for first and last DAY of month (NOT exact timestamp, i.e.
    % NOT midnight).
    dv1 = datevec(datetime([yearMonth(1), yearMonth(2),   1]));
    dv2 = datevec(datetime([yearMonth(1), yearMonth(2)+1, 0]));
    
    yearMonthStr = sprintf('%04i%02i%02i-%04i%02i%02i', dv1(1:3), dv2(1:3));
    
    filename = sprintf('solo_L3_rpw-bia-vht_%s_V%02i.cdf', yearMonthStr, verNbr);
end