%
% Create filename to be used for OFFICIAL summary plot.
%
%
% NOTE: For exact filenaming convention, see solo.sp.create_SP_basename().
%
%
% ARGUMENTS
% =========
% versionNbr
%       Should be dataset version number (not plot file version).
%
%
% RETURN VALUES
% =============
% filename
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-08-13.
%
function filename = create_SP_filename(srcDatasetId, dateVec3, versionNbr)
    % YK 2020-0x-xx: Use "png" as image format.
    %
    % PROPOSAL: Merge with solo.sp.create_SP_basename().
    
    spBasename = solo.sp.create_SP_basename(...
        srcDatasetId, dateVec3, versionNbr);
    filename   = sprintf('%s.png', spBasename);
end
