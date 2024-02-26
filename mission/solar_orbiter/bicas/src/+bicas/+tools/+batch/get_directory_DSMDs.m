%
% Given a path to the reference directory, return DSMD array for it.
%
%
% ARGUMENTS
% =========
% dirPathsCa
%       Cell array of paths to directories.
%
%
% RETURN VALUES
% =============
% DsmdArray
%       1D array of DSMDs to datasets under the specified directories,
%       recursively.
% OiArray
%       The combined result of the underlying calls to dir(), but only for the
%       found datasets.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function [DsmdArray, OiArray] = get_directory_DSMDs(dirPathsCa)
% PROPOSAL: Make "generic": Move to solo.adm.

    t = tic();

    % Create empty struct array of the right size.
    % IMPLEMENTATION NOTE: Calling dir() just to make sure that the code uses
    % the exact struct which it produces. This should make the code more
    % future-proof.
    OiArray = dir('~');
    OiArray = OiArray([], 1);    % Column array.

    DsmdArray = solo.adm.DSMD.empty(0,1);
    for i = 1:numel(dirPathsCa)
        % Assert that directory exists.
        % IMPLEMENTATION NOTE: dir(fullfile(dirPathsCa{i}, '**')) will NOT raise
        % error for non-existing directory.
        irf.assert.dir_exists(dirPathsCa{i})

        DirOiArray = dir(fullfile(dirPathsCa{i}, '**'));
        DirOiArray = DirOiArray(~[DirOiArray.isdir]);
        DirOiArray = DirOiArray(:);
        % CASE: DirOiArray is a column array.
        dirFilesPathsCa = arrayfun(@(Oi) (fullfile(Oi.folder, Oi.name)), DirOiArray, 'UniformOutput', false);

        [DirDsmdArray, bIsDataSetArray] = solo.adm.paths_to_DSMD_array(dirFilesPathsCa);
        DirOiArray = DirOiArray(bIsDataSetArray);

        OiArray = [...
            OiArray; ...
            DirOiArray];
        DsmdArray = [...
            DsmdArray; ...
            DirDsmdArray];
    end

    wallTimeSec = toc(t);
    fprintf('SPEED: Wall time to obtain DSMDs for %g paths: %.1f [s]\n', ...
        numel(dirPathsCa), wallTimeSec);
end
