%
% Autocreate ONE BPCI from DSMDs, including default output filenames while only
% considering one SWM.
%
%
% NOTE: Not obvious how to derive output file time range automatically.
%   ** If one uses the UNION of ALL input dataset time intervals, then the
%      interval becomes very large when current datasets are used.
%   ** If one uses the INTERSECTION of ALL input dataset time intervals, then
%      the interval can become negative if the CURRENT DSMD was derived using
%      zVar Epoch and begins after the other datasets.
%
%
% ARGUMENTS
% =========
% Swm
%       bicas.swm.SoftwareMode object.
% DsmdArray
%       DSMD array. Must contain exactly one DSI for each input DSI required by
%       the SWM.
% createOutputPathFh
%       Function for obtaining the path to one output dataset.
%           path = func(outputDsi, InputDsmdArray)
%       Arguments describe the SWM inputs (all) & output (one).
%       NOTE: bicas.tools.batch.default_get_BPCI_output_filename() is
%       meant to be used for helping to construct such functions by default.
%       That function does however only construct the filename (not specify the
%       parent directory).
%
%
% RETURN VALUE
% ============
% Bpci
%       bicas.tools.batch.BicasProcessingCallInfo object.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-05-27.
%
function Bpci = autocreate_one_SWM_BPCI(Swm, DsmdArray, createOutputPathFh)

% ASSERTIONS
assert(isa(Swm, 'bicas.swm.SoftwareMode') && isscalar(Swm))
irf.assert.castring_set({Swm.inputsList.dsi })
irf.assert.castring_set({Swm.outputsList.dsi})
assert(isa(createOutputPathFh, 'function_handle'))
assert(numel(Swm.inputsList) == numel(DsmdArray))

% t = tic();

% cohbCa      = cell(0, 1);
inputsArray = bicas.tools.batch.BpciInput.empty(0, 1);
for i = 1:numel(Swm.inputsList)
  % Iterate over inputs
  % ===================
  % NOTE: Better to iterate over SWM inputs than DSMDs, since
  % (1) Easier to permit irrelevant DSMDs (not implemented)
  % (2) Check that there is exactly one matching DSMD (by DSI).

  dsiToSearchFor = Swm.inputsList(i).dsi;
  iDsmd          = find(strcmp(dsiToSearchFor, {DsmdArray.datasetId}));

  % ASSERTIONS
  if numel(iDsmd) == 0
    error('Can not find any dataset with DSI="%s".', dsiToSearchFor)
  elseif numel(iDsmd) >= 2
    error('Found more than one dataset with DSI="%s".', dsiToSearchFor)
  end

  inputsArray(i, 1) = bicas.tools.batch.BpciInput(...
    Swm.inputsList(i).cliOptionHeaderBody, ...
    Swm.inputsList(i).dsi, ...
    DsmdArray(iDsmd).path);

  % Prepare for calling output filenaming function.
  % cohbCa{i}     = Swm.inputsList(iDsmd).cliOptionHeaderBody;
  iDsmdArray(i) = iDsmd;
end

% Re-arrange order so that it is equal to that of cohbCa.
% DsmdArray = DsmdArray(iDsmdArray);

%fprintf('SPEED: autocreate_one_SWM_BPCI(): t=%.1f [s]\n', toc(t));

outputsArray = bicas.tools.batch.BpciOutput.empty(0, 1);
for i = 1:numel(Swm.outputsList)
  filePath = createOutputPathFh(Swm.outputsList(i).dsi, DsmdArray);
  outputsArray(i, 1) = bicas.tools.batch.BpciOutput(...
    Swm.outputsList(i).cliOptionHeaderBody, ...
    Swm.outputsList(i).dsi, ...
    filePath...
    );
end
%fprintf('SPEED: autocreate_one_SWM_BPCI(): t=%.1f [s]\n', toc(t));

Bpci = bicas.tools.batch.BicasProcessingCallInfo(Swm.cliOption, inputsArray, outputsArray);
end
