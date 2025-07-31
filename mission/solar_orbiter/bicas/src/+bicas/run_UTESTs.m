%
% Run all of BICAS's automatic tests. This function is not called by BICAS
% proper. It is only intended to be called manually during development to
% trigger all automated tests related to BICAS.
%
% NOTE: Will ALWAYS only run in the CURRENT irfu-matlab git repo, not
% necessarily the irfu-matlab git repo used for running BICAS.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-19.
%
function run_UTESTs()
tic

% IMPLEMENTATION NOTE: runtests('bicas', 'IncludeSubpackages', true); used to
% work, but stopped working after creating a file src/bicas.m, presumably
% because it clashes with src/+bicas/. Must therefore use another command.
runtests(bicas.utils.get_BICAS_source_dir(), 'IncludeSubfolders', true);

toc
end
