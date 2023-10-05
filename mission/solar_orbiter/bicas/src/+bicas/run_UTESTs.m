%
% Run all of BICAS's automatic tests. This function is not called by BICAS
% proper. It is only intended to be called manually during development to easily
% check that changes have not obviously broken anything.
%
% NOTE: Will ALWAYS only run in the CURRENT irfu-matlab git repo, not
% necessarily the irfu-matlab git repo used for running BICAS.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2021-08-19.
%
function run_UTESTs()
    tic
    runtests('bicas', 'IncludeSubpackages', true);
    toc
end
