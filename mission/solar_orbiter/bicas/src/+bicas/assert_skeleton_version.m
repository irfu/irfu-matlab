%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2019-08-02
%
function assert_skeleton_version(skeletonVersionStr)

EJ_library.utils.assert.castring_regexp(skeletonVersionStr, '[0-9][0-9]')
end