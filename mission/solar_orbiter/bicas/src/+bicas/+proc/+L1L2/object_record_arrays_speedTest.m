%
% Performance tests to test the feasibility of different forms of potential
% future refactoring. /2024-09-19
%
% NOTE: This code is not part of the BICAS implementation, but rather serves as
% a type of documentation.
%
% NOTE: The largest L1R LFR CWF dataset found 2024-09-16
% (solo_L1R_rpw-lfr-surv-cwf-e-cdag_20210809_V18.cdf) has 22117536 CDF records.
% This can be used as an estimate for how large arrays needs to be tested for.
%
function object_record_arrays_speedTest()

N = 22117536;
f = 1;
N2 = round(N*f);

%test_ASID_arrays(N2)
test_uint8_arrays(N=N2)
end



function test_ASID_arrays(N)
N
tic
AsidArray = repmat(bicas.proc.L1L2.AntennaSignalId.C.DC_V1, N, 1);
toc
whos AsidArray

toc
AsidArray2 = AsidArray;
% Requires eq() method for non-handle objects.
bArray = [AsidArray == AsidArray2];
assert(all(bArray))
assert(isequal(size(bArray), size(AsidArray)))
toc



%===============================================================================
% NOTE: bicas.proc.L1L2.AntennaSignalId as NOT handle object.
% N =
%     22117536
% Elapsed time is 27.210939 seconds.
%   Name                  Size                Bytes  Class                              Attributes
%
%   AsidArray      22117536x1             575055936  bicas.proc.L1L2.AntennaSignalId
%
% Elapsed time is 56.935381 seconds.
% Elapsed time is 80.516625 seconds.
% ==> 548 MiB
% irony becomes noticebly slow. MATLAB hangs (temporarily?).


%===============================================================================
% NOTE: bicas.proc.L1L2.AntennaSignalId as HANDLE object without eq() method.
% NOTE: Is slower with eq() method.
% N =
%     22117536
% Elapsed time is 1.769454 seconds.
%   Name                  Size                Bytes  Class                              Attributes
%
%   AsidArray      22117536x1             176940288  bicas.proc.L1L2.AntennaSignalId
%
% Elapsed time is 1.969811 seconds.
% Elapsed time is 2.410135 seconds.
% ==> 169 MiB

end



function test_uint8_arrays(N)
N
tic
uint8Array = repmat(uint8(3), N, 1);
toc
whos uint8Array

tic
bArray = [uint8Array == 3];
assert(isequal(size(bArray), size(uint8Array)))
toc
whos bArray

% N =
%     22117536
% Elapsed time is 0.005099 seconds.
%   Name                   Size               Bytes  Class    Attributes
%
%   uint8Array      22117536x1             22117536  uint8
%
% Elapsed time is 0.015606 seconds.
%   Name               Size               Bytes  Class      Attributes
%
%   bArray      22117536x1             22117536  logical


% ==> 21 MiB

end
