% Speed test of dictionary lookup which sometimes seems suspiciously slow for
% SSID/SDID/ASID dictionaries and large arrays.
%
%
function dict_lookup___speedTest()
% Size [843 2048 5] has been observed in actual processing.
SIZE = [843 2048 5];

ssidAr          = repmat(uint8(101), SIZE);
ssidAr(1:2:end) = uint8(102);

assert(isequal(size(ssidAr), SIZE))

tic
valueAr1 = bicas.utils.dict_lookup(bicas.proc.L1L2.const.C.SSID_ASID_DICT, ssidAr);
toc

tic
valueAr2 = bicas.proc.L1L2.const.C.SSID_ASID_DICT(ssidAr);
toc

assert(isequaln(valueAr1, valueAr2))
end