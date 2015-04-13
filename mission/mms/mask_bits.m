function data = mask_bits(data, bitmask, bit)
%MASK_BIT  mask a particular bit (set to NaN)
%
% DATA = MASK_BIT(DATA, BITMASK, BIT)
%
% Set DATA to NaN where BIT is set in BITMASK. This setting is "in place",
% i.e. the input data is being modified.

narginchk(3,3)

if ~all( size(data) == size(bitmask))
  error('DATA and BITMASK')
end

data(logical(bitand(bitmask, bit))) = NaN;

