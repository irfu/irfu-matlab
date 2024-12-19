%
% Functions for working with mixed radix numbers.
%
% Ex: Can interpret [dayNbr, hours, minutes, second] as digits in mixed radix
%     number (if not using leap seconds).
% Ex: Can interpret [seconds, milliseconds, microseconds, nanoseconds] as used
%     by the spdf* TT2000 functions.
%
% See e.g. https://en.wikipedia.org/wiki/Mixed_radix .
%
%
% NOTES
% =====
% Code is currently very permissive in the accepted input values, as long as the
% algorithm is a well defined and consistent extrapolation of the canonical
% behaviour. Therefore:
% () Integer as an input argument may be negative.
% () MRD as an input argument may go outside range 0 to (base-1).
% () Base>=1 is permitted. Digits have the same value as the next "less
%    significant" digit (which technically is equally significant). Could be
%    extended to permit base=<0.
%
%
% CONVENTIONS
% ===========
% MRD = Mixed Radix Digits, i.e. the digits in a mixed radix number.
% --
% iDigit=1 <==> Least significant digit.
%            NOTE: This means that digits are in the opposite order to what is
%            convention when printing out digits as a literal at the MTLAB
%            command-line. It is also the opposite convention compared to
%            MATLAB's date vectors. Therefore often necessary to use "fliplr()".
% iNbr     : Index to separate number, to be converted/treated separately.
% --
% NOTE: Arguments and return values follow (almost) same conventions, including
%       sizes.
% bases : (iDigit).
% mrd   : (iNbr, iDigit).
%         NOTE: This implies that an mrd which represents a single integer is a
%               ROW vector, not column.
% n     : (iNbr). May be negative as input argument.
% --
% NOTE: All arguments will be rounded to integers (int64) internally. The caller
% must make sure to not go outside of range used by internal data type.
% NOTE: Return values are always int64.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-10-21
%
classdef mixed_radix
  % PROPOSAL: Shorter class name.
  % PROPOSAL: Move outside of utils package.

  % PROPOSAL: Permit infinitely high highest base.
  %   CON: Caller can simply set the highest base to a high value.
  %   CON: Ill-defined result for negative integers.
  % PROPOSAL: Permit negative integers (not MRD).
  %   CON: Caller can do mod before calling instead.
  %       CON: Extra work for caller.
  %
  % PROPOSAL: "mrd" should be column vector (not row).
  %   PRO: Consistency with "bases", general convention.
  %   CON: Only for mrd that represents single integer. Size is
  %        (nNbrs,nDigits). Is otherwise not consistent with other "1D vector
  %        of 1D vector" matrices(?)
  %       Ex: datevec
  %
  % PROPOSAL: Optional input argument assertions.
  %   Base >=1
  %   Base >=2
  %   N >= 0
  %   MRD >= 0
  %   MRD <= base-1
  %   NOTE: All of these can be easily implemented by caller, except
  %         MRD<=base-1 for multiple integers (which is probably still a
  %         one-liner).
  %   PROPOSAL: Policy argument.
  %       TODO-DEC: Format?
  %   PROPOSAL: irf.utils.interpret_settings_args().
  %       'minBase' : Number: -Inf, 1, 2,
  %       'minN'    : Number: -Inf, 0
  %       'minMrd'  : Number: -Inf, 0
  %       'maxMrd'  : 'base cap' (base-1), Inf (number?)
  % PROPOSAL: Implement strict assertions until relaxation has prove useful.



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % NOTE: Permits n to be outside of range defined by bases. Result can be
    % interpreted as using argument n := mod(n, prod(bases)).
    % RATIONALE: Useful for in particular negative values where one can not
    % emulate the behaviour by setting a very large highest base.
    %   Ex: Negative day numbers+time of day.
    %
    function mrd = integer_to_MRD(n, bases)
      % ASSERTIONS
      [nNbrs, nBases] = irf.assert.sizes(...
        n,     -1, ...
        bases, -2);
      assert(all(bases >= 1))

      n     = int64(n);
      bases = int64(bases);

      mrd = int64(zeros(nNbrs, nBases));
      B   = int64(1);
      for i = 1:nBases
        mrd(:, i) = mod(idivide(int64(n), B, 'floor'), bases(i));
        B = B*bases(i);
      end

      % ASSERTION: Check argument (sic!).
      % Done last since requires B.
      %assert(all(n<B), 'n is too large for the specified bases.')
    end



    % NOTE: Permits mrd components to be outside of range as defined by
    % bases, including negative numbers.
    % RATIONALE: This is useful for in particular the highest component
    % where it is often useful to think of the highest base as
    % "infinite", or for negative values which can not be emulated with a
    % very high highest base.
    %   Ex: Day numbers+time of day.
    %
    function n = MRD_to_integer(mrd, bases)

      % ASSERTIONS
      [nNbrs, nBases] = irf.assert.sizes(...
        mrd,   [-1, -2], ...
        bases, -2);
      assert(all(bases >= 1))
      %             % NOTE: Does not check for negative MRDs.
      %             for i = 1:nBases
      %                 assert(all(mrd(:, i) < bases(i)), ...
      %                     'Mixed radix digits must be lower than digit base.')
      %             end

      mrd   = int64(mrd);
      bases = int64(bases);

      n = int64(zeros(nNbrs, 1));
      B = int64(1);
      for i = 1:nBases
        n = n + mrd(:, i) * B;
        B = B * bases(i);
      end
    end



  end



end
