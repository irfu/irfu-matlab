%
% Class which instances encapsulate and store one IMMUTABLE copy of a "Laplace/Fourier transform" on the format
%
%     nc(1)*s^0 + ... + nc(nNc)*s^(nNc-1)
% Z = -----------------------------------,   s = i*omega
%     dc(1)*s^0 + ... + dc(nDc)*s^(nDc-1)
% .
%
%
% IMPLEMENTATION NOTE
% ===================
% This class does NOT implement assertions for all conceivable sanity checks on transforms. Reasons:
% ** It is ambiguous which sanity checks one should use.
% ** Want to keep the code generic, e.g.
%    ** be able to store "backward" TFs,
%    ** be useful for numeric experiments,
%    ** be able to use it for both Laplace transforms and Fourier transforms(?).
% Instead it provides extra methods to make relevant assertions easy to implement by the user.
%
% This class deliberately requires (via assertions):
%   ** The TF can always be evaluated to finite values (not NaN, Inf) by requiring sufficient non-zero coefficients.
% This class deliberately permits the following overlapping categories of TFs:
%   ** Non-real coefficients
%      (only real coefficients <=> real impulse response <=> real input signal produces real output signal)
%   ** "Backwards" transfer functions (<=> "non-causal" TFs <=> non-zero impulse response before t=0 <=>? unstable TFs)
%   ** TFs with poles with non-negative real value (oscillations and divergent impulse responses).
%   ** Z does not converge to zero when omega goes to infinity.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-11-12
%
classdef rational_func_transform

  % PROPOSAL: Implement subclasses with more strict requirements on permitted transforms.
  % PROPOSAL: Change name to specify subset of analytical transforms.
  %       "analytical_rational_transform".
  %       "rational_func_transform"
  %           NOTE: Implicit that has factor i^k for every coefficient.
  % PROPOSAL: Capitalize initial in class name.



  % PRIVATE, IMMUTABLE
  properties(SetAccess=immutable, GetAccess=public)
    numeratorCoeffs
    denominatorCoeffs
  end



  methods(Access=public)



    function obj = rational_func_transform(numeratorCoeffs, denominatorCoeffs)
      % IMPLEMENTATION NOTE: See top-level file comments.

      % ASSERTION: Require 1D vectors
      irf.assert.vector(numeratorCoeffs)
      irf.assert.vector(denominatorCoeffs)
      % ASSERTION: Exclude NaN, +-Inf
      assert(all(isfinite(  numeratorCoeffs)))
      assert(all(isfinite(denominatorCoeffs)))
      % ASSERTION: Require at least non-zero value in each polynomial separately.
      % PROPOSITION: Bad, non-generic assertions?
      assert(max(abs(  numeratorCoeffs)) > 0)
      assert(max(abs(denominatorCoeffs)) > 0)

      % ENFORCE: Column vectors
      % NOTE: Column vectors are required for using "flipud" at evaluation.
      obj.numeratorCoeffs   = numeratorCoeffs(:);
      obj.denominatorCoeffs = denominatorCoeffs(:);
    end



    function Tf = inverse(obj)
      Tf = irf.utils.rational_func_transform(obj.denominatorCoeffs, obj.numeratorCoeffs);
    end



    % Check if transfer function converges toward zero in the limit of high frequencies.
    % This is a basic sanity check and should be a necessary(?) (but not sufficient?) criterion for stable
    % transfer functions.
    %
    % Useful for implementing separate assertions.
    %
    function zeroInHfLimit = zero_in_high_freq_limit(obj)
      % PROPOSAL: Better name
      %     high_freq_toward_zero
      %     Z_zero_high_freq
      %     Z_toward_zero_as_omega_toward_inf
      %     amplitude_toward_zero_as_freq_toward_inf
      %     Z_zero_omega_inf
      %     Z_zero_at_omega_inf
      %     Z_zero_at_high_freq
      %     Z_zero_at_inf_freq
      %     smaller_at_higher_freq
      %     decreasing_at_higher_freq
      %     high_freq_limit_zero
      %     zero_at_high_freq_limit

      % Denominator polynomial should have at least as high a degree as the numerator polynomial.
      % <==> Z should not diverge as omega-->inf.
      % Detect highest-order non-zero coefficients.
      nNc = find(obj.numeratorCoeffs,   1, 'last');    % NC = Numerator Coefficients
      nDc = find(obj.denominatorCoeffs, 1, 'last');    % DC = Denominator Coefficients
      zeroInHfLimit = (nDc > nNc);
    end



    % Return true if-and-only-if all coefficients are real.
    %
    % This criterion is equivalent with (or at least sufficient for)
    % (1) the corresponding impulse response is real, and
    % (2) that a real input signal fed through this transfer function produces a real output.
    %
    % Useful for implementing separate assertions.
    %
    %
    % NOTES
    % =====
    % x(t), y(t) are real
    % <=> X(omega)=X(*(-omega), Y(omega)=Y(*(-omega)
    % <=> H(omega)=H*(-omega) (<=> h(t) is real)
    % <== All coefficients are real (likely equivalent).
    %
    function isReal = is_real(obj)
      % NOTE: "isreal" returns a scalar value for matrices.
      isReal = isreal(obj.numeratorCoeffs) && isreal(obj.denominatorCoeffs);
    end

    function hasRealImpulseResponse = has_real_impulse_response(obj)
      hasRealImpulseResponse = obj.is_real();
    end



    % Evaluate the transfer function for specified frequencies.
    %
    function Z = eval(obj, omegaRps)
      %assert(iscolumn(obj.numeratorCoeffs  ))
      %assert(iscolumn(obj.denominatorCoeffs))

      % Calculate Z
      % NOTE: MATLAB's "polyval" uses the opposite convention for the order/index of coefficients.
      % NOTE: Column coefficient vectors are required for using "flipud".
      s = 1i * omegaRps;
      Z = polyval(flipud(obj.numeratorCoeffs), s) ./ polyval(flipud(obj.denominatorCoeffs), s);
    end



  end

end
