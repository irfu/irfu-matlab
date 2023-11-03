%
% Class which instances store and encapsulate one IMMUTABLE transfer function in the form of tabulated values: 1D
% vectors of Z and omega values.
%
%
% IMPLEMENTATION NOTES
% ====================
% This class does NOT implement assertions for all conceivable sanity checks on transfer functions. Reasons:
% ** It is ambiguous which sanity checks one should use.
% ** Want to keep the code generic, e.g.
%    ** be able to store "backward" TFs,
%    ** useful for numeric experiments,
%    ** use it for both Laplace transforms and Fourier transforms(?).
% Instead it provides extra methods to make relevant assertions easy to implement by the user.
% --
% This class deliberately requires (via assertions):
%   ** The TF can be evaluated to finite values (not NaN, Inf) by requiring tabulated Z values to be finite.
%   ** The TF is not extrapolated beyond the tabulated omega values.
%   ** Real omega values
% This class deliberately permits the following overlapping categories of TFs:
%   ** Separately tabulated positive and negative frequencies <==> Not using relationship between Z(omega) and
%      Z(-omega), e.g. not enforcing real impulse response (<=> Z(omega) = Z*(-omega)).
%   ** "Backwards" transfer functions ("non-causal" TFs <=> non-zero impulse response before t=0; unstable TFs)
%   ** TFs with poles with non-negative real value (oscillations and divergent impulse responses).
%   ** Z does not converge to zero when omega goes to infinity (could define value
%      for omega=inf).
% --
% This class deliberately does NOT implement an "eval" function at arbitrary frequencies since it is ambiguous.
%   ** How extrapolate beyond the frequency limits of the table?
%       ** Assert that frequency is within table (forbid extrapolation)?
%       ** Extrapolate, in particular to 0 Hz?
%   ** How interpolate (linear interpolation? quadratic interpolation? splines?)
%   ** It is easy for the caller to use e.g. interp1 for interpolation.
%
%
% RATIONALE
% =========
% This class may seem to not "do much", but it is still useful since
%   ** It is a "standard struct" for tabulated TFs. Does not need to manually keep frequencies and Z together.
%   ** Documentation
%   ** It is immutable.
%   ** The constructor can initialize using two common formats.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-11-12
%
classdef tabulated_transform
  % PROBLEM: Evaluating a TF is not unambiguous.
  %   Ex: How interpolate? Linear, quadratic, spline etc. Interpolate amplitude & phase separately?!
  %   Ex: How handle values outside of table?
  %           Assertion against extrapolation
  %           Extrapolate below min. frequency: How?
  %           Extrapolate above max. frequency: How?
  % PROPOSAL: Redefine as class for only the data, not the evaluation.
  % PROPOSAL: Redefine as class for general-purpose table.
  % PROPOSAL: Implement subclasses with more strict requirements on permitted transfer functions.
  %
  % PROPOSAL: New class name(s), without "transform".
  %   PROPOSAL: tabulated_TF, rational_func_TF
  %   CON: "transform" refers to Fourier/Laplace transform.
  % PROPOSAL: Capitalize initial in class name.



  % PRIVATE, IMMUTABLE
  properties(SetAccess=immutable, GetAccess=public)
    omegaRps    % Rps = Radians Per Second
    Z
  end



  methods(Access=public)



    % ARGUMENTS
    % =========
    % SYNTAX 1: (omegaRps, Z)
    % SYNTAX 2: (omegaRps, amplitude, phaseRad)
    % --
    % omegaRps  : 1D array of frequencies (positive, real; radians/s).
    % Z         : 1D array of complex transfer function values.
    % amplitude : 1D array of amplification values (positive, real; not dB or similar).
    % phaseRad  : 1D array of phase shifts (radians).
    % --
    % NOTE: Argument arrays must have same size.
    % NOTE: There is no automatic extrapolation to e.g. 0 rad/s or infinite rad/s. The caller has to supply this.
    %
    function obj = tabulated_transform(omegaRps, varargin)
      % PROPOSAL: Assertions for isfinite, non-length zero
      % PROPOSAL: Constructor string constants for extra assertions, extrapolations
      %   "Extrapolate to omega=0"
      %   "Assert non-negative omega"
      % PROPOSAL: Use interpret_settings_args. Require key+value.
      %   NOTE: Must rename string constant.
      %       PROPOSAL: "extrapolatePositiveFreqZtoZ0"

      %             DEFAULT_SETTINGS = [];
      %             DEFAULT_SETTINGS.extrapolatePositiveFreqZtoZero = false;

      %==========================
      % Handle omegaRps argument
      %==========================
      % ASSERTIONS: omegaRps
      irf.assert.vector(omegaRps)
      assert(all(~isnan(omegaRps)))   % NOTE: Permit +-Inf, but not NaN.
      % Do not really limit the type of TFs, but still good sanity check on argument (e.g. when confusing
      % arguments).
      assert(isreal(omegaRps))
      assert(issorted(omegaRps))

      %=====================================
      % Handle TF amplitude/phase arguments
      %=====================================
      if (numel(varargin) >= 2) && isnumeric(varargin{1}) && isnumeric(varargin{2})
        %=======================================
        % CASE: Arguments for amplitude + phase
        %=======================================
        amplitude = varargin{1};
        phaseRad  = varargin{2};
        varargin  = varargin(3:end);

        % ASSERTIONS
        irf.assert.vector(amplitude)
        irf.assert.vector(phaseRad)
        irf.assert.all_equal([...
          numel(omegaRps), ...
          numel(amplitude), ...
          numel(phaseRad)])
        % NOTE: isfinite assertion implemented via assertions on Z (later).

        Z = amplitude .* exp(1i*phaseRad);
      elseif (numel(varargin) >= 1) && isnumeric(varargin{1})
        %======================
        % CASE: Argument for Z
        %======================
        Z = varargin{1};
        varargin = varargin(2:end);

        % ASSERTIONS
        irf.assert.vector(Z)
        assert( numel(omegaRps) == numel(Z) )
      else
        error('Illegal combination of arguments.')
      end

      %===========================
      % Handle settings arguments
      %===========================
      %             Settings = irf.utils.interpret_settings_args(DEFAULT_SETTINGS, varargin);
      %             if Settings.extrapolatePositiveFreqZtoZero
      %                 % IMPLEMENTATION NOTE: Require only tabulated positive frequencies, since
      %                 % (1) it is ambiguous how one would interpolate from Z(omega<0) and Z(omega>0) (and it is interpolation,
      %                 %     not extrapolation as in the settings key name string),
      %                 % (2) "eval" method would interpolate anyway.
      %                 % NOTE: Could also implement extrapolation from Z(omega<0) (alone).
      %                 % NOTE: Settings key name itself does not indicate the requirement of positive frequencies.
      %                 assert(omegaRps(1) > 0)
      %
      %                 % NOTE: Can not just use the lowest-frequency Z value for 0 Hz since it has to be real (not complex).
      %                 Z1 = Z(1);
      %                 signZ0 = sign(real(Z1));
      %                 assert(signZ0 ~= 0, 'Can not extrapolate due to ambiguity. real(Z(1)) = 0.')
      %                 Z0 = abs(Z1) * signZ0;
      %
      %                 omegaRps = [0;  omegaRps(:)];
      %                 Z        = [Z0; Z(:)       ];
      %             end



      % ASSERTIONS: Z
      assert(all(isfinite(Z)))

      %===========================
      % Assign instance variables
      %===========================
      obj.omegaRps = omegaRps;
      obj.Z        = Z;
    end



    function Tf = inverse(obj)
      Tf = irf.utils.tabulated_transform(obj.omegaRps, 1./obj.Z);
    end



    % "Sanity check" one might want to use for assertions. Determine, "approximately", whether the tabulated abs(Z)
    % decreases toward zero for high frequencies (or is zero).
    %
    % NOTE: There is no unambiguous correct algorithm for this, except when there is a defined Z(inf)=0.
    % Might fail if TF is noisy.
    %
    function towardZeroHighFreq = toward_zero_at_high_freq(obj)
      if obj.Z(end) == 0
        towardZeroHighFreq = true;
      else
        towardZeroHighFreq = abs(obj.Z(end-1)) > abs(obj.Z(end));
      end
    end



    %         % NOTE: Interpolates Z; not amplitude and phase separately.
    %         function Z = eval_linear(obj, omegaRps)
    %             % NOTE: interp1 return NaN for values outside range.
    %             Z = interp1(obj.omegaRps, obj.Z, omegaRps, 'linear');
    %
    %             if ~all(isfinite(Z))
    %                 % IMPLEMENTATION NOTE: Experience shows that it is useful to have an extended error message confirming
    %                 % that the requested frequence range is outside the tabulated one, and by how much.
    %                 errorMsg = sprintf(...
    %                     ['Can not evaluate tabulated transfer function for frequencies outside of the range of tabulated frequencies.\n', ...
    %                     'Range of frequencies for which there are tabulated Z values:\n', ...
    %                     '    min(obj.omegaRps) = %g\n', ...
    %                     '    max(obj.omegaRps) = %g\n', ...
    %                     'Range of frequencies for which evaluation (interpolation) of Z was attempted:\n', ...
    %                     '    min(omegaRps)     = %g\n', ...
    %                     '    max(omegaRps)     = %g\n'], ...
    %                     min(obj.omegaRps), ...
    %                     max(obj.omegaRps), ...
    %                     min(omegaRps), ...
    %                     max(omegaRps));
    %
    %                 error('Assertion', errorMsg)
    %             end
    %         end



  end    % methods(Access=public)




end
