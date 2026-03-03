%
% Code for applying a TF to samples in the form of a kernel. In practice, one
% main function plus helper functions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef kernel
  % PROPOSAL: Name "kern".
  %   PRO: Shorter
  %   PRO: Same length as "freq" and "time".



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Apply a transfer function in the form of a specified kernel to a series
    % of samples, i.e. convolution. The function has functionality for how to
    % treat edges of signal.
    %
    % The function can be seens a wrapper around MATLAB's conv() function.
    %
    % NOTE: The function does not have any concept of time more than element
    % indices (there is no argument "dt" etc.).
    %
    % NOTE: The terminology of the function is to apply an impulse response,
    % but in practice the function can be used for any 1D "point-spread"
    % function, i.e. both "backward" and "forward" impulse responses. The
    % function is in reality intended for reversing/inverting the effects of a
    % transfer function.
    %
    %
    % ARGUMENTS
    % =========
    % y1
    %       Real column vector. Input time series.
    %       NOTE: May contain NaN.
    % yKernel
    %       Real column vector. Kernel/impulse response.
    %       NOTE: May contain NaN.
    % iKernelOrigin
    %       Scalar index into yKernel. Determines what is to be regarded as the
    %       "coordinate origin" or center of the kernel, i.e. e.g. "x=0", or
    %       t=0 in a impulse response). For example, if the kernel is non-zero
    %       only at this index, then applying the kernel to the signal is
    %       equivalent to just scaling/multiplying the signal.
    %       NOTE: See IMPLEMENTATION NOTES.
    % edgePolicy
    %       String constant. Specifies how to handle edges.
    %       'ZEROS'
    %       'CYCLIC'
    %       'MIRROR'
    % --
    % NOTE: No argument for sampling frequency or period (dt) is required.
    %
    %
    % NOT-A-NUMBER
    % ============
    % This function is currently designed to permit NaN and allow values to
    % propagate so that it is ~analogous with bicas.tf.freq.apply_TF().
    %
    %
    % IMPLEMENTATION NOTES
    % ====================
    % IMPLEMENTATION NOTE: In practice, a user probably wants iKernelOrigin to
    % be the middle of the vector (~length/2). However, it is useful to be able
    % to specify this arbitrarily in an argument since:
    %     PRO: It makes the function more generic.
    %     PRO: Middle index is ambiguous for even-length kernels (depends on
    %          rounding). The caller must need to know the center index anyway
    %          since it creates the kernel. Thus having an argument for it
    %          prevents the code from encoding the same design decision twice
    %          (i.e. code could be inconsistent ==> bug).
    %     PRO: Having it as an explicit variable clarifies the code/algorithm
    %          in the same way as having a global constant instead of a
    %          hardcoded value. In this particular case(?), the hardcoded
    %          value would probably have been further obscured by simplifying
    %          expressions with it.
    %     CON: Requires more test cases.
    %     CON: One more argument.
    %
    %
    % RETURN VALUES
    % =============
    % y2
    %       Real column vector.
    %
    %
    % Author: Erik P G Johansson, Uppsala, Sweden
    % First created 2021-08-08.
    %
    function y2 = apply_kernel(y1, yKernel, iKernelOrigin, edgePolicy)
      % TODO-NI: Speed up?
      %   PROPOSAL: Permit multiple input & output signals.
      %     Ex: Snapshots.
      %   PROPOSAL: Do not use plain conv(). Use "overlap save method" as in
      %             c_efw_invert_tf.m:block_conv()?
      %
      % PROPOSAL: Upon detection of NaN in signal or kernel, make entire output
      %           signal NaN.
      %   PRO: More consistent with bicas.tf.freq.apply_TF().
      %   CON: Bad for CWF.
      %   CON-PROPOSAL: Do in wrapper.
      %
      % TODO-DEC: Does argument iKernelOrigin make sense?
      %   PROPOSAL: Always assume iKernelOrigin=1.
      %   PROPOSAL: Assume iKernelOrigin = lengthKernel/2
      %     PROPOSAL: Require odd-length kernel where iKernelOrigin is always
      %               the middle element.
      %     PROPOSAL: Support both odd and even-length kernels.
      % PROPOSAL: Internal implementation pads yKernel to be symmetric
      %           (odd-valued length; iKernelOrigin in the middle).
      % ~BUG: Imposes kernel length limit, but the real conceptual motivation
      %       for it is that it should be possible to pad the signal with the
      %       side of the kernel before/after the kernel origin (the "radius"
      %       of the kernel), not the kernel length in total. The kernel could
      %       thus be longer.

      lenKernel = length(yKernel);

      %============
      % ASSERTIONS
      %============
      bicas.tf.kernel.assert_y(y1,      'y1')
      bicas.tf.kernel.assert_y(yKernel, 'yKernel')
      %
      % IMPLEMENTATION NOTE: Requiring 1 <= lenKernel <= length(y1), but in
      % principle, the upper limit is not entirely fundamental. However, (1)
      % imposing this requirement is "natural", and (2) using the *same*
      % requirement globally (independent of edgePolicy) simplifies the
      % interface (makes it easier to understand).
      %
      % The "real" constraints on kernel length really refer to number of
      % samples before and after iKernelOrigin:
      % edgePolicy=ZEROS
      %       There is no conceptual nor implementation constraint on max
      %       kernel length.
      % edgePolicy=CYCLIC
      %       There is no conceptual constraint on max kernel length.
      %       Implementation (padding) is a constraint.
      % edgePolicy=MIRROR
      %       There is a conceptual constraint on max kernel length, maybe.
      %       Implementation (padding) is a constraint.
      assert(~isempty(yKernel), ...
        'BICAS:Assertion:IllegalArgument', ...   % NOTE: EMID is required in test.
        'Argument yKernel is empty.')
      assert(lenKernel <= length(y1))
      %
      assert(isscalar(iKernelOrigin) & isnumeric(iKernelOrigin))
      assert((1 <= iKernelOrigin) & (iKernelOrigin <= lenKernel))



      %-----------------------------------------------------
      % Lengths of minimum necessary padding before & after
      %-----------------------------------------------------
      % Padding length BEFORE signal == Length of kernel AFTER  origin.
      nPadA = lenKernel - iKernelOrigin;
      % Padding length AFTER  signal == Length of kernel BEFORE origin.
      nPadB = iKernelOrigin - 1;

      %====================================
      % Pad signal y1 depending on setting
      %====================================
      y1padded = bicas.tf.kernel.add_padding(y1, nPadA, nPadB, edgePolicy);

      %=====================================================
      % CONVOLVE PADDED SIGNAL USING MATLAB FUNCTION conv()
      %=====================================================
      y2padded = conv(y1padded, yKernel);

      %================
      % Remove padding
      %================
      i0 = nPadA + iKernelOrigin - 1;    % Index before first index to keep.
      y2 = y2padded(i0 + [1:length(y1)]);
    end



    function y2 = add_padding(y1, nPadA, nPadB, edgePolicy)
      % PROPOSAL: Add tests.
      %   NOTE: Currently only tested indirectly.

      lenY1 = length(y1);

      switch(edgePolicy)
        case 'ZEROS'
          %================
          % Pad with zeros
          %================
          % NOTE: Due to how conv() and the algorithm works, padding with
          % zeros is equivalent to not padding at all, IF THERE ARE NO
          % NOT-A-NUMBER in the signal. Therefore padding with zeros anyway.

          yPadA = zeros(nPadA, 1);
          yPadB = zeros(nPadB, 1);

        case 'CYCLIC'
          %==============================================
          % Pad with signal itself, as if it were cyclic
          %==============================================
          % NOTE: This mode is implemented to make it possible to get the exact
          % same result with bicas.tf.time.apply_TF() as with
          % bicas.tf.freq.apply_TF().

          % IMPLEMENTATION NOTE: Could (?) be implemented with MATLAB's
          % cconv(), but that would defeat the purpose of having this case for
          % testing (to test other code).

          % ASSERTION
          % NOTE: Could in principle update implementation to eliminate this
          % constraint.
          assert(max(nPadA, nPadB) <= lenY1,...
            ['CYCLIC: Can not pad with more signal samples', ...
            ' than thera are samples available.'])

          yPadA = y1(end-nPadA+1 : end,   1);
          yPadB = y1(1           : nPadB, 1);

        case 'MIRROR'
          %=============================================================
          % Pad edges with mirrored signals (mirrored around the edges)
          %=============================================================
          % NOTE: The implementation uses mirror symmetry axes located at
          % "array indices" 0.5 and end+0.5, i.e. the very first and last
          % samples are themselves mirrored (duplicated).
          % NOTE: One could also use symmetry around indices "1" and "end" and
          % not mirror the very first and last samples.

          % ASSERTION
          assert(max(nPadA, nPadB) <= lenY1,...
            ['MIRROR: Can not pad with more mirrored signal samples', ...
            ' than thera are samples available.'])

          % NOTE: Y = wrev(X) reverses the 1D vector X.
          yPadA = wrev(y1(1           : nPadA, 1));
          yPadB = wrev(y1(end-nPadB+1 : end,   1));

        otherwise
          error('Illegal argument edgePolicy="%s".', edgePolicy)

      end

      % Pad signal.
      y2 = [yPadA; y1; yPadB];
    end



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function assert_y(y, argName)
      assert(iscolumn(y),  'Argument %s is not a column vector.', argName)
      assert(isnumeric(y), 'Argument %s is not numeric.', argName)
      assert(isreal(y),    'Argument %s is not real.', argName)
    end



  end    % methods(Static, Access=private)



end
