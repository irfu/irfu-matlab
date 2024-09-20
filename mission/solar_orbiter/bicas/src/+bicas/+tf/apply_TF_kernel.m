%
% Apply a transfer function in the form of a kernel to a time series, i.e.
% convolution. Has functionality for how to treat edges of signal.
%
% Can be seens a wrapper around MATLAB's conv() function.
%
% NOTE: The terminology of the function is to apply an impulse response, but
% in practice the function can be used for any 1D "point-spread" function, i.e.
% both "backward" and "forward" impulse responses. The function is in reality
% intended for reversing/inverting the effects of a transfer function.
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
%       "coordinate origin" or center of the kernel, i.e. e.g. "x=0", or t=0 in
%       a impulse response). For example, if the kernel is non-zero only at this
%       index, then applying the kernel to the signal is equivalent to just
%       scaling/multiplying the signal.
%       NOTE: See IMPLEMENTATION NOTES.
% edgePolicy
%       String constant. Specifies how to handle edges.
%           NOTE: Edges are never handled as if the signal was cyclic.
%       'ZEROS'
%       'CYCLIC'
%       'MIRROR'
%           NOTE: Kernel length determines the amount of padding, which is taken
%           from the signal, which implies that there is a smallest number of
%           signal samples required.
%           length(yKernel)/2 > length(y1) (approx.) ==> Error
% --
% NOTE: No argument for sampling frequency is required.
%
%
% NOT-A-NUMBER
% ============
% Function is currently designed to permit NaN and allow values to propagate so
% that it is ~analogous with bicas.tf.apply_TF_freq().
%
%
% IMPLEMENTATION NOTES
% ====================
% IMPLEMENTATION NOTE: In practice, a user probably wants iKernelOrigin to be
% the middle of the vector (~length/2). However, it is useful to be able to
% specify this arbitrarily in an argument since:
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
%     CON: More arguments.
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
function y2 = apply_TF_kernel(y1, yKernel, iKernelOrigin, edgePolicy)
%
% PROPOSAL: Permit multiple input & output signals.
%   PRO: More efficient(?)
%       Ex: Snapshots.
%
% TODO-DEC: Which functionality should be placed in
%   (1) this function,
%   (2) wrapper(s)?
%   PROPOSAL: Edge handling.
%   PROPOSAL: Modify impulse response for stability.
%       Ex: Hann window.
%       NOTE: Includes modifying for scaling, offsets?
%   PROPOSAL: Interpolating kernel/change sampling frequency?!
%   PROPOSAL: De- & re-trending?!!
%   PROPOSAL: Somehow reuse the de- & re-trending now in
%             bicas.tf.apply_TF.m.
%
% PROPOSAL: Permit zero length kernel.
%   CON: iKernelOrigin is then ambiguous.
%       CON-PROPOSAL: iKernelOrigin does not conceptually have to be within
%                     the index range of the actual kernel variable. It
%                     could be outside too.
%   PROPOSAL: Implement by normalizing empty kernel to [0].
%       NOTE: conv() seems inconsistent for empty vectors.
%           length(conv(zeros(1,0), zeros(1,0))) == 0
%           length(conv(zeros(1,1), zeros(1,0))) == 1
%           length(conv(zeros(1,0), zeros(1,1))) == 1
%           length(conv(zeros(1,1), zeros(1,1))) == 1
%
% PROPOSAL: Option to pad with signal itself cyclically (to mimic circular
%       convolution).
%   PRO: Can be used for automatic tests that can be applied to both
%        bicas.tf.apply_TF() and bicas.tf.apply_TF_time().
%
% TODO-NI: Speed up?
%   PROPOSAL: Do not use plain conv(). Use "overlap save method" as in
%             c_efw_invert_tf.m:block_conv()?
%
% PROBLEM: Uses naming convention where 1 & 2 represent both
%   (a) signal before and after applying kernel.
%   (b) padding before and after input signal.
%   PROPOSAL: Input/output (signal), before/after (signal or kernel), A/B
%
% PROPOSAL: Upon detection of NaN in signal or kernel, make entire output
%           signal NaN.
%   PRO: More consistent with bicas.tf.apply_TF_freq().
%   CON: Bad for CWF.
%   CON-PROPOSAL: Do in wrapper.

EMID = 'BICAS:Assertion:IllegalArgument';

lenKernel = length(yKernel);
lenY1     = length(y1);

%============
% ASSERTIONS
%============
assert_y(y1,      'y1',      EMID)
assert_y(yKernel, 'yKernel', EMID)
assert(~isempty(yKernel), ...
  'BICAS:Assertion:IllegalArgument', ...
  'Argument yKernel is empty.')
assert(isscalar(iKernelOrigin) & isnumeric(iKernelOrigin))
assert((1 <= iKernelOrigin) & (iKernelOrigin <= lenKernel))



%-----------------------------------------------------
% Lengths of minimum necessary padding before & after
%-----------------------------------------------------
% Padding length BEFORE signal == Length of kernel AFTER origin.
nPad1 = lenKernel - iKernelOrigin;
% Padding length AFTER signal == Length of kernel BEFORE origin.
nPad2 = iKernelOrigin - 1;

%====================================
% Pad signal y1 depending on setting
%====================================
switch(edgePolicy)
  case 'ZEROS'
    %================
    % Pad with zeros
    %================
    % NOTE: Due to how conv() and the algorithm works, padding with
    % zeros is equivalent to not padding at all, IF THERE ARE NO
    % NOT-A-NUMBER in the signal. Therefore paddign with zeros anyway.

    yPad1 = zeros(nPad1, 1);
    yPad2 = zeros(nPad2, 1);

  case 'CYCLIC'
    %==============================================
    % Pad with signal itself, as if it were cyclic
    %==============================================
    % NOTE: This mode is implemented to make it possible to get the
    % exact same result with bicas.tf.apply_TF_time() as with
    % bicas.tf.apply_TF_freq().
    % IMPLEMENTATION NOTE: Could (?) be implemented with MATLAB's
    % cconv(), but that would defeat the purpose of having this case for
    % testing (to test other code).

    % ASSERTION
    % NOTE: Could update implementation to eliminate this constraint.
    assert(max(nPad1, nPad2) <= lenY1,...
      EMID, ...
      ['Kernel length implies padding with more mirrored signal', ...
      ' samples than thera are samples available.'])

    yPad1 = y1(end-nPad1+1 : end,   1);
    yPad2 = y1(1           : nPad2, 1);

  case 'MIRROR'
    %=============================================================
    % Pad edges with mirrored signals (mirrored around the edges)
    %=============================================================
    % NOTE: The implementation uses mirror symmetry axes located at
    % "indices" 0.5 and end+0.5, i.e. the very first and last samples
    % are mirrored (duplicated).
    % NOTE: One could also use symmetry around indices "1" and
    % "end" and not mirror the very first and last samples.

    % ASSERTION
    assert(max(nPad1, nPad2) <= lenY1,...
      EMID, ...
      ['Kernel length implies padding with more mirrored signal', ...
      ' samples than thera are samples available.'])

    % NOTE: Y = wrev(X) reverses the 1D vector X.
    yPad1 = wrev(y1(1           : nPad1, 1));
    yPad2 = wrev(y1(end-nPad2+1 : end,   1));

  otherwise
    error(EMID, 'Illegal argument edgePolicy="%s".', edgePolicy)

end

% Pad signal.
y1b   = [yPad1; y1; yPad2];

%=====================================================
% CONVOLVE PADDED SIGNAL USING MATLAB FUNCTION conv()
%=====================================================
y2b = conv(y1b, yKernel);

%================
% Remove padding
%================
y2  = y2b(nPad1 + iKernelOrigin-1 + [1:lenY1]);

end



function assert_y(y, argName, EMID)

if ~iscolumn(y)
  error(EMID, 'Argument %s is not a column vector.', argName)

elseif ~isnumeric(y)
  error(EMID, 'Argument %s is not numeric.', argName)

elseif ~isreal(y)
  error(EMID, '%s is not real.', argName)

end
end
