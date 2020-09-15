%
% Take multiple Specrec structs as returned from irf_powerfft and merge them into one Specrec struct. The input structs
% may have different sets of frequencies, but should not overlap in time.
%
% NOTE: This function is A BIT OF A HACK to make it possible to apply irf_powerfft to data with (wildly) changing
% sampling frequencies, yet call irf_spectrogram only once.
% The intended use is to
% (1) Split the data into time intervals with one approximately constant sampling frequency per time interval.
% (2) For each time interval, call for irf_powerfft.
% (3) Merge the resulting Specrec sstructs using this function.
% (4) Call irf_spectrogram once using the result from step (3).
%
% The resulting Specrec struct describes a larger spectrum, with potentially large parts set to NaN.
% It has the union of the source structs' frequencies and timestamps, plus some. The extra frequencies and timestamps
% are there to be filled with NaN so that irf_spectrogram plots correctly, and does not (due to how MATLAB's plotting
% works) "inapropriately bind together" areas of the spectrum with each other.
% 
%
% ARGUMENTS
% =========
% SpecrecCa : 1D cell array of "Specrec" structs as returned by irf_powerfft.
%             NOTE: Must only contain one sample per timestamp.
%
%
% RETURN VALUES
% =============
% Specrec : "Specrec" struct as returned by irf_powerfft. Consists of the merger of structs in SpecrecCa. p=NaN for
%           values not assigned by any argument.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-08-14.
%
function Specrec = merge_Specrec(SpecrecCa)
    % PROPOSAL: Automatic test code.
    %
    % PROBLEM: irf_spectrogram connects areas with data.
    %   PROPOSAL: Insert NaN in between.
    %       NOTE: Needs to sort Specrecs in time.
    %   PROPOSAL: Add NaN before and after each specrec.
    %   PROBLEM: Need to consider spacing between spectras.
    %
    % INCOMPLETE: Does not complement NaN between non-NaN values in the frequency coordinate.
    %   PROPOSAL: Generic function x,y-->x,y such that NaN values are replaced by nearest value, if not farther away
    %   than threeshold.
    %       CON: When function is applied, coordianates have already been merged, and it is hard to calculate the max
    %            nearest distance for replaceing NaN to use.
    
    % NOTE: Memory use could be a potential problem since internal data size should be on the same order as data set
    % zVars. Has not yet observed to be a problem though. /2020-08-16
    
    
    
    N = numel(SpecrecCa);
    assert(N >= 1, 'SpecrecCa is empty.')
    
    
    
    % Specrec which will grow with the content of other Specrecs. Its content
    % will be overwritten.
    S = SpecrecCa{1};
    
    for iS = 1:N
        S2 = SpecrecCa{iS};
        
        % ASSERTIONS
        EJ_library.assert.struct(S2, {'t', 'p', 'f'}, {})
        assert(isscalar(S2.p))
        EJ_library.assert.sizes(S2.t, [-1], S2.p{1}, [-1, -2], S2.f, [-2]);
        
            
        if iS >= 2
            S2 = pad_NaN(S2);
        
            [tf, pArray] = EJ_library.utils.merge_coordinated_arrays(NaN, {S.t, S.f}, S.p{1}, {S2.t, S2.f}, S2.p{1});
            
            S.t = tf{1};
            S.f = tf{2};
            S.p = {pArray};
        end
        
    end
    
    for iTime = 1:numel(S.t)
        %d = mode(diff(S.f(~isnan(S.p{1}(iTime, :))))) / 2;
        
        %S.p{1}(iTime, :) = EJ_library.utils.fill_NaN(S.f, S.p{1}(iTime, :), d);
        
        S.p{1}(iTime, :) = use_nearest_nonNaN(S.f, S.p{1}(iTime, :));
    end
    
    Specrec = S;
end



function Specrec = pad_NaN(Specrec)
    
    t = Specrec.t;
    % NOTE: Requires numel(t) >= 2.
    t = [2*t(1) - t(2); t; 2*t(end) - t(end-1)];
    Specrec.t = t;
    
    Specrec.p{1} = padarray(Specrec.p{1}, [1, 0], NaN, 'both');
    % NOTE: .f unchanged
end



% Replace NaN values with nearest non-NaN value, unless outside range of non-NaN
% values
function y = use_nearest_nonNaN(x,y)
    bFinite = ~isnan(y);
    x2 = x(bFinite);
    y2 = y(bFinite);
    if any(bFinite)
        y = interp1(x2, y2, x, 'nearest');   % No extrapolation.
    end
end
