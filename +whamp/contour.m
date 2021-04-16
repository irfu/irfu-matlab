function contour(p,z,f)
% WHAMP.CONTOUR(p,z,f)
%
% plots contours with some options that are often used
%
q1 = input('Contour or Surface, c/s: ','s');
switch lower(q1)
  case 'c'
    % Contour
    fToPlot = f; % Default f to be plotted
    disp('Levels')
    disp('1 Automaticaly linear')
    disp('2 Automaticaly log')
    disp('3 [0 0.1 0.2 ... 1.0]')
    disp('4 [0 0.01 0.02 0.03 0.04 0.05 0.06]')
    disp('5 [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1]')
    disp('9 Specify levels')
    q2 = input('Your selection, [1-9]: ');
    switch q2
      case {1, 2}
        if q2==2
          fToPlot = log10(f); % f to be plotted should be log, per user input.
        end
        q2 = input('Number of levels?');
        v = q2;
      case 3
        v = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
      case 4
        v = [0 0.01 0.02 0.03 0.04 0.05 0.06];
      case 5
        v = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1];
      case 9
        q2 = input('level1 level2 ...=','s');
        eval(['v=[' q2 ']']);
      otherwise
        % Unknown input, return error
        error('Unkown input, levels should be a number "1", "2", "3", "4", "5" or "9".');
    end
    cs = contour(p, z, fToPlot, v); % Now plot as per user input
    if strcmpi(input('Label whamp.contours? y/n','s'), 'y')
      clabel(cs);
    end
    xlabel('k_{perp}');
    ylabel('k_{par}');
  case 's'
    % Surface
    pcolor(p,z,f);
    colorbar;
    q2 = input('Color axis (if nothing automatically) cmin cmax =','s');
    if (~isempty(q2))
      caxis(eval(['[' q2 ']']));
    else
      caxis('auto')
    end
    colorbar;
    q2 = 'flat'; % set default..
    q2 = irf_ask('shading (flat,faceted,interp,no) [%] >', 'q2','flat'); % NOTE: This function is part of irfu-matlab
    if ismember(lower(q2), {'flat', 'faceted', 'interp', 'no'})
      if ~strcmpi(q2, 'no')
        eval(['shading ' lower(q2)]);
      end
    else
      warning('Unexpected shading repsonse ignored, should be one of "flat", "faceted", "interp" or "no".');
    end
  otherwise
    % Unknown input, return error
    error('Unknown input, should be "c" or "s" for contour or surface respectivly.');
end % end switch lower(q1)

end