%
% Simple script for testing performance of function bicas.tf.time.apply_TF().
%
% CONCLUSIONS
% ===========
% (~ = "proportional to")
% lenKernel ~ lenY1    ==> tSec ~ lenY1^2
% lenKernel = constant ==> tSec ~ lenY1
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function time_apply_TF___speedTest()
% TODO-NI: What is typical long sample sequence?

close all

xAr    = [];
tSecAr = [];
for x = logspace(log10(1e2), log10(1e5), 10)
  %for x = logspace(log10(1e2), log10(1e8), 10)
  x = round(x);

  lenY1     = x;
  lenKernel = x;
  %lenKernel = min(x, 100);
  dt = 0.1;
  t  = [0:(lenY1-1)]' * dt;
  tf = @(omegaRps) (2+0.1*omegaRps);    % Has no deep interpretation.
  y1 = sin(t);                          % Has no deep interpretation.
  edgePolicy        = 'MIRROR';
  hannWindowEnabled = true;

  % Log
  fprintf('lengthY1 = %g\n', lenY1);

  xAr(end+1)    = lenY1;
  tSec          = test(dt, y1, tf, lenKernel, edgePolicy, hannWindowEnabled);
  tSecAr(end+1) = tSec;

  % Log
  fprintf('------------------------\n')
  fprintf('lengthY1          = %g\n', lenY1);
  fprintf('lenKernel         = %g\n', lenKernel);
  fprintf('tSec              = %g\n', tSec);
end

figure('WindowState','maximized')
for iSubplot = 1:2
  subplot(2, 1, iSubplot)
  switch(iSubplot)
    case 1
      semilogy(xAr, tSecAr, 'o-')
      title('Linear x, logarithmic y')
    case 2
      loglog(  xAr, tSecAr, 'o-')
      title('Log-log')
  end
  xlabel('x')
  ylabel('Time [s]')
  grid on
end
end



function tSec = test(dt, y1, tf, lenKernel, edgePolicy, hannWindowEnabled)
tTicToc = tic();

[~] = bicas.tf.time.apply_TF(dt, y1, tf, lenKernel, edgePolicy, hannWindowEnabled);

tSec = toc(tTicToc);
end
