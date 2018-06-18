function QR = thor_QR(origR,Rout)
% Quality factor for bow shock crossing.
%   R<Rout: QR = 1
%   R>Rout: QR = (Rout/R)^3
%
%   QR = THOR_QR(R,Rout);

if isa(origR,'TSeries') 
  Time = origR.time;
  R = origR.data;   
else 
  R = origR;
end

QR = (Rout./R).^(3);
QR(R<Rout) = 1;

if isa(origR,'TSeries')   
  QR = irf.ts_scalar(origR.time,QR);
  QR.name = 'Quality factor QR';
end