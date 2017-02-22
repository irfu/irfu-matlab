function res = ts2mat(Ts)
%IRF.TS2MAT  Convert TSeries to MAT
%
%   res = IRF.TS2MAT(Ts)

res = [];
if isempty(Ts), return, end

if ~isa(Ts,'TSeries') 
  error('TS must be TSeries object')
elseif Ts.tensorOrder>1
  error('TS must be scalar or vector')
end

res = [Ts.time.epochUnix double(Ts.data)];