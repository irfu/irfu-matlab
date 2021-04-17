function [out]=irf_convert_fac(inp,B0,r)
%IRF_CONVERT_FAC transforms to a field-aligned coordinate (FAC) system
%
%  out = irf_convert_fac(inp,B0,[r])
%
%  Transforms to a field-aligned coordinate (FAC) system defined as:
%  z = aligned with the background magnetic field
%  y = z cross product with the position vector r of the spacecraft 
%      (nominally eastward at the equator), if r not given used r=[1 0 0]
%  x = y x z
%  If inp is single column vector then assume that it is amplitude of vector 
%  along r direction. Output out has 2 columns corresponding to inp[perp, para] projections
%
%  rotMatrix = irf_convert_fac([],B0,[r])
%
%  Return only the rotation matrix to the FAC system
%
%
%  out = irf_convert_fac(inp,rotMatrix)
%
%  Perform the roration defined by rotMatrix
%
%
%  Input:
%    r   = position vector of spacecraft, columns (t x y z)
%    B0  = background magnetic field, columns (t x y z)
%    inp = vector that is to be transformed to FAC, columns (t x y z)
%			     or inp can be cell array of such vectors
%    out = output in the same form as inp
%
% Example:
%   Efac  = irf_convert_fac(Egse, Bgse, [1, 0, 0]);
%   Emfac = irf_convert_fac(Em, Bgse, Mdir);
%
% Note: all input parameters must be in the same coordinate system

% Begin temporary fix to convert TS format to older format (Must include spacecraft position)

isinpTS = isa(inp,'TSeries');
if isinpTS
  inptemp = inp;
  ttemp = inp.time.epochUnix;
  datatemp = double(inp.data);
  inp = [ttemp, double(datatemp)];
end
if isa(B0,'TSeries')
  ttemp = B0.time.epochUnix;
  datatemp = double(B0.data);
  B0 = [ttemp, datatemp];
end
% End of temporary fix

if nargin<3, r=[1 0 0];
end
if isa(r,'TSeries')
  ttemp = r.time.epochUnix;
  datatemp = double(r.data);
  r = [ttemp, datatemp];
end
Rpar = []; Rperpy = []; Rperpx = [];
if size(r, 2) == 3, r = [B0(1, 1), r];end

if isempty(inp) % return the transformation matrix instead
  r = irf_resamp(r,B0);
else
  if isstruct(B0) % tr matrix given as input
    rotMatrix = B0;
    Rperpx = zeros(length(rotMatrix.t),4); Rperpx(:,1) = rotMatrix.t;
    Rperpy = Rperpx; Rpar = Rperpx;
    Rperpx(:,2:4) = rotMatrix.rotMatrix(:,1,:);
    Rperpy(:,2:4) = rotMatrix.rotMatrix(:,2,:);
    Rpar(:,2:4) = rotMatrix.rotMatrix(:,3,:);
    if size(inp,1) ~= length(rotMatrix.t)
      Rperpx = irf_resamp(Rperpx,inp);
      Rperpy = irf_resamp(Rperpy,inp);
      Rpar = irf_resamp(Rpar,inp);
    end
  else
    if size(inp,1) ~= size(B0,1), B0 = irf_resamp(B0,inp); end
    if size(inp,1) ~= size(r,1),   r = irf_resamp(r,inp);  end
  end
end

if isempty(Rpar)
  % the direction of background magnetic field
  bn=irf_norm(B0);
  [Rpar,Rperpy,Rperpx]=define_ref_syst();
end

if iscell(inp)
  out=inp;
  for j = 1:numel(inp)
    if j>1 && (size(inp{j},1) ~= size(inp{1},1) ...
        || ~all(inp{j}(:,1) == inp{1}(:,1)) ) % the same time axis
      error('all inputs must have the same time axis')
    end
    out{j}=calculate_out(inp{j});
    if isinpTS
      inptemp.data =  out(:, 2:4);
      out = inptemp;
    end
  end
else
  if isempty(inp)  % return the transformation matrix instead
    out = struct('t',B0(:,1),'rotMatrix',zeros(length(B0(:,1)),3,3),...
      'b',B0(:,2:4),'r',r(:,2:4));
    out.rotMatrix(:,1,:) = Rperpx(:,2:4);
    out.rotMatrix(:,2,:) = Rperpy(:,2:4);
    out.rotMatrix(:,3,:) = Rpar(:,2:4);
  else
    out=calculate_out(inp);
    if isinpTS
      if size(out, 2) == 4
        out = irf.ts_vec_xyz(inptemp.time, out(:, 2: 4));
      elseif size(out, 2) == 3
        out = irf.ts_vec_xy(inptemp.time, out(:, 2: 3));
      end
      %inptemp.data =  out(:,[2:4]);
      %out = inptemp;
    end
  end
end

  function [Rpar,Rperpy,Rperpx]=define_ref_syst
    Rpar=bn;
    Rperpy=irf_norm(irf_cross(Rpar, r));
    Rperpx=irf_norm(irf_cross(Rperpy, B0));
  end
  function out=calculate_out(inp)
    %A_mfa=A;
    ndata=size(inp,1);
    ndim = size(inp, 2);
    if ndim == 4
      out=zeros(ndata,4);
      out(:,1)=inp(:,1);
      out(:,4)=  Rpar(:,2).*inp(:,2)+  Rpar(:,3).*inp(:,3)+  Rpar(:,4).*inp(:,4);
      out(:,2)=Rperpx(:,2).*inp(:,2)+Rperpx(:,3).*inp(:,3)+Rperpx(:,4).*inp(:,4);
      out(:,3)=Rperpy(:,2).*inp(:,2)+Rperpy(:,3).*inp(:,3)+Rperpy(:,4).*inp(:,4);
    elseif ndim ==2
      out=zeros(ndata,3);
      out(:, 1)=inp(:,1);
      out(:, 2)=inp(:, 2) .*(Rperpx(:,2) .* r(:, 2) + Rperpx(:,3) .* r(:, 3) + Rperpx(:,4) .* r(:, 4));
      out(:, 3)=inp(:, 2) .*(Rpar(:,2) .* r(:, 2) + Rpar(:,3) .* r(:, 3) + Rpar(:,4) .* r(:, 4));
    end
  end
end