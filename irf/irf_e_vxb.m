function res = irf_e_vxb(v,b,flag)
%IRF_E_VXB   Compute VxB and ExB/B^2
%
% [E] = IRF_E_VXB(V,B)
%        Calculate electric field E = -VxB
%
% [VExB] = IRF_E_VXB(E,B,-1)
%     Calculate convection velocity VExB = ExB/B^2
%
%     Units: v[km/s], B[nT], E[mV/m]
%     V,B can be TSeries or arrays. If arrays, then:
%     V = [T VX VY VZ] (T..VZ column vectors)
%     B = [T BX BY BZ]
%     V and B can be at different sampling
%     Resulting E is at B sampling, VExB is at E sampling

%% Default flags false
estimateExB    = false;
estimateVxB    = false;
inputTSeries   = false;
inputNumeric   = false;
inputVConstant = false;

%% Check input type

if (nargin ==3) && (flag == -1)
  estimateExB = true;
else
  estimateVxB = true;
end
if isa(v,'TSeries') || isa(b,'TSeries')
  inputTSeries = true;
elseif isnumeric(v) && isnumeric(b)
  inputNumeric = true;
else
  errStr = 'irf_e_vxb: input neither TSeries or numeric.';
  irf.log('critical',errStr);error(errStr);
end
if isnumeric(v) && all( size(v) == [1 3])
  inputVConstant = true;
end

%% do calculation
if inputNumeric
  if estimateExB
    e = v;
    if size(b,1) ~= size(e,1)
      b = irf_resamp(b,e);
    end
    res = irf_vec_x_scal( irf_cross(e,b), [b(:,1) irf_abs(b,1)], -2 );
    res = irf_tappl(res,'*1e3');

  elseif estimateVxB
    if inputVConstant
      v = [b(:,1) repmat(v,[size(b,1) 1])];
    end
    if size(v,2) <= 3
      error('irf_e_vxb: v has too few components');
    end
    if size(v,1) == 1
      v = irf_resamp(v,b);
    else
      b = irf_resamp(b,v);
    end
    res = irf_tappl( irf_cross(v,b), '*1e3*1e-9*1e3*(-1)' );
  end
end

if inputTSeries
  if estimateExB
    e = v;
    if e.length ~= b.length
      b = b.resample(e.time);
    end
    res = cross(e,b);
    res.data = res.data./(b.abs.data.^2*[1 1 1])*1e3;
    if strcmp(e.units,'mV/m') && strcmp (b.units,'nT')
      res.units = 'km/s';
    else
      irf.log('warning','irf_e_vxb: units not correct!'); % TODO: implement units
    end
    res.name = 'Velocity';
    res.userData.LABLAXIS = 'V';
  elseif estimateVxB
    if inputVConstant
      res = b;
      res.data = cross(repmat(v,res.length, 1),b.data)*(-1)*1e-3;
      v = [];
      v.units = 'km/s';
    else
      res = cross(v,b)*(-1)*1e-3;
    end
    if (strcmp(v.units,'km/s') || strcmp(v.units,'km s^-1')) ...
        && strcmp (b.units,'nT')
      res.units = 'mV/m';
    else
      irf.log('warning','irf_e_vxb: units not correct!'); % TODO: implement units
    end
    res.units = 'mV/m';
    res.name = 'Electric field';
    res.userData.LABLAXIS = 'E';

  end
end

