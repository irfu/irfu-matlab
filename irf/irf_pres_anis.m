function res = irf_pres_anis(Pi, B)
%IRF_PRES_ANIS  Compute ion pressure anisotropy factor: (Pi_para - Pi_perp) * mu0 / B^2
%
%   alpha = irf_pres_anis(Pi, B)
%       Units: Pi [nPa], B [nT]
%       Pi = [T, Ppar P12 P13; P12 Pperp1 P23; P13 P23 Pperp2]
%       B = [T, B1, B2, B3]
%       resampling B to Pi

%%  1. default flags false
    inputTSeries   = false;
    inputNumeric   = false;

%%  2. check input type
    if isa(Pi,'TSeries') || isa(B,'TSeries')
        inputTSeries = true;
    elseif isnumeric(Pi) && isnumeric(B)
        inputNumeric = true;
    else
        errStr = 'irf_pres_anis: input neither TSeries or numeric.';
        irf.log('critical',errStr);error(errStr);
    end

%%  3. compute
    if inputTSeries
        Bres = B.resample(Pi);
        Pi_para = Pi.xx;
        Pi_perp = (Pi.yy + Pi.zz)/2;
        res = (Pi_para - Pi_perp) / Bres.abs2 * 4 * 3.14159e2;
        res = irf.ts_scalar(res.time, res.data);        % tensorOrder from '2' to '0'        
        if strcmp(Pi.units, 'nPa') && strcmp (B.units,'nT')
			res.units = '1';
        else
			irf.log('warning','irf_pres_anis: units not correct!'); % TODO: implement units
        end
		res.name = 'ion pressure anisotropy';
		res.userData.LABLAXIS = 'alpha';
    end
    
end
