function Output = run(PlasmaModel,InputParameters)
% WHAMP.RUN run WHAMP code 
%
% PlasmaModel is a structure:
%    PlasmaModel.B        - B field in nT
%    PlasmaModel.Species  - cell array of Species
%
% Species is a structure:
%    Species.m      - mass/mass_proton (0 means electron)
%    Species.n      - density in cm^-3
%    Species.t      - temperature in eV
%    Species.a      - anisotropy, Tperp/Tpar
%    Species.d      - loss cone parameter 1
%    Species.b      - loss cone parameter 2
%    Species.vd     - Vpar_drift/Vthermal 
%
% InputParameters is a structure: 
%     InputParameters.fstart - start frequency
%     InputParameters.kp     - kp (perpendicular wave vector) as scalar or
%                              [kp_start kp_step kp_end]
%     InputParameters.kp     - kz (parallel wave vector) or 
%                              [kz_start kz_step kz_end]
%     
%     InputParameters.useLog - 1 input log10(p) and log10(z),(default)
%                              0 input given as p and z
%     InputParameters.fstart - start frequency
%
% Output.InputParameters
% Output.PlasmaModel
% Output.p
% Output.z
% Output.f
% Output.E
% Output.Ex
% Output.Ey
% Output.Ez
% Output.B
% Output.Bx
% Output.By
% Output.Bz
% Output.EB
% Output.VGP
% Output.VGZ
% Output.S
% Output.Sx
% Output.Sy
% Output.Sz
% Output.u
% Output.SGP
% Output.SGZ
%
% Examples:
%    Inp=struct('fstart',0.4,'kp',ad);Out=whamp.run(Plasma,Inp);











