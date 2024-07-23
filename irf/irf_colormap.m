function cmap1=irf_colormap(varargin)
% IRF_COLORMAP return colormap by name or apply and freeze the colormap
%
% CMAP = IRF_COLORMAP(colormap_name)
%  Colormap_names:
%       'standard'  - (default), same as 'space','cmap' (commonly used showing space data)
%       'poynting'  - white in center and blue/green for negative and red/black for positive values
%       'poynting_gray'  - gray in center and blue/green for negative and red/black for positive values
%       'solo'
%       'bluered'
%       'waterfall' - fancy-schmancy
%       'batlow' - scientific colour map, (v8.0.1, DOI: 10.5281/zenodo.8409685)
%
% IRF_COLORMAP(AX,colormap_name) - apply colormap to axis AX


[ax,args,nargs] = axescheck(varargin{:});

if nargs == 0 % show only help
  help irf_colormap;
  return
end

% check which axis to apply
if isempty(ax)
  axes(gca);
else
  axes(ax(1));
end

colormap_name=args{1};

load caa/cmap.mat % default map
if nargs > 0
  switch lower(colormap_name)
    case 'poynting'
      it=0:.02:1;it=it(:);
      cmap=[ [0*it flipud(it) it];[it it it*0+1];[it*0+1 flipud(it) flipud(it)]; [flipud(it) 0*it 0*it]]; clear it;
    case {'poynting_grey','poynting_gray'}
      it=0:.02:1;it=it(:);
      cmap=[ [0*it flipud(it) it];...
        [it*.8 it*.8 it*0+1];...
        [it*0+1 flipud(it*.8) flipud(it*.8)];...
        [flipud(it) 0*it 0*it]];
      clear it;
    case 'solo'
      it=0:.02:1;it=it(:);
      cmap=[ [it it it*0+1];[it*0+1 flipud(it) flipud(it)]; [flipud(it) 0*it 0*it]]; clear it;
    case {'parula','jet','hsv','hot','cool','spring','summer','autumn',...
        'winter','gray','bone','copper','pink','lines','colorcube','prism','flag','white'}
      cmap = colormap(colormap_name);
    case {'bluered'}
      rr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
      gg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
      bb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
      cmap = [rr' gg' bb'];
    case 'waterfall' % fancy-schmancy
      c = [55,137,187;...
        106,193,165;...
        172,220,166;...
        230,244,157;...
        255,254,194;...
        253,223,144;...
        251,173,104;...
        242,109,074;...
        211,064,082]/255;
      cmap = interp1(linspace(1,64,size(c,1)),c,1:64);
    case 'cubehelix' % stolen from matplotlib
      cmap = [0, 0, 0;
        1, 0, 1;
        3, 1, 3;
        4, 1, 4;
        6, 2, 6;
        8, 2, 8;
        9, 3, 9;
        10, 4, 11;
        12, 4, 13;
        13, 5, 15;
        14, 6, 17;
        15, 6, 19;
        17, 7, 21;
        18, 8, 23;
        19, 9, 25;
        20, 10, 27;
        20, 11, 29;
        21, 11, 31;
        22, 12, 33;
        23, 13, 35;
        23, 14, 37;
        24, 15, 39;
        24, 17, 41;
        25, 18, 43;
        25, 19, 45;
        25, 20, 47;
        26, 21, 48;
        26, 22, 50;
        26, 24, 52;
        26, 25, 54;
        26, 26, 56;
        26, 28, 57;
        26, 29, 59;
        26, 31, 60;
        26, 32, 62;
        26, 34, 63;
        26, 35, 65;
        25, 37, 66;
        25, 38, 67;
        25, 40, 69;
        25, 41, 70;
        24, 43, 71;
        24, 45, 72;
        24, 46, 73;
        23, 48, 74;
        23, 50, 74;
        23, 52, 75;
        23, 53, 76;
        22, 55, 76;
        22, 57, 77;
        22, 58, 77;
        21, 60, 77;
        21, 62, 78;
        21, 64, 78;
        21, 66, 78;
        21, 67, 78;
        21, 69, 78;
        20, 71, 78;
        20, 73, 78;
        20, 74, 77;
        21, 76, 77;
        21, 78, 77;
        21, 79, 76;
        21, 81, 76;
        21, 83, 75;
        22, 84, 75;
        22, 86, 74;
        22, 88, 73;
        23, 89, 73;
        23, 91, 72;
        24, 92, 71;
        25, 94, 70;
        26, 95, 69;
        27, 97, 68;
        27, 98, 67;
        28, 99, 66;
        30, 101, 66;
        31, 102, 65;
        32, 103, 64;
        33, 104, 63;
        35, 106, 61;
        36, 107, 60;
        38, 108, 59;
        39, 109, 58;
        41, 110, 58;
        43, 111, 57;
        45, 112, 56;
        47, 113, 55;
        49, 114, 54;
        51, 114, 53;
        53, 115, 52;
        55, 116, 51;
        57, 116, 51;
        60, 117, 50;
        62, 118, 49;
        65, 118, 49;
        67, 119, 48;
        70, 119, 48;
        72, 120, 47;
        75, 120, 47;
        78, 120, 47;
        81, 121, 46;
        83, 121, 46;
        86, 121, 46;
        89, 121, 46;
        92, 122, 46;
        95, 122, 47;
        98, 122, 47;
        101, 122, 47;
        104, 122, 48;
        107, 122, 48;
        110, 122, 49;
        113, 122, 50;
        116, 122, 50;
        120, 122, 51;
        123, 122, 52;
        126, 122, 53;
        129, 122, 55;
        132, 122, 56;
        135, 122, 57;
        138, 121, 59;
        141, 121, 60;
        144, 121, 62;
        147, 121, 64;
        150, 121, 65;
        153, 121, 67;
        155, 121, 69;
        158, 121, 71;
        161, 121, 74;
        164, 120, 76;
        166, 120, 78;
        169, 120, 81;
        171, 120, 83;
        174, 120, 86;
        176, 120, 88;
        178, 120, 91;
        181, 120, 94;
        183, 120, 96;
        185, 120, 99;
        187, 121, 102;
        189, 121, 105;
        191, 121, 108;
        193, 121, 111;
        194, 121, 114;
        196, 122, 117;
        198, 122, 120;
        199, 122, 124;
        201, 123, 127;
        202, 123, 130;
        203, 124, 133;
        204, 124, 136;
        205, 125, 140;
        206, 125, 143;
        207, 126, 146;
        208, 127, 149;
        209, 127, 153;
        209, 128, 156;
        210, 129, 159;
        211, 130, 162;
        211, 131, 165;
        211, 131, 169;
        212, 132, 172;
        212, 133, 175;
        212, 135, 178;
        212, 136, 181;
        212, 137, 184;
        212, 138, 186;
        212, 139, 189;
        212, 140, 192;
        211, 142, 195;
        211, 143, 197;
        211, 144, 200;
        210, 146, 203;
        210, 147, 205;
        210, 149, 207;
        209, 150, 210;
        208, 152, 212;
        208, 154, 214;
        207, 155, 216;
        207, 157, 218;
        206, 158, 220;
        205, 160, 222;
        205, 162, 224;
        204, 164, 226;
        203, 165, 227;
        203, 167, 229;
        202, 169, 230;
        201, 171, 231;
        201, 172, 233;
        200, 174, 234;
        199, 176, 235;
        199, 178, 236;
        198, 180, 237;
        197, 182, 238;
        197, 183, 239;
        196, 185, 239;
        196, 187, 240;
        195, 189, 241;
        195, 191, 241;
        194, 193, 242;
        194, 194, 242;
        194, 196, 242;
        193, 198, 243;
        193, 200, 243;
        193, 202, 243;
        193, 203, 243;
        193, 205, 243;
        193, 207, 243;
        193, 208, 243;
        193, 210, 243;
        193, 212, 243;
        193, 213, 243;
        194, 215, 242;
        194, 216, 242;
        195, 218, 242;
        195, 219, 242;
        196, 221, 241;
        196, 222, 241;
        197, 224, 241;
        198, 225, 241;
        199, 226, 240;
        200, 228, 240;
        200, 229, 240;
        202, 230, 239;
        203, 231, 239;
        204, 232, 239;
        205, 233, 239;
        206, 235, 239;
        208, 236, 238;
        209, 237, 238;
        210, 238, 238;
        212, 239, 238;
        213, 240, 238;
        215, 240, 238;
        217, 241, 238;
        218, 242, 238;
        220, 243, 239;
        222, 244, 239;
        223, 244, 239;
        225, 245, 240;
        227, 246, 240;
        229, 247, 240;
        231, 247, 241;
        232, 248, 242;
        234, 248, 242;
        236, 249, 243;
        238, 250, 244;
        240, 250, 245;
        242, 251, 246;
        244, 251, 247;
        245, 252, 248;
        247, 252, 249;
        249, 253, 250;
        251, 253, 252;
        253, 254, 253;
        255, 255, 255] / 255;
    case 'batlow'
      % scientific colour map from Fabio Crameri, v8.0.1, DOI: 10.5281/zenodo.8409685
      % SPDX-License-Identifier: MIT
      % LICENSE
      % Copyright (c) 2023, Fabio Crameri
      % Permission is hereby granted, free of charge, to any person ob-
      % taining a copy of this software and associated documentation
      % files (the ”Software”), to deal in the Software without restric-
      % tion, including without limitation the rights to use, copy, modify,
      % merge, publish, distribute, sublicense, and/or sell copies of the
      % Software, and to permit persons to whom the Software is fur-
      % nished to do so, subject to the following conditions:
      % The above copyright notice and this permission notice shall be
      % included in all copies or substantial portions of the Software.
      % THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF
      % ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
      % TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
      % PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
      % SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
      % ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN AC-
      % TION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
      % OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
      % OTHER DEALINGS IN THE SOFTWARE.
      % END OF LICENSE
      cmap = [0.005193, 0.098238, 0.349842;
        0.009065, 0.104487, 0.350933;
        0.012963, 0.110779, 0.351992;
        0.016530, 0.116913, 0.353070;
        0.019936, 0.122985, 0.354120;
        0.023189, 0.129035, 0.355182;
        0.026291, 0.135044, 0.356210;
        0.029245, 0.140964, 0.357239;
        0.032053, 0.146774, 0.358239;
        0.034853, 0.152558, 0.359233;
        0.037449, 0.158313, 0.360216;
        0.039845, 0.163978, 0.361187;
        0.042104, 0.169557, 0.362151;
        0.044069, 0.175053, 0.363084;
        0.045905, 0.180460, 0.364007;
        0.047665, 0.185844, 0.364915;
        0.049378, 0.191076, 0.365810;
        0.050795, 0.196274, 0.366684;
        0.052164, 0.201323, 0.367524;
        0.053471, 0.206357, 0.368370;
        0.054721, 0.211234, 0.369184;
        0.055928, 0.216046, 0.369974;
        0.057033, 0.220754, 0.370750;
        0.058032, 0.225340, 0.371509;
        0.059164, 0.229842, 0.372252;
        0.060167, 0.234299, 0.372978;
        0.061052, 0.238625, 0.373691;
        0.062060, 0.242888, 0.374386;
        0.063071, 0.247085, 0.375050;
        0.063982, 0.251213, 0.375709;
        0.064936, 0.255264, 0.376362;
        0.065903, 0.259257, 0.376987;
        0.066899, 0.263188, 0.377594;
        0.067921, 0.267056, 0.378191;
        0.069002, 0.270922, 0.378774;
        0.070001, 0.274713, 0.379342;
        0.071115, 0.278497, 0.379895;
        0.072192, 0.282249, 0.380434;
        0.073440, 0.285942, 0.380957;
        0.074595, 0.289653, 0.381452;
        0.075833, 0.293321, 0.381922;
        0.077136, 0.296996, 0.382376;
        0.078517, 0.300622, 0.382814;
        0.079984, 0.304252, 0.383224;
        0.081553, 0.307858, 0.383598;
        0.083082, 0.311461, 0.383936;
        0.084778, 0.315043, 0.384240;
        0.086503, 0.318615, 0.384506;
        0.088353, 0.322167, 0.384731;
        0.090281, 0.325685, 0.384910;
        0.092304, 0.329220, 0.385040;
        0.094462, 0.332712, 0.385116;
        0.096618, 0.336161, 0.385134;
        0.099015, 0.339621, 0.385090;
        0.101481, 0.343036, 0.384981;
        0.104078, 0.346410, 0.384801;
        0.106842, 0.349774, 0.384548;
        0.109695, 0.353098, 0.384217;
        0.112655, 0.356391, 0.383807;
        0.115748, 0.359638, 0.383310;
        0.118992, 0.362849, 0.382713;
        0.122320, 0.366030, 0.382026;
        0.125889, 0.369160, 0.381259;
        0.129519, 0.372238, 0.380378;
        0.133298, 0.375282, 0.379395;
        0.137212, 0.378282, 0.378315;
        0.141260, 0.381240, 0.377135;
        0.145432, 0.384130, 0.375840;
        0.149706, 0.386975, 0.374449;
        0.154073, 0.389777, 0.372934;
        0.158620, 0.392531, 0.371320;
        0.163246, 0.395237, 0.369609;
        0.167952, 0.397889, 0.367784;
        0.172788, 0.400496, 0.365867;
        0.177752, 0.403041, 0.363833;
        0.182732, 0.405551, 0.361714;
        0.187886, 0.408003, 0.359484;
        0.193050, 0.410427, 0.357177;
        0.198310, 0.412798, 0.354767;
        0.203676, 0.415116, 0.352253;
        0.209075, 0.417412, 0.349677;
        0.214555, 0.419661, 0.347019;
        0.220112, 0.421864, 0.344261;
        0.225707, 0.424049, 0.341459;
        0.231362, 0.426197, 0.338572;
        0.237075, 0.428325, 0.335634;
        0.242795, 0.430418, 0.332635;
        0.248617, 0.432493, 0.329571;
        0.254452, 0.434529, 0.326434;
        0.260320, 0.436556, 0.323285;
        0.266241, 0.438555, 0.320085;
        0.272168, 0.440541, 0.316831;
        0.278171, 0.442524, 0.313552;
        0.284175, 0.444484, 0.310243;
        0.290214, 0.446420, 0.306889;
        0.296294, 0.448357, 0.303509;
        0.302379, 0.450282, 0.300122;
        0.308517, 0.452205, 0.296721;
        0.314648, 0.454107, 0.293279;
        0.320834, 0.456006, 0.289841;
        0.327007, 0.457900, 0.286377;
        0.333235, 0.459794, 0.282937;
        0.339469, 0.461685, 0.279468;
        0.345703, 0.463563, 0.275998;
        0.351976, 0.465440, 0.272492;
        0.358277, 0.467331, 0.269037;
        0.364589, 0.469213, 0.265543;
        0.370922, 0.471085, 0.262064;
        0.377291, 0.472952, 0.258588;
        0.383675, 0.474842, 0.255131;
        0.390070, 0.476711, 0.251665;
        0.396505, 0.478587, 0.248212;
        0.402968, 0.480466, 0.244731;
        0.409455, 0.482351, 0.241314;
        0.415967, 0.484225, 0.237895;
        0.422507, 0.486113, 0.234493;
        0.429094, 0.488011, 0.231096;
        0.435714, 0.489890, 0.227728;
        0.442365, 0.491795, 0.224354;
        0.449052, 0.493684, 0.221074;
        0.455774, 0.495585, 0.217774;
        0.462539, 0.497497, 0.214518;
        0.469368, 0.499393, 0.211318;
        0.476221, 0.501314, 0.208148;
        0.483123, 0.503216, 0.205037;
        0.490081, 0.505137, 0.201976;
        0.497089, 0.507058, 0.198994;
        0.504153, 0.508984, 0.196118;
        0.511253, 0.510898, 0.193296;
        0.518425, 0.512822, 0.190566;
        0.525637, 0.514746, 0.187990;
        0.532907, 0.516662, 0.185497;
        0.540225, 0.518584, 0.183099;
        0.547599, 0.520486, 0.180884;
        0.555024, 0.522391, 0.178854;
        0.562506, 0.524293, 0.176964;
        0.570016, 0.526186, 0.175273;
        0.577582, 0.528058, 0.173775;
        0.585199, 0.529927, 0.172493;
        0.592846, 0.531777, 0.171449;
        0.600520, 0.533605, 0.170648;
        0.608240, 0.535423, 0.170104;
        0.615972, 0.537231, 0.169826;
        0.623739, 0.539002, 0.169814;
        0.631513, 0.540752, 0.170075;
        0.639301, 0.542484, 0.170622;
        0.647098, 0.544183, 0.171465;
        0.654889, 0.545863, 0.172603;
        0.662691, 0.547503, 0.174044;
        0.670477, 0.549127, 0.175747;
        0.678244, 0.550712, 0.177803;
        0.685995, 0.552274, 0.180056;
        0.693720, 0.553797, 0.182610;
        0.701421, 0.555294, 0.185478;
        0.709098, 0.556772, 0.188546;
        0.716731, 0.558205, 0.191851;
        0.724322, 0.559628, 0.195408;
        0.731878, 0.561011, 0.199174;
        0.739393, 0.562386, 0.203179;
        0.746850, 0.563725, 0.207375;
        0.754268, 0.565033, 0.211761;
        0.761629, 0.566344, 0.216322;
        0.768942, 0.567630, 0.221045;
        0.776208, 0.568899, 0.225930;
        0.783416, 0.570162, 0.230962;
        0.790568, 0.571421, 0.236160;
        0.797665, 0.572682, 0.241490;
        0.804709, 0.573928, 0.246955;
        0.811692, 0.575187, 0.252572;
        0.818610, 0.576462, 0.258303;
        0.825472, 0.577725, 0.264197;
        0.832272, 0.579026, 0.270211;
        0.838999, 0.580339, 0.276353;
        0.845657, 0.581672, 0.282631;
        0.852247, 0.583037, 0.289036;
        0.858747, 0.584440, 0.295572;
        0.865168, 0.585882, 0.302255;
        0.871505, 0.587352, 0.309112;
        0.877741, 0.588873, 0.316081;
        0.883878, 0.590450, 0.323195;
        0.889900, 0.592087, 0.330454;
        0.895809, 0.593765, 0.337865;
        0.901590, 0.595507, 0.345429;
        0.907242, 0.597319, 0.353142;
        0.912746, 0.599191, 0.360986;
        0.918103, 0.601126, 0.368999;
        0.923300, 0.603137, 0.377139;
        0.928323, 0.605212, 0.385404;
        0.933176, 0.607369, 0.393817;
        0.937850, 0.609582, 0.402345;
        0.942332, 0.611867, 0.411006;
        0.946612, 0.614218, 0.419767;
        0.950697, 0.616649, 0.428624;
        0.954574, 0.619137, 0.437582;
        0.958244, 0.621671, 0.446604;
        0.961696, 0.624282, 0.455702;
        0.964943, 0.626934, 0.464860;
        0.967983, 0.629639, 0.474057;
        0.970804, 0.632394, 0.483290;
        0.973424, 0.635183, 0.492547;
        0.975835, 0.638012, 0.501826;
        0.978052, 0.640868, 0.511090;
        0.980079, 0.643752, 0.520350;
        0.981918, 0.646664, 0.529602;
        0.983574, 0.649590, 0.538819;
        0.985066, 0.652522, 0.547998;
        0.986392, 0.655470, 0.557142;
        0.987567, 0.658422, 0.566226;
        0.988596, 0.661378, 0.575265;
        0.989496, 0.664329, 0.584246;
        0.990268, 0.667280, 0.593174;
        0.990926, 0.670230, 0.602031;
        0.991479, 0.673165, 0.610835;
        0.991935, 0.676091, 0.619575;
        0.992305, 0.679007, 0.628251;
        0.992595, 0.681914, 0.636869;
        0.992813, 0.684815, 0.645423;
        0.992967, 0.687705, 0.653934;
        0.993064, 0.690579, 0.662398;
        0.993111, 0.693451, 0.670810;
        0.993112, 0.696314, 0.679177;
        0.993074, 0.699161, 0.687519;
        0.993002, 0.702006, 0.695831;
        0.992900, 0.704852, 0.704114;
        0.992771, 0.707689, 0.712380;
        0.992619, 0.710530, 0.720639;
        0.992447, 0.713366, 0.728892;
        0.992258, 0.716210, 0.737146;
        0.992054, 0.719049, 0.745403;
        0.991837, 0.721893, 0.753673;
        0.991607, 0.724754, 0.761959;
        0.991367, 0.727614, 0.770270;
        0.991116, 0.730489, 0.778606;
        0.990855, 0.733373, 0.786976;
        0.990586, 0.736265, 0.795371;
        0.990307, 0.739184, 0.803810;
        0.990018, 0.742102, 0.812285;
        0.989720, 0.745039, 0.820804;
        0.989411, 0.747997, 0.829372;
        0.989089, 0.750968, 0.837979;
        0.988754, 0.753949, 0.846627;
        0.988406, 0.756949, 0.855332;
        0.988046, 0.759964, 0.864078;
        0.987672, 0.762996, 0.872864;
        0.987280, 0.766047, 0.881699;
        0.986868, 0.769105, 0.890573;
        0.986435, 0.772184, 0.899493;
        0.985980, 0.775272, 0.908448;
        0.985503, 0.778378, 0.917444;
        0.985002, 0.781495, 0.926468;
        0.984473, 0.784624, 0.935531;
        0.983913, 0.787757, 0.944626;
        0.983322, 0.790905, 0.953748;
        0.982703, 0.794068, 0.962895;
        0.982048, 0.797228, 0.972070;
        0.981354, 0.800406, 0.981267];
  end
end

if nargout == 0 % apply the colormap and freeze
  colormap(cmap);
  freezeColors;
  hcb = cbhandle;
  if hcb % workaround cbfreeze bug that cbfreeze removes cblabel
    hy=get(hcb,'ylabel');
    ylabel_string=get(hy,'string');
    ylabel_fontsize=get(hy,'fontsize');
    new_hcb = cbfreeze(hcb);
    new_hy=get(new_hcb,'ylabel');
    set(new_hy,'string',ylabel_string,'fontsize',ylabel_fontsize);
  end
  %    cbfreeze;
elseif nargout == 1 % only return colormap
  cmap1=cmap;
end

