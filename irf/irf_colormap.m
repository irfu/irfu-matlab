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
%       'batlow' - scientific colour map
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
    case 'batlow' % scientific colour map from Fabio Crameri
      cmap = [0.0052    0.0982    0.3498;
      0.0091    0.1045    0.3509;
      0.0130    0.1108    0.3520;
      0.0165    0.1169    0.3531;
      0.0199    0.1230    0.3541;
      0.0232    0.1290    0.3552;
      0.0263    0.1350    0.3562;
      0.0292    0.1410    0.3572;
      0.0321    0.1468    0.3582;
      0.0349    0.1526    0.3592;
      0.0374    0.1583    0.3602;
      0.0398    0.1640    0.3612;
      0.0421    0.1696    0.3622;
      0.0441    0.1751    0.3631;
      0.0459    0.1805    0.3640;
      0.0477    0.1858    0.3649;
      0.0494    0.1911    0.3658;
      0.0508    0.1963    0.3667;
      0.0522    0.2013    0.3675;
      0.0535    0.2064    0.3684;
      0.0547    0.2112    0.3692;
      0.0559    0.2160    0.3700;
      0.0570    0.2208    0.3707;
      0.0580    0.2253    0.3715;
      0.0592    0.2298    0.3723;
      0.0602    0.2343    0.3730;
      0.0611    0.2386    0.3737;
      0.0621    0.2429    0.3744;
      0.0631    0.2471    0.3750;
      0.0640    0.2512    0.3757;
      0.0649    0.2553    0.3764;
      0.0659    0.2593    0.3770;
      0.0669    0.2632    0.3776;
      0.0679    0.2671    0.3782;
      0.0690    0.2709    0.3788;
      0.0700    0.2747    0.3793;
      0.0711    0.2785    0.3799;
      0.0722    0.2822    0.3804;
      0.0734    0.2859    0.3810;
      0.0746    0.2897    0.3815;
      0.0758    0.2933    0.3819;
      0.0771    0.2970    0.3824;
      0.0785    0.3006    0.3828;
      0.0800    0.3043    0.3832;
      0.0816    0.3079    0.3836;
      0.0831    0.3115    0.3839;
      0.0848    0.3150    0.3842;
      0.0865    0.3186    0.3845;
      0.0884    0.3222    0.3847;
      0.0903    0.3257    0.3849;
      0.0923    0.3292    0.3850;
      0.0945    0.3327    0.3851;
      0.0966    0.3362    0.3851;
      0.0990    0.3396    0.3851;
      0.1015    0.3430    0.3850;
      0.1041    0.3464    0.3848;
      0.1068    0.3498    0.3845;
      0.1097    0.3531    0.3842;
      0.1127    0.3564    0.3838;
      0.1157    0.3596    0.3833;
      0.1190    0.3628    0.3827;
      0.1223    0.3660    0.3820;
      0.1259    0.3692    0.3813;
      0.1295    0.3722    0.3804;
      0.1333    0.3753    0.3794;
      0.1372    0.3783    0.3783;
      0.1413    0.3812    0.3771;
      0.1454    0.3841    0.3758;
      0.1497    0.3870    0.3744;
      0.1541    0.3898    0.3729;
      0.1586    0.3925    0.3713;
      0.1632    0.3952    0.3696;
      0.1680    0.3979    0.3678;
      0.1728    0.4005    0.3659;
      0.1778    0.4030    0.3638;
      0.1827    0.4056    0.3617;
      0.1879    0.4080    0.3595;
      0.1931    0.4104    0.3572;
      0.1983    0.4128    0.3548;
      0.2037    0.4151    0.3523;
      0.2091    0.4174    0.3497;
      0.2146    0.4197    0.3470;
      0.2201    0.4219    0.3443;
      0.2257    0.4240    0.3415;
      0.2314    0.4262    0.3386;
      0.2371    0.4283    0.3356;
      0.2428    0.4304    0.3326;
      0.2486    0.4325    0.3296;
      0.2545    0.4345    0.3264;
      0.2603    0.4366    0.3233;
      0.2662    0.4386    0.3201;
      0.2722    0.4405    0.3168;
      0.2782    0.4425    0.3136;
      0.2842    0.4445    0.3102;
      0.2902    0.4464    0.3069;
      0.2963    0.4484    0.3035;
      0.3024    0.4503    0.3001;
      0.3085    0.4522    0.2967;
      0.3146    0.4541    0.2933;
      0.3208    0.4560    0.2898;
      0.3270    0.4579    0.2864;
      0.3332    0.4598    0.2829;
      0.3395    0.4617    0.2795;
      0.3457    0.4636    0.2760;
      0.3520    0.4654    0.2725;
      0.3583    0.4673    0.2690;
      0.3646    0.4692    0.2655;
      0.3709    0.4711    0.2621;
      0.3773    0.4730    0.2586;
      0.3837    0.4748    0.2551;
      0.3901    0.4767    0.2517;
      0.3965    0.4786    0.2482;
      0.4030    0.4805    0.2447;
      0.4095    0.4824    0.2413;
      0.4160    0.4842    0.2379;
      0.4225    0.4861    0.2345;
      0.4291    0.4880    0.2311;
      0.4357    0.4899    0.2277;
      0.4424    0.4918    0.2244;
      0.4491    0.4937    0.2211;
      0.4558    0.4956    0.2178;
      0.4625    0.4975    0.2145;
      0.4694    0.4994    0.2113;
      0.4762    0.5013    0.2081;
      0.4831    0.5032    0.2050;
      0.4901    0.5051    0.2020;
      0.4971    0.5071    0.1990;
      0.5042    0.5090    0.1961;
      0.5113    0.5109    0.1933;
      0.5184    0.5128    0.1906;
      0.5256    0.5147    0.1880;
      0.5329    0.5167    0.1855;
      0.5402    0.5186    0.1831;
      0.5476    0.5205    0.1809;
      0.5550    0.5224    0.1789;
      0.5625    0.5243    0.1770;
      0.5700    0.5262    0.1753;
      0.5776    0.5281    0.1738;
      0.5852    0.5299    0.1725;
      0.5928    0.5318    0.1714;
      0.6005    0.5336    0.1706;
      0.6082    0.5354    0.1701;
      0.6160    0.5372    0.1698;
      0.6237    0.5390    0.1698;
      0.6315    0.5408    0.1701;
      0.6393    0.5425    0.1706;
      0.6471    0.5442    0.1715;
      0.6549    0.5459    0.1726;
      0.6627    0.5475    0.1740;
      0.6705    0.5491    0.1757;
      0.6782    0.5507    0.1778;
      0.6860    0.5523    0.1801;
      0.6937    0.5538    0.1826;
      0.7014    0.5553    0.1855;
      0.7091    0.5568    0.1885;
      0.7167    0.5582    0.1919;
      0.7243    0.5596    0.1954;
      0.7319    0.5610    0.1992;
      0.7394    0.5624    0.2032;
      0.7469    0.5637    0.2074;
      0.7543    0.5650    0.2118;
      0.7616    0.5663    0.2163;
      0.7689    0.5676    0.2210;
      0.7762    0.5689    0.2259;
      0.7834    0.5702    0.2310;
      0.7906    0.5714    0.2362;
      0.7977    0.5727    0.2415;
      0.8047    0.5739    0.2470;
      0.8117    0.5752    0.2526;
      0.8186    0.5765    0.2583;
      0.8255    0.5777    0.2642;
      0.8323    0.5790    0.2702;
      0.8390    0.5803    0.2764;
      0.8457    0.5817    0.2826;
      0.8522    0.5830    0.2890;
      0.8587    0.5844    0.2956;
      0.8652    0.5859    0.3023;
      0.8715    0.5874    0.3091;
      0.8777    0.5889    0.3161;
      0.8839    0.5904    0.3232;
      0.8899    0.5921    0.3305;
      0.8958    0.5938    0.3379;
      0.9016    0.5955    0.3454;
      0.9072    0.5973    0.3531;
      0.9127    0.5992    0.3610;
      0.9181    0.6011    0.3690;
      0.9233    0.6031    0.3771;
      0.9283    0.6052    0.3854;
      0.9332    0.6074    0.3938;
      0.9378    0.6096    0.4023;
      0.9423    0.6119    0.4110;
      0.9466    0.6142    0.4198;
      0.9507    0.6166    0.4286;
      0.9546    0.6191    0.4376;
      0.9582    0.6217    0.4466;
      0.9617    0.6243    0.4557;
      0.9649    0.6269    0.4649;
      0.9680    0.6296    0.4741;
      0.9708    0.6324    0.4833;
      0.9734    0.6352    0.4925;
      0.9758    0.6380    0.5018;
      0.9781    0.6409    0.5111;
      0.9801    0.6438    0.5203;
      0.9819    0.6467    0.5296;
      0.9836    0.6496    0.5388;
      0.9851    0.6525    0.5480;
      0.9864    0.6555    0.5571;
      0.9876    0.6584    0.5662;
      0.9886    0.6614    0.5753;
      0.9895    0.6643    0.5842;
      0.9903    0.6673    0.5932;
      0.9909    0.6702    0.6020;
      0.9915    0.6732    0.6108;
      0.9919    0.6761    0.6196;
      0.9923    0.6790    0.6283;
      0.9926    0.6819    0.6369;
      0.9928    0.6848    0.6454;
      0.9930    0.6877    0.6539;
      0.9931    0.6906    0.6624;
      0.9931    0.6935    0.6708;
      0.9931    0.6963    0.6792;
      0.9931    0.6992    0.6875;
      0.9930    0.7020    0.6958;
      0.9929    0.7049    0.7041;
      0.9928    0.7077    0.7124;
      0.9926    0.7105    0.7206;
      0.9924    0.7134    0.7289;
      0.9923    0.7162    0.7371;
      0.9921    0.7190    0.7454;
      0.9918    0.7219    0.7537;
      0.9916    0.7248    0.7620;
      0.9914    0.7276    0.7703;
      0.9911    0.7305    0.7786;
      0.9909    0.7334    0.7870;
      0.9906    0.7363    0.7954;
      0.9903    0.7392    0.8038;
      0.9900    0.7421    0.8123;
      0.9897    0.7450    0.8208;
      0.9894    0.7480    0.8294;
      0.9891    0.7510    0.8380;
      0.9888    0.7539    0.8466;
      0.9884    0.7569    0.8553;
      0.9880    0.7600    0.8641;
      0.9877    0.7630    0.8729;
      0.9873    0.7660    0.8817;
      0.9869    0.7691    0.8906;
      0.9864    0.7722    0.8995;
      0.9860    0.7753    0.9084;
      0.9855    0.7784    0.9174;
      0.9850    0.7815    0.9265;
      0.9845    0.7846    0.9355;
      0.9839    0.7878    0.9446;
      0.9833    0.7909    0.9537;
      0.9827    0.7941    0.9629;
      0.9820    0.7972    0.9721;
      0.9814    0.8004    0.9813];
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

