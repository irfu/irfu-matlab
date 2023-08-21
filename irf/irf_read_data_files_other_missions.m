function [res]=irf_read_data_files_other_missions(file_name,data_source)
% [res]=irf_read_data_other_missions(file_name,data_source);
%
% data_source:
%   DMSP_SSIES - data from http://cindispace.utdallas.edu/DMSP/dmsp_data_at_utdallas.html
%   DMSP_SSJ - data from SSJ exported using DMSP Tool http://sd-www.jhuapl.edu/Aurora/spectrogram/index.html
%   IAGA2002 -  e.g. AE index from WDC http://wdc.kugi.kyoto-u.ac.jp/aeasy/index.html

switch lower(data_source)
  case 'dmsp_ssies'
    fid = fopen(file_name);
    C = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',3);
    fclose(fid);
    year=1900+floor(C{1}/1000); %
    zero=zeros(size(year));
    one=zero+1;
    res.t=irf_time([year one one zero zero zero])+(mod(C{1},1000)-1)*24*3600+C{2}; % time in isdat epoch
    res.Vx=C{10};res.Vx(res.Vx==-9999.0)=NaN;
    res.Vy=C{11};res.Vx(res.Vy==-9999.0)=NaN;
    res.Vz=C{12};res.Vx(res.Vz==-9999.0)=NaN;
    res.n=C{16};
  case 'dmsp_ssj'
    C = txt2mat(file_name);
    res.t=irf_time([C(:,1) -C(:,2) -C(:,3) C(:,4) C(:,5) C(:,6)]); % time in isdat epoch
    res.GLAT=C(:,7);
    res.GLON=C(:,8);
    res.MLAT=C(:,9);
    res.MLT=C(:,10);
    res.JEe=C(:,11);
    res.AvgEe=C(:,12);
    res.JEi=C(:,13);
    res.AvgEi=C(:,14);
    ii=15;
    res.eCts=C(:,ii:ii+18);ii=ii+19;
    res.eDef=C(:,ii:ii+18);ii=ii+19;
    res.eErr=C(:,ii:ii+18);ii=ii+19;
    res.iCts=C(:,ii:ii+18);ii=ii+19;
    res.iDef=C(:,ii:ii+18);ii=ii+19;
    res.iErr=C(:,ii:ii+18);ii=ii+19; %#ok<NASGU>
    res.eEn=[3.100000E+01 4.500000E+01 6.600000E+01 9.700000E+01 1.420000E+02 2.080000E+02 3.040000E+02 4.450000E+02 6.520000E+02 9.540000E+02 1.396000E+03 2.057000E+03 3.030000E+03 4.465000E+03 6.578000E+03 9.693000E+03 1.428100E+04 2.104100E+04 3.100100E+04]';  % F14
    res.iEn=[ 3.200000E+01 4.600000E+01 6.800000E+01 9.900000E+01 1.450000E+02 2.130000E+02 3.120000E+02 4.570000E+02 6.690000E+02 9.800000E+02 1.398000E+03 2.056000E+03 3.023000E+03 4.445000E+03 6.537000E+03 9.612000E+03 1.413400E+04 2.078400E+04 3.056200E+04]'; % F14
    res.Def_unit='eV/cm^2 s sr eV';
    res.JE_unit='eV/cm^2 s sr';
    res.En_unit='eV';
  case 'dmsp_ssm'
    fid = fopen(file_name);
    C = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',28);
    fclose(fid);
    res.t=irf_time([C{1} C{2} C{3} C{4} C{5} C{6}]); % time in isdat epoch
    res.lat=C{7};
    res.long=C{8};
    res.alt=C{9};
    res.mdiplat=C{10};
    res.mlat=C{11};
    res.mlt=C{12};
    res.Bx=C{13};res.By=C{14};res.Bz=C{15};res.B=C{16};
    res.dBx=C{17};res.dBy=C{18};res.dBz=C{19};res.dB=C{20};
  case 'iaga2002'
    fid = fopen(file_name);
    for j=1:3, fgetl(fid);end
    iagacode=textscan(fgetl(fid),'%s');
    for j=5:14, fgetl(fid);end
    textscan(fgetl(fid),'%s'); % header
    %        C = textscan(fid,'%f-%f-%f %f:%f:%f %*f %f %f %f %f','headerlines',15);
    switch lower(iagacode{1}{3})
      case 'ae'
        C = textscan(fid,'%f-%f-%f %f:%f:%f %*f %f %f %f %f');
        fclose(fid);
        res.t=irf_time([C{1} C{2} C{3} C{4} C{5} C{6}]); % time in isdat epoch
        res.AE=C{7};
        res.AU=C{8};
        res.AL=C{9};
        res.AO=C{10};
      otherwise % assume magnetometer data
        res.station=iagacode{1}{3}; % magnetometer station name
        C = textscan(fid,'%f-%f-%f %f:%f:%f %*f %f %f %f %f');
        fclose(fid);
        res.t=irf_time([C{1} C{2} C{3} C{4} C{5} C{6}]); % time in isdat epoch
        res.D=C{7};
        res.H=C{8};
        res.Z=C{9};
        res.F=C{10};
    end

  otherwise
    disp('Data source not recognized');
    disp('Reading assuming first row is variable, second units and then comes data');
end

