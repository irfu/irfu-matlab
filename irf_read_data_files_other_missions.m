function [res]=irf_read_data_files_other_missions(file_name,data_source)
% [res]=irf_read_data_other_missions(file_name,data_source);
%
% data_source:
%   DMSP_SSIES - data from http://cindispace.utdallas.edu/DMSP/dmsp_data_at_utdallas.html
%   IAGA2002 -  e.g. AE index from WDC http://wdc.kugi.kyoto-u.ac.jp/aeasy/index.html

switch lower(data_source)
    case 'dmsp_ssies'
        fid = fopen(file_name);
        C = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',3);
        fclose(fid);
        year=1900+floor(C{1}/1000); %
        zero=zeros(size(year));
        one=zero+1;
        res.t=toepoch([year one one zero zero zero])+mod(C{1},1000)*24*3600+C{2}; % time in isdat epoch
        res.Vx=C{10};res.Vx(res.Vx==-9999.0)=NaN;
        res.Vy=C{11};res.Vx(res.Vy==-9999.0)=NaN;
        res.Vz=C{12};res.Vx(res.Vz==-9999.0)=NaN;
        res.n=C{16};
    case 'iaga2002'
        fid = fopen(file_name);
        for j=1:3, fgetl(fid);end
        iagacode=textscan(fgetl(fid),'%s');
        for j=5:14, fgetl(fid);end
        header=textscan(fgetl(fid),'%s');
        %        C = textscan(fid,'%f-%f-%f %f:%f:%f %*f %f %f %f %f','headerlines',15);
        switch lower(iagacode{1}{3})
            case 'ae'
                C = textscan(fid,'%f-%f-%f %f:%f:%f %*f %f %f %f %f');
                fclose(fid);
                res.t=toepoch([C{1} C{2} C{3} C{4} C{5} C{6}]); % time in isdat epoch
                res.AE=C{7};
                res.AU=C{8};
                res.AL=C{9};
                res.AO=C{10};
            otherwise % assume magnetometer data
                res.station=iagacode{1}{3}; % magnetometer station name
                C = textscan(fid,'%f-%f-%f %f:%f:%f %*f %f %f %f %f');
                fclose(fid);                
                res.t=toepoch([C{1} C{2} C{3} C{4} C{5} C{6}]); % time in isdat epoch
                res.D=C{7};
                res.H=C{8};
                res.Z=C{9};
                res.F=C{10};
        end
        
    otherwise
        disp('Data source not recognized');
        disp('Reading assuming first row is variable, second units and then comes data');
end

