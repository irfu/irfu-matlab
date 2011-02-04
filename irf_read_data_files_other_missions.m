function [res]=irf_read_data_files_other_missions(file_name,data_source)
% [res]=irf_read_data_other_missions(file_name,data_source);
%
% data_source:
%   DMSP_SSIES - data from http://cindispace.utdallas.edu/DMSP/dmsp_data_at_utdallas.html
%

switch data_source
    case 'DMSP_SSIES'
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
    otherwise
        disp('Data source not recognized');
        disp('Reading assuming first row is variable, second units and then comes data');
end

