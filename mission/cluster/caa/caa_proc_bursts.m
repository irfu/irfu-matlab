function caa_proc_bursts(burstfile, plotFlag)
% Run caa_get_bursts on a list file with internal burst file names from
% /data/cluster/burst.
%
% burstfile : burst list file name
% ibplot    : plot flag 0=no 1=yes
%
% Failed files are logged in the bflog file
%
% Create an iburst list file like this in a shell:
% cd /data/cluster/burst
% ls -1 09*we.0? >~/ib2009.list
% The list file is in the home directory ~/.
%
% Ex. data file:
% 020930090251we.04
% 020930103730we.02
% :
% :
%
    narginchk(1,2);
    if nargin < 2
        plotFlag = 0;
    end

    cd([getenv('HOME') '/matlab']);
    if ismac
        bflog=['/Volumes/caa/log/' burstfile '.ibfail.log'];
    else
       bflog=['/data/caa/log/' burstfile '.ibfail.log'];
%       bflog=[getenv('HOME') '/iburstfail.log'];
    end
    fid = fopen(burstfile,'r'); %Open the text file that contain the information about the burst.
    if fid==-1
        error(['Can not find burst list file ' burstfile]);
    end 
    fout = fopen(bflog,'a');
    if fout==-1
        fclose(fid);
        error(['Can not open/write file ' bflog]);
    end 
    while ~feof(fid)
        tline = fgetl(fid);
        if length(tline)>=17
            fn = tline(1:17)
            ret = caa_get_bursts(fn,plotFlag);
            if ret
                fprintf(fout,'%s ret:%d\n',fn,ret);
                irf_log('proc',['Burst file ' fn ' failed!']);
            end
        end
    end
    fclose(fid);
    fclose(fout);
end
