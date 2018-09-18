function eDist_new = remove_edist_background(varargin)
% MMS.REMOVE_EDIST_BACKGROUND: remove secondary photoelectrons from electron distribution function
%
% [eDist_new] = MMS..REMOVE_EDIST_BACKGROUND(eDist, Tint)
%
% Input: 
%       eDist - Electron distribution as a TSeries
%
% Options:
%   'tint' - time interval
% 
% References:
%   1. MMS FPI Data Users Guide: 
%       https://lasp.colorado.edu/galaxy/display/MFDPG/DES+Photoelectrons+-+further+details
% 
% Next step: move ion distribution background noise
% 
%   Example:     
%     eDist1_bgremoved = remove_edist_background(eDist1, 'tint', Tint);
%     eDist1_bgremoved = remove_edist_background(eDist1);
%
% 2018-09-18, wyli   

    %%  1. basic
    [~, args, nargs] = axescheck(varargin{:});
    if nargs == 1 
        eDist = args{1};
        Tint = irf.tint(eDist.time(1), eDist.time(end));
    elseif nargs ==3
        if strcmp(args{2}, 'tint')  %
        eDist = args{1};
        Tint = args{3};
        end
    end
    ic = eDist.name(4);
    eDist_tmp = eDist.tlim(Tint);
    eDist_new = eDist_tmp;
    
    %%  2. load data 
    time_mid = EpochUnix((Tint.start.epochUnix + Tint.stop.epochUnix)/2);
    % find the photoelectron model & scale (photo electron density)
    c_eval('emoms_fn = mms.get_filepath([''mms?_fpi_brst_l2_des-moms''], time_mid);', ic);
    c_eval('edist_fn = mms.get_filepath([''mms?_fpi_brst_l2_des-dist''], time_mid);', ic);
    emoms_obj = dataobj(emoms_fn);
    edist_obj = dataobj(edist_fn);
    estartdelphi_count = get_ts(edist_obj, 'mms1_des_startdelphi_count_brst');
    estartdelphi_angle = get_ts(edist_obj, 'mms1_des_startdelphi_angle_brst');
    ebgdist_fn = emoms_obj.GlobalAttributes.Photoelectron_model_filenames{1};
    ephoto_scale = str2double(emoms_obj.GlobalAttributes.Photoelectron_model_scaling_factor{1});
    estartdelphi_count_Tint = estartdelphi_count.tlim(Tint);
    estartdelphi_angle_Tint = estartdelphi_angle.tlim(Tint);
    
    %%  3. load ebgdist data 
    localFileDbRoot = datastore('mms_db', 'local_file_db_root');
    ebgdist_obj = dataobj([localFileDbRoot '/models/fpi/' ebgdist_fn]);
    ebgdist0 = get_variable(ebgdist_obj, 'mms_des_bgdist_p0_brst');
    ebgdist1 = get_variable(ebgdist_obj, 'mms_des_bgdist_p1_brst');
    energy0 = get_variable(ebgdist_obj, 'mms_des_energy0_brst');
    energy1 = get_variable(ebgdist_obj, 'mms_des_energy1_brst');
    phi = get_variable(ebgdist_obj, 'mms_des_phi_brst');
    theta = get_variable(ebgdist_obj, 'mms_des_phi_brst');
    
    %%  4. remove ebgdist
    ntime = length(eDist_tmp.time);
    edist_size = size(eDist_new.data);
    eDist_new_tmp = zeros(1, edist_size(2), edist_size(3), edist_size(4));      % preallocating for speed
    tic;
    for it = 1: ntime
        iestartdelphi_count = double(estartdelphi_count_Tint.data(it));
        iestartdelphi_angle = double(estartdelphi_angle_Tint.data(it));
        iebgdist = fix(iestartdelphi_count/16);
        if iebgdist == 0;       iebgdist = 360;     end
%         if eDist_new.ancillary.esteptable(it) == 0
%             eDist_new_tmp = eDist_new.data(it, :, :, :) - ephoto_scale * ebgdist0.data(iebgdist, :, :, :);
%         elseif eDist_new.ancillary.esteptable(it) == 1
%             eDist_new_tmp = eDist_new.data(it, :, :, :) - ephoto_scale * ebgdist1.data(iebgdist, :, :, :);
%         end
        c_eval('eDist_new_tmp = eDist_new.data(it, :, :, :) - ephoto_scale * ebgdist?.data(iebgdist, :, :, :);', eDist_new.ancillary.esteptable(it)); 
        eDist_new_tmp(eDist_new_tmp <= 0) = 0;              % may consider 'NaN';
        eDist_new.data(it, :, :, :) = eDist_new_tmp;
    end
    toc;
    
end
%%