function [eDist_new, eDist_bg, ephoto_scale] = remove_edist_background(varargin)
% MMS.REMOVE_EDIST_BACKGROUND: remove secondary photoelectrons from electron distribution function
%
% [eDist_new] = MMS.REMOVE_EDIST_BACKGROUND(eDist, Tint)
%
% Input: 
%       eDist - Electron distribution as a TSeries
%
% Options:
%   'tint' - time interval; DEFAULT: time interval of input eDist;
%   'Nphotoe_art' - artificial photo electron density; DEFAULT: ephoto_scale from des-emoms GlobalAttributes data;
%   'ZeroNaN' - set the negative; DEFAULT: 0; [0/NaN]
% 
% References:
%   1. MMS FPI Data Users Guide: 
%       https://lasp.colorado.edu/galaxy/display/MFDPG/DES+Photoelectrons+-+further+details
% 
% Next step: move ion distribution background noise
% 
%   Example:     
%     [eDist1_bgremoved, eDist1_bg, ephoto_scale] = remove_edist_background(eDist1, 'tint', Tint, 'Nphotoe_art', 1.0, 'ZeroNaN', NaN);
%     [eDist1_bgremoved, eDist1_bg, ephoto_scale] = remove_edist_background(eDist1, 'tint', Tint);
%     [eDist1_bgremoved, eDist1_bg, ephoto_scale] = remove_edist_background(eDist1);
%
% 2018-09-28, wyli;

    %%  1. basic
    Nphotoe_art = -1;
    Nphoto = 0;
    ZeroNaN = 0;
    [~, args, nargs] = axescheck(varargin{:});
    eDist = args{1};
    Tint = irf.tint(eDist.time(1), eDist.time(end));    
    args_tmp = args(2: end);
    if nargs > 1, have_options = 1; end
    while have_options 
        l = 1;
       switch(lower(args_tmp{1}))
            case 'tint'
                l = 2;
                Tint = args_tmp{2};
            case 'nphotoe_art'
                l = 2;
                Nphotoe_art = args_tmp{2};
            case 'zeronan'
                l = 2;
                ZeroNaN = args_tmp{2};                
        end
        args_tmp = args_tmp(l+1:end);
        if isempty(args_tmp), break, end  
    end
    ic = eDist.name(4);
    eDist_tmp = eDist.tlim(Tint);
    eDist_new = eDist_tmp;
    eDist_bg = eDist_tmp;
    
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
    energy0 = get_variable(ebgdist_obj, 'mms_des_energy0_brst');        % not used
    energy1 = get_variable(ebgdist_obj, 'mms_des_energy1_brst');        % not used
    phi = get_variable(ebgdist_obj, 'mms_des_phi_brst');
    theta = get_variable(ebgdist_obj, 'mms_des_phi_brst');
    
    %%  4. remove ebgdist
    % 4.1. set photo electron density
    if Nphotoe_art > 0
        Nphoto = Nphotoe_art;
    else
        Nphoto = ephoto_scale;
    end    
    % 4.2. make bg removed eDist & bg edist
    ntime = length(eDist_tmp.time);
    edist_size = size(eDist_new.data);
    eDist_new_tmp = zeros(1, edist_size(2), edist_size(3), edist_size(4));          % preallocating for speed
    eDist_bg_tmp = zeros(1, edist_size(2), edist_size(3), edist_size(4));           % preallocating for speed
    tic;
    for it = 1: ntime
        iestartdelphi_count = double(estartdelphi_count_Tint.data(it));
        iestartdelphi_angle = double(estartdelphi_angle_Tint.data(it));
        iebgdist = fix(iestartdelphi_count/16);
        if iebgdist == 0;       iebgdist = 360;     end  
        c_eval('eDist_new_tmp = eDist_new.data(it, :, :, :) - Nphoto * ebgdist?.data(iebgdist, :, :, :);', eDist_new.ancillary.esteptable(it)); 
        c_eval('eDist_bg_tmp = Nphoto * ebgdist?.data(iebgdist, :, :, :);', eDist_new.ancillary.esteptable(it)); 
        %c_eval('eDist_new_tmp = eDist_new.data(it, :, :, :) - 1 * ebgdist?.data(iebgdist, :, :, :);', eDist_new.ancillary.esteptable(it)); 
        eDist_new_tmp(eDist_new_tmp <= 0) = ZeroNaN;                  % may consider 'NaN';
        eDist_bg_tmp(eDist_bg_tmp <= 0) = ZeroNaN;                    % may consider 'NaN';
        eDist_new.data(it, :, :, :) = eDist_new_tmp;
        eDist_bg.data(it, :, :, :) = eDist_bg_tmp;        
    end
    toc;
    % 4.3. disp information
    irf.log('warning', ['* Secondary photoelectron model used: ''../mms/models/fpi/' ebgdist_fn '''.']);
    irf.log('warning', ['* Secondary photoelectron scale (number density): ' num2str(ephoto_scale)]);
    if Nphotoe_art > 0, rf.log('warning', ['* Artificial Secondary photoelectron scale (number density): ' num2str(Nphotoe_art)]); end
    irf.log('warning', ['* Negative phase-space density set to : ' num2str(ZeroNaN)]);

end
%%