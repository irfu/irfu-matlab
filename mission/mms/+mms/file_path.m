
function filepath = file_path(varargin)
% FILE_PATH make the filepath of FPI/FGM/EDP data
%
%   filepath = file_path(ic, interval, 'harddrive', harddrive, 'fgm_mode', 'srvy', 'edp_mode', 'brst');
%
%   Input:
%   ic - spacecraft #
%   interval - data file time: "YYYYMMDDhhmmss"
%   Optional Inputs:
%   'harddrive' - data harddrive; DEFAULT: '/Volumes/mms';
%   'fpi_mode' - FPI mode: brst [default]/fast; [under construction];
%   'fgm_mode' - B field data mode: brst [default]/srvy;
%   'edp_mode' - E field data mode: brst [default]/fast;
%   'scpot_mode' - sc potential data mode: brst [default]/fast
%   Output: 
%   filepath - structure containing filepath;
%   Potential bug: fake version selection.

    %   0. input check
    if (nargin < 2),
        nargin
        help filepath;
        return;
    end
    ic = varargin{1};
    interval = varargin{2};
    % default setting;
    harddrive = '/Volumes/mms';
    fpi_mode = 'brst';
    fgm_mode = 'brst';
    edp_mode = 'brst';
    scpot_mode = 'brst';
    % optional input
    args = varargin(3: end);
    if numel(args)> 0, 
        options = 1;
    else
        options = 0;
    end
    while options
    l = 2;
    switch(lower(args{1}))
        case 'harddrive'
            if numel(args) > 1 && ~isempty(args{2}),
                harddrive = args{2};
            end
        case 'fpi_mode'
            if numel(args) > 1 && ~isempty(args{2}),
                fpi_mode = args{2};
                if strcmp(fpi_mode, 'fast'), error('Under construction!'); end
            end
        case 'fgm_mode'
            if numel(args) > 1 && ~isempty(args{2}),
                fgm_mode = args{2};
            end
        case 'edp_mode'
            if numel(args) > 1 && ~isempty(args{2}),
                edp_mode = args{2};
            end            
        case 'scpot_mode'
            if numel(args) > 1 && ~isempty(args{2}),
                scpot_mode = args{2};
            end
        otherwise
            irf.log('critical',['Unknown flag: ' args{1}]);
            l=1;
            break
    end
    args = args(l+1:end);
    if isempty(args), options=0; end
    end
    
    %   1. time variables
    yyyy = interval(1:4);
    mm = interval(5:6);
    dd = interval(7:8);

    %   2. FPI direction and filename, NOW ONLY brst
    str_fpi = {'idist', 'edist', 'imoms', 'emoms'};
    c_eval('idist_dr = [harddrive,''/mms?/fpi/brst/l2/dis-dist/'' yyyy ''/'' mm ''/'' dd ''/''];', ic);
    c_eval('idist_fn = [''mms?_fpi_brst_l2_dis-dist_'' interval ''_v*.cdf''];', ic);
    c_eval('edist_dr = [harddrive,''/mms?/fpi/brst/l2/des-dist/'' yyyy ''/'' mm ''/'' dd ''/''];', ic);
    c_eval('edist_fn = [''mms?_fpi_brst_l2_des-dist_'' interval ''_v*.cdf''];', ic);    
    c_eval('imoms_dr = [harddrive,''/mms?/fpi/brst/l2/dis-moms/'' yyyy ''/'' mm ''/'' dd ''/''];', ic);
    c_eval('imoms_fn = [''mms?_fpi_brst_l2_dis-moms_'' interval ''_v*.cdf''];', ic);
    c_eval('emoms_dr = [harddrive,''/mms?/fpi/brst/l2/des-moms/'' yyyy ''/'' mm ''/'' dd ''/''];', ic);
    c_eval('emoms_fn = [''mms?_fpi_brst_l2_des-moms_'' interval ''_v*.cdf''];', ic);
    c_eval('str_? = dir([?_dr ?_fn]);', str_fpi);
    if isempty(str_idist) || isempty(str_edist) || isempty(str_imoms) || isempty(str_emoms), ...
            error('Check connection of internet and harddrive and existence of this data'); end
    c_eval('?_fn = str_?(end).name;', str_fpi);
    for i=1: 4
        c_eval('tmp = str_?;', str_fpi(i));   
        if length(tmp) > 1
            disp(tmp(:).name);
            c_eval('disp(?_fn);', str_fpi(i));
            
        end
    end

    % 3. FGM direction and filename: brst [default]/srvy
    if strcmp(fgm_mode, 'brst')
        c_eval('fgm_dr = [harddrive,''/mms?/fgm/brst/l2/'' yyyy ''/'' mm ''/'' dd ''/''];', ic);
        c_eval('fgm_fn = [''mms?_fgm_brst_l2_'' interval ''_v*.cdf''];', ic);            
        str_fgm = dir([fgm_dr fgm_fn]);
        if isempty(str_fgm), error('Check connection of internet and harddrive and existence of this data'); end
        fgm_fn = str_fgm(end).name;
        if length(str_fgm) > 1
            disp(str_fgm(:).name);
            disp(fgm_fn);
        end
    else if strcmp(fgm_mode, 'srvy')
        c_eval('fgm_dr = [harddrive,''/mms?/fgm/srvy/l2/'' yyyy ''/'' mm ''/''];', ic);
        c_eval('fgm_fn = [''mms?_fgm_srvy_l2_'' interval(1:8) ''_v*.cdf''];', ic);            
        str_fgm = dir([fgm_dr fgm_fn]);
        if isempty(str_fgm), error('Check connection of internet and harddrive and existence of this data'); end
        fgm_fn = str_fgm(end).name;
        if length(str_fgm) > 1
            disp(str_fgm(:).name);
            disp(fgm_fn);
        end
        else
            error('FGM mode error');
        end
    end
    
    % 4. edp direction and filename: brst [default]/fast;
    if strcmp(edp_mode, 'brst')
        c_eval('edp_dr = [harddrive,''/mms?/edp/brst/l2/dce/'' yyyy ''/'' mm ''/'' dd ''/''];', ic);
        c_eval('edp_fn = [''mms?_edp_brst_l2_dce_'' interval ''_v*.cdf''];', ic);            
        str_edp = dir([edp_dr edp_fn]);
        if isempty(str_edp), error('Check connection of internet and harddrive and existence of this data'); end
        edp_fn = str_edp(end).name;
        if length(str_edp) > 1
            disp(str_edp(:).name);
            disp(edp_fn);
        end
    else if strcmp(edp_mode, 'fast')
        c_eval('edp_dr = [harddrive,''/mms?/edp/fast/l2/dce/'' yyyy ''/'' mm ''/''];', ic);
        c_eval('edp_fn = [''mms?_edp_fast_l2_dce_'' interval(1:8) ''_v*.cdf''];', ic);            
        str_edp = dir([edp_dr edp_fn]);
        if isempty(str_edp), error('Check connection of internet and harddrive and existence of this data'); end
        edp_fn = str_edp(end).name;
        if length(str_edp) > 1
            disp(str_edp(:).name);
            disp(edp_fn);
        end
        else
            error('EDP mode error');
        end
    end

    % 5. scpot direction and filename: brst [default]/fast;    
    if strcmp(scpot_mode, 'brst')
        c_eval('scpot_dr = [harddrive,''/mms?/edp/brst/l2/scpot/'' yyyy ''/'' mm ''/'' dd ''/''];', ic);
        c_eval('scpot_fn = [''mms?_edp_brst_l2_scpot_'' interval ''_v*.cdf''];', ic);            
        str_scpot = dir([scpot_dr scpot_fn]);
        if isempty(str_scpot), error('Check connection of internet and harddrive and existence of this data'); end
        scpot_fn = str_scpot(end).name;
        if length(str_scpot) > 1
            disp(str_scpot(:).name);
            disp(scpot_fn);
        end
    else if strcmp(scpot_mode, 'fast')
        c_eval('scpot_dr = [harddrive,''/mms?/edp/fast/l2/scpot/'' yyyy ''/'' mm ''/''];', ic);
        c_eval('scpot_fn = [''mms?_edp_fast_l2_scpot_'' interval(1:8) ''000000_v*.cdf''];', ic);            
        str_scpot = dir([scpot_dr scpot_fn]);
        if isempty(str_scpot), error('Check connection of internet and harddrive and existence of this data'); end
        scpot_fn = str_scpot(end).name;
        if length(str_scpot) > 1
            disp(str_scpot(:).name);
            disp(scpot_fn);
        end
        else
            error('EDP/scpot mode error');
        end
    end    
    
    % 6. output structure
    filepath =struct('idist', [idist_dr idist_fn], 'edist', [edist_dr edist_fn], ...
        'imoms', [imoms_dr imoms_fn], 'emoms', [emoms_dr emoms_fn], ...
        'fgm', [fgm_dr fgm_fn], 'edp', [edp_dr edp_fn], 'scpot',[scpot_dr scpot_fn]);
    
end
%% updated on 2016-05-03; irfu.