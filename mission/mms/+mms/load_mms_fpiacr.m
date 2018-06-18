
    function fpiacr = load_mms_fpiacr(ic, fpiacr_dr)
% * loading 37.5/7.5 ms FPI ACR data;
% How to use: 
%    ic = 1;
%    fpiacr_dr = '/.../';
%   fpiacr = load_mms_fpiacr(ic, fpiacr_dr);
%    idist = fpiacr.idist;       edist = fpiacr.edist;
%    ni = fpiacr.ni;             nee = fpiacr.nee;
%    gsevi = fpiacr.gsevi;       gseve = fpiacr.gseve;
%    tiperp = fpiacr.tiperp;     tipara = fpiacr.tipara;
%    teperp = fpiacr.teperp;     tepara = fpiacr.tepara;
% History: 
%   1. version 0.1 created on 2017-07-03; 
%

    % 1. basic
    c_eval('edist_fn = [''mms?_fpi_brst_acr_des-dist_*''];', ic); 
    c_eval('emoms_fn = [''mms?_fpi_brst_acr_des-moms_*''];', ic); 
    c_eval('idist_fn = [''mms?_fpi_brst_acr_dis-dist_*''];', ic); 
    c_eval('imoms_fn = [''mms?_fpi_brst_acr_dis-moms_*''];', ic); 
    str_edist = dir([fpiacr_dr edist_fn]);
    str_emoms = dir([fpiacr_dr emoms_fn]);
    str_idist = dir([fpiacr_dr idist_fn]);
    str_imoms = dir([fpiacr_dr imoms_fn]);    
    emoms_obj = dataobj([fpiacr_dr str_emoms(1).name]);
    imoms_obj = dataobj([fpiacr_dr str_imoms(1).name]);
    
    % 2. load dist data
    [edist, edisterr] = mms.make_pdist([fpiacr_dr str_edist(1).name]);
    [idist, idisterr] = mms.make_pdist([fpiacr_dr str_idist(1).name]);
    idist.data(idist.data <0) = 0;
    edist.data(edist.data <0) = 0;
    %edist.data = abs(edist.data);
    
    % 3. load ion moms data
    c_eval('ni = get_ts(imoms_obj, ''mms?_dis_numberdensity_brst'');', ic);
    c_eval('nierr = get_ts(imoms_obj, ''mms?_dis_numberdensity_err_brst'');', ic);
    %c_eval('niexlow = get_ts(imoms_obj, ''mms?_dis_densityextrapolation_low_brst'');', ic);
    %c_eval('niexhigh = get_ts(imoms_obj, ''mms?_dis_densityextrapolation_high_brst'');', ic);
    c_eval('gsevi = get_ts(imoms_obj, ''mms?_dis_bulkv_gse_brst'');', ic);
    c_eval('dbcsvi = get_ts(imoms_obj, ''mms?_dis_bulkv_dbcs_brst'');', ic);
    c_eval('tipara = get_ts(imoms_obj, ''mms?_dis_temppara_brst'');', ic);
    c_eval('tiperp = get_ts(imoms_obj, ''mms?_dis_tempperp_brst'');', ic);
    c_eval('dbcsPi33 = get_ts(imoms_obj, ''mms?_dis_prestensor_dbcs_brst'');', ic);
        
    % 4. load electron moms data
    c_eval('nee = get_ts(emoms_obj, ''mms?_des_numberdensity_brst'');', ic);
    c_eval('neerr = get_ts(emoms_obj, ''mms?_des_numberdensity_err_brst'');', ic);
    c_eval('gseve = get_ts(emoms_obj, ''mms?_des_bulkv_gse_brst'');', ic);
    c_eval('dbcsve = get_ts(emoms_obj, ''mms?_des_bulkv_dbcs_brst'');', ic);
    c_eval('tepara = get_ts(emoms_obj, ''mms?_des_temppara_brst'');', ic);
    c_eval('teperp = get_ts(emoms_obj, ''mms?_des_tempperp_brst'');', ic);   
    c_eval('dbcsPe33 = get_ts(emoms_obj, ''mms?_des_prestensor_dbcs_brst'');', ic);
    
    % 5. make output
    fpiacr = struct('idist', idist, 'edist', edist, 'idisterr', idisterr, ...
        'edisterr', edisterr, 'ni', ni, 'nierr', nierr, 'gsevi', gsevi, ...
        'dbcsvi', dbcsvi, 'tipara', tipara, 'tiperp', tiperp, 'dbcsPi33', dbcsPi33, ...
        'nee', nee, 'neerr', neerr, 'gseve', gseve, 'dbcsve', dbcsve, 'tepara', ...
        tepara, 'teperp', teperp, 'dbcsPe33', dbcsPe33);
    
    end
%%    