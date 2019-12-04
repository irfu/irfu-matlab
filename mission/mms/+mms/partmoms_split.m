function res = partmoms_split(partvar)
% split FPI partmoms data into a struct
%
%   Example:
%       1. PartNi? = mms.partmoms_split(partNi);
%   History:
%       1. + on 20190910; wyli

%%  1. check
    tmp = partvar;
    
%%  2. make data 
    for ienergy = 1: 32
        if strfind(partvar.name, 'numberdensity')
            c_eval('restmp? = irf.ts_scalar(tmp.time, tmp.data(:, ienergy));', ienergy);
        elseif strfind(partvar.name, 'bulkv')
            c_eval('restmp? = irf.ts_vec_xyz(tmp.time, squeeze(tmp.data(:, ienergy, :)));', ienergy);
        elseif or(~isempty(strfind(partvar.name, 'temptensor')), ~isempty(strfind(partvar.name, 'prestensor')))
            c_eval('restmp? = irf.ts_tensor_xyz(tmp.time, squeeze(tmp.data(:, :, ienergy, :)));', ienergy);
        else
            error('Not part-moms data. ');
        end
            c_eval('restmp?.userData = tmp.userData;', ienergy);
            c_eval('restmp?.name = tmp.name;', ienergy);
            c_eval('restmp?.units = tmp.units;', ienergy);
            c_eval('restmp?.siConversion = tmp.siConversion;', ienergy);            
    end 
    
%%  2. make struct    
    res = struct('part1', restmp1, 'part2', restmp2, 'part3', restmp3, 'part4', restmp4, ...
        'part5', restmp5, 'part6', restmp6, 'part7', restmp7, 'part8', restmp8, ...
        'part9', restmp9, 'part10', restmp10, 'part11', restmp11, 'part12', restmp12, ...
        'part13', restmp13, 'part14', restmp14, 'part15', restmp15, 'part16', restmp16, ...
        'part17', restmp17, 'part18', restmp18, 'part19', restmp19, 'part20', restmp20, ...
        'part21', restmp21, 'part22', restmp22, 'part23', restmp23, 'part24', restmp24, ...
        'part25', restmp25, 'part26', restmp26, 'part27', restmp27, 'part28', restmp28, ...
        'part29', restmp29, 'part30', restmp30, 'part31', restmp31, 'part32', restmp32);
   
end