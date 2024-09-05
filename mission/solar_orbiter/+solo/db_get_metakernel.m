function db_get_metakernel(flown_or_predicted)
% SOLO.DB_GET_METAKERNEL  Load a local SolO SPICE metakernel
%
% solo.db_get_metakernel(flown_or_predicted)
%
% Function to load a locally adapted SPICE metakernel for SolO (assuming
% folder "SPICE" is a subfolder of the database root).
%
% flown_or_predicted = 'flown' return metakernel of actually
%                       flown orbit only
% flown_or_predicted = 'predicted' return metakernel of latest orbit
%                      prediction (while still using reconstructed flown
%                      orbit information as much as possible).
narginchk(1,2)
global SOLO_DB; if isempty(SOLO_DB), solo.db_init(), end

SOLO_DB.get_metakernel(flown_or_predicted);
end
