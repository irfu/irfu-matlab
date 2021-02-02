# Mirroring public PSP data locally

Goto directory where you want to mirror FIELDS data into directory "fields" and SWEAP data into directory "sweap".

% FIELDS data (the example downloads 2020 data) 
wget -c -N  --http-user=xxxx --http-password=xxxxxxxx -e robots=off --recursive --wait=2 -A "*_2020*_v0x*" --no-parent --level=inf --reject-regex '[\?]'  --cut-dirs=4 -nH  http://research.ssl.berkeley.edu/data/psp/data/sci/fields/l2

% SWEAP data (the example downloads 2020 data)
wget -c -N -e robots=off --recursive --wait=2 -A "*_2020*_v0x*" --no-parent --level=inf --reject-regex '[\?]'  --cut-dirs=3 -nH http://sweap.cfa.harvard.edu/pub/data/sci/

See also https://gist.github.com/pulupa/fd1b6810091a347fccaf1657cadefe39

% In both commands above, "_v0x" is the file's version and "x" should always be substituted by the latest number version (x = 1, 2 or 3) available for the cdf file. Otherwise all the versions will be fetched even though psp_load will always load the latest version.

% This works only if you are downloading an specific set of data, let's say 2020's dfb_wf_scm (v02) or dfb_ac_spec (v01), because different data products can be available in different versions. If you are making a full dowload of 2020 data, just delete the "_v0x*" part, keeping only "*_2020*". 

% The "x'es" in "--http-user=xxxx" and  "--http-password=xxxxxxx" must be substituted by the username and password for the data server. If you do not have an user and password, and the server is public, wget will override this flag and dowload the data anyway.

% For a full description of wget flags:
wget -h


