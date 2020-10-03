# Mirroring public PSP data locally

% FIELDS data (the example downloads 2020 data)
wget -c -N -e robots=off --recursive --wait=2 -A "*_2020*" --no-parent --level=inf --reject-regex '[\?]' http://research.ssl.berkeley.edu/data/psp/data/sci/fields/l2

% SWEAP data (the example downloads 2020 data)
wget -c -N -e robots=off --recursive --wait=2 -A "*_2020*" --no-parent --level=inf --reject-regex '[\?]' http://sweap.cfa.harvard.edu/pub/data/sci/

See also https://gist.github.com/pulupa/fd1b6810091a347fccaf1657cadefe39

% Creating directory with links to all files
mkdir allcdfs
rm allcdfs/*
cd allcdfs
find .. -type f -name "*.cdf" -exec ln -s "{}" ./ ';' 

