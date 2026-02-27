%% Example initalizing the solo database
% The data from nas24 is mounted under /Volumes/solo

solo.db_init('local_file_db','/Volumes/solo/');
solo.db_init('local_file_db','/Volumes/solo/data_irfu');
solo.db_init('db_cache_size_max',4096)
solo.db_cache('on','save')

