%Documentation for c_ri
%
%
%
%
%
%
%
% see also c_ri_run_all
%
%##### Top function #####
%c_ri_run_all
%This function 7 of the 9 main functions. Goes all the way from
%finding MP-crossings to generating plots of the events. 
%
%using:
% c_ri_auto_event_search
% c_ri_run_get_B
% c_ri_four_B_files_2_one
% c_ri_load_B_for_preprocess
% c_ri_run_calc_angles_w_pre
% c_ri_run_class_angle_as_event
% c_ri_run_events_into_pictures
%
%##### main functions #####
% c_ri_auto_event_search
% c_ri_run_get_B
% c_ri_four_B_files_2_one
% c_ri_load_B_for_preprocess
% c_ri_run_calc_angles_w_pre
% c_ri_run_class_angle_as_event
% c_ri_run_events_into_pictures
%
% download_B -(using ISDAT)
% calc_angles.m -(can calculate the angles for assymetric B)
% view_and_filter_angles 
%
%#### lower level functions #####
%
% used by: c_ri_auto_event_search
% isGetContentLite
% toepoch
% search_events
%
% used by: search_events
% get_pos_from_ISDAT
% distance_to_MP
% gse2gsm
%
%used by: get_B_from_ISDAT
%Mat_DbOpen
%database UNIX:99 !not a function!
%isGetDataLite
%R_c_despin
%R_c_gse2dsc
%Mat_DbClose
%create_downloadblocks
%
%used by: get_pos_from_ISDAT
%Mat_DbOpen
%isGetDataLite
%Mat_DbClose
%database UNIX:99 !not a function!
%
%used by: search_events
%get_pos_from_ISDAT
%distance_to_MP
%gse2gsm
%
%function c_ri_auto_event_search
%isGetContentLite
%toepoch
%search_events
%save /share/robert/mp_crossing_(from date)_to_(to_date) 
%
%used by: download_B
%isGetContetLite
%get_sample_fq_for_period
%t_and_dt
%download_B_4_cl
%R_datestring
%load (filename and path to a file that contains the crossing of the MP)
%save /share/robert/B_data/B_[from]_to_[to]
%
%used by: download_B_intervall
%get_B_from_ISDAT
%add_A2M
%
%used by: calc_angles
%fromepoch
%Ut2number
%toepoch
%time_synch
%four_vector_angles
%saves /share/robert/angle/A_[from]_to_[to]
%
%used by:view_and_filter_angles
%fromepoch
%toepoch
%time2row
%find_max_angles
%plot_max_angles 
%class_angle_as_event
%
%function c_riget_B
%using	-ddscut
%		-fgmtel
%		-fgmcal
%		-fgmhrt
%		-fgmvec
%
%function c_ri_get_many_B
%using:	-create_timetable
%		-c_ri_get_B
%		-datestring
%
%function c_ri_four_B_files_2_one
%using:	-find_str
%		-load_file
%
%function load_file
%using:	-timestr2epoch
%		-hhmmss2epoch
%
%function load_B_for_preprocess
%using:	-preprocess
%
%function c_ri_calc_angle_w_preprocess
%using:	-c_ri_angles
%
%function c_ri_run_get_B
%using:	-c_ri_get_many_B
%
%function c_ri_run_class_angle_as_event
%Using:
% c_ri_timestr_within_intervall
% class_angle_as_event
%
%fucntion c_ri_run_event_to_picture
%using: c_ri_event_to_picture
%
%##### List of functions #####
%These functions are active in the program and are used by one of the
%main functions
%c_ri_angles.m
%c_ri_angles_and_ampl.m
%c_ri_auto_event_search.m
%c_ri_calc_angle_w_preprocess.m
%c_ri_calc_dt.m
%c_ri_comple.m
%c_ri_events_into_pictures.m
%c_ri_four_B_files_2_one.m
%c_ri_get_B.m
%c_ri_get_many_B.m
%c_ri_load_B_for_preprocess.m
%c_ri_run_all.m
%c_ri_run_calc_angles.m
%c_ri_run_calc_angles_w_pre.m
%c_ri_run_class_angle_as_event.m
%c_ri_run_events_into_pictures.m
%c_ri_run_get_B.m
%c_ri_timestr_within_intervall.m
%c_ri_timestr_within_intervall_E.m
%c_ri_timestr_within_intervall_MP.m
%c_run_events_into_pictures.m
%calc_angles.m
%class_angle_as_event.m
%create_download_blocks.m
%create_file.m
%create_timetable.m
%dist_to_MP_shue.m
%distance_to_MP.m
%download_B.m
%download_B_4_cl.m
%download_intervall.m
%find_max_angles.m
%find_row.m
%find_str.m
%find_vector_to_time.m
%four_vector_angles.m
%get_B_from_ISDAT.m
%get_event_data.m
%get_pos_from_ISDAT.m
%get_sample_fq_for_period.m
%get_timetable.m
%hhmmss2epoch.m
%ind2nr.m
%load_file.m
%plot_angles.m
%plot_max_angles.m
%preprocess_B.m
%search_events.m
%t_and_dt.m
%time2row.m
%time_synch.m
%timestr2epoch.m
%vector_angles.m
%view_and_filter_angles.m
