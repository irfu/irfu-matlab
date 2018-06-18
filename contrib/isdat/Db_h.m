%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure like implementation in matlab
% Db_h
% include the Dbdef_h.m file and Isdef_h.m ( definition of DbMAX_DIMS and others)
% Warning: these structures fit the ISDAT 1.6 release.
% PG 96/04/26
 

Dbdef_h
Isdef_h


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbDataSpec Index

data_spec_project	= 1;
data_spec_member	= 2;
data_spec_instrument	= 3;
data_spec_sensor	= 4;
data_spec_signal	= 5;
data_spec_channel	= 6;
data_spec_parameter	= 7;

DbDataSpecIndex=[data_spec_project data_spec_member data_spec_instrument...
 data_spec_sensor data_spec_signal data_spec_channel data_spec_parameter];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbSpecName Index

spec_name_project	= 1:DbSPEC_NAME_DIM;
spec_name_member	= max(spec_name_project)	+ 1:DbSPEC_NAME_DIM;
spec_name_instrument	= max(spec_name_member)		+ 1:DbSPEC_NAME_DIM;
spec_name_sensor	= max(spec_name_instrument)	+ 1:DbSPEC_NAME_DIM;
spec_name_signal	= max(spec_name_sensor)		+ 1:DbSPEC_NAME_DIM;
spec_name_channel	= max(spec_name_signal)		+ 1:DbSPEC_NAME_DIM;
spec_name_parameter	= max(spec_name_channel)	+ 1:DbSPEC_NAME_DIM;

DbSpecNameIndex = [spec_name_project spec_name_member spec_name_instrument...
 spec_name_sensor spec_name_signal spec_name_channel spec_name_parameter];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbData Index

data_RF0_re	= 1;			DbDataRF0Index	= data_RF0_re;
data_RF1_re	= 1:3;		DbDataRF1Index	= data_RF1_re;
data_RF2_re	= 1:6;		DbDataRF2Index	= data_RF2_re;
data_RF3_re	= 1:9;		DbDataRF3Index	= data_RF3_re;

data_RF2D_re	= 1:2;		DbDataRF2DIndex	= data_RF2D_re;

data_RD0_re	= 1;			DbDataRD0Index	= data_RD0_re;
data_RD1_re	= 1:3;		DbDataRD1Index	= data_RD1_re;
data_RD2_re	= 1:6;		DbDataRD2Index	= data_RD2_re;

data_CF0_re	= 1;
data_CF0_im	= 2;			DbDataCF0Index	= [data_CF0_re data_CF0_im];
data_CF1_re	= 1:3;
data_CF1_im	= 4:6;		DbDataCF1Index	= [data_CF1_re data_CF1_im];
data_CF2_re	= 1:6;
data_CF2_im	= 7:12;		DbDataCF2Index	= [data_CF2_re data_CF2_im];
data_CF3_re	= 1:9;
data_CF3_im	= 10:18;		DbDataCF3Index	= [data_CF3_re data_CF3_im];

data_CF2D_re	= 1:6;
data_CF2D_im	= 7:12;		DbDataCF2DIndex	= [data_CF2D_re data_CF2D_im];
data_CD0_re	= 1;
data_CD0_im	= 2;			DbDataCD0Index	= [data_CD0_re data_CD0_im];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbDataInfo Index

data_info_units		 = 1;
data_info_quantity	 = max (data_info_units		)	+ 1;
data_info_scaleType	 = max (data_info_quantity	)	+ 1;
data_info_scaleMin	 = max (data_info_scaleType	)	+ 1;
data_info_scaleMax	 = max (data_info_scaleMin	)	+ 1;
data_info_samplingFreq	 = max (data_info_scaleMax	)	+ 1;
data_info_filterFreq	 = max (data_info_samplingFreq	)	+ 1;
data_info_unitString	 = max (data_info_filterFreq	)	+ 1:32;
data_info_quantityString = max (data_info_unitString	)	+ 1:32;
data_info_conversion	 = max (data_info_quantityString)	+ 1:80;

DbDataInfoIndex = [ data_info_units data_info_quantity data_info_scaleType ...
		    data_info_scaleMin data_info_scaleMax data_info_samplingFreq ...
		    data_info_filterFreq data_info_unitString ...
		    data_info_quantityString data_info_conversion ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbCoordinate Index

coordinate_reference	= 1;
coordinate_system	= 2;
coordinate_rot		= max(coordinate_system)	+ DbDataRD2Index;
coordinate_location	= max(coordinate_rot)		+ DbDataRD1Index;
coordinate_direction	= max(coordinate_location)	+ DbDataRD1Index;

DbCoordinateIndex =[ coordinate_reference coordinate_system coordinate_rot	...
			coordinate_location coordinate_direction ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbDataRequest Index

data_request_start	= IsTimeIndex;
data_request_interval	= max(data_request_start)	+ IsTimeIndex;
data_request_spec	= max(data_request_interval)	+ DbDataSpecIndex;
data_request_units	= max(data_request_spec)	+ 1;
data_request_reduction	= max(data_request_units)	+ 1;
data_request_samples	= max(data_request_reduction)	+ 1;
data_request_gapFill	= max(data_request_samples)	+ 1;
data_request_pack	= max(data_request_gapFill)	+ 1;
data_request_dataVersion= max(data_request_pack)	+ 1;

DbDataRequestIndex = [ data_request_start data_request_interval data_request_spec ...
data_request_units data_request_reduction data_request_samples data_request_gapFill ...
data_request_pack data_request_dataVersion ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbQuery_desc Index

query_desc_mode		= 1;
query_desc_category	= max(query_desc_mode)		+ 1;
query_desc_level	= max(query_desc_category)	+ 1;
query_desc_time		= max(query_desc_level)		+ IsTimeIndex;
query_desc_spec		= max(query_desc_time)		+ DbDataSpecIndex;

DbQueryDescIndex = [ query_desc_mode query_desc_category query_desc_level	...
			query_desc_time query_desc_spec ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbQuery_data Index

query_data_value	= 1;
query_data_groupId	= max(query_data_value)		+ 1;
query_data_event	= max(query_data_groupId)	+ 1;
query_data_name		= max(query_data_event)		+ 1:32;
query_data_rank		= max(query_data_name)		+ 1;
query_data_complete	= max(query_data_rank)		+ 1;
query_data_dataType	= max(query_data_complete)	+ 1;
query_data_dimension	= max(query_data_dataType)	+ 1;
query_data_category	= max(query_data_dimension)	+ 1;
query_data_mapType	= max(query_data_category)	+ 1;

DbQueryDataIndex = [ query_data_value query_data_groupId query_data_event query_data_name ...
query_data_rank	query_data_complete query_data_dataType query_data_dimension ...
query_data_category query_data_mapType ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbContent_desc Index

content_desc_spec		= DbDataSpecIndex;
content_desc_dataVersion	= max(content_desc_spec)	+ 1;
content_desc_sections		= max(content_desc_dataVersion)	+ 1;

DbContentDescIndex = [ content_desc_spec content_desc_dataVersion content_desc_sections ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbOverview_desc Index

overview_desc_start		= IsTimeIndex;
overview_desc_interval		= max(overview_desc_start)		+ IsTimeIndex;
overview_desc_spec		= max(overview_desc_interval)		+ DbDataSpecIndex;
overview_desc_event		= max(overview_desc_spec)		+ 1;
overview_desc_dataVersion	= max(overview_desc_event)		+ 1;
overview_desc_sections		= max(overview_desc_dataVersion)	+ 1;

DbOverviewDescIndex = [ overview_desc_start overview_desc_interval overview_desc_spec ...
		overview_desc_event overview_desc_dataVersion overview_desc_sections ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbPrepareDataDesc Index

prepare_desc_start	= IsTimeIndex;
prepare_desc_interval	= max(prepare_desc_start)	+ IsTimeIndex;
prepare_desc_spec	= max(prepare_desc_interval)	+ DbDataSpecIndex;

DbPrepareDescIndex = [ prepare_desc_start prepare_desc_interval prepare_desc_spec ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbControlDesc Index

control_desc_function		= 1;
control_desc_spec		= max(control_desc_function)	+ DbDataSpecIndex;

DbControlDescIndex = [ control_desc_function control_desc_spec ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DbInfoDesc Index

info_desc_spec			= DbDataSpecIndex;

DbInfoDescIndex = info_desc_spec;
