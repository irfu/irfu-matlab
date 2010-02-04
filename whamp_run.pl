#! perl -w
#
print(" \n Running Whamp in a bit stupid way, anybody wants to make proper solution?:) \n\nInput model filename\n");
$whamp_home="./";
open(WHAMP,"$whamp_home/whamp_cygwin |");
$file_save="$whamp_home/temp";
$file_wh="$whamp_home/wh";
@wh_columns = (1, 2, 3, 4);    # 1,2,3,4 should be first to use irf_whamp.m in matlab
#$flag_save=0;
#$flag_input=0;
@result=();
$iresult = 0;
$ifinal =0;
do{
	while(<WHAMP>){
		/^INPUT/ && do{
		$ifinal=$iresult-1;
		$iresult=0; 		# reset memory	
		if ($ifinal > 10) {save();};
		};
		print $_;
		$result[$iresult]=$_;
		$iresult ++;
	}
	print("hej\n");
} until 1;

sub save
{
print("--> ",$ifinal," lines saved to file: ",$file_save,"\n");
open(TEMP,">$file_save");
open(WH,">$file_wh");
for ($i=1;$i<=$ifinal;$i++) {
	print TEMP ($result[$i]);
	@res = split(' ',$result[$i]);	
	if ($#res>=$#wh_columns){
		foreach $wh_col (@wh_columns)
		{
			print WH $res[$wh_col-1],"\t";
		}
		print WH "\n";
	}
}
close(TEMP);
close(WH)

}


