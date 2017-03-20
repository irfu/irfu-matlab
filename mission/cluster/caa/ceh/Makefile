
clean:
	rm -f *~ C[1-4]_CH_EFW_L2_HK.ceh C[1-4]_CH_EFW_INST.ceh \
	C[1-4]_CH_EFW_L[1-3]_P.ceh C[1-4]_CH_EFW_L[1-3]_E.ceh \
	C[1-4]_CH_EFW_L3_SFIT.ceh C[1-4]_CH_EFW_L3_DER.ceh \
	C[1-4]_CH_EFW_L1_P[1-4].ceh C[1-4]_CH_EFW_L1_P[1,3][2,4].ceh \
	C[1-4]_CQ_EFW_INST.ceh C[1-4]_CH_EFW_L1_IB.ceh C[1-4]_CH_EFW_L2_PB.ceh \
	C[1-4]_CH_EFW_L2_EB.ceh C[1-4]_CH_EFW_L2_BB.ceh \
	C[1-4]_CH_AUX_MAARBLE_ULF_PC1.ceh \
	C[1-4]_CH_AUX_MAARBLE_ULF_PC12.ceh C[1-4]_CH_AUX_MAARBLE_ULF_PC35.ceh \
	C[1-4]_CH_AUX_MAARBLE_ULF_FACMATR.ceh C[1-4]_CH_AUX_MAARBLE_SA_VLF.ceh \
        CC_CH_AUX_MAARBLE_TH[A-E]_ULF_PC12.ceh \
        CC_CH_AUX_MAARBLE_TH[A-E]_ULF_PC35.ceh \
        CC_CH_AUX_MAARBLE_TH[A-E]_ULF_FACMATR.ceh \
        CC_CH_AUX_MAARBLE_G1[1-2]_ULF_PC12.ceh \
        CC_CH_AUX_MAARBLE_G1[1-2]_ULF_PC35.ceh \
        CC_CH_AUX_MAARBLE_G1[1-2]_ULF_FACMATR.ceh \
	C[1-4]_CP_EFW_*_V00.cef

upload-maarble-irf:
	ls *.ceh| grep -v CX | grep -v THZ | grep -v GZZ | xargs -I % cp % /data/caa/MAARBLE/Upload/irf/HEADERS/

upload-maarble-irf-scp:
	ls *.ceh| grep -v CX | grep -v THZ | grep -v GZZ | xargs -I % scp % yuri@db.irfu.se:/data/caa/MAARBLE/Upload/irf/HEADERS/

upload-maarble-iap:
	ls *.ceh| grep -v CX | grep -v ULF_PC12 | grep -v ULF_PC35 | grep -v _TH | grep -v EFW | grep -v FACMATR | xargs -I % cp % /data/caa/MAARBLE/Upload/iap/HEADERS/

upload-maarble-caa:
	scp C[1-4]_CH_AUX_MAARBLE_*.ceh maarble@caa-delivery.estec.esa.int:HEADERS/ 

update-dropbox:
	ls Makefile make_maarble_headers.sh *.ceh| xargs -I % cp % /Users/yuri/Dropbox/Projects/MAARBLE/WP3/HeadersAndExamples/Headers/
