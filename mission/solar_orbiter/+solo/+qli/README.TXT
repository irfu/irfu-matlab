
QLI = QuickLooks IRFU ("IRFU quicklooks")

Code for generating SolO quicklooks available at https://www.space.irfu.se/soloql/.

NOTE: Not to be confused with LESIA's summary plots or BIA_QL3.
NOTE: These quicklooks use both SolO/RPW datasets AND datasets from other instruments.



=================================
 How to generate quicklook files
=================================
For an example of how to generate quicklooks locally, e.g. in bulk, look at below file:
    generate_quicklooks_demo.m



===========================================
 How to edit the content of the quicklooks
===========================================
generate_quicklooks_24h_6h_2h.m
generate_quicklook_7days.m
    Files which do the actual plotting.

generate_quicklooks_24h_6h_2h_using_DB_SPICE.m
generate_quicklook_7days_using_DB_SPICE.m
    Files which call generate_quicklooks_24h_6h_2h.m and generate_quicklook_7days.m
    and supply them with bulk data. If you want the quicklooks to have access
    to more data e.g. science data (zVariables), then you need to edit these.
