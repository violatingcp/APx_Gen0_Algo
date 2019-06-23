## Set the top level module
set_top algo_layer2_v3
##
#### Add source code
add_files src/algo_layer2_v3.cpp 
add_files src/sorting_network_axi.cpp
add_files src/firmware/tau_nn.cpp
add_files src/GlobalCorrelator_HLS/firmware/simple_fullpfalgo.cpp -cflags "-DTESTMP7 -DHLS_pipeline_II=3"
add_files src/GlobalCorrelator_HLS/puppi/firmware/simple_puppi.cpp -cflags "-DTESTMP7 -DHLS_pipeline_II=3"
#
### Add testbed files
add_files -tb src/algo_layer2_v3_tb.cpp 


