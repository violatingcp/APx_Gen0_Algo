## Set the top level module
set_top algo_unpacked
##
#### Add source code
add_files ${PROJ_DIR}/src/algo_unpacked.cpp
#
### Add testbed files
add_files -tb ${PROJ_DIR}/src/algo_unpacked_tb.cpp -cflags ${CFLAGS}

### Add test input files
add_files -tb ${PROJ_DIR}/data/test1_inp.txt
add_files -tb ${PROJ_DIR}/data/test1_out_ref.txt
