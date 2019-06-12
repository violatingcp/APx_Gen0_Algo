#include "bitonic_sort.hpp"
#include "sorting_network.hpp"


void input(PFChargedObj datas[DATA_SIZE], hls::stream<PFChargedObj> axis_in[DATA_SIZE]) {
	// input data from stream
	for (int i = 0; i < DATA_SIZE; ++i) {
	#pragma HLS unroll
		datas[i] = axis_in[i].read();
	}
}


void output(PFChargedObj datas[DATA_SIZE], hls::stream<PFChargedObj > axis_out[DATA_SIZE]) {
	// output data to stream
	for (int i = 0; i < DATA_SIZE; ++i) {
	#pragma HLS unroll
		axis_out[i].write(datas[i]);
	}
}

void bitonic_sort(hls::stream<PFChargedObj> data_in[DATA_SIZE],hls::stream<PFChargedObj> data_out[DATA_SIZE]) { 
  //void bitonic_sort(hls::stream<ap_int<DATA_W> > axis_in[DATA_SIZE], hls::stream<ap_int<DATA_W> > axis_out[DATA_SIZE]) {
#pragma HLS interface ap_ctrl_hs port=return
#pragma HLS interface axis port=data_in 
#pragma HLS interface axis port=data_out 
#pragma HLS PIPELINE II=25

  PFChargedObj datas[DATA_SIZE];
  #pragma HLS ARRAY_RESHAPE variable=datas complete dim=1
	
  for(int i0 = 0; i0 < 5; i0++) { 
    input(datas, data_in);
    sorting_network_128(datas);
    output(datas, data_out);
  }
}
