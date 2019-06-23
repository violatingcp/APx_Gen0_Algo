#ifndef __SORTNING_NETWORK_H__
#define __SORTNING_NETWORK_H__

//#include "GlobalCorrelator_HLS/firmware/simple_fullpfalgo.h"
#include "hls_stream.h"
#include <ap_int.h>
//#include "data.h"

typedef ap_uint<64> axi2_t;
#define NTAU 6
#define DATA_SIZE 128
#define NTAUPARTS 10

void swap1(axi2_t &data1,axi2_t &data2);
void swap2(axi2_t &data1,axi2_t &data2);
void sorting_network_64_in(axi2_t datas[DATA_SIZE]);
void sorting_network_64(hls::stream<axi2_t > data_in[], hls::stream<axi2_t > data_out[]);
void sorting_network_64_nn(axi2_t data_in[DATA_SIZE],axi2_t data_out[NTAU]);
void sorting_network_128_in(axi2_t datas[DATA_SIZE]);
void sorting_network_128(hls::stream<axi2_t > data_in[], hls::stream<axi2_t > data_out[]);
void sorting_network_128_nn(hls::stream<axi2_t> data_in[DATA_SIZE],hls::stream<axi2_t> &data_out);
#endif
