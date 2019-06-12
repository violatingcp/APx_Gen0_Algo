#ifndef __SORTNING_NETWORK_H__
#define __SORTNING_NETWORK_H__

#include "GlobalCorrelator_HLS/firmware/simple_fullpfalgo.h"
#include "hls_stream.h"

#define DATA_SIZE 128
#define NTAU 3
#define NTAUPARTS 10

void swap1(PFChargedObj &data1,PFChargedObj &data2);
void swap2(PFChargedObj &data1,PFChargedObj &data2);
void sorting_network_64_in(PFChargedObj datas[DATA_SIZE]);
void sorting_network_64(hls::stream<PFChargedObj > data_in[], hls::stream<PFChargedObj > data_out[]);
void sorting_network_128_in(PFChargedObj datas[DATA_SIZE]);
void sorting_network_128(hls::stream<PFChargedObj > data_in[], hls::stream<PFChargedObj > data_out[]);
void sorting_network_128_nn(hls::stream<PFChargedObj> data_in[DATA_SIZE],hls::stream<PFChargedObj> &data_out);
#endif
