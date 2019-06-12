#ifndef ALGO_INPUTS_LAYER2_V2_H
#define ALGO_INPUTS_LAYER2_V2_H

#include <ap_int.h>
#include "GlobalCorrelator_HLS/firmware/simple_fullpfalgo.h"
#include "GlobalCorrelator_HLS/puppi/firmware/simple_puppi.h"

#define DATA_SIZE 32
#define NTAU  6
#define NREGIONS 36
#define NPART 25
#define DEPTH NREGIONS*2
#define NTAUPARTS  2
#define DRCONE 8410
#define DR2MAX 10000
#define EMOFFS 2*NTRACK
#define HAOFFS 2*NEMCALO+EMOFFS
#define MUOFFS 2*NCALO+HAOFFS
#define MP7_NCHANN 72
typedef ap_uint<32> MP7DataWord;

static float PT_SCALE = 4.0;     // quantize in units of 0.25 GeV (can be changed)
static float ETAPHI_FACTOR = 4;  // size of an ecal crystal in phi in integer units (our choice)
static float ETAPHI_SCALE = ETAPHI_FACTOR*(180./M_PI);  // M_PI/180 is the size of an ECal crystal; we make a grid that is 4 times that size
static int16_t PHI_WRAP = 360*ETAPHI_FACTOR;            // what is 3.14 in integer

void algo_inputs_layer2_v3(MP7DataWord input[MP7_NCHANN],MP7DataWord output[MP7_NCHANN]);

#endif
