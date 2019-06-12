#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "ap_int.h"
#include "algo_inputs_layer2_v3.h"

using namespace std;

typedef ap_int<64>            tmpaxi_t;
int main(int argc, char ** argv) {
  float dphi  = 0.02; 
  float deta  = 0.02; 
  float pt    = 1.;
  float pterr = 0.1;
  tmpaxi_t chparts[DEPTH][NTRACK];  
  tmpaxi_t emparts[DEPTH][NEMCALO];  
  tmpaxi_t neparts[DEPTH][NCALO];  
  tmpaxi_t muparts[DEPTH][NMU];  
  for(int idepth = 0; idepth < DEPTH; idepth++) { 
    int   pEta  = idepth % 4;
    int   pPhi  = idepth/3;
    for(int i0 = 0; i0 < NTRACK; i0++) { 
      float pTEta = 0. + pEta*0.7+dphi*int(i0/2);
      float pTPhi = 0. + pPhi*0.7+deta*(1+int(i0/2));
      tmpaxi_t pTmp;
      pTmp.range(15, 0) = pt*PT_SCALE;
      pTmp.range(24,16) = pTEta*ETAPHI_SCALE;
      pTmp.range(34,25) = pTPhi*ETAPHI_SCALE;
      pTmp.range(38,35) = 0;
      pTmp.range(48,39) = 0;
      chparts[idepth][i0] = pTmp;
    }
    for(int i0 = 0; i0 < NEMCALO; i0++) { 
      float pTEta = 0. + pEta*0.7+dphi*(1+int(i0/2));
      float pTPhi = 0. + pPhi*0.7+deta*(3+int(i0/2));
      tmpaxi_t pTmp;
      pTmp.range(15, 0) = pt*PT_SCALE;
      pTmp.range(24,16) = pTEta*ETAPHI_SCALE;
      pTmp.range(34,25) = pTPhi*ETAPHI_SCALE;
      pTmp.range(38,35) = 2;
      pTmp.range(55,39) = pt*PT_SCALE;
      emparts[idepth][i0] = pTmp;
    }
    for(int i0 = 0; i0 < NSELCALO; i0++) { 
      float pTEta = 0. + pEta*0.7+dphi*(2+int(i0/2));
      float pTPhi = 0. + pPhi*0.7+deta*(6+int(i0/2));
      tmpaxi_t pTmp;
      pTmp.range(15, 0) = pt*PT_SCALE;
      pTmp.range(24,16) = pTEta*ETAPHI_SCALE;
      pTmp.range(34,25) = pTPhi*ETAPHI_SCALE;
      pTmp.range(38,35) = 2;
      pTmp.range(55,39) = pt*PT_SCALE;
      neparts[idepth][i0] = pTmp;
    }
    for(int i0 = 0; i0 < NMU; i0++) { 
      float pTEta = 0. + pEta*0.7+dphi*(5+int(i0/2));
      float pTPhi = 0. + pPhi*0.7+deta*(1+int(i0/2));
      tmpaxi_t pTmp;
      pTmp.range(15, 0) = pt*PT_SCALE;
      pTmp.range(24,16) = pTEta*ETAPHI_SCALE;
      pTmp.range(34,25) = pTPhi*ETAPHI_SCALE;
      pTmp.range(38,35) = 0;
      pTmp.range(48,39) = 0;
      muparts[idepth][i0] = pTmp;
    }
  } 
  //void algo_inputs_layer2_v3(
  MP7DataWord input[DEPTH][MP7_NCHANN];
  MP7DataWord output[DEPTH][MP7_NCHANN];
  for(int idepth = 0; idepth < DEPTH; idepth++) { 
    for(int i0 = 0; i0 < NTRACK;  i0++) {
      input[idepth][i0*2+0] = chparts[idepth][i0].range(31,0);  
      input[idepth][i0*2+1] = chparts[idepth][i0].range(63,32);  
    }
    for(int i0 = 0; i0 < NEMCALO; i0++) { 
      input[idepth][i0*2+0+EMOFFS] = emparts[idepth][i0].range(31,0);  
      input[idepth][i0*2+1+EMOFFS] = emparts[idepth][i0].range(63,32);  
    }
    for(int i0 = 0; i0 < NCALO;   i0++) { 
      input[idepth][i0*2+0+HAOFFS] = neparts[idepth][i0].range(31,0);  
      input[idepth][i0*2+1+HAOFFS] = neparts[idepth][i0].range(63,32);
    }  
    for(int i0 = 0; i0 < NMU;     i0++) { 
      input[idepth][i0*2+0+MUOFFS]  = muparts[idepth][i0].range(31,0);  
      input[idepth][i0*2+0+MUOFFS]  = muparts[idepth][i0].range(63,32);  
    }
    algo_inputs_layer2_v3(input[idepth],output[idepth]);
  }
  for(int idepth = 0; idepth < NTAU; idepth++) {
    for(int ipart = 0; ipart < DATA_SIZE; ipart++) { 
      PFChargedObj pTmp;
      float pPt  = output[idepth][ipart].range(15,0);  pPt/=PT_SCALE;
      float pEta = output[idepth][ipart].range(31,16); pEta/=ETAPHI_SCALE;
      float pPhi = output[idepth][ipart].range(15,0);  pPhi/=ETAPHI_SCALE;
      int   pId  = output[idepth][ipart].range(31,16); 
      std::cout << "===> tau part " << idepth << " -- part " << ipart << " vector " << pPt << "-- " << pEta << " -- " << pPhi << " - Id - " << pId << std::endl;
    }
  }
}

