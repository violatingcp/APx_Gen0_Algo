#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "ap_int.h"
#include "algo_layer2_v2.h"

using namespace std;

typedef ap_int<64>            tmpaxi_t;
typedef ap_axis <64*NPART,1,1,1> axi_t;
typedef hls::stream<axi_t> stream_t;

int main(int argc, char ** argv) {
  hls::stream<axi_t>     ch_link_in;
  hls::stream<axi_t>     ne_link_in;
  hls::stream<axi_t>     em_link_in;
  hls::stream<axi_t>     mu_link_in;
  hls::stream<PFChargedObj> link_out;

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
  for(int idepth = 0; idepth < DEPTH; idepth++) { 
     axi_t tmpch; 
     axi_t tmpne; 
     axi_t tmpem; 
     axi_t tmpmu; 
     for(int i0 = 0; i0 < NTRACK;  i0++) tmpch.data.range(63*(i0+1),64*(i0)) = chparts[idepth][i0];  
     for(int i0 = 0; i0 < NEMCALO; i0++) tmpem.data.range(63*(i0+1),64*(i0)) = emparts[idepth][i0];  
     for(int i0 = 0; i0 < NCALO;   i0++) tmpne.data.range(63*(i0+1),64*(i0)) = neparts[idepth][i0];  
     for(int i0 = 0; i0 < NMU;     i0++) tmpmu.data.range(63*(i0+1),64*(i0)) = muparts[idepth][i0];  
     //tmp.
     ch_link_in.write(tmpch);
     ne_link_in.write(tmpne);
     em_link_in.write(tmpem);
     mu_link_in.write(tmpmu);
  }
  algo_layer2_v2(ch_link_in,ne_link_in, em_link_in, mu_link_in,  link_out);
  for(int idepth = 0; idepth < NTAU; idepth++) {
    PFChargedObj pTmp; link_out.read(pTmp);
    float pPt  = pTmp.hwPt;pPt/=PT_SCALE;
    float pEta = pTmp.hwEta;pEta/=ETAPHI_SCALE;
    float pPhi = pTmp.hwPhi;pPhi/=ETAPHI_SCALE;
    std::cout << "===> depth " << idepth << " -- part " << " vector " << pPt << "-- " << pEta << " -- " << pPhi << std::endl;
  }
}

