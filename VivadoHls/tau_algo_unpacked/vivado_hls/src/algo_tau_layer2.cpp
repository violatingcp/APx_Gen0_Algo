#include "algo_tau_layer2.h"
#include "firmware/tau_nn.h"
#include "sorting_network.h"

//16+10+10+3+10 for pt,eta,phi,particleid,z0
//typedef ap_axis <64*NPART,1,1,1> axi_t;
//typedef hls::stream<axi_t> hls::stream<axi_t>;
//stream depth is TM6 and 240 MHz gives 24 clocks

template<typename T, int NIn, int NOut>
void ptsort_hwopt_ind(T in[NIn], T out[NOut]) { 
    #pragma HLS PIPELINE
    T tmp[NOut];
    #pragma HLS ARRAY_PARTITION variable=tmp complete
    for (int iout = 0; iout < NOut; ++iout) {
        #pragma HLS unroll
        tmp[iout].hwPt = 0;
    }

    for (int it = 0; it < NIn; ++it) {
        for (int iout = NOut-1; iout >= 0; --iout) {
            if (tmp[iout].hwPt <= in[it].hwPt) {
                if (iout == 0 || tmp[iout-1].hwPt > in[it].hwPt) {
                    tmp[iout] = in[it];
                } else {
                    tmp[iout] = tmp[iout-1];
                }
            }
        }

    }
    for (int iout = 0; iout < NOut; ++iout) {
        out[iout] = tmp[iout];
    }

}
template<unsigned int N, unsigned int OFFS> 
inline void mp7_pack(PFChargedObj obj[N], axi_t &data) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
      #pragma HLS UNROLL
      int pOffset = i*64;
      data.data.range(15+pOffset,0 +pOffset) = obj[i].hwPt;
      data.data.range(24+pOffset,16+pOffset) = obj[i].hwEta;
      data.data.range(34+pOffset,25+pOffset) = obj[i].hwPhi;
      data.data.range(38+pOffset,35+pOffset) = obj[i].hwId;
      data.data.range(48+pOffset,39+pOffset) = obj[i].hwZ0;
    }
}
template<unsigned int N, unsigned int OFFS> 
inline void mp7_pack(hls::stream<PFChargedObj> obj[], axi_t &data) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
      #pragma HLS UNROLL
      int pOffset = i*64;
      PFChargedObj tmpobj = obj[i].read();
      data.data.range(15+pOffset,0 +pOffset) = tmpobj.hwPt;
      data.data.range(24+pOffset,16+pOffset) = tmpobj.hwEta;
      data.data.range(34+pOffset,25+pOffset) = tmpobj.hwPhi;
      data.data.range(38+pOffset,35+pOffset) = tmpobj.hwId;
      data.data.range(48+pOffset,39+pOffset) = tmpobj.hwZ0;
    }
}
inline void mp7_pack(hls::stream<PFChargedObj> &obj, axi_t &data) {
    #pragma HLS inline
  PFChargedObj tmpobj = obj.read();
  data.data.range(15,0 ) = tmpobj.hwPt;
  data.data.range(24,16) = tmpobj.hwEta;
  data.data.range(34,25) = tmpobj.hwPhi;
  data.data.range(38,35) = tmpobj.hwId;
  data.data.range(48,39) = tmpobj.hwZ0;
}
template<unsigned int N>
void make_inputs(input_t nn_data[N*8], hls::stream<PFChargedObj> pf[N]) {
  //#pragma HLS inline
  #pragma HLS PIPELINE
  for (int i = 0; i < N; i++) {
        #pragma HLS PIPELINE II=1
    PFChargedObj tmpobj = pf[i].read();
    nn_data[i*8+0] = input_t(tmpobj.hwPt);
    nn_data[i*8+1] = input_t(tmpobj.hwEta);
    nn_data[i*8+2] = input_t(tmpobj.hwPhi);
    nn_data[i*8+3] = input_t(tmpobj.hwId == 2 ? 1 : 0);
    nn_data[i*8+4] = input_t(tmpobj.hwId == 3 ? 1 : 0);
    nn_data[i*8+5] = input_t(tmpobj.hwId == 4 ? 1 : 0);
    nn_data[i*8+6] = input_t(tmpobj.hwId == 1 ? 1 : 0);
    nn_data[i*8+7] = input_t(tmpobj.hwId == 0 ? 1 : 0);
  }
} 
void algo_tau_layer2(hls::stream<PFChargedObj > allparts_in [DATA_SIZE],hls::stream<PFChargedObj > &allparts_out) { //hls::stream<axi_t> &link_out) {
  #pragma HLS PIPELINE II=12
  #pragma HLS INTERFACE axis port=allparts_in
  #pragma HLS INTERFACE axis port=allparts_out
  #pragma HLS INTERFACE s_axilite port=return bundle=ctrl

  //hls::stream<PFChargedObj > allparts_out;
  #pragma HLS stream variable=allparts_in  depth=6
  #pragma HLS stream variable=allparts_out depth=6  
  //#pragma HLS stream variable=link_out     depth=5

  sorting_network_64_nn(allparts_in,allparts_out);
  //axi_t data_out;
  //mp7_pack(allparts_out,data_out);
  //link_out.write(data_out);
  /*
  for(int itau = 0; itau < NTAU; itau++) {
    PFChargedObj dummyc; 
    result_t taunn[N_OUTPUTS];
    input_t  nn_data[NTAUPARTS*8];
    make_inputs<NTAUPARTS>(nn_data,allparts_out);            
    tau_nn(nn_data,taunn);
    dummyc.hwPt = taunn[0]*100; dummyc.hwEta = nn_data[1]; dummyc.hwPhi = nn_data[2]; dummyc.hwId = 0; dummyc.hwZ0 = 0;
    //taus[itau] = dummyc;
    axi_t data_out;
    mp7_pack<1,0>(&dummyc,data_out);
    link_out.write(data_out);
  }
  //PFChargedObj tausout[NTAU]; 
  //#pragma HLS ARRAY_PARTITION variable=tausout complete
  //ptsort_hwopt_ind<PFChargedObj,NTAU,NTAU>(taus, tausout);
  */
}

