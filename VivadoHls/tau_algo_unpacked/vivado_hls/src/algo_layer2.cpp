#include "algo_layer2.h"
#include "firmware/tau_nn.h"
#include "sorting_network.h"

//16+10+10+3+10 for pt,eta,phi,particleid,z0
//typedef ap_axis <64*NPART,1,1,1> axi_t;
//typedef hls::stream<axi_t> hls::stream<axi_t>;
//stream depth is TM6 and 240 MHz gives 24 clocks

template<int NB>
ap_uint<NB> dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<NB> max) {
    #pragma HLS inline
    #pragma HLS PIPELINE
    //hardcode for etaphi size
    int tmpe = eta1-eta2;
    ap_uint<NB> deta = (tmpe > 0 ? tmpe : -tmpe);
    int tmpp = phi1-phi2;
    ap_uint<NB> dphi = (tmpp > 0 ? tmpp : -tmpp);
    int dr2 = max;
    if ((deta >> (NB/2))==0 && (dphi >> (NB/2))==0) {
        ap_uint<NB> deta2 = deta*deta;
        ap_uint<NB> dphi2 = dphi*dphi;
        dr2 = deta2 + dphi2;
    }
    return (dr2 < int(max) ? ap_uint<NB>(dr2) : max);
}
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
template<unsigned int N> 
inline void mp7_unpack(axi_t data, PFChargedObj emcalo[N]) {
    #pragma HLS inline
    #pragma HLS PIPELINE
    for (unsigned int i = 0; i < N; ++i) {
      int pOffset = (i)*64;
      emcalo[i].hwPt       = data.data(15+pOffset, 0+pOffset);
      emcalo[i].hwEta      = data.data(24+pOffset,16+pOffset);
      emcalo[i].hwPhi      = data.data(34+pOffset,25+pOffset);
      emcalo[i].hwId       = data.data(38+pOffset,35+pOffset);
      emcalo[i].hwZ0       = data.data(48+pOffset,39+pOffset);
    }
}

template<unsigned int N> 
inline void mp7_unpack_seed(axi_t data, PFChargedObj &seed) {
    #pragma HLS inline
    #pragma HLS PIPELINE
    seed.hwPt          = data.data(15, 0);
    seed.hwEta         = data.data(24,16);
    seed.hwPhi         = data.data(34,25);
    seed.hwId          = data.data(38,35);
    seed.hwZ0          = data.data(48,39);
}

template<unsigned int N,unsigned int NR> 
inline void mp7_unpack(axi_t data, PFChargedObj emcalo[N*NR],int iBase) {
    #pragma HLS inline
    #pragma HLS PIPELINE
    for (unsigned int i = 0; i < N; ++i) {
      int pOffset = (i)*64;
      emcalo[iBase+i].hwPt       = data.data(15+pOffset, 0+pOffset);
      emcalo[iBase+i].hwEta      = data.data(24+pOffset,16+pOffset);
      emcalo[iBase+i].hwPhi      = data.data(34+pOffset,25+pOffset);
      emcalo[iBase+i].hwId       = data.data(38+pOffset,35+pOffset);
      emcalo[iBase+i].hwZ0       = data.data(48+pOffset,39+pOffset);
    }
}

template<unsigned int N> 
inline void mp7_unpack(axi_t data, PFNeutralObj emcalo[N]) { //,int iBase) {
    #pragma HLS inline
    #pragma HLS PIPELINE
    for (unsigned int i = 0; i < N; ++i) {
      int pOffset = (i)*64;
      emcalo[i].hwPt       = data.data(15+pOffset, 0+pOffset);
      emcalo[i].hwEta      = data.data(24+pOffset,16+pOffset);
      emcalo[i].hwPhi      = data.data(34+pOffset,25+pOffset);
      emcalo[i].hwId       = data.data(38+pOffset,35+pOffset);
      emcalo[i].hwPtPuppi  = data.data(55+pOffset,39+pOffset);
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
template<unsigned int N> void sumparts(pt_t &ipt,PFChargedObj iCol[4*N]) {
  #pragma HLS inline
  pt_t pt1 = 0; 
  pt_t pt2 = 0; 
  pt_t pt3 = 0; 
  pt_t pt4 = 0; 
  for(int i = 0; i < N; i++) { 
    pt1 = pt1 + iCol[i].hwPt;
  }
  for(int i = 0; i < N; i++) { 
    pt2 = pt2 + iCol[i+N].hwPt;
  }
  for(int i = 0; i < N; i++) { 
    pt3 = pt3 + iCol[i+N*2].hwPt;
  }
  for(int i = 0; i < N; i++) { 
    pt4 = pt4 + iCol[i+N*3].hwPt;
  }
  ipt = pt1 + pt2 + pt3 + pt4;
} 
void sumpt(pt_t &taupt, PFChargedObj pfch[NTRACK*4], PFChargedObj pfpho[NPHOTON*4], PFChargedObj pfne[NSELCALO*4], PFChargedObj pfmu[NMU*4]) {
  #pragma HLS inline
  #pragma HLS PIPELINE
  pt_t tauptch   = 0;
  pt_t tauptpho  = 0;
  pt_t tauptne   = 0;
  pt_t tauptmu   = 0;
  sumparts<NTRACK>  (tauptch ,pfch);
  sumparts<NPHOTON> (tauptpho,pfpho);
  sumparts<NSELCALO>(tauptne ,pfne);
  sumparts<NMU>     (tauptmu ,pfmu);
  taupt = tauptch + tauptpho + tauptne + tauptmu;
}
template<int N,int NMAX> 
void deltaR(int iOffSet, etaphi_t seedeta,etaphi_t seedphi,PFChargedObj pfch[N],PFChargedObj pfout[NMAX*4]) { 
  #pragma HLS inline
  #pragma HLS PIPELINE 
  const ap_int<16> eDR2MAX = DR2MAX;
  PFChargedObj dummyc; dummyc.hwPt = 0; dummyc.hwEta = 0; dummyc.hwPhi = 0; dummyc.hwId = 0; dummyc.hwZ0 = 0;
  for (int i = 0; i < NMAX; i++) {
   #pragma HLS UNROLL
    int drcheck = dr2_int_cap<12>(seedeta,seedphi,pfch[i].hwEta,pfch[i].hwPhi,eDR2MAX);
    if(drcheck < DRCONE) { 
     pfout[iOffSet+i] = pfch[i];
    } else { 
     pfout[iOffSet+i] = dummyc;
    }
  }
}
template<int N> 
void deltaR(int iOffSet, etaphi_t seedeta,etaphi_t seedphi,PFNeutralObj pfne[N],PFNeutralObj pfne_sel[N*4]) { 
  #pragma HLS inline 
  #pragma HLS PIPELINE 
  const ap_int<16> eDR2MAX = DR2MAX;
  PFNeutralObj dummyn; dummyn.hwPt = 0; dummyn.hwEta = 0; dummyn.hwPhi = 0; dummyn.hwId = 0; dummyn.hwPtPuppi = 0;
  for (int i = 0; i < N; i++) {
   #pragma HLS UNROLL
   int drcheck = dr2_int_cap<12>(seedeta,seedphi,pfne[i].hwEta,pfne[i].hwPhi,eDR2MAX);
   if(drcheck < DRCONE) { 
    pfne_sel[iOffSet+i] = pfne[i];
   } else { 
    pfne_sel[iOffSet+i] = dummyn;
   }
  }
}

template<int N_TABLE>
void init_map_table(ap_int<8>  arr[N_TABLE], int iId){
  //4 blocks across*2 8 blocks                                                                                                                                                                                                                                                 
  for(int iregion = 0; iregion < N_TABLE; iregion++) {
    bool pUp    = true;
    bool pLeft  = true;
    bool pEdgeLeft  = false;
    bool pEdgeRight = false;
    bool pEdgeTop   = false;
    bool pEdgeBot   = false;
    int  pRow    = iregion/8;
    int  pX      = (iregion % 8)/2;
    int  pY      = iregion/16;
    int  pRegion = pY*4+pX;
    if(iregion % 2 == 1) pLeft        = false;
    if(pRow    % 2 == 1) pUp          = false;
    if(iregion % 8 == 0) pEdgeLeft    = true;
    if(iregion % 8 == 7) pEdgeRight   = true;
    if(iregion < 8     ) pEdgeTop     = true;
    if(iregion > N_TABLE-9) pEdgeBot  = true;
    if(pUp && pLeft) {
      if(iId == 0) arr[iregion] = pRegion;
      if(iId == 1) arr[iregion] = pRegion-4;
      if(iId == 2) arr[iregion] = pRegion-5;
      if(iId == 3) arr[iregion] = pRegion-1;
    }
    if(pUp && !pLeft) {
      if(iId == 0) arr[iregion] = pRegion;
      if(iId == 1) arr[iregion] = pRegion-4;
      if(iId == 2) arr[iregion] = pRegion-3;
      if(iId == 3) arr[iregion] = pRegion+1;
    }
    if(!pUp && pLeft) {
      if(iId == 0) arr[iregion] = pRegion;
      if(iId == 1) arr[iregion] = pRegion+4;
      if(iId == 2) arr[iregion] = pRegion+3;
      if(iId == 3) arr[iregion] = pRegion-1;
    }
    if(!pUp && !pLeft) {
      if(iId == 0) arr[iregion] = pRegion;
      if(iId == 1) arr[iregion] = pRegion+4;
      if(iId == 2) arr[iregion] = pRegion+5;
      if(iId == 3) arr[iregion] = pRegion+1;
    }
    if(pEdgeTop && iId == 1) arr[iregion] += N_TABLE/4;
    if(pEdgeTop && iId == 2) arr[iregion] += N_TABLE/4;
    if(pEdgeBot && iId == 1) arr[iregion] -= N_TABLE/4;
    if(pEdgeBot && iId == 2) arr[iregion] -= N_TABLE/4;
    if(pEdgeLeft || pEdgeRight) {
      if(iId == 2) arr[iregion] = 0;
      if(iId == 3) arr[iregion] = 0;
    }
  }
}

template<int N_TABLE>
void init_map_table0(ap_int<8>  arr[N_TABLE]){
  //4 blocks across*2 8 blocks                                                                                                                                                                                                                                                 
  for(int iregion = 0; iregion < N_TABLE; iregion++) {
    bool pUp    = true;
    bool pLeft  = true;
    bool pEdgeLeft  = false;
    bool pEdgeRight = false;
    bool pEdgeTop   = false;
    bool pEdgeBot   = false;
    int  pRow    = iregion/8;
    int  pX      = (iregion % 8)/2;
    int  pY      = iregion/16;
    int  pRegion = pY*4+pX;
    if(iregion % 2 == 1) pLeft        = false;
    if(pRow    % 2 == 1) pUp          = false;
    if(iregion % 8 == 0) pEdgeLeft    = true;
    if(iregion % 8 == 7) pEdgeRight   = true;
    if(iregion < 8     ) pEdgeTop     = true;
    if(iregion > N_TABLE-9) pEdgeBot  = true;
    if(pUp && pLeft) {
      arr[iregion] = pRegion;
    }
    if(pUp && !pLeft) {
      arr[iregion] = pRegion;
    }
    if(!pUp && pLeft) {
      arr[iregion] = pRegion;      
    }
    if(!pUp && !pLeft) {
      arr[iregion] = pRegion;
    }
  }
}
template<int N_TABLE>
void init_map_table1(ap_int<8>  arr[N_TABLE]){
  //4 blocks across*2 8 blocks                                                                                                                                                                                                                                                 
  for(int iregion = 0; iregion < N_TABLE; iregion++) {
    bool pUp    = true;
    bool pLeft  = true;
    bool pEdgeLeft  = false;
    bool pEdgeRight = false;
    bool pEdgeTop   = false;
    bool pEdgeBot   = false;
    int  pRow    = iregion/8;
    int  pX      = (iregion % 8)/2;
    int  pY      = iregion/16;
    int  pRegion = pY*4+pX;
    if(iregion % 2 == 1) pLeft        = false;
    if(pRow    % 2 == 1) pUp          = false;
    if(iregion % 8 == 0) pEdgeLeft    = true;
    if(iregion % 8 == 7) pEdgeRight   = true;
    if(iregion < 8     ) pEdgeTop     = true;
    if(iregion > N_TABLE-9) pEdgeBot  = true;
    if(pUp && pLeft) {
      arr[iregion] = pRegion-4;
    }
    if(pUp && !pLeft) {
      arr[iregion] = pRegion-4;
    }
    if(!pUp && pLeft) {
      arr[iregion] = pRegion+4;
    }
    if(!pUp && !pLeft) {
      arr[iregion] = pRegion+4;
    }
    if(pEdgeTop) arr[iregion] += N_TABLE/4;
    if(pEdgeBot) arr[iregion] -= N_TABLE/4;
  }
}

template<int N_TABLE>
void init_map_table2(ap_int<8>  arr[N_TABLE]){
  //4 blocks across*2 8 blocks                                                                                                                                                                                                                                                 
  for(int iregion = 0; iregion < N_TABLE; iregion++) {
    bool pUp    = true;
    bool pLeft  = true;
    bool pEdgeLeft  = false;
    bool pEdgeRight = false;
    bool pEdgeTop   = false;
    bool pEdgeBot   = false;
    int  pRow    = iregion/8;
    int  pX      = (iregion % 8)/2;
    int  pY      = iregion/16;
    int  pRegion = pY*4+pX;
    if(iregion % 2 == 1) pLeft        = false;
    if(pRow    % 2 == 1) pUp          = false;
    if(iregion % 8 == 0) pEdgeLeft    = true;
    if(iregion % 8 == 7) pEdgeRight   = true;
    if(iregion < 8     ) pEdgeTop     = true;
    if(iregion > N_TABLE-9) pEdgeBot  = true;
    if(pUp && pLeft) {
      arr[iregion] = pRegion-5;
    }
    if(pUp && !pLeft) {
      arr[iregion] = pRegion-3;
    }
    if(!pUp && pLeft) {
      arr[iregion] = pRegion+3;
    }
    if(!pUp && !pLeft) {
      arr[iregion] = pRegion+5;
    }
    if(pEdgeTop) arr[iregion] += N_TABLE/4;
    if(pEdgeBot) arr[iregion] -= N_TABLE/4;
    if(pEdgeLeft || pEdgeRight) {
      arr[iregion] = 0;
    }
  }
}

template<int N_TABLE>
void init_map_table3(ap_int<8>  arr[N_TABLE]){
  //4 blocks across*2 8 blocks                                                                                                                                                                                                                                                 
  for(int iregion = 0; iregion < N_TABLE; iregion++) {
    bool pUp    = true;
    bool pLeft  = true;
    bool pEdgeLeft  = false;
    bool pEdgeRight = false;
    bool pEdgeTop   = false;
    bool pEdgeBot   = false;
    int  pRow    = iregion/8;
    int  pX      = (iregion % 8)/2;
    int  pY      = iregion/16;
    int  pRegion = pY*4+pX;
    if(iregion % 2 == 1) pLeft        = false;
    if(pRow    % 2 == 1) pUp          = false;
    if(iregion % 8 == 0) pEdgeLeft    = true;
    if(iregion % 8 == 7) pEdgeRight   = true;
    if(iregion < 8     ) pEdgeTop     = true;
    if(iregion > N_TABLE-9) pEdgeBot  = true;
    if(pUp && pLeft) {
      arr[iregion] = pRegion-1;
    }
    if(pUp && !pLeft) {
      arr[iregion] = pRegion+1;
    }
    if(!pUp && pLeft) {
      arr[iregion] = pRegion-1;
    }
    if(!pUp && !pLeft) {
      arr[iregion] = pRegion+1;
    }
    if(pEdgeLeft || pEdgeRight) {
      arr[iregion] = 0;
    }
  }
}

void arraymap(etaphi_t iEta,etaphi_t iPhi,ap_int<8> arrs[4]) {
 #pragma HLS inline
 #pragma HLS PIPELINE
#ifdef __HLS_SYN__
  bool initialized = false;
  ap_int<8> arr0[NREGIONS*4];
  ap_int<8> arr1[NREGIONS*4];
  ap_int<8> arr2[NREGIONS*4];
  ap_int<8> arr3[NREGIONS*4];
#else
  static bool initialized = false;
  static ap_int<8> arr0[NREGIONS*4];
  static ap_int<8> arr1[NREGIONS*4];
  static ap_int<8> arr2[NREGIONS*4];
  static ap_int<8> arr3[NREGIONS*4];
#endif
  if (!initialized) {
    init_map_table0<NREGIONS*4>(arr0);
    init_map_table1<NREGIONS*4>(arr1);
    init_map_table2<NREGIONS*4>(arr2);
    init_map_table3<NREGIONS*4>(arr3);
    initialized = true;
  }
  ap_int<8> pEta = ap_int<16>(iEta*0.0125) % 8;
  ap_int<8> pPhi = ap_int<16>(iPhi*0.0125) % (NREGIONS/4);
  arrs[0] = arr0[pPhi*8+pEta];
  arrs[1] = arr1[pPhi*8+pEta];
  arrs[2] = arr2[pPhi*8+pEta];
  arrs[3] = arr3[pPhi*8+pEta];
}
void convert_out(PFChargedObj datas[NTAUPARTS], hls::stream<PFChargedObj> axis_out[DATA_SIZE]) {
	// input data from stream
	for (int i = 0; i < NTAUPARTS; ++i) {
	#pragma HLS UNROLL
	  datas[i] = axis_out[i].read();
	}
	/*  Below is needed to remove some errors when NTAUPARTS != DATA_SIZE
	for (int i = NTAUPARTS; i < DATA_SIZE; ++i) {
	  #pragma HLS unroll
	  axis_out[i].read();
	 }*/
}
void convert_in(PFChargedObj datas[DATA_SIZE], hls::stream<PFChargedObj > axis_in[DATA_SIZE]) {
	// output data to stream
	for (int i = 0; i < DATA_SIZE; ++i) {
	#pragma HLS unroll
		axis_in[i].write(datas[i]);
	}
}
template<unsigned int N0,unsigned int N1, unsigned int N2, unsigned int N3,unsigned int NOUT>
void tausort(PFChargedObj allparts[NOUT], PFChargedObj pfch[N0], PFChargedObj pfpho[N1], PFChargedObj pfne[N2], PFChargedObj pfmu[N3]) {
  #pragma HLS PIPELINE
  
  #pragma HLS ARRAY_PARTITION variable=allparts complete
  for(int i = 0; i < N0; i++) { 
    #pragma HLS UNROLL
    allparts[i] = pfch[i];
  }
  for(int i = 0; i < N1; i++) { 
    #pragma HLS UNROLL
    allparts[i+N0] = pfpho[i];
  }
  for(int i = 0; i < N2; i++) { 
    #pragma HLS UNROLL
    allparts[i+N0+N1] = pfne[i];
  }
  for(int i = 0; i < N3; i++) { 
    #pragma HLS UNROLL
    allparts[i+N0+N1+N2] = pfmu[i];
  }
}
template<unsigned int N>
void make_inputs(input_t nn_data[N*8], hls::stream<PFChargedObj> pf[N]) {
  #pragma HLS inline
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
void algo_layer2(hls::stream<axi_t> &ch_link_in,hls::stream<axi_t> &ne_link_in, hls::stream<axi_t> &em_link_in, hls::stream<axi_t> &mu_link_in,
		 hls::stream<axi_t> &link_out) {
  #pragma HLS PIPELINE
  #pragma HLS INTERFACE axis port=ch_link_in
  #pragma HLS INTERFACE axis port=ne_link_in
  #pragma HLS INTERFACE axis port=em_link_in
  #pragma HLS INTERFACE axis port=mu_link_in
  #pragma HLS INTERFACE axis port=link_out
  #pragma HLS INTERFACE s_axilite port=return bundle=ctrl
  
  #pragma HLS stream variable=ch_link_in depth=36
  #pragma HLS stream variable=ne_link_in depth=36
  #pragma HLS stream variable=em_link_in depth=36
  #pragma HLS stream variable=mu_link_in depth=36
  #pragma HLS stream variable=link_out   depth=36

  PFChargedObj pfch[NREGIONS][NTRACK];
  PFChargedObj pfpho[NREGIONS][NPHOTON];
  PFChargedObj pfne[NREGIONS][NSELCALO];
  PFChargedObj pfmu[NREGIONS][NMU];
  #pragma HLS ARRAY_RESHAPE variable=pfch  block factor=25   dim=2
  #pragma HLS ARRAY_RESHAPE variable=pfpho block factor=20   dim=2
  #pragma HLS ARRAY_RESHAPE variable=pfne  block factor=20   dim=2
  #pragma HLS ARRAY_RESHAPE variable=pfmu  block factor=2    dim=2
  #pragma HLS RESOURCE      variable=pfch  core=RAM_2P_BRAM
  #pragma HLS RESOURCE      variable=pfpho core=RAM_2P_BRAM
  #pragma HLS RESOURCE      variable=pfne  core=RAM_2P_BRAM
  #pragma HLS RESOURCE      variable=pfmu  core=RAM_2P_BRAM
  
  //First take the seed from each region
  PFChargedObj chseed[NREGIONS];
  #pragma HLS ARRAY_PARTITION variable=chseed complete
  LoopA:
  for(int idepth =0; idepth < NREGIONS; idepth++) { 
    #pragma HLS PIPELINE II=1
    axi_t ch_in = ch_link_in.read();
    mp7_unpack_seed<NTRACK> (ch_in, chseed[idepth]);
    mp7_unpack<NTRACK> (ch_link_in.read(), pfch [idepth]);
    mp7_unpack<NEMCALO>(em_link_in.read(), pfpho[idepth]);
    mp7_unpack<NCALO>  (ne_link_in.read(), pfne [idepth]);
    mp7_unpack<NMU>    (mu_link_in.read(), pfmu [idepth]);
  }
  //Next iterate through seeds and compute DR
  PFChargedObj seeds[NTAU];
  #pragma HLS ARRAY_PARTITION variable=seeds complete
  ptsort_hwopt_ind<PFChargedObj,NREGIONS,NTAU>(chseed, seeds);
  const ap_int<16> eDR2MAX = DR2MAX;
  ap_uint<6> seedsort2[NTAU];
  #pragma HLS ARRAY_PARTITION variable=seedsort2 complete
  LoopB:
  for(int iseed0=0; iseed0 < NTAU; iseed0++) {
   #pragma HLS UNROLL
    seedsort2[iseed0]=iseed0;
  }
  LoopC:
  for(int iseed0=0; iseed0 < NTAU; iseed0++) {
    #pragma HLS UNROLL
    for(int iseed1 = 1; iseed1 < NTAU; iseed1++) {
      if(iseed1 < iseed0+1) continue;                                                                                                                                                                                                                                         
      int drcheck = dr2_int_cap<16>(seeds[iseed0].hwEta,seeds[iseed0].hwPhi,seeds[iseed1].hwEta,seeds[iseed0].hwPhi,eDR2MAX);
      if(drcheck < DRCONE && iseed1 > iseed0) {
        seedsort2[iseed1] = 0;
      }
    }
  }
  //Now get all the candidates in each region => Note that a max of 4 regoins can fill array                                                                                                                                                                                  
  PFChargedObj pfchreg [NTAU][NTAUPARTS*4];
  PFChargedObj pfemreg [NTAU][NTAUPARTS*4];
  PFChargedObj pfnereg [NTAU][NTAUPARTS*4];
  PFChargedObj pfmureg [NTAU][NMU*4];
  #pragma HLS ARRAY_RESHAPE variable=pfchreg  dim=0 complete
  #pragma HLS ARRAY_RESHAPE variable=pfemreg  dim=0 complete
  #pragma HLS ARRAY_RESHAPE variable=pfnereg  dim=0 complete
  #pragma HLS ARRAY_RESHAPE variable=pfmureg  dim=0 complete
  PFChargedObj dummyc; dummyc.hwPt = 0; dummyc.hwEta = 0; dummyc.hwPhi = 0; dummyc.hwId = 0; dummyc.hwZ0 = 0;
  LoopD:
  for(int itau = 0; itau < NTAU; itau++) {
    #pragma HLS UNROLL
    etaphi_t seedeta = seeds[itau].hwEta;
    etaphi_t seedphi = seeds[itau].hwPhi;
    ap_int<8> arr[4];
    #pragma HLS ARRAY_PARTITION variable=arr  complete dim=0
    arraymap(seedeta,seedphi,arr);
    for(int iregion = 0; iregion < 4; iregion++) {
      #pragma HLS UNROLL
      int pRegion = arr[iregion];
      int NTRACKOFFSET   = iregion*NTAUPARTS;//NTRACK;
      int NPHOTONOFFSET  = iregion*NTAUPARTS;//NPHOTON;
      int NSELCALOOFFSET = iregion*NTAUPARTS;//NSELCALO;
      int NMUOFFSET      = iregion*NMU;
      deltaR<NTRACK,NTAUPARTS>   (NTRACKOFFSET,  seedeta,seedphi,pfch [pRegion],pfchreg[itau]);
      deltaR<NPHOTON,NTAUPARTS>  (NPHOTONOFFSET, seedeta,seedphi,pfpho[pRegion],pfemreg[itau]);
      deltaR<NSELCALO,NTAUPARTS> (NSELCALOOFFSET,seedeta,seedphi,pfne [pRegion],pfnereg[itau]);
      deltaR<NMU,NMU>            (NMUOFFSET,     seedeta,seedphi,pfmu [pRegion],pfmureg[itau]);
    }
  }
  hls::stream<PFChargedObj > allparts_in [DATA_SIZE];
  hls::stream<PFChargedObj > allparts_out[DATA_SIZE];
  #pragma HLS stream variable=allparts_in  depth=5
  #pragma HLS stream variable=allparts_out depth=5
  for(int itau = 0; itau < NTAU; itau++) {
    PFChargedObj allparts[DATA_SIZE];
    #pragma HLS ARRAY_PARTITION variable=allparts complete
    tausort<4*NTAUPARTS,4*NTAUPARTS,4*NTAUPARTS,4*NMU,DATA_SIZE>(allparts,pfchreg[itau],pfnereg[itau],pfemreg[itau],pfmureg[itau]);  
    convert_in(allparts,allparts_in);
  }
  sorting_network_128(allparts_in,allparts_out);
  PFChargedObj taus[NTAU];
  #pragma HLS ARRAY_PARTITION variable=taus complete
  for(int itau = 0; itau < NTAU; itau++) {
    input_t nn_data[NTAUPARTS*8];
    //#pragma HLS ARRAY_PARTITION variable=nn_data complete
    make_inputs<NTAUPARTS>(nn_data,allparts_out);            
    result_t taunn[N_OUTPUTS];
    tau_nn(nn_data,taunn);
    PFChargedObj dummyc; 
    dummyc.hwPt = taunn[0]*100;; dummyc.hwEta = nn_data[1]; dummyc.hwPhi = nn_data[2]; dummyc.hwId = 0; dummyc.hwZ0 = 0;
    taus[itau] = dummyc;
  }
  PFChargedObj tausout[NTAU]; 
  #pragma HLS ARRAY_PARTITION variable=tausout complete
  ptsort_hwopt_ind<PFChargedObj,NTAU,NTAU>(taus, tausout);
  axi_t data_out;
  mp7_pack<NTAU,0>(tausout,data_out);
  link_out.write(data_out);
}
