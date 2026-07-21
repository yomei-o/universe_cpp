// Universe OS — Chemical spiral (Barkley reaction-diffusion)  (C++/WASM)
// Excitable-medium PDE with an activator u and inhibitor v (Barkley model),
// integrated with a tiny dt and 12 sub-steps/frame for stability. A cross-block
// initial seed spawns a rotating spiral core; clicking injects fresh
// activator/inhibitor to nucleate new spirals. State + physics + rendering in C++.
#define OLIVEC_IMPLEMENTATION
#include "olive.c"
#include <vector>
#include <cstdint>
#include <cmath>
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define KEEP EMSCRIPTEN_KEEPALIVE
#else
#define KEEP
#endif

static const int CELL=4, COLS=200, ROWS=140;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 800 x 560
static std::vector<double> uMap, vMap, nU, nV;
static std::vector<uint32_t> px;
static long total_time=0;

// Barkley parameters
static double g_a=0.75, g_b=0.02, g_eps=0.02, g_Du=0.20, g_dt=0.002;
static const int SUBSTEPS=12;

static uint32_t rng=2463534242u;
static inline double rnd(){ rng^=rng<<13; rng^=rng>>17; rng^=rng<<5; return (rng&0xFFFFFF)/(double)0x1000000; }
static inline int IDX(int i,int j){ return i*ROWS+j; }
static inline uint32_t rgb(int r,int g,int b){ return 0xFF000000u|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r; }

extern "C" {
KEEP int sim_w(){ return FW; }
KEEP int sim_h(){ return FH; }

KEEP void sim_reset(){
    uMap.assign((size_t)COLS*ROWS,0.0); vMap.assign((size_t)COLS*ROWS,0.0);
    nU.assign((size_t)COLS*ROWS,0.0);   nV.assign((size_t)COLS*ROWS,0.0);
    total_time=0;
    for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){ uMap[IDX(i,j)]=rnd()*0.01; vMap[IDX(i,j)]=rnd()*0.01; }
    // cross-block seed -> guaranteed free end -> spiral
    int cx=COLS/2, cy=ROWS/2, bs=30;
    for(int i=cx-bs;i<=cx;++i) for(int j=cy-bs;j<=cy;++j) if(i>0&&j>0&&i<COLS&&j<ROWS) uMap[IDX(i,j)]=1.0;
    for(int i=cx-bs;i<=cx;++i) for(int j=cy+1;j<=cy+bs;++j) if(i>0&&j>0&&i<COLS&&j<ROWS) vMap[IDX(i,j)]=0.6;
}
KEEP int sim_init(int,int){ px.assign((size_t)FW*FH,0); sim_reset(); return 1; }

KEEP void sim_set(int id,double v){
    if(id==0) g_eps=1.0/v;      // reaction speed (raw 20..80)
    else if(id==1) g_a=v;       // excitation threshold
}
KEEP void sim_action(int id){ if(id==0) sim_reset(); }

// inject: right half of the brush = activator, left half = inhibitor (makes a free end)
KEEP void sim_click(double nx,double ny){
    int gi=(int)(nx*COLS), gj=(int)(ny*ROWS), R=10;
    for(int di=-R;di<=R;++di) for(int dj=-R;dj<=R;++dj){
        int ii=gi+di, jj=gj+dj;
        if(ii>0&&ii<COLS-1&&jj>0&&jj<ROWS-1 && di*di+dj*dj<=R*R){
            if(di>0) uMap[IDX(ii,jj)]=1.0; else vMap[IDX(ii,jj)]=0.6;
        }
    }
}

KEEP void sim_step(int frames){
    for(int f=0;f<frames;++f){
        for(int s=0;s<SUBSTEPS;++s){
            for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
                double u=uMap[IDX(i,j)], v=vMap[IDX(i,j)];
                double lapU=uMap[IDX(i+1,j)]+uMap[IDX(i-1,j)]+uMap[IDX(i,j+1)]+uMap[IDX(i,j-1)]-4.0*u;
                double fu=(1.0/g_eps)*u*(1.0-u)*(u-(v+g_b)/g_a);
                double gv=u-v;
                nU[IDX(i,j)]=u+g_dt*(fu+g_Du*lapU);
                nV[IDX(i,j)]=v+g_dt*gv;
            }
            for(int i=0;i<COLS;++i){ nU[IDX(i,0)]=0;nU[IDX(i,ROWS-1)]=0;nV[IDX(i,0)]=0;nV[IDX(i,ROWS-1)]=0; }
            for(int j=0;j<ROWS;++j){ nU[IDX(0,j)]=0;nU[IDX(COLS-1,j)]=0;nV[IDX(0,j)]=0;nV[IDX(COLS-1,j)]=0; }
            std::swap(uMap,nU); std::swap(vMap,nV);
        }
        total_time++;
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc=olivec_canvas(px.data(),FW,FH,FW);
    olivec_fill(oc, rgb(1,2,10));
    for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
        double u=uMap[IDX(i,j)], v=vMap[IDX(i,j)];
        if(u>0.01||v>0.01){
            double iu=std::min(1.0,u)*255.0, iv=std::min(1.0,v)*180.0;
            uint32_t c = iu>20 ? rgb((int)(iu*0.3),255,235)
                               : rgb((int)(35+iv*0.8),0,(int)(55+iv*1.0));
            olivec_rect(oc, i*CELL, j*CELL, CELL, CELL, c);
        }
    }
    olivec_text(oc, "BARKLEY REACTION-DIFFUSION SPIRAL", 16, FH-22, olivec_default_font, 2, rgb(0,255,204));
    return (uint8_t*)px.data();
}
}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
int main(int argc,char**argv){
    sim_init(0,0);
    int steps=argc>1?atoi(argv[1]):300; sim_step(steps);
    uint8_t*p=sim_render();
    int nb=0; for(int k=0;k<FW*FH;++k) if(px[k]!=rgb(1,2,10)) nb++;
    printf("spiral_os native: %dx%d after %d frames, non-bg=%d\n",FW,FH,steps,nb);
    stbi_write_png("spiral_os_preview.png",FW,FH,4,p,FW*4);
    return 0;
}
#endif
