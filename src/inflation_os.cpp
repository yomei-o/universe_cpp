// Universe OS — Cosmic Inflation Kernel  (C++/WASM port)
// A 3-phase state machine on a 2D wave grid. PRE_INFLATION: a tiny central
// Planck seed jitters with high-frequency quantum noise. INFLATING: the seed is
// inverse-mapped (stretched) across the whole grid faster than the L3 propagation
// speed, so the waveform is "frozen". POST_CMB: a normal damped wave equation
// (5-point Laplacian + leapfrog) evolves the frozen fluctuations into the CMB /
// cosmic-web pattern. All state, update and rendering (olive.c into an RGBA
// framebuffer) live in C++; JS only blits + forwards UI events.
#define OLIVEC_IMPLEMENTATION
#include "olive.c"
#include <vector>
#include <cstdint>
#include <cmath>
#include <random>
#include <cstring>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define KEEP EMSCRIPTEN_KEEPALIVE
#else
#define KEEP
#endif

static const int COLS=120, ROWS=75, CELL=6;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 720 x 450

static std::vector<float> prevX, currX, nextX;
static std::vector<uint32_t> px;

enum Phase { PRE_INFLATION=0, INFLATING=1, POST_CMB=2 };
static int bootPhase=PRE_INFLATION;
static double inflationScale=1.0;
static long elapsedFrames=0;

static std::mt19937 g_rng(12345);
static inline double frand(){ return (double)g_rng()/(double)0xFFFFFFFFu; } // [0,1]

static inline int IDX(int i,int j){ return i*ROWS+j; }

static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.0f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}

static const int MIDX=COLS/2, MIDY=ROWS/2;

extern "C" {

KEEP int  sim_w(){ return FW; }
KEEP int  sim_h(){ return FH; }

KEEP void sim_reset(){
    prevX.assign((size_t)COLS*ROWS,0.f);
    currX.assign((size_t)COLS*ROWS,0.f);
    nextX.assign((size_t)COLS*ROWS,0.f);
    bootPhase=PRE_INFLATION;
    inflationScale=1.0;
    elapsedFrames=0;
    g_rng.seed(12345);
}

KEEP int sim_init(int hintW,int hintH){
    (void)hintW;(void)hintH;
    px.assign((size_t)FW*FH,0);
    sim_reset();
    return 1;
}

KEEP void sim_set(int,double){}

// button: trigger the inflation boot sequence
KEEP void sim_action(int id){
    if(id==0 && bootPhase==PRE_INFLATION){
        bootPhase=INFLATING;
        elapsedFrames=0;
    }
}
KEEP void sim_click(double,double){}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        elapsedFrames++;
        const double OMEGA_SQ=0.24;

        if(bootPhase==PRE_INFLATION){
            for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
                double dx=i-MIDX, dy=j-MIDY; double dist=std::sqrt(dx*dx+dy*dy);
                if(dist<5.0){ float v=(float)((frand()-0.5)*35.0); currX[IDX(i,j)]=v; prevX[IDX(i,j)]=v; }
                else { currX[IDX(i,j)]=0; prevX[IDX(i,j)]=0; }
            }
        }
        else if(bootPhase==INFLATING){
            inflationScale += 0.8;
            std::vector<float> tmp((size_t)COLS*ROWS,0.f);
            for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
                double origI=MIDX+(i-MIDX)/inflationScale;
                double origJ=MIDY+(j-MIDY)/inflationScale;
                int srcI=(int)std::floor(origI), srcJ=(int)std::floor(origJ);
                if(srcI>=0&&srcI<COLS&&srcJ>=0&&srcJ<ROWS){
                    double v=(frand()-0.5)*4.0;
                    if(inflationScale<15){
                        double di=srcI-MIDX, dj=srcJ-MIDY;
                        if(std::sqrt(di*di+dj*dj)<5.0) v+=(frand()-0.5)*20.0;
                    }
                    tmp[IDX(i,j)]=(float)v;
                }
            }
            for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
                currX[IDX(i,j)]=tmp[IDX(i,j)];
                prevX[IDX(i,j)]=tmp[IDX(i,j)]*0.95f;
            }
            if(elapsedFrames>45) bootPhase=POST_CMB;
        }
        else { // POST_CMB — standard damped wave
            for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
                double lap=currX[IDX(i+1,j)]+currX[IDX(i-1,j)]+currX[IDX(i,j+1)]+currX[IDX(i,j-1)]-4.0*currX[IDX(i,j)];
                double nv=2.0*currX[IDX(i,j)]-prevX[IDX(i,j)]+OMEGA_SQ*lap;
                nv*=0.992;
                nextX[IDX(i,j)]=(float)nv;
            }
            for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
                prevX[IDX(i,j)]=currX[IDX(i,j)];
                currX[IDX(i,j)]=nextX[IDX(i,j)];
            }
        }
    }
}

static void draw_ring(Olivec_Canvas oc,int cx,int cy,int r,uint32_t c){
    const int N=48; int px0=0,py0=0;
    for(int k=0;k<=N;++k){ double a=2.0*M_PI*k/N;
        int x=cx+(int)(r*std::cos(a)), y=cy+(int)(r*std::sin(a));
        if(k>0) olivec_line(oc,px0,py0,x,y,c);
        px0=x; py0=y;
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(2,2,5,1.f));                    // #020205
    for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
        float amp=currX[IDX(i,j)];
        double energy=std::fabs(amp)*12.0;
        if(energy>1.5){
            float alpha=(float)std::min(1.0, energy/255.0);
            uint32_t c;
            if(bootPhase==PRE_INFLATION) c=rgba(0,255,200,alpha);
            else if(bootPhase==INFLATING) c=rgba(212,0,255,alpha*0.7f);
            else c = amp>0 ? rgba(255,100,0,alpha*0.85f) : rgba(130,0,255,alpha*0.6f);
            olivec_rect(oc, i*CELL, j*CELL, CELL-1, CELL-1, c);
        }
    }
    if(bootPhase==PRE_INFLATION){
        draw_ring(oc, MIDX*CELL, MIDY*CELL, 6*CELL, rgba(0,255,200,0.3f));
        olivec_text(oc, "PLANCK SEED", MIDX*CELL-40, MIDY*CELL-52, olivec_default_font, 2, rgba(0,255,200,0.6f));
    }
    const char* ph = bootPhase==PRE_INFLATION ? "QUANTUM JITTER"
                    : bootPhase==INFLATING ? "EXPONENTIAL ALLOCATION" : "FROZEN CMB MATRIX";
    char buf[128];
    std::snprintf(buf,sizeof(buf),"OS KERNEL PHASE : %s", ph);
    olivec_text(oc, buf, 20, FH-40, olivec_default_font, 2, rgba(0,255,204,1.f));
    const char* log = bootPhase==PRE_INFLATION ? "PRE-BOOT: quantum noise flickers in the Planck seed"
                    : bootPhase==INFLATING ? "INFLATION ACTIVE: exceeding L3 propagation speed C"
                    : "CMB fluctuation frozen. Reheating complete.";
    olivec_text(oc, log, 20, FH-22, olivec_default_font, 1, rgba(0,255,204,1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
int main(int argc,char**argv){
    sim_init(0,0);
    sim_action(0);                        // trigger inflation
    int steps = argc>1? atoi(argv[1]) : 80;
    sim_step(steps);
    uint8_t* p=sim_render();
    int nonbg=0; uint32_t bg=rgba(2,2,5,1.f);
    for(int k=0;k<FW*FH;++k) if(px[k]!=bg) nonbg++;
    printf("inflation_os native: %dx%d after %d steps, phase=%d, non-bg px=%d\n", FW,FH,steps,bootPhase,nonbg);
    stbi_write_png("inflation_os_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
