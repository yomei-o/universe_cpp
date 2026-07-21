// Universe OS — Wave-Particle Duality Kernel  (C++/WASM port)
// A 2D finite-difference wave grid with a double-slit wall. Normally the field
// propagates as a continuous probability wave (float). When the detector is
// enabled (observe mode), the wave that reaches the screen column collapses:
// one row is sampled with probability weight = amplitude^2 and tallied into a
// histogram, so the interference pattern emerges dot-by-dot. All state, the
// PDE step AND the rendering (olive.c into RGBA) live in C++; JS only blits.
//
// Common WASM ABI shared by every Universe OS sim:
//   sim_init(hintW,hintH) -> 1 ok ; sim_w()/sim_h() give the framebuffer size
//   sim_reset() ; sim_step(n) ; sim_render() -> RGBA* ; sim_click(nx,ny)
//   sim_set(id,val) ; sim_action(id)
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

static const int COLS=110, ROWS=56, CELL=7;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 770 x 392
static const int SLIT_COL=35;
static const int SLIT_OPEN1=ROWS/2-5;          // 23  (openings 23,24)
static const int SLIT_OPEN2=ROWS/2+5;          // 33  (openings 33,34)
static const int SCREEN_COL=COLS-10;           // 100

static std::vector<float> prevX, currX, nextX;
static std::vector<int> screenBuffer;
static std::vector<uint32_t> px;
static long pulseTimer=0;
static bool isObserving=false;
static int lastHitRow=-1, lastHitTimer=0;

static double g_omega=0.25, g_damp=0.993;

// simple xorshift PRNG -> uniform [0,1)
static uint32_t rngState=0x2545F491u;
static inline double frand(){
    rngState^=rngState<<13; rngState^=rngState>>17; rngState^=rngState<<5;
    return (rngState & 0xFFFFFFu)/16777216.0;
}

static inline int IDX(int i,int j){ return i*ROWS+j; }

static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.0f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}

static inline bool isSlitOpen(int j){
    return j==SLIT_OPEN1||j==SLIT_OPEN1+1||j==SLIT_OPEN2||j==SLIT_OPEN2+1;
}

extern "C" {

KEEP int  sim_w(){ return FW; }
KEEP int  sim_h(){ return FH; }

KEEP void sim_reset(){
    prevX.assign((size_t)COLS*ROWS,0.f);
    currX.assign((size_t)COLS*ROWS,0.f);
    nextX.assign((size_t)COLS*ROWS,0.f);
    screenBuffer.assign((size_t)ROWS,0);
    pulseTimer=0; lastHitRow=-1; lastHitTimer=0;
    rngState=0x2545F491u;
}

KEEP int sim_init(int hintW,int hintH){
    (void)hintW;(void)hintH;
    px.assign((size_t)FW*FH,0);
    sim_reset();
    return 1;
}

KEEP void sim_set(int id,double v){
    if(id==0) g_omega=v;
    else if(id==1) g_damp=v;
}

// action 0: toggle observe (collapse) mode
KEEP void sim_action(int id){
    if(id==0){
        isObserving=!isObserving;
        if(isObserving){ screenBuffer.assign((size_t)ROWS,0); lastHitRow=-1; lastHitTimer=0; }
    }
}

KEEP void sim_click(double nx,double ny){
    int i=(int)(nx*COLS), j=(int)(ny*ROWS);
    if(i>1&&i<COLS-1&&j>1&&j<ROWS-1){ currX[IDX(i,j)]=40.f; prevX[IDX(i,j)]=40.f; }
}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        pulseTimer++;
        // plane-wave injection from the left, every 60 frames
        if(pulseTimer%60==1){
            int midY=ROWS/2;
            for(int dy=-3;dy<=3;++dy){
                float v=70.0f*std::exp(-(float)(dy*dy)/4.0f);
                currX[IDX(2,midY+dy)]=v; prevX[IDX(2,midY+dy)]=v;
            }
        }
        const double OMEGA=g_omega;
        for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
            if(i==SLIT_COL && !isSlitOpen(j)){ currX[IDX(i,j)]=0.f; continue; }
            double lap = currX[IDX(i+1,j)]+currX[IDX(i-1,j)]+currX[IDX(i,j+1)]+currX[IDX(i,j-1)]-4.0*currX[IDX(i,j)];
            double nv = 2.0*currX[IDX(i,j)] - prevX[IDX(i,j)] + OMEGA*lap;
            nv *= g_damp;
            nextX[IDX(i,j)] = (float)nv;
        }
        // wavefront reaches the screen at frame 46 of each cycle -> collapse
        if(pulseTimer%60==46 && isObserving){
            double total=0.0;
            for(int j=2;j<ROWS-2;++j){ double e=(double)currX[IDX(SCREEN_COL,j)]; total+=e*e; }
            if(total>0.01){
                double r=frand()*total, sum=0.0;
                for(int j=2;j<ROWS-2;++j){
                    double e=(double)currX[IDX(SCREEN_COL,j)]; sum+=e*e;
                    if(r<=sum){ screenBuffer[j]++; lastHitRow=j; lastHitTimer=20; break; }
                }
            }
        }
        if(lastHitTimer>0) lastHitTimer--;
        std::swap(prevX,currX);
        std::swap(currX,nextX);
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(5,5,10,1.f));
    for(int x=0;x<FW;x+=40) olivec_rect(oc, x,0, 1,FH, rgba(0,17,5,1.f));

    // wave field (culled right of the screen while observing)
    for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
        if(isObserving && i>=SCREEN_COL) continue;
        float energy=std::fabs(currX[IDX(i,j)])*250.f;
        if(energy>5.f){ float a=energy/255.f; if(a>1.f)a=1.f;
            olivec_rect(oc, i*CELL, j*CELL, CELL-1, CELL-1, rgba(0,255,200,a)); }
    }

    // double-slit wall
    olivec_rect(oc, SLIT_COL*CELL, 0, CELL, FH, rgba(0x33,0x44,0x55,1.f));
    olivec_rect(oc, SLIT_COL*CELL, SLIT_OPEN1*CELL, CELL, CELL*2, rgba(5,5,10,1.f));
    olivec_rect(oc, SLIT_COL*CELL, SLIT_OPEN2*CELL, CELL, CELL*2, rgba(5,5,10,1.f));

    // detector screen line
    uint32_t scrCol = isObserving ? rgba(255,170,0,1.f) : rgba(0,170,255,1.f);
    olivec_frame(oc, SCREEN_COL*CELL, 0, 10, FH-1, 1, scrCol);

    // screen readout
    for(int j=0;j<ROWS;++j){
        if(!isObserving){
            float energy=std::fabs(currX[IDX(SCREEN_COL,j)])*200.f;
            if(energy>2.f){ float a=energy/255.f; if(a>1.f)a=1.f;
                olivec_rect(oc, SCREEN_COL*CELL+14, j*CELL, 40, CELL-1, rgba(0,255,170,a)); }
        } else {
            int hits=screenBuffer[j];
            if(hits>0){
                float a=std::min(255,hits*12)/255.f;
                int w=std::min(150,hits*6);
                olivec_rect(oc, SCREEN_COL*CELL+14, j*CELL, w, CELL-1, rgba(255,180,0,a));
                olivec_rect(oc, SCREEN_COL*CELL+3, j*CELL+CELL/2-1, 4, 2, rgba(255,255,255,0.4f));
            }
        }
    }

    // impact flash
    if(isObserving && lastHitRow!=-1 && lastHitTimer>0){
        float fa=lastHitTimer/20.f;
        olivec_rect(oc, SCREEN_COL*CELL-4, lastHitRow*CELL-2, 18, CELL+4, rgba(255,255,255,fa));
        olivec_rect(oc, SCREEN_COL*CELL+14, lastHitRow*CELL, 40, CELL-1, rgba(255,215,0,fa));
    }

    olivec_text(oc, isObserving?"SCREEN: PARTICLE DOTS (QUANTIZED)":"SCREEN: WAVE INTERFERENCE (CONTINUOUS)",
                SCREEN_COL*CELL-190, 8, olivec_default_font, 1, scrCol);
    const char* log = isObserving ? "INTERRUPT: Detector active. Sampling float->int."
                                  : "STATUS: COHERENT WAVE (unobserved probability field)";
    olivec_text(oc, log, 16, FH-16, olivec_default_font, 1, rgba(0,255,200,1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
int main(int argc,char**argv){
    sim_init(0,0);
    sim_action(0);        // observe mode ON to exercise histogram
    int steps = argc>1? atoi(argv[1]) : 400;
    sim_step(steps);
    uint8_t* p=sim_render();
    printf("wave_particle_os native: %dx%d after %d steps\n", FW,FH,steps);
    stbi_write_png("wave_particle_os_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
