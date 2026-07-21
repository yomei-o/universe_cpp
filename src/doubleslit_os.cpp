// Universe OS — Quantum Measurement / Double-Slit Kernel  (C++/WASM port)
// A wave packet is pulsed from the left edge and passes through a two-slit wall
// solved with an explicit finite-difference wave equation. With the observation
// sensor OFF the wave crosses both slits and builds an interference fringe on the
// right screen; turning it ON forces a "memory-sampling interrupt" at the slit
// cells (phase snapped to a particle amplitude), collapsing the fringe pattern.
// All state, the PDE step, AND rendering (olive.c) live in C++; JS only blits.
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
#include <cstdio>
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define KEEP EMSCRIPTEN_KEEPALIVE
#else
#define KEEP
#endif

static const int COLS=120, ROWS=60, CELL=7;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 840 x 420
static std::vector<float> prevX, currX, nextX;
static std::vector<double> screenPat;
static std::vector<uint32_t> px;

static const int WALL_COL=40, SLIT_1_ROW=22, SLIT_2_ROW=38, SLIT_WIDTH=3;
static const int SCREEN_COL=COLS-5;
static long spawnTimer=0;
static bool isObserving=false;

static double g_omega=0.15;

static inline int IDX(int i,int j){ return i*ROWS+j; }
static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.0f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}
static inline bool inSlit(int j){
    return (j>=SLIT_1_ROW&&j<SLIT_1_ROW+SLIT_WIDTH)||(j>=SLIT_2_ROW&&j<SLIT_2_ROW+SLIT_WIDTH);
}

extern "C" {

KEEP int  sim_w(){ return FW; }
KEEP int  sim_h(){ return FH; }

KEEP void sim_reset(){
    prevX.assign((size_t)COLS*ROWS,0.f);
    currX.assign((size_t)COLS*ROWS,0.f);
    nextX.assign((size_t)COLS*ROWS,0.f);
    screenPat.assign(ROWS,0.0);
    spawnTimer=0;
}

KEEP int sim_init(int hintW,int hintH){
    (void)hintW;(void)hintH;
    px.assign((size_t)FW*FH,0);
    sim_reset();
    return 1;
}

KEEP void sim_set(int id,double v){ if(id==0) g_omega=v; }
// toggle observation sensor (resets the accumulated fringe)
KEEP void sim_action(int id){ if(id==0){ isObserving=!isObserving; screenPat.assign(ROWS,0.0); } }

KEEP void sim_click(double nx,double ny){
    int i=(int)(nx*COLS), j=(int)(ny*ROWS);
    if(i>0&&i<COLS-1&&j>0&&j<ROWS-1){ currX[IDX(i,j)]+=15.f; prevX[IDX(i,j)]+=15.f; }
}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        // 1. electron gun pulse from left edge
        spawnTimer++;
        if(spawnTimer%25==0){ int cr=ROWS/2; currX[IDX(2,cr)]=15.f; prevX[IDX(2,cr)]=15.f; }

        // 2. L3 spring difference update
        const double OMEGA_SQ=g_omega;
        for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
            if(i==WALL_COL && !inSlit(j)){ nextX[IDX(i,j)]=0.f; continue; }
            double lap=currX[IDX(i+1,j)]+currX[IDX(i-1,j)]+currX[IDX(i,j+1)]+currX[IDX(i,j-1)]-4.0*currX[IDX(i,j)];
            double nv=2.0*currX[IDX(i,j)]-prevX[IDX(i,j)]+OMEGA_SQ*lap;
            nextX[IDX(i,j)]=(float)nv;
            // observation: memory-sampling interrupt at slit cells
            if(isObserving && i==WALL_COL && inSlit(j)){
                if(std::fabs(nextX[IDX(i,j)])>0.1f) nextX[IDX(i,j)]= nextX[IDX(i,j)]>0?3.0f:-3.0f;
            }
        }

        // 3. accumulate screen energy (interference fringe)
        for(int j=0;j<ROWS;++j){
            double e=(double)currX[IDX(SCREEN_COL,j)]*currX[IDX(SCREEN_COL,j)];
            if(e>0.05) screenPat[j]+=e*0.02;
        }

        // 4. buffer shift
        std::swap(prevX,currX);
        std::swap(currX,nextX);
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(5,5,10,1.f));                               // #05050a

    // field brightness
    for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
        float val=std::fabs(currX[IDX(i,j)])*35.f;
        if(val>5.f){ float a=val/255.f; if(a>1.f)a=1.f;
            olivec_rect(oc, i*CELL, j*CELL, CELL-1, CELL-1, rgba(0,255,200,a)); }
    }

    // mid wall (#555566) with two slit gaps
    int wx=WALL_COL*CELL;
    olivec_rect(oc, wx,0, CELL, SLIT_1_ROW*CELL, rgba(85,85,102,1.f));
    olivec_rect(oc, wx,(SLIT_1_ROW+SLIT_WIDTH)*CELL, CELL, (SLIT_2_ROW-(SLIT_1_ROW+SLIT_WIDTH))*CELL, rgba(85,85,102,1.f));
    olivec_rect(oc, wx,(SLIT_2_ROW+SLIT_WIDTH)*CELL, CELL, FH-(SLIT_2_ROW+SLIT_WIDTH)*CELL, rgba(85,85,102,1.f));

    // observation sensor eyes (red)
    if(isObserving){
        olivec_circle(oc, wx+CELL/2, (int)((SLIT_1_ROW+SLIT_WIDTH/2.0)*CELL), 8, rgba(255,51,102,1.f));
        olivec_circle(oc, wx+CELL/2, (int)((SLIT_2_ROW+SLIT_WIDTH/2.0)*CELL), 8, rgba(255,51,102,1.f));
    }

    // electron gun
    olivec_rect(oc, 0,(ROWS/2-2)*CELL, 15, 4*CELL, rgba(153,153,170,1.f));

    // screen receiving line (#005544)
    olivec_rect(oc, SCREEN_COL*CELL,0, 2,FH, rgba(0,85,68,1.f));

    // accumulated fringe bars (#00ffcc)
    for(int j=0;j<ROWS;++j){
        int gw=(int)(screenPat[j]*5.0);
        if(gw>0) olivec_rect(oc, SCREEN_COL*CELL, j*CELL+1, gw, CELL-2, rgba(0,255,200,0.4f));
    }

    olivec_text(oc, isObserving?"OBSERVING: fringe collapsed (particle)":"WAVE: interference fringe forming",
                12, FH-16, olivec_default_font, 2, rgba(0,255,204,1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
int main(int argc,char**argv){
    sim_init(0,0);
    int steps = argc>1? atoi(argv[1]) : 900;
    sim_step(steps);
    uint8_t* p=sim_render();
    long nb=0; uint32_t bg=rgba(5,5,10,1.f);
    for(size_t k=0;k<(size_t)FW*FH;++k) if(px[k]!=bg) nb++;
    printf("doubleslit_os native: %dx%d after %d steps, non-bg px=%ld\n", FW,FH,steps,nb);
    stbi_write_png("doubleslit_os_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
