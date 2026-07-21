// Universe OS — Superconductivity & Meissner Kernel  (C++/WASM port)
// A 2D finite-difference wave grid modelling electron packets in a conductor.
// At high temperature, thermal jitter (noise) is injected and the field is
// damped (resistance). Cooling to T=0 phase-locks the central layer: waves
// there propagate loss-free (persistent current) while the surroundings still
// decay, and vertical magnetic field lines are expelled upward from the layer
// (Meissner effect). All state, PDE step AND rendering live in C++.
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

static const int COLS=100, ROWS=56, CELL=8;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 800 x 448
static const int SC_Y_START=ROWS/2-2;          // 26
static const int SC_Y_END  =ROWS/2+8;          // 36

static std::vector<float> prevX, currX, nextX;
static std::vector<uint32_t> px;
static double temperature=1.0;
static bool isSuperconducting=false;
static double g_omega=0.25;

static uint32_t rngState=0x9E3779B9u;
static inline double frand(){
    rngState^=rngState<<13; rngState^=rngState>>17; rngState^=rngState<<5;
    return (rngState & 0xFFFFFFu)/16777216.0;
}

static inline int IDX(int i,int j){ return i*ROWS+j; }

static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.0f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}

static void inject_packet(){
    int midY=ROWS/2+4;  // 32
    for(int i=10;i<20;++i){
        float v=20.f*std::sin((i-10)*0.3f);
        currX[IDX(i,midY)]=v; prevX[IDX(i,midY)]=v;
    }
}

extern "C" {

KEEP int  sim_w(){ return FW; }
KEEP int  sim_h(){ return FH; }

KEEP void sim_reset(){
    prevX.assign((size_t)COLS*ROWS,0.f);
    currX.assign((size_t)COLS*ROWS,0.f);
    nextX.assign((size_t)COLS*ROWS,0.f);
    temperature=1.0; isSuperconducting=false;
    rngState=0x9E3779B9u;
}

KEEP int sim_init(int hintW,int hintH){
    (void)hintW;(void)hintH;
    px.assign((size_t)FW*FH,0);
    sim_reset();
    return 1;
}

// id 0: temperature slider (0 -> superconducting)
KEEP void sim_set(int id,double v){
    if(id==0){ temperature=v; isSuperconducting=(temperature<=0.0); }
    else if(id==1) g_omega=v;
}

// action 0: instant cool to T=0 ; action 1: inject electron packet
KEEP void sim_action(int id){
    if(id==0){ temperature=0.0; isSuperconducting=true; }
    else if(id==1) inject_packet();
}

KEEP void sim_click(double nx,double ny){
    int i=(int)(nx*COLS), j=(int)(ny*ROWS);
    if(i>1&&i<COLS-1&&j>1&&j<ROWS-1){ currX[IDX(i,j)]=20.f; prevX[IDX(i,j)]=20.f; }
}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        const double OMEGA=g_omega;
        for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
            double lap = currX[IDX(i+1,j)]+currX[IDX(i-1,j)]+currX[IDX(i,j+1)]+currX[IDX(i,j-1)]-4.0*currX[IDX(i,j)];
            double nv = 2.0*currX[IDX(i,j)] - prevX[IDX(i,j)] + OMEGA*lap;
            if(temperature>0){
                nv += (frand()-0.5)*temperature*2.5;
                nv *= 0.95;                 // resistive damping
            } else {
                if(j>=SC_Y_START && j<=SC_Y_END) nv *= 1.0;   // loss-free layer
                else                              nv *= 0.92;  // outside conductor
            }
            nextX[IDX(i,j)] = (float)nv;
        }
        std::swap(prevX,currX);
        std::swap(currX,nextX);
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(5,5,10,1.f));

    // superconducting layer plate
    uint32_t plateFill = isSuperconducting ? rgba(0,255,200,0.08f) : rgba(100,100,120,0.10f);
    olivec_rect(oc, 0, SC_Y_START*CELL, FW, (SC_Y_END-SC_Y_START)*CELL, plateFill);
    uint32_t plateEdge = isSuperconducting ? rgba(0,255,200,1.f) : rgba(0x55,0x55,0x66,1.f);
    olivec_frame(oc, 0, SC_Y_START*CELL, FW-1, (SC_Y_END-SC_Y_START)*CELL, 1, plateEdge);
    olivec_text(oc, "SUPERCONDUCTING LAYER", 12, SC_Y_START*CELL+8, olivec_default_font, 1, plateEdge);

    // electron wave (persistent current) brightness
    for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
        float energy=std::fabs(currX[IDX(i,j)])*15.f;
        if(energy>5.f){ float a=energy/255.f; if(a>1.f)a=1.f;
            uint32_t c = isSuperconducting ? rgba(0,255,200,a) : rgba(150,160,180,a);
            olivec_rect(oc, i*CELL, j*CELL, CELL-1, CELL-1, c); }
    }

    // magnetic field lines (COLS/4 vertical lines), expelled upward when SC
    float pushTop = 0.f;
    if(isSuperconducting) pushTop = (SC_Y_END-SC_Y_START+1)*1.8f;
    int jmax=(int)(ROWS*0.7);
    for(int i=0;i<COLS;i+=4){
        int px0=i*CELL;
        int prevY=0;
        for(int j=0;j<jmax;++j){
            float y;
            if(isSuperconducting && j>=SC_Y_START-8 && j<=SC_Y_END)
                y = j*CELL - pushTop*(1.f-(j-SC_Y_START+8)/10.f);
            else
                y = (float)(j*CELL);
            if(j>0) olivec_line(oc, px0, prevY, px0, (int)y, rgba(0,150,255,0.7f));
            prevY=(int)y;
        }
    }

    // external magnet bar
    olivec_rect(oc, FW/2-60, 0, 120, 15, rgba(255,51,51,1.f));
    olivec_text(oc, "EXTERNAL MAGNET (N)", FW/2-50, 3, olivec_default_font, 1, rgba(255,255,255,1.f));

    // status log
    const char* log = isSuperconducting ? "LOG: Phase-lock established. Resistance 0.0 Ohm"
                                         : "STATUS: HIGH TEMPERATURE (thermal jitter)";
    olivec_text(oc, log, 16, FH-16, olivec_default_font, 1, rgba(0,255,200,1.f));
    char tbuf[64];
    snprintf(tbuf,sizeof(tbuf),"SYSTEM TEMP (T): %.2f %s", temperature, temperature==0?"(CRITICAL ZERO)":"JITTER_LOAD");
    olivec_text(oc, tbuf, 16, FH-32, olivec_default_font, 1, rgba(0,255,200,1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
int main(int argc,char**argv){
    sim_init(0,0);
    sim_action(0);        // cool to superconducting
    sim_action(1);        // inject packet
    int steps = argc>1? atoi(argv[1]) : 200;
    sim_step(steps);
    uint8_t* p=sim_render();
    printf("superconductivity_os native: %dx%d after %d steps\n", FW,FH,steps);
    stbi_write_png("superconductivity_os_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
