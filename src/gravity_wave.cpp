// Universe OS — Gravity Wave Kernel  (C++/WASM port)
// A rotating binary (quadrupole source) radiates gravitational waves across a
// 2D spacetime grid solved with an explicit finite-difference wave equation
// (5-point Laplacian + leapfrog), with absorbing boundaries. All state, the
// PDE step, AND the rendering (olive.c into an RGBA framebuffer) are in C++;
// JS only blits the buffer and forwards UI events.
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

static const int COLS=160, ROWS=100, CELL=5;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 800 x 500
static std::vector<float> prevX, currX, nextX;
static std::vector<uint32_t> px;
static long pulseTimer=0;

// tunable params
static double g_speed=0.18, g_mass=10.0, g_omega=0.25, g_damp=0.996;

static inline int IDX(int i,int j){ return i*ROWS+j; }

// pack an olive color (0xAABBGGRR in memory = RGBA bytes) from rgb + alpha[0..1]
static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.0f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}

extern "C" {

KEEP int  sim_w(){ return FW; }
KEEP int  sim_h(){ return FH; }

KEEP void sim_reset(){
    prevX.assign((size_t)COLS*ROWS,0.f);
    currX.assign((size_t)COLS*ROWS,0.f);
    nextX.assign((size_t)COLS*ROWS,0.f);
    pulseTimer=0;
}

KEEP int sim_init(int hintW,int hintH){
    (void)hintW;(void)hintH;
    px.assign((size_t)FW*FH,0);
    sim_reset();
    return 1;
}

KEEP void sim_set(int id,double v){
    if(id==0) g_speed=v;
    else if(id==1) g_damp=v;
    else if(id==2) g_mass=v;
}
KEEP void sim_action(int){}

// click injects a positive pulse at the grid cell under the cursor (nx,ny in [0,1])
KEEP void sim_click(double nx,double ny){
    int i=(int)(nx*COLS), j=(int)(ny*ROWS);
    if(i>0&&i<COLS-1&&j>0&&j<ROWS-1) currX[IDX(i,j)] += g_mass*2.0;
}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        pulseTimer++;
        // 1. rotating binary injection (two anti-phase point masses)
        int cx=COLS/2, cy=ROWS/2; double r=4.0;
        int x1=(int)(cx+r*std::cos(pulseTimer*g_speed)), y1=(int)(cy+r*std::sin(pulseTimer*g_speed));
        int x2=(int)(cx-r*std::cos(pulseTimer*g_speed)), y2=(int)(cy-r*std::sin(pulseTimer*g_speed));
        if(x1>0&&x1<COLS&&y1>0&&y1<ROWS) currX[IDX(x1,y1)]= (float)g_mass;
        if(x2>0&&x2<COLS&&y2>0&&y2<ROWS) currX[IDX(x2,y2)]=-(float)g_mass;

        const double OMEGA=g_omega;
        // 2. leapfrog wave update with damping + absorbing boundary
        for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
            double lap = currX[IDX(i+1,j)]+currX[IDX(i-1,j)]+currX[IDX(i,j+1)]+currX[IDX(i,j-1)]-4.0*currX[IDX(i,j)];
            double nv = 2.0*currX[IDX(i,j)] - prevX[IDX(i,j)] + OMEGA*lap;
            nv *= g_damp;
            int de = std::min(std::min(i,COLS-1-i), std::min(j,ROWS-1-j));
            if(de<8) nv *= (de/8.0);
            nextX[IDX(i,j)] = (float)nv;
        }
        // 3. buffer rotation
        std::swap(prevX,currX);
        std::swap(currX,nextX);
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(5,5,10,1.f));                       // backdrop #05050a
    for(int x=0;x<FW;x+=50) olivec_rect(oc, x,0, 1,FH, rgba(0,17,5,1.f)); // faint grid
    for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
        float val=currX[IDX(i,j)], amp=std::fabs(val);
        if(amp>0.01f){ float a=amp*6.f; if(a>1.f)a=1.f;
            uint32_t c = val>0 ? rgba(0,255,200,a) : rgba(255,0,150,a);
            olivec_rect(oc, i*CELL, j*CELL, CELL, CELL, c); }
    }
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
// native self-test: run a while, dump a PNG to eyeball the spiral
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
int main(int argc,char**argv){
    sim_init(0,0);
    int steps = argc>1? atoi(argv[1]) : 400;
    sim_step(steps);
    uint8_t* p=sim_render();
    printf("gravity_wave native: %dx%d after %d steps\n", FW,FH,steps);
    stbi_write_png("gw_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
