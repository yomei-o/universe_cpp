// Universe OS — Wormhole Routing Kernel  (C++/WASM port)
// A 2D finite-difference wave grid split into two disconnected "universes"
// (A on top, B on bottom) by forcing the middle row to zero every frame.
// Opening the wormhole directly couples two distant cells (pointer routing):
// the entrance amplitude is copied to the exit and drained from the entrance,
// letting a wave packet teleport with zero propagation latency. All state,
// the PDE step AND the rendering (olive.c into RGBA) live in C++; JS only
// blits and forwards UI.
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

static const int COLS=110, ROWS=50, CELL=7;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 770 x 350
static const int MID_Y=ROWS/2;                 // 25 (boundary row)
static const int ENTRANCE_X=75, ENTRANCE_Y=12; // gate A (upper universe)
static const int EXIT_X=25,     EXIT_Y=37;     // gate B (lower universe)

static std::vector<float> prevX, currX, nextX;
static std::vector<uint32_t> px;
static bool wormholeOpen=false;

static double g_omega=0.23, g_damp=0.985;

static inline int IDX(int i,int j){ return i*ROWS+j; }

static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.0f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}

static void inject_pulse(int cx,int cy){
    for(int dy=-2;dy<=2;++dy){
        int y=cy+dy;
        if(cx>0&&cx<COLS&&y>0&&y<ROWS){
            float v=25.0f*std::exp(-(float)(dy*dy)/2.0f);
            currX[IDX(cx,y)]=v; prevX[IDX(cx,y)]=v;
        }
    }
}

extern "C" {

KEEP int  sim_w(){ return FW; }
KEEP int  sim_h(){ return FH; }

KEEP void sim_reset(){
    prevX.assign((size_t)COLS*ROWS,0.f);
    currX.assign((size_t)COLS*ROWS,0.f);
    nextX.assign((size_t)COLS*ROWS,0.f);
    wormholeOpen=false;
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

// action 0: toggle wormhole ; action 1: fire pulse at A-side entrance
KEEP void sim_action(int id){
    if(id==0) wormholeOpen=!wormholeOpen;
    else if(id==1) inject_pulse(3, ENTRANCE_Y);
}

KEEP void sim_click(double nx,double ny){
    int i=(int)(nx*COLS), j=(int)(ny*ROWS);
    if(j==MID_Y) j++;
    inject_pulse(i,j);
}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        const double OMEGA=g_omega;
        // 1. vacuum-spring wave update (each half propagates independently)
        for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
            if(j==MID_Y){ currX[IDX(i,j)]=0.f; continue; }  // sever the two spaces
            double lap = currX[IDX(i+1,j)]+currX[IDX(i-1,j)]+currX[IDX(i,j+1)]+currX[IDX(i,j-1)]-4.0*currX[IDX(i,j)];
            double nv = 2.0*currX[IDX(i,j)] - prevX[IDX(i,j)] + OMEGA*lap;
            nv *= g_damp;
            nextX[IDX(i,j)] = (float)nv;
        }
        // 2. wormhole pointer routing (direct cell coupling)
        if(wormholeOpen){
            float tmp = currX[IDX(ENTRANCE_X,ENTRANCE_Y)];
            currX[IDX(EXIT_X,EXIT_Y)] += tmp*0.95f;
            currX[IDX(ENTRANCE_X,ENTRANCE_Y)] *= 0.1f;
        }
        // 3. buffer rotation (prev<-curr, curr<-next)
        std::swap(prevX,currX);
        std::swap(currX,nextX);
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(5,5,10,1.f));
    // faint vertical grid (every 5 cells ~ 35px)
    for(int i=0;i<COLS;i+=5) olivec_rect(oc, i*CELL,0, 1,FH, rgba(0,17,5,1.f));
    // opaque kernel boundary between universe A and B
    olivec_rect(oc, 0, MID_Y*CELL-2, FW, 4, rgba(0x22,0x22,0x33,1.f));
    // wave (L3 packet) render, 250x gain
    for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
        float amp=std::fabs(currX[IDX(i,j)]);
        float energy=amp*250.f;
        if(energy>5.f){ float a=energy/255.f; if(a>1.f)a=1.f;
            olivec_rect(oc, i*CELL, j*CELL, CELL-1, CELL-1, rgba(0,255,200,a)); }
    }
    // space labels
    olivec_text(oc, "SPACE LAYER A", 12, 8, olivec_default_font, 2, rgba(0x55,0x55,0x77,1.f));
    olivec_text(oc, "SPACE LAYER B", 12, MID_Y*CELL+8, olivec_default_font, 2, rgba(0x55,0x55,0x77,1.f));

    int ex=ENTRANCE_X*CELL+CELL/2, ey=ENTRANCE_Y*CELL+CELL/2;
    int xx=EXIT_X*CELL+CELL/2,      xy=EXIT_Y*CELL+CELL/2;
    uint32_t gateCol = wormholeOpen ? rgba(0,170,255,1.f) : rgba(0x55,0x55,0x66,1.f);
    if(wormholeOpen){
        // dashed index-link between the gates
        float dx=(float)(xx-ex), dy=(float)(xy-ey);
        float len=std::sqrt(dx*dx+dy*dy); int segs=(int)(len/8.f);
        for(int k=0;k<segs;k+=2){
            float t0=k/(float)segs, t1=(k+1)/(float)segs;
            olivec_line(oc, (int)(ex+dx*t0),(int)(ey+dy*t0), (int)(ex+dx*t1),(int)(ey+dy*t1), rgba(0,170,255,0.8f));
        }
    }
    // entrance / exit gate disks
    olivec_circle(oc, ex,ey, 9, wormholeOpen?rgba(0,170,255,0.35f):rgba(100,100,100,0.15f));
    olivec_circle(oc, xx,xy, 9, wormholeOpen?rgba(0,170,255,0.35f):rgba(100,100,100,0.15f));
    olivec_text(oc, "ENTRANCE", ex-30, ey-22, olivec_default_font, 1, gateCol);
    olivec_text(oc, "EXIT",     xx-14, xy-22, olivec_default_font, 1, gateCol);
    // status log
    const char* log = wormholeOpen ? "ROUTE: pointer bridge (75,12)<->(25,37)"
                                    : "STATUS: ISOLATED MATRICES";
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
    sim_action(0);        // open wormhole
    sim_action(1);        // fire a pulse into entrance
    int steps = argc>1? atoi(argv[1]) : 300;
    sim_step(steps);
    uint8_t* p=sim_render();
    printf("wormhole_os native: %dx%d after %d steps\n", FW,FH,steps);
    stbi_write_png("wormhole_os_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
