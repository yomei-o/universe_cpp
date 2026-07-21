// Universe OS — Anti-Matter Annihilation Kernel  (C++/WASM port)
// Matter (+A) and anti-matter (-A) register indices march toward the center of a
// 2D vacuum spring grid, each injecting a signed static amplitude. When the two
// indices collide, the cell is zeroed ("garbage-collected") and its mass load is
// refactored into a violent high-frequency L3 wave (a 3x3 explosive kick) that
// propagates via an explicit finite-difference wave equation. All state, the PDE
// step, AND rendering (olive.c into an RGBA framebuffer) live in C++; JS only
// blits the buffer and forwards UI events.
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

static const int COLS=110, ROWS=56, CELL=7;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 770 x 392
static std::vector<float> prevX, currX, nextX;
static std::vector<uint32_t> px;

// particle register indices
static int matterX=18, matterY=ROWS/2;
static int antiX=92,   antiY=ROWS/2;
static bool isRunning=true, isAnnihilated=false;
static int flashTimer=0;
static long clockT=0;

// tunable params
static double g_omega=0.26, g_damp=0.992;

static inline int IDX(int i,int j){ return i*ROWS+j; }

// pack an olive color (0xAABBGGRR in memory = RGBA bytes) from rgb + alpha[0..1]
static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.0f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}

static void trigger(){
    prevX.assign((size_t)COLS*ROWS,0.f);
    currX.assign((size_t)COLS*ROWS,0.f);
    nextX.assign((size_t)COLS*ROWS,0.f);
    matterX=18; matterY=ROWS/2;
    antiX=92;   antiY=ROWS/2;
    isRunning=true; isAnnihilated=false; flashTimer=0; clockT=0;
}

extern "C" {

KEEP int  sim_w(){ return FW; }
KEEP int  sim_h(){ return FH; }

KEEP void sim_reset(){ trigger(); }

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
// re-run the annihilation crash
KEEP void sim_action(int id){ if(id==0) trigger(); }

KEEP void sim_click(double nx,double ny){
    int i=(int)(nx*COLS), j=(int)(ny*ROWS);
    if(i>0&&i<COLS-1&&j>0&&j<ROWS-1){ currX[IDX(i,j)]+=60.0; prevX[IDX(i,j)]+=60.0; }
}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        clockT++;
        // 1. movement / injection phase (before annihilation)
        if(isRunning && !isAnnihilated){
            if(clockT%2==0){ matterX+=1; antiX-=1; }
            currX[IDX(matterX,matterY)]= 40.f; prevX[IDX(matterX,matterY)]= 40.f;
            currX[IDX(antiX,antiY)]   =-40.f; prevX[IDX(antiX,antiY)]   =-40.f;
            if(matterX>=antiX){
                isAnnihilated=true; flashTimer=25;
                int cx=matterX, cy=matterY;
                currX[IDX(cx,cy)]=0.f; prevX[IDX(cx,cy)]=0.f;
                for(int dx=-1;dx<=1;++dx) for(int dy=-1;dy<=1;++dy){
                    int i=cx+dx, j=cy+dy;
                    if(i>=0&&i<COLS&&j>=0&&j<ROWS){
                        currX[IDX(i,j)]=95.f*((dx==0&&dy==0)?1.0f:-0.5f);
                        prevX[IDX(i,j)]=0.f;
                    }
                }
            }
        }
        // 2. L3 vacuum spring wave update (fixed boundaries)
        const double OMEGA_SQ=g_omega;
        for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
            double lap=currX[IDX(i+1,j)]+currX[IDX(i-1,j)]+currX[IDX(i,j+1)]+currX[IDX(i,j-1)]-4.0*currX[IDX(i,j)];
            double nv=2.0*currX[IDX(i,j)]-prevX[IDX(i,j)]+OMEGA_SQ*lap;
            nv*=g_damp;
            nextX[IDX(i,j)]=(float)nv;
        }
        // 3. buffer shift
        std::swap(prevX,currX);
        std::swap(currX,nextX);
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(3,3,8,1.f));                                 // #030308
    for(int x=0;x<FW;x+=40) olivec_rect(oc, x,0, 1,FH, rgba(5,17,0,1.f)); // #051100 grid

    // released pure-energy wave (10x gain photon green)
    for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
        if(!isAnnihilated && ((i==matterX&&j==matterY)||(i==antiX&&j==antiY))) continue;
        float amp=std::fabs(currX[IDX(i,j)]);
        float energy=amp*250.f;
        if(energy>4.f){ float a=energy/255.f; if(a>1.f)a=1.f;
            olivec_rect(oc, i*CELL, j*CELL, CELL-1, CELL-1, rgba(0,255,150,a)); }
    }

    // matter (cobalt blue +) and anti-matter (crimson -) buffers
    if(!isAnnihilated && isRunning){
        int mcx=matterX*CELL+CELL/2, mcy=matterY*CELL+CELL/2;
        olivec_circle(oc, mcx,mcy, 8, rgba(0,170,255,1.f));
        olivec_text(oc, "+", mcx-2, mcy-CELL-6, olivec_default_font, 2, rgba(255,255,255,1.f));
        int acx=antiX*CELL+CELL/2, acy=antiY*CELL+CELL/2;
        olivec_circle(oc, acx,acy, 8, rgba(255,51,102,1.f));
        olivec_text(oc, "-", acx-2, acy-CELL-6, olivec_default_font, 2, rgba(255,255,255,1.f));
    }

    // annihilation refactoring flash
    if(isAnnihilated && flashTimer>0){
        float a=flashTimer/25.f;
        int cx=matterX*CELL+CELL/2, cy=matterY*CELL+CELL/2;
        olivec_circle(oc, cx,cy, (25-flashTimer)*4, rgba(255,255,255,a));
        flashTimer--;
    }

    // system log
    const char* st = isAnnihilated ? "CRITICAL: +A + (-A) = 0. Refactoring to L3 Wave..."
                                   : "KERNEL LOG: Injecting signed particles -> merge index...";
    olivec_text(oc, st, 12, FH-20, olivec_default_font, 2, rgba(0,255,204,1.f));
    olivec_text(oc, isAnnihilated?"L1 MASS REGISTER : 0.00 (MEMORY CLEANED)"
                                 :"L1 MASS REGISTER : ALLOCATED (+A / -A)",
                12, FH-40, olivec_default_font, 2,
                isAnnihilated?rgba(204,102,255,1.f):rgba(0,170,255,1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
int main(int argc,char**argv){
    sim_init(0,0);
    int steps = argc>1? atoi(argv[1]) : 400;
    sim_step(steps);
    uint8_t* p=sim_render();
    long nb=0; uint32_t bg=rgba(3,3,8,1.f);
    for(size_t k=0;k<(size_t)FW*FH;++k) if(px[k]!=bg) nb++;
    printf("annihilation_os native: %dx%d after %d steps, non-bg px=%ld\n", FW,FH,steps,nb);
    stbi_write_png("annihilation_os_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
