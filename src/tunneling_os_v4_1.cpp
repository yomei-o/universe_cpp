// Universe OS — Quantum Tunneling Kernel V4.1  (C++/WASM port)
// A Gaussian wave packet is fired periodically from the left edge across a 2D
// vacuum spring grid solved with an explicit finite-difference wave equation. A
// barrier column attenuates any passing amplitude (*0.15 + small jitter); the
// faint leaked field beyond the barrier is boosted 10x for display. A detector
// column near the right counts successful tunnel events. All state, the PDE step,
// AND rendering (olive.c) live in C++; JS only blits.
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
#include <cstdlib>
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define KEEP EMSCRIPTEN_KEEPALIVE
#else
#define KEEP
#endif

static const int COLS=100, ROWS=40, CELL=8;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 800 x 320
static std::vector<float> prevX, currX, nextX;
static std::vector<uint32_t> px;

static const int BARRIER_COL=45;
static long tunnelCount=0, totalShot=0, pulseTimer=0;

static double g_omega=0.25, g_barrier=0.15;

static inline int IDX(int i,int j){ return i*ROWS+j; }
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
    tunnelCount=0; totalShot=0; pulseTimer=0;
}

KEEP int sim_init(int hintW,int hintH){
    (void)hintW;(void)hintH;
    px.assign((size_t)FW*FH,0);
    sim_reset();
    return 1;
}

KEEP void sim_set(int id,double v){
    if(id==0) g_omega=v;
    else if(id==1) g_barrier=v;
}
KEEP void sim_action(int){}

KEEP void sim_click(double nx,double ny){
    int i=(int)(nx*COLS), j=(int)(ny*ROWS);
    if(i>0&&i<COLS-1&&j>0&&j<ROWS-1){ currX[IDX(i,j)]+=25.f; prevX[IDX(i,j)]+=25.f; }
}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        pulseTimer++;
        // fire a Gaussian packet from the left edge
        if(pulseTimer%45==1){
            totalShot++; int midY=ROWS/2;
            for(int dy=-2;dy<=2;++dy){ int j=midY+dy;
                float amp=25.f*std::exp(-(dy*dy)/2.f);
                currX[IDX(3,j)]=amp; prevX[IDX(3,j)]=amp; }
        }

        const double OMEGA_SQ=g_omega;
        for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
            double lap=currX[IDX(i+1,j)]+currX[IDX(i-1,j)]+currX[IDX(i,j+1)]+currX[IDX(i,j-1)]-4.0*currX[IDX(i,j)];
            double nv=2.0*currX[IDX(i,j)]-prevX[IDX(i,j)]+OMEGA_SQ*lap;
            if(i==BARRIER_COL){ nv*=g_barrier; nv+=((double)rand()/RAND_MAX-0.5)*0.12; }
            nextX[IDX(i,j)]=(float)nv;
        }

        // detector column near the right
        int DETECT_COL=COLS-15, midY=ROWS/2;
        if(pulseTimer%45==24){
            double mx=0;
            for(int dy=-5;dy<=5;++dy) mx=std::max(mx,(double)std::fabs(currX[IDX(DETECT_COL,midY+dy)]));
            if(mx>0.08) tunnelCount++;
        }

        std::swap(prevX,currX);
        std::swap(currX,nextX);
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(5,5,10,1.f));                               // #05050a
    for(int x=0;x<FW;x+=40) olivec_rect(oc, x,0, 1,FH, rgba(0,17,5,1.f)); // #001105 grid

    // field brightness (10x boost beyond the barrier)
    for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
        float amp=std::fabs(currX[IDX(i,j)]);
        float energy = (i>BARRIER_COL) ? amp*250.f : amp*35.f;
        if(energy>5.f){ float a=energy/255.f; if(a>1.f)a=1.f;
            olivec_rect(oc, i*CELL, j*CELL, CELL-1, CELL-1, rgba(0,255,200,a)); }
    }

    // barrier wall
    int bx=BARRIER_COL*CELL;
    olivec_rect(oc, bx,0, CELL,FH, rgba(255,30,80,0.4f));
    olivec_rect(oc, bx,0, CELL,2,  rgba(255,30,80,1.f));
    olivec_rect(oc, bx,FH-2, CELL,2, rgba(255,30,80,1.f));

    int midY=ROWS/2, DETECT_COL=COLS-15;
    // electron gun
    olivec_rect(oc, 0,(midY-2)*CELL, 15, 4*CELL, rgba(136,136,153,1.f));
    // detector box (#00aaff)
    olivec_rect(oc, DETECT_COL*CELL,(midY-5)*CELL, 15,2, rgba(0,170,255,1.f));
    olivec_rect(oc, DETECT_COL*CELL,(midY+5)*CELL, 15,2, rgba(0,170,255,1.f));
    olivec_rect(oc, DETECT_COL*CELL,(midY-5)*CELL, 2,10*CELL, rgba(0,170,255,1.f));
    olivec_rect(oc, DETECT_COL*CELL+13,(midY-5)*CELL, 2,10*CELL, rgba(0,170,255,1.f));

    char buf[64];
    snprintf(buf,sizeof buf,"TOTAL SHOTS   : %ld", totalShot);
    olivec_text(oc, buf, 20, FH-52, olivec_default_font, 2, rgba(0,255,204,1.f));
    snprintf(buf,sizeof buf,"TUNNELED (*) : %ld", tunnelCount);
    olivec_text(oc, buf, 20, FH-34, olivec_default_font, 2, rgba(0,255,204,1.f));
    if(totalShot>0){ snprintf(buf,sizeof buf,"TUNNEL RATE   : %.1f %%", 100.0*tunnelCount/totalShot);
        olivec_text(oc, buf, 20, FH-16, olivec_default_font, 2, rgba(0,255,204,1.f)); }
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
int main(int argc,char**argv){
    sim_init(0,0);
    int steps = argc>1? atoi(argv[1]) : 200;
    sim_step(steps);
    uint8_t* p=sim_render();
    long nb=0; uint32_t bg=rgba(5,5,10,1.f);
    for(size_t k=0;k<(size_t)FW*FH;++k) if(px[k]!=bg) nb++;
    printf("tunneling_os_v4_1 native: %dx%d after %d steps, non-bg px=%ld\n", FW,FH,steps,nb);
    stbi_write_png("tunneling_os_v4_1_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
