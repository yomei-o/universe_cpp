// Universe OS — Quantum Entanglement Kernel  (C++/WASM port)
// Two particles are emitted from a central source and drift apart. While
// unmeasured they share a common phase pointer, so anti-phase amplitudes are
// injected at both particle columns in perfect sync (sin of the system clock).
// Measuring particle A collapses a shared spin (+/-1) once, and B instantly
// mirrors the opposite fixed amplitude. The surrounding vacuum spring grid is
// solved with an explicit finite-difference wave equation. All state, the PDE
// step, AND rendering (olive.c) live in C++; JS only blits.
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

static const int COLS=120, ROWS=40, CELL=7;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 840 x 280
static std::vector<float> prevX, currX, nextX;
static std::vector<uint32_t> px;

static bool isMeasured=false;
static int  spinValue=0;
static double phaseSeed=0.0;
static double aCol=60.0, bCol=60.0;   // particle indices (A right, B left)
static long systemClock=0;

static double g_omega=0.2, g_damp=0.96;

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
    isMeasured=false; spinValue=0; phaseSeed=0.0;
    aCol=60.0; bCol=60.0; systemClock=0;
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
// measure particle A -> collapse shared spin
KEEP void sim_action(int id){
    if(id==0 && !isMeasured){ isMeasured=true; spinValue = ((double)rand()/RAND_MAX)>0.5?1:-1; }
}

KEEP void sim_click(double nx,double ny){
    int i=(int)(nx*COLS), j=(int)(ny*ROWS);
    if(i>0&&i<COLS-1&&j>0&&j<ROWS-1){ currX[IDX(i,j)]+=15.f; prevX[IDX(i,j)]+=15.f; }
}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        systemClock++;
        if(!isMeasured) phaseSeed=std::sin(systemClock*0.4);
        if(aCol<COLS-5) aCol+=0.6;
        if(bCol>5)      bCol+=-0.6;

        const double OMEGA_SQ=g_omega;
        for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
            double lap=currX[IDX(i+1,j)]+currX[IDX(i-1,j)]+currX[IDX(i,j+1)]+currX[IDX(i,j-1)]-4.0*currX[IDX(i,j)];
            double nv=2.0*currX[IDX(i,j)]-prevX[IDX(i,j)]+OMEGA_SQ*lap;
            nv*=g_damp;
            nextX[IDX(i,j)]=(float)nv;
        }

        // inject shared-state amplitude at both particle cells (into curr, pre-shift)
        int idxA=(int)aCol, idxB=(int)bCol, midY=ROWS/2;
        if(!isMeasured){
            currX[IDX(idxA,midY)]=(float)( phaseSeed*15.0);
            currX[IDX(idxB,midY)]=(float)(-phaseSeed*15.0);
        }else{
            currX[IDX(idxA,midY)]=(float)( spinValue*12.0);
            currX[IDX(idxB,midY)]=(float)(-spinValue*12.0);
        }

        std::swap(prevX,currX);
        std::swap(currX,nextX);
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(5,5,10,1.f));                               // #05050a
    for(int x=0;x<FW;x+=30) olivec_rect(oc, x,0, 1,FH, rgba(0,17,5,1.f)); // #001105 grid

    // L3 wavefront brightness (phase sign coloured)
    for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
        float amp=currX[IDX(i,j)]; float energy=std::fabs(amp)*20.f;
        if(energy>2.f){ float a=energy/255.f; if(a>1.f)a=1.f;
            uint32_t c = amp>0 ? rgba(0,255,200,a) : rgba(0,150,255,a);
            olivec_rect(oc, i*CELL, j*CELL, CELL-1, CELL-1, c); }
    }

    int midY=ROWS/2;
    int ax=(int)(aCol*CELL), bx=(int)(bCol*CELL), y=midY*CELL+CELL/2;

    // source marker (col 60)
    olivec_rect(oc, 60*CELL-10, midY*CELL-10, 20, 20, rgba(85,85,102,1.f));
    olivec_text(oc, "SOURCE", 60*CELL-18, midY*CELL-24, olivec_default_font, 1, rgba(85,85,102,1.f));

    // particle A (blue) & B (red)
    olivec_circle(oc, ax,y, 8, rgba(51,153,255,1.f));
    olivec_circle(oc, bx,y, 8, rgba(255,51,102,1.f));

    if(!isMeasured){
        olivec_text(oc, "Spin ?", ax-24, y-26, olivec_default_font, 2, rgba(51,153,255,1.f));
        olivec_text(oc, "Spin ?", bx-24, y-26, olivec_default_font, 2, rgba(255,51,102,1.f));
    }else{
        olivec_text(oc, spinValue==1?"Spin UP":"Spin DN", ax-30, y-26, olivec_default_font, 2, rgba(51,153,255,1.f));
        olivec_text(oc, spinValue==1?"Spin DN":"Spin UP", bx-30, y-26, olivec_default_font, 2, rgba(255,51,102,1.f));
    }

    char buf[96];
    olivec_text(oc, isMeasured?"INTERRUPT: A measured -> B synchronized instantly!":"STATUS: WAVE SUPERPOSITION",
                12, FH-16, olivec_default_font, 2, rgba(0,255,204,1.f));
    snprintf(buf,sizeof buf,"DISTANCE (n) : %d nodes", (int)(aCol-bCol));
    olivec_text(oc, buf, 12, FH-34, olivec_default_font, 2, rgba(0,255,204,1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
int main(int argc,char**argv){
    sim_init(0,0);
    int steps = argc>1? atoi(argv[1]) : 150;
    sim_step(steps);
    uint8_t* p=sim_render();
    long nb=0; uint32_t bg=rgba(5,5,10,1.f);
    for(size_t k=0;k<(size_t)FW*FH;++k) if(px[k]!=bg) nb++;
    printf("entanglement_os native: %dx%d after %d steps, non-bg px=%ld\n", FW,FH,steps,nb);
    stbi_write_png("entanglement_os_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
