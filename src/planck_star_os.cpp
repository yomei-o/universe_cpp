// Universe OS — Planck Star Safety Kernel  (C++/WASM port)
// A 2D finite-difference wave grid. Triggering COLLAPSE injects a spherical
// shell of mass and applies a pseudo-gravity nudge pulling every cell toward
// the center. As density concentrates, a Planck safety clamp caps the central
// cells at +/-45 and inverts their phase, so instead of a singularity the core
// bounces energy back outward (a Planck star). All state, PDE step AND
// rendering (olive.c into RGBA) live in C++; JS only blits.
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

static const int COLS=100, ROWS=62, CELL=8;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 800 x 496
static const int midX=COLS/2, midY=ROWS/2;     // 50, 31
static const double PLANCK_LIMIT=45.0;
static const double PI=3.14159265358979323846;

static std::vector<float> prevX, currX, nextX;
static std::vector<uint32_t> px;
static bool isCollapsing=false, planckCritical=false;
static long collapseTimer=0;
static double g_omega=0.25, g_damp=0.985;

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
    isCollapsing=false; planckCritical=false; collapseTimer=0;
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

// action 0: trigger gravity collapse sequence
KEEP void sim_action(int id){
    if(id==0 && !isCollapsing){
        isCollapsing=true; collapseTimer=0; planckCritical=false;
        for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
            double dist=std::sqrt((double)((i-midX)*(i-midX)+(j-midY)*(j-midY)));
            if(dist>15 && dist<25){
                float v=12.0f*std::sin(dist*0.5);
                currX[IDX(i,j)]=v; prevX[IDX(i,j)]=v;
            }
        }
    }
}

KEEP void sim_click(double nx,double ny){
    int i=(int)(nx*COLS), j=(int)(ny*ROWS);
    if(i>1&&i<COLS-1&&j>1&&j<ROWS-1){ currX[IDX(i,j)]=15.f; prevX[IDX(i,j)]=15.f; }
}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        if(isCollapsing){
            collapseTimer++;
            if(collapseTimer<150){
                for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
                    double dist=std::sqrt((double)((i-midX)*(i-midX)+(j-midY)*(j-midY)));
                    if(dist>1){
                        double push=0.18/dist;
                        currX[IDX(i,j)] += (midX-i)*push*0.05;
                        currX[IDX(i,j)] += (midY-j)*push*0.05;
                    }
                }
            }
        }
        planckCritical=false;
        const double OMEGA=g_omega;
        for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
            double lap = currX[IDX(i+1,j)]+currX[IDX(i-1,j)]+currX[IDX(i,j+1)]+currX[IDX(i,j-1)]-4.0*currX[IDX(i,j)];
            double nv = 2.0*currX[IDX(i,j)] - prevX[IDX(i,j)] + OMEGA*lap;
            nv *= g_damp;
            double dc=std::sqrt((double)((i-midX)*(i-midX)+(j-midY)*(j-midY)));
            if(dc<=3){
                if(std::fabs(nv)>=PLANCK_LIMIT){
                    planckCritical=true;
                    nv = (nv>0?PLANCK_LIMIT:-PLANCK_LIMIT)*-0.8;  // clip + phase invert
                }
            }
            nextX[IDX(i,j)] = (float)nv;
        }
        std::swap(prevX,currX);
        std::swap(currX,nextX);
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(3,3,7,1.f));
    for(int x=0;x<FW;x+=40) olivec_rect(oc, x,0, 1,FH, rgba(0x11,0x05,0x05,1.f));

    for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
        float energy=std::fabs(currX[IDX(i,j)])*280.f;
        if(energy>4.f){ float a=energy/255.f; if(a>1.f)a=1.f;
            double dist=std::sqrt((double)((i-midX)*(i-midX)+(j-midY)*(j-midY)));
            uint32_t c;
            if(dist<=3 && planckCritical)      c=rgba(255,255,255,a);
            else if(dist<=3)                   c=rgba(255,50,50,a);
            else                               c=rgba(255,140,0,a*0.8f);
            olivec_rect(oc, i*CELL, j*CELL, CELL-1, CELL-1, c); }
    }

    // event horizon prediction ring (dashed), radius 12 cells
    if(isCollapsing){
        int cx=midX*CELL+CELL/2, cy=midY*CELL+CELL/2, R=12*CELL;
        uint32_t hc = planckCritical ? rgba(255,50,50,0.4f) : rgba(150,50,50,0.2f);
        int N=96;
        for(int k=0;k<N;k+=2){
            double a0=2*PI*k/N, a1=2*PI*(k+1)/N;
            olivec_line(oc, cx+(int)(R*std::cos(a0)), cy+(int)(R*std::sin(a0)),
                            cx+(int)(R*std::cos(a1)), cy+(int)(R*std::sin(a1)), hc);
        }
        olivec_text(oc, "EVENT HORIZON", cx-45, cy-R-14, olivec_default_font, 1, rgba(255,100,100,0.8f));
    }

    const char* status = planckCritical ? "HARDWARE PATCH: PLANCK CRITICAL! Energy bounced."
                        : (isCollapsing && collapseTimer>=150) ? "LOG: Core collapse halted. Planck Star stable."
                        : "SYSTEM IDLE (vacuum register stable)";
    olivec_text(oc, status, 16, FH-16, olivec_default_font, 1, rgba(0,255,200,1.f));
    olivec_text(oc, planckCritical?"CORE DENSITY: 1.00 PLANCK LIMIT (MAX)":"CORE DENSITY: SUB-CRITICAL STATE",
                16, FH-32, olivec_default_font, 1, planckCritical?rgba(255,51,51,1.f):rgba(0,255,200,1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
int main(int argc,char**argv){
    sim_init(0,0);
    sim_action(0);        // trigger collapse
    int steps = argc>1? atoi(argv[1]) : 120;
    sim_step(steps);
    uint8_t* p=sim_render();
    printf("planck_star_os native: %dx%d after %d steps\n", FW,FH,steps);
    stbi_write_png("planck_star_os_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
