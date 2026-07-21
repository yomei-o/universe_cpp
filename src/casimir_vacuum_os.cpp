// Universe OS — Pure Discrete Casimir Engine  (C++/WASM port)
// No external force formula: space is a finite rational lattice. Every frame each
// cell samples a quantum-jitter zero-point energy. Cells BETWEEN the two plates
// are filtered by a standing-wave eigenvalue mask (long wavelengths excluded),
// so their average energy drops. The lattice sums outer vs inner energy density;
// the (outer - inner) asymmetry is injected directly into the plate velocity,
// making an attractive pressure emerge natively. Sub-cell accumulators convert
// that velocity into discrete plate-index steps.
// All state, update and rendering (olive.c) live in C++; JS only blits + UI.
#define OLIVEC_IMPLEMENTATION
#include "olive.c"
#include <vector>
#include <cstdint>
#include <cmath>
#include <random>
#include <cstdio>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define KEEP EMSCRIPTEN_KEEPALIVE
#else
#define KEEP
#endif

static const int COLS=85, ROWS=50, CELL=10;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 850 x 500
static const int centerY=FH/2;
static const int plateWidth=15;

static std::vector<float> field;               // last-computed cell energy
static std::vector<uint32_t> px;

static int plateL=32, plateR=53;
static double plateVelocity=0.0, timeT=0.0;
static double bufL=0.0, bufR=0.0;
static double lastOuter=0.0, lastInner=0.0, lastDelta=0.0;

static std::mt19937 g_rng(9001);
static inline double frand(){ return (double)g_rng()/(double)0xFFFFFFFFu; } // [0,1]

static inline int IDX(int i,int j){ return i*ROWS+j; }
static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.0f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}

extern "C" {

KEEP int  sim_w(){ return FW; }
KEEP int  sim_h(){ return FH; }

KEEP void sim_reset(){
    field.assign((size_t)COLS*ROWS,0.f);
    plateL=32; plateR=53;
    plateVelocity=0.0; timeT=0.0; bufL=0.0; bufR=0.0;
    lastOuter=lastInner=lastDelta=0.0;
    g_rng.seed(9001);
}

KEEP int sim_init(int hintW,int hintH){
    (void)hintW;(void)hintH;
    px.assign((size_t)FW*FH,0);
    sim_reset();
    return 1;
}

KEEP void sim_set(int,double){}
KEEP void sim_action(int id){ if(id==0) sim_reset(); }   // reset & re-place plates
KEEP void sim_click(double,double){}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        timeT+=0.3;
        double outerSum=0, innerSum=0; int outerN=0, innerN=0;
        int cellDistance=plateR-plateL-1; if(cellDistance<1) cellDistance=1;

        for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
            double jitter=std::sin(i*0.12+timeT)*std::cos(j*0.15-timeT*0.4)*0.5+0.5;
            jitter+=frand()*0.3;
            double e=jitter;
            bool isInner=(i>plateL && i<plateR);
            if(isInner){
                double f=std::sin((i-plateL)*(M_PI/(cellDistance+1)));
                e*=std::fabs(f)*0.35;
                innerSum+=e; innerN++;
            } else {
                if(std::abs(i-plateL)==1 || std::abs(i-plateR)==1) e*=1.25;
                outerSum+=e; outerN++;
            }
            field[IDX(i,j)]=(float)e;
        }

        double avgOuter = outerN>0 ? outerSum/outerN : 0.0;
        double avgInner = innerN>0 ? innerSum/innerN : 0.0;
        lastOuter=avgOuter; lastInner=avgInner;
        double net=avgOuter-avgInner; if(net<0) net=0; lastDelta=net;

        if(cellDistance>1) plateVelocity += net*0.012;
        else plateVelocity=0;

        bufL+=plateVelocity; bufR-=plateVelocity;   // faithful to source (see note)
        if(bufL>=1){ if(plateL<plateR-1) plateL++; bufL=0; }
        if(bufR>=1){ if(plateR>plateL+1) plateR--; bufR=0; }

        if(plateL>=plateR-1){ plateL=(plateL+plateR)/2; plateR=plateL+1; plateVelocity=0; }
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(2,2,6,1.f));                    // #020206
    for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
        float alpha=(float)std::min(1.0, std::max(0.04, field[IDX(i,j)]*0.45));
        olivec_rect(oc, i*CELL+1, j*CELL+1, CELL-2, CELL-2, rgba(0,255,200,alpha));
    }

    int cellDistance=plateR-plateL-1; if(cellDistance<1) cellDistance=1;
    int xL = plateL*CELL+CELL;     // right edge of left plate
    int xR = plateR*CELL;          // left edge of right plate
    int pTop=50, pH=FH-100;

    // left plate: gradient blue(0,170,255) -> white
    for(int k=0;k<plateWidth;++k){ double t=(double)k/(plateWidth-1);
        int r=(int)(0+(255-0)*t), g=(int)(170+(255-170)*t), b=255;
        olivec_rect(oc, xL-plateWidth+k, pTop, 1, pH, rgba(r,g,b,1.f)); }
    // right plate: gradient white -> blue
    for(int k=0;k<plateWidth;++k){ double t=(double)k/(plateWidth-1);
        int r=(int)(255+(0-255)*t), g=(int)(255+(170-255)*t), b=255;
        olivec_rect(oc, xR+k, pTop, 1, pH, rgba(r,g,b,1.f)); }

    // emergent pressure arrows (pushing plates inward from outside)
    if(cellDistance>1){
        uint32_t ac=rgba(255,170,0,1.f);
        int len=(int)std::min(60.0, 5.0+lastDelta*120.0);
        int lx=xL-plateWidth;
        olivec_line(oc, lx-len,centerY, lx-2,centerY, ac);
        olivec_line(oc, lx-2,centerY, lx-10,centerY-8, ac);
        olivec_line(oc, lx-2,centerY, lx-10,centerY+8, ac);
        int rx=xR+plateWidth;
        olivec_line(oc, rx+len,centerY, rx+2,centerY, ac);
        olivec_line(oc, rx+2,centerY, rx+10,centerY-8, ac);
        olivec_line(oc, rx+2,centerY, rx+10,centerY+8, ac);
    }

    char buf[128];
    std::snprintf(buf,sizeof(buf),"GRID DISTANCE d : %d cells", cellDistance);
    olivec_text(oc, buf, 24, 24, olivec_default_font, 1, rgba(0,255,200,0.8f));
    std::snprintf(buf,sizeof(buf),"OUTER A : %.5f", lastOuter);
    olivec_text(oc, buf, 24, 36, olivec_default_font, 1, rgba(0,255,200,0.8f));
    std::snprintf(buf,sizeof(buf),"INNER B : %.5f", lastInner);
    olivec_text(oc, buf, 24, 48, olivec_default_font, 1, rgba(0,255,200,0.8f));
    std::snprintf(buf,sizeof(buf),"FORCE A-B : %.5f", lastDelta);
    olivec_text(oc, buf, 24, 60, olivec_default_font, 1, rgba(0,255,200,0.8f));
    std::snprintf(buf,sizeof(buf),"VELOCITY : %.6f step/frame", plateVelocity);
    olivec_text(oc, buf, 24, 72, olivec_default_font, 1, rgba(0,255,200,0.8f));

    const char* st = cellDistance<=1 ? "CORE EQUILIBRIUM: minimum grid resolution"
                    : (plateVelocity>0 ? "REAL SUM DELTA: outer > inner pressure" : "GRID LOCKED");
    olivec_text(oc, st, 24, FH-24, olivec_default_font, 2, rgba(0,255,204,1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
int main(int argc,char**argv){
    sim_init(0,0);
    int steps = argc>1? atoi(argv[1]) : 120;
    sim_step(steps);
    uint8_t* p=sim_render();
    int nonbg=0; uint32_t bg=rgba(2,2,6,1.f);
    for(int k=0;k<FW*FH;++k) if(px[k]!=bg) nonbg++;
    printf("casimir_vacuum_os native: %dx%d after %d steps, plateL=%d plateR=%d, non-bg px=%d\n",
        FW,FH,steps,plateL,plateR,nonbg);
    stbi_write_png("casimir_vacuum_os_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
