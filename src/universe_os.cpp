// Universe OS — 3-Layer Kernel Mainframe  (C++/WASM port)
// A small grid of INDEPENDENT damped oscillators (NOT a Laplacian):
//   X_{n+1} = (2 - omega^2 * dtau^2) X_n - X_{n-1}
// A clicked "black hole" (high-load sector) slows the local sampling rate dtau
// via a distance-based compute load, freezing nearby cells. The grid is rendered
// as a deformable MESH: vertical + horizontal polylines whose vertices are
// displaced by the oscillator value and pulled toward the black hole.
// All state, update and rendering (olive.c) live in C++; JS only blits + UI.
#define OLIVEC_IMPLEMENTATION
#include "olive.c"
#include <vector>
#include <cstdint>
#include <cmath>
#include <random>
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define KEEP EMSCRIPTEN_KEEPALIVE
#else
#define KEEP
#endif

static const int COLS=40, ROWS=30, CELL=20;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 800 x 600
static const double SPX=CELL, SPY=CELL;
static const double C_BASE=137.0, OMEGA_0=0.2;

static std::vector<float> prevX, currX, nextX;
static std::vector<uint32_t> px;
static long total_time=1;
static bool bhOn=false;
static int bhCol=0, bhRow=0;
static const double BH_ENERGY=15.0;

static std::mt19937 g_rng(2024);
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
    prevX.assign((size_t)COLS*ROWS,0.f);
    currX.assign((size_t)COLS*ROWS,0.f);
    nextX.assign((size_t)COLS*ROWS,0.f);
    total_time=1; bhOn=false;
    g_rng.seed(2024);
}

KEEP int sim_init(int hintW,int hintH){
    (void)hintW;(void)hintH;
    px.assign((size_t)FW*FH,0);
    sim_reset();
    return 1;
}

KEEP void sim_set(int,double){}
KEEP void sim_action(int){}

// click places the black hole on the nearest grid cell
KEEP void sim_click(double nx,double ny){
    int c=(int)(nx*COLS), r=(int)(ny*ROWS);
    if(c<0)c=0; if(c>=COLS)c=COLS-1; if(r<0)r=0; if(r>=ROWS)r=ROWS-1;
    bhCol=c; bhRow=r; bhOn=true;
}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        double current_c = C_BASE/(1.0+total_time*0.0001);
        double omega_sq = OMEGA_0*OMEGA_0;
        for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
            double dtau=1.0;
            if(bhOn){
                double di=i-bhCol, dj=j-bhRow;
                double dn=std::sqrt(di*di+dj*dj)+0.5;
                double load=BH_ENERGY/std::pow(current_c*dn,1.5);
                dtau=std::max(0.0,1.0-load);
            }
            double nv=(2.0-omega_sq*dtau*dtau)*currX[IDX(i,j)]-prevX[IDX(i,j)];
            if(frand()<0.001) nv+=(frand()-0.5)*5.0;
            nextX[IDX(i,j)]=(float)nv;
        }
        for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
            prevX[IDX(i,j)]=currX[IDX(i,j)];
            currX[IDX(i,j)]=nextX[IDX(i,j)];
        }
        total_time++;
    }
}

// map grid vertex (i,j) to displaced mesh coordinate
static inline void vertex(int i,int j,double& dx,double& dy){
    double baseX=i*SPX+SPX/2.0;
    double baseY=j*SPY+SPY/2.0 + currX[IDX(i,j)];
    dx=baseX; dy=baseY;
    if(bhOn){
        double bhX=bhCol*SPX+SPX/2.0, bhY=bhRow*SPY+SPY/2.0;
        double d=std::sqrt((dx-bhX)*(dx-bhX)+(dy-bhY)*(dy-bhY));
        double pull=(BH_ENERGY*50.0)/(d+20.0);
        dx += ((bhX-dx)/(d+1.0))*pull;
        dy += ((bhY-dy)/(d+1.0))*pull;
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(5,5,10,1.f));                   // #05050a
    uint32_t mesh=rgba(0,85,68,1.f);                     // #005544
    // vertical polylines
    for(int i=0;i<COLS;++i){
        double px0=0,py0=0;
        for(int j=0;j<ROWS;++j){ double x,y; vertex(i,j,x,y);
            if(j>0) olivec_line(oc,(int)px0,(int)py0,(int)x,(int)y,mesh);
            px0=x; py0=y; }
    }
    // horizontal polylines
    for(int j=0;j<ROWS;++j){
        double px0=0,py0=0;
        for(int i=0;i<COLS;++i){ double x,y; vertex(i,j,x,y);
            if(i>0) olivec_line(oc,(int)px0,(int)py0,(int)x,(int)y,mesh);
            px0=x; py0=y; }
    }
    // black hole core
    if(bhOn){
        int bx=(int)(bhCol*SPX+SPX/2.0), by=(int)(bhRow*SPY+SPY/2.0);
        olivec_circle(oc, bx,by, 14, rgba(255,51,102,0.25f)); // glow
        olivec_circle(oc, bx,by, 10, rgba(255,51,102,1.f));   // #ff3366
    }
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
int main(int argc,char**argv){
    sim_init(0,0);
    sim_click(0.5,0.5);                    // drop a black hole in the center
    int steps = argc>1? atoi(argv[1]) : 200;
    sim_step(steps);
    uint8_t* p=sim_render();
    int nonbg=0; uint32_t bg=rgba(5,5,10,1.f);
    for(int k=0;k<FW*FH;++k) if(px[k]!=bg) nonbg++;
    printf("universe_os native: %dx%d after %d steps, non-bg px=%d\n", FW,FH,steps,nonbg);
    stbi_write_png("universe_os_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
