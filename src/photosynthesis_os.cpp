// Universe OS — Photosynthesis / quantum coherent energy transport (C++/WASM)
// A pure 2D wave equation on a grid whose per-cell damping is raised along
// "chlorophyll channels" (waveguides) around randomly placed nodes, so a sine
// wave injected at the antenna (Source, left) is funnelled with low loss to the
// reaction centre (Sink, right), which absorbs the arriving amplitude.
// All state, physics AND rendering (olive.c -> RGBA) are C++; JS only blits + UI.
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

static const int CELL=8, COLS=100, ROWS=63;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 800 x 504
static std::vector<double> curX, prvX, nxtX, damp;
static std::vector<uint32_t> px;

struct Node{ double x,y; int kind; };           // kind: 0 node, 1 source, 2 sink
static std::vector<Node> nodes;
static long total_time=0;
static double energy=0.0;
static double g_coherence=0.99;

static uint32_t rng=88172645u;
static inline double rnd(){ rng^=rng<<13; rng^=rng>>17; rng^=rng<<5; return (rng&0xFFFFFF)/(double)0x1000000; }

static inline int IDX(int i,int j){ return i*ROWS+j; }
static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}

static void buildDamping(){
    for(auto&d:damp) d=0.84;
    for(const Node&n:nodes){
        int gi=(int)(n.x/CELL), gj=(int)(n.y/CELL);
        int radius = (n.kind==2)?8:6;
        for(int di=-radius;di<=radius;++di) for(int dj=-radius;dj<=radius;++dj){
            double dist=std::sqrt((double)di*di+dj*dj);
            int ii=gi+di, jj=gj+dj;
            if(ii>0&&ii<COLS-1&&jj>0&&jj<ROWS-1&&dist<=radius) damp[IDX(ii,jj)]=g_coherence;
        }
    }
}

extern "C" {
KEEP int sim_w(){ return FW; }
KEEP int sim_h(){ return FH; }

KEEP void sim_reset(){
    curX.assign((size_t)COLS*ROWS,0.0); prvX.assign((size_t)COLS*ROWS,0.0);
    nxtX.assign((size_t)COLS*ROWS,0.0); damp.assign((size_t)COLS*ROWS,0.84);
    nodes.clear(); total_time=0; energy=0.0;
    nodes.push_back({80.0, FH/2.0, 1});            // Source (antenna)
    nodes.push_back({FW-80.0, FH/2.0, 2});         // Sink (reaction centre)
    int num=(int)((double)FW*FH/6500.0);
    for(int i=0;i<num;++i) nodes.push_back({130.0+rnd()*(FW-260), 50.0+rnd()*(FH-100), 0});
    buildDamping();
}
KEEP int sim_init(int,int){ px.assign((size_t)FW*FH,0); sim_reset(); return 1; }

KEEP void sim_set(int id,double v){ if(id==0){ g_coherence=v; buildDamping(); } }
KEEP void sim_action(int id){ if(id==0) sim_reset(); }      // regenerate network
KEEP void sim_click(double,double){}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        int sgi=(int)(80.0/CELL), sgj=(int)((FH/2.0)/CELL);
        int egi=(int)((FW-80.0)/CELL), egj=(int)((FH/2.0)/CELL);
        curX[IDX(sgi,sgj)] = std::sin(total_time*0.4)*15.0;    // inject
        const double OMEGA=0.24;
        for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
            double lap=curX[IDX(i+1,j)]+curX[IDX(i-1,j)]+curX[IDX(i,j+1)]+curX[IDX(i,j-1)]-4.0*curX[IDX(i,j)];
            double nv=(2.0*curX[IDX(i,j)]-prvX[IDX(i,j)])+OMEGA*lap;
            nxtX[IDX(i,j)]=nv*damp[IDX(i,j)];
        }
        double arrived=std::fabs(curX[IDX(egi,egj)]);
        if(arrived>0.05) energy += arrived*0.02;
        std::swap(prvX,curX); std::swap(curX,nxtX);
        total_time++;
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc=olivec_canvas(px.data(),FW,FH,FW);
    olivec_fill(oc, rgba(2,2,6,1.f));
    // wave field (min brightness boosted to 0.15)
    for(int i=1;i<COLS;++i) for(int j=1;j<ROWS;++j){
        double val=curX[IDX(i,j)], a=std::fabs(val);
        if(a>0.01){ float in=0.15f+(float)(std::min(1.0,a/4.0)*0.85);
            uint32_t c = val>0 ? rgba(255,230,0,in) : rgba(0,255,180,in);
            olivec_rect(oc, i*CELL-CELL/2, j*CELL-CELL/2, CELL, CELL, c); }
    }
    // faint waveguide topology
    for(int i=1;i<COLS;++i) for(int j=1;j<ROWS;++j)
        if(damp[IDX(i,j)]>0.9) olivec_rect(oc, i*CELL-CELL/2, j*CELL-CELL/2, CELL, CELL, rgba(0,255,100,0.02f));
    // nodes
    for(const Node&n:nodes){
        if(n.kind==1) olivec_circle(oc,(int)n.x,(int)n.y,8, rgba(0,170,255,1.f));
        else if(n.kind==2) olivec_circle(oc,(int)n.x,(int)n.y,12, rgba(255,51,102,1.f));
        else olivec_circle(oc,(int)n.x,(int)n.y,3, rgba(0,255,100,0.4f));
    }
    char buf[96]; snprintf(buf,sizeof(buf),"WAVE ENERGY ABSORBED  %.2f",energy);
    olivec_text(oc, buf, 16, FH-22, olivec_default_font, 2, rgba(0,255,204,1.f));
    return (uint8_t*)px.data();
}
}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
int main(int argc,char**argv){
    sim_init(0,0);
    int steps=argc>1?atoi(argv[1]):600; sim_step(steps);
    uint8_t*p=sim_render();
    printf("photosynthesis_os native: %dx%d after %d steps, energy=%.2f\n",FW,FH,steps,energy);
    stbi_write_png("photosynthesis_os_preview.png",FW,FH,4,p,FW*4);
    return 0;
}
#endif
