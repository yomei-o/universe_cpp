// Universe OS — Neutron Decay Kernel (Wave-Particle Patch)  (C++/WASM port)
// Neutrons scattered in a central disk undergo random beta decay -> protons.
// Two render engines toggled by a button:
//   WAVE mode: each decay injects a pulse into a damped 2D wave field
//     (5-point Laplacian + leapfrog, absorbing border).
//   PARTICLE mode: each decay ejects an electron (pink) + antineutrino (green)
//     as free particles advanced by Euler integration with a finite lifetime.
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

static const int COLS=140, ROWS=65, CELL=6;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 840 x 390
static const double RADIUS = CELL*0.4;         // neutron/proton radius (px)

static std::vector<float> prevW, currW, nextW;
static std::vector<int>   grid;                // 0 vac / 1 neutron / 2 proton
static std::vector<uint32_t> px;

struct Particle { double x,y,vx,vy; int type; int life; }; // type 0=electron 1=neutrino
static std::vector<Particle> parts;

static bool particleMode=false;
static long pulseTimer=0, decayCount=0, neutronCount=0;
static const double DECAY_PROB=0.0003;

static std::mt19937 g_rng(777);
static inline double frand(){ return (double)g_rng()/(double)0xFFFFFFFFu; } // [0,1]

static inline int IDX(int i,int j){ return i*ROWS+j; }
static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.0f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}

static const int CX=COLS/2, CY=ROWS/2;

static void seed_neutrons(){
    grid.assign((size_t)COLS*ROWS,0);
    neutronCount=0;
    for(int i=20;i<COLS-20;++i) for(int j=10;j<ROWS-10;++j){
        double di=i-CX, dj=j-CY;
        if(std::sqrt(di*di+dj*dj)<25.0 && frand()<0.25){ grid[IDX(i,j)]=1; neutronCount++; }
    }
}

extern "C" {

KEEP int  sim_w(){ return FW; }
KEEP int  sim_h(){ return FH; }

KEEP void sim_reset(){
    prevW.assign((size_t)COLS*ROWS,0.f);
    currW.assign((size_t)COLS*ROWS,0.f);
    nextW.assign((size_t)COLS*ROWS,0.f);
    parts.clear();
    particleMode=false;
    pulseTimer=0; decayCount=0;
    g_rng.seed(777);
    seed_neutrons();
}

KEEP int sim_init(int hintW,int hintH){
    (void)hintW;(void)hintH;
    px.assign((size_t)FW*FH,0);
    sim_reset();
    return 1;
}

KEEP void sim_set(int,double){}

// toggle wave <-> particle mode
KEEP void sim_action(int id){
    if(id==0){
        particleMode=!particleMode;
        if(particleMode){ std::fill(currW.begin(),currW.end(),0.f); std::fill(prevW.begin(),prevW.end(),0.f); }
        else parts.clear();
    }
}
KEEP void sim_click(double,double){}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        pulseTimer++;
        const double OMEGA_SQ=0.25;
        for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
            if(grid[IDX(i,j)]==1){
                if(frand()<DECAY_PROB){
                    grid[IDX(i,j)]=2; decayCount++; neutronCount--;
                    if(!particleMode){
                        currW[IDX(i,j)]=20.0f;
                    } else {
                        double pX=i*CELL+CELL/2.0, pY=j*CELL+CELL/2.0;
                        double ang=frand()*2.0*M_PI;
                        double sE=3.5+frand()*1.5, sN=5.0+frand()*1.5;
                        parts.push_back({pX,pY, std::cos(ang)*sE, std::sin(ang)*sE, 0,120});
                        parts.push_back({pX,pY,-std::cos(ang)*sN,-std::sin(ang)*sN, 1,120});
                    }
                }
            }
            if(!particleMode){
                double lap=currW[IDX(i+1,j)]+currW[IDX(i-1,j)]+currW[IDX(i,j+1)]+currW[IDX(i,j-1)]-4.0*currW[IDX(i,j)];
                double nv=2.0*currW[IDX(i,j)]-prevW[IDX(i,j)]+OMEGA_SQ*lap;
                nv*=0.992;
                int de=std::min(std::min(i,COLS-1-i),std::min(j,ROWS-1-j));
                if(de<6) nv*=(de/6.0);
                nextW[IDX(i,j)]=(float)nv;
            }
        }
        if(particleMode){
            for(auto& p:parts){ p.x+=p.vx; p.y+=p.vy; p.life--; }
            std::vector<Particle> keep; keep.reserve(parts.size());
            for(auto& p:parts) if(p.life>0 && p.x>0 && p.x<FW && p.y>0 && p.y<FH) keep.push_back(p);
            parts.swap(keep);
        } else {
            for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
                prevW[IDX(i,j)]=currW[IDX(i,j)];
                currW[IDX(i,j)]=nextW[IDX(i,j)];
            }
        }
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(5,5,10,1.f));                       // #05050a
    for(int x=0;x<FW;x+=50) olivec_rect(oc, x,0, 1,FH, rgba(0,17,5,1.f));  // #001105

    if(!particleMode){
        for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
            float val=currW[IDX(i,j)], amp=std::fabs(val);
            if(amp>0.01f){
                float alpha=(float)std::min(0.25, (double)amp*1.5);
                uint32_t c = val>0 ? rgba(0,255,200,alpha) : rgba(255,0,150,alpha);
                olivec_rect(oc, i*CELL, j*CELL, CELL, CELL, c);
            }
        }
    }
    // neutrons / protons
    for(int i=0;i<COLS;++i) for(int j=0;j<ROWS;++j){
        int t=grid[IDX(i,j)]; if(t==0) continue;
        int pX=(int)(i*CELL+CELL/2.0), pY=(int)(j*CELL+CELL/2.0);
        if(t==1) olivec_circle(oc, pX,pY,(int)RADIUS, rgba(255,255,255,0.95f));
        else     olivec_circle(oc, pX,pY,(int)RADIUS, rgba(0,150,255,1.f));
    }
    if(particleMode){
        for(auto& p:parts){
            float alpha=(float)std::min(1.0, p.life/30.0);
            uint32_t c = p.type==0 ? rgba(255,0,150,alpha) : rgba(0,255,200,alpha);
            int r=(int)(RADIUS*0.5); if(r<1) r=1;
            olivec_circle(oc, (int)p.x,(int)p.y, r, c);
        }
    }
    char buf[160];
    std::snprintf(buf,sizeof(buf),"NEUTRON:%ld  PROTON:%ld  ENGINE:%s  t:%ld",
        neutronCount, decayCount, particleMode?"PARTICLE":"WAVE", pulseTimer);
    olivec_text(oc, buf, 20, FH-22, olivec_default_font, 2, rgba(0,255,204,1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
int main(int argc,char**argv){
    sim_init(0,0);
    int steps = argc>1? atoi(argv[1]) : 2500;   // wave mode: let some neutrons decay
    sim_step(steps);
    uint8_t* p=sim_render();
    int nonbg=0; uint32_t bg=rgba(5,5,10,1.f);
    for(int k=0;k<FW*FH;++k) if(px[k]!=bg) nonbg++;
    printf("neutron_decay_os native: %dx%d after %d steps, decays=%ld, non-bg px=%d\n", FW,FH,steps,decayCount,nonbg);
    stbi_write_png("neutron_decay_os_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
