// Universe OS — Graviton Emission Kernel  (C++/WASM port)
// A central binary (quadrupole moment) rotates; each frame it probabilistically
// emits pairs of quantized gravitons tangent to the orbit — neon-green (positive
// strain) and pink (negative) — that stream outward as fading spiral sparks and
// are culled at the transparent boundary. State, update AND rendering (olive.c
// -> RGBA) are all in C++; JS only blits.
#define OLIVEC_IMPLEMENTATION
#include "olive.c"
#include <vector>
#include <cstdint>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define KEEP EMSCRIPTEN_KEEPALIVE
#else
#define KEEP
#endif

static const int FW=800, FH=520;
static const double CX=FW/2.0, CY=FH/2.0;
static const double orbitRadius=35.0, orbitSpeed=0.08;
static const int COLS=140;
static const int ROWS=(int)(COLS*(double)FH/FW);            // 91
static const double DX=(double)FW/COLS, DY=(double)FH/ROWS;
static const double PR = (DX<DY?DX:DY)*0.4;                  // particle radius base

struct Grav { double x,y,vx,vy; int life; bool positive; };
static std::vector<Grav> stack;
static std::vector<uint32_t> px;
static long pulseTimer=0, totalGravitons=0;

static uint32_t rng=1234567u;
static inline double rnd(){ rng^=rng<<13; rng^=rng>>17; rng^=rng<<5; return (rng&0xFFFFFF)/(double)0x1000000; }

static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}
static void ring(Olivec_Canvas oc,double cx,double cy,double r,uint32_t c){
    int seg=(int)(r*0.8)+16; double px_=cx+r,py_=cy;
    for(int i=1;i<=seg;i++){ double a=2*M_PI*i/seg,x=cx+std::cos(a)*r,y=cy+std::sin(a)*r;
        olivec_line(oc,(int)px_,(int)py_,(int)x,(int)y,c); px_=x; py_=y; }
}

extern "C" {
KEEP int sim_w(){ return FW; }
KEEP int sim_h(){ return FH; }

KEEP void sim_reset(){ stack.clear(); pulseTimer=0; totalGravitons=0; rng=1234567u; }
KEEP int sim_init(int,int){ px.assign((size_t)FW*FH,0); sim_reset(); return 1; }

KEEP void sim_set(int,double){}
KEEP void sim_action(int){}
KEEP void sim_click(double,double){}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        pulseTimer++;
        double x1=CX+orbitRadius*std::cos(pulseTimer*orbitSpeed);
        double y1=CY+orbitRadius*std::sin(pulseTimer*orbitSpeed);
        double x2=CX-orbitRadius*std::cos(pulseTimer*orbitSpeed);
        double y2=CY-orbitRadius*std::sin(pulseTimer*orbitSpeed);
        if(rnd()<0.65){
            double a1=pulseTimer*orbitSpeed + M_PI/2;
            double speed=4.5+rnd()*1.0;
            stack.push_back({x1,y1,std::cos(a1)*speed,std::sin(a1)*speed,180,true});
            double a2=pulseTimer*orbitSpeed - M_PI/2;
            stack.push_back({x2,y2,std::cos(a2)*speed,std::sin(a2)*speed,180,false});
            totalGravitons+=2;
        }
        for(auto &p:stack){ p.x+=p.vx; p.y+=p.vy; p.life--; }
        std::vector<Grav> keep; keep.reserve(stack.size());
        for(auto &p:stack) if(p.life>0 && p.x>0 && p.x<FW && p.y>0 && p.y<FH) keep.push_back(p);
        stack.swap(keep);
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc=olivec_canvas(px.data(),FW,FH,FW);
    olivec_fill(oc, rgba(5,5,10,1.f));                         // #05050a
    for(int x=0;x<FW;x+=50) olivec_rect(oc,x,0,1,FH, rgba(0,17,5,1.f)); // #001105

    int gr=(int)(PR*0.6); if(gr<1)gr=1;
    for(auto &p:stack){
        float alpha=(float)p.life/40.f; if(alpha>1)alpha=1;
        uint32_t c = p.positive? rgba(0,255,200,alpha) : rgba(255,0,150,alpha);
        olivec_circle(oc,(int)p.x,(int)p.y,gr,c);
    }
    // central binary
    double x1=CX+orbitRadius*std::cos(pulseTimer*orbitSpeed);
    double y1=CY+orbitRadius*std::sin(pulseTimer*orbitSpeed);
    double x2=CX-orbitRadius*std::cos(pulseTimer*orbitSpeed);
    double y2=CY-orbitRadius*std::sin(pulseTimer*orbitSpeed);
    int mr=(int)(PR*1.8); if(mr<3)mr=3;
    olivec_circle(oc,(int)x1,(int)y1,mr, rgba(255,255,255,1.f));
    ring(oc,x1,y1,mr, rgba(0,255,200,0.6f));
    olivec_circle(oc,(int)x2,(int)y2,mr, rgba(255,255,255,1.f));
    ring(oc,x2,y2,mr, rgba(255,0,150,0.6f));

    char buf[160];
    snprintf(buf,sizeof buf,"GRAVITON ACTIVE REGISTERS %zu   TOTAL EMITTED %ld   CLOCK T %ld",
             stack.size(), totalGravitons, pulseTimer);
    olivec_text(oc,buf,20,FH-24,olivec_default_font,2, rgba(0,255,204,1.f));
    return (uint8_t*)px.data();
}
}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
int main(int argc,char**argv){
    sim_init(0,0);
    int steps=argc>1?atoi(argv[1]):120; sim_step(steps);
    uint8_t*p=sim_render();
    printf("graviton_spiral_os native: %dx%d, active=%zu total=%ld\n",FW,FH,stack.size(),totalGravitons);
    stbi_write_png("graviton_spiral_os_preview.png",FW,FH,4,p,FW*4);
    return 0;
}
#endif
