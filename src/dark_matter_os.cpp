// Universe OS — Galactic Strain Kernel (Dark Matter)  (C++/WASM port)
// 160 stars orbit a galactic core on Keplerian curves (v = 0.8/sqrt(r)). In the
// unpatched/"buggy" mode the outer stars lose gravity and fly apart; toggling the
// patch pre-allocates a vacuum "strain buffer" that flattens the rotation curve
// (the dark-matter analogue) and locks the outer stars into fast, flat rotation.
// All state, update AND rendering (olive.c -> RGBA) are in C++; JS only blits.
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

static const int FW=800, FH=500;
static const double CX=FW/2.0, CY=FH/2.0;
static const int numStars=160;

struct Star { double r,angle,size,baseSpeed; };
static std::vector<Star> stars;
static std::vector<uint32_t> px;
static bool patched=false;

static uint32_t rng=88172645u;
static inline double rnd(){ rng^=rng<<13; rng^=rng>>17; rng^=rng<<5; return (rng&0xFFFFFF)/(double)0x1000000; }

static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}
// hollow ring outline (olivec_circle only fills disks)
static void ring(Olivec_Canvas oc,double cx,double cy,double r,uint32_t c){
    int seg=(int)(r*0.6)+16;
    double px_=cx+r, py_=cy;
    for(int i=1;i<=seg;i++){ double a=2*M_PI*i/seg; double x=cx+std::cos(a)*r, y=cy+std::sin(a)*r;
        olivec_line(oc,(int)px_,(int)py_,(int)x,(int)y,c); px_=x; py_=y; }
}
static void dashed_line(Olivec_Canvas oc,double x1,double y1,double x2,double y2,double dash,double gap,uint32_t c){
    double dx=x2-x1,dy=y2-y1,len=std::sqrt(dx*dx+dy*dy); if(len<1e-6)return;
    double ux=dx/len,uy=dy/len,pos=0; bool on=true;
    while(pos<len){ double e=std::min(len,pos+(on?dash:gap));
        if(on) olivec_line(oc,(int)(x1+ux*pos),(int)(y1+uy*pos),(int)(x1+ux*e),(int)(y1+uy*e),c);
        pos=e; on=!on; }
}
static void thick_line(Olivec_Canvas oc,double x1,double y1,double x2,double y2,double w,uint32_t c){
    double dx=x2-x1,dy=y2-y1,len=std::sqrt(dx*dx+dy*dy);
    if(len<1e-6)return; double nx=-dy/len,ny=dx/len; int half=(int)(w/2); if(half<1)half=1;
    for(int o=-half;o<=half;o++) olivec_line(oc,(int)(x1+nx*o),(int)(y1+ny*o),(int)(x2+nx*o),(int)(y2+ny*o),c);
}
static void quad(Olivec_Canvas oc,double x0,double y0,double cxp,double cyp,double x1,double y1,uint32_t c){
    double pxp=x0,pyp=y0;
    for(int s=1;s<=24;s++){ double t=s/24.0,u=1-t;
        double bx=u*u*x0+2*u*t*cxp+t*t*x1, by=u*u*y0+2*u*t*cyp+t*t*y1;
        thick_line(oc,pxp,pyp,bx,by,2,c); pxp=bx; pyp=by; }
}

extern "C" {
KEEP int sim_w(){ return FW; }
KEEP int sim_h(){ return FH; }

KEEP void sim_reset(){
    rng=88172645u; stars.clear();
    for(int i=0;i<numStars;i++){
        double r=rnd()*180+20;
        stars.push_back({ r, rnd()*2*M_PI, rnd()*2+1, 0.8/std::sqrt(r) });
    }
    patched=false;
}
KEEP int sim_init(int,int){ px.assign((size_t)FW*FH,0); sim_reset(); return 1; }

KEEP void sim_set(int,double){}
KEEP void sim_action(int id){ if(id==0) patched=!patched; }
KEEP void sim_click(double,double){}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        for(auto &st:stars){
            double cur=st.baseSpeed;
            if(patched){
                if(st.r>60) cur=0.045+rnd()*0.005;
            } else {
                if(st.r>100) st.r+=0.2;
            }
            st.angle+=cur;
            if(!patched && st.r>300) st.r=rnd()*80+30;
        }
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc=olivec_canvas(px.data(),FW,FH,FW);
    olivec_fill(oc, rgba(2,2,5,1.f));                          // #020205

    // 1. concentric spacetime mesh
    uint32_t mesh = patched? rgba(0,255,150,0.15f) : rgba(0,150,255,0.05f);
    for(double r=40;r<240;r+=40) ring(oc,CX,CY,r,mesh);

    // strain buffer radial spokes when patched
    if(patched){
        uint32_t spoke=rgba(0,255,150,0.2f);
        for(double a=0;a<2*M_PI-1e-6;a+=M_PI/6)
            dashed_line(oc,CX,CY,CX+std::cos(a)*240,CY+std::sin(a)*240,4,12,spoke);
    }

    // 2. galactic core: radial-gradient approximation (transparent cyan -> white)
    for(int rr=25;rr>=1;--rr){
        double f=rr/25.0; int cr,cg,cb; float a;
        if(f<=0.4){ double t=f/0.4; cr=(int)(255*(1-t)); cg=255; cb=255; a=1.f; }
        else { double t=(f-0.4)/0.6; cr=0; cg=255; cb=255; a=(float)(1.0-t); }
        olivec_circle(oc,(int)CX,(int)CY,rr, rgba(cr,cg,cb,a));
    }

    // 3. stars
    for(auto &st:stars){
        double sx=CX+std::cos(st.angle)*st.r, sy=CY+std::sin(st.angle)*st.r;
        uint32_t c;
        if(patched) c = st.r>100? rgba(100,255,200,0.9f) : rgba(255,255,255,1.f);
        else        c = st.r>120? rgba(255,100,100,0.7f) : rgba(100,200,255,0.9f);
        int sz=(int)st.size; if(sz<1)sz=1;
        olivec_rect(oc,(int)sx,(int)sy,sz,sz,c);
    }

    // 4. rotation-curve graph (bottom-left)
    int gx=30, gy=FH-160, gw=220, gh=110;
    olivec_rect(oc,gx,gy,gw,gh, rgba(0,10,20,0.8f));
    olivec_frame(oc,gx,gy,gw,gh,1, rgba(0,255,204,1.f));
    olivec_text(oc,"ROTATION CURVE V VS R", gx+10, gy+6, olivec_default_font,1, rgba(0,255,204,1.f));
    uint32_t axc=rgba(85,85,85,1.f);
    olivec_line(oc,50,FH-70,230,FH-70,axc);   // R axis
    olivec_line(oc,50,FH-70,50,FH-140,axc);   // V axis
    if(patched) quad(oc,50,FH-120,100,FH-115,230,FH-115, rgba(51,255,51,1.f));   // flat
    else        quad(oc,50,FH-120, 90,FH-110,230,FH-75,  rgba(255,51,51,1.f));   // drop-off

    // 5. system log
    const char* msg = patched
        ? "PATCH 10 APPLIED. STRAIN REGISTER PRE-ALLOCATED. COHERENCE 1.00"
        : "CRITICAL WARNING. GALACTIC OUTER REGISTER DISRUPTION. TENSION LOST.";
    olivec_text(oc,msg,20,FH-24,olivec_default_font,2, rgba(0,255,204,1.f));
    return (uint8_t*)px.data();
}
}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
int main(int argc,char**argv){
    sim_init(0,0);
    int steps=argc>1?atoi(argv[1]):200; sim_step(steps);
    uint8_t*p=sim_render();
    printf("dark_matter_os native: %dx%d, patched=%d\n",FW,FH,(int)patched);
    stbi_write_png("dark_matter_os_preview.png",FW,FH,4,p,FW*4);
    return 0;
}
#endif
