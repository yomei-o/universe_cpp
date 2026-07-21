// Universe OS — Strong Interaction Kernel (Quark Confinement)  (C++/WASM port)
// Adjacent quark pairs are bound by a non-linear "flux string": the spring
// constant ω² stiffens with distance (confinement). Beyond BREAK_LIMIT the
// string snaps and its energy materialises a fresh quark/anti-quark pair
// (pair production). Drag a quark with the pointer to stretch the string.
// All state, physics AND rendering (olive.c -> RGBA) live in C++; JS only
// blits the framebuffer and forwards pointer/UI events.
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

static const int FW=800, FH=600;

// physics constants (match original)
static const double OMEGA_0=0.05, GS=0.0002, COULOMB=20.0, BREAK_LIMIT=220.0;

struct Quark { double x,y,vx,vy; uint32_t color; bool dragging; };
static std::vector<Quark> quarks;
static std::vector<uint32_t> px;

static int   g_grab=-1;      // index of grabbed quark, -1 none
static int   g_since=999;    // frames since last click (auto-release)

// xorshift PRNG in place of Math.random()
static uint32_t rng=2463534242u;
static inline double rnd(){ rng^=rng<<13; rng^=rng>>17; rng^=rng<<5; return (rng&0xFFFFFF)/(double)0x1000000; }

// pack olive color (memory 0xAABBGGRR) from rgb + alpha[0..1]
static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}
// thick line via perpendicular offset copies
static void thick_line(Olivec_Canvas oc,double x1,double y1,double x2,double y2,double w,uint32_t c){
    double dx=x2-x1,dy=y2-y1,len=std::sqrt(dx*dx+dy*dy);
    if(len<1e-6){ olivec_circle(oc,(int)x1,(int)y1,(int)(w/2),c); return; }
    double nx=-dy/len, ny=dx/len; int half=(int)(w/2); if(half<1)half=1;
    for(int o=-half;o<=half;o++)
        olivec_line(oc,(int)(x1+nx*o),(int)(y1+ny*o),(int)(x2+nx*o),(int)(y2+ny*o),c);
}

extern "C" {
KEEP int sim_w(){ return FW; }
KEEP int sim_h(){ return FH; }

KEEP void sim_reset(){
    quarks.clear();
    quarks.push_back({350,300,0,0, rgba(255,51,102,1.f), false}); // red  #ff3366
    quarks.push_back({450,300,0,0, rgba(51,102,255,1.f), false}); // blue #3366ff
    g_grab=-1; g_since=999;
}
KEEP int sim_init(int,int){ px.assign((size_t)FW*FH,0); sim_reset(); return 1; }

KEEP void sim_set(int,double){}
KEEP void sim_action(int){}

// pointer down/move: grab nearest quark within 25px then drag it to cursor
KEEP void sim_click(double nx,double ny){
    double mx=nx*FW, my=ny*FH;
    if(mx<0)mx=0; if(mx>FW)mx=FW; if(my<0)my=0; if(my>FH)my=FH;  // keep the dragged quark on-canvas
    g_since=0;
    if(g_grab<0){
        for(size_t i=0;i<quarks.size();++i){
            double d=std::hypot(quarks[i].x-mx, quarks[i].y-my);
            if(d<25){ g_grab=(int)i; break; }
        }
    }
    if(g_grab>=0 && g_grab<(int)quarks.size()){
        quarks[g_grab].x=mx; quarks[g_grab].y=my;
        quarks[g_grab].vx=0; quarks[g_grab].vy=0;
    }
}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        // auto-release grab a few frames after the last pointer event
        if(++g_since>3) g_grab=-1;
        for(size_t i=0;i<quarks.size();++i) quarks[i].dragging = ((int)i==g_grab);

        // pairwise interaction over adjacent pairs
        for(size_t i=0;i+1<quarks.size();i+=2){
            Quark &q1=quarks[i], &q2=quarks[i+1];
            double dx=q2.x-q1.x, dy=q2.y-q1.y;
            double dist=std::sqrt(dx*dx+dy*dy); if(dist<1) dist=1;
            double omega_sq=OMEGA_0+GS*dist*dist;
            double mag=omega_sq*dist + COULOMB/(dist*dist);
            double fx=(dx/dist)*mag, fy=(dy/dist)*mag;
            if(!q1.dragging){ q1.vx+=fx*0.1; q1.vy+=fy*0.1; }
            if(!q2.dragging){ q2.vx-=fx*0.1; q2.vy-=fy*0.1; }

            if(dist>BREAK_LIMIT){
                Quark a=q1, b=q2;               // copy before we rebuild
                double mx=(a.x+b.x)/2, my=(a.y+b.y)/2;
                g_grab=-1;
                Quark n1{mx-10,my,-2,0, rgba(51,102,255,1.f), false}; // partner for q1 (blue)
                Quark n2{mx+10,my, 2,0, rgba(255,51,102,1.f), false}; // partner for q2 (red)
                // faithful to original: string break refactors to exactly [q1,n1,n2,q2].
                // The whole register is reassigned, which keeps the system bounded/stable
                // (the confinement force is purely attractive and energetically unstable,
                //  so preserving all quarks would pair-produce without bound).
                quarks = { a, n1, n2, b };
                break; // recompute next frame
            }
        }

        // integrate + wall bounce.
        // The confinement spring stiffens as dist^2, so hard yanks can blow the
        // velocity/position up to huge values; feeding those into the line
        // rasteriser makes it walk billions of pixels and freezes the browser.
        // Clamp speed AND clamp the position back inside the canvas to stay bounded.
        const double VMAX=45.0;
        for(auto &q:quarks){
            if(!q.dragging){
                double sp=std::sqrt(q.vx*q.vx+q.vy*q.vy);
                if(!(sp<VMAX)){ if(sp>1e-9){ q.vx=q.vx/sp*VMAX; q.vy=q.vy/sp*VMAX; } else { q.vx=q.vy=0; } }
                q.x+=q.vx; q.y+=q.vy; q.vx*=0.95; q.vy*=0.95;
                if(q.x<20){ q.x=20; q.vx=-q.vx; } else if(q.x>FW-20){ q.x=FW-20; q.vx=-q.vx; }
                if(q.y<20){ q.y=20; q.vy=-q.vy; } else if(q.y>FH-20){ q.y=FH-20; q.vy=-q.vy; }
            }
        }
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc=olivec_canvas(px.data(),FW,FH,FW);
    olivec_fill(oc, rgba(5,5,10,1.f));                 // #05050a
    // background spacetime grid #002211
    for(int x=0;x<FW;x+=25) olivec_rect(oc,x,0,1,FH,rgba(0,34,17,1.f));
    for(int y=0;y<FH;y+=25) olivec_rect(oc,0,y,FW,1,rgba(0,34,17,1.f));

    // flux strings
    for(size_t i=0;i+1<quarks.size();i+=2){
        Quark &q1=quarks[i], &q2=quarks[i+1];
        double dist=std::hypot(q2.x-q1.x,q2.y-q1.y);
        double heat=dist/BREAK_LIMIT; if(heat>1)heat=1;
        uint32_t c=rgba((int)(heat*255),(int)((1-heat)*255),200,1.f);
        thick_line(oc,q1.x,q1.y,q2.x,q2.y, 2+heat*8, c);
    }
    // quarks
    for(auto &q:quarks){
        if(q.dragging){ uint32_t g=q.color; // dim halo
            olivec_circle(oc,(int)q.x,(int)q.y,18, (g&0x00FFFFFF)|0x50000000); }
        olivec_circle(oc,(int)q.x,(int)q.y,12,q.color);
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
    // grab a quark and yank it far to snap the string (pair production)
    sim_click(350.0/FW,300.0/FH);   // grab left quark
    sim_click( 80.0/FW,300.0/FH);   // yank it far left -> dist>BREAK_LIMIT
    sim_step(2);                    // string snaps, new pair materialises
    int steps=argc>1?atoi(argv[1]):260; sim_step(steps); // release + relax
    uint8_t*p=sim_render();
    printf("quark_os native: %dx%d, quarks=%zu\n",FW,FH,quarks.size());
    stbi_write_png("quark_os_preview.png",FW,FH,4,p,FW*4);
    return 0;
}
#endif
