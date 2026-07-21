// Universe OS — Weak Interaction / Beta-Decay Kernel  (C++/WASM port)
// ~50 particles drift with digital friction and wall bounces. When a neutron's
// kinetic energy (|v|) exceeds V_WEAK it flashes white for 10 frames, then
// decays into a proton (heavier, slowed) and ejects a fresh electron at speed
// 8 in a random direction. Pairwise proximity collisions push particles apart
// and swap velocities. Clicking injects a radial energy shockwave. When protons
// exceed 40 the vacuum resets. All state, update AND rendering (olive.c into an
// RGBA framebuffer) live in C++; JS only blits the buffer and forwards events.
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
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define KEEP EMSCRIPTEN_KEEPALIVE
#else
#define KEEP
#endif

static const int FW=800, FH=600;

// tunable params
static double g_vweak=3.5;    // decay energy (speed) threshold
static double g_friction=0.98;

enum PType { NEUTRON=0, PROTON=1, ELECTRON=2 };
struct Particle { double x,y,vx,vy; int type; double radius; int flash; };

static std::vector<Particle> parts;
static std::vector<uint32_t> px;

static uint32_t rng_state=0x9e3779b9u;
static inline double frand(){
    rng_state^=rng_state<<13; rng_state^=rng_state>>17; rng_state^=rng_state<<5;
    return (rng_state&0xFFFFFF)/(double)0x1000000;
}

static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.0f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}

static void initParticles(){
    parts.clear(); parts.reserve(64);
    for(int i=0;i<50;++i){
        Particle p;
        p.x=frand()*(FW-100.0)+50.0;
        p.y=frand()*(FH-100.0)+50.0;
        p.vx=(frand()-0.5)*2.0;
        p.vy=(frand()-0.5)*2.0;
        p.type=NEUTRON; p.radius=10.0; p.flash=0;
        parts.push_back(p);
    }
}

extern "C" {

KEEP int  sim_w(){ return FW; }
KEEP int  sim_h(){ return FH; }

KEEP void sim_reset(){ rng_state=0x9e3779b9u; initParticles(); }

KEEP int sim_init(int,int){
    px.assign((size_t)FW*FH,0);
    sim_reset();
    return 1;
}

KEEP void sim_set(int id,double v){
    if(id==0) g_vweak=v;
    else if(id==1) g_friction=v;
}
KEEP void sim_action(int){}

// click injects a radial energy shockwave
KEEP void sim_click(double nx,double ny){
    double mx=nx*FW, my=ny*FH;
    for(auto &p: parts){
        double dx=p.x-mx, dy=p.y-my;
        double dist=std::sqrt(dx*dx+dy*dy);
        if(dist<150){
            double force=(150.0-dist)*0.08;
            p.vx+=(dx/(dist+1.0))*force;
            p.vy+=(dy/(dist+1.0))*force;
        }
    }
}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        std::vector<Particle> spawned;
        size_t n=parts.size();
        for(size_t i=0;i<n;++i){
            Particle &p=parts[i];
            double energy=std::sqrt(p.vx*p.vx+p.vy*p.vy);
            // weak-force step function: arm a decay flash
            if(p.type==NEUTRON && energy>g_vweak && p.flash==0) p.flash=10;
            // decay completes on the last flash frame
            if(p.flash==1){
                p.type=PROTON; p.radius=15.0;
                p.vx*=0.1; p.vy*=0.1;
                double ang=frand()*M_PI*2.0;
                Particle e; e.x=p.x; e.y=p.y;
                e.vx=std::cos(ang)*8.0; e.vy=std::sin(ang)*8.0;
                e.type=ELECTRON; e.radius=4.0; e.flash=0;
                spawned.push_back(e);
            }
            if(p.flash>0) p.flash--;
            p.x+=p.vx; p.y+=p.vy;
            p.vx*=g_friction; p.vy*=g_friction;
            if(p.x<p.radius||p.x>FW-p.radius) p.vx*=-1;
            if(p.y<p.radius||p.y>FH-p.radius) p.vy*=-1;
            // pairwise proximity collision: push apart + swap velocities
            for(size_t j=0;j<n;++j){
                if(j==i) continue;
                Particle &o=parts[j];
                double dx=o.x-p.x, dy=o.y-p.y;
                double dist=std::sqrt(dx*dx+dy*dy);
                double minDist=p.radius+o.radius;
                if(dist<minDist){
                    double overlap=minDist-dist;
                    p.x-=(dx/(dist+1.0))*overlap*0.5;
                    p.y-=(dy/(dist+1.0))*overlap*0.5;
                    double tvx=p.vx, tvy=p.vy;
                    p.vx=o.vx*0.9; p.vy=o.vy*0.9;
                    o.vx=tvx*0.9;  o.vy=tvy*0.9;
                }
            }
        }
        for(auto &e: spawned) parts.push_back(e);
        int protons=0; for(auto &p: parts) if(p.type==PROTON) protons++;
        if(protons>40) initParticles();
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(5,5,10,1.f));
    for(int x=0;x<FW;x+=40) olivec_rect(oc, x,0, 1,FH, rgba(0,17,34,1.f));
    for(int y=0;y<FH;y+=40) olivec_rect(oc, 0,y, FW,1, rgba(0,17,34,1.f));
    for(auto &p: parts){
        int r,g,b;
        if(p.type==NEUTRON){ r=51;g=204;b=102; }
        else if(p.type==PROTON){ r=255;g=153;b=0; }
        else { r=0;g=255;b=255; }
        int rad=(int)p.radius;
        if(p.flash>0){
            // white flash with strong glow at the moment of decay
            olivec_circle(oc,(int)p.x,(int)p.y,rad+8, rgba(255,255,255,0.35f));
            olivec_circle(oc,(int)p.x,(int)p.y,rad,   rgba(255,255,255,1.0f));
        } else {
            // glow halo + solid core
            olivec_circle(oc,(int)p.x,(int)p.y,rad+ (p.type==ELECTRON?6:3), rgba(r,g,b,0.30f));
            olivec_circle(oc,(int)p.x,(int)p.y,rad, rgba(r,g,b,1.0f));
        }
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
    sim_click(0.5,0.5);   // inject a shockwave to trigger decays
    sim_click(0.3,0.6);
    int steps = argc>1? atoi(argv[1]) : 120;
    sim_step(steps);
    uint8_t* p=sim_render();
    int protons=0,electrons=0; for(auto &q:parts){ if(q.type==PROTON)protons++; if(q.type==ELECTRON)electrons++; }
    printf("weak_os native: %dx%d after %d steps  parts=%d protons=%d electrons=%d\n",
           FW,FH,steps,(int)parts.size(),protons,electrons);
    stbi_write_png("weak_os_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
