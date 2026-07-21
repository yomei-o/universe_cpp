// Universe OS — Electron Orbital / Electron Cloud Kernel  (C++/WASM port)
// NUM_ELECTRONS point-electrons feel a Coulomb-impedance attraction from 1-2
// nuclei (F = charge/(C_BASE*dist*0.01)^2), a short-range electron-electron
// repulsion, and near-nucleus quantization jitter. Positions accumulate an 8
// point trail rendered as a faint cloud. Clicking drops a 2nd nucleus and the
// s-orbital deforms into a covalent "peanut" molecular orbital. All state,
// update AND rendering (olive.c into an RGBA framebuffer) live in C++; JS only
// blits the buffer and forwards UI events.
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
#include <algorithm>
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
static const int NUM_ELECTRONS=120;
static const double C_BASE=137.0;     // your constant C=137
static const double NUCLEUS_CHARGE=40.0;

// tunable params
static double g_damp=0.95;    // velocity damping (thermal equilibration)
static double g_noise=1.0;    // quantization-noise multiplier
static double g_repel=0.5;    // electron-electron repulsion strength

struct Electron { double x,y,vx,vy; double hx[8],hy[8]; int hlen; };
struct Nucleus  { double x,y,charge; };

static std::vector<Electron> els;
static std::vector<Nucleus>  nuclei;
static std::vector<uint32_t> px;

// simple xorshift PRNG in [0,1)
static uint32_t rng_state=0x1234567u;
static inline double frand(){
    rng_state^=rng_state<<13; rng_state^=rng_state>>17; rng_state^=rng_state<<5;
    return (rng_state&0xFFFFFF)/(double)0x1000000;
}

// pack an olive color (0xAABBGGRR in memory = RGBA bytes) from rgb + alpha[0..1]
static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.0f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}

extern "C" {

KEEP int  sim_w(){ return FW; }
KEEP int  sim_h(){ return FH; }

KEEP void sim_reset(){
    rng_state=0x1234567u;
    nuclei.clear();
    nuclei.push_back({FW/2.0, FH/2.0, NUCLEUS_CHARGE});
    els.clear(); els.reserve(NUM_ELECTRONS);
    for(int i=0;i<NUM_ELECTRONS;++i){
        double ang=frand()*M_PI*2.0;
        double r=frand()*100.0+20.0;
        Electron e;
        e.x=FW/2.0+std::cos(ang)*r;
        e.y=FH/2.0+std::sin(ang)*r;
        e.vx=(frand()-0.5)*5.0;
        e.vy=(frand()-0.5)*5.0;
        e.hlen=0;
        els.push_back(e);
    }
    // start the persistent trail buffer at opaque background
    if(!px.empty()) std::fill(px.begin(), px.end(), rgba(5,5,10,1.f));
}

KEEP int sim_init(int,int){
    px.assign((size_t)FW*FH,0);
    sim_reset();
    return 1;
}

KEEP void sim_set(int id,double v){
    if(id==0) g_damp=v;
    else if(id==1) g_noise=v;
    else if(id==2) g_repel=v;
}
KEEP void sim_action(int){}

// click drops a 2nd nucleus (or overwrites it) -> covalent molecular orbital
KEEP void sim_click(double nx,double ny){
    double mx=nx*FW, my=ny*FH;
    if(nuclei.size()<2) nuclei.push_back({mx,my,NUCLEUS_CHARGE});
    else nuclei[1]={mx,my,NUCLEUS_CHARGE};
}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        for(auto &e: els){
            double fx=0, fy=0;
            // layer-2 attraction from every nucleus: F = charge/(C*dist*0.01)^2
            for(auto &n: nuclei){
                double dx=n.x-e.x, dy=n.y-e.y;
                double dist=std::sqrt(dx*dx+dy*dy);
                if(dist<5) dist=5;                       // Planck-distance guard
                double denom=C_BASE*(dist*0.01);
                double force=n.charge/(denom*denom);
                fx+=(dx/dist)*force;
                fy+=(dy/dist)*force;
            }
            // electron-electron exclusion (short-range repulsion)
            for(auto &o: els){
                if(&o==&e) continue;
                double dx=o.x-e.x, dy=o.y-e.y;
                double dist=std::sqrt(dx*dx+dy*dy);
                if(dist<30){
                    fx-=(dx/(dist+1.0))*(g_repel/(dist+1.0));
                    fy-=(dy/(dist+1.0))*(g_repel/(dist+1.0));
                }
            }
            e.vx+=fx; e.vy+=fy;
            // near-nucleus quantization jitter (sampling fold-back noise)
            for(auto &n: nuclei){
                double dx=n.x-e.x, dy=n.y-e.y;
                double dist=std::sqrt(dx*dx+dy*dy);
                if(dist<80){
                    double ns=(80.0-dist)*0.15*g_noise;
                    e.vx+=(frand()-0.5)*ns;
                    e.vy+=(frand()-0.5)*ns;
                }
            }
            e.vx*=g_damp; e.vy*=g_damp;
            e.x+=e.vx; e.y+=e.vy;
            if(e.x<10||e.x>FW-10) e.vx*=-1;
            if(e.y<10||e.y>FH-10) e.vy*=-1;
            // record trail (max 8)
            if(e.hlen<8){ e.hx[e.hlen]=e.x; e.hy[e.hlen]=e.y; e.hlen++; }
            else { for(int i=0;i<7;++i){e.hx[i]=e.hx[i+1]; e.hy[i]=e.hy[i+1];} e.hx[7]=e.x; e.hy[7]=e.y; }
        }
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    // fade previous frame slightly to leave a cloud trail (rgba(5,5,10,0.2))
    olivec_rect(oc, 0,0, FW,FH, rgba(5,5,10,0.2f));
    // faint spatial grid
    for(int x=0;x<FW;x+=50) olivec_rect(oc, x,0, 1,FH, rgba(0,17,8,1.f));
    for(int y=0;y<FH;y+=50) olivec_rect(oc, 0,y, FW,1, rgba(0,17,8,1.f));
    // electrons: trail polyline + core dot
    for(auto &e: els){
        for(int i=1;i<e.hlen;++i)
            olivec_line(oc,(int)e.hx[i-1],(int)e.hy[i-1],(int)e.hx[i],(int)e.hy[i], rgba(0,255,200,0.14f));
        olivec_circle(oc,(int)e.x,(int)e.y,2, rgba(0,255,200,0.8f));
    }
    // nuclei: glowing positive cores
    for(auto &n: nuclei){
        olivec_circle(oc,(int)n.x,(int)n.y,14, rgba(255,170,0,0.20f));
        olivec_circle(oc,(int)n.x,(int)n.y,8,  rgba(255,170,0,1.0f));
    }
    char buf[64];
    snprintf(buf,sizeof(buf),"e-:%d  nuclei:%d", (int)els.size(), (int)nuclei.size());
    olivec_text(oc, buf, 10, 10, olivec_default_font, 2, rgba(0,255,200,0.8f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
int main(int argc,char**argv){
    sim_init(0,0);
    // drop a 2nd nucleus to exercise the covalent path
    sim_click(0.62, 0.5);
    int steps = argc>1? atoi(argv[1]) : 400;
    sim_step(steps);
    uint8_t* p=sim_render();
    printf("electron_os native: %dx%d after %d steps\n", FW,FH,steps);
    stbi_write_png("electron_os_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
