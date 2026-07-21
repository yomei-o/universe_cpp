// Universe OS — Three-Body Chaos (Figure-8) Kernel  (C++/WASM port)
// Three equal masses launched on the exact figure-8 choreography initial
// conditions integrate under pure Newtonian gravity (G=1) with a symplectic
// Euler solver (20 substeps of dt=0.001 per frame, eps=0). Each body keeps a
// 500-point trail. "Inject bug" perturbs the launch velocities ~10%, and the
// butterfly effect blows the coherent orbit into chaos after a few passes.
// All state, integration AND rendering (olive.c into an RGBA framebuffer) live
// in C++; JS only blits the buffer and forwards UI events.
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
#include <cstring>
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define KEEP EMSCRIPTEN_KEEPALIVE
#else
#define KEEP
#endif

static const int FW=850, FH=500;
static const double CX=FW/2.0, CY=FH/2.0;

// tunable params
static double g_scale=180.0;   // world -> pixel scale

// strict initial-value register (figure-8 choreography)
static const double IX1=-0.97000436, IY1=0.24308753;
static const double IVX1=0.46620531, IVY1=0.43236573;

struct Body { double x,y,vx,vy,mass; int r,g,b; std::vector<double> px_,py_; };
static Body bodies[3];
static bool g_bug=false;
static char logMsg[96];

static std::vector<uint32_t> px;

static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.0f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}

static void initSimulation(bool injectBug){
    bodies[0]={ IX1, IY1, IVX1, IVY1, 1.0, 0,255,204,{},{} };   // #00ffcc
    bodies[1]={-IX1,-IY1, IVX1, IVY1, 1.0, 255,51,102,{},{} };  // #ff3366
    bodies[2]={ 0,  0, -2*IVX1, -2*IVY1, 1.0, 255,170,0,{},{} }; // #ffaa00
    if(injectBug){
        bodies[0].vx+=0.045;
        bodies[1].vy-=0.020;
        std::strcpy(logMsg,"ERROR INJECTED. Butterfly effect growing exponentially...");
    } else {
        std::strcpy(logMsg,"PERFECT COHERENCE. Running infinite stable choreography.");
    }
    g_bug=injectBug;
}

extern "C" {

KEEP int  sim_w(){ return FW; }
KEEP int  sim_h(){ return FH; }

KEEP void sim_reset(){ initSimulation(false); }

KEEP int sim_init(int,int){
    px.assign((size_t)FW*FH,0);
    sim_reset();
    return 1;
}

KEEP void sim_set(int id,double v){
    if(id==0) g_scale=v;
}
// action 0 = figure-8 (no bug), action 1 = inject chaos bug
KEEP void sim_action(int id){
    if(id==1) initSimulation(true);
    else initSimulation(false);
}
KEEP void sim_click(double,double){}

KEEP void sim_step(int steps){
    const double G=1.0;
    const int substeps=20;
    const double dt=0.02/substeps;   // 0.001
    for(int s=0;s<steps;++s){
        for(int step=0;step<substeps;++step){
            double ax[3]={0,0,0}, ay[3]={0,0,0};
            for(int i=0;i<3;++i) for(int j=0;j<3;++j){
                if(i==j) continue;
                double dx=bodies[j].x-bodies[i].x;
                double dy=bodies[j].y-bodies[i].y;
                double distSq=dx*dx+dy*dy;
                double dist=std::sqrt(distSq);
                double force=(G*bodies[j].mass)/(distSq*dist);   // pure Newtonian, eps=0
                ax[i]+=force*dx; ay[i]+=force*dy;
            }
            for(int i=0;i<3;++i){
                Body &b=bodies[i];
                b.vx+=ax[i]*dt; b.vy+=ay[i]*dt;
                b.x+=b.vx*dt;   b.y+=b.vy*dt;
            }
        }
        for(int i=0;i<3;++i){
            Body &b=bodies[i];
            b.px_.push_back(b.x); b.py_.push_back(b.y);
            if((int)b.px_.size()>500){ b.px_.erase(b.px_.begin()); b.py_.erase(b.py_.begin()); }
        }
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(2,2,6,1.f));                       // backdrop #020206
    for(int x=0;x<FW;x+=40) olivec_rect(oc, x,0, 1,FH, rgba(5,17,21,1.f)); // grid #051115
    // trails as alpha polylines
    for(int i=0;i<3;++i){
        Body &b=bodies[i];
        int n=(int)b.px_.size();
        for(int p=1;p<n;++p){
            int x0=(int)(CX+b.px_[p-1]*g_scale), y0=(int)(CY+b.py_[p-1]*g_scale);
            int x1p=(int)(CX+b.px_[p]*g_scale),  y1p=(int)(CY+b.py_[p]*g_scale);
            olivec_line(oc,x0,y0,x1p,y1p, rgba(b.r,b.g,b.b,0.35f));
        }
    }
    // body cores (glow halo -> color -> white center)
    for(int i=0;i<3;++i){
        Body &b=bodies[i];
        int sx=(int)(CX+b.x*g_scale), sy=(int)(CY+b.y*g_scale);
        olivec_circle(oc,sx,sy,10, rgba(b.r,b.g,b.b,0.30f));
        olivec_circle(oc,sx,sy,6,  rgba(b.r,b.g,b.b,0.9f));
        olivec_circle(oc,sx,sy,3,  rgba(255,255,255,1.0f));
    }
    olivec_text(oc, "G = 1.0", 20, 14, olivec_default_font, 2, rgba(0,255,200,0.6f));
    olivec_text(oc, logMsg, 20, FH-22, olivec_default_font, 2,
                g_bug? rgba(255,51,102,0.9f) : rgba(0,255,200,0.9f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
int main(int argc,char**argv){
    sim_init(0,0);
    int steps = argc>1? atoi(argv[1]) : 300;
    sim_step(steps);
    uint8_t* p=sim_render();
    printf("three_body_chaos_working native: %dx%d after %d frames  b0=(%.3f,%.3f)\n",
           FW,FH,steps,bodies[0].x,bodies[0].y);
    stbi_write_png("three_body_chaos_working_preview.png", FW,FH, 4, p, FW*4);
    return 0;
}
#endif
