// Universe OS — Open superconductor system / Cooper-pair transport (C++/WASM)
// Electrons are injected on the left and pushed right by a voltage. A phonon
// field (2D wave grid) is dented under each electron; the resulting attraction,
// plus Pauli repulsion and a velocity-match test, lets slow same-speed electrons
// form Cooper pairs (red, near-frictionless) that stream across as a supercurrent.
// Raising temperature adds thermal noise that breaks pairs (ohmic scattering).
// State + physics + rendering (olive.c mesh + particles) all in C++.
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

static const int CELL=15, COLS=53, ROWS=35;
static const int FW=COLS*CELL, FH=ROWS*CELL;   // 795 x 525
static const int TARGET=30;

static std::vector<double> curX, prvX, nxtX;
static std::vector<uint32_t> px;

struct Elec{ int id; double x,y,vx,vy; bool paired; int partner; };
static std::vector<Elec> es;
static int idc=0;
static double g_temp=0.20, g_volt=0.50;
static long total_time=0;

static uint32_t rng=123456789u;
static inline double rnd(){ rng^=rng<<13; rng^=rng>>17; rng^=rng<<5; return (rng&0xFFFFFF)/(double)0x1000000; }
static inline int IDX(int i,int j){ return i*ROWS+j; }
static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}
static int findById(int id){ for(size_t k=0;k<es.size();++k) if(es[k].id==id) return (int)k; return -1; }

extern "C" {
KEEP int sim_w(){ return FW; }
KEEP int sim_h(){ return FH; }

KEEP void sim_reset(){
    curX.assign((size_t)COLS*ROWS,0.0); prvX.assign((size_t)COLS*ROWS,0.0); nxtX.assign((size_t)COLS*ROWS,0.0);
    es.clear(); idc=0; total_time=0;
    for(int i=0;i<15;++i) es.push_back({idc++, rnd()*(FW-100)+50, rnd()*FH, 0,0, false, -1});
}
KEEP int sim_init(int,int){ px.assign((size_t)FW*FH,0); sim_reset(); return 1; }

KEEP void sim_set(int id,double v){ if(id==0) g_temp=v; else if(id==1) g_volt=v; }
KEEP void sim_action(int){}
KEEP void sim_click(double,double){}

KEEP void sim_step(int frames){
    for(int f=0;f<frames;++f){
        // inject from the left
        if(((int)es.size()<TARGET || rnd()<0.15) && (int)es.size()<45)
            es.push_back({idc++, 5.0+rnd()*15.0, rnd()*FH, g_volt*2.0, (rnd()-0.5)*2.0, false, -1});

        // phonon wave field
        const double OMEGA=0.15;
        for(int i=1;i<COLS-1;++i) for(int j=1;j<ROWS-1;++j){
            double lap=curX[IDX(i+1,j)]+curX[IDX(i-1,j)]+curX[IDX(i,j+1)]+curX[IDX(i,j-1)]-4.0*curX[IDX(i,j)];
            double nv=(2.0*curX[IDX(i,j)]-prvX[IDX(i,j)])+OMEGA*lap;
            nv *= (i<5||i>COLS-6)?0.75:0.92;
            nxtX[IDX(i,j)]=nv;
        }
        for(auto&e:es) e.paired=false;

        for(int i=(int)es.size()-1;i>=0;--i){
            Elec&e1=es[i];
            if(e1.x>FW-10){ if(e1.partner>=0){ int p=findById(e1.partner); if(p>=0) es[p].partner=-1; }
                es.erase(es.begin()+i); continue; }
            double fx=g_volt, fy=0.0;
            for(int j=0;j<(int)es.size();++j){
                if(i==j) continue; Elec&e2=es[j];
                double dx=e1.x-e2.x, dy=e1.y-e2.y, dist=std::sqrt(dx*dx+dy*dy); if(dist<5)dist=5;
                if(e1.partner==e2.id){
                    e1.paired=true;
                    double pf=(dist-25.0)*0.06; fx-=(dx/dist)*pf; fy-=(dy/dist)*pf;
                    if(g_temp>0.35||dist>65){ e1.partner=-1; e2.partner=-1; }
                } else {
                    double rep=350.0/(dist*dist); if(dist<32) rep+=(32.0-dist)*12.0;
                    fx+=(dx/dist)*rep; fy+=(dy/dist)*rep;
                    if(e1.partner<0&&e2.partner<0&&e1.x>60&&e2.x>60&&g_temp<0.35&&dist>32&&dist<55){
                        if(std::fabs(e1.vx-e2.vx)<1.2&&std::fabs(e1.vy-e2.vy)<1.2){ e1.partner=e2.id; e2.partner=e1.id; }
                    }
                }
            }
            int gi=(int)(e1.x/CELL), gj=(int)(e1.y/CELL);
            if(gi>0&&gi<COLS-1&&gj>0&&gj<ROWS-1){
                double gX=curX[IDX(gi+1,gj)]-curX[IDX(gi-1,gj)], gY=curX[IDX(gi,gj+1)]-curX[IDX(gi,gj-1)];
                fx+=gX*2.5; fy+=gY*2.5;
                nxtX[IDX(gi,gj)]+=0.7; if(nxtX[IDX(gi,gj)]>4.0) nxtX[IDX(gi,gj)]=4.0;
            }
            fx+=(rnd()-0.5)*g_temp*10.0; fy+=(rnd()-0.5)*g_temp*12.0;
            e1.vx+=fx*0.1; e1.vy+=fy*0.1;
            double fr = e1.paired?0.996:(0.85-g_temp*0.1); e1.vx*=fr; e1.vy*=fr;
            if(e1.y<10||e1.y>FH-10) e1.vy*=-1;
            double sp=std::sqrt(e1.vx*e1.vx+e1.vy*e1.vy), mx=e1.paired?4.5:2.0;
            if(sp>mx){ e1.vx=e1.vx/sp*mx; e1.vy=e1.vy/sp*mx; }
            e1.x+=e1.vx; e1.y+=e1.vy;
        }
        std::swap(prvX,curX); std::swap(curX,nxtX);
        total_time++;
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc=olivec_canvas(px.data(),FW,FH,FW);
    olivec_fill(oc, rgba(2,2,6,1.f));
    // phonon lattice mesh (displaced by the field)
    uint32_t mesh=rgba(0,255,200,0.10f);
    for(int i=1;i<COLS;++i) for(int j=1;j<ROWS;++j){
        double val=curX[IDX(i,j)]; int pxx=i*CELL, pyy=(int)(j*CELL+val*2);
        if(val>1.0){ float a=(float)std::min(val*0.1,0.3); olivec_rect(oc,pxx-CELL/2,pyy-CELL/2,CELL,CELL, rgba(0,150,255,a)); }
        if(i<COLS-1) olivec_line(oc,pxx,pyy,(i+1)*CELL,(int)(j*CELL+curX[IDX(i+1,j)]*2), mesh);
        if(j<ROWS-1) olivec_line(oc,pxx,pyy,i*CELL,(int)((j+1)*CELL+curX[IDX(i,j+1)]*2), mesh);
    }
    // Cooper-pair bonds
    for(auto&e:es) if(e.partner>=0 && e.id<e.partner){ int p=findById(e.partner);
        if(p>=0) olivec_line(oc,(int)e.x,(int)e.y,(int)es[p].x,(int)es[p].y, rgba(255,51,102,0.95f)); }
    // electrons
    for(auto&e:es) olivec_circle(oc,(int)e.x,(int)e.y,4, e.paired?rgba(255,51,102,1.f):rgba(0,255,255,1.f));
    int pairs=0; for(auto&e:es) if(e.paired) pairs++;
    double res = es.size()? (1.0-(double)pairs/es.size())*100.0 : 100.0;
    char buf[96]; snprintf(buf,sizeof(buf),"ELECTRONS %d  RESISTANCE %.0f%%", (int)es.size(), res);
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
    int pairs=0; for(auto&e:es) if(e.paired) pairs++;
    printf("superconductor_os native: %dx%d after %d steps, electrons=%zu pairs=%d\n",FW,FH,steps,es.size(),pairs);
    stbi_write_png("superconductor_os_preview.png",FW,FH,4,p,FW*4);
    return 0;
}
#endif
