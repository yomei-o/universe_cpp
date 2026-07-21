// Universe OS — Time Dilation Kernel (Twin Paradox)  (C++/WASM port)
// The stationary Earth twin's clock advances at full rate; the ship twin flies a
// scripted round-trip (accelerate / cruise / U-turn / return / land). High
// velocity raises the coordinate-rewrite "task load", so the ship's proper-time
// step deltaTau = 1 - load is throttled (frame-skipped) — accumulating the twin
// paradox. State, update AND rendering (olive.c -> RGBA) are all in C++; JS blits.
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

static const int FW=800, FH=450;

static double earthTime, shipTime, shipX, shipVelocity;
static int shipY=280, moveDirection, missionStep;
static bool isMoving;
static char logMsg[128];
static std::vector<uint32_t> px;

static inline uint32_t rgba(int r,int g,int b,float a){
    int A=(int)(a*255.f); if(A<0)A=0; if(A>255)A=255;
    return ((uint32_t)A<<24)|((uint32_t)b<<16)|((uint32_t)g<<8)|(uint32_t)r;
}
static void thick_line(Olivec_Canvas oc,double x1,double y1,double x2,double y2,double w,uint32_t c){
    double dx=x2-x1,dy=y2-y1,len=std::sqrt(dx*dx+dy*dy);
    if(len<1e-6)return; double nx=-dy/len,ny=dx/len; int half=(int)(w/2); if(half<1)half=1;
    for(int o=-half;o<=half;o++) olivec_line(oc,(int)(x1+nx*o),(int)(y1+ny*o),(int)(x2+nx*o),(int)(y2+ny*o),c);
}
static void ring(Olivec_Canvas oc,double cx,double cy,double r,double w,uint32_t c){
    int seg=(int)(r*0.8)+16; double px_=cx+r,py_=cy;
    for(int i=1;i<=seg;i++){ double a=2*M_PI*i/seg,x=cx+std::cos(a)*r,y=cy+std::sin(a)*r;
        thick_line(oc,px_,py_,x,y,w,c); px_=x; py_=y; }
}
static void dashed_line(Olivec_Canvas oc,double x1,double y1,double x2,double y2,double dash,double gap,uint32_t c){
    double dx=x2-x1,dy=y2-y1,len=std::sqrt(dx*dx+dy*dy); if(len<1e-6)return;
    double ux=dx/len,uy=dy/len,pos=0; bool on=true;
    while(pos<len){ double e=std::min(len,pos+(on?dash:gap));
        if(on) olivec_line(oc,(int)(x1+ux*pos),(int)(y1+uy*pos),(int)(x1+ux*e),(int)(y1+uy*e),c);
        pos=e; on=!on; }
}

extern "C" {
KEEP int sim_w(){ return FW; }
KEEP int sim_h(){ return FH; }

KEEP void sim_reset(){
    earthTime=0; shipTime=0; shipX=150; shipVelocity=0;
    isMoving=false; moveDirection=1; missionStep=0;
    snprintf(logMsg,sizeof logMsg,"SYSTEM STATUS. SYSTEM IDLE. CLOCKS SYNCHRONISED.");
}
KEEP int sim_init(int,int){ px.assign((size_t)FW*FH,0); sim_reset(); return 1; }

KEEP void sim_set(int,double){}
KEEP void sim_action(int id){
    if(id==0 && !isMoving){
        isMoving=true; shipVelocity=0; moveDirection=1; missionStep=0;
        snprintf(logMsg,sizeof logMsg,"MISSION LOG. SHIP LAUNCHED. HIGH VELOCITY REWRITE ACTIVE.");
    }
}
KEEP void sim_click(double,double){}

KEEP void sim_step(int steps){
    for(int s=0;s<steps;++s){
        earthTime+=1.0;
        double deltaTau=1.0;
        if(isMoving){
            missionStep++;
            if(missionStep<120) shipVelocity=std::min(6.5, shipVelocity+0.1);
            else if(missionStep<180){ shipVelocity=std::max(0.2, shipVelocity-0.2); if(shipVelocity<=0.2) moveDirection=-1; }
            else if(missionStep<300) shipVelocity=std::min(6.5, shipVelocity+0.1);
            else { shipVelocity=std::max(0.0, shipVelocity-0.15);
                   if(shipVelocity==0.0){ isMoving=false;
                       snprintf(logMsg,sizeof logMsg,"MISSION COMPLETE. EARTH %d STEPS. SHIP %d STEPS. PARADOX SOLVED.",
                                (int)std::floor(earthTime),(int)std::floor(shipTime)); } }
            shipX+=shipVelocity*moveDirection;
            double taskLoad=(shipVelocity/7.0)*0.85;
            deltaTau=1.0-taskLoad;
            if(isMoving)
                snprintf(logMsg,sizeof logMsg,"SHIP VELOCITY %.0f KM/S. TASK LOAD %.0f PCT. CLOCK RATE %.2f",
                         shipVelocity*40000.0, taskLoad*100.0, deltaTau);
        }
        shipTime+=deltaTau;
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc=olivec_canvas(px.data(),FW,FH,FW);
    olivec_fill(oc, rgba(5,5,10,1.f));                        // #05050a
    for(int x=0;x<FW;x+=50) olivec_rect(oc,x,0,1,FH, rgba(0,21,10,1.f)); // #00150a

    // guideline + turn point
    dashed_line(oc,150,280,750,280,5,5, rgba(255,255,255,0.1f));
    olivec_frame(oc,750,270,5,20,1, rgba(255,51,51,1.f));
    olivec_text(oc,"TURN POINT",725,300,olivec_default_font,1, rgba(255,51,51,1.f));

    // Earth clock
    int ex=150,ey=140; uint32_t ec=rgba(0,170,255,1.f);
    ring(oc,ex,ey,45,3,ec);
    double ea=std::fmod(earthTime*0.05,2*M_PI);
    thick_line(oc,ex,ey, ex+std::cos(ea)*35, ey+std::sin(ea)*35, 4, ec);
    char b[80];
    snprintf(b,sizeof b,"EARTH CLOCK %d",(int)std::floor(earthTime));
    olivec_text(oc,b,ex+60,ey-8,olivec_default_font,2,ec);
    olivec_text(oc,"STATIONARY LOAD 0 PCT",ex+60,ey+10,olivec_default_font,1,ec);

    // Earth planet
    olivec_circle(oc,150,280,15, rgba(0,85,255,1.f));
    olivec_text(oc,"EARTH",135,308,olivec_default_font,1, rgba(255,255,255,1.f));

    // Ship clock (follows ship)
    int sx=(int)shipX, sy=shipY-80; uint32_t sc=isMoving?rgba(0,255,102,1.f):rgba(136,136,136,1.f);
    ring(oc,sx,sy,45,3,sc);
    double sa=std::fmod(shipTime*0.05,2*M_PI);
    thick_line(oc,sx,sy, sx+std::cos(sa)*35, sy+std::sin(sa)*35, 4, sc);
    snprintf(b,sizeof b,"SHIP CLOCK %d",(int)std::floor(shipTime));
    olivec_text(oc,b,sx-60,sy-58,olivec_default_font,2,sc);

    // Ship (triangle)
    uint32_t st=isMoving?rgba(0,255,102,1.f):rgba(102,102,102,1.f);
    if(moveDirection==1)
        olivec_triangle(oc,(int)shipX+15,shipY,(int)shipX-15,shipY-8,(int)shipX-15,shipY+8,st);
    else
        olivec_triangle(oc,(int)shipX-15,shipY,(int)shipX+15,shipY-8,(int)shipX+15,shipY+8,st);

    olivec_text(oc,logMsg,30,FH-24,olivec_default_font,2, rgba(0,255,204,1.f));
    return (uint8_t*)px.data();
}
}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
int main(int argc,char**argv){
    sim_init(0,0);
    sim_action(0);                       // launch mission
    int steps=argc>1?atoi(argv[1]):90; sim_step(steps);
    uint8_t*p=sim_render();
    printf("time_dilation_os native: %dx%d, earth=%.0f ship=%.0f x=%.0f\n",
           FW,FH,earthTime,shipTime,shipX);
    stbi_write_png("time_dilation_os_preview.png",FW,FH,4,p,FW*4);
    return 0;
}
#endif
