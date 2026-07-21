// Universe OS — H2O Vibrational Mode Kernel  (C++/WASM port)
// A 3-atom water molecule (O + 2H) held by classical bond-stretch springs and an
// angle-bend restoring force, integrated with semi-implicit Euler. Buttons inject
// the three normal modes (v1 symmetric stretch, v2 bend, v3 asymmetric stretch);
// energy is auto re-injected when the vibration decays. The purple arrow is the
// live dipole moment vector. State, physics AND rendering (olive.c -> RGBA) are
// all in C++; JS only blits.
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
static const double CX=FW/2.0, CY=FH/2.0-30.0;
static const double R0=110.0;
static const double THETA0=104.5*M_PI/180.0;

struct Atom { double x,y,vx,vy,mass,charge; };
static Atom O,H1,H2;
static int currentMode=1;
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
// radial-gradient disk approximation: inner color -> mid color (@0.8) -> transparent
static void grad_disk(Olivec_Canvas oc,double cx,double cy,int R,
                      int r0,int g0,int b0,int r1,int g1,int b1){
    for(int rr=R;rr>=1;--rr){
        double f=(double)rr/R; int cr,cg,cb; float a;
        if(f>0.8){ double t=(f-0.8)/0.2; cr=r1; cg=g1; cb=b1; a=(float)(1.0-t); }
        else { double t=f/0.8; cr=(int)(r0+(r1-r0)*t); cg=(int)(g0+(g1-g0)*t); cb=(int)(b0+(b1-b0)*t); a=1.f; }
        olivec_circle(oc,(int)cx,(int)cy,rr, rgba(cr,cg,cb,a));
    }
}

static void resetToEquilibrium(){
    O={0,-30,0,0,16,-0.66};
    double h=THETA0/2;
    H1={ -R0*std::sin(h), -30+R0*std::cos(h),0,0,1,0.33};
    H2={  R0*std::sin(h), -30+R0*std::cos(h),0,0,1,0.33};
}
static void injectMode(int mode){
    resetToEquilibrium();
    double h=THETA0/2, s=std::sin(h), c=std::cos(h);
    if(mode==1){ double p=2.5;
        H1.vx=s*p; H1.vy=-c*p; H2.vx=-s*p; H2.vy=-c*p;
        O.vx=-(H1.vx+H2.vx)/O.mass; O.vy=-(H1.vy+H2.vy)/O.mass;
    } else if(mode==2){ double p=3.0;
        H1.vx=c*p; H1.vy=s*p; H2.vx=-c*p; H2.vy=s*p;
        O.vy=-(H1.vy+H2.vy)/O.mass;
    } else if(mode==3){ double p=3.0;
        H1.vx=s*p; H1.vy=-c*p; H2.vx=s*p; H2.vy=-c*p;
        O.vx=-(H1.vx+H2.vx)/O.mass; O.vy=-(H1.vy+H2.vy)/O.mass;
    }
}

extern "C" {
KEEP int sim_w(){ return FW; }
KEEP int sim_h(){ return FH; }

KEEP void sim_reset(){ currentMode=1; injectMode(1); }
KEEP int sim_init(int,int){ px.assign((size_t)FW*FH,0); sim_reset(); return 1; }

KEEP void sim_set(int,double){}
KEEP void sim_action(int id){ if(id>=0&&id<=2){ currentMode=id+1; injectMode(currentMode); } }
KEEP void sim_click(double,double){}

KEEP void sim_step(int steps){
    const double kr=0.12, ktheta=45.0, dt=0.5;
    for(int s=0;s<steps;++s){
        double dx1=H1.x-O.x, dy1=H1.y-O.y, r1=std::sqrt(dx1*dx1+dy1*dy1);
        double f_r1=-kr*(r1-R0); double fx_r1=(dx1/r1)*f_r1, fy_r1=(dy1/r1)*f_r1;
        double dx2=H2.x-O.x, dy2=H2.y-O.y, r2=std::sqrt(dx2*dx2+dy2*dy2);
        double f_r2=-kr*(r2-R0); double fx_r2=(dx2/r2)*f_r2, fy_r2=(dy2/r2)*f_r2;

        double cosang=(dx1*dx2+dy1*dy2)/(r1*r2);
        if(cosang>1)cosang=1; if(cosang<-1)cosang=-1;
        double theta=std::acos(cosang);
        double f_theta=-ktheta*(theta-THETA0);
        double fx_t1=(-dy1/r1)*f_theta/r1, fy_t1=(dx1/r1)*f_theta/r1;
        double fx_t2=(dy2/r2)*f_theta/r2, fy_t2=(-dx2/r2)*f_theta/r2;

        H1.vx+=(fx_r1+fx_t1)/H1.mass*dt; H1.vy+=(fy_r1+fy_t1)/H1.mass*dt;
        H2.vx+=(fx_r2+fx_t2)/H2.mass*dt; H2.vy+=(fy_r2+fy_t2)/H2.mass*dt;
        O.vx+=-(fx_r1+fx_r2+fx_t1+fx_t2)/O.mass*dt;
        O.vy+=-(fy_r1+fy_r2+fy_t1+fy_t2)/O.mass*dt;

        for(Atom*p:{&O,&H1,&H2}){ p->vx*=0.995; p->vy*=0.995; p->x+=p->vx*dt; p->y+=p->vy*dt; }

        double velSum=std::fabs(H1.vx)+std::fabs(H1.vy)+std::fabs(H2.vx)+std::fabs(H2.vy);
        if(velSum<0.15) injectMode(currentMode);
    }
}

KEEP uint8_t* sim_render(){
    Olivec_Canvas oc=olivec_canvas(px.data(),FW,FH,FW);
    olivec_fill(oc, rgba(3,3,12,1.f));                         // #03030c
    for(int x=0;x<FW;x+=40) olivec_rect(oc,x,0,1,FH, rgba(5,17,17,1.f)); // #051111

    // bonds (world coords)
    uint32_t bond=rgba(85,102,136,1.f);
    thick_line(oc, CX+O.x,CY+O.y, CX+H1.x,CY+H1.y, 6, bond);
    thick_line(oc, CX+O.x,CY+O.y, CX+H2.x,CY+H2.y, 6, bond);

    // dipole vector arrow
    double dipX=(O.x*O.charge+H1.x*H1.charge+H2.x*H2.charge)*4.0;
    double dipY=(O.y*O.charge+H1.y*H1.charge+H2.y*H2.charge)*4.0;
    double ax=CX+O.x, ay=CY+O.y+10, bx=CX+O.x+dipX, by=CY+O.y+10+dipY;
    uint32_t dip=rgba(204,51,255,1.f);
    thick_line(oc,ax,ay,bx,by,4,dip);
    { double dx=bx-ax,dy=by-ay,l=std::sqrt(dx*dx+dy*dy);
      if(l>1){ double ux=dx/l,uy=dy/l; double h=10;
        thick_line(oc,bx,by, bx-ux*h-uy*h*0.5, by-uy*h+ux*h*0.5, 3, dip);
        thick_line(oc,bx,by, bx-ux*h+uy*h*0.5, by-uy*h-ux*h*0.5, 3, dip); } }

    // atoms
    grad_disk(oc, CX+O.x, CY+O.y, 24, 255,85,85, 170,34,34);   // O red
    olivec_text(oc,"O",(int)(CX+O.x)-3,(int)(CY+O.y)-3,olivec_default_font,1, rgba(255,255,255,1.f));
    grad_disk(oc, CX+H1.x, CY+H1.y, 14, 85,170,255, 34,102,170); // H blue
    olivec_text(oc,"H",(int)(CX+H1.x)-3,(int)(CY+H1.y)-3,olivec_default_font,1, rgba(255,255,255,1.f));
    grad_disk(oc, CX+H2.x, CY+H2.y, 14, 85,170,255, 34,102,170);
    olivec_text(oc,"H",(int)(CX+H2.x)-3,(int)(CY+H2.y)-3,olivec_default_font,1, rgba(255,255,255,1.f));

    // status readouts
    double dx1=H1.x-O.x,dy1=H1.y-O.y,r1=std::sqrt(dx1*dx1+dy1*dy1);
    double dx2=H2.x-O.x,dy2=H2.y-O.y,r2=std::sqrt(dx2*dx2+dy2*dy2);
    double cosang=(dx1*dx2+dy1*dy2)/(r1*r2); if(cosang>1)cosang=1; if(cosang<-1)cosang=-1;
    double deg=std::acos(cosang)*180.0/M_PI;
    double dipMag=std::sqrt(dipX*dipX+dipY*dipY);
    char b[96]; uint32_t fg=rgba(0,255,204,1.f);
    snprintf(b,sizeof b,"BOND LENGTH R1  %.1f PM",r1);           olivec_text(oc,b,30,FH-100,olivec_default_font,2,fg);
    snprintf(b,sizeof b,"BOND ANGLE      %.1f DEG (REF 104.5)",deg); olivec_text(oc,b,30,FH-78,olivec_default_font,2,fg);
    snprintf(b,sizeof b,"DIPOLE MOMENT   %.1f",dipMag);          olivec_text(oc,b,30,FH-56,olivec_default_font,2,fg);
    const char* mn = currentMode==1?"MODE V1 SYMMETRIC STRETCH":currentMode==2?"MODE V2 BEND":"MODE V3 ASYMMETRIC STRETCH";
    olivec_text(oc,mn,30,20,olivec_default_font,2, rgba(0,170,255,1.f));
    return (uint8_t*)px.data();
}
}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
int main(int argc,char**argv){
    sim_init(0,0);
    int steps=argc>1?atoi(argv[1]):20; sim_step(steps);
    uint8_t*p=sim_render();
    printf("water_vibration_os native: %dx%d, mode=%d\n",FW,FH,currentMode);
    stbi_write_png("water_vibration_os_preview.png",FW,FH,4,p,FW*4);
    return 0;
}
#endif
