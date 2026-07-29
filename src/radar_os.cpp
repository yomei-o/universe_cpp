// Universe OS — Radar: Two Charts, One Answer  (C++/WASM)
//
// A numerical relativity code and a "c*t = const" calculation, running side by
// side on the same physical experiment: bounce a light pulse off a mirror near
// a black hole and time the round trip at the dish.
//
// Both panels integrate the SAME 3+1 formula for the coordinate light speed,
//
//     dx^i / dt  =  -beta^i  +-  alpha / sqrt(gamma_ii)
//
// and differ only in what (alpha, beta, gamma) they are handed:
//
//   pg   (numerical relativity's usual choice: horizon penetrating)
//        alpha = 1,  beta^r = sqrt(2m/r),  gamma_rr = 1
//        ->  dr/dt = +-1 - sqrt(2m/r)         space flat, free fallers slide
//                                             through a grid that stays put
//
//   iso  (the c*t = const chart: flat grid, light speed varies)
//        alpha = (1-u)/(1+u),  beta = 0,  gamma_rr = (1+u)^4,  u = m/2rho
//        ->  drho/dt = +- (1-u)/(1+u)^3       grid still, c varies
//
// The grid never falls in either of them -- a grid point is a coordinate label.
// alpha says whether the observers normal to the slices are in free fall; beta
// says whether the grid goes along with them. In pg, alpha = 1 means they are
// falling (n^mu = (1, -sqrt(2m/r), 0, 0), rain from rest at infinity) and the
// shift is there precisely so the grid does NOT follow: the radial label stays
// the areal radius. Set beta = 0 instead and the grid really does ride the
// free fallers -- that is geodesic slicing, and from the momentarily static
// slice the throat reaches r = 0 at proper time exactly pi*m = 3.141592654,
// which is the entire lifetime of such a run.
//
// The point of the demo is that the agreement is exactly as wide as what you
// can measure and not one bit wider:
//
//   when the pulse turned around   pg 35.867448   iso 54.742402   DIFFER by 19
//   when the pulse got back        agree to 2e-12
//
// "Where is the pulse now" is a chart's private bookkeeping, because "now" is a
// choice of slicing. The round trip is one clock reading minus another clock
// reading at the same place, so it cannot depend on the choice -- and it does
// not, to twelve digits, from two integrals that share no integrand:
//
//   pg    integral of 2 dr / (1 - 2m/r)          over the areal radius
//   iso   integral of 2 drho / c(rho)            over the isotropic radius
//   both  = 2[(r1-r2) + 4m ln((r1-2m)/(r2-2m))] = 109.484804043632 for 3m..50m
#define OLIVEC_IMPLEMENTATION
#include "olive.c"
#include <vector>
#include <cstdint>
#include <cmath>
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

static const int FW = 960, FH = 620;
static std::vector<uint32_t> px;

static inline uint32_t rgba(int r, int g, int b, float a) {
    int A = (int)(a * 255.0f); if (A < 0) A = 0; if (A > 255) A = 255;
    return ((uint32_t)A << 24) | ((uint32_t)b << 16) | ((uint32_t)g << 8) | (uint32_t)r;
}

// ----------------------------------------------- chart 1: Painleve-Gullstrand
// what a numerical relativity code actually uses: flat spatial metric, unit
// lapse, all the gravity in the shift. horizon penetrating.
static double pg_alpha(double) { return 1.0; }
static double pg_beta(double r) { return std::sqrt(2.0 / r); }
static double pg_gam(double) { return 1.0; }
static double pg_areal(double r) { return r; }
static double pg_x_of_areal(double r) { return r; }

// ------------------------------------------------- chart 2: isotropic, c*t=const
// flat grid, no shift at all, and the whole of gravity sitting in a light speed
// that depends on where you are.
static inline double iso_u(double p) { return 0.5 / p; }
static double iso_alpha(double p) { double u = iso_u(p); return (1.0 - u) / (1.0 + u); }
static double iso_beta(double) { return 0.0; }
static double iso_gam(double p) { double B = 1.0 + iso_u(p); return B * B * B * B; }
static double iso_areal(double p) { double B = 1.0 + iso_u(p); return p * B * B; }
static double iso_x_of_areal(double r) {
    double a = r - 1.0;                                  // rho + m + m^2/4rho = r
    double d = a * a - 1.0; if (d < 0.0) d = 0.0;
    return 0.5 * (a + std::sqrt(d));
}

// the one formula both charts are fed
struct Chart {
    double (*alpha)(double);
    double (*beta)(double);
    double (*gam)(double);
    double (*areal)(double);
    double (*x_of_areal)(double);
    // grid state
    int N = 0;
    std::vector<double> x, uin, uout, a1, a2, b1, b2, s1, s2;
    double dx = 0, dt = 0, t = 0;
    int idet = 0;
    double tarr = -1.0, pk = 0.0;                        // measured return time
    double p0 = 0, p1 = 0, p2 = 0;                       // last 3 dish samples
    double tp0 = 0, tp1 = 0, tp2 = 0;
    std::vector<double> rin, rout, rt;                        // dish record
    double xpk_in = 0, xpk_out = 0;                      // pulse peaks
    double amp_in = 0, amp_out = 0;
    double quad = 0, half_in = 0;                        // exact quadrature results
};
static Chart PG, ISO;

static inline double lam_out(const Chart& c, double x) {
    return -c.beta(x) + c.alpha(x) / std::sqrt(c.gam(x));
}
static inline double lam_in(const Chart& c, double x) {
    return -c.beta(x) - c.alpha(x) / std::sqrt(c.gam(x));
}

// ---------------------------------------------------------------- quadrature
// dt = dx / |coordinate light speed|, integrated in each chart's own variable
static double simpson_leg(const Chart& c, double xa, double xb, int n, int which) {
    if (n % 2) ++n;
    double h = (xb - xa) / n, s = 0.0;
    for (int i = 0; i <= n; ++i) {
        double x = xa + i * h;
        double f = 0.0;
        if (which != 1) f += 1.0 / std::fabs(lam_in(c, x));
        if (which != 0) f += 1.0 / lam_out(c, x);
        double w = (i == 0 || i == n) ? 1.0 : ((i & 1) ? 4.0 : 2.0);
        s += w * f;
    }
    return s * h / 3.0;
}

static const double R_PAD = 12.0;                        // grid extends past the dish
static double g_mirror = 6.0, g_dish = 50.0, g_speed = 1.0, g_width = 2.5;
static int g_N = 1600;
static bool dirty = true;
static double T_EXACT = 0.0;
static const double REC_DT = 0.1;
static const int REC_N = 2200;

// ---------------------------------------------------------------- grid solver
static void ddx(const std::vector<double>& u, double dx, std::vector<double>& d) {
    int N = (int)u.size();
    double c = 1.0 / (12.0 * dx);
    for (int i = 2; i < N - 2; ++i)
        d[i] = (-u[i + 2] + 8.0 * u[i + 1] - 8.0 * u[i - 1] + u[i - 2]) * c;
    // 4th order one sided at the edges: the mirror sits at i = 0 and a sloppy
    // boundary there would land straight on the measurement
    d[0] = (-25.0 * u[0] + 48.0 * u[1] - 36.0 * u[2] + 16.0 * u[3] - 3.0 * u[4]) * c;
    d[1] = (-3.0 * u[0] - 10.0 * u[1] + 18.0 * u[2] - 6.0 * u[3] + u[4]) * c;
    d[N - 1] = (25.0 * u[N - 1] - 48.0 * u[N - 2] + 36.0 * u[N - 3]
                - 16.0 * u[N - 4] + 3.0 * u[N - 5]) * c;
    d[N - 2] = (3.0 * u[N - 1] + 10.0 * u[N - 2] - 18.0 * u[N - 3]
                + 6.0 * u[N - 4] - u[N - 5]) * c;
}
// mirror at the inner edge: the out-going family's in-flow boundary is fed by
// the in-going family. outer edge: nothing comes in from outside.
static void bcs(Chart& c, std::vector<double>& ui, std::vector<double>& uo) {
    uo[0] = ui[0];
    ui[c.N - 1] = 0.0;
}
static void rhs(Chart& c, std::vector<double>& ui, std::vector<double>& uo,
                std::vector<double>& di, std::vector<double>& dou) {
    bcs(c, ui, uo);
    ddx(ui, c.dx, di);
    ddx(uo, c.dx, dou);
    for (int i = 0; i < c.N; ++i) {
        di[i] = -lam_in(c, c.x[i]) * di[i];
        dou[i] = -lam_out(c, c.x[i]) * dou[i];
    }
}
static void kreiss_oliger(std::vector<double>& u, double sig) {
    int N = (int)u.size();
    static std::vector<double> t;
    t = u;
    for (int i = 2; i < N - 2; ++i)
        u[i] = t[i] - sig * (t[i - 2] - 4.0 * t[i - 1] + 6.0 * t[i] - 4.0 * t[i + 1] + t[i + 2]);
}
static void step(Chart& c) {
    const int N = c.N;
    double h = c.dt;
    static std::vector<double> yi, yo, ki, ko, li, lo;
    yi = c.uin; yo = c.uout;
    ki.assign(N, 0.0); ko.assign(N, 0.0);
    li.assign(N, 0.0); lo.assign(N, 0.0);
    // classic RK4, accumulating into (li, lo)
    rhs(c, yi, yo, ki, ko);
    for (int i = 0; i < N; ++i) { li[i] = ki[i]; lo[i] = ko[i]; }
    c.a1 = c.uin; c.a2 = c.uout;
    for (int i = 0; i < N; ++i) { c.a1[i] += 0.5 * h * ki[i]; c.a2[i] += 0.5 * h * ko[i]; }
    rhs(c, c.a1, c.a2, ki, ko);
    for (int i = 0; i < N; ++i) { li[i] += 2.0 * ki[i]; lo[i] += 2.0 * ko[i]; }
    c.a1 = c.uin; c.a2 = c.uout;
    for (int i = 0; i < N; ++i) { c.a1[i] += 0.5 * h * ki[i]; c.a2[i] += 0.5 * h * ko[i]; }
    rhs(c, c.a1, c.a2, ki, ko);
    for (int i = 0; i < N; ++i) { li[i] += 2.0 * ki[i]; lo[i] += 2.0 * ko[i]; }
    c.a1 = c.uin; c.a2 = c.uout;
    for (int i = 0; i < N; ++i) { c.a1[i] += h * ki[i]; c.a2[i] += h * ko[i]; }
    rhs(c, c.a1, c.a2, ki, ko);
    for (int i = 0; i < N; ++i) { li[i] += ki[i]; lo[i] += ko[i]; }
    for (int i = 0; i < N; ++i) {
        c.uin[i] += h / 6.0 * li[i];
        c.uout[i] += h / 6.0 * lo[i];
    }
    kreiss_oliger(c.uin, 0.01);
    kreiss_oliger(c.uout, 0.01);
    bcs(c, c.uin, c.uout);
    c.t += h;

    // dish record and return-time measurement
    int k = (int)(c.t / REC_DT);
    if (k >= 0 && k < REC_N && c.rt[k] < 0.0) {          // first sample in the bin,
        c.rt[k] = c.t;                                   // and keep its real time
        c.rin[k] = c.uin[c.idet];
        c.rout[k] = c.uout[c.idet];
    }
}
// The value of u is constant along a characteristic, so the peak of the
// returning pulse rides the characteristic that started at the dish -- its
// arrival time IS the round trip. Fit a quadratic to the top of the record
// (the record is sampled every REC_DT, coarse enough for a stable fit).
static void measure(Chart& c) {
    int km = -1; double mx = 0.0;
    for (int k = 0; k < REC_N; ++k) if (c.rout[k] > mx) { mx = c.rout[k]; km = k; }
    c.pk = mx;
    if (km <= 0 || mx < 0.05) { c.tarr = -1.0; return; }
    // three point parabola on the record. narrow on purpose: the arriving pulse
     // is not symmetric, so a wide fit buys precision at the cost of a bias that
     // does not go away when the grid is refined.
    double y0 = c.rout[km - 1], y1 = c.rout[km], y2 = c.rout[km + 1];
    double A = c.rt[km - 1] - c.rt[km], B = c.rt[km + 1] - c.rt[km];
    double d0 = y0 - y1, d2 = y2 - y1;
    if (A == 0.0 || B == 0.0 || A == B) { c.tarr = c.rt[km]; return; }
    double a = (d0 / A - d2 / B) / (A - B);
    double b = d0 / A - a * A;
    if (a >= 0.0) { c.tarr = c.rt[km]; return; }
    double sh = -0.5 * b / a;
    if (sh > B) sh = B; if (sh < A) sh = A;
    c.tarr = c.rt[km] + sh;
}
static void setup(Chart& c, int Nwant) {
    // Put BOTH the mirror and the dish exactly on grid points, in both charts.
    // Otherwise each chart's dish sits up to dx/2 away from areal radius 50 and
    // that sub-cell offset lands directly on the measured travel time -- an
    // O(dx) error whose sign wanders with N, which looks exactly like a code
    // that does not converge.
    double x0 = c.x_of_areal(g_mirror), xd = c.x_of_areal(g_dish);
    double x1 = c.x_of_areal(g_dish + R_PAD);
    int J = (int)(Nwant * (xd - x0) / (x1 - x0)); if (J < 8) J = 8;
    c.dx = (xd - x0) / J;
    int M = (int)std::ceil((x1 - xd) / c.dx);
    int N = J + M + 1;
    c.N = N; c.idet = J;
    c.x.assign(N, 0.0); c.uin.assign(N, 0.0); c.uout.assign(N, 0.0);
    c.a1.assign(N, 0.0); c.a2.assign(N, 0.0);
    c.s1.assign(N, 0.0); c.s2.assign(N, 0.0);
    double mx = 0.0;
    for (int i = 0; i < N; ++i) {
        c.x[i] = x0 + i * c.dx;
        double d = (c.areal(c.x[i]) - g_dish) / g_width;
        c.uin[i] = std::exp(-d * d);
        double a = std::fabs(lam_in(c, c.x[i])), b = std::fabs(lam_out(c, c.x[i]));
        if (a > mx) mx = a; if (b > mx) mx = b;
    }
    c.dt = 0.4 * c.dx / mx;
    c.t = 0.0; c.tarr = -1.0; c.pk = 0.0; c.p0 = c.p1 = c.p2 = 0.0;
    c.rin.assign(REC_N, 0.0); c.rout.assign(REC_N, 0.0); c.rt.assign(REC_N, -1.0);
    // exact quadrature for this pair of radii, in this chart's own variable
    c.quad = simpson_leg(c, x0, xd, 200000, 2);
    c.half_in = simpson_leg(c, x0, xd, 200000, 0);
}
static void peaks(Chart& c) {
    c.amp_in = 0; c.amp_out = 0; c.xpk_in = c.x[0]; c.xpk_out = c.x[0];
    for (int i = 0; i < c.N; ++i) {
        if (c.uin[i] > c.amp_in) { c.amp_in = c.uin[i]; c.xpk_in = c.x[i]; }
        if (c.uout[i] > c.amp_out) { c.amp_out = c.uout[i]; c.xpk_out = c.x[i]; }
    }
}

static void rebuild() {
    PG.alpha = pg_alpha; PG.beta = pg_beta; PG.gam = pg_gam;
    PG.areal = pg_areal; PG.x_of_areal = pg_x_of_areal;
    ISO.alpha = iso_alpha; ISO.beta = iso_beta; ISO.gam = iso_gam;
    ISO.areal = iso_areal; ISO.x_of_areal = iso_x_of_areal;
    setup(PG, g_N);
    setup(ISO, g_N);
    T_EXACT = 2.0 * ((g_dish - g_mirror) + 2.0 * std::log((g_dish - 2.0) / (g_mirror - 2.0)));
    dirty = false;
}

// ------------------------------------------------------------------- drawing
static void vline(Olivec_Canvas oc, int x, int y0, int y1, uint32_t col) {
    if (x < 0 || x >= FW) return;
    olivec_line(oc, x, y0, x, y1, col);
}
static void trace(Olivec_Canvas oc, const Chart& c, const std::vector<double>& u,
                  bool by_areal, double ax0, double ax1,
                  int px0, int px1, int base, int hgt, uint32_t col) {
    int lx = -1, ly = -1;
    for (int i = 0; i < c.N; ++i) {
        double xv = by_areal ? c.areal(c.x[i]) : c.x[i];
        double f = (xv - ax0) / (ax1 - ax0);
        if (f < 0.0 || f > 1.0) { lx = -1; continue; }
        int sx = px0 + (int)(f * (px1 - px0));
        double a = u[i]; if (a > 1.2) a = 1.2;
        int sy = base - (int)(a * hgt);
        if (lx >= 0) olivec_line(oc, lx, ly, sx, sy, col);
        lx = sx; ly = sy;
    }
}

extern "C" {

KEEP int sim_w() { return FW; }
KEEP int sim_h() { return FH; }
KEEP void sim_reset() {
    g_mirror = 6.0; g_dish = 50.0; g_speed = 1.0; g_N = 1600; dirty = true;
}
KEEP int sim_init(int, int) { px.assign((size_t)FW * FH, 0); sim_reset(); rebuild(); return 1; }
KEEP void sim_set(int id, double v) {
    if (id == 0) { g_mirror = v; dirty = true; }
    else if (id == 1) { g_dish = v; dirty = true; }
    else if (id == 2) { g_N = (int)v; dirty = true; }
    else if (id == 3) g_speed = v;
}
KEEP void sim_action(int a) {
    if (a == 0) { g_mirror = 3.0; dirty = true; }        // the photon sphere
    else if (a == 1) { g_mirror = 6.0; dirty = true; }
    else if (a == 2) { g_mirror = 2.2; dirty = true; }   // just outside the horizon
    else if (a == 3) { g_N = 400; dirty = true; }
    else if (a == 4) { g_N = 3200; dirty = true; }
}
KEEP void sim_click(double nx, double) {
    double r = g_mirror + nx * (g_dish + R_PAD - g_mirror);
    if (r < 2.05) r = 2.05; if (r > g_dish - 4.0) r = g_dish - 4.0;
    g_mirror = r; dirty = true;
}
KEEP void sim_step(int steps) {
    if (dirty) rebuild();
    double want = 0.55 * g_speed * (steps > 0 ? steps : 1);
    double tgt = PG.t + want;
    for (int k = 0; k < 4000 && PG.t < tgt; ++k) step(PG);
    tgt = ISO.t + want;
    for (int k = 0; k < 4000 && ISO.t < tgt; ++k) step(ISO);
    peaks(PG); peaks(ISO);
    measure(PG); measure(ISO);
}

KEEP uint8_t* sim_render() {
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(3, 5, 11, 1.f));
    Olivec_Font ft = olivec_default_font;
    char buf[220];
    const uint32_t C_PG = rgba(80, 235, 255, 1.f), C_PGO = rgba(40, 150, 190, 1.f);
    const uint32_t C_IS = rgba(255, 130, 230, 1.f), C_ISO2 = rgba(180, 70, 160, 1.f);
    const uint32_t C_AX = rgba(40, 60, 72, 1.f), C_MARK = rgba(230, 180, 55, 1.f);

    olivec_text(oc, "radar. two charts, one answer. same 3 plus 1 formula, different alpha beta gamma.",
                12, 6, ft, 1, rgba(150, 230, 220, 1.f));
    olivec_text(oc, "pg   alpha 1, beta sqrt 2m over r, gamma 1        light speed  plus minus 1 - beta",
                12, 20, ft, 1, C_PG);
    olivec_text(oc, "iso  alpha a, beta 0, gamma b to the 4th          light speed  plus minus alpha over b squared",
                12, 32, ft, 1, C_IS);

    // ---------------- panel 1: each chart on its own grid ----------------
    {
        const int TOP = 52, BOT = 190, HW = FW / 2;
        olivec_text(oc, "each chart on its own coordinate. different grids, different light speeds.",
                    12, TOP - 8, ft, 1, rgba(120, 150, 160, 1.f));
        for (int s = 0; s < 2; ++s) {
            Chart& c = s ? ISO : PG;
            int px0 = s ? HW + 8 : 8, px1 = s ? FW - 8 : HW - 8;
            double ax0 = c.x[0], ax1 = c.x[c.N - 1];
            olivec_line(oc, px0, BOT, px1, BOT, C_AX);
            for (int k = 0; k <= 10; ++k) {              // grid spacing ticks
                int sx = px0 + k * (px1 - px0) / 10;
                olivec_line(oc, sx, BOT, sx, BOT + 4, C_AX);
            }
            double f = (c.x_of_areal(g_dish) - ax0) / (ax1 - ax0);
            vline(oc, px0 + (int)(f * (px1 - px0)), TOP + 8, BOT, rgba(120, 120, 60, 1.f));
            vline(oc, px0, TOP + 8, BOT, C_MARK);        // the mirror is the inner edge
            trace(oc, c, c.uin, false, ax0, ax1, px0, px1, BOT, 118, s ? C_IS : C_PG);
            trace(oc, c, c.uout, false, ax0, ax1, px0, px1, BOT, 118, s ? C_ISO2 : C_PGO);
            std::snprintf(buf, sizeof buf, s ? "iso   rho from %.4f to %.4f, dx %.5f, dt %.5f"
                                             : "pg    r   from %.4f to %.4f, dx %.5f, dt %.5f",
                          ax0, ax1, c.dx, c.dt);
            olivec_text(oc, buf, px0 + 2, TOP - 8 + 12, ft, 1, s ? C_IS : C_PG);
            std::snprintf(buf, sizeof buf, "t %.3f", c.t);
            olivec_text(oc, buf, px1 - 70, BOT - 132, ft, 1, rgba(150, 170, 180, 1.f));
        }
    }

    // ---------------- panel 2: both against the areal radius ----------------
    {
        const int TOP = 214, BOT = 352;
        olivec_text(oc, "both against the areal radius. mid flight they disagree, because now is a choice.",
                    12, TOP - 10, ft, 1, rgba(120, 150, 160, 1.f));
        double ax0 = g_mirror, ax1 = g_dish + R_PAD;
        int px0 = 12, px1 = FW - 12;
        olivec_line(oc, px0, BOT, px1, BOT, C_AX);
        for (int k = 0; k <= 10; ++k) {
            int sx = px0 + k * (px1 - px0) / 10;
            olivec_line(oc, sx, BOT, sx, BOT + 4, C_AX);
            std::snprintf(buf, sizeof buf, "%.0f", ax0 + k * (ax1 - ax0) / 10);
            olivec_text(oc, buf, sx - 6, BOT + 6, ft, 1, rgba(90, 110, 120, 1.f));
        }
        vline(oc, px0, TOP, BOT, C_MARK);
        {
            double f = (g_dish - ax0) / (ax1 - ax0);
            vline(oc, px0 + (int)(f * (px1 - px0)), TOP, BOT, rgba(120, 120, 60, 1.f));
        }
        trace(oc, PG, PG.uin, true, ax0, ax1, px0, px1, BOT, 116, C_PG);
        trace(oc, PG, PG.uout, true, ax0, ax1, px0, px1, BOT, 116, C_PGO);
        trace(oc, ISO, ISO.uin, true, ax0, ax1, px0, px1, BOT, 116, C_IS);
        trace(oc, ISO, ISO.uout, true, ax0, ax1, px0, px1, BOT, 116, C_ISO2);
        // the gap between the two charts' pulse peaks, drawn
        double ap = (PG.amp_in > PG.amp_out) ? PG.areal(PG.xpk_in) : PG.areal(PG.xpk_out);
        double ai = (ISO.amp_in > ISO.amp_out) ? ISO.areal(ISO.xpk_in) : ISO.areal(ISO.xpk_out);
        int sxp = px0 + (int)((ap - ax0) / (ax1 - ax0) * (px1 - px0));
        int sxi = px0 + (int)((ai - ax0) / (ax1 - ax0) * (px1 - px0));
        olivec_line(oc, sxp, BOT - 128, sxi, BOT - 128, rgba(210, 190, 90, 1.f));
        std::snprintf(buf, sizeof buf, "same t, areal radius  pg %.3f   iso %.3f   apart by %.3f",
                      ap, ai, std::fabs(ap - ai));
        olivec_text(oc, buf, px0 + 4, BOT - 126, ft, 1, rgba(210, 190, 90, 1.f));
        olivec_text(oc, "bright  on the way in.   dim  on the way back.",
                    px0 + 4, BOT - 114, ft, 1, rgba(110, 130, 140, 1.f));
    }

    // ---------------- panel 3: the dish record ----------------
    {
        const int TOP = 382, BOT = 486;
        olivec_text(oc, "what the dish records. one clock, two readings. this is the observable.",
                    12, TOP - 12, ft, 1, rgba(120, 150, 160, 1.f));
        int px0 = 12, px1 = FW - 12;
        double t1 = T_EXACT * 1.25 + 6.0;
        olivec_line(oc, px0, BOT, px1, BOT, C_AX);
        for (int k = 0; k <= 10; ++k) {
            int sx = px0 + k * (px1 - px0) / 10;
            olivec_line(oc, sx, BOT, sx, BOT + 4, C_AX);
            if (k == 10) continue;
            std::snprintf(buf, sizeof buf, "%.0f", k * t1 / 10);
            olivec_text(oc, buf, sx - 6, BOT + 6, ft, 1, rgba(90, 110, 120, 1.f));
        }
        {
            int sx = px0 + (int)(T_EXACT / t1 * (px1 - px0));
            vline(oc, sx, TOP, BOT, rgba(120, 120, 60, 1.f));
            olivec_text(oc, "exact", sx - 16, TOP - 2, ft, 1, rgba(150, 150, 80, 1.f));
        }
        for (int s = 0; s < 2; ++s) {
            Chart& c = s ? ISO : PG;
            uint32_t col = s ? C_IS : C_PG, cold = s ? C_ISO2 : C_PGO;
            int lx = -1, ly = -1, lx2 = -1, ly2 = -1;
            for (int k = 0; k < REC_N; ++k) {
                double tt = k * REC_DT; if (tt > t1) break;
                int sx = px0 + (int)(tt / t1 * (px1 - px0));
                int sy = BOT - (int)(c.rout[k] * 88.0);
                int sy2 = BOT - (int)(c.rin[k] * 88.0);
                if (lx >= 0) olivec_line(oc, lx, ly, sx, sy, col);
                if (lx2 >= 0) olivec_line(oc, lx2, ly2, sx, sy2, cold);
                lx = sx; ly = sy; lx2 = sx; ly2 = sy2;
            }
        }
        olivec_text(oc, "bright  the pulse coming back.  dim  the pulse leaving.",
                    px0 + 4, BOT - 96, ft, 1, rgba(110, 130, 140, 1.f));
    }

    // ---------------- numbers ----------------
    std::snprintf(buf, sizeof buf, "round trip   pg %.12f   iso %.12f", PG.quad, ISO.quad);
    olivec_text(oc, buf, 12, FH - 116, ft, 2, rgba(180, 210, 150, 1.f));
    std::snprintf(buf, sizeof buf, "exact %.12f      pg vs iso %.2e",
                  T_EXACT, std::fabs(PG.quad - ISO.quad));
    olivec_text(oc, buf, 12, FH - 92, ft, 2, rgba(180, 210, 150, 1.f));
    std::snprintf(buf, sizeof buf,
                  "turned around at  pg %.6f   iso %.6f   apart by %.6f",
                  PG.half_in, ISO.half_in, std::fabs(PG.half_in - ISO.half_in));
    olivec_text(oc, buf, 12, FH - 66, ft, 2, rgba(235, 190, 90, 1.f));
    if (PG.tarr < 0.0 || ISO.tarr < 0.0)
        std::snprintf(buf, sizeof buf, "grid measured return   still on its way out.   n %d", g_N);
    else
        std::snprintf(buf, sizeof buf, "grid measured return   pg %.6f   iso %.6f   apart %.1e",
                      PG.tarr, ISO.tarr, std::fabs(PG.tarr - ISO.tarr));
    olivec_text(oc, buf, 12, FH - 40, ft, 2, rgba(140, 210, 230, 1.f));
    std::snprintf(buf, sizeof buf,
                  "dish clock %.6f   mirror areal %.3f   dish areal %.3f",
                  T_EXACT * std::sqrt(1.0 - 2.0 / g_dish), g_mirror, g_dish);
    olivec_text(oc, buf, 12, FH - 16, ft, 2, rgba(150, 170, 185, 1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
int main(int argc, char** argv) {
    sim_init(0, 0);
    if (argc > 2) { g_mirror = atof(argv[2]); dirty = true; rebuild(); }
    printf("mirror at areal %.3f, dish at areal %.3f\n", g_mirror, g_dish);
    printf("round trip coordinate time, each chart in its own variable\n");
    printf("  pg   alpha 1, beta sqrt(2m/r), gamma 1     %.12f\n", PG.quad);
    printf("  iso  alpha A, beta 0, gamma B^4            %.12f\n", ISO.quad);
    printf("  exact 2[(r1-r2) + 4m ln((r1-2m)/(r2-2m))]  %.12f\n", T_EXACT);
    printf("  pg vs iso %.3e    pg vs exact %.3e    iso vs exact %.3e\n",
           std::fabs(PG.quad - ISO.quad), std::fabs(PG.quad - T_EXACT),
           std::fabs(ISO.quad - T_EXACT));
    printf("turn-around, the halfway event the two charts do NOT share\n");
    printf("  pg  %.9f   iso %.9f   apart by %.9f\n",
           PG.half_in, ISO.half_in, std::fabs(PG.half_in - ISO.half_in));
    // The honest test for two different discretisations: do they converge to
     // the same number, and at the order the scheme promises?
    printf("convergence of the grid measurement, same experiment on both grids\n");
    double ep0 = 0, ei0 = 0;
    for (int k = 0; k < 4; ++k) {
        g_N = 400 << k; dirty = true; rebuild();
        while (PG.t < T_EXACT + 6.0) step(PG);
        while (ISO.t < T_EXACT + 6.0) step(ISO);
        measure(PG); measure(ISO);
        double ep = std::fabs(PG.tarr - T_EXACT), ei = std::fabs(ISO.tarr - T_EXACT);
        printf("  n %5d   pg %.6f err %.3e    iso %.6f err %.3e    pg vs iso %.3e",
               g_N, PG.tarr, ep, ISO.tarr, ei, std::fabs(PG.tarr - ISO.tarr));
        if (k) printf("   order pg %.2f iso %.2f", std::log2(ep0 / ep), std::log2(ei0 / ei));
        printf("\n");
        ep0 = ep; ei0 = ei;
    }
    printf("  returning peak amplitude  pg %.6f   iso %.6f   (1 means no loss)\n",
           PG.pk, ISO.pk);
    g_N = 1600; dirty = true; rebuild();
    int frames = argc > 1 ? atoi(argv[1]) : 100;
    for (int i = 0; i < frames; ++i) sim_step(1);
    printf("preview frame at t = %.4f, mid flight\n", PG.t);
    printf("  pulse peak areal radius  pg %.4f   iso %.4f   apart by %.4f\n",
           PG.amp_in > PG.amp_out ? PG.areal(PG.xpk_in) : PG.areal(PG.xpk_out),
           ISO.amp_in > ISO.amp_out ? ISO.areal(ISO.xpk_in) : ISO.areal(ISO.xpk_out),
           std::fabs((PG.amp_in > PG.amp_out ? PG.areal(PG.xpk_in) : PG.areal(PG.xpk_out))
                     - (ISO.amp_in > ISO.amp_out ? ISO.areal(ISO.xpk_in) : ISO.areal(ISO.xpk_out))));
    uint8_t* p = sim_render();
    stbi_write_png("radar_os_preview.png", FW, FH, 4, p, FW * 4);
    printf("wrote radar_os_preview.png\n");
    return 0;
}
#endif
