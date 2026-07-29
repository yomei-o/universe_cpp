// Universe OS — Radar on the AMR machine  (C++/WASM)
//
// The radar round trip of radar_os, put on the Berger-Oliger machinery of
// amr_os. Two charts of the same Schwarzschild spacetime, each with its own
// grid hierarchy, both timing the same experiment: a light pulse from the dish
// at areal radius 50m, off a mirror parked just outside the horizon, and back.
//
// Why refinement is needed here has nothing to do with where the pulse is. It
// is about where the BLACK HOLE is. The out-going coordinate light speed
// collapses near the horizon,
//
//     pg    1 - sqrt(2m/r)      0.800 at the dish, 0.0123 at a 2.05m mirror
//     iso   (1-u)/(1+u)^3       0.960 at the dish, 0.0522 at the same place
//
// so an out-going pulse gets squeezed by 65x (pg) or 18x (iso) on its way in,
// and a uniform grid would have to use the smallest spacing everywhere. The
// hole does not move, so the boxes do not either: this is fixed mesh
// refinement, which is what a real run for a single stationary hole uses. The
// moving-box machinery lives in amr_os.
//
// Each level shares the INNER edge with the physical boundary, because that is
// where the mirror is. So the mirror condition (the out-going family's inflow
// boundary is fed by the in-going family) is applied on every level, and only
// the OUTER edge of a refined level is an AMR interface. Refinement touching a
// physical boundary is ordinary in real codes.
//
// Ground truth is exact: the round trip is
//     2[(r1-r2) + 4m ln((r1-2m)/(r2-2m))] = 123.367733137848
// and both charts' own quadratures agree with it to 1e-12 from integrands that
// share nothing.
//
// The result worth reading is at the bottom of the screen: the two charts need
// DIFFERENT hierarchies for the same physics -- and the c*t = const chart is
// the cheaper of the two here. Not a law, just this problem with this mirror.
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

// ------------------------------------------------------------------ the charts
static double pg_alpha(double) { return 1.0; }
static double pg_beta(double r) { return std::sqrt(2.0 / r); }
static double pg_gam(double) { return 1.0; }
static double pg_areal(double r) { return r; }
static double pg_x_of_areal(double r) { return r; }

static inline double iso_u(double p) { return 0.5 / p; }
static double iso_alpha(double p) { double u = iso_u(p); return (1.0 - u) / (1.0 + u); }
static double iso_beta(double) { return 0.0; }
static double iso_gam(double p) { double B = 1.0 + iso_u(p); return B * B * B * B; }
static double iso_areal(double p) { double B = 1.0 + iso_u(p); return p * B * B; }
static double iso_x_of_areal(double r) {
    double a = r - 1.0, d = a * a - 1.0; if (d < 0.0) d = 0.0;
    return 0.5 * (a + std::sqrt(d));
}

static const int NG = 3;
static const double CFL = 0.4;
static double KO = 0.005;
static const double REC_DT = 0.05;
static const int REC_N = 4200;

static double g_mirror = 2.05, g_dish = 50.0, g_width = 4.0, g_speed = 1.0;
static const double R_PAD = 12.0;
static int g_base = 150;
static int g_extra = 0;                        // levels beyond what the squeeze asks for
static bool dirty = true;
static double T_EXACT = 0.0;

// ------------------------------------------------------------------- one level
struct Lev {
    double h = 0, dt = 0, t = 0;
    int n = 0;                                 // interior points, x(i) = x_in + i*h
    std::vector<double> ui, uo;
    std::vector<double> ui1, uo1, ui2, uo2;
    double t1 = 0, t2 = 0;
    bool have2 = false;
    long steps = 0;
};
struct Chart {
    const char* tag;
    double (*alpha)(double); double (*beta)(double); double (*gam)(double);
    double (*areal)(double); double (*x_of_areal)(double);
    int NL = 1;
    std::vector<Lev> lv;
    double x_in = 0, x_out = 0, h0 = 0, lam_max = 1.0;
    int idet = 0;                              // dish index on level 0
    std::vector<double> rin, rout, rt;
    double tarr = -1.0, quad = 0.0, squeeze = 1.0;
};
static Chart PG, ISO;

static inline double lam_out(const Chart& c, double x) {
    return -c.beta(x) + c.alpha(x) / std::sqrt(c.gam(x));
}
static inline double lam_in(const Chart& c, double x) {
    return -c.beta(x) - c.alpha(x) / std::sqrt(c.gam(x));
}
static inline double xof(const Chart& c, int l, int i) { return c.x_in + i * c.lv[l].h; }

// ---------------------------------------------------------------- quadrature
static double simpson_rt(const Chart& c, double xa, double xb, int n) {
    if (n % 2) ++n;
    double h = (xb - xa) / n, s = 0.0;
    for (int i = 0; i <= n; ++i) {
        double x = xa + i * h;
        double f = 1.0 / std::fabs(lam_in(c, x)) + 1.0 / lam_out(c, x);
        s += ((i == 0 || i == n) ? 1.0 : ((i & 1) ? 4.0 : 2.0)) * f;
    }
    return s * h / 3.0;
}

// ---------------------------------------------------------------- derivatives
// The inner end is the mirror: a physical boundary on every level, so it takes
// one sided differences and never reads ghosts. The outer end of a REFINED level
// is an AMR interface, and there the stencil must reach into the ghost cells --
// that is the only path by which the parent's solution enters the child. Using
// one sided differences there instead (the easy mistake) turns the interface
// into an outflow boundary and the pulse can never get into the fine grid.
static void ddx(const std::vector<double>& u, int n, double h, std::vector<double>& d,
                bool outer_ghosts) {
    double c = 1.0 / (12.0 * h);
    auto K = [&](int i) { return (size_t)(i + NG); };
    int hi = outer_ghosts ? n : n - 2;
    for (int i = 2; i < hi; ++i) {
        size_t k = K(i);
        d[k] = (-u[k + 2] + 8.0 * u[k + 1] - 8.0 * u[k - 1] + u[k - 2]) * c;
    }
    d[K(0)] = (-25.0 * u[K(0)] + 48.0 * u[K(1)] - 36.0 * u[K(2)]
               + 16.0 * u[K(3)] - 3.0 * u[K(4)]) * c;
    d[K(1)] = (-3.0 * u[K(0)] - 10.0 * u[K(1)] + 18.0 * u[K(2)]
               - 6.0 * u[K(3)] + u[K(4)]) * c;
    if (!outer_ghosts) {
        d[K(n - 1)] = (25.0 * u[K(n - 1)] - 48.0 * u[K(n - 2)] + 36.0 * u[K(n - 3)]
                       - 16.0 * u[K(n - 4)] + 3.0 * u[K(n - 5)]) * c;
        d[K(n - 2)] = (3.0 * u[K(n - 1)] + 10.0 * u[K(n - 2)] - 18.0 * u[K(n - 3)]
                       + 6.0 * u[K(n - 4)] - u[K(n - 5)]) * c;
    }
}
static void ko_filter(std::vector<double>& u, int n, double sig) {
    static std::vector<double> t; t = u;
    for (int i = 3; i < n - 3; ++i) {
        size_t k = (size_t)(i + NG);
        u[k] = t[k] + sig * (t[k - 3] - 6.0 * t[k - 2] + 15.0 * t[k - 1] - 20.0 * t[k]
                             + 15.0 * t[k + 1] - 6.0 * t[k + 2] + t[k + 3]);
    }
}

// ---------------------------------------------- prolongation from the parent
static double tinterp(double a2, double a1, double a0,
                      double t2, double t1, double t0, double t, bool have2) {
    double d21 = t1 - t2, d10 = t0 - t1;
    if (!have2 || std::fabs(d21) < 1e-14 || std::fabs(d10) < 1e-14
        || std::fabs(t0 - t2) < 1e-14) {
        double w = (std::fabs(d10) < 1e-14) ? 0.0 : (t - t1) / d10;
        return a1 + w * (a0 - a1);
    }
    double L2 = (t - t1) * (t - t0) / ((t2 - t1) * (t2 - t0));
    double L1 = (t - t2) * (t - t0) / ((t1 - t2) * (t1 - t0));
    double L0 = (t - t2) * (t - t1) / ((t0 - t2) * (t0 - t1));
    return a2 * L2 + a1 * L1 + a0 * L0;
}
// level l's point i sits at parent index i/2 when i is even, and halfway between
// two parent points when it is odd
static double parent_at(const Lev& pa, const std::vector<double>& cur,
                        const std::vector<double>& h1, const std::vector<double>& h2,
                        int i, double t) {
    auto val = [&](int J) -> double {
        if (J < -NG) J = -NG; if (J > pa.n - 1 + NG) J = pa.n - 1 + NG;
        size_t k = (size_t)(J + NG);
        return tinterp(h2[k], h1[k], cur[k], pa.t2, pa.t1, pa.t, t, pa.have2);
    };
    if ((i & 1) == 0) return val(i >> 1);
    int Jl = (i - 1) >> 1;
    return (-val(Jl - 1) + 9.0 * val(Jl) + 9.0 * val(Jl + 1) - val(Jl + 2)) / 16.0;
}
// the mirror: the out-going family's inflow edge is fed by the in-going family.
// applied on EVERY level, since every level reaches the physical inner edge.
static void bcs(Chart& c, int l, std::vector<double>& ui, std::vector<double>& uo, double t) {
    Lev& L = c.lv[l];
    uo[(size_t)(0 + NG)] = ui[(size_t)(0 + NG)];
    if (l == 0) {
        ui[(size_t)(L.n - 1 + NG)] = 0.0;      // nothing comes in from outside
        for (int k = 1; k <= NG; ++k) {        // outflow for the out-going family
            ui[(size_t)(L.n - 1 + k + NG)] = 0.0;
            uo[(size_t)(L.n - 1 + k + NG)] = uo[(size_t)(L.n - 1 + NG)];
        }
        return;
    }
    Lev& pa = c.lv[l - 1];                     // outer edge is an amr interface
    for (int k = 1; k <= NG; ++k) {
        int i = L.n - 1 + k;
        ui[(size_t)(i + NG)] = parent_at(pa, pa.ui, pa.ui1, pa.ui2, i, t);
        uo[(size_t)(i + NG)] = parent_at(pa, pa.uo, pa.uo1, pa.uo2, i, t);
    }
}
static void rk4(Chart& c, int l) {
    Lev& L = c.lv[l];
    const int n = L.n; const double h = L.h, dt = L.dt, t = L.t;
    static std::vector<double> yi, yo, ki, ko, ai, ao, di, dou;
    size_t m = (size_t)n + 2 * NG;
    yi = L.ui; yo = L.uo;
    ki.assign(m, 0.0); ko.assign(m, 0.0);
    ai.assign(m, 0.0); ao.assign(m, 0.0);
    di.assign(m, 0.0); dou.assign(m, 0.0);
    auto stage = [&](std::vector<double>& a, std::vector<double>& b, double tt, double w) {
        bcs(c, l, a, b, tt);
        ddx(a, n, h, di, l > 0); ddx(b, n, h, dou, l > 0);
        for (int i = 0; i < n; ++i) {
            size_t k = (size_t)(i + NG);
            double x = xof(c, l, i);
            di[k] = -lam_in(c, x) * di[k];
            dou[k] = -lam_out(c, x) * dou[k];
            ki[k] += w * di[k]; ko[k] += w * dou[k];
        }
    };
    stage(yi, yo, t, 1.0);
    for (size_t k = 0; k < m; ++k) { ai[k] = yi[k] + 0.5 * dt * di[k]; ao[k] = yo[k] + 0.5 * dt * dou[k]; }
    stage(ai, ao, t + 0.5 * dt, 2.0);
    for (size_t k = 0; k < m; ++k) { ai[k] = yi[k] + 0.5 * dt * di[k]; ao[k] = yo[k] + 0.5 * dt * dou[k]; }
    stage(ai, ao, t + 0.5 * dt, 2.0);
    for (size_t k = 0; k < m; ++k) { ai[k] = yi[k] + dt * di[k]; ao[k] = yo[k] + dt * dou[k]; }
    stage(ai, ao, t + dt, 1.0);
    for (int i = 0; i < n; ++i) {
        size_t k = (size_t)(i + NG);
        L.ui[k] = yi[k] + dt / 6.0 * ki[k];
        L.uo[k] = yo[k] + dt / 6.0 * ko[k];
    }
    L.t = t + dt;
    bcs(c, l, L.ui, L.uo, L.t);
    ko_filter(L.ui, n, KO); ko_filter(L.uo, n, KO);
    ++L.steps;
    if (l == 0) {                              // the dish record lives on level 0
        int k = (int)(L.t / REC_DT);
        if (k >= 0 && k < REC_N && c.rt[k] < 0.0) {
            c.rt[k] = L.t;
            c.rin[k] = L.ui[(size_t)(c.idet + NG)];
            c.rout[k] = L.uo[(size_t)(c.idet + NG)];
        }
    }
}
static void restrict_up(Chart& c, int l) {
    Lev& L = c.lv[l]; Lev& pa = c.lv[l - 1];
    for (int i = 0; i < L.n; i += 2) {
        int j = i >> 1;
        if (j >= pa.n) break;
        pa.ui[(size_t)(j + NG)] = L.ui[(size_t)(i + NG)];
        pa.uo[(size_t)(j + NG)] = L.uo[(size_t)(i + NG)];
    }
}
static void advance(Chart& c, int l) {
    Lev& L = c.lv[l];
    L.ui2 = L.ui1; L.uo2 = L.uo1; L.t2 = L.t1;
    L.ui1 = L.ui; L.uo1 = L.uo; L.t1 = L.t;
    rk4(c, l);
    L.have2 = (L.t1 > L.t2 + 1e-14);
    if (l + 1 < c.NL) {
        advance(c, l + 1);
        advance(c, l + 1);
        restrict_up(c, l + 1);
    }
}

// -------------------------------------------------------------------- setup
static void setup(Chart& c, int nl_extra) {
    c.x_in = c.x_of_areal(g_mirror);
    c.x_out = c.x_of_areal(g_dish + R_PAD);
    double xd = c.x_of_areal(g_dish);
    // how badly the out-going speed collapses on the way in decides the depth
    c.squeeze = lam_out(c, xd) / lam_out(c, c.x_in);
    int nl = (int)std::ceil(std::log2(c.squeeze)) + 1 + nl_extra;
    if (nl < 1) nl = 1; if (nl > 10) nl = 10;
    c.NL = nl;
    // level 0: put the dish exactly on a grid point, and every finer level
    // shares the inner edge, so its points are a subset refined by 2
    int J = g_base;
    c.h0 = (xd - c.x_in) / J;
    int n0 = J + (int)std::ceil((c.x_out - xd) / c.h0) + 1;
    c.idet = J;
    c.lam_max = 0.0;
    for (int i = 0; i < n0; ++i) {
        double x = c.x_in + i * c.h0;
        c.lam_max = std::fmax(c.lam_max, std::fmax(std::fabs(lam_in(c, x)), std::fabs(lam_out(c, x))));
    }
    c.lv.assign(c.NL, Lev());
    double W0 = (double)(n0 - 1) * c.h0;
    for (int l = 0; l < c.NL; ++l) {
        Lev& L = c.lv[l];
        L.h = c.h0 / (double)(1 << l);
        L.dt = CFL * L.h / c.lam_max;          // one lam_max, so subcycling stays 2:1
        L.t = 0.0; L.t1 = 0.0; L.t2 = 0.0; L.have2 = false; L.steps = 0;
        double W = (l == 0) ? W0 : W0 / (double)(1 << l);
        L.n = (int)std::llround(W / L.h) + 1;
        if (L.n > 4000) L.n = 4000;
        size_t m = (size_t)L.n + 2 * NG;
        L.ui.assign(m, 0.0); L.uo.assign(m, 0.0);
        for (int i = 0; i < L.n; ++i) {         // the pulse, as an areal profile
            double d = (c.areal(c.x_in + i * L.h) - g_dish) / g_width;
            L.ui[(size_t)(i + NG)] = std::exp(-d * d);
        }
        bcs(c, l, L.ui, L.uo, 0.0);
        L.ui1 = L.ui; L.uo1 = L.uo; L.ui2 = L.ui; L.uo2 = L.uo;
    }
    c.rin.assign(REC_N, 0.0); c.rout.assign(REC_N, 0.0); c.rt.assign(REC_N, -1.0);
    c.tarr = -1.0;
    c.quad = simpson_rt(c, c.x_in, xd, 200000);
}
static void measure(Chart& c) {
    // the two charts have different dt, so one of them skips record bins.
    // walk outward to the nearest bins that actually got a sample.
    int km = -1; double mx = 0.0;
    for (int k = 0; k < REC_N; ++k)
        if (c.rt[k] >= 0.0 && c.rout[k] > mx) { mx = c.rout[k]; km = k; }
    if (km < 0 || mx < 0.05) { c.tarr = -1.0; return; }
    int ka = -1, kb = -1;
    for (int k = km - 1; k >= 0; --k) if (c.rt[k] >= 0.0) { ka = k; break; }
    for (int k = km + 1; k < REC_N; ++k) if (c.rt[k] >= 0.0) { kb = k; break; }
    if (ka < 0 || kb < 0) { c.tarr = c.rt[km]; return; }
    double y0 = c.rout[ka], y1 = c.rout[km], y2 = c.rout[kb];
    double A = c.rt[ka] - c.rt[km], B = c.rt[kb] - c.rt[km];
    double d0 = y0 - y1, d2 = y2 - y1;
    if (A == 0.0 || B == 0.0 || A == B) { c.tarr = c.rt[km]; return; }
    double a = (d0 / A - d2 / B) / (A - B), b = d0 / A - a * A;
    if (a >= 0.0) { c.tarr = c.rt[km]; return; }
    double sh = -0.5 * b / a;
    if (sh > B) sh = B; if (sh < A) sh = A;
    c.tarr = c.rt[km] + sh;
}
static long points(const Chart& c) { long s = 0; for (auto& L : c.lv) s += L.n; return s; }
static long uniform_points(const Chart& c) {
    return (long)std::llround((c.x_out - c.x_in) / c.lv.back().h) + 1;
}
static double work_amr(const Chart& c) {
    double w = 0.0;
    for (int l = 0; l < c.NL; ++l) w += (double)c.lv[l].n * (double)(1 << l);
    return w;
}
static double work_uniform(const Chart& c) {
    return (double)uniform_points(c) * (double)(1 << (c.NL - 1));
}

static void rebuild() {
    PG.tag = "pg"; PG.alpha = pg_alpha; PG.beta = pg_beta; PG.gam = pg_gam;
    PG.areal = pg_areal; PG.x_of_areal = pg_x_of_areal;
    ISO.tag = "iso"; ISO.alpha = iso_alpha; ISO.beta = iso_beta; ISO.gam = iso_gam;
    ISO.areal = iso_areal; ISO.x_of_areal = iso_x_of_areal;
    setup(PG, g_extra); setup(ISO, g_extra);
    T_EXACT = 2.0 * ((g_dish - g_mirror) + 2.0 * std::log((g_dish - 2.0) / (g_mirror - 2.0)));
    dirty = false;
}

extern "C" {

KEEP int sim_w() { return FW; }
KEEP int sim_h() { return FH; }
KEEP void sim_reset() { g_mirror = 2.05; g_dish = 50.0; g_extra = 0; g_speed = 1.0; g_base = 150; dirty = true; }
KEEP int sim_init(int, int) { px.assign((size_t)FW * FH, 0); sim_reset(); rebuild(); return 1; }
KEEP void sim_set(int id, double v) {
    if (id == 0) { g_mirror = v; dirty = true; }
    else if (id == 1) { g_dish = v; dirty = true; }
    else if (id == 2) { g_extra = (int)v; dirty = true; }
    else if (id == 4) { g_base = (int)v; dirty = true; }
    else if (id == 3) g_speed = v;
}
KEEP void sim_action(int a) {
    if (a == 0) { g_mirror = 2.05; dirty = true; }
    else if (a == 1) { g_mirror = 3.0; dirty = true; }
    else if (a == 2) { g_mirror = 6.0; dirty = true; }
    else if (a == 3) { g_extra = 0; dirty = true; }
    else if (a == 4) dirty = true;
}
KEEP void sim_click(double, double) { }
KEEP void sim_step(int steps) {
    if (dirty) rebuild();
    double want = 0.12 * g_speed * (steps > 0 ? steps : 1);
    if (PG.lv[0].t < T_EXACT * 1.15) {
        double tgt = PG.lv[0].t + want;
        for (int k = 0; k < 20000 && PG.lv[0].t < tgt; ++k) advance(PG, 0);
        tgt = ISO.lv[0].t + want;
        for (int k = 0; k < 20000 && ISO.lv[0].t < tgt; ++k) advance(ISO, 0);
    }
    measure(PG); measure(ISO);
}

// ------------------------------------------------------------------ rendering
// The refined boxes live between areal radius 2.05 and 2.52 while the domain
// runs to 62, so on a linear axis the entire hierarchy is under one percent of
// the width and simply cannot be seen. Plot against log10(r - 2m) instead: the
// geometric nesting then comes out evenly spaced, which is also how anyone
// working on this actually looks at it.
static inline double axmap(double areal) { return std::log10(std::fmax(areal - 2.0, 1e-6)); }
static void plot(Olivec_Canvas oc, const Chart& c, int l, bool outgoing,
                 double a0, double a1, int px0, int px1, int base, int hgt, uint32_t col) {
    const Lev& L = c.lv[l];
    int lx = -1, ly = -1;
    for (int i = 0; i < L.n; ++i) {
        double av = axmap(c.areal(c.x_in + i * L.h));
        double fr = (av - a0) / (a1 - a0);
        if (fr < 0.0 || fr > 1.0) { lx = -1; continue; }
        int sx = px0 + (int)(fr * (px1 - px0));
        double v = outgoing ? L.uo[(size_t)(i + NG)] : L.ui[(size_t)(i + NG)];
        if (v > 1.3) v = 1.3; if (v < -0.3) v = -0.3;
        int sy = base - (int)(v * hgt);
        if (lx >= 0) olivec_line(oc, lx, ly, sx, sy, col);
        lx = sx; ly = sy;
    }
}

KEEP uint8_t* sim_render() {
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(3, 5, 11, 1.f));
    Olivec_Font ft = olivec_default_font;
    char buf[240];
    const uint32_t C_AX = rgba(38, 56, 68, 1.f);

    olivec_text(oc, "radar on the amr machine. one experiment, two charts, two different hierarchies.",
                12, 6, ft, 1, rgba(150, 230, 220, 1.f));
    olivec_text(oc, "refinement is not tracking the pulse. it is sitting on the hole, where the",
                12, 18, ft, 1, rgba(120, 150, 160, 1.f));
    olivec_text(oc, "out going light speed collapses. the axis is log of r minus 2m, or the",
                12, 30, ft, 1, rgba(120, 150, 160, 1.f));

    for (int s = 0; s < 2; ++s) {
        Chart& c = s ? ISO : PG;
        const int TOP = 48 + s * 168, BOT = TOP + 118;
        uint32_t base = s ? rgba(255, 130, 230, 1.f) : rgba(80, 235, 255, 1.f);
        auto lcol = [&](int l) {
            float u = c.NL > 1 ? (float)l / (float)(c.NL - 1) : 0.f;
            return s ? rgba((int)(120 + 135 * u), (int)(40 + 40 * u), (int)(150 + 80 * u), 1.f)
                     : rgba((int)(20 + 60 * u), (int)(110 + 140 * u), (int)(200 - 20 * u), 1.f);
        };
        std::snprintf(buf, sizeof buf,
                      "%s   squeeze %.1f times, so %d levels.  h from %.5f down to %.7f",
                      s ? "iso  flat grid, light speed varies" : "pg   flat space, grid stands still",
                      c.squeeze, c.NL, c.h0, c.lv.back().h);
        olivec_text(oc, buf, 12, TOP - 10, ft, 1, base);
        double a0 = axmap(g_mirror), a1 = axmap(g_dish + R_PAD);
        int px0 = 12, px1 = FW - 12;
        olivec_line(oc, px0, BOT, px1, BOT, C_AX);
        for (int k = 0; k < 5; ++k) {            // decade ticks in areal r minus 2m
            double v = -1.0 + k;
            double fr = (v - a0) / (a1 - a0);
            if (fr < 0.0 || fr > 1.0) continue;
            int sx = px0 + (int)(fr * (px1 - px0));
            olivec_line(oc, sx, BOT, sx, BOT + 3, C_AX);
            std::snprintf(buf, sizeof buf, "r minus 2m is 1e%.0f", v);
            if (s == 1) olivec_text(oc, buf, sx + 2, BOT + 5, ft, 1, rgba(80, 100, 112, 1.f));
        }
        // the boxes
        for (int l = 0; l < c.NL; ++l) {
            double aend = axmap(c.areal(c.x_in + (c.lv[l].n - 1) * c.lv[l].h));
            double fr = (aend - a0) / (a1 - a0); if (fr > 1.0) fr = 1.0;
            int y = BOT + 8 + l * 5;
            int b = px0 + (int)(fr * (px1 - px0));
            olivec_line(oc, px0, y, b, y, lcol(l));
            olivec_line(oc, b, y - 2, b, y + 2, lcol(l));
        }
        for (int l = 0; l < c.NL; ++l) {
            plot(oc, c, l, false, a0, a1, px0, px1, BOT, 96, lcol(l));
            plot(oc, c, l, true, a0, a1, px0, px1, BOT, 96,
                 s ? rgba(255, 190, 120, 1.f) : rgba(230, 190, 70, 1.f));
        }
        std::snprintf(buf, sizeof buf, "t %.3f   points %ld   a flat grid at the finest h needs %ld",
                      c.lv[0].t, points(c), uniform_points(c));
        olivec_text(oc, buf, px0 + 4, TOP + 2, ft, 1, rgba(120, 140, 155, 1.f));
    }

    // the dish record
    {
        const int TOP = 396, BOT = 470, px0 = 12, px1 = FW - 12;
        olivec_text(oc, "what the dish records. bright coming back, dim leaving.",
                    12, TOP - 10, ft, 1, rgba(120, 150, 160, 1.f));
        double t1 = T_EXACT * 1.15;
        olivec_line(oc, px0, BOT, px1, BOT, C_AX);
        int sx = px0 + (int)(T_EXACT / t1 * (px1 - px0));
        olivec_line(oc, sx, TOP, sx, BOT, rgba(120, 120, 60, 1.f));
        olivec_text(oc, "exact", sx - 16, TOP - 2, ft, 1, rgba(150, 150, 80, 1.f));
        for (int s = 0; s < 2; ++s) {
            Chart& c = s ? ISO : PG;
            uint32_t col = s ? rgba(255, 130, 230, 1.f) : rgba(80, 235, 255, 1.f);
            uint32_t dim = s ? rgba(150, 60, 130, 1.f) : rgba(40, 130, 160, 1.f);
            int lx = -1, ly = -1, lx2 = -1, ly2 = -1;
            for (int k = 0; k < REC_N; ++k) {
                if (c.rt[k] < 0.0) continue;
                double tt = c.rt[k]; if (tt > t1) break;
                int qx = px0 + (int)(tt / t1 * (px1 - px0));
                int qy = BOT - (int)(c.rout[k] * 64.0);
                int qy2 = BOT - (int)(c.rin[k] * 64.0);
                if (lx >= 0) olivec_line(oc, lx, ly, qx, qy, col);
                if (lx2 >= 0) olivec_line(oc, lx2, ly2, qx, qy2, dim);
                lx = qx; ly = qy; lx2 = qx; ly2 = qy2;
            }
        }
    }

    // numbers
    {
        int Y = FH - 136;
        std::snprintf(buf, sizeof buf, "round trip exact %.12f   mirror areal %.3f", T_EXACT, g_mirror);
        olivec_text(oc, buf, 12, Y, ft, 2, rgba(180, 210, 150, 1.f));
        std::snprintf(buf, sizeof buf, "quadrature   pg %.9f   iso %.9f", PG.quad, ISO.quad);
        olivec_text(oc, buf, 12, Y + 24, ft, 2, rgba(180, 210, 150, 1.f));
        if (PG.tarr > 0.0 && ISO.tarr > 0.0)
            std::snprintf(buf, sizeof buf, "on the grid  pg %.6f   iso %.6f   apart %.2e",
                          PG.tarr, ISO.tarr, std::fabs(PG.tarr - ISO.tarr));
        else
            std::snprintf(buf, sizeof buf, "on the grid  still on its way back");
        olivec_text(oc, buf, 12, Y + 48, ft, 2, rgba(140, 210, 230, 1.f));
        if (PG.tarr > 0.0 && ISO.tarr > 0.0)
            std::snprintf(buf, sizeof buf, "error vs exact   pg %.3e   iso %.3e",
                          std::fabs(PG.tarr - T_EXACT), std::fabs(ISO.tarr - T_EXACT));
        else
            std::snprintf(buf, sizeof buf, "error vs exact   waiting");
        olivec_text(oc, buf, 12, Y + 72, ft, 2, rgba(140, 210, 230, 1.f));
        std::snprintf(buf, sizeof buf, "work saved by amr   pg %.1f times   iso %.1f times",
                      work_uniform(PG) / work_amr(PG), work_uniform(ISO) / work_amr(ISO));
        olivec_text(oc, buf, 12, Y + 96, ft, 2, rgba(235, 190, 90, 1.f));
        std::snprintf(buf, sizeof buf,
                      "same physics, and iso wants %d levels where pg wants %d. this problem only.",
                      ISO.NL, PG.NL);
        olivec_text(oc, buf, 12, Y + 120, ft, 2, rgba(200, 200, 210, 1.f));
    }
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
int main(int argc, char** argv) {
    sim_init(0, 0);
    if (argc > 2) { g_mirror = atof(argv[2]); dirty = true; rebuild(); }
    if (argc > 3) { g_extra = atoi(argv[3]); dirty = true; rebuild(); }
    if (argc > 4) { g_base = atoi(argv[4]); dirty = true; rebuild(); }
    printf("radar on the amr machine. mirror areal %.3f, dish areal %.3f\n", g_mirror, g_dish);
    printf("round trip exact  %.12f\n", T_EXACT);
    printf("  pg  quadrature  %.12f   iso quadrature  %.12f\n", PG.quad, ISO.quad);
    for (int s = 0; s < 2; ++s) {
        Chart& c = s ? ISO : PG;
        printf("%s  squeeze %.2f  levels %d  h %.6f to %.8f  points %ld  flat needs %ld"
               "  work saved %.1f x\n",
               c.tag, c.squeeze, c.NL, c.h0, c.lv.back().h, points(c), uniform_points(c),
               work_uniform(c) / work_amr(c));
        for (int l = 0; l < c.NL; ++l)
            printf("    level %d  n %5d  h %.8f  dt %.8f  out to areal %.4f\n",
                   l, c.lv[l].n, c.lv[l].h, c.lv[l].dt,
                   c.areal(c.x_in + (c.lv[l].n - 1) * c.lv[l].h));
    }
    int frames = argc > 1 ? atoi(argv[1]) : 400;
    for (int i = 0; i < frames; ++i) sim_step(1);
    printf("\nt now  pg %.4f  iso %.4f\n", PG.lv[0].t, ISO.lv[0].t);
    printf("steps per level, pg:");
    for (int l = 0; l < PG.NL; ++l) printf("  %ld", PG.lv[l].steps);
    printf("\n  ratio to level 0:");
    for (int l = 0; l < PG.NL; ++l)
        printf("  %.2f", (double)PG.lv[l].steps / (double)PG.lv[0].steps);
    printf("\nmeasured round trip\n");
    printf("  pg  %.6f  error %.3e\n", PG.tarr, std::fabs(PG.tarr - T_EXACT));
    printf("  iso %.6f  error %.3e\n", ISO.tarr, std::fabs(ISO.tarr - T_EXACT));
    for (int s2 = 0; s2 < 2; ++s2) {
        Chart& c = s2 ? ISO : PG;
        double mx = 0.0; int km = -1;
        for (int k = 0; k < REC_N; ++k) if (c.rout[k] > mx) { mx = c.rout[k]; km = k; }
        printf("  %s returning peak %.6f at t %.4f   in going peak seen %.6f\n", c.tag, mx,
               km >= 0 ? c.rt[km] : -1.0,
               [&] { double q = 0; for (int k = 0; k < REC_N; ++k) q = std::fmax(q, c.rin[k]);
                     return q; }());
    }
    printf("  the two grids differ by %.3e\n", std::fabs(PG.tarr - ISO.tarr));
    uint8_t* p = sim_render();
    stbi_write_png("radaramr_os_preview.png", FW, FH, 4, p, FW * 4);
    printf("wrote radaramr_os_preview.png\n");
    return 0;
}
#endif
