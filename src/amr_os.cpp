// Universe OS — AMR: how numerical relativity actually computes  (C++/WASM)
//
// This demo is not about Einstein's equations. It is about the MACHINE that
// numerical relativity codes are built out of, which is independent of what
// equations you feed it. Swap the right hand side for BSSN and this is a
// numerical relativity code.
//
// The equation here is deliberately the simplest hyperbolic system there is,
// the 1D wave equation in first order form:
//
//     dt f = p ,   dt p = dx s ,   dt s = dx p          characteristic speeds +-1
//
// with a closed form solution to check against: f = g(x-t), p = -g'(x-t),
// s = g'(x-t). And note that s - dx f = 0 is a CONSTRAINT -- true analytically,
// not enforced by the evolution, drifting numerically. That is exactly the
// standing of the hamiltonian constraint H = 0 in a real code, so we monitor it
// the same way.
//
// The machinery, all of which is what real codes do:
//
//   method of lines        4th order centred differences in space, RK4 in time,
//                          Kreiss-Oliger dissipation, CFL from the char. speeds
//   berger-oliger amr      a hierarchy of levels, refinement factor 2, and
//                          SUBCYCLING IN TIME: level l takes 2^l steps per
//                          coarse step, because its dt is 2^l times smaller
//   prolongation           a fine level's ghost zones are filled from its
//                          parent, interpolated 4th order in SPACE and 2nd
//                          order in TIME (the fine level needs parent values at
//                          instants the parent never visited). this is the part
//                          real codes most often get wrong
//   restriction            when a fine level catches up, its interior is
//                          injected back into the parent
//   nesting and buffers    each box sits strictly inside its parent with a
//                          buffer, so interface error cannot reach the feature
//                          before the next regrid
//   regridding             boxes move to follow the feature, and newly exposed
//                          cells are filled by prolongation. this is the same
//                          mechanism that lets punctures move across the grid
//
// What it costs, and why 3D changes the answer: in 1D this hierarchy saves a
// factor of tens over a uniform grid at the finest spacing. The SAME hierarchy
// in 3D saves a factor of ~10^5, which is the difference between a run that
// finishes and one that never starts.
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

// ----------------------------------------------------------------- the problem
static const double XMIN = 0.0, XMAX = 240.0;
static const double H0 = 1.0;                  // coarsest spacing
static const double CFL = 0.4;
static double PW = 0.4;                 // pulse width: 0.3 cells on level 0
static const double X0 = 12.0;                  // where it starts
static const int NG = 3;                       // ghost points per side
static double KO = 0.005;                      // kreiss-oliger strength

static inline double g_of(double z) { double q = z / PW; return std::exp(-q * q); }
static inline double gp_of(double z) { double q = z / PW; return -2.0 * q / PW * std::exp(-q * q); }
static inline double ex_f(double x, double t) { return g_of(x - X0 - t); }
static inline double ex_p(double x, double t) { return -gp_of(x - X0 - t); }
static inline double ex_s(double x, double t) { return gp_of(x - X0 - t); }

// ------------------------------------------------------------------- one level
// A level is a uniform patch. Its interior points are indexed 0..n-1 and its
// position on the global lattice of spacing h is Ilo (an integer), so
// x(i) = XMIN + (Ilo + i) * h. Ilo is kept EVEN and n ODD, which puts both ends
// of the patch on the parent's grid -- the nesting condition, in one line.
struct Lev {
    int l = 0;
    double h = 0, dt = 0, t = 0;
    long Ilo = 0;
    int n = 0;
    std::vector<double> f, p, s;
    std::vector<double> f1, p1, s1, f2, p2, s2;   // history for time interpolation
    double t1 = 0, t2 = 0;
    bool have2 = false;
    long steps = 0;

    void alloc(int nn) {
        n = nn;
        size_t m = (size_t)n + 2 * NG;
        f.assign(m, 0.0); p.assign(m, 0.0); s.assign(m, 0.0);
        f1 = f; p1 = p; s1 = s; f2 = f; p2 = p; s2 = s;
        have2 = false;
    }
    inline double x(int i) const { return XMIN + (double)(Ilo + i) * h; }
    inline double& F(int i) { return f[(size_t)(i + NG)]; }
    inline double& P(int i) { return p[(size_t)(i + NG)]; }
    inline double& S(int i) { return s[(size_t)(i + NG)]; }
    inline double F(int i) const { return f[(size_t)(i + NG)]; }
    inline double P(int i) const { return p[(size_t)(i + NG)]; }
    inline double S(int i) const { return s[(size_t)(i + NG)]; }
};

static std::vector<Lev> LEV;                   // LEV[0] is the coarsest
static int g_nlev = 8;
static double g_speed = 1.0;
static int g_regrid = 1;                       // regrid every this many coarse steps
static int g_buf = 12;                         // buffer cells beyond the pulse
static bool g_show_coarse = true;
static bool dirty = true;

// the single level coarse run, for comparison
static Lev COARSE;
// and a record of who stepped when, for the berger-oliger diagram
static const int EVN = 520;
static int8_t ev[EVN];
static int evhead = 0;

// --------------------------------------------------------- spatial derivative
static void ddx(const std::vector<double>& u, int n, double h, std::vector<double>& d) {
    double c = 1.0 / (12.0 * h);
    for (int i = 0; i < n; ++i) {
        size_t k = (size_t)(i + NG);
        d[k] = (-u[k + 2] + 8.0 * u[k + 1] - 8.0 * u[k - 1] + u[k - 2]) * c;
    }
}

// -------------------------------------------------- prolongation from a parent
// 4th order in space (the classic -1/16, 9/16, 9/16, -1/16 midpoint stencil)
// and 2nd order in time through the parent's three stored instants.
static double tinterp(double a2, double a1, double a0, double t2, double t1, double t0,
                      double t, bool have2) {
    // Only three DISTINCT instants can carry a quadratic. Right after a level's
    // first step its history still has t2 == t1, and reaching for the quadratic
    // there divides by zero and returns nan -- which is exactly what happened,
    // and only showed up through the regrid, because that is the one place that
    // reads a parent's history after the flag has been set.
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
// value of the parent field at global fine index I and time t
static double parent_at(const Lev& pa, const std::vector<double>& cur,
                        const std::vector<double>& h1, const std::vector<double>& h2,
                        long I, double t) {
    auto val = [&](long J) -> double {          // parent point J, time interpolated
        int i = (int)(J - pa.Ilo);
        if (i < -NG) i = -NG; if (i > pa.n - 1 + NG) i = pa.n - 1 + NG;
        size_t k = (size_t)(i + NG);
        return tinterp(h2[k], h1[k], cur[k], pa.t2, pa.t1, pa.t, t, pa.have2);
    };
    if ((I & 1) == 0) return val(I >> 1);       // sits on a parent point
    long Jl = (I - 1) >> 1;                     // between Jl and Jl+1
    double m1 = val(Jl - 1), a0 = val(Jl), a1 = val(Jl + 1), p1 = val(Jl + 2);
    return (-m1 + 9.0 * a0 + 9.0 * a1 - p1) / 16.0;
}
static void fill_ghosts(int l, std::vector<double>& f, std::vector<double>& p,
                        std::vector<double>& s, double t) {
    Lev& L = LEV[l];
    if (l == 0) {                               // testbed outer boundary: exact
        for (int k = 1; k <= NG; ++k) {
            int i = -k;
            f[(size_t)(i + NG)] = ex_f(L.x(i), t);
            p[(size_t)(i + NG)] = ex_p(L.x(i), t);
            s[(size_t)(i + NG)] = ex_s(L.x(i), t);
            i = L.n - 1 + k;
            f[(size_t)(i + NG)] = ex_f(L.x(i), t);
            p[(size_t)(i + NG)] = ex_p(L.x(i), t);
            s[(size_t)(i + NG)] = ex_s(L.x(i), t);
        }
        return;
    }
    Lev& pa = LEV[l - 1];
    for (int k = 1; k <= NG; ++k) {
        for (int side = 0; side < 2; ++side) {
            int i = side ? (L.n - 1 + k) : (-k);
            long I = L.Ilo + i;
            f[(size_t)(i + NG)] = parent_at(pa, pa.f, pa.f1, pa.f2, I, t);
            p[(size_t)(i + NG)] = parent_at(pa, pa.p, pa.p1, pa.p2, I, t);
            s[(size_t)(i + NG)] = parent_at(pa, pa.s, pa.s1, pa.s2, I, t);
        }
    }
}

// ------------------------------------------------------------------ one rk4 step
// Kreiss-Oliger dissipation, and the ORDER of it matters. The 5 point (p=2)
// filter is O(h^3), which is bigger than a 4th order scheme's truncation error,
// so it dominates: its per-step factor on mode k is 1 - 16 sig sin^4(kh/2), and
// on the finest level that accumulates over 2^L times as many steps. With
// 10 points across the pulse it ate 24 percent of the amplitude.
// Real 4th order codes use the p=3, 7 point operator instead: O(h^5), factor
// 1 - 64 sig sin^6(kh/2), which leaves 4th order convergence intact.
static void ko_filter(std::vector<double>& u, int n, double sig) {
    static std::vector<double> t;
    t = u;
    for (int i = 0; i < n; ++i) {
        size_t k = (size_t)(i + NG);
        u[k] = t[k] + sig * (t[k - 3] - 6.0 * t[k - 2] + 15.0 * t[k - 1] - 20.0 * t[k]
                             + 15.0 * t[k + 1] - 6.0 * t[k + 2] + t[k + 3]);
    }
}
static void rk4(int l) {
    Lev& L = LEV[l];
    const int n = L.n;
    const double h = L.h, dt = L.dt, t = L.t;
    static std::vector<double> yf, yp, ys, kf, kp, ks, af, ap, as, df, dp, ds;
    size_t m = (size_t)n + 2 * NG;
    yf = L.f; yp = L.p; ys = L.s;
    kf.assign(m, 0.0); kp.assign(m, 0.0); ks.assign(m, 0.0);
    af.assign(m, 0.0); ap.assign(m, 0.0); as.assign(m, 0.0);
    df.assign(m, 0.0); dp.assign(m, 0.0); ds.assign(m, 0.0);

    auto stage = [&](std::vector<double>& uf, std::vector<double>& up,
                     std::vector<double>& us, double tt, double w) {
        fill_ghosts(l, uf, up, us, tt);
        ddx(us, n, h, dp);                      // dt p = dx s
        ddx(up, n, h, ds);                      // dt s = dx p
        for (int i = 0; i < n; ++i) {
            size_t k = (size_t)(i + NG);
            df[k] = up[k];                      // dt f = p
        }
        for (int i = 0; i < n; ++i) {
            size_t k = (size_t)(i + NG);
            kf[k] += w * df[k]; kp[k] += w * dp[k]; ks[k] += w * ds[k];
        }
    };
    stage(yf, yp, ys, t, 1.0);
    for (int i = -NG; i < n + NG; ++i) {
        size_t k = (size_t)(i + NG);
        af[k] = yf[k] + 0.5 * dt * df[k]; ap[k] = yp[k] + 0.5 * dt * dp[k];
        as[k] = ys[k] + 0.5 * dt * ds[k];
    }
    stage(af, ap, as, t + 0.5 * dt, 2.0);
    for (int i = -NG; i < n + NG; ++i) {
        size_t k = (size_t)(i + NG);
        af[k] = yf[k] + 0.5 * dt * df[k]; ap[k] = yp[k] + 0.5 * dt * dp[k];
        as[k] = ys[k] + 0.5 * dt * ds[k];
    }
    stage(af, ap, as, t + 0.5 * dt, 2.0);
    for (int i = -NG; i < n + NG; ++i) {
        size_t k = (size_t)(i + NG);
        af[k] = yf[k] + dt * df[k]; ap[k] = yp[k] + dt * dp[k]; as[k] = ys[k] + dt * ds[k];
    }
    stage(af, ap, as, t + dt, 1.0);
    for (int i = 0; i < n; ++i) {
        size_t k = (size_t)(i + NG);
        L.f[k] = yf[k] + dt / 6.0 * kf[k];
        L.p[k] = yp[k] + dt / 6.0 * kp[k];
        L.s[k] = ys[k] + dt / 6.0 * ks[k];
    }
    L.t = t + dt;
    fill_ghosts(l, L.f, L.p, L.s, L.t);
    ko_filter(L.f, n, KO); ko_filter(L.p, n, KO); ko_filter(L.s, n, KO);
    ++L.steps;
    ev[evhead] = (int8_t)l; evhead = (evhead + 1) % EVN;
}

// ----------------------------------------------------------------- restriction
static void restrict_to_parent(int l) {
    Lev& L = LEV[l];
    Lev& pa = LEV[l - 1];
    for (int i = 0; i < L.n; ++i) {
        long I = L.Ilo + i;
        if (I & 1) continue;                    // only points shared with the parent
        int j = (int)((I >> 1) - pa.Ilo);
        if (j < 0 || j >= pa.n) continue;
        pa.F(j) = L.F(i); pa.P(j) = L.P(i); pa.S(j) = L.S(i);
    }
}

#ifndef __EMSCRIPTEN__
static int nan_said = 0;
static void nan_check(int l, const char* where) {
    if (nan_said) return;
    Lev& L = LEV[l];
    for (int i = -NG; i < L.n + NG; ++i) {
        double v = L.F(i);
        if (!(v == v)) {
            nan_said = 1;
            printf("FIRST NAN  level %d  %s  i %d of n %d  x %.6f  t %.6f  step %ld\n",
                   l, where, i, L.n, L.x(i), L.t, L.steps);
            return;
        }
    }
}
#else
static void nan_check(int, const char*) {}
#endif

// -------------------------------------------------- berger-oliger recursive step
static void advance(int l) {
    Lev& L = LEV[l];
    L.f2 = L.f1; L.p2 = L.p1; L.s2 = L.s1; L.t2 = L.t1;
    L.f1 = L.f; L.p1 = L.p; L.s1 = L.s; L.t1 = L.t;
    rk4(l);
    nan_check(l, "after rk4");
    // the history is only usable for a quadratic once two steps have been taken,
    // so let the times say so rather than a flag that gets set too early
    L.have2 = (L.t1 > L.t2 + 1e-14);
    if (l + 1 < (int)LEV.size()) {
        advance(l + 1);                         // subcycling: two fine steps
        advance(l + 1);                         // per one coarse step
        restrict_to_parent(l + 1);
        nan_check(l, "after restriction from the child");
    }
}

// ----------------------------------------------------------------- regridding
static double track_center() {
    Lev& L = LEV.back();
    int im = 0; double mx = -1.0;
    for (int i = 0; i < L.n; ++i) { double a = std::fabs(L.F(i)); if (a > mx) { mx = a; im = i; } }
    if (im > 0 && im < L.n - 1) {
        double y0 = std::fabs(L.F(im - 1)), y1 = mx, y2 = std::fabs(L.F(im + 1));
        double den = y0 - 2.0 * y1 + y2;
        double sh = (den != 0.0) ? 0.5 * (y0 - y2) / den : 0.0;
        if (sh > 1.0) sh = 1.0; if (sh < -1.0) sh = -1.0;
        return L.x(im) + sh * L.h;
    }
    return L.x(im);
}
// half width of level l's box, in x. innermost holds the pulse plus a buffer,
// and each coarser level doubles it -- so every level costs the same number of
// points, which is the whole trick.
static double box_hw(int l, int nl) {
    double inner = g_buf * (H0 / (double)(1 << (nl - 1))) + 4.0 * PW;
    return inner * (double)(1 << (nl - 1 - l));
}
static void make_level(int l, int nl, double xc, bool from_parent) {
    Lev& L = LEV[l];
    L.l = l;
    L.h = H0 / (double)(1 << l);
    L.dt = CFL * L.h;
    if (l == 0) {
        L.Ilo = 0;
        int nn = (int)std::llround((XMAX - XMIN) / L.h) + 1;
        if (!(nn & 1)) ++nn;
        if (!from_parent) { L.alloc(nn); L.t = 0.0; L.t1 = 0.0; L.t2 = 0.0; }
        return;
    }
    double hw = box_hw(l, nl);
    long Ic = (long)std::llround((xc - XMIN) / L.h);
    long lo = Ic - (long)std::llround(hw / L.h);
    long hi = Ic + (long)std::llround(hw / L.h);
    // keep the box strictly inside the parent's interior, and aligned
    Lev& pa = LEV[l - 1];
    long plo = 2 * pa.Ilo + 2 * NG, phi = 2 * (pa.Ilo + pa.n - 1) - 2 * NG;
    if (lo < plo) lo = plo;
    if (hi > phi) hi = phi;
    if (lo & 1) ++lo;
    if (!((hi - lo) & 1)) { /* n odd wants hi-lo even */ } else --hi;
    if (hi <= lo + 8) { hi = lo + 8; }
    int nn = (int)(hi - lo) + 1;

    long oIlo = L.Ilo; int on = L.n;
    bool had = (on > 0);
    size_t m = (size_t)nn + 2 * NG;
    // Carry the TIME HISTORY through the regrid as well. Dropping it (which is
    // the easy thing to do) forces this level's children back onto linear time
    // interpolation at every regrid, quietly costing an order.
    // keep what was already refined; prolongate only the newly exposed cells.
    // passing the same array in all three history slots makes parent_at a pure
    // spatial interpolation, which is what a history slot wants.
    auto regrid_one = [&](const std::vector<double>& old, const std::vector<double>& pcur) {
        std::vector<double> out(m, 0.0);
        for (int i = 0; i < nn; ++i) {
            long I = lo + i;
            size_t k = (size_t)(i + NG);
            int oi = (int)(I - oIlo);
            if (had && oi >= 0 && oi < on) out[k] = old[(size_t)(oi + NG)];
            else out[k] = parent_at(pa, pcur, pcur, pcur, I, pa.t);
        }
        return out;
    };
    std::vector<double> nf = regrid_one(L.f, pa.f);
    std::vector<double> np = regrid_one(L.p, pa.p);
    std::vector<double> ns = regrid_one(L.s, pa.s);
    std::vector<double> nf1 = regrid_one(L.f1, pa.f1);
    std::vector<double> np1 = regrid_one(L.p1, pa.p1);
    std::vector<double> ns1 = regrid_one(L.s1, pa.s1);
    std::vector<double> nf2 = regrid_one(L.f2, pa.f2);
    std::vector<double> np2 = regrid_one(L.p2, pa.p2);
    std::vector<double> ns2 = regrid_one(L.s2, pa.s2);
    L.Ilo = lo;
    L.n = nn;
    L.f = nf; L.p = np; L.s = ns;
    L.f1 = nf1; L.p1 = np1; L.s1 = ns1;
    L.f2 = nf2; L.p2 = np2; L.s2 = ns2;
    if (!had) { L.t = pa.t; L.t1 = L.t; L.t2 = L.t; L.have2 = false; }
}

static void init_all() {
    LEV.assign(g_nlev, Lev());
    for (int l = 0; l < g_nlev; ++l) {
        Lev& L = LEV[l];
        L.l = l; L.h = H0 / (double)(1 << l); L.dt = CFL * L.h;
        L.n = 0; L.steps = 0;
    }
    make_level(0, g_nlev, X0, false);
    for (int i = 0; i < LEV[0].n; ++i) {
        LEV[0].F(i) = ex_f(LEV[0].x(i), 0.0);
        LEV[0].P(i) = ex_p(LEV[0].x(i), 0.0);
        LEV[0].S(i) = ex_s(LEV[0].x(i), 0.0);
    }
    fill_ghosts(0, LEV[0].f, LEV[0].p, LEV[0].s, 0.0);
    LEV[0].f1 = LEV[0].f; LEV[0].p1 = LEV[0].p; LEV[0].s1 = LEV[0].s;
    LEV[0].f2 = LEV[0].f; LEV[0].p2 = LEV[0].p; LEV[0].s2 = LEV[0].s;
    for (int l = 1; l < g_nlev; ++l) {
        make_level(l, g_nlev, X0, true);
        for (int i = 0; i < LEV[l].n; ++i) {     // exact initial data on every level
            LEV[l].F(i) = ex_f(LEV[l].x(i), 0.0);
            LEV[l].P(i) = ex_p(LEV[l].x(i), 0.0);
            LEV[l].S(i) = ex_s(LEV[l].x(i), 0.0);
        }
        fill_ghosts(l, LEV[l].f, LEV[l].p, LEV[l].s, 0.0);
        LEV[l].f1 = LEV[l].f; LEV[l].p1 = LEV[l].p; LEV[l].s1 = LEV[l].s;
        LEV[l].f2 = LEV[l].f; LEV[l].p2 = LEV[l].p; LEV[l].s2 = LEV[l].s;
    }
    // the single level coarse run, same scheme, no refinement anywhere
    COARSE.l = 0; COARSE.h = H0; COARSE.dt = CFL * H0; COARSE.Ilo = 0; COARSE.t = 0;
    COARSE.alloc(LEV[0].n);
    for (int i = 0; i < COARSE.n; ++i) {
        COARSE.F(i) = ex_f(COARSE.x(i), 0.0);
        COARSE.P(i) = ex_p(COARSE.x(i), 0.0);
        COARSE.S(i) = ex_s(COARSE.x(i), 0.0);
    }
    for (int i = 0; i < EVN; ++i) ev[i] = -1;
    evhead = 0;
    dirty = false;
}
// the comparison run gets its own rk4, so the amr path stays untouched
static void rk4_single(Lev& L) {
    const int n = L.n; const double h = L.h, dt = L.dt, t = L.t;
    static std::vector<double> yf, yp, ys, kf, kp, ks, af, ap, as, df, dp, ds;
    size_t m = (size_t)n + 2 * NG;
    yf = L.f; yp = L.p; ys = L.s;
    kf.assign(m, 0.0); kp.assign(m, 0.0); ks.assign(m, 0.0);
    af.assign(m, 0.0); ap.assign(m, 0.0); as.assign(m, 0.0);
    df.assign(m, 0.0); dp.assign(m, 0.0); ds.assign(m, 0.0);
    auto bc = [&](std::vector<double>& uf, std::vector<double>& up,
                  std::vector<double>& us, double tt) {
        for (int k = 1; k <= NG; ++k) {
            int i = -k;
            uf[(size_t)(i + NG)] = ex_f(L.x(i), tt);
            up[(size_t)(i + NG)] = ex_p(L.x(i), tt);
            us[(size_t)(i + NG)] = ex_s(L.x(i), tt);
            i = n - 1 + k;
            uf[(size_t)(i + NG)] = ex_f(L.x(i), tt);
            up[(size_t)(i + NG)] = ex_p(L.x(i), tt);
            us[(size_t)(i + NG)] = ex_s(L.x(i), tt);
        }
    };
    auto stage = [&](std::vector<double>& uf, std::vector<double>& up,
                     std::vector<double>& us, double tt, double w) {
        bc(uf, up, us, tt);
        ddx(us, n, h, dp); ddx(up, n, h, ds);
        for (int i = 0; i < n; ++i) { size_t k = (size_t)(i + NG); df[k] = up[k]; }
        for (int i = 0; i < n; ++i) {
            size_t k = (size_t)(i + NG);
            kf[k] += w * df[k]; kp[k] += w * dp[k]; ks[k] += w * ds[k];
        }
    };
    stage(yf, yp, ys, t, 1.0);
    for (int i = -NG; i < n + NG; ++i) { size_t k = (size_t)(i + NG);
        af[k] = yf[k] + 0.5 * dt * df[k]; ap[k] = yp[k] + 0.5 * dt * dp[k]; as[k] = ys[k] + 0.5 * dt * ds[k]; }
    stage(af, ap, as, t + 0.5 * dt, 2.0);
    for (int i = -NG; i < n + NG; ++i) { size_t k = (size_t)(i + NG);
        af[k] = yf[k] + 0.5 * dt * df[k]; ap[k] = yp[k] + 0.5 * dt * dp[k]; as[k] = ys[k] + 0.5 * dt * ds[k]; }
    stage(af, ap, as, t + 0.5 * dt, 2.0);
    for (int i = -NG; i < n + NG; ++i) { size_t k = (size_t)(i + NG);
        af[k] = yf[k] + dt * df[k]; ap[k] = yp[k] + dt * dp[k]; as[k] = ys[k] + dt * ds[k]; }
    stage(af, ap, as, t + dt, 1.0);
    for (int i = 0; i < n; ++i) { size_t k = (size_t)(i + NG);
        L.f[k] = yf[k] + dt / 6.0 * kf[k];
        L.p[k] = yp[k] + dt / 6.0 * kp[k];
        L.s[k] = ys[k] + dt / 6.0 * ks[k]; }
    L.t = t + dt;
    bc(L.f, L.p, L.s, L.t);
    ko_filter(L.f, n, KO); ko_filter(L.p, n, KO); ko_filter(L.s, n, KO);
}

// -------------------------------------------------------------------- measures
// Max error on the finest level -- but only meaningful if the pulse is actually
// inside the finest box. If the boxes have lost it, both the solution and the
// exact answer are zero there and the "error" comes out beautifully small while
// meaning nothing. Ask first.
static bool amr_covers_pulse() {
    Lev& L = LEV.back();
    double xe = X0 + LEV[0].t;
    return xe > L.x(0) + 2.0 * PW && xe < L.x(L.n - 1) - 2.0 * PW;
}
static double err_amr() {
    Lev& L = LEV.back();
    if (!amr_covers_pulse()) return -1.0;
    double e = 0.0;
    for (int i = 0; i < L.n; ++i) {
        double d = std::fabs(L.F(i) - ex_f(L.x(i), L.t));
        if (d > e) e = d;
    }
    return e;
}
static double err_coarse() {
    double e = 0.0;
    for (int i = 0; i < COARSE.n; ++i) {
        double d = std::fabs(COARSE.F(i) - ex_f(COARSE.x(i), COARSE.t));
        if (d > e) e = d;
    }
    return e;
}
// the constraint s - dx f = 0, level by level. true analytically, not enforced.
static double constraint(int l, int& where) {
    Lev& L = LEV[l];
    static std::vector<double> d;
    d.assign(L.f.size(), 0.0);
    ddx(L.f, L.n, L.h, d);
    double e = 0.0; where = 0;
    for (int i = 2; i < L.n - 2; ++i) {
        double c = std::fabs(L.S(i) - d[(size_t)(i + NG)]);
        if (c > e) { e = c; where = i; }
    }
    return e;
}
static long total_points() {
    long s = 0; for (auto& L : LEV) s += L.n; return s;
}
static long uniform_points() {                 // a flat grid at the finest spacing
    return (long)std::llround((XMAX - XMIN) / LEV.back().h) + 1;
}
static double work_amr() {                     // points x steps, per coarse step
    double w = 0.0;
    for (int l = 0; l < (int)LEV.size(); ++l) w += (double)LEV[l].n * (double)(1 << l);
    return w;
}
static double work_uniform() {
    return (double)uniform_points() * (double)(1 << (LEV.size() - 1));
}

extern "C" {

KEEP int sim_w() { return FW; }
KEEP int sim_h() { return FH; }
KEEP void sim_reset() { dirty = true; }
KEEP int sim_init(int, int) { px.assign((size_t)FW * FH, 0); init_all(); return 1; }
KEEP void sim_set(int id, double v) {
    if (id == 0) { int nl = (int)v; if (nl < 2) nl = 2; if (nl > 8) nl = 8;
                   if (nl != g_nlev) { g_nlev = nl; dirty = true; } }
    else if (id == 1) g_speed = v;
    else if (id == 2) { int r = (int)v; if (r < 1) r = 1; g_regrid = r; }
    else if (id == 3) { int b = (int)v; if (b < 4) b = 4; if (b != g_buf) { g_buf = b; dirty = true; } }
}
KEEP void sim_action(int a) {
    if (a == 0) { g_nlev = 2; dirty = true; }
    else if (a == 1) { g_nlev = 6; dirty = true; }
    else if (a == 2) { g_nlev = 8; dirty = true; }
    else if (a == 3) g_show_coarse = !g_show_coarse;
    else if (a == 4) dirty = true;
}
KEEP void sim_click(double, double) { }

static long g_coarse_steps = 0;
KEEP void sim_step(int steps) {
    if (dirty) { init_all(); g_coarse_steps = 0; }
    int want = (int)(g_speed * (steps > 0 ? steps : 1) + 0.5); if (want < 1) want = 1;
    for (int k = 0; k < want; ++k) {
        if (LEV[0].t > (XMAX - 6.0) - X0) break;           // pulse near the far edge
        advance(0);
        rk4_single(COARSE);
        ++g_coarse_steps;
        if (g_coarse_steps % g_regrid == 0) {
            double xc = track_center();
            for (int l = 1; l < (int)LEV.size(); ++l) {
                make_level(l, g_nlev, xc, true);
                nan_check(l, "after the regrid");
            }
        }
    }
}

// ------------------------------------------------------------------- rendering
static void plot_lev(Olivec_Canvas oc, const Lev& L, double ax0, double ax1,
                     int px0, int px1, int base, int hgt, uint32_t col, bool dots) {
    int lx = -1, ly = -1;
    for (int i = 0; i < L.n; ++i) {
        double xv = L.x(i);
        double fr = (xv - ax0) / (ax1 - ax0);
        if (fr < 0.0 || fr > 1.0) { lx = -1; continue; }
        int sx = px0 + (int)(fr * (px1 - px0));
        double a = L.F(i); if (a > 1.3) a = 1.3; if (a < -0.4) a = -0.4;
        int sy = base - (int)(a * hgt);
        if (lx >= 0) olivec_line(oc, lx, ly, sx, sy, col);
        if (dots && sx >= px0 && sx <= px1) OLIVEC_PIXEL(oc, sx, base + 3) = col;
        lx = sx; ly = sy;
    }
}

KEEP uint8_t* sim_render() {
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(3, 5, 11, 1.f));
    Olivec_Font ft = olivec_default_font;
    char buf[240];
    const uint32_t C_AX = rgba(38, 56, 68, 1.f);
    const uint32_t C_EX = rgba(190, 200, 215, 1.f);
    const uint32_t C_CO = rgba(255, 90, 100, 1.f);
    const int NL = (int)LEV.size();
    auto lcol = [&](int l) {
        float u = NL > 1 ? (float)l / (float)(NL - 1) : 0.f;
        return rgba((int)(30 + 60 * u), (int)(120 + 135 * u), (int)(220 - 40 * u), 1.f);
    };

    olivec_text(oc, "amr. how numerical relativity actually computes. the equation is only a wave,",
                12, 6, ft, 1, rgba(150, 230, 220, 1.f));
    olivec_text(oc, "the machinery is the real thing. berger-oliger levels, subcycling in time,",
                12, 18, ft, 1, rgba(150, 230, 220, 1.f));
    olivec_text(oc, "prolongation in space and time, restriction, nesting with buffers, regridding.",
                12, 30, ft, 1, rgba(150, 230, 220, 1.f));

    // ---------------- panel 1: the hierarchy over the whole domain ------------
    {
        const int TOP = 52, ROW = 15, px0 = 74, px1 = FW - 14;
        olivec_text(oc, "the boxes, following the pulse. each level halves h and halves its box,",
                    12, TOP - 6, ft, 1, rgba(120, 150, 160, 1.f));
        olivec_text(oc, "so every level costs about the same number of points.",
                    12, TOP + 5, ft, 1, rgba(120, 150, 160, 1.f));
        for (int l = 0; l < NL; ++l) {
            int y = TOP + 20 + l * ROW;
            Lev& L = LEV[l];
            double f0 = (L.x(0) - XMIN) / (XMAX - XMIN);
            double f1 = (L.x(L.n - 1) - XMIN) / (XMAX - XMIN);
            int a = px0 + (int)(f0 * (px1 - px0)), b = px0 + (int)(f1 * (px1 - px0));
            olivec_line(oc, px0, y, px1, y, C_AX);
            uint32_t c = lcol(l);
            olivec_line(oc, a, y, b, y, c);
            olivec_line(oc, a, y - 4, a, y + 4, c);
            olivec_line(oc, b, y - 4, b, y + 4, c);
            int step = (b - a) / 40; if (step < 1) step = 1;
            for (int q = a; q <= b; q += step) OLIVEC_PIXEL(oc, q, y - 5) = c;
            std::snprintf(buf, sizeof buf, "l %d", l);
            olivec_text(oc, buf, 14, y - 3, ft, 1, c);
            std::snprintf(buf, sizeof buf, "h %.5f  n %d  dt %.5f  steps %ld",
                          L.h, L.n, L.dt, L.steps);
            olivec_text(oc, buf, 42, y - 3, ft, 1, rgba(110, 130, 145, 1.f));
        }
    }

    // ---------------- panel 2: zoom on the pulse ------------------------------
    {
        const int TOP = 210, BOT = 366, px0 = 14, px1 = FW - 14;
        double xc = track_center();
        double zw = 8.0 * PW + 6.0 * H0 / (double)(1 << (NL - 1)) * 4.0;
        double ax0 = xc - zw, ax1 = xc + zw;
        olivec_text(oc, "zoomed on the pulse. dots under the axis are the finest level grid points.",
                    12, TOP - 10, ft, 1, rgba(120, 150, 160, 1.f));
        olivec_line(oc, px0, BOT, px1, BOT, C_AX);
        // exact
        int lx = -1, ly = -1;
        for (int q = px0; q <= px1; ++q) {
            double xv = ax0 + (ax1 - ax0) * (double)(q - px0) / (double)(px1 - px0);
            double a = ex_f(xv, LEV[0].t);
            int sy = BOT - (int)(a * 130.0);
            if (lx >= 0) olivec_line(oc, lx, ly, q, sy, C_EX);
            lx = q; ly = sy;
        }
        if (g_show_coarse) plot_lev(oc, COARSE, ax0, ax1, px0, px1, BOT, 130, C_CO, false);
        for (int l = 0; l < NL; ++l)
            plot_lev(oc, LEV[l], ax0, ax1, px0, px1, BOT, 130, lcol(l), l == NL - 1);
        olivec_text(oc, "white  exact.   blue shades  the levels.   red  one coarse grid alone.",
                    px0 + 4, TOP + 2, ft, 1, rgba(110, 130, 145, 1.f));
        std::snprintf(buf, sizeof buf, "window %.4f wide, centred on the tracked peak %.4f",
                      ax1 - ax0, xc);
        olivec_text(oc, buf, px0 + 4, TOP + 14, ft, 1, rgba(110, 130, 145, 1.f));
    }

    // ---------------- panel 3: who stepped when -------------------------------
    {
        const int TOP = 386, ROW = 9, px0 = 74, px1 = FW - 14;
        olivec_text(oc, "the subcycling, as it happened. one column per step taken.",
                    12, TOP - 10, ft, 1, rgba(120, 150, 160, 1.f));
        int shown = px1 - px0; if (shown > EVN) shown = EVN;
        for (int l = 0; l < NL; ++l) {
            int y = TOP + 4 + l * ROW;
            olivec_line(oc, px0, y, px1, y, rgba(22, 30, 38, 1.f));
            std::snprintf(buf, sizeof buf, "l %d", l);
            olivec_text(oc, buf, 40, y - 3, ft, 1, lcol(l));
        }
        for (int k = 0; k < shown; ++k) {
            int idx = (evhead - shown + k + 2 * EVN) % EVN;
            int l = ev[idx];
            if (l < 0 || l >= NL) continue;
            int y = TOP + 4 + l * ROW;
            olivec_line(oc, px0 + k, y - 3, px0 + k, y + 3, lcol(l));
        }
    }

    // ---------------- numbers -------------------------------------------------
    {
        int w0 = 0, wf = 0;
        double c0 = constraint(0, w0), cf = constraint(NL - 1, wf);
        long tp = total_points(), up = uniform_points();
        double wa = work_amr(), wu = work_uniform();
        double s3 = 0.0;
        for (int l = 0; l < NL; ++l)
            s3 += std::pow((double)LEV[l].n, 3.0) * (double)(1 << l);
        double u3 = std::pow((double)up, 3.0) * (double)(1 << (NL - 1));
        int Y = FH - 112;
        // olive's font has no plus sign, so keep exponents negative or spell the
        // number out -- otherwise 1.1e+05 silently renders as 1.1e 05
        std::snprintf(buf, sizeof buf, "t %.1f   levels %d   points %ld   flat would need %ld",
                      LEV[0].t, NL, tp, up);
        olivec_text(oc, buf, 12, Y, ft, 2, rgba(180, 210, 150, 1.f));
        std::snprintf(buf, sizeof buf, "work per coarse step  amr %.0f  flat %.0f  saved %.1f times",
                      wa, wu, wu / wa);
        olivec_text(oc, buf, 12, Y + 24, ft, 2, rgba(180, 210, 150, 1.f));
        std::snprintf(buf, sizeof buf, "the same hierarchy in 3d would save %.0f times", u3 / s3);
        olivec_text(oc, buf, 12, Y + 48, ft, 2, rgba(235, 190, 90, 1.f));
        double ea = err_amr();
        if (ea < 0.0)
            std::snprintf(buf, sizeof buf, "error  the boxes lost the pulse.  one coarse grid %.3e",
                          err_coarse());
        else
            std::snprintf(buf, sizeof buf, "error  amr %.3e   one coarse grid %.3e   better by %.0f",
                          ea, err_coarse(), err_coarse() / ea);
        olivec_text(oc, buf, 12, Y + 72, ft, 2, rgba(140, 210, 230, 1.f));
        std::snprintf(buf, sizeof buf, "constraint %.2f on level 0, %.2e on the finest, regrid %d",
                      c0, cf, g_regrid);
        olivec_text(oc, buf, 12, Y + 96, ft, 2, rgba(150, 170, 185, 1.f));
    }
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
int main(int argc, char** argv) {
    sim_init(0, 0);
    int frames = argc > 1 ? atoi(argv[1]) : 300;
    if (argc > 2) { PW = atof(argv[2]); dirty = true; }
    int want_nlev = (argc > 3) ? atoi(argv[3]) : 8;
    if (argc > 4) g_regrid = atoi(argv[4]);
    printf("berger-oliger amr on the first order wave equation\n");
    printf("domain %.1f to %.1f, coarse h %.4f, pulse width %.3f\n", XMIN, XMAX, H0, PW);
    for (int nl = 2; nl <= 8; ++nl) {
        g_nlev = nl; dirty = true; sim_step(0);
        long tp = total_points(), up = uniform_points();
        double wa = work_amr(), wu = work_uniform();
        double s3 = 0.0;
        for (int l = 0; l < nl; ++l) s3 += std::pow((double)LEV[l].n, 3.0) * (double)(1 << l);
        double u3 = std::pow((double)up, 3.0) * (double)(1 << (nl - 1));
        printf("  levels %d  finest h %.6f  points %5ld  flat would need %6ld"
               "   work saved 1d %6.1f x   3d %9.3g x\n",
               nl, LEV.back().h, tp, up, wu / wa, u3 / s3);
    }
    // The test that says the machinery is right: add a level, so the finest
    // spacing halves, and the error must fall by 16. Everything is on while this
    // runs -- subcycling, prolongation in space and time, restriction, and a
    // regrid on every coarse step.
    printf("\nconvergence with the whole machine running\n");
    double eprev = 0.0;
    for (int nl = 5; nl <= 8; ++nl) {
        g_nlev = nl; dirty = true; g_regrid = 1;
        sim_step(0);
        for (int i = 0; i < 60; ++i) sim_step(1);
        double e = err_amr();
        printf("  levels %d  finest h %.7f  error %.4e", nl, LEV.back().h, e);
        if (eprev > 0.0 && e > 0.0)
            printf("   fell by %.2f, order %.2f", eprev / e, std::log2(eprev / e));
        printf("%s\n", e < 0.0 ? "   (boxes lost the pulse, number meaningless)" : "");
        eprev = e;
    }
    g_nlev = want_nlev; if (argc > 4) g_regrid = atoi(argv[4]); dirty = true; sim_step(0);
    bool told = false;
    for (int i = 0; i < frames; ++i) {
        sim_step(1);
        for (int l = 0; l < (int)LEV.size() && !told; ++l) {
            for (int i = -NG; i < LEV[l].n + NG; ++i) {
                double v = LEV[l].F(i);
                if (!(v == v) || std::fabs(v) > 1e3) {
                    told = true;
                    printf("BLEW UP on level %d at i %d of %d, value %g, coarse step %ld, t %.4f\n",
                           l, i, LEV[l].n, v, LEV[l].steps, LEV[l].t);
                    printf("  that point is at x %.6f, the pulse is at %.6f\n",
                           LEV[l].x(i), X0 + LEV[0].t);
                    printf("  level %d box  x %.4f .. %.4f\n",
                           l, LEV[l].x(0), LEV[l].x(LEV[l].n - 1));
                    break;
                }
            }
        }
        double xc = track_center(), xe = X0 + LEV[0].t;
        if (!told && std::fabs(xc - xe) > 0.2) {
            told = true;
            printf("TRACKER LOST THE PULSE at coarse step %ld, t %.4f\n", LEV[0].steps, LEV[0].t);
            printf("  tracked %.4f but the peak is at %.4f\n", xc, xe);
            for (int l = 0; l < (int)LEV.size(); ++l)
                printf("  level %d  x %.4f .. %.4f  n %d\n",
                       l, LEV[l].x(0), LEV[l].x(LEV[l].n - 1), LEV[l].n);
        }
    }
    printf("\nafter %d coarse steps, t = %.6f\n", (int)LEV[0].steps, LEV[0].t);
    printf("steps taken per level:");
    for (int l = 0; l < (int)LEV.size(); ++l) printf("  l%d %ld", l, LEV[l].steps);
    printf("\n  ratio to level 0:");
    for (int l = 0; l < (int)LEV.size(); ++l)
        printf("  %.3f", (double)LEV[l].steps / (double)LEV[0].steps);
    printf("   (should be 1 2 4 8 ...)\n");
    printf("where the boxes actually are (exact peak is at x %.4f)\n", X0 + LEV[0].t);
    for (int l = 0; l < (int)LEV.size(); ++l)
        printf("  level %d  x %.4f .. %.4f   Ilo %ld  n %d\n",
               l, LEV[l].x(0), LEV[l].x(LEV[l].n - 1), LEV[l].Ilo, LEV[l].n);
    printf("  tracked peak %.4f, peak value on the finest level %.6f\n",
           track_center(), [] { double m = 0; for (int i = 0; i < LEV.back().n; ++i)
                                   m = std::fmax(m, std::fabs(LEV.back().F(i))); return m; }());
    printf("error against the closed form solution\n");
    printf("  amr finest level %.4e     one coarse grid alone %.4e     ratio %.4g\n",
           err_amr(), err_coarse(), err_coarse() / err_amr());
    printf("constraint s - dx f, level by level\n");
    for (int l = 0; l < (int)LEV.size(); ++l) {
        int w; double c = constraint(l, w);
        printf("  level %d  %.4e   at i %d of %d\n", l, c, w, LEV[l].n);
    }
    uint8_t* p = sim_render();
    stbi_write_png("amr_os_preview.png", FW, FH, 4, p, FW * 4);
    printf("wrote amr_os_preview.png\n");
    return 0;
}
#endif
