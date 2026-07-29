// Universe OS — Flat Grid / Photon Ring  (C++/WASM)
//
// THE POINT: the grid drawn on screen is EXACTLY a Cartesian lattice. It never
// bends. What varies is the local speed of light:
//
//     u = m/(2rho),   A^2 = ((1-u)/(1+u))^2,   B^4 = (1+u)^4
//     c(rho) = A/B^2 = (1-u)/(1+u)^3        <- goes to 0 at rho = m/2
//
// This is Schwarzschild in isotropic coordinates, read as "flat space with a
// position-dependent light speed". Light is traced with the canonical equations
// of the super-Hamiltonian
//
//     H = 1/2 ( -E^2 F(rho) + |p|^2 G(rho) ) = 0     F = 1/A^2, G = 1/B^4
//     dx_i/dl = p_i G,   dp_i/dl = -1/2 ( -E^2 F' + |p|^2 G' ) x_i/rho
//
// NO Christoffel symbols, NO curvature tensor, NO geodesic equation.
//
// The turning point of a ray is where rho/c(rho) = b, so the shadow radius is
//
//     b_crit = min over rho of  rho/c(rho)  =  sqrt(27) m = 5.196152423 m
//
// and as b -> b_crit the deflection diverges logarithmically, producing PHOTON
// RINGS at deflection = pi, 3pi, 5pi, ...  Their distance from the shadow edge
// shrinks by exactly e^(2pi) = 535.4917 per turn -- which is why only the n=1
// ring is ever visible on a linear image (the n=2 ring is sub-pixel).
//
// Verified against pure-Python reference (miharashi/flatgrid2.py):
//   b_crit  = 5.196152423 m   (analytic sqrt(27), 12 digits)
//   rho_ph  = 1.866025 m      (= areal radius 3m)
//   rings at b/b_crit-1 = 3.094717e-2, 5.409356e-5, 1.009951e-7, 1.885839e-10
//   Sun's limb deflection 1.751201" vs VLBI 1.75119"
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

static const int FW = 900, FH = 620;
static std::vector<uint32_t> px;

static inline uint32_t rgba(int r, int g, int b, float a) {
    int A = (int)(a * 255.0f); if (A < 0) A = 0; if (A > 255) A = 255;
    return ((uint32_t)A << 24) | ((uint32_t)b << 16) | ((uint32_t)g << 8) | (uint32_t)r;
}

// ------------------------------------------------------------ the two functions
// m = 1 throughout. u = 1/(2 rho).
static inline void FG(double rho, double& F, double& G, double& dF, double& dG) {
    double u = 0.5 / rho, du = -u / rho;
    double om = 1.0 - u, op = 1.0 + u;
    F = (op * op) / (om * om);
    G = 1.0 / (op * op * op * op);
    double dF_du = 4.0 * op / (om * om * om);      // d/du [ (1+u)^2/(1-u)^2 ]
    double dG_du = -4.0 / (op * op * op * op * op);
    dF = dF_du * du; dG = dG_du * du;
}
static inline double c_local(double rho) {
    double F, G, dF, dG; FG(rho, F, G, dF, dG);
    return std::sqrt(G / F);                        // = |1-u|/(1+u)^3
}
static inline double Vpot(double rho) { return rho / c_local(rho); }

static double B_CRIT = 0.0, RHO_PH = 0.0;
static void find_bcrit() {
    double lo = 0.5000001, hi = 60.0;
    for (int i = 0; i < 300; ++i) {
        double a = lo + (hi - lo) / 3.0, b = hi - (hi - lo) / 3.0;
        if (Vpot(a) < Vpot(b)) hi = b; else lo = a;
    }
    RHO_PH = 0.5 * (lo + hi);
    B_CRIT = Vpot(RHO_PH);
}



// olivec_circle FILLS. We need outlines, so draw them with the midpoint
// algorithm. (Using the filled version turned the shadow into a solid disc and
// the outermost equal-speed ring painted the whole frame.)
static void ring(Olivec_Canvas oc, int cx, int cy, int r, uint32_t col) {
    if (r < 1) return;
    int x = r, y = 0, err = 1 - r;
    auto put = [&](int a, int b) {
        if (a >= 0 && a < (int)oc.width && b >= 0 && b < (int)oc.height)
            OLIVEC_PIXEL(oc, a, b) = col;
    };
    while (x >= y) {
        put(cx + x, cy + y); put(cx + y, cy + x);
        put(cx - y, cy + x); put(cx - x, cy + y);
        put(cx - x, cy - y); put(cx - y, cy - x);
        put(cx + y, cy - x); put(cx + x, cy - y);
        ++y;
        if (err < 0) err += 2 * y + 1;
        else { --x; err += 2 * (y - x) + 1; }
    }
}

// ---------------------------------------------------- exact deflection (1-D integral)
// Turning point: V(rho) = rho/c(rho) = b (outer root).  Then
//   dphi/drho = (b/rho^2) / sqrt( 1/c^2 - b^2/rho^2 )
// Substituting rho = rho0/w, w = 1-v^2 removes the 1/sqrt singularity at the
// turning point.  Both endpoints need their analytic limits (v=0 is 0/0; putting
// 0 there makes the error fall as 1/N instead of 1/N^4 -- that bug ate the
// weak-field 4m/b entirely when this was first written in Python).
static double turning_point(double b) {
    double lo = RHO_PH, hi = (2.0 * b > 4.0 ? 2.0 * b : 4.0);
    while (Vpot(hi) < b) hi *= 2.0;
    for (int i = 0; i < 200; ++i) {
        double mid = 0.5 * (lo + hi);
        if (Vpot(mid) < b) lo = mid; else hi = mid;
        if (hi - lo < 1e-16 * hi) break;
    }
    return 0.5 * (lo + hi);
}
static double deflection_exact(double b, int N = 4000) {
    if (b <= B_CRIT) return -1.0;                       // captured: no turning point
    double r0 = turning_point(b);
    const double dd = 1e-7, w = 1.0 - dd;
    double cc = c_local(r0 / w);
    double Q = 1.0 / (cc * cc) - b * b * w * w / (r0 * r0);
    double f0 = (Q > 0.0) ? 2.0 * b / (r0 * std::sqrt(Q / dd)) : 0.0;
    if (N % 2) ++N;
    double h = 1.0 / N, sum = 0.0;
    for (int i = 0; i <= N; ++i) {
        double v = i * h, val;
        if (i == 0) val = f0;
        else {
            double ww = 1.0 - v * v;
            if (ww <= 0.0) val = 2.0 * v * b / r0;       // rho = infinity, Q -> 1
            else {
                double rho = r0 / ww, c2 = c_local(rho);
                double QQ = 1.0 / (c2 * c2) - b * b / (rho * rho);
                val = (QQ > 0.0) ? 2.0 * v * b / (r0 * std::sqrt(QQ)) : 0.0;
            }
        }
        double wt = (i == 0 || i == N) ? 1.0 : ((i & 1) ? 4.0 : 2.0);
        sum += wt * val;
    }
    return 2.0 * sum * h / 3.0 - M_PI;
}

// ------------------------------------------------------------ ray tracing
struct Pt { float x, y; };
static std::vector<Pt> path;
static bool captured = false;
static double winding = 0.0, defl = 0.0;

static void rhs(const double s[4], double d[4]) {
    double x = s[0], y = s[1], pxx = s[2], pyy = s[3];
    double rho = std::sqrt(x * x + y * y);
    double F, G, dF, dG; FG(rho, F, G, dF, dG);
    double p2 = pxx * pxx + pyy * pyy;
    double k = -0.5 * (-dF + p2 * dG) / rho;        // E = 1
    d[0] = pxx * G; d[1] = pyy * G; d[2] = k * x; d[3] = k * y;
}

// trace a ray coming from x=-X0 at height b; store the polyline
static void trace(double b, double X0) {
    path.clear(); captured = false; winding = 0.0; defl = 0.0;
    // The impact parameter that matters is b = L/E, and L = y0 * p_x with
    // p_x = sqrt(F/G) at the launch radius. Launching at y0 = b makes L = b*sqrt(F/G),
    // which at rho = 200 is 1% off -- larger than the ring-2 offset of 5.4e-5, so the
    // ray misses the photon sphere entirely. Launch at y0 = b*sqrt(G/F) instead, which
    // makes L = b exactly.
    double y0 = b, F, G, dF, dG, rho = 0.0;
    for (int it = 0; it < 4; ++it) {
        rho = std::sqrt(X0 * X0 + y0 * y0);
        FG(rho, F, G, dF, dG);
        y0 = b * std::sqrt(G / F);
    }
    rho = std::sqrt(X0 * X0 + y0 * y0);
    FG(rho, F, G, dF, dG);
    double s[4] = { -X0, y0, std::sqrt(F / G), 0.0 };
    double th_prev = std::atan2(s[1], s[0]), th = 0.0;
    path.push_back({ (float)s[0], (float)s[1] });
    const int MAXN = 400000;
    const double ETA = 0.0015;
    for (int n = 0; n < MAXN; ++n) {
        rho = std::sqrt(s[0] * s[0] + s[1] * s[1]);
        if (rho < 0.52) { captured = true; break; }             // fell to the throat
        // Escape means "far away and heading outward". Testing x > X0 is wrong:
        // a ray deflected by pi goes back the way it came, so x never exceeds X0
        // and the loop just spins until the step cap.
        if (rho > X0 && (s[0] * s[2] + s[1] * s[3]) > 0.0) break;
        // step must shrink toward the throat: c -> 0 there and dF diverges,
        // so a fixed step overshoots past rho = m/2 and the state explodes.
        double d0 = rho - 0.5;
        double h = ETA * (d0 < 1.0 ? (d0 > 0.02 ? d0 : 0.02) : rho);
        double a[4], bb[4], cc[4], dd[4], t[4];
        rhs(s, a);
        for (int i = 0; i < 4; ++i) t[i] = s[i] + 0.5 * h * a[i]; rhs(t, bb);
        for (int i = 0; i < 4; ++i) t[i] = s[i] + 0.5 * h * bb[i]; rhs(t, cc);
        for (int i = 0; i < 4; ++i) t[i] = s[i] + h * cc[i]; rhs(t, dd);
        for (int i = 0; i < 4; ++i) s[i] += h / 6.0 * (a[i] + 2 * bb[i] + 2 * cc[i] + dd[i]);
        if (!(s[0] == s[0]) || !(s[1] == s[1]) || !(s[2] == s[2])) {
            captured = true; break;      // catch nan BEFORE it pollutes theta
        }
        // accumulate the azimuth CONTINUOUSLY (folding it into +-pi is the bug
        // that produced -5.4e8 arcsec/century in the Mercury run)
        double raw = std::atan2(s[1], s[0]);
        double d = raw - th_prev;
        while (d > M_PI) d -= 2.0 * M_PI;
        while (d < -M_PI) d += 2.0 * M_PI;
        th += d; th_prev = raw;
        if ((n & 3) == 0) path.push_back({ (float)s[0], (float)s[1] });
        // The photon sphere is an UNSTABLE orbit, so a ray that is numerically
        // a hair off critical can get stuck circling it for thousands of turns
        // until the state overflows to nan. Four turns is plenty to look at.
        if (std::fabs(th) > 8.0 * M_PI) break;
    }
    path.push_back({ (float)s[0], (float)s[1] });
    winding = std::fabs(th) / (2.0 * M_PI);
    defl = std::fabs(th) - M_PI;                    // total deflection
}

// ------------------------------------------------------------ state
static double g_view = 14.0;       // half-width of the view, in units of m
static double g_leps = -1.5;       // log10( b/b_crit - 1 )
static int g_marker = 0;
static bool dirty = true;
static bool g_fan = true;

static const double RING_EPS[4] = { 3.094717e-2, 5.409356e-5, 1.009951e-7, 1.885839e-10 };

extern "C" {

KEEP int sim_w() { return FW; }
KEEP int sim_h() { return FH; }

KEEP void sim_reset() { g_view = 14.0; g_leps = -1.5; g_marker = 0; g_fan = true; dirty = true; }

KEEP int sim_init(int, int) {
    px.assign((size_t)FW * FH, 0);
    find_bcrit();
    sim_reset();
    return 1;
}

KEEP void sim_set(int id, double v) {
    if (id == 0) { g_leps = v; dirty = true; }
    else if (id == 1) { g_view = v; }
    else if (id == 2) { g_fan = (v > 0.5); }
}

KEEP void sim_action(int a) {
    if (a >= 0 && a <= 3) { g_leps = std::log10(RING_EPS[a]); dirty = true; g_marker = 0; }
    else if (a == 4) { g_leps = -12.0; dirty = true; g_marker = 0; }   // essentially the edge
    else if (a == 5) { g_fan = !g_fan; }
}

// clicking picks b from the vertical distance to the axis
KEEP void sim_click(double nx, double ny) {
    (void)nx;
    double y = (0.5 - ny) * 2.0 * g_view * (double)FH / (double)FW;
    double b = std::fabs(y);
    double e = b / B_CRIT - 1.0;
    if (e < 1e-12) e = 1e-12;
    g_leps = std::log10(e); if (g_leps > 0.7) g_leps = 0.7;
    dirty = true; g_marker = 0;
}

KEEP void sim_step(int steps) {
    if (dirty) {
        double b = B_CRIT * (1.0 + std::pow(10.0, g_leps));
        trace(b, 6.0 * g_view + 40.0);
        dirty = false; g_marker = 0;
    }
    g_marker += 6 * (steps > 0 ? steps : 1);
    if (path.size() && g_marker >= (int)path.size()) g_marker = 0;
}

// ------------------------------------------------------------ render
static double SC;                      // pixels per m
static inline int SX(double x) { return (int)(FW * 0.5 + x * SC); }
static inline int SY(double y) { return (int)(FH * 0.5 - y * SC); }

static void draw_ray(Olivec_Canvas oc, double b, uint32_t col, bool thick) {
    // trace into a scratch buffer without disturbing the main path
    std::vector<Pt> keep; keep.swap(path);
    bool kc = captured; double kw = winding, kd = defl;
    trace(b, 6.0 * g_view + 40.0);
    for (size_t i = 1; i < path.size(); ++i) {
        olivec_line(oc, SX(path[i - 1].x), SY(path[i - 1].y), SX(path[i].x), SY(path[i].y), col);
        if (thick) olivec_line(oc, SX(path[i - 1].x), SY(path[i - 1].y) + 1,
                               SX(path[i].x), SY(path[i].y) + 1, col);
    }
    path.swap(keep); captured = kc; winding = kw; defl = kd;
}

KEEP uint8_t* sim_render() {
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(4, 5, 10, 1.f));
    SC = (FW * 0.5) / g_view;

    // ---- the grid. It is EXACTLY Cartesian. Nothing here bends. ----
    uint32_t gcol = rgba(0, 70, 60, 1.f);
    for (double gx = -60; gx <= 60; gx += 1.0) {
        int sx = SX(gx); if (sx < 0 || sx >= FW) continue;
        olivec_line(oc, sx, 0, sx, FH - 1, gcol);
    }
    for (double gy = -60; gy <= 60; gy += 1.0) {
        int sy = SY(gy); if (sy < 0 || sy >= FH) continue;
        olivec_line(oc, 0, sy, FW - 1, sy, gcol);
    }

    // ---- rings where the local light speed has dropped to a round fraction ----
    uint32_t ccol = rgba(0, 120, 150, 1.f);
    for (int k = 1; k <= 9; ++k) {
        double frac = 0.1 * k;                       // c/c_inf = frac
        // solve |1-u|/(1+u)^3 = frac for u in (0,1)
        double lo = 0.0, hi = 1.0;
        for (int i = 0; i < 60; ++i) {
            double mid = 0.5 * (lo + hi);
            double val = (1.0 - mid) / std::pow(1.0 + mid, 3);
            if (val > frac) lo = mid; else hi = mid;
        }
        double u = 0.5 * (lo + hi), rho = 0.5 / u;
        int r = (int)(rho * SC);
        if (r > 3 && r < FW) ring(oc, FW / 2, FH / 2, r, rgba(0, 62, 78, 1.f));
        (void)ccol;
    }

    // ---- shadow (b_crit), photon sphere, throat ----
    ring(oc, FW / 2, FH / 2, (int)(B_CRIT * SC), rgba(150, 80, 190, 1.f));
    ring(oc, FW / 2, FH / 2, (int)(RHO_PH * SC), rgba(220, 170, 50, 1.f));
    olivec_circle(oc, FW / 2, FH / 2, (int)(0.5 * SC) + 2, rgba(0, 0, 0, 1.f));
    ring(oc, FW / 2, FH / 2, (int)(0.5 * SC) + 2, rgba(120, 60, 150, 1.f));

    // ---- a fan of ordinary rays, to show plain bending ----
    if (g_fan) {
        for (int k = 1; k <= 6; ++k) {
            double b = B_CRIT * (1.0 + 0.35 * k);
            draw_ray(oc, b, rgba(20, 130, 115, 1.f), false);
            draw_ray(oc, -b, rgba(20, 130, 115, 1.f), false);
        }
    }

    // ---- the selected ray ----
    uint32_t rc = captured ? rgba(255, 70, 90, 1.f) : rgba(0, 255, 204, 1.f);
    for (size_t i = 1; i < path.size(); ++i) {
        olivec_line(oc, SX(path[i - 1].x), SY(path[i - 1].y), SX(path[i].x), SY(path[i].y), rc);
        olivec_line(oc, SX(path[i - 1].x) + 1, SY(path[i - 1].y), SX(path[i].x) + 1, SY(path[i].y), rc);
    }
    // travelling photon
    if (path.size() > 2) {
        int i = g_marker % (int)path.size();
        olivec_circle(oc, SX(path[i].x), SY(path[i].y), 4, rgba(255, 255, 200, 1.f));
    }

    // ---- HUD ----
    char buf[160];
    Olivec_Font f = olivec_default_font;
    uint32_t tc = rgba(120, 220, 200, 1.f), tw = rgba(255, 210, 120, 1.f);
    olivec_text(oc, "flat grid - the lattice never bends. c of rho does.", 14, 12, f, 2, tc);
    std::snprintf(buf, sizeof buf, "b crit  min rho over c  %.9f m    sqrt27 %.9f",
                  B_CRIT, std::sqrt(27.0));
    olivec_text(oc, buf, 14, 36, f, 2, tc);
    std::snprintf(buf, sizeof buf, "b over b crit minus 1   1e%+.2f", g_leps);
    olivec_text(oc, buf, 14, 64, f, 2, tw);
    if (captured) {
        olivec_text(oc, "captured. inside the shadow.", 14, 92, f, 2, rgba(255, 90, 110, 1.f));
    } else {
        std::snprintf(buf, sizeof buf, "deflection %.4f rad    winding %.3f turns",
                      defl, winding);
        olivec_text(oc, buf, 14, 92, f, 2, tw);
    }
    std::snprintf(buf, sizeof buf, "photon sphere rho %.6f m, areal 3m.  throat rho 0.5 m where c goes to 0",
                  RHO_PH);
    olivec_text(oc, buf, 14, FH - 26, f, 2, rgba(150, 130, 90, 1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
int main(int argc, char** argv) {
    sim_init(0, 0);
    printf("b_crit  = %.12f   sqrt(27) = %.12f   ratio = %.12f\n",
           B_CRIT, std::sqrt(27.0), B_CRIT / std::sqrt(27.0));
    printf("rho_ph  = %.9f   1+sqrt3/2 = %.9f\n", RHO_PH, 1.0 + std::sqrt(3.0) / 2.0);
    printf("areal r = %.9f   (should be 3)\n", RHO_PH * std::pow(1.0 + 0.5 / RHO_PH, 2));
    // deflection at each ring b: should be pi, 3pi, 5pi, 7pi  (exact 1-D integral)
    for (int k = 0; k < 4; ++k) {
        double b = B_CRIT * (1.0 + RING_EPS[k]);
        double de = deflection_exact(b, 40000);
        printf("ring %d: b/b_c-1=%.4e  exact=%9.6f  target=%9.6f  ratio=%.6f\n",
               k + 1, RING_EPS[k], de, (2 * k + 1) * M_PI, de / ((2 * k + 1) * M_PI));
    }
    // weak field: 4m/b
    for (double b : {1e6, 1e4, 1e2, 30.0}) {
        double de = deflection_exact(b, 20000);
        printf("weak b=%10.0f m: exact=%.6e   4m/b=%.6e   ratio=%.6f\n",
               b, de, 4.0 / b, de / (4.0 / b));
    }
    // the tracer is for DRAWING only; report how far it gets
    for (int k = 0; k < 3; ++k) {
        trace(B_CRIT * (1.0 + RING_EPS[k]), 200.0);
        printf("tracer ring %d: turns=%.4f captured=%d points=%d\n",
               k + 1, winding, (int)captured, (int)path.size());
    }
    int steps = argc > 1 ? atoi(argv[1]) : 1;
    g_leps = std::log10(RING_EPS[1]); dirty = true;
    sim_step(steps);
    uint8_t* p = sim_render();
    stbi_write_png("flatgrid_os_preview.png", FW, FH, 4, p, FW * 4);
    printf("wrote flatgrid_os_preview.png\n");
    return 0;
}
#endif
