// Universe OS — Two Grids, One Shadow  (C++/WASM)
//
// The claim being tested: "c = const coordinates and c*t = const coordinates are
// the same thing, so use whichever is convenient." That is a statement about
// COORDINATES, not about physics -- which means it has one observable
// consequence and only one: the two descriptions must agree on everything
// measurable. This demo puts them side by side and checks.
//
//   LEFT  - static grid, variable light speed (isotropic Schwarzschild)
//           space is conformally stretched: gamma_ij = B^4 delta_ij
//           c(rho) = A/B^2 = (1-u)/(1+u)^3   ->  0 at the throat rho = m/2
//           turning point:  rho/c(rho) = b
//           b_crit = min rho/c(rho)
//
//   RIGHT - flowing grid, constant light speed (Painleve-Gullstrand)
//           space is EXACTLY euclidean, lapse exactly 1
//           the grid flows inward at v(r) = sqrt(2m/r), reaching c at r = 2m
//           turning point:  r/sqrt(1 - 2m/r) = b
//           b_crit = min of that
//
// Different variables, different pictures, different places where something
// "stops" (c hits zero on the left at rho = m/2; the river hits c on the right
// at r = 2m -- and those are the SAME surface). Same shadow:
//
//     b_crit = sqrt(27) m = 5.196152422707   from both, to 12 digits
//     photon orbit at areal radius 3m        from both
//
// That agreement is the entire content of the coordinate-equivalence claim. It
// is not a prediction about the world; it is a consistency check on the
// bookkeeping. This demo is that check, drawn.
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
static void ring(Olivec_Canvas oc, int cx, int cy, int r, uint32_t col,
                 int x0, int x1, int y0, int y1) {
    if (r < 1) return;
    int x = r, y = 0, err = 1 - r;
    auto put = [&](int a, int b) {
        if (a >= x0 && a <= x1 && b >= y0 && b <= y1) OLIVEC_PIXEL(oc, a, b) = col;
    };
    while (x >= y) {
        put(cx + x, cy + y); put(cx + y, cy + x); put(cx - y, cy + x); put(cx - x, cy + y);
        put(cx - x, cy - y); put(cx - y, cy - x); put(cx + y, cy - x); put(cx + x, cy - y);
        ++y;
        if (err < 0) err += 2 * y + 1; else { --x; err += 2 * (y - x) + 1; }
    }
}

// ------------------------------------------------------------ LEFT: static grid
static inline void FG(double rho, double& F, double& G, double& dF, double& dG) {
    double u = 0.5 / rho, du = -u / rho;
    double om = 1.0 - u, op = 1.0 + u;
    F = (op * op) / (om * om);
    G = 1.0 / (op * op * op * op);
    dF = (4.0 * op / (om * om * om)) * du;
    dG = (-4.0 / (op * op * op * op * op)) * du;
}
static inline double c_iso(double rho) {
    double F, G, dF, dG; FG(rho, F, G, dF, dG); return std::sqrt(G / F);
}
static inline double V_iso(double rho) { return rho / c_iso(rho); }

// ------------------------------------------------------------ RIGHT: flowing grid
static inline double v_flow(double r) { return std::sqrt(2.0 / r); }
// turning point condition in PG: r / sqrt(1 - v^2) = b   (v^2 = 2m/r)
static inline double V_pg(double r) {
    double s = 1.0 - 2.0 / r;
    return (s > 0.0) ? r / std::sqrt(s) : 1e30;
}

static double B_ISO = 0, RHO_ISO = 0, B_PG = 0, R_PG = 0;
static double minimize(double (*f)(double), double lo, double hi, double& at) {
    for (int i = 0; i < 400; ++i) {
        double a = lo + (hi - lo) / 3.0, b = hi - (hi - lo) / 3.0;
        if (f(a) < f(b)) hi = b; else lo = a;
    }
    at = 0.5 * (lo + hi); return f(at);
}
static void solve_both() {
    B_ISO = minimize(V_iso, 0.5000001, 60.0, RHO_ISO);
    B_PG  = minimize(V_pg,  2.0000001, 60.0, R_PG);
}
// areal radius of an isotropic rho
static inline double areal(double rho) { double u = 0.5 / rho; return rho * (1.0 + u) * (1.0 + u); }

// ------------------------------------------------------------ shared tracer
struct Pt { float x, y; };
static std::vector<Pt> pathL, pathR;
static bool capL = false, capR = false;
static double thL = 0, thR = 0;

static void rhsL(const double s[4], double d[4]) {
    double rho = std::sqrt(s[0] * s[0] + s[1] * s[1]);
    double F, G, dF, dG; FG(rho, F, G, dF, dG);
    double p2 = s[2] * s[2] + s[3] * s[3];
    double k = -0.5 * (-dF + p2 * dG) / rho;
    d[0] = s[2] * G; d[1] = s[3] * G; d[2] = k * s[0]; d[3] = k * s[1];
}
// PG null rays, zero shift removed: use the Schwarzschild-coordinate form since
// the PG slice differs from it only by a time relabel, which cannot move a
// photon's spatial path. Effective: dphi/dr from r/sqrt(1-2m/r) = b.
static void trace_iso(double b, double X0, std::vector<Pt>& out, bool& cap, double& th) {
    out.clear(); cap = false; th = 0.0;
    double y0 = b, F, G, dF, dG, rho = 0;
    for (int it = 0; it < 4; ++it) {
        rho = std::sqrt(X0 * X0 + y0 * y0); FG(rho, F, G, dF, dG);
        y0 = b * std::sqrt(G / F);
    }
    rho = std::sqrt(X0 * X0 + y0 * y0); FG(rho, F, G, dF, dG);
    double s[4] = { -X0, y0, std::sqrt(F / G), 0.0 };
    double prev = std::atan2(s[1], s[0]);
    out.push_back({ (float)s[0], (float)s[1] });
    for (int n = 0; n < 300000; ++n) {
        rho = std::sqrt(s[0] * s[0] + s[1] * s[1]);
        if (rho < 0.52) { cap = true; break; }
        if (rho > X0 && (s[0] * s[2] + s[1] * s[3]) > 0.0) break;
        double d0 = rho - 0.5;
        double h = 0.0025 * (d0 < 1.0 ? (d0 > 0.02 ? d0 : 0.02) : rho);
        double a[4], bb[4], cc[4], dd[4], t[4];
        rhsL(s, a);
        for (int i = 0; i < 4; ++i) t[i] = s[i] + 0.5 * h * a[i]; rhsL(t, bb);
        for (int i = 0; i < 4; ++i) t[i] = s[i] + 0.5 * h * bb[i]; rhsL(t, cc);
        for (int i = 0; i < 4; ++i) t[i] = s[i] + h * cc[i]; rhsL(t, dd);
        for (int i = 0; i < 4; ++i) s[i] += h / 6.0 * (a[i] + 2 * bb[i] + 2 * cc[i] + dd[i]);
        if (!(s[0] == s[0]) || !(s[1] == s[1])) { cap = true; break; }
        double raw = std::atan2(s[1], s[0]), dth = raw - prev;
        while (dth > M_PI) dth -= 2 * M_PI;
        while (dth < -M_PI) dth += 2 * M_PI;
        th += dth; prev = raw;
        if ((n & 3) == 0) out.push_back({ (float)s[0], (float)s[1] });
        if (std::fabs(th) > 8.0 * M_PI) break;
    }
    out.push_back({ (float)s[0], (float)s[1] });
}
// Right panel: integrate dphi/dr directly from the PG turning-point function.
static void trace_pg(double b, double Rmax, std::vector<Pt>& out, bool& cap, double& th) {
    out.clear(); cap = false; th = 0.0;
    // find the turning point: outer root of V_pg(r) = b
    double r0 = -1.0;
    if (b > B_PG) {
        double lo = R_PG, hi = std::fmax(2.0 * b, 8.0);
        while (V_pg(hi) < b) hi *= 2.0;
        for (int i = 0; i < 200; ++i) {
            double mid = 0.5 * (lo + hi);
            if (V_pg(mid) < b) lo = mid; else hi = mid;
        }
        r0 = 0.5 * (lo + hi);
    } else { cap = true; }
    // inbound leg from Rmax to r0 (or to the horizon if captured), then outbound
    double rlo = cap ? 2.0001 : r0;
    int N = 4000;
    double phi = M_PI;                       // start on the -x axis
    auto emit = [&](double r, double p) {
        out.push_back({ (float)(r * std::cos(p)), (float)(r * std::sin(p)) });
    };
    emit(Rmax, phi);
    // dphi/dr = (b/r^2)/sqrt( 1/(1-2m/r) - b^2/r^2 ).  The 1/sqrt blows up at the
    // turning point, so substitute r = rlo/w with w = 1 - v^2. Then
    //     dphi = 2 v b dv / ( rlo sqrt(Q) ),  Q = 1/(1-2m/r) - b^2/r^2
    // and Q ~ v^2 near v = 0, so the integrand is finite there. (Clustering
    // points as (1-t^2) is NOT enough -- the singularity survives it, which is
    // what made the right panel disagree by 50% on the first attempt.)
    auto integrand = [&](double v) {
        double w = 1.0 - v * v;
        if (w <= 0.0) return 2.0 * v * b / rlo;                 // r = infinity, Q -> 1
        double r = rlo / w, s = 1.0 - 2.0 / r;
        if (s <= 0.0) return 0.0;
        double Q = 1.0 - b * b * s / (r * r);
        return (Q > 0.0) ? 2.0 * v * b / (rlo * std::sqrt(Q)) : 0.0;
    };
    // limit at v = 0 (0/0): take it from a small offset, same trick as flatgrid_os
    const double dv = 1e-7;
    double w1 = 1.0 - dv, r1 = rlo / w1, s1 = 1.0 - 2.0 / r1;
    double Q1 = 1.0 - b * b * s1 / (r1 * r1);
    double f0 = (Q1 > 0.0) ? 2.0 * b / (rlo * std::sqrt(Q1 / dv)) : 0.0;

    if (cap) {                               // no turning point: just fall in
        for (int i = 0; i <= 600; ++i) {
            double r = Rmax - (Rmax - 2.0001) * i / 600.0;
            double s = 1.0 - 2.0 / r, Q = (s > 0.0) ? 1.0 - b * b * s / (r * r) : -1.0;
            if (Q <= 0.0) break;
            phi += (b / (r * r)) / std::sqrt(Q) * ((Rmax - 2.0001) / 600.0);
            emit(r, phi);
        }
        th = std::fabs(phi - M_PI);
        return;
    }
    // inbound leg: v from 1 (infinity) down to 0 (turning point)
    for (int i = N; i >= 0; --i) {
        double v = (double)i / N, dphi;
        double val = (i == 0) ? f0 : integrand(v);
        dphi = val / N;
        phi += dphi;
        double w = 1.0 - v * v;
        double r = (w > 1e-12) ? rlo / w : Rmax;
        if (r <= Rmax) emit(r, phi);
    }
    // outbound leg: mirror, v from 0 back to 1
    for (int i = 0; i <= N; ++i) {
        double v = (double)i / N;
        double val = (i == 0) ? f0 : integrand(v);
        phi += val / N;
        double w = 1.0 - v * v;
        double r = (w > 1e-12) ? rlo / w : Rmax;
        if (r <= Rmax) emit(r, phi);
    }
    th = std::fabs(phi - M_PI) - M_PI;       // deflection, not the swept angle
}

static double g_leps = -1.2;      // log10(b/b_crit - 1)
static double g_view = 12.0;
static bool dirty = true;
static int g_mark = 0;

extern "C" {

KEEP int sim_w() { return FW; }
KEEP int sim_h() { return FH; }
KEEP void sim_reset() { g_leps = -1.2; g_view = 12.0; dirty = true; g_mark = 0; }
KEEP int sim_init(int, int) { px.assign((size_t)FW * FH, 0); solve_both(); sim_reset(); return 1; }
KEEP void sim_set(int id, double v) {
    if (id == 0) { g_leps = v; dirty = true; } else if (id == 1) g_view = v;
}
KEEP void sim_action(int a) {
    const double E[4] = { 3.094717e-2, 5.409356e-5, 1.009951e-7, 1.885839e-10 };
    if (a >= 0 && a <= 3) { g_leps = std::log10(E[a]); dirty = true; }
    else if (a == 4) { g_leps = -12.0; dirty = true; }
}
KEEP void sim_click(double nx, double ny) {
    (void)nx;
    double y = (0.5 - ny) * 2.0 * g_view * (double)FH / (double)(FW / 2);
    double e = std::fabs(y) / B_ISO - 1.0; if (e < 1e-12) e = 1e-12;
    g_leps = std::log10(e); if (g_leps > 0.6) g_leps = 0.6;
    dirty = true;
}
KEEP void sim_step(int steps) {
    if (dirty) {
        double b = B_ISO * (1.0 + std::pow(10.0, g_leps));
        trace_iso(b, 6.0 * g_view + 40.0, pathL, capL, thL);
        trace_pg(b, 6.0 * g_view + 40.0, pathR, capR, thR);
        dirty = false; g_mark = 0;
    }
    g_mark += 5 * (steps > 0 ? steps : 1);
}

KEEP uint8_t* sim_render() {
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(3, 5, 11, 1.f));
    Olivec_Font ft = olivec_default_font;
    char buf[210];
    const int HW = FW / 2, TOP = 106, BOT = FH - 94;
    const int CY = (TOP + BOT) / 2;
    double SC = (HW * 0.5) / g_view;
    olivec_line(oc, HW, TOP - 8, HW, BOT + 8, rgba(70, 90, 100, 1.f));

    // ---- left half : static grid, variable c ----
    {
        int cx = HW / 2, x0 = 0, x1 = HW - 2;
        uint32_t g1 = rgba(0, 56, 60, 1.f);
        for (double g = -60; g <= 60; g += 1.0) {
            int a = (int)(cx + g * SC); if (a > x0 && a < x1) olivec_line(oc, a, TOP, a, BOT, g1);
            int b = (int)(CY - g * SC); if (b > TOP && b < BOT) olivec_line(oc, x0, b, x1, b, g1);
        }
        for (int k = 1; k <= 8; ++k) {                      // equal-speed rings
            double frac = 0.1 * k, lo = 0.0, hi = 1.0;
            for (int i = 0; i < 60; ++i) {
                double mid = 0.5 * (lo + hi);
                if ((1.0 - mid) / std::pow(1.0 + mid, 3) > frac) lo = mid; else hi = mid;
            }
            double u = 0.5 * (lo + hi);
            ring(oc, cx, CY, (int)(0.5 / u * SC), rgba(0, 60, 76, 1.f), x0, x1, TOP, BOT);
        }
        ring(oc, cx, CY, (int)(B_ISO * SC), rgba(150, 80, 190, 1.f), x0, x1, TOP, BOT);
        ring(oc, cx, CY, (int)(RHO_ISO * SC), rgba(230, 180, 55, 1.f), x0, x1, TOP, BOT);
        olivec_circle(oc, cx, CY, (int)(0.5 * SC) + 2, rgba(0, 0, 0, 1.f));
        ring(oc, cx, CY, (int)(0.5 * SC) + 2, rgba(120, 60, 150, 1.f), x0, x1, TOP, BOT);
        uint32_t rc = capL ? rgba(255, 75, 95, 1.f) : rgba(0, 255, 204, 1.f);
        for (size_t i = 1; i < pathL.size(); ++i) {
            int ax = (int)(cx + pathL[i - 1].x * SC), ay = (int)(CY - pathL[i - 1].y * SC);
            int bx = (int)(cx + pathL[i].x * SC), by = (int)(CY - pathL[i].y * SC);
            if (ax < x0 || ax > x1 || bx < x0 || bx > x1) continue;
            if (ay < TOP || ay > BOT || by < TOP || by > BOT) continue;
            olivec_line(oc, ax, ay, bx, by, rc);
        }
    }
    // ---- right half : flowing grid, constant c ----
    {
        int cx = HW + HW / 2, x0 = HW + 2, x1 = FW - 1;
        uint32_t g1 = rgba(0, 52, 62, 1.f);
        for (double g = -60; g <= 60; g += 1.0) {
            int a = (int)(cx + g * SC); if (a > x0 && a < x1) olivec_line(oc, a, TOP, a, BOT, g1);
            int b = (int)(CY - g * SC); if (b > TOP && b < BOT) olivec_line(oc, x0, b, x1, b, g1);
        }
        for (int k = 1; k <= 8; ++k) {                       // equal-flow rings
            double frac = 0.1 * k;                           // v = frac  ->  r = 2m/frac^2
            double r = 2.0 / (frac * frac);
            if (r < 60.0) ring(oc, cx, CY, (int)(r * SC), rgba(0, 56, 74, 1.f), x0, x1, TOP, BOT);
        }
        ring(oc, cx, CY, (int)(B_PG * SC), rgba(150, 80, 190, 1.f), x0, x1, TOP, BOT);
        ring(oc, cx, CY, (int)(R_PG * SC), rgba(230, 180, 55, 1.f), x0, x1, TOP, BOT);
        ring(oc, cx, CY, (int)(2.0 * SC), rgba(255, 80, 100, 1.f), x0, x1, TOP, BOT);
        olivec_circle(oc, cx, CY, 3, rgba(0, 0, 0, 1.f));
        uint32_t rc = capR ? rgba(255, 75, 95, 1.f) : rgba(0, 255, 204, 1.f);
        for (size_t i = 1; i < pathR.size(); ++i) {
            int ax = (int)(cx + pathR[i - 1].x * SC), ay = (int)(CY - pathR[i - 1].y * SC);
            int bx = (int)(cx + pathR[i].x * SC), by = (int)(CY - pathR[i].y * SC);
            if (ax < x0 || ax > x1 || bx < x0 || bx > x1) continue;
            if (ay < TOP || ay > BOT || by < TOP || by > BOT) continue;
            olivec_line(oc, ax, ay, bx, by, rc);
        }
    }

    // ---- HUD ----
    olivec_text(oc, "two grids, one shadow. coordinates differ, observables cannot.",
                14, 8, ft, 2, rgba(140, 230, 210, 1.f));
    olivec_text(oc, "left  static grid, light speed varies, c goes to 0 at rho m over 2",
                14, 34, ft, 2, rgba(110, 200, 190, 1.f));
    olivec_text(oc, "right flowing grid, light speed exactly c, river reaches c at r 2m",
                14, 56, ft, 2, rgba(110, 180, 210, 1.f));
    std::snprintf(buf, sizeof buf, "b crit   left %.12f     right %.12f", B_ISO, B_PG);
    olivec_text(oc, buf, 14, 80, ft, 2, rgba(255, 220, 120, 1.f));
    std::snprintf(buf, sizeof buf, "sqrt27 %.12f   difference %.2e   same to machine precision",
                  std::sqrt(27.0), std::fabs(B_ISO - B_PG));
    olivec_text(oc, buf, 14, FH - 78, ft, 2, rgba(180, 200, 140, 1.f));
    std::snprintf(buf, sizeof buf, "photon orbit areal r   left %.6f   right %.6f   both 3m",
                  areal(RHO_ISO), R_PG);
    olivec_text(oc, buf, 14, FH - 52, ft, 2, rgba(180, 200, 140, 1.f));
    std::snprintf(buf, sizeof buf, "b over b crit minus 1 1e%+.2f   deflection left %.5f right %.5f rad",
                  g_leps, std::fabs(thL) - M_PI, thR);
    olivec_text(oc, buf, 14, FH - 26, ft, 2, rgba(120, 210, 190, 1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
int main(int argc, char** argv) {
    sim_init(0, 0);
    printf("static grid  : b_crit = %.12f  at rho = %.9f  (areal %.9f)\n",
           B_ISO, RHO_ISO, areal(RHO_ISO));
    printf("flowing grid : b_crit = %.12f  at r   = %.9f\n", B_PG, R_PG);
    printf("sqrt(27)     =          %.12f\n", std::sqrt(27.0));
    printf("difference between the two grids = %.3e   ratio = %.12f\n",
           std::fabs(B_ISO - B_PG), B_ISO / B_PG);
    printf("areal radius of the photon orbit: left %.9f   right %.9f   diff %.2e\n",
           areal(RHO_ISO), R_PG, std::fabs(areal(RHO_ISO) - R_PG));
    int steps = argc > 1 ? atoi(argv[1]) : 1;
    g_leps = std::log10(3.094717e-2); dirty = true;
    sim_step(steps);
    printf("at ring 1: deflection left %.6f rad, right %.6f rad, target pi = %.6f\n",
           std::fabs(thL) - M_PI, thR, M_PI);
    uint8_t* p = sim_render();
    stbi_write_png("twogrid_os_preview.png", FW, FH, 4, p, FW * 4);
    printf("wrote twogrid_os_preview.png\n");
    return 0;
}
#endif
