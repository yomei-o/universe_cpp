// Universe OS — Photon Sphere  (C++/WASM)
//
// At areal radius r = 3m (isotropic rho = 1 + sqrt3/2 = 1.866025 m) light can
// orbit in a circle. On the flat grid that circle has a one-line reason:
//
//     turning point of a ray  <=>  rho / c(rho) = b        with  c = (1-u)/(1+u)^3
//
// V(rho) = rho/c(rho) goes to infinity at BOTH ends -- at large rho because rho
// grows, and at rho -> m/2 because c -> 0. So it has a MINIMUM, and the minimum
// is the circular photon orbit:
//
//     rho_ph = 1.866025404 m   (areal 3m),   b_crit = V_min = sqrt(27) m
//
// But it is a MAXIMUM of the effective potential for the radial motion, so the
// orbit is UNSTABLE. This demo makes that concrete: launch a photon exactly
// tangentially at radius rho and watch it
//   * spiral in      if rho < rho_ph
//   * orbit and peel if rho = rho_ph (numerically it always peels: the error
//     grows by e^(2pi) = 535.49 per turn, which is the Lyapunov factor)
//   * spiral out     if rho > rho_ph
//
// The instability is not a numerical artefact -- it is the SAME e^(2pi) that
// sets the spacing of the photon rings in flatgrid_os. One number does both:
// how fast a photon leaves the sphere, and how tightly the rings crowd the
// shadow edge.
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

static const int FW = 900, FH = 660;
static std::vector<uint32_t> px;

static inline uint32_t rgba(int r, int g, int b, float a) {
    int A = (int)(a * 255.0f); if (A < 0) A = 0; if (A > 255) A = 255;
    return ((uint32_t)A << 24) | ((uint32_t)b << 16) | ((uint32_t)g << 8) | (uint32_t)r;
}
static void ring(Olivec_Canvas oc, int cx, int cy, int r, uint32_t col) {
    if (r < 1) return;
    int x = r, y = 0, err = 1 - r;
    auto put = [&](int a, int b) {
        if (a >= 0 && a < (int)oc.width && b >= 0 && b < (int)oc.height)
            OLIVEC_PIXEL(oc, a, b) = col;
    };
    while (x >= y) {
        put(cx + x, cy + y); put(cx + y, cy + x); put(cx - y, cy + x); put(cx - x, cy + y);
        put(cx - x, cy - y); put(cx - y, cy - x); put(cx + y, cy - x); put(cx + x, cy - y);
        ++y;
        if (err < 0) err += 2 * y + 1; else { --x; err += 2 * (y - x) + 1; }
    }
}

// ---------------------------------------------------------------- the grid
static inline void FG(double rho, double& F, double& G, double& dF, double& dG) {
    double u = 0.5 / rho, du = -u / rho;
    double om = 1.0 - u, op = 1.0 + u;
    F = (op * op) / (om * om);
    G = 1.0 / (op * op * op * op);
    dF = (4.0 * op / (om * om * om)) * du;
    dG = (-4.0 / (op * op * op * op * op)) * du;
}
static inline double c_local(double rho) {
    double F, G, dF, dG; FG(rho, F, G, dF, dG); return std::sqrt(G / F);
}
static inline double Vpot(double rho) { return rho / c_local(rho); }

static double RHO_PH = 0.0, B_CRIT = 0.0;
static void find_ph() {
    double lo = 0.5000001, hi = 60.0;
    for (int i = 0; i < 300; ++i) {
        double a = lo + (hi - lo) / 3.0, b = hi - (hi - lo) / 3.0;
        if (Vpot(a) < Vpot(b)) hi = b; else lo = a;
    }
    RHO_PH = 0.5 * (lo + hi); B_CRIT = Vpot(RHO_PH);
}

// ---------------------------------------------------------------- one photon
struct Ph { double x, y, px_, py_; bool alive; double th; int fate; };  // fate 0 live 1 in 2 out
static std::vector<Ph> phs;
static std::vector<std::vector<float>> trails;      // x,y pairs

static void rhs(const double s[4], double d[4]) {
    double rho = std::sqrt(s[0] * s[0] + s[1] * s[1]);
    double F, G, dF, dG; FG(rho, F, G, dF, dG);
    double p2 = s[2] * s[2] + s[3] * s[3];
    double k = -0.5 * (-dF + p2 * dG) / rho;
    d[0] = s[2] * G; d[1] = s[3] * G; d[2] = k * s[0]; d[3] = k * s[1];
}

// launch tangentially at radius r0, angle a0
static void launch(double r0, double a0) {
    double F, G, dF, dG; FG(r0, F, G, dF, dG);
    double pm = std::sqrt(F / G);                    // |p| for a null ray, E = 1
    double cx = std::cos(a0), sy = std::sin(a0);
    Ph p;
    p.x = r0 * cx; p.y = r0 * sy;
    p.px_ = -pm * sy; p.py_ = pm * cx;               // tangential
    p.alive = true; p.th = 0.0; p.fate = 0;
    phs.push_back(p);
    trails.push_back({ (float)p.x, (float)p.y });
}

static double g_r0 = 1.866025404;
static double g_view = 7.0;
static double g_dt = 1.0;               // step multiplier

extern "C" {

KEEP int sim_w() { return FW; }
KEEP int sim_h() { return FH; }

KEEP void sim_reset() {
    phs.clear(); trails.clear();
    g_r0 = RHO_PH; g_view = 7.0; g_dt = 1.0;
    launch(g_r0, 0.0);
}
KEEP int sim_init(int, int) { px.assign((size_t)FW * FH, 0); find_ph(); sim_reset(); return 1; }

KEEP void sim_set(int id, double v) {
    if (id == 0) { g_r0 = v; }
    else if (id == 1) { g_view = v; }
    else if (id == 2) { g_dt = v; }
}
KEEP void sim_action(int a) {
    if (a == 0) { phs.clear(); trails.clear(); launch(g_r0, 0.0); }
    else if (a == 1) { phs.clear(); trails.clear(); launch(RHO_PH, 0.0); g_r0 = RHO_PH; }
    else if (a == 2) { phs.clear(); trails.clear(); launch(RHO_PH * 0.97, 0.0); g_r0 = RHO_PH * 0.97; }
    else if (a == 3) { phs.clear(); trails.clear(); launch(RHO_PH * 1.03, 0.0); g_r0 = RHO_PH * 1.03; }
    else if (a == 4) {                       // a whole family at once
        phs.clear(); trails.clear();
        for (int k = -3; k <= 3; ++k) launch(RHO_PH * (1.0 + 0.02 * k), 0.0);
    }
}
KEEP void sim_click(double nx, double ny) {
    double x = (nx - 0.5) * 2.0 * g_view;
    double y = (0.5 - ny) * 2.0 * g_view * (double)FH / (double)FW;
    double r = std::sqrt(x * x + y * y);
    if (r > 0.55) { g_r0 = r; phs.clear(); trails.clear(); launch(r, std::atan2(y, x)); }
}

KEEP void sim_step(int steps) {
    for (int s = 0; s < (steps > 0 ? steps : 1); ++s) {
        for (size_t i = 0; i < phs.size(); ++i) {
            Ph& p = phs[i];
            if (!p.alive) continue;
            for (int sub = 0; sub < 4; ++sub) {
                double st[4] = { p.x, p.y, p.px_, p.py_ };
                double rho = std::sqrt(st[0] * st[0] + st[1] * st[1]);
                double d0 = rho - 0.5;
                double h = 0.0025 * g_dt * (d0 < 1.0 ? (d0 > 0.02 ? d0 : 0.02) : rho);
                double a[4], b[4], c[4], d[4], t[4];
                rhs(st, a);
                for (int k = 0; k < 4; ++k) t[k] = st[k] + 0.5 * h * a[k]; rhs(t, b);
                for (int k = 0; k < 4; ++k) t[k] = st[k] + 0.5 * h * b[k]; rhs(t, c);
                for (int k = 0; k < 4; ++k) t[k] = st[k] + h * c[k]; rhs(t, d);
                double before = std::atan2(st[1], st[0]);
                for (int k = 0; k < 4; ++k) st[k] += h / 6.0 * (a[k] + 2 * b[k] + 2 * c[k] + d[k]);
                if (!(st[0] == st[0]) || !(st[1] == st[1])) { p.alive = false; p.fate = 1; break; }
                double after = std::atan2(st[1], st[0]);
                double dth = after - before;
                while (dth > M_PI) dth -= 2 * M_PI;
                while (dth < -M_PI) dth += 2 * M_PI;
                p.th += dth;
                p.x = st[0]; p.y = st[1]; p.px_ = st[2]; p.py_ = st[3];
                double rr = std::sqrt(p.x * p.x + p.y * p.y);
                if (rr < 0.53) { p.alive = false; p.fate = 1; break; }
                if (rr > 3.0 * g_view + 8.0) { p.alive = false; p.fate = 2; break; }
            }
            if (trails[i].size() < 200000) {
                trails[i].push_back((float)p.x); trails[i].push_back((float)p.y);
            }
        }
    }
}

KEEP uint8_t* sim_render() {
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(3, 5, 11, 1.f));
    Olivec_Font ft = olivec_default_font;
    char buf[190];
    const int PY0 = FH - 168, PY1 = FH - 48, PL = 74, PR = FW - 40;
    const int TOP0 = 96, TOP1 = PY0 - 30;
    const int CX = FW / 2, CY = (TOP0 + TOP1) / 2;
    double SC = (FW * 0.5) / g_view;
    auto SX = [&](double x) { return (int)(CX + x * SC); };
    auto SY = [&](double y) { return (int)(CY - y * SC); };

    uint32_t gcol = rgba(0, 54, 58, 1.f);
    for (double g = -40; g <= 40; g += 1.0) {
        int a = SX(g); if (a >= 0 && a < FW) olivec_line(oc, a, TOP0, a, TOP1, gcol);
        int b = SY(g); if (b >= TOP0 && b <= TOP1) olivec_line(oc, 0, b, FW - 1, b, gcol);
    }
    ring(oc, CX, CY, (int)(B_CRIT * SC), rgba(140, 75, 180, 1.f));
    ring(oc, CX, CY, (int)(RHO_PH * SC), rgba(230, 180, 55, 1.f));
    ring(oc, CX, CY, (int)(RHO_PH * SC) + 1, rgba(120, 92, 25, 1.f));
    olivec_circle(oc, CX, CY, (int)(0.5 * SC) + 2, rgba(0, 0, 0, 1.f));
    ring(oc, CX, CY, (int)(0.5 * SC) + 2, rgba(110, 55, 140, 1.f));

    // trails
    for (size_t i = 0; i < trails.size(); ++i) {
        const auto& tr = trails[i];
        uint32_t col = (phs[i].fate == 1) ? rgba(255, 80, 100, 1.f)
                     : (phs[i].fate == 2) ? rgba(120, 200, 255, 1.f)
                                          : rgba(0, 255, 204, 1.f);
        for (size_t k = 2; k + 1 < tr.size(); k += 2) {
            int ax = SX(tr[k - 2]), ay = SY(tr[k - 1]), bx = SX(tr[k]), by = SY(tr[k + 1]);
            if (ay < TOP0 || ay > TOP1 || by < TOP0 || by > TOP1) continue;
            olivec_line(oc, ax, ay, bx, by, col);
        }
        if (phs[i].alive) {
            int sy = SY(phs[i].y);
            if (sy >= TOP0 && sy <= TOP1)
                olivec_circle(oc, SX(phs[i].x), sy, 3, rgba(255, 255, 210, 1.f));
        }
    }

    // ---- V(rho) = rho/c(rho) profile: the minimum IS the photon sphere ----
    olivec_rect(oc, 0, PY0 - 26, FW, FH - (PY0 - 26), rgba(2, 4, 9, 1.f));
    auto QX = [&](double r) { return (int)(PL + (r / 8.0) * (PR - PL)); };
    auto QY = [&](double v) { return (int)(PY1 - ((v - 4.0) / 6.0) * (PY1 - PY0)); };
    olivec_line(oc, PL, PY1, PR, PY1, rgba(60, 90, 100, 1.f));
    olivec_line(oc, PL, PY0, PL, PY1, rgba(60, 90, 100, 1.f));
    olivec_line(oc, PL, QY(B_CRIT), PR, QY(B_CRIT), rgba(110, 60, 140, 1.f));
    int ox = 0, oy = 0; bool first = true;
    for (int i = 1; i <= 700; ++i) {
        double r = 0.52 + (8.0 - 0.52) * i / 700.0, v = Vpot(r);
        if (v > 10.0 || v < 4.0) { first = true; continue; }
        int cx = QX(r), cy = QY(v);
        if (!first) olivec_line(oc, ox, oy, cx, cy, rgba(0, 210, 250, 1.f));
        ox = cx; oy = cy; first = false;
    }
    olivec_circle(oc, QX(RHO_PH), QY(B_CRIT), 4, rgba(230, 180, 55, 1.f));
    olivec_line(oc, QX(g_r0), PY0, QX(g_r0), PY1, rgba(255, 220, 120, 1.f));

    // ---- HUD (lowercase only) ----
    olivec_text(oc, "photon sphere - light can circle, but the orbit is unstable", 14, 8, ft, 2,
                rgba(140, 230, 210, 1.f));
    std::snprintf(buf, sizeof buf, "rho ph %.9f m, areal %.6f m.  b crit %.9f m",
                  RHO_PH, RHO_PH * std::pow(1.0 + 0.5 / RHO_PH, 2), B_CRIT);
    olivec_text(oc, buf, 14, 36, ft, 2, rgba(230, 180, 55, 1.f));
    double turns = phs.empty() ? 0.0 : std::fabs(phs[0].th) / (2.0 * M_PI);
    const char* fate = phs.empty() ? "" : (phs[0].fate == 1 ? "fell in" :
                       phs[0].fate == 2 ? "escaped" : "still circling");
    std::snprintf(buf, sizeof buf, "launch radius %.9f m,  offset %+.3e   %s",
                  g_r0, g_r0 / RHO_PH - 1.0, fate);
    olivec_text(oc, buf, 14, 64, ft, 2, rgba(255, 220, 120, 1.f));
    std::snprintf(buf, sizeof buf, "turns so far %.3f     one turn multiplies any error by 535.49",
                  turns);
    olivec_text(oc, buf, 14, 92, ft, 2, rgba(110, 190, 175, 1.f));
    olivec_text(oc, "v of rho  rho over c of rho . the minimum is the photon sphere",
                PL + 4, PY0 - 18, ft, 2, rgba(0, 190, 225, 1.f));
    olivec_text(oc, "cyan still circling, red fell in, blue escaped", 14, FH - 24, ft, 2,
                rgba(120, 170, 190, 1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
int main(int argc, char** argv) {
    sim_init(0, 0);
    printf("rho_ph  = %.9f   1+sqrt3/2 = %.9f\n", RHO_PH, 1.0 + std::sqrt(3.0) / 2.0);
    printf("areal r = %.9f   (should be 3)\n", RHO_PH * std::pow(1.0 + 0.5 / RHO_PH, 2));
    printf("b_crit  = %.12f   sqrt27 = %.12f\n", B_CRIT, std::sqrt(27.0));
    printf("V is a MINIMUM here: V(0.98 rho) %.6f  V(rho) %.6f  V(1.02 rho) %.6f\n",
           Vpot(0.98 * RHO_PH), Vpot(RHO_PH), Vpot(1.02 * RHO_PH));
    // The instability rate IS e^(2pi). Cutting the launch offset by 100 should
    // buy exactly ln(100)/(2pi) = 0.73299 extra turns, because that is how long
    // the exponential takes to eat two more decades. It is the same number that
    // sets the photon-ring spacing in flatgrid_os: one constant does both.
    printf("\ninstability: turns survived vs launch offset\n");
    double prev = -1.0;
    for (double off : {1e-2, 1e-4, 1e-6, 1e-8}) {
        phs.clear(); trails.clear(); launch(RHO_PH * (1.0 + off), 0.0);
        for (int i = 0; i < 60000 && phs[0].alive; ++i) sim_step(1);
        double t = std::fabs(phs[0].th) / (2.0 * M_PI);
        printf("  offset %+8.0e -> %-8s after %.3f turns", off,
               phs[0].fate == 1 ? "fell in" : (phs[0].fate == 2 ? "escaped" : "going"), t);
        if (prev >= 0.0)
            printf("   gain %.3f   ln100/2pi = %.5f", t - prev, std::log(100.0) / (2.0 * M_PI));
        printf("\n");
        prev = t;
    }
    printf("  last row is error-dominated: the intended offset is below the\n");
    printf("  integrator noise, so it survives longer than the law predicts.\n");
    int steps = argc > 1 ? atoi(argv[1]) : 4000;
    phs.clear(); trails.clear(); launch(RHO_PH, 0.0); g_r0 = RHO_PH;
    sim_step(steps);
    uint8_t* p = sim_render();
    stbi_write_png("photonsphere_os_preview.png", FW, FH, 4, p, FW * 4);
    printf("wrote photonsphere_os_preview.png (%.3f turns)\n",
           std::fabs(phs[0].th) / (2.0 * M_PI));
    return 0;
}
#endif
