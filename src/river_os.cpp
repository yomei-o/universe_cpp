// Universe OS — River (Painleve-Gullstrand)  (C++/WASM)
//
// The other way to write the same black hole. Here space is EXACTLY Euclidean --
// not conformally stretched, not curved, no factor in front of dx^2 at all --
// and the lapse is EXACTLY 1, so the speed of light with respect to the water is
// EXACTLY c. What moves is the grid itself:
//
//     ds^2 = -dt^2 + ( dr + v(r) dt )^2 + r^2 dOmega^2 ,   v(r) = sqrt(2m/r)
//
// v(r) is Newton's escape velocity. The horizon r = 2m is not "where light
// stops"; it is where THE RIVER REACHES THE SPEED OF LIGHT. Inside, the river is
// supersonic: an outgoing photon still moves at c through the water and is still
// carried downstream.
//
//     radial light:   dr/dt = +-1 - v(r)
//     free fall from rest at infinity:  dr/dt = -v(r)      <- Newton's own law
//
// Nothing negative ever appears. Compare flatgrid_os, where the same spacetime is
// written with a static grid and a variable light speed that hits zero at the
// throat: two flat grids, same physics, and they agree on the shadow radius
// sqrt(27) m to 10 digits.
//
// Proper time from horizon to singularity (verified against miharashi/flatgrid5.py):
//     from rest at infinity   tau = 4m/3   = 1.333333332 m
//     longest possible        tau = pi m   = 3.141592654 m
// For a solar mass that is 6.567 us; for M87* it is 711.5 minutes.
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
// olivec_circle fills; we need outlines
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

// m = 1 throughout
static inline double vflow(double r) { return (r > 1e-6) ? std::sqrt(2.0 / r) : 1e3; }

// ------------------------------------------------------------ swimmers
// Each tracer is a radial photon or a faller. Advanced in coordinate time.
struct Swim { double r, th; int kind; bool alive; };   // kind 0 out, 1 in, 2 faller
static std::vector<Swim> sw;
static double g_view = 8.0;         // view half-width in m
static double g_dt = 0.012;
static int    g_show = 1;           // 1 = show the flow arrows
static long   g_tick = 0;

static void seed() {
    sw.clear();
    // rings of outgoing photons at several radii, plus infallers
    const double rs[] = { 1.2, 1.6, 2.0, 2.6, 3.4, 4.4, 5.6 };
    for (double r : rs)
        for (int k = 0; k < 16; ++k)
            sw.push_back({ r, 2.0 * M_PI * k / 16.0, 0, true });
    for (int k = 0; k < 10; ++k) sw.push_back({ 7.4, 2.0 * M_PI * k / 10.0 + 0.31, 2, true });
    for (int k = 0; k < 10; ++k) sw.push_back({ 7.0, 2.0 * M_PI * k / 10.0 + 0.10, 1, true });
}

extern "C" {

KEEP int sim_w() { return FW; }
KEEP int sim_h() { return FH; }
KEEP void sim_reset() { g_view = 8.0; g_dt = 0.012; g_show = 1; g_tick = 0; seed(); }
KEEP int sim_init(int, int) { px.assign((size_t)FW * FH, 0); sim_reset(); return 1; }

KEEP void sim_set(int id, double v) {
    if (id == 0) g_dt = v; else if (id == 1) g_view = v; else if (id == 2) g_show = (int)(v + 0.5);
}
KEEP void sim_action(int a) {
    if (a == 0) seed();
    else if (a == 1) {                       // one outgoing photon exactly at the horizon
        sw.clear();
        for (int k = 0; k < 24; ++k) sw.push_back({ 2.0, 2.0 * M_PI * k / 24.0, 0, true });
    } else if (a == 2) {                     // outgoing photons started INSIDE
        sw.clear();
        for (int k = 0; k < 24; ++k) sw.push_back({ 1.4, 2.0 * M_PI * k / 24.0, 0, true });
    } else if (a == 3) { g_show = !g_show; }
}
// click drops an outgoing photon there
KEEP void sim_click(double nx, double ny) {
    double x = (nx - 0.5) * 2.0 * g_view;
    double y = (0.5 - ny) * 2.0 * g_view * (double)FH / (double)FW;
    double r = std::sqrt(x * x + y * y);
    if (r > 0.05) sw.push_back({ r, std::atan2(y, x), 0, true });
}

KEEP void sim_step(int steps) {
    for (int s = 0; s < (steps > 0 ? steps : 1); ++s) {
        for (auto& p : sw) {
            if (!p.alive) continue;
            double v = vflow(p.r);
            double dr = (p.kind == 0 ? (1.0 - v) : (p.kind == 1 ? (-1.0 - v) : -v));
            p.r += dr * g_dt;
            if (p.r <= 0.02) { p.alive = false; }
            if (p.r > 1.4 * g_view + 6.0) p.alive = false;
        }
        ++g_tick;
        if ((g_tick % 900) == 0) seed();
    }
}

KEEP uint8_t* sim_render() {
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(3, 6, 12, 1.f));
    Olivec_Font ft = olivec_default_font;
    char buf[190];
    const int PY0 = FH - 158, PY1 = FH - 46, PL = 70, PR = FW - 40;
    const int TOP0 = 92, TOP1 = PY0 - 30;          // upper panel bounds
    const int CX = FW / 2, CY = (TOP0 + TOP1) / 2;
    double SC = (FW * 0.5) / g_view;
    auto SX = [&](double x) { return (int)(CX + x * SC); };
    auto SY = [&](double y) { return (int)(CY - y * SC); };

    // ---- the grid: EXACTLY Euclidean. No conformal factor anywhere. ----
    uint32_t gcol = rgba(0, 58, 62, 1.f);
    for (double g = -40; g <= 40; g += 1.0) {
        int a = SX(g); if (a >= 0 && a < FW) olivec_line(oc, a, TOP0, a, TOP1, gcol);
        int b = SY(g); if (b >= TOP0 && b <= TOP1) olivec_line(oc, 0, b, FW - 1, b, gcol);
    }

    // ---- the flow: inward arrows whose length is v(r) ----
    if (g_show) {
        for (int i = 0; i < 26; ++i) {
            double th = 2.0 * M_PI * i / 26.0;
            for (double r = 0.55; r < g_view * 1.35; r += 0.62) {
                double v = vflow(r); if (v > 3.2) continue;
                double x = r * std::cos(th), y = r * std::sin(th);
                double L = 0.30 * v;                       // arrow length tracks the flow
                int x0 = SX(x), y0 = SY(y);
                int x1 = SX(x - L * std::cos(th)), y1 = SY(y - L * std::sin(th));
                if (y0 < TOP0 || y0 > TOP1 || y1 < TOP0 || y1 > TOP1) continue;
                uint32_t col = (v >= 1.0) ? rgba(200, 60, 90, 1.f) : rgba(40, 120, 170, 1.f);
                olivec_line(oc, x0, y0, x1, y1, col);
            }
        }
    }

    // ---- horizon (river = c) and the photon sphere ----
    ring(oc, CX, CY, (int)(2.0 * SC), rgba(255, 80, 100, 1.f));
    ring(oc, CX, CY, (int)(2.0 * SC) - 1, rgba(150, 40, 60, 1.f));
    ring(oc, CX, CY, (int)(3.0 * SC), rgba(210, 165, 50, 1.f));
    olivec_circle(oc, CX, CY, (int)(0.06 * SC) + 2, rgba(0, 0, 0, 1.f));

    // ---- swimmers ----
    for (auto& p : sw) {
        if (!p.alive) continue;
        double x = p.r * std::cos(p.th), y = p.r * std::sin(p.th);
        uint32_t col = (p.kind == 0) ? rgba(0, 255, 204, 1.f)
                     : (p.kind == 1) ? rgba(120, 150, 255, 1.f)
                                     : rgba(255, 220, 120, 1.f);
        int sy = SY(y); if (sy < TOP0 + 3 || sy > TOP1 - 3) continue;
        olivec_circle(oc, SX(x), sy, (p.kind == 0 ? 3 : 2), col);
    }

    // ---- v(r) profile strip along the bottom (own opaque panel) ----
    olivec_rect(oc, 0, PY0 - 26, FW, FH - (PY0 - 26), rgba(2, 4, 9, 1.f));
    for (int k = 0; k <= 10; ++k) {
        int gx = PL + k * (PR - PL) / 10;
        olivec_line(oc, gx, PY0, gx, PY1, rgba(0, 44, 48, 1.f));
    }
    olivec_line(oc, PL, PY1, PR, PY1, rgba(60, 90, 100, 1.f));
    olivec_line(oc, PL, PY0, PL, PY1, rgba(60, 90, 100, 1.f));
    auto QX = [&](double r) { return (int)(PL + (r / 10.0) * (PR - PL)); };
    auto QY = [&](double v) { return (int)(PY1 - (v / 3.0) * (PY1 - PY0)); };
    olivec_line(oc, PL, QY(1.0), PR, QY(1.0), rgba(90, 70, 40, 1.f));   // v = c
    int ox = 0, oy = 0; bool first = true;
    for (int i = 1; i <= 400; ++i) {
        double r = 10.0 * i / 400.0, v = vflow(r); if (v > 3.0) { first = true; continue; }
        int cx = QX(r), cy = QY(v);
        if (!first) olivec_line(oc, ox, oy, cx, cy, rgba(40, 140, 200, 1.f));
        ox = cx; oy = cy; first = false;
    }
    // outgoing light speed 1 - v : crosses zero exactly at r = 2m
    first = true;
    for (int i = 1; i <= 400; ++i) {
        double r = 10.0 * i / 400.0, o = 1.0 - vflow(r);
        if (o < -2.0) { first = true; continue; }
        int cx = QX(r), cy = QY(o);
        if (!first) olivec_line(oc, ox, oy, cx, cy, rgba(0, 255, 204, 1.f));
        ox = cx; oy = cy; first = false;
    }
    olivec_line(oc, QX(2.0), PY0, QX(2.0), PY1, rgba(255, 80, 100, 1.f));

    // ---- HUD (lowercase only: olive's font has no capitals) ----
    olivec_text(oc, "river model - space is exactly euclidean. the grid flows.", 14, 8, ft, 2,
                rgba(140, 230, 210, 1.f));
    olivec_text(oc, "lapse is 1, so light does exactly c through the water", 14, 34, ft, 2,
                rgba(110, 190, 175, 1.f));
    std::snprintf(buf, sizeof buf, "v of r  sqrt of 2m over r , newtons escape speed. v %.3f at r 2m",
                  vflow(2.0));
    olivec_text(oc, buf, 14, 60, ft, 2, rgba(255, 200, 90, 1.f));
    olivec_text(oc, "red arrows mean the river is faster than light", 14, FH - 24, ft, 2,
                rgba(255, 110, 130, 1.f));
    olivec_text(oc, "cyan  1 minus v , outgoing light. zero at r 2m", PL + 250, PY0 - 18, ft, 2,
                rgba(0, 220, 180, 1.f));
    olivec_text(oc, "blue  v of r", PL + 6, PY0 - 18, ft, 2, rgba(60, 160, 220, 1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
static double tau_from_infinity() {           // int_0^{2m} sqrt(r/2m) dr = 4m/3
    const int N = 200000; double h = 2.0 / N, s = 0.0;
    for (int i = 0; i <= N; ++i) {
        double r = i * h, f = std::sqrt(r / 2.0);
        s += ((i == 0 || i == N) ? 1.0 : ((i & 1) ? 4.0 : 2.0)) * f;
    }
    return s * h / 3.0;
}
static double tau_from_horizon() {            // r = 2m sin^2(eta) -> int 4m sin^2 = pi m
    const int N = 400000; double h = (M_PI / 2) / N, s = 0.0;
    for (int i = 0; i <= N; ++i) {
        double e = i * h, f = 4.0 * std::sin(e) * std::sin(e);
        s += ((i == 0 || i == N) ? 1.0 : ((i & 1) ? 4.0 : 2.0)) * f;
    }
    return s * h / 3.0;
}
int main(int argc, char** argv) {
    sim_init(0, 0);
    printf("v(2m) = %.12f   (should be exactly 1)\n", vflow(2.0));
    printf("outgoing dr/dt at r=2m  = %.3e  (should be 0)\n", 1.0 - vflow(2.0));
    for (double r : {10.0, 4.0, 2.0, 1.0, 0.25}) {
        printf("r=%5.2f m: v=%.6f  out=%+.6f  in=%+.6f\n",
               r, vflow(r), 1.0 - vflow(r), -1.0 - vflow(r));
    }
    printf("\nproper time horizon to singularity:\n");
    printf("  from rest at infinity  %.9f m   (4/3 = %.9f)\n", tau_from_infinity(), 4.0 / 3.0);
    printf("  longest possible       %.9f m   (pi  = %.9f)\n", tau_from_horizon(), M_PI);
    const double MSUN = 1476.6250, C = 299792458.0;
    printf("  solar mass: %.3f us      M87*: %.1f min\n",
           4.0 / 3.0 * MSUN / C * 1e6, 4.0 / 3.0 * 6.5e9 * MSUN / C / 60.0);
    int steps = argc > 1 ? atoi(argv[1]) : 40;
    sim_step(steps);
    uint8_t* p = sim_render();
    stbi_write_png("river_os_preview.png", FW, FH, 4, p, FW * 4);
    printf("wrote river_os_preview.png\n");
    return 0;
}
#endif
