// Universe OS — Gauge Shock  (C++/WASM)
//
// A bug inside perfect physics.
//
// Numerical relativity does not just evolve the geometry; it also has to CHOOSE
// the coordinates as it goes. The lapse alpha (how much proper time passes per
// step of the time coordinate) is evolved by its own equation. The standard
// family is Bona-Masso:
//
//     d(alpha)/dt = -alpha^2 f(alpha) K
//
// In 1+1 relativity the gauge sector is hyperbolic and its characteristic speed
// is
//     lambda(alpha) = alpha * sqrt( f(alpha) )
//
// The speed depends on the very thing that is propagating. So fast parts of the
// profile overtake slow parts, the profile steepens, and THE CHARACTERISTICS
// CROSS. That is a genuine shock -- but a shock in the COORDINATES, not in the
// spacetime. The geometry is untouched; the chart tears. Codes crash, waveforms
// go bad, and no physics is wrong anywhere.
//
// This demo solves the simple wave EXACTLY by characteristics, so nothing here
// can be blamed on a finite-difference scheme:
//
//     alpha is constant along  x(xi,t) = xi + lambda(alpha_0(xi)) t
//     dx/dxi = 1 + t lambda'(alpha_0) alpha_0'(xi)
//     shock at  t_s = min over xi of  -1 / ( lambda'(alpha_0) alpha_0'(xi) )
//
// With alpha_0 = 1 + A sin(2 pi x) and A = 0.2 the exact shock times are
//     harmonic      f = 1            lambda = alpha            t_s = 0.7958
//     1+log         f = 2/alpha      lambda = sqrt(2 alpha)    t_s = 1.1254
//     shock-avoid   f = 1 + k/a^2    lambda = sqrt(a^2+k)      t_s = 2.6393 (k=10)
//                                                              t_s = 7.9974 (k=100)
//
// HONEST LIMIT OF THIS MODEL: here lambda' = 0 needs alpha^2 f = const, i.e.
// f = k/alpha^2. Alcubierre's actual shock-avoiding condition for the full 1+1
// system is  alpha f'(alpha) + 2(f-1) = 0,  whose solution is f = 1 + k/alpha^2
// -- the "+1" comes from the two-variable system this reduction does not carry.
// The mechanism is the same; the constant term is not reproduced. And whether
// these shocks actually bite in 3D runs is contested (Reimann et al., "Are gauge
// shocks really shocks?", CQG 22 (2005) 4215).
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

static const int FW = 900, FH = 640;
static std::vector<uint32_t> px;

static inline uint32_t rgba(int r, int g, int b, float a) {
    int A = (int)(a * 255.0f); if (A < 0) A = 0; if (A > 255) A = 255;
    return ((uint32_t)A << 24) | ((uint32_t)b << 16) | ((uint32_t)g << 8) | (uint32_t)r;
}

// ------------------------------------------------------------ the slicing choice
// 0 = harmonic f=1 ; 1 = 1+log f=2/alpha ; 2 = shock-avoiding f=1+kappa/alpha^2
static int   g_slice = 1;
static double g_kappa = 10.0;
static double g_amp   = 0.20;
static double g_t     = 0.0;
static double g_dt    = 0.004;

static inline double f_of(double a) {
    if (g_slice == 0) return 1.0;
    if (g_slice == 1) return 2.0 / a;
    return 1.0 + g_kappa / (a * a);
}
static inline double lam(double a) { return a * std::sqrt(f_of(a)); }
static inline double dlam(double a) {                     // d lambda / d alpha
    const double h = 1e-7;
    return (lam(a + h) - lam(a - h)) / (2.0 * h);
}
static inline double a0(double xi)  { return 1.0 + g_amp * std::sin(2.0 * M_PI * xi); }
static inline double da0(double xi) { return 2.0 * M_PI * g_amp * std::cos(2.0 * M_PI * xi); }

// dx/dxi along the characteristic map. Shock when this first reaches zero.
static inline double stretch(double xi, double t) { return 1.0 + t * dlam(a0(xi)) * da0(xi); }

static double shock_time() {
    double best = 1e30;
    for (int i = 0; i < 4000; ++i) {
        double xi = (double)i / 4000.0;
        double d = dlam(a0(xi)) * da0(xi);
        if (d < -1e-12) { double t = -1.0 / d; if (t < best) best = t; }
    }
    return best;
}
static double T_SHOCK = 1.0;

extern "C" {

KEEP int sim_w() { return FW; }
KEEP int sim_h() { return FH; }

KEEP void sim_reset() {
    g_t = 0.0; g_slice = 1; g_kappa = 10.0; g_amp = 0.20; g_dt = 0.004;
    T_SHOCK = shock_time();
}
KEEP int sim_init(int, int) { px.assign((size_t)FW * FH, 0); sim_reset(); return 1; }

KEEP void sim_set(int id, double v) {
    if (id == 0) g_dt = v;
    else if (id == 1) { g_amp = v; T_SHOCK = shock_time(); }
    else if (id == 2) { g_kappa = v; T_SHOCK = shock_time(); }
}
KEEP void sim_action(int a) {
    if (a >= 0 && a <= 2) { g_slice = a; g_t = 0.0; T_SHOCK = shock_time(); }
    else if (a == 3) g_t = 0.0;
    else if (a == 4) g_t = T_SHOCK;             // jump exactly to the shock
}
KEEP void sim_click(double nx, double ny) { (void)ny; g_t = nx * 2.6 * T_SHOCK; }

KEEP void sim_step(int steps) {
    for (int s = 0; s < (steps > 0 ? steps : 1); ++s) {
        g_t += g_dt;
        if (g_t > 2.6 * T_SHOCK) g_t = 0.0;
    }
}

// ------------------------------------------------------------ render
KEEP uint8_t* sim_render() {
    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(4, 5, 10, 1.f));
    Olivec_Font ft = olivec_default_font;
    char buf[200];

    const int L = 70, R = FW - 30;                 // x range in pixels
    const int TOP_Y0 = 150, TOP_Y1 = 340;           // alpha panel
    const int BOT_Y0 = 392, BOT_Y1 = FH - 46;      // x-t panel
    const double XW = 1.0;                         // x in [0,1)
    const double TVIEW = 2.6 * T_SHOCK;
    const double LAM_BAR = lam(1.0);   // co-moving frame for the x-t panel

    auto PX = [&](double x) { return (int)(L + (x / XW) * (R - L)); };
    auto PY_A = [&](double a) {                    // alpha in [1-1.6A, 1+1.6A]
        double lo = 1.0 - 1.7 * g_amp, hi = 1.0 + 1.7 * g_amp;
        return (int)(TOP_Y1 - (a - lo) / (hi - lo) * (TOP_Y1 - TOP_Y0));
    };
    auto PY_T = [&](double t) { return (int)(BOT_Y1 - (t / TVIEW) * (BOT_Y1 - BOT_Y0)); };

    uint32_t axis = rgba(50, 80, 90, 1.f), grid = rgba(0, 45, 45, 1.f);
    uint32_t cyan = rgba(0, 255, 204, 1.f), amber = rgba(255, 200, 90, 1.f);
    uint32_t red = rgba(255, 70, 90, 1.f), dim = rgba(0, 110, 100, 1.f);

    // ---------------- panel frames + gridlines ----------------
    for (int k = 0; k <= 10; ++k) {
        int x = PX(k / 10.0);
        olivec_line(oc, x, TOP_Y0, x, TOP_Y1, grid);
        olivec_line(oc, x, BOT_Y0, x, BOT_Y1, grid);
    }
    olivec_line(oc, L, TOP_Y1, R, TOP_Y1, axis);
    olivec_line(oc, L, TOP_Y0, L, TOP_Y1, axis);
    olivec_line(oc, L, BOT_Y1, R, BOT_Y1, axis);
    olivec_line(oc, L, BOT_Y0, L, BOT_Y1, axis);
    // alpha = 1 line
    olivec_line(oc, L, PY_A(1.0), R, PY_A(1.0), grid);

    // ---------------- bottom: the characteristic fan ----------------
    // Each seed xi gives a STRAIGHT line x = xi + lambda(alpha0(xi)) t.
    // Where two of them meet, the chart has torn.
    const int NC = 150;
    for (int i = 0; i < NC; ++i) {
        double xi = (double)i / NC;
        double v = lam(a0(xi));
        // Plot in the frame that moves with the MEAN characteristic speed.
        // Without this every line wraps the periodic box several times and the
        // panel turns into a hash in which no crossing can be seen.
        int prevx = PX(xi), prevy = PY_T(0.0);
        for (int j = 1; j <= 80; ++j) {
            double t = TVIEW * j / 80.0;
            double x = xi + (v - LAM_BAR) * t;
            x -= std::floor(x);                     // periodic
            int cx = PX(x), cy = PY_T(t);
            if (std::fabs(cx - prevx) < (R - L) / 2)
                olivec_line(oc, prevx, prevy, cx, cy, (i % 5 == 0) ? dim : grid);
            prevx = cx; prevy = cy;
        }
    }
    // shock time line
    olivec_line(oc, L, PY_T(T_SHOCK), R, PY_T(T_SHOCK), red);
    olivec_text(oc, "t shock", L + 6, PY_T(T_SHOCK) - 20, ft, 2, red);
    // current time line
    olivec_line(oc, L, PY_T(g_t), R, PY_T(g_t), amber);

    // ---------------- top: alpha profile, drawn parametrically ----------------
    // (x(xi,t), alpha0(xi)) -- after t_shock this curve FOLDS OVER and the
    // profile becomes multivalued. The fold is the shock.
    // faint initial profile
    for (int i = 1; i <= 600; ++i) {
        double x0 = (double)(i - 1) / 600.0, x1 = (double)i / 600.0;
        olivec_line(oc, PX(x0), PY_A(a0(x0)), PX(x1), PY_A(a0(x1)), grid);
    }
    const int NS = 3000;
    int px_prev = 0, py_prev = 0; bool have = false;
    bool folded = false;
    for (int i = 0; i <= NS; ++i) {
        double xi = (double)i / NS;
        double a = a0(xi);
        double x = xi + (lam(a) - LAM_BAR) * g_t;
        x -= std::floor(x);
        double st = stretch(xi, g_t);
        if (st < 0.0) folded = true;
        int cx = PX(x), cy = PY_A(a);
        uint32_t col = (st < 0.0) ? red : cyan;
        if (have && std::fabs(cx - px_prev) < (R - L) / 2) {
            olivec_line(oc, px_prev, py_prev, cx, cy, col);
            olivec_line(oc, px_prev, py_prev + 1, cx, cy + 1, col);
        }
        px_prev = cx; py_prev = cy; have = true;
    }

    // ---------------- HUD ----------------
    const char* names[3] = { "harmonic   f = 1",
                             "1+log      f = 2/alpha",
                             "shock-avoid f = 1 + k/alpha^2" };
    // olive.c default font has ONLY lowercase, digits and , - .  Everything
    // else (uppercase, = + / parentheses) renders as blank. Labels obey that.
    olivec_text(oc, "gauge shock - spacetime is fine. the chart tears.", 14, 8, ft, 2,
                rgba(140, 230, 210, 1.f));
    std::snprintf(buf, sizeof buf, "slicing %s", names[g_slice]);
    olivec_text(oc, buf, 14, 36, ft, 2, cyan);
    std::snprintf(buf, sizeof buf, "amp %.3f   kappa %.1f   dlambda-dalpha %.6f",
                  g_amp, g_kappa, dlam(1.0));
    olivec_text(oc, buf, 14, 64, ft, 2, rgba(110, 190, 175, 1.f));
    std::snprintf(buf, sizeof buf, "t %.4f    t shock %.4f    ratio %.3f",
                  g_t, T_SHOCK, g_t / T_SHOCK);
    olivec_text(oc, buf, 14, 92, ft, 2, amber);
    olivec_text(oc, "alpha of x     red means folded, multivalued", L, TOP_Y0 - 18, ft, 2,
                rgba(90, 150, 140, 1.f));
    olivec_text(oc, "characteristics x  xi plus lambda t     they cross at t shock",
                L, BOT_Y0 - 18, ft, 2, rgba(90, 150, 140, 1.f));
    if (folded)
        olivec_text(oc, "shock. characteristics crossed. alpha is multivalued now.",
                    L, FH - 24, ft, 2, red);
    else
        olivec_text(oc, "smooth. no crossing yet.", L, FH - 24, ft, 2, rgba(0, 200, 160, 1.f));
    return (uint8_t*)px.data();
}

}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
int main(int argc, char** argv) {
    sim_init(0, 0);
    g_amp = 0.20;
    printf("alpha0 = 1 + %.2f sin(2 pi x)\n\n", g_amp);
    printf("%-30s %-14s %-14s %s\n", "slicing", "lambda'(1)", "t_shock", "analytic");
    struct { int s; double k; const char* n; double an; } T[] = {
        {0, 0.0,   "harmonic      f=1",        0.795775},
        {1, 0.0,   "1+log         f=2/alpha",  1.125395},
        {2, 1.0,   "shock-avoid   k=1",        1.125395},
        {2, 10.0,  "shock-avoid   k=10",       2.639323},
        {2, 100.0, "shock-avoid   k=100",      7.997435},
    };
    for (auto& e : T) {
        g_slice = e.s; g_kappa = e.k;
        double ts = shock_time();
        printf("%-30s %-14.6f %-14.6f %.6f  ratio=%.6f\n", e.n, dlam(1.0), ts, e.an, ts / e.an);
    }
    // verify: at t = t_shock the stretch dx/dxi first touches zero
    g_slice = 1; g_kappa = 10.0; T_SHOCK = shock_time();
    double mn = 1e30; for (int i = 0; i < 20000; ++i) {
        double st = stretch((double)i / 20000.0, T_SHOCK); if (st < mn) mn = st;
    }
    printf("\n1+log: min dx/dxi at t=t_shock = %.3e   (should be ~0)\n", mn);
    double mn2 = 1e30; for (int i = 0; i < 20000; ++i) {
        double st = stretch((double)i / 20000.0, 0.98 * T_SHOCK); if (st < mn2) mn2 = st;
    }
    printf("1+log: min dx/dxi at t=0.98 t_shock = %.5f  (should be > 0)\n", mn2);
    int steps = argc > 1 ? atoi(argv[1]) : 0;
    g_t = (steps ? 1.12 : 1.0) * T_SHOCK;
    uint8_t* p = sim_render();
    stbi_write_png("gaugeshock_os_preview.png", FW, FH, 4, p, FW * 4);
    printf("wrote gaugeshock_os_preview.png  (t = %.4f = %.2f t_shock)\n", g_t, g_t / T_SHOCK);
    return 0;
}
#endif
