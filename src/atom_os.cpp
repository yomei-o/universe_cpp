// ============================================================================
// Universe OS - Atom / Nucleus / Quark  (C++ / WASM)
//
// Three physical levels, each integrated by its OWN derived difference equation.
// Nothing is scripted, nothing is drawn from a canned animation: every dot on
// screen is a degree of freedom whose update rule is written out below.
//
//   L1  electrons around the nuclei      Variational Monte Carlo, |Psi_T|^2
//   L2  protons + neutrons in a nucleus  Quantum Molecular Dynamics
//   L3  up/down quarks inside a nucleon  relativistic Cornell + Y-string
//
// Why L1 is Monte Carlo and not "classical electrons in a Coulomb field":
// a classical point electron in a 1/r potential has NO shell structure, no
// Pauli exclusion and no s/p lobes; it spirals in.  The honest particle-based
// way to get the real orbitals is to let point walkers sample the many-body
// probability density |Psi(r_1..r_N)|^2.  That is what quantum chemistry
// actually does (VMC/DMC), and the antisymmetry of the Slater determinant is
// what produces shells, Hund's rule and the p-lobes.
//
// ---------------------------------------------------------------------------
// DERIVATION 1 - electrons: the Metropolis drift-diffusion difference equation
// ---------------------------------------------------------------------------
// Trial function (atomic units hbar = m_e = e = 1):
//     Psi_T(R) = det[ Phi_j(r_i) ]_up * det[ Phi_j(r_i) ]_dn * exp(J)
//     Phi_j    = sum_k c_jk chi_k          (LCAO molecular orbital)
//     chi      = hydrogenic Slater orbital on a nucleus, exponent zeta
//     J        = sum_{i<j} a_ij u(r_ij),  u(r) = r/(1+b r)
//     a_ij     = 1/2 (antiparallel), 1/4 (parallel)   <- e-e cusp conditions
//
// We want walkers distributed as f(R) = |Psi_T(R)|^2.  Take the Fokker-Planck
// equation whose stationary solution is exactly that:
//     df/dt = D grad.( grad f - F f ),      F = grad ln|Psi_T|^2 = 2 grad Psi/Psi
// The Ito discretisation of the matching Langevin equation, for ONE electron i
// moved per step, is the difference equation
//
//     r_i^{n+1} = r_i^n + D tau_i F_i(R^n) + sqrt(2 D tau_i) xi ,   xi ~ N(0,1)^3
//
// with D = 1/2.  A finite tau makes that biased, so we correct it exactly with
// a Metropolis-Hastings accept/reject on the proposal density
//     G(r->r') = (4 pi D tau)^{-3/2} exp( -|r' - r - D tau F(r)|^2 / (4 D tau) )
//     A = min( 1, |Psi(R')|^2 G(r'->r) / ( |Psi(R)|^2 G(r->r') ) )
// which makes the sampled distribution |Psi_T|^2 for ANY tau (zero time-step
// error in the distribution).  tau_i is allowed to depend on r_i (small steps
// near a nucleus, large ones in the valence region); MH stays exact as long as
// the reverse move uses tau(r_i').  That is why both normalisations are kept.
//
// The Psi ratio for a single-electron move is O(N) using the inverse Slater
// matrix Ainv (A_ij = Phi_j(r_i)):
//     ratio_det = sum_j Phi_j(r_i') * Ainv[j][i]
//     grad ln det |_{r_i'} = ( sum_j Ainv[j][i] grad Phi_j(r_i') ) / ratio_det
//
// Energy (the observable that proves the sampling is right):
//     E_L(R) = -1/2 sum_i (lap Psi/Psi)_i + V(R)
//     (lap Psi/Psi)_i = (lap det/det)_i + lap_i J + |grad_i J|^2
//                       + 2 (grad_i ln det).(grad_i J)
//     (lap det/det)_i = sum_j Ainv[j][i] lap Phi_j(r_i)
//     V = -sum_A Z_A/|r_i-R_A| + sum_{i<j} 1/r_ij + sum_{A<B} Z_A Z_B/R_AB
// <E_L> over the walkers is the variational energy; it is printed next to the
// exact non-relativistic value so you can see how good the trial state is.
//
// ---------------------------------------------------------------------------
// DERIVATION 2 - nucleons: QMD difference equation
// ---------------------------------------------------------------------------
// Each nucleon is a Gaussian wave packet of fixed width L, centred at R_i with
// mean momentum P_i.  Taking the expectation of the effective NN Hamiltonian in
// that product state gives a classical Hamiltonian H(R,P) (all terms below are
// the analytic Gaussian folds):
//     rho_ij  = (4 pi L)^{-3/2} exp( -r_ij^2 / (4 L) )        (overlap density)
//     H = sum_i P_i^2/(2M)
//       + (alpha/(2 rho0)) sum_{i!=j} rho_ij                   Skyrme 2-body
//       + (beta/((g+1) rho0^g)) sum_i ( sum_{j!=i} rho_ij )^g  Skyrme density
//       + (Cs/(2 rho0)) sum_{i!=j} tau_i tau_j rho_ij          symmetry energy
//       + (e^2/2) sum_{i!=j, p p} erf(r_ij/(2 sqrt L))/r_ij    Coulomb
//       + Cp sum_{i<j} d_tau d_spin exp(-r_ij^2/(2 q0^2) - p_ij^2/(2 p0^2))
//                                                              Pauli potential
// The Pauli term is what replaces antisymmetry in a classical nucleon MD: it
// keeps identical nucleons out of the same phase-space cell and supplies the
// Fermi motion, so the nucleus neither collapses nor freezes.
// Because H depends on P beyond the kinetic term, dR/dt = dH/dP is NOT just
// P/M, so the usual kick-drift-kick leapfrog is not applicable.  The derived
// update is the explicit midpoint (2nd-order Runge-Kutta) difference equation
//     R* = R^n + (dt/2) dH/dP(R^n,P^n) ,   P* = P^n - (dt/2) dH/dR(R^n,P^n)
//     R^{n+1} = R^n + dt dH/dP(R*,P*) ,    P^{n+1} = P^n - dt dH/dR(R*,P*)
// Analytic derivatives used (with d rho_ij/dR_i = -rho_ij (R_i-R_j)/(2L)):
//     dE_sk2/dR_i = (alpha/rho0) sum_j d rho_ij/dR_i
//     dE_sk3/dR_i = (beta g/((g+1) rho0^g)) sum_j (rb_i^{g-1}+rb_j^{g-1}) drho/dR_i
//     d/dr [erf(r/a)/r] = 2 exp(-r^2/a^2)/(a sqrt(pi) r) - erf(r/a)/r^2
//     dVp/dR_i = -Vp (R_i-R_j)/q0^2 ,  dVp/dP_i = -Vp (P_i-P_j)/p0^2
// The ground state is prepared by frictional cooling of the same difference
// equation (P <- P (1-lambda) plus random kicks), i.e. by the model itself.
//
// ---------------------------------------------------------------------------
// DERIVATION 3 - quarks: relativistic leapfrog with a Y-string
// ---------------------------------------------------------------------------
//     H = sum_i sqrt(p_i^2 + m_i^2) + V(x)
//     V = -kappa sum_{i<j} erf(r_ij/r0)/r_ij + sigma * L_Y(x) + V0
//     kappa = (2/3) alpha_s hbar c   (colour factor of a 3q singlet)
//     L_Y   = minimal total string length = Steiner (Fermat) tree of 3 points
// Because H = T(p) + V(x) is separable, the exactly symplectic kick-drift-kick
// difference equation applies and conserves H to O(dt^2) with no drift:
//     p^{n+1/2} = p^n       - (dt/2) grad V(x^n)
//     x^{n+1}   = x^n       + dt * p^{n+1/2} / sqrt((p^{n+1/2})^2 + m^2)
//     p^{n+1}   = p^{n+1/2} - (dt/2) grad V(x^{n+1})
// (dx/dt = p/E is exactly the relativistic velocity v/c.)
// The string gradient uses the envelope theorem: because the junction S is at
// the minimum of L, dL/dx_i = (x_i - S)/|x_i - S|, a unit vector.  Confinement
// is therefore a constant inward force sigma no matter how far you pull: the
// quarks can never leave.  That is not imposed, it is what sigma*L means.
//
// Common WASM ABI (identical to every other Universe OS sim):
//   sim_init(hintW,hintH) ; sim_w() ; sim_h() ; sim_reset() ; sim_step(n)
//   sim_render() -> RGBA* ; sim_click(nx,ny) ; sim_set(id,val) ; sim_action(id)
// ============================================================================
#define OLIVEC_IMPLEMENTATION
#include "olive.c"
#include <vector>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <algorithm>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define KEEP EMSCRIPTEN_KEEPALIVE
#else
#define KEEP
#endif

// ---------------------------------------------------------------------------
// framebuffer geometry
// ---------------------------------------------------------------------------
static const int FW = 1000, FH = 680;

// ---------------------------------------------------------------------------
// random numbers (xorshift64* + Box-Muller), fully deterministic per reset
// ---------------------------------------------------------------------------
static uint64_t g_rs = 88172645463325252ull;
static inline double urand()
{
    g_rs ^= g_rs >> 12; g_rs ^= g_rs << 25; g_rs ^= g_rs >> 27;
    return (double)((g_rs * 2685821657736338717ull) >> 11) * (1.0 / 9007199254740992.0);
}
static double g_gcache = 0.0; static bool g_ghas = false;
static inline double nrand()
{
    if (g_ghas) { g_ghas = false; return g_gcache; }
    double u1 = urand(), u2 = urand();
    if (u1 < 1e-300) u1 = 1e-300;
    double r = std::sqrt(-2.0 * std::log(u1)), th = 2.0 * M_PI * u2;
    g_gcache = r * std::sin(th); g_ghas = true;
    return r * std::cos(th);
}

// ===========================================================================
//                        S Y S T E M   T A B L E
// ===========================================================================
// Element data. Eref = exact non-relativistic total electronic energy (Ha)
// [Chakravorty et al. PRA 47, 3649 (1993)] - used only as a printed reference.
// BA = experimental binding energy per nucleon (MeV), Rch = charge radius (fm).
struct ElemInfo {
    const char* sym; int A; double Eref; double BA; double Rch;
};
static const ElemInfo ELEM[19] = {
    {"",  0,    0.0,      0.000, 0.000},
    {"H",  1,   -0.5,     0.000, 0.8783},
    {"He", 4,   -2.9037,  7.074, 1.6755},
    {"Li", 7,   -7.4781,  5.606, 2.4440},
    {"Be", 9,  -14.6674,  6.463, 2.5190},
    {"B",  11, -24.6539,  6.928, 2.4060},
    {"C",  12, -37.8450,  7.680, 2.4702},
    {"N",  14, -54.5892,  7.476, 2.5582},
    {"O",  16, -75.0673,  7.976, 2.6991},
    {"F",  19, -99.7339,  7.779, 2.8976},
    {"Ne", 20,-128.9376,  8.032, 3.0055},
    {"Na", 23,-162.2546,  8.111, 2.9936},
    {"Mg", 24,-200.0530,  8.261, 3.0570},
    {"Al", 27,-242.3460,  8.332, 3.0610},
    {"Si", 28,-289.3590,  8.448, 3.1224},
    {"P",  31,-341.2590,  8.481, 3.1889},
    {"S",  32,-398.1100,  8.493, 3.2611},
    {"Cl", 35,-460.1480,  8.520, 3.3654},
    {"Ar", 40,-527.5400,  8.595, 3.4274},
};

static const double ANG = 1.8897259886;   // 1 Angstrom in bohr

// A "system" is one atom or one molecule at its experimental geometry.
struct SysDef {
    const char* name;      // ascii label (drawn with the built-in font)
    int  nat;
    int  Z[8];
    double R[8][3];        // bohr, centre of mass at origin-ish
    double Eref;           // reference total energy in Ha (0 = unknown)
};
static std::vector<SysDef> SYS;

static void push_atom(int Z)
{
    SysDef s; std::memset(&s, 0, sizeof(s));
    s.name = ELEM[Z].sym; s.nat = 1; s.Z[0] = Z;
    s.R[0][0] = s.R[0][1] = s.R[0][2] = 0.0;
    s.Eref = ELEM[Z].Eref;
    SYS.push_back(s);
}
static void push_mol(const char* nm, double Eref, int n,
                     const int* Zs, const double* xyz_ang)
{
    SysDef s; std::memset(&s, 0, sizeof(s));
    s.name = nm; s.nat = n; s.Eref = Eref;
    double cx = 0, cy = 0, cz = 0, w = 0;
    for (int i = 0; i < n; ++i) {
        s.Z[i] = Zs[i];
        for (int k = 0; k < 3; ++k) s.R[i][k] = xyz_ang[i * 3 + k] * ANG;
        double m = ELEM[Zs[i]].A;
        cx += m * s.R[i][0]; cy += m * s.R[i][1]; cz += m * s.R[i][2]; w += m;
    }
    for (int i = 0; i < n; ++i) { s.R[i][0] -= cx / w; s.R[i][1] -= cy / w; s.R[i][2] -= cz / w; }
    SYS.push_back(s);
}

static void build_systems()
{
    if (!SYS.empty()) return;
    for (int Z = 1; Z <= 18; ++Z) push_atom(Z);

    // Experimental geometries (bond lengths in Angstrom).  Eref = exact
    // non-relativistic total energy where a reliable value is available.
    { const int z[2] = {1,1};
      const double p[6] = {0,0,-0.3705, 0,0,0.3705};
      push_mol("H2", -1.1745, 2, z, p); }                       // r = 0.7414
    { const int z[3] = {8,1,1};
      // r(OH)=0.9575, angle=104.51 deg, molecule in the yz plane
      const double a = 104.51 * M_PI / 180.0 * 0.5, r = 0.9575;
      const double p[9] = {0,0,0,
                           0, r*std::sin(a),  r*std::cos(a),
                           0,-r*std::sin(a),  r*std::cos(a)};
      push_mol("H2O", -76.438, 3, z, p); }
    { const int z[2] = {8,8};
      const double p[6] = {0,0,-0.604, 0,0,0.604};
      push_mol("O2", -150.327, 2, z, p); }                      // r = 1.2075
    { const int z[2] = {7,7};
      const double p[6] = {0,0,-0.5488, 0,0,0.5488};
      push_mol("N2", -109.542, 2, z, p); }                      // r = 1.0977
    { const int z[5] = {6,1,1,1,1};
      const double d = 1.0870 / std::sqrt(3.0);
      const double p[15] = {0,0,0,  d,d,d,  d,-d,-d,  -d,d,-d,  -d,-d,d};
      push_mol("CH4", -40.515, 5, z, p); }
    { const int z[4] = {7,1,1,1};
      // r(NH)=1.0124, HNH=106.67 deg  -> polar angle from C3 axis = 111.7 deg
      const double r = 1.0124, th = 111.70 * M_PI / 180.0;
      const double st = r * std::sin(th), ct = r * std::cos(th);
      const double p[12] = {0,0,0,
                            st, 0, ct,
                            st * -0.5,  st * 0.8660254, ct,
                            st * -0.5,  st * -0.8660254, ct};
      push_mol("NH3", -56.564, 4, z, p); }
    { const int z[2] = {1,9};
      const double p[6] = {0,0,-0.8168, 0,0,0.0999};
      push_mol("HF", -100.459, 2, z, p); }                      // r = 0.9168
    { const int z[2] = {6,8};
      const double p[6] = {0,0,-0.6446, 0,0,0.4834};
      push_mol("CO", -113.317, 2, z, p); }                      // r = 1.128
}

// ===========================================================================
//         L E V E L   1  :   e l e c t r o n s   ( V M C )
// ===========================================================================
static const int MAXAT = 8;
static const int MAXAO = 40;
static const int MAXMO = 24;
static const int MAXEL = 24;
static const int MAXW  = 64;
static const int NTRAIL = 26;

struct Atom { int Z, A; double R[3]; };

// atomic orbital: hydrogenic, exponent zeta, on centre `at`
//   l=0 : s          l=1 : m=0,1,2 -> x,y,z
struct AO { int at, n, l, m; double zeta; };

struct MO { double c[MAXAO]; };   // coefficients over the AO list

static std::vector<Atom> g_at;
static AO  g_ao[MAXAO];  static int g_nao = 0;
static MO  g_mo[MAXMO];  static int g_nmo = 0;
static int g_occup[MAXMO], g_occdn[MAXMO];  // MO index per column
static int g_nup = 0, g_ndn = 0, g_nel = 0;
static double g_jb = 1.0;                   // Jastrow b
static double g_zscale = 1.0;               // legacy global scale (kept for the view code)
static double g_zs[4]  = {1.0, 1.0, 1.0, 1.0};  // per-shell variational scale on zeta
static double g_zp     = 1.0;               // extra scale on the p orbitals
static double g_zh     = 1.0;               // extra scale on hydrogen 1s
static double g_hyb    = 1.0;               // scale on the hybrid s weight
static double g_enuc = 0.0;                 // nuclear repulsion (Ha)
static int    g_sysid = 7;                  // default: oxygen atom

// -- Slater screening (Slater's rules) gives the starting exponents ---------
static void shell_zetas(int Z, double z[4])
{
    int n1 = std::min(Z, 2);
    int n2 = std::max(0, std::min(Z - 2, 8));
    int n3 = std::max(0, std::min(Z - 10, 8));
    z[1] = Z - 0.30 * (n1 - 1);
    z[2] = Z - 0.85 * n1 - 0.35 * (n2 > 0 ? n2 - 1 : 0);
    z[3] = Z - 1.00 * (n1 + n2) - 0.35 * (n3 > 0 ? n3 - 1 : 0);
    for (int i = 1; i <= 3; ++i) if (z[i] < 0.5) z[i] = 0.5;
}

// -- hydrogenic radial polynomials -----------------------------------------
// s-type:  chi = P(r) exp(-lam r)            P as below, lam = zeta/n
// p-type:  chi = x_m * Q(r) exp(-lam r)
struct AOVal { double v, g[3], lap; };

static inline void ao_eval(const AO& a, const double* r, AOVal& o)
{
    double dx = r[0] - g_at[a.at].R[0];
    double dy = r[1] - g_at[a.at].R[1];
    double dz = r[2] - g_at[a.at].R[2];
    double rr = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (rr < 1e-8) rr = 1e-8;
    double lam = a.zeta / a.n;
    double u = std::exp(-lam * rr);
    if (a.l == 0) {
        double P, P1, P2;
        if (a.n == 1)      { P = 1.0;                                   P1 = 0.0;                    P2 = 0.0; }
        else if (a.n == 2) { P = 1.0 - lam * rr;                        P1 = -lam;                   P2 = 0.0; }
        else               { double t = lam * rr;
                             P = 1.0 - 2.0 * t + (2.0 / 3.0) * t * t;   P1 = lam * (-2.0 + (4.0 / 3.0) * t);
                             P2 = (4.0 / 3.0) * lam * lam; }
        double g0  = P * u;
        double g1  = (P1 - lam * P) * u;
        double g2  = (P2 - 2.0 * lam * P1 + lam * lam * P) * u;
        o.v = g0;
        o.g[0] = g1 * dx / rr; o.g[1] = g1 * dy / rr; o.g[2] = g1 * dz / rr;
        o.lap = g2 + 2.0 * g1 / rr;
    } else {
        double Q, Q1, Q2;
        if (a.n == 2) { Q = 1.0;                Q1 = 0.0;        Q2 = 0.0; }
        else          { Q = 1.0 - 0.5 * lam * rr; Q1 = -0.5 * lam; Q2 = 0.0; }
        double h0 = Q * u;
        double h1 = (Q1 - lam * Q) * u;
        double h2 = (Q2 - 2.0 * lam * Q1 + lam * lam * Q) * u;
        double d[3] = {dx, dy, dz};
        double c = d[a.m];
        o.v = c * h0;
        for (int k = 0; k < 3; ++k)
            o.g[k] = (k == a.m ? h0 : 0.0) + c * h1 * d[k] / rr;
        o.lap = c * (h2 + 4.0 * h1 / rr);
    }
}

// molecular orbital = sum of AOs
static inline void mo_eval(int j, const double* r, AOVal& o)
{
    o.v = 0.0; o.g[0] = o.g[1] = o.g[2] = 0.0; o.lap = 0.0;
    const MO& m = g_mo[j];
    for (int k = 0; k < g_nao; ++k) {
        double c = m.c[k];
        if (c == 0.0) continue;
        AOVal t; ao_eval(g_ao[k], r, t);
        o.v += c * t.v; o.lap += c * t.lap;
        o.g[0] += c * t.g[0]; o.g[1] += c * t.g[1]; o.g[2] += c * t.g[2];
    }
}

// -- small dense linear algebra (n <= MAXMO) -------------------------------
// Gauss-Jordan inverse with partial pivoting; returns false if singular.
static bool mat_inv(const double* A, double* Ai, int n)
{
    double M[MAXMO][2 * MAXMO];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) { M[i][j] = A[i * MAXMO + j]; M[i][n + j] = (i == j) ? 1.0 : 0.0; }
    }
    for (int c = 0; c < n; ++c) {
        int p = c; double best = std::fabs(M[c][c]);
        for (int i = c + 1; i < n; ++i) { double v = std::fabs(M[i][c]); if (v > best) { best = v; p = i; } }
        if (best < 1e-300) return false;
        if (p != c) for (int j = 0; j < 2 * n; ++j) std::swap(M[c][j], M[p][j]);
        double iv = 1.0 / M[c][c];
        for (int j = 0; j < 2 * n; ++j) M[c][j] *= iv;
        for (int i = 0; i < n; ++i) {
            if (i == c) continue;
            double f = M[i][c];
            if (f == 0.0) continue;
            for (int j = 0; j < 2 * n; ++j) M[i][j] -= f * M[c][j];
        }
    }
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) Ai[i * MAXMO + j] = M[i][n + j];
    return true;
}

struct Walker {
    double r[MAXEL][3];
    double Au[MAXMO * MAXMO], Aui[MAXMO * MAXMO];   // A[i*MAXMO+j] = Phi_j(r_i)
    double Ad[MAXMO * MAXMO], Adi[MAXMO * MAXMO];
    bool   ok;
    double eloc;
    float  tx[MAXEL][NTRAIL], ty[MAXEL][NTRAIL], tz[MAXEL][NTRAIL];
    int    tn;
};
static Walker g_w[MAXW];
static int    g_nw = 40;

// ---------------------------------------------------------------------------
// energy statistics.  Successive VMC samples are correlated, so the naive
// sqrt(var/N) badly underestimates the error.  We therefore average in blocks
// and take the scatter OF THE BLOCK MEANS, which is what gets displayed.
// ---------------------------------------------------------------------------
static const long EBLK = 240;
static double g_esum = 0.0, g_esum2 = 0.0; static long g_ecnt = 0;
static double g_eblk = 0.0; static long g_eblkn = 0;
static double g_bsum = 0.0, g_bsum2 = 0.0; static long g_bcnt = 0;
static double g_acc = 0.0, g_att = 0.0;

static inline void stat_add(double e)
{
    g_esum += e; g_esum2 += e * e; g_ecnt++;
    g_eblk += e; g_eblkn++;
    if (g_eblkn >= EBLK) {
        double m = g_eblk / g_eblkn;
        g_bsum += m; g_bsum2 += m * m; g_bcnt++;
        g_eblk = 0.0; g_eblkn = 0;
    }
}
static inline void stat_clear()
{
    g_esum = g_esum2 = 0.0; g_ecnt = 0;
    g_eblk = 0.0; g_eblkn = 0;
    g_bsum = g_bsum2 = 0.0; g_bcnt = 0;
    g_acc = g_att = 0.0;
}
static inline double stat_mean() { return g_ecnt ? g_esum / g_ecnt : 0.0; }
static inline double stat_err()
{
    if (g_bcnt < 3) return 0.0;
    double m = g_bsum / g_bcnt;
    double v = g_bsum2 / g_bcnt - m * m;
    return v > 0 ? std::sqrt(v / (g_bcnt - 1)) : 0.0;
}

static inline int spin_of(int i) { return i < g_nup ? 0 : 1; }

static void build_slater(Walker& w)
{
    w.ok = true;
    for (int i = 0; i < g_nup; ++i)
        for (int j = 0; j < g_nup; ++j) { AOVal t; mo_eval(g_occup[j], w.r[i], t); w.Au[i * MAXMO + j] = t.v; }
    for (int i = 0; i < g_ndn; ++i)
        for (int j = 0; j < g_ndn; ++j) { AOVal t; mo_eval(g_occdn[j], w.r[g_nup + i], t); w.Ad[i * MAXMO + j] = t.v; }
    if (g_nup && !mat_inv(w.Au, w.Aui, g_nup)) w.ok = false;
    if (g_ndn && !mat_inv(w.Ad, w.Adi, g_ndn)) w.ok = false;
}

// Jastrow pair coefficient: 1/2 antiparallel, 1/4 parallel (cusp conditions)
static inline double jcoef(int i, int j) { return spin_of(i) == spin_of(j) ? 0.25 : 0.5; }
static inline double jU(double r) { return r / (1.0 + g_jb * r); }

// sum over pairs involving electron i, with r_i replaced by rp
static double jastrow_part(const Walker& w, int i, const double* rp)
{
    double s = 0.0;
    for (int j = 0; j < g_nel; ++j) {
        if (j == i) continue;
        double dx = rp[0] - w.r[j][0], dy = rp[1] - w.r[j][1], dz = rp[2] - w.r[j][2];
        s += jcoef(i, j) * jU(std::sqrt(dx * dx + dy * dy + dz * dz));
    }
    return s;
}

// grad_i J and lap_i J at position rp
static void jastrow_deriv(const Walker& w, int i, const double* rp, double* gr, double* lp)
{
    gr[0] = gr[1] = gr[2] = 0.0; double L = 0.0;
    for (int j = 0; j < g_nel; ++j) {
        if (j == i) continue;
        double d[3] = {rp[0] - w.r[j][0], rp[1] - w.r[j][1], rp[2] - w.r[j][2]};
        double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (r < 1e-10) r = 1e-10;
        double a = jcoef(i, j), q = 1.0 + g_jb * r;
        double u1 = 1.0 / (q * q);              // u'
        double u2 = -2.0 * g_jb / (q * q * q);  // u''
        for (int k = 0; k < 3; ++k) gr[k] += a * u1 * d[k] / r;
        L += a * (u2 + 2.0 * u1 / r);
    }
    *lp = L;
}

// position dependent time step: small near a nucleus, tau0 far away
static double g_tau0 = 0.25;
static inline double tau_at(const double* r)
{
    double best = 1e30; int ib = 0;
    for (size_t A = 0; A < g_at.size(); ++A) {
        double dx = r[0] - g_at[A].R[0], dy = r[1] - g_at[A].R[1], dz = r[2] - g_at[A].R[2];
        double d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < best) { best = d2; ib = (int)A; }
    }
    double a = 1.0 / std::max(1.0, (double)g_at[ib].Z);   // 1s length scale
    return g_tau0 * best / (a * a + best);
}

// determinant ratio + drift for a single-electron trial move
struct MoveInfo { double ratio; double gdet[3]; };
static bool det_ratio(const Walker& w, int i, const double* rp, MoveInfo& mi)
{
    int sp = spin_of(i);
    int n  = sp == 0 ? g_nup : g_ndn;
    int li = sp == 0 ? i : i - g_nup;
    const double* Ai = sp == 0 ? w.Aui : w.Adi;
    const int* occ = sp == 0 ? g_occup : g_occdn;
    double rt = 0.0, gd[3] = {0, 0, 0};
    for (int j = 0; j < n; ++j) {
        AOVal t; mo_eval(occ[j], rp, t);
        double a = Ai[j * MAXMO + li];
        rt += a * t.v;
        gd[0] += a * t.g[0]; gd[1] += a * t.g[1]; gd[2] += a * t.g[2];
    }
    if (std::fabs(rt) < 1e-300) return false;
    mi.ratio = rt;
    for (int k = 0; k < 3; ++k) mi.gdet[k] = gd[k] / rt;
    return true;
}

// one Metropolis drift-diffusion sweep (all electrons of one walker)
static void vmc_sweep(Walker& w)
{
    if (!w.ok) { build_slater(w); if (!w.ok) return; }
    const double D = 0.5;
    for (int i = 0; i < g_nel; ++i) {
        double* r = w.r[i];
        MoveInfo m0;
        if (!det_ratio(w, i, r, m0)) { build_slater(w); if (!w.ok) return; if (!det_ratio(w, i, r, m0)) return; }
        double gj0[3], lj0; jastrow_deriv(w, i, r, gj0, &lj0);
        double F0[3];  // grad ln|Psi|^2 = 2 grad ln Psi
        for (int k = 0; k < 3; ++k) F0[k] = 2.0 * (m0.gdet[k] + gj0[k]);
        double t0 = tau_at(r), s0 = std::sqrt(2.0 * D * t0);
        double mu0[3], step = 0.0;
        for (int k = 0; k < 3; ++k) { mu0[k] = D * t0 * F0[k]; step += mu0[k] * mu0[k]; }
        double cap = 2.0 * s0;   // deterministic drift cap (keeps MH exact)
        step = std::sqrt(step);
        if (step > cap) for (int k = 0; k < 3; ++k) mu0[k] *= cap / step;
        double rp[3];
        for (int k = 0; k < 3; ++k) rp[k] = r[k] + mu0[k] + s0 * nrand();

        MoveInfo m1;
        g_att += 1.0;
        if (!det_ratio(w, i, rp, m1)) continue;
        double gj1[3], lj1; jastrow_deriv(w, i, rp, gj1, &lj1);
        double dJ = jastrow_part(w, i, rp) - jastrow_part(w, i, r);
        double psir = m1.ratio * std::exp(dJ);

        double F1[3];
        for (int k = 0; k < 3; ++k) F1[k] = 2.0 * (m1.gdet[k] + gj1[k]);
        double t1 = tau_at(rp), s1 = std::sqrt(2.0 * D * t1);
        double mu1[3], st1 = 0.0;
        for (int k = 0; k < 3; ++k) { mu1[k] = D * t1 * F1[k]; st1 += mu1[k] * mu1[k]; }
        double cap1 = 2.0 * s1; st1 = std::sqrt(st1);
        if (st1 > cap1) for (int k = 0; k < 3; ++k) mu1[k] *= cap1 / st1;

        // ln G(r'->r) - ln G(r->r')   (normalisations kept: tau differs)
        double q0 = 0.0, q1 = 0.0;
        for (int k = 0; k < 3; ++k) {
            double a = rp[k] - r[k] - mu0[k];      q0 += a * a;   // forward
            double b = r[k] - rp[k] - mu1[k];      q1 += b * b;   // reverse
        }
        double lg = -3.0 * std::log(s1 / s0) - q1 / (2.0 * s1 * s1) + q0 / (2.0 * s0 * s0);
        double lnA = 2.0 * std::log(std::fabs(psir)) + lg;
        if (lnA >= 0.0 || urand() < std::exp(lnA)) {
            for (int k = 0; k < 3; ++k) r[k] = rp[k];
            // refresh the Slater matrix row and its inverse
            int sp = spin_of(i), n = sp == 0 ? g_nup : g_ndn, li = sp == 0 ? i : i - g_nup;
            double* A  = sp == 0 ? w.Au  : w.Ad;
            double* Ai = sp == 0 ? w.Aui : w.Adi;
            const int* occ = sp == 0 ? g_occup : g_occdn;
            for (int j = 0; j < n; ++j) { AOVal t; mo_eval(occ[j], r, t); A[li * MAXMO + j] = t.v; }
            if (!mat_inv(A, Ai, n)) { build_slater(w); }
            g_acc += 1.0;
        }
    }
}

// local energy of one walker
static double local_energy(Walker& w)
{
    if (!w.ok) return 0.0;
    double kin = 0.0, ven = 0.0, vee = 0.0;
    for (int i = 0; i < g_nel; ++i) {
        int sp = spin_of(i), n = sp == 0 ? g_nup : g_ndn, li = sp == 0 ? i : i - g_nup;
        const double* Ai = sp == 0 ? w.Aui : w.Adi;
        const int* occ = sp == 0 ? g_occup : g_occdn;
        double lapd = 0.0, gd[3] = {0, 0, 0};
        for (int j = 0; j < n; ++j) {
            AOVal t; mo_eval(occ[j], w.r[i], t);
            double a = Ai[j * MAXMO + li];
            lapd += a * t.lap;
            for (int k = 0; k < 3; ++k) gd[k] += a * t.g[k];
        }
        double gj[3], lj; jastrow_deriv(w, i, w.r[i], gj, &lj);
        double gj2 = gj[0] * gj[0] + gj[1] * gj[1] + gj[2] * gj[2];
        double cross = 2.0 * (gd[0] * gj[0] + gd[1] * gj[1] + gd[2] * gj[2]);
        kin += -0.5 * (lapd + lj + gj2 + cross);
        for (size_t A = 0; A < g_at.size(); ++A) {
            double dx = w.r[i][0] - g_at[A].R[0], dy = w.r[i][1] - g_at[A].R[1], dz = w.r[i][2] - g_at[A].R[2];
            double d = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (d < 1e-10) d = 1e-10;
            ven += -(double)g_at[A].Z / d;
        }
        for (int j = i + 1; j < g_nel; ++j) {
            double dx = w.r[i][0] - w.r[j][0], dy = w.r[i][1] - w.r[j][1], dz = w.r[i][2] - w.r[j][2];
            double d = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (d < 1e-10) d = 1e-10;
            vee += 1.0 / d;
        }
    }
    return kin + ven + vee + g_enuc;
}

// ---------------------------------------------------------------------------
// build the AO list / MOs / occupations for the selected system
// ---------------------------------------------------------------------------
// subshell layout inside one atom's AO block: 1s 2s 2px 2py 2pz 3s 3px 3py 3pz
static const int SS_START[5] = {0, 1, 2, 5, 6};
static const int SS_COUNT[5] = {1, 1, 3, 1, 3};
static const int SS_N[5]     = {1, 2, 2, 3, 3};
static const int SS_L[5]     = {0, 0, 1, 0, 1};

static int  g_aobase[MAXAT];    // first AO index of atom a
static int  g_naoat[MAXAT];     // how many AOs atom a contributes
static double g_zt[MAXAT][4];
static double g_lam = 1.0;      // bond polarity, variational
static bool g_molmo = false;    // true when genuine LCAO molecular orbitals are in use
static bool build_mos_molecule();
static const char* g_moname[MAXMO];   // label per MO, for the readout

// determinant of bare atomic orbitals: correct for a single atom, and the
// fallback if the molecular builder cannot handle a geometry
static void build_mos_atomic()
{
    g_nmo = 0;
    for (int i = 0; i < g_nao; ++i) {
        std::memset(g_mo[g_nmo].c, 0, sizeof(double) * MAXAO);
        g_mo[g_nmo].c[i] = 1.0;
        g_moname[g_nmo] = "atomic";
        g_nmo++;
    }
    g_nup = g_ndn = 0;
    for (size_t a = 0; a < g_at.size(); ++a) {
        int rem = g_at[a].Z, base = g_aobase[a];
        for (int s = 0; s < 5 && rem > 0; ++s) {
            int k = std::min(rem, SS_COUNT[s]);
            for (int c = 0; c < k; ++c) g_occup[g_nup++] = base + SS_START[s] + c;
            rem -= k;
            int k2 = std::min(rem, SS_COUNT[s]);
            for (int c = 0; c < k2; ++c) g_occdn[g_ndn++] = base + SS_START[s] + c;
            rem -= k2;
        }
    }
    g_nel = g_nup + g_ndn;
}

static void setup_electrons(int sysid)
{
    build_systems();
    if (sysid < 0) sysid = 0;
    if (sysid >= (int)SYS.size()) sysid = (int)SYS.size() - 1;
    g_sysid = sysid;
    const SysDef& S = SYS[sysid];

    g_at.clear();
    for (int i = 0; i < S.nat; ++i) {
        Atom a; a.Z = S.Z[i]; a.A = ELEM[S.Z[i]].A;
        for (int k = 0; k < 3; ++k) a.R[k] = S.R[i][k];
        g_at.push_back(a);
    }
    // nuclear repulsion
    g_enuc = 0.0;
    for (size_t i = 0; i < g_at.size(); ++i)
        for (size_t j = i + 1; j < g_at.size(); ++j) {
            double dx = g_at[i].R[0] - g_at[j].R[0], dy = g_at[i].R[1] - g_at[j].R[1], dz = g_at[i].R[2] - g_at[j].R[2];
            g_enuc += (double)g_at[i].Z * g_at[j].Z / std::sqrt(dx * dx + dy * dy + dz * dz);
        }

    // ---- AO list: enough shells to hold each atom's own electrons ----
    g_nao = 0;
    for (size_t a = 0; a < g_at.size(); ++a) {
        int Z = g_at[a].Z;
        shell_zetas(Z, g_zt[a]);
        int nsh = (Z <= 2) ? 1 : (Z <= 10 ? 3 : 5);   // 1s | +2s2p | +3s3p
        g_aobase[a] = g_nao;
        for (int s = 0; s < nsh; ++s)
            for (int c = 0; c < SS_COUNT[s]; ++c) {
                AO& o = g_ao[g_nao++];
                o.at = (int)a; o.n = SS_N[s]; o.l = SS_L[s];
                o.m = (SS_L[s] == 1) ? c : 0;
                o.zeta = g_zt[a][SS_N[s]] * g_zs[SS_N[s]] * (SS_L[s] == 1 ? g_zp : 1.0)
                       * (Z == 1 ? g_zh : 1.0);
            }
        g_naoat[a] = g_nao - g_aobase[a];
    }

    build_mos_atomic();
    // a molecule needs LCAO orbitals; if the builder cannot handle the shape,
    // fall back to the atomic determinant (and say so on screen)
    g_molmo = false;
    if (g_at.size() > 1) {
        if (build_mos_molecule()) g_molmo = true;
        else                      build_mos_atomic();
    }
}

// ===========================================================================
//        M O L E C U L A R   O R B I T A L S
// ===========================================================================
// The determinant needs one-electron orbitals.  For a molecule those must be
// LCAO combinations, or there is no chemical bond at all: a determinant of bare
// atomic orbitals is just two atoms sitting next to each other, and its energy
// comes out ABOVE the separated atoms because of the nuclear repulsion.
//
// The combinations here are not fitted to anything - they are the ones symmetry
// forces, built from the geometry:
//   diatomic     sigma/pi bonding and antibonding pairs along the bond axis,
//                filled by Aufbau + Hund over the degenerate pi levels (which
//                is what makes O2 come out as a spin triplet, on its own)
//   AX_n         sp^m hybrids on the central atom pointing at each ligand, the
//                hybrid mixing fixed by the measured bond angle through the
//                orthogonality condition a^2 + b^2 cos(theta) = 0, plus the
//                orthogonal complement of those hybrids as the lone pairs
// The single free number is the bond polarity lambda, and that is optimised
// variationally like every other trial-function parameter.
// ===========================================================================
static inline void mo_clear(int j)
{
    std::memset(g_mo[j].c, 0, sizeof(double) * MAXAO);
    g_moname[j] = "mo";
}
static inline void mo_normalise(int j)
{
    double s = 0;
    for (int k = 0; k < g_nao; ++k) s += g_mo[j].c[k] * g_mo[j].c[k];
    if (s > 1e-300) { s = 1.0 / std::sqrt(s); for (int k = 0; k < g_nao; ++k) g_mo[j].c[k] *= s; }
}
// local AO offsets inside one atom's block
enum { L1S = 0, L2S = 1, L2PX = 2, L2PY = 3, L2PZ = 4, L3S = 5, L3PX = 6 };

// fill `ngroup` degeneracy groups with `nel` electrons, Hund's rule inside each
static bool fill_groups(const int (*grp)[3], const int* gsz, int ngroup, int nel)
{
    g_nup = g_ndn = 0;
    int rem = nel;
    for (int g = 0; g < ngroup && rem > 0; ++g) {
        int k = std::min(rem, gsz[g]);
        for (int c = 0; c < k; ++c) g_occup[g_nup++] = grp[g][c];
        rem -= k;
        int k2 = std::min(rem, gsz[g]);
        for (int c = 0; c < k2; ++c) g_occdn[g_ndn++] = grp[g][c];
        rem -= k2;
    }
    g_nel = g_nup + g_ndn;
    return rem == 0;
}

static bool build_mos_molecule()
{
    int nat = (int)g_at.size();
    int nel = 0;
    for (int a = 0; a < nat; ++a) nel += g_at[a].Z;
    double lam = g_lam;

    // ---------------- diatomic ----------------
    if (nat == 2) {
        int A = 0, B = 1;
        // put the atom with the more negative z first so the sigma/pi signs below hold
        if (g_at[0].R[2] > g_at[1].R[2]) { A = 1; B = 0; }
        int ba = g_aobase[A], bb = g_aobase[B];
        int na = g_naoat[A], nb = g_naoat[B];
        int grp[10][3]; int gsz[10], ng = 0, nm = 0;

        if (na == 1 && nb == 1) {                 // H2 : one sigma_g
            mo_clear(nm); g_mo[nm].c[ba + L1S] = 1.0; g_mo[nm].c[bb + L1S] = 1.0;
            g_moname[nm] = "sigma g"; mo_normalise(nm);
            grp[ng][0] = nm; gsz[ng] = 1; ng++; nm++;
        } else if (na == 1 || nb == 1) {          // H-X : core, 2s, sigma, 2 pi lone pairs
            int H = (na == 1) ? A : B, X = (na == 1) ? B : A;
            int bh = g_aobase[H], bx = g_aobase[X];
            double sgn = (g_at[H].R[2] > g_at[X].R[2]) ? 1.0 : -1.0;
            mo_clear(nm); g_mo[nm].c[bx + L1S] = 1.0; g_moname[nm] = "core 1s";
            grp[ng][0] = nm; gsz[ng] = 1; ng++; nm++;
            mo_clear(nm); g_mo[nm].c[bx + L2S] = 1.0; g_moname[nm] = "2s lone pair";
            grp[ng][0] = nm; gsz[ng] = 1; ng++; nm++;
            mo_clear(nm); g_mo[nm].c[bx + L2PZ] = sgn; g_mo[nm].c[bh + L1S] = lam;
            g_moname[nm] = "sigma bond"; mo_normalise(nm);
            grp[ng][0] = nm; gsz[ng] = 1; ng++; nm++;
            mo_clear(nm); g_mo[nm].c[bx + L2PX] = 1.0; g_moname[nm] = "pi lone pair";
            mo_clear(nm + 1); g_mo[nm + 1].c[bx + L2PY] = 1.0; g_moname[nm + 1] = "pi lone pair";
            grp[ng][0] = nm; grp[ng][1] = nm + 1; gsz[ng] = 2; ng++; nm += 2;
        } else {                                  // X-Y : the full sigma/pi ladder
            mo_clear(nm); g_mo[nm].c[ba + L1S] = 1.0; g_moname[nm] = "core 1s a";
            grp[ng][0] = nm; gsz[ng] = 1; ng++; nm++;
            mo_clear(nm); g_mo[nm].c[bb + L1S] = 1.0; g_moname[nm] = "core 1s b";
            grp[ng][0] = nm; gsz[ng] = 1; ng++; nm++;
            mo_clear(nm); g_mo[nm].c[ba + L2S] = 1.0; g_mo[nm].c[bb + L2S] = lam;
            g_moname[nm] = "sigma 2s"; mo_normalise(nm);
            grp[ng][0] = nm; gsz[ng] = 1; ng++; nm++;
            mo_clear(nm); g_mo[nm].c[ba + L2S] = 1.0; g_mo[nm].c[bb + L2S] = -lam;
            g_moname[nm] = "sigma 2s star"; mo_normalise(nm);
            grp[ng][0] = nm; gsz[ng] = 1; ng++; nm++;
            mo_clear(nm);     g_mo[nm].c[ba + L2PX] = 1.0;     g_mo[nm].c[bb + L2PX] = lam;
            mo_clear(nm + 1); g_mo[nm + 1].c[ba + L2PY] = 1.0; g_mo[nm + 1].c[bb + L2PY] = lam;
            g_moname[nm] = "pi bonding"; g_moname[nm + 1] = "pi bonding";
            mo_normalise(nm); mo_normalise(nm + 1);
            grp[ng][0] = nm; grp[ng][1] = nm + 1; gsz[ng] = 2; ng++; nm += 2;
            mo_clear(nm); g_mo[nm].c[ba + L2PZ] = 1.0; g_mo[nm].c[bb + L2PZ] = -lam;
            g_moname[nm] = "sigma 2p"; mo_normalise(nm);
            grp[ng][0] = nm; gsz[ng] = 1; ng++; nm++;
            mo_clear(nm);     g_mo[nm].c[ba + L2PX] = 1.0;     g_mo[nm].c[bb + L2PX] = -lam;
            mo_clear(nm + 1); g_mo[nm + 1].c[ba + L2PY] = 1.0; g_mo[nm + 1].c[bb + L2PY] = -lam;
            g_moname[nm] = "pi star"; g_moname[nm + 1] = "pi star";
            mo_normalise(nm); mo_normalise(nm + 1);
            grp[ng][0] = nm; grp[ng][1] = nm + 1; gsz[ng] = 2; ng++; nm += 2;
            mo_clear(nm); g_mo[nm].c[ba + L2PZ] = 1.0; g_mo[nm].c[bb + L2PZ] = lam;
            g_moname[nm] = "sigma 2p star"; mo_normalise(nm);
            grp[ng][0] = nm; gsz[ng] = 1; ng++; nm++;
        }
        g_nmo = nm;
        return fill_groups(grp, gsz, ng, nel);
    }

    // ---------------- central atom + n ligands ----------------
    int C = 0;
    for (int a = 1; a < nat; ++a) if (g_at[a].Z > g_at[C].Z) C = a;
    if (g_naoat[C] < 5) return false;
    int nb = 0, lig[6];
    for (int a = 0; a < nat; ++a) if (a != C) { if (nb >= 6) return false; lig[nb++] = a; }
    if (nb > 4) return false;
    for (int i = 0; i < nb; ++i) if (g_naoat[lig[i]] != 1) return false;   // H ligands only

    // bond unit vectors and the hybrid mixing from the measured bond angle
    double u[4][3];
    for (int i = 0; i < nb; ++i) {
        double d[3], n = 0;
        for (int k = 0; k < 3; ++k) { d[k] = g_at[lig[i]].R[k] - g_at[C].R[k]; n += d[k] * d[k]; }
        n = std::sqrt(n); if (n < 1e-9) return false;
        for (int k = 0; k < 3; ++k) u[i][k] = d[k] / n;
    }
    double ct = -1.0 / 3.0;
    if (nb >= 2) {
        double s = 0; int c = 0;
        for (int i = 0; i < nb; ++i)
            for (int j = i + 1; j < nb; ++j) {
                double d = 0; for (int k = 0; k < 3; ++k) d += u[i][k] * u[j][k];
                s += d; c++;
            }
        ct = s / c;
    }
    if (ct > -0.02) ct = -1.0 / 3.0;
    double a2 = -ct / (1.0 - ct) * g_hyb;         // s weight, from orthogonality
    if (a2 < 0.02) a2 = 0.02; if (a2 > 0.9) a2 = 0.9;
    double ah = std::sqrt(a2), bh = std::sqrt(1.0 - a2);

    // hybrids in the 4-dim space (s, px, py, pz), then their complement
    double V[4][4]; int nv = 0;
    for (int i = 0; i < nb; ++i) {
        V[nv][0] = ah;
        for (int k = 0; k < 3; ++k) V[nv][1 + k] = bh * u[i][k];
        nv++;
    }
    // Gram-Schmidt the hybrids, then grow the basis to 4 with the largest residuals
    for (int i = 0; i < nv; ++i) {
        for (int j = 0; j < i; ++j) {
            double d = 0; for (int k = 0; k < 4; ++k) d += V[i][k] * V[j][k];
            for (int k = 0; k < 4; ++k) V[i][k] -= d * V[j][k];
        }
        double n = 0; for (int k = 0; k < 4; ++k) n += V[i][k] * V[i][k];
        n = std::sqrt(n); if (n < 1e-9) return false;
        for (int k = 0; k < 4; ++k) V[i][k] /= n;
    }
    while (nv < 4) {
        double bestv[4] = {0, 0, 0, 0}, bestn = -1.0;
        for (int e = 0; e < 4; ++e) {
            double t[4] = {0, 0, 0, 0}; t[e] = 1.0;
            for (int j = 0; j < nv; ++j) {
                double d = 0; for (int k = 0; k < 4; ++k) d += t[k] * V[j][k];
                for (int k = 0; k < 4; ++k) t[k] -= d * V[j][k];
            }
            double n = 0; for (int k = 0; k < 4; ++k) n += t[k] * t[k];
            if (n > bestn) { bestn = n; for (int k = 0; k < 4; ++k) bestv[k] = t[k]; }
        }
        if (bestn < 1e-9) return false;
        double n = std::sqrt(bestn);
        for (int k = 0; k < 4; ++k) V[nv][k] = bestv[k] / n;
        nv++;
    }

    int bc = g_aobase[C];
    int grp[8][3]; int gsz[8], ng = 0, nm = 0;
    mo_clear(nm); g_mo[nm].c[bc + L1S] = 1.0; g_moname[nm] = "core 1s";
    grp[ng][0] = nm; gsz[ng] = 1; ng++; nm++;
    for (int i = 0; i < nb; ++i) {               // bonds
        mo_clear(nm);
        g_mo[nm].c[bc + L2S]  = V[i][0];
        g_mo[nm].c[bc + L2PX] = V[i][1];
        g_mo[nm].c[bc + L2PY] = V[i][2];
        g_mo[nm].c[bc + L2PZ] = V[i][3];
        g_mo[nm].c[g_aobase[lig[i]] + L1S] = lam;
        g_moname[nm] = "sigma bond"; mo_normalise(nm);
        grp[ng][0] = nm; gsz[ng] = 1; ng++; nm++;
    }
    for (int i = nb; i < 4; ++i) {               // lone pairs
        mo_clear(nm);
        g_mo[nm].c[bc + L2S]  = V[i][0];
        g_mo[nm].c[bc + L2PX] = V[i][1];
        g_mo[nm].c[bc + L2PY] = V[i][2];
        g_mo[nm].c[bc + L2PZ] = V[i][3];
        g_moname[nm] = "lone pair"; mo_normalise(nm);
        grp[ng][0] = nm; gsz[ng] = 1; ng++; nm++;
    }
    g_nmo = nm;
    return fill_groups(grp, gsz, ng, nel);
}

// which atom / AO dominates an MO (used for the initial walker placement and
// for labelling the orbital on screen)
static int mo_main_ao(int mo)
{
    int best = 0; double bc = -1.0;
    for (int k = 0; k < g_nao; ++k)
        if (std::fabs(g_mo[mo].c[k]) > bc) { bc = std::fabs(g_mo[mo].c[k]); best = k; }
    return best;
}

static void vmc_reset_walkers()
{
    for (int wi = 0; wi < g_nw; ++wi) {
        Walker& w = g_w[wi];
        for (int i = 0; i < g_nel; ++i) {
            // place electron i near the nucleus whose block owns its MO
            int mo = (i < g_nup) ? g_occup[i] : g_occdn[i - g_nup];
            int ao = mo_main_ao(mo);
            int at = g_ao[ao].at;
            double sc = 1.5 / std::max(1.0, (double)g_at[at].Z) * (g_ao[ao].n * g_ao[ao].n);
            for (int k = 0; k < 3; ++k) w.r[i][k] = g_at[at].R[k] + nrand() * sc;
        }
        build_slater(w);
        w.tn = 0;
        for (int i = 0; i < g_nel; ++i)
            for (int t = 0; t < NTRAIL; ++t) { w.tx[i][t] = (float)w.r[i][0]; w.ty[i][t] = (float)w.r[i][1]; w.tz[i][t] = (float)w.r[i][2]; }
        w.eloc = 0.0;
    }
    stat_clear();
}

// short VMC run used only by the variational scan below
static double vmc_trial_energy(int equil, int samp, int nw)
{
    int save = g_nw; g_nw = nw;
    vmc_reset_walkers();
    for (int s = 0; s < equil; ++s) for (int i = 0; i < nw; ++i) vmc_sweep(g_w[i]);
    double s1 = 0.0; long c = 0;
    for (int s = 0; s < samp; ++s)
        for (int i = 0; i < nw; ++i) {
            vmc_sweep(g_w[i]);
            if (g_w[i].ok) { double e = local_energy(g_w[i]); if (std::isfinite(e)) { s1 += e; c++; } }
        }
    g_nw = save;
    return c ? s1 / c : 0.0;
}

// ---------------------------------------------------------------------------
// Variational optimisation.  Coordinate descent on the trial-function
// parameters: one exponent scale per shell, an extra scale on the p orbitals,
// and the Jastrow length b.  The objective is the model's OWN energy <E_L>,
// never an experimental number, so this is the variational principle doing the
// work: whatever it lands on is an upper bound to the true ground state.
// Common random numbers are used so the comparisons are not noise-dominated.
// ---------------------------------------------------------------------------
static void vmc_optimise()
{
    static const double SC[5] = {0.90, 0.95, 1.00, 1.05, 1.11};
    static const double JB[5] = {0.3, 0.5, 0.8, 1.2, 1.8};
    int Z0 = 0;
    for (size_t a = 0; a < g_at.size(); ++a) Z0 = std::max(Z0, g_at[a].Z);
    int nsh = (Z0 <= 2) ? 1 : (Z0 <= 10 ? 2 : 3);
    uint64_t seed = g_rs;
    for (int i = 0; i < 4; ++i) g_zs[i] = 1.0;
    g_zp = 1.0; g_jb = 0.8; g_lam = 1.0; g_zh = 1.0; g_hyb = 1.0;
    setup_electrons(g_sysid);
    g_rs = seed; g_ghas = false;
    double best = vmc_trial_energy(400, 1600, 6);

    struct P { double* p; const double* v; };
    static const double LM[5] = {0.55, 0.75, 0.95, 1.15, 1.40};
    static const double HB[5] = {0.55, 0.78, 1.00, 1.30, 1.70};
    static const double ZH[5] = {0.80, 0.90, 1.00, 1.12, 1.28};
    P plist[8] = {{&g_zs[1], SC}, {&g_zs[2], SC}, {&g_zs[3], SC}, {&g_zp, SC},
                  {&g_jb, JB}, {&g_lam, LM}, {&g_zh, ZH}, {&g_hyb, HB}};
    int npar = 1 + (nsh >= 2 ? 1 : 0) + (nsh >= 3 ? 1 : 0);
    bool hasH = false, hasHeavy = false, hybrid = (int)g_at.size() > 2;
    for (size_t a = 0; a < g_at.size(); ++a) { if (g_at[a].Z == 1) hasH = true; else hasHeavy = true; }
    int order[10]; int no = 0;
    if (g_molmo) order[no++] = 5;       // bond polarity first: it matters most
    if (hybrid)  order[no++] = 7;
    order[no++] = 0;
    if (nsh >= 2) order[no++] = 1;
    if (nsh >= 3) order[no++] = 2;
    if (nsh >= 2) order[no++] = 3;      // p scale only matters if there is a p shell
    if (hasH && hasHeavy) order[no++] = 6;
    order[no++] = 4;
    if (g_molmo) order[no++] = 5;
    (void)npar;
    for (int sweep = 0; sweep < 2; ++sweep) {
        bool moved = false;
        for (int oi = 0; oi < no; ++oi) {
            P& pp = plist[order[oi]];
            double keep = *pp.p, bv = keep;
            for (int i = 0; i < 5; ++i) {
                *pp.p = pp.v[i];
                setup_electrons(g_sysid);
                g_rs = seed; g_ghas = false;
                double e = vmc_trial_energy(400, 1600, 6);
                if (std::isfinite(e) && e < best - 1e-6) { best = e; bv = pp.v[i]; moved = true; }
            }
            *pp.p = bv;
        }
        if (!moved) break;
    }
    g_zscale = g_zs[nsh];
    g_rs = seed;
    setup_electrons(g_sysid);
}

// ---------------------------------------------------------------------------
// Cached optimiser output, one row per system, in the order of the SYS table.
// Regenerate with:  atom_os opt      (mode "opt" prints exactly this block)
// This is a cache of the variational optimiser's own result, not a fit to any
// measured energy - the sim can re-run the optimiser live from the panel.
// ---------------------------------------------------------------------------
struct TrialParam { double z1, z2, z3, zp, jb, lam, zh, hyb; };
static const TrialParam TRIAL_DEFAULT = {1, 1, 1, 1, 0.8, 1, 1, 1};
static const TrialParam* trial_table(int* n);

static void trial_apply(const TrialParam& t)
{
    g_zs[1] = t.z1; g_zs[2] = t.z2; g_zs[3] = t.z3;
    g_zp = t.zp; g_jb = t.jb; g_lam = t.lam; g_zh = t.zh; g_hyb = t.hyb;
    int Z0 = 0;
    for (size_t a = 0; a < g_at.size(); ++a) Z0 = std::max(Z0, g_at[a].Z);
    g_zscale = g_zs[(Z0 <= 2) ? 1 : (Z0 <= 10 ? 2 : 3)];
    setup_electrons(g_sysid);
}

static void trial_load(int sysid)
{
    int n = 0;
    const TrialParam* T = trial_table(&n);
    trial_apply((sysid >= 0 && sysid < n) ? T[sysid] : TRIAL_DEFAULT);
}

// ===========================================================================
//        L E V E L   2  :   n u c l e o n s   ( Q M D )
// ===========================================================================
static const double MN    = 938.918;   // MeV   average nucleon mass
static const double HBARC = 197.327;   // MeV fm
static const double E2    = 1.43996;   // MeV fm   (e^2 = alpha hbar c)

// ---------------------------------------------------------------------------
// QMD parameters.  The Fermi (zero-point) kinetic energy is carried by an
// explicit Thomas-Fermi density functional rather than by the classical
// momenta, because a classical particle has no zero-point motion and a nucleus
// without it collapses.  For symmetric matter
//     T/A = (3/5) E_F = (3/5)(hbar^2/2M)(3 pi^2/2)^{2/3} rho^{2/3}
//         = C_TF rho^{2/3},   C_TF = 0.6 * 20.735 * 6.036 = 75.09 MeV fm^2
// so the model energy per nucleon in uniform matter is
//     E/A(x) = C_TF rho0^{2/3} x^{2/3} + (alpha/2) x + (beta/(gamma+1)) x^gamma
// with x = rho/rho0.  Imposing the two measured saturation properties
//     E/A(1) = -16 MeV     and     dE/dx|_1 = 0        (gamma = 2)
// gives alpha = -124.9 MeV, beta = 70.8 MeV - i.e. exactly the standard soft
// Skyrme QMD set below.  Those two numbers are therefore not free: they are
// the solution of the saturation conditions given C_TF.
//
// Consequence: the classical momenta P now describe only motion ABOVE the
// ground state, so the relaxed nucleus is static and the visible motion is a
// stated excitation energy E* - which is what gets printed on screen.
// ---------------------------------------------------------------------------
static double QL     = 2.0;      // fm^2   Gaussian width parameter
static double QRHO0  = 0.168;    // fm^-3
static double QALPHA = -124.0;   // MeV
static double QBETA  = 70.5;     // MeV
static double QGAMMA = 2.0;
static double QCTF   = 0.0;      // MeV fm^2   optional Thomas-Fermi kinetic term
static double QRBMIN = 5.0e-3;   // fm^-3      floor keeping rho^{-1/3} finite
static double QCS    = 69.984;   // MeV   symmetry (this normalisation, see below)
static double QCP    = 25.0;     // MeV   momentum-space Pauli potential
static double QQ0    = 1.80;     // fm
static double QP0    = 148.80;   // MeV/c
static double QCZ    = 41.30;    // MeV   short-range NN repulsive core
static double QAZ    = 0.65;     // fm
static double g_ndt  = 0.20;     // fm/c  nucleus time step

struct Nucleon { double R[3], P[3]; int tau, spin; float tr[NTRAIL][3]; };
static std::vector<Nucleon> g_nuc;
static int g_A = 16, g_Zn = 8;
static double g_nE = 0.0, g_nE0 = 0.0, g_nRrms = 0.0;
static int g_selnuc = 0;
static double g_omega = 0.0;     // collective rotation put in at preparation
static const int NUC_TRSTRIDE = 14;
static int    g_nuctr = 0;

struct NForce { double dR[MAXEL * 4][3]; };

// dH/dR (into fR) and dH/dP (into fP) for the whole nucleus
static void qmd_forces(const std::vector<Nucleon>& s,
                       std::vector<double>& fR, std::vector<double>& fP,
                       double* Eout)
{
    int n = (int)s.size();
    fR.assign(n * 3, 0.0); fP.assign(n * 3, 0.0);
    double pref = std::pow(4.0 * M_PI * QL, -1.5);
    std::vector<double> rb(n, 0.0);
    // pass 1: overlap densities
    std::vector<double> rho((size_t)n * n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            double d[3] = {s[i].R[0] - s[j].R[0], s[i].R[1] - s[j].R[1], s[i].R[2] - s[j].R[2]};
            double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            double v = pref * std::exp(-r2 / (4.0 * QL));
            rho[(size_t)i * n + j] = rho[(size_t)j * n + i] = v;
            rb[i] += v; rb[j] += v;
        }
    double E = 0.0;
    for (int i = 0; i < n; ++i) {
        double p2 = s[i].P[0] * s[i].P[0] + s[i].P[1] * s[i].P[1] + s[i].P[2] * s[i].P[2];
        E += p2 / (2.0 * MN);
        for (int k = 0; k < 3; ++k) fP[i * 3 + k] += s[i].P[k] / MN;   // dH/dP kinetic
        E += QBETA / ((QGAMMA + 1.0) * std::pow(QRHO0, QGAMMA)) * std::pow(rb[i], QGAMMA);
        E += QCTF * std::pow(rb[i] + QRBMIN, 2.0 / 3.0);   // Thomas-Fermi (quantum) T
    }
    double c2 = QALPHA / QRHO0;                                             // per pair
    double c3 = QBETA * QGAMMA / ((QGAMMA + 1.0) * std::pow(QRHO0, QGAMMA));
    double ct = QCTF * (2.0 / 3.0);
    double cs = QCS / QRHO0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            double d[3] = {s[i].R[0] - s[j].R[0], s[i].R[1] - s[j].R[1], s[i].R[2] - s[j].R[2]};
            double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            double r  = std::sqrt(r2); if (r < 1e-8) r = 1e-8;
            double rj = rho[(size_t)i * n + j];
            // drho/dR_i = -rho (R_i-R_j)/(2L)
            double dcoef = -rj / (2.0 * QL);
            double ti = (double)s[i].tau, tj = (double)s[j].tau;
            // energies
            E += c2 * rj;
            E += cs * ti * tj * rj;
            // gradient coefficient multiplying d[]
            double gc = c2 * dcoef
                      + c3 * (std::pow(rb[i], QGAMMA - 1.0) + std::pow(rb[j], QGAMMA - 1.0)) * dcoef
                      + ct * (std::pow(rb[i] + QRBMIN, -1.0 / 3.0) + std::pow(rb[j] + QRBMIN, -1.0 / 3.0)) * dcoef
                      + cs * ti * tj * dcoef;
            // Coulomb (protons only)
            if (s[i].tau > 0 && s[j].tau > 0) {
                double a = 2.0 * std::sqrt(QL);
                double er = std::erf(r / a);
                E += E2 * er / r;
                double dv = E2 * (2.0 * std::exp(-r2 / (a * a)) / (a * std::sqrt(M_PI) * r) - er / r2);
                gc += dv / r;
            }
            // short-range repulsive core, all pairs.  In the real NN interaction
            // this is omega exchange / quark Pauli; in this centroid model it
            // also stands in for the zero-point repulsion that a classical
            // particle does not have (without it a 4-nucleon system, which has
            // no identical (tau,sigma) pair at all, collapses to a point).
            if (QCZ != 0.0) {
                double vz = QCZ * std::exp(-r2 / (2.0 * QAZ * QAZ));
                E += vz;
                gc += -vz / (QAZ * QAZ);
            }
            // Pauli potential (identical fermions only)
            if (s[i].tau == s[j].tau && s[i].spin == s[j].spin) {
                double dp[3] = {s[i].P[0] - s[j].P[0], s[i].P[1] - s[j].P[1], s[i].P[2] - s[j].P[2]};
                double p2 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];
                double vp = QCP * std::exp(-r2 / (2.0 * QQ0 * QQ0) - p2 / (2.0 * QP0 * QP0));
                E += vp;
                gc += -vp / (QQ0 * QQ0);
                for (int k = 0; k < 3; ++k) {
                    fP[i * 3 + k] += -vp * dp[k] / (QP0 * QP0);
                    fP[j * 3 + k] -= -vp * dp[k] / (QP0 * QP0);
                }
            }
            for (int k = 0; k < 3; ++k) { fR[i * 3 + k] += gc * d[k]; fR[j * 3 + k] -= gc * d[k]; }
        }
    if (Eout) *Eout = E;
}

// explicit midpoint (RK2) step of the QMD difference equation
static void qmd_step(double dt)
{
    int n = (int)g_nuc.size(); if (!n) return;
    std::vector<double> fR, fP, gR, gP;
    double E;
    qmd_forces(g_nuc, fR, fP, &E);
    std::vector<Nucleon> mid = g_nuc;
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) {
            mid[i].R[k] = g_nuc[i].R[k] + 0.5 * dt * fP[i * 3 + k];
            mid[i].P[k] = g_nuc[i].P[k] - 0.5 * dt * fR[i * 3 + k];
        }
    qmd_forces(mid, gR, gP, nullptr);
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) {
            g_nuc[i].R[k] += dt * gP[i * 3 + k];
            g_nuc[i].P[k] -= dt * gR[i * 3 + k];
        }
    qmd_forces(g_nuc, fR, fP, &g_nE);
    // rms radius
    double cx = 0, cy = 0, cz = 0;
    for (int i = 0; i < n; ++i) { cx += g_nuc[i].R[0]; cy += g_nuc[i].R[1]; cz += g_nuc[i].R[2]; }
    cx /= n; cy /= n; cz /= n;
    double s2 = 0.0;
    for (int i = 0; i < n; ++i) {
        double dx = g_nuc[i].R[0] - cx, dy = g_nuc[i].R[1] - cy, dz = g_nuc[i].R[2] - cz;
        s2 += dx * dx + dy * dy + dz * dz;
    }
    g_nRrms = std::sqrt(s2 / n);   // rms radius of the packet centroids
    // Trail sampling stride.  qmd_step runs many times per rendered frame, so
    // recording every step gives a trail that spans a fraction of one frame and
    // is therefore invisible.  With this stride the trail covers
    // NTRAIL * NUC_TRSTRIDE * dt = 26 * 14 * 0.2 = 73 fm/c of real time.
    if (++g_nuctr >= NUC_TRSTRIDE) {
        g_nuctr = 0;
        for (int i = 0; i < n; ++i) {
            for (int t = NTRAIL - 1; t > 0; --t)
                for (int k = 0; k < 3; ++k) g_nuc[i].tr[t][k] = g_nuc[i].tr[t - 1][k];
            for (int k = 0; k < 3; ++k) g_nuc[i].tr[0][k] = (float)g_nuc[i].R[k];
        }
    }
}

static void qmd_recentre()
{
    int n = (int)g_nuc.size(); if (!n) return;
    double px = 0, py = 0, pz = 0, cx = 0, cy = 0, cz = 0;
    for (int i = 0; i < n; ++i) {
        px += g_nuc[i].P[0]; py += g_nuc[i].P[1]; pz += g_nuc[i].P[2];
        cx += g_nuc[i].R[0]; cy += g_nuc[i].R[1]; cz += g_nuc[i].R[2];
    }
    for (int i = 0; i < n; ++i) {
        g_nuc[i].P[0] -= px / n; g_nuc[i].P[1] -= py / n; g_nuc[i].P[2] -= pz / n;
        g_nuc[i].R[0] -= cx / n; g_nuc[i].R[1] -= cy / n; g_nuc[i].R[2] -= cz / n;
    }
}

// ---------------------------------------------------------------------------
// Ground state = steepest descent on the SAME H, in phase space.  The descent
// must run in R *and* P: with a Pauli potential the minimum of H has non-zero
// momenta (identical nucleons are pushed apart in momentum space too), so
// simply damping P towards zero would converge to the wrong state.
//     R <- R - muR dH/dR ,   P <- P - muP dH/dP
// A decaying random perturbation in the first phase is annealing, not physics
// input: the end point is still a true local minimum of H.
//
// Consequence worth stating plainly: a fully relaxed CLASSICAL nucleus is
// static (dR/dt = dH/dP = 0 at a minimum).  Real nuclei have ~20 MeV/nucleon
// of quantum zero-point Fermi motion, which a classical model cannot carry as
// motion.  The visible motion here is therefore a real, dialled-in excitation
// energy E* plus an optional collective rotation - both printed on screen.
// ---------------------------------------------------------------------------
static double qmd_relax(int iter)
{
    int n = (int)g_nuc.size(); if (n < 2) return 0.0;
    std::vector<double> fR, fP;
    double muR = 2.0e-4, muP = 30.0, Eprev = 1e300, E = 0.0;
    const double MAXDR = 0.06, MAXDP = 12.0;   // per-step trust region
    for (int it = 0; it < iter; ++it) {
        qmd_forces(g_nuc, fR, fP, &E);
        if (!std::isfinite(E)) return 1e30;
        if (E > Eprev) { muR *= 0.75; muP *= 0.75; } else { muR *= 1.03; muP *= 1.03; }
        if (muR > 4e-3) muR = 4e-3; if (muR < 1e-9) muR = 1e-9;
        if (muP > 4e2) muP = 4e2;   if (muP < 1e-4) muP = 1e-4;
        Eprev = E;
        for (int i = 0; i < n; ++i) {
            double dr = 0, dp = 0;
            for (int k = 0; k < 3; ++k) { dr += fR[i * 3 + k] * fR[i * 3 + k]; dp += fP[i * 3 + k] * fP[i * 3 + k]; }
            double sr = muR * std::sqrt(dr), sp = muP * std::sqrt(dp);
            double kr = (sr > MAXDR) ? MAXDR / sr : 1.0;
            double kp = (sp > MAXDP) ? MAXDP / sp : 1.0;
            for (int k = 0; k < 3; ++k) {
                g_nuc[i].R[k] -= kr * muR * fR[i * 3 + k];
                g_nuc[i].P[k] -= kp * muP * fP[i * 3 + k];
            }
        }
        double frac = (double)it / (0.45 * iter);
        if (frac < 1.0) {
            double s = 1.0 - frac;
            for (int i = 0; i < n; ++i)
                for (int k = 0; k < 3; ++k) {
                    g_nuc[i].R[k] += nrand() * 0.12 * s;
                    g_nuc[i].P[k] += nrand() * 18.0 * s;
                }
            Eprev = 1e300;
        }
        qmd_recentre();
    }
    qmd_forces(g_nuc, fR, fP, &E);
    return E;
}

// ---------------------------------------------------------------------------
// Excite the relaxed nucleus to a target E*.  Draw ONE random momentum
// direction set xi (with zero net momentum), set P = P_gs + lambda xi, and
// bisect on lambda until H reaches E_gs + E*.  H(lambda) increases
// monotonically, so this always converges.
// (Repeatedly adding fresh random kicks does not work: random kicks can only
// raise the energy, so a single overshoot runs away without limit.)
// ---------------------------------------------------------------------------
static void qmd_excite(double Estar)
{
    int n = (int)g_nuc.size(); if (n < 2 || Estar <= 0.0) return;
    std::vector<double> fR, fP, P0((size_t)n * 3), xi((size_t)n * 3);
    double E = 0.0;
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) { P0[i * 3 + k] = g_nuc[i].P[k]; xi[i * 3 + k] = nrand(); }
    for (int k = 0; k < 3; ++k) {                  // zero net momentum in the kick
        double m = 0;
        for (int i = 0; i < n; ++i) m += xi[i * 3 + k];
        m /= n;
        for (int i = 0; i < n; ++i) xi[i * 3 + k] -= m;
    }
    double target = g_nE0 + Estar;
    auto energy_at = [&](double lam) {
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < 3; ++k) g_nuc[i].P[k] = P0[i * 3 + k] + lam * xi[i * 3 + k];
        qmd_forces(g_nuc, fR, fP, &E);
        return E;
    };
    double lo = 0.0, hi = 20.0;
    int guard = 0;
    while (energy_at(hi) < target && guard++ < 40) hi *= 1.7;
    for (int it = 0; it < 70; ++it) {
        double mid = 0.5 * (lo + hi);
        if (energy_at(mid) < target) lo = mid; else hi = mid;
    }
    energy_at(0.5 * (lo + hi));
}

// Rigid collective rotation about z: P += M (omega x R).  Returns J in hbar.
static double qmd_spin(double omega)
{
    int n = (int)g_nuc.size(); if (n < 2 || omega == 0.0) return 0.0;
    for (auto& x : g_nuc) { x.P[0] += -omega * MN * x.R[1]; x.P[1] += omega * MN * x.R[0]; }
    double J = 0.0;
    for (auto& x : g_nuc) J += x.R[0] * x.P[1] - x.R[1] * x.P[0];
    return J / HBARC;
}

static int    g_nucIter = 2200;  // relaxation iterations for the ground state
static double g_nEstar = 1.2;    // MeV per nucleon of excitation energy
static double g_nJ = 0.0;
static double g_nR0gs = 0.0;
static double g_nEref = 0.0;

static void nucleus_setup(int Z, int A, bool cool)
{
    g_Zn = Z; g_A = A;
    g_nuc.clear();
    double R0 = 1.15 * std::pow((double)A, 1.0 / 3.0);
    int np = 0, nn = 0;
    for (int i = 0; i < A; ++i) {
        Nucleon x; std::memset(&x, 0, sizeof(x));
        bool proton = (np < Z) && (nn >= A - Z || (i % 2 == 0));
        if (np >= Z) proton = false;
        if (nn >= A - Z) proton = true;
        x.tau = proton ? 1 : -1; if (proton) np++; else nn++;
        // spin: pair them up so the Pauli potential sees at most 2 per state
        x.spin = ((proton ? np : nn) % 2 == 1) ? 1 : -1;
        double u = urand() * 2.0 - 1.0, ph = urand() * 2.0 * M_PI;
        double rr = R0 * std::pow(urand(), 1.0 / 3.0);
        double st = std::sqrt(std::max(0.0, 1.0 - u * u));
        x.R[0] = rr * st * std::cos(ph); x.R[1] = rr * st * std::sin(ph); x.R[2] = rr * u;
        double pf = 200.0 * std::pow(urand(), 1.0 / 3.0);
        u = urand() * 2.0 - 1.0; ph = urand() * 2.0 * M_PI; st = std::sqrt(std::max(0.0, 1.0 - u * u));
        x.P[0] = pf * st * std::cos(ph); x.P[1] = pf * st * std::sin(ph); x.P[2] = pf * u;
        for (int t = 0; t < NTRAIL; ++t) for (int k = 0; k < 3; ++k) x.tr[t][k] = (float)x.R[k];
        g_nuc.push_back(x);
    }
    std::vector<double> fR, fP;
    if (A == 1) {
        std::memset(g_nuc[0].R, 0, sizeof(double) * 3);
        std::memset(g_nuc[0].P, 0, sizeof(double) * 3);
        g_nE0 = g_nE = 0.0; g_nRrms = g_nR0gs = 0.0; g_nJ = 0.0; g_selnuc = 0;
        for (int t = 0; t < NTRAIL; ++t) for (int k = 0; k < 3; ++k) g_nuc[0].tr[t][k] = 0.0f;
        return;
    }
    if (cool) g_nE0 = qmd_relax(g_nucIter);
    else      qmd_forces(g_nuc, fR, fP, &g_nE0);
    {   // ground-state rms radius: this is the number to compare with experiment
        int n = (int)g_nuc.size(); double s2 = 0;
        for (int i = 0; i < n; ++i)
            s2 += g_nuc[i].R[0] * g_nuc[i].R[0] + g_nuc[i].R[1] * g_nuc[i].R[1] + g_nuc[i].R[2] * g_nuc[i].R[2];
        g_nR0gs = g_nRrms = std::sqrt(s2 / n);
    }
    qmd_excite(g_nEstar * A);
    g_nJ = qmd_spin(g_omega);
    qmd_forces(g_nuc, fR, fP, &g_nE);
    g_nEref = g_nE;                       // integrator reference, after excitation
    for (int i = 0; i < (int)g_nuc.size(); ++i)
        for (int t = 0; t < NTRAIL; ++t) for (int k = 0; k < 3; ++k) g_nuc[i].tr[t][k] = (float)g_nuc[i].R[k];
    g_selnuc = 0;
}

// ===========================================================================
//        L E V E L   3  :   q u a r k s   ( C o r n e l l  +  Y - s t r i n g )
// ===========================================================================
static double QK_M[2]   = {336.0, 336.0};   // MeV  constituent u, d (isospin symmetric)
static double QK_ALPHAS = 0.60;
static double QK_SIGMA  = 900.0;            // MeV/fm   string tension
static double QK_R0     = 0.08;             // fm       OGE smearing
static double QK_V0     = 0.0;              // MeV      additive constant (fitted)
static double g_qdt     = 0.0012;           // fm/c

struct Quark { double x[3], p[3]; double m; int flavour; int color; float tr[NTRAIL][3]; };
static Quark g_q[3];
static double g_qS[3];              // string junction
static double g_qE = 0.0, g_qE0 = 0.0, g_qL = 0.0, g_qRrms = 0.0;
static int    g_qIsProton = 1;
static const int QK_TRSTRIDE = 26;   // trail sampling stride (see quark_step)
static int    g_qtr = 0;

static inline double qk_kappa() { return (2.0 / 3.0) * QK_ALPHAS * HBARC; }   // MeV fm

// Steiner / Y-string geometry.  Returns the obtuse vertex index if one of the
// three angles is >= 120 deg (the junction then degenerates onto that vertex),
// otherwise -1 and S = Fermat point found by Weiszfeld iteration.
static int y_junction(const double* P[3], double S[3])
{
    for (int i = 0; i < 3; ++i) {
        const double* o1 = P[(i + 1) % 3]; const double* o2 = P[(i + 2) % 3];
        double u[3], v[3], nu = 0, nv = 0, d = 0;
        for (int k = 0; k < 3; ++k) { u[k] = o1[k] - P[i][k]; v[k] = o2[k] - P[i][k]; nu += u[k] * u[k]; nv += v[k] * v[k]; }
        nu = std::sqrt(nu); nv = std::sqrt(nv);
        if (nu < 1e-9 || nv < 1e-9) continue;
        for (int k = 0; k < 3; ++k) d += (u[k] / nu) * (v[k] / nv);
        if (d <= -0.5) { for (int k = 0; k < 3; ++k) S[k] = P[i][k]; return i; }   // >= 120 deg
    }
    for (int k = 0; k < 3; ++k) S[k] = (P[0][k] + P[1][k] + P[2][k]) / 3.0;
    for (int it = 0; it < 64; ++it) {
        double num[3] = {0, 0, 0}, den = 0.0;
        for (int i = 0; i < 3; ++i) {
            double d2 = 0; for (int k = 0; k < 3; ++k) { double t = P[i][k] - S[k]; d2 += t * t; }
            double d = std::sqrt(d2); if (d < 1e-9) d = 1e-9;
            for (int k = 0; k < 3; ++k) num[k] += P[i][k] / d;
            den += 1.0 / d;
        }
        double nx[3], mv = 0;
        for (int k = 0; k < 3; ++k) nx[k] = num[k] / den;
        for (int k = 0; k < 3; ++k) { double t = nx[k] - S[k]; mv += t * t; S[k] = nx[k]; }
        if (mv < 1e-18) break;
    }
    return -1;
}

// V(x) and -dV/dx for the three quarks
static double quark_pot(const Quark* q, double F[3][3], double* Lout)
{
    const double* P[3] = {q[0].x, q[1].x, q[2].x};
    double S[3];
    int ob = y_junction(P, S);
    for (int k = 0; k < 3; ++k) g_qS[k] = S[k];
    double L = 0.0;
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) F[i][k] = 0.0;
    if (ob >= 0) {
        // degenerate junction: two straight segments meeting at quark `ob`.
        // L = |x_ob - x_j| + |x_ob - x_k|, gradients are the segment unit vectors.
        for (int t = 1; t <= 2; ++t) {
            int j = (ob + t) % 3;
            double d[3], dn = 0;
            for (int k = 0; k < 3; ++k) { d[k] = q[ob].x[k] - q[j].x[k]; dn += d[k] * d[k]; }
            dn = std::sqrt(dn); if (dn < 1e-9) dn = 1e-9;
            L += dn;
            for (int k = 0; k < 3; ++k) {
                F[ob][k] -= QK_SIGMA * d[k] / dn;
                F[j][k]  += QK_SIGMA * d[k] / dn;
            }
        }
    } else {
        // dL/dx_i = (x_i - S)/|x_i - S| by the envelope theorem (S is optimal)
        for (int i = 0; i < 3; ++i) {
            double d[3], dn = 0;
            for (int k = 0; k < 3; ++k) { d[k] = q[i].x[k] - S[k]; dn += d[k] * d[k]; }
            dn = std::sqrt(dn); L += dn;
            if (dn > 1e-9) for (int k = 0; k < 3; ++k) F[i][k] -= QK_SIGMA * d[k] / dn;
        }
    }
    double V = QK_SIGMA * L + QK_V0;
    double kap = qk_kappa();
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) {
            double d[3], r2 = 0;
            for (int k = 0; k < 3; ++k) { d[k] = q[i].x[k] - q[j].x[k]; r2 += d[k] * d[k]; }
            double r = std::sqrt(r2); if (r < 1e-7) r = 1e-7;
            double er = std::erf(r / QK_R0);
            V += -kap * er / r;
            // d/dr[-kap erf(r/r0)/r] = -kap( 2 e^{-r^2/r0^2}/(r0 sqrt(pi) r) - erf/r^2 )
            double dv = -kap * (2.0 * std::exp(-r2 / (QK_R0 * QK_R0)) / (QK_R0 * std::sqrt(M_PI) * r) - er / r2);
            for (int k = 0; k < 3; ++k) { F[i][k] -= dv * d[k] / r; F[j][k] += dv * d[k] / r; }
        }
    if (Lout) *Lout = L;
    return V;
}

static double quark_energy(const Quark* q, double* Vout, double* Lout)
{
    double F[3][3], L;
    double V = quark_pot(q, F, &L);
    double T = 0.0;
    for (int i = 0; i < 3; ++i) {
        double p2 = q[i].p[0] * q[i].p[0] + q[i].p[1] * q[i].p[1] + q[i].p[2] * q[i].p[2];
        T += std::sqrt(p2 + q[i].m * q[i].m);
    }
    if (Vout) *Vout = V; if (Lout) *Lout = L;
    return T + V;
}

// exactly symplectic kick-drift-kick for H = sum sqrt(p^2+m^2) + V(x)
static void quark_step(double dt)
{
    double F[3][3], L;
    quark_pot(g_q, F, &L);
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) g_q[i].p[k] += 0.5 * dt * F[i][k];
    for (int i = 0; i < 3; ++i) {
        double p2 = g_q[i].p[0] * g_q[i].p[0] + g_q[i].p[1] * g_q[i].p[1] + g_q[i].p[2] * g_q[i].p[2];
        double E = std::sqrt(p2 + g_q[i].m * g_q[i].m);
        for (int k = 0; k < 3; ++k) g_q[i].x[k] += dt * g_q[i].p[k] / E;
    }
    quark_pot(g_q, F, &L);
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) g_q[i].p[k] += 0.5 * dt * F[i][k];
    double V;
    g_qE = quark_energy(g_q, &V, &g_qL);
    double cx = 0, cy = 0, cz = 0;
    for (int i = 0; i < 3; ++i) { cx += g_q[i].x[0]; cy += g_q[i].x[1]; cz += g_q[i].x[2]; }
    cx /= 3; cy /= 3; cz /= 3;
    double s2 = 0;
    for (int i = 0; i < 3; ++i) {
        double dx = g_q[i].x[0] - cx, dy = g_q[i].x[1] - cy, dz = g_q[i].x[2] - cz;
        s2 += dx * dx + dy * dy + dz * dz;
    }
    g_qRrms = std::sqrt(s2 / 3.0);
    // same stride reasoning as the nucleons: 26 * 26 * 0.0012 = 0.81 fm/c,
    // about a third of the 2.7 fm/c orbital period
    if (++g_qtr >= QK_TRSTRIDE) {
        g_qtr = 0;
        for (int i = 0; i < 3; ++i) {
            for (int t = NTRAIL - 1; t > 0; --t) for (int k = 0; k < 3; ++k) g_q[i].tr[t][k] = g_q[i].tr[t - 1][k];
            for (int k = 0; k < 3; ++k) g_q[i].tr[0][k] = (float)g_q[i].x[k];
        }
    }
}

// ---------------------------------------------------------------------------
// Ground configuration of the 3-quark system - and where its SIZE comes from.
//
// A classical relative equilibrium exists: an equilateral triangle of circum-
// radius R rotating rigidly.  Put the three quarks on a circle of radius R and
// let V_eq(R) be the model potential of that configuration.  Each quark needs
// an inward force (1/3) dV_eq/dR, and a relativistic circular orbit supplies
//      |dp/dt| = p v / R = p^2 / (E_i R),      E_i = sqrt(p^2 + m^2)
// so force balance is
//      (A)   p^2 / sqrt(p^2 + m^2) = (R/3) dV_eq/dR
//
// That single equation does NOT fix the size: it has a solution for every R,
// and the classical energy keeps falling as R -> 0, i.e. a purely classical
// nucleon collapses.  The nucleon's size is quantum, so we close the system
// with the same Bohr-Sommerfeld condition that fixes the Bohr radius of
// hydrogen - quantised angular momentum on the closed orbit:
//      (B)   p R = n hbar          (n = 1 for the ground pinwheel)
// Substituting p = n hbar c / R into (A) leaves one equation in R, solved here
// by bisection.  Nothing is dialled in: R and p come out of the model.
//
// V0 is the usual additive constant of the constituent quark model.  It is
// fitted ONCE, at n = 1, to the measured nucleon mass and then never touched,
// so the n = 2, 3 ... orbits come out genuinely heavier, like real N* states.
//
// Honest caveat: this pinwheel carries 3 n hbar of orbital angular momentum,
// while the real nucleon ground state is J = 1/2 with L = 0.  A classical
// trajectory cannot be an s-wave; this is the classical state of the same
// Hamiltonian at the same size and momentum scale.
// ---------------------------------------------------------------------------
static double g_qR0 = 0.0, g_qP0v = 0.0;
static int    g_qn  = 1;          // Bohr-Sommerfeld quantum number
static double g_qwob = 0.0;       // shape perturbation of the pinwheel

static double quark_Veq(double R)
{
    Quark t[3];
    for (int i = 0; i < 3; ++i) {
        double th = 2.0 * M_PI * i / 3.0;
        t[i] = g_q[i];
        t[i].x[0] = R * std::cos(th); t[i].x[1] = R * std::sin(th); t[i].x[2] = 0.0;
    }
    double F[3][3], L, saved = QK_V0; QK_V0 = 0.0;
    double V = quark_pot(t, F, &L);
    QK_V0 = saved;
    return V;
}

// f(R) = p^2/E_i - (R/3) dV_eq/dR   with p = n hbar c / R.  Root = the orbit.
static double quark_balance(double R, int n)
{
    double m = QK_M[0];
    double p = n * HBARC / R;
    double h = 1e-5 * std::max(0.01, R);
    double dV = (quark_Veq(R + h) - quark_Veq(R - h)) / (2.0 * h);
    return p * p / std::sqrt(p * p + m * m) - (R / 3.0) * dV;
}

static void quark_setup(bool proton, bool refit)
{
    g_qIsProton = proton ? 1 : 0;
    int fl[3];
    if (proton) { fl[0] = 0; fl[1] = 0; fl[2] = 1; } else { fl[0] = 0; fl[1] = 1; fl[2] = 1; }
    for (int i = 0; i < 3; ++i) {
        g_q[i].flavour = fl[i]; g_q[i].m = QK_M[fl[i]]; g_q[i].color = i;
        for (int k = 0; k < 3; ++k) { g_q[i].x[k] = 0; g_q[i].p[k] = 0; }
    }
    // bisect f(R) = 0 : f > 0 at small R (kinetic wins), f < 0 at large R
    double lo = 0.02, hi = 6.0;
    for (int it = 0; it < 90; ++it) {
        double mid = 0.5 * (lo + hi);
        if (quark_balance(mid, g_qn) > 0.0) lo = mid; else hi = mid;
    }
    double R = 0.5 * (lo + hi);
    double p = g_qn * HBARC / R;
    g_qR0 = R; g_qP0v = p;
    for (int i = 0; i < 3; ++i) {
        double th = 2.0 * M_PI * i / 3.0;
        g_q[i].x[0] = R * std::cos(th); g_q[i].x[1] = R * std::sin(th); g_q[i].x[2] = 0.0;
        g_q[i].p[0] = -p * std::sin(th); g_q[i].p[1] = p * std::cos(th); g_q[i].p[2] = 0.0;
    }
    if (refit) {
        QK_V0 = 0.0;
        double V, L;
        double E = quark_energy(g_q, &V, &L);
        QK_V0 = 938.918 - E;   // isospin symmetric: p/n split is EM + m_d-m_u
    }
    // a small perturbation so the motion is a genuine quasi-periodic 3-body
    // tumble instead of a perfectly rigid pinwheel (net momentum kept at zero)
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k) g_q[i].p[k] += p * g_qwob * nrand();
    double s[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) s[k] += g_q[i].p[k];
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) g_q[i].p[k] -= s[k] / 3.0;
    // fit V0 to the state that is actually going to be integrated, so the mass
    // read off the screen is the mass of THIS state
    if (refit) {
        double V, L;
        QK_V0 = 0.0;
        double E = quark_energy(g_q, &V, &L);
        QK_V0 = 938.918 - E;
    }
    double V;
    g_qE = g_qE0 = quark_energy(g_q, &V, &g_qL);
    for (int i = 0; i < 3; ++i)
        for (int t = 0; t < NTRAIL; ++t) for (int k = 0; k < 3; ++k) g_q[i].tr[t][k] = (float)g_q[i].x[k];
}

// ===========================================================================
//                        R E N D E R I N G
// ===========================================================================
static std::vector<uint32_t> g_px;
static std::vector<float>    g_dens;      // electron-cloud accumulation buffer
static std::vector<float>    g_densU;     // spin-up only, for the spin colouring
static std::vector<float>    g_dpos;      // selected orbital, phi > 0 lobes
static std::vector<float>    g_dneg;      // selected orbital, phi < 0 lobes

static inline uint32_t rgba(int r, int g, int b, float a)
{
    int A = (int)(a * 255.0f); if (A < 0) A = 0; if (A > 255) A = 255;
    if (r < 0) r = 0; if (r > 255) r = 255;
    if (g < 0) g = 0; if (g > 255) g = 255;
    if (b < 0) b = 0; if (b > 255) b = 255;
    return ((uint32_t)A << 24) | ((uint32_t)b << 16) | ((uint32_t)g << 8) | (uint32_t)r;
}

// ---------------------------------------------------------------------------
// text.  olive's built-in font has lowercase, digits and , - . only; every
// other code point is a blank that still advances.  These 5x6 glyphs fill in
// the capitals and the punctuation the readouts need.
// ---------------------------------------------------------------------------
static uint8_t g_gly[128][6];
static bool    g_glyhas[128];
static bool    g_glyinit = false;

static void font_init()
{
    if (g_glyinit) return;
    g_glyinit = true;
    std::memset(g_gly, 0, sizeof(g_gly));
    std::memset(g_glyhas, 0, sizeof(g_glyhas));
    struct E { char c; const char* h; };
    static const E T[] = {
        {'A',"0E11111F1111"},{'B',"1E111E11111E"},{'C',"0E111010110E"},
        {'D',"1E111111111E"},{'E',"1F101E10101F"},{'F',"1F101E101010"},
        {'G',"0E111017110E"},{'H',"11111F111111"},{'I',"1F040404041F"},
        {'J',"07020202120C"},{'K',"11121C121111"},{'L',"10101010101F"},
        {'M',"111B15111111"},{'N',"111915131111"},{'O',"0E111111110E"},
        {'P',"1E11111E1010"},{'Q',"0E1111150E03"},{'R',"1E11111E1211"},
        {'S',"0F100E01110E"},{'T',"1F0404040404"},{'U',"11111111110E"},
        {'V',"111111110A04"},{'W',"111111151B11"},{'X',"110A04040A11"},
        {'Y',"11110A040404"},{'Z',"1F020408101F"},
        {':',"000400000400"},{'=',"00001F001F00"},{'+',"0004041F0404"},
        {'/',"010204040810"},{'(',"020408080402"},{')',"080402020408"},
        {'%',"110204081100"},{'*',"000A040A0000"},{'<',"000408100804"},
        {'>',"000402010204"},{'?',"0E1102040004"},{'!',"040404040004"},
        {';',"000400000408"},{0x27,"040400000000"},{'"',"0A0A00000000"},
        {'_',"00000000001F"},{'^',"040A11000000"},{'|',"040404040404"},
    };
    for (const E& e : T) {
        for (int r = 0; r < 6; ++r) {
            char a = e.h[r * 2], b = e.h[r * 2 + 1];
            int hi = (a >= '0' && a <= '9') ? a - '0' : (a >= 'A' && a <= 'F' ? a - 'A' + 10 : 0);
            int lo = (b >= '0' && b <= '9') ? b - '0' : (b >= 'A' && b <= 'F' ? b - 'A' + 10 : 0);
            g_gly[(int)(unsigned char)e.c][r] = (uint8_t)(hi * 16 + lo);
        }
        g_glyhas[(int)(unsigned char)e.c] = true;
    }
}

static void txt(Olivec_Canvas oc, const char* s, int x, int y, int sz, uint32_t col)
{
    font_init();
    for (int i = 0; s[i]; ++i) {
        int c = (unsigned char)s[i];
        int gx = x + i * 6 * sz;
        if (c < 128 && g_glyhas[c]) {
            for (int r = 0; r < 6; ++r)
                for (int b = 0; b < 5; ++b)
                    if (g_gly[c][r] & (1 << (4 - b)))
                        olivec_rect(oc, gx + b * sz, y + r * sz, sz, sz, col);
        } else if (c > 32) {
            char one[2] = {(char)c, 0};
            olivec_text(oc, one, gx, y, olivec_default_font, sz, col);
        }
    }
}
static inline int txtw(const char* s, int sz) { return (int)std::strlen(s) * 6 * sz; }

// ---------------------------------------------------------------------------
// camera: yaw/pitch rotation then orthographic projection
// ---------------------------------------------------------------------------
static double g_yaw = 0.6, g_pitch = -0.30, g_rotspd = 0.0;

struct Proj { double sx, sy, dz; };
static inline Proj project(const double* p, double cx, double cy, double sc)
{
    double cy1 = std::cos(g_yaw), sy1 = std::sin(g_yaw);
    double cp = std::cos(g_pitch), sp = std::sin(g_pitch);
    double x1 =  p[0] * cy1 + p[2] * sy1;
    double z1 = -p[0] * sy1 + p[2] * cy1;
    double y2 =  p[1] * cp - z1 * sp;
    double z2 =  p[1] * sp + z1 * cp;
    Proj o; o.sx = cx + x1 * sc; o.sy = cy - y2 * sc; o.dz = z2;
    return o;
}
static inline Proj projectf(const float* p, double cx, double cy, double sc)
{
    double d[3] = {p[0], p[1], p[2]};
    return project(d, cx, cy, sc);
}

// ---------------------------------------------------------------------------
// view state
// ---------------------------------------------------------------------------
static int    g_focus   = 0;      // 0 all, 1 atom, 2 nucleus, 3 nucleon
static int    g_orbview = 0;      // 0 total cloud, else 1+MO index
static double g_zoom    = 1.0;
static double g_persist = 0.994;
static int    g_sweeps  = 2;
static bool   g_bonds   = true;
static bool   g_trails  = true;
static int    g_selatom = 0;
static double g_viewR   = 4.0;    // half width of the atom panel, in bohr
static bool   g_busy    = false;

struct Rect { int x, y, w, h; };
static Rect r_atom, r_nuc, r_qk;

static void layout()
{
    const int HT = 26;
    if (g_focus == 1)      { r_atom = {0, HT, FW, FH - HT}; r_nuc = {0,0,0,0}; r_qk = {0,0,0,0}; }
    else if (g_focus == 2) { r_nuc  = {0, HT, FW, FH - HT}; r_atom = {0,0,0,0}; r_qk = {0,0,0,0}; }
    else if (g_focus == 3) { r_qk   = {0, HT, FW, FH - HT}; r_atom = {0,0,0,0}; r_nuc = {0,0,0,0}; }
    else {
        r_atom = {0, HT, 600, FH - HT};
        r_nuc  = {600, HT, FW - 600, (FH - HT) / 2};
        r_qk   = {600, HT + (FH - HT) / 2, FW - 600, FH - HT - (FH - HT) / 2};
    }
}

// splat one point into the density accumulator
static inline void splat(std::vector<float>& buf, int x, int y, float w, const Rect& R)
{
    if (x < R.x + 2 || x >= R.x + R.w - 2 || y < R.y + 2 || y >= R.y + R.h - 2) return;
    static const float K[5][5] = {
        {0.018f, 0.082f, 0.135f, 0.082f, 0.018f},
        {0.082f, 0.368f, 0.607f, 0.368f, 0.082f},
        {0.135f, 0.607f, 1.000f, 0.607f, 0.135f},
        {0.082f, 0.368f, 0.607f, 0.368f, 0.082f},
        {0.018f, 0.082f, 0.135f, 0.082f, 0.018f}};
    for (int dy = -2; dy <= 2; ++dy) {
        float* row = &buf[(size_t)(y + dy) * FW];
        for (int dx = -2; dx <= 2; ++dx) row[x + dx] += w * K[dy + 2][dx + 2];
    }
}

static void panel_frame(Olivec_Canvas oc, const Rect& R, const char* title, uint32_t col)
{
    if (!R.w) return;
    olivec_frame(oc, R.x, R.y, R.w - 1, R.h - 1, 1, rgba(30, 60, 55, 1.f));
    txt(oc, title, R.x + 8, R.y + 6, 2, col);
}

// scale bar: pick a round length that is about a quarter of the panel
static void scale_bar(Olivec_Canvas oc, const Rect& R, double px_per_unit,
                      const char* unit, uint32_t col)
{
    if (!R.w) return;
    double want = R.w * 0.22 / px_per_unit;
    double mag = std::pow(10.0, std::floor(std::log10(want)));
    double n = want / mag;
    double v = (n >= 5.0) ? 5.0 * mag : (n >= 2.0 ? 2.0 * mag : mag);
    int len = (int)(v * px_per_unit);
    if (len < 12) len = 12;
    int bx = R.x + 12, by = R.y + R.h - 20;
    olivec_rect(oc, bx, by, len, 2, col);
    olivec_rect(oc, bx, by - 3, 2, 8, col);
    olivec_rect(oc, bx + len - 2, by - 3, 2, 8, col);
    char b[48];
    if (v >= 1.0) std::snprintf(b, sizeof(b), "%.0f %s", v, unit);
    else          std::snprintf(b, sizeof(b), "%.2f %s", v, unit);
    txt(oc, b, bx + len + 6, by - 4, 1, col);
}

// ---------------------------------------------------------------------------
// orbital probe walkers: Metropolis sampling of |phi_j(r)|^2 for ONE orbital,
// so the s / p / bonding shape can be seen on its own.  Same difference
// equation as the many-body sampler, one particle and no Jastrow:
//     r' = r + s xi ,  accept with min(1, |phi(r')|^2 / |phi(r)|^2)
// ---------------------------------------------------------------------------
static const int NPROBE = 2600;
static double g_pr[NPROBE][3];
static float  g_prs[NPROBE];      // sign of phi at that point
static int    g_prorb = -1;

static void probe_reset(int mo)
{
    g_prorb = mo;
    for (int i = 0; i < NPROBE; ++i) {
        int at = 0;
        if (mo >= 0 && mo < g_nmo) at = g_ao[mo_main_ao(mo)].at;
        for (int k = 0; k < 3; ++k) g_pr[i][k] = g_at[at].R[k] + nrand() * 1.2;
        g_prs[i] = 1.0f;
    }
}

static void probe_step(int mo)
{
    if (mo < 0 || mo >= g_nmo) return;
    if (mo != g_prorb) probe_reset(mo);
    double s = 0.16 * g_viewR;
    for (int i = 0; i < NPROBE; ++i) {
        AOVal a; mo_eval(mo, g_pr[i], a);
        double rp[3];
        for (int k = 0; k < 3; ++k) rp[k] = g_pr[i][k] + nrand() * s;
        AOVal b; mo_eval(mo, rp, b);
        double p = (a.v * a.v > 1e-300) ? (b.v * b.v) / (a.v * a.v) : 1.0;
        if (p >= 1.0 || urand() < p) {
            for (int k = 0; k < 3; ++k) g_pr[i][k] = rp[k];
            g_prs[i] = b.v >= 0.0 ? 1.0f : -1.0f;
        } else {
            g_prs[i] = a.v >= 0.0 ? 1.0f : -1.0f;
        }
    }
}

// ---------------------------------------------------------------------------
// electron configuration string, e.g. "1s2 2s2 2p4"
// ---------------------------------------------------------------------------
static const char* SS_NAME[5] = {"1s", "2s", "2p", "3s", "3p"};
static void config_str(char* out, int n)
{
    out[0] = 0;
    if (g_at.size() != 1) { std::snprintf(out, n, "%d MO occupied", (g_nup + g_ndn + 1) / 2); return; }
    int Z = g_at[0].Z, rem = Z, p = 0;
    for (int s = 0; s < 5 && rem > 0; ++s) {
        int cap = 2 * SS_COUNT[s], k = std::min(rem, cap);
        p += std::snprintf(out + p, n - p, "%s%d ", SS_NAME[s], k);
        rem -= k;
        if (p > n - 8) break;
    }
}

// ---------------------------------------------------------------------------
// panel 1 : the atom / molecule (electron cloud + nuclei + orbital)
// ---------------------------------------------------------------------------
static void draw_atom(Olivec_Canvas oc)
{
    const Rect R = r_atom; if (!R.w) return;
    double cx = R.x + R.w * 0.5, cy = R.y + R.h * 0.5;
    double sc = (std::min(R.w, R.h) * 0.42) / g_viewR * g_zoom;

    // decay the accumulated cloud, then add this frame's walkers
    float k = (float)g_persist;
    for (int y = R.y; y < R.y + R.h; ++y) {
        float* d = &g_dens[(size_t)y * FW]; float* u = &g_densU[(size_t)y * FW];
        for (int x = R.x; x < R.x + R.w; ++x) { d[x] *= k; u[x] *= k; }
    }
    for (int wi = 0; wi < g_nw; ++wi) {
        if (!g_w[wi].ok) continue;
        for (int i = 0; i < g_nel; ++i) {
            Proj p = project(g_w[wi].r[i], cx, cy, sc);
            splat(g_dens, (int)p.sx, (int)p.sy, 1.0f, R);
            if (i < g_nup) splat(g_densU, (int)p.sx, (int)p.sy, 1.0f, R);
        }
    }

    // auto exposure: the accumulated counts depend on walker number, persistence
    // and zoom, so normalise by a smoothed peak instead of a magic gain
    float peak = 0.0f;
    for (int y = R.y; y < R.y + R.h; ++y) {
        const float* d = &g_dens[(size_t)y * FW];
        for (int x = R.x; x < R.x + R.w; ++x) if (d[x] > peak) peak = d[x];
    }
    static float s_peak = 0.0f;
    s_peak = (s_peak <= 0.0f) ? peak : (0.94f * s_peak + 0.06f * peak);
    // the 1s density near a heavy nucleus is orders of magnitude above the
    // valence density, so a linear map shows either the core or the cloud but
    // never both: use a log map, auto-referenced to the running peak
    float vlow = (s_peak > 1e-9f) ? s_peak * 2.0e-3f : 1.0f;
    float lnorm = 1.0f / std::log(1.0f + s_peak / vlow);
    float cut = vlow * 0.6f;

    // colour map: cyan -> white, tinted towards magenta where spin-down wins
    for (int y = R.y; y < R.y + R.h; ++y) {
        uint32_t* px = &g_px[(size_t)y * FW];
        const float* d = &g_dens[(size_t)y * FW]; const float* u = &g_densU[(size_t)y * FW];
        for (int x = R.x; x < R.x + R.w; ++x) {
            float v = d[x];
            if (v <= cut) continue;
            float t = std::log(1.0f + v / vlow) * lnorm;
            if (t > 1.0f) t = 1.0f;
            if (g_orbview > 0) t *= 0.30f;      // let the selected orbital dominate
            float up = v > 1e-6f ? u[x] / v : 0.5f;
            float mg = 1.0f - up;                       // spin-down fraction
            int r = (int)(255.0f * t * (0.20f + 0.80f * mg * mg));
            int g = (int)(255.0f * t * (0.95f - 0.35f * mg));
            int b = (int)(255.0f * t * (0.80f + 0.20f * mg));
            uint32_t c = rgba(r, g, b, 1.f);
            olivec_blend_color(&px[x], c);
        }
    }

    // selected orbital: accumulate |phi|^2 into two buffers, one per sign of
    // phi, so the lobes come out as smooth density with the phase visible
    if (g_orbview > 0) {
        int mo = g_orbview - 1;
        probe_step(mo);
        float kk = 0.985f;
        for (int y = R.y; y < R.y + R.h; ++y) {
            float* a = &g_dpos[(size_t)y * FW]; float* b = &g_dneg[(size_t)y * FW];
            for (int x = R.x; x < R.x + R.w; ++x) { a[x] *= kk; b[x] *= kk; }
        }
        for (int i = 0; i < NPROBE; ++i) {
            Proj p = project(g_pr[i], cx, cy, sc);
            splat(g_prs[i] > 0 ? g_dpos : g_dneg, (int)p.sx, (int)p.sy, 1.0f, R);
        }
        float pk = 0.0f;
        for (int y = R.y; y < R.y + R.h; ++y) {
            const float* a = &g_dpos[(size_t)y * FW]; const float* b = &g_dneg[(size_t)y * FW];
            for (int x = R.x; x < R.x + R.w; ++x) { if (a[x] > pk) pk = a[x]; if (b[x] > pk) pk = b[x]; }
        }
        static float s_pk = 0.0f;
        s_pk = (s_pk <= 0.0f) ? pk : (0.93f * s_pk + 0.07f * pk);
        float lo2 = (s_pk > 1e-9f) ? s_pk * 4.0e-3f : 1.0f;
        float ln2 = 1.0f / std::log(1.0f + s_pk / lo2);
        for (int y = R.y; y < R.y + R.h; ++y) {
            uint32_t* px = &g_px[(size_t)y * FW];
            const float* a = &g_dpos[(size_t)y * FW]; const float* b = &g_dneg[(size_t)y * FW];
            for (int x = R.x; x < R.x + R.w; ++x) {
                float tp = (a[x] > lo2 * 0.6f) ? std::log(1.0f + a[x] / lo2) * ln2 : 0.0f;
                float tn = (b[x] > lo2 * 0.6f) ? std::log(1.0f + b[x] / lo2) * ln2 : 0.0f;
                if (tp <= 0.0f && tn <= 0.0f) continue;
                if (tp > 1.0f) tp = 1.0f; if (tn > 1.0f) tn = 1.0f;
                int r = (int)(255.0f * (0.42f * tp + 1.00f * tn));
                int g = (int)(255.0f * (0.92f * tp + 0.62f * tn));
                int bl = (int)(255.0f * (1.00f * tp + 0.16f * tn));
                olivec_blend_color(&px[x], rgba(r, g, bl, 1.f));
            }
        }
    }

    // Snapshot of where the point electrons are right now, for every walker.
    // These are samples of |Psi|^2, NOT trajectories - a Metropolis step is a
    // jump, not a path - so they are deliberately drawn as dots and never
    // joined up with lines.
    int ndot = (g_orbview > 0) ? 2 : std::min(g_nw, 10);
    for (int wi = 0; wi < ndot; ++wi) {
        const Walker& w = g_w[wi];
        if (!w.ok) continue;
        for (int i = 0; i < g_nel; ++i) {
            bool up = (i < g_nup);
            Proj p = project(w.r[i], cx, cy, sc);
            if (p.sx < R.x || p.sx >= R.x + R.w || p.sy < R.y || p.sy >= R.y + R.h) continue;
            olivec_circle(oc, (int)p.sx, (int)p.sy, 1,
                          up ? rgba(220, 255, 250, 0.42f) : rgba(255, 215, 255, 0.42f));
        }
    }

    // nuclei
    for (size_t a = 0; a < g_at.size(); ++a) {
        Proj p = project(g_at[a].R, cx, cy, sc);
        int rr = 5 + (int)(std::pow((double)g_at[a].Z, 1.0 / 3.0) * 3.0);
        bool sel = ((int)a == g_selatom);
        olivec_circle(oc, (int)p.sx, (int)p.sy, rr + 7, rgba(255, 170, 40, 0.10f));
        olivec_circle(oc, (int)p.sx, (int)p.sy, rr, rgba(255, 190, 70, 0.95f));
        if (sel) olivec_frame(oc, (int)p.sx - rr - 5, (int)p.sy - rr - 5, 2 * rr + 10, 2 * rr + 10, 1,
                              rgba(255, 255, 255, 0.8f));
        const char* sym = ELEM[g_at[a].Z].sym;
        txt(oc, sym, (int)p.sx - txtw(sym, 2) / 2, (int)p.sy - rr - 18, 2, rgba(255, 220, 140, 1.f));
    }
    // bonds (drawn only as a guide to the eye: nuclei are fixed, Born-Oppenheimer)
    for (size_t a = 0; a < g_at.size(); ++a)
        for (size_t b = a + 1; b < g_at.size(); ++b) {
            double dx = g_at[a].R[0] - g_at[b].R[0], dy = g_at[a].R[1] - g_at[b].R[1], dz = g_at[a].R[2] - g_at[b].R[2];
            double d = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (d > 3.6) continue;
            Proj p = project(g_at[a].R, cx, cy, sc), q = project(g_at[b].R, cx, cy, sc);
            olivec_line(oc, (int)p.sx, (int)p.sy, (int)q.sx, (int)q.sy, rgba(140, 120, 70, 0.35f));
        }

    // readouts (on a dimmed strip so they stay readable over the cloud)
    char b[96];
    int ty = R.y + 26, dy = 13;
    olivec_rect(oc, R.x + 1, R.y + 22, std::min(R.w - 2, 380), 6 * dy + 4, rgba(2, 6, 10, 0.75f));
    uint32_t cl = rgba(120, 230, 210, 0.95f);
    uint32_t cd = rgba(90, 160, 150, 0.9f);
    double em = stat_mean(), ev = stat_err();
    config_str(b, sizeof(b));
    txt(oc, b, R.x + 8, ty, 1, cd); ty += dy;
    std::snprintf(b, sizeof(b), "E VMC   %10.4f Ha  +/- %.4f", em, ev);
    txt(oc, b, R.x + 8, ty, 1, cl); ty += dy;
    if (SYS[g_sysid].Eref != 0.0) {
        std::snprintf(b, sizeof(b), "E exact %10.4f Ha   err %+.2f%%", SYS[g_sysid].Eref,
                      100.0 * (em - SYS[g_sysid].Eref) / std::fabs(SYS[g_sysid].Eref));
        txt(oc, b, R.x + 8, ty, 1, cd); ty += dy;
    }
    std::snprintf(b, sizeof(b), "walkers %d   accept %.0f%%   zeta x%.2f   jastrow b %.2f",
                  g_nw, 100.0 * (g_att > 0 ? g_acc / g_att : 0.0), g_zscale, g_jb);
    txt(oc, b, R.x + 8, ty, 1, cd); ty += dy;
    if (g_orbview > 0) {
        int mo = g_orbview - 1;
        int best = mo_main_ao(mo);
        int at = g_ao[best].at;
        const char* nm = "s";
        if (g_ao[best].l == 1) nm = (g_ao[best].m == 0 ? "px" : (g_ao[best].m == 1 ? "py" : "pz"));
        if (g_molmo)
            std::snprintf(b, sizeof(b), "orbital %d of %d  %s  mostly %s %d%s",
                          mo + 1, g_nmo, g_moname[mo], ELEM[g_at[at].Z].sym, g_ao[best].n, nm);
        else
            std::snprintf(b, sizeof(b), "orbital  %s %d%s", ELEM[g_at[at].Z].sym, g_ao[best].n, nm);
        txt(oc, b, R.x + 8, ty, 1, rgba(255, 190, 90, 0.95f)); ty += dy;
        txt(oc, "cyan and orange are the two signs of phi", R.x + 8, ty, 1, rgba(205, 155, 75, 0.9f));
    }
    scale_bar(oc, R, sc / 52.917721, "pm", cd);
    std::snprintf(b, sizeof(b), "1 px = %.1f pm", 52.917721 / sc);
    txt(oc, b, R.x + 12, R.y + R.h - 34, 1, cd);
}

// ---------------------------------------------------------------------------
// panel 2 : the nucleus (QMD)
// ---------------------------------------------------------------------------
static void draw_nucleus(Olivec_Canvas oc)
{
    const Rect R = r_nuc; if (!R.w) return;
    double cx = R.x + R.w * 0.5, cy = R.y + R.h * 0.52;
    double span = std::max(2.0, g_nR0gs * 3.4);
    double sc = (std::min(R.w, R.h) * 0.40) / span;
    int n = (int)g_nuc.size();

    // strong-interaction bonds: brightness is the actual overlap density that
    // enters the Hamiltonian, so what you see is the interaction being computed
    if (g_bonds) {
        double pref = std::pow(4.0 * M_PI * QL, -1.5);
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) {
                double d[3] = {g_nuc[i].R[0] - g_nuc[j].R[0], g_nuc[i].R[1] - g_nuc[j].R[1], g_nuc[i].R[2] - g_nuc[j].R[2]};
                double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                double rho = pref * std::exp(-r2 / (4.0 * QL));
                float a = (float)(rho / pref);
                if (a < 0.16f) continue;
                Proj p = project(g_nuc[i].R, cx, cy, sc), q = project(g_nuc[j].R, cx, cy, sc);
                float w = (a - 0.16f) / 0.84f;
                olivec_line(oc, (int)p.sx, (int)p.sy, (int)q.sx, (int)q.sy,
                            rgba(90, 255, 165, 0.16f + 0.72f * w));
                if (w > 0.45f)      // strong pairs get a thicker tube
                    olivec_line(oc, (int)p.sx, (int)p.sy + 1, (int)q.sx, (int)q.sy + 1,
                                rgba(140, 255, 190, 0.40f * w));
            }
    }
    // trails
    if (g_trails)
        for (int i = 0; i < n; ++i) {
            bool pr = g_nuc[i].tau > 0;
            for (int t = 1; t < NTRAIL; ++t) {
                Proj p = projectf(g_nuc[i].tr[t - 1], cx, cy, sc), q = projectf(g_nuc[i].tr[t], cx, cy, sc);
                float a = 0.55f * (1.0f - (float)t / NTRAIL);
                olivec_line(oc, (int)p.sx, (int)p.sy, (int)q.sx, (int)q.sy,
                            pr ? rgba(255, 110, 90, a) : rgba(110, 170, 255, a));
            }
        }
    // nucleons, back to front
    int order[64]; for (int i = 0; i < n && i < 64; ++i) order[i] = i;
    double dz[64];
    for (int i = 0; i < n && i < 64; ++i) dz[i] = project(g_nuc[i].R, cx, cy, sc).dz;
    for (int i = 0; i < n && i < 64; ++i)
        for (int j = i + 1; j < n && j < 64; ++j)
            if (dz[order[j]] < dz[order[i]]) std::swap(order[i], order[j]);
    int rr = std::max(3, (int)(0.40 * sc));
    for (int oi = 0; oi < n && oi < 64; ++oi) {
        int i = order[oi];
        Proj p = project(g_nuc[i].R, cx, cy, sc);
        bool pr = g_nuc[i].tau > 0;
        float sh = (float)(0.55 + 0.45 / (1.0 + std::exp(-dz[i] * 0.5)));
        olivec_circle(oc, (int)p.sx, (int)p.sy, rr + 3, pr ? rgba(255, 90, 70, 0.10f) : rgba(90, 150, 255, 0.10f));
        olivec_circle(oc, (int)p.sx, (int)p.sy, rr,
                      pr ? rgba((int)(255 * sh), (int)(95 * sh), (int)(75 * sh), 0.95f)
                         : rgba((int)(100 * sh), (int)(165 * sh), (int)(255 * sh), 0.95f));
        if (i == g_selnuc)
            olivec_frame(oc, (int)p.sx - rr - 4, (int)p.sy - rr - 4, 2 * rr + 8, 2 * rr + 8, 1, rgba(255, 255, 255, 0.9f));
    }

    char b[112];
    int ty = R.y + 26, dyy = 12;
    uint32_t cl = rgba(140, 235, 200, 0.95f), cd = rgba(95, 165, 150, 0.9f);
    olivec_rect(oc, R.x + 1, R.y + 22, std::min(R.w - 2, 390), 5 * dyy + 6, rgba(2, 6, 10, 0.72f));
    int Z = g_Zn, A = g_A;
    std::snprintf(b, sizeof(b), "%s-%d   %d p   %d n", ELEM[Z].sym, A, Z, A - Z);
    txt(oc, b, R.x + 8, ty, 1, cl); ty += dyy;
    if (A > 1) {
        std::snprintf(b, sizeof(b), "E/A %7.3f MeV  exp %7.3f  err %+.1f%%",
                      g_nE0 / A, -ELEM[Z].BA, 100.0 * (g_nE0 / A + ELEM[Z].BA) / ELEM[Z].BA);
        txt(oc, b, R.x + 8, ty, 1, cd); ty += dyy;
        double re = std::sqrt(std::max(0.0, ELEM[Z].Rch * ELEM[Z].Rch - 0.769));
        std::snprintf(b, sizeof(b), "r rms %5.3f fm   exp %5.3f   err %+.1f%%", g_nR0gs, re,
                      re > 0 ? 100.0 * (g_nR0gs - re) / re : 0.0);
        txt(oc, b, R.x + 8, ty, 1, cd); ty += dyy;
        std::snprintf(b, sizeof(b), "E* %5.2f MeV/A   J %5.1f hbar   dE %+.1e",
                      (g_nE - g_nE0) / A, g_nJ,
                      g_nEref != 0 ? (g_nE - g_nEref) / std::fabs(g_nEref) : 0.0);
        txt(oc, b, R.x + 8, ty, 1, cd); ty += dyy;
    } else {
        txt(oc, "single proton, no nuclear force to solve", R.x + 8, ty, 1, cd); ty += dyy;
    }
    txt(oc, "red p   blue n   green lines are the overlap density", R.x + 8, ty, 1, rgba(80, 140, 130, 0.85f));
    scale_bar(oc, R, sc, "fm", cd);
}

// ---------------------------------------------------------------------------
// panel 3 : one nucleon (3 quarks + Y string)
// ---------------------------------------------------------------------------
static void draw_nucleon(Olivec_Canvas oc)
{
    const Rect R = r_qk; if (!R.w) return;
    double cx = R.x + R.w * 0.5, cy = R.y + R.h * 0.52;
    double span = std::max(0.5, g_qR0 * 3.2);
    double sc = (std::min(R.w, R.h) * 0.40) / span;

    // the Y string: three flux tubes meeting at the Steiner junction
    Proj S = project(g_qS, cx, cy, sc);
    for (int i = 0; i < 3; ++i) {
        Proj p = project(g_q[i].x, cx, cy, sc);
        for (int w = 4; w >= 1; --w)
            olivec_line(oc, (int)p.sx, (int)p.sy + (w - 2), (int)S.sx, (int)S.sy + (w - 2),
                        rgba(255, 230, 120, 0.05f * w));
        olivec_line(oc, (int)p.sx, (int)p.sy, (int)S.sx, (int)S.sy, rgba(255, 240, 170, 0.55f));
    }
    olivec_circle(oc, (int)S.sx, (int)S.sy, 3, rgba(255, 255, 200, 0.8f));

    if (g_trails)
        for (int i = 0; i < 3; ++i)
            for (int t = 1; t < NTRAIL; ++t) {
                Proj p = projectf(g_q[i].tr[t - 1], cx, cy, sc), q = projectf(g_q[i].tr[t], cx, cy, sc);
                olivec_line(oc, (int)p.sx, (int)p.sy, (int)q.sx, (int)q.sy, rgba(210, 255, 255, 0.45f * (1.0f - (float)t / NTRAIL)));
            }

    static const int CC[3][3] = {{255, 80, 80}, {90, 255, 120}, {110, 150, 255}};   // colour charge
    for (int i = 0; i < 3; ++i) {
        Proj p = project(g_q[i].x, cx, cy, sc);
        int c = g_q[i].color % 3;
        bool isu = (g_q[i].flavour == 0);
        // outer halo = colour charge, inner disc = flavour, so switching between
        // a proton (u u d) and a neutron (u d d) is visible at a glance
        olivec_circle(oc, (int)p.sx, (int)p.sy, 16, rgba(CC[c][0], CC[c][1], CC[c][2], 0.10f));
        olivec_circle(oc, (int)p.sx, (int)p.sy, 12, rgba(CC[c][0], CC[c][1], CC[c][2], 0.78f));
        olivec_circle(oc, (int)p.sx, (int)p.sy, 8,
                      isu ? rgba(255, 250, 235, 0.96f) : rgba(38, 38, 56, 0.96f));
        txt(oc, isu ? "u" : "d", (int)p.sx - 3, (int)p.sy - 3, 1,
            isu ? rgba(20, 20, 30, 1.f) : rgba(240, 240, 255, 1.f));
        const char* q = isu ? "+2/3" : "-1/3";
        txt(oc, q, (int)p.sx - txtw(q, 1) / 2, (int)p.sy + 19, 1,
            isu ? rgba(255, 240, 200, 0.9f) : rgba(180, 200, 255, 0.9f));
    }

    char b[112];
    int ty = R.y + 26, dyy = 12;
    uint32_t cl = rgba(255, 220, 140, 0.95f), cd = rgba(170, 150, 110, 0.9f);
    olivec_rect(oc, R.x + 1, R.y + 22, std::min(R.w - 2, 390), 5 * dyy + 22, rgba(2, 6, 10, 0.72f));
    std::snprintf(b, sizeof(b), "%s  %s", g_qIsProton ? "PROTON" : "NEUTRON", g_qIsProton ? "u u d" : "u d d");
    txt(oc, b, R.x + 8, ty, 2, cl); ty += dyy + 5;
    std::snprintf(b, sizeof(b), "mass %8.2f MeV   exp %8.2f", g_qE, g_qIsProton ? 938.272 : 939.565);
    txt(oc, b, R.x + 8, ty, 1, cd); ty += dyy;
    double v = g_qP0v / std::sqrt(g_qP0v * g_qP0v + QK_M[0] * QK_M[0]);
    std::snprintf(b, sizeof(b), "n %d   r0 %.3f fm   p %.0f MeV/c   v %.3f c", g_qn, g_qR0, g_qP0v, v);
    txt(oc, b, R.x + 8, ty, 1, cd); ty += dyy;
    std::snprintf(b, sizeof(b), "string L %.3f fm   dE %+.2e", g_qL, g_qE0 != 0 ? (g_qE - g_qE0) / std::fabs(g_qE0) : 0.0);
    txt(oc, b, R.x + 8, ty, 1, cd); ty += dyy;
    txt(oc, "rgb rings are colour charge, yellow is the flux tube", R.x + 8, ty, 1, rgba(150, 130, 95, 0.85f));
    scale_bar(oc, R, sc, "fm", cd);
}

// ===========================================================================
//                        A B I
// ===========================================================================
// generated by `atom_os opt`; see trial_table() above
//        z1      z2      z3      zp      jb      lam     zh      hyb
static const TrialParam TRIAL[] = {
    {1.0000, 1.0000, 1.0000, 1.0000, 0.8000, 1.0000, 1.0000, 1.0000},   // H
    {0.9500, 1.0000, 1.0000, 1.0000, 0.8000, 1.0000, 1.0000, 1.0000},   // He
    {1.0000, 1.0000, 1.0000, 1.0000, 0.8000, 1.0000, 1.0000, 1.0000},   // Li
    {1.0000, 1.0000, 1.0000, 1.0000, 0.8000, 1.0000, 1.0000, 1.0000},   // Be
    {1.0000, 1.0000, 1.0000, 1.0000, 0.8000, 1.0000, 1.0000, 1.0000},   // B
    {0.9500, 1.1100, 1.0000, 1.0000, 0.5000, 1.0000, 1.0000, 1.0000},   // C
    {0.9500, 1.0500, 1.0000, 1.1100, 0.8000, 1.0000, 1.0000, 1.0000},   // N
    {1.0000, 1.0500, 1.0000, 1.0000, 1.2000, 1.0000, 1.0000, 1.0000},   // O
    {0.9500, 1.0000, 1.0000, 1.0000, 0.8000, 1.0000, 1.0000, 1.0000},   // F
    {1.1100, 1.0500, 1.0000, 1.0000, 0.8000, 1.0000, 1.0000, 1.0000},   // Ne
    {1.1100, 1.0500, 0.9000, 1.0000, 1.8000, 1.0000, 1.0000, 1.0000},   // Na
    {1.1100, 1.1100, 1.0000, 1.1100, 0.8000, 1.0000, 1.0000, 1.0000},   // Mg
    {0.9500, 1.1100, 1.0000, 1.0500, 0.3000, 1.0000, 1.0000, 1.0000},   // Al
    {1.1100, 1.0500, 0.9000, 1.0000, 0.8000, 1.0000, 1.0000, 1.0000},   // Si
    {1.0000, 1.0000, 1.0000, 1.0000, 0.8000, 1.0000, 1.0000, 1.0000},   // P
    {1.1100, 1.1100, 0.9000, 1.0000, 1.8000, 1.0000, 1.0000, 1.0000},   // S
    {0.9000, 1.1100, 0.9000, 1.1100, 0.8000, 1.0000, 1.0000, 1.0000},   // Cl
    {1.1100, 1.1100, 1.0500, 0.9500, 0.8000, 1.0000, 1.0000, 1.0000},   // Ar
    {1.0500, 1.0000, 1.0000, 1.0000, 0.8000, 1.0000, 1.0000, 1.0000},   // H2
    {1.1100, 1.0500, 1.0000, 1.1100, 1.2000, 0.5500, 1.1200, 1.0000},   // H2O
    {0.9500, 1.0500, 1.0000, 1.0500, 1.8000, 1.1500, 1.0000, 1.0000},   // O2
    {1.0000, 1.1100, 1.0000, 1.0500, 1.8000, 0.7500, 1.0000, 1.0000},   // N2
    {0.9500, 0.9500, 1.0000, 1.0000, 1.8000, 1.1500, 1.0000, 1.3000},   // CH4
    {1.0000, 1.0000, 1.0000, 1.0000, 0.8000, 0.5500, 1.0000, 1.7000},   // NH3
    {1.1100, 1.0000, 1.0000, 1.0000, 0.8000, 1.0000, 1.0000, 1.0000},   // HF
    {0.9500, 1.1100, 1.0000, 1.0500, 1.8000, 1.1500, 1.0000, 1.0000},   // CO
};
static const TrialParam* trial_table(int* n)
{
    *n = (int)(sizeof(TRIAL) / sizeof(TRIAL[0]));
    return TRIAL;
}

static void full_setup()
{
    g_busy = true;
    setup_electrons(g_sysid);
    trial_load(g_sysid);
    vmc_reset_walkers();
    for (int s = 0; s < 400; ++s) for (int i = 0; i < g_nw; ++i) vmc_sweep(g_w[i]);
    stat_clear();
    // the nucleus shown is the one of the selected (by default heaviest) atom
    g_selatom = 0;
    for (size_t a = 1; a < g_at.size(); ++a) if (g_at[a].Z > g_at[g_selatom].Z) g_selatom = (int)a;
    int Z = g_at[g_selatom].Z;
    nucleus_setup(Z, ELEM[Z].A, true);
    quark_setup(true, true);
    g_viewR = 0.0;
    for (size_t a = 0; a < g_at.size(); ++a) {
        double rn = std::sqrt(g_at[a].R[0] * g_at[a].R[0] + g_at[a].R[1] * g_at[a].R[1] + g_at[a].R[2] * g_at[a].R[2]);
        int Za = g_at[a].Z;
        int nmax = (Za <= 2) ? 1 : (Za <= 10 ? 2 : 3);
        double rout = 3.0 * nmax * nmax / std::max(0.8, g_zt[a][nmax] * g_zscale);
        g_viewR = std::max(g_viewR, rn + rout);
    }
    if (g_viewR < 1.2) g_viewR = 1.2;
    std::fill(g_dens.begin(), g_dens.end(), 0.0f);
    std::fill(g_densU.begin(), g_densU.end(), 0.0f);
    std::fill(g_dpos.begin(), g_dpos.end(), 0.0f);
    std::fill(g_dneg.begin(), g_dneg.end(), 0.0f);
    g_prorb = -1;
    g_busy = false;
}

extern "C" {
KEEP int sim_w() { return FW; }
KEEP int sim_h() { return FH; }

KEEP void sim_reset()
{
    g_rs = 88172645463325252ull; g_ghas = false;
    full_setup();
}

KEEP int sim_init(int, int)
{
    g_px.assign((size_t)FW * FH, 0);
    g_dens.assign((size_t)FW * FH, 0.0f);
    g_densU.assign((size_t)FW * FH, 0.0f);
    g_dpos.assign((size_t)FW * FH, 0.0f);
    g_dneg.assign((size_t)FW * FH, 0.0f);
    build_systems();
    font_init();
    layout();
    sim_reset();
    return 1;
}

KEEP void sim_set(int id, double v)
{
    switch (id) {
    case 0: { int s = (int)(v + 0.5);
              if (s != g_sysid && s >= 0 && s < (int)SYS.size()) { g_sysid = s; full_setup(); } } break;
    case 1: { int o = (int)(v + 0.5); if (o < 0) o = 0; if (o > g_nmo) o = g_nmo;
              if (o != g_orbview) { g_orbview = o;
                  std::fill(g_dpos.begin(), g_dpos.end(), 0.0f);
                  std::fill(g_dneg.begin(), g_dneg.end(), 0.0f); } } break;
    case 2: g_zoom = v; break;
    case 3: g_rotspd = v; break;
    case 4: if (std::fabs(v - g_nEstar) > 1e-9) { g_nEstar = v; int Z = g_at[g_selatom].Z; nucleus_setup(Z, ELEM[Z].A, true); } break;
    case 5: if (std::fabs(v - g_omega) > 1e-12) { g_omega = v; int Z = g_at[g_selatom].Z; nucleus_setup(Z, ELEM[Z].A, true); } break;
    case 6: { int n = (int)(v + 0.5); if (n < 1) n = 1; if (n > 5) n = 5;
              if (n != g_qn) { g_qn = n; quark_setup(g_qIsProton != 0, n == 1); } } break;
    case 7: g_persist = v; break;
    case 8: g_sweeps = std::max(1, std::min(6, (int)(v + 0.5))); break;
    case 9: if (std::fabs(v - g_qwob) > 1e-9) { g_qwob = v; quark_setup(g_qIsProton != 0, true); } break;
    default: break;
    }
}

KEEP void sim_action(int id)
{
    switch (id) {
    case 0: g_focus = (g_focus + 1) % 4; layout();
            std::fill(g_dens.begin(), g_dens.end(), 0.0f);
            std::fill(g_densU.begin(), g_densU.end(), 0.0f);
            std::fill(g_dpos.begin(), g_dpos.end(), 0.0f);
            std::fill(g_dneg.begin(), g_dneg.end(), 0.0f); break;
    case 1: g_bonds = !g_bonds; break;
    case 2: g_trails = !g_trails; break;
    case 3: quark_setup(g_qIsProton == 0, false); break;
    case 4: { int Z = g_at[g_selatom].Z; nucleus_setup(Z, ELEM[Z].A, true); } break;
    case 6: g_busy = true; vmc_optimise(); vmc_reset_walkers();
            for (int s = 0; s < 400; ++s) for (int i = 0; i < g_nw; ++i) vmc_sweep(g_w[i]);
            stat_clear(); g_busy = false; break;
    case 5: g_yaw = 0.6; g_pitch = -0.30;
            std::fill(g_dens.begin(), g_dens.end(), 0.0f);
            std::fill(g_densU.begin(), g_densU.end(), 0.0f);
            std::fill(g_dpos.begin(), g_dpos.end(), 0.0f);
            std::fill(g_dneg.begin(), g_dneg.end(), 0.0f); break;
    default: break;
    }
}

static double g_lastnx = -1, g_lastny = -1;

KEEP void sim_click(double nx, double ny)
{
    if (nx < 0.0 || ny < 0.0) { g_lastnx = g_lastny = -1.0; return; }   // pointer up
    double px = nx * FW, py = ny * FH;
    // drag inside the atom panel rotates the camera
    if (r_atom.w && px >= r_atom.x && px < r_atom.x + r_atom.w && py >= r_atom.y) {
        if (g_lastnx >= 0) {
            g_yaw   += (px - g_lastnx) * 0.008;
            g_pitch += (py - g_lastny) * 0.006;
            if (g_pitch > 1.4) g_pitch = 1.4; if (g_pitch < -1.4) g_pitch = -1.4;
            std::fill(g_dens.begin(), g_dens.end(), 0.0f);
            std::fill(g_densU.begin(), g_densU.end(), 0.0f);
            std::fill(g_dpos.begin(), g_dpos.end(), 0.0f);
            std::fill(g_dneg.begin(), g_dneg.end(), 0.0f);
        }
        g_lastnx = px; g_lastny = py;
        // clicking near a nucleus selects which atom's nucleus is shown
        double cx = r_atom.x + r_atom.w * 0.5, cy = r_atom.y + r_atom.h * 0.5;
        double sc = (std::min(r_atom.w, r_atom.h) * 0.42) / g_viewR * g_zoom;
        for (size_t a = 0; a < g_at.size(); ++a) {
            Proj p = project(g_at[a].R, cx, cy, sc);
            if ((p.sx - px) * (p.sx - px) + (p.sy - py) * (p.sy - py) < 400.0 && (int)a != g_selatom) {
                g_selatom = (int)a;
                int Z = g_at[a].Z; nucleus_setup(Z, ELEM[Z].A, true);
                break;
            }
        }
        return;
    }
    // clicking a nucleon selects it and shows its quarks
    if (r_nuc.w && px >= r_nuc.x && px < r_nuc.x + r_nuc.w && py >= r_nuc.y && py < r_nuc.y + r_nuc.h) {
        double cx = r_nuc.x + r_nuc.w * 0.5, cy = r_nuc.y + r_nuc.h * 0.52;
        double span = std::max(2.0, g_nR0gs * 3.4);
        double sc = (std::min(r_nuc.w, r_nuc.h) * 0.40) / span;
        int bi = -1; double bd = 1e30;
        for (int i = 0; i < (int)g_nuc.size(); ++i) {
            Proj p = project(g_nuc[i].R, cx, cy, sc);
            double d = (p.sx - px) * (p.sx - px) + (p.sy - py) * (p.sy - py);
            if (d < bd) { bd = d; bi = i; }
        }
        if (bi >= 0 && bd < 900.0) {
            g_selnuc = bi;
            bool pr = g_nuc[bi].tau > 0;
            if ((g_qIsProton != 0) != pr) quark_setup(pr, false);
        }
    }
}

KEEP void sim_step(int n)
{
    if (g_busy) return;
    for (int s = 0; s < n; ++s) {
        for (int k = 0; k < g_sweeps; ++k)
            for (int i = 0; i < g_nw; ++i) {
                vmc_sweep(g_w[i]);
                if (g_w[i].ok) {
                    double e = local_energy(g_w[i]);
                    if (std::isfinite(e)) { g_w[i].eloc = e; stat_add(e); }
                }
            }
        for (int k = 0; k < 8; ++k) qmd_step(g_ndt);
        for (int k = 0; k < 34; ++k) quark_step(g_qdt);
        g_yaw += g_rotspd * 0.012;
        if (g_rotspd > 1e-6) {
            for (size_t i = 0; i < g_dens.size(); ++i) { g_dens[i] *= 0.90f; g_densU[i] *= 0.90f; }
        }
    }
}

KEEP uint8_t* sim_render()
{
    Olivec_Canvas oc = olivec_canvas(g_px.data(), FW, FH, FW);
    olivec_fill(oc, rgba(5, 6, 12, 1.f));
    layout();
    draw_atom(oc);
    draw_nucleus(oc);
    draw_nucleon(oc);
    // header
    olivec_rect(oc, 0, 0, FW, 24, rgba(9, 14, 18, 1.f));
    olivec_rect(oc, 0, 24, FW, 1, rgba(30, 70, 62, 1.f));
    char b[160];
    std::snprintf(b, sizeof(b), "%s   %d electrons   %d nuclei", SYS[g_sysid].name, g_nel, (int)g_at.size());
    txt(oc, b, 8, 8, 2, rgba(150, 255, 225, 1.f));
    const char* fn[4] = {"all three scales", "atom", "nucleus", "nucleon"};
    std::snprintf(b, sizeof(b), "view %s", fn[g_focus]);
    txt(oc, b, FW - txtw(b, 1) - 10, 10, 1, rgba(110, 190, 175, 0.9f));
    // panel titles and the honest zoom ladder
    double satom = (r_atom.w ? (std::min(r_atom.w, r_atom.h) * 0.42) / g_viewR * g_zoom : 1.0);
    double snuc  = (r_nuc.w  ? (std::min(r_nuc.w, r_nuc.h) * 0.40) / std::max(2.0, g_nR0gs * 3.4) : 1.0);
    double sqk   = (r_qk.w   ? (std::min(r_qk.w, r_qk.h) * 0.40) / std::max(0.5, g_qR0 * 3.2) : 1.0);
    panel_frame(oc, r_atom, "ELECTRONS  VMC", rgba(120, 255, 220, 1.f));
    std::snprintf(b, sizeof(b), "NUCLEUS  QMD   zoom x%.1e", (snuc / 1e-15) / (satom / 0.529177e-10));
    panel_frame(oc, r_nuc, b, rgba(140, 255, 200, 1.f));
    std::snprintf(b, sizeof(b), "NUCLEON  3 QUARKS   zoom x%.1e", (sqk / 1e-15) / (satom / 0.529177e-10));
    panel_frame(oc, r_qk, b, rgba(255, 220, 140, 1.f));
    return (uint8_t*)g_px.data();
}
}  // extern "C"

#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
// ---------------------------------------------------------------------------
// native self test.  Modes:  e  (electrons)   n (nuclei)   q (quarks)
//                            png (render a frame)
//                            cal (QMD parameter calibration)   all (default)
// ---------------------------------------------------------------------------
// experimental point-nucleon rms radius = sqrt(Rch^2 - <r_p^2>), <r_p^2>=0.769 fm^2
static double exp_point_rms(int Z)
{
    double rc = ELEM[Z].Rch;
    double v = rc * rc - 0.769;
    return v > 0 ? std::sqrt(v) : 0.0;
}

// run the built-in variational optimiser for every system and print it as a C
// table.  The table is only a cache: it is the optimiser's own output, so
// shipping it changes nothing except how long a reset takes.  The live sim can
// re-run the optimiser at any time (action 6).
static void emit_param_table()
{
    build_systems();
    std::printf("static const TrialParam TRIAL[] = {\n");
    for (int id = 0; id < (int)SYS.size(); ++id) {
        g_rs = 88172645463325252ull; g_ghas = false;
        g_sysid = id;
        setup_electrons(id);
        vmc_optimise();
        std::printf("    {%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f},   // %s\n",
                    g_zs[1], g_zs[2], g_zs[3], g_zp, g_jb, g_lam, g_zh, g_hyb, SYS[id].name);
        std::fflush(stdout);
    }
    std::printf("};\n");
}

static void test_electrons(int NS)
{
    std::printf("== level 1: VMC energies (Ha) ==\n");
    std::printf("%-6s %4s %11s %11s %7s %6s\n", "sys", "Nel", "E_vmc", "E_exact", "err%", "acc%");
    int list[] = {0, 1, 7, 9, 17, 18, 19, 20, 21, 22, 23, 24, 25};
    for (int li = 0; li < (int)(sizeof(list) / sizeof(int)); ++li) {
        int id = list[li];
        if (id >= (int)SYS.size()) continue;
        g_rs = 88172645463325252ull; g_ghas = false;
        g_sysid = id;
        setup_electrons(id);
        trial_load(id);
        g_nw = 16; vmc_reset_walkers();
        for (int s = 0; s < 3000; ++s) for (int i = 0; i < g_nw; ++i) vmc_sweep(g_w[i]);
        stat_clear();
        for (int s = 0; s < NS; ++s)
            for (int i = 0; i < g_nw; ++i) {
                vmc_sweep(g_w[i]);
                double e = local_energy(g_w[i]);
                if (std::isfinite(e)) stat_add(e);
            }
        double m = stat_mean();
        double sd = stat_err();
        double er = SYS[id].Eref != 0.0 ? 100.0 * (m - SYS[id].Eref) / std::fabs(SYS[id].Eref) : 0.0;
        std::printf("%-6s %4d %11.4f %11.4f %7.2f %6.1f   se %.4f\n",
                    SYS[id].name, g_nel, m, SYS[id].Eref, er, 100.0 * g_acc / g_att, sd);
    }
}

// integrator stability: relax, excite, then evolve and watch H
static void test_nuc_dynamics()
{
    std::printf("\n== level 2: dynamics stability ==\n");
    for (double dt : {0.05, 0.10, 0.20, 0.30, 0.50}) {
        g_ndt = dt;
        for (int Z : {8, 18}) {
            int A = ELEM[Z].A;
            g_rs = 4242ull; g_ghas = false;
            g_nEstar = 1.2; g_omega = 0.0; g_nucIter = 2200;
            nucleus_setup(Z, A, true);
            double E0 = g_nE, emin = E0, emax = E0, rmax = 0;
            bool bad = false;
            for (int s = 0; s < 4000; ++s) {
                qmd_step(dt);
                if (!std::isfinite(g_nE)) { bad = true; break; }
                if (g_nE < emin) emin = g_nE;
                if (g_nE > emax) emax = g_nE;
                if (g_nRrms > rmax) rmax = g_nRrms;
            }
            std::printf("dt %4.2f  %-3s A %2d  E0/A %8.3f  drift %+.3e  band %.3f MeV/A  rmax %6.2f fm %s\n",
                        dt, ELEM[Z].sym, A, E0 / A, (g_nE - E0) / std::fabs(E0), (emax - emin) / A, rmax,
                        bad ? "DIVERGED" : "");
        }
    }
    g_ndt = 0.20;
}

static void test_nuclei()
{
    std::printf("\n== level 2: QMD ground states ==\n");
    std::printf("%-4s %3s %3s %9s %9s %8s %8s\n", "nuc", "Z", "A", "E/A", "exp", "Rrms", "expRpt");
    int zz[] = {2, 4, 6, 8, 10, 14, 18};
    for (int i = 0; i < 7; ++i) {
        int Z = zz[i], A = ELEM[Z].A;
        g_rs = 12345678901ull; g_ghas = false;
        double sE = g_nEstar, sO = g_omega; g_nEstar = 0.0; g_omega = 0.0;
        nucleus_setup(Z, A, true);
        g_nEstar = sE; g_omega = sO;
        std::printf("%-4s %3d %3d %9.3f %9.3f %8.3f %8.3f\n",
                    ELEM[Z].sym, Z, A, g_nE0 / A, -ELEM[Z].BA, g_nR0gs, exp_point_rms(Z));
    }
}

// ground state of one nuclide, averaged over `nseed` independent anneals so a
// lucky (or unlucky) random start cannot decide the calibration
static bool qmd_ground(int Z, int A, int nseed, int iter, double* eaOut, double* rrOut)
{
    double ea = 0.0, rr = 0.0;
    double sE = g_nEstar, sO = g_omega; g_nEstar = 0.0; g_omega = 0.0;
    for (int s = 0; s < nseed; ++s) {
        g_rs = 999331ull + (uint64_t)s * 7771ull; g_ghas = false;
        g_nucIter = iter;
        nucleus_setup(Z, A, true);
        if (!std::isfinite(g_nE0) || !std::isfinite(g_nR0gs) || g_nR0gs > 14.0 || g_nR0gs < 0.25) {
            g_nEstar = sE; g_omega = sO; return false;
        }
        ea += g_nE0 / A; rr += g_nR0gs;
    }
    g_nEstar = sE; g_omega = sO;
    *eaOut = ea / nseed; *rrOut = rr / nseed;
    return true;
}

static const int CAL_Z[9] = {2, 4, 6, 8, 10, 12, 14, 16, 18};

static double qmd_score(bool verbose, int nseed, int iter)
{
    double score = 0.0;
    for (int t = 0; t < 9; ++t) {
        int Z = CAL_Z[t], A = ELEM[Z].A;
        double ea, rr;
        if (!qmd_ground(Z, A, nseed, iter, &ea, &rr)) return 1e9;
        double re = exp_point_rms(Z);
        double de = (ea + ELEM[Z].BA) / ELEM[Z].BA, dr = (rr - re) / re;
        score += de * de + dr * dr;
        if (verbose)
            std::printf("   %-3s A %2d  E/A %8.3f (%7.3f)  err %6.2f%%   rrms %6.3f (%6.3f)  err %6.2f%%\n",
                        ELEM[Z].sym, A, ea, -ELEM[Z].BA,
                        100.0 * (ea + ELEM[Z].BA) / ELEM[Z].BA, rr, re, 100.0 * (rr - re) / re);
    }
    return std::sqrt(score / 18.0);   // rms relative error over both observables
}

// ---------------------------------------------------------------------------
// Coordinate descent over the model parameters that are NOT fixed by nuclear
// matter saturation.  alpha, beta, gamma stay at the standard soft-Skyrme
// values throughout.  Everything is fitted to two measured observables only:
// binding energy per nucleon and point-nucleon rms radius, for nine nuclides
// from He-4 to Ar-40.
// ---------------------------------------------------------------------------
static void calibrate_qmd()
{
    struct Knob { const char* name; double* p; int n; double v[9]; };
    Knob knobs[] = {
        {"QL",  &QL,  8, {1.0, 1.4, 1.8, 2.2, 2.7, 3.3, 4.0, 5.0}},
        {"QCP", &QCP, 8, {0, 8, 15, 25, 40, 60, 90, 140}},
        {"QQ0", &QQ0, 7, {0.9, 1.2, 1.5, 1.8, 2.2, 2.6, 3.1}},
        {"QP0", &QP0, 7, {80, 120, 160, 210, 270, 340, 430}},
        {"QCZ", &QCZ, 9, {0, 25, 60, 120, 240, 450, 800, 1400, 2400}},
        {"QAZ", &QAZ, 8, {0.3, 0.4, 0.5, 0.65, 0.8, 1.0, 1.25, 1.55}},
        {"QCS", &QCS, 5, {0, 12, 25, 40, 60}},
    };
    const int NK = (int)(sizeof(knobs) / sizeof(knobs[0]));
    const int NSEED = 2, ITER = 2200;
    double best = qmd_score(false, NSEED, ITER);
    std::printf("start score %.5f\n", best);
    for (int sweep = 0; sweep < 7; ++sweep) {
        bool moved = false;
        for (int k = 0; k < NK; ++k) {
            double keep = *knobs[k].p, bv = keep;
            for (int i = 0; i < knobs[k].n; ++i) {
                *knobs[k].p = knobs[k].v[i];
                double sc = qmd_score(false, NSEED, ITER);
                if (sc < best - 1e-9) { best = sc; bv = knobs[k].v[i]; moved = true; }
            }
            *knobs[k].p = bv;
            std::printf("  sweep %d  %-4s -> %8.3f   score %.5f\n", sweep, knobs[k].name, bv, best);
        }
        if (!moved) break;
    }
    std::printf("coarse  QL %.2f QCP %.1f QQ0 %.2f QP0 %.1f QCZ %.1f QAZ %.2f QCS %.1f  score %.5f\n",
                QL, QCP, QQ0, QP0, QCZ, QAZ, QCS, best);

    // ---- refinement: multiplicative steps around the coarse winner, with a
    // longer relaxation and more seeds so the score itself is more precise ----
    const int NSEED2 = 3, ITER2 = 2800;
    static const double MUL[6] = {0.72, 0.85, 0.93, 1.08, 1.18, 1.40};
    best = qmd_score(false, NSEED2, ITER2);
    std::printf("refine start score %.5f\n", best);
    for (int sweep = 0; sweep < 5; ++sweep) {
        bool moved = false;
        for (int k = 0; k < NK; ++k) {
            double keep = *knobs[k].p, bv = keep;
            if (keep == 0.0) continue;
            for (int i = 0; i < 6; ++i) {
                *knobs[k].p = keep * MUL[i];
                double sc = qmd_score(false, NSEED2, ITER2);
                if (sc < best - 1e-9) { best = sc; bv = keep * MUL[i]; moved = true; }
            }
            *knobs[k].p = bv;
            std::printf("  refine %d  %-4s -> %8.3f   score %.5f\n", sweep, knobs[k].name, bv, best);
        }
        if (!moved) break;
    }
    std::printf("\nCALIBRATED\n");
    std::printf("static double QL     = %.4f;\nstatic double QCP    = %.4f;\n"
                "static double QQ0    = %.4f;\nstatic double QP0    = %.4f;\n"
                "static double QCZ    = %.4f;\nstatic double QAZ    = %.4f;\n"
                "static double QCS    = %.4f;\n", QL, QCP, QQ0, QP0, QCZ, QAZ, QCS);
    std::printf("rms relative error %.5f\n", best);
    qmd_score(true, 4, 3200);
}

static void test_quarks()
{
    std::printf("\n== level 3: quarks ==\n");
    g_rs = 5551212ull; g_ghas = false;
    quark_setup(true, true);
    double V, L, E0 = quark_energy(g_q, &V, &L);
    double rmin = 1e9, rmax = 0, lmin = 1e9, lmax = 0, emin = 1e9, emax = -1e9;
    for (int s = 0; s < 400000; ++s) {
        quark_step(g_qdt);
        if (g_qRrms < rmin) rmin = g_qRrms;
        if (g_qRrms > rmax) rmax = g_qRrms;
        if (g_qL < lmin) lmin = g_qL;
        if (g_qL > lmax) lmax = g_qL;
        if (g_qE < emin) emin = g_qE;
        if (g_qE > emax) emax = g_qE;
    }
    double E1 = quark_energy(g_q, &V, &L);
    std::printf("nucleon mass  %.3f MeV        V0 %.1f MeV\n", E0, QK_V0);
    std::printf("orbit r0 %.4f fm   p0 %.1f MeV per c   v per c %.3f\n",
                g_qR0, g_qP0v, g_qP0v / std::sqrt(g_qP0v * g_qP0v + QK_M[0] * QK_M[0]));
    std::printf("after 400k steps E %.4f  drift %.3e   E band %.3f..%.3f\n",
                E1, (E1 - E0) / E0, emin, emax);
    std::printf("string L %.3f fm (%.3f..%.3f)   rrms %.3f fm (%.3f..%.3f)\n",
                L, lmin, lmax, g_qRrms, rmin, rmax);
}

int main(int argc, char** argv)
{
    build_systems();
    g_px.assign((size_t)FW * FH, 0);
    const char* mode = argc > 1 ? argv[1] : "all";
    int NS = argc > 2 ? atoi(argv[2]) : 1200;
    if (!std::strcmp(mode, "png")) {
        int sys = argc > 2 ? atoi(argv[2]) : 7;
        int nst = argc > 3 ? atoi(argv[3]) : 400;
        int orb = argc > 4 ? atoi(argv[4]) : 0;
        int foc = argc > 5 ? atoi(argv[5]) : 0;
        g_sysid = sys;
        sim_init(0, 0);
        sim_set(1, orb);
        g_focus = foc; layout();
        // render every frame, like the browser does, so the accumulated
        // electron-density buffer actually builds up
        uint8_t* px = nullptr;
        for (int f = 0; f < nst; ++f) { sim_step(1); px = sim_render(); }
        char fn[128]; std::snprintf(fn, sizeof(fn), "atom_os_%s_o%d_f%d.png", SYS[sys].name, orb, foc);
        stbi_write_png(fn, FW, FH, 4, px, FW * 4);
        std::printf("wrote %s   %s  Nel %d  E %.4f (exact %.4f)\n", fn, SYS[sys].name, g_nel,
                    g_ecnt ? g_esum / g_ecnt : 0.0, SYS[sys].Eref);
        return 0;
    }
    if (!std::strcmp(mode, "opt")) { emit_param_table(); return 0; }
    if (!std::strcmp(mode, "diag")) {
        // one hydrogen 1s orbital with a KNOWN exponent: the sampler must give
        // <1/r> = zeta and <r> = 3/(2 zeta) exactly.  Any deviation is a bug in
        // the Metropolis drift-diffusion step, not statistics.
        double zs = argc > 2 ? atof(argv[2]) : 1.08;
        long NS = argc > 3 ? atol(argv[3]) : 400000;
        g_rs = 88172645463325252ull; g_ghas = false;
        g_sysid = 0; g_zs[1] = zs; g_zscale = zs; g_jb = 0.8;
        setup_electrons(0);
        g_nw = 8; vmc_reset_walkers();
        for (int s = 0; s < 4000; ++s) for (int i = 0; i < g_nw; ++i) vmc_sweep(g_w[i]);
        double s_ir = 0, s_r = 0, s_e = 0; long c = 0;
        for (long s = 0; s < NS; ++s)
            for (int i = 0; i < g_nw; ++i) {
                vmc_sweep(g_w[i]);
                double r = std::sqrt(g_w[i].r[0][0] * g_w[i].r[0][0] + g_w[i].r[0][1] * g_w[i].r[0][1]
                                   + g_w[i].r[0][2] * g_w[i].r[0][2]);
                s_ir += 1.0 / r; s_r += r; s_e += local_energy(g_w[i]); c++;
            }
        double z = zs;
        std::printf("zeta %.4f   tau0 %.3f\n", z, g_tau0);
        std::printf("  <1/r>  %10.5f   exact %10.5f   err %+.3f%%\n", s_ir / c, z, 100.0 * (s_ir / c - z) / z);
        std::printf("  <r>    %10.5f   exact %10.5f   err %+.3f%%\n", s_r / c, 1.5 / z, 100.0 * (s_r / c - 1.5 / z) / (1.5 / z));
        std::printf("  <E_L>  %10.5f   exact %10.5f   accept %.1f%%\n", s_e / c, 0.5 * z * z - z,
                    100.0 * g_acc / g_att);
        return 0;
    }
    if (!std::strcmp(mode, "scan")) {
        // show how E/A and R respond to the packet width, at QCP = 0
        for (double L : {0.7, 0.9, 1.2, 1.6, 2.0, 2.5, 3.0, 3.6}) {
            QL = L; QCP = 0.0;
            std::printf("--- QL %.2f ---\n", L);
            test_nuclei();
        }
        return 0;
    }
    if (!std::strcmp(mode, "e"))        test_electrons(NS);
    else if (!std::strcmp(mode, "n"))   test_nuclei();
    else if (!std::strcmp(mode, "dyn")) test_nuc_dynamics();
    else if (!std::strcmp(mode, "q"))   test_quarks();
    else if (!std::strcmp(mode, "cal")) calibrate_qmd();
    else { test_electrons(NS); test_nuclei(); test_quarks(); }
    return 0;
}
#endif
