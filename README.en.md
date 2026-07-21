# 🌌 Universe OS — physics simulations (C++ / WASM)

[日本語](README.md) | **English**

21 physics simulations where **both the computation and the rendering are written in C++** and compiled to WebAssembly. JavaScript only does what the browser requires (blitting the pixel buffer to a canvas, forwarding UI events). Rendering uses the single-header 2D graphics library [olive.c](https://github.com/tsoding/olive.c) (MIT); the C++ side draws straight into an RGBA framebuffer.

The original JavaScript versions live at [yomei-o.github.io/universe](https://github.com/yomei-o/yomei-o.github.io/tree/main/universe).

## 🎮 Live demos

**[▶ Gallery (all demos)](https://yomei-o.github.io/universe_wasm/)**

| # | Demo | What it shows |
|---|------|---------------|
| 🌌 | [Gravity waves](https://yomei-o.github.io/universe_wasm/gravity_wave/) | Rotating binary (quadrupole) radiating spacetime waves; absorbing boundary |
| 💥 | [Annihilation](https://yomei-o.github.io/universe_wasm/annihilation_os/) | Matter/antimatter annihilation into an energy wave |
| 🌊 | [Double slit](https://yomei-o.github.io/universe_wasm/doubleslit_os/) | Interference fringes collapsing under observation |
| 🔗 | [Entanglement](https://yomei-o.github.io/universe_wasm/entanglement_os/) | Anti-phase pair; measuring one fixes both spins |
| 🚧 | [Tunneling](https://yomei-o.github.io/universe_wasm/tunneling_os_v4_1/) | Wave packet through a barrier; transmission counter |
| 🕳️ | [Wormhole](https://yomei-o.github.io/universe_wasm/wormhole_os/) | Amplitude transfer linking two isolated regions |
| ⚛️ | [Wave–particle](https://yomei-o.github.io/universe_wasm/wave_particle_os/) | Wave interference vs. probabilistic particle dots |
| 🧲 | [Superconductivity](https://yomei-o.github.io/universe_wasm/superconductivity_os/) | Zero resistance + Meissner field expulsion on cooling |
| ⭐ | [Planck star](https://yomei-o.github.io/universe_wasm/planck_star_os/) | Gravitational collapse with a singularity-avoiding bounce |
| 🎆 | [Inflation](https://yomei-o.github.io/universe_wasm/inflation_os/) | Exponential stretch of quantum noise, freezing into CMB |
| ☢️ | [Neutron decay](https://yomei-o.github.io/universe_wasm/neutron_decay_os/) | Beta decay (wave mode / particle mode) |
| 🕸️ | [Spacetime mesh](https://yomei-o.github.io/universe_wasm/universe_os/) | Click a black hole to warp the lattice |
| 📐 | [Casimir effect](https://yomei-o.github.io/universe_wasm/casimir_vacuum_os/) | Vacuum-fluctuation pressure difference moving plates |
| 🌫️ | [Electron cloud](https://yomei-o.github.io/universe_wasm/electron_os/) | Electron orbitals; click to add a nucleus (covalent) |
| ⚡ | [Weak interaction](https://yomei-o.github.io/universe_wasm/weak_os/) | Neutron → proton + electron (beta decay) |
| ♾️ | [Three-body](https://yomei-o.github.io/universe_wasm/three_body_chaos_working/) | Figure-8 choreography and chaos (error injection) |
| 🎗️ | [Quark confinement](https://yomei-o.github.io/universe_wasm/quark_os/) | Distance-stiffening strong force; string snap = pair production |
| 🌀 | [Dark matter](https://yomei-o.github.io/universe_wasm/dark_matter_os/) | Galaxy rotation curve (Keplerian vs flat) |
| 🌠 | [Graviton spiral](https://yomei-o.github.io/universe_wasm/graviton_spiral_os/) | Gravitons emitted in a spiral from a binary |
| 💧 | [H₂O vibration](https://yomei-o.github.io/universe_wasm/water_vibration_os/) | The three normal vibrational modes of water |
| ⏱️ | [Time dilation](https://yomei-o.github.io/universe_wasm/time_dilation_os/) | Twin paradox — the fast twin's clock runs slow |

## 💡 Architecture

- **Compute core**: every sim's state update (wave PDE via 5-point Laplacian + leapfrog, N-body integration, particle systems, …) is in C++.
- **Rendering**: C++ draws with olive.c straight into an RGBA (uint32) framebuffer (fills, lines, circles, text). olive colors are `0xAABBGGRR` (bytes R,G,B,A in memory), which map directly onto a canvas `ImageData`.
- **JS harness**: one thin, shared harness. It loads `sim.js` (Emscripten output) and every frame runs `sim_step()` → `sim_render()` → blit to canvas, forwarding sliders/buttons/clicks to C++.
- **Common ABI** (implemented by every sim):
  `sim_init` / `sim_w` / `sim_h` / `sim_reset` / `sim_step` / `sim_render` / `sim_click` / `sim_set` / `sim_action`

Most wave sims share one finite-difference kernel `next = 2·curr − prev + Ω²·∇²curr` (plus per-sim damping / boundary / injection).

## 🛠 Build

Requires [Emscripten](https://emscripten.org/).

```bash
# build one (e.g. gravity waves)
EMSDK=/path/to/emsdk ./build.sh gravity_wave
# -> wasmdist/gravity_wave/sim.js (+.wasm)

# build all
for f in src/*.cpp; do ./build.sh "$(basename "$f" .cpp)"; done
```

Each `src/<name>.cpp` also carries a native self-test under `#ifndef __EMSCRIPTEN__` (renders one frame to a PNG), so you can verify it standalone with `clang++ -I. src/<name>.cpp` (use clang — olive.c's C99 designated initializers don't compile under g++).

## 📁 Layout

| Path | Role |
|------|------|
| `src/<name>.cpp` | Per-sim compute + olive.c rendering (common ABI) |
| `olive.c` | Single-header 2D graphics library (MIT, [tsoding/olive.c](https://github.com/tsoding/olive.c)) |
| `stb_image_write.h` | PNG output for the native self-tests (public domain) |
| `build.sh` | Build one sim with Emscripten |
| `wasmdist/<name>/` | Prebuilt demo (`index.html` + `sim.js` + `sim.wasm`) |
| `wasmdist/index.html` | Gallery of all demos |

## 📝 License / note

- The ported code is original; the physics models stay faithful to the original JS versions ([universe](https://github.com/yomei-o/yomei-o.github.io/tree/main/universe)).
- olive.c is MIT; stb is public domain — see each header for its license.
- On-canvas overlay text uses olive's built-in (ASCII) font, so labels are in English; the Japanese copy lives in each `index.html`.
