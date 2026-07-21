# 🌌 Universe OS — 物理シミュレーション集 (C++ / WASM)

**日本語** | [English](README.en.md)

24 個の物理シミュレーションを、**計算も描画もすべて C++ で実装して WebAssembly にコンパイル**したものです。
JavaScript はブラウザに必須の部分（ピクセルバッファの canvas への転送・UI イベント）だけを担当します。描画には単一ヘッダの 2D グラフィックスライブラリ [olive.c](https://github.com/tsoding/olive.c)（MIT）を使い、C++ 側が RGBA フレームバッファへ直接描いています。

元になった JavaScript 版は [yomei-o.github.io/universe](https://github.com/yomei-o/yomei-o.github.io/tree/main/universe) にあります。

## 🎮 オンラインデモ

**[▶ ギャラリー（全デモ一覧）](https://yomei-o.github.io/universe_wasm/)**

| # | デモ | 内容 |
|---|------|------|
| 🌌 | [重力波](https://yomei-o.github.io/universe_wasm/gravity_wave/) | 回転連星（四重極）が放つ時空の波・吸収境界 |
| 💥 | [反物質消滅](https://yomei-o.github.io/universe_wasm/annihilation_os/) | 物質・反物質の対消滅とエネルギー波 |
| 🌊 | [二重スリット](https://yomei-o.github.io/universe_wasm/doubleslit_os/) | 干渉縞と観測による収束 |
| 🔗 | [量子もつれ](https://yomei-o.github.io/universe_wasm/entanglement_os/) | 反位相の対・測定でスピン確定 |
| 🚧 | [トンネル効果](https://yomei-o.github.io/universe_wasm/tunneling_os_v4_1/) | 障壁を透過する波束・透過カウンタ |
| 🕳️ | [ワームホール](https://yomei-o.github.io/universe_wasm/wormhole_os/) | 2 領域を直結する振幅転送 |
| ⚛️ | [波と粒](https://yomei-o.github.io/universe_wasm/wave_particle_os/) | 波の干渉と観測での粒子ドット化 |
| 🧲 | [超伝導・マイスナー](https://yomei-o.github.io/universe_wasm/superconductivity_os/) | 冷却で無抵抗・磁場排除 |
| ⭐ | [プランク星](https://yomei-o.github.io/universe_wasm/planck_star_os/) | 重力崩壊と特異点回避のバウンス |
| 🎆 | [インフレーション](https://yomei-o.github.io/universe_wasm/inflation_os/) | 量子ゆらぎの指数膨張と CMB 凍結 |
| ☢️ | [中性子崩壊](https://yomei-o.github.io/universe_wasm/neutron_decay_os/) | β崩壊（波モード／粒子モード） |
| 🕸️ | [時空メッシュ](https://yomei-o.github.io/universe_wasm/universe_os/) | クリックしたブラックホールで格子が歪む |
| 📐 | [カシミール効果](https://yomei-o.github.io/universe_wasm/casimir_vacuum_os/) | 真空ゆらぎの圧力差で板が動く |
| 🌫️ | [電子雲](https://yomei-o.github.io/universe_wasm/electron_os/) | 核に束縛された電子軌道（クリックで共有結合） |
| ⚡ | [弱い相互作用](https://yomei-o.github.io/universe_wasm/weak_os/) | 中性子→陽子＋電子（β崩壊） |
| ♾️ | [三体問題](https://yomei-o.github.io/universe_wasm/three_body_chaos_working/) | 8 の字周期軌道とカオス（誤差注入） |
| 🎗️ | [クォーク閉じ込め](https://yomei-o.github.io/universe_wasm/quark_os/) | 距離で硬くなる強い力・弦の破断で対生成 |
| 🌀 | [ダークマター](https://yomei-o.github.io/universe_wasm/dark_matter_os/) | 銀河回転曲線（ケプラー vs 平坦） |
| 🌠 | [グラビトン放出](https://yomei-o.github.io/universe_wasm/graviton_spiral_os/) | 連星から螺旋状に放出される重力子 |
| 💧 | [水分子振動](https://yomei-o.github.io/universe_wasm/water_vibration_os/) | H₂O の 3 基準振動モード |
| ⏱️ | [時間遅延・双子](https://yomei-o.github.io/universe_wasm/time_dilation_os/) | 高速移動する双子の時計の遅れ |
| 🌿 | [光合成・量子輸送](https://yomei-o.github.io/universe_wasm/photosynthesis_os/) | 導波路を低損失で伝わり反応中心へ集約する波 |
| 🌀 | [化学スパイラル](https://yomei-o.github.io/universe_wasm/spiral_os/) | Barkley 反応拡散の回転スパイラル（クリックで注入） |
| 🧊 | [開放系超伝導](https://yomei-o.github.io/universe_wasm/superconductor_os/) | クーパー対の輸送・温度で対が壊れ抵抗が出る |

## 💡 構成

- **計算コア**：各シミュレーションの状態更新（波動 PDE の 5 点ラプラシアン＋leapfrog、N 体積分、粒子系など）はすべて C++。
- **描画**：C++ が olive.c を使って RGBA(uint32) フレームバッファへ直接描画（塗り・線・円・文字）。olive の色は `0xAABBGGRR`（メモリ上は R,G,B,A バイト）で、canvas の `ImageData` にそのまま転送できます。
- **JS ハーネス**：全デモで共通の薄い JS。`sim.js`（Emscripten 出力）を読み込み、毎フレーム `sim_step()`→`sim_render()`→canvas へ blit し、スライダー／ボタン／クリックを C++ に転送するだけ。
- **共通 ABI**（全 sim が実装）:
  `sim_init` / `sim_w` / `sim_h` / `sim_reset` / `sim_step` / `sim_render` / `sim_click` / `sim_set` / `sim_action`

多くの波動系は同一の差分カーネル `next = 2·curr − prev + Ω²·∇²curr`（＋各 sim 固有の減衰・境界・注入）を使っています。

## 🛠 ビルド

[Emscripten](https://emscripten.org/) が必要です。

```bash
# 1 つビルド（例: 重力波）
EMSDK=/path/to/emsdk ./build.sh gravity_wave
# 出力: wasmdist/gravity_wave/sim.js (+.wasm)

# 全部ビルド
for f in src/*.cpp; do ./build.sh "$(basename "$f" .cpp)"; done
```

各 `src/<name>.cpp` は `#ifndef __EMSCRIPTEN__` のネイティブ自己テスト（1 フレーム描いて PNG 出力）も持っており、`clang++ -I. src/<name>.cpp` で単体検証できます（※ olive.c の C99 指定イニシャライザは g++ では通らないため clang 系を使用）。

## 📁 ファイル構成

| パス | 役割 |
|------|------|
| `src/<name>.cpp` | 各シミュレーションの計算＋olive.c 描画（共通 ABI）|
| `olive.c` | 単一ヘッダ 2D グラフィックスライブラリ（MIT, [tsoding/olive.c](https://github.com/tsoding/olive.c)）|
| `stb_image_write.h` | ネイティブ自己テストの PNG 出力用（public domain）|
| `build.sh` | 1 つの sim を Emscripten でビルド |
| `wasmdist/<name>/` | ビルド済みデモ（`index.html` + `sim.js` + `sim.wasm`）|
| `wasmdist/index.html` | 全デモへのギャラリー |

## 📝 ライセンス / 注意

- 移植コードはオリジナル。物理モデルは元 JS 版（[universe](https://github.com/yomei-o/yomei-o.github.io/tree/main/universe)）に忠実。
- olive.c は MIT、stb は public domain。各ヘッダのライセンス表記を参照。
- canvas 描画フォントは olive の内蔵フォント（ASCII）を使用するため、オーバーレイ文字は英字です（説明文は各 index.html に日本語で記載）。
