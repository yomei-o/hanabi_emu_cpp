// Hanabi OS — 打ち上げ花火シミュレータ / 菊 (Chrysanthemum shell)   (C++ / WASM)
//
// 本物の花火玉の中身をそのまま数値計算する。
//   1) 打揚:  打揚薬で玉が上昇。玉は球体として二次空気抵抗(F = -1/2 rho Cd A |v| v)+重力で減速し、
//             昇り曲導の火の粉を曳きながら頂点へ。初速は目標到達高度から解析解で逆算する。
//               v0 = sqrt( (g/k) * (exp(2 k H) - 1) ),  k = 1/2 rho Cd A / m
//   2) 開発:  頂点付近で割薬が働き、玉皮の内側に一層(または芯入りで多層)に詰められた「星」が
//             球面状に飛散する。方向はフィボナッチ球で等方に配り、初速に数%のばらつきを与える。
//   3) 星:    星は1個ずつ独立に積分する。燃焼で半径 r が線形に減る(表面燃焼)ので、
//             質量 m = rho_s (4/3)pi r^3、断面積 A = pi r^2 が時々刻々変化し、
//             弾道係数 k = 1/2 rho Cd A / m = 0.8636 / (4 rho_s r) は r とともに増大する。
//             = 小さくなるほど急減速する → これが菊の「開いて止まって垂れる」形を作る。
//             輝度は燃焼表面積 (r/r0)^2 に比例させ、消え際に自然に暗くなる。
//   4) 引先:  菊の尾は星から剥がれ落ちる燃えかす(火の粉)そのもの。半径 0.3mm 程度の粒は
//             k ~ 0.3 と抗力が桁違いに大きく、放出直後にほぼ空気に対して静止し、
//             その場で燃えながら重力でゆっくり垂れる。これを実粒子として撃ち出している。
//
// 計算も描画もすべて C++。粒子は float の加算バッファ(HDR)に加算合成し、
// x/(1+x) のトーンマップで夜空にコンポジットする。テキストのみ olive.c。
#define OLIVEC_IMPLEMENTATION
#include "olive.c"
#include <vector>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <algorithm>
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define KEEP EMSCRIPTEN_KEEPALIVE
#else
#define KEEP
#endif

// ---------------------------------------------------------------- 定数
static const int FW = 960, FH = 600;
static const int SUBSTEPS = 4;                 // 1frame = 4 substep
static const double DT = 1.0 / (60.0 * SUBSTEPS);
static const double G_ACC = 9.80665;           // 重力加速度 [m/s^2]
static const double RHO_AIR = 1.225;           // 空気密度 [kg/m^3]
static const double CD = 0.47;                 // 球の抗力係数
static const double RHO_STAR = 1750.0;         // 星の密度 [kg/m^3]
static const double RHO_SHELL = 800.0;         // 玉の平均密度 [kg/m^3]
static const size_t MAX_SPARKS = 560000;   // 連打(スターマイン)を捌くための粒子予算
static const size_t MAX_STARS = 120000;

// ---------------------------------------------------------------- 乱数
static uint32_t rng = 88675123u;
static inline double rnd() { rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5; return (rng & 0xFFFFFF) / (double)0x1000000; }
static inline double rnds() { return rnd() * 2.0 - 1.0; }
static inline float smoothstep(float a, float b, float x) {
    float t = (x - a) / (b - a);
    if (t < 0) t = 0; if (t > 1) t = 1;
    return t * t * (3.0f - 2.0f * t);
}
static inline uint32_t rgb(int r, int g, int b) {
    return 0xFF000000u | ((uint32_t)b << 16) | ((uint32_t)g << 8) | (uint32_t)r;
}

// ---------------------------------------------------------------- 玉の諸元
// 号 = 寸。3号 = 90mm、10号 = 尺玉 = 300mm。
//
// 入力は「火薬が与える初速」と「導火線の秒時」だけ。到達高度も開花直径も
// 積分(差分式)の結果として出てくる — 逆算ソルバや解析解は一切使わない。
struct ShellSpec {
    double go;        // 号数
    double dia;       // 玉直径 [m]
    double mass;      // 玉質量 [kg]
    double liftV;     // 打揚薬が与える初速 [m/s]
    double fuseT;     // 時限導火線の秒時 [s]
    double starV;     // 割薬が星に与える初速 [m/s]
    double starR;     // 星の半径 [m]
    double burnT;     // 星の燃焼時間 [s]
    int    nstar;     // 一層あたりの星数
    double nomH;      // 画角決めに使う想定高度 [m] (表示や物理には使わない)
    double nomR;      // 画角決めに使う想定半径 [m]
};
static ShellSpec spec_of(double go) {
    ShellSpec s;
    s.go = go;
    s.dia = go * 0.03;                          // 3号=90mm … 10号(尺)=300mm
    s.mass = RHO_SHELL * (4.0 / 3.0) * M_PI * pow(s.dia * 0.5, 3.0);
    s.liftV = 100.0 + 1.6 * go;                 // 臼砲の初速 105〜116 m/s
    s.fuseT = 3.90 + 0.37 * go;                 // 秒時 5.0〜7.6s (頂点付近で開発)
    s.starV = 41.0 + 6.3 * go;                  // 割薬 60〜104 m/s
    s.starR = (1.4 * go + 0.3) * 0.001;         // 3号 4.5mm … 10号 14.3mm
    s.burnT = 2.15 + 0.485 * go;                // 3号 3.6s … 10号 7.0s
    s.nstar = (int)(14.0 * pow(go, 1.7) + 0.5); // 3号 91 … 10号 700
    s.nomH = 100.0 + 24.0 * go;
    s.nomR = 16.0 + 14.3 * go;
    return s;
}

// 星の弾道係数 k = 1/2 rho Cd (pi r^2) / (rho_s 4/3 pi r^3) = 0.8636 / (4 rho_s r)
static inline double star_k(double r) { return (0.5 * RHO_AIR * CD * 3.0) / (4.0 * RHO_STAR * r); }

// ---------------------------------------------------------------- カメラ
// ピンホール。地平線を画面下端(1.02*FH)、花の頂点を画面上端(0.06*FH)に置くよう
// 焦点距離 f とピッチ角 p を決める。二式の比から t = tan(p) の二次方程式になるので
// 反復せず一発で解ける(物理ではなく画角の話なので毎フレームは要らない)。
//     q*ry*t^2 + (1+q)*D*t - ry = 0,   q = 0.44/0.52,  ry = 花の頂点高さ - 視点高さ
static double camD = 500, camPitch = 0.22, camF = 1400, camCP = 1, camSP = 0;
static const double CAM_EYE = 1.6;

static inline bool project(double x, double y, double z, double& sx, double& sy, double& sc) {
    double ry = y - CAM_EYE, rz = z + camD;
    double yp = ry * camCP - rz * camSP;
    double zp = ry * camSP + rz * camCP;
    if (zp < 1.0) return false;
    sc = camF / zp;
    sx = FW * 0.5 + x * sc;
    sy = FH * 0.5 - yp * sc;
    return true;
}
static void setup_camera(const ShellSpec& s) {
    camD = 2.5 * s.nomH;                             // 観客の距離 = 想定高度の2.5倍
    double ry = s.nomH + s.nomR * 1.18 - CAM_EYE;
    const double A = 0.52, B = 0.44, q = B / A;
    double a = q * ry, b = (1.0 + q) * camD, c = -ry;
    double t = (-b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
    camPitch = atan(t);
    camF = (A * FH) / t;
    camCP = cos(camPitch); camSP = sin(camPitch);
}

// ---------------------------------------------------------------- 粒子
struct Vec3 { double x, y, z; };
struct Col { float r, g, b; };

struct Star {
    double x, y, z, vx, vy, vz;
    float bx, by, bz;      // 開発点(開花直径の実測用)
    double r, r0, dr;      // 現在半径 / 初期半径 / 燃焼速度 [m/s]
    double shed;           // 火の粉放出アキュムレータ
    Col c0, c1;            // 変化前 / 変化後の炎色
    float chg;             // 色が変わる燃焼進行度 (1.0 なら変化なし)
    Col tail;              // 引先(尾)の色
    float shedRate;        // [個/s]
    // --- 蜂(ハチ): ノズルの推力ベクトルが軸まわりに回る ---
    float ux, uy, uz;      // 推力が回る平面の基底 u
    float wx, wy, wz;      //                        w  (u ⊥ w ⊥ 回転軸)
    float thrust;          // ノズルの推力/質量 [m/s^2] (0 なら普通の星)
    float spin;            // 回転角速度 [rad/s] — 空気抵抗でだんだん落ちる
    float spinB;           // 回転の抵抗係数 [1/rad]  dω/dt = -B|ω|ω
    float sphase;          // 現在の位相 [rad]
    float spinDelay;       // ノズルに火が回るまでの時間 [s]。それまでは普通の星として飛ぶ
    float exhaust;         // 噴射速度 [m/s]。火の粉はノズルの反対向きにこの速さで出る
    float delay;           // 暗飛行(点火までの時間)[s] — 飛んではいるが光らない
    bool  crackle;         // 消え際に分砲(パチパチ)する星か
    bool  alive;
};
// 蜂・渦蜂 — ノズル付きの星。彗星のように噴射し、その反作用で飛ぶ。
//   機体が自転しているのでノズルの向き n̂ が回り、推力の向きも一緒に回る。
//     n̂(θ) = u cosθ + w sinθ,   θ ← θ + ω dt,   dω/dt = -B|ω|ω
//     a  += (T r/r₀) n̂                      ← 軌跡は積分の結果。押しつけない
//     噴射(火の粉)は反対向き: v_spark = v_star - v_ex n̂
//   結果として速度は振幅 a/ω で振れ、位置は a/ω² の輪を描く。
//   ω が速い(渦蜂)ほど輪は小さくなり、ほぼ直進しながら噴射だけが回る。
struct Star;
static void nozzle_dir(const Star& st, double phase, double& nx, double& ny, double& nz);

struct Spark {
    float x, y, z, vx, vy, vz;
    float k;               // 弾道係数
    float life, life0;
    Col c;
    bool alive;
};
struct Shell {
    double x, y, z, vx, vy, vz;
    double k;              // 玉の弾道係数(一定)
    double shed;
    double fuse;           // 時限導火線の残り [s] — 毎ステップ dt だけ焼ける
    ShellSpec sp;
    int layers;
    int scheme;
    int type;              // 玉の種類(下の enum)
    int gen;               // 世代 (0=親玉, 1=千輪の子玉)
    int fixc;              // テーマで炎色を固定するとき FIXCOL の添字。-1 で自由
    bool alive;
};

// 玉の種類。見た目を直接描くのではなく、星に与える物理量を変えているだけ。
enum {
    T_KIKU = 0,   // 菊     — 引先を長く曳く。基本形
    T_BOTAN,      // 牡丹   — 尾なし。点の集まりで開く
    T_YANAGI,     // 柳     — 低速・長燃焼で垂れ下がる
    T_HACHI,      // 蜂     — 星が自分で推力を出し、その向きが回る = 思い思いにクルクル飛ぶ
    T_UZU,        // 渦蜂   — 回転軸を玉全体で揃えた蜂。強い推力＋ゆっくりの回転で大きな渦を巻く
    T_SENRIN,     // 千輪   — 子玉を撒き、少し遅れて一斉に小さく開く
    T_HEART,      // 型物   — ハート
    T_RING,       // 型物   — 輪
    T_SATURN,     // 型物   — 土星(芯＋傾いた輪)
    T_KAMURO,     // 錦冠   — 柳の親玉。10秒級の長燃焼で地面近くまで垂れ、消え際だけ紅くなる
    T_COUNT
};
static inline bool is_shaped(int t) { return t == T_HEART || t == T_RING || t == T_SATURN; }
struct Flash { double x, y, z; double t, t0, pw; bool alive; };

// 予約発射 — t 秒後に上げる。大玉と小玉は秒時が違うので、同時に「開かせる」には
// 開発までの時間の差だけ小玉を遅らせて撃つ必要がある(長岡のフェニックスがこれ)
struct Pending { double t, nx, go, fuseScale, fuseJit; int type, fixc; };

// 位相 phase のときのノズルの向き(=推力の向き)。噴射はこの反対に出る
static void nozzle_dir(const Star& st, double phase, double& nx, double& ny, double& nz) {
    double cs = cos(phase), sn = sin(phase);
    nx = st.ux * cs + st.wx * sn;
    ny = st.uy * cs + st.wy * sn;
    nz = st.uz * cs + st.wz * sn;
}

static std::vector<Star>  stars;
static std::vector<Spark> sparks;
static std::vector<Shell> shells;
static std::vector<Flash> flashes;
static std::vector<Pending> pending;

static std::vector<float>    acc;   // HDR 加算バッファ (RGB)
static std::vector<uint32_t> px;    // 出力 RGBA
static std::vector<uint32_t> bg;    // 夜空(静的)

// ---------------------------------------------------------------- パラメータ
static double p_go = 5.0;        // 号数
static double p_density = 1.0;   // 星数スケール
static double p_interval = 2.2;  // 打ち上げ間隔 [s]
static double p_wind = 1.5;      // 風速 [m/s] (+x)
static double p_shed = 1.0;      // 火の粉の量
static double p_layers = 1.0;    // 芯の数 (0=菊, 1=芯入, 2=八重芯)
static double p_glow = 0.82;     // 残光(加算バッファの減衰)
static bool   p_auto = true;

static double p_shots = 100.0;   // スターマインの発数
static double p_rapid = 0.22;    // 連打間隔 [s]
static double p_quiet = 6.0;     // 連打後の余韻(静寂)の長さ [s]
static double p_starInt = 120.0; // 自動スターマインの間隔 [s] (0 で自動なし)
static double p_hud = 1.0;       // 画面内の文字表示 (0 で花火だけ)
static double p_type = -1.0;     // 単発・自動連発の玉の種類 (-1 = おまかせ)
static double p_theme = 0.0;     // スターマインのテーマ (下の THEMES の添字)

static double simTime = 0.0, launchTimer = 0.5;
static ShellSpec curSpec;
static long   frameNo = 0;

// スターマイン(早打ち)の状態 — 全部インクリメンタルに減っていく
// barPhase: 0=通常 1=連打中 2=余韻(空が空になるまで待ち、さらに静寂の時間だけ何も上げない)
static int    barLeft = 0, barTotal = 0, barPhase = 0;
static double barTimer = 0.0, quietTimer = 0.0;
static double starTimer = 120.0; // 次の自動スターマインまで [s] — 通常時だけ減っていく
static int    curTheme = 1;      // この連打で実際に使うテーマ(おまかせはここで確定させる)

// 粒子予算に対する負荷から引先の量を自動調整する係数(一次遅れ = 差分式)
static double qual = 1.0;
// 実測値(毎ステップ running max で更新する = これも差分式)
static double measH = 0.0, measD = 0.0;

// ---------------------------------------------------------------- 炎色
// 実際の炎色反応。Sr=紅, Ba=緑, Cu=青, Na=黄, Mg/Al=銀, 木炭=金
static const Col C_GOLD   = { 1.00f, 0.62f, 0.20f };
static const Col C_AMBER  = { 1.00f, 0.42f, 0.08f };
static const Col C_RED    = { 1.00f, 0.16f, 0.14f };
static const Col C_GREEN  = { 0.26f, 1.00f, 0.34f };
static const Col C_BLUE   = { 0.28f, 0.52f, 1.00f };
static const Col C_SILVER = { 1.00f, 0.97f, 0.90f };
static const Col C_LEMON  = { 1.00f, 0.93f, 0.30f };
static const Col C_PURPLE = { 0.72f, 0.36f, 1.00f };
static const Col C_PINK   = { 1.00f, 0.48f, 0.66f };   // 桜色
// テーマで色を固定するときの表
static const Col FIXCOL[] = { C_PINK, C_RED, C_GREEN, C_BLUE, C_LEMON, C_SILVER, C_PURPLE, C_GOLD };

static inline Col mixc(Col a, Col b, float t) {
    return { a.r + (b.r - a.r) * t, a.g + (b.g - a.g) * t, a.b + (b.b - a.b) * t };
}
// 色替わり(変化菊)を含む配色表
static void scheme_colors(int s, int layer, Col& c0, Col& c1, float& chg, Col& tail) {
    static const Col tbl[] = { C_RED, C_GREEN, C_BLUE, C_LEMON, C_SILVER, C_PURPLE };
    switch (s) {
    case 0: // 錦冠 — 木炭の金一色。菊の基本形
        c0 = C_GOLD; c1 = C_AMBER; chg = 0.55f; tail = C_GOLD; break;
    case 1: // 単色 + 金の引先
        c0 = tbl[(int)(rnd() * 6)]; c1 = c0; chg = 1.0f; tail = mixc(c0, C_GOLD, 0.38f); break;
    case 2: // 変化菊 — 途中で色が変わる
        c0 = tbl[(int)(rnd() * 6)]; c1 = tbl[(int)(rnd() * 6)]; chg = 0.5f;
        tail = mixc(c0, C_GOLD, 0.35f); break;
    default: // 芯替わり — 層ごとに別の色
        c0 = tbl[(layer * 2 + 1) % 6]; c1 = c0; chg = 1.0f;
        tail = mixc(c0, C_GOLD, 0.35f); break;
    }
    if (s >= 3 && layer == 0) { c0 = C_GOLD; c1 = C_AMBER; chg = 0.6f; tail = C_GOLD; }
}

// ---------------------------------------------------------------- 夜空
static void build_bg() {
    bg.assign((size_t)FW * FH, 0);
    // 地平線の画面 y (無限遠・地面)
    double horiz = FH * 0.5 + camF * (camSP / camCP);
    uint32_t seed = rng; rng = 20240807u;
    for (int j = 0; j < FH; ++j) {
        double t = (double)j / FH;
        int r = (int)(4 + 10 * t * t), g = (int)(6 + 13 * t * t), b = (int)(16 + 26 * t * t);
        uint32_t c = rgb(r, g, b);
        for (int i = 0; i < FW; ++i) bg[(size_t)j * FW + i] = c;
    }
    // 星
    for (int n = 0; n < 700; ++n) {
        int i = (int)(rnd() * FW), j = (int)(rnd() * (horiz < FH ? horiz : FH));
        if (j < 0 || j >= FH) continue;
        double br = rnd(); br = br * br * br;
        int v = (int)(40 + 215 * br);
        uint32_t& d = bg[(size_t)j * FW + i];
        int r = std::min(255, (int)(d & 255) + v), g = std::min(255, (int)((d >> 8) & 255) + v);
        int b = std::min(255, (int)((d >> 16) & 255) + (int)(v * 1.05));
        d = rgb(r, g, b);
    }
    // 地面と稜線
    for (int i = 0; i < FW; ++i) {
        double u = i * 0.017;
        double h = 14.0 * sin(u * 0.31) + 9.0 * sin(u * 0.73 + 1.7) + 5.0 * sin(u * 1.9 + 0.4)
                 + 3.0 * sin(u * 4.1 + 2.2);
        int top = (int)(horiz - 16.0 - h);
        for (int j = std::max(0, top); j < FH; ++j) {
            double d = (j - top) / 40.0; if (d > 1) d = 1;
            bg[(size_t)j * FW + i] = rgb((int)(3 + 3 * (1 - d)), (int)(4 + 4 * (1 - d)), (int)(9 + 6 * (1 - d)));
        }
        // 街あかり
        if (top >= 0 && top < FH && rnd() < 0.02) {
            int v = 60 + (int)(rnd() * 120);
            for (int j = top; j < std::min(FH, top + 2); ++j)
                bg[(size_t)j * FW + i] = rgb(v, (int)(v * 0.75), (int)(v * 0.35));
        }
    }
    rng = seed;
}

// ---------------------------------------------------------------- 加算描画
static inline void add_px(int i, int j, float r, float g, float b) {
    if (i < 0 || j < 0 || i >= FW || j >= FH) return;
    float* p = &acc[((size_t)j * FW + i) * 3];
    p[0] += r; p[1] += g; p[2] += b;
}
// サブピクセル(バイリニア)で1点を加算 — 火の粉用
static inline void splat_point(double sx, double sy, Col c, float amp) {
    int ix = (int)floor(sx), iy = (int)floor(sy);
    if (ix < 0 || iy < 0 || ix >= FW - 1 || iy >= FH - 1) return;
    float fx = (float)(sx - ix), fy = (float)(sy - iy);
    float w[4] = { (1 - fx) * (1 - fy), fx * (1 - fy), (1 - fx) * fy, fx * fy };
    const int ox[4] = { 0,1,0,1 }, oy[4] = { 0,0,1,1 };
    for (int k = 0; k < 4; ++k) {
        float a = w[k] * amp;
        add_px(ix + ox[k], iy + oy[k], c.r * a, c.g * a, c.b * a);
    }
}
// ガウス円盤 — 星頭・閃光用
static inline void splat_disc(double sx, double sy, float rad, Col c, float amp) {
    if (rad < 0.8f) { splat_point(sx, sy, c, amp); return; }
    int R = (int)ceil(rad * 2.0f);
    if (R > 64) R = 64;
    int cx = (int)floor(sx), cy = (int)floor(sy);
    if (cx + R < 0 || cy + R < 0 || cx - R >= FW || cy - R >= FH) return;
    float inv = 1.0f / (rad * rad);
    for (int dy = -R; dy <= R; ++dy) for (int dx = -R; dx <= R; ++dx) {
        float px_ = (float)(cx + dx + 0.5 - sx), py_ = (float)(cy + dy + 0.5 - sy);
        float d2 = (px_ * px_ + py_ * py_) * inv;
        if (d2 > 4.0f) continue;
        float w = expf(-d2 * 1.6f) * amp;
        add_px(cx + dx, cy + dy, c.r * w, c.g * w, c.b * w);
    }
}

// ---------------------------------------------------------------- 生成
static void emit_spark(double x, double y, double z, double vx, double vy, double vz, Col c, double sz) {
    if (sparks.size() >= MAX_SPARKS) return;
    Spark s;
    s.x = (float)x; s.y = (float)y; s.z = (float)z;
    s.vx = (float)vx; s.vy = (float)vy; s.vz = (float)vz;
    double r = sz * (0.6 + 0.8 * rnd());
    s.k = (float)((0.5 * RHO_AIR * CD * 3.0) / (4.0 * RHO_STAR * r));
    s.life0 = s.life = (float)(0.5 + 1.2 * rnd());
    s.c = c;
    s.alive = true;
    sparks.push_back(s);
}

// 型物の輪郭。u∈[0,1) を受けて、観客を向いた平面上の点(半径1に正規化)を返す
static void shape_point(int type, double u, double& sx, double& sy) {
    double t = u * 6.2831853;
    if (type == T_HEART) {
        double s3 = sin(t);
        sx = 16.0 * s3 * s3 * s3;
        sy = 13.0 * cos(t) - 5.0 * cos(2 * t) - 2.0 * cos(3 * t) - cos(4 * t);
        sx /= 17.0; sy /= 17.0;
    } else if (type == T_RING) {
        sx = cos(t); sy = sin(t);
    } else {                       // 土星の輪 — 寝かせて傾ける(芯は burst 側で3Dの球にする)
        double rx = cos(t), ry = sin(t) * 0.26;
        const double ti = 0.30;
        sx = rx * cos(ti) - ry * sin(ti);
        sy = rx * sin(ti) + ry * cos(ti);
    }
}

static void spawn_child(const Shell& p, double dx, double dy, double dz, double v,
                        double go, double fuse, int type);

static void burst(Shell& sh) {
    const ShellSpec& s = sh.sp;

    // 千輪 — 星ではなく子玉を撒く。子玉は短い秒時で、少し遅れて一斉に小さく開く
    if (sh.type == T_SENRIN) {
        int nc = (int)(9.0 + 4.0 * s.go);
        double cgo = std::max(1.4, s.go * 0.26);
        for (int i = 0; i < nc; ++i) {
            double yy = 1.0 - 2.0 * (i + 0.5) / nc;
            double rr = sqrt(std::max(0.0, 1.0 - yy * yy));
            double ph = i * 2.39996322973 + rnd();
            spawn_child(sh, rr * cos(ph), yy, rr * sin(ph),
                        s.starV * 0.85 * (1.0 + rnds() * 0.06), cgo,
                        1.45 + rnd() * 0.55, rnd() < 0.5 ? T_BOTAN : T_KIKU);
        }
        if (flashes.size() < 64) {
            Flash f; f.x = sh.x; f.y = sh.y; f.z = sh.z;
            f.t0 = f.t = 0.10; f.pw = 3.2 * s.go; f.alive = true;
            flashes.push_back(f);
        }
        sh.alive = false;
        return;
    }

    // 玉種による味付け。物理量(初速・燃焼時間・剥離量・推力)を変えるだけで見た目は勝手に変わる
    double tv = 1.0, tb = 1.0, tshed = 1.0, tthrust = 0.0;
    if (sh.type == T_BOTAN)  { tv = 1.10; tb = 0.75; tshed = 0.0; }   // 牡丹 — 尾を引かない
    if (sh.type == T_YANAGI) { tv = 0.48; tb = 1.95; tshed = 1.7; }   // 柳   — 低速で長く燃え垂れる
    // 蜂の類は推進薬を推力と回転に使い切るので、燃焼は短い(1.5〜2秒)。光る時間も短い
    if (sh.type == T_HACHI)  { tv = 0.50; tb = 1.45; tshed = 0.85; tthrust = 1.0; }  // 蜂
    // 渦蜂は「まっすぐ飛ぶ 1〜2秒」＋「回る 1.5〜2秒」なので、燃焼はその合計ぶん要る
    if (sh.type == T_UZU)    { tv = 0.90; tb = 0.80; tshed = 3.60; tthrust = 2.0; }  // 渦蜂
    if (is_shaped(sh.type))  { tv = 1.00; tb = 0.68; tshed = 0.24; } // 型物 — 形が読めるよう尾は短く
    // 錦冠(金冠) — 柳をうんと強くしたもの。低速・超長燃焼で地面近くまで垂れる
    if (sh.type == T_KAMURO) { tv = 0.40; tb = 2.30; tshed = 2.20; }
    int layers = 1 + ((sh.type == T_YANAGI || sh.type == T_HACHI || sh.type == T_UZU
                       || sh.type == T_KAMURO || is_shaped(sh.type)) ? 0 : sh.layers);
    bool crackle = (rnd() < 0.20);          // 消え際に分砲(パチパチ)する玉か

    // 型物を並べる平面の向き。玉は飛びながら回っているので、開く向きは制御しきれない。
    // 観客向き(e1=X, e2=Y)を基準に、玉ごとにランダムな軸まわりにランダムな角度だけ傾ける。
    // (真横を向くと線にしか見えないので、傾きは最大60度までにとどめてある)
    double e1x = 1, e1y = 0, e1z = 0, e2x = 0, e2y = 1, e2z = 0, enx = 0, eny = 0, enz = 1;
    if (is_shaped(sh.type)) {
        double roll = rnd() * 6.2831853;                    // 面内の回転(絵の傾き)
        double cr = cos(roll), sr = sin(roll);
        e1x = cr;  e1y = sr; e1z = 0;
        e2x = -sr; e2y = cr; e2z = 0;
        // 傾ける軸(面内のどこか)と傾き角
        double al = rnd() * 6.2831853, ph = (10.0 + 38.0 * rnd()) * M_PI / 180.0;
        double kx = cos(al), ky = sin(al), kz = 0.0;        // 回転軸(単位ベクトル)
        double cp = cos(ph), sp2 = sin(ph);
        // ロドリゲスの回転公式で 3 本まとめて回す
        auto rot = [&](double& vx, double& vy, double& vz) {
            double dotk = kx * vx + ky * vy + kz * vz;
            double crx = ky * vz - kz * vy, cry = kz * vx - kx * vz, crz = kx * vy - ky * vx;
            double nx = vx * cp + crx * sp2 + kx * dotk * (1 - cp);
            double ny = vy * cp + cry * sp2 + ky * dotk * (1 - cp);
            double nz = vz * cp + crz * sp2 + kz * dotk * (1 - cp);
            vx = nx; vy = ny; vz = nz;
        };
        rot(e1x, e1y, e1z); rot(e2x, e2y, e2z);
        enx = e1y * e2z - e1z * e2y;                        // 面の法線
        eny = e1z * e2x - e1x * e2z;
        enz = e1x * e2y - e1y * e2x;
    }
    // 玉の中の星は同じ組成なので、回る速さ・落ち方・点火の時期はだいたい揃う。
    // ただし玉に詰められたときの向きはバラバラなので、回転軸と回る向きは星ごとにランダム。
    double gomega = 20.0 + 12.0 * rnd();    // [rad/s] 点火直後の自転(3〜5回転/秒)
    // 回転の落ち方。R = 推力/ω^2 なので、強く落とすと後半で輪が膨らんでしまう。
    // ほぼ「星がその場で回っている」ように見せたいので、ω は高いまま保つ
    double gdamp  = 0.020 + 0.030 * rnd();
    double gdelay = 1.0 + 0.9 * rnd();      // [s] 割れてからノズルに火が回るまで
    for (int L = 0; L < layers; ++L) {
        double vscale = (L == 0) ? 1.0 : (0.62 - 0.20 * (L - 1));   // 芯は内側 = 遅い
        // 蜂は1個1個の軌跡を見せたいので星数を減らす(実物も星数は少ない)
        double tn = (sh.type == T_HACHI) ? 0.38 : (sh.type == T_UZU ? 0.42 : 1.0);
        int n = (int)(s.nstar * p_density * tn * (L == 0 ? 1.0 : 0.45 + 0.1 * L));
        if (n < 8) n = 8;
        if (stars.size() + n > MAX_STARS) n = (int)(MAX_STARS - std::min(MAX_STARS, stars.size()));
        if (n <= 0) break;

        Col c0, c1, tail; float chg;
        scheme_colors(sh.scheme, L, c0, c1, chg, tail);

        // 玉ごとにランダムに回した フィボナッチ球 = 星の配置
        double a1 = rnd() * 6.2831853, a2 = rnd() * 6.2831853;
        double ca1 = cos(a1), sa1 = sin(a1), ca2 = cos(a2), sa2 = sin(a2);
        double rstar = s.starR * (L == 0 ? 1.0 : 0.82);
        double burnT = s.burnT * (L == 0 ? 1.0 : 0.8) * tb;
        if (sh.type == T_YANAGI) { c0 = C_GOLD; c1 = C_AMBER; chg = 0.35f; tail = C_GOLD; }
        if (sh.type == T_HACHI)  { c0 = C_SILVER; c1 = C_GOLD; chg = 0.4f; tail = C_GOLD; }
        if (sh.type == T_UZU)    { c0 = C_SILVER; c1 = C_GOLD; chg = 0.3f; tail = C_GOLD; }
        // 菊先紅 — 燃焼の最後の1割だけ紅に変わる。錦冠の見せ場
        if (sh.type == T_KAMURO) { c0 = C_GOLD; c1 = C_RED; chg = 0.90f; tail = C_GOLD; }
        if (is_shaped(sh.type))  { c1 = c0; chg = 1.0f; tail = mixc(c0, C_GOLD, 0.25f); }
        // テーマで炎色を固定する場合(例: 桜 = ピンクの菊だけ)
        if (sh.fixc >= 0) {
            c0 = c1 = FIXCOL[sh.fixc]; chg = 1.0f;
            tail = mixc(c0, C_GOLD, 0.30f);
        }

        for (int i = 0; i < n; ++i) {
            double dx, dy, dz, v;
            int ncore = (sh.type == T_SATURN) ? (n * 38) / 100 : 0;   // 土星の芯の星数
            if (is_shaped(sh.type) && i >= ncore) {
                // 型物 — 玉の中の平面に形どおり並べ、中心から距離に比例した速度で押し出す。
                // 形はそのまま相似に拡大していく(外側ほど速いので抗力で少しだけ崩れる = 実物と同じ)。
                // 平面の向きは玉ごとにランダムに傾けてある(下の e1/e2/en)
                double px2, py2;
                shape_point(sh.type, (double)(i - ncore) / (n - ncore), px2, py2);
                double rad = sqrt(px2 * px2 + py2 * py2);
                if (rad < 1e-6) rad = 1e-6;
                double a = px2 / rad, b = py2 / rad, c = rnds() * 0.06;   // 板厚のゆらぎ
                dx = e1x * a + e2x * b + enx * c;
                dy = e1y * a + e2y * b + eny * c;
                dz = e1z * a + e2z * b + enz * c;
                v = s.starV * 0.92 * rad * (1.0 + rnds() * 0.012);
            } else if (ncore > 0 && i < ncore) {
                // 土星の芯 — こちらは3Dの小さな球
                double yy = 1.0 - 2.0 * (i + 0.5) / ncore;
                double rr = sqrt(std::max(0.0, 1.0 - yy * yy));
                double ph = i * 2.39996322973;
                dx = rr * cos(ph); dy = yy; dz = rr * sin(ph);
                v = s.starV * 0.26 * (1.0 + rnds() * 0.05);
            } else {
                double yy = 1.0 - 2.0 * (i + 0.5) / n;
                double rr = sqrt(std::max(0.0, 1.0 - yy * yy));
                double ph = i * 2.39996322973;
                dx = rr * cos(ph); dy = yy; dz = rr * sin(ph);
                // 微小な向きの揺らぎ(玉込めの誤差)
                dx += rnds() * 0.02; dy += rnds() * 0.02; dz += rnds() * 0.02;
                double t;
                t = dy * ca1 - dz * sa1; dz = dy * sa1 + dz * ca1; dy = t;   // X 回転
                t = dx * ca2 - dz * sa2; dz = dx * sa2 + dz * ca2; dx = t;   // Y 回転
                double dl = sqrt(dx * dx + dy * dy + dz * dz); dx /= dl; dy /= dl; dz /= dl;
                v = s.starV * vscale * tv * (1.0 + rnds() * 0.035);
            }

            Star st;
            st.x = sh.x; st.y = sh.y; st.z = sh.z;
            st.bx = (float)sh.x; st.by = (float)sh.y; st.bz = (float)sh.z;
            st.vx = sh.vx + dx * v; st.vy = sh.vy + dy * v; st.vz = sh.vz + dz * v;
            st.r0 = st.r = rstar * (1.0 + rnds() * 0.04);
            st.dr = st.r0 / (burnT * (1.0 + rnds() * 0.06));
            st.shed = rnd();
            if (ncore > 0 && i < ncore) {        // 土星の芯だけ銀色にして輪と区別する
                st.c0 = C_SILVER; st.c1 = C_SILVER; st.chg = 1.0f;
                st.tail = mixc(C_SILVER, C_GOLD, 0.4f);
            } else {
                st.c0 = c0; st.c1 = c1; st.chg = chg; st.tail = tail;
            }
            st.shedRate = (float)(115.0 * p_shed * tshed);
            st.delay = (L == 0) ? 0.0f : (float)(0.10 + 0.16 * L);   // 芯は少し遅れて点火(時差開発)
            st.crackle = crackle;
            st.thrust = 0.0f; st.spin = 0.0f; st.spinB = 0.0f; st.sphase = 0.0f;
            st.spinDelay = 0.0f; st.exhaust = 0.0f;
            st.ux = st.uy = st.uz = st.wx = st.wy = st.wz = 0.0f;

            if (tthrust > 0.0) {
                // 蜂 — 星が自分で推力を出し、その向きが軸まわりに回る。
                // 推力方向 = u cos(θ) + w sin(θ),  θ += spin*dt
                // 回転が速いと推力が打ち消し合ってブルブル、遅いと大きな輪を描く。
                // 渦の半径はだいたい (推力)/(角速度)^2 になる
                bool uzu = (sh.type == T_UZU);
                // 回転軸は星ごとに 3D でランダム(玉に詰められたときの向きがバラバラだから)
                double axx = rnds(), axy = rnds(), axz = rnds();
                double al = sqrt(axx * axx + axy * axy + axz * axz);
                if (al < 1e-6) { axx = 0; axy = 1; axz = 0; al = 1; }
                axx /= al; axy /= al; axz /= al;
                double tx = (fabs(axx) < 0.9) ? 1.0 : 0.0, ty = (fabs(axx) < 0.9) ? 0.0 : 1.0;
                double ux = axy * 0.0 - axz * ty, uy = axz * tx - axx * 0.0, uz = axx * ty - axy * tx;
                double ul = sqrt(ux * ux + uy * uy + uz * uz);
                ux /= ul; uy /= ul; uz /= ul;
                double wx = axy * uz - axz * uy, wy = axz * ux - axx * uz, wz = axx * uy - axy * ux;
                st.ux = (float)ux; st.uy = (float)uy; st.uz = (float)uz;
                st.wx = (float)wx; st.wy = (float)wy; st.wz = (float)wz;
                if (uzu) {
                    // 渦蜂 — 自転が速いぶん軌跡の輪 a/ω^2 は小さく、ほぼ直進しながら
                    // 噴射だけがぐるぐる回る。割れてしばらくしてからノズルに火が回る
                    st.thrust  = (float)(900.0 + 600.0 * rnd());
                    // 速さと落ち方は組成が同じなのでほぼ揃う。回る向き(左右)は星ごとにランダム
                    st.spin    = (float)((rnd() < 0.5 ? -1.0 : 1.0) * gomega * (1.0 + rnds() * 0.14));
                    st.spinB   = (float)(gdamp * (1.0 + rnds() * 0.15));
                    st.spinDelay = (float)(gdelay * (1.0 + rnds() * 0.10));
                    st.exhaust = (float)(26.0 + 16.0 * rnd());   // 噴射速度 [m/s]
                } else {
                    // 蜂 — 自転が遅いので推力が軌跡そのものを大きく振り回す
                    st.thrust  = (float)(45.0 + 60.0 * rnd());
                    st.spin    = (float)((0.7 + 5.5 * rnd()) * (rnd() < 0.5 ? -1.0 : 1.0));
                    st.spinB   = 0.0f;
                    st.spinDelay = 0.0f;
                    st.exhaust = (float)(6.0 + 6.0 * rnd());
                }
                st.sphase = (float)(rnd() * 6.2831853);
            }
            st.alive = true;
            stars.push_back(st);
        }
    }
    if (flashes.size() < 64) {
        Flash f; f.x = sh.x; f.y = sh.y; f.z = sh.z;
        f.t0 = f.t = 0.10; f.pw = 3.2 * s.go; f.alive = true;
        flashes.push_back(f);
    }
    sh.alive = false;
}

// 千輪の子玉 — 親の開発点から飛び出し、短い秒時で小さく開く
static void spawn_child(const Shell& p, double dx, double dy, double dz, double v,
                        double go, double fuse, int type) {
    if (shells.size() >= 400) return;
    (void)0;
    ShellSpec s = spec_of(go);
    Shell c;
    c.x = p.x; c.y = p.y; c.z = p.z;
    c.vx = p.vx + dx * v; c.vy = p.vy + dy * v; c.vz = p.vz + dz * v;
    c.k = 0.5 * RHO_AIR * CD * (M_PI * s.dia * s.dia * 0.25) / s.mass;
    c.shed = 0;
    c.fuse = fuse;
    c.sp = s;
    c.type = type;
    c.layers = 0;
    c.scheme = (rnd() < 0.5) ? 0 : 1;
    c.fixc = p.fixc;                 // 千輪の子玉はテーマの色を引き継ぐ
    c.gen = p.gen + 1;
    c.alive = true;
    shells.push_back(c);
}

// nx: 画面内の水平位置 [-1,1] (-1000 でランダム) / go: 号数 / type: 玉の種類(T_*)
// fixc: 炎色の固定(FIXCOL の添字。-1 で自由)
static void launch_ex(double nx, double go, int type, double fuseJit = 0.03, int fixc = -1,
                      double fuseScale = 1.0) {
    if (shells.size() >= 90) return;
    if (go < 3.0) go = 3.0;
    ShellSpec s = spec_of(go);
    Shell sh;
    double halfW = (FW * 0.5) * camD / camF;
    sh.x = (nx > -900) ? (nx * halfW * 0.80) : rnds() * halfW * 0.70;
    sh.z = rnds() * camD * 0.16;
    sh.y = 0.0;
    sh.k = 0.5 * RHO_AIR * CD * (M_PI * s.dia * s.dia * 0.25) / s.mass;
    // 打揚薬が与える初速。ここから先は積分するだけ — 到達高度は結果として出る
    sh.vx = rnds() * 1.2; sh.vy = s.liftV * (1.0 + rnds() * 0.02); sh.vz = rnds() * 1.2;
    sh.shed = 0;
    // 時限導火線に点火。秒時のばらつき=開発高度のばらつき。
    // fuseScale を縮めると頂点まで待たずに低い位置で開く(フェニックスの低い小玉)
    sh.fuse = s.fuseT * fuseScale * (1.0 + rnds() * fuseJit);
    sh.sp = s;
    sh.type = type;
    sh.gen = 0;
    sh.fixc = fixc;
    sh.layers = (int)(p_layers + 0.5);
    double q = rnd();
    sh.scheme = (sh.layers > 0 && q < 0.30) ? 3 : (q < 0.55 ? 0 : (q < 0.85 ? 1 : 2));
    if (is_shaped(type)) sh.scheme = 1;              // 型物は形が読めるよう単色に
    sh.alive = true;
    shells.push_back(sh);
}
// p_type が -1(おまかせ)なら重み付きでランダムに選ぶ
static int pick_type() {
    if (p_type >= 0) return (int)(p_type + 0.5);
    double q = rnd();
    if (q < 0.30) return T_KIKU;
    if (q < 0.42) return T_BOTAN;
    if (q < 0.53) return T_YANAGI;
    if (q < 0.64) return T_HACHI;
    if (q < 0.74) return T_UZU;
    if (q < 0.85) return T_SENRIN;
    if (q < 0.90) return T_HEART;
    if (q < 0.95) return T_RING;
    return T_SATURN;
}
// スターマインのテーマ。1テーマ = 1種類に揃えるので、菊とハートが混ざって濁らない。
// 添字 0 は「おまかせ」= 打ち上げ開始時に 1 以降からランダムに1つ選ぶ(全種類を混ぜはしない)。
// type: -2=型物からランダム / それ以外は T_* をそのまま
// fixc: FIXCOL の添字で炎色を固定(-1 で自由)
// layout: 0=通常(掃引しながら1〜3発) / 1=フェニックス(大玉を横一列＋低い小玉を同時開花)
struct Theme { int type; int fixc; int layout; };
static const Theme THEMES[] = {
    { T_KIKU,   -1, 0 },  // 0 おまかせ(実際には下から1つ選ばれる。ここは保険)
    { T_KIKU,   -1, 0 },  // 1 菊 — 色とりどり
    { T_KIKU,    0, 0 },  // 2 桜 — ピンクの菊だけ
    { T_KIKU,    5, 0 },  // 3 銀世界 — 銀の菊だけ
    { T_BOTAN,  -1, 0 },  // 4 牡丹
    { T_YANAGI, -1, 0 },  // 5 金の柳
    { T_HACHI,  -1, 0 },  // 6 蜂
    { T_UZU,    -1, 0 },  // 7 渦蜂
    { T_SENRIN, -1, 0 },  // 8 千輪
    { -2,       -1, 0 },  // 9 型物
    { T_KAMURO, -1, 0 },  // 10 錦冠(金冠)
    { T_KAMURO,  5, 1 },  // 11 フェニックス — 銀冠の大玉を横一列＋低い小玉を同時に
};
static const int THEME_COUNT = (int)(sizeof(THEMES) / sizeof(THEMES[0]));

static int theme_layout() {
    int ti = curTheme;
    if (ti < 1) ti = 1; if (ti >= THEME_COUNT) ti = THEME_COUNT - 1;
    return THEMES[ti].layout;
}
static void theme_pick(int& ty, int& fx) {
    int ti = curTheme;
    if (ti < 1) ti = 1; if (ti >= THEME_COUNT) ti = THEME_COUNT - 1;
    const Theme& th = THEMES[ti];
    fx = th.fixc;
    if (th.type == -2) { const int s[3] = { T_HEART, T_RING, T_SATURN }; ty = s[(int)(rnd() * 3)]; }
    else               ty = th.type;
}
static void launch(double nx) { launch_ex(nx, p_go, pick_type()); }

// スターマイン(早打ち) — 発数だけ装填して、あとは連打間隔を差分で消化していく
static void schedule_launch(double delay, double nx, double go, int type, int fixc,
                            double fuseScale, double fuseJit) {
    if (pending.size() >= 600) return;
    Pending p;
    p.t = delay; p.nx = nx; p.go = go; p.type = type; p.fixc = fixc;
    p.fuseScale = fuseScale; p.fuseJit = fuseJit;
    pending.push_back(p);
}

// フェニックス(長岡) — 大玉を横一列に、その手前の低い位置に小玉を大量に。
// 号数が違えば秒時も違うので、開発までの時間差だけ小玉を遅らせて撃ち、空で同時に開かせる。
static void phoenix_wave(int& budget) {
    int ty, fx; theme_pick(ty, fx);

    // --- 大玉: 横一列に等間隔。開くタイミングを少しだけずらして「完全な同時」を避ける
    int nbig = 4 + (int)(rnd() * 3.0);                       // 4〜6発
    double bigGo = p_go;
    double bigFuse = spec_of(bigGo).fuseT;
    for (int i = 0; i < nbig && budget > 0; ++i, --budget) {
        double x = -0.98 + 1.96 * (i + 0.5) / nbig + rnds() * 0.04;
        // 号数と秒時を少しずつ変えて、高さと大きさを揃えすぎない
        double g = (rnd() < 0.32) ? std::max(3.0, bigGo - 2.0) : bigGo;
        schedule_launch(rnd() * 0.28, x, g, ty, fx, 0.90 + rnd() * 0.18, 0.015);
    }

    // --- 小玉: うんと低い位置で、大玉と同時に開くように遅らせて予約
    double smallGo = std::max(3.0, floor(p_go * 0.30 + 0.5));
    double smallScale = 0.45;                                 // 頂点よりだいぶ手前で開く = 低い
    double smallFuse = spec_of(smallGo).fuseT * smallScale;
    double delay = std::max(0.0, bigFuse - smallFuse);        // これで空での開花が揃う
    int nsmall = 13 + (int)(rnd() * 7.0);                     // 13〜19発
    for (int i = 0; i < nsmall && budget > 0; ++i, --budget) {
        double x = -0.97 + 1.94 * (i + 0.5) / nsmall + rnds() * 0.05;
        int sty = (rnd() < 0.65) ? T_BOTAN : T_KIKU;          // 小玉は尾の短いものを主体に
        schedule_launch(delay + rnds() * 0.12, x, smallGo, sty, fx, smallScale, 0.02);
    }
}

static void start_barrage() {
    barTotal = barLeft = (int)(p_shots + 0.5);
    barTimer = 0.0;
    barPhase = 1;
    quietTimer = 0.0;
    starTimer = p_starInt;       // 次の自動打ち上げまでを測り直す
    // おまかせ = 全種類を混ぜるのではなく、テーマを1つ引く(1発ごとではなく連打ごとに確定)
    int ti = (int)(p_theme + 0.5);
    curTheme = (ti <= 0) ? (1 + (int)(rnd() * (THEME_COUNT - 1))) : ti;
    if (curTheme >= THEME_COUNT) curTheme = THEME_COUNT - 1;
}

// ---------------------------------------------------------------- ABI
extern "C" {

KEEP int sim_w() { return FW; }
KEEP int sim_h() { return FH; }

KEEP void sim_reset() {
    stars.clear(); sparks.clear(); shells.clear(); flashes.clear(); pending.clear();
    stars.reserve(20000); sparks.reserve(60000);
    acc.assign((size_t)FW * FH * 3, 0.0f);
    simTime = 0; launchTimer = 0.4; frameNo = 0;
    measH = 0; measD = 0;
    barLeft = barTotal = barPhase = 0; barTimer = 0; quietTimer = 0; qual = 1.0;
    starTimer = p_starInt;
    curSpec = spec_of(p_go);
    setup_camera(curSpec);
    build_bg();
}

KEEP int sim_init(int, int) {
    px.assign((size_t)FW * FH, 0);
    sim_reset();
    return 1;
}

KEEP void sim_set(int id, double v) {
    switch (id) {
    case 0: {
        double g = floor(v + 0.5);
        if (g != p_go) {
            p_go = g;
            curSpec = spec_of(p_go);
            measH = 0; measD = 0;
            setup_camera(curSpec);
            build_bg();
        }
        break;
    }
    case 1: p_density = v; break;
    case 2: p_interval = v; break;
    case 3: p_wind = v; break;
    case 4: p_shed = v; break;
    case 5: p_layers = floor(v + 0.5); break;
    case 6: p_glow = v; break;
    case 7: p_shots = v; break;
    case 8: p_rapid = v; break;
    case 9: p_quiet = v; break;
    case 10: p_type = v; break;
    case 11: p_theme = v; break;
    case 12: {                               // 自動スターマインの間隔
        double old = p_starInt;
        p_starInt = v;
        if (v > 0.5 && (old <= 0.5 || starTimer > v)) starTimer = v;
        break;
    }
    case 13: p_hud = v; break;
    }
}

// UI に状態を返す(universe_cpp の共通 ABI への追加分)
KEEP double sim_get(int id) {
    switch (id) {
    case 0: return (double)barLeft;
    case 1: return (double)barTotal;
    case 2: return (double)stars.size();
    case 3: return (double)sparks.size();
    case 4: return (double)shells.size();
    case 5: return p_auto ? 1.0 : 0.0;
    case 6: return measH;
    case 7: return measD * 2.0;
    case 8: return qual;
    case 9: return simTime;
    case 10: return (double)barPhase;
    case 11: return quietTimer;
    case 12: return p_theme;
    case 13: return (double)curTheme;    // この連打で実際に使われているテーマ
    case 14: return starTimer;           // 次の自動スターマインまで [s]
    case 15: return p_starInt;
    }
    return 0.0;
}

KEEP void sim_action(int id) {
    if (id == 0) launch(-1000.0);            // 手動で1発
    else if (id == 1) sim_reset();
    else if (id == 2) p_auto = !p_auto;
    else if (id == 3) start_barrage();       // スターマイン
    else if (id == 4) barLeft = 0;           // 連打中止
}

KEEP void sim_click(double nx, double ny) {
    (void)ny;
    launch(nx * 2.0 - 1.0);
}

KEEP void sim_step(int frames) {
    for (int f = 0; f < frames; ++f) {
        for (int ss = 0; ss < SUBSTEPS; ++ss) {
            const double dt = DT;
            simTime += dt;

            // --- 粒子予算の負荷から引先の量を追従調整(一次遅れ)
            {
                double load = (double)sparks.size() / (double)MAX_SPARKS;
                double tgt = (load < 0.55) ? 1.0 : std::max(0.10, 1.0 - (load - 0.55) / 0.45);
                qual += (tgt - qual) * 0.05;
            }

            // --- スターマイン(早打ち連打)
            // --- 予約発射(フェニックスの小玉など)。これも差分で焼いていくだけ
            for (size_t i = 0; i < pending.size();) {
                pending[i].t -= dt;
                if (pending[i].t <= 0.0) {
                    Pending p = pending[i];
                    pending[i] = pending.back(); pending.pop_back();
                    launch_ex(p.nx, p.go, p.type, p.fuseJit, p.fixc, p.fuseScale);
                } else ++i;
            }

            if (barLeft > 0 && theme_layout() == 1) {
                // フェニックス — ウェーブ単位で大玉＋小玉をまとめて
                barTimer -= dt;
                if (barTimer <= 0.0) {
                    phoenix_wave(barLeft);
                    barTimer = 2.2 + rnd() * 1.0;
                }
            } else if (barLeft > 0) {
                barTimer -= dt;
                if (barTimer <= 0.0) {
                    int volley = 1 + (int)(rnd() * 2.7);          // 1〜3発同時に上がる
                    for (int v = 0; v < volley && barLeft > 0; ++v) {
                        bool last = (barLeft == 1);
                        int ty, fx; theme_pick(ty, fx);
                        double go = (last || rnd() < 0.12)
                                  ? p_go
                                  : floor(p_go * (0.45 + 0.50 * rnd()) + 0.5);   // 小玉主体の早打ち
                        // 左右に掃引しながら上げる(実際のスターマインの並び)
                        double prog = 1.0 - (double)barLeft / (double)barTotal;
                        double sweep = sin(prog * 15.7) * 0.82 + rnds() * 0.30;   // 左右2.5往復
                        if (last) sweep = rnds() * 0.15;
                        launch_ex(sweep, go, ty, last ? 0.02 : 0.13, fx);
                        barLeft--;
                    }
                    barTimer = p_rapid * (0.55 + 0.9 * rnd());
                    if (barLeft == 1) barTimer += 0.9;            // 締めの前に一拍おく
                }
            }

            // --- 連打が終わったら余韻。最後の玉が開いて星が燃え尽きるまで待ち、
            //     空が完全に空になってから静寂の時間を差分で消化する
            if (barPhase == 1 && barLeft == 0) { barPhase = 2; quietTimer = p_quiet; }
            if (barPhase == 2) {
                if (shells.empty() && stars.empty() && pending.empty()) {
                    quietTimer -= dt;
                    if (quietTimer <= 0.0) { barPhase = 0; launchTimer = p_interval; }
                } else {
                    quietTimer = p_quiet;   // まだ空に何か残っている = 余韻はまだ始まらない
                }
            }

            // --- 自動スターマイン。単発が上がっている通常時だけ時計が進む
            //     (連打中と余韻のあいだは止まるので、間隔は「静かな時間」で測られる)
            if (barPhase == 0 && p_starInt > 0.5) {
                starTimer -= dt;
                if (starTimer <= 0.0) start_barrage();
            }

            // --- 通常の打ち上げ間隔
            if (p_auto && barPhase == 0) {
                launchTimer -= dt;
                if (launchTimer <= 0.0) {
                    launch_ex(-1000.0, p_go, pick_type());
                    launchTimer = p_interval * (0.75 + 0.5 * rnd());
                }
            }

            // --- 玉(上昇中)
            for (size_t i = 0; i < shells.size(); ++i) {
                Shell& sh = shells[i];
                if (!sh.alive) continue;
                double vrx = sh.vx - p_wind, vry = sh.vy, vrz = sh.vz;
                double sp = sqrt(vrx * vrx + vry * vry + vrz * vrz);
                sh.vx += (-sh.k * sp * vrx) * dt;
                sh.vy += (-sh.k * sp * vry - G_ACC) * dt;
                sh.vz += (-sh.k * sp * vrz) * dt;
                sh.x += sh.vx * dt; sh.y += sh.vy * dt; sh.z += sh.vz * dt;
                // 昇り曲導
                sh.shed += 150.0 * p_shed * qual * (sh.gen ? 0.30 : 1.0) * dt;
                while (sh.shed >= 1.0) {
                    sh.shed -= 1.0;
                    emit_spark(sh.x, sh.y, sh.z,
                               sh.vx * 0.05 + rnds() * 1.5, sh.vy * 0.05 + rnds() * 1.5, sh.vz * 0.05 + rnds() * 1.5,
                               mixc(C_GOLD, C_SILVER, (float)rnd() * 0.5f), 0.00035);
                }
                if (sh.y > measH) measH = sh.y;   // 到達高度は積分の結果
                sh.fuse -= dt;                    // 時限導火線が焼ける
                // 親玉は秒時が狂っても落下し始めたら開く。子玉(千輪)は下向きに撒かれるので秒時のみ
                if (sh.fuse <= 0.0 || (sh.gen == 0 && sh.vy < -12.0)) burst(sh);
            }
            for (size_t i = 0; i < shells.size();) {
                if (!shells[i].alive) { shells[i] = shells.back(); shells.pop_back(); }
                else ++i;
            }

            // --- 星
            for (size_t i = 0; i < stars.size(); ++i) {
                Star& st = stars[i];
                double k = star_k(st.r);
                double vrx = st.vx - p_wind, vry = st.vy, vrz = st.vz;
                double sp = sqrt(vrx * vrx + vry * vry + vrz * vrz);
                double ax = -k * sp * vrx, ay = -k * sp * vry - G_ACC, az = -k * sp * vrz;
                double sph0 = st.sphase;      // このサブステップの開始位相(火の粉の補間に使う)
                if (st.thrust > 0.0f && st.spinDelay > 0.0f) {
                    // まだノズルに火が回っていない。割薬で飛び出したまま普通の星として飛ぶ
                    st.spinDelay -= (float)dt;
                } else if (st.thrust > 0.0f) {
                    // 機体の自転。空気抵抗で落ちていく: dω/dt = -B|ω|ω (これも差分)
                    st.spin -= st.spinB * fabsf(st.spin) * st.spin * (float)dt;
                    st.sphase += st.spin * (float)dt;
                    // 噴射の反作用。ノズルの向きに加速する。軌跡はこの積分の結果として出る
                    double nx, ny, nz;
                    nozzle_dir(st, st.sphase, nx, ny, nz);
                    double th = st.thrust * (st.r / st.r0);   // 燃え尽きるにつれ推力も落ちる
                    ax += th * nx; ay += th * ny; az += th * nz;
                }
                st.vx += ax * dt; st.vy += ay * dt; st.vz += az * dt;
                st.x += st.vx * dt; st.y += st.vy * dt; st.z += st.vz * dt;
                if (st.delay > 0.0f) { st.delay -= (float)dt; continue; }   // 暗飛行中は燃えない
                st.r -= st.dr * dt;
                if (st.r <= 1e-4 || st.y < 0.0) {
                    // 消え際の分砲(パチパチ)
                    if (st.crackle && st.y > 0.0)
                        for (int c = 0; c < 10; ++c)
                            emit_spark(st.x, st.y, st.z,
                                       st.vx + rnds() * 9.0, st.vy + rnds() * 9.0, st.vz + rnds() * 9.0,
                                       C_SILVER, 0.00022);
                    st.alive = false; continue;
                }
                // 開花直径も積分の結果 — 開発点からの水平距離の running max
                double ddx = st.x - st.bx, ddz = st.z - st.bz;
                double d2 = ddx * ddx + ddz * ddz;
                if (d2 > measD * measD) measD = sqrt(d2);
                // 引先(尾) — 燃えかすを実際に放出する
                double prog = 1.0 - st.r / st.r0;
                Col hc = (prog < st.chg) ? st.c0 : st.c1;
                (void)hc;
                st.shed += st.shedRate * qual * dt;
                if (st.shed >= 1.0) {
                    int nsh = (int)st.shed;
                    bool jet = (st.thrust > 0.0f && st.spinDelay <= 0.0f);
                    for (int q = 0; q < nsh; ++q) {
                        st.shed -= 1.0;
                        if (jet) {
                            // 噴射。ノズルの反対向きに吹き出す(これが彗星の尾になる)。
                            // 回転が速いと 1 サブステップで位相が進むので、
                            // このサブステップ内の円弧上に位相を補間して撒く
                            double f = (q + 0.5) / nsh;
                            double nx, ny, nz;
                            nozzle_dir(st, sph0 + (st.sphase - sph0) * f, nx, ny, nz);
                            double ve = st.exhaust * (0.75 + 0.5 * rnd());
                            emit_spark(st.x, st.y, st.z,
                                       st.vx * 0.10 - ve * nx + rnds() * 2.0,
                                       st.vy * 0.10 - ve * ny + rnds() * 2.0,
                                       st.vz * 0.10 - ve * nz + rnds() * 2.0,
                                       st.tail, 0.00032);
                        } else {
                            emit_spark(st.x, st.y, st.z,
                                       st.vx * 0.12 + rnds() * 2.0,
                                       st.vy * 0.12 + rnds() * 2.0,
                                       st.vz * 0.12 + rnds() * 2.0,
                                       st.tail, 0.00032);
                        }
                    }
                }
            }
            for (size_t i = 0; i < stars.size();) {
                if (!stars[i].alive) { stars[i] = stars.back(); stars.pop_back(); }
                else ++i;
            }

            // --- 火の粉
            for (size_t i = 0; i < sparks.size(); ++i) {
                Spark& s = sparks[i];
                float vrx = s.vx - (float)p_wind, vry = s.vy, vrz = s.vz;
                float sp = sqrtf(vrx * vrx + vry * vry + vrz * vrz);
                float kk = s.k * sp * (float)dt;
                s.vx -= kk * vrx; s.vy -= kk * vry + (float)(G_ACC * dt); s.vz -= kk * vrz;
                s.x += s.vx * (float)dt; s.y += s.vy * (float)dt; s.z += s.vz * (float)dt;
                s.life -= (float)dt;
                if (s.life <= 0.0f || s.y < 0.0f) s.alive = false;
            }
            for (size_t i = 0; i < sparks.size();) {
                if (!sparks[i].alive) { sparks[i] = sparks.back(); sparks.pop_back(); }
                else ++i;
            }

            // --- 開発の閃光
            for (size_t i = 0; i < flashes.size();) {
                flashes[i].t -= dt;
                if (flashes[i].t <= 0) { flashes[i] = flashes.back(); flashes.pop_back(); }
                else ++i;
            }
        }
        frameNo++;
    }
}

KEEP uint8_t* sim_render() {
    // 残光の減衰
    float dec = (float)p_glow;
    size_t n3 = acc.size();
    for (size_t i = 0; i < n3; ++i) acc[i] *= dec;

    // 火の粉 — 燃え尽きるにつれて暗く、赤く冷える
    for (size_t i = 0; i < sparks.size(); ++i) {
        const Spark& s = sparks[i];
        double sx, sy, sc;
        if (!project(s.x, s.y, s.z, sx, sy, sc)) continue;
        float lf = s.life / s.life0;
        Col c = mixc({ 0.95f, 0.16f, 0.03f }, s.c, lf);
        float amp = 0.34f * lf * lf * (float)(sc / camF * 520.0);
        splat_point(sx, sy, c, amp);
    }

    // 星頭 — 輝度は燃焼表面積 (r/r0)^2 に比例
    for (size_t i = 0; i < stars.size(); ++i) {
        const Star& st = stars[i];
        if (st.delay > 0.0f) continue;          // 暗飛行中は光らない
        double sx, sy, sc;
        if (!project(st.x, st.y, st.z, sx, sy, sc)) continue;
        double prog = 1.0 - st.r / st.r0;
        Col c = (prog < st.chg) ? st.c0 : st.c1;
        float f = (float)((st.r / st.r0) * (st.r / st.r0));
        float amp = (0.9f + 2.6f * f) * (float)(sc / camF * 520.0);
        splat_disc(sx, sy, 1.0f + 1.4f * f, c, amp);
    }

    // 開発の閃光
    for (size_t i = 0; i < flashes.size(); ++i) {
        const Flash& fl = flashes[i];
        double sx, sy, sc;
        if (!project(fl.x, fl.y, fl.z, sx, sy, sc)) continue;
        float t = (float)(fl.t / fl.t0);
        splat_disc(sx, sy, (float)(fl.pw * t), { 1.0f, 0.92f, 0.72f }, 2.2f * t * t);
    }

    // トーンマップして夜空に合成
    //
    // **既定は成分ごとの r/(1+r)。** 明るいところから順に色が抜けるので、明るい芯は
    // クリーム色〜白へ寄り、色は暗い先端に残る。フィルムやセンサが飽和したときの
    // 写り方に近く、花火ではこれが自然に見える。
    //
    // `-DLIGHT_VIVID` で、明るさ(最大成分)だけ圧縮して色の向きを保つ式に切り替わる。
    // 色は最後まで残るが、**花火では鮮やかすぎて嘘っぽく見える**(2026-08-14 にユーザ判断で
    // 既定から外した)。点光源の粒がはっきり分かれている drone_emu_cpp では、こちらのほうが
    // 実物に近かった。同じ式でも題材で結論が変わるので、両方ビルドして
    // compare.html で見比べられるようにしてある。
    for (size_t i = 0, np = (size_t)FW * FH; i < np; ++i) {
        float r = acc[i * 3], g = acc[i * 3 + 1], b = acc[i * 3 + 2];
        if (r + g + b < 0.004f) { px[i] = bg[i]; continue; }
#ifdef LIGHT_VIVID
        float mx = r > g ? (r > b ? r : b) : (g > b ? g : b);
        float t = mx / (1.0f + mx);           // 圧縮した明るさ
        float s = t / mx;                     // 色の向きを保ったまま合わせる倍率
        float wht = smoothstep(4.0f, 16.0f, mx);   // ここから上は飛んで白くなる
        int R = (int)(255.0f * (r * s + (t - r * s) * wht));
        int G = (int)(255.0f * (g * s + (t - g * s) * wht));
        int B = (int)(255.0f * (b * s + (t - b * s) * wht));
#else
        int R = (int)(255.0f * (r / (1.0f + r)));
        int G = (int)(255.0f * (g / (1.0f + g)));
        int B = (int)(255.0f * (b / (1.0f + b)));
#endif
        uint32_t d = bg[i];
        R += (int)(d & 255); G += (int)((d >> 8) & 255); B += (int)((d >> 16) & 255);
        px[i] = rgb(R > 255 ? 255 : R, G > 255 ? 255 : G, B > 255 ? 255 : B);
    }

    if (p_hud < 0.5) return (uint8_t*)px.data();     // 文字なし(フルスクリーン鑑賞モード)

    Olivec_Canvas oc = olivec_canvas(px.data(), FW, FH, FW);
    char buf[160];
    snprintf(buf, sizeof(buf), "KIKU %d-GO %.0fmm  lift %.0f  burst %.0f m/s -> H %.0fm  flower %.0fm",
             (int)p_go, curSpec.dia * 1000.0, curSpec.liftV, curSpec.starV, measH, measD * 2.0);
    olivec_text(oc, buf, 14, 12, olivec_default_font, 2, rgb(255, 190, 90));
    if (barPhase == 1)
        snprintf(buf, sizeof(buf), "STARMINE %d/%d   stars %d  sparks %d  shells %d  q%.2f",
                 barTotal - barLeft, barTotal, (int)stars.size(), (int)sparks.size(),
                 (int)shells.size(), qual);
    else if (barPhase == 2)
        snprintf(buf, sizeof(buf), "... afterglow %.1fs   stars %d  sparks %d  shells %d",
                 quietTimer, (int)stars.size(), (int)sparks.size(), (int)shells.size());
    else
        snprintf(buf, sizeof(buf), "stars %5d   sparks %6d   shells %d   t=%.1fs",
                 (int)stars.size(), (int)sparks.size(), (int)shells.size(), simTime);
    olivec_text(oc, buf, 14, FH - 26, olivec_default_font, 2, rgb(150, 160, 190));
    return (uint8_t*)px.data();
}

}  // extern "C"

// ---------------------------------------------------------------- native self-test
#ifndef __EMSCRIPTEN__
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <cstdio>
#include <cstdlib>
int main(int argc, char** argv) {
    sim_init(0, 0);
    int steps = argc > 1 ? atoi(argv[1]) : 480;
    const char* out = argc > 2 ? argv[2] : "hanabi_preview.png";
    int barrage = argc > 3 ? atoi(argv[3]) : 0;   // 4番目の引数 = スターマインの発数
    if (argc > 4) sim_set(10, atof(argv[4]));     // 5番目の引数 = 玉の種類(-1でおまかせ)
    if (argc > 5) sim_set(11, atof(argv[5]));     // 6番目の引数 = スターマインのテーマ
    if (argc > 6) sim_set(12, atof(argv[6]));     // 7番目の引数 = 自動スターマインの間隔
    if (argc > 7) sim_set(13, atof(argv[7]));     // 8番目の引数 = 0 で HUD の文字を消す
    if (argc > 8) sim_set(0,  atof(argv[8]));     // 9番目の引数 = 号数
    if (barrage) { sim_set(7, barrage); sim_action(3); }
    for (int i = 0; i < steps; ++i) {
        sim_step(1); sim_render();                                   // 毎フレーム描画 = 実負荷
        if ((barrage || argc > 6) && i % 30 == 0)
            printf("  t=%5.1fs phase=%d theme=%d bar=%3d/%-3d quiet=%4.1f next=%5.1f "
                   "shells=%2d stars=%5d sparks=%6d q=%.2f\n",
                   sim_get(9), (int)sim_get(10), (int)sim_get(13), barTotal - barLeft, barTotal,
                   sim_get(11), sim_get(14), (int)sim_get(4), (int)sim_get(2), (int)sim_get(3),
                   sim_get(8));
    }
    uint8_t* p = sim_render();
    int nb = 0; for (int k = 0; k < FW * FH; ++k) if (px[k] != bg[k]) nb++;
    printf("hanabi_os native: %dx%d  frames=%d  stars=%d sparks=%d  lit=%d\n",
           FW, FH, steps, (int)stars.size(), (int)sparks.size(), nb);
    printf("  spec : %d-go dia=%.0fmm mass=%.2fkg lift=%.0fm/s fuse=%.1fs starV=%.0fm/s starR=%.1fmm burn=%.1fs n=%d\n",
           (int)p_go, curSpec.dia * 1000, curSpec.mass, curSpec.liftV, curSpec.fuseT,
           curSpec.starV, curSpec.starR * 1000, curSpec.burnT, curSpec.nstar);
    printf("  meas : apogee=%.0fm  flower dia=%.0fm   (both emerge from the integration)\n", measH, measD * 2.0);
    stbi_write_png(out, FW, FH, 4, p, FW * 4);
    return 0;
}
#endif
