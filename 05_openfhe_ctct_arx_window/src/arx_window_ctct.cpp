// 05 — OpenFHE CKKS CT-CT: delayed-window ARX control with bootstrapped inputs.
//
// u(t) is computed every step from input/output sequences only (no encrypted state):
//
//   block k covers t = kK .. kK+K-1, window tau_k = kK - d, m = t - tau_k - 1
//   u(t) = sum_j alpha[m][j] u_bt(tau_k - j) + sum_j beta[m][j] y(tau_k - j)
//          + sum_{l<m} gamma[l] y(t-1-l)
//   alpha = H F^m Fb^j R,  beta = H F^m Fb^j G,  gamma = H F^l G,  Fb = F - R H (Fb^n = 0)
//
// The same window is reused for K steps, so each bootstrapped window input has d sampling
// periods to finish.  Only window inputs are bootstrapped (min(n,K)/K per step).  Before the
// first window exists, u(t) = H F^t xc0 + sum_l gamma[l] y(t-1-l) with an encrypted xc0.
// All coefficients are encrypted (CT-CT); every u(t) has multiplicative depth one.
// See tex/formulation.tex for the derivation.
#include "openfhe.h"
#include "../../include/controller_profiles.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace lbcrypto;
using namespace tutorial_profiles;
using Cipher = Ciphertext<DCRTPoly>;
using Clock = std::chrono::steady_clock;

namespace {
constexpr std::size_t kN = 4;
constexpr uint32_t kSlots = 4;
constexpr uint32_t kCorrectionFactor = 13;

struct Options {
    std::string profile = "stable";
    std::size_t K = 4;
    std::size_t d = 0;  // 0 -> d = K
    std::size_t steps = 200;
    std::string resultsDir = "results";
    bool keepImag = false;
    uint32_t bootIters = 2;
    uint32_t levelsAfterBootstrap = 2;
    bool packY = false;
    bool packU = false;
    bool secure128 = false;      // HEStd_128_classic (OpenFHE picks N = 65536 for this depth)
    bool compactInputs = false;  // encrypt coefficients / y / xc0 with only the levels they need
};

Options ParseOptions(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) throw std::invalid_argument("missing value for " + a);
            return argv[++i];
        };
        if (a == "--profile") o.profile = next();
        else if (a == "--K") o.K = std::stoul(next());
        else if (a == "--d") o.d = std::stoul(next());
        else if (a == "--steps") o.steps = std::stoul(next());
        else if (a == "--results-dir") o.resultsDir = next();
        else if (a == "--keep-imag") o.keepImag = true;
        else if (a == "--boot-iters") o.bootIters = std::stoul(next());
        else if (a == "--levels-after-bootstrap") o.levelsAfterBootstrap = std::stoul(next());
        else if (a == "--pack-y") o.packY = true;
        else if (a == "--pack-u") o.packU = true;
        else if (a == "--secure128") o.secure128 = true;
        else if (a == "--compact-inputs") o.compactInputs = true;
        else
            throw std::invalid_argument(
                "usage: arx_window_ctct [--profile stable|unstable_low|unstable_high] [--K N] [--d N] "
                "[--steps N] [--results-dir PATH] [--keep-imag] [--boot-iters 1|2] "
                "[--levels-after-bootstrap N>=2] [--pack-y] [--pack-u] [--secure128] [--compact-inputs]");
    }
    // Default delay d = K - n + 1: the window then lies inside one propagation segment, so the
    // block-reset jump is not fed into the next window (see tex/formulation.tex, sec:dstab).
    // d > K - n + 1 makes the jump recursion J_k = -F^{d-1} Fb^{K-d+1} J_{k-1}, which diverges
    // for large K when rho(F) > 1.
    if (o.d == 0) o.d = o.K + 1 > kN ? o.K + 1 - kN : 1;
    if (o.K == 0 || o.steps == 0) throw std::invalid_argument("--K and --steps must be positive");
    if (o.bootIters < 1 || o.bootIters > 2) throw std::invalid_argument("--boot-iters must be 1 or 2");
    if (o.levelsAfterBootstrap < 2) throw std::invalid_argument("--levels-after-bootstrap must be >= 2");
    return o;
}

// ---- plain linear algebra -------------------------------------------------------------
double Dot(const Vec4& a, const Vec4& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < kN; ++i) s += a[i] * b[i];
    return s;
}
Vec4 MatVec(const Mat4& M, const Vec4& v) {
    Vec4 r{};
    for (std::size_t i = 0; i < kN; ++i) r[i] = Dot(M[i], v);
    return r;
}
Vec4 RowMat(const Vec4& row, const Mat4& M) {
    Vec4 r{};
    for (std::size_t j = 0; j < kN; ++j)
        for (std::size_t i = 0; i < kN; ++i) r[j] += row[i] * M[i][j];
    return r;
}
Vec4 AddScaled(const Vec4& a, const Vec4& b, double s) {
    Vec4 r{};
    for (std::size_t i = 0; i < kN; ++i) r[i] = a[i] + s * b[i];
    return r;
}
double InfDiff(const Vec4& a, const Vec4& b) {
    double r = 0.0;
    for (std::size_t i = 0; i < kN; ++i) r = std::max(r, std::abs(a[i] - b[i]));
    return r;
}
double Ms(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}

// ---- coefficient tables ---------------------------------------------------------------
struct Coefficients {
    Vec4 R{};
    Mat4 Fb{};
    std::vector<Vec4> HFm;               // H F^m
    std::vector<Vec4> FbjR, FbjG;        // Fb^j R, Fb^j G
    double alpha(std::size_t m, std::size_t j) const { return Dot(HFm[m], FbjR[j]); }
    double beta(std::size_t m, std::size_t j) const { return Dot(HFm[m], FbjG[j]); }
    double gamma(std::size_t l, const Vec4& G) const { return Dot(HFm[l], G); }
};

Coefficients BuildCoefficients(const ControllerProfile& p, std::size_t maxPower) {
    Coefficients c;
    // Deadbeat R = F^n O^{-1} e_n, which is the last column of M_u = F^n O^{-1}.
    for (std::size_t i = 0; i < kN; ++i) c.R[i] = p.M_u[i][kN - 1];
    for (std::size_t i = 0; i < kN; ++i)
        for (std::size_t j = 0; j < kN; ++j) c.Fb[i][j] = p.F[i][j] - c.R[i] * p.H[j];
    c.HFm.push_back(p.H);
    for (std::size_t m = 1; m <= maxPower; ++m) c.HFm.push_back(RowMat(c.HFm.back(), p.F));
    c.FbjR.push_back(c.R);
    c.FbjG.push_back(p.G);
    for (std::size_t j = 1; j < kN; ++j) {
        c.FbjR.push_back(MatVec(c.Fb, c.FbjR.back()));
        c.FbjG.push_back(MatVec(c.Fb, c.FbjG.back()));
    }
    // Sanity: Fb^n must vanish.
    Vec4 probe{1.0, -0.5, 0.25, 2.0};
    for (std::size_t j = 0; j < kN; ++j) probe = MatVec(c.Fb, probe);
    double mx = 0.0;
    for (double v : probe) mx = std::max(mx, std::abs(v));
    if (mx > 1e-9) throw std::runtime_error("deadbeat check failed: Fb^n != 0");
    return c;
}

// ---- CKKS helpers ---------------------------------------------------------------------
// OpenFHE's encoder throws "Scaling factor too small" for nonzero |v| * 2^59 < 1 (e.g. a
// structurally-zero coefficient computed as 1e-19, or a converged plant output).  Such values
// are far below CKKS precision, so they are flushed to zero.
constexpr double kFlushToZero = 1e-15;
// Level at which fresh inputs are encrypted.  0 = top level.  With --compact-inputs every input
// is encrypted with just the levels its circuit needs (at most two multiplications before the
// bootstrap input level), which shrinks ciphertexts and speeds up every operation on them.
uint32_t gInputLevel = 0;  // inputs that pass through two operations (y, packed y coefficients)
uint32_t gCoefLevel = 0;   // inputs that pass through one multiplication (alpha, eta, xc0)
Cipher EncryptVector(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk,
                     std::vector<double> values, uint32_t slots, uint32_t level) {
    for (double& v : values)
        if (std::abs(v) < kFlushToZero) v = 0.0;
    values.resize(slots, 0.0);
    return cc->Encrypt(pk, cc->MakeCKKSPackedPlaintext(values, 1, level, nullptr, slots));
}
Cipher EncryptVector(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk,
                     std::vector<double> values, uint32_t slots) {
    return EncryptVector(cc, pk, std::move(values), slots, gInputLevel);
}
Cipher EncryptScalar(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, double v,
                     uint32_t slots = kSlots) {
    return EncryptVector(cc, pk, {v}, slots);
}
Cipher EncryptCoef(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, double v) {
    return EncryptVector(cc, pk, {v}, kSlots, gCoefLevel);
}

struct Decoded {
    Vec4 slots{};
    double logError = 0.0;  // OpenFHE's imaginary-part based log2 error estimate
};
Decoded Decode(const CryptoContext<DCRTPoly>& cc, const PrivateKey<DCRTPoly>& sk, const Cipher& c) {
    Plaintext p;
    cc->Decrypt(sk, c, &p);
    p->SetLength(kSlots);
    const auto v = p->GetRealPackedValue();
    return {{v[0], v[1], v[2], v[3]}, p->GetLogError()};
}

Cipher NormalizeScaleDegree(const CryptoContext<DCRTPoly>& cc, Cipher c) {
    while (c->GetNoiseScaleDeg() > 1) c = cc->GetScheme()->ModReduceInternal(c, 1);
    return c;
}
Cipher ReduceToLevel(const CryptoContext<DCRTPoly>& cc, Cipher c, std::size_t level) {
    if (c->GetLevel() > level) throw std::runtime_error("ciphertext already below bootstrap input level");
    if (c->GetLevel() < level) c = cc->GetScheme()->LevelReduceInternal(c, level - c->GetLevel());
    return c;
}

// y-part of u(t) with packed outputs: <C, Y>, Y = [y(t-1), y(t-2), ..., 0, ...] (S slots).
// One ct-ct product and log2(S) rotate-and-add steps put the sum in every slot; a public
// period-4 mask [1,0,0,0,1,0,0,0,...] then turns it into the 4-slot sparse scalar [s,0,0,0].
// Every input is a fresh (top-level) ciphertext, so the mask level does not touch the
// post-bootstrap level budget, and the bootstrapped u keeps the sparse layout.
Cipher PackedInnerProduct(const CryptoContext<DCRTPoly>& cc, const Cipher& coefVec, const Cipher& reg,
                          uint32_t S, std::size_t& rotations) {
    Cipher p = NormalizeScaleDegree(cc, cc->EvalMult(coefVec, reg));
    for (uint32_t s = 1; s < S; s <<= 1) {
        p = cc->EvalAdd(p, cc->EvalRotate(p, static_cast<int32_t>(s)));
        ++rotations;
    }
    std::vector<double> mask(S, 0.0);
    for (uint32_t i = 0; i < S; i += kSlots) mask[i] = 1.0;
    auto pt = cc->MakeCKKSPackedPlaintext(mask, p->GetNoiseScaleDeg(), p->GetLevel(), nullptr, S);
    p = cc->GetScheme()->ModReduceInternal(cc->EvalMult(p, pt), 1);
    p->SetSlots(kSlots);  // content is period-4 now: [s,0,0,0] in the 4-slot view
    return p;
}

uint32_t EstimateBootstrapPrecisionBits(const CryptoContext<DCRTPoly>& cc, const KeyPair<DCRTPoly>& kp,
                                        std::size_t inputLevel) {
    Cipher probe = ReduceToLevel(cc, EncryptScalar(cc, kp.publicKey, 0.5), inputLevel);
    const double before = Decode(cc, kp.secretKey, probe).slots[0];
    const double after = Decode(cc, kp.secretKey, cc->EvalBootstrap(probe)).slots[0];
    const double bits = -std::log2(std::max(std::abs(after - before), 1e-16));
    return static_cast<uint32_t>(std::clamp(static_cast<int>(std::floor(bits - 5.0)), 1, 40));
}

// ---- real-time schedule model ---------------------------------------------------------
// Server computes u(t) right after y(t) arrives at t*Ts (duration uMs[t]).  A separate
// bootstrap worker processes window inputs FIFO, each released when its u is computed.
// Input s must be ready by the start of the first block that uses it (deadline[s]*Ts).
struct BootJob {
    std::size_t s, deadlineStep;
    double ms;
};
bool Feasible(double Ts, const std::vector<double>& uMs, const std::vector<BootJob>& jobs) {
    for (double u : uMs)
        if (u > Ts) return false;
    double worker = 0.0;
    for (const auto& job : jobs) {
        const double release = job.s * Ts + uMs[job.s];
        worker = std::max(worker, release) + job.ms;
        if (worker > job.deadlineStep * Ts) return false;
    }
    return true;
}
double MinFeasibleTs(const std::vector<double>& uMs, const std::vector<BootJob>& jobs) {
    double lo = 0.0, hi = 1.0;
    while (!Feasible(hi, uMs, jobs)) hi *= 2.0;
    for (int it = 0; it < 60; ++it) {
        const double mid = 0.5 * (lo + hi);
        (Feasible(mid, uMs, jobs) ? hi : lo) = mid;
    }
    return hi;
}
}  // namespace

int main(int argc, char** argv) {
    Options o;
    std::ofstream failure;
    std::size_t currentStep = 0;
    std::string stage = "setup";
    try {
        o = ParseOptions(argc, argv);
        const auto& p = GetControllerProfile(o.profile);
        const std::size_t K = o.K, d = o.d, n = kN;

        // First block whose window [tau_k - n + 1, tau_k] is entirely at t >= 0.
        const std::size_t kStar = (n - 1 + d + K - 1) / K;
        const std::size_t tFirst = kStar * K;
        const std::size_t mMin = d - 1, mMax = d + K - 2;
        const std::size_t gammaMax = std::max(mMax, tFirst);  // gamma[l] for l < max(m, t)
        const Coefficients coef = BuildCoefficients(p, std::max(mMax, tFirst));

        // Which inputs enter a window, and when they are first needed.
        std::vector<long> firstUse(o.steps, -1);
        for (std::size_t k = kStar; k * K < o.steps; ++k) {
            const std::size_t tau = k * K - d;
            for (std::size_t j = 0; j < n; ++j) {
                const std::size_t s = tau - j;
                if (firstUse[s] < 0) firstUse[s] = static_cast<long>(k * K);
            }
        }

        std::filesystem::create_directories(o.resultsDir);
        std::ofstream trace(o.resultsDir + "/sample_trace.csv"), timing(o.resultsDir + "/timing.csv"),
            boot(o.resultsDir + "/bootstrap.csv"), coefCsv(o.resultsDir + "/coefficients.csv");
        failure.open(o.resultsDir + "/failure.csv");
        failure << "step,stage,message\n";

        // ---- CKKS context ---------------------------------------------------------------
        const std::vector<uint32_t> levelBudget{1, 1}, bsgsDim{0, 0};
        const SecretKeyDist dist = UNIFORM_TERNARY;
        const uint32_t bootDepth = FHECKKSRNS::GetBootstrapDepth(levelBudget, dist);
        const uint32_t depth = bootDepth + o.levelsAfterBootstrap + (o.bootIters - 1);
        const uint32_t bootInputLevel = depth - 1;
        CCParams<CryptoContextCKKSRNS> params;
        params.SetSecretKeyDist(dist);
        if (o.secure128) {
            params.SetSecurityLevel(HEStd_128_classic);
        } else {
            params.SetSecurityLevel(HEStd_NotSet);
            params.SetRingDim(4096);
        }
        params.SetBatchSize(kSlots);
        params.SetNumLargeDigits(3);
        params.SetKeySwitchTechnique(HYBRID);
        params.SetScalingTechnique(FLEXIBLEAUTO);
        params.SetScalingModSize(59);
        params.SetFirstModSize(60);
        params.SetMultiplicativeDepth(depth);
        auto cc = GenCryptoContext(params);
        for (auto f : {PKE, KEYSWITCH, LEVELEDSHE, ADVANCEDSHE, FHE}) cc->Enable(f);
        cc->EvalBootstrapSetup(levelBudget, bsgsDim, kSlots, kCorrectionFactor);
        auto kp = cc->KeyGen();
        cc->EvalMultKeyGen(kp.secretKey);
        cc->EvalBootstrapKeyGen(kp.secretKey, kSlots);  // also generates the conjugation key (M-1)
        const uint32_t conjIndex = cc->GetCyclotomicOrder() - 1;
        // y-part: mult + mask (2 levels); alpha/eta parts: one mult.  Inputs at depth-3 end at or
        // below the bootstrap input level (depth-1).
        if (o.compactInputs) {
            gInputLevel = depth - 3;
            gCoefLevel = depth - 2;
        }

        stage = "bootstrap_warmup";
        (void)cc->EvalBootstrap(ReduceToLevel(cc, EncryptScalar(cc, kp.publicKey, 0.0125), bootInputLevel));
        const uint32_t precBits = o.bootIters > 1 ? EstimateBootstrapPrecisionBits(cc, kp, bootInputLevel) : 0;

        std::cout << std::setprecision(6)
                  << "profile=" << p.name << ", rho_F=" << p.rho_F << ", rho_Acl=" << p.rho_Acl << '\n'
                  << "K=" << K << ", d=" << d << ", n=" << n << ", steps=" << o.steps
                  << ", first_window_step=" << tFirst << ", m_range=[" << mMin << "," << mMax << "]"
                  << (d + n > K + 1 ? "  WARNING: d > K-n+1, block-reset jumps can accumulate" : "") << '\n'
                  << "ring_dimension=" << cc->GetRingDimension() << ", depth=" << depth
                  << ", bootstrap_depth=" << bootDepth << ", levels_after_bootstrap=" << o.levelsAfterBootstrap
                  << ", boot_iters=" << o.bootIters << ", iterative_precision_bits=" << precBits
                  << ", clear_imag=" << (o.keepImag ? "off" : "on")
                  << ", security=" << (o.secure128 ? "HEStd_128_classic" : "HEStd_NotSet")
                  << ", input_level=" << gInputLevel << '\n';
#ifdef _OPENMP
        std::cout << "OpenMP_max_threads=" << omp_get_max_threads() << '\n';
#endif

        // ---- encrypted coefficients (client setup) ----------------------------------------
        stage = "encrypt_coefficients";
        const double imagFold = o.keepImag ? 1.0 : 0.5;  // window inputs are 2*Re(u) when cleared
        coefCsv << "kind,m_or_t,j,value\n" << std::setprecision(17);
        std::vector<std::array<Cipher, kN>> alphaCt(K), betaCt(K);
        double maxAlpha = 0.0, maxBeta = 0.0;
        for (std::size_t m = mMin; m <= mMax; ++m)
            for (std::size_t j = 0; j < n; ++j) {
                const double a = coef.alpha(m, j), b = coef.beta(m, j);
                maxAlpha = std::max(maxAlpha, std::abs(a));
                maxBeta = std::max(maxBeta, std::abs(b));
                alphaCt[m - mMin][j] = EncryptCoef(cc, kp.publicKey, imagFold * a);
                if (!o.packY) betaCt[m - mMin][j] = EncryptScalar(cc, kp.publicKey, b);
                coefCsv << "alpha," << m << ',' << j << ',' << a << "\nbeta," << m << ',' << j << ',' << b << '\n';
            }
        std::vector<Cipher> gammaCt;
        double maxGamma = 0.0;
        for (std::size_t l = 0; l < gammaMax; ++l) {
            const double g = coef.gamma(l, p.G);
            maxGamma = std::max(maxGamma, std::abs(g));
            if (!o.packY) gammaCt.push_back(EncryptScalar(cc, kp.publicKey, g));
            coefCsv << "gamma," << l << ",0," << g << '\n';
        }
        // Initial-phase coefficients H F^t e_i are each used once; the client-side encryption is
        // done lazily right before step t (outside the server timing) to bound memory at N = 65536.
        for (std::size_t t = 0; t < tFirst && t < o.steps; ++t)
            for (std::size_t i = 0; i < n; ++i) coefCsv << "eta," << t << ',' << i << ',' << coef.HFm[t][i] << '\n';
        std::array<Cipher, kN> x0Ct;
        for (std::size_t i = 0; i < n; ++i) x0Ct[i] = EncryptCoef(cc, kp.publicKey, p.xc0[i]);
        std::cout << "max|alpha|=" << maxAlpha << ", max|beta|=" << maxBeta << ", max|gamma|=" << maxGamma << '\n';

        // ---- packed y (optional) ------------------------------------------------------------
        // Register for block k holds y(kK-d-n+1 .. t-1) as [y(t-1), y(t-2), ...]; it is rebuilt
        // for every block, so no ciphertext lives longer than d+n-1+K steps.
        uint32_t S = kSlots;
        std::vector<Cipher> coefVecCt;
        {
            std::vector<int32_t> rot;
            if (o.packU)
                for (int32_t r : {1, 2, 3, -1, -2, -3}) rot.push_back(r);  // window pack / unpack
            if (o.packY) {
                while (S < std::max<std::size_t>(d + n + K - 2, tFirst)) S <<= 1;
                rot.push_back(-1);
                for (uint32_t s = 1; s < S; s <<= 1) rot.push_back(static_cast<int32_t>(s));
            }
            std::sort(rot.begin(), rot.end());
            rot.erase(std::unique(rot.begin(), rot.end()), rot.end());
            if (!rot.empty()) cc->EvalRotateKeyGen(kp.secretKey, rot);
        }
        if (o.packY) {
            for (std::size_t m = mMin; m <= mMax; ++m) {
                std::vector<double> c(S, 0.0);
                for (std::size_t l = 0; l < m; ++l) c[l] = coef.gamma(l, p.G);
                for (std::size_t j = 0; j < n; ++j) c[m + j] = coef.beta(m, j);
                coefVecCt.push_back(EncryptVector(cc, kp.publicKey, c, S));
            }
        }
        std::cout << "y_packing=" << (o.packY ? "on" : "off") << ", y_slots=" << S
                  << ", window_packing=" << (o.packU ? "on (one bootstrap per window)" : "off") << '\n';

        trace << "global_sample,block,within_block,phase,m,y_plain,y_physical,u_plain,u_ckks,u_formula,"
                 "e_u_abs,e_local_abs,u_inactive_leakage,u_level,u_log_error,bootstrapped";
        for (auto nm : {"plant_plain", "plant_ckks"})
            for (std::size_t i = 0; i < n; ++i) trace << ',' << nm << '_' << i;
        trace << '\n' << std::setprecision(17);
        timing << "global_sample,y_encrypt_ms,u_eval_ms,register_ms,u_decode_ms,ctct_mults,rotations,bootstrap_ms\n";
        boot << "global_sample,deadline_step,slack_steps,u_applied,pre_active_error,post_active_error,"
                "post_imag_log_error,bootstrap_ms\n"
             << std::setprecision(17);

        // ---- closed-loop simulation -----------------------------------------------------
        Vec4 xpRef = p.xp0, xcRef = p.xc0, xpEnc = p.xp0;
        std::map<std::size_t, Cipher> cy, bu;       // encrypted y history, bootstrapped window inputs
        std::map<std::size_t, Cipher> reg;          // packed-y registers, one per block
        Cipher regInit;                             // packed-y register for the initial phase
        std::map<std::size_t, Cipher> preU;         // pre-bootstrap window inputs (pack-u)
        std::map<std::size_t, std::array<Cipher, kN>> rw;  // Rot_j(bootstrapped window), per block
        std::vector<double> ySent(o.steps), uApplied(o.steps), uMsAll(o.steps);
        std::vector<BootJob> jobs;
        double maxEu = 0.0, maxPlant = 0.0, maxLocal = 0.0, maxLeak = 0.0;

        for (std::size_t t = 0; t < o.steps; ++t) {
            currentStep = t;
            const std::size_t k = t / K;
            const bool init = k < kStar;
            const double yRef = Dot(p.C, xpRef), uRef = Dot(p.H, xcRef);

            stage = "encrypt_y";
            ySent[t] = Dot(p.C, xpEnc);
            auto t0 = Clock::now();
            cy[t] = EncryptScalar(cc, kp.publicKey, ySent[t], o.packY ? S : kSlots);
            const double encMs = Ms(t0, Clock::now());

            // client-side, lazily encrypted initial-phase coefficients (not part of server timing)
            std::array<Cipher, kN> etaCt;
            Cipher initVecCt;
            if (init) {
                stage = "encrypt_init_coefficients";
                for (std::size_t i = 0; i < n; ++i) etaCt[i] = EncryptCoef(cc, kp.publicKey, coef.HFm[t][i]);
                if (o.packY && t > 0) {
                    std::vector<double> c(S, 0.0);
                    for (std::size_t l = 0; l < t; ++l) c[l] = coef.gamma(l, p.G);
                    initVecCt = EncryptVector(cc, kp.publicKey, c, S);
                }
            }

            stage = "u_eval";
            std::size_t mults = 0, rotations = 0, m = 0;
            Cipher acc, yPart;
            auto term = [&](const Cipher& a, const Cipher& b) {
                Cipher prod = cc->EvalMult(a, b);
                acc = acc ? cc->EvalAdd(acc, prod) : prod;
                ++mults;
            };
            double uFormula = 0.0;  // same formula in plain arithmetic on the applied u / sent y
            t0 = Clock::now();
            if (init) {
                for (std::size_t i = 0; i < n; ++i) {
                    term(etaCt[i], x0Ct[i]);
                    uFormula += coef.HFm[t][i] * p.xc0[i];
                }
                for (std::size_t l = 0; l < t; ++l) {
                    if (!o.packY) term(gammaCt[l], cy.at(t - 1 - l));
                    uFormula += coef.gamma(l, p.G) * ySent[t - 1 - l];
                }
                if (o.packY && t > 0) {
                    yPart = PackedInnerProduct(cc, initVecCt, regInit, S, rotations);
                    ++mults;
                }
            } else {
                const std::size_t tau = k * K - d;
                m = t - tau - 1;
                for (std::size_t j = 0; j < n; ++j) {
                    // pack-u: Rot_j(W) has u(tau-j) in slot 0; the sparse alpha keeps [.,0,0,0]
                    term(alphaCt[m - mMin][j], o.packU ? rw.at(k)[j] : bu.at(tau - j));
                    if (!o.packY) term(betaCt[m - mMin][j], cy.at(tau - j));
                    uFormula += coef.alpha(m, j) * uApplied[tau - j] + coef.beta(m, j) * ySent[tau - j];
                }
                for (std::size_t l = 0; l < m; ++l) {
                    if (!o.packY) term(gammaCt[l], cy.at(t - 1 - l));
                    uFormula += coef.gamma(l, p.G) * ySent[t - 1 - l];
                }
                if (o.packY) {
                    yPart = PackedInnerProduct(cc, coefVecCt[m - mMin], reg.at(k), S, rotations);
                    ++mults;
                }
            }
            acc = NormalizeScaleDegree(cc, acc);
            if (yPart) acc = cc->EvalAdd(acc, yPart);
            const double uMs = Ms(t0, Clock::now());

            stage = "u_decode";
            t0 = Clock::now();
            const Decoded du = Decode(cc, kp.secretKey, acc);
            const double decMs = Ms(t0, Clock::now());
            const double uCkks = du.slots[0];
            uApplied[t] = uCkks;
            uMsAll[t] = uMs;

            auto clearImag = [&](Cipher c) {
                if (o.keepImag) return c;
                return cc->EvalAdd(c, cc->EvalAutomorphism(c, conjIndex, cc->GetEvalAutomorphismKeyMap(c->GetKeyTag())));
            };
            double bootMs = 0.0;
            if (o.packU) {
                // Keep the sparse pre-bootstrap u; once the newest element of a window exists,
                // pack W = sum_j Rot_{-j}(u(tau-j)) = [u(tau), u(tau-1), u(tau-2), u(tau-3)] and
                // bootstrap it once.  Rotations cost no level; the 4-slot bootstrap costs the same
                // whether one or four slots are used.
                if (firstUse[t] >= 0) preU[t] = acc;
                const std::size_t kNext = (t + d) / K;
                if ((t + d) % K == 0 && kNext >= kStar && kNext * K < o.steps) {
                    stage = "window_bootstrap";
                    t0 = Clock::now();
                    Cipher w = preU.at(t);
                    for (std::size_t j = 1; j < n; ++j)
                        w = cc->EvalAdd(w, cc->EvalRotate(preU.at(t - j), -static_cast<int32_t>(j)));
                    const Decoded pre = Decode(cc, kp.secretKey, w);
                    Cipher in = ReduceToLevel(cc, w, bootInputLevel);
                    Cipher out = clearImag(o.bootIters > 1 ? cc->EvalBootstrap(in, o.bootIters, precBits)
                                                           : cc->EvalBootstrap(in));
                    std::array<Cipher, kN> rots;
                    rots[0] = out;
                    for (std::size_t j = 1; j < n; ++j) rots[j] = cc->EvalRotate(out, static_cast<int32_t>(j));
                    rw[kNext] = rots;
                    bootMs = Ms(t0, Clock::now());
                    const Decoded post = Decode(cc, kp.secretKey, out);
                    double preErr = 0.0, postErr = 0.0;
                    for (std::size_t j = 0; j < n; ++j) {
                        preErr = std::max(preErr, std::abs(pre.slots[j] - uApplied[t - j]));
                        postErr = std::max(postErr, std::abs(post.slots[j] * imagFold - uApplied[t - j]));
                    }
                    boot << t << ',' << kNext * K << ',' << kNext * K - t << ',' << uCkks << ',' << preErr << ','
                         << postErr << ',' << post.logError << ',' << bootMs << '\n';
                    jobs.push_back({t, kNext * K, bootMs});
                }
            } else if (firstUse[t] >= 0) {
                stage = "bootstrap";
                t0 = Clock::now();
                Cipher in = ReduceToLevel(cc, acc, bootInputLevel);
                Cipher out = o.bootIters > 1 ? cc->EvalBootstrap(in, o.bootIters, precBits) : cc->EvalBootstrap(in);
                if (!o.keepImag)
                    out = cc->EvalAdd(out, cc->EvalAutomorphism(out, conjIndex,
                                                                cc->GetEvalAutomorphismKeyMap(out->GetKeyTag())));
                bootMs = Ms(t0, Clock::now());
                bu[t] = out;
                const Decoded post = Decode(cc, kp.secretKey, out);
                const double postVal = post.slots[0] * imagFold;
                boot << t << ',' << firstUse[t] << ',' << firstUse[t] - static_cast<long>(t) << ',' << uCkks
                     << ',' << std::abs(du.slots[0] - uFormula) << ',' << std::abs(postVal - uCkks) << ','
                     << post.logError << ',' << bootMs << '\n';
                jobs.push_back({t, static_cast<std::size_t>(firstUse[t]), bootMs});
            }

            // shift y(t) into every packed register that still needs it (after u(t) was computed)
            double regMs = 0.0;
            if (o.packY) {
                stage = "register_update";
                t0 = Clock::now();
                auto shiftIn = [&](Cipher& r) {
                    r = cc->EvalAdd(cc->EvalRotate(r, -1), cy[t]);
                    ++rotations;
                };
                if (t + 2 <= tFirst) {
                    if (t == 0) regInit = cy[t];
                    else shiftIn(regInit);
                }
                for (std::size_t kk = std::max(kStar, k); kk * K <= t + d + n - 1; ++kk) {
                    const std::size_t start = kk * K - d - n + 1;
                    if (t < start || t + 2 > kk * K + K) continue;
                    if (t == start) reg[kk] = cy[t];
                    else shiftIn(reg.at(kk));
                }
                if (t + 1 == k * K + K) reg.erase(k);  // block k is finished
                if (t + 1 == tFirst) regInit.reset();
                cy.erase(t);
                regMs = Ms(t0, Clock::now());
                uMsAll[t] += regMs;
            }

            const double eU = std::abs(uCkks - uRef), eLocal = std::abs(uCkks - uFormula);
            const double leak = std::max({std::abs(du.slots[1]), std::abs(du.slots[2]), std::abs(du.slots[3])});
            maxEu = std::max(maxEu, eU);
            maxLocal = std::max(maxLocal, eLocal);
            maxLeak = std::max(maxLeak, leak);
            maxPlant = std::max(maxPlant, InfDiff(xpEnc, xpRef));
            trace << t << ',' << k << ',' << t % K << ',' << (init ? "init" : "window") << ',' << m << ','
                  << yRef << ',' << ySent[t] << ',' << uRef << ',' << uCkks << ',' << uFormula << ',' << eU << ','
                  << eLocal << ',' << leak << ',' << acc->GetLevel() << ',' << du.logError << ','
                  << (firstUse[t] >= 0);
            for (double v : xpRef) trace << ',' << v;
            for (double v : xpEnc) trace << ',' << v;
            trace << '\n';
            timing << t << ',' << encMs << ',' << uMs << ',' << regMs << ',' << decMs << ',' << mults << ','
                   << rotations << ',' << bootMs << '\n';

            // plant / reference controller updates
            xpRef = AddScaled(MatVec(p.A, xpRef), p.B, uRef);
            xcRef = AddScaled(MatVec(p.F, xcRef), p.G, yRef);
            xpEnc = AddScaled(MatVec(p.A, xpEnc), p.B, uCkks);

            // drop ciphertexts no future step can reference
            if (t + 1 >= tFirst) {
                const std::size_t kn = (t + 1) / K;
                const std::size_t oldest = kn * K - d - (n - 1);
                cy.erase(cy.begin(), cy.lower_bound(oldest));
                bu.erase(bu.begin(), bu.lower_bound(oldest));
                preU.erase(preU.begin(), preU.lower_bound(oldest));
                rw.erase(rw.begin(), rw.lower_bound(kn));
            }
        }

        // ---- summary ----------------------------------------------------------------------
        double uSum = 0.0, uMax = 0.0, bSum = 0.0;
        std::size_t steady = 0;
        for (std::size_t t = tFirst; t < o.steps; ++t) { uSum += uMsAll[t]; uMax = std::max(uMax, uMsAll[t]); ++steady; }
        for (const auto& j : jobs) bSum += j.ms;
        const double minTs = MinFeasibleTs(uMsAll, jobs);
        std::cout << std::setprecision(6) << "completed_steps=" << o.steps << '\n'
                  << "max |u_ckks - u_plain|=" << maxEu << '\n'
                  << "max plant error inf=" << maxPlant << '\n'
                  << "max |u_ckks - u_formula| (per-step CKKS error)=" << maxLocal << '\n'
                  << "max u inactive-slot leakage=" << maxLeak << '\n'
                  << "bootstraps=" << jobs.size() << " (" << double(jobs.size()) / o.steps << " per step)"
                  << ", mean_bootstrap_ms=" << (jobs.empty() ? 0.0 : bSum / jobs.size()) << '\n'
                  << "steady_u_eval_ms mean=" << (steady ? uSum / steady : 0.0) << " max=" << uMax << '\n'
                  << "min_feasible_sampling_period_ms=" << minTs << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAILED at step " << currentStep << " (" << stage << "): " << e.what() << '\n';
        if (failure.is_open()) failure << currentStep << ',' << stage << ",\"" << e.what() << "\"\n";
        std::cout << "completed_steps=" << currentStep << '\n';
        return 1;
    }
}
