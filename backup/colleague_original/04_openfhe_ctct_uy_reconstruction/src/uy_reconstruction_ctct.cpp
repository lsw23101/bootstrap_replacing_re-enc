#include "openfhe.h"
#include "../../include/controller_profiles.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace lbcrypto;
using namespace tutorial_profiles;

namespace {
Mat4 A{}, F{}, O4{}, T4{}, M_u{}, M_y{}, P{}, q{};
Vec4 B{}, C{}, G{}, H{}, Xp0{}, Xc0{};
double rho_F=0.0, rho_Acl=0.0, cond_O4=0.0, M_u_inf=0.0, M_y_inf=0.0, Gamma_B=0.0;
int rank_O4=4;

void LoadProfile(const ControllerProfile& p) {
    A=p.A; B=p.B; C=p.C; F=p.F; G=p.G; H=p.H; Xp0=p.xp0; Xc0=p.xc0;
    O4=p.O4; T4=p.T4; M_u=p.M_u; M_y=p.M_y; P=p.P; q=p.q;
    rho_F=p.rho_F; rho_Acl=p.rho_Acl; cond_O4=p.cond_O4;
    M_u_inf=p.M_u_inf; M_y_inf=p.M_y_inf; Gamma_B=p.Gamma_B;
}

constexpr uint32_t kSlots = 4;
constexpr uint32_t kCorrectionFactor = 13;
constexpr uint32_t kBootstrapIterations = 2;
constexpr double kBootstrapTarget = 0.50;
constexpr double kApplicationResolution = 1.0e-6;
constexpr double kMaxBootstrapScale = 1.0e5;
constexpr double kScaledBootstrapGuard = 0.90;
// The production fix uses two-iteration (META-BTS) bootstrapping. OpenFHE's
// iterative bootstrap requires one extra level for the second iteration.
//
// The recursive ct-ct path itself needs three levels available after bootstrap.
// If public pre-bootstrap dynamic scaling is enabled, multiplying c_u by the
// block-dependent scalar and normalizing its scale degree consumes one additional
// level before the next bootstrap. That extra level must be present in the modulus
// chain; otherwise block 1 tries to move a level-20 ciphertext back to level 19.
constexpr uint32_t kDefaultLevelsAfterBootstrap = 3;
constexpr uint32_t kDynamicScalingExtraLevel = 1;
using Cipher = Ciphertext<DCRTPoly>;
using Clock = std::chrono::steady_clock;

struct Options {
    std::string profile = "stable";
    std::size_t blocks = 3;
    std::string resultsDir = "results/ctct_smoke";
    bool freshBootstrapProbe = false;
    bool singleBootstrap = false;
    bool disableDynamicScaling = true;
    uint32_t levelsAfterBootstrap = kDefaultLevelsAfterBootstrap;
};

struct CipherMetadata {
    std::size_t level = 0;
    std::size_t towers = 0;
    std::size_t noiseScaleDegree = 0;
    double scalingFactor = 0.0;
};

struct FailureState {
    std::size_t block = 0;
    std::size_t completedBlocks = 0;
    std::size_t completedSamples = 0;
    long sample = -1;
    std::string stage = "setup";
    std::string operation = "none";
    std::string branch = "production";
    CipherMetadata metadata{};
    bool hasMetadata = false;
} gFailure;

struct EncryptedCoefficients {
    // pMasks[sample][state_index] encrypts [P[sample][state_index],0,0,0].
    std::array<std::array<Cipher, 4>, 4> pMasks;
    // qMasks[sample][history] encrypts [q[sample][history],0,0,0].
    std::array<std::array<Cipher, 4>, 4> qMasks;
    // Reconstruction masks: [0,...,M_[row][column],...,0].
    std::array<std::array<Cipher, 4>, 4> muMasks;
    std::array<std::array<Cipher, 4>, 4> myMasks;
};

struct CsvFiles {
    std::ofstream sample;
    std::ofstream block;
    std::ofstream bootstrap;
    std::ofstream timing;
    std::ofstream ctct;
    std::ofstream failure;
    std::ofstream mode;
};

Options ParseOptions(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--profile" && i + 1 < argc)
            options.profile = argv[++i];
        else if (arg == "--blocks" && i + 1 < argc)
            options.blocks = std::stoul(argv[++i]);
        else if (arg == "--results-dir" && i + 1 < argc)
            options.resultsDir = argv[++i];
        else if (arg == "--fresh-bootstrap-probe")
            options.freshBootstrapProbe = true;
        else if (arg == "--single-bootstrap")
            options.singleBootstrap = true;
        else if (arg == "--disable-dynamic-scaling")
            options.disableDynamicScaling = true;
        else if (arg == "--enable-dynamic-scaling")
            options.disableDynamicScaling = false;
        else if (arg == "--levels-after-bootstrap" && i + 1 < argc)
            options.levelsAfterBootstrap = static_cast<uint32_t>(std::stoul(argv[++i]));
        else
            throw std::invalid_argument(
                "usage: uy_reconstruction_ctct [--profile stable|unstable_low|unstable_high] [--blocks N] "
                "[--results-dir PATH] [--fresh-bootstrap-probe] [--single-bootstrap] "
                "[--disable-dynamic-scaling] [--enable-dynamic-scaling] "
                "[--levels-after-bootstrap N]");
    }
    if (options.blocks == 0)
        throw std::invalid_argument("--blocks must be positive");
    if (options.levelsAfterBootstrap < 3)
        throw std::invalid_argument("--levels-after-bootstrap must be at least 3 for this ct-ct circuit");
    return options;
}

double Milliseconds(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

Vec4 Add(const Vec4& a, const Vec4& b) {
    Vec4 result{};
    for (std::size_t i = 0; i < 4; ++i)
        result[i] = a[i] + b[i];
    return result;
}

Vec4 Subtract(const Vec4& a, const Vec4& b) {
    Vec4 result{};
    for (std::size_t i = 0; i < 4; ++i)
        result[i] = a[i] - b[i];
    return result;
}

Vec4 Scale(const Vec4& a, double scalar) {
    Vec4 result{};
    for (std::size_t i = 0; i < 4; ++i)
        result[i] = a[i] * scalar;
    return result;
}

Vec4 MatrixVector(const Mat4& matrix, const Vec4& vector) {
    Vec4 result{};
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            result[i] += matrix[i][j] * vector[j];
    return result;
}

Mat4 MatrixMultiply(const Mat4& lhs, const Mat4& rhs) {
    Mat4 result{};
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t k = 0; k < 4; ++k)
            for (std::size_t j = 0; j < 4; ++j)
                result[i][j] += lhs[i][k] * rhs[k][j];
    return result;
}

double Dot(const Vec4& a, const Vec4& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < 4; ++i)
        result += a[i] * b[i];
    return result;
}

double InfinityNorm(const Vec4& value) {
    double result = 0.0;
    for (double element : value)
        result = std::max(result, std::abs(element));
    return result;
}

double ActiveSlotError(const Vec4& value, double expected) {
    return std::abs(value[0] - expected);
}

double InactiveSlotLeakage(const Vec4& value) {
    double result = 0.0;
    for (std::size_t slot = 1; slot < 4; ++slot)
        result = std::max(result, std::abs(value[slot]));
    return result;
}

double SparseScalarErrorInf(const Vec4& value, double expected) {
    return std::max(ActiveSlotError(value, expected), InactiveSlotLeakage(value));
}

CipherMetadata GetMetadata(const Cipher& ciphertext) {
    return {ciphertext->GetLevel(), ciphertext->GetElements()[0].GetNumOfElements(),
            ciphertext->GetNoiseScaleDeg(), ciphertext->GetScalingFactor()};
}

Cipher ReduceToLevel(const CryptoContext<DCRTPoly>& cc, Cipher ciphertext, std::size_t level) {
    if (ciphertext->GetLevel() > level)
        throw std::runtime_error("cannot raise a ciphertext to a lower level index");
    if (ciphertext->GetLevel() < level)
        ciphertext = cc->GetScheme()->LevelReduceInternal(ciphertext, level - ciphertext->GetLevel());
    return ciphertext;
}

Cipher NormalizeScaleDegree(const CryptoContext<DCRTPoly>& cc, Cipher ciphertext) {
    while (ciphertext->GetNoiseScaleDeg() > 1)
        ciphertext = cc->GetScheme()->ModReduceInternal(ciphertext, 1);
    return ciphertext;
}

Vec4 Decode(const CryptoContext<DCRTPoly>& cc, const PrivateKey<DCRTPoly>& secretKey,
            const Cipher& ciphertext) {
    Plaintext plaintext;
    cc->Decrypt(secretKey, ciphertext, &plaintext);
    plaintext->SetLength(kSlots);
    const auto values = plaintext->GetRealPackedValue();
    return {values[0], values[1], values[2], values[3]};
}

Cipher EncryptVector(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& publicKey,
                     const Vec4& value) {
    const std::vector<double> values{value[0], value[1], value[2], value[3]};
    return cc->Encrypt(publicKey, cc->MakeCKKSPackedPlaintext(values, 1, 0, nullptr, kSlots));
}

Cipher EncryptSparseScalar(const CryptoContext<DCRTPoly>& cc,
                           const PublicKey<DCRTPoly>& publicKey, double value) {
    return EncryptVector(cc, publicKey, {value, 0.0, 0.0, 0.0});
}

Vec4 OneHot(std::size_t slot, double value) {
    Vec4 result{};
    result.at(slot) = value;
    return result;
}

void WriteMetadataHeader(std::ostream& output, const std::string& prefix) {
    output << ',' << prefix << "_level," << prefix << "_towers," << prefix
           << "_noiseScaleDegree," << prefix << "_scalingFactor";
}

void WriteMetadata(std::ostream& output, const CipherMetadata& metadata) {
    output << ',' << metadata.level << ',' << metadata.towers << ','
           << metadata.noiseScaleDegree << ',' << metadata.scalingFactor;
}

void WriteVecHeader(std::ostream& output, const std::string& prefix) {
    for (std::size_t i = 0; i < 4; ++i)
        output << ',' << prefix << '_' << i;
}

void WriteVector(std::ostream& output, const Vec4& value) {
    for (double element : value)
        output << ',' << element;
}

void OpenCsvFiles(const Options& options, CsvFiles& files) {
    std::filesystem::create_directories(options.resultsDir);
    files.sample.open(options.resultsDir + "/sample_trace.csv");
    files.block.open(options.resultsDir + "/block_metrics.csv");
    files.bootstrap.open(options.resultsDir + "/bootstrap_diagnostics.csv");
    files.timing.open(options.resultsDir + "/timing.csv");
    files.ctct.open(options.resultsDir + "/ctct_metadata.csv");
    files.failure.open(options.resultsDir + "/failure.csv");
    files.mode.open(options.resultsDir + "/error_mode_diagnostics.csv");
    for (auto* stream : {&files.sample, &files.block, &files.bootstrap, &files.timing,
                         &files.ctct, &files.failure, &files.mode})
        *stream << std::scientific << std::setprecision(17);

    files.sample << "block,within_block_sample,global_sample,y_nominal,y_physical,y_semantic,"
                    "u_nominal,u_semantic,u_ckks,e_u_semantic_abs,e_u_nominal_abs,"
                    "y_active_error,y_inactive_leakage,u_active_error,u_inactive_leakage";
    for (const auto* name : {"plant_nominal", "plant_ckks", "plant_error", "y_slots", "u_slots"})
        WriteVecHeader(files.sample, name);
    WriteMetadataHeader(files.sample, "u");
    files.sample << '\n';

    files.block << "block";
    for (const auto* name : {"controller_nominal", "controller_semantic", "controller_ckks",
                             "controller_error", "plant_error"})
        WriteVecHeader(files.block, name);
    files.block << ",controller_error_inf,plant_error_inf,semantic_identity_inf";
    WriteMetadataHeader(files.block, "x_next");
    files.block << '\n';

    files.bootstrap << "block,within_block_sample,u_semantic,bootstrap_scale,scaled_input_abs_max,num_iterations,iterative_precision_bits";
    WriteVecHeader(files.bootstrap, "u_pre");
    WriteVecHeader(files.bootstrap, "u_post");
    files.bootstrap << ",u_active_error_pre,u_inactive_leakage_pre,e_u_boot_inf,"
                       "u_active_error_post,u_inactive_leakage_post";
    WriteMetadataHeader(files.bootstrap, "before_reduction");
    WriteMetadataHeader(files.bootstrap, "before_bootstrap");
    WriteMetadataHeader(files.bootstrap, "after_bootstrap");
    files.bootstrap << ",bootstrap_ms,fresh_probe_enabled,e_B_prod,e_B_fresh,"
                       "fresh_bootstrap_ms\n";

    files.timing << "block,u_evaluation_ms,u4_bootstrap_ms,reconstruction_ms,full_block_ms\n";

    files.ctct << "block,sample,operation";
    for (const auto* name : {"lhs_input", "rhs_input", "lhs_aligned", "rhs_aligned",
                             "after_mult", "after_relinearize", "after_rescale"})
        WriteMetadataHeader(files.ctct, name);
    files.ctct << '\n';

    files.failure << "requested_blocks,completed_blocks,completed_samples,block,sample,stage,"
                     "operation,branch,message";
    WriteMetadataHeader(files.failure, "last");
    files.failure << '\n';

    files.mode << "block";
    for (const auto* name : {"same_y_reference_start", "decoded_start", "delta_start",
                             "same_y_reference_end", "decoded_end", "delta_end",
                             "F4_delta_start", "innovation"})
        WriteVecHeader(files.mode, name);
    files.mode << ",delta_start_inf,delta_end_inf,block_growth_inf,"
                  "unstable_start_abs,unstable_end_abs,unstable_growth_abs,"
                  "theoretical_unstable_block_growth,innovation_inf\n";
}

Cipher MultiplyCiphertexts(const CryptoContext<DCRTPoly>& cc, const Cipher& lhsInput,
                           const Cipher& rhsInput, CsvFiles& files, std::size_t block,
                           long sample, const std::string& operation) {
    gFailure.stage = "ct_ct_multiply";
    gFailure.operation = operation;
    gFailure.metadata = GetMetadata(lhsInput);
    gFailure.hasMetadata = true;

    Cipher lhs = NormalizeScaleDegree(cc, lhsInput);
    Cipher rhs = NormalizeScaleDegree(cc, rhsInput);
    const std::size_t targetLevel = std::max(lhs->GetLevel(), rhs->GetLevel());
    lhs = ReduceToLevel(cc, lhs, targetLevel);
    rhs = ReduceToLevel(cc, rhs, targetLevel);
    if (lhs->GetLevel() != rhs->GetLevel() || lhs->GetNoiseScaleDeg() != 1 ||
        rhs->GetNoiseScaleDeg() != 1)
        throw std::runtime_error("explicit ct-ct input alignment failed");

    const Cipher raw = cc->EvalMultNoRelin(lhs, rhs);
    const Cipher relinearized = cc->Relinearize(raw);
    const Cipher rescaled = cc->GetScheme()->ModReduceInternal(relinearized, 1);

    files.ctct << block << ',' << sample << ',' << operation;
    for (const auto& metadata : {GetMetadata(lhsInput), GetMetadata(rhsInput), GetMetadata(lhs),
                                GetMetadata(rhs), GetMetadata(raw), GetMetadata(relinearized),
                                GetMetadata(rescaled)})
        WriteMetadata(files.ctct, metadata);
    files.ctct << '\n';
    return rescaled;
}

Cipher AddAligned(const CryptoContext<DCRTPoly>& cc, Cipher lhs, Cipher rhs) {
    const std::size_t targetLevel = std::max(lhs->GetLevel(), rhs->GetLevel());
    lhs = ReduceToLevel(cc, lhs, targetLevel);
    rhs = ReduceToLevel(cc, rhs, targetLevel);
    if (lhs->GetNoiseScaleDeg() != rhs->GetNoiseScaleDeg())
        throw std::runtime_error("addition noise-scale degrees do not match");
    return cc->EvalAdd(lhs, rhs);
}

Cipher SparseInnerEncrypted(
    const CryptoContext<DCRTPoly>& cc, const Cipher& state,
    const std::array<Cipher, 4>& encryptedMasks, CsvFiles& files,
    std::size_t block, long sample, const std::string& operationPrefix) {
    Cipher sum;
    for (std::size_t stateIndex = 0; stateIndex < 4; ++stateIndex) {
        Cipher rotated = stateIndex == 0
                             ? state
                             : cc->EvalRotate(state, static_cast<int32_t>(stateIndex));
        Cipher term = MultiplyCiphertexts(
            cc, rotated, encryptedMasks[stateIndex], files, block, sample,
            operationPrefix + "_x" + std::to_string(stateIndex));
        sum = sum ? AddAligned(cc, sum, term) : term;
    }
    return sum;
}

Cipher ScatterEncrypted(
    const CryptoContext<DCRTPoly>& cc, const Cipher& sparseScalar,
    const std::array<Cipher, 4>& encryptedMasks, CsvFiles& files,
    std::size_t block, long sample, const std::string& operationPrefix,
    double coefficientScale = 1.0) {
    Cipher sum;
    for (std::size_t row = 0; row < 4; ++row) {
        Cipher placed = row == 0
                            ? sparseScalar
                            : cc->EvalRotate(sparseScalar, -static_cast<int32_t>(row));
        Cipher mask = encryptedMasks[row];
        if (coefficientScale != 1.0)
            mask = cc->EvalMult(mask, coefficientScale);
        Cipher term = MultiplyCiphertexts(
            cc, placed, mask, files, block, sample,
            operationPrefix + "_row" + std::to_string(row));
        sum = sum ? AddAligned(cc, sum, term) : term;
    }
    return sum;
}

EncryptedCoefficients EncryptCoefficients(const CryptoContext<DCRTPoly>& cc,
                                          const PublicKey<DCRTPoly>& publicKey) {
    EncryptedCoefficients encrypted;
    for (std::size_t sample = 0; sample < 4; ++sample) {
        for (std::size_t stateIndex = 0; stateIndex < 4; ++stateIndex)
            encrypted.pMasks[sample][stateIndex] =
                EncryptVector(cc, publicKey, OneHot(0, P[sample][stateIndex]));
        for (std::size_t history = 0; history < 4; ++history)
            encrypted.qMasks[sample][history] =
                EncryptVector(cc, publicKey, OneHot(0, q[sample][history]));
        for (std::size_t row = 0; row < 4; ++row) {
            encrypted.muMasks[sample][row] =
                EncryptVector(cc, publicKey, OneHot(row, M_u[row][sample]));
            encrypted.myMasks[sample][row] =
                EncryptVector(cc, publicKey, OneHot(row, M_y[row][sample]));
        }
    }
    return encrypted;
}

std::vector<double> BuildBootstrapScales(std::size_t blocks, bool disableDynamicScaling) {
    std::vector<double> scales(blocks, 1.0);
    if (disableDynamicScaling)
        return scales;

    Vec4 plant=Xp0;
    Vec4 controller=Xc0;
    for (std::size_t block = 0; block < blocks; ++block) {
        double blockMaxU = 0.0;
        for (std::size_t sample = 0; sample < 4; ++sample) {
            const double y = Dot(C, plant);
            const double u = Dot(H, controller);
            blockMaxU = std::max(blockMaxU, std::abs(u));
            plant = Add(MatrixVector(A, plant), Scale(B, u));
            controller = Add(MatrixVector(F, controller), Scale(G, y));
        }
        const double publicEnvelope = std::max(blockMaxU, kApplicationResolution);
        const double desiredScale = kBootstrapTarget / publicEnvelope;
        scales[block] = std::clamp(desiredScale, 1.0, kMaxBootstrapScale);
    }
    return scales;
}

double MaxAbs(const Vec4& value) {
    double result = 0.0;
    for (double element : value)
        result = std::max(result, std::abs(element));
    return result;
}

uint32_t EstimateBootstrapPrecisionBits(const CryptoContext<DCRTPoly>& cc,
                                        const PublicKey<DCRTPoly>& publicKey,
                                        const PrivateKey<DCRTPoly>& secretKey,
                                        std::size_t bootstrapInputLevel) {
    const Vec4 probeValue{0.50, 0.0, 0.0, 0.0};
    Cipher probe = EncryptVector(cc, publicKey, probeValue);
    probe = ReduceToLevel(cc, probe, bootstrapInputLevel);
    const Vec4 before = Decode(cc, secretKey, probe);
    const Cipher afterCipher = cc->EvalBootstrap(probe);
    const Vec4 after = Decode(cc, secretKey, afterCipher);
    const double error = std::max(InfinityNorm(Subtract(after, before)), 1.0e-16);
    const double measuredBits = -std::log2(error);
    const int bufferedBits = static_cast<int>(std::floor(measuredBits - 5.0));
    return static_cast<uint32_t>(std::clamp(bufferedBits, 1, 40));
}

void RecordFailure(CsvFiles& files, const Options& options, const std::string& message) {
    std::string escaped = message;
    for (std::size_t position = 0;
         (position = escaped.find('"', position)) != std::string::npos; position += 2)
        escaped.insert(position, 1, '"');
    files.failure << options.blocks << ',' << gFailure.completedBlocks << ','
                  << gFailure.completedSamples << ',' << gFailure.block << ',' << gFailure.sample
                  << ',' << gFailure.stage << ',' << gFailure.operation << ',' << gFailure.branch
                  << ",\"" << escaped << '"';
    WriteMetadata(files.failure, gFailure.hasMetadata ? gFailure.metadata : CipherMetadata{});
    files.failure << '\n';
    files.failure.flush();
}
}  // namespace

int main(int argc, char** argv) {
    Options options;
    CsvFiles files;
    try {
        options = ParseOptions(argc, argv);
        const auto& selectedProfile = GetControllerProfile(options.profile);
        LoadProfile(selectedProfile);
        OpenCsvFiles(options, files);

        const std::vector<uint32_t> levelBudget{1, 1};
        const std::vector<uint32_t> bsgsDim{0, 0};
        const SecretKeyDist secretKeyDistribution = UNIFORM_TERNARY;
        const uint32_t bootstrapDepth =
            FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDistribution);
        const uint32_t bootstrapIterations = options.singleBootstrap ? 1 : kBootstrapIterations;
        const uint32_t iterativeExtraDepth = bootstrapIterations > 1 ? (bootstrapIterations - 1) : 0;
        const uint32_t preBootstrapNormalizationLevels =
            options.disableDynamicScaling ? 0 : kDynamicScalingExtraLevel;
        const uint32_t levelsAvailableAfterBootstrap =
            options.levelsAfterBootstrap + preBootstrapNormalizationLevels;
        const uint32_t depth =
            bootstrapDepth + levelsAvailableAfterBootstrap + iterativeExtraDepth;
        const uint32_t bootstrapInputLevel = depth - 1;

        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetSecretKeyDist(secretKeyDistribution);
        parameters.SetSecurityLevel(HEStd_NotSet);
        parameters.SetRingDim(4096);
        parameters.SetBatchSize(kSlots);
        parameters.SetNumLargeDigits(3);
        parameters.SetKeySwitchTechnique(HYBRID);
        parameters.SetScalingTechnique(FLEXIBLEAUTO);
        parameters.SetScalingModSize(59);  // 20 effective bits required by application; 59-bit CKKS scale leaves margin for ct-ct + bootstrap.
        parameters.SetFirstModSize(60);
        parameters.SetMultiplicativeDepth(depth);

        auto cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
        cc->Enable(FHE);
        cc->EvalBootstrapSetup(levelBudget, bsgsDim, kSlots, kCorrectionFactor);
        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);
        cc->EvalRotateKeyGen(keys.secretKey, {1, 2, 3, -1, -2, -3});
        cc->EvalBootstrapKeyGen(keys.secretKey, kSlots);

        gFailure.stage = "coefficient_encryption";
        const auto encrypted = EncryptCoefficients(cc, keys.publicKey);

        std::cout << std::setprecision(17)
                  << "OpenFHE_version=v1.5.1, source_commit="
                     "1306d14f8c26bb6150d3e6ad54f28dfe1007689e, library_version_string="
                  << GetOPENFHEVersion() << ", ring_dimension=" << cc->GetRingDimension()
                  << ", logical_slots=4, depth=" << depth
                  << ", bootstrap_depth=" << bootstrapDepth
                  << ", levels_after_bootstrap=" << levelsAvailableAfterBootstrap
                  << ", requested_postbootstrap_levels=" << options.levelsAfterBootstrap
                  << ", prebootstrap_normalization_levels="
                  << preBootstrapNormalizationLevels << '\n';
        std::cout << "experiment=ctct_uy_reconstruction, profile=" << selectedProfile.name
                  << ", label=" << selectedProfile.label << ", blocks=" << options.blocks
                  << ", rho_F=" << rho_F << ", rho_Acl=" << rho_Acl
                  << ", rank_O4=" << rank_O4 << ", cond_O4=" << cond_O4
                  << ", M_u_inf=" << M_u_inf << ", M_y_inf=" << M_y_inf
                  << ", Gamma_B=" << Gamma_B
                  << ", g_F_4=" << selectedProfile.g_F_4
                  << ", g_HF_4=" << selectedProfile.g_HF_4 << '\n';
        std::cout << "precision_fix: bootstrap_iterations=" << bootstrapIterations
                  << ", dynamic_scaling=" << (options.disableDynamicScaling ? "disabled" : "enabled")
                  << ", bootstrap_target=" << kBootstrapTarget
                  << ", application_resolution=" << kApplicationResolution
                  << ", max_bootstrap_scale=" << kMaxBootstrapScale
                  << ", scaled_bootstrap_guard=" << kScaledBootstrapGuard << '\n';
        std::cout << "profile_note: for unstable_high, the ideal u/plant trajectory matches unstable_low, "
                     "but internal state/coefficient magnitudes are intentionally larger. "
                     "Do not interpret an earlier CKKS failure as physical closed-loop instability.\n";
        std::cout << "signal_representation: each y(t+i) and u(t+i) is encrypted separately "
                     "as [value,0,0,0]; controller state remains [x1,x2,x3,x4]\n";
        std::cout << "per_block_operations: ct_ct_mult=54, relinearizations=54, "
                     "product_rescales=54, controller_coeff_ct_pt_mult=0, "
                     "public_bootstrap_normalization_ct_pt_mult<=8, rotations=36, "
                     "additions=49, production_bootstraps=4\n";
        std::cout << "fresh_bootstrap_probe="
                  << (options.freshBootstrapProbe ? "enabled" : "disabled")
                  << ", diagnostic_bootstraps_per_block="
                  << (options.freshBootstrapProbe ? 4 : 0) << '\n';
        std::cout << "encrypted_once_at_setup: P_sparse_masks=16, q_sparse_masks=16, "
                     "M_u_sparse_masks=16, M_y_sparse_masks=16\n";
        std::cout << "compiler=" << __VERSION__
                  << ", hardware_concurrency=" << std::thread::hardware_concurrency();
#ifdef _OPENMP
        std::cout << ", OpenMP_max_threads=" << omp_get_max_threads();
#endif
        std::cout << ", OMP_NUM_THREADS="
                  << (std::getenv("OMP_NUM_THREADS") ? std::getenv("OMP_NUM_THREADS") : "not set")
                  << '\n';

        { // Verify the rotation convention assumed by sparse control/reconstruction.
            const Vec4 probeState{1.0, 2.0, 3.0, 4.0};
            Cipher probe = EncryptVector(cc, keys.publicKey, probeState);
            const Vec4 left = Decode(cc, keys.secretKey, cc->EvalRotate(probe, 1));
            Cipher sparseProbe = EncryptSparseScalar(cc, keys.publicKey, 1.0);
            const Vec4 right = Decode(cc, keys.secretKey, cc->EvalRotate(sparseProbe, -1));
            if (std::abs(left[0] - 2.0) > 1e-6 || std::abs(right[1] - 1.0) > 1e-6)
                throw std::runtime_error("rotation convention smoke check failed");
            std::cout << "rotation_smoke: Rot(+1) brings x2 to slot0; Rot(-1) moves "
                         "sparse slot0 to slot1: PASS\n";
        }

        gFailure.stage = "bootstrap_warmup";
        Cipher warm = EncryptSparseScalar(cc, keys.publicKey, 0.0125);
        warm = ReduceToLevel(cc, warm, bootstrapInputLevel);
        (void)cc->EvalBootstrap(warm);

        const uint32_t iterativePrecisionBits =
            bootstrapIterations > 1
                ? EstimateBootstrapPrecisionBits(cc, keys.publicKey, keys.secretKey,
                                                 bootstrapInputLevel)
                : 0;
        const auto bootstrapScales =
            BuildBootstrapScales(options.blocks, options.disableDynamicScaling);
        std::cout << "measured_single_bootstrap_precision_bits_for_meta_bts="
                  << iterativePrecisionBits << '\n';
        std::cout << "bootstrap_scale_schedule:";
        for (std::size_t i = 0; i < std::min<std::size_t>(bootstrapScales.size(), 20); ++i)
            std::cout << (i == 0 ? " " : ",") << bootstrapScales[i];
        if (bootstrapScales.size() > 20)
            std::cout << ",...";
        std::cout << '\n';

        Vec4 nominalPlant=Xp0;
        Vec4 nominalController=Xc0;
        Vec4 encryptedPlant = nominalPlant;
        Cipher controllerCipher = EncryptVector(cc, keys.publicKey, nominalController);
        Vec4 controllerSemantic = Decode(cc, keys.secretKey, controllerCipher);

        // Plaintext controller driven by the SAME physical measurement y as the encrypted
        // loop.  Hence delta_x = Dec(c_x) - x_same_y removes plant-trajectory mismatch and
        // isolates the controller-side HE numerical error.  In exact arithmetic, between
        // block boundaries its homogeneous part obeys delta_{j+1} = F^4 delta_j.
        Vec4 sameYReferenceController = nominalController;
        const Mat4 F2 = MatrixMultiply(F, F);
        const Mat4 F4 = MatrixMultiply(F2, F2);
        const double theoreticalUnstableBlockGrowth = std::pow(rho_F, 4.0);

        bool smokeChecked = false;

        for (std::size_t block = 0; block < options.blocks; ++block) {
            const auto blockBegin = Clock::now();
            gFailure.block = block;
            gFailure.sample = -1;
            gFailure.stage = "block_start";
            gFailure.operation = "none";
            const Vec4 semanticStart = controllerSemantic;
            const Vec4 sameYReferenceStart = sameYReferenceController;
            const Vec4 decodedStart = controllerSemantic;
            const Vec4 deltaStart = Subtract(decodedStart, sameYReferenceStart);
            Vec4 semanticState = semanticStart;
            Vec4 semanticU{};
            Vec4 semanticY{};
            std::array<Cipher, 4> controlCiphertexts;
            std::array<Cipher, 4> measurementCiphertexts;
            double controlEvaluationMs = 0.0;

            for (std::size_t sample = 0; sample < 4; ++sample) {
                gFailure.sample = static_cast<long>(sample);
                gFailure.stage = "sample";
                const Vec4 nominalPlantBefore = nominalPlant;
                const Vec4 encryptedPlantBefore = encryptedPlant;
                const double nominalY = Dot(C, nominalPlant);
                const double nominalU = Dot(H, nominalController);
                const double physicalY = Dot(C, encryptedPlant);
                // Each measurement is a separate scalar ciphertext: [y(t+i),0,0,0].
                measurementCiphertexts[sample] =
                    EncryptSparseScalar(cc, keys.publicKey, physicalY);
                const Vec4 decodedY = Decode(cc, keys.secretKey, measurementCiphertexts[sample]);
                semanticY[sample] = decodedY[0];
                semanticU[sample] = Dot(H, semanticState);

                const auto evaluationBegin = Clock::now();
                // Sparse scalar control ciphertext:
                // [u(t+i),0,0,0] = sum_j Enc([P_i[j],0,0,0]) * Rot^j(c_x)
                //                    + sum_h Enc([q_ih,0,0,0]) * c_y(t+h).
                controlCiphertexts[sample] = SparseInnerEncrypted(
                    cc, controllerCipher, encrypted.pMasks[sample], files, block, sample,
                    "P_" + std::to_string(sample));
                for (std::size_t history = 0; history < sample; ++history) {
                    Cipher term = MultiplyCiphertexts(
                        cc, measurementCiphertexts[history],
                        encrypted.qMasks[sample][history], files, block, sample,
                        "q_" + std::to_string(sample) + "_" + std::to_string(history));
                    controlCiphertexts[sample] =
                        AddAligned(cc, controlCiphertexts[sample], term);
                }
                controlEvaluationMs += Milliseconds(evaluationBegin, Clock::now());

                const Vec4 decodedU = Decode(cc, keys.secretKey, controlCiphertexts[sample]);
                const double ckksU = decodedU[0];
                if (!smokeChecked) {
                    const double activeError = ActiveSlotError(decodedU, semanticU[sample]);
                    const double inactiveLeakage = InactiveSlotLeakage(decodedU);
                    std::cout << "sparse_scalar_smoke_u" << sample
                              << "_expected=" << semanticU[sample]
                              << ", active_error=" << activeError
                              << ", inactive_leakage=" << inactiveLeakage << '\n';
                    if (std::max(activeError, inactiveLeakage) > 1e-5)
                        throw std::runtime_error("sparse ct-ct scalar smoke check failed");
                }

                files.sample << block << ',' << sample << ',' << block * 4 + sample << ','
                             << nominalY << ',' << physicalY << ',' << semanticY[sample] << ','
                             << nominalU << ',' << semanticU[sample] << ',' << ckksU << ','
                             << std::abs(ckksU - semanticU[sample]) << ','
                             << std::abs(ckksU - nominalU) << ','
                             << ActiveSlotError(decodedY, physicalY) << ','
                             << InactiveSlotLeakage(decodedY) << ','
                             << ActiveSlotError(decodedU, semanticU[sample]) << ','
                             << InactiveSlotLeakage(decodedU);
                WriteVector(files.sample, nominalPlantBefore);
                WriteVector(files.sample, encryptedPlantBefore);
                WriteVector(files.sample, Subtract(encryptedPlantBefore, nominalPlantBefore));
                WriteVector(files.sample, decodedY);
                WriteVector(files.sample, decodedU);
                WriteMetadata(files.sample, GetMetadata(controlCiphertexts[sample]));
                files.sample << '\n';
                files.sample.flush();
                ++gFailure.completedSamples;

                nominalPlant = Add(MatrixVector(A, nominalPlant), Scale(B, nominalU));
                nominalController =
                    Add(MatrixVector(F, nominalController), Scale(G, nominalY));
                encryptedPlant = Add(MatrixVector(A, encryptedPlant), Scale(B, ckksU));
                semanticState =
                    Add(MatrixVector(F, semanticState), Scale(G, semanticY[sample]));
                sameYReferenceController =
                    Add(MatrixVector(F, sameYReferenceController), Scale(G, physicalY));
            }
            smokeChecked = true;

            std::array<Cipher, 4> bootstrapped;
            double bootstrapMs = 0.0;
            const double bootstrapScale = bootstrapScales[block];
            for (std::size_t sample = 0; sample < 4; ++sample) {
                gFailure.sample = static_cast<long>(sample);
                gFailure.stage = "bootstrap";
                gFailure.operation = "EvalBootstrap_u" + std::to_string(sample);

                const CipherMetadata before = GetMetadata(controlCiphertexts[sample]);
                const Vec4 rawPre = Decode(cc, keys.secretKey, controlCiphertexts[sample]);

                // Public precision normalization. The controller coefficient path remains ct-ct;
                // this scalar multiplication only uses public block-dependent normalization.
                Cipher scaled = controlCiphertexts[sample];
                if (bootstrapScale != 1.0) {
                    scaled = cc->EvalMult(scaled, bootstrapScale);
                    scaled = NormalizeScaleDegree(cc, scaled);
                }
                Cipher input = ReduceToLevel(cc, scaled, bootstrapInputLevel);
                const CipherMetadata immediatelyBefore = GetMetadata(input);
                const Vec4 preScaled = Decode(cc, keys.secretKey, input);
                const double scaledInputMax = MaxAbs(preScaled);
                if (scaledInputMax >= kScaledBootstrapGuard)
                    throw std::runtime_error(
                        "scaled bootstrap input exceeded guard: max_abs=" +
                        std::to_string(scaledInputMax) + ", scale=" +
                        std::to_string(bootstrapScale));

                const auto begin = Clock::now();
                Cipher bootScaled =
                    bootstrapIterations > 1
                        ? cc->EvalBootstrap(input, bootstrapIterations, iterativePrecisionBits)
                        : cc->EvalBootstrap(input);
                const double elapsed = Milliseconds(begin, Clock::now());
                bootstrapMs += elapsed;

                // Keep the bootstrapped controls in the scaled domain. Reconstruction applies
                // the inverse public scale to the encrypted M_u columns, avoiding an extra
                // post-bootstrap plaintext multiply on the recursive signal ciphertext.
                bootstrapped[sample] = bootScaled;
                const CipherMetadata after = GetMetadata(bootstrapped[sample]);
                gFailure.metadata = after;
                gFailure.hasMetadata = true;
                const Vec4 postScaled = Decode(cc, keys.secretKey, bootstrapped[sample]);

                Vec4 preEquivalent{};
                Vec4 post{};
                for (std::size_t slot = 0; slot < 4; ++slot) {
                    preEquivalent[slot] = preScaled[slot] / bootstrapScale;
                    post[slot] = postScaled[slot] / bootstrapScale;
                }
                const double productionBootstrapError =
                    InfinityNorm(Subtract(post, preEquivalent));

                double freshBootstrapError = std::numeric_limits<double>::quiet_NaN();
                double freshBootstrapMs = 0.0;
                if (options.freshBootstrapProbe) {
                    gFailure.branch = "fresh_bootstrap_probe";
                    gFailure.stage = "fresh_bootstrap_probe";
                    gFailure.operation = "EvalBootstrap_fresh_u" + std::to_string(sample);
                    Cipher fresh = EncryptSparseScalar(cc, keys.publicKey, rawPre[0]);
                    if (bootstrapScale != 1.0) {
                        fresh = cc->EvalMult(fresh, bootstrapScale);
                        fresh = NormalizeScaleDegree(cc, fresh);
                    }
                    fresh = ReduceToLevel(cc, fresh, bootstrapInputLevel);
                    const Vec4 freshPreScaled = Decode(cc, keys.secretKey, fresh);
                    const auto freshBegin = Clock::now();
                    Cipher freshBootScaled =
                        bootstrapIterations > 1
                            ? cc->EvalBootstrap(fresh, bootstrapIterations, iterativePrecisionBits)
                            : cc->EvalBootstrap(fresh);
                    freshBootstrapMs = Milliseconds(freshBegin, Clock::now());
                    const Vec4 freshPostScaled = Decode(cc, keys.secretKey, freshBootScaled);
                    Vec4 freshPre{};
                    Vec4 freshPost{};
                    for (std::size_t slot = 0; slot < 4; ++slot) {
                        freshPre[slot] = freshPreScaled[slot] / bootstrapScale;
                        freshPost[slot] = freshPostScaled[slot] / bootstrapScale;
                    }
                    freshBootstrapError = InfinityNorm(Subtract(freshPost, freshPre));
                    gFailure.branch = "production";
                    gFailure.stage = "bootstrap";
                }

                files.bootstrap << block << ',' << sample << ',' << semanticU[sample] << ','
                                << bootstrapScale << ',' << scaledInputMax << ','
                                << bootstrapIterations << ',' << iterativePrecisionBits;
                WriteVector(files.bootstrap, preEquivalent);
                WriteVector(files.bootstrap, post);
                files.bootstrap << ',' << ActiveSlotError(preEquivalent, semanticU[sample]) << ','
                                << InactiveSlotLeakage(preEquivalent) << ','
                                << productionBootstrapError << ','
                                << ActiveSlotError(post, semanticU[sample]) << ','
                                << InactiveSlotLeakage(post);
                WriteMetadata(files.bootstrap, before);
                WriteMetadata(files.bootstrap, immediatelyBefore);
                WriteMetadata(files.bootstrap, after);
                files.bootstrap << ',' << elapsed << ',' << options.freshBootstrapProbe << ','
                                << productionBootstrapError << ',' << freshBootstrapError << ','
                                << freshBootstrapMs << '\n';
                files.bootstrap.flush();
            }

            gFailure.sample = 4;
            gFailure.stage = "reconstruction";
            const auto reconstructionBegin = Clock::now();
            Cipher reconstructedU;
            Cipher reconstructedY;
            for (std::size_t column = 0; column < 4; ++column) {
                // u(t+column) and y(t+column) are separate sparse scalar ciphertexts.
                // Scatter each active scalar slot into the four controller-state slots with
                // one-hot encrypted reconstruction coefficients.
                Cipher uTerm = ScatterEncrypted(
                    cc, bootstrapped[column], encrypted.muMasks[column], files, block, 4,
                    "M_u_col_" + std::to_string(column), 1.0 / bootstrapScale);
                Cipher yTerm = ScatterEncrypted(
                    cc, measurementCiphertexts[column], encrypted.myMasks[column], files, block,
                    4, "M_y_col_" + std::to_string(column));
                reconstructedU = reconstructedU ? AddAligned(cc, reconstructedU, uTerm) : uTerm;
                reconstructedY = reconstructedY ? AddAligned(cc, reconstructedY, yTerm) : yTerm;
            }
            controllerCipher = AddAligned(cc, reconstructedU, reconstructedY);
            const double reconstructionMs =
                Milliseconds(reconstructionBegin, Clock::now());
            const CipherMetadata controllerMetadata = GetMetadata(controllerCipher);
            const Vec4 decodedController = Decode(cc, keys.secretKey, controllerCipher);
            controllerSemantic = decodedController;
            const Vec4 exactSemanticReconstruction =
                Add(MatrixVector(M_u, semanticU), MatrixVector(M_y, semanticY));
            const Vec4 controllerError = Subtract(decodedController, nominalController);
            const Vec4 plantError = Subtract(encryptedPlant, nominalPlant);

            // Diagnostic only: same-y controller numerical mode.  This is NOT the success
            // metric; u and plant trajectory remain the primary performance metrics.
            const Vec4 deltaEnd = Subtract(decodedController, sameYReferenceController);
            const Vec4 homogeneousPrediction = MatrixVector(F4, deltaStart);
            const Vec4 innovation = Subtract(deltaEnd, homogeneousPrediction);
            const double deltaStartInf = InfinityNorm(deltaStart);
            const double deltaEndInf = InfinityNorm(deltaEnd);
            const double blockGrowthInf =
                deltaStartInf > 1e-30 ? deltaEndInf / deltaStartInf
                                      : std::numeric_limits<double>::quiet_NaN();
            const double unstableStartAbs = std::abs(deltaStart[0]);
            const double unstableEndAbs = std::abs(deltaEnd[0]);
            const double unstableGrowthAbs =
                unstableStartAbs > 1e-30 ? unstableEndAbs / unstableStartAbs
                                         : std::numeric_limits<double>::quiet_NaN();

            files.mode << block;
            WriteVector(files.mode, sameYReferenceStart);
            WriteVector(files.mode, decodedStart);
            WriteVector(files.mode, deltaStart);
            WriteVector(files.mode, sameYReferenceController);
            WriteVector(files.mode, decodedController);
            WriteVector(files.mode, deltaEnd);
            WriteVector(files.mode, homogeneousPrediction);
            WriteVector(files.mode, innovation);
            files.mode << ',' << deltaStartInf << ',' << deltaEndInf << ',' << blockGrowthInf
                       << ',' << unstableStartAbs << ',' << unstableEndAbs << ','
                       << unstableGrowthAbs << ',' << theoreticalUnstableBlockGrowth << ','
                       << InfinityNorm(innovation) << '\n';
            files.mode.flush();

            files.block << block;
            WriteVector(files.block, nominalController);
            WriteVector(files.block, exactSemanticReconstruction);
            WriteVector(files.block, decodedController);
            WriteVector(files.block, controllerError);
            WriteVector(files.block, plantError);
            files.block << ',' << InfinityNorm(controllerError) << ',' << InfinityNorm(plantError)
                        << ',' << InfinityNorm(Subtract(exactSemanticReconstruction, semanticState));
            WriteMetadata(files.block, controllerMetadata);
            files.block << '\n';
            files.block.flush();

            files.timing << block << ',' << controlEvaluationMs << ',' << bootstrapMs << ','
                         << reconstructionMs << ',' << Milliseconds(blockBegin, Clock::now())
                         << '\n';
            files.timing.flush();
            ++gFailure.completedBlocks;
            std::cout << "block=" << block
                      << ", control_error_inf=" << InfinityNorm(controllerError)
                      << ", plant_error_inf=" << InfinityNorm(plantError)
                      << ", semantic_identity_inf="
                      << InfinityNorm(Subtract(exactSemanticReconstruction, semanticState))
                      << ", bootstrap_scale=" << bootstrapScale
                      << ", same_y_delta_inf=" << deltaEndInf
                      << ", unstable_delta_abs=" << unstableEndAbs
                      << ", unstable_block_growth=" << unstableGrowthAbs
                      << ", expected_F4_growth=" << theoreticalUnstableBlockGrowth
                      << ", u4_bootstrap_ms=" << bootstrapMs
                      << ", x_level=" << controllerMetadata.level
                      << ", x_towers=" << controllerMetadata.towers << '\n';
        }

        std::cout << "completed_blocks=" << gFailure.completedBlocks << '\n';
        return 0;
    }
    catch (const std::exception& error) {
        if (files.failure.is_open())
            RecordFailure(files, options, error.what());
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
