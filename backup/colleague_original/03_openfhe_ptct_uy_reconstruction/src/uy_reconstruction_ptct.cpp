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
double rho_F=0.0, rho_Acl=0.0, cond_O4=0.0;

void LoadProfile(const ControllerProfile& p) {
    A=p.A; B=p.B; C=p.C; F=p.F; G=p.G; H=p.H; Xp0=p.xp0; Xc0=p.xc0;
    O4=p.O4; T4=p.T4; M_u=p.M_u; M_y=p.M_y; P=p.P; q=p.q;
    rho_F=p.rho_F; rho_Acl=p.rho_Acl; cond_O4=p.cond_O4;
}

constexpr uint32_t kSlots = 4;
constexpr uint32_t kCorrectionFactor = 13;
constexpr uint32_t kLevelsAfterBootstrap = 3;  // two circuit levels plus one spare
using Cipher = Ciphertext<DCRTPoly>;
using Clock = std::chrono::steady_clock;
using Vec8 = std::array<double, 8>;
using Mat8 = std::array<Vec8, 8>;

struct Options {
    std::string mode = "ckks_u4_bootstrap";
    std::string profile = "stable";
    std::string coeffMode = "plaintext";
    std::size_t blocks = 3;
    std::string resultsDir = "results";
    bool clearImag = false;
};
struct CipherMetadata {
    std::size_t level = 0, towers = 0, noiseScaleDegree = 0;
    double scalingFactor = 0.0;
};
struct FailureState {
    std::size_t block = 0, completedBlocks = 0, completedSamples = 0;
    long sample = -1;
    std::string stage = "setup";
    CipherMetadata metadata{};
    bool hasMetadata = false;
} gFailure;

Options ParseOptions(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--mode" && i + 1 < argc) o.mode = argv[++i];
        else if (a == "--profile" && i + 1 < argc) o.profile = argv[++i];
        else if (a == "--coeff-mode" && i + 1 < argc) o.coeffMode = argv[++i];
        else if (a == "--blocks" && i + 1 < argc) o.blocks = std::stoul(argv[++i]);
        else if (a == "--results-dir" && i + 1 < argc) o.resultsDir = argv[++i];
        else if (a == "--clear-imag") o.clearImag = true;
        else throw std::invalid_argument("usage: uy_reconstruction_ptct [--profile stable|unstable_low|unstable_high] [--mode plaintext_reference|ckks_no_bootstrap|ckks_u4_bootstrap] [--blocks N] [--results-dir PATH] [--clear-imag]");
    }
    if (o.mode != "plaintext_reference" && o.mode != "ckks_no_bootstrap" && o.mode != "ckks_u4_bootstrap")
        throw std::invalid_argument("invalid --mode");
    if (o.coeffMode != "plaintext")
        throw std::invalid_argument("encrypted coefficient mode is deferred until the plaintext-coefficient long run passes");
    if (!o.blocks) throw std::invalid_argument("--blocks must be positive");
    return o;
}

double Ms(Clock::time_point a, Clock::time_point b) { return std::chrono::duration<double, std::milli>(b-a).count(); }
Vec4 Add(const Vec4& a, const Vec4& b) { Vec4 r{}; for(size_t i=0;i<4;++i) r[i]=a[i]+b[i]; return r; }
Vec4 Sub(const Vec4& a, const Vec4& b) { Vec4 r{}; for(size_t i=0;i<4;++i) r[i]=a[i]-b[i]; return r; }
Vec4 Scale(const Vec4& a,double s) { Vec4 r{}; for(size_t i=0;i<4;++i) r[i]=a[i]*s; return r; }
Vec4 Mv(const Mat4& m,const Vec4& v) { Vec4 r{}; for(size_t i=0;i<4;++i) for(size_t j=0;j<4;++j) r[i]+=m[i][j]*v[j]; return r; }
double Dot(const Vec4& a,const Vec4& b) { double r=0; for(size_t i=0;i<4;++i) r+=a[i]*b[i]; return r; }
double Inf(const Vec4& v) { double r=0; for(double x:v) r=std::max(r,std::abs(x)); return r; }
double ActiveError(const Vec4& v,double expected) { return std::abs(v[0]-expected); }
double InactiveLeakage(const Vec4& v) { double r=0; for(size_t i=1;i<4;++i) r=std::max(r,std::abs(v[i])); return r; }
double SparseErrorInf(const Vec4& v,double expected) { return std::max(ActiveError(v,expected),InactiveLeakage(v)); }
double Inf8(const Vec8& v) { double r=0; for(double x:v) r=std::max(r,std::abs(x)); return r; }
Vec8 Join(const Vec4&a,const Vec4&b) { Vec8 r{}; std::copy(a.begin(),a.end(),r.begin()); std::copy(b.begin(),b.end(),r.begin()+4); return r; }
Vec8 Add8(const Vec8&a,const Vec8&b) { Vec8 r{}; for(size_t i=0;i<8;++i)r[i]=a[i]+b[i]; return r; }
Vec8 Sub8(const Vec8&a,const Vec8&b) { Vec8 r{}; for(size_t i=0;i<8;++i)r[i]=a[i]-b[i]; return r; }
Vec8 Mv8(const Mat8&m,const Vec8&v) { Vec8 r{}; for(size_t i=0;i<8;++i)for(size_t j=0;j<8;++j)r[i]+=m[i][j]*v[j]; return r; }
Mat8 Mm8(const Mat8&a,const Mat8&b) { Mat8 r{}; for(size_t i=0;i<8;++i)for(size_t k=0;k<8;++k)for(size_t j=0;j<8;++j)r[i][j]+=a[i][k]*b[k][j]; return r; }
Mat8 Eye8(){ Mat8 r{}; for(size_t i=0;i<8;++i)r[i][i]=1; return r; }
Mat8 Pow8(const Mat8&a,size_t n){ Mat8 r=Eye8(); for(size_t i=0;i<n;++i)r=Mm8(r,a); return r; }
Mat8 Acl(){ Mat8 r{}; for(size_t i=0;i<4;++i)for(size_t j=0;j<4;++j){r[i][j]=A[i][j];r[i][j+4]=B[i]*H[j];r[i+4][j]=G[i]*C[j];r[i+4][j+4]=F[i][j];}return r; }

CipherMetadata Meta(const Cipher& c){ return {c->GetLevel(),c->GetElements()[0].GetNumOfElements(),c->GetNoiseScaleDeg(),c->GetScalingFactor()}; }
Vec4 Decode(const CryptoContext<DCRTPoly>&cc,const PrivateKey<DCRTPoly>&sk,const Cipher&c){ Plaintext p;cc->Decrypt(sk,c,&p);p->SetLength(4);const auto v=p->GetRealPackedValue();return {v[0],v[1],v[2],v[3]}; }
Cipher EncryptSparseScalar(const CryptoContext<DCRTPoly>&cc,const PublicKey<DCRTPoly>&pk,double x){ const std::vector<double>v{x,0.0,0.0,0.0};return cc->Encrypt(pk,cc->MakeCKKSPackedPlaintext(v,1,0,nullptr,4)); }
Cipher EncryptVec(const CryptoContext<DCRTPoly>&cc,const PublicKey<DCRTPoly>&pk,const Vec4&x){ const std::vector<double>v{x[0],x[1],x[2],x[3]};return cc->Encrypt(pk,cc->MakeCKKSPackedPlaintext(v,1,0,nullptr,4)); }
Cipher MulPlain(const CryptoContext<DCRTPoly>&cc,const Cipher&c,const Vec4&v){ const std::vector<double>values{v[0],v[1],v[2],v[3]};auto p=cc->MakeCKKSPackedPlaintext(values,c->GetNoiseScaleDeg(),c->GetLevel(),nullptr,4);auto product=cc->EvalMult(c,p);return cc->GetScheme()->ModReduceInternal(product,1); }
Cipher AddAligned(const CryptoContext<DCRTPoly>&cc,Cipher lhs,Cipher rhs){ const auto target=std::max(lhs->GetLevel(),rhs->GetLevel()); if(lhs->GetLevel()<target)lhs=cc->GetScheme()->LevelReduceInternal(lhs,target-lhs->GetLevel()); if(rhs->GetLevel()<target)rhs=cc->GetScheme()->LevelReduceInternal(rhs,target-rhs->GetLevel()); return cc->EvalAdd(lhs,rhs); }
Cipher SparseInnerPlain(const CryptoContext<DCRTPoly>&cc,const Cipher&state,const Vec4&row){
 Cipher sum;
 for(size_t j=0;j<4;++j){
  Cipher rotated=(j==0)?state:cc->EvalRotate(state,static_cast<int32_t>(j));
  Vec4 mask{}; mask[0]=row[j];
  Cipher term=MulPlain(cc,rotated,mask);
  sum=sum?AddAligned(cc,sum,term):term;
 }
 return sum;
}
Cipher SparseScalePlain(const CryptoContext<DCRTPoly>&cc,const Cipher&scalar,double coefficient){return MulPlain(cc,scalar,{coefficient,0.0,0.0,0.0});}
Cipher ScatterPlain(const CryptoContext<DCRTPoly>&cc,const Cipher&scalar,const Vec4&column){
 Cipher sum;
 for(size_t row=0;row<4;++row){
  Cipher placed=(row==0)?scalar:cc->EvalRotate(scalar,-static_cast<int32_t>(row));
  Vec4 mask{}; mask[row]=column[row];
  Cipher term=MulPlain(cc,placed,mask);
  sum=sum?AddAligned(cc,sum,term):term;
 }
 return sum;
}
// 2*Re(c) = c + conj(c).  Plant only sees the real part of the decoded u, so any
// imaginary component fed back into the reconstruction evolves open-loop through F^K.
// Uses the conjugation key (index M-1) generated by EvalBootstrapKeyGen.  The 1/2 is
// folded into the M_u scatter masks so that no extra level is consumed.
Cipher DoubleRealPart(const CryptoContext<DCRTPoly>&cc,const Cipher&c){
 auto conj=cc->EvalAutomorphism(c,cc->GetCyclotomicOrder()-1,cc->GetEvalAutomorphismKeyMap(c->GetKeyTag()));
 return cc->EvalAdd(c,conj);
}
Cipher Reduce(const CryptoContext<DCRTPoly>&cc,Cipher c,uint32_t level){if(c->GetLevel()<level)c=cc->GetScheme()->LevelReduceInternal(c,level-c->GetLevel());return c;}

void H4(std::ostream&o,const std::string&n){for(size_t i=0;i<4;++i)o<<','<<n<<'_'<<i;}
void H8(std::ostream&o,const std::string&n){for(size_t i=0;i<8;++i)o<<','<<n<<'_'<<i;}
template<size_t N>void W(std::ostream&o,const std::array<double,N>&v){for(double x:v)o<<','<<x;}
void HM(std::ostream&o,const std::string&n){o<<','<<n<<"_level,"<<n<<"_towers,"<<n<<"_noiseScaleDegree,"<<n<<"_scalingFactor";}
void WM(std::ostream&o,const CipherMetadata&m){o<<','<<m.level<<','<<m.towers<<','<<m.noiseScaleDegree<<','<<m.scalingFactor;}
void Print4(const std::string&n,const Vec4&v){std::cout<<n<<"=[";for(size_t i=0;i<4;++i)std::cout<<(i?", ":"")<<v[i];std::cout<<"]\n";}

double PlainReference(size_t samples){Vec4 xp=Xp0,xc=Xc0;double e=0;for(size_t b=0;b<samples/4;++b){Vec4 u{},y{};for(size_t i=0;i<4;++i){y[i]=Dot(C,xp);u[i]=Dot(H,xc);xp=Add(Mv(A,xp),Scale(B,u[i]));xc=Add(Mv(F,xc),Scale(G,y[i]));}e=std::max(e,Inf(Sub(Add(Mv(M_u,u),Mv(M_y,y)),xc)));}return e;}

struct Csv { std::ofstream sample,block,boot,time,failure; };
void OpenCsv(const Options&o,Csv&f){
 std::filesystem::create_directories(o.resultsDir);f.sample.open(o.resultsDir+"/sample_trace.csv");f.block.open(o.resultsDir+"/block_metrics.csv");f.boot.open(o.resultsDir+"/bootstrap_diagnostics.csv");f.time.open(o.resultsDir+"/timing.csv");f.failure.open(o.resultsDir+"/failure.csv");
 for(auto*s:{&f.sample,&f.block,&f.boot,&f.time,&f.failure})*s<<std::scientific<<std::setprecision(17);
 f.sample<<"block,within_block_sample,global_sample,y_plain,y_physical,y_semantic,epsilon_y,u_plain,u_semantic_exact,u_ckks,epsilon_u,epsilon_y_abs,epsilon_u_abs,y_active_error,y_inactive_leakage,u_active_error,u_inactive_leakage";for(auto n:{"plant_plain","plant_ckks","plant_error","y_decoded_slots","u_decoded_slots"})H4(f.sample,n);HM(f.sample,"y");HM(f.sample,"u");f.sample<<'\n';
 f.block<<"block";for(auto n:{"start_state","z_end","x_semantic_next","x_reconstructed","epsilon_x"})H4(f.block,n);for(auto n:{"E_start","E_end","Acl4_E_start","W_u","W_y","W_x","W_total","closure"})H8(f.block,n);f.block<<",epsilon_x_inf,E_start_inf,E_end_inf,W_u_inf,W_y_inf,W_x_inf,W_total_inf,closure_inf,semantic_identity_inf";HM(f.block,"x_next");f.block<<'\n';
 f.boot<<"block,within_block_sample,u_semantic";H4(f.boot,"u_pre");H4(f.boot,"u_post");f.boot<<",u_active_error_pre,u_inactive_leakage_pre,e_u_boot_inf,u_active_error_post,u_inactive_leakage_post";HM(f.boot,"before_reduction");HM(f.boot,"immediately_before_bootstrap");HM(f.boot,"after_bootstrap");f.boot<<",bootstrap_ms\n";
 f.time<<"block,u_evaluation_ms,u4_bootstrap_ms,u0_bootstrap_ms,u1_bootstrap_ms,u2_bootstrap_ms,u3_bootstrap_ms,reconstruction_ms,full_block_ms\n";
 f.failure<<"requested_mode,requested_blocks,completed_blocks,completed_samples,block,sample,stage,message";HM(f.failure,"last");f.failure<<'\n';
}
void Fail(Csv&f,const Options&o,const std::string&m){std::string e=m;for(size_t p=0;(p=e.find('"',p))!=std::string::npos;p+=2)e.insert(p,1,'"');f.failure<<o.mode<<','<<o.blocks<<','<<gFailure.completedBlocks<<','<<gFailure.completedSamples<<','<<gFailure.block<<','<<gFailure.sample<<','<<gFailure.stage<<",\""<<e<<'"';WM(f.failure,gFailure.hasMetadata?gFailure.metadata:CipherMetadata{});f.failure<<'\n';f.failure.flush();}
}

int main(int argc,char**argv){
 Options o;Csv csv;
 try{
  o=ParseOptions(argc,argv);
  const auto& selectedProfile = GetControllerProfile(o.profile);
  LoadProfile(selectedProfile);
  if(o.mode=="plaintext_reference"){double e=PlainReference(std::max<size_t>(1000,o.blocks*4));std::cout<<std::scientific<<std::setprecision(17)<<"plaintext_reconstruction_max_abs_error="<<e<<'\n';return e<1e-11?0:1;}
  OpenCsv(o,csv);
  const std::vector<uint32_t>budget{1,1},bsgs{0,0};const SecretKeyDist dist=UNIFORM_TERNARY;
  uint32_t bootDepth=FHECKKSRNS::GetBootstrapDepth(budget,dist),depth=bootDepth+kLevelsAfterBootstrap,inputLevel=depth-1;
  CCParams<CryptoContextCKKSRNS>p;p.SetSecretKeyDist(dist);p.SetSecurityLevel(HEStd_NotSet);p.SetRingDim(4096);p.SetBatchSize(4);p.SetNumLargeDigits(3);p.SetKeySwitchTechnique(HYBRID);p.SetScalingTechnique(FLEXIBLEAUTO);p.SetScalingModSize(59);p.SetFirstModSize(60);p.SetMultiplicativeDepth(depth);
  auto cc=GenCryptoContext(p);cc->Enable(PKE);cc->Enable(KEYSWITCH);cc->Enable(LEVELEDSHE);cc->Enable(ADVANCEDSHE);cc->Enable(FHE);
  bool doBoot=o.mode=="ckks_u4_bootstrap";if(doBoot)cc->EvalBootstrapSetup(budget,bsgs,4,kCorrectionFactor);auto kp=cc->KeyGen();cc->EvalMultKeyGen(kp.secretKey);cc->EvalRotateKeyGen(kp.secretKey,{1,2,3,-1,-2,-3});if(doBoot)cc->EvalBootstrapKeyGen(kp.secretKey,4);
  std::cout<<std::setprecision(17)<<"OpenFHE_version=v1.5.1, source_commit=1306d14f8c26bb6150d3e6ad54f28dfe1007689e, library_version_string="<<GetOPENFHEVersion()<<", ring_dimension="<<cc->GetRingDimension()<<", logical_slots=4, depth="<<depth<<", bootstrap_depth="<<bootDepth<<", levels_after_bootstrap="<<kLevelsAfterBootstrap<<'\n';
  std::cout<<"profile="<<selectedProfile.name<<", label="<<selectedProfile.label<<", rho_F="<<rho_F<<", rho_Acl="<<rho_Acl<<", cond_O4="<<cond_O4<<"\n";
  std::cout<<"mode="<<o.mode<<", coeff_mode="<<o.coeffMode<<", blocks="<<o.blocks<<", clear_imag="<<(o.clearImag?"on":"off")<<", scaling=FLEXIBLEAUTO, scaling_modulus_size=59, first_modulus_size=60, security=HEStd_NotSet\n";
  std::cout<<"signal_representation: each y(t+i) and u(t+i) is a separate sparse scalar ciphertext [value,0,0,0]\n";
  std::cout<<"per_block_bootstraps: u="<<(doBoot?4:0)<<", x=0, y=0, r=0\nper_block_baseline_operations: plaintext_mult=54, rotations=36, additions=49, ct_ct_mult=0, relinearizations=0, explicit_rescales=54, bootstrap_input_normalization_calls="<<(doBoot?4:0)<<"\n";
  std::cout<<"compiler="<<__VERSION__<<", hardware_concurrency="<<std::thread::hardware_concurrency();
#ifdef _OPENMP
  std::cout<<", OpenMP_max_threads="<<omp_get_max_threads();
#endif
  std::cout<<", OMP_NUM_THREADS="<<(std::getenv("OMP_NUM_THREADS")?std::getenv("OMP_NUM_THREADS"):"not set")<<'\n';
  { // Verify the OpenFHE rotation convention used by the sparse layout.
   const Vec4 probeState{1.0,2.0,3.0,4.0};
   auto probe=EncryptVec(cc,kp.publicKey,probeState);
   auto left=Decode(cc,kp.secretKey,cc->EvalRotate(probe,1));
   auto sparse=EncryptSparseScalar(cc,kp.publicKey,1.0);
   auto right=Decode(cc,kp.secretKey,cc->EvalRotate(sparse,-1));
   if(std::abs(left[0]-2.0)>1e-6 || std::abs(right[1]-1.0)>1e-6)
    throw std::runtime_error("rotation convention smoke check failed");
   std::cout<<"rotation_smoke: Rot(+1) brings x2 to slot0; Rot(-1) moves sparse slot0 to slot1: PASS\n";
  }
  if(doBoot){gFailure.stage="bootstrap_warmup";auto c=Reduce(cc,EncryptSparseScalar(cc,kp.publicKey,.0125),inputLevel);auto d=cc->EvalBootstrap(c);(void)d;}

  Vec4 xp=Xp0,xc=Xc0,xpC=xp;auto anchor=EncryptVec(cc,kp.publicKey,xc);Vec4 anchorSem=Decode(cc,kp.secretKey,anchor);Mat8 acl=Acl(),acl4=Pow8(acl,4);bool smoke=false;
  for(size_t b=0;b<o.blocks;++b){
   auto blockBegin=Clock::now();gFailure.block=b;gFailure.sample=-1;gFailure.stage="block_start";double uMs=0,reconMs=0;std::array<double,4>bootMs{};Vec4 start=anchorSem,z=start,uSem{},ySem{};Vec8 eStart=Join(Sub(xpC,xp),Sub(start,xc)),wU{},wY{};std::array<Cipher,4>cu,cy;std::array<Vec4,4>du{};
   for(size_t i=0;i<4;++i){
    gFailure.sample=i;gFailure.stage="sample";Vec4 xp0=xp,xpC0=xpC;double yp=Dot(C,xp),up=Dot(H,xc),yPhysical=Dot(C,xpC);cy[i]=EncryptSparseScalar(cc,kp.publicKey,yPhysical);auto ym=Meta(cy[i]);Vec4 dy=Decode(cc,kp.secretKey,cy[i]);ySem[i]=dy[0];double ey=ySem[i]-yPhysical;uSem[i]=Dot(H,z);
    auto t=Clock::now();cu[i]=SparseInnerPlain(cc,anchor,P[i]);for(size_t j=0;j<i;++j)cu[i]=AddAligned(cc,cu[i],SparseScalePlain(cc,cy[j],q[i][j]));uMs+=Ms(t,Clock::now());auto um=Meta(cu[i]);du[i]=Decode(cc,kp.secretKey,cu[i]);double ua=du[i][0],eu=ua-uSem[i];
    if(!smoke){double ae=ActiveError(du[i],uSem[i]),leak=InactiveLeakage(du[i]);Print4("sparse_scalar_smoke_u"+std::to_string(i),du[i]);std::cout<<"sparse_scalar_smoke_u"<<i<<"_expected="<<uSem[i]<<", active_error="<<ae<<", inactive_leakage="<<leak<<'\n';if(std::max(ae,leak)>1e-5)throw std::runtime_error("sparse scalar u smoke diagnostic failed");}
    wU=Add8(wU,Mv8(Pow8(acl,3-i),Join(Scale(B,eu),Vec4{})));wY=Add8(wY,Mv8(Pow8(acl,3-i),Join(Vec4{},Scale(G,ey))));
    csv.sample<<b<<','<<i<<','<<b*4+i<<','<<yp<<','<<yPhysical<<','<<ySem[i]<<','<<ey<<','<<up<<','<<uSem[i]<<','<<ua<<','<<eu<<','<<std::abs(ey)<<','<<std::abs(eu)<<','<<ActiveError(dy,yPhysical)<<','<<InactiveLeakage(dy)<<','<<ActiveError(du[i],uSem[i])<<','<<InactiveLeakage(du[i]);W(csv.sample,xp0);W(csv.sample,xpC0);W(csv.sample,Sub(xpC0,xp0));W(csv.sample,dy);W(csv.sample,du[i]);WM(csv.sample,ym);WM(csv.sample,um);csv.sample<<'\n';csv.sample.flush();++gFailure.completedSamples;
    xp=Add(Mv(A,xp),Scale(B,up));xc=Add(Mv(F,xc),Scale(G,yp));xpC=Add(Mv(A,xpC),Scale(B,ua));z=Add(Mv(F,z),Scale(G,ySem[i]));
   } smoke=true;
   std::array<Cipher,4>bu;
   for(size_t i=0;i<4;++i){gFailure.sample=i;gFailure.stage="u_bootstrap_"+std::to_string(i);auto before=Meta(cu[i]);Cipher in=doBoot?Reduce(cc,cu[i],inputLevel):cu[i];auto immediate=Meta(in);Vec4 pre=Decode(cc,kp.secretKey,in);if(doBoot){auto t=Clock::now();bu[i]=cc->EvalBootstrap(in);if(o.clearImag)bu[i]=DoubleRealPart(cc,bu[i]);bootMs[i]=Ms(t,Clock::now());}else bu[i]=in;auto after=Meta(bu[i]);gFailure.metadata=after;gFailure.hasMetadata=true;Vec4 post=Decode(cc,kp.secretKey,bu[i]);if(o.clearImag)post=Scale(post,0.5);csv.boot<<b<<','<<i<<','<<uSem[i];W(csv.boot,pre);W(csv.boot,post);csv.boot<<','<<ActiveError(pre,uSem[i])<<','<<InactiveLeakage(pre)<<','<<Inf(Sub(post,pre))<<','<<ActiveError(post,uSem[i])<<','<<InactiveLeakage(post);WM(csv.boot,before);WM(csv.boot,immediate);WM(csv.boot,after);csv.boot<<','<<bootMs[i]<<'\n';csv.boot.flush();}
   gFailure.sample=4;gFailure.stage="column_reconstruction";auto rt=Clock::now();Cipher xu;
   for(size_t j=0;j<4;++j){Vec4 col{};for(size_t i=0;i<4;++i)col[i]=M_u[i][j]*(o.clearImag?0.5:1.0);auto term=ScatterPlain(cc,bu[j],col);xu=xu?AddAligned(cc,xu,term):term;}
   Cipher xy;for(size_t j=0;j<4;++j){Vec4 col{};for(size_t i=0;i<4;++i)col[i]=M_y[i][j];auto term=ScatterPlain(cc,cy[j],col);xy=xy?AddAligned(cc,xy,term):term;}anchor=AddAligned(cc,xu,xy);reconMs=Ms(rt,Clock::now());auto xm=Meta(anchor);Vec4 decoded=Decode(cc,kp.secretKey,anchor);anchorSem=decoded;Vec4 xSemantic=Add(Mv(M_u,uSem),Mv(M_y,ySem)),epsX=Sub(decoded,xSemantic);Vec8 eEnd=Join(Sub(xpC,xp),Sub(decoded,xc)),prop=Mv8(acl4,eStart),wX=Join(Vec4{},epsX),wTot=Add8(Add8(wU,wY),wX),closure=Sub8(Sub8(eEnd,prop),wTot);
   csv.block<<b;W(csv.block,start);W(csv.block,z);W(csv.block,xSemantic);W(csv.block,decoded);W(csv.block,epsX);W(csv.block,eStart);W(csv.block,eEnd);W(csv.block,prop);W(csv.block,wU);W(csv.block,wY);W(csv.block,wX);W(csv.block,wTot);W(csv.block,closure);csv.block<<','<<Inf(epsX)<<','<<Inf8(eStart)<<','<<Inf8(eEnd)<<','<<Inf8(wU)<<','<<Inf8(wY)<<','<<Inf8(wX)<<','<<Inf8(wTot)<<','<<Inf8(closure)<<','<<Inf(Sub(xSemantic,z));WM(csv.block,xm);csv.block<<'\n';csv.block.flush();
   double bt=bootMs[0]+bootMs[1]+bootMs[2]+bootMs[3];csv.time<<b<<','<<uMs<<','<<bt;for(double x:bootMs)csv.time<<','<<x;csv.time<<','<<reconMs<<','<<Ms(blockBegin,Clock::now())<<'\n';csv.time.flush();++gFailure.completedBlocks;
   std::cout<<"block="<<b<<", epsilon_x_inf="<<Inf(epsX)<<", closure_inf="<<Inf8(closure)<<", semantic_identity_inf="<<Inf(Sub(xSemantic,z))<<", u4_bootstrap_ms="<<bt<<", x_level="<<xm.level<<", x_towers="<<xm.towers<<'\n';
  }
  std::cout<<"completed_blocks="<<gFailure.completedBlocks<<'\n';return 0;
 }catch(const std::exception&e){if(csv.failure.is_open())Fail(csv,o,e.what());std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}
