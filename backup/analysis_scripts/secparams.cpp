#include "openfhe.h"
using namespace lbcrypto;
int main(int argc,char**argv){
  bool std128 = argc>1;
  const std::vector<uint32_t>budget{1,1},bsgs{0,0};const SecretKeyDist dist=UNIFORM_TERNARY;
  uint32_t bootDepth=FHECKKSRNS::GetBootstrapDepth(budget,dist),depth=bootDepth+3;
  CCParams<CryptoContextCKKSRNS>p;p.SetSecretKeyDist(dist);
  if(std128){p.SetSecurityLevel(HEStd_128_classic);} else {p.SetSecurityLevel(HEStd_NotSet);p.SetRingDim(4096);}
  p.SetBatchSize(4);p.SetNumLargeDigits(3);p.SetKeySwitchTechnique(HYBRID);p.SetScalingTechnique(FLEXIBLEAUTO);p.SetScalingModSize(59);p.SetFirstModSize(60);p.SetMultiplicativeDepth(depth);
  auto cc=GenCryptoContext(p);
  auto cp=std::dynamic_pointer_cast<CryptoParametersRNS>(cc->GetCryptoParameters());
  double logQ=cc->GetElementParams()->GetModulus().GetMSB();
  double logP=cp->GetParamsP()->GetModulus().GetMSB();
  std::cout<<(std128?"HEStd_128_classic":"current(NotSet)")<<": N="<<cc->GetRingDimension()<<" depth="<<depth<<" towers="<<cc->GetElementParams()->GetParams().size()<<" logQ~"<<logQ<<" logP~"<<logP<<" logQP~"<<logQ+logP<<"\n";
}
