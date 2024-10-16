#!/bin/bash

root -l<<EOC
.L code_dedx.C
TChain *chain_dedx = new TChain("stage/ttree");
chain_dedx->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodJune2024_CMSSW_14_0_7/ZeroBias/Run2024C_v1/*.root");
chain_dedx->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodJune2024_CMSSW_14_0_7/ZeroBias/Run2024D_v1/*.root");
run2analysis t2(chain_dedx);
t2.Loop();
EOC

