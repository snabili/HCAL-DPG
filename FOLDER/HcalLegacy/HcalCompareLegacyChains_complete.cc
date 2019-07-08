// -*- C++ -*-
//
// Package:    HcalCompareLegacyChains
// Class:      HcalCompareLegacyChains
// 
/**\class HcalCompareLegacyChains HcalCompareLegacyChains.cc HcalDebug/CompareChans/src/HcalCompareLegacyChains.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  matthias wolf
//         Created:  Fri Aug 26 11:37:21 CDT 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"

#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1Upgrade.h"

#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"

#include "DataFormats/METReco/interface/METFwd.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMetDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoJetDataFormat.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableVertexReconstructor.h"
#include "RecoVertex/TrimmedKalmanVertexFinder/interface/KalmanTrimmedVertexFinder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"

#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"

#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMetFilterDataFormat.h"

#include "FastSimulation/CalorimeterProperties/interface/Calorimeter.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TTree.h"
#include <cmath>
//#include <iostream>
#include <vector>
// class declaration
//

class HcalCompareLegacyChains : public edm::EDAnalyzer {
   public:
      explicit HcalCompareLegacyChains(const edm::ParameterSet&);
      ~HcalCompareLegacyChains();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;

      double get_cosh(const HcalDetId&);
      double get_cosh_ee(const EEDetId&);
      double get_cosh_eb(const EBDetId&);
      double setupGeometry(const CaloGeometry&);
      // ----------member data ---------------------------
      bool first_;

      std::vector<edm::InputTag> frames_;
      std::vector<edm::InputTag> rechits_;
      edm::InputTag digis_;

      edm::ESHandle<CaloGeometry> gen_geo_;
      edm::ESHandle<CaloGeometry> EcalBarrelGeometry_;
      edm::ESHandle<CaloGeometry> EcalEndcapGeometry_;
      edm::EDGetTokenT<l1t::EtSumBxCollection> sumToken_;
      //edm::EDGetTokenT<edm::SortedCollection<EcalRecHit>> barrelEcalHits_token_;
      //edm::EDGetTokenT<edm::SortedCollection<EcalRecHit>> endcapEcalHits_token_;

      edm::EDGetTokenT<reco::PFMETCollection>      metToken_;
      edm::EDGetTokenT<reco::CaloMETCollection>   caloMetToken_;
      edm::EDGetTokenT<reco::CaloMETCollection>   caloMetBEToken_;

      bool swap_iphi_;
      int max_severity_;

      L1Analysis::L1AnalysisRecoMetDataFormat* met_data;
      L1Analysis::L1AnalysisRecoMetFilterDataFormat* metFilter_data;

      TH2D *df_multiplicity_;
      TH2D *tp_multiplicity_;

      TH2D *h_Online_Offline_MET;
      TH2D *h_On_Off_MET;
      TH2D *h_L1MET_CaloMET;
      TH2D *h_Online_Offline_SumEt;
      TH2D *h_OfflineMET_Phi;
      TH2D *h_CaloMET_Phi;
      TH2D *h_OnlineMET_Phi;
      TH2D *h_OnlineMET_OnlineSumEt;
      TH2D *h_OfflineMET_OfflineSumEt;

      TH2D *h_PFMET_RechitM0_Et;
      TH2D *h_PFMET_RechitM0_Et_HF;
      TH2D *h_PFMET_RechitM0_Et_BE;

      TH2D *h_L1MET_RechitM0_Et;

      TH2D *h_CaloMET_RechitM0_Et;
      TH2D *h_CaloMET_RechitM0_Et_HF;
      TH2D *h_CaloMET_RechitM0_Et_BE;
      TH2D *h_CaloMET_BE_RechitM0_Et_BE;

      TH2D *h_OfflineMET_TPEt;
      TH2D *h_OnlineMET_TPEt;

      TH2D *h_OfflineMET_Calo_PF;

      TH2D *h_OfflineCaloMET_RechitM0_Et;

      TH2D *h_OfflineMET_NewRechitM0_Et;
      TH2D *h_OfflineMET_NewTP;
      TH2D *h_OnlineMET_NewTP;

      TH2D * h_Calo_PF_MET;
      TH2D * h_Calo_PF_SumEt;

      TH1D *h_et_kTotalEt;

      TH1D *h_et_offline;
      TH1D *h_PFSumEt;

      TH1D *h_et_kMissingEtHF;
      TH1D *h_et_kMissingEt;

      //TH1D *h_nVtx;

      TH1D *h_RechitM0_HF;
      TH1D *h_RechitM0_BE;
      TH1D *h_RechitM0;

      TH1D *h_CaloMETBE;
      TH1D *h_CaloMET;
      TH1D *h_CaloSumEt;

      TH2D *h_CaloMET_NewRechitM0_BE;

      TH1D *h_NewRechitM0_BE;

      TH2D *h_L1_SumEt_RechitM0_Et;

      TH2D *h_Calo_SumEt_RechitM0_Et_HF;
      TH2D *h_Calo_SumEt_RechitM0_Et;

      TTree *tps_;
      TTree *tpsplit_;
      TTree *events_;
      TTree *matches_;

      TTree * tree_;

      double tp_energy_;
      int tp_ieta_;
      int tp_iphi_;
      int tp_soi_;

      double tpsplit_energy_;
      double tpsplit_oot_;
      int tpsplit_ieta_;
      int tpsplit_iphi_;
      int tpsplit_depth_;
      double tpsplit_ettot_;
      double tpsplit_rise_avg_;
      double tpsplit_rise_rms_;
      double tpsplit_fall_avg_;
      double tpsplit_fall_rms_;

      double ev_rh_energy0_;
      double ev_rh_energy2_;
      double ev_rh_energy3_;
      double ev_tp_energy_;
      int ev_rh_unmatched_;
      int ev_tp_unmatched_;

      double mt_rh_energy0_;
      double mt_rh_energy2_;
      double mt_rh_energy3_;
      double mt_tp_energy_;

      int mt_ieta_;
      int mt_iphi_;
      int mt_version_;
      int mt_tp_soi_;

      //int max_severity_;

      const HcalChannelQuality* status_;
      const HcalSeverityLevelComputer* comp_;
      edm::EDGetTokenT<EcalRecHitCollection>   tok_EB_;
      edm::EDGetTokenT<EcalRecHitCollection>   tok_EE_;
};

HcalCompareLegacyChains::HcalCompareLegacyChains(const edm::ParameterSet& config) :
   edm::EDAnalyzer(),
   first_(true),
   frames_(config.getParameter<std::vector<edm::InputTag>>("dataFrames")),
   rechits_(config.getParameter<std::vector<edm::InputTag>>("recHits")),
   digis_(config.getParameter<edm::InputTag>("triggerPrimitives")),
   swap_iphi_(config.getParameter<bool>("swapIphi")),
   max_severity_(config.getParameter<int>("maxSeverity")) 
{
   consumes<HcalTrigPrimDigiCollection>(digis_);
   consumes<HBHEDigiCollection>(frames_[0]);
   consumes<HFDigiCollection>(frames_[1]);
   consumes<edm::SortedCollection<HBHERecHit>>(rechits_[0]);
   consumes<edm::SortedCollection<HFRecHit>>(rechits_[1]);

   tok_EB_     = consumes<EcalRecHitCollection>(edm::InputTag("ecalRecHit","EcalRecHitsEB"));
   tok_EE_     = consumes<EcalRecHitCollection>(edm::InputTag("ecalRecHit","EcalRecHitsEE"));

   consumes<reco::PFMETCollection> (edm::InputTag("pfMet"));

   sumToken_ = consumes<l1t::EtSumBxCollection>(config.getUntrackedParameter<edm::InputTag>("sumToken"));

   //barrelEcalHits_token_ = consumes<edm::SortedCollection<EcalRecHit>>(config.getParameter<edm::InputTag>("barrelEcalHits"));
   //endcapEcalHits_token_ = consumes<edm::SortedCollection<EcalRecHit>>(config.getParameter<edm::InputTag>("endcapEcalHits"));

   metToken_ = consumes<reco::PFMETCollection>(config.getUntrackedParameter<edm::InputTag>("metToken"));
   caloMetToken_ = consumes<reco::CaloMETCollection>(config.getUntrackedParameter<edm::InputTag>("caloMetToken"));
   caloMetBEToken_ = consumes<reco::CaloMETCollection>(config.getUntrackedParameter("caloMetBEToken",edm::InputTag("caloMetBE")));

   edm::Service<TFileService> fs;
   met_data = new L1Analysis::L1AnalysisRecoMetDataFormat();
   metFilter_data = new L1Analysis::L1AnalysisRecoMetFilterDataFormat();


   df_multiplicity_ = fs->make<TH2D>("df_multiplicity", "DataFrame multiplicity;ieta;iphi", 65, -32.5, 32.5, 72, 0.5, 72.5);
   tp_multiplicity_ = fs->make<TH2D>("tp_multiplicity", "TrigPrim multiplicity;ieta;iphi", 65, -32.5, 32.5, 72, 0.5, 72.5);
   h_et_kMissingEtHF = fs->make<TH1D>("h_et_kMissingEtHF", "BX = 0 and sumType = kMissingEtHF",100,-20,700);   
   h_et_kMissingEt = fs->make<TH1D>("h_et_kMissingEt", "BX = 0 and sumType = kMissingEt",100,-20,3000);
   h_et_kTotalEt = fs->make<TH1D>("h_et_kTotalEt", "sumType = kTotalEt",100,-20,3000);

   h_et_offline = fs->make<TH1D>("h_et_offline", "PF MET",100,-20,2500);
   h_PFSumEt = fs->make<TH1D>("h_PFSumEt","PF SumEt",100,-20,2500);

   h_Online_Offline_MET = fs->make<TH2D>("h_Online_Offline_MET","Scattered Plot for Online Vs Offline MET",100,0,3000,100,0,1800);
   h_On_Off_MET = fs->make<TH2D>("h_On_Off_MET","Scattered Plot for L1MET Vs PF MET",100,0,2200,100,0,1800);
   h_L1MET_CaloMET = fs->make<TH2D>("h_L1MET_CaloMET","Scattered Plot for L1MET Vs Calo MET",100,0,2200,100,0,1800);

   h_Online_Offline_SumEt = fs->make<TH2D>("h_Online_Offline_SumEt","Scattered Plot for L1 Vs PF SumEt",100,0,2500,100,0,2500);

   h_CaloMET_Phi = fs->make<TH2D>("h_CaloMET_Phi","Calo MET Vs Phi",100,-3.5,3.5,100,0,1800);
   h_OnlineMET_Phi = fs->make<TH2D>("h_OnlineMET_Phi","L1MET Vs Phi",100,-3.5,3.5,100,0,2200);
   h_OfflineMET_Phi = fs->make<TH2D>("h_OfflineMET_Phi","PF MET Vs Phi",100,-3.5,3.5,100,0,1800);

   h_OnlineMET_OnlineSumEt = fs->make<TH2D>("h_OnlineMET_OnlineSumEt","L1 MET vs SumEt",100, 0,2500,100,0,2500);
   h_OfflineMET_OfflineSumEt = fs->make<TH2D>("h_OfflineMET_OfflineSumEt", "PF MET vs SumEt",100, 0,1800,100,0,2500);

   h_PFMET_RechitM0_Et = fs->make<TH2D>("h_PFMET_RechitM0_Et","PF MET vs Rechit Method Zero",100,0,1800,100,0,1800);
   h_PFMET_RechitM0_Et_HF = fs->make<TH2D>("h_PFMET_RechitM0_Et_HF","PF MET vs Rechit Method Zero in HF",100,0,1800,100,0,1800);
   h_PFMET_RechitM0_Et_BE = fs->make<TH2D>("h_PFMET_RechitM0_Et_BE","PF MET vs Rechit Method Zero in HBHE",100,0,1800,100,0,1800);

   h_L1MET_RechitM0_Et = fs->make<TH2D>("h_L1MET_RechitM0_Et","L1MET vs Rechit Method Zero",100,0,2500,100,0,1800);

   h_CaloMET_RechitM0_Et = fs->make<TH2D>("h_CaloMET_RechitM0_Et","Calo eMET vs Rechit Method Zero",100,0,600,100,0,600);
   h_CaloMET_RechitM0_Et_HF = fs->make<TH2D>("h_CaloMET_RechitM0_Et_HF","Calo MET vs Rechit Method Zero in HF",100,0,1800,100,0,1800);
   h_CaloMET_RechitM0_Et_BE = fs->make<TH2D>("h_CaloMET_RechitM0_Et_BE","Calo MET vs Rechit Method Zero in HBHE",100,0,600,100,0,600);
   h_CaloMET_BE_RechitM0_Et_BE = fs->make<TH2D>("h_CaloMET_BE_RechitM0_Et_BE","Calo METBE vs Rechit Method Zero in HBHE",100,0,1800,100,0,1800);

   h_OfflineMET_TPEt = fs->make<TH2D>("h_OfflineMET_TPEt","PF MET vs Trigger Primitive Et",100,0,1800,100,0,1800);
   h_OnlineMET_TPEt = fs->make<TH2D>("h_OnlineMET_TPEt","L1MET vs Trigger Primitive Et",100,0,2200,100,0,1800);

   h_OfflineMET_Calo_PF = fs->make<TH2D>("h_OfflineMET_Calo_PF","OfflineMET vs Trigger Primitive Et",100,0,1800,100,0,1800);


   h_OfflineCaloMET_RechitM0_Et = fs->make<TH2D>("h_OfflineCaloMET_RechitM0_Et","CaloMET vs Trigger Primitive Et",100,0,1800,100,0,1800);

   h_OfflineMET_NewRechitM0_Et = fs->make<TH2D>("h_OfflineMET_NewRechitM0_Et","PFMET vs RechitM0 New Geometry Et",100,0,1800,100,0,1800);

   h_OfflineMET_NewTP = fs->make<TH2D>("h_OfflineMET_NewTP","PFMET vs TP New Geometry Et",100,0,1800,100,0,1800);
   h_OnlineMET_NewTP = fs->make<TH2D>("h_OnlineMET_NewTP","L1MET vs TP New Geometry Et",100,0,2200,100,0,1800);


   h_Calo_PF_MET = fs->make<TH2D>("h_Calo_PF_MET","Calo vs PF MET",100,0,1800,100,0,1800);
   h_Calo_PF_SumEt = fs->make<TH2D>("h_Calo_PF_SumEt","Calo vs PF SumEt",100,0,2500,100,0,2500);


   h_RechitM0_HF = fs->make<TH1D>("h_RechitM0_HF","RechitM0 HF",100,-20,2000);
   h_RechitM0_BE = fs->make<TH1D>("h_RechitM0_BE","RechitM0 BE",100,-20,2000);
   h_RechitM0 = fs->make<TH1D>("h_RechitM0","RechitM0",100,-20,700);

   h_CaloMETBE = fs->make<TH1D>("h_CaloMETBE","CaloMET for HB and HE",100,-10,2500);
   h_CaloMET = fs->make<TH1D>("h_CaloMET","CaloMET for HBHE and HF",100,-10,700);
   h_CaloSumEt = fs->make<TH1D>("h_CaloSumEt","Calo SumEt",100,-10,2500);

   h_CaloMET_NewRechitM0_BE = fs->make<TH2D>("h_CaloMET_NewRechitM0_BE","CaloMET vs RechitM0 New Geometry HB+HE",100,0,1800,100,0,1800);

   h_NewRechitM0_BE = fs->make<TH1D>("h_NewRechitM0_BE","New Geometry RechitM0 HB+HE",100,-20,2000);

   h_L1_SumEt_RechitM0_Et = fs->make<TH2D>("h_L1_SumEt_RechitM0_Et","L1SumEt vs RechitM0",100,0,1800,100,0,1800);

   h_Calo_SumEt_RechitM0_Et_HF = fs->make<TH2D>("h_Calo_SumEt_RechitM0_Et_HF","Calo SumEt vs Rechit Method Zero in HF",100,0,3000,100,0,1800);
   h_Calo_SumEt_RechitM0_Et = fs->make<TH2D>("h_Calo_SumEt_RechitM0_Et","Calo SumEt vs Rechit Method Zero",100,0,3000,100,0,1800);

   void SetHistFillStyle(Style_t styl = 3);
   void SetHistFillSize(Size_t size = 1);

   tps_ = fs->make<TTree>("tps", "Trigger primitives");
   tps_->Branch("et", &tp_energy_);
   tps_->Branch("ieta", &tp_ieta_);
   tps_->Branch("iphi", &tp_iphi_);
   tps_->Branch("soi", &tp_soi_);

   tpsplit_ = fs->make<TTree>("tpsplit", "Trigger primitives");
   tpsplit_->Branch("et", &tpsplit_energy_);
   tpsplit_->Branch("oot", &tpsplit_oot_);
   tpsplit_->Branch("ieta", &tpsplit_ieta_);
   tpsplit_->Branch("iphi", &tpsplit_iphi_);
   tpsplit_->Branch("depth", &tpsplit_depth_);
   tpsplit_->Branch("etsum", &tpsplit_ettot_);
   tpsplit_->Branch("rise_avg", &tpsplit_rise_avg_);
   tpsplit_->Branch("rise_rms", &tpsplit_rise_rms_);
   tpsplit_->Branch("fall_avg", &tpsplit_fall_avg_);
   tpsplit_->Branch("fall_rms", &tpsplit_fall_rms_);

   events_ = fs->make<TTree>("events", "Event quantities");
   events_->Branch("RH_energyM0", &ev_rh_energy0_);
   events_->Branch("RH_energyM2", &ev_rh_energy2_);
   events_->Branch("RH_energyM3", &ev_rh_energy3_);
   events_->Branch("TP_energy", &ev_tp_energy_);
   events_->Branch("RH_unmatched", &ev_rh_unmatched_);
   events_->Branch("TP_unmatched", &ev_tp_unmatched_);

   matches_ = fs->make<TTree>("matches", "Matched RH and TP");
   matches_->Branch("RH_energyM0", &mt_rh_energy0_);
   matches_->Branch("RH_energyM2", &mt_rh_energy2_);
   matches_->Branch("RH_energyM3", &mt_rh_energy3_);
   matches_->Branch("TP_energy", &mt_tp_energy_);
   matches_->Branch("ieta", &mt_ieta_);
   matches_->Branch("iphi", &mt_iphi_);
   matches_->Branch("tp_version", &mt_version_);
   matches_->Branch("tp_soi", &mt_tp_soi_);
}

HcalCompareLegacyChains::~HcalCompareLegacyChains() {}

double
HcalCompareLegacyChains::get_cosh(const HcalDetId& id)
{
   const auto *sub_geo = dynamic_cast<const HcalGeometry*>(gen_geo_->getSubdetectorGeometry(id));
   auto eta = sub_geo->getPosition(id).eta();
   return cosh(eta);
}

double
HcalCompareLegacyChains::get_cosh_eb(const EBDetId& id)
{
   const auto *sub_geo_eb = dynamic_cast<const EcalBarrelGeometry*>(EcalBarrelGeometry_->getSubdetectorGeometry(id));
   auto eta_eb = sub_geo_eb->getGeometry(id)->getPosition().eta();
   return cosh(eta_eb);
}

double
HcalCompareLegacyChains::get_cosh_ee(const EEDetId& id)
{
   const auto *sub_geo_ee = dynamic_cast<const EcalEndcapGeometry*>(EcalEndcapGeometry_->getSubdetectorGeometry(id));
   auto eta_ee = sub_geo_ee->getGeometry(id)->getPosition().eta();
   return cosh(eta_ee);
}

void
HcalCompareLegacyChains::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& setup)
{
   edm::ESHandle<HcalChannelQuality> status;
   setup.get<HcalChannelQualityRcd>().get("withTopo", status);
   status_ = status.product();
   edm::ESHandle<HcalSeverityLevelComputer> comp;
   setup.get<HcalSeverityLevelComputerRcd>().get(comp);
   comp_ = comp.product();
}

unsigned int ievt = 0;

void
HcalCompareLegacyChains::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
   using namespace edm;
   using namespace reco;

   ++ievt;

   double OnlineMET = 0;
   double OnlineSumEt = 0;
   double OnlinePhi = 0;
   double OfflineMET = 0;
   double OnMET = 0;
   //double OfflineSumEt = 0;
   //double METPhi = 0;

   //============Offline PFMET ===================================================================== 

   edm::Handle<reco::PFMETCollection> metLabel_;
   event.getByToken(metToken_, metLabel_);
   const reco::PFMETCollection *metCol = metLabel_.product();
   const reco::PFMET theMet = metCol->front();

   met_data->met     = theMet.et();
   met_data->metPhi  = theMet.phi();
   met_data->sumEt   = theMet.sumEt();


   //============Online MET =====================================================================  
  
   edm::Handle<l1t::EtSumBxCollection> sums;
   event.getByToken(sumToken_, sums);
   if(sums.isValid()){// && MET==true){
     for (int ibx = sums->getFirstBX(); ibx <= sums->getLastBX(); ++ibx) {
       for (l1t::EtSumBxCollection::const_iterator it=sums->begin(ibx); it!=sums->end(ibx); it++) {
         int type = static_cast<int>( it->getType() );
         
	 
         if (type == 8 && ibx ==0){
	   if(it->et()<1000){
             h_et_kMissingEtHF->Fill(it->et());
             OnlineMET = it->et();
             OnlinePhi = it->phi();
           }
           OnMET = it->et();
         }
 
         if (type == 2 && ibx ==0){
	   if(it->et()<1000){
             h_et_kMissingEt->Fill(it->et());
           }
         }
         if (type == 0 && ibx ==0){
	   if(it->et()<1000){
             h_et_kTotalEt->Fill(it->et());
	     OnlineSumEt = it->et();
	   }

         }
 
       }
     }
   }
   else{
      std::cout<<" invalid sums"<<std::endl;
   }

   //unsigned int y = 0;

   //============Offline CaloMET =====================================================================  
   edm::Handle<reco::CaloMETCollection> caloMet;
   event.getByToken(caloMetToken_, caloMet);
   const reco::CaloMETCollection *CalometCol = caloMet.product();
   const reco::CaloMET CaloMet = CalometCol->front();
   met_data->caloMet     = CaloMet.et();
   met_data->caloMetPhi  = CaloMet.phi();
   met_data->caloSumEt   = CaloMet.sumEt();

   edm::Handle<reco::CaloMETCollection> caloMetBE;
   event.getByToken(caloMetBEToken_, caloMetBE);
   const reco::CaloMETCollection *BEmetCol = caloMetBE.product();
   const reco::CaloMET CaloMetBE = BEmetCol->front();
   met_data->caloMetBE   = CaloMetBE.et(); 
  
  
   //==========Offline and Online Histograms=======================================================
   if(OnMET<1000){
     h_Online_Offline_SumEt->Fill(OnlineSumEt,theMet.sumEt());
     h_Online_Offline_SumEt->SetMarkerStyle(3);
     h_Online_Offline_SumEt->SetMarkerSize(1);
     h_Online_Offline_SumEt->GetXaxis()->SetTitle("L1 Sum Et (GeV)");
     h_Online_Offline_SumEt->GetYaxis()->SetTitle("PF Sum Et (GeV)");
   
     h_OfflineMET_Phi->Fill(theMet.phi(),theMet.et());
     h_OfflineMET_Phi->SetMarkerStyle(3);
     h_OfflineMET_Phi->SetMarkerSize(1);
     h_OfflineMET_Phi->GetYaxis()->SetTitle("Phi ");
     h_OfflineMET_Phi->GetXaxis()->SetTitle("PF MET (GeV)");


     h_CaloMET_Phi->Fill(CaloMet.phi(),CaloMet.et());
     h_CaloMET_Phi->SetMarkerStyle(3);
     h_CaloMET_Phi->SetMarkerSize(1);
     h_CaloMET_Phi->GetYaxis()->SetTitle("Phi ");
     h_CaloMET_Phi->GetXaxis()->SetTitle("Calo MET (GeV)");


     h_OnlineMET_Phi->Fill(OnlinePhi,OnlineMET);
     h_OnlineMET_Phi->SetMarkerStyle(3);
     h_OnlineMET_Phi->SetMarkerSize(1);
     h_OnlineMET_Phi->GetYaxis()->SetTitle("Phi ");
     h_OnlineMET_Phi->GetXaxis()->SetTitle("L1MET (GeV)");
   
     h_et_offline->Fill(theMet.et());
     h_et_offline->GetYaxis()->SetTitle("PF MET (GeV)");

     h_On_Off_MET->Fill(OnlineMET,theMet.et());
     h_On_Off_MET->SetMarkerStyle(3);
     h_On_Off_MET->SetMarkerSize(1);
     h_On_Off_MET->GetXaxis()->SetTitle("L1MET (GeV)");
     h_On_Off_MET->GetYaxis()->SetTitle("PF MET (GeV)");

     h_L1MET_CaloMET->Fill(OnlineMET,CaloMet.et());
     h_L1MET_CaloMET->SetMarkerStyle(3);
     h_L1MET_CaloMET->SetMarkerSize(1);
     h_L1MET_CaloMET->GetXaxis()->SetTitle("L1MET (GeV)");
     h_L1MET_CaloMET->GetYaxis()->SetTitle("Calo MET (GeV)");
   
     h_OfflineMET_OfflineSumEt->Fill(theMet.et(),theMet.sumEt());
     h_OfflineMET_OfflineSumEt->SetMarkerStyle(3);
     h_OfflineMET_OfflineSumEt->SetMarkerSize(1);
     h_OfflineMET_OfflineSumEt->GetXaxis()->SetTitle("PF MET (GeV)");
     h_OfflineMET_OfflineSumEt->GetYaxis()->SetTitle("PF Sum Et (GeV)");
   
     h_OnlineMET_OnlineSumEt->Fill(OnlineMET,OnlineSumEt);
     h_OnlineMET_OnlineSumEt->SetMarkerStyle(3);
     h_OnlineMET_OnlineSumEt->SetMarkerSize(1);
     h_OnlineMET_OnlineSumEt->GetXaxis()->SetTitle("L1MET (GeV)");
     h_OnlineMET_OnlineSumEt->GetYaxis()->SetTitle("L1 Sum Et (GeV)");   

     h_PFSumEt->Fill(theMet.sumEt());
     h_PFSumEt->GetYaxis()->SetTitle("PF SumEt (GeV)");

     h_CaloMETBE->Fill(CaloMetBE.et());
     h_CaloMET->Fill(CaloMet.et());
     h_CaloSumEt->Fill(CaloMet.sumEt());


     h_Calo_PF_MET->Fill(CaloMet.et(),theMet.et());
     h_Calo_PF_MET->SetMarkerStyle(3);
     h_Calo_PF_MET->SetMarkerSize(1);
     h_Calo_PF_MET->GetXaxis()->SetTitle("CaloMET (GeV)");
     h_Calo_PF_MET->GetYaxis()->SetTitle("PFMET (GeV)");

     h_Calo_PF_SumEt->Fill(CaloMet.sumEt(),theMet.sumEt());
     h_Calo_PF_SumEt->SetMarkerStyle(3);
     h_Calo_PF_SumEt->SetMarkerSize(1);
     h_Calo_PF_SumEt->GetXaxis()->SetTitle("Calo SumEt (GeV)");
     h_Calo_PF_SumEt->GetYaxis()->SetTitle("PF SumEt (GeV)");
   }
  //-------------------------------------------------------------------------------------------------------------

   edm::ESHandle<HcalTrigTowerGeometry> tpd_geo_h;
   setup.get<CaloGeometryRecord>().get(tpd_geo_h);
   edm::ESHandle<HcalDbService> conditions;
   setup.get<HcalDbRecord>().get(conditions);
   const HcalTrigTowerGeometry& tpd_geo = *tpd_geo_h;

   // ==========
   // Dataframes
   // ==========

   Handle<HBHEDigiCollection> frames;
   Handle<HFDigiCollection> hfframes;
   if (frames_.size() == 2) {
      if (first_ && event.getByLabel(frames_[0], frames) && event.getByLabel(frames_[1], hfframes)) {
         std::set<HcalTrigTowerDetId> ids;

         for (const auto& frame: *(frames.product())) {
            auto mapped = tpd_geo_h->towerIds(frame.id());

            for (const auto& id: mapped) {
               df_multiplicity_->Fill(id.ieta(), id.iphi());
               ids.insert(id);
            }
         }

         for (const auto& frame: *(hfframes.product())) {
            auto mapped = tpd_geo_h->towerIds(frame.id());

            for (const auto& id: mapped) {
               df_multiplicity_->Fill(id.ieta(), id.iphi());
               ids.insert(id);
            }
         }

         for (const auto& id: ids) {
            tp_multiplicity_->Fill(id.ieta(), id.iphi());
         }

         first_ = false;
      }
   }

   // ==============
   // Matching stuff
   // ==============

   ev_rh_energy0_ = 0.;
   ev_rh_energy2_ = 0.;
   ev_rh_energy3_ = 0.;
   ev_rh_unmatched_ = 0.;
   ev_tp_energy_ = 0.;
   ev_tp_unmatched_ = 0.;

   std::map<HcalTrigTowerDetId, std::vector<HBHERecHit>> rhits;
   std::map<HcalTrigTowerDetId, std::vector<HFRecHit>> fhits;
   std::map<EcalTrigTowerDetId, std::vector<EcalRecHit>> ehits;
   std::map<HcalTrigTowerDetId, std::vector<HcalTriggerPrimitiveDigi>> tpdigis;



   Handle<HcalTrigPrimDigiCollection> digis;
   if (!event.getByLabel(digis_, digis)) {
      LogError("HcalTrigPrimDigiCleaner") <<
         "Can't find hcal trigger primitive digi collection with tag '" <<
         digis_ << "'" << std::endl;
      return;
   }

   edm::Handle< edm::SortedCollection<HBHERecHit> > hits;
   if (!event.getByLabel(rechits_[0], hits)) {
      edm::LogError("HcalCompareLegacyChains") <<
         "Can't find rec hit collection with tag '" << rechits_[0] << "'" << std::endl;
      /* return; */
   }

   edm::Handle< edm::SortedCollection<HFRecHit> > hfhits;
   if (!event.getByLabel(rechits_[1], hfhits)) {
      edm::LogError("HcalCompareLegacyChains") <<
         "Can't find rec hit collection with tag '" << rechits_[1] << "'" << std::endl;
      /* return; */
   }


   //std::cout<<"SCE getting handle"<<std::endl;
   edm::Handle<EcalRecHitCollection> ecalhits;

   if (!event.getByToken(tok_EB_, ecalhits)) {
     edm::LogError("HcalCompareLegacyChains") <<
       "Can't find rec hit collection" << std::endl;
      /* return; */
   }



   if (!event.getByToken(tok_EE_, ecalhits)) {
     edm::LogError("HcalCompareLegacyChains") <<
       "Can't find rec hit collection" << std::endl;
      /* return; */
   }


   /*edm::Handle< edm::SortedCollection<EcalRecHit> > eecalhits;
   //edm::Handle<EcalRecHit> eecalhits;
   if (!event.getByToken(endcapEcalHits_token_, eecalhits)) {
      edm::LogError("HcalCompareLegacyChains") <<
         "Can't find rec hit collection with tag '" << endcapEcalHits_token_ << "'" << std::endl;
       return; 
   }

   edm::Handle< edm::SortedCollection<EcalRecHit> > ebcalhits;
   //edm::Handle<EcalRecHit> ebcalhits;
   if (!event.getByToken(barrelEcalHits_token_, ebcalhits)) {
      edm::LogError("HcalCompareLegacyChains") <<
         "Can't find rec hit collection with tag '" << barrelEcalHits_token_ << "'" << std::endl;
       return; 
   }*/

   setup.get<CaloGeometryRecord>().get(gen_geo_);
   //setup.get<CaloGeometryRecord>().get(EcalBarrelGeometry_);
   //setup.get<CaloGeometryRecord>().get(EcalEndcapGeometry_);

   
   //std::cout<<" pain number 1"<<std::endl;
   auto isValid = [&](const auto& hit) {
      HcalDetId id(hit.id());
      auto s = status_->getValues(id);
      int level = comp_->getSeverityLevel(id, 0, s->getValue());
      return level <= max_severity_;
   };

   /*std::cout<<" pain number 2"<<std::endl;
   auto isValideb = [&](const auto& hiteb) {
      EBDetId ideb(hiteb.id());
      auto s = status_->getValues(ideb);
      int level = comp_->getSeverityLevel(ideb, 0, s->getValue());
      return level <= max_severity_;
   };
   */

   /*std::cout<<" pain number 3"<<std::endl;
   auto isValidee = [&](const auto& hitee) {
      EEDetId idee(hitee.id());
      auto s = status_->getValues(idee);
      int level = comp_->getSeverityLevel(idee, 0, s->getValue());
      return level <= max_severity_;
   };*/

   // Get EcalRecHits  taken from RecoEgamma/PhotonIdentification/src/PhotonMIPHaloTagger.cc

   //edm::Handle<EcalRecHitCollection> barrelHitHandle;
   //EcalRecHitCollection barrelRecHits;
   //event.getByToken(barrelEcalHits_token_, barrelHitHandle);
   //edm::SortedCollection<EcalRecHit> endcapHitHandle;
   //edm::Handle<EcalRecHitCollection> endcapHitHandle;
   //event.getByToken(endcapEcalHits_token_, endcapHitHandle);
   //EcalRecHitCollection endcapRecHits;

   double RechitM0_Et = 0;
   double RechitM0_Et_HF = 0;
   double RechitM0_Et_BE = 0;
   double Rechit0Ethbhe = 0;
   double Rechit0Ethf = 0;
   //double Rechit0Eteb = 0;
   //double Rechit0Etee = 0;

   double ev_phi_hbhe = 0;
   double ev_phi_hf = 0;
   //double ev_phi_ee = 0;
   //double ev_phi_eb = 0;

   double Ex = 0;
   double Ey = 0;

   double hbhe_Ex = 0;
   double hbhe_Ey = 0;
   double hf_Ex = 0;
   double hf_Ey = 0;
   double eb_Ex = 0;
   //double eb_Ey = 0;
   double ee_Ex = 0;
   double ee_Ey = 0;
  
   if (hits.isValid()){
     OfflineMET =theMet.et();
     for (auto& hit: *(hits.product())) {
        HcalDetId id(hit.id());
        if (not isValid(hit))
        continue;
        ev_rh_energy0_ += hit.eraw() / get_cosh(id);
        ev_rh_energy2_ += hit.energy() / get_cosh(id);
        ev_rh_energy3_ += hit.eaux() / get_cosh(id);

        ev_phi_hbhe = id.iphi()*5*M_PI/180;
        Rechit0Ethbhe = hit.energy() / get_cosh(id);
        hbhe_Ex += Rechit0Ethbhe*cos(ev_phi_hbhe);
  	hbhe_Ey += Rechit0Ethbhe*sin(ev_phi_hbhe);
        auto tower_ids = tpd_geo.towerIds(id);
        for (auto& tower_id: tower_ids) {
          tower_id = HcalTrigTowerDetId(tower_id.ieta(), tower_id.iphi(), 1);
          rhits[tower_id].push_back(hit);
        }  
     }
   }


   if (hfhits.isValid()) {
      for (auto& hit: *(hfhits.product())) {
         HcalDetId id(hit.id());
         if (not isValid(hit))
            continue;
         ev_phi_hf = id.iphi()*5*M_PI/180;
         Rechit0Ethf = hit.energy() / get_cosh(id);
         hf_Ex += Rechit0Ethf*cos(ev_phi_hf);
         hf_Ey += Rechit0Ethf*sin(ev_phi_hf);
         ev_rh_energy0_ += hit.energy() / get_cosh(id);
         ev_rh_energy2_ += hit.energy() / get_cosh(id);
         ev_rh_energy3_ += hit.energy() / get_cosh(id);

         auto tower_ids = tpd_geo.towerIds(id);
         for (auto& tower_id: tower_ids) {
            tower_id = HcalTrigTowerDetId(tower_id.ieta(), tower_id.iphi(), 1, tower_id.version());
            fhits[tower_id].push_back(hit);
         }
      }
   }
    double ev_phi_eb = 0;
    if (ecalhits.isValid()) {
      for (auto& ebhit: *(ecalhits.product())) {
         EBDetId id(ebhit.id());
         //const auto *sub_geo_eb = dynamic_cast<const EcalBarrelGeometry*>(EcalBarrelGeometry_->getSubdetectorGeometry(id));
         //auto eta_eb = sub_geo_eb->getGeometry(id)->getPosition().eta();
	 ev_phi_eb = id.iphi()*5*M_PI/180;
         //Rechit0Eteb = ebhit.energy() / cosh(eta_eb);
         //eb_Ex += Rechit0Eteb*cos(ev_phi_eb);
         //eb_Ey += Rechit0Eteb*sin(ev_phi_eb);
       }
    }


    /*if (ecalhits.isValid()) {
      for (auto& eehit: *(ecalhits.product())) {
         EEDetId id(eehit.id());
         const auto *sub_geo_ee = dynamic_cast<const EcalEndcapGeometry*>(EcalEndcapGeometry_->getSubdetectorGeometry(id));
         auto eta_ee = sub_geo_ee->getGeometry(id)->getPosition().eta();
         ev_phi_ee = id.iphi()*5*M_PI/180;
         Rechit0Etee = eehit.energy() / cosh(eta_ee);
         ee_Ex += Rechit0Etee*cos(ev_phi_ee);
         ee_Ey += Rechit0Etee*sin(ev_phi_ee);
       }
    }*/


  /* if(barrelHitHandle.isValid()){
     for (auto& hit: *(barrelHitHandle.product())) {
        EBDetId ideb(hit.id());
        if (not isValid(hit))
        continue;
        ev_phi_ee = ideb.iphi()*5*M_PI/180;
        Rechit0Etee = hit.energy() / get_cosh_ee(ideb);
        ee_Ex += Rechit0Etee*cos(ev_phi_ee);
        ee_Ey += Rechit0Etee*sin(ev_phi_ee);
        auto tower_ids = tpd_geo.towerIds(ideb);
        for (auto& tower_ideb: tower_ids) {
          tower_ideb = EcalTrigTowerDetId( 1, EcalEndcap, tower_ideb.ieta(), tower_ideb.iphi(), 1);
          ehits[tower_ideb].push_back(hit);
        }
     }
   }*/
   //std::cout<<" pain number 4"<<std::endl;
   /*if(endcapHitHandle.isValid()){
   std::cout<<" pain number 4a"<<std::endl;
     for (auto& hitee: *(endcapHitHandle.product())) {
   std::cout<<" pain number 4b"<<std::endl;
        EBDetId idee(hitee.id());
   std::cout<<" pain number 4c"<<std::endl;
        if (not isValidee(hitee))
        continue;
   std::cout<<" pain number 4d"<<std::endl;
        ev_phi_eb = idee.iphi()*5*M_PI/180;
   std::cout<<" pain number 4e"<<std::endl;
        Rechit0Eteb = hitee.energy() / get_cosh_eb(idee);
   std::cout<<" pain number 4f"<<std::endl;
        eb_Ex += Rechit0Eteb*cos(ev_phi_eb);
   std::cout<<" pain number 4g"<<std::endl;
        eb_Ey += Rechit0Eteb*sin(ev_phi_eb);
   std::cout<<" pain number 4h"<<std::endl;
        auto tower_ids = tpd_geo.towerIds(idee);
   std::cout<<" pain number 4i"<<std::endl;
        for (auto& tower_idee: tower_ids) {
   std::cout<<" pain number 4j"<<std::endl;
          tower_idee = EcalTrigTowerDetId( 1, EcalBarrel, tower_idee.ieta(), tower_idee.iphi(), 1);
   std::cout<<" pain number 4k"<<std::endl;
          ehits[tower_idee].push_back(hitee);
   std::cout<<" pain number 4l"<<std::endl;
        }
     }
   }*/
   //std::cout<<" pain number 5"<<std::endl;

   if(OnMET<1000){
     Ex = hf_Ex + hbhe_Ex + ee_Ex + eb_Ex;
     Ey = hf_Ey + hbhe_Ey + ee_Ey + eb_Ex;

     RechitM0_Et = sqrt(Ex*Ex+Ey*Ey);
     RechitM0_Et_HF = sqrt(hf_Ex*hf_Ex+hf_Ey*hf_Ey);
     RechitM0_Et_BE = sqrt(hbhe_Ex*hbhe_Ex+hbhe_Ey*hbhe_Ey);

     h_PFMET_RechitM0_Et->Fill(OfflineMET,RechitM0_Et);
     h_PFMET_RechitM0_Et->SetMarkerStyle(3);
     h_PFMET_RechitM0_Et->SetMarkerSize(1);
     h_PFMET_RechitM0_Et->SetMarkerColor(8);
     h_PFMET_RechitM0_Et->GetXaxis()->SetTitle("PF MET (GeV)");
     h_PFMET_RechitM0_Et->GetYaxis()->SetTitle("RechitM0 (GeV)");
     h_PFMET_RechitM0_Et->GetYaxis()->SetTitleOffset(1.5);

     h_RechitM0->Fill(RechitM0_Et);

     h_CaloMET_RechitM0_Et->Fill(CaloMet.et(),RechitM0_Et);
     h_CaloMET_RechitM0_Et->SetMarkerStyle(3);
     h_CaloMET_RechitM0_Et->SetMarkerSize(1);
     h_CaloMET_RechitM0_Et->SetMarkerColor(4);
     h_CaloMET_RechitM0_Et->GetXaxis()->SetTitle("Calo MET (GeV)");
     h_CaloMET_RechitM0_Et->GetYaxis()->SetTitle("RechitM0 (GeV)");
     h_CaloMET_RechitM0_Et->GetYaxis()->SetTitleOffset(1.5);

     h_L1MET_RechitM0_Et->Fill(OnlineMET,RechitM0_Et);
     h_L1MET_RechitM0_Et->SetMarkerStyle(3);
     h_L1MET_RechitM0_Et->SetMarkerSize(1);
     h_L1MET_RechitM0_Et->SetMarkerColor(2);
     h_L1MET_RechitM0_Et->GetXaxis()->SetTitle("L1MET (GeV)");
     h_L1MET_RechitM0_Et->GetYaxis()->SetTitle("RechitM0 (GeV)");
     h_L1MET_RechitM0_Et->GetYaxis()->SetTitleOffset(1.5);

     h_L1_SumEt_RechitM0_Et->Fill(OnlineSumEt,RechitM0_Et);
     h_L1_SumEt_RechitM0_Et->SetMarkerStyle(3);
     h_L1_SumEt_RechitM0_Et->SetMarkerSize(1);
     h_L1_SumEt_RechitM0_Et->SetMarkerColor(2);
     h_L1_SumEt_RechitM0_Et->GetXaxis()->SetTitle("L1SumEt (GeV)");
     h_L1_SumEt_RechitM0_Et->GetYaxis()->SetTitle("RechitM0 (GeV)");
     h_L1_SumEt_RechitM0_Et->GetYaxis()->SetTitleOffset(1.5);

     h_CaloMET_RechitM0_Et_HF->Fill(CaloMet.et(),RechitM0_Et_HF);
     h_CaloMET_RechitM0_Et_HF->SetMarkerStyle(3);
     h_CaloMET_RechitM0_Et_HF->SetMarkerSize(1);
     h_CaloMET_RechitM0_Et_HF->SetMarkerColor(9);
     h_CaloMET_RechitM0_Et_HF->GetXaxis()->SetTitle("Calo MET (GeV)");
     h_CaloMET_RechitM0_Et_HF->GetYaxis()->SetTitle("RechitM0 (GeV)");
     h_CaloMET_RechitM0_Et_HF->GetYaxis()->SetTitleOffset(1.5);

     h_Calo_SumEt_RechitM0_Et_HF->Fill(CaloMet.sumEt(),RechitM0_Et_HF);
     h_Calo_SumEt_RechitM0_Et_HF->SetMarkerStyle(3);
     h_Calo_SumEt_RechitM0_Et_HF->SetMarkerSize(1);
     h_Calo_SumEt_RechitM0_Et_HF->SetMarkerColor(9);
     h_Calo_SumEt_RechitM0_Et_HF->GetXaxis()->SetTitle("Calo SumEt (GeV)");
     h_Calo_SumEt_RechitM0_Et_HF->GetYaxis()->SetTitle("RechitM0 (GeV)");
     h_Calo_SumEt_RechitM0_Et_HF->GetYaxis()->SetTitleOffset(1.5);

     h_Calo_SumEt_RechitM0_Et->Fill(CaloMet.sumEt(),RechitM0_Et);
     h_Calo_SumEt_RechitM0_Et->SetMarkerStyle(3);
     h_Calo_SumEt_RechitM0_Et->SetMarkerSize(1);
     h_Calo_SumEt_RechitM0_Et->SetMarkerColor(9);
     h_Calo_SumEt_RechitM0_Et->GetXaxis()->SetTitle("Calo SumEt (GeV)");
     h_Calo_SumEt_RechitM0_Et->GetYaxis()->SetTitle("RechitM0 (GeV)");
     h_Calo_SumEt_RechitM0_Et->GetYaxis()->SetTitleOffset(1.5);

     h_CaloMET_BE_RechitM0_Et_BE->Fill(CaloMetBE.et(),RechitM0_Et_BE);
     h_CaloMET_BE_RechitM0_Et_BE->SetMarkerStyle(3);
     h_CaloMET_BE_RechitM0_Et_BE->SetMarkerSize(1);
     h_CaloMET_BE_RechitM0_Et_BE->SetMarkerColor(6);
     h_CaloMET_BE_RechitM0_Et_BE->GetXaxis()->SetTitle("Calo MET BE (GeV)");
     h_CaloMET_BE_RechitM0_Et_BE->GetYaxis()->SetTitle("RechitM0 (GeV)");
     h_CaloMET_BE_RechitM0_Et_BE->GetYaxis()->SetTitleOffset(1.5);

     h_CaloMET_RechitM0_Et_BE->Fill(CaloMet.et(),RechitM0_Et_BE);
     h_CaloMET_RechitM0_Et_BE->SetMarkerStyle(3);
     h_CaloMET_RechitM0_Et_BE->SetMarkerSize(1);
     h_CaloMET_RechitM0_Et_BE->SetMarkerColor(6);
     h_CaloMET_RechitM0_Et_BE->GetXaxis()->SetTitle("Calo MET (GeV)");
     h_CaloMET_RechitM0_Et_BE->GetYaxis()->SetTitle("RechitM0 (GeV)");
     h_CaloMET_RechitM0_Et_BE->GetYaxis()->SetTitleOffset(1.5);

     h_PFMET_RechitM0_Et_HF->Fill(theMet.et(),RechitM0_Et_HF);
     h_PFMET_RechitM0_Et_HF->SetMarkerStyle(3);
     h_PFMET_RechitM0_Et_HF->SetMarkerSize(1);
     h_PFMET_RechitM0_Et_HF->SetMarkerColor(9);
     h_PFMET_RechitM0_Et_HF->GetXaxis()->SetTitle("Calo MET (GeV)");
     h_PFMET_RechitM0_Et_HF->GetYaxis()->SetTitle("RechitM0 (GeV)");
     h_PFMET_RechitM0_Et_HF->GetYaxis()->SetTitleOffset(1.5);


     h_PFMET_RechitM0_Et_BE->Fill(theMet.et(),RechitM0_Et_BE);
     h_PFMET_RechitM0_Et_BE->SetMarkerStyle(3);
     h_PFMET_RechitM0_Et_BE->SetMarkerSize(1);
     h_PFMET_RechitM0_Et_BE->SetMarkerColor(6);
     h_PFMET_RechitM0_Et_BE->GetXaxis()->SetTitle("Calo MET (GeV)");
     h_PFMET_RechitM0_Et_BE->GetYaxis()->SetTitle("RechitM0 (GeV)");
     h_PFMET_RechitM0_Et_BE->GetYaxis()->SetTitleOffset(1.5);

     h_RechitM0_HF->Fill(RechitM0_Et_HF);
     h_RechitM0_BE->Fill(RechitM0_Et_BE);
   }

   /*std::cout<<"SCE setup geo"<<std::endl;

   setup.get<CaloGeometryRecord>().get(EcalBarrelGeometry_);

   std::cout<<"SCE loop"<<std::endl;
   if (ecalhits.isValid()) {
     std::cout<<"1"<<std::endl;
      for (auto& ebhit: *(ecalhits.product())) {
     std::cout<<"2"<<std::endl;
         EBDetId id(ebhit.id());
     std::cout<<"3"<<std::endl;
         std::cout<<"ecal energy is "<<ebhit.energy()<<std::endl;
         const auto *sub_geo_eb = dynamic_cast<const EcalBarrelGeometry*>(EcalBarrelGeometry_->getSubdetectorGeometry(id));
         auto eta_eb = sub_geo_eb->getGeometry(id)->getPosition().eta();
         std::cout<<"eta is"<<eta_eb<<std::endl;
	 
     
       }
   }*/


   ESHandle<CaloTPGTranscoder> decoder;
   setup.get<CaloTPGRecord>().get(decoder);

   double tp_phi = 0;
   double tp_Ex = 0;
   double tp_Ey = 0;
   double tp_Et = 0;

   for (const auto& digi: *digis){
     HcalTrigTowerDetId id = digi.id();
     id = HcalTrigTowerDetId(id.ieta(), id.iphi(), 1, id.version());
     ev_tp_energy_ += decoder->hcaletValue(id, digi.t0());
     tpdigis[id].push_back(digi);
     tp_energy_ = decoder->hcaletValue(id, digi.t0());
     tp_ieta_ = id.ieta();
     tp_iphi_ = id.iphi();
     tp_phi = tp_iphi_*5*M_PI/180;
     tp_Ex += tp_energy_*cos(tp_phi);
     tp_Ey += tp_energy_*sin(tp_phi);
     tp_soi_ = digi.SOI_compressedEt();
 
     tps_->Fill();    
    
   }
   tp_Et = sqrt(tp_Ex*tp_Ex+tp_Ey*tp_Ey);
   h_OfflineMET_TPEt->Fill(theMet.et(),tp_Et);
   h_OfflineMET_TPEt->SetMarkerStyle(3);
   h_OfflineMET_TPEt->SetMarkerSize(1);
   h_OfflineMET_TPEt->GetXaxis()->SetTitle("PF MET (GeV)");
   h_OfflineMET_TPEt->GetYaxis()->SetTitle("Trigger Primitive (GeV)");

   h_OnlineMET_TPEt->Fill(OnlineMET,tp_Et);
   h_OnlineMET_TPEt->SetMarkerStyle(3);
   h_OnlineMET_TPEt->SetMarkerSize(1);
   h_OnlineMET_TPEt->GetXaxis()->SetTitle("L1MET (GeV)");
   h_OnlineMET_TPEt->GetYaxis()->SetTitle("Trigger Primitive (GeV)");

   double Extp = 0;
   double Eytp = 0;
   double NewRechitM0_Et = 0;
   double NewRechitM0_Et_BE = 0;
   double new_Et = 0;
   double new_tp_Ex = 0;
   double new_tp_Ey = 0;

   double rh_Ex = 0;
   double rh_Ey = 0;

   for (const auto& pair: tpdigis) {
      auto id = pair.first;

      auto new_id(id);
      if (swap_iphi_ and id.version() == 1 and id.ieta() > 28 and id.ieta() < 40) {
         if (id.iphi() % 4 == 1)
            new_id = HcalTrigTowerDetId(id.ieta(), (id.iphi() + 70) % 72, id.depth(), id.version());
         else
            new_id = HcalTrigTowerDetId(id.ieta(), (id.iphi() + 2) % 72 , id.depth(), id.version());
      }

      auto rh = rhits.find(new_id);
      auto fh = fhits.find(new_id);

      if (rh != rhits.end() and fh != fhits.end()) {
         assert(0);
      }

      double newgeo_phi = 0;
      mt_ieta_ = new_id.ieta();
      mt_iphi_ = new_id.iphi();
      newgeo_phi = new_id.iphi()*5*M_PI/180;
      mt_version_ = new_id.version();
      mt_tp_energy_ = 0;
      mt_tp_soi_ = 0;

      double newgeotpet = 0;

      for (const auto& tp: pair.second) {
         mt_tp_energy_ += decoder->hcaletValue(new_id, tp.t0());
         newgeotpet = decoder->hcaletValue(new_id, tp.t0());
	 new_tp_Ex += newgeotpet*cos(newgeo_phi);
         new_tp_Ey += newgeotpet*sin(newgeo_phi);
         
         mt_tp_soi_ = tp.SOI_compressedEt();
      }
           
      mt_rh_energy0_ = 0.;
      mt_rh_energy2_ = 0.;
      mt_rh_energy3_ = 0.;

      double fh_Ex = 0;
      double fh_Ey = 0;

      //double NewRechitM0_rh = 0;
      //double NewRechitM0_fh = 0;

      //double newIDphi = mt_iphi_*5*M_PI/180;

      /*if (rh != rhits.end()) {
         for (const auto& hit: rh->second) {
            HcalDetId id(hit.id());
            auto tower_ids = tpd_geo.towerIds(id);
            auto count = std::count_if(std::begin(tower_ids), std::end(tower_ids),
                  [&](const auto& t) { return t.version() == new_id.version(); });
            mt_rh_energy0_ += hit.eraw() / get_cosh(id) / count;
            mt_rh_energy2_ += hit.energy() / get_cosh(id) / count;
            mt_rh_energy3_ += hit.eaux() / get_cosh(id) / count;

            NewRechitM0_rh = hit.eraw() / get_cosh(id) / count;
            rh_Ex += NewRechitM0_rh*cos(newIDphi);
	    rh_Ey += NewRechitM0_rh*sin(newIDphi);
            
         }
         matches_->Fill();
         rhits.erase(rh);
      } else if (fh != fhits.end()) {
         for (const auto& hit: fh->second) {
            HcalDetId id(hit.id());
            auto tower_ids = tpd_geo.towerIds(id);
            auto count = std::count_if(std::begin(tower_ids), std::end(tower_ids),
                  [&](const auto& t) { return t.version() == new_id.version(); });
            mt_rh_energy0_ += hit.energy() / get_cosh(id) / count;
            mt_rh_energy2_ += hit.energy() / get_cosh(id) / count;
            mt_rh_energy3_ += hit.energy() / get_cosh(id) / count;

	    NewRechitM0_rh = hit.energy() / get_cosh(id) / count;
            fh_Ex += NewRechitM0_fh*cos(newIDphi);
            fh_Ey += NewRechitM0_fh*sin(newIDphi);
            //std::cout<<"Ex fh: "<<fh_Ex<<" Ey fh: "<<fh_Ey<<std::endl;
         }
         
         matches_->Fill();
         fhits.erase(fh);
      } else {
         ++ev_tp_unmatched_;
      }*/

      Extp += rh_Ex+fh_Ex;
      Eytp += rh_Ey+fh_Ey;

   }
   if(OnlineMET<1000){
     new_Et = sqrt(new_tp_Ex*new_tp_Ex+new_tp_Ey*new_tp_Ey);
     h_OnlineMET_NewTP->Fill(theMet.et(),new_Et);
     h_OnlineMET_NewTP->SetMarkerStyle(3);
     h_OnlineMET_NewTP->Fill(OnlineMET,new_Et);
     h_OnlineMET_NewTP->GetXaxis()->SetTitle("L1MET (GeV)");
     h_OnlineMET_NewTP->GetYaxis()->SetTitle("TP (GeV)");

     NewRechitM0_Et = sqrt(Extp*Extp+Eytp*Eytp);
     h_OfflineMET_NewRechitM0_Et->Fill(theMet.et(),NewRechitM0_Et);
     h_OfflineMET_NewRechitM0_Et->SetMarkerStyle(3);
     h_OfflineMET_NewRechitM0_Et->SetMarkerSize(1);
     h_OfflineMET_NewRechitM0_Et->GetXaxis()->SetTitle("OfflineMET (GeV)");
     h_OfflineMET_NewRechitM0_Et->GetYaxis()->SetTitle("New Geometry RechitM0 (GeV)");

     NewRechitM0_Et_BE = sqrt(rh_Ex*rh_Ex+rh_Ey*rh_Ey);
     h_CaloMET_NewRechitM0_BE->Fill(CaloMetBE.et(),NewRechitM0_Et_BE);
     h_CaloMET_NewRechitM0_BE->SetMarkerStyle(3);
     h_CaloMET_NewRechitM0_BE->SetMarkerSize(1);
     h_CaloMET_NewRechitM0_BE->GetXaxis()->SetTitle("CaloMET (GeV)");
     h_CaloMET_NewRechitM0_BE->GetYaxis()->SetTitle("New Geometry RechitM0 BE (GeV)");

     h_NewRechitM0_BE->Fill(NewRechitM0_Et_BE);
   }

   for (const auto& pair: rhits) {
      ev_rh_unmatched_ += pair.second.size();
   }

   events_->Fill();
}

void
HcalCompareLegacyChains::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalCompareLegacyChains);
