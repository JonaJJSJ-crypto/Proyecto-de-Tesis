// -*- C++ -*-
//
// Package:    AcausalAnalyzer
// Class:      AcausalAnalyzer
// 
/**\class AcausalAnalyzer AcausalAnalyzer.cc Proyecto-de-Tesis/AcausalAnalyzer/src/AcausalAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Wed May 19 05:35:26 CEST 2021
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "math.h"

#include "TFile.h"
#include "TTree.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/ParameterSet/interface/Registry.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Utilities/interface/RegexMatch.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/print.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include <iostream>
#include <cassert>
#include <vector>
#include <string>
#include <iomanip>
#include <boost/foreach.hpp>

const static std::vector<std::string> interestingTriggers = {
/*      "HLT_HT400_Ele60_CaloIdL_TkrIdt",
      "HLT_HT450_Ele60_CaloIdL_TkrIdt",*/
      "HLT_DoubleEle45_CaloIdL",
      "HLT_Ele8",
};

template <typename T>
void subtractInvisible(T g, reco::Candidate::LorentzVector& p4) {
  auto daughters = (*g).daughterRefVector();
  for (auto d = daughters.begin(); d != daughters.end(); d++) {
    const auto pdgId = (*d)->pdgId();
    if (std::abs(pdgId) == 12 || std::abs(pdgId) == 14 ||
        std::abs(pdgId) == 16 || std::abs(pdgId) == 18) {
      p4 = p4 - (*d)->p4();
    }
    subtractInvisible(*d, p4);
  }
}

template <typename T>
int findBestVisibleMatch(T& gens, reco::Candidate::LorentzVector& p4) {
  float minDeltaR = 999.0;
  int idx = -1;
  for (auto g = gens.begin(); g != gens.end(); g++) {
    auto tmp_p4 = g->p4();
    subtractInvisible(g, tmp_p4);
    const auto tmp = deltaR(tmp_p4, p4);
    if (tmp < minDeltaR) {
      minDeltaR = tmp;
      idx = g - gens.begin();
    }
  }
  return idx;
}

template <typename T>
int findBestMatch(T& gens, reco::Candidate::LorentzVector& p4) {
  float minDeltaR = 999.0;
  int idx = -1;
  for (auto g = gens.begin(); g != gens.end(); g++) {
    const auto tmp = deltaR(g->p4(), p4);
    if (tmp < minDeltaR) {
      minDeltaR = tmp;
      idx = g - gens.begin();
    }
  }
  return idx;
}
//
// class declaration
//

class AcausalAnalyzer : public edm::EDAnalyzer {
   public:
      explicit AcausalAnalyzer(const edm::ParameterSet&);
      ~AcausalAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
//tracks
      edm::InputTag trackTags_;


      // from HLTEventAnalyzerAOD.h
      /// module config parameters
      std::string   processName_;
      edm::InputTag triggerResultsTag_;
      edm::InputTag triggerEventTag_;
    /// HLT trigger names
    edm::ParameterSetID triggerNamesID_;

      // additional class data memebers
      // these are actually the containers where we will store
      // the trigger information
      edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
      edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
      HLTConfigProvider hltConfig_;

      //inspired by https://github.com/cms-sw/cmssw/blob/CMSSW_5_3_X/HLTrigger/HLTfilters/interface/HLTHighLevel.h
      // input patterns that will be expanded into trigger names
      std::vector<std::string>  HLTPatterns_;

    /// list of required HLT triggers by HLT name
      std::vector<std::string>  HLTPathsByName_;


      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
	
     bool cmsStandardCuts(const edm::Event&, const edm::EventSetup&);
     bool matchingCuts( bool , double , int , int , double, double, double);
     double deltaR(double , double , double, double);
     double conePt(int , int , double , double , int ,const edm::Event& , const edm::EventSetup& );
     double mCos(double , double , double , double  );
     double mTheta(double , double , double , double );
     double invMass(double , double , double , double  , double ,  double );
     double dotProduct(double , double , double , double );
     bool impactParameterCut(reco::TrackCollection::const_iterator, reco::TrackCollection::const_iterator, reco::BeamSpot );
      bool providesGoodLumisection(const edm::Event &iEvent);
      bool isData;

      TTree *tree;
      TFile *fs;
      // ----------member data ---------------------------
   // Event information
  Int_t value_run;
  UInt_t value_lumi_block;
  ULong64_t value_event;

  // Trigger
  const static int max_trig = 1000;
  int value_trig[max_trig];
  float prescale_trig[max_trig];

  // Vertices
  int value_ve_n;
  float value_ve_x;
  float value_ve_y;
  float value_ve_z;


  // Electrons
  const static int max_el = 1000;
  UInt_t value_el_n;
  float value_el_pt[max_el];
  float value_el_eta[max_el];
  float value_el_phi[max_el];
  float value_el_mass[max_el];
  int value_el_charge[max_el];
  float value_el_pfreliso03all[max_el];
  float value_el_dxy[max_el];
  float value_el_dxyErr[max_el];
  float value_el_dz[max_el];
  float value_el_dzErr[max_el];
  bool value_el_cutbasedid[max_el];
  bool value_el_pfid[max_el];
  int value_el_genpartidx[max_el];
  int value_el_jetidx[max_el];

    // Jets
  const static int max_jet = 1000;
  UInt_t value_jet_n;
  float value_jet_pt[max_jet];
  float value_jet_eta[max_jet];
  float value_jet_phi[max_jet];
  float value_jet_mass[max_jet];
  bool value_jet_puid[max_jet];
  float value_jet_btag[max_jet];


// Generator particles
  const static int max_gen = 1000;
  UInt_t value_gen_n;
  float value_gen_pt[max_gen];
  float value_gen_eta[max_gen];
  float value_gen_phi[max_gen];
  float value_gen_mass[max_gen];
  int value_gen_pdgid[max_gen];
  int value_gen_motherid[max_gen];
  int value_gen_status[max_gen];
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
AcausalAnalyzer::AcausalAnalyzer(const edm::ParameterSet &iConfig):
trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks")),
processName_(iConfig.getParameter<std::string>("processName")),
triggerResultsTag_(iConfig.getParameter<edm::InputTag>("triggerResults")),
triggerEventTag_(iConfig.getParameter<edm::InputTag>("triggerEvent")),
triggerNamesID_(),
HLTPatterns_(iConfig.getParameter<std::vector<std::string> >("triggerPatterns")),
HLTPathsByName_(),
isData(iConfig.getParameter<bool>("isData"))
{
   //now do what ever initialization is needed
  fs = new TFile("EleInfo.root","RECREATE");
  tree = new TTree("tree","EleInfo");


  // Event information
  tree->Branch("run", &value_run);
  tree->Branch("luminosityBlock", &value_lumi_block);
  tree->Branch("event", &value_event);

  // Trigger
  for(size_t i = 0; i < interestingTriggers.size(); i++) {
    tree->Branch(interestingTriggers[i].c_str(), value_trig + i, (interestingTriggers[i] + "/i").c_str());
    tree->Branch((interestingTriggers[i]+"_ps").c_str(), prescale_trig + i, (interestingTriggers[i]+"_ps"+"/F").c_str());
  }

  // Vertices
  tree->Branch("PV_npvs", &value_ve_n, "PV_npvs/I");
  tree->Branch("PV_x", &value_ve_x, "PV_x/F");
  tree->Branch("PV_y", &value_ve_y, "PV_y/F");
  tree->Branch("PV_z", &value_ve_z, "PV_z/F");

 // Electrons
  tree->Branch("nElectron", &value_el_n, "nElectron/i");
  tree->Branch("Electron_pt", value_el_pt, "Electron_pt[nElectron]/F");
  tree->Branch("Electron_eta", value_el_eta, "Electron_eta[nElectron]/F");
  tree->Branch("Electron_phi", value_el_phi, "Electron_phi[nElectron]/F");
  tree->Branch("Electron_mass", value_el_mass, "Electron_mass[nElectron]/F");
  tree->Branch("Electron_charge", value_el_charge, "Electron_charge[nElectron]/I");
  tree->Branch("Electron_pfRelIso03_all", value_el_pfreliso03all, "Electron_pfRelIso03_all[nElectron]/F");
  tree->Branch("Electron_dxy", value_el_dxy, "Electron_dxy[nElectron]/F");
  tree->Branch("Electron_dxyErr", value_el_dxyErr, "Electron_dxyErr[nElectron]/F");
  tree->Branch("Electron_dz", value_el_dz, "Electron_dz[nElectron]/F");
  tree->Branch("Electron_dzErr", value_el_dzErr, "Electron_dzErr[nElectron]/F");
  tree->Branch("Electron_cutBasedId", value_el_cutbasedid, "Electron_cutBasedId[nElectron]/O");
  tree->Branch("Electron_pfId", value_el_pfid, "Electron_pfId[nElectron]/O");
  tree->Branch("Electron_jetIdx", value_el_jetidx, "Electron_jetIdx[nElectron]/I");
  tree->Branch("Electron_genPartIdx", value_el_genpartidx, "Electron_genPartIdx[nElectron]/I");

  // Jets
  tree->Branch("nJet", &value_jet_n, "nJet/i");
  tree->Branch("Jet_pt", value_jet_pt, "Jet_pt[nJet]/F");
  tree->Branch("Jet_eta", value_jet_eta, "Jet_eta[nJet]/F");
  tree->Branch("Jet_phi", value_jet_phi, "Jet_phi[nJet]/F");
  tree->Branch("Jet_mass", value_jet_mass, "Jet_mass[nJet]/F");
  tree->Branch("Jet_puId", value_jet_puid, "Jet_puId[nJet]/O");
  tree->Branch("Jet_btag", value_jet_btag, "Jet_btag[nJet]/F");

  // Generator particles
  if (!isData) {
    tree->Branch("nGenPart", &value_gen_n, "nGenPart/i");
    tree->Branch("GenPart_pt", value_gen_pt, "GenPart_pt[nGenPart]/F");
    tree->Branch("GenPart_eta", value_gen_eta, "GenPart_eta[nGenPart]/F");
    tree->Branch("GenPart_phi", value_gen_phi, "GenPart_phi[nGenPart]/F");
    tree->Branch("GenPart_mass", value_gen_mass, "GenPart_mass[nGenPart]/F");
    tree->Branch("GenPart_pdgId", value_gen_pdgid, "GenPart_pdgId[nGenPart]/I");
    tree->Branch("GenPart_motherId", value_gen_motherid, "GenPart_motherId[nGenPart]/I");
    tree->Branch("GenPart_status", value_gen_status, "GenPart_status[nGenPart]/I");
  }

}


AcausalAnalyzer::~AcausalAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void AcausalAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;


  // Event information
  value_run = iEvent.run();
  value_lumi_block = iEvent.luminosityBlock();
  value_event = iEvent.id().event();

  // Trigger results
  Handle<TriggerResults> trigger;
  iEvent.getByLabel(InputTag("TriggerResults", "", "HLT"), trigger);
  auto psetRegistry = edm::pset::Registry::instance();
  auto triggerParams = psetRegistry->getMapped(trigger->parameterSetID());
  TriggerNames triggerNames(*triggerParams);
  TriggerResultsByName triggerByName(&(*trigger), &triggerNames);


  for (size_t i = 0; i < interestingTriggers.size(); i++) {
    value_trig[i] = 0;
    prescale_trig[i] = 0;
  }
  const auto names = triggerByName.triggerNames();
  for (size_t i = 0; i < names.size(); i++) {
    const auto name = names[i];
    for (size_t j = 0; j < interestingTriggers.size(); j++) {
      const auto interest = interestingTriggers[j];
      if (name.find(interest) == 0) {
        const auto substr = name.substr(interest.length(), 2);
        if (substr.compare("_v") == 0) {
          const auto status = triggerByName.state(name);

            cout<<"status "<<status<<endl;
            cout<<"nombre "<<names[i]<<endl;
            auto pathName = names[i];
            const std::pair<int,int> prescales(hltConfig_.prescaleValues(iEvent,iSetup,pathName));
            cout<<"pre  "<<prescales.first<<","<<prescales.second<<endl;      
            prescale_trig[j] = prescales.first/prescales.second;
            value_trig[j] = status;

          if (status == 1) {
            value_trig[j] = 1;
            cout<<"status "<<status<<endl;
            cout<<"nombre "<<names[i]<<endl;
            break;
          }
        }
      }
    }
  }

//Tracks
  Handle<TrackCollection> tracks;
  iEvent.getByLabel(trackTags_,tracks);
/*  int matchedTrack[tracks->size()];
  for (int i = 0; i<(int)tracks->size(); i++)
  {
	matchedTrack[i] = 0;
  }*/
  cout<<"number tracks 	"<<tracks->size()<<endl; 
  for(TrackCollection::const_iterator itTrack = tracks->begin();
      itTrack != tracks->end();
      ++itTrack)
  {
    if(itTrack->pt() > 5)
    {
     cout<<"vert "<<itTrack->pt()<<endl;
    }
  }

  // Vertex
  Handle<VertexCollection> vertices;
  iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);
  value_ve_n = vertices->size();
  value_ve_x = vertices->begin()->x();
  value_ve_y = vertices->begin()->y();
  value_ve_z = vertices->begin()->z();
  math::XYZPoint pv(vertices->begin()->position());

  // Electrons
  Handle<GsfElectronCollection> electrons;
  iEvent.getByLabel(InputTag("gsfElectrons"), electrons);

  value_el_n = 0;
  const float el_min_pt = 5;
  std::vector<GsfElectron> selectedElectrons;
  for (auto it = electrons->begin(); it != electrons->end(); it++) {
    if (it->pt() > el_min_pt) {
      std::cout<<it->pt()<<"================================AQUI===========================\n ";
      selectedElectrons.emplace_back(*it);
      value_el_pt[value_el_n] = it->pt();
      value_el_eta[value_el_n] = it->eta();
      value_el_phi[value_el_n] = it->phi();
      value_el_charge[value_el_n] = it->charge();
      value_el_mass[value_el_n] = it->mass();
      value_el_cutbasedid[value_el_n] = it->passingCutBasedPreselection();
      value_el_pfid[value_el_n] = it->passingPflowPreselection();
      if (it->passingPflowPreselection()) {
        auto iso03 = it->pfIsolationVariables();
        value_el_pfreliso03all[value_el_n] =
            (iso03.chargedHadronIso + iso03.neutralHadronIso + iso03.photonIso)/it->pt();
      } else {
        value_el_pfreliso03all[value_el_n] = -999;
      }
      auto trk = it->gsfTrack();
      value_el_dxy[value_el_n] = trk->dxy(pv);
      value_el_dz[value_el_n] = trk->dz(pv);
      value_el_dxyErr[value_el_n] = trk->d0Error();
      value_el_dzErr[value_el_n] = trk->dzError();
      value_el_jetidx[value_el_n] = -1;
      value_el_genpartidx[value_el_n] = -1;
      value_el_n++;
    }
  }

  // Jets
  // Jet ID recommendations:
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_8_TeV_data_a
  // B-tag recommendations:
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation53XReReco
  Handle<CaloJetCollection> jets;
  iEvent.getByLabel(InputTag("ak5CaloJets"), jets);
  Handle<JetTagCollection> btags;
  iEvent.getByLabel(InputTag("combinedSecondaryVertexBJetTags"), btags);

  const float jet_min_pt = 15;
  value_jet_n = 0;
  std::vector<CaloJet> selectedJets;
  for (auto it = jets->begin(); it != jets->end(); it++) {
    if (it->pt() > jet_min_pt) {
      selectedJets.emplace_back(*it);
      value_jet_pt[value_jet_n] = it->pt();
      value_jet_eta[value_jet_n] = it->eta();
      value_jet_phi[value_jet_n] = it->phi();
      value_jet_mass[value_jet_n] = it->mass();
      value_jet_puid[value_jet_n] = it->emEnergyFraction() > 0.01 && it->n90() > 1;
      value_jet_btag[value_jet_n] = btags->operator[](it - jets->begin()).second;
      value_jet_n++;
    }
  }

  if (!isData) {
    Handle<GenParticleCollection> gens;
    iEvent.getByLabel(InputTag("genParticles"), gens);

    value_gen_n = 0;
    std::vector<GenParticle> interestingGenParticles;
    for (auto it = gens->begin(); it != gens->end(); it++) {
      const auto status = it->status();
      const auto pdgId = std::abs(it->pdgId());
      if (status == 1 && pdgId == 13) { // muon
        interestingGenParticles.emplace_back(*it);
      }
      if (status == 1 && pdgId == 11) { // electron
        interestingGenParticles.emplace_back(*it);
      }
      /*
      if (status == 1 && pdgId == 22) { // photon
        interestingGenParticles.emplace_back(*it);
      }
      */
      if (status == 2 && pdgId == 15) { // tau
        interestingGenParticles.emplace_back(*it);
      }
    }


    // Match electrons with gen particles
    for (auto p = selectedElectrons.begin(); p != selectedElectrons.end(); p++) {
      // Gen particle matching
      auto p4 = p->p4();
      auto idx = findBestVisibleMatch(interestingGenParticles, p4);
      if (idx != -1) {
        auto g = interestingGenParticles.begin() + idx;
        value_gen_pt[value_gen_n] = g->pt();
        value_gen_eta[value_gen_n] = g->eta();
        value_gen_phi[value_gen_n] = g->phi();
        value_gen_mass[value_gen_n] = g->mass();
        value_gen_pdgid[value_gen_n] = g->pdgId();
        value_gen_motherid[value_gen_n] = g->mother()->pdgId();
	//Mother
	auto motherRefs = g->motherRefVector();
        for(reco::GenParticleRefVector::const_iterator imr = motherRefs.begin(); imr!= motherRefs.end(); ++imr) {
	 cout<<"   - Mother "<<(*imr).key()<<" "<<(*imr)->pdgId()<<endl;
       }
        value_gen_status[value_gen_n] = g->status();
        value_el_genpartidx[p - selectedElectrons.begin()] = value_gen_n;
        value_gen_n++;
      }
    }
  } // !isData

  // Fill event
  tree->Fill();
  return;
}


// ------------ method called once each job just before starting event loop  ------------
void 
AcausalAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AcausalAnalyzer::endJob() 
{
    fs->Write();
}

// ------------ method called when starting to processes a run  ------------

bool

AcausalAnalyzer::cmsStandardCuts(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 bool ret = false;
 int disp = 0;
 int numHighPurity = 0;
 int numTracks = 0;

 double purityFrac = 0.0;
 using namespace edm;
 using namespace reco;
   Handle<reco::VertexCollection> vertHand;
   iEvent.getByLabel( "offlinePrimaryVertices",vertHand);
   
   Handle<TrackCollection> tracks;
   iEvent.getByLabel(trackTags_,tracks);
   
   for(reco::VertexCollection::const_iterator itVert = vertHand->begin();
       itVert != vertHand->begin()+1 && itVert != vertHand->end();
       ++itVert)
       {
		  
		    for(reco::Vertex::trackRef_iterator itTrack = itVert->tracks_begin();
       itTrack != itVert->tracks_end();
       ++itTrack)
           {
		        
		 double tx = (**itTrack).vx();
	     double ty = (**itTrack).vy();
		 double tz = (**itTrack).vz();
		 double tr = sqrt(tx*tx +ty*ty);
		 if ( (tr < 2 )&&(tz < 24) )
		 {
			 disp ++;
		 }
 
    	   }
		    
	   }
	 
	   
	   for(TrackCollection::const_iterator itTrack = tracks->begin();
       itTrack != tracks->end();                      
       ++itTrack) 
       {   
		  
		   numTracks ++;
		   
		   if (itTrack->quality(reco::Track::highPurity))
		   {
			   numHighPurity++;
		   }
		   
	   }
	   
	   purityFrac = (double)numHighPurity/numTracks;
	   if (numTracks<10)
	   {
		   purityFrac = 1;
	   }
	   
	   if(disp >=4 && purityFrac >= 0.25)
	   {
		   ret = true;
	   }
	   
	
}	
	return ret;
	

}




bool 
{
	bool ret = false;
	
	
		
	  if(purity && pt > 41 && hits >= 6   && eta < 2  && hits3D >1)
	  if(true)
	  
	  {
		  ret = true;

	  }	
	 
	  

	
	return ret;
}
double 
AcausalAnalyzer::deltaR(double obj1Phi, double obj1Eta, double obj2Phi, double obj2Eta)
{
	double dPhi = obj1Phi - obj2Phi;
	if (abs(dPhi)>3.1415/2)
	{
		dPhi =  3.1415 -dPhi;
	}
	double dEta = obj1Eta - obj2Eta;
	double dR = sqrt(dPhi*dPhi + dEta*dEta);
	return dR;
}
double 
AcausalAnalyzer::conePt(int forbiddenIndex1, int forbiddenIndex2, double eta, double phi, int numTracks,const edm::Event& iEvent, const edm::EventSetup& iSetup )
{  
	using namespace edm;
	using namespace reco;
    Handle<TrackCollection> tracks;
   iEvent.getByLabel(trackTags_,tracks);
    int i =0;
	double sumPt = 0.0;
	for(TrackCollection::const_iterator itTrack = tracks->begin();
       itTrack != tracks->end();                      
       ++itTrack) 
	{
		if(deltaR(phi, eta, itTrack->phi(), itTrack->eta()) < 0.3 && deltaR(phi, eta, itTrack->phi(), itTrack->eta()) > 0.03 && i != forbiddenIndex1 && i != forbiddenIndex2)
		{
			sumPt = itTrack->pt() +sumPt;
			i++;
		}
	}
	
	return sumPt;
}
double 
AcausalAnalyzer::mCos(double phi1, double eta1, double phi2, double eta2 )
{double cosAlpha = cos(eta1)*cos(phi1)*cos(eta2)*cos(phi2) + cos(eta1)*sin(phi1)*cos(eta2)*sin(phi2) + sin(eta1)*cos(eta2);
	
	return cosAlpha;
}
double 
AcausalAnalyzer::mTheta(double ax, double ay, double bx, double by)
{
	double cosAlpha = ax*bx + ay*by;
	double theta;
	cosAlpha = cosAlpha/(sqrt(ax*ax+ay*ay)*sqrt(bx*bx+by*by));
	//std::cout<<"cosAlpha: "<<cosAlpha<<std::endl;
	theta  = acos(cosAlpha);
	if (theta > 3.1514)
	{
		theta =  2*3.1514-theta;
	}
	return theta;
}
double 
AcausalAnalyzer::invMass(double px1, double py1, double pz1, double px2 , double py2,  double pz2)
{
  double E1 =  sqrt(px1*px1 + py1*py1 + pz1*pz1);  // asummes rest mass energy to be negligible
  double E2 =  sqrt(px2*px2 + py2*py2 + pz2*pz2);
  
  double E = E1+E2;
  double px = px1 + px2;
  double py = py1 + py2;
  double pz = pz1 + pz2;
  double mass = sqrt(E*E - px*px - py*py - pz*pz);
  
  return mass;
  
  

	
}
double
AcausalAnalyzer::dotProduct(double x1, double y1 , double x2, double y2)
{
	double ret;
	ret = x1*x2+y1*y2;
	return ret;
}
bool
AcausalAnalyzer::impactParameterCut(reco::TrackCollection::const_iterator it1, reco::TrackCollection::const_iterator it2, reco::BeamSpot beamSpot)
{    
	double dxy1, dxy1Err, dxy2, dxy2Err, sig1, sig2;
	bool ret = false;
	dxy1 = it1->dxy(beamSpot);
	dxy2 = it2->dxy(beamSpot);
	dxy1Err = it1->dxyError();
	dxy2Err = it2->dxyError();
	
	sig1 = dxy1/dxy1Err;
	sig2 = dxy2/dxy2Err;
	
	if (sig1 > 2 && sig2 > 2)
	{
		ret = true;
	}
	
	
	return ret;
}


void 
AcausalAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{

    using namespace std;
    using namespace edm;
    processName_="HLT";
    bool changed(true);
    hltConfig_.init(iRun,iSetup,processName_,changed);
    
    if (changed)
    {
        cout<<"HLTConfig has changed . . . "<<endl;
        
    }

}

// ------------ method called when ending the processing of a run  ------------
void 
AcausalAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
AcausalAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
AcausalAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AcausalAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AcausalAnalyzer);
