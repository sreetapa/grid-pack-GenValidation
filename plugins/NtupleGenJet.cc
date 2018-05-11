// -*- C++ -*-
//
// Package:    validation/NtupleGenJet
// Class:      NtupleGenJet
// 
/**\class NtupleGenJet NtupleGenJet.cc validation/NtupleGenJet/plugins/NtupleGenJet.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lata Panwar
//         Created:  Wed, 14 Mar 2018 12:41:16 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "TTree.h"
#include "TH1.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class NtupleGenJet : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit NtupleGenJet(const edm::ParameterSet&);
      ~NtupleGenJet();
      
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      // ----------member data ---------------------------
       edm::EDGetTokenT <reco::GenJetCollection> genJetsToken;
	
        int genjet_n_;
        TH1D *njet_histo;
        TH1F *pt_histo_lead;
        TH1F *eta_histo_lead;
  	TH1F *phi_histo_lead;
	TH1F *pt_histo_sublead;
        TH1F *eta_histo_sublead;
        TH1F *phi_histo_sublead;
	TH1F *pt_histo_3rdlead;
        TH1F *eta_histo_3rdlead;
        TH1F *phi_histo_3rdlead;
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
NtupleGenJet::NtupleGenJet(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   genJetsToken	   = consumes <reco::GenJetCollection> (std::string("ak4GenJets"));
   edm::Service<TFileService> fs;
   njet_histo = fs->make<TH1D>("Njet" , ";Njet;Events;;" , 50 , 0 , 50 );
   pt_histo_lead = fs->make<TH1F>("pT1" , ";p_{T}leading AK4Jet[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo_lead=fs->make<TH1F>("eta1" , ";#eta leading AK4Jet;Events;;" , 50 , -5 , 5 );
   phi_histo_lead=fs->make<TH1F>("phi1" , ";#phi leading AK4Jet;Events;;" , 50 , -5 , 5 );
   pt_histo_sublead = fs->make<TH1F>("pT2" , ";p_{T}subleading AK4Jet[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo_sublead=fs->make<TH1F>("eta2" , ";#eta subleading AK4Jet;Events;;" , 50 , -5 , 5 );
   phi_histo_sublead=fs->make<TH1F>("phi2" , ";#phi subleading AK4Jet;Events;;" , 50 , -5 , 5 );
   pt_histo_3rdlead = fs->make<TH1F>("pT3" , ";p_{T}3rd leading AK4Jet[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo_3rdlead=fs->make<TH1F>("eta3" , ";#eta 3rd leading AK4Jet;Events;;" , 50 , -5 , 5 );
   phi_histo_3rdlead=fs->make<TH1F>("phi3" , ";#phi 3rd leading AK4Jet;Events;;" , 50 , -5 , 5 );

}



NtupleGenJet::~NtupleGenJet()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
NtupleGenJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    edm::Handle< reco::GenJetCollection > genjets_h;
    iEvent.getByToken(genJetsToken, genjets_h);
    const reco::GenJetCollection& genjets = *genjets_h;
   

    genjet_n_ = genjets.size();
    
    njet_histo->Fill(genjet_n_);

    pt_histo_lead->Fill(genjets.at(0).pt());
    eta_histo_lead->Fill(genjets.at(0).eta());
    phi_histo_lead->Fill(genjets.at(0).phi());
    
    pt_histo_sublead->Fill(genjets.at(1).pt());
    eta_histo_sublead->Fill(genjets.at(1).eta());
    phi_histo_sublead->Fill(genjets.at(1).phi());
  
    pt_histo_3rdlead->Fill(genjets.at(2).pt());
    eta_histo_3rdlead->Fill(genjets.at(2).eta());
    phi_histo_3rdlead->Fill(genjets.at(2).phi());

    //std::cout << genjet_n_<<std::endl;
    //for(const auto& jet : genjets.at(0)){
       //std::cout <<"pt="<< jet.pt();
      // pt_histo->Fill(jet.pt());
      // eta_histo->Fill(jet.eta());
   //}
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NtupleGenJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NtupleGenJet);
