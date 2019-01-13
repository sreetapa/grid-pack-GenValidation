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
       edm::EDGetTokenT <reco::GenParticleCollection> genparticlesToken;
	
        int genHiggs_n_=0;
        TH1D *nHiggs_histo;
        TH1F *pt_histo;
        TH1F *eta_histo;
	TH1F *phi_histo;  
	TH1F *mass_histo;
	TH1F *pt_histo_lead_b;
        TH1F *eta_histo_lead_b;
        TH1F *phi_histo_lead_b;
	TH1F *pt_histo_sublead_b;
        TH1F *eta_histo_sublead_b;
        TH1F *phi_histo_sublead_b;
	TH1F *pt_histo_add_b;
        TH1F *eta_histo_add_b;
        TH1F *phi_histo_add_b;
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
   genparticlesToken	   = consumes <reco::GenParticleCollection> (std::string("genParticles"));
   edm::Service<TFileService> fs;
   nHiggs_histo = fs->make<TH1D>("N_higgs" , ";N_{H};Events;;" , 50 , 0 , 50 );
   pt_histo = fs->make<TH1F>("pT_H" , ";p_{T} of Higgs[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo=fs->make<TH1F>("eta_H" , ";#eta of Higgs;Events;;" , 50 , -5 , 5 );
   phi_histo=fs->make<TH1F>("phi_H" , ";#phi of Higgs;Events;;" , 50 , -5 , 5 );
   mass_histo=fs->make<TH1F>("mass_H" , ";mass of Higgs;Events;;" , 100,0,500);
   pt_histo_lead_b = fs->make<TH1F>("pT_b1" , ";p_{T} of leading b[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo_lead_b=fs->make<TH1F>("eta_b1" , ";#eta of leading b;Events;;" , 50 , -5 , 5 );
   phi_histo_lead_b=fs->make<TH1F>("phi_b1" , ";#phi of leading b;Events;;" , 50 , -5 , 5 );
   pt_histo_sublead_b = fs->make<TH1F>("pT_b2" , ";p_{T} of sub-leading b[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo_sublead_b=fs->make<TH1F>("eta_b2" , ";#eta of sub-leading b;Events;;" , 50 , -5 , 5 );
   phi_histo_sublead_b=fs->make<TH1F>("phi_b2" , ";#phi of sub-leading b;Events;;" , 50 , -5 , 5 );
   pt_histo_add_b = fs->make<TH1F>("pT_b" , ";p_{T} of additional b[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo_add_b=fs->make<TH1F>("eta_b" , ";#eta of additional b;Events;;" , 50 , -5 , 5 );
   phi_histo_add_b=fs->make<TH1F>("phi_b" , ";#phi of additional b;Events;;" , 50 , -5 , 5 );
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

    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(genparticlesToken, genParticles);
//for(reco::GenParticle jet : *(gen_h.product())){
 // for(const auto& jet : genparticles){ 
   for(size_t i = 0; i < genParticles->size(); ++ i) {
     const reco::GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
     //int st = p.status();  
     //     double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
     //std::cout << "pdg id =  "<< id << std::endl;
     if (id == 25){ // check if it is H
       int n = p.numberOfDaughters();
       if(n < 2 ) continue;
       //std::cout << "number of daughter:  " << n << std::endl;
       const reco::Candidate * d1 = p.daughter( 0 );
       const reco::Candidate * d2 = p.daughter( 1 );
       //std::cout << "pdg id of d1=  " << d1->pdgId() << " pdg id of d2=  " << d2->pdgId() << std::endl;
       
       if (std::abs(d1->pdgId())==15 && std::abs(d2->pdgId())==15){ // check when H decays to tau tau
	 ++genHiggs_n_;
	 nHiggs_histo->Fill(genHiggs_n_); // mind it will give a flat histogram at "1" for every entries
	 pt_histo->Fill(p.pt());
	 eta_histo->Fill(p.eta());
	 phi_histo ->Fill(p.phi());
	 mass_histo->Fill(p.mass());
	 //// plotting for b's
	 pt_histo_lead_b->Fill(d1->pt());
	 eta_histo_lead_b->Fill(d1->eta());
	 phi_histo_lead_b->Fill(d1->phi());
	 
	 pt_histo_sublead_b->Fill(d2->pt());
	 eta_histo_sublead_b->Fill(d2->eta());
	 phi_histo_sublead_b->Fill(d2->phi());
       }
     }
     
     if(id == 5 || id == -5){
       const reco::Candidate * mom = p.mother();
       if (mom->pdgId()!=25){
	 pt_histo_add_b->Fill(p.pt());
	 eta_histo_add_b->Fill(p.eta());
	 phi_histo_add_b->Fill(p.phi());
       }
     }
   }
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
