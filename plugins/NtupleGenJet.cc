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


#include "Math/VectorUtil.h"
#include "Math/Point3D.h"
#include "Math/GenVector/CoordinateSystemTags.h"
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> > Point3D;

#include "Math/PxPyPzE4D.h"
#include "Math/LorentzVector.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > MyLorentzVector;

using namespace std;
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
	
        int genHiggs_n_=1;
  int noOfGenParticle = 0;
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
  TH1D *nTPrimeId;
  MyLorentzVector genLead, genSublead;
  TH1F *invMass;
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
   nTPrimeId = fs->make<TH1D>("N_TPrimeId" , ";T/T' Id;Events;;" , 20 , -10 , 10 );
   nHiggs_histo = fs->make<TH1D>("N_Tprime" , ";N_{T/T'};Events;;" , 5 , 0 , 5 );
   pt_histo = fs->make<TH1F>("pT_Tprime" , ";p_{T} of Tprime[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo=fs->make<TH1F>("eta_Tprime" , ";#eta of Tprime;Events;;" , 50 , -5 , 5 );
   phi_histo=fs->make<TH1F>("phi_Tprime" , ";#phi of Tprime;Events;;" , 50 , -5 , 5 );
   mass_histo=fs->make<TH1F>("mass_Tprime" , ";mass of Tprime;Events;;" , 100,200,700);
   pt_histo_lead_b = fs->make<TH1F>("pT_b1" , ";p_{T} of leading b[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo_lead_b=fs->make<TH1F>("eta_b1" , ";#eta of leading b;Events;;" , 50 , -5 , 5 );
   phi_histo_lead_b=fs->make<TH1F>("phi_b1" , ";#phi of leading b;Events;;" , 50 , -5 , 5 );
   pt_histo_sublead_b = fs->make<TH1F>("pT_b2" , ";p_{T} of sub-leading b[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo_sublead_b=fs->make<TH1F>("eta_b2" , ";#eta of sub-leading b;Events;;" , 50 , -5 , 5 );
   phi_histo_sublead_b=fs->make<TH1F>("phi_b2" , ";#phi of sub-leading b;Events;;" , 50 , -5 , 5 );
   pt_histo_add_b = fs->make<TH1F>("pT_b" , ";p_{T} of additional b[GeV];Events;;" , 100 , 0 , 500 );
   eta_histo_add_b=fs->make<TH1F>("eta_b" , ";#eta of additional b;Events;;" , 50 , -5 , 5 );
   phi_histo_add_b=fs->make<TH1F>("phi_b" , ";#phi of additional b;Events;;" , 50 , -5 , 5 );
   invMass=fs->make<TH1F>("invMass",";mass;Events;;",500,200.,700.);
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
  noOfGenParticle++;
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genparticlesToken, genParticles);
  //for(reco::GenParticle jet : *(gen_h.product())){
  // for(const auto& jet : genparticles){ 
  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const reco::GenParticle & p = (*genParticles)[i];
    int id = p.pdgId();
    //     noOfGenParticle++;
    //int st = p.status();  
    //     double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
    //std::cout << "pdg id =  "<< id << std::endl;
    if (id == 8 || id == -8){ // check if it is T/Tprime
      //nTPrimeId->Fill(id);
      int n = p.numberOfDaughters();
      if(n < 2 ) continue;
      nTPrimeId->Fill(id);
      //std::cout << "number of daughter:  " << n << std::endl;
      const reco::Candidate * d1 = p.daughter( 0 );
      const reco::Candidate * d2 = p.daughter( 1 );
      //std::cout << "pdg id of d1=  " << d1->pdgId() << " pdg id of d2=  " << d2->pdgId() << std::endl;
      
      if ((std::abs(d1->pdgId())==6 && std::abs(d2->pdgId())==21) || (std::abs(d1->pdgId())==21 && std::abs(d2->pdgId())==6)) // 
	{
	  const reco::Candidate * mom = p.mother();
	  if (std::abs(mom->pdgId()) != 8) continue;

	  nHiggs_histo->Fill(genHiggs_n_); // mind it will give a histogram at "1" for every entries
	  pt_histo->Fill(p.pt());
	  eta_histo->Fill(p.eta());
	  phi_histo ->Fill(p.phi());
	  mass_histo->Fill(p.mass());
	  //// plotting for all daughters perperties
	  
	  pt_histo_lead_b->Fill(d1->pt()); 
	  eta_histo_lead_b->Fill(d1->eta());
	  phi_histo_lead_b->Fill(d1->phi());
	  
	  pt_histo_sublead_b->Fill(d2->pt());
	  eta_histo_sublead_b->Fill(d2->eta());
	  phi_histo_sublead_b->Fill(d2->phi());
	  
	  // calculate invariant mass
	  genLead.SetPxPyPzE(d1->px(),d1->py(),d1->pz(),d1->energy());
	  genSublead.SetPxPyPzE(d2->px(),d2->py(),d2->pz(),d2->energy());
	  MyLorentzVector diElectron = genLead + genSublead;
	  invMass->Fill(diElectron.M());
	  
	}
    }
     /*
     if(id == 5 || id == -5){
       const reco::Candidate * mom = p.mother();
       if (mom->pdgId()!=25){
	 pt_histo_add_b->Fill(p.pt());
	 eta_histo_add_b->Fill(p.eta());
	 phi_histo_add_b->Fill(p.phi());
       }
     }
     */
   }
  //   std::cout << "noOfGenParticle: " << noOfGenParticle << endl;
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
