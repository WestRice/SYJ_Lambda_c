#include "SYJ_Lambda_c/Sigma0PionEta.h"
#include "PDG.h"
#include "Cut.h"
#include <cmath>

#include "McDecayModeSvc/McDecayModeSvc.h"
#include "BestDTagSvc/BestDTagSvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EvtRecEvent/EvtRecDTag.h"
#include "EvtRecEvent/EvtRecVeeVertex.h"
#include "EvtRecEvent/EvtRecPi0.h"
#include "DstEvent/TofHitStatus.h"
#include "MdcRecEvent/RecMdcKalTrack.h"

#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"

#include "VertexFit/Helix.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KalmanKinematicFit.h"
//#include "SimplePIDSvc/ISimplePIDSvc.h"

#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "ParticleID/ParticleID.h"

#include "McTruth/McParticle.h"
#include "EvTimeEvent/RecEsTime.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;

#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include <vector>
typedef std::vector<int> Vint;
typedef std::vector<double> Vdou;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<WTrackParameter> VWTrkPara;
enum PARTICLE{
	ELETRON,
	MUON,
	PION,
	KAON,
	PROTON
};
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
typedef std::vector<int> Vint;
typedef std::vector<double> VDouble;
typedef std::vector<HepLorentzVector> Vp4;

static long m_cout_all(0), m_cout_ngood(0), /*m_cout_pkpi(0),*/m_cout_skim(0);


Sigma0PionEta::Sigma0PionEta(const std::string& name, ISvcLocator* pSvcLocator):Algorithm(name,pSvcLocator){
		declareProperty("Vr0cut",   m_vr0cut=1.0);
		declareProperty("Vz0cut",   m_vz0cut=10.0);
		declareProperty("Probcut",   m_prob_cut=0.0);
		declareProperty("BestCandidate",   m_BestCandidate=true);
		declareProperty("UseOverallTof",   m_use_Total_TOF=false);
		declareProperty("SkimFlag", m_skim=false);
		declareProperty("beamE", m_beamE=2.300);
		declareProperty("ecms", m_ecms=4.600);
		declareProperty("debug", m_debug=true);
		declareProperty("use1c_sigma0", m_use1c_sigma0=true);
		declareProperty("CosThetaCut", m_CosThetaCut=0.93);
		declareProperty("CheckTotal", m_checktotal=true);

		//good shower
		declareProperty("PhotonMinEnergy", m_minEnergy = 0.025);
		declareProperty("GammaAngleCut", m_gammaAngleCut=20.0);
		declareProperty("GammatlCut",   m_gammatlCut=0.0);
		declareProperty("GammathCut",  m_gammathCut=14.0);
		declareProperty("PhotonMaxCosThetaBarrel", m_maxCosThetaBarrel = 0.8);
		declareProperty("PhotonMinCosThetaEndcap", m_minCosThetaEndcap = 0.84);
		declareProperty("PhotonMaxCosThetaEndcap", m_maxCosThetaEndcap = 0.92);
		declareProperty("PhotonMinEndcapEnergy",   m_minEndcapEnergy   = 0.050);

}

Sigma0PionEta::~Sigma0PionEta(){

	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter Deconstructor"<<"!!!!!!!!!!"<<std::endl;
		//add your code for deconstructor

}
void Sigma0PionEta::addItem()
{
	StatusCode status;
        status = m_tuple->addItem ("run", m_run);
		status = m_tuple->addItem ("event", m_event);
		status = m_tuple->addItem ("flag", m_flag);
        status = m_tuple->addItem ("index", m_index,0,150);

		status = m_tuple->addIndexedItem ("charge", m_index, m_charge);
        status = m_tuple->addItem ("nGood", m_ngood);
        status = m_tuple->addIndexedItem ("Ebeam", m_index, m_ebeam);
        status = m_tuple->addIndexedItem ("mBC", m_index, m_mbc);
        status = m_tuple->addIndexedItem ("mass", m_index, m_mass);
        status = m_tuple->addIndexedItem ("De", m_index, m_deltaE);
        status = m_tuple->addIndexedItem ("De_lmd_pip_eta", m_index, m_deltaE_lmd_pip_eta);

        //status = m_tuple->addItem ("lambda_p4",m_index, m_lambda_p4);
        status = m_tuple->addIndexedItem ("lambda_mass",m_index, m_lambda_mass);
        status = m_tuple->addIndexedItem ("lamda_chis1",m_index, m_lambda_chis1);
        status = m_tuple->addIndexedItem ("lambda_lchue",m_index, m_lambda_lchue);


        status = m_tuple->addIndexedItem ("eta_mass",m_index, m_eta_mass);
        status = m_tuple->addIndexedItem ("eta_gam1_p4",m_index, 4, m_eta_gam1_p4);
        status = m_tuple->addIndexedItem ("eta_gam2_p4",m_index, 4, m_eta_gam2_p4);
        status = m_tuple->addIndexedItem ("sigma0_gam_p4",m_index, 4, m_sigma0_gam_p4);

        status = m_tuple->addIndexedItem ("sigma0_mass",m_index, m_sigma0_mass);
        status = m_tuple->addIndexedItem ("sigma0_chis",m_index, m_sigma0_chis);

        status = m_tuple->addIndexedItem ("eta_chis",m_index, m_eta_chis);
        //status = m_tuple->addItem ("lmdc_chis",m_lmdc_chis);

		if(m_checktotal)
		{
			status = m_tuple->addItem ("mode1", m_mode1);
			status = m_tuple->addItem ("mode2", m_mode2);
			status = m_tuple->addItem ("mode3", m_mode3);
			status = m_tuple->addItem ("ndaughterAp", m_ndaughterAp, 0, 15);
			status = m_tuple->addIndexedItem ("Ap_id", m_ndaughterAp, m_Ap_id);
			status = m_tuple->addIndexedItem ("Ap_ptruth", m_ndaughterAp, 4, m_Ap_ptruth);
			status = m_tuple->addItem ("ndaughterAm",       m_ndaughterAm, 0, 15);
			status = m_tuple->addIndexedItem ("Am_id", m_ndaughterAm, m_Am_id);
			status = m_tuple->addIndexedItem ("Am_ptruth", m_ndaughterAm,  4, m_Am_ptruth);
		}
}

StatusCode Sigma0PionEta::initialize(){
		MsgStream log(msgSvc(), name());
		log<<MSG::INFO<<"Sigma0PionEta::initialize()"<<endreq;
		if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!Enter Initialize!!!!!!!!!"<<std::endl;
		//add your code here
		StatusCode status;
		NTuplePtr nt1(ntupleSvc(), "FILE1/Lmdc");
		if ( nt1 ) m_tuple = nt1;
		else {
				m_tuple = ntupleSvc()->book ("FILE1/Lmdc", CLID_ColumnWiseTuple, "exam N-Tuple example");
				if ( m_tuple )    
				{
				addItem();
				}
				else    {
						log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple) << endmsg;
						return StatusCode::FAILURE;
				}
		}

		if(m_checktotal)
		{
		NTuplePtr nt_mc(ntupleSvc(), "FILE1/mc_truth");
		if ( nt_mc ) m_tuple_mc = nt_mc;
		else {
				m_tuple_mc = ntupleSvc()->book ("FILE1/mc_truth", CLID_ColumnWiseTuple, "exam N-Tuple example");
				if ( m_tuple_mc )    
				{
					status = m_tuple_mc->addItem ("run", mc_run);
					status = m_tuple_mc->addItem ("event", mc_event);
					status = m_tuple_mc->addItem ("mode1", mc_mode1);
					status = m_tuple_mc->addItem ("mode2", mc_mode2);
					status = m_tuple_mc->addItem ("mode3", mc_mode3);
					status = m_tuple_mc->addItem ("ndaughterAp", mc_ndaughterAp, 0, 15);
					status = m_tuple_mc->addIndexedItem ("Ap_id", mc_ndaughterAp, mc_Ap_id);
					status = m_tuple_mc->addIndexedItem ("Ap_ptruth", mc_ndaughterAp, 4, mc_Ap_ptruth);
					status = m_tuple_mc->addItem ("ndaughterAm",       mc_ndaughterAm, 0, 15);
					status = m_tuple_mc->addIndexedItem ("Am_id", mc_ndaughterAm, mc_Am_id);
					status = m_tuple_mc->addIndexedItem ("Am_ptruth", mc_ndaughterAm,  4, mc_Am_ptruth);
				
				}
				else    {
						log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple_mc) << endmsg;
						return StatusCode::FAILURE;
				}
		}
		}

		return StatusCode::SUCCESS;

}

StatusCode Sigma0PionEta::beginRun(){
	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter beginRun"<<"!!!!!!!!!!"<<std::endl;
		MsgStream log(msgSvc(), name());
		log<<MSG::INFO<<"Sigma0PionEta::beginRun()"<<endreq;
		//add your code here
		return StatusCode::SUCCESS;

}
StatusCode Sigma0PionEta::execute(){

	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter execute"<<"!!!!!!!!!!"<<std::endl;
		MsgStream log(msgSvc(), name());
		log<<MSG::INFO<<"Sigma0PionEta::execute()"<<endreq;
		SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
		m_cout_all++;

		//initi
		const int LAMBDAC_SIGMA0_PION_PI0 = 11;

		m_flag= -1;

		int runNo = eventHeader->runNumber();
		int eventNo = eventHeader->eventNumber();
		
		
		mc_run = runNo;
		mc_event = eventNo;


	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Run number: "<<runNo<<"!!!!!!!!!!"<<std::endl;
		log << MSG::DEBUG <<"run, evtnum = "
				<< runNo << " , "
				<< eventNo <<endreq;

	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	log << MSG::DEBUG <<"ncharg, nneu, tottks = "
			<< evtRecEvent->totalCharged() << " , "
			<< evtRecEvent->totalNeutral() << " , "
			<< evtRecEvent->totalTracks() <<endreq;

	int mode1= ((eventHeader->flag1()/1000000))%1000; //QUESTION: What is the flag1?
	int mode2= (eventHeader->flag1()/1000)%1000 ;
	int mode3= eventHeader->flag1()%1000;

	Vint pdgid; pdgid.clear();
	Vint motherindex; motherindex.clear();

	int pdgid_p[100];
	int motheridx_p[100];

	int pdgid_m[100];
	int motheridx_m[100];

	int numParticle = 0;
	int numParticle_p = 0;
	int numParticle_m = 0;

	int ndaughterAp=0; double	Ap_ptruth[15][4]; int Ap_id[15];
	for ( int aa = 0; aa < 15; aa++ ) 
		for ( int ll = 0; ll < 4; ll++ ) 
			Ap_ptruth[aa][ll]=0;

	for ( int aa = 0; aa < 15; aa++ ) 
			Ap_id[aa]=0;

	int	ndaughterAm=0; double	Am_ptruth[15][4]; int Am_id[15];
	for ( int aa = 0; aa < 15; aa++ ) 
		for ( int ll = 0; ll < 4; ll++ ) 
			Am_ptruth[aa][ll]=0;

	for ( int aa = 0; aa < 15; aa++ ) 
			Am_id[aa]=0;

	if(m_debug) std::cerr<<"Enter McParticleCol"<<std::endl;

	if (eventHeader->runNumber()<0 && m_checktotal)
	{
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");

		if (!mcParticleCol)
		{
			std::cout << "Could not retrieve McParticelCol" << std::endl;
			return StatusCode::FAILURE;
		}
		else
		{
			Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
	
			for(; iter_mc != mcParticleCol->end(); ++iter_mc)
			{
				if ((*iter_mc)->primaryParticle()) continue; //QUESTION: What are primary particles?
				if (!(*iter_mc)->decayFromGenerator()) continue;
				int pdg = (*iter_mc)->particleProperty();
				int motherpdg = ((*iter_mc)->mother()).particleProperty();
				int mmotherpdg = (((*iter_mc)->mother()).mother()).particleProperty();

				if(pdg == 4122 && mode1 == 0 /*&& mode2 == 11*/) //lambda_c+
				{
					const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
					for(unsigned int ii = 0; ii < gc.size(); ii++) 
					{
						if( gc[ii]->particleProperty() == -22) continue; 
						Ap_id[ndaughterAp]=gc[ii]->particleProperty();
						
						for( int ll = 0; ll < 4; ll++ ) 
							Ap_ptruth[ndaughterAp][ll]=gc[ii]->initialFourMomentum()[ll];
						
						ndaughterAp++;
					}// End of "gc.size() > 0" IF

					m_flag = LAMBDAC_SIGMA0_PION_PI0;
				}

				
				if (pdg == -4122 && mode1 == 1 /*&& mode3 == 11*/) //lambda_c-
				{
					const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
					for(unsigned int ii = 0; ii < gc.size(); ii++) 
					{
						if( gc[ii]->particleProperty() == -22) continue;
						
						Am_id[ndaughterAm]=gc[ii]->particleProperty();
						for( int ll = 0; ll < 4; ll++ ) 
							Am_ptruth[ndaughterAm][ll]=gc[ii]->initialFourMomentum()[ll];
						
						ndaughterAm++;
					}// End of "gc.size() > 0" IF
					
					m_flag = LAMBDAC_SIGMA0_PION_PI0;
				}

				//if(mode1 ==0 && mode2 == 6 && pdg == 4122) 
				//{
					//const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
					//for(unsigned int ii = 0; ii < gc.size(); ii++) 
					//{
						//if( gc[ii]->particleProperty() == -22) continue;
					
						//Ap_id[ndaughterAp]=gc[ii]->particleProperty();
						//for( int ll = 0; ll < 4; ll++ ) 
							//Ap_ptruth[ndaughterAp][ll]=gc[ii]->initialFourMomentum()[ll];
					
						//ndaughterAp++;
					//}// End of "gc.size() > 0" IF

				//}

				//if (pdg == -4122 && mode1 == 0 && mode3 == 6) //lambda_c-
				//{
					//const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
					//for(unsigned int ii = 0; ii < gc.size(); ii++) 
					//{
						//if( gc[ii]->particleProperty() == -22) continue;
						
						//Am_id[ndaughterAm]=gc[ii]->particleProperty();
						//for( int ll = 0; ll < 4; ll++ ) 
							//Am_ptruth[ndaughterAm][ll]=gc[ii]->initialFourMomentum()[ll];
						
						//ndaughterAm++;
					//}// End of "gc.size() > 0" IF
					
				//}

				//mode2
				if(mode1 ==0 /*&& mode2 == 11*/ && pdg == 3212 && motherpdg == 4122) //sigma0 -> lmd gam
				{
					const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
					for(unsigned int ii = 0; ii < gc.size(); ii++) 
					{
						if( gc[ii]->particleProperty() == -22) continue;
					
						Ap_id[ndaughterAp]=gc[ii]->particleProperty();
						for( int ll = 0; ll < 4; ll++ ) 
							Ap_ptruth[ndaughterAp][ll]=gc[ii]->initialFourMomentum()[ll];
					
						ndaughterAp++;
					}// End of "gc.size() > 0" IF
			
				}

				if(mode1 == 0 /*&& mode2 == 11*/ && pdg == 221 && motherpdg == 4122) //eta -> gam gam
				{
					const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
					for(unsigned int ii = 0; ii < gc.size(); ii++) 
					{
						if( gc[ii]->particleProperty() == -22) continue;
					
						Ap_id[ndaughterAp]=gc[ii]->particleProperty();
						for( int ll = 0; ll < 4; ll++ ) 
							Ap_ptruth[ndaughterAp][ll]=gc[ii]->initialFourMomentum()[ll];
					
						ndaughterAp++;
					}// End of "gc.size() > 0" IF
			
				}

				if(mode1 == 0 /*&& mode2 == 11*/ && pdg == 3122 && motherpdg == 3212 && mmotherpdg == 4122) //lmd -> p pi
				{
					const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
					for(unsigned int ii = 0; ii < gc.size(); ii++) 
					{
						if( gc[ii]->particleProperty() == -22) continue;
					
						Ap_id[ndaughterAp]=gc[ii]->particleProperty();
						for( int ll = 0; ll < 4; ll++ ) 
							Ap_ptruth[ndaughterAp][ll]=gc[ii]->initialFourMomentum()[ll];
					
						ndaughterAp++;
					}// End of "gc.size() > 0" IF
			
				}

				//mode2
				if(mode1 == 1 /*&& mode3 == 11*/ && pdg == -3212 && motherpdg == -4122) //sigma0 -> lmd gam
				{
					const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
					for(unsigned int ii = 0; ii < gc.size(); ii++) 
					{
						if( gc[ii]->particleProperty() == -22) continue;
					
						Am_id[ndaughterAm]=gc[ii]->particleProperty();
						for( int ll = 0; ll < 4; ll++ ) 
							Am_ptruth[ndaughterAm][ll]=gc[ii]->initialFourMomentum()[ll];
					
						ndaughterAm++;
					}// End of "gc.size() > 0" IF
			
				}

				if(mode1 == 1 /*&& mode3 == 11*/ && pdg == 221 && motherpdg == -4122) //eta -> gam gam
				{
					const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
					for(unsigned int ii = 0; ii < gc.size(); ii++) 
					{
						if( gc[ii]->particleProperty() == -22) continue;
					
						Am_id[ndaughterAm]=gc[ii]->particleProperty();
						for( int ll = 0; ll < 4; ll++ ) 
							Am_ptruth[ndaughterAm][ll]=gc[ii]->initialFourMomentum()[ll];
					
						ndaughterAm++;
					}// End of "gc.size() > 0" IF
			
				}

				if(mode1 == 1 /*&& mode3 == 11*/  && pdg == -3122 && motherpdg == -3212 && mmotherpdg == -4122) //lmd -> p pi
				{
					const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
					for(unsigned int ii = 0; ii < gc.size(); ii++) 
					{
						if( gc[ii]->particleProperty() == -22) continue;
					
						Am_id[ndaughterAm]=gc[ii]->particleProperty();
						for( int ll = 0; ll < 4; ll++ ) 
							Am_ptruth[ndaughterAm][ll]=gc[ii]->initialFourMomentum()[ll];
					
						ndaughterAm++;
					}// End of "gc.size() > 0" IF
			
				}
			} 
		}
	}	 

	if(m_debug) std::cerr<<"!!!!! End McParticleCol !!!!"<<std::endl;
	if(m_checktotal)
	{
		mc_mode1 = mode1;
		mc_mode2 = mode2;
		mc_mode3 = mode3;

		mc_ndaughterAp = ndaughterAp;

		for ( int aa = 0; aa < ndaughterAp; aa++ ) 
		{
			mc_Ap_id[aa] = Ap_id[aa];
		}
		
		for ( int aa = 0; aa < ndaughterAp; aa++ ) 
			for ( int ll = 0; ll < 4; ll++ ) 
			{	
				mc_Ap_ptruth[aa][ll]=Ap_ptruth[aa][ll];
			}

		mc_ndaughterAm = ndaughterAm;
		
		for ( int aa = 0; aa < ndaughterAm; aa++ ) 
		{
			mc_Am_id[aa]=Am_id[aa];
		}

		
		for ( int aa = 0; aa < ndaughterAm; aa++ ) 
			for ( int ll = 0; ll < 4; ll++ ) 
			{		
				mc_Am_ptruth[aa][ll] = Am_ptruth[aa][ll];
			}

		m_tuple_mc->write();
	}
	
	if(m_debug) std::cerr<<"!!!! Begin analysis !!!!!"<<std::endl;

	//Begin analysis
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);

	Vint iGood;iGood.clear(); 
	Vint iProton; iProton.clear(); 
	Vint iPion; iPion.clear(); 
	Vint iKaon; iKaon.clear();
		//VDouble Rxy;Rxy.clear(); VDouble Rz; Rz.clear();
	Vint iProtonm; iProtonm.clear(); 
	Vint iProtonp; iProtonp.clear();
	Vint iProtonp_loose; iProtonp_loose.clear();
	Vint iProtonm_loose; iProtonm_loose.clear();
	Vint iPionm; iPionm.clear();  
	Vint iPionp; iPionp.clear();
	Vint iPionm_loose; iPionm_loose.clear();
	Vint iPionp_loose; iPionp_loose.clear();

		//Vint iKaonm; iKaonm.clear();  Vint iKaonp; iKaonp.clear();
	Vint iGam; iGam.clear();

	int nProton = 0, nPion = 0; //note: we dont require proton or pion to good!

		//VDouble Rxy_pm,Rz_pm; Rxy_pm.clear(); Rz_pm.clear();
		//VDouble Rxy_km,Rz_km; Rxy_km.clear(); Rz_km.clear();
		//VDouble Rxy_pim,Rz_pim; Rxy_pim.clear(); Rz_pim.clear();
		//VDouble Rxy_pp,Rz_pp; Rxy_pp.clear(); Rz_pp.clear();
		//VDouble Rxy_kp,Rz_kp; Rxy_kp.clear(); Rz_kp.clear();
		//VDouble Rxy_pip,Rz_pip; Rxy_pip.clear(); Rz_pip.clear();

	if(m_debug) std::cerr<<"ngood is "<<iGood.size()<<endl;

	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter GoodTrack"<<"!!!!!!!!!!"<<std::endl;
	if(!evtRecTrkCol)
	{
		std::cerr<<"!!!!!!!!!!!!!!!!"<<"Empty evtRecTrkCol"<<"!!!!!!"<<std::endl;
		return StatusCode::SUCCESS;
	}
	for(int i = 0; i < evtRecEvent->totalCharged(); i++){
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		RecMdcTrack *mdcTrack = (*itTrk)->mdcTrack();
		int charge = mdcTrack->charge();
		double vz=-10, vxy=-10;
				
		if(isGoodTrackForLambda(*itTrk))
		{
			if(charge > 0)
			{
				if(isProton(itTrk))
				{
					++nProton;
					iProtonp_loose.push_back(i);
				}
				else if(isPion(itTrk))
				{
					++nPion;
					iPionp_loose.push_back(i);
				}
			}
			else
			{
				if(isProton(itTrk))
				{
					++nProton;
					iProtonm_loose.push_back(i);
				}
				else if(isPion(itTrk))
				{
					++nPion;
					iPionm_loose.push_back(i);
				}
			}

		}


		if(!(isGoodTrk(itTrk, vz, vxy)))  
			continue;
		else
			iGood.push_back(i);

		if(m_debug) std::cerr<<"vz, vxy "<<vz << " , " << vxy <<endl;
	}

	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"End GoodTrack"<<"!!!!!!!!!!"<<std::endl;
	if(m_debug) std::cerr<<"ngood is "<<iGood.size()<<endl;

	if(iGood.size()<1) return StatusCode::SUCCESS;
	m_ngood = iGood.size();
	m_cout_ngood++;

	//int const Ng=iGood.size();

	//double rxy[Ng][3]={-10}, rz[Ng][3]={-10};
	//double rxy_m[Ng][3]={-10}, rz_m[Ng][3]={-10};

	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter PID"<<"!!!!!!!!!!"<<std::endl;
	for(int i=0; i<iGood.size(); i++)
	{
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + iGood[i];
		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
		int charge =  mdcTrk->charge();
		if(isPion(itTrk)){
			iPion.push_back(iGood[i]); 
			if(charge>0) {
				iPionp.push_back(iGood[i]);
			} else { 
				iPionm.push_back(iGood[i]);
			} 
		}

		//if(IsKaon(itTrk)){iKaon.push_back(iGood[i]); 
		//	if(charge>0) { 
		//		iKaonp.push_back(iGood[i]); 
		//		Rxy_kp.push_back(Rxy[i]);
		//		Rz_kp.push_back(Rz[i]);
		//	} 
		//	else { 
		//		iKaonm.push_back(iGood[i]);
		//		Rxy_km.push_back(Rxy[i]);
		//		Rz_km.push_back(Rz[i]); 
		//	} 
		//}

		if(isProton(itTrk))
		{
			iProton.push_back(iGood[i]); 
			if(charge>0) {
				iProtonp.push_back(iGood[i]); 
				//		Rxy_pp.push_back(Rxy[i]);Rz_pp.push_back(Rz[i]);
			}

			else { 
				iProtonm.push_back(iGood[i]); 
				//			Rxy_pm.push_back(Rxy[i]);
				//			Rz_pm.push_back(Rz[i]); 
			}
		}

	}

	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter Good Shower"<<"!!!!!!!!!!"<<std::endl;
		for(int i=evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++){
			EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
			double vz=-10, vxy=-10;
			if(!(isGoodShower(*itTrk)))  continue;
			
			iGam.push_back(i);
			//m_vz[cnt_good] = vz;
			//vz_mean+=vz;
			//m_vxy[cnt_good] = vxy;
			//cnt_good++;
		}
	
	if(iGam.size() < 3 || iGam.size() > 10) 
		return StatusCode::SUCCESS;
	
	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter  Vector Size"<<"!!!!!!!!!!"<<std::endl;
	if(iPion.size()<1) return StatusCode::SUCCESS;

	//loose pion and proton
	if(nPion < 2) return StatusCode::SUCCESS;
	if(nProton < 1) return StatusCode::SUCCESS;

	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"End Vector Size"<<"!!!!!!!!!!"<<std::endl;

		//************* Looking for eta **********
		Vint iEtagam1,iEtagam2;iEtagam1.clear(),iEtagam2.clear();
		Vdou massEta; massEta.clear();
		Vdou chisEta; chisEta.clear();
		Vp4 p4Eta; p4Eta.clear();
		Vp4 p4Eta1c; p4Eta1c.clear();
		Vp4 p4Etagam1, p4Etagam2; p4Etagam1.clear(), p4Etagam2.clear();
		vector<WTrackParameter> wtpEta1c;
		int nEta;
	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter Good Eta"<<"!!!!!!!!!!"<<std::endl;
		for(int i = 0; i < iGam.size()-1; i++) 
		{
			for(int j = i+1; j < iGam.size(); j++) 
			{
				EvtRecTrackIterator itTrki = evtRecTrkCol->begin() + iGam[i];
				RecEmcShower* shr1 = (*itTrki)->emcShower();

				EvtRecTrackIterator itTrkj = evtRecTrkCol->begin() + iGam[j];
				RecEmcShower* shr2 = (*itTrkj)->emcShower();

				HepLorentzVector p4_eta(0,0,0,0),p4_eta_1c(0,0,0,0);
				HepLorentzVector p4_gam1(0,0,0,0),p4_gam2(0,0,0,0);
				double eta_mass;
				double eta_chis;
				WTrackParameter wtp_eta_1c;
				if(isGoodEta(shr1,shr2,eta_mass,p4_eta,eta_chis,p4_eta_1c, p4_gam1, p4_gam2, wtp_eta_1c)){
					iEtagam1.push_back(iGam[i]);
					iEtagam2.push_back(iGam[j]);
					p4Etagam1.push_back(p4_gam1);
					p4Etagam2.push_back(p4_gam2);

					massEta.push_back(eta_mass);
					chisEta.push_back(eta_chis);
					p4Eta.push_back(p4_eta);
					p4Eta1c.push_back(p4_eta_1c);
					wtpEta1c.push_back(wtp_eta_1c);
				}

			}
		}

      nEta = iEtagam1.size();
      if(nEta==0) return StatusCode::SUCCESS;
		
	
	//************** Looking for lambda *************
	Vint iLmdpp,iLmdpim;iLmdpp.clear(),iLmdpim.clear();
	Vdou massLmd; massLmd.clear();
	Vdou chis1Lmd; chis1Lmd.clear();
	Vdou chis2Lmd; chis2Lmd.clear();
	Vdou lchueLmd; lchueLmd.clear();
	vector<WTrackParameter> wtpLmd1s; wtpLmd1s.clear();
	int nLmd;

	Vint iLmdpm,iLmdpip;iLmdpm.clear(),iLmdpip.clear();
	Vdou massALmd; massALmd.clear();
	Vdou chis1ALmd; chis1ALmd.clear();
	Vdou chis2ALmd; chis2ALmd.clear();
	Vdou lchueALmd; lchueALmd.clear();
	vector<WTrackParameter> wtpALmd1s; wtpALmd1s.clear();
	int nALmd;

	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter Good Lambda"<<"!!!!!!!!!!"<<std::endl;
	for(int i = 0; i < iProtonp_loose.size(); i++) 
	{
		for(int j = 0; j < iPionm_loose.size(); j++) 
		{
			EvtRecTrackIterator itTrki = evtRecTrkCol->begin() + iProtonp_loose[i];
			RecMdcKalTrack *ppKalTrk = (*itTrki)->mdcKalTrack();

			EvtRecTrackIterator itTrkj = evtRecTrkCol->begin() + iPionm_loose[j];
			RecMdcKalTrack *pimKalTrk = (*itTrkj)->mdcKalTrack();

			WTrackParameter wtp_lmd_1s;
			double lmd_1chis,lmd_2chis,lmd_lchue,lmd_mass;

			if(isGoodLambda(ppKalTrk, pimKalTrk, lmd_1chis, lmd_2chis, lmd_lchue, wtp_lmd_1s, lmd_mass)){
				iLmdpp.push_back(iProtonp_loose[i]);
				iLmdpim.push_back(iPionm_loose[j]);
				 massLmd.push_back(lmd_mass);
				chis1Lmd.push_back(lmd_1chis);
				chis2Lmd.push_back(lmd_2chis);
				lchueLmd.push_back(lmd_lchue);
				wtpLmd1s.push_back(wtp_lmd_1s);
			}

		}
	}

	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter Good anti Lambda"<<"!!!!!!!!!!"<<std::endl;
				// Loop each pbar pi+ pair, check if it is a Anti-lambda
	for(int i = 0; i < iProtonm_loose.size(); i++) 
	{
		for(int j = 0; j < iPionp_loose.size(); j++) 
		{
			EvtRecTrackIterator itTrki = evtRecTrkCol->begin() + iProtonm_loose[i];
			RecMdcKalTrack *pmKalTrk = (*itTrki)->mdcKalTrack();

			EvtRecTrackIterator itTrkj = evtRecTrkCol->begin() + iPionp_loose[j];
			RecMdcKalTrack *pipKalTrk = (*itTrkj)->mdcKalTrack();

			WTrackParameter wtp_lmd_1s;
			double lmd_1chis,lmd_2chis,lmd_lchue,lmd_mass;

			if(isGoodLambda(pmKalTrk, pipKalTrk, lmd_1chis, lmd_2chis, lmd_lchue, wtp_lmd_1s, lmd_mass)){
				iLmdpm.push_back(iProtonm_loose[i]);
				iLmdpip.push_back(iPionp_loose[j]);
				massALmd.push_back(lmd_mass);
				chis1ALmd.push_back(lmd_1chis);
				chis2ALmd.push_back(lmd_2chis);
				lchueALmd.push_back(lmd_lchue);
				wtpALmd1s.push_back(wtp_lmd_1s);
			}

		}
	}

	nLmd = iLmdpp.size();
	nALmd = iLmdpm.size();
	if(nLmd==0&&nALmd==0) return StatusCode::SUCCESS;
		
	//********* Looking for sigma0 ****************
	Vdou massSigma0; massSigma0.clear();
	Vdou chisSigma0; chisSigma0.clear();
	Vint iSigma0gam; iSigma0gam.clear();
	Vint iSigma0lmd; iSigma0lmd.clear();
	Vp4  p4Sigma01c; p4Sigma01c.clear();
	Vp4  p4Sigma0gam; p4Sigma0gam.clear();
	vector<WTrackParameter> wtpSigma01c; 
	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter Good Sigma0"<<"!!!!!!!!!!"<<std::endl;
	for(int iLmd = 0; iLmd < iLmdpp.size(); ++iLmd)
		for( int iShr = 0; iShr < iGam.size(); ++iShr)
		{
			EvtRecTrackIterator itTrki = evtRecTrkCol->begin() + iGam[iShr];
			RecEmcShower* gamShr = (*itTrki)->emcShower();
			
			double sigma0_mass, sigma0_chis;
			HepLorentzVector p4_sigma0, p4_sigma01c, p4_gam;
			WTrackParameter wtp_sigma0_1c;
			if(isGoodSigma0(&wtpLmd1s[iLmd], gamShr, sigma0_mass, p4_sigma0, sigma0_chis, p4_sigma01c, p4_gam, wtp_sigma0_1c))
			{
				massSigma0.push_back(sigma0_mass);
				chisSigma0.push_back(sigma0_chis);
				iSigma0gam.push_back(iGam[iShr]);
				p4Sigma0gam.push_back(p4_gam);
				p4Sigma01c.push_back(p4_sigma01c);
				iSigma0lmd.push_back(iLmd);
				wtpSigma01c.push_back(wtp_sigma0_1c);
			}
		}	
	
	// Anti-sigma0
	Vdou massASigma0; massASigma0.clear();
	Vdou chisASigma0; chisASigma0.clear();
	Vint iASigma0gam; iASigma0gam.clear();
	Vint iASigma0lmd; iASigma0lmd.clear();
	Vp4  p4ASigma01c; p4ASigma01c.clear();
	Vp4  p4ASigma0gam; p4ASigma0gam.clear();
	vector<WTrackParameter> wtpASigma01c; 

	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter Good anti Sigma0"<<"!!!!!!!!!!"<<std::endl;
	for(int iALmd = 0; iALmd < iLmdpm.size(); ++iALmd)
		for( int iShr = 0; iShr < iGam.size(); ++iShr)
		{
			EvtRecTrackIterator itTrki = evtRecTrkCol->begin() + iGam[iShr];
			RecEmcShower* gamShr = (*itTrki)->emcShower();
			
			double Asigma0_mass, Asigma0_chis;
			HepLorentzVector p4_Asigma0, p4_Asigma01c, p4_gam;
			WTrackParameter wtp_sigma0_1c;
			if(isGoodSigma0(&wtpALmd1s[iALmd], gamShr, Asigma0_mass, p4_Asigma0, Asigma0_chis, p4_Asigma01c, p4_gam, wtp_sigma0_1c))
			{
				massASigma0.push_back(Asigma0_mass);
				chisASigma0.push_back(Asigma0_chis);
				iASigma0gam.push_back(iGam[iShr]);
				p4ASigma0gam.push_back(p4_gam);
				p4ASigma01c.push_back(p4_Asigma01c);
				iASigma0lmd.push_back(iALmd);
				wtpASigma01c.push_back(wtp_sigma0_1c);
			}
		}	

	//************* Analysis of lambda_c -> sigma0 pion eta ************* 

	HepLorentzVector proton_p4, kaon_p4, pion_p4;

	int cnt_Lp=0, cnt_Lm=0; 
	int index = 0;
	double tmp_chisq_lp = INFINITY, tmp_chisq_lm = INFINITY;
	
	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter Lambda_c +"<<"!!!!!!!!!!"<<std::endl;
	//lambda_c + --> sigma0 pi+ eta
	bool have_best_lambda_cp = false;
	m_run=runNo;
	m_event=eventNo;
	if(m_checktotal)
	{
		m_mode1 = mode1;
		m_mode2 = mode2;
		m_mode3 = mode3;

		m_ndaughterAp = ndaughterAp;

		for ( int aa = 0; aa < ndaughterAp; aa++ ) 
		{
			m_Ap_id[aa] = Ap_id[aa];
		}
		
		for ( int aa = 0; aa < ndaughterAp; aa++ ) 
			for ( int ll = 0; ll < 4; ll++ ) 
			{	
				m_Ap_ptruth[aa][ll]=Ap_ptruth[aa][ll];
			}

		m_ndaughterAm = ndaughterAm;
		
		for ( int aa = 0; aa < ndaughterAm; aa++ ) 
		{
			m_Am_id[aa]=Am_id[aa];
		}

		
		for ( int aa = 0; aa < ndaughterAm; aa++ ) 
			for ( int ll = 0; ll < 4; ll++ ) 
			{		
				m_Am_ptruth[aa][ll] = Am_ptruth[aa][ll];
			}

	}

	if(m_debug) cerr<<"nPionp= "<<iPionp.size()<<endl;
	for(int i_pionp=0;i_pionp<iPionp.size();i_pionp++)
	{
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + iPionp[i_pionp];
		RecMdcKalTrack::setPidType  (RecMdcKalTrack::pion);
		RecMdcKalTrack *pipTrk = (*itTrk)->mdcKalTrack();
		HepLorentzVector p4_pionp = pipTrk->p4(PDG::Pion);
		WTrackParameter wtp_pip(PDG::Pion, pipTrk->getZHelix(),pipTrk->getZError());

		for(int i_sigma0 = 0; i_sigma0 < wtpSigma01c.size(); i_sigma0++)
		{

			WTrackParameter wtp_sigma0 = wtpSigma01c[i_sigma0];
			HepLorentzVector p4_sigma0_1c = p4Sigma01c[i_sigma0];

			for(int i_eta = 0; i_eta < wtpEta1c.size(); i_eta++)
			{
				if(m_debug) cerr<<"Enter Loop!!!!!!!!!"<<endl;
				if(iEtagam1[i_eta] == iSigma0gam[i_sigma0] || 
				   iEtagam2[i_eta] == iSigma0gam[i_sigma0]  )
					continue;

				int iLmd = iSigma0lmd[i_sigma0];
				if( iLmdpp[iLmd] == iPionp[i_pionp] )
					continue;

				WTrackParameter wtp_eta = wtpEta1c[i_eta];	
				HepLorentzVector p4_eta_1c = p4Eta1c[i_eta];

				//boost
				double lmdc_chisq = -999;
				double mbc = -10;
				HepLorentzVector tot_p4, tot_lmd;
				//if(!Mass1c_Lmdc(wtp_sigma0, wtp_pip, wtp_eta,lmdc_chisq, tot_p4))
				//	continue;
			
				tot_p4 =  p4_sigma0_1c + p4_pionp + p4_eta_1c;  
				tot_lmd = wtpLmd1s[iLmd].p() + p4_pionp + p4_eta_1c;

				HepLorentzVector cms(0,0,0,m_ecms);
				HepLorentzVector tot_p4_boost, tot_lmd_boost;
				tot_p4_boost = tot_p4.boost(-0.011,0.,0.);
				tot_lmd_boost = tot_lmd.boost(-0.011,0.,0.);
				double mbc2 = m_beamE*m_beamE- tot_p4_boost.v().mag2();
				if(mbc2 > 0)
					mbc = sqrt(mbc2);
				else 
					continue;

				double deltaE = tot_p4_boost.t() - m_beamE;
				double deltaE_lmd_pip_eta = tot_lmd_boost.t() - m_beamE;
				
				if( !(2.25 < mbc && mbc < 2.3 && -0.1 < deltaE && deltaE < 0.1) )
					continue;

				if(m_debug) cerr<<"Enter Selection!!!!!!!!!!!"<<endl;
				if(m_BestCandidate)
				{
					if(fabs(lmdc_chisq)<tmp_chisq_lp)
					{
						have_best_lambda_cp = true;
						cnt_Lp++;
						tmp_chisq_lp=fabs(lmdc_chisq);

					//	m_lmdc_chis[index] = lmdc_chisq;
						
						m_mass[index] = tot_p4.m();
				//		HepLorentzVector p4_lambda = wtpLmd1s[iSigma0lmd[i_sigma0]].p(); 
						for(int pp=0;pp<4;pp++)
						{
							m_eta_gam1_p4[index][pp] = p4Etagam1[i_eta][pp];
							m_eta_gam2_p4[index][pp] = p4Etagam2[i_eta][pp];
							m_sigma0_gam_p4[index][pp] =  p4Sigma0gam[i_sigma0][pp];
						}
						m_charge[index]=1;
						m_ebeam[index] = m_beamE;
						m_mbc[index] =mbc;
						m_deltaE[index] =deltaE;
						m_deltaE_lmd_pip_eta[index] =deltaE_lmd_pip_eta;
						m_index = 1;
						
						m_eta_mass[index] = massEta[i_eta];
						
						m_lambda_mass[index] = massLmd[iSigma0lmd[i_sigma0]];
						m_lambda_chis1[index] = chis1Lmd[iSigma0lmd[i_sigma0]];
						m_lambda_lchue[index] = lchueLmd[iSigma0lmd[i_sigma0]];

						m_sigma0_mass[index] = massSigma0[i_sigma0];
						m_sigma0_chis[index] = chisSigma0[i_sigma0];

						m_eta_chis[index] = chisEta[i_eta];

					}
				}
				else 
				{
					++cnt_Lp;
					//m_lmdc_chis[index] = lmdc_chisq;

					if(m_debug) cerr<<"Before!!!!!!!!!!!"<<endl;
					m_mass[index] = tot_p4.m();
					if(m_debug) cerr<<"After!!!!!!!!!!!"<<endl;
					HepLorentzVector p4_lambda = wtpLmd1s[iSigma0lmd[i_sigma0]].p(); 

					for(int pp=0;pp<4;pp++)
					{
						m_eta_gam1_p4[index][pp] = p4Etagam1[i_eta][pp];
						m_eta_gam2_p4[index][pp] = p4Etagam2[i_eta][pp];
						m_sigma0_gam_p4[index][pp] =  p4Sigma0gam[i_sigma0][pp];
					}
					m_charge[index]=1;
					m_ebeam[index] = m_beamE;
					m_mbc[index] =mbc;
					m_deltaE[index] =deltaE;
					m_deltaE_lmd_pip_eta[index] =deltaE_lmd_pip_eta;

					m_lambda_mass[index] = massLmd[iSigma0lmd[i_sigma0]];
					m_lambda_chis1[index] = chis1Lmd[iSigma0lmd[i_sigma0]];
					m_lambda_lchue[index] = lchueLmd[iSigma0lmd[i_sigma0]];

					m_eta_mass[index] = massEta[i_eta];
					m_sigma0_mass[index] = massSigma0[i_sigma0];
					m_sigma0_chis[index] = chisSigma0[i_sigma0];

					m_eta_chis[index] = chisEta[i_eta];

					++index;
					m_index = index;

				}
			}

		}

	}

	if(m_BestCandidate)
	{	
		if(have_best_lambda_cp)
			m_tuple->write();
	}
	else if(index > 0)
		m_tuple->write();

	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<cnt_Lp<<"!!!!!!!!!!"<<std::endl;

	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter Good Lambda_c -"<<"!!!!!!!!!!"<<std::endl;
	//lamdac- -- > Asigma0 pi- eta
	bool have_best_lambda_cm = false;
	index = 0;
	m_run=runNo;
	m_event=eventNo;
	if(m_checktotal)
	{
		m_mode1 = mode1;
		m_mode2 = mode2;
		m_mode3 = mode3;

		m_ndaughterAp = ndaughterAp;

		for ( int aa = 0; aa < ndaughterAp; aa++ ) 
		{
			m_Ap_id[aa] = Ap_id[aa];
		}
		
		for ( int aa = 0; aa < ndaughterAp; aa++ ) 
			for ( int ll = 0; ll < 4; ll++ ) 
			{	
				m_Ap_ptruth[aa][ll]=Ap_ptruth[aa][ll];
			}

		m_ndaughterAm = ndaughterAm;
		
		for ( int aa = 0; aa < ndaughterAm; aa++ ) 
		{
			m_Am_id[aa]=Am_id[aa];
		}

		
		for ( int aa = 0; aa < ndaughterAm; aa++ ) 
			for ( int ll = 0; ll < 4; ll++ ) 
			{		
				m_Am_ptruth[aa][ll] = Am_ptruth[aa][ll];
			}

	}

	
	if(m_debug) cerr<<"nPionm: "<<iPionm.size();
	for(int i_pionm=0;i_pionm<iPionm.size();i_pionm++)
	{
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + iPionm[i_pionm];
		RecMdcKalTrack::setPidType  (RecMdcKalTrack::pion);
		RecMdcKalTrack *pimTrk = (*itTrk)->mdcKalTrack();
		HepLorentzVector p4_pionm = pimTrk->p4(PDG::Pion);
		WTrackParameter wtp_pim(PDG::Pion, pimTrk->getZHelix(),pimTrk->getZError());

		for(int i_sigma0 = 0; i_sigma0 < wtpASigma01c.size(); i_sigma0++)
		{
			WTrackParameter wtp_sigma0 = wtpASigma01c[i_sigma0];
			HepLorentzVector p4_sigma0_1c = p4ASigma01c[i_sigma0];

			for(int i_eta = 0; i_eta < wtpEta1c.size(); i_eta++)
			{
				if(m_debug) cerr<<"Enter Loop!!!"<<endl;
				if(iEtagam1[i_eta] == iASigma0gam[i_sigma0] || 
						iEtagam2[i_eta] == iASigma0gam[i_sigma0]  )
					continue;

				int iLmd = iASigma0lmd[i_sigma0];
				if( iLmdpm[iLmd] == iPionm[i_pionm] )
					continue;

				WTrackParameter wtp_eta = wtpEta1c[i_eta];	
				HepLorentzVector p4_eta_1c = p4Eta1c[i_eta];

				//boost
				double lmdc_chisq = -999;
				double mbc = -10;
				HepLorentzVector tot_p4, tot_lmd;
//				if(!Mass1c_Lmdc(wtp_sigma0, wtp_pim, wtp_eta,lmdc_chisq, tot_p4))
//					continue;

				tot_p4 =  p4_sigma0_1c + p4_pionm + p4_eta_1c;  
				tot_lmd = wtpALmd1s[iLmd].p() + p4_pionm + p4_eta_1c;

				HepLorentzVector cms(0,0,0,m_ecms);
				HepLorentzVector tot_p4_boost,tot_lmd_boost;
				tot_p4_boost = tot_p4.boost(-0.011,0.,0.);
				tot_lmd_boost = tot_lmd.boost(-0.011,0.,0.);
				double mbc2 = m_beamE*m_beamE- tot_p4_boost.v().mag2();
				if(mbc2 > 0)
					mbc = sqrt(mbc2);
				else 
					continue;
				double deltaE = tot_p4_boost.t() - m_beamE;
				double deltaE_lmd_pip_eta = tot_lmd_boost.t() - m_beamE;

				if( !(2.25 < mbc && mbc < 2.3 && -0.1 < deltaE && deltaE < 0.1) )
					continue;

				if(m_BestCandidate)
				{
					if(fabs(lmdc_chisq)<tmp_chisq_lm)
					{
						have_best_lambda_cm = true;
						cnt_Lm++;
						//m_lmdc_chis[index] = lmdc_chisq;
						m_mass[index] = tot_p4.m();
		//				HepLorentzVector p4_lambda = wtpALmd1s[iASigma0lmd[i_sigma0]].p(); 
						for(int pp=0;pp<4;pp++)
						{
							m_eta_gam1_p4[index][pp] = p4Etagam1[i_eta][pp];
							m_eta_gam2_p4[index][pp] = p4Etagam2[i_eta][pp];
							m_sigma0_gam_p4[index][pp] =  p4ASigma0gam[i_sigma0][pp];
						}
						m_charge[index] = -1;
						m_ebeam[index] = m_beamE;
						m_mbc[index] = mbc;
						m_deltaE[index] = deltaE;
						m_deltaE_lmd_pip_eta[index] = deltaE_lmd_pip_eta;

						m_lambda_mass[index] = massALmd[iASigma0lmd[i_sigma0]];
						m_lambda_chis1[index] = chis1ALmd[iASigma0lmd[i_sigma0]];
						m_lambda_lchue[index] = lchueALmd[iASigma0lmd[i_sigma0]];

						m_eta_mass[index] = massEta[i_eta];
						m_sigma0_mass[index] = massASigma0[i_sigma0];
						m_sigma0_chis[index] = chisASigma0[i_sigma0];

						m_eta_chis[index] = chisEta[i_eta];

						m_index = 1;
					}
				}
				else 
				{
					++cnt_Lm;
					m_mass[index] = tot_p4.m();
				//	m_lmdc_chis[index] = lmdc_chisq;
					HepLorentzVector p4_lambda = wtpALmd1s[iASigma0lmd[i_sigma0]].p(); 

					for(int pp=0;pp<4;pp++)
					{
						m_eta_gam1_p4[index][pp] = p4Etagam1[i_eta][pp];
						m_eta_gam2_p4[index][pp] = p4Etagam2[i_eta][pp];
						m_sigma0_gam_p4[index][pp] =  p4ASigma0gam[i_sigma0][pp];
					}
					m_charge[index] = -1;
					m_ebeam[index] = m_beamE;
					m_mbc[index] = mbc;
					m_deltaE[index] =deltaE;
					m_deltaE_lmd_pip_eta[index] = deltaE_lmd_pip_eta;

					m_lambda_mass[index] = massALmd[iASigma0lmd[i_sigma0]];
					m_lambda_chis1[index] = chis1ALmd[iASigma0lmd[i_sigma0]];
					m_lambda_lchue[index] = lchueALmd[iASigma0lmd[i_sigma0]];

					m_eta_mass[index] = massEta[i_eta];
					m_sigma0_mass[index] = massASigma0[i_sigma0];
					m_sigma0_chis[index] = chisASigma0[i_sigma0];

					m_eta_chis[index] = chisEta[i_eta];

					++index;
					m_index = index;

				}
			}

		}

	}

	if(m_BestCandidate)
	{	
		if(have_best_lambda_cm)
			m_tuple->write();
	}
	else if(index > 0)
		m_tuple->write();


	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<cnt_Lm<<"!!!!!!!!!!"<<std::endl;
	return StatusCode::SUCCESS;

}
StatusCode Sigma0PionEta::endRun(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"Sigma0PionEta::endRun()"<<endreq;
	//add your code here
	if(m_debug) std::cerr<<"!!!!!!!!!!!!!!!!!"<<"Enter endRun"<<"!!!!!!!!!!"<<std::endl;
	return StatusCode::SUCCESS;

}
StatusCode Sigma0PionEta::finalize(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"Sigma0PionEta::finalize()"<<endreq;
	cout<<"total events "<< m_cout_all <<endl;
	cout<<"pass ngood  "<< m_cout_ngood <<endl;
	//cout<<"pass p k pi  "<< m_cout_pkpi <<endl;

	//add your code here
	return StatusCode::SUCCESS;
}

bool Sigma0PionEta::isGoodTrk(EvtRecTrackIterator itTrk, double &vz , double &vxy) {
	if ( !(*itTrk)->isMdcTrackValid() ) return false;
	if ( !(*itTrk)->isMdcKalTrackValid() ) return false;
	RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();

	Hep3Vector xorigin(0,0,0);
	IVertexDbSvc*  vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if(vtxsvc->isVertexValid()){
		double* dbv = vtxsvc->PrimaryVertex();
		double* vv = vtxsvc->SigmaPrimaryVertex();
		//pretect
		if(vtxsvc->PrimaryVertex()[0]>100||vtxsvc->PrimaryVertex()[1]>100 || vtxsvc->PrimaryVertex()[2]>100)
		{
			cout<<"Vertex is abnormal! check your jobOption "<<endl;
			dbv[0]=0;
			dbv[1]=0;
			dbv[2]=0;
		}       
		if(m_debug)cout<<dbv[0]<< " , "<< dbv[1] <<" , "<< dbv[2] <<endl;
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
	}
	else { cout<<"warning !!! IVertexDbSvc is inValid!"<<endl;}
	if(m_debug) cout<<" xorigin : "<< xorigin[0] << " , "<< xorigin[1]<< " , "<< xorigin[2]<<endl;

	HepVector a = mdcTrk->helix();
	HepSymMatrix Ea = mdcTrk->err();
	HepPoint3D point0(0.,0.,0.);
	HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
	VFHelix helixip3(point0,a,Ea);
	helixip3.pivot(IP);
	HepVector  vecipa = helixip3.a();
	double dr=(vecipa[0]);
		double dz=(vecipa[3]);
		double costheta=cos(mdcTrk->theta());
		if (  fabs(dr)>= m_vr0cut) return false;
		if (  fabs(dz)>= m_vz0cut ) return false;
		if ( fabs(costheta) >= m_CosThetaCut ) return false;
		vz  = dz ;
		vxy = dr;

		return true;
}


//add your code here,for other member-functions
bool Sigma0PionEta::isProton(EvtRecTrackIterator itTrk)
{

		//double m_prob_cut =0.001;
		if(!(*itTrk)->isMdcKalTrackValid()) return false;
		ParticleID *pid = ParticleID::instance();
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		if(m_use_Total_TOF)
		{pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2()| pid->useTof());}
		else { pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2());}
		pid->identify(pid->onlyProton()|pid->onlyPion()|pid->onlyKaon());
		pid->calculate();
		if(!(pid->IsPidInfoValid())) return false;

		if(pid->prob(4)<m_prob_cut) return  false;

		if( pid->prob(4)< pid->prob(3) || pid->prob(4)< pid->prob(2)) return false;
		return true;

}

bool Sigma0PionEta::isPion( EvtRecTrackIterator itTrk )
{
		//double m_prob_cut =0.001;
		if(!(*itTrk)->isMdcKalTrackValid()) return false;
		ParticleID *pid = ParticleID::instance();
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		if(m_use_Total_TOF)
		{pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2()| pid->useTof());}
		else { pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2());}
		pid->identify(pid->onlyPion()|pid->onlyKaon());
		pid->calculate();
		if(!(pid->IsPidInfoValid())) return false;
		//prob_k = pid->prob(3);
		//prob_pi=pid->prob(2);

		if(pid->prob(2)<m_prob_cut) return  false;

		//if( pid->prob(2)< pid->prob(3) || pid->prob(2)< pid->prob(4)) return false;
		if( pid->prob(2)< pid->prob(3)) return false;
		return true;

}
bool Sigma0PionEta::isKaon( EvtRecTrackIterator itTrk )
{
		//double m_prob_cut =0.001;
		if(!(*itTrk)->isMdcKalTrackValid()) return false;
		ParticleID *pid = ParticleID::instance();
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		if(m_use_Total_TOF)
		{pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2()| pid->useTof());}
		else { pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2());}
		pid->identify(pid->onlyPion()|pid->onlyKaon());
		pid->calculate();
		if(!(pid->IsPidInfoValid())) return false;

		if(pid->prob(3)<m_prob_cut) return  false;

		//if( pid->prob(3)< pid->prob(2) || pid->prob(3)< pid->prob(4)) return false;
		if( pid->prob(3)< pid->prob(2)) return false;
		return true;

}


bool Sigma0PionEta::isGoodShower(EvtRecTrack* trk){

	if(!trk->isEmcShowerValid()){
		return false;
	}


	Hep3Vector xorigin(0,0,0);
	IVertexDbSvc*  vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if(vtxsvc->isVertexValid()){
		double* dbv = vtxsvc->PrimaryVertex();
		double*  vv = vtxsvc->SigmaPrimaryVertex();
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
	}

	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);

	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);


	RecEmcShower *emcTrk = trk->emcShower();

	Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
	HepLorentzVector shP4 = getP4(emcTrk,xorigin);
	double cosThetaSh = shP4.vect().cosTheta();
	double eraw = emcTrk->energy();
	double getTime = emcTrk->time();
	if(getTime>m_gammathCut||getTime<m_gammatlCut) return false;
	if(!((fabs(cosThetaSh) < m_maxCosThetaBarrel&&eraw > m_minEnergy)||((fabs(cosThetaSh) > m_minCosThetaEndcap)&&(fabs(cosThetaSh) < m_maxCosThetaEndcap)&&(eraw > m_minEndcapEnergy))))  return false;

	double dang = 200.; 
	for(int j = 0; j < evtRecEvent->totalCharged(); j++)
	{
		EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j; 
		if(!(*jtTrk)->isExtTrackValid()) continue;
		RecExtTrack *extTrk = (*jtTrk)->extTrack();
		if(extTrk->emcVolumeNumber() == -1) continue;
		Hep3Vector extpos = extTrk->emcPosition();
		double angd = extpos.angle(emcpos);
		if(angd < dang) dang = angd;
	}
		//if(dang>=200) return false;
	dang = dang * 180 / (CLHEP::pi);
	if(fabs(dang) < m_gammaAngleCut) return false;

	return true;
}

bool Sigma0PionEta::isGoodLambda(RecMdcKalTrack* ppTrk, RecMdcKalTrack* pimTrk, double& lmd_1chis, double& lmd_2chis, double& lmd_lchue, WTrackParameter &swvlmd, double& lmd_mass){
		lmd_1chis=-100;
		lmd_2chis=-100;
		lmd_lchue =-100;
		lmd_mass =-100;

		HepPoint3D vx(0., 0., 0.);
		HepSymMatrix Evx(3, 0);
		double bx = 1E+6;
		double by = 1E+6;
		double bz = 1E+6;
		Evx[0][0] = bx*bx;
		Evx[1][1] = by*by;
		Evx[2][2] = bz*bz;
		VertexParameter vxpar;
		vxpar.setVx(vx);
		vxpar.setEvx(Evx);

		VertexFit *vtxfit_s = VertexFit::instance();
		SecondVertexFit *vtxfit2 = SecondVertexFit::instance();

		WTrackParameter  wvpplmdTrk, wvpimlmdTrk;

		RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);
		wvpplmdTrk = WTrackParameter(xmass[4], ppTrk->getZHelixP(),ppTrk->getZErrorP());

		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
		wvpimlmdTrk = WTrackParameter(xmass[2], pimTrk->getZHelix(),pimTrk->getZError());

		//******primary vertex fit
		vtxfit_s->init();
		vtxfit_s->setChisqCut(200);
		vtxfit_s->AddTrack(0, wvpplmdTrk);
		vtxfit_s->AddTrack(1, wvpimlmdTrk);
		vtxfit_s->AddVertex(0, vxpar, 0, 1);
		bool okvs=vtxfit_s->Fit(0);

		if(!okvs) return false;

		vtxfit_s->Swim(0);
		vtxfit_s->BuildVirtualParticle(0);
		lmd_1chis = vtxfit_s->chisq(0);
		//Cut
		//if(lmd_1chis >= 100)
		//	return false;
	
		WTrackParameter  wvlmd = vtxfit_s->wVirtualTrack(0);
		HepLorentzVector p4_lmd_1s=wvlmd.p();

		lmd_mass = p4_lmd_1s.m();

		if(!Cut::Lambda(lmd_mass)) return false;

		WTrackParameter wppFlmd = vtxfit_s->wtrk(0);
		WTrackParameter wpimFlmd = vtxfit_s->wtrk(1);
		//******second vertex fit
		HepPoint3D newvx(0., 0., 0.);
		HepSymMatrix newEvx(3, 0);
		VertexParameter primaryVpar;
		IVertexDbSvc*  vtxsvc;
		Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
		if(!(vtxsvc->isVertexValid())) {
				cout<<"Attention --!(vtxsvc->isVertexValid())"<<endl;
		}
		double* db_vx = vtxsvc->PrimaryVertex();
		double* db_vx_err = vtxsvc->SigmaPrimaryVertex();
		newvx.setX(db_vx[0]);
		newvx.setY(db_vx[1]);
		newvx.setZ(db_vx[2]);
		newEvx[0][0] = db_vx_err[0]*db_vx_err[0];
		newEvx[1][1] = db_vx_err[1]*db_vx_err[1];
		newEvx[2][2] = db_vx_err[2]*db_vx_err[2];
		primaryVpar.setVx(newvx);
		primaryVpar.setEvx(newEvx);
		vtxfit2->init();
		vtxfit2->setChisqCut(200);
		vtxfit2->setPrimaryVertex(primaryVpar);
		vtxfit2->AddTrack(0, wvlmd);
		vtxfit2->setVpar(vtxfit_s->vpar(0));
		bool okv2=vtxfit2->Fit();
		if(!okv2) return false;

		lmd_2chis= vtxfit2->chisq();
		swvlmd = vtxfit2->wpar();
		//p4_lmd_2s=wlmd.p();

		double lmd_dl  = vtxfit2->decayLength();
		double lmd_dle = vtxfit2->decayLengthError();
		lmd_lchue = lmd_dl/lmd_dle;
		
		if(lmd_lchue > 2)
			return true;
		else
			return false;
}


bool Sigma0PionEta::isGoodEta(RecEmcShower *shr1,RecEmcShower *shr2,double& eta_mass,HepLorentzVector& p4_eta,double& eta_chis,HepLorentzVector& p4_eta_1c,
		HepLorentzVector &p4_gam1, HepLorentzVector &p4_gam2, WTrackParameter &etaTrk){

		Hep3Vector xorigin(0,0,0);
		IVertexDbSvc*  vtxsvc;
		Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
		if(vtxsvc->isVertexValid()){
				double* dbv = vtxsvc->PrimaryVertex();
				double*  vv = vtxsvc->SigmaPrimaryVertex();
				xorigin.setX(dbv[0]);
				xorigin.setY(dbv[1]);
				xorigin.setZ(dbv[2]);
		}

		eta_mass=-100;
		eta_chis=-100;

		HepLorentzVector g1P4 = getP4(shr1,xorigin);
		HepLorentzVector g2P4 = getP4(shr2,xorigin);
		p4_eta = g1P4 + g2P4;
		eta_mass = p4_eta.m();

		if(!Cut::Eta(eta_mass)) return false;

		//double xmeta=0.134976;

		KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
		kmfit->init();
		kmfit->setChisqCut(200);
		kmfit->setIterNumber(5);
		kmfit->AddTrack(0, 0.0, shr1);
		kmfit->AddTrack(1, 0.0, shr2);
		kmfit->AddResonance(0, PDG::Eta, 0, 1);

		bool olmdq =kmfit->Fit(0);
		if(!olmdq) return false;

		kmfit->BuildVirtualParticle(0);
		eta_chis = kmfit->chisq(0);
//		if(eta_chis>100) return false;

		p4_gam1 = kmfit->pfit(0);
		p4_gam2 = kmfit->pfit(1);

		p4_eta_1c=kmfit->pfit(0)+kmfit->pfit(1);
		etaTrk = kmfit->wVirtualTrack(0);

		return true;
}

bool Sigma0PionEta::isGoodSigma0(WTrackParameter *lmd,RecEmcShower *gammaShr,double& sigma0_mass,HepLorentzVector& p4_sigma0,double& sigma0_chis,HepLorentzVector& p4_sigma0_1c, HepLorentzVector &p4_gam, WTrackParameter &wtp_sigma0_1c){

		Hep3Vector xorigin(0,0,0);
		IVertexDbSvc*  vtxsvc;
		Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
		if(vtxsvc->isVertexValid()){
				double* dbv = vtxsvc->PrimaryVertex();
				double*  vv = vtxsvc->SigmaPrimaryVertex();
				xorigin.setX(dbv[0]);
				xorigin.setY(dbv[1]);
				xorigin.setZ(dbv[2]);
		}

		sigma0_mass=-100;
		sigma0_chis=-100;

		HepLorentzVector gamP4 = getP4(gammaShr,xorigin);
		p4_sigma0 = lmd->p() + gamP4;
		sigma0_mass = p4_sigma0.m();

		if(!m_use1c_sigma0)
		{
			if(Cut::Sigma0(sigma0_mass))
			{
				p4_sigma0_1c = p4_sigma0;
				sigma0_chis = -1;
				return true;
			}
			else
				return false;
		}
		
		if(!Cut::Sigma0(sigma0_mass)) return false;


		KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
		kmfit->init();
		kmfit->setChisqCut(200);
		kmfit->setIterNumber(5);
		kmfit->AddTrack(0, *lmd);
		kmfit->AddTrack(1, 0.0, gammaShr);
		kmfit->AddResonance(0, PDG::Sigma0, 0, 1);

		bool olmdq =kmfit->Fit(0);
		if(!olmdq) return false;

		kmfit->BuildVirtualParticle(0);
		sigma0_chis = kmfit->chisq(0);
		//if(sigma0_chis>100) return false;
		
		p4_gam = kmfit->pfit(1);
		p4_sigma0_1c = kmfit->pfit(0)+kmfit->pfit(1);
		wtp_sigma0_1c = kmfit->wVirtualTrack(0);

		return true;
}

bool Sigma0PionEta::Mass1c_Lmdc(WTrackParameter &sigma0, WTrackParameter &piTrk, WTrackParameter &etaTrk,double& lmdc_chis,HepLorentzVector& p4_lmdc_1c){
		lmdc_chis=-100;

		HepLorentzVector p4_Lab(0.011*m_ecms,0,0,m_ecms);
		HepLorentzVector p4_lmdc = sigma0.p() + piTrk.p() + etaTrk.p();
		double lmdc_mass = p4_lmdc.m();

		if(!Cut::Lambda_c(lmdc_mass)) return false;
		else 
		{
			p4_lmdc_1c = p4_lmdc;
			return true;
		}

		KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
		kmfit->init();
		kmfit->setChisqCut(INFINITY);
		kmfit->setIterNumber(5);
		kmfit->AddTrack(0, sigma0);
		kmfit->AddTrack(1, piTrk);
		kmfit->AddTrack(2, etaTrk);
		kmfit->AddMissTrack(3, PDG::Lambda_c);
		kmfit->AddFourMomentum(0, p4_Lab);

		bool olmdq =kmfit->Fit(0);
		if(!olmdq) return false;

		kmfit->BuildVirtualParticle(0);
		lmdc_chis = kmfit->chisq(0);
		
		p4_lmdc_1c = kmfit->pfit(0)+kmfit->pfit(1) + kmfit->pfit(2);

		return true;
}

HepLorentzVector Sigma0PionEta::getP4(RecEmcShower* gTrk, Hep3Vector origin){

		//double eraw = gTrk->energy();
		//double phi =  gTrk->phi();
		//double the =  gTrk->theta();

		//return HepLorentzVector( eraw * sin(the) * cos(phi),	eraw * sin(the) * sin(phi), 	eraw * cos(the),	eraw );


		Hep3Vector Gm_Vec(gTrk->x(), gTrk->y(), gTrk->z());
		Hep3Vector Gm_Mom = Gm_Vec - origin;
		Gm_Mom.setMag(gTrk->energy());
		HepLorentzVector pGm(Gm_Mom, gTrk->energy());
		return pGm; 
}


bool Sigma0PionEta::isGoodTrackForLambda(EvtRecTrack* trk){

		if( !trk->isMdcKalTrackValid()) {
				return false;
		}

		Hep3Vector xorigin(0,0,0);
		IVertexDbSvc*  vtxsvc;
		Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
		if(vtxsvc->isVertexValid()){
				double* dbv = vtxsvc->PrimaryVertex();
				double*  vv = vtxsvc->SigmaPrimaryVertex();
				xorigin.setX(dbv[0]);
				xorigin.setY(dbv[1]);
				xorigin.setZ(dbv[2]);
		}


		RecMdcTrack *mdcTrk = trk->mdcTrack();

		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();

		HepPoint3D point0(0.,0.,0.); // the initial point for MDC recosntruction
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
		VFHelix helixip(point0,a,Ea);
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
		double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
		double costheta=cos(mdcTrk->theta());

		if(fabs(Rvz0) < 20 && fabs(costheta)<m_CosThetaCut)  return true;

		return false;
}
