#ifndef SigmapPi0Pi0_Header
#define SigmapPi0Pi0_Header

#include "GaudiKernel/Algorithm.h"
//you can add oher necessary header files
#include "GaudiKernel/NTuple.h"
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "BestDTagSvc/BestDTagSvc.h"
#include "GaudiKernel/Algorithm.h"
#include "VertexFit/WTrackParameter.h"

//
//namespace
//
using CLHEP::HepLorentzVector;

class SigmapPi0Pi0:public Algorithm {
  public:
    SigmapPi0Pi0(const std::string& name, ISvcLocator* pSvcLocator);
    ~SigmapPi0Pi0();
    StatusCode initialize();
    StatusCode beginRun();   
    StatusCode execute();
    StatusCode endRun();
    StatusCode finalize();

    bool isGoodTrk(EvtRecTrackIterator itTrk, double &vz, double &vxy);
    bool isProton(EvtRecTrackIterator itTrk);
    bool isPion (EvtRecTrackIterator itTrk);
    bool isKaon (EvtRecTrackIterator itTrk);
    bool isGoodShower(EvtRecTrack* trk, double &);

    // Not Lambda_c !!!!
    bool isGoodLambda(RecMdcKalTrack* ppTrk, RecMdcKalTrack* pimTrk, 
                double& lmd_1chis, double& lmd_2chis, double& lmd_lchue,
                WTrackParameter &wvlmd, double& lmd_mass);
    //Notice that the first parameter is from member fuction isGoodLambda
    bool isGoodSigma0(WTrackParameter* virtual_lmd_trk, RecEmcShower* gammaShr, double& sigma0_mass,HepLorentzVector& p4_sigma0,double& sigma0_chis,HepLorentzVector& p4_sigma0_1c, HepLorentzVector&,WTrackParameter &);

    bool isGoodPi0(RecEmcShower *shr1,RecEmcShower *shr2,double& pi0_mass,HepLorentzVector& p4_pi0,double& pi0_chis,HepLorentzVector& p4_pi0_1c, HepLorentzVector&, HepLorentzVector&,WTrackParameter &);

	bool    isGoodTrackForLambda(EvtRecTrack* trk);

	bool Mass1c_Lmdc(WTrackParameter &sigma0, WTrackParameter &piTrk, WTrackParameter &pi0Trk,double& lmdc_chis,HepLorentzVector& p4_lmdc_1c);

  private:

     int m_irun;
     bool m_debug;
     bool m_use_Total_TOF;
     bool m_BestCandidate;
     bool m_use1c_sigma0;
	 bool m_checktotal;

    double m_vr0cut;
    double m_vz0cut;
    double m_skim;
    double m_beamE;
    double m_ecms;
    double m_prob_cut;
    double m_CosThetaCut;

    double m_minEnergy;
    double m_gammaAngleCut;
    double m_gammatlCut;  
    double m_gammathCut;  
    double m_maxCosThetaBarrel;
    double m_minCosThetaEndcap;
    double m_maxCosThetaEndcap;
    double m_minEndcapEnergy;

    //NTuple::Item<int> m_idxmc;
    //NTuple::Array<int>m_pdgid;
    //NTuple::Array<int>m_motheridx;
    
//Because there are Lambda_c+ and Lambda_c- thus we make it twice( 0 for lambda_c+ and 1 for lambda_c-)
    NTuple::Tuple* m_tuple;   
	NTuple::Tuple *m_tuple_mc;
    NTuple::Item<int> m_run;
    NTuple::Item<int> m_event;
    NTuple::Item<int>m_flag;
    NTuple::Item<int>m_index;

	NTuple::Item<int> m_mode1;
	NTuple::Item<int> m_mode2;
	NTuple::Item<int> m_mode3;
	NTuple::Item<int> m_ndaughterAp;
	NTuple::Array<int> m_Ap_id;
	NTuple::Matrix<double> m_Ap_ptruth;
	NTuple::Item<int> m_ndaughterAm;
	NTuple::Array<int> m_Am_id;
	NTuple::Matrix<double> m_Am_ptruth;
	NTuple::Item<int> m_idxmc_p;
	NTuple::Item<int> m_idxmc_m;
	NTuple::Array<int> m_pdgid_p;
	NTuple::Array<int> m_motheridx_p;
	NTuple::Array<int> m_pdgid_m;
	NTuple::Array<int> m_motheridx_m;

	NTuple::Item<int> mc_run;
	NTuple::Item<int> mc_event;
	NTuple::Item<int> mc_mode1;
	NTuple::Item<int> mc_mode2;
	NTuple::Item<int> mc_mode3;
	NTuple::Item<int> mc_ndaughterAp;
	NTuple::Array<int> mc_Ap_id;
	NTuple::Matrix<double> mc_Ap_ptruth;
	NTuple::Item<int> mc_ndaughterAm;
	NTuple::Array<int> mc_Am_id;
	NTuple::Matrix<double> mc_Am_ptruth;

	NTuple::Array<int>m_charge;
	NTuple::Item<int>m_ngood;
	NTuple::Array<double>m_mass;
	NTuple::Array<double>m_ebeam;
	NTuple::Array<double>m_mbc;
	NTuple::Array<double>m_deltaE;
	NTuple::Array<double>m_deltaE_lmd_pip_pi0;

	NTuple::Array<double> m_lambda_p4;
	NTuple::Array<double>  m_lambda_mass;
	NTuple::Array<double>  m_lambda_chis1;
	NTuple::Array<double>  m_lambda_lchue;

	NTuple::Array<double>  m_lmdc_chis;

	NTuple::Matrix<double> m_gam0_p4;
	NTuple::Matrix<double> m_gam1_p4;
	NTuple::Matrix<double> m_gam2_p4;
	NTuple::Matrix<double> m_gam3_p4;
	NTuple::Matrix<double> m_gam4_p4;
	NTuple::Matrix<double> m_gam5_p4;
	NTuple::Matrix<double> m_gam0_p4_before;
	NTuple::Matrix<double> m_gam1_p4_before;
	NTuple::Matrix<double> m_gam2_p4_before;
	NTuple::Matrix<double> m_gam3_p4_before;
	NTuple::Matrix<double> m_gam4_p4_before;
	NTuple::Matrix<double> m_gam5_p4_before;
	NTuple::Array<double>  m_ang_gam_antiproton;

	void addItem();

	HepLorentzVector getP4(RecEmcShower *gTrk, Hep3Vector origin);
bool isGoodSigmap(RecMdcKalTrack *ppTrk, RecEmcShower *shr1, RecEmcShower *shr2);
  protected:

};

typedef struct{
	int index;
	double ang_antiproton;
}Gam_info;

//add your inline methods

//

#endif//SigmapPi0Pi0_Header
