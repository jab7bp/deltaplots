#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>
#include <TF1.h>
#include "TStopwatch.h"

using namespace std::chrono;
#include "/w/halla-scshelf2102/sbs/jboyd/include/include_files.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/GEM_lookups.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/beam_variables.h"

Double_t fit_gaus(Double_t * x, Double_t *par){

	Double_t g = 0.0;

	g = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	
	return g;
}

template<typename T>
double VectorMean(std::vector<T> const& v){
	if(v.empty()){
		return 0;
	}
	return std::accumulate(v.begin(), v.end(), 0.0)/v.size();
}

bool single_run = false;
bool multi_run = true;

bool calc_W = true;
bool plot_dxdy = true;
bool dxdy_only = false;

bool use_heavy_cut = false;

//Run info and lookups
TString run_target = "LD2";
int kine = 9;
int sbsfieldscale = 70;

int runnum = lookup_parsed_runnums(run_target.Data(), kine, sbsfieldscale, 0);
// vector<int> runnum_vec = {13585, 13586, 13587, 13581, 13582, 13583, 13584};
vector<int> runnum_vec;
TString runs_string;
TString experiment = "gmn";
int pass = 1;


double E_beam = lookup_beam_energy_from_kine(kine); //Electron beam energy (electron energy) in GeV.q

double SBS_field = sbsfieldscale/100;

bool theta_pq_cut = true;
double theta_pq_p_thresh = 0.04;
double theta_pq_n_thresh = 0.04;

double p_recon, nu_recon, E_loss, E_corr, theta_pq_n, theta_pq_p;
double p_Beam, E_loss_outgoing;
double Eloss, E_beam_final, p_corr;

TString rootfile_dir;
// TString input_rootfile;

TFile *outfile;
TChain *TC = new TChain("T");
vector<TString> master_cut_vec;
TString master_cut_string;

TString elastic_yield_str = "";
TCut master_cut = "";

//Experimental Constants, Thresholds, cuts, etc DEFINITIONS
const double pi = TMath::Pi();
const double Mp = 0.938272; //Mass of proton [GeV]
const double Mn = 0.939565; //Mass of neutron [GeV]
double W2_mean; //Invariant Mass-squared (mean val) {With perfect optics W2 = Mp. Can be calculated run-by-run}
double W2_sigma; //Invariant Mass-squared sigma {Reasonable default/guess. Can be calculated run-by-run from W plot}

//SBS Magnet
const Double_t Dgap = 48.0*2.54/100.0; //about 1.22 m
const Double_t maxSBSfield = 1.26; //Tesla
const Double_t SBSdist = 2.25; //m
const Double_t dipGap = 1.22; //m
const Double_t sbsmaxfield = 3.1 * atan( 0.85/(11.0 - 2.25 - 1.22/2.0 ))/0.3/1.22/0.7;

//HCal constants and stuff
double tdiff = 510;		//Time difference between BBCal and HCal signals
double tdiff_max = 10;	//Maximum time difference from coincidences through tdctrig cut
double HCal_dist; 	//Distace from HCal face to target chamber
double HCal_theta;		//Theta angle for HCal from downstream beamline
const double hcalheight =  0.2897; //-0.2897; // Height of HCal above beamline

const Double_t sampfrac = 0.077; 	//Estimate of the sampling fraction from MC
const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t kNtrack = 100; // Reasonable max number of tracks per event
const Int_t kNtdc = 1000; // Reasonable max number of tdc signals per event
const Double_t Xi = -2.20; // Distance from beam center to top of HCal in m
const Double_t Xf = 1.47; // Distance from beam center to bottom of HCal in m
const Double_t Yi = -0.853; // Distance from beam center to opposite-beam side of HCal in m
const Double_t Yf = 0.853; // Distance from beam center to beam side of HCal in m

//Static Target Parameters
const double l_tgt = 0.15; // Length of the target (m)
const double rho_tgt = 0.0723; // Density of target (g/cc)
const double rho_Al = 2.7; // Density of aluminum windows (g/cc)
const double cell_diameter = 1.6*2.54; //cm, right now this is a guess
const double Ztgt = 1.0;
const double Atgt = 1.0;
const double Mmol_tgt = 1.008; //g/mol

//For energy-loss correction to beam energy:
const double dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const double dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const double uwallthick_LH2 = 0.0145; //cm
const double dwallthick_LH2 = 0.015; //cm
const double cellthick_LH2 = 0.02; //cm, this is a guess;
const double Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch

//Declare vars
Double_t atime[kNcell], row[kNcell], col[kNcell], tdctime[kNcell], cblkid[kNcell], cblke[kNcell];
Double_t nblk, nclus, SH_nclus, PS_nclus, hcal_x, hcal_y, hcal_e;

Double_t par[3];

double BB_dist, BB_theta;

Double_t bb_tr_p[kNtrack], bb_tr_px[kNtrack], bb_tr_py[kNtrack], bb_tr_pz[kNtrack];
Double_t bb_tr_vz[kNtrack];
Double_t bb_tr_n, bb_ps_x, bb_ps_y, bb_ps_e, bb_sh_x, bb_sh_y, bb_sh_e;

Double_t TDCT_id[kNtdc], TDCT_tdc[kNtdc], hodo_tmean[kNtdc]; 
Int_t TDCTndata;

Long64_t Nevents;

//INITIALIZE ALL HISTOGRAMS:
TH1D *h_Ep, *h_PS, *h_HCal_e, *h_SHPS, *h_W2, *h_W2recon, *h_vert, *h_W, *h_Wrecon;
TH1D *hin_Ep, *hin_PS, *hin_HCal_e, *hin_SHPS, *hin_W2, *hin_W2recon, *hin_W, *hin_dx;
TH1D *hin_tr_p, *hin_tr_p_wcut;
TH1D *hin_Wrecon, *hin_dxdy_wcut, *hin_dxdy_all;
TH1D *h_tr_p, *h_tr_p_wcut;

TH1D *h_dx, *h_dy, *h_Y;
TH2D *h_dxdy_all, *h_dxdy_wcut;

double Ep_center, Ep_sigma, PS_center, PS_sigma, PS_min, SHPS_center, SHPS_sigma, HCal_e;
double Ep_mean_pre, Ep_sig_pre, PS_mean_pre, PS_sig_pre, SHPS_mean_pre, SHPS_sig_pre, W2_mean_pre, W2_sig_pre, W_mean_pre, W_sig_pre;

double HCal_e_min, W2_fit_center, W2_fit_sigma, W_fit_center, W_fit_sigma;
double bb_tr_p_min = 0.0, bb_tr_p_wcut_min = 0.0, bb_tr_p_cut = 0.0;
vector<double> cuts_from_fits;
int target_int;

// TH1D *h_atime, *h_E, *h_E_cut, *h_E_exp, *h_vert, *h_W2, *h_W2recon, *h_dx, *h_dy, *h_KE_p, *h_KElow, *h_diff, *h_X, *h_Y;

TString pq_cut_String = "";

void calibrate_parsed_files(){
	auto total_time_start = high_resolution_clock::now();

	TStopwatch *StopWatch = new TStopwatch();

	gStyle->SetPalette(55);
	cout << "--------------------------------------" << endl;
	cout << "Analysis started. " << endl;
	cout << "--------------------------------------" << endl;

	if( run_target == "LH2" ){ target_int = 0; }
	if( run_target == "LD2" ){ target_int = 1; }

	for(int i = 0; i < lookup_parsed_runs_cnt(run_target.Data(), kine, sbsfieldscale); i++){
		runnum_vec.push_back(lookup_parsed_runnums(run_target.Data(), kine, sbsfieldscale, i));
	}

	cout << "Cut booleans: " << endl;
	cout << "Calc W: " << calc_W << endl;
	cout << "Heavy cut: " << use_heavy_cut << endl;

	if( theta_pq_cut ){
		pq_cut_String = "_pqCut";
	}


	if( !calc_W ){
		outfile = new TFile(Form("rootfiles/%s_SBS%i_mag%i%s_prime_calibration_parsed.root", run_target.Data(), kine, sbsfieldscale, pq_cut_String.Data()), "RECREATE");
	}
	if( calc_W ){
		outfile = new TFile(Form("rootfiles/%s_SBS%i_mag%i%s_full_calibration_parsed.root", run_target.Data(), kine, sbsfieldscale, pq_cut_String.Data()), "RECREATE");
	}		


	//Define Histograms:
	h_Ep = new TH1D("h_Ep", Form("E/p - SBS%i = %i%%, %s", kine, sbsfieldscale, run_target.Data()), 200, 0, 2);
	h_PS = new TH1D("h_PS", Form("Pre-Shower Clus. E - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0, 3);
	h_HCal_e = new TH1D("h_HCal_e", Form("HCal Clus. E - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 200, 0, 0.4);
	h_SHPS = new TH1D("h_SHPS", Form("SH + PS Clus. E - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0, 5);
	h_W = new TH1D( "h_W", Form("Invariant Mass, W (No Cuts) - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 200, 0.0, 2.0 );
	h_Wrecon = new TH1D( "h_Wrecon", Form("Recon. Invariant Mass, W  - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0.0, 3.0 );
	h_W2 = new TH1D( "h_W2", Form("Invariant Mass, W^{2} - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 400, 0.0, 4.0 );
	h_W2recon = new TH1D( "h_W2recon", Form("Reconstructed W^{2} - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 400, 0.0, 4.0 );
	h_vert = new TH1D( "h_vert", "Vertex Position; m", 200, -.1, .1 );	
	h_tr_p = new TH1D("h_tr_p", Form("Scattered electron track momentum - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0.0, 5.0);
	h_tr_p_wcut = new TH1D("h_tr_p_wcut", Form("Scattered electron track momentum (with W2 cut)- SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0.0, 5.0);

	BB_dist = lookup_BB_dist_by_kine(kine);
	BB_theta = (pi/180.0)*lookup_BB_angle_by_kine(kine, "deg");
	HCal_dist = lookup_HCal_dist_by_kine(kine);
	HCal_theta = (pi/180.0)*lookup_HCal_angle_by_kine(kine, "deg");

	if( plot_dxdy ){
		h_dxdy_wcut = new TH2D("h_dxdy_wcut",Form("HCal dxdy (W cuts) - SBS%i = %i%%, %s;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 400, -2, 2, 400, -2.0, 2.0 );
		h_dxdy_all = new TH2D("h_dxdy_all",Form("HCal dxdy (NO CUTS) - SBS%i = %i%%, %s;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()),125,-2,2,125,-4,6);
		h_dx = new TH1D( "h_dx", Form("HCal dx - SBS%i = %i%%, %s; m", kine, sbsfieldscale, run_target.Data()), 200, -4.0, 2.0 );
		h_dy = new TH1D( "h_dy", Form("HCal dy - SBS%i = %i%%, %s; m", kine, sbsfieldscale, run_target.Data()), 100, -1.2, 1.2 );
		h_Y = new TH1D( "h_Y", Form("HCal Y - SBS%i = %i%%, %s; m", kine, sbsfieldscale, run_target.Data()), 100, -1.2, 1.2 );
	}

	HCal_dist = lookup_HCal_dist_by_kine(kine); 	//Distace from HCal face to target chamber
	HCal_theta = (pi/180.0)*lookup_HCal_angle_by_kine(kine, "deg");		//Theta angle for HCal from downstream beamline
	cout << endl << "-----------------------------------" << endl;
	cout << "Run parameters: " << endl;
	cout << "Run: " << runnum << endl;
	cout << "Kinematic: SBS" << kine << endl;
	cout << "Target: " << run_target.Data() << endl;
	cout << "Beam Energy: " << E_beam << endl;
	cout << "SBS Field: " << sbsfieldscale << "%" << endl;
	cout << "-----------------------------------" << endl;
	cout << "BB angle: " << lookup_BB_angle_by_kine(kine, "deg") << endl;
	cout << "BB distance: " << lookup_BB_dist_by_kine(kine) << endl;
	cout << "SBS angle: " << lookup_SBS_angle_by_kine(kine, "deg") << endl;
	cout << "HCal angle: " << HCal_theta << endl;
	cout << "HCal distance: " << HCal_dist << endl;
	cout << "-----------------------------------" << endl << endl;;

	// rootfile_dir = Form("/lustre19/expphy/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass%i/SBS%i/%s/rootfiles", pass, kine, run_target.Data());
	if( kine == 11 ){
		rootfile_dir = "/volatile/halla/sbs/adr/Rootfiles/gmn_parsed/SBS11/pass1/";
	}
	else{
		rootfile_dir = Form("/volatile/halla/sbs/adr/Rootfiles/gmn_parsed/SBS%i/pass%i/", kine, pass);
		// rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/analysis_rootfiles/jboyd_parsed";
	}


	if( single_run ){
		cout << "Running in single run mode: " << runnum << endl;
		cout << "--------------------------------------" << endl;
		cout << "Adding files to TChain from: " << rootfile_dir.Data() << endl;
		// TC->Add(Form("%s/*%i*.root", rootfile_dir.Data(), runnum));
		if( kine == 11 ){
			TC->Add(Form("%s/gmn_parsed_SBS%i_targ%s_sbsmagscale%i.root", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale) );
		}
		else{
			TC->Add(Form("%s/gmn_parsed_SBS%i_targ%s_sbsmagscale%i.root", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale) );
			// TC->Add(Form("%s/gmn_parsed_%s_SBS%i_mag%i.root", rootfile_dir.Data(), run_target.Data(), kine, sbsfieldscale));		
		}
		
	}
	cout << "--------------------------------------" << endl;
	if( multi_run ){
		runnum = runnum_vec[0];
		cout << "Running in multi-run mode for runs: "  << endl;
		for(size_t run = 0; run < runnum_vec.size(); run++){
			cout << runnum_vec[run] << " ";
		}
		cout << endl;
		cout << "--------------------------------------" << endl;
		cout << "Adding files to TChain from: " << rootfile_dir.Data() << endl;
		for(size_t run = 0; run < runnum_vec.size(); run++){
			TC->Add(Form("%s/*%i*.root", rootfile_dir.Data(), runnum_vec[run]));

		}
		// runs_string = Form("Runs: %i thru %i", *min_element(runnum_vec.begin(), runnum_vec.end()), *max_element(runnum_vec.begin(), runnum_vec.end()));
	}
	cout << "Finished adding files to TChain. " << endl;
	cout << "--------------------------------------" << endl;

	// if( calc_W || plot_dxdy ){
	W2_mean = lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "W2_mean");
	W2_sigma = lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "W2_sigma");

	E_loss_outgoing = cell_diameter/2.0/sin(BB_theta)*rho_tgt*dEdx_tgt; //Should be about 1 MeV

	// }
	// else{
	// 	W2_mean = 0.92;
	// 	W2_sigma = 0.3;
	// }

	//DEFINE CUTS

	if( !calc_W ){
		master_cut_vec = {
			"sbs.hcal.nclus>0",
			"bb.ps.nclus>0",
			"bb.sh.nclus>0",
			"abs(bb.tr.vz[0])<=0.075",
			"bb.gem.track.nhits[0]>3",
			"bb.tr.n==1",
		};
	}
	if( calc_W ){
		master_cut_vec = {
			"sbs.hcal.nclus>0",
			"bb.ps.nclus>0",
			"bb.sh.nclus>0",
			"abs(bb.tr.vz)<=0.075",
			"bb.gem.track.nhits[0]>3",
			"bb.tr.n==1",
			// Form("bb.ps.e>%f", lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "PS_min")),
			// Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))>(%f)&&((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))<(%f)", lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "Ep_min"), lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "Ep_max")),
			// // Form("((abs(((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))))-%f)<%f", lookup_parsed_cut(runnum, "Ep"), lookup_parsed_cut(runnum, "Ep_sigma")),
			// Form("sbs.hcal.e>%f", lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "HCal_min")),
			// Form("(bb.sh.e+bb.ps.e)>%f", lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "SH_PS_min") - lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "SH_PS_max"))
		};
	}
	if( use_heavy_cut ){
		master_cut_vec = {
			"sbs.hcal.nclus>0",
			"bb.ps.nclus>0",
			"bb.sh.nclus>0",
			"abs(bb.tr.vz[0])<=0.075",
			"bb.gem.track.nhits[0]>4",
			"bb.tr.n==1",
			Form("bb.ps.e>%f", 0.3),
			Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))>(%f)&&((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))<(%f)", lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "Ep_min") + 0.15, lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "Ep_max") - .15),
			Form("sbs.hcal.e>%f", 0.03),
			Form("(bb.sh.e+bb.ps.e)>%f", lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "SH_PS_min") + 0.25)

		};
	}
	for(size_t cut = 0; cut < master_cut_vec.size(); cut++){
		if(cut == master_cut_vec.size() - 1){
			master_cut_string.Append(Form("%s", master_cut_vec[cut].Data()));
		}
		else{
			master_cut_string.Append(Form("%s%s", master_cut_vec[cut].Data(), "&&"));
		}
	}

	master_cut = Form("%s", master_cut_string.Data());

	cout << "--------------------------------------" << endl;
	cout << "Setting up branches... ";
	//Setup Branches
	TC->SetMakeClass(1); //Allows for viewing of general detector params from Event_Branch
	TC->SetBranchStatus( "*", 0 );
	TC->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );  
	TC->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
	TC->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
	TC->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
	TC->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
	TC->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
	TC->SetBranchStatus( "sbs.hcal.x", 1 );
	TC->SetBranchStatus( "sbs.hcal.y", 1 );
	TC->SetBranchStatus( "sbs.hcal.e", 1 );
	TC->SetBranchStatus( "sbs.hcal.nblk", 1 );
	TC->SetBranchStatus( "sbs.hcal.nclus", 1 );
	TC->SetBranchStatus( "bb.sh.nclus", 1 );
	TC->SetBranchStatus( "bb.ps.nclus", 1 );
	TC->SetBranchStatus( "bb.tr.n", 1 );
	TC->SetBranchStatus( "bb.tr.px", 1 );
	TC->SetBranchStatus( "bb.tr.py", 1 );
	TC->SetBranchStatus( "bb.tr.pz", 1 );
	TC->SetBranchStatus( "bb.tr.p", 1 );
	TC->SetBranchStatus( "bb.tr.vz", 1 );
	TC->SetBranchStatus( "bb.ps.e", 1 );
	TC->SetBranchStatus( "bb.ps.x", 1 );
	TC->SetBranchStatus( "bb.ps.y", 1 );
	TC->SetBranchStatus( "bb.sh.e", 1 );
	TC->SetBranchStatus( "bb.sh.x", 1 );
	TC->SetBranchStatus( "bb.sh.y", 1 );
	TC->SetBranchStatus( "bb.tdctrig.tdc", 1 );
	TC->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
	TC->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
	TC->SetBranchStatus( "bb.gem.track.nhits", 1);

	TC->SetBranchAddress( "sbs.hcal.clus_blk.atime", atime );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.row", row );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.col", col );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", tdctime );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke );
	TC->SetBranchAddress( "sbs.hcal.x", &hcal_x );
	TC->SetBranchAddress( "sbs.hcal.y", &hcal_y );
	TC->SetBranchAddress( "sbs.hcal.e", &hcal_e );
	TC->SetBranchAddress( "sbs.hcal.nblk", &nblk );
	TC->SetBranchAddress( "sbs.hcal.nclus", &nclus );
	TC->SetBranchAddress( "bb.sh.nclus", &SH_nclus );
	TC->SetBranchAddress( "bb.ps.nclus", &PS_nclus );
	TC->SetBranchAddress( "bb.tr.n", &bb_tr_n );
	TC->SetBranchAddress( "bb.tr.px", bb_tr_px );
	TC->SetBranchAddress( "bb.tr.py", bb_tr_py );
	TC->SetBranchAddress( "bb.tr.pz", bb_tr_pz );
	TC->SetBranchAddress( "bb.tr.vz", bb_tr_vz );
	TC->SetBranchAddress( "bb.tr.p", bb_tr_p );
	TC->SetBranchAddress( "bb.ps.e", &bb_ps_e );
	TC->SetBranchAddress( "bb.ps.x", &bb_ps_x );
	TC->SetBranchAddress( "bb.ps.y", &bb_ps_y );
	TC->SetBranchAddress( "bb.sh.e", &bb_sh_e );
	TC->SetBranchAddress( "bb.sh.x", &bb_sh_x );
	TC->SetBranchAddress( "bb.sh.y", &bb_sh_y );
	TC->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
	TC->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
	TC->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );
	cout << " done. " << endl;
	cout << "--------------------------------------" << endl;
	cout << "Applying master cut: " << endl;
	cout << master_cut_string.Data() << endl;
	cout << "--------------------------------------" << endl;

	TEventList *ev_list = new TEventList("ev_list", "Elastic Event List");

//Event list with cuts:
	TC->Draw(">>ev_list", master_cut);
	cout << "--------------------------------------" << endl;
	cout << "--------------------------------------" << endl;
	
	Nevents = ev_list->GetN();

	cout << "Number of events to analyze: " << Nevents << endl;
	cout << "--------------------------------------" << endl;
	cout << "--------------------------------------" << endl;
	cout << "Starting analysis loop on events..... " << endl;

	//Beam and detector vars including definition of HCal coordinate axis
	TLorentzVector Pbeam( 0, 0, E_beam, E_beam );
	TLorentzVector Ptarg( 0, 0, 0, 0.5*(Mp+Mn) ); //Average of proton and neutron rest mass
	TVector3 hcal_origin( -HCal_dist*sin(HCal_theta), 0, HCal_dist*cos(HCal_theta) );
	TVector3 hcal_zaxis = hcal_origin.Unit();
	TVector3 hcal_xaxis(0,-1,0);
	TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();

	Int_t elastic_yield = 0;

	int watch_cnt = 0;
	int five_percent = int(0.05*Nevents);
	vector<double> time_for_five;
	double average_time = 0.0, time_remaining;
	StopWatch->Start();

	for(Long64_t nevent = 0; nevent < Nevents; nevent++){
		TC->GetEntry( ev_list->GetEntry( nevent ));
		if( plot_dxdy ){
			elastic_yield_str = Form("Elastic yield = %i.", elastic_yield);
		}

		if( nevent%five_percent == 0){
			StopWatch->Stop();

			if( watch_cnt == 0){
				cout << "Evt: " << nevent <<"/" << Nevents << Form("(%.0f/100%%)", 100.0*double(1.0*nevent/Nevents)) << ". " << elastic_yield_str.Data() << endl;
			}

			if( watch_cnt > 0 ){
				time_for_five.push_back(StopWatch->RealTime());	
				average_time = VectorMean(time_for_five)/60;
				time_remaining = average_time*100*( 1.0 - double(nevent)/double(Nevents));
				cout << "Evt: " << nevent <<"/" << Nevents << Form("(%.0f/100%%)", 100.0*double(1.0*nevent/Nevents)) << ". " << elastic_yield_str.Data() << " Time left: " << time_remaining << " minutes." << endl;
			}
			StopWatch->Reset();
			StopWatch->Continue();
			watch_cnt++;
		}
		cout.flush();
	
		double bbcal_time=0.0, hcal_time=0.0;

		for(int ihit=0; ihit<TDCTndata; ihit++){
			if(TDCT_id[ihit]==5){
				bbcal_time=TDCT_tdc[ihit];
			}
			if(TDCT_id[ihit]==0){
				hcal_time=TDCT_tdc[ihit];
			}
		}

		double diff = hcal_time - bbcal_time; 
		if( fabs(diff - tdiff)>tdiff_max ){
			continue;
		}

		//Sanity check
      	if( (int)bb_tr_n!=1 ){
      		cout << "**************************************************************" << endl;
      		cout << "--------------------------------------------------------------" << endl;
      		cout << endl << endl << "WARNING: Total tracks not as expected from global cut. Check globalcut for errors." << endl << endl;
      		cout << "--------------------------------------------------------------" << endl;
      		cout << "**************************************************************" << endl;
      	} 
	
		Double_t Ep = (bb_ps_e + bb_sh_e)/(bb_tr_p[0]);
		Double_t PS = bb_ps_e;
		Double_t SHPS = bb_ps_e + bb_sh_e;
		Double_t HCal_e = hcal_e;

		h_Ep->Fill(Ep);
		h_PS->Fill(PS);
		h_SHPS->Fill(SHPS);
		h_HCal_e->Fill(hcal_e);
		h_tr_p->Fill(bb_tr_p[0]);

//----------------------------------------------------
//---------------INVARIANT MASS -----------------------
		if( calc_W ){

	      	TVector3 vertex( 0, 0, bb_tr_vz[0] ); // z location of vertex in hall coordinates
			TLorentzVector P_beam( 0, 0, E_beam_final, E_beam_final ); //Mass of e negligable
			TLorentzVector k_prime( bb_tr_px[0], bb_tr_py[0], bb_tr_pz[0], bb_tr_p[0] );
			TLorentzVector P_targ( 0, 0, 0, Mp );

			TLorentzVector kprime( bb_tr_px[0], bb_tr_py[0], bb_tr_pz[0], bb_tr_p[0] );
			TLorentzVector q = Pbeam - kprime; //Standard q-vector
			TVector3 qunit = q.Vect().Unit(); //q-vector direction

			TVector3 HCAL_zaxis( sin(-HCal_theta ), 0, cos(-HCal_theta) );
			TVector3 HCAL_xaxis( 0, -1, 0 );
			TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();

			TVector3 HCAL_origin = HCal_dist*HCAL_zaxis + hcalheight*HCAL_xaxis;

			TVector3 HCAL_pos = HCAL_origin + (hcal_x*HCAL_xaxis) + (hcal_y*HCAL_yaxis);

			Double_t sintersect = (hcal_origin-vertex).Dot( hcal_zaxis )/qunit.Dot( hcal_zaxis );
			TVector3 hcal_intersect = vertex + sintersect*qunit;
			Double_t x_expect = hcal_intersect.Dot( hcal_xaxis );
			Double_t y_expect = hcal_intersect.Dot( hcal_yaxis );

			Double_t W2recon = (Ptarg + q).M2();
			Double_t Wrecon = (Ptarg + q).M();
			Double_t E_ep = bb_tr_p[0]; // Obtain the scattered electron energy, neglect mass e
			Double_t p_ep = bb_tr_p[0]; // Obtain the magnitude of scattered electron momentum
			Double_t Q2 = 2*E_beam*E_ep*( 1-(bb_tr_pz[0]/p_ep) );
			Double_t nu = E_beam-E_ep; // Obtain energy transfer
			Double_t W2 = pow( Mp,2 )+2*Mp*nu-Q2; // Obtain W2 from Q2 and nu
			Double_t W = sqrt(W2);

		    h_vert->Fill(bb_tr_vz[0]);
		    h_W2->Fill(W2);
		    h_W2recon->Fill(W2recon);

		    h_W->Fill(W);
		    h_Wrecon->Fill(Wrecon);


			if( plot_dxdy ){
				h_dxdy_all->Fill( hcal_y - y_expect, hcal_x - x_expect );
			}

			if( plot_dxdy ){
				if( fabs(W2 - W2_mean)>W2_sigma ){
					continue;
				}


				// cout << "--------------------------------" << endl;
				// cout << "dx: " << hcal_y - y_expect << ", dy: " << hcal_x - x_expect << endl;
				// cout << "--------------------------------" << endl;
				h_dxdy_wcut->Fill( hcal_y - y_expect, hcal_x - x_expect );
				h_dx->Fill( hcal_x - x_expect );
				h_dy->Fill( hcal_y - y_expect );
				h_Y->Fill( hcal_y );
				h_tr_p_wcut->Fill(bb_tr_p[0]);
			}
	//--------------------------------------------------------------
	//Calculate theta pq variables
	      	//Reconstructed momentum, corrected for mean loss exiting the target
			p_recon = bb_tr_p[0] + E_loss_outgoing; 

			TLorentzVector k_prime_recon(p_recon*bb_tr_px[0]/bb_tr_p[0], p_recon*bb_tr_py[0]/bb_tr_p[0], p_recon*bb_tr_pz[0]/bb_tr_p[0], p_recon);
			TLorentzVector q_recon = P_beam - k_prime_recon;
			TVector3 qvec_recon = q_recon.Vect();

		//Calculate q vector as beam momentum - scattered k
			// TLorentzVector q = P_beam - k_prime;

		//Expected neutron direction
			TVector3 Neutron_Direction = (HCAL_pos - vertex).Unit();

		//Expected proton direction
			//Need to incorporate deflection due to SBS magnet
			double Bdl = sbsfieldscale*maxSBSfield*Dgap;
			double Proton_Deflection = tan( 0.3*Bdl/qvec_recon.Mag() )*(HCal_dist - (SBSdist + Dgap/2.0) ); 

			TVector3 Proton_Direction = (HCAL_pos + Proton_Deflection*HCAL_xaxis - vertex).Unit();

			theta_pq_n = acos( Neutron_Direction.Dot( qvec_recon.Unit() ) );
			theta_pq_p = acos( Proton_Direction.Dot( qvec_recon.Unit() ) );

			if( theta_pq_cut ){
				if( run_target == "LD2" ){
					if( theta_pq_p < theta_pq_p_thresh && theta_pq_n < theta_pq_n_thresh){
						continue;
					}
				}
			}

	//--------------------------------------------------------------
			
			elastic_yield++;

		}

	}

	cout << "---------------------------------------" << endl;
	cout << "-----Finished going through events-----" << endl;
	cout << "---------------------------------------" << endl;
	cout << "--- Writing to output file and saving ---" << endl;
	cout << "---------------------------------------" << endl;

	outfile->Write();
	outfile->Close();

	cout << "---------------------------------------" << endl;
	cout << "       Opening written file and plotting" << endl;
	cout << "---------------------------------------" << endl;

	TFile *infile;
	if( !calc_W ){
		infile = new TFile(Form("rootfiles/%s_SBS%i_mag%i_prime_calibration_parsed.root", run_target.Data(), kine, sbsfieldscale), "READ");
	}

	if( calc_W ){
		infile = new TFile(Form("rootfiles/%s_SBS%i_mag%i_full_calibration_parsed.root", run_target.Data(), kine, sbsfieldscale), "READ");
		hin_W2recon = static_cast<TH1D*>(infile->Get("h_W2recon"));
		hin_W2 = static_cast<TH1D*>(infile->Get("h_W2"));
		hin_Wrecon = static_cast<TH1D*>(infile->Get("h_Wrecon"));
		hin_W = static_cast<TH1D*>(infile->Get("h_W"));
		hin_tr_p = static_cast<TH1D*>(infile->Get("h_tr_p"));
		hin_tr_p_wcut = static_cast<TH1D*>(infile->Get("h_tr_p_wcut"));
		hin_dx = static_cast<TH1D*>(infile->Get("h_dx"));
	}
	if( plot_dxdy){
		hin_dxdy_wcut = static_cast<TH1D*>(infile->Get("h_dxdy_wcut"));
		hin_dxdy_all = static_cast<TH1D*>(infile->Get("h_dxdy_all"));
	}


	hin_Ep = static_cast<TH1D*>(infile->Get("h_Ep"));
	hin_PS = static_cast<TH1D*>(infile->Get("h_PS"));
	hin_SHPS = static_cast<TH1D*>(infile->Get("h_SHPS"));
	hin_HCal_e = static_cast<TH1D*>(infile->Get("h_HCal_e"));

//---------------------------------------
//---------------- Ep -----------------------

	TCanvas *c_Ep = new TCanvas("c_Ep", "E/p", 600, 500);
	hin_Ep->Draw();

	Ep_mean_pre = lookup_pre_parsed_cut(run_target, kine, "Ep");
	Ep_sig_pre = lookup_pre_parsed_cut(run_target, kine, "Ep_sigma");

	if( !calc_W ){
		TF1 *fit_Ep = new TF1("fit_Ep", fit_gaus, 0.1, 1.7, 3);
		fit_Ep->SetParName(0, "Ep Norm");
		fit_Ep->SetParName(1, "Ep Center");
		fit_Ep->SetParName(2, "Ep Sigma");

		fit_Ep->SetParLimits(0, 0.85*hin_Ep->GetMaximum(), hin_Ep->GetMaximum());
		fit_Ep->SetParLimits(1, (0.75)*Ep_mean_pre, (1.25)*Ep_mean_pre);
		fit_Ep->SetParLimits(2, (0.75)*Ep_sig_pre, 1.25*Ep_sig_pre);			

		hin_Ep->Fit("fit_Ep", "R");
		fit_Ep->Draw("same");
		Ep_center = fit_Ep->GetParameter(1);
		Ep_sigma = fit_Ep->GetParameter(2);

		par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;
	}


//---------------------------------------
//---------------- PS -----------------------

	TCanvas *c_PS = new TCanvas("c_PS", "PS", 600, 500);
	hin_PS->Draw();

	PS_mean_pre = lookup_pre_parsed_cut(run_target, kine, "PS");
	PS_sig_pre = lookup_pre_parsed_cut(run_target, kine, "PS_sigma");

	if( !calc_W ){
		TF1 *fit_PS = new TF1("fit_PS", fit_gaus, 0.15, 3, 3);
		fit_PS->SetParName(0, "PS Norm");
		fit_PS->SetParName(1, "PS Center");
		fit_PS->SetParName(2, "PS Sigma");

		fit_PS->SetParLimits(0, 0, hin_PS->GetMaximum());
		fit_PS->SetParLimits(1, (0.75)*PS_mean_pre, (1.25)*PS_mean_pre);
		fit_PS->SetParLimits(2, (0.75)*PS_sig_pre, (1.25)*PS_sig_pre);	
		
		hin_PS->Fit("fit_PS", "R");
		// fit_PS->Draw("same");
		PS_center = fit_PS->GetParameter(1);
		PS_sigma = fit_PS->GetParameter(2);

		par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;
	

		// cout << "Setting range of PS to : " << (3.0/300.0)*hin_PS->FindFirstBinAbove((0.75)*(hin_PS->GetMaximum())) << " - " << PS_center << endl;
		// hin_PS->GetXaxis()->SetRangeUser((3.0/300.0)*hin_PS->FindFirstBinAbove(0), PS_center);
		// cout << "---------------------------------------" << endl;
		// PS_min = (3.0/300.0)*hin_PS->GetMinimumBin();
		// cout << "PS_min before rounding: " << PS_min << endl;
		// PS_min = ceil(100*PS_min)/100.0;
		// cout << "PS_min set to: " << PS_min << endl;
		// cout << "---------------------------------------" << endl;

		PS_min = (3.0/300.0)*hin_PS->FindFirstBinAbove(0);
		cout << "PS_min set to: " << PS_min << endl;
		cout << "---------------------------------------" << endl;

		hin_PS->GetXaxis()->SetRangeUser(0, 3);
		TLine *tl_PS = new TLine(PS_min, 0, PS_min, hin_PS->GetMaximum());
		tl_PS->SetLineColor(6);
		tl_PS->Draw("same");
	}


//---------------------------------------
//---------------- SHPS -----------------------
	TCanvas *c_SHPS = new TCanvas("c_SHPS", "SHPS", 600, 500);
	hin_SHPS->Draw();

	SHPS_mean_pre = lookup_pre_parsed_cut(run_target, kine, "SH_PS_clus_e_cut");
	SHPS_sig_pre = lookup_pre_parsed_cut(run_target, kine, "SH_PS_sigma");

	if( !calc_W ){
		TF1 *fit_SHPS = new TF1("fit_SHPS", fit_gaus, 0.5, 5, 3);
		fit_SHPS->SetParName(0, "SHPS Norm");
		fit_SHPS->SetParName(1, "SHPS Center");
		fit_SHPS->SetParName(2, "SHPS Sigma");

		fit_SHPS->SetParLimits(0, 0, hin_SHPS->GetMaximum());
		fit_SHPS->SetParLimits(1, (0.75)*SHPS_mean_pre, (1.25)*SHPS_mean_pre);
		fit_SHPS->SetParLimits(2, (0.75)*SHPS_sig_pre, (1.25)*SHPS_sig_pre);

		hin_SHPS->Fit("fit_SHPS", "R");
		fit_SHPS->Draw("same");
		SHPS_center = fit_SHPS->GetParameter(1);
		SHPS_sigma = fit_SHPS->GetParameter(2);

		par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;
	}

//---------------------------------------
//---------------- HCal_e -----------------------

	TCanvas *c_HCal_e = new TCanvas("c_HCal_e", "HCal_e", 600, 500);
	hin_HCal_e->Draw();
	if( !calc_W ){
		cout << "---------------------------------------" << endl;
		HCal_e_min = (0.4/200.0)*hin_HCal_e->FindFirstBinAbove(0);
		cout << "HCal_e_min before rounding: " << HCal_e_min << endl;
		// HCal_e_min = ceil(100*HCal_e_min)/100.0;
		cout << "HCal_e_min set to: " << HCal_e_min << endl;
		cout << "---------------------------------------" << endl;

		TLine *tl_HCal_e = new TLine(HCal_e_min, 0, HCal_e_min, hin_HCal_e->GetMaximum());
		tl_HCal_e->SetLineColor(6);
		tl_HCal_e->Draw("same");

		par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;
	}
//---------------------------------------
//---------------- W2 Recon -----------------------	

	W2_mean_pre = lookup_pre_parsed_cut(run_target, kine, "W2");
	W2_sig_pre = lookup_pre_parsed_cut(run_target, kine, "W2_sigma");

	if( calc_W ){
		TCanvas *c_W2recon = new TCanvas("c_W2recon", "W2recon", 600, 500);
		hin_W2recon->Draw();
		TF1 *fit_W2recon = new TF1("fit_W2recon", fit_gaus, 0.5, 1.25, 3);
		fit_W2recon->SetParName(0, "W2recon Norm");
		fit_W2recon->SetParName(1, "W2recon Center");
		fit_W2recon->SetParName(2, "W2recon Sigma");

		fit_W2recon->SetParLimits(0, 0, hin_W2recon->GetMaximum());
		fit_W2recon->SetParLimits(1, (0.75)*W2_mean_pre, (1.25)*W2_mean_pre);
		fit_W2recon->SetParLimits(2, (0.75)*W2_sig_pre, (1.25)*W_sig_pre);

		hin_W2recon->Fit("fit_W2recon", "R");
		fit_W2recon->Draw("same");
		W2_fit_center = fit_W2recon->GetParameter(1);
		W2_fit_sigma = fit_W2recon->GetParameter(2);

		hin_W2->SetLineColor(6);
		hin_W2->SetLineStyle(6);
		hin_W2->Draw("same");
		par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;
	}

//---------------------------------------
//---------------- W Recon -----------------------	
	W_mean_pre = lookup_pre_parsed_cut(run_target, kine, "W");
	W_sig_pre = lookup_pre_parsed_cut(run_target, kine, "W_sigma");

	if( calc_W ){
		TCanvas *c_Wrecon = new TCanvas("c_Wrecon", "Wrecon", 600, 500);
		hin_Wrecon->Draw();
		TF1 *fit_Wrecon = new TF1("fit_Wrecon", fit_gaus, 0.5, 1.25, 3);
		fit_Wrecon->SetParName(0, "Wrecon Norm");
		fit_Wrecon->SetParName(1, "Wrecon Center");
		fit_Wrecon->SetParName(2, "Wrecon Sigma");

		fit_Wrecon->SetParLimits(0, 0, hin_Wrecon->GetMaximum());
		fit_Wrecon->SetParLimits(1, (0.75)*W_mean_pre, (1.25)*W_mean_pre);
		fit_Wrecon->SetParLimits(2, (0.75)*W_sig_pre, (1.25)*W_sig_pre);

		hin_Wrecon->Fit("fit_Wrecon", "R");
		fit_Wrecon->Draw("same");
		W_fit_center = fit_Wrecon->GetParameter(1);
		W_fit_sigma = fit_Wrecon->GetParameter(2);

		hin_W->SetLineColor(6);
		hin_W->SetLineStyle(6);
		hin_W->Draw("same");
		par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;
	}

//---------------------------------------
//---------------- dx -----------------------	
	if( calc_W ){
		TCanvas *c_dx = new TCanvas("c_dx", "dx", 600, 500);
		hin_dx->Draw();

	}

//TRACK MOMENTUM CUTS
	// if( calc_W	){

	// 	bb_tr_p_min = (0.01)*hin_tr_p->FindFirstBinAbove(100);
	// 	bb_tr_p_wcut_min = (0.01)*hin_tr_p_wcut->FindFirstBinAbove(100);
	// 	bb_tr_p_cut = min( bb_tr_p_min, bb_tr_p_wcut_min );
	// }


//---------------------------------------
//---------------- dxdy -----------------------	
	if( plot_dxdy ){
		TCanvas *c_dxdy_wcut = new TCanvas("c_dxdy_wcut", "dxdy_wcut", 600, 500);
		hin_dxdy_wcut->Draw("colz");
	}	

//---------------------------------------
//---------------- CUT VECTOR -----------------------	

	cuts_from_fits.push_back(target_int);
	cuts_from_fits.push_back(kine);
	cuts_from_fits.push_back(sbsfieldscale);
	cuts_from_fits.push_back(PS_min);
	cuts_from_fits.push_back(SHPS_center);
	cuts_from_fits.push_back(SHPS_sigma);
	cuts_from_fits.push_back(HCal_e_min);
	cuts_from_fits.push_back(Ep_center);
	cuts_from_fits.push_back(Ep_sigma);


	if( calc_W ){
		cuts_from_fits.push_back(W2_fit_center);
		cuts_from_fits.push_back(W2_fit_sigma);
		cuts_from_fits.push_back(W_fit_center);
		cuts_from_fits.push_back(W_fit_sigma);
		// cuts_from_fits.push_back(bb_tr_p_cut);
	}


	cout << "---------------------------------------" << endl;
	cout << "---------------------------------------" << endl << endl;	
	cout << "Vector with cuts: " << endl;
	cout << "{";
	for(size_t cut = 0; cut < cuts_from_fits.size(); cut++){
		if( cut < cuts_from_fits.size() -1 ){
			cout << cuts_from_fits[cut] << ", ";
		}
		if( cut == cuts_from_fits.size() -1 ){
			cout << cuts_from_fits[cut] << " }" << endl;
		}
	}


	// TFile *inputfile = new TFile(Form("rootfiles/%i_dxdy", runnum), "READ");



	// TH1D *h_atime = new TH1D( "atime", Form("HCal ADC Time, All Channels - SBS = %i%%, %s, Run %i; ns", sbsfieldscale, run_target.Data(), runnum), 160, 0, 160 );
	// TH1D *h_E = new TH1D( "E", Form("HCal Cluster E, All Channels - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 500, 0, 0.5 );
	// TH1D *h_E_cut = new TH1D( "E_cut", Form("HCal Cluster E All Cuts, All Channels - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 100, 0, 0.5 );
	// TH1D *h_E_exp = new TH1D( "E_exp", Form("Expected Energy Dep in HCal - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 100, 0, 0.2 );
	// TH1D *h_vert = new TH1D( "vert", Form("Vertex Position - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 200, -1.0, 1.0 );
	// TH1D *hdx = new TH1D( "dx", Form("HCal dx - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 200, -4.0, 2.0 );
	// TH1D *hdy = new TH1D( "dy", Form("HCal dy - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 100, -1.2, 1.2 );
	// TH1D *hKE_p = new TH1D( "KE_p", Form("Scattered Proton Kinetic Energy - SBS = %i%%, %s, Run %i", sbsfieldscale, run_target.Data(), runnum), 500, 0.0, 5.0 );
	// TH1D *hKElow = new TH1D( "KElow", Form("Lowest Elastic E Sampled in HCal - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 500, 0.0, 0.2 );
	// TH1D *hDiff = new TH1D( "hDiff",Form("HCal time - BBCal time - SBS = %i%%, %s, Run %i; ns", sbsfieldscale, run_target.Data(), runnum), 1300, -500, 800 );
	// TH1D *hX = new TH1D( "X", Form("HCal X - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 100, -4.0, 2.0 );
	// TH1D *hY = new TH1D( "Y", Form("HCal Y - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 100, -1.2, 1.2 );
	// TH2D *h2_CvCh = new TH2D( "CvCh", Form("HCal Coeff Single Block Clusters - SBS = %i%%, %s, Run %i; channel, GeV", sbsfieldscale, run_target.Data(), runnum), kNcell, 0, kNcell, 200, 0, 1.0 );
	// TH2D *h2_EvCh = new TH2D( "EvCh", Form("HCal Cluster E Single Block Clusters - SBS = %i%%, %s, Run %i; channel, GeV", sbsfieldscale, run_target.Data(), runnum), kNcell, 0, kNcell, 50, 0, 0.5 );
	// TH2D *h2_dxdy = new TH2D("dxdy",Form("Delta Plot (q-proj. dxdy) - SBS = %i%%, %s, Run %i;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", sbsfieldscale, run_target.Data(), runnum), 50, -2, 2, 100, -4.0, 2.0 );
	// TH2D *h2_dxdy_all = new TH2D("dxdy_all",Form("Delta Plot (NO CUT) - SBS = %i%%, %s, Run %i;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", sbsfieldscale, run_target.Data(), runnum),125,-2,2,125,-4,6);
	// TH2D *h2_dxVE = new TH2D("dxVE",Form("dxVE - SBS = %i%%, %s, Run %i;x_{HCAL}-x_{expect} (m); GeV", sbsfieldscale, run_target.Data(), runnum), 100, -4.0, 2.0, 100, 0, 0.5 );
	// TH2D *h2_rowcol = new TH2D( "hrowcol", Form("HCal Block Position Elastics, HCal - SBS = %i%%, %s, Run %i; Col; Row", sbsfieldscale, run_target.Data(), runnum), kNcols, 0, kNcols, kNrows, -kNrows, 0 );
cout << "------------------------------------------------------------------"<< endl;
cout << "                       ANALYSIS FINISHED" << endl;
cout << "------------------------------------------------------------------"<< endl;
	cout << "Run parameters: " << endl;
	cout << "Run: " << runnum << endl;
	cout << "Kinematic: SBS" << kine << endl;
	cout << "Target: " << run_target.Data() << endl;
	cout << "Beam Energy: " << E_beam << endl;
	cout << "SBS Field: " << sbsfieldscale << "%" << endl;
	cout << "-----------------------------------" << endl;
	cout << "BB angle: " << lookup_BB_angle(runnum) << endl;
	cout << "SBS angle: " << lookup_SBS_angle(runnum) << endl;
	cout << "HCal angle: " << HCal_theta << endl;
	cout << "HCal distance: " << HCal_dist << endl;
	cout << "-----------------------------------" << endl << endl;
	cout << "Elastic yield: " << elastic_yield << endl << endl;
	cout << "---------------------------------------" << endl << endl;	
	cout << "Vector with cuts: " << endl;

	if( !calc_W ){
		cout << "{ Target, Kine, mag, PS_cut, SH_PS_mean, SH_PS_sigma, HCal_cut, Ep_mean, Ep_sigma }" << endl << endl;
		cout << "{ ";
		for(size_t cut = 0; cut < cuts_from_fits.size(); cut++){
			if( cut < cuts_from_fits.size() -1 ){
				cout << cuts_from_fits[cut] << ", ";
			}
			if( cut == cuts_from_fits.size() -1 ){
				cout << cuts_from_fits[cut] << " }" << endl;
			}
		}
	}
	if( calc_W ){
		cout << "{ Target, Kine, mag, PS_cut, SH_PS_mean, SH_PS_sigma, HCal_cut, Ep_mean, Ep_sigmaW2, W2_sigma, W, W_sigma, bb_tr_p }" << endl << endl;
		cout << "{";
		cout << target_int << ", " << kine << ", " << sbsfieldscale << ", "<< lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "Ep_mean") << ", " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "Ep_sigma") << ", " ;
		cout << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "PS_min") << ", " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "SH_PS_mean") << ", ";
		cout << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "SH_PS_sigma") << ", " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "HCal_min") << ", ";
		cout  << W2_fit_center <<", " << W2_fit_sigma << ", ";
		cout << W_fit_center << ", " << W_fit_sigma << ", " << bb_tr_p_cut << " }" << endl;
	}

	cout << endl << "Output file: " << endl;
	cout << outfile->GetName() << endl << endl;

	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;


}

