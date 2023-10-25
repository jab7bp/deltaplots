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
#include "/w/halla-scshelf2102/sbs/jboyd/include/utility_functions.h"
// #include "/w/halla-scshelf2102/sbs/jboyd/include/MC_lookups.h"

bool print_extra = false;

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
bool use_parsed = false;

bool calc_W = true;
bool use_heavy_cut = false;

bool correct_beam_energy = true;
bool use_acceptance_cut = true;
bool fiducial_cut = true;

bool calibrate = true;

//Bool that will switch depending on particle location
bool apply_fcut = false;

bool apply_adc_diff_time_cut = false;
bool apply_ADC_cut = false;
TString ADC_timing_string = "";

//Run info and lookups
TString run_target = "LD2";
int kine = 8;
int sbsfieldscale = 70;
double sbsfieldscale_Frac = sbsfieldscale/100.0;

//*-*-*-*-*-*-*-*-*-*-*
//BEST CLUSTER BOOLEANS
//hcal_cluster_minimize
//Which method to use for "Best cluster selection"
//Choices are: "coin_time_and_maxE", "coin_time"
TString hcal_cluster_minimize = "coin_time_and_maxE";
TString best_cluster_OnOff_indicated = "";

vector<double> best_cluster_OrigIndex_Energy_Index_Timing_Score = {0};
bool sort_hcal_cluster_energy = true;
bool use_best_cluster = true;
bool use_scoring = false;
int num_best_cluster_types = 3; ////DON"T TOUCH!!
//Choose the best cluster type to prioritize:
// 1 --> Highest energy cluster
// 2 --> Cluster to minimize "dx"
// 3 --> Cluster to minimize ADC time different

int priority_best_cluster_type = 2;

bool theta_pq_cut = false;
double theta_pq_p_thresh = 0.18; //0.05;
double theta_pq_n_thresh = 0.18; //0.04;


//*-*-*-*-*-*-*-*-*-*-*

int runnum = lookup_parsed_runnums(run_target.Data(), kine, sbsfieldscale, 0);

vector<int> runnum_vec;

TString experiment = "gmn";
int pass;

//Experimental Lookup Parameters
double E_beam = lookup_beam_energy_from_kine(kine); //Electron beam energy (electron energy) in GeV
double SBS_field = sbsfieldscale; //Strength (in percentage) of SBS magnet

double BB_dist, BB_theta, W_mean, W_sigma;
double dx_p, dx_p_sigma, dy_p, dy_p_sigma, dx_n, dx_n_sigma, dy_n, dy_n_sigma, dx_pn_max, dx_pn_mean_center;

TString rootfile_dir;
// TString input_rootfile;

TFile *outfile;
TChain *TC = new TChain("T");
vector<TString> master_cut_vec;
TString master_cut_string;

TString elastic_yield_str = "";
TCut master_cut = "";

TTreeFormula *master_cut_formula;

double dx_p_scale = 1.0;
double dx_n_scale = 1.0;
double dy_scale = 1.0;

//Experimental Constants, Thresholds, cuts, etc DEFINITIONS
const double pi = TMath::Pi();
const double Mp = 0.938272; //Mass of proton [GeV]
const double Mn = 0.939565; //Mass of neutron [GeV]
const double Me = 0.00051; //Mass of electron [GeV]

//SBS Magnet
const Double_t Dgap = 48.0*2.54/100.0; //about 1.22 m
const Double_t maxSBSfield = 1.26; //Tesla
const Double_t SBSdist = 2.25; //m
const Double_t dipGap = 1.22; //m
const Double_t sbsmaxfield = 3.1 * atan( 0.85/(11.0 - 2.25 - 1.22/2.0 ))/0.3/1.22/0.7;

double W2_mean; //Invariant Mass-squared (mean val) {With perfect optics W2 = Mp. Can be calculated run-by-run}
double W2_sigma; //Invariant Mass-squared sigma {Reasonable default/guess. Can be calculated run-by-run from W plot}

//HCal constants and stuff
double tdiff = 510;		//Time difference between BBCal and HCal signals
double tdiff_max = 10;	//Maximum time difference from coincidences through tdctrig cut
double HCal_dist; 	//Distace from HCal face to target chamber
double HCal_theta;		//Theta angle for HCal from downstream beamline
double scint_intersect, x_expected_HCal, y_expected_HCal;
double ADC_time_min, ADC_time_max, ADC_diff_time_min, ADC_diff_time_max, ADC_time_mean;

//Scattered kinematics
double e_prime_theta; //Scattered electron theta angle
double e_prime_phi; //Scattered electron phi angle
double p_el, nu, pp, nucleon_theta, nucleon_phi, E_ep, p_ep, Q2, W, W2, E_pp, E_nucleon, KE_p, dx, dy;

//Static Detector Parameters
const int maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const int maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
// const double hcal_height = -0.2897; // Height of HCal above beamline
double hcal_height;

const Double_t sampfrac = 0.077; 	//Estimate of the sampling fraction from MC
const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t kNtrack = 100; // Reasonable max number of tracks per event
const Int_t kNtdc = 1000; // Reasonable max number of tdc signals per event
const Int_t max_clus = 10;

//Values have ben updated -- Known as of Oct 1, 2023......
const Double_t Xi = -2.655; //Distance from beam center to top of HCal in meters, from database
const Double_t Xf = 1.155; //Distance from beam center to bottom of HCal in meters, from database
const Double_t Yi = -0.92964; //Distance from beam center to opposite-beam side of HCal in meters, from MC database
const Double_t Yf = 0.92964; //Distance from beam center to beam side of HCal in meters, from MC database

//PRevious values:
// const Double_t Xi = -2.20; // Distance from beam center to top of HCal in m
// const Double_t Xf = 1.47; // Distance from beam center to bottom of HCal in m
// const Double_t Yi = -0.853; // Distance from beam center to opposite-beam side of HCal in m
// const Double_t Yf = 0.853; // Distance from beam center to beam side of HCal in m

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

double p_recon, nu_recon, E_loss, E_corr, theta_pq_n, theta_pq_p;

int useAlshield = 0;

//Declare vars
Double_t atime[kNcell], row[kNcell], col[kNcell], tdctime[kNcell], cblkid[kNcell], cblke[kNcell];
Double_t nblk, nclus, SH_nclus, PS_nclus, hcal_x, hcal_y, hcal_e;
Double_t hcal_clus_e[max_clus], hcal_clus_x[max_clus], hcal_clus_y[max_clus], hcal_clus_atime[max_clus], hcal_clus_tdctime[max_clus];
Array1DValueWithIndex hcal_clus_e_sorted[maxTracks];
// Array1DScoredWithIndex hcal_clus_scored[maxTracks];
vector<vector<int>> hcal_clus_scored;
//{ SCORE, HIGHEST_E_CLUS_INDEX, BEST_TIMING_INDEX, BEST_DXDY_INDEX, BEST_THETA_PQ_INDEX}
vector<int> hcal_clus_row = {0, -1, -1, -1};


Double_t dx_bestcluster, dy_bestcluster, HCal_ADC_time_bestcluster;
Double_t dx_cluster_final, dy_cluster_final, HCal_ADC_time_final;
Double_t hcal_clus_id, hcal_clus_nblk;

Double_t par[3];

double bb_tr_p[maxTracks], bb_tr_px[maxTracks], bb_tr_py[maxTracks], bb_tr_pz[maxTracks];
double bb_tr_vx[maxTracks], bb_tr_vy[maxTracks], bb_tr_vz[maxTracks], bb_tr_chi2[maxTracks];
double bb_fp_x[maxTracks], bb_fp_y[maxTracks], bb_fp_th[maxTracks], bb_fp_ph[maxTracks];
double bb_tgt_x[maxTracks], bb_tgt_y[maxTracks], bb_tgt_th[maxTracks], bb_tgt_ph[maxTracks];
double hcal_clusblk_ADC_time[max_clus]; //Maximum number of blocks in a cluster is 15 as per S. Seeds
double bb_sh_atimeblk;
double hcal_clus_ADCtime_diff = 0.0;
double bb_tr_n, bb_ps_x, bb_ps_y, bb_ps_e, bb_sh_x, bb_sh_y, bb_sh_e;

Double_t TDCT_id[kNtdc], TDCT_tdc[kNtdc], hodo_tmean[kNtdc]; 
Int_t TDCTndata;

Long64_t Nevents;

//INITIALIZE ALL HISTOGRAMS:
TH1D *h_atime, *h_W, *h_W2, *h_W2recon, *h_KE_p, *h_KE_low, *h_Diff, *h_X, *h_Y, *h_E_eloss, *h_hcal_clusblk_ADC_time, *h_hcal_clusblk_ADC_time_diff, *h_bb_sh_atimeblk, *h_HCal_e;
TH1D *h_hcal_clusblk_ADC_time_cut, *h_hcal_clusblk_ADC_time_diff_cut;
TH1D *h_W_cut, *h_W_fcut, *h_vz_cut, *h_Ep, *h_PS, *h_SHPS;
TH1D *h_SHPS_wcut, *h_Ep_wcut;
TH1D *h_theta_pq_n, *h_theta_pq_p, *h_theta_pq_n_cut, *h_theta_pq_p_cut, *h_theta_pq_n_anticut, *h_theta_pq_p_anticut;

TH1D *h_Q2, *h_E_ep, *h_E_pp;
TH1D *h_dy, *h_dy_cut, *h_dy_wcut, *h_dx, *h_dx_cut, *h_dx_wcut, *h_dx_fcut, *h_dx_wcut_fcut, *h_dy_wcut_fcut;
TH1D *h_tr_p, *h_tr_p_wcut;
TH1D *h_Nevents, *h_eclus_sel;
TH1D *h_hcal_e, *h_hcal_clus_e, *h_hcal_clus_e_sorted;
TH1D *h_hcal_clus_e_indexed, *h_hcal_clus_e_sorted_indexed;

TH2D *h2_hcal_e_V_clus_e, *h2_hcal_e_V_clus_e_sorted, *h2_hcal_clus_e_V_hcal_clus_e_sorted;
 
TH2D *h_E_ecorr_vs_vert;
TH2D *h_dxdy, *h_dxdy_cut, *h_dxdy_wcut, *h_dxdy_ncut, *h_dxdy_pcut, *h_dxdy_fcut, *h_dxdy_wcut_fcut, *h_dxdy_wcut_fcut_ADCtiming, *h_dxdy_wcut_fcut_AntiADCtiming;
TH2D *h_dxdy_wcut_2multfcut, *h_dxdy_wcut_15multfcut;
TH2D *h_xy, *h_xy_cut, *h_xy_fcut, *h_xy_cut_p, *h_xy_cut_n, *h_PAngleCorr_theta, *h_PAngleCorr_phi;

double n_integral, p_integral, n_center, n_sigma, p_center, p_sigma;
double p_Beam, E_loss_outgoing;
double Eloss, E_beam_final, p_corr;
int n_counts, p_counts, elastic_yield, dxdy_cnt, dxdy_wcut_cnt;
int n_hcal_clusblk_atime_cut = 0, n_adc_diff_time_cnt = 0;


double Ep_sig_mult, SH_PS_sig_mult, Ep_min, Ep_max, SH_PS_min;
double bb_tr_p_min, bb_tr_p_wcut_min, bb_tr_p_cut;

TString pq_cut_String = "";

int Nclusters = 0, badHcalESort = 0;
Int_t Nhcal_clus_id;

double hcal_e_ADC_passed_maxElement;
int hcal_e_ADC_passed_maxElementIndex;

double fiducial_active_area_xmax, fiducial_active_area_xmin, fiducial_active_area_ymax, fiducial_active_area_ymin;


void dxdy_parsed_files_byKine(){

	auto total_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	gStyle->SetPalette(55);
	cout << "--------------------------------------" << endl;
	cout << "Analysis started. " << endl;
	cout << "--------------------------------------" << endl;
	for(int i = 0; i < lookup_parsed_runs_cnt(run_target.Data(), kine, sbsfieldscale); i++){
		runnum_vec.push_back(lookup_parsed_runnums(run_target.Data(), kine, sbsfieldscale, i));
	}

	if( kine == 4 ){
		pass = 0;
	}
	if( kine == 8 ){
		pass = 1;
	}

	if( theta_pq_cut ){
		pq_cut_String = "_pqCut";
	}
	TString parsed_sel_string = "";

	if( use_parsed ){
		parsed_sel_string = "_parsed";
	}
	if( !use_parsed){
		parsed_sel_string = "_NOTparsed";
	}
	if( apply_ADC_cut || apply_adc_diff_time_cut ){
		ADC_timing_string = "_ADCtiming";
	}

	if( use_best_cluster ){
		best_cluster_OnOff_indicated = Form("_%s", hcal_cluster_minimize.Data() );
	}

	outfile = new TFile(Form("rootfiles/%s_SBS%i_mag%i%s_dxdy%s_elastics_only_trPfact_100_bestCluster_%i%s%s_08_10_2023.root", run_target.Data(), kine, sbsfieldscale, pq_cut_String.Data(), parsed_sel_string.Data(), use_best_cluster,  best_cluster_OnOff_indicated.Data(), ADC_timing_string.Data() ), "RECREATE");		

Int_t nBins_x_dxdy = 300;
Double_t xmin_dxdy = -1.5;
Double_t xmax_dxdy = 1.5;

Int_t nBins_y_dxdy = 500;
Double_t ymin_dxdy = -2.5;
Double_t ymax_dxdy = 2.5;

	h_E_eloss = new TH1D("E_eloss", Form("Scattered Electron Energy Loss in Target - SBS%i %i, %s", kine, sbsfieldscale, run_target.Data()), 500, 0.0, (0.1)*E_beam);
	h_E_ecorr_vs_vert = new TH2D("h_E_ecorr_vs_vert", Form("Corrected Beam Energy vs Vertex - SBS%i %i, %s; E_{e} (GeV); Z_{vertex} (m)", kine, sbsfieldscale, run_target.Data()), 250, -0.125, 0.125, 500, 0, 0.001);
	h_Q2 = new TH1D("h_Q2", Form("Momentum Transfer Q^2 - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0, 6.0);
	h_E_ep = new TH1D("h_E_ep", Form("Scattered Electron Energy - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0.0, 1.5*E_beam);
	h_E_pp = new TH1D("h_E_pp", Form("Scattered Proton Energy - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0.0, 1.5*E_beam);
	h_W = new TH1D("h_W", Form("Invariant Mass W - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0.0, 3.0);
	h_W_cut = new TH1D("h_W_cut", Form("Invariant Mass W (Coin & Vert Cuts) - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0.0, 3.0);
	h_W_fcut = new TH1D("h_W_fcut", Form("Invariant Mass W (Fiduc. Cuts) - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0.0, 3.0);
	h_W2 = new TH1D("h_W2", Form("Invariant Mass Squared W^{2} - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0.0, 3.0);
	h_W2recon = new TH1D("h_W2recon", Form("Invariant Mass Squared W^{2} Recon - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0.0, 3.0);
	h_KE_p = new TH1D("h_KE_p", Form("Scattered Proton Kinetic Energy - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0.0, 1.5*E_beam);
	h_tr_p = new TH1D("h_tr_p", Form("Scattered electron track momentum - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0.0, 5.0);
	h_tr_p_wcut = new TH1D("h_tr_p_wcut", Form("Scattered electron track momentum (with W2 cut)- SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0.0, 5.0);
	h_HCal_e = new TH1D("h_HCal_e", Form("HCal Clus. E (HCal_E > 0) - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0, 1.0);

	h_Ep = new TH1D("h_Ep", Form("E/p - SBS%i = %i%%, %s", kine, sbsfieldscale, run_target.Data()), 200, 0, 2);
	h_Ep_wcut = new TH1D("h_Ep_wcut", Form("E/p (wcut)- SBS%i = %i%%, %s", kine, sbsfieldscale, run_target.Data()), 200, 0, 2);	
	h_PS = new TH1D("h_PS", Form("Pre-Shower Clus. E - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0, 3);
	h_HCal_e = new TH1D("h_HCal_e", Form("HCal E (hcal_e)- SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 200, 0, 0.4);
	h_hcal_clus_e_sorted = new TH1D("h_hcal_clus_e_sorted", Form("HCal Clus E Sorted in Desc. Order of E (hcal_e_sorted)- SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 200, 0, 0.4);
	h_hcal_clus_e = new TH1D("h_hcal_clus_e", Form("HCal Clus. E (hcal_clus_e) - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 200, 0, 0.4);
	h2_hcal_e_V_clus_e = new TH2D("h2_hcal_e_V_clus_e", "h2_hcal_e_V_clus_e", 200, 0.0, 0.4, 200, 0.0, 0.4);
	h2_hcal_e_V_clus_e_sorted = new TH2D("h2_hcal_e_V_clus_e_sorted", "h2_hcal_e_V_clus_e_sorted", 200, 0.0, 0.4, 200, 0.0, 0.4);
	h2_hcal_clus_e_V_hcal_clus_e_sorted = new TH2D("h2_hcal_clus_e_V_hcal_clus_e_sorted", "h2_hcal_clus_e_V_hcal_clus_e_sorted", 200, 0.0, 0.4, 200, 0.0, 0.4);
	h_hcal_clus_e_indexed = new TH1D("h_hcal_clus_e_indexed", "h_hcal_clus_e_indexed", 10, 0, 10);
	h_hcal_clus_e_sorted_indexed = new TH1D("h_hcal_clus_e_sorted_indexed", "h_hcal_clus_e_sorted_indexed", 10, 0, 10);

	h_SHPS = new TH1D("h_SHPS", Form("SH + PS Clus. E - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0, 5);
	h_SHPS_wcut = new TH1D("h_SHPS_wcut", Form("SH + PS Clus. E (wcut)- SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0, 5);

	h_dx = new TH1D("h_dx",Form("dx (NO CUTS) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), nBins_y_dxdy, ymin_dxdy, ymax_dxdy);
	h_dx_cut = new TH1D("h_dx_cut",Form("dx (Basic CUTS) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), nBins_y_dxdy, ymin_dxdy, ymax_dxdy);
	h_dx_wcut = new TH1D("h_dx_wcut",Form("dx (W cut) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), nBins_y_dxdy, ymin_dxdy, ymax_dxdy);
	h_dx_fcut = new TH1D("h_dx_fcut",Form("dx (f cut) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), nBins_y_dxdy, ymin_dxdy, ymax_dxdy);
	h_dx_wcut_fcut = new TH1D("h_dx_wcut_fcut",Form("dx (W & Fiduc. Cuts) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), nBins_y_dxdy, ymin_dxdy, ymax_dxdy);
	h_dy = new TH1D("h_dy",Form("dy (NO CUTS) - SBS%i %i, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -1.25,1.25);
	h_dy_cut = new TH1D("h_dy_cut",Form("dy (Basic Cuts) - SBS%i %i, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -1.25,1.25);  
	h_dy_wcut = new TH1D("h_dy_wcut",Form("dy (W Cuts) - SBS%i %i, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -1.25,1.25);
	h_dy_wcut_fcut = new TH1D("h_dy_wcut_fcut",Form("dy (W & Fiduc. Cuts) - SBS%i %i, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -1.25,1.25);    

	h_dxdy = new TH2D("h_dxdy", Form("Hadron Spot(s) on HCal (NO CUTS) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), nBins_x_dxdy, xmin_dxdy, xmax_dxdy, nBins_y_dxdy, ymin_dxdy, ymax_dxdy );
	h_dxdy_wcut = new TH2D("h_dxdy_wcut", Form("Hadron Spot(s) on HCal (W cut) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), nBins_x_dxdy, xmin_dxdy, xmax_dxdy, nBins_y_dxdy, ymin_dxdy, ymax_dxdy );
	h_dxdy_wcut_fcut = new TH2D("h_dxdy_wcut_fcut", Form("Hadron Spot(s) on HCal (W & Fiduc. Cuts) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), nBins_x_dxdy, xmin_dxdy, xmax_dxdy, nBins_y_dxdy, ymin_dxdy, ymax_dxdy );
	h_dxdy_wcut_2multfcut = new TH2D("h_dxdy_wcut_2multfcut", Form("Hadron Spot(s) on HCal (W & 2*Mult_Fiduc. Cuts) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), nBins_x_dxdy, xmin_dxdy, xmax_dxdy, nBins_y_dxdy, ymin_dxdy, ymax_dxdy );
	h_dxdy_wcut_15multfcut = new TH2D("h_dxdy_wcut_15multfcut", Form("Hadron Spot(s) on HCal (W & 1.5*Mult_Fiduc. Cuts) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), nBins_x_dxdy, xmin_dxdy, xmax_dxdy, nBins_y_dxdy, ymin_dxdy, ymax_dxdy );
	h_dxdy_wcut_fcut_ADCtiming = new TH2D("h_dxdy_wcut_fcut_ADCtiming", Form("Hadron Spot(s) on HCal (W, Fiduc. & ADC timing Cuts) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), nBins_x_dxdy, xmin_dxdy, xmax_dxdy, nBins_y_dxdy, ymin_dxdy, ymax_dxdy );
	h_dxdy_wcut_fcut_AntiADCtiming = new TH2D("h_dxdy_wcut_fcut_AntiADCtiming", Form("Hadron Spot(s) on HCal (W, Fiduc. & Anti-ADC timing Cuts) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), nBins_x_dxdy, xmin_dxdy, xmax_dxdy, nBins_y_dxdy, ymin_dxdy, ymax_dxdy );

	h_dxdy_cut = new TH2D("h_dxdy_cut", Form("Hadron Spot(s) on HCal (Basic cuts) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), nBins_x_dxdy, xmin_dxdy, xmax_dxdy, nBins_y_dxdy, ymin_dxdy, ymax_dxdy );

	h_dxdy_ncut = new TH2D("h_dxdy_ncut", Form("Hadron Spot(s) on HCal (n cut) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 250, -2.5, 2.5 );
	h_dxdy_pcut = new TH2D("h_dxdy_pcut", Form("Hadron Spot(s) on HCal (p cut) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 250, -2.5, 2.5 );
	h_dxdy_fcut = new TH2D("h_dxdy_fcut", Form("Hadron Spot(s) on HCal (f cut) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 250, -2.5, 2.5 );
	h_xy = new TH2D("h_xy",Form("HCal Hadron Spots (x, y) (NO CUTS) - SBS%i %i, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_cut = new TH2D("h_xy_cut", Form("HCal Hadron Spots (x, y) (BASIC CUTS) - SBS%i %i, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_fcut = new TH2D("h_xy_fcut", Form("HCal Hadron Spots (x, y) (Fiduc. CUTS) - SBS%i %i, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_cut_p = new TH2D("h_xy_cut_p", Form("HCal Hadron Spots (x, y) (p CUT) - SBS%i %i, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_cut_n = new TH2D("h_xy_cut_n", Form("HCal Hadron Spots (x, y) (n CUT) - SBS%i %i, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);

	h_hcal_clusblk_ADC_time = new TH1D("h_hcal_clusblk_ADC_time", Form("ADC time of the highest energy block in the largest cluster - SBS%i %i, %s; ADC Time (ns)", kine, sbsfieldscale, run_target.Data()), 300, -100, 200);
	h_hcal_clusblk_ADC_time_cut = new TH1D("h_hcal_clusblk_ADC_time_cut", Form("ADC time of the highest energy block in the largest cluster (Cut window)- SBS%i %i, %s; ADC Time (ns)", kine, sbsfieldscale, run_target.Data()), 300, -100, 200);

	h_hcal_clusblk_ADC_time_diff = new TH1D("h_hcal_clusblk_ADC_time_diff", Form("HCal ADC Time - BBCal Shower ADC Time - SBS%i %i, %s; ADC Time (ns)", kine, sbsfieldscale, run_target.Data()), 100, 0, 100);
	h_hcal_clusblk_ADC_time_diff_cut = new TH1D("h_hcal_clusblk_ADC_time_diff_cut", Form("HCal ADC Time - BBCal Shower ADC Time (cut window)- SBS%i %i, %s; ADC Time (ns)", kine, sbsfieldscale, run_target.Data()), 100, 0, 100);

	h_bb_sh_atimeblk = new TH1D("h_bb_sh_atimeblk", Form("ADC Shower Timing for BB Block - SBS%i %i, %s; ADC Time (ns)", kine, sbsfieldscale, run_target.Data()), 300, -100, 200);
	
	h_theta_pq_n = new TH1D("h_theta_pq_n", Form("Theta pq for neutron - SBS%i %i, %s; theta_pq_n (rad);", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 0.60);
	h_theta_pq_p = new TH1D("h_theta_pq_p", Form("Theta pq for proton - SBS%i %i, %s; theta_pq_p (rad);", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 0.60);

	h_theta_pq_n_cut = new TH1D("h_theta_pq_n_cut", Form("dx for theta_pq neutron (events cut by theta_pq_cut) - SBS%i %i, %s; theta_pq_n (rad);", kine, sbsfieldscale, run_target.Data()), nBins_y_dxdy, ymin_dxdy, ymax_dxdy);
	h_theta_pq_n_anticut = new TH1D("h_theta_pq_n_anticut", Form("dx for theta_pq neutron (events cut by theta_pq_anticut) - SBS%i %i, %s; theta_pq_n (rad);", kine, sbsfieldscale, run_target.Data()), nBins_y_dxdy, ymin_dxdy, ymax_dxdy);
	h_theta_pq_p_cut = new TH1D("h_theta_pq_p_cut", Form("dx for theta_pq proton (events cut by theta_pq_cut) - SBS%i %i, %s; theta_pq_p (rad);", kine, sbsfieldscale, run_target.Data()), nBins_y_dxdy, ymin_dxdy, ymax_dxdy);
	h_theta_pq_p_anticut = new TH1D("h_theta_pq_p_anticut", Form("dx for theta_pq proton (events cut by theta_pq_anticut) - SBS%i %i, %s; theta_pq_p (rad);", kine, sbsfieldscale, run_target.Data()), nBins_y_dxdy, ymin_dxdy, ymax_dxdy);

	h_PAngleCorr_theta = new TH2D( "h_PAngCorr_theta",Form("BB theta vs HCal theta - SBS%i %i, %s", kine, sbsfieldscale, run_target.Data()), 200, 0.55, 0.75, 300, 0.35, 0.65 );
	h_PAngleCorr_phi = new TH2D( "h_PAngCorr_phi",Form("BB phi vs HCal phi - SBS%i %i, %s", kine, sbsfieldscale, run_target.Data()), 500, -0.4, 0.1, 500, 2.7, 3.2 );
	h_vz_cut = new TH1D("h_vz_cut",Form("BB phi vs HCal phi - SBS%i %i, %s; vertex z (m);", kine, sbsfieldscale, run_target.Data()), 250,-0.125,0.125);

	h_Nevents = new TH1D("h_Nevents", "Number of passing events", 1, 0, 1);
	h_eclus_sel = new TH1D("h_eclus_sel", "Index of hcal cluster selected", 10, 0, 10);
	cout << "Pulling experimental, fit, and other variables..." << endl;

	if( kine == 4 ){
		hcal_height = -0.312479; // Height of HCal above beamline
		//offset: -1.37521e-01
	}
	if( kine == 8 ){
		//From MC, dx_n for Mag70 and LD2 is: -2.39171e-02
		//From offset with -0.450 on data the dx_n is: -7.22499e-02
		//So, new offset is: -0.450 + (7.22499e-02 - 2.39171e-02) = -0.40166720 --> off by 0.0033489 --> -0.4016672+0.003349 = -0.39831820
		// hcal_height = -0.39831820; // Height of HCal above beamline --> dx_n = -0.0170160
		// hcal_height = -0.40010038000; //16930;
		//more negative puts the dx_n more negative
		hcal_height = -0.37064142; //-0.37733888; //-0.38368348; //16930; -0.3954, -0.36061770, -0.36765899
	}
	if( kine == 9 ){
		hcal_height = -0.36191616; //-0.36307351; // -0.36249082 //-0.35863270 //more negative puts the dx_n more negative
	}


	BB_dist = lookup_BB_dist_by_kine(kine);
	BB_theta = lookup_BB_angle_by_kine(kine, "rad");
	HCal_dist = lookup_HCal_dist_by_kine(kine);
	HCal_theta = lookup_HCal_angle_by_kine(kine, "rad");
	W_mean = lookup_parsed_cut(run_target, kine, sbsfieldscale, "W");
	W_sigma = lookup_parsed_cut(run_target, kine, sbsfieldscale, "W_sigma");

	ADC_time_min = lookup_ADC_time_cut(run_target, kine, sbsfieldscale, "ADC_time_min");
	ADC_time_max = lookup_ADC_time_cut(run_target, kine, sbsfieldscale, "ADC_time_max");
	ADC_time_mean = lookup_ADC_time_cut(run_target, kine, sbsfieldscale, "ADC_time_mean");

	ADC_diff_time_min = lookup_ADC_diff_time(run_target, kine, sbsfieldscale, "ADC_diff_time_min");
	ADC_diff_time_max = lookup_ADC_diff_time(run_target, kine, sbsfieldscale, "ADC_diff_time_max");

	dx_p = lookup_dxdy_by_kine_and_mag(run_target, kine, sbsfieldscale, "dx_p");
	if( dx_p == -1 ){ std::cout << "Lookup for dx_p may not exist. Returned value of -1. " << std::endl;}

	dx_p_sigma = (dx_p_scale)*lookup_dxdy_by_kine_and_mag(run_target, kine, sbsfieldscale, "dx_p_sigma");
	if( dx_p_sigma == -1 ){ std::cout << "Lookup for dx_p_sigma may not exist. Returned value of -1. " << std::endl;}

	dy_p = lookup_dxdy_by_kine_and_mag(run_target, kine, sbsfieldscale, "dy");
	if( dy_p == -1 ){ std::cout << "Lookup for dy_p may not exist. Returned value of -1. " << std::endl;}

	dy_p_sigma = (dy_scale)*lookup_dxdy_by_kine_and_mag(run_target, kine, sbsfieldscale, "dy_sigma");
	if( dy_p_sigma == -1 ){ std::cout << "Lookup for dy_p_sigma may not exist. Returned value of -1. " << std::endl;}

	if( run_target == "LD2" ){
		dx_n = lookup_dxdy_by_kine_and_mag(run_target, kine, sbsfieldscale, "dx_n");
		if( dx_n == -1 ){ std::cout << "Lookup for dx_n may not exist. Returned value of -1. " << std::endl;}

		dx_n_sigma = (dx_n_scale)*lookup_dxdy_by_kine_and_mag(run_target, kine, sbsfieldscale, "dx_n_sigma");
		if( dx_n_sigma == -1 ){ std::cout << "Lookup for dx_n_sigma may not exist. Returned value of -1. " << std::endl;}

		dx_pn_mean_center = (dx_p + dx_n)/2.0;
	}
	if( run_target == "LH2" ){
		dx_n = dx_p;
		dx_n_sigma = dx_p_sigma;
	}

	dy_n = dy;
	dy_n_sigma = dy_p_sigma;

	if( run_target == "LH2" ){
		dx_pn_max = abs(dx_p);
	}
	if( run_target == "LD2" ){
		dx_pn_max = abs( dx_p - dx_n );
		// dx_pn_max = 2.0;		
	}


	if( (dx_p == -1) || (dx_p_sigma == -1) || (dy_p == -1) || (dy_p_sigma == -1) || (dx_n == -1) || (dx_n_sigma == -1) || (dy_n == -1) || (dy_n_sigma == -1) || (dx_pn_max == -1)){
		cout << endl << "Press any key to continue...." << endl;
		std::cin.ignore(1000000000, '\n');
	}
	cout << "Finished pulling variables. " << endl << endl;

	cout << "Run parameters: " << endl;
	cout << "Kinematic: SBS" << kine << endl;
	cout << "Target: " << run_target.Data() << endl;
	cout << "Beam Energy: " << E_beam << endl;
	cout << "SBS Field: " << sbsfieldscale << "%" << endl;
	cout << "-----------------------------------" << endl;
	cout << "BB angle [deg]: " << (180/pi)*BB_theta << endl;
	cout << "SBS angle: " << lookup_SBS_angle_by_kine(kine, "deg") << endl;
	cout << "HCal angle [deg]: " << (180/pi)*HCal_theta << endl;
	cout << "HCal distance: " << HCal_dist << endl;
	cout << "dx_p = : " << dx_p << "; dx_p_sigma = " << dx_p_sigma << endl;
	cout << "dx_n = : " << dx_n << "; dx_n_sigma = " << dx_n_sigma << endl;
	cout << "dx_pn_max = " << dx_pn_max << " abs(dx_p - dx_n): " << abs(dx_p - dx_n) << endl;
	cout << "dy_p = " << dy_p << "; dy_p_sigma = " << dy_p_sigma << endl;
	cout << "-----------------------------------" << endl << endl;
	cout << "Lookups: " << endl;
	cout << "W_mean: " << W_mean << endl;
	cout << "W_sigma: " << W_sigma << endl;
	cout << "-----------------------------------" << endl << endl;


	if( kine == 11 ){
		rootfile_dir = "/volatile/halla/sbs/adr/Rootfiles/gmn_parsed/SBS11/pass1/";
	}
	if( kine == 4 ){
		rootfile_dir = Form("/volatile/halla/sbs/adr/Rootfiles/gmn_parsed/SBS%i/pass%i/", kine, pass);
	}
	if( kine == 8 ){
 		if( multi_run ){
 			rootfile_dir = Form("/work/halla/sbs/sbs-gmn/pass%i/SBS%i/%s/rootfiles/", 1, kine, run_target.Data());
 		}
 		if( !multi_run ){
  			rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/analysis_rootfiles/jboyd_parsed";			
 		}		
 		 if( use_parsed ){
 			rootfile_dir = "/volatile/halla/sbs/jboyd/analysis_rootfiles/jboyd_parsed";
 		}
	}
	if( kine == 9 ){
		 if( multi_run ){
 			rootfile_dir = Form("/work/halla/sbs/sbs-gmn/pass%i/SBS%i/%s/rootfiles/", 1, kine, run_target.Data());
 		}
 		if( !multi_run ){
  			rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/analysis_rootfiles/jboyd_parsed";			
 		}
 		if( use_parsed ){
 			rootfile_dir = "/volatile/halla/sbs/jboyd/analysis_rootfiles/jboyd_parsed";
 		}
	}
	else{
		// rootfile_dir = Form("/volatile/halla/sbs/seeds/parse/sbs%i_%s_mag%i/", kine, run_target.Data(), sbsfieldscale);
 		// rootfile_dir = Form("/volatile/halla/sbs/adr/Rootfiles/gmn_parsed/SBS%i/pass%i/", kine, pass);
 		if( multi_run ){
 			rootfile_dir = Form("/work/halla/sbs/sbs-gmn/pass%i/SBS%i/%s/rootfiles/", 1, kine, run_target.Data());
 		}
 		if( !multi_run ){
  			rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/analysis_rootfiles/jboyd_parsed";			
 		}

 	}



	if( single_run ){

		cout << "--------------------------------------" << endl;
		cout << "Adding files to TChain from: " << rootfile_dir.Data() << endl;
			if( kine == 11 ){
				TC->Add(Form("%s/gmn_parsed_SBS%i_targ%s_sbsmagscale%i.root", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale) );
			}
			if( kine == 4 ){
				TC->Add(Form("%s/gmn_parsed_SBS%i_targ%s_sbsmagscale%i.root", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale) );
			}
			else{
				// TC->Add(Form("%s/gmn_parsed_fulltree_SBS%i_%s_mag%i*",rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale));
				// TC->Add(Form("%s/gmn_parsed_SBS%i_targ%s_sbsmagscale%i.root", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale) );
				cout << "Adding: " << Form("%s/gmn_parsed_%s_SBS%i_mag%i.root", rootfile_dir.Data(), run_target.Data(), kine, sbsfieldscale) << endl;
				TC->Add(Form("%s/gmn_parsed_%s_SBS%i_mag%i.root", rootfile_dir.Data(), run_target.Data(), kine, sbsfieldscale));
			}
		// TC->Add(Form("%s/gmn_parsed_fulltree_SBS%i_%s_mag%i.root", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale));
		// TC->Add(Form("%s/*%i*.root", rootfile_dir.Data(), runnum));
	}
	cout << "--------------------------------------" << endl;
	if( multi_run && !use_parsed ){
		runnum = runnum_vec[0];
		cout << "Running in multi-run mode for runs: "  << endl;
		for(size_t run = 0; run < runnum_vec.size(); run++){
			cout << runnum_vec[run] << " ";
		}
		cout << endl;
		cout << "--------------------------------------" << endl;
		cout << "Adding files to TChain from: " << rootfile_dir.Data() << endl;
		// for(size_t run = 0; run < runnum_vec.size(); run++){
		// 	TC->Add(Form("%s/*%i*.root", rootfile_dir.Data(), runnum_vec[run]));
		// }
		TC->Add(Form("%s/*.root", rootfile_dir.Data()));
	}

	if( use_parsed ){
		rootfile_dir = "/volatile/halla/sbs/jboyd/analysis_rootfiles/jboyd_parsed";
		cout << "--------------------------------------" << endl;
		cout << "            USING PARSED FILES " << endl;
		cout << "--------------------------------------" << endl;
		cout << "Rootfile dir: " << rootfile_dir.Data() << endl;

		// cout << "Adding: " << Form("%s/gmn_parsed_%s_SBS%i_mag%i.root", rootfile_dir.Data(), run_target.Data(), kine, sbsfieldscale ) << endl;
		// TC->Add(Form("%s/gmn_parsed_%s_SBS%i_mag%i.root", rootfile_dir.Data(), run_target.Data(), kine, sbsfieldscale ) );
		
		// cout << "Adding: /lustre19/expphy/volatile/halla/sbs/adr/Rootfiles/gmn_parsed/SBS8/pass1/gmn_parsed_SBS8_targLD2_sbsmagscale70.root" << endl;
		// TC->Add("/lustre19/expphy/volatile/halla/sbs/adr/Rootfiles/gmn_parsed/SBS8/pass1/gmn_parsed_SBS8_targLD2_sbsmagscale70.root");

		TC->Add(Form("/lustre19/expphy/volatile/halla/sbs/jboyd/analysis_rootfiles/jboyd_parsed/SBS%i/%s/mag%i/*.root", kine, run_target.Data(), sbsfieldscale) );

		// cout << "Adding: /volatile/halla/sbs/jboyd/analysis_rootfiles/jboyd_parsed/gmn_parsed_LD2_SBS8_mag70_05_10_2023.root" << endl;
		// TC->Add("/volatile/halla/sbs/jboyd/analysis_rootfiles/jboyd_parsed/gmn_parsed_LD2_SBS8_mag70_05_10_2023.root");
	}

	cout << "--------------------------------------" << endl;
	cout << "Setting up branches... ";
	// Declare root tree variables and set values to memory locations in root file
	// Switch them on
	TC->SetBranchStatus( "*", 0 );

	// HCal
	TC->SetBranchStatus( "sbs.hcal.x", 1 );
	TC->SetBranchStatus( "sbs.hcal.y", 1 );
	TC->SetBranchStatus( "sbs.hcal.e", 1 );
	TC->SetBranchStatus( "sbs.hcal.nclus", 1);

	// HClas Cluster Tree Variables
	TC->SetBranchStatus( "sbs.hcal.clus.e", 1);
	TC->SetBranchStatus( "sbs.hcal.clus.x", 1);
	TC->SetBranchStatus( "sbs.hcal.clus.y", 1);
	TC->SetBranchStatus( "sbs.hcal.clus.atime", 1);
	TC->SetBranchStatus( "sbs.hcal.clus.tdctime", 1);
	TC->SetBranchStatus( "sbs.hcal.clus.id", 1);
	TC->SetBranchStatus( "sbs.hcal.clus.nblk", 1);

	// BB track
	TC->SetBranchStatus( "bb.tr.chi2", 1 );
	TC->SetBranchStatus( "bb.tr.n", 1 );
	TC->SetBranchStatus( "bb.tr.px", 1 );
	TC->SetBranchStatus( "bb.tr.py", 1 );
	TC->SetBranchStatus( "bb.tr.pz", 1 );    
	TC->SetBranchStatus( "bb.tr.p", 1 );
	TC->SetBranchStatus( "bb.tr.vx", 1 );
	TC->SetBranchStatus( "bb.tr.vy", 1 );
	TC->SetBranchStatus( "bb.tr.vz", 1 );
	TC->SetBranchStatus( "bb.tr.r_x", 1 );
	TC->SetBranchStatus( "bb.tr.r_y", 1 );
	TC->SetBranchStatus( "bb.tr.r_th", 1 );
	TC->SetBranchStatus( "bb.tr.r_ph", 1 );
	TC->SetBranchStatus( "bb.tr.tg_x", 1 );
	TC->SetBranchStatus( "bb.tr.tg_y", 1 );
	TC->SetBranchStatus( "bb.tr.tg_th", 1 );
	TC->SetBranchStatus( "bb.tr.tg_ph", 1 );
	TC->SetBranchStatus( "bb.gem.track.nhits", 1);
	TC->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1);
	TC->SetBranchStatus( "bb.sh.atimeblk", 1);

	// BBCal shower preshower
	TC->SetBranchStatus( "bb.ps.e", 1 );
	TC->SetBranchStatus( "bb.ps.x", 1 );
	TC->SetBranchStatus( "bb.ps.y", 1 );
	TC->SetBranchStatus( "bb.sh.e", 1 );
	TC->SetBranchStatus( "bb.sh.x", 1 );
	TC->SetBranchStatus( "bb.sh.y", 1 );
	TC->SetBranchStatus( "bb.sh.nclus", 1 );
	TC->SetBranchStatus( "bb.ps.nclus", 1 );
	TC->SetBranchStatus( "e.kine.W2", 1);

	// Trigger TDC
	TC->SetBranchStatus( "bb.tdctrig.tdc", 1 );
	TC->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
	TC->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );

// Set BRANCH ADDRESSES
	// HCal
	TC->SetBranchAddress( "sbs.hcal.x", &hcal_x );
	TC->SetBranchAddress( "sbs.hcal.y", &hcal_y );
	TC->SetBranchAddress( "sbs.hcal.e", &hcal_e );
	TC->SetBranchAddress( "sbs.hcal.nclus", &nclus );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.atime", &hcal_clusblk_ADC_time );
	TC->SetBranchAddress( "Ndata.sbs.hcal.clus.id", &Nhcal_clus_id );

	// HCal Cluster Tree Variables
	TC->SetBranchAddress( "sbs.hcal.clus.e", &hcal_clus_e);
	TC->SetBranchAddress( "sbs.hcal.clus.x", &hcal_clus_x);
	TC->SetBranchAddress( "sbs.hcal.clus.y", &hcal_clus_y);
	TC->SetBranchAddress( "sbs.hcal.clus.atime", &hcal_clus_atime);
	TC->SetBranchAddress( "sbs.hcal.clus.tdctime", &hcal_clus_tdctime);
	TC->SetBranchAddress( "sbs.hcal.clus.id", &hcal_clus_id);
	TC->SetBranchAddress( "sbs.hcal.clus.nblk", &hcal_clus_nblk);

	// BB track
	TC->SetBranchAddress( "bb.tr.chi2", bb_tr_chi2 );
	TC->SetBranchAddress( "bb.tr.n", &bb_tr_n );
	TC->SetBranchAddress( "bb.tr.px", bb_tr_px );
	TC->SetBranchAddress( "bb.tr.py", bb_tr_py );
	TC->SetBranchAddress( "bb.tr.pz", bb_tr_pz );
	TC->SetBranchAddress( "bb.tr.p", bb_tr_p );
	TC->SetBranchAddress( "bb.tr.vx", bb_tr_vx );
	TC->SetBranchAddress( "bb.tr.vy", bb_tr_vy );
	TC->SetBranchAddress( "bb.tr.vz", bb_tr_vz );
	TC->SetBranchAddress( "bb.tr.r_x", bb_fp_x );
	TC->SetBranchAddress( "bb.tr.r_y", bb_fp_y );
	TC->SetBranchAddress( "bb.tr.r_th", bb_fp_th );
	TC->SetBranchAddress( "bb.tr.r_ph", bb_fp_ph );
	TC->SetBranchAddress( "bb.tr.tg_x", bb_tgt_x );
	TC->SetBranchAddress( "bb.tr.tg_y", bb_tgt_y );
	TC->SetBranchAddress( "bb.tr.tg_th", bb_tgt_th );
	TC->SetBranchAddress( "bb.tr.tg_ph", bb_tgt_ph );

	// BBCal shower preshower
	TC->SetBranchAddress( "bb.ps.e", &bb_ps_e );
	TC->SetBranchAddress( "bb.ps.x", &bb_ps_x );
	TC->SetBranchAddress( "bb.ps.y", &bb_ps_y );
	TC->SetBranchAddress( "bb.sh.e", &bb_sh_e );
	TC->SetBranchAddress( "bb.sh.x", &bb_sh_x );
	TC->SetBranchAddress( "bb.sh.y", &bb_sh_y );
	TC->SetBranchAddress( "bb.sh.nclus", &SH_nclus );
	TC->SetBranchAddress( "bb.ps.nclus", &PS_nclus );
	TC->SetBranchAddress( "bb.sh.atimeblk", &bb_sh_atimeblk );

	// Trigger TDC
	TC->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
	TC->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
	TC->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );
	cout << " done. " << endl;
	cout << "--------------------------------------" << endl;

	p_Beam = E_beam/(1.0 + E_beam/Mp*(1.0 - cos(BB_theta)));

	E_loss_outgoing = cell_diameter/2.0/sin(BB_theta)*rho_tgt*dEdx_tgt; //Should be about 1 MeV
	if( useAlshield !=0 ) E_loss_outgoing += Alshieldthick*rho_Al*dEdx_Al;

//FIDUCIAL CUT:

//Now getting fiducial cut values from beam_variables.h
	// hcal_y_fmin = -0.75;
	// hcal_y_fmax = 0.75;
	// hcal_x_fmin = -2.015;
	// hcal_x_fmax = 1.285;

	if( kine == 8 ){
		Ep_sig_mult = 2.0;
		SH_PS_sig_mult = 2.0;		
	}
	if( kine == 9 ){
		Ep_sig_mult = 2.0;
		SH_PS_sig_mult = 2.0;
	}

	Ep_min = lookup_parsed_cut(run_target, kine, sbsfieldscale, "Ep") - Ep_sig_mult*lookup_parsed_cut(run_target, kine, sbsfieldscale, "Ep_sigma");
	Ep_max = lookup_parsed_cut(run_target, kine, sbsfieldscale, "Ep") + Ep_sig_mult*lookup_parsed_cut(run_target, kine, sbsfieldscale, "Ep_sigma");
	SH_PS_min = lookup_parsed_cut(run_target, kine, sbsfieldscale, "SH_PS_mean") - SH_PS_sig_mult*lookup_parsed_cut(run_target, kine, sbsfieldscale, "SH_PS_sigma")/2.0;

	if( kine == 8 ){
		if( !use_parsed ){
			master_cut_vec = {
				"sbs.hcal.nclus>0",
				"bb.ps.nclus>0",
				"bb.sh.nclus>0",
				"abs(bb.tr.vz[0])<0.075",
				"bb.gem.track.nhits[0]>2",
				"bb.tr.n==1",
				"bb.ps.e>0.15",
				
				// "sbs.hcal.nclus>0",
				// "bb.ps.nclus>0",
				// "bb.sh.nclus>0",
				// "abs(bb.tr.vz[0])<0.075",
				// "bb.gem.track.nhits[0]>3",
				// "bb.tr.n==1",
				// "bb.ps.e>0.2",
				// // Form("bb.tr.p[0]>%f", 1.10*lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "SH_PS_mean") ),
				// "sbs.hcal.e>0.005",
				// // // "((abs(((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))))-0.98)<0.15",
				// // "bb.ps.e+bb.sh.e>2.5"; -->bb.tr.p[0]
				// // // Form("bb.ps.e>%f", lookup_parsed_cut(run_target, kine, sbsfieldscale, "PS_min")),
				// // Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))>(%f)", 0.8),
				// // Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))<(%f)", 1.25),
				// // // // // Form("((abs(((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))))-%f)<%f", lookup_parsed_cut(runnum, "Ep"), lookup_parsed_cut(runnum, "Ep_sigma")),
				// // Form("sbs.hcal.e>%f",lookup_parsed_cut(run_target, kine, sbsfieldscale, "HCal_clus_e_cut")),
				// // Form("(bb.sh.e+bb.ps.e)>%f", SH_PS_min),
				// // "bb.sh.e>2.6",
				// // "bb.tr.p[0]>3.0",
				// Form("sbs.hcal.clus_blk.atime[0]>%f", 0.75*ADC_time_min),
				// Form("sbs.hcal.clus_blk.atime[0]<%f", 1.25*ADC_time_max),
				// // "bb.sh.atimeblk>-12",
				// // "bb.sh.atimeblk<12", 		
				// Form("(bb.sh.e+bb.ps.e)>%f", SH_PS_min),
				// Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))>(%f)&&((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))<(%f)", Ep_min, Ep_max),
				};
			}
		if( use_parsed ){
			master_cut_vec = {
				"sbs.hcal.nclus>0",
				"bb.ps.nclus>0",
				"bb.sh.nclus>0",
				"abs(bb.tr.vz[0])<0.075",
				"bb.gem.track.nhits[0]>2",
				"bb.tr.n==1",
				"bb.ps.e>0.15",
				};
			}
	}
	if( kine == 9 ){
		master_cut_vec = {
				"sbs.hcal.nclus>0",
				"bb.ps.nclus>0",
				"bb.sh.nclus>0",
				"abs(bb.tr.vz[0])<0.075",
				"bb.gem.track.nhits[0]>2",
				"bb.tr.n==1",
				"bb.ps.e>0.15",
				// Form("bb.tr.p[0]>%f", 1.10*lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "SH_PS_mean") ),
				// "sbs.hcal.e>0.05",
				// // "((abs(((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))))-0.98)<0.15",
				// "bb.ps.e+bb.sh.e>2.5"; -->bb.tr.p[0]
				// // Form("bb.ps.e>%f", lookup_parsed_cut(run_target, kine, sbsfieldscale, "PS_min")),
				Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))>(%f)&&((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))<(%f)", Ep_min, Ep_max),
				// // // // Form("((abs(((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))))-%f)<%f", lookup_parsed_cut(runnum, "Ep"), lookup_parsed_cut(runnum, "Ep_sigma")),
				// Form("sbs.hcal.e>%f",lookup_parsed_cut(run_target, kine, sbsfieldscale, "HCal_clus_e_cut")),
				Form("(bb.sh.e+bb.ps.e)>%f", SH_PS_min)
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
	cout << "--------------------------------------" << endl << endl;
	cout << "Number of RAW entries: " << TC->GetEntries() << endl << endl;
	cout << "--------------------------------------" << endl;
	cout << "--------------------------------------" << endl;
	cout << "Applying Master Cut: " << endl;
	cout << master_cut << endl;

//-----EVENT LIST WITH CUTS METHOD------
	// TEventList *ev_list = new TEventList("ev_list", "Elastic Events List");
	// TC->Draw(">>ev_list", master_cut);
	// Long64_t Nevents = ev_list->GetN();
	// h_Nevents->SetBinContent(1, Nevents);

	// cout << "--------------------------------------" << endl;
	// cout << "Number of events to analyze: " << Nevents << endl;
	// cout << "--------------------------------------" << endl;
//-----EVENT LIST WITH CUTS METHOD------

//-----TTREEFORMULA Method
	master_cut_formula = new TTreeFormula( "master_cut_formula", master_cut, TC );
	Long64_t Nevents = TC->GetEntries();
	Int_t treenum = 0, currenttreenum = 0;

//-----TTREEFORMULA Method

	cout << "--------------------------------------" << endl;
	cout << "Starting analysis loop on events..... " << endl;

	elastic_yield = 0;
	dxdy_cnt = 0;
	dxdy_wcut_cnt = 0;

	int watch_cnt = 0;	
	int five_percent = int(0.05*Nevents);
	vector<double> time_for_five;
	double average_time = 0.0, time_remaining;
	StopWatch->Start();
	
//-----EVENT LIST WITH CUTS METHOD------
	// for(Long64_t nevent = 0; nevent < Nevents; nevent++){

	// 	TC->GetEntry( ev_list->GetEntry( nevent ));
//-----EVENT LIST WITH CUTS METHOD------

//-----TTREEFORMULA Method
	long nevent = 0;
	while( TC->GetEntry(nevent++)){

		currenttreenum = TC->GetTreeNumber();
		if( nevent == 1 || currenttreenum != treenum ){
			treenum = currenttreenum;
			master_cut_formula->UpdateFormulaLeaves();
		}

		bool failedMasterCut = master_cut_formula->EvalInstance(0) == 0;

		if( failedMasterCut ){ continue; }
//-----TTREEFORMULA Method

		if( nevent%five_percent == 0){		
			StopWatch->Stop();

			if( watch_cnt == 0){
				cout << "Evt: " << nevent <<"/" << Nevents << Form("(%.0f/100%%)", 100.0*double(1.0*nevent/Nevents)) << ". Elastic yield = " << elastic_yield << ". dxdy_cnt: " << dxdy_cnt << ", dxdy_wcut_cnt: " << dxdy_wcut_cnt << ". " << endl;
			}

			if( watch_cnt > 0 ){
				time_for_five.push_back(StopWatch->RealTime());	
				average_time = VectorMean(time_for_five);
				// cout << "average time for 5 = " << average_time << endl;
				time_remaining = average_time*( 1.0 - double(nevent)/double(Nevents));
				cout << "Evt: " << nevent <<"/" << Nevents << Form("(%.0f/100%%)", 100.0*double(1.0*nevent/Nevents)) << ". Elastic yield = " << elastic_yield << ". dxdy_cnt: " << dxdy_cnt << ", dxdy_wcut_cnt: " << dxdy_wcut_cnt << ". Time left: " << time_remaining << endl;
			}
			StopWatch->Reset();
			StopWatch->Continue();
			watch_cnt++;
		}

		if( sort_hcal_cluster_energy ){
		//SORT THE HCAL CLUSTER ENERGY ARRAY BY DESCENDING ORDER
			for(int clus = 0; clus < max_clus; clus++){
				hcal_clus_e_sorted[clus].ArrayValue = hcal_clus_e[clus];
				hcal_clus_e_sorted[clus].index = clus;
			}

			sort( hcal_clus_e_sorted, hcal_clus_e_sorted + max_clus, CompareArrayValueWithIndex );

			hcal_x = hcal_clus_x[hcal_clus_e_sorted[0].index];
			hcal_y = hcal_clus_y[hcal_clus_e_sorted[0].index];			
		}

		h_HCal_e->Fill(hcal_e);
		h2_hcal_e_V_clus_e->Fill(hcal_clus_e[0], hcal_e);
		if( sort_hcal_cluster_energy ){
			h_hcal_clus_e_sorted->Fill(hcal_clus_e_sorted[0].ArrayValue);
			h2_hcal_e_V_clus_e_sorted->Fill(hcal_clus_e_sorted[0].ArrayValue, hcal_e);
			h2_hcal_clus_e_V_hcal_clus_e_sorted->Fill(hcal_clus_e[0], hcal_clus_e_sorted[0].ArrayValue);
		}
		h_hcal_clus_e->Fill(hcal_clus_e[0]);

		for( int index = 0; index < 10; index++ ){
			h_hcal_clus_e_indexed->SetBinContent(index, h_hcal_clus_e_indexed->GetBinContent(index) + hcal_clus_e[index] );
			h_hcal_clus_e_sorted_indexed->SetBinContent(index, h_hcal_clus_e_sorted_indexed->GetBinContent(index) + hcal_clus_e_sorted[index].ArrayValue);
		}

	//Timing stuff
		h_hcal_clusblk_ADC_time->Fill(hcal_clusblk_ADC_time[0]);
		
		hcal_clus_ADCtime_diff = hcal_clusblk_ADC_time[0] - bb_sh_atimeblk;
		h_bb_sh_atimeblk->Fill(bb_sh_atimeblk);
		h_hcal_clusblk_ADC_time_diff->Fill(hcal_clus_ADCtime_diff);

		bool b_hcal_clusblk_ADC_cut = false;
		bool b_adc_diff_time_cut = false;

		if( apply_ADC_cut ){
		// if( !use_best_cluster ){
			if( (hcal_clusblk_ADC_time[0] < 0.9*ADC_time_min) || (hcal_clusblk_ADC_time[0] > 1.1*ADC_time_max) ){
				b_hcal_clusblk_ADC_cut = true;
				n_hcal_clusblk_atime_cut++;
				h_hcal_clusblk_ADC_time_cut->Fill(hcal_clusblk_ADC_time[0]);
				continue;
			}			
		}


		if( apply_adc_diff_time_cut ){
			if( (hcal_clus_ADCtime_diff < 0.9*ADC_diff_time_min) || (hcal_clus_ADCtime_diff > 1.1*ADC_diff_time_max) ){
				b_adc_diff_time_cut = true;
				n_adc_diff_time_cnt++;
				h_hcal_clusblk_ADC_time_diff_cut->Fill(hcal_clus_ADCtime_diff);
				continue;
			}			
		}


		// double bbcal_time=0.0, hcal_time=0.0;

		// for(int ihit=0; ihit<TDCTndata; ihit++){
		// 	if(TDCT_id[ihit]==5){
		// 		bbcal_time=TDCT_tdc[ihit];
		// 	}
		// 	if(TDCT_id[ihit]==0){
		// 		hcal_time=TDCT_tdc[ihit];
		// 	}
		// }

		// double diff = hcal_time - bbcal_time; 

		// if( fabs(diff - tdiff)>tdiff_max ){
		// 	continue;
		// }

		//Sanity check
      	// if( (int)bb_tr_n!=1 ){
      	// 	cout << "**************************************************************" << endl;
      	// 	cout << "--------------------------------------------------------------" << endl;
      	// 	cout << endl << endl << "WARNING: Total tracks not as expected from global cut. Check globalcut for errors." << endl << endl;
      	// 	cout << "--------------------------------------------------------------" << endl;
      	// 	cout << "**************************************************************" << endl;
      	// } 


      	if( correct_beam_energy ){
      		Eloss = (bb_tr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      		h_E_eloss->Fill( Eloss );

      		E_beam_final = E_beam - Eloss;
      		h_E_ecorr_vs_vert->Fill( bb_tr_vz[0], E_beam_final);
      	}
      	if( !correct_beam_energy){
      		Eloss = (bb_tr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      		h_E_eloss->Fill( Eloss );

      		E_beam_final = E_beam;
      		h_E_ecorr_vs_vert->Fill( bb_tr_vz[0], E_beam_final);
      	}

      	////Corrections
	    //Correct the beam energy with energy loss in target using vertex position
	    Double_t E_loss = (bb_tr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
	    Double_t E_corr = E_beam - Eloss;	

      	p_corr = bb_tr_p[0] - E_loss_outgoing; //Neglecting mass of e'


    //Proceed only if one track exists in BB arm - lowest chi2 track always first element
      	// if( bb_tr_n > 1){
      	// 	continue;
      	// }
		
		p_Beam = E_beam/(1.0 + E_beam/Mp*(1.0 - cos(BB_theta)));
      	
      	e_prime_theta = acos( bb_tr_pz[0]/bb_tr_p[0] ); //Ucorrected track momenutm to reconstruct e' theta
      	e_prime_phi = atan2( bb_tr_py[0], bb_tr_px[0]);

      	TVector3 vertex( 0, 0, bb_tr_vz[0] ); // z location of vertex in hall coordinates
		TLorentzVector P_beam( 0, 0, E_beam_final, E_beam_final ); //Mass of e negligable
		TLorentzVector k_prime( bb_tr_px[0], bb_tr_py[0], bb_tr_pz[0], bb_tr_p[0] );
		TLorentzVector P_targ( 0, 0, 0, Mp );

		p_el = E_beam_final/( 1.0+E_beam_final/Mp*( 1.0-cos(e_prime_theta) ) );
		nu = E_beam_final - bb_tr_p[0];
		pp = sqrt( pow(nu,2)+2.*Mp*nu );
		nucleon_phi = e_prime_phi + pi; //assume coplanarity
		//double thetanucleon = acos( (E_corr - BBtr_p[0]*cos(etheta))/pp ); //use elastic constraint on nucleon kinematics
		nucleon_theta = acos( (E_beam_final - bb_tr_pz[0])/pp ); //use elastic constraint on nucleon kinematics

		TVector3 pNhat( sin(nucleon_theta)*cos(nucleon_phi), sin(nucleon_theta)*sin(nucleon_phi), cos(nucleon_theta) );

		//Define HCal coordinate system
		TVector3 HCAL_zaxis( sin(-HCal_theta ), 0, cos(-HCal_theta) );
		TVector3 HCAL_xaxis( 0, -1, 0 );
		TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();

		TVector3 HCAL_origin = HCal_dist*HCAL_zaxis + hcal_height*HCAL_xaxis;

		TVector3 HCAL_pos = HCAL_origin + (hcal_x*HCAL_xaxis) + (hcal_y*HCAL_yaxis);

		//Define intersection points for hadron vector
		scint_intersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) );
		TVector3 HCAL_intersect = vertex + scint_intersect * pNhat;

		//Define the expected position of hadron on HCal from BB track
		x_expected_HCal = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
		y_expected_HCal = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );

//--------------------------------------------------------------
//Calculate theta pq variables
      	//Reconstructed momentum, corrected for mean loss exiting the target
		p_recon = bb_tr_p[0] + E_loss_outgoing; 

		TLorentzVector k_prime_recon(p_recon*bb_tr_px[0]/bb_tr_p[0], p_recon*bb_tr_py[0]/bb_tr_p[0], p_recon*bb_tr_pz[0]/bb_tr_p[0], p_recon);
		TLorentzVector q_recon = P_beam - k_prime_recon;
		TVector3 qvec_recon = q_recon.Vect();

	//Calculate q vector as beam momentum - scattered k
		TLorentzVector q = P_beam - k_prime;

	//Expected neutron direction
		TVector3 Neutron_Direction = (HCAL_pos - vertex).Unit();

	//Expected proton direction
		//Need to incorporate deflection due to SBS magnet
		double Bdl = sbsfieldscale_Frac*sbsmaxfield*Dgap;//maxSBSfield
		double Proton_Deflection = tan( 0.3*Bdl/qvec_recon.Mag() )*(HCal_dist - (SBSdist + Dgap/2.0) ); 
		// double Proton_Deflection = tan( 0.3*Bdl/qvec_recon.Mag() )*(HCal_dist - (SBSdist + Dgap/2.0) ); 
		// double Proton_Deflection = dx_pn_max;

		TVector3 Proton_Direction = (HCAL_pos + Proton_Deflection*HCAL_xaxis - vertex).Unit();

		// theta_pq_n = acos( Neutron_Direction.Dot( q.Vect() ) );
		// theta_pq_p = acos( Proton_Direction.Dot( q.Vect() ) );

		theta_pq_n = acos(Neutron_Direction.Dot( qvec_recon.Unit() ));
		theta_pq_p = acos(Proton_Direction.Dot( qvec_recon.Unit() ));

		// h_theta_pq_n->Fill(theta_pq_n);
		// h_theta_pq_p->Fill(theta_pq_p);

		h_theta_pq_n->Fill( theta_pq_n );
		h_theta_pq_p->Fill( theta_pq_p );

	//---------------------------
		if( theta_pq_cut ){
			if( run_target == "LD2" ){
				if( theta_pq_n > theta_pq_n_thresh ){
					continue;
				}
				if( theta_pq_p > theta_pq_p_thresh){
					continue;
				}

			}
		}

//-----------------------------------------------------
// Search for best cluster on HCal (Not just highest energy cluster)
//-----------------------------------------------------
	
	if( use_best_cluster ){
		Int_t clus_sel_dx_atime = 0, clus_sel_dx_atimediff_simple = 0, clus_sel_dx_atimediff = 0, clus_sel_dx_atime_and_maxE = 0;
		Double_t clus_sel_atimediff = 1000.0;
		Double_t clus_sel_atime = 1000.0;

		Int_t clus_sel_dx_theta_pq = 0;
		Double_t clus_sel_theta_pq_diff = 1000.0;

		Int_t clus_sel_dx_dxdy = 0;
		Double_t clus_sel_dxdy_diff = 1000.0;

		Double_t clus_sel_theta_pq;

		//Add number of elastic clusters to: Nclusters
		Nclusters+=Nhcal_clus_id;
		// cout << "Nhcal_clus_id to search through: " << Nhcal_clus_id << endl;

	//START OF LOOP FOR BEST CLUSTER SELECTION:
		vector<Double_t> hcal_e_ADC_coin_passed;
		vector<Int_t> hcal_e_index_ADC_coin_passed;
		vector<Double_t> ADC_diff_time_coin_passed;
		int ADC_coin_pass_cnt = 0;

//------------------------------------
	//We first need to build the vector of the clusters which pass the coincidence window:
	/// ADC Coincidence timing check. Save the clusters that pass and then assess the clusters' energies
		//Here we save clusters that pass through the timing window
		for( Int_t clus_sel = 0; clus_sel < Nhcal_clus_id; clus_sel++){
		
		//Timing variables and minimization
			//Get time values from tree
			double hcal_clus_sel_ADCtime = hcal_clus_atime[clus_sel];
			double hcal_clus_sel_ADCtime_diff = hcal_clus_sel_ADCtime - bb_sh_atimeblk;

			bool passed_ADC_timing = false;
			passed_ADC_timing = ( hcal_clus_sel_ADCtime_diff >= ADC_diff_time_min ) && ( hcal_clus_sel_ADCtime_diff <= ADC_diff_time_max );
	
			if( passed_ADC_timing ){
				hcal_e_ADC_coin_passed.push_back(hcal_clus_e[clus_sel]);
				hcal_e_index_ADC_coin_passed.push_back(clus_sel);
				ADC_diff_time_coin_passed.push_back(hcal_clus_sel_ADCtime_diff);
				ADC_coin_pass_cnt++;
			}
		}

	//Now we have the vector of clusters passing the timing window: 
		//hcal_e_ADC_coin_passed && hcal_e_index_ADC_coin_passed & ADC_diff_time_coin_passed
		//Let us send this to our best cluster selection function:
		if( ADC_coin_pass_cnt > 0 ){
			if( !use_scoring ){
				bestClusterSelectionHcalE_CoinPassAndMaxE( hcal_e_ADC_coin_passed, hcal_e_index_ADC_coin_passed, best_cluster_OrigIndex_Energy_Index_Timing_Score );
			}			
			if( use_scoring ){
				bestClusterSelectionHcalE_CoinTimingMaxE_NormalizedScoring( hcal_e_ADC_coin_passed, hcal_e_index_ADC_coin_passed, ADC_diff_time_coin_passed, best_cluster_OrigIndex_Energy_Index_Timing_Score );
				// bestClusterSelectionHcalE_CoinTimingMaxE_EuclideanScoring( hcal_e_ADC_coin_passed, hcal_e_index_ADC_coin_passed, ADC_diff_time_coin_passed, best_cluster_OrigIndex_Energy_Index_Timing_Score );
			}
		}

 		if( hcal_cluster_minimize == "coin_time_and_maxE" ){
 			int best_energy_index = -1;
			best_energy_index = static_cast<Int_t>(best_cluster_OrigIndex_Energy_Index_Timing_Score[2]);
			dx_bestcluster = hcal_clus_x[best_energy_index] - x_expected_HCal;
			dy_bestcluster = hcal_clus_y[best_energy_index] - y_expected_HCal;
			HCal_ADC_time_bestcluster = hcal_clus_atime[best_energy_index];	

			// dx_bestcluster = hcal_clus_x[clus_sel_dx_atime_and_maxE] - x_expected_HCal;
			// dy_bestcluster = hcal_clus_y[clus_sel_dx_atime_and_maxE] - y_expected_HCal;
			// HCal_ADC_time_bestcluster = hcal_clus_atime[clus_sel_dx_atime_and_maxE];	
		} 
//SELECTION: Only looking at passing the coincidence timing window:
		// if( hcal_cluster_minimize == "coin_time" ){
		// 	dx_bestcluster = hcal_clus_x[clus_sel_dx_atime] - x_expected_HCal;
		// 	dy_bestcluster = hcal_clus_y[clus_sel_dx_atime] - y_expected_HCal;
		// 	HCal_ADC_time_bestcluster = hcal_clus_atime[clus_sel_dx_atime];				
		// }

	}
// --------------------------------------


	//--------------------------------------------------------------

		TLorentzVector P_gammaN = P_targ + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)
		TLorentzVector P_gammaN_second = P_targ + q_recon;
		double W_recon = 0.0;

		E_ep = sqrt( pow(Me,2) + pow(bb_tr_p[0],2) ); // Obtain the scattered electron energy
		h_E_ep->Fill( E_ep );

		p_ep = bb_tr_p[0];

		h_tr_p->Fill(bb_tr_p[0]);

		Q2 = 2*E_beam_final*E_ep*( 1-(bb_tr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
		h_Q2->Fill( Q2 );

		//Get invariant mass transfer W from the four-momentum of the scattered nucleon
		W = P_gammaN.M();
		W_recon = P_gammaN_second.M();
		W2 = pow(W, 2);
		h_W->Fill( W );
		h_W2->Fill(W2);
		h_W2recon->Fill( pow(W_recon, 2) );

		//Use the electron kinematics to predict the proton momedntum assuming elastic scattering on free proton at rest (will need to correct for fermi motion):
		E_pp = nu + Mp; // Get energy of the proton
		E_nucleon = sqrt(pow(pp,2)+pow(Mp,2)); // Check on E_pp, same
		h_E_pp->Fill( E_pp ); // Fill histogram

		KE_p = nu; // For elastics
		h_KE_p->Fill( KE_p );

		Double_t Ep = (bb_ps_e + bb_sh_e)/(bb_tr_p[0]);
		Double_t PS = bb_ps_e;
		Double_t SHPS = bb_ps_e + bb_sh_e;
		Double_t HCal_e = hcal_e;

		h_Ep->Fill(Ep);
		h_PS->Fill(PS);
		h_SHPS->Fill(SHPS);
		h_HCal_e->Fill(hcal_e);


		if( use_best_cluster ){
			dx = dx_bestcluster;
			dy = dy_bestcluster;
		}

		if( !use_best_cluster ){
			dx = hcal_x - x_expected_HCal;
			dy = hcal_y - y_expected_HCal;				
		}

	//If using the scoring method for best cluster we need to sort through the scores and matches and pick a best cluster
	//This is all based on the original priority_best_cluster.
	//I think we should work chronogically from score 0 to score 3

		//Resolve the hadron spots without cuts
		h_dx->Fill( dx );
		h_dy->Fill( dy );
		h_dxdy->Fill( dy, dx );
		h_xy->Fill( hcal_y, hcal_x );

		dxdy_cnt++;

	// Coincidence timing cut and vertex cut to resolve W well
		// if( fabs(diff - tdiff)<tdiff_max && fabs(bb_tr_vz[0])<=0.075 ){
		// 	h_W_cut->Fill( W );
		// } 

		double w2_mult;
		if( kine == 8 ){
			w2_mult = 1.5;
		}
		if( kine == 9 ){
			w2_mult = 2.0;
		}

		// Preliminary HCal projections with single cut on W
		if( fabs(W - W_mean) < w2_mult*W_sigma ){
		// if( ( W > (W_mean - 0.75*W_sigma) ) && ( W < (W_mean + 0.75*W_sigma) ) ){
		// if( W2 < 1.05 && W2 > 0.85){
			h_dx_wcut->Fill( dx );
			h_dy_wcut->Fill ( dy );
			h_dxdy_wcut->Fill( dy, dx );
			h_W_cut->Fill( W );
			h_tr_p_wcut->Fill( bb_tr_p[0] );
			h_SHPS_wcut->Fill( SHPS );
			h_Ep_wcut->Fill( Ep );

			if( true ){
				if( theta_pq_n < theta_pq_n_thresh){
					h_theta_pq_n_cut->Fill( dx );
				}
				if( theta_pq_n > theta_pq_n_thresh){
					h_theta_pq_n_anticut->Fill( dx );
				}

				if( theta_pq_p < theta_pq_p_thresh){
					h_theta_pq_p_cut->Fill( dx );
				}
				if( theta_pq_p > theta_pq_p_thresh){
					h_theta_pq_p_anticut->Fill( dx );
				}		
			}

			dxdy_wcut_cnt++;
		}

		//Populate position histograms with cuts
		h_dxdy_cut->Fill( dy, dx );
		h_dx_cut->Fill( dx );
		h_dy_cut->Fill( dy );

		//Populate BB/HCal correlation histograms from elastics
		h_PAngleCorr_phi->Fill( e_prime_phi, nucleon_phi );
		h_PAngleCorr_theta->Fill( e_prime_theta, nucleon_theta );

		//Fill vertex position histogram for cut on tracks
    	h_vz_cut->Fill( bb_tr_vz[0] );



		//Check "elastic" events on center HCal for id with spot checks
		bool HCal_on = false;
		bool is_p = false;
		bool is_n = false;

		double fcut_mult = 1;

		if( use_acceptance_cut ){
			//Acceptance cut
			if( hcal_y>hcal_y_fmin && hcal_y<hcal_y_fmax && hcal_x>hcal_x_fmin && hcal_x<hcal_x_fmax ){
				HCal_on = true;
			}			
		}
		if( !use_acceptance_cut ){
			HCal_on = true;
		}

	//FIDUCIAL Cut
		if( fiducial_cut ){
			fiducial_active_area_xmax = hcal_x_fmax - dx_pn_max - fcut_mult*dx_p_sigma;
			fiducial_active_area_xmin = hcal_x_fmin + dx_pn_max + fcut_mult*dx_p_sigma;
			fiducial_active_area_ymax = hcal_y_fmax - fcut_mult*dy_p_sigma;
			fiducial_active_area_ymin = hcal_y_fmin + fcut_mult*dy_p_sigma;
			// if( ((y_expected_HCal - fcut_mult*dy_p_sigma) > hcal_y_fmin) && ((y_expected_HCal + fcut_mult*dy_p_sigma) < hcal_y_fmax) && ((x_expected_HCal - dx_pn_max - fcut_mult*dx_p_sigma) > hcal_x_fmin) && ((x_expected_HCal + dx_pn_max + fcut_mult*dx_p_sigma) < hcal_x_fmax) ){
			// 	apply_fcut = true;
			// }
			if( (hcal_x < fiducial_active_area_xmax) && (hcal_x > fiducial_active_area_xmin) && (hcal_y < fiducial_active_area_ymax) && (hcal_y > fiducial_active_area_ymin) ){
				apply_fcut = true;
			}

		//For testing out the effect of fcut_mult.... let's have a section to always look at 2*fcut_mult and 3*fcut_mult
			if( ((y_expected_HCal - 2.0*fcut_mult*dy_p_sigma) > hcal_y_fmin) && ((y_expected_HCal + 2.0*fcut_mult*dy_p_sigma) < hcal_y_fmax) && ((x_expected_HCal - dx_pn_max - 2.0*fcut_mult*dx_p_sigma) > hcal_x_fmin) && ((x_expected_HCal + dx_pn_max + 2.0*fcut_mult*dx_p_sigma) < hcal_x_fmax) ){
				if( fabs(W - W_mean) < w2_mult*W_sigma ){
					h_dxdy_wcut_2multfcut->Fill( dy, dx );
				}
			}	
			if( ((y_expected_HCal - 1.5*fcut_mult*dy_p_sigma) > hcal_y_fmin) && ((y_expected_HCal + 1.5*fcut_mult*dy_p_sigma) < hcal_y_fmax) && ((x_expected_HCal - dx_pn_max - 1.5*fcut_mult*dx_p_sigma) > hcal_x_fmin) && ((x_expected_HCal + dx_pn_max + 1.5*fcut_mult*dx_p_sigma) < hcal_x_fmax) ){
				if( fabs(W - W_mean) < w2_mult*W_sigma ){
					h_dxdy_wcut_15multfcut->Fill( dy, dx );
				}
			}	
			if( ( pow( (hcal_x - x_expected_HCal - dx_p)/dx_p_sigma,2) + pow( (hcal_y - y_expected_HCal - dy_p)/dy_p_sigma,2) ) <= pow(1.5,2) ){
				is_p = true;
				is_n = false;
			}

			if( ( pow( (hcal_x - x_expected_HCal - dx_n)/dx_n_sigma,2) + pow( (hcal_y - y_expected_HCal - dy_n)/dy_n_sigma,2) ) <= pow(1.5,2) ){
				is_p = false;
				is_n = true;
			}

	//Fill respective histograms for these checks.
			if( HCal_on && is_n && apply_fcut ) h_dxdy_ncut->Fill( dy, dx );
			if( HCal_on && is_p && apply_fcut ) h_dxdy_pcut->Fill( dy, dx );

	//----------ALL HADRONS
			if( HCal_on && apply_fcut ){
				h_dxdy_fcut->Fill( dy, dx );
				h_dx_fcut->Fill( dx );
				h_W_fcut->Fill( W );
				h_xy_fcut->Fill( hcal_y, hcal_x );
				//including work cut
				if( fabs(W - W_mean) < w2_mult*W_sigma ){
				// if( fabs(W - 0.94) < 0.375 ){
				// if( fabs(W - W_mean) < 1.0*W_sigma ){
					h_dx_wcut_fcut->Fill( dx );
					h_dy_wcut_fcut->Fill ( dy );
					h_dxdy_wcut_fcut->Fill( dy, dx );
					if( (hcal_clusblk_ADC_time[0] > 0.75*ADC_time_min) || (hcal_clusblk_ADC_time[0] < 1.25*ADC_time_max) ){
						h_dxdy_wcut_fcut_ADCtiming->Fill( dy, dx );
					}
					if( (hcal_clusblk_ADC_time[0] < 0.75*ADC_time_min) || (hcal_clusblk_ADC_time[0] > 1.25*ADC_time_max) ){
						h_dxdy_wcut_fcut_AntiADCtiming->Fill( dy, dx );
					}
				}

				elastic_yield++;				
			}

	//----------neutron
			if( HCal_on && is_n && apply_fcut ){
				// if( (hcal_x - dx_pn_max )>hcal_x_fmin ){
					h_xy_cut_n->Fill( hcal_y, hcal_x );
				// }
			}
	//----------proton
			else if( HCal_on && is_p && apply_fcut ){
				// if( (hcal_x + dx_pn_max)<hcal_x_fmax ){
					h_xy_cut_p->Fill( hcal_y, hcal_x );
				// }
			}
		}
	//END OF FIDUCIAL Cut
		//Still should count elastic yields if we got this far.....
		if( !fiducial_cut ){
			elastic_yield++;
		}
//end of events loop  
    }
	
	bb_tr_p_min = (0.01)*h_tr_p->FindFirstBinAbove(100);
	bb_tr_p_wcut_min = (0.01)*h_tr_p_wcut->FindFirstBinAbove(100);
	bb_tr_p_cut = min( bb_tr_p_min, bb_tr_p_wcut_min );

	cout << "---------------------------------------" << endl;
	cout << "-----Finished going through events-----" << endl;
	cout << "---------------------------------------" << endl;

	outfile->Write();
	outfile->Close();

    if( calibrate ){
    	TFile *infile = new TFile(outfile->GetName(), "READ");

    	TString h_plot_select = "wcut";
    	cout << "Selected plot type: " << h_plot_select.Data() << endl << endl;
    	TH1D *hin_dx_select;
    	TH1D *hin_dy_select;
    	TH1D *hin_dxdy_select;
    	TH1D *hin_dx = static_cast<TH1D*>(infile->Get("h_dx"));

    	hin_dx_select = static_cast<TH1D*>(infile->Get(Form("h_dx_%s", h_plot_select.Data())));
    	hin_dy_select = static_cast<TH1D*>(infile->Get(Form("h_dy_%s", h_plot_select.Data())));
    	hin_dxdy_select = static_cast<TH1D*>(infile->Get(Form("h_dxdy_%s", h_plot_select.Data())));

    	TH1D *hin_hcal_clusblk_ADC_time = static_cast<TH1D*>(infile->Get("h_hcal_clusblk_ADC_time"));

    	// TCanvas *c_hcal_clusblk_ADC_time = new TCanvas("c_hcal_clusblk_ADC_time", "c_hcal_clusblk_ADC_time", 600, 500);
    	// h_hcal_clusblk_ADC_time->Draw();

    	TCanvas *c_dx = new TCanvas("c_dx", "c_dx", 600, 500);
    	hin_dx_select->Draw();
  	
  	//------ p -------
    	TF1 *fit_dx_p = new TF1("fit_dx_p", fit_gaus, -1.5, -0.2, 3);
  
    	fit_dx_p->SetParName(0, "dx_p Norm");
		fit_dx_p->SetParName(1, "dx_p Center");
		fit_dx_p->SetParName(2, "dx_p Sigma");
		fit_dx_p->SetLineColor(2);

		if( kine == 4 && sbsfieldscale == 30){
			fit_dx_p->SetParLimits(0, 0, hin_dx_select->GetMaximum());
			fit_dx_p->SetParLimits(1, -0.7, -0.5);
			fit_dx_p->SetParLimits(2, 0.1, 0.19);
		}

		if( kine == 4 && sbsfieldscale == 50){
			fit_dx_p->SetParLimits(0, 0, hin_dx_select->GetMaximum());
			fit_dx_p->SetParLimits(1, -1.1, -1.05);
			fit_dx_p->SetParLimits(2, 0.1, 0.14);
		}

		if( kine == 8 && sbsfieldscale == 70){
			fit_dx_p->SetParLimits(0, 0, hin_dx_select->GetMaximum());
			fit_dx_p->SetParLimits(1, -1.0, -0.85);
			fit_dx_p->SetParLimits(2, 0.1, 0.19);			
		}
		if( kine == 9 && sbsfieldscale == 70){
			fit_dx_p->SetParLimits(0, 0, hin_dx_select->GetMaximum());
			fit_dx_p->SetParLimits(1, -0.95, -0.85);
			fit_dx_p->SetParLimits(2, 0.1, 0.16);			
		}	
	
		hin_dx_select->Fit("fit_dx_p", "R+");
		dx_p = fit_dx_p->GetParameter(1);
		dx_p_sigma = fit_dx_p->GetParameter(2);	

	//------ n -------
    	TF1 *fit_dx_n = new TF1("fit_dx_n", fit_gaus, -0.2, 0.5, 3);
  
    	fit_dx_n->SetParName(0, "dx_n Norm");
		fit_dx_n->SetParName(1, "dx_n Center");
		fit_dx_n->SetParName(2, "dx_n Sigma");
		fit_dx_n->SetLineColor(3);

		if( kine == 4 ){
			fit_dx_n->SetParLimits(0, 0, (0.35)*hin_dx_select->GetMaximum());
			fit_dx_n->SetParLimits(1, 0.0, 0.1);
			fit_dx_n->SetParLimits(2, 0.1, 0.19);
		}	

		if( kine == 8 && sbsfieldscale == 70){
			hin_dx_select->GetXaxis()->SetRangeUser(-dx_n_sigma, dx_n_sigma);
			cout << "Max value for n: " << hin_dx_select->GetMaximum() << endl;
			fit_dx_n->SetParLimits(0, 0, hin_dx_select->GetMaximum());
			fit_dx_n->SetParLimits(1, -0.08, 0.08);
			fit_dx_n->SetParLimits(2, 0.1, 0.18);
			hin_dx_select->GetXaxis()->SetRangeUser(-2.5, 2.5);
		}
		if( kine == 9 && sbsfieldscale == 70){
			fit_dx_n->SetParLimits(0, 0, hin_dx_select->GetMaximum());
			fit_dx_n->SetParLimits(1, -0.05, 0.05);
			fit_dx_n->SetParLimits(2, 0.1, 0.165);
			hin_dx_select->GetXaxis()->SetRangeUser(-2.5, 2.5);
		}	
	
		hin_dx_select->Fit("fit_dx_n", "R+");
		dx_n = fit_dx_n->GetParameter(1);
		dx_n_sigma = fit_dx_n->GetParameter(2);	

	//------- dy -------
    	TCanvas *c_dy = new TCanvas("c_dy", "c_dy", 600, 500);
    	hin_dy_select->Draw();

    	TF1 *fit_dy = new TF1("fit_dy", fit_gaus, -1.5, 1.5, 3);
  
    	fit_dy->SetParName(0, "dy Norm");
		fit_dy->SetParName(1, "dy Center");
		fit_dy->SetParName(2, "dy Sigma");
		fit_dy->SetLineColor(3);

		if( kine == 9 ){
			fit_dy->SetParLimits(0, 0, hin_dy_select->GetMaximum());
			fit_dy->SetParLimits(1, -0.08, 0.08);
			fit_dy->SetParLimits(2, 0.1, 0.19);				
		}
		else{
			fit_dy->SetParLimits(0, 0, hin_dy_select->GetMaximum());
			fit_dy->SetParLimits(1, -0.15, 0.15);
			fit_dy->SetParLimits(2, 0.1, 0.25);				
		}

	
		hin_dy_select->Fit("fit_dy", "R+");
		dy_p = fit_dy->GetParameter(1);
		dy_p_sigma = fit_dy->GetParameter(2);	

	//------- dxdy -------
    	TCanvas *c_dxdy = new TCanvas("c_dxdy", "c_dxdy", 600, 500);
    	hin_dxdy_select->Draw("colz");

    	TCanvas *c_h_dx = new TCanvas("c_h_dx", "c_dx", 600, 500);
    	hin_dx->Draw();

    }

	cout << "------------------------------------------------------------------"<< endl;
	cout << "                       ANALYSIS FINISHED" << endl;
	cout << "------------------------------------------------------------------"<< endl;
	cout << "Run parameters: " << endl;
	cout << "Kinematic: SBS" << kine << endl;
	cout << "Target: " << run_target.Data() << endl;
	cout << "Beam Energy: " << E_beam << endl;
	cout << "SBS Field: " << sbsfieldscale << "%" << endl;
	cout << "-----------------------------------" << endl;
	cout << "BB angle: " << lookup_BB_angle_by_kine(kine, "deg") << endl;
	cout << "SBS angle: " << lookup_SBS_angle_by_kine(kine, "deg") << endl;
	cout << "HCal angle [deg]: " << (180/pi)*HCal_theta << endl;
	cout << "HCal distance: " << HCal_dist << endl;
	cout << "-----------------------------------" << endl << endl;
	cout << "Elastic yield: " << elastic_yield << endl << endl;
	cout << "---------------------------------------" << endl << endl;	
	cout << "------------------------------------------------------------------"<< endl;
	cout << "Cut info:" << endl;
	cout << "Master cut: " << master_cut << endl << endl;
	cout << "Basic cuts applied: " << endl;
	cout << "Pre-shower & Shower clusters each > 0" << endl;
	cout << "Shower + Pre-Shower cut: "<< SH_PS_min << endl;
	cout << "Vertex <= 0.075" << endl;
	cout << "Number of GEM planes hit > 3" << endl;
	cout << "Number of tracks per event = 1" << endl;
	cout << "E/p: center = " << lookup_parsed_cut(run_target, kine, sbsfieldscale, "Ep") << endl;
	cout << "HCal cluster energy: " << lookup_parsed_cut(run_target, kine, sbsfieldscale, "HCal_clus_e_cut");
	cout << "------------------------------------------------------------------"<< endl;
	cout << "Cut lookups: " << endl;
	cout << "PS_min: " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "PS_clus_e_cut") << endl;
	cout << "SH_PS: " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "SH_PS") << "; " << endl;
	cout << "HCal_clus_e: " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "HCal_clus_e_cut") << endl;
	cout << "Ep: " << Ep_min << " min, " << Ep_max << " max; " << endl;
	cout << "W2: " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "W2_min") << "; " << lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "W2_max") << endl;
	cout << "------------------------------------------------------------------"<< endl << endl;
	if( calibrate ){
		cout << "dx_p = : " << dx_p << "; dx_p_sigma = " << dx_p_sigma << endl;
		cout << "dx_n = : " << dx_n << "; dx_n_sigma = " << dx_n_sigma << endl;
		cout << "dy_p = " << dy_p << "; dy_p_sigma = " << dy_p_sigma << endl;
	}
	cout << "------------------------------------------------------------------"<< endl << endl;
	cout << "Output file: " << outfile->GetName() << endl << endl;
	cout << "------------------------------------------------------------------"<< endl;
	
	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;

	TCanvas *c_dx_dycut = new TCanvas("c_dx_dycut", "c_dx_dycut", 600, 500);

	TFile *TF = new TFile(outfile->GetName());
	TH2D *h2 = static_cast<TH2D*>(TF->Get("h_dxdy_wcut_fcut"));
	TH2D *h2copy = static_cast<TH2D*>(TF->Get("h_dxdy_wcut_fcut"));
	h2->GetXaxis()->SetRangeUser(-0.3, 0.3);
	h2->GetYaxis()->SetRangeUser(-2.0, 1.0);
	h2->ProjectionY()->Draw();

	TCanvas *c_dx_dycut_tight = new TCanvas("c_dx_dycut_tight", "c_dx_dycut_tight", 600, 500);
	h2copy->GetXaxis()->SetRangeUser(-0.2, 0.2);
	h2copy->GetYaxis()->SetRangeUser(-2.0, 1.0);
	h2copy->ProjectionY()->Draw();

}