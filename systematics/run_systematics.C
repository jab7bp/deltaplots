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

int kine_preset = 9;
int kine;
int sbsfieldscale = 70;
TString run_target = "LD2";

TString this_systematic_var_str;
double loop_val_min, loop_val_max, loop_val_step, max_steps;

#include "/w/halla-scshelf2102/sbs/jboyd/include/SYST_classes.h"
#include "/w/halla-scshelf2102/sbs/jboyd/analysis/gmn/deltaplots/systematics/dxdy_systematics.C"

void run_systematics(){

}

//possible systematic_var_str:
// inv_mass

void run_inv_mass_fixed_mean_vary_sigma_mult(int kine_select = -1){

	if( kine_select != -1 ){
		kine = kine_select;
	}
	else{ kine = kine_preset; }

	this_systematic_var_str = "inv_mass_fixed_mean_vary_sigma_mult";

	auto inv_mass_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	if( kine == 4 ){
		loop_val_min = 0.1;
		loop_val_max = 4.0;		
	}
	if( kine == 8 ){
		loop_val_min = 0.1;
		loop_val_max = 4.0;		
	}
	if( kine == 9 ){
		loop_val_min = 0.1;
		loop_val_max = 4.0;		
	}
	loop_val_step  = 0.1;
	max_steps = (loop_val_max - loop_val_min)/loop_val_step;

	cout << "Running systematics on :" << endl;
	cout << this_systematic_var_str.Data() << endl;
	cout << "Min: " << loop_val_min << ", max: " << loop_val_max << ", step: " << loop_val_step << endl;
	cout << "Max steps: " << max_steps << endl;

	int step_cnt = 0;

	cout << "Step: ";
	vector<double> loop_val_vec = {0};

	for( double value = loop_val_min; value <= loop_val_max; value += loop_val_step){
		cout << step_cnt << " ";
		loop_val_vec = {0};
		loop_val_vec[0] = value;
		dxdy_systematics(kine, sbsfieldscale, run_target, "inv_mass_fixed_mean_vary_sigma_mult", value, value );
		
		step_cnt++;

	}
	cout << endl;

	auto inv_mass_time_end = high_resolution_clock::now();
	auto inv_mass_time_duration = duration_cast<minutes>(inv_mass_time_end - inv_mass_time_start);

	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
	cout << "Total time for " << this_systematic_var_str.Data() << "analysis: " << inv_mass_time_duration.count() << " minutes. " << endl;
	cout << endl;
	cout << "Starting inv. mass val: " << loop_val_min << ", end value: " << loop_val_max << ", total steps: " << max_steps << endl;
	cout << endl;
	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;

}

void run_inv_mass_fixed_min_vary_max(int kine_select = -1){

	if( kine_select != -1 ){
		kine = kine_select;
	}
	else{ kine = kine_preset; }

	this_systematic_var_str = "run_inv_mass_fixed_min_vary_max";

	auto inv_mass_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	if( kine_select == 4 ){
		loop_val_min = 0.4;
		loop_val_max = 2.3;		
	}
	else{
		loop_val_min = 0.4;
		loop_val_max = 2.6;		
	}
	loop_val_step  = 0.1;
	
	max_steps = (loop_val_max - loop_val_min)/loop_val_step;

	cout << "Running systematics on :" << endl;
	cout << this_systematic_var_str.Data() << endl;
	cout << "Min: " << loop_val_min << ", max: " << loop_val_max << ", step: " << loop_val_step << endl;
	cout << "Max steps: " << max_steps << endl;

	int step_cnt = 0;

	cout << "Step: ";

	vector<double> loop_val_vec = {0, 0};

	for( double value = loop_val_min + loop_val_step; value <= loop_val_max; value += loop_val_step){
		cout << step_cnt << " ";

		loop_val_vec = {0, 0};
		loop_val_vec = {loop_val_min, value};

		dxdy_systematics(kine, sbsfieldscale, run_target, "inv_mass_fixed_min_vary_max", loop_val_min, value );
		
		step_cnt++;

	}
	cout << endl;

	auto inv_mass_time_end = high_resolution_clock::now();
	auto inv_mass_time_duration = duration_cast<minutes>(inv_mass_time_end - inv_mass_time_start);

	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
	cout << "Total time for " << this_systematic_var_str.Data() << " analysis: " << inv_mass_time_duration.count() << " minutes. " << endl;
	cout << endl;
	cout << "Starting inv. mass val: " << loop_val_min << ", end value: " << loop_val_max << ", total steps: " << max_steps << endl;
	cout << endl;
	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;

}

void inv_mass_fixed_moving_range(int kine_select = -1){

	if( kine_select != -1 ){
		kine = kine_select;
	}
	else{ kine = kine_preset; }

	this_systematic_var_str = "inv_mass_fixed_moving_range";

	auto inv_mass_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	loop_val_min = 0.4;
	loop_val_max = 1.5;
	loop_val_step  = 0.1;
	max_steps = (loop_val_max - loop_val_min)/loop_val_step;

	cout << "Running systematics on :" << endl;
	cout << this_systematic_var_str.Data() << endl;
	cout << "Min: " << loop_val_min << ", max: " << loop_val_max << ", step: " << loop_val_step << endl;
	cout << "Max steps: " << max_steps << endl;

	int step_cnt = 0;

	cout << "Step: ";

	vector<double> loop_val_vec = {0, 0};

	for( double value = loop_val_min; value < loop_val_max; value += loop_val_step){
		cout << step_cnt << " ";

		loop_val_vec = {0, 0};
		loop_val_vec = {value, value + loop_val_step};

		dxdy_systematics(kine, sbsfieldscale, run_target, "inv_mass_fixed_moving_range", value, value + loop_val_step );
		
		step_cnt++;

	}

	cout << endl;

	auto inv_mass_time_end = high_resolution_clock::now();
	auto inv_mass_time_duration = duration_cast<minutes>(inv_mass_time_end - inv_mass_time_start);

	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
	cout << "Total time for " << this_systematic_var_str.Data() << " analysis: " << inv_mass_time_duration.count() << " minutes. " << endl;
	cout << endl;
	cout << "Starting inv. mass val: " << loop_val_min << ", end value: " << loop_val_max << ", total steps: " << max_steps << endl;
	cout << endl;
	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;

}

void run_SHPS_sigma_mult(int kine_select = -1){

	if( kine_select != -1 ){
		kine = kine_select;
	}
	else{ kine = kine_preset; }

	this_systematic_var_str = "SH_PS_sigma_mult";

	auto SHPS_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	if( kine == 4 ){
		loop_val_min = 0.8;
		loop_val_max = 2.6;	
	}
	if( kine == 8 ){
		loop_val_min = 1.7;
		loop_val_max = 4.5;	
	}
	if( kine == 9 ){
		loop_val_min = 0.7;
		loop_val_max = 2.3;

	}
	loop_val_step  = 0.1;	
	max_steps = (loop_val_max - loop_val_min)/loop_val_step;

	cout << "Running systematics on :" << endl;
	cout << this_systematic_var_str.Data() << endl;
	cout << "Min: " << loop_val_min << ", max: " << loop_val_max << ", step: " << loop_val_step << endl;
	cout << "Max steps: " << max_steps << endl;

	int step_cnt = 0;

	cout << "Step: ";
	vector<double> loop_val_vec = {0};

	for( double value = loop_val_min; value <= loop_val_max; value += loop_val_step){
		cout << step_cnt << " ";
		loop_val_vec = {0};
		loop_val_vec[0] = value;
		dxdy_systematics(kine, sbsfieldscale, run_target, "SH_PS_sigma_mult", value, value );
		
		step_cnt++;

	}
	cout << endl;

	auto SHPS_time_end = high_resolution_clock::now();
	auto SHPS_time_duration = duration_cast<minutes>(SHPS_time_end - SHPS_time_start);

	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
	cout << "Total time for " << this_systematic_var_str.Data() << "analysis: " << SHPS_time_duration.count() << " minutes. " << endl;
	cout << endl;
	cout << "Starting SHPS val: " << loop_val_min << ", end value: " << loop_val_max << ", total steps: " << max_steps << endl;
	cout << endl;
	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;

}

void run_Ep_min(int kine_select = -1){

	if( kine_select != -1 ){
		kine = kine_select;
	}
	else{ kine = kine_preset; }

	this_systematic_var_str = "Ep_min";

	auto Ep_min_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	if( kine == 4 ){
		loop_val_min = 0.6;
		loop_val_max = 1.3;	
	}
	if( kine == 8 ){
		loop_val_min = 0.6;
		loop_val_max = 1.4;	
	}
	if( kine == 9 ){
		loop_val_min = 0.75;
		loop_val_max = 1.5;
	}
	loop_val_step  = 0.025;
	max_steps = (loop_val_max - loop_val_min)/loop_val_step;

	cout << "Running systematics on :" << endl;
	cout << this_systematic_var_str.Data() << endl;
	cout << "Min: " << loop_val_min << ", max: " << loop_val_max << ", step: " << loop_val_step << endl;
	cout << "Max steps: " << max_steps << endl;

	int step_cnt = 0;

	cout << "Step: ";
	vector<double> loop_val_vec = {0};

	for( double value = loop_val_min; value <= loop_val_max; value += loop_val_step){
		cout << step_cnt << " ";
		loop_val_vec = {0};
		loop_val_vec[0] = value;
		dxdy_systematics(kine, sbsfieldscale, run_target, "Ep_min", value, value );
		
		step_cnt++;

	}
	cout << endl;

	auto Ep_min_time_end = high_resolution_clock::now();
	auto Ep_min_time_duration = duration_cast<minutes>(Ep_min_time_end - Ep_min_time_start);

	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
	cout << "Total time for " << this_systematic_var_str.Data() << "analysis: " << Ep_min_time_duration.count() << " minutes. " << endl;
	cout << endl;
	cout << "Starting Ep_min val: " << loop_val_min << ", end value: " << loop_val_max << ", total steps: " << max_steps << endl;
	cout << endl;
	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;

}

void run_PS_min(int kine_select = -1){

	if( kine_select != -1 ){
		kine = kine_select;
	}
	else{ kine = kine_preset; }

	this_systematic_var_str = "PS_min";

	auto PS_min_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	if( kine == 4 ){
		loop_val_min = 0.0;
		loop_val_max = 1.5;	
	}
	if( kine == 8 ){
		loop_val_min = 0.0;
		loop_val_max = 2.0;	
	}
	if( kine == 9 ){
		loop_val_min = 0.0;
		loop_val_max = 1.2;
	}
	loop_val_step  = 0.05;
	max_steps = (loop_val_max - loop_val_min)/loop_val_step;

	cout << "Running systematics on :" << endl;
	cout << this_systematic_var_str.Data() << endl;
	cout << "Min: " << loop_val_min << ", max: " << loop_val_max << ", step: " << loop_val_step << endl;
	cout << "Max steps: " << max_steps << endl;

	int step_cnt = 0;

	cout << "Step: ";
	vector<double> loop_val_vec = {0};

	for( double value = loop_val_min; value <= loop_val_max; value += loop_val_step){
		cout << step_cnt << " ";
		loop_val_vec = {0};
		loop_val_vec[0] = value;
		dxdy_systematics(kine, sbsfieldscale, run_target, "PS_min", value, value );
		
		step_cnt++;

	}
	cout << endl;

	auto PS_min_time_end = high_resolution_clock::now();
	auto PS_min_time_duration = duration_cast<minutes>(PS_min_time_end - PS_min_time_start);

	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
	cout << "Total time for " << this_systematic_var_str.Data() << "analysis: " << PS_min_time_duration.count() << " minutes. " << endl;
	cout << endl;
	cout << "Starting PS_min val: " << loop_val_min << ", end value: " << loop_val_max << ", total steps: " << max_steps << endl;
	cout << endl;
	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;

}

void run_HCal_clus_e_min(int kine_select = -1){

	if( kine_select != -1 ){
		kine = kine_select;
	}
	else{ kine = kine_preset; }

	this_systematic_var_str = "HCal_clus_e_min";

	auto HCal_clus_e_min_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	if( kine == 4 ){
		loop_val_min = 0.0;
		loop_val_max = 0.25;	
	}
	if( kine == 8 ){
		loop_val_min = 0.0;
		loop_val_max = 0.2;	
	}
	if( kine == 9 ){
		loop_val_min = 0.0;
		loop_val_max = 0.2;
	}
	loop_val_step  = 0.005;
	max_steps = (loop_val_max - loop_val_min)/loop_val_step;

	cout << "Running systematics on :" << endl;
	cout << this_systematic_var_str.Data() << endl;
	cout << "Min: " << loop_val_min << ", max: " << loop_val_max << ", step: " << loop_val_step << endl;
	cout << "Max steps: " << max_steps << endl;

	int step_cnt = 0;

	cout << "Step: ";
	vector<double> loop_val_vec = {0};

	for( double value = loop_val_min; value <= loop_val_max; value += loop_val_step){
		cout << step_cnt << " ";
		loop_val_vec = {0};
		loop_val_vec[0] = value;
		dxdy_systematics(kine, sbsfieldscale, run_target, "HCal_clus_e_min", value, value );
		
		step_cnt++;

	}
	cout << endl;

	auto HCal_clus_e_min_time_end = high_resolution_clock::now();
	auto HCal_clus_e_min_time_duration = duration_cast<minutes>(HCal_clus_e_min_time_end - HCal_clus_e_min_time_start);

	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
	cout << "Total time for " << this_systematic_var_str.Data() << "analysis: " << HCal_clus_e_min_time_duration.count() << " minutes. " << endl;
	cout << endl;
	cout << "Starting HCal_clus_e_min val: " << loop_val_min << ", end value: " << loop_val_max << ", total steps: " << max_steps << endl;
	cout << endl;
	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;

}

void run_gem_track_nhits_min(int kine_select = -1){

	if( kine_select != -1 ){
		kine = kine_select;
	}
	else{ kine = kine_preset; }

	this_systematic_var_str = "gem_track_nhits_min";

	auto gem_track_nhits_min_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	loop_val_min = 1;
	loop_val_max = 3;
	loop_val_step  = 1;
	max_steps = (loop_val_max - loop_val_min - 1)/loop_val_step;

	cout << "Running systematics on :" << endl;
	cout << this_systematic_var_str.Data() << endl;
	cout << "Min: " << loop_val_min << ", max: " << loop_val_max << ", step: " << loop_val_step << endl;
	cout << "Max steps: " << max_steps << endl;

	int step_cnt = 0;

	cout << "Step: ";
	vector<double> loop_val_vec = {0};

	for( double value = loop_val_min; value <= loop_val_max; value += loop_val_step){
		cout << step_cnt << " ";
		loop_val_vec = {0};
		loop_val_vec[0] = value;
		dxdy_systematics(kine, sbsfieldscale, run_target, "gem_track_nhits_min", value, value );
		
		step_cnt++;

	}
	cout << endl;

	auto gem_track_nhits_min_time_end = high_resolution_clock::now();
	auto gem_track_nhits_min_time_duration = duration_cast<minutes>(gem_track_nhits_min_time_end - gem_track_nhits_min_time_start);

	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
	cout << "Total time for " << this_systematic_var_str.Data() << "analysis: " << gem_track_nhits_min_time_duration.count() << " minutes. " << endl;
	cout << endl;
	cout << "Starting gem_track_nhits_min val: " << loop_val_min << ", end value: " << loop_val_max << ", total steps: " << max_steps << endl;
	cout << endl;
	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;

}

void run_fcut_mult(int kine_select = -1){

	if( kine_select != -1 ){
		kine = kine_select;
	}
	else{ kine = kine_preset; }

	this_systematic_var_str = "fcut_mult";

	auto fcut_mult_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	loop_val_min = 1.5;
	loop_val_max = 2.0;
	loop_val_step  = 0.1;
	max_steps = (loop_val_max - loop_val_min)/loop_val_step;

	cout << "Running systematics on :" << endl;
	cout << this_systematic_var_str.Data() << endl;
	cout << "Min: " << loop_val_min << ", max: " << loop_val_max << ", step: " << loop_val_step << endl;
	cout << "Max steps: " << max_steps << endl;

	int step_cnt = 0;

	cout << "Step: ";
	vector<double> loop_val_vec = {0};

	for( double value = loop_val_min; value <= loop_val_max; value += loop_val_step){
		cout << step_cnt << " ";
		loop_val_vec = {0};
		loop_val_vec[0] = value;
		dxdy_systematics(kine, sbsfieldscale, run_target, "fcut_mult", value, value );
		
		step_cnt++;

	}
	cout << endl;

	auto fcut_mult_time_end = high_resolution_clock::now();
	auto fcut_mult_time_duration = duration_cast<minutes>(fcut_mult_time_end - fcut_mult_time_start);

	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
	cout << "Total time for " << this_systematic_var_str.Data() << "analysis: " << fcut_mult_time_duration.count() << " minutes. " << endl;
	cout << endl;
	cout << "Starting fcut_mult val: " << loop_val_min << ", end value: " << loop_val_max << ", total steps: " << max_steps << endl;
	cout << endl;
	cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;

}
