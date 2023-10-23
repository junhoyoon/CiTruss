//#include <metis.h>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <string>
#include <string.h>
#include <vector>
#include <chrono>
#include "mega_scggm.h"
#include "util.h"

#ifdef _OPENMP
	#include <omp.h>
#endif

using namespace std;

void exit_with_help() {
	printf(
		"Usage: ./citruss_sum [options] "
		"num_samples num_genes num_SNPs Ysum_file Xm_file Xp_file reg_V reg_F\n"
    "X,Y: (SNP or gene)-by-samples. Each row is a SNP or gene."
		
		"options:\n"
		"		-o output prefix(./): \n"
		"		--init-V V0_filename(none): filename with initial V\n"
		"		--init-F F0_filename(none): filename with initial F\n"
		"		-v verbose(1): show information or not (0 or 1)\n"
		"		-i max_iters(10): max number of outer iterations\n"
		"		-s sigma(1e-4): backtracking termination criterion\n"
		"		-q tol(1e-2): tolerance for terminating outer loop\n"
		"		-j obj_tol(1e-13): CG tolerance for calculating objective function\n"
		"		-g grad_tol(1e-10): CG tolerance for calculating gradient\n"
		"		-h hess_tol(1e-8): CG tolerance for calculating hessian\n"
		"		-m memory_usage(32000): memory capacity in MB\n"	
		"		-n threads(16) : set the max number of threads\n"		
		"		-r refit(false): update (Lambda0,Theta0) without adding edges\n"
	);
	exit(1);
}

bool strsame_ignorecase(string s1, string s2) {
  transform(s1.begin(), s1.end(), s1.begin(), ::tolower);
  transform(s2.begin(), s2.end(), s2.begin(), ::tolower);
  
  if (s1.compare(s2) != 0)
    return false;
  else
    return true;
}

int main(int argc, char **argv) {

	int num_reqd_args = 8;
	if (argc < 1 + num_reqd_args) {
		fprintf(stderr,"not enough arguments\n");
		exit_with_help();
	}

	if ( (argc - num_reqd_args - 1) % 2 != 0) {
		fprintf(stderr,"option-value pair is incorrectly given\n");
		exit_with_help();
	}

	long n = atol(argv[1]);
	long q = atol(argv[2]);
	long p = atol(argv[3]);

	string Ysum_filename = argv[4];
	string Xm_filename = argv[5];
	string Xp_filename = argv[6];

	double reg_V = atof(argv[7]);
	double reg_F = atof(argv[8]);

  int is_output_given = 0;
  ofstream output_V;
  ofstream output_F;

	for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-o") == 0) {
			string output_prefix(argv[i+1]);
			output_V.open(output_prefix + "V.txt");
			output_F.open(output_prefix + "F.txt");
						
			is_output_given = 1;
			break;
		}
	}
	if (is_output_given == 0) {
		output_V.open("./V.txt");
		output_F.open("./F.txt");
	}

	string V0_filename = "";
	string F0_filename = "";

	for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"--init-V") == 0) {
			V0_filename = argv[i+1];
			break;
		}
	}
	for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"--init-F") == 0) {
			F0_filename = argv[i+1];
			break;
		}
	}
 
  CGGMOptions options_sum;
  options_sum.refit = false;

 
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-v") == 0) {
			options_sum.quiet = (atoi(argv[i+1]) == 0);
			break;
		}
	}

  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-i") == 0) {
			options_sum.max_outer_iters = atoi(argv[i+1]);
			break;
		}
	}
	
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-s") == 0) {
			options_sum.sigma = atof(argv[i+1]);
			break;
		}
	}

	for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-q") == 0) {
			options_sum.tol = atof(argv[i+1]);
			break;
		}
	}
 
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-j") == 0) {
			options_sum.obj_tol = atof(argv[i+1]);
			break;
		}
	}
 
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-g") == 0) {
			options_sum.grad_tol = atof(argv[i+1]);
			break;
		}
	}
 
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-h") == 0) {
			options_sum.hess_tol = atof(argv[i+1]);
			break;
		}
	}
  
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-m") == 0) {
			options_sum.memory_usage = atol(argv[i+1]);
			break;
		}
	}
  
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-n") == 0) {
			options_sum.max_threads = atoi(argv[i+1]);
			break;
		}
	}
  
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-r") == 0) {
			options_sum.refit = atoi(argv[i+1]) != 0;
			break;
		}
	}
  
  auto time_begin = chrono::high_resolution_clock::now();
  //
  // MODEL START
  //
  
  double val;
  string svalm;
  string svalp;  

	// Read input data from file
	vector< vector<double> > Ysum(q, vector<double>(n, 0));
  
	
	ifstream input_Ysum(Ysum_filename);
	for (long j = 0; j < q; j++) {
    for (long i = 0; i < n; i++) {
			if (!input_Ysum.good()) {
				fprintf(stderr, "line %ld column %ld \n", i, j);
				fprintf(stderr, "error reading Ysum_file\n");
				exit_with_help();
			}
			input_Ysum >> val;
			Ysum[j][i] = val;
		}
	}
	input_Ysum.close();

  
  vector< vector<double> > Xsum(p, vector<double>(n, 0));
  { // new scope where Xm, Xp temporarily exist
   
  ifstream input_Xm(Xm_filename);
  ifstream input_Xp(Xp_filename);
	for (long j = 0; j < p; j++) {
		for (long i = 0; i < n; i++) {
      
      double Xm_ji;
      double Xp_ji;
      
			if (!input_Xm.good()) {
				fprintf(stderr, "line %ld column %ld \n", i, j);
				fprintf(stderr, "error reading Xm_file\n");
				exit_with_help();
			}
			input_Xm >> Xm_ji;
      
      if (!input_Xp.good()) {
				fprintf(stderr, "line %ld column %ld \n", i, j);
				fprintf(stderr, "error reading Xp_file\n");
				exit_with_help();
			}
      input_Xp >> Xp_ji;
      Xsum[j][i]  = Xm_ji + Xp_ji;
      //Xdiff[j][i] = Xm[j][i] - Xp[j][i];
		}
	}
	input_Xm.close();
  input_Xp.close(); 
  
  } // Xm,Xp scope ends
  
  if (!options_sum.quiet) {
		fprintf(stdout, "finished reading data\n");
	}
  
  // normalize (center) data
	for (long i = 0; i < q; i++) {
		double mean = 0;
		for (long nn = 0; nn < n; nn++) {
			mean += Ysum[i][nn];
		}
		mean = mean/n;
    
    double std = 0;
    for (long nn = 0; nn < n; nn++) {
			std += (Ysum[i][nn] - mean) * (Ysum[i][nn] - mean);
		}
    std = sqrt(std/(n-1));
       
		for (long nn = 0; nn < n; nn++) {
			Ysum[i][nn] = (Ysum[i][nn] - mean)/std;
		}
	}
  
  for (long i = 0; i < p; i++) {
		double mean = 0;
		for (long nn = 0; nn < n; nn++) {
			mean += Xsum[i][nn];
		}
		mean = mean/n;
    
    double std = 0;
    for (long nn = 0; nn < n; nn++) {
			std += (Xsum[i][nn] - mean) * (Xsum[i][nn] - mean);
		}
    std = sqrt(std/(n-1));
    
		for (long nn = 0; nn < n; nn++) {
			Xsum[i][nn] = (Xsum[i][nn] - mean)/std;
		}
	}
  
	if (!options_sum.quiet) {
		fprintf(stdout, "finished normalizing data for sum\n");
	}
  
  auto time_end = chrono::high_resolution_clock::now();
	double time_elapsed = chrono::duration<double>(time_end - time_begin).count();
	fprintf(stdout, "Reading: %f sec.\n", time_elapsed);
  time_begin = chrono::high_resolution_clock::now();

	// Initialize parameters
	smat_t V;
	vector<Triplet> triplets;
	for (long i = 0; i < q; i++) {
		double tmp = 0;
		for (long nn = 0; nn < n; nn++) {
			tmp += Ysum[i][nn] * Ysum[i][nn];
		}
		if (tmp <= 1e-8) {
			fprintf(stderr, "Ysum variable %li has variance <= 1e-8\n", i);
			exit(1);
		}
		triplets.push_back(Triplet(i, i, 1.0/(0.01 + (1.0/n)*tmp)));
	}
	V.setTriplets(q, triplets);
	sparse_t F(p, q);
	for (long i = 0; i < p; i++) {
		double variance = 0;
		for (long nn = 0; nn < n; nn++) {
			variance += Xsum[i][nn] * Xsum[i][nn];
		}
		if (variance <= 1e-8) {
			fprintf(stderr, "Xsum variable %li has variance <= 1e-8\n", i);
			//exit(1);
		}
	}

	// Initialize V0 if specified by user
	if (!V0_filename.empty()) {
		ifstream ifV(V0_filename.c_str(), ifstream::in);
		long V0_p, V0_q, V0_nnz;
		ifV >> V0_p >> V0_q >> V0_nnz;
		if (V0_p != q || V0_q != q) {
			fprintf(stderr, "error reading V0_file\n");
			exit(1);
		}
		vector<Triplet> triplets;
		long i, j;
		double val;
		for (long nn = 0; nn < V0_nnz; nn++) {
			ifV >> i >> j >> val;
			if (!ifV.good()) {
				fprintf(stderr, "error reading V0_file\n");
				exit(1);
			}
			if (i >= j) {
				triplets.push_back(Triplet(i-1, j-1, val));
			}
		}
		V.setTriplets(q, triplets);
		ifV.close();
	}

	// Initialize F0 if specified by user
	if (!F0_filename.empty()) {
    ifstream ifF(F0_filename.c_str(), ifstream::in);
		long F0_p, F0_q, F0_nnz;
		ifF >> F0_p >> F0_q >> F0_nnz;
		if (F0_p != p || F0_q != q) {
			fprintf(stderr, "error reading F0_file\n");
			exit(1);
		}
		vector<Triplet> triplets;
		long i, j;
		double val;
		for (long nn = 0; nn < F0_nnz; nn++) {
			ifF >> i >> j >> val;
			if (!ifF.good()) {
				fprintf(stderr, "error reading F0_file\n");
				exit(1);
			}
			triplets.push_back(Triplet(i-1, j-1, val));
		}
		F.setTriplets(p, q, triplets);
		ifF.close();
	}

  // Run optimization for sum
  CGGMStats stats_sum;
	mega_scggm(Ysum, Xsum, reg_V, reg_F, options_sum, V, F, stats_sum);
  
  time_end = chrono::high_resolution_clock::now();
	time_elapsed = chrono::duration<double>(time_end - time_begin).count();
	fprintf(stdout, "Whole model: %f sec.\n", time_elapsed);
  time_begin = chrono::high_resolution_clock::now();
  
  // Output results for sum model
  // Output sparse V
	long true_V_nnz = 0;
	for (long idx = 0; idx < V.nnz; idx++) {
		if (V.values[idx] != 0) {
			true_V_nnz++;
		}
	}
	output_V.precision(12);
	//output_V << q << " " << q << " " << true_V_nnz << endl;
	for (long i = 0; i < q; i++) {
		for (long idx = V.row_ptr[i]; idx < V.row_ptr[i+1]; idx++) {
			if (V.values[idx] != 0) {
				output_V << i+1 << " " << V.col_idx[idx]+1 << " " 
					<< V.values[idx] << endl;
			}
		}
	}
	output_V.close();

	// Output sparse F
	long true_F_nnz = 0;
	for (long idx = 0; idx < F.nnz; idx++) {
		if (F.values[idx] != 0) {
			true_F_nnz++;
		}
	}
	output_F.precision(12);
	//output_F << p << " " << q << " " << true_F_nnz << endl;
	
  //fprintf(stdout, "F: nnz:%d\n", F.nnz);
  
  for (long i = 0; i < p; i++) {
		for (long idx = F.row_ptr[i]; idx < F.row_ptr[i+1]; idx++) {
			if (F.values[idx] != 0) {
				output_F << i+1 << " " << F.col_idx[idx]+1 << " " 
					<< F.values[idx] << endl;
			}
      //fprintf(stdout, "(%d,%d): %f\n", i+1, F.col_idx[idx]+1, F.values[idx]);
		}
	}
	output_F.close();
  
  time_end = chrono::high_resolution_clock::now();
	time_elapsed = chrono::duration<double>(time_end - time_begin).count();
	fprintf(stdout, "Saving: %f sec.\n", time_elapsed);
  
  



/*
	// Output stats
	ofstream fS(stats_filename.c_str(), ofstream::out);
	fS.precision(13);
	
 
  fS << "Trans stats" << endl;
	fS << "objval ";
	for (int i = 0; i < stats.objval.size(); i++) {
		fS << stats.objval[i] << " ";
	}
	fS << endl;
	fS << "time ";
	for (int i = 0; i < stats.time.size(); i++) {
		fS << stats.time[i] << " ";
	}
	fS << endl;
	fS << "active_set_size ";
	for (int i = 0; i < stats.active_set_size.size(); i++) {
		fS << stats.active_set_size[i] << " ";
	}
	fS << endl;
	fS << "active_theta ";
	for (int i = 0; i < stats.active_theta.size(); i++) {
		fS << stats.active_theta[i] << " ";
	}
	fS << endl;
	fS << "active_lambda ";
	for (int i = 0; i < stats.active_lambda.size(); i++) {
		fS << stats.active_lambda[i] << " ";
	}
	fS << endl;
	fS << "subgrad ";
	for (int i = 0; i < stats.subgrad.size(); i++) {
		fS << stats.subgrad[i] << " ";
	}
	fS << endl;
	fS << "l1norm ";
	for (int i = 0; i < stats.l1norm.size(); i++) {
		fS << stats.l1norm[i] << " ";
	}
	fS << endl;
	fS << "time_objective ";
	for (int i = 0; i < stats.time_objective.size(); i++) {
		fS << stats.time_objective[i] << " ";
	}
	fS << endl;
	fS << "time_clustering ";
	for (int i = 0; i < stats.time_clustering.size(); i++) {
		fS << stats.time_clustering[i] << " ";
	}
	fS << endl;
	fS << "time_theta_active ";
	for (int i = 0; i < stats.time_theta_active.size(); i++) {
		fS << stats.time_theta_active[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_active ";
	for (int i = 0; i < stats.time_lambda_active.size(); i++) {
		fS << stats.time_lambda_active[i] << " ";
	}
	fS << endl;
	fS << "time_theta_cd ";
	for (int i = 0; i < stats.time_theta_cd.size(); i++) {
		fS << stats.time_theta_cd[i] << " ";
	}
	fS << endl;
	fS << "time_theta_cd_prep ";
	for (int i = 0; i < stats.time_theta_cd_prep.size(); i++) {
		fS << stats.time_theta_cd_prep[i] << " ";
	}
	fS << endl;
	fS << "time_theta_cd_cd ";
	for (int i = 0; i < stats.time_theta_cd_cd.size(); i++) {
		fS << stats.time_theta_cd_cd[i] << " ";
	}
	fS << endl;
	fS << "time_theta_cd_qr ";
	for (int i = 0; i < stats.time_theta_cd_qr.size(); i++) {
		fS << stats.time_theta_cd_qr[i] << " ";
	}
	fS << endl;
	fS << "time_theta_cd_inv ";
	for (int i = 0; i < stats.time_theta_cd_inv.size(); i++) {
		fS << stats.time_theta_cd_inv[i] << " ";
	}
	fS << endl;
	fS << "time_theta_cd_s ";
	for (int i = 0; i < stats.time_theta_cd_s.size(); i++) {
		fS << stats.time_theta_cd_s[i] << " ";
	}
	fS << endl;
	fS << "time_theta_cd_vp ";
	for (int i = 0; i < stats.time_theta_cd_vp.size(); i++) {
		fS << stats.time_theta_cd_vp[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd ";
	for (int i = 0; i < stats.time_lambda_cd.size(); i++) {
		fS << stats.time_lambda_cd[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_prep ";
	for (int i = 0; i < stats.time_lambda_cd_prep.size(); i++) {
		fS << stats.time_lambda_cd_prep[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_prep_z ";
	for (int i = 0; i < stats.time_lambda_cd_prep_z.size(); i++) {
		fS << stats.time_lambda_cd_prep_z[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_prep_q ";
	for (int i = 0; i < stats.time_lambda_cd_prep_q.size(); i++) {
		fS << stats.time_lambda_cd_prep_q[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_prep_d ";
	for (int i = 0; i < stats.time_lambda_cd_prep_d.size(); i++) {
		fS << stats.time_lambda_cd_prep_d[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_cd ";
	for (int i = 0; i < stats.time_lambda_cd_cd.size(); i++) {
		fS << stats.time_lambda_cd_cd[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_cd_prep ";
	for (int i = 0; i < stats.time_lambda_cd_cd_prep.size(); i++) {
		fS << stats.time_lambda_cd_cd_prep[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_cd_apply ";
	for (int i = 0; i < stats.time_lambda_cd_cd_apply.size(); i++) {
		fS << stats.time_lambda_cd_cd_apply[i] << " ";
	}
	fS << endl;
	fS << "time_lambda_cd_ls ";
	for (int i = 0; i < stats.time_lambda_cd_ls.size(); i++) {
		fS << stats.time_lambda_cd_ls[i] << " ";
	}
	fS << endl;
	fS << "boundary_lambda ";
	for (int i = 0; i < stats.boundary_lambda.size(); i++) {
		fS << stats.boundary_lambda[i] << " ";
	}
	fS << endl;
	fS << "boundary_theta ";
	for (int i = 0; i < stats.boundary_theta.size(); i++) {
		fS << stats.boundary_theta[i] << " ";
	}
	fS << endl;
	fS << "blocks_lambda ";
	for (int i = 0; i < stats.blocks_lambda.size(); i++) {
		fS << stats.blocks_lambda[i] << " ";
	}
	fS << endl;
	fS << "blocks_theta ";
	for (int i = 0; i < stats.blocks_theta.size(); i++) {
		fS << stats.blocks_theta[i] << " ";
	}
	fS << endl;

	fS.close();
*/
	return 0;
}

