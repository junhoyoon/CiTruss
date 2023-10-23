//#include <metis.h>
#include <unordered_map>
#include <iostream>
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
		"Usage: ./citruss_diff [options] "
		"num_samples num_genes num_SNPs Ym_file Yp_file Xm_file Xp_file gene_info snp_info reg_Psi\n"
    "X,Y: (SNP or gene)-by-samples. Each row is a SNP or gene."
		
		"options:\n"
		"		-o output prefix(./): \n"
		"		--init-Gamma Gamma0_filename(none): filename with initial Gamma\n"
		"		--init-Psi Psi0_filename(none): filename with initial Psi\n"
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
	int num_reqd_args = 10;
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

	string Ym_filename = argv[4];
	string Yp_filename = argv[5];
	string Xm_filename = argv[6];
	string Xp_filename = argv[7];

  string gene_info_filename = argv[8];
  string snp_info_filename = argv[9];

	double reg_Gamma = 0;
	double reg_Psi = atof(argv[10]);

  int is_output_given = 0;
  ofstream output_Gamma;
  ofstream output_Psi;

	for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-o") == 0) {
			string output_prefix(argv[i+1]);
			output_Gamma.open(output_prefix + "Gamma.txt");
			output_Psi.open(output_prefix + "Psi.txt");
						
			is_output_given = 1;
			break;
		}
	}
	if (is_output_given == 0) {
		output_Gamma.open("./Gamma.txt");
		output_Psi.open("./Psi.txt");
	}
  
	string Gamma0_filename = "";
	string Psi0_filename = "";

	for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"--init-Gamma") == 0) {
			Gamma0_filename = argv[i+1];
			break;
		}
	}
	for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"--init-Psi") == 0) {
			Psi0_filename = argv[i+1];
			break;
		}
	}
  CGGMOptions options_diff;
  options_diff.refit = false;

  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-v") == 0) {
      options_diff.quiet = (atoi(argv[i+1]) == 0);
			break;
		}
	}

  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-i") == 0) {
      options_diff.max_outer_iters = atoi(argv[i+1]);
			break;
		}
	}
	
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-s") == 0) {
      options_diff.sigma = atof(argv[i+1]);
			break;
		}
	}

	for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-q") == 0) {
      options_diff.tol = atof(argv[i+1]);
			break;
		}
	}
 
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-j") == 0) {
      options_diff.obj_tol = atof(argv[i+1]);
			break;
		}
	}
 
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-g") == 0) {
      options_diff.grad_tol = atof(argv[i+1]);
			break;
		}
	}
 
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-h") == 0) {
      options_diff.hess_tol = atof(argv[i+1]);
			break;
		}
	}
  
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-m") == 0) {
      options_diff.memory_usage = atol(argv[i+1]);
			break;
		}
	}
  
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-n") == 0) {
      options_diff.max_threads = atoi(argv[i+1]);
			break;
		}
	}
  
  for (int i = num_reqd_args+1; i < argc-1; i++) {
		if (strcmp(argv[i],"-r") == 0) {
      options_diff.refit = atoi(argv[i+1]) != 0;
			break;
		}
	}
  
  auto time_begin = chrono::high_resolution_clock::now();
  
  double val;
  string svalm;
  string svalp;
  
  vector< vector<double> > Ydiff(q);
  vector< vector<double> > Xdiff(p, vector<double>(n, 0));
  
  vector< vector<bool> > Ynan(q, vector<bool>(n, false));

  // Read Y's and compute Ydiff
  
  ifstream input_Ym(Ym_filename);
	ifstream input_Yp(Yp_filename);
	for (long j = 0; j < q; j++) {
    vector<double> tmp_input;
    for (long i = 0; i < n; i++) {
			if (!input_Ym.good()) {
				fprintf(stderr, "line %ld column %ld \n", i, j);
				fprintf(stderr, "error reading Ym_file\n");
				exit_with_help();
			}
      if (!input_Yp.good()) {
				fprintf(stderr, "line %ld column %ld \n", i, j);
				fprintf(stderr, "error reading Yp_file\n");
				exit_with_help();
			}
			input_Ym >> svalm;
      input_Yp >> svalp;
      if (!(strsame_ignorecase(svalm, "nan") || strsame_ignorecase(svalp, "nan")))
        tmp_input.push_back(stod(svalm) - stod(svalp));
      else
        Ynan[j][i] = true;
		}
    Ydiff[j] = tmp_input;
	}
	input_Ym.close();
  input_Yp.close();

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
      //Xsum[j][i]  = Xm_ji + Xp_ji;
      Xdiff[j][i] = Xm_ji - Xp_ji;
		}
	}
	input_Xm.close();
  input_Xp.close();

  // read gene and snp chr
  vector<int> gene_chr(q);
  vector<int> snp_chr; // saves start snp_idx for chromosomes
  
  string gene_start;
  string gene_end;
  
  ifstream input_gene_info(gene_info_filename);
  for (long i = 0; i < q; i++) {
    if (!input_gene_info.good()) {
			fprintf(stderr, "line %ld\n", i+1);
			fprintf(stderr, "error reading gene_info\n");
			exit_with_help();
		}
    input_gene_info >> gene_chr[i];
    input_gene_info >> gene_start;
    input_gene_info >> gene_end;
  }
  input_gene_info.close();

  string snp_pos;
  int prev_snp_chr = 0;
  ifstream input_snp_info(snp_info_filename);
  for (long i = 0; i < p; i++) {
    if (!input_snp_info.good()) {
			fprintf(stderr, "line %ld\n", i+1);
			fprintf(stderr, "error reading snp_info\n");
			exit_with_help();
		}
    int curr_snp_chr;
    input_snp_info >> curr_snp_chr;
    input_snp_info >> snp_pos;
    
    if (curr_snp_chr != prev_snp_chr) {
      snp_chr.push_back(i);
      prev_snp_chr = curr_snp_chr;
    }
  }
  input_snp_info.close();
  
  vector<CGGMStats> stats_diff(q);
  
  output_Gamma.precision(12);
  output_Psi.precision(12);
  
  for (long i = 0; i < q; i++) {
    long ni = Ydiff[i].size();
    vector< vector<double> > Ydiffi(1, vector<double>(ni, 0));
    
    // decide target SNPs
    int curr_gene_chr = gene_chr[i];
    long target_start = snp_chr[curr_gene_chr-1];
    long target_end;
    if (curr_gene_chr == snp_chr.size()) { // last chromosome
      target_end = p; // exclusive
    } else {
      target_end = snp_chr[curr_gene_chr]; // exclusive
    }
    vector< vector<double> > Xdiffi(target_end-target_start, vector<double>(ni, 0));

    double mean = 0;
		for (long nn = 0; nn < ni; nn++) {
			mean += Ydiff[i][nn];
		}
		mean = mean/ni;
    
    double std = 0;
    for (long nn = 0; nn < ni; nn++) {
			std += (Ydiff[i][nn] - mean) * (Ydiff[i][nn] - mean);
		}
    std = sqrt(std/(ni-1));
    
		for (long nn = 0; nn < ni; nn++) {
			Ydiff[i][nn] = (Ydiff[i][nn] - mean)/std;
		}
    Ydiffi[0] = Ydiff[i];
    
    for (long j = target_start; j < target_end; j++) {
  		mean = 0;
      std = 0;
  		for (long nn = 0; nn < n; nn++) {
        if (!Ynan[i][nn])
  			  mean += Xdiff[j][nn];
  		}
  		mean = mean/ni;
      
      for (long nn = 0; nn < n; nn++) {
        if (!Ynan[i][nn])
  			  std += (Xdiff[j][nn] - mean) * (Xdiff[j][nn] - mean);
  		}
      std = sqrt(std/(ni-1));
      
      int sample_count = 0;
  		for (long nn = 0; nn < n; nn++) {
        if (!Ynan[i][nn]) {
          Xdiffi[j-target_start][sample_count] = (Xdiff[j][nn] - mean)/std;
          sample_count++;
        }
  		}
  	}
    
    if (!options_diff.quiet) {
  		fprintf(stdout, "finished normalizing data for diff %li\n", i);
  	}
    fprintf(stdout, "target start=%li target_end=%li\n", target_start, target_end);
    smat_t Gamma(1);
    //vector<Triplet> triplets;
    //triplets.push_back(Triplet(0,0,1.0));
    //Gamma.setTriplets(1, triplets);
    Gamma.reset();
	  sparse_t Psi(target_end-target_start,1);
    //Psi.setZeros(target_end-target_start,1);
    
    if (Gamma0_filename.empty()) {
      //smat_t Gammai;
    	vector<Triplet> triplets;
  		double tmp = 0;
  		for (long nn = 0; nn < ni; nn++) {
  			tmp += Ydiffi[0][nn] * Ydiffi[0][nn];
  		}
  		if (tmp <= 1e-8) {
  			fprintf(stderr, "Ydiff variable %li has variance <= 1e-8\n", i);
  			exit(1);
  		}
  		triplets.push_back(Triplet(0, 0, 1.0/(0.01 + (1.0/ni)*tmp)));
    	Gamma.setTriplets(1, triplets);
    }/* else {
      vector<Triplet> triplets;
      triplets.push_back(Triplet(0, 0, Gamma0[i]));
      Gamma.setTriplets(1, triplets);
    }*/
    //Gamma.print(stdout);
    
    if (!Psi0_filename.empty()) {
      ifstream ifPsi(Psi0_filename.c_str(), ifstream::in);
  		long Psi0_p, Psi0_q, Psi0_nnz;
  		ifPsi >> Psi0_p >> Psi0_q >> Psi0_nnz;
  		if (Psi0_p != p || Psi0_q != q) {
  			fprintf(stderr, "error reading Psi0_file\n");
  			exit(1);
  		}
  		vector<Triplet> triplets;
  		long ii, jj;
  		double val;
  		for (long nn = 0; nn < Psi0_nnz; nn++) {
  			ifPsi >> ii >> jj >> val;
  			if (!ifPsi.good()) {
  				fprintf(stderr, "error reading Psi0_file\n");
  				exit(1);
  			}
        if (ii-1 > target_end || ii-1 < target_start)
          continue;
        if (jj-1 == i) {
  			  triplets.push_back(Triplet(ii-1-target_start, 0, val));
        }
  		}
      print(triplets);
  		Psi.setTriplets(target_end-target_start, 1, triplets);
  		ifPsi.close();
    }
    //fprintf(stdout, "outside4 %ld %ld\n", Xdiffi.size(), Xdiffi[0].size());
    fprintf(stdout,"DIFF %d\n",i);
    fflush(stdout);
  	mega_scggm(Ydiffi, Xdiffi, reg_Gamma, reg_Psi, options_diff, Gamma, Psi, stats_diff[i]);
    
    
    
    
    //fprintf(stdout, "Gamma: %d, %f\n", Gamma.nnz, Gamma.values[0]);
	  output_Gamma << Gamma.values[0] << "\n";
     
    //fprintf(stdout, "Psi: nnz:%d\n", Psi.nnz);
    for (long ii = 0; ii < target_end-target_start; ii++) {  
        for (long idx = Psi.row_ptr[ii]; idx < Psi.row_ptr[ii+1]; idx++) {
    			if (Psi.values[idx] != 0) {
    				output_Psi << ii+1+target_start << " " << i+1 << " " 
    					<< Psi.values[idx] << endl;
    			}
          //fprintf(stdout, "(%d,%d): %f\n", ii+1+target_start, i+1, Psi.values[idx]);
    		}
    }
	}
  output_Gamma.close();
  output_Psi.close();
  
  auto time_end = chrono::high_resolution_clock::now();
	double time_elapsed = chrono::duration<double>(time_end - time_begin).count();
	fprintf(stdout, "Elapsed time is %f sec.\n", time_elapsed);

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
	// TODO- use hdf5 format
	return 0;
}

