// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;

// [[Rcpp::export]]
arma::vec regress_m(arma::vec z, arma::vec x, arma::vec alpha, String family){
    
    arma::vec mu;
    arma::mat zx_vec = join_cols(ones<vec>(1), z, x);
	double eta = as_scalar(zx_vec.t() * alpha);
	if(family == "gaussian"){
		mu = eta;
	} else if(family == "binomial"){
		mu = exp(eta)/(1 + exp(eta));
	} else {Rcpp::stop("Unknown distribution family\n");}
	return mu;
}

// [[Rcpp::export]]
arma::vec regress_m_div(arma::vec z, arma::vec x, arma::vec alpha, String family){
    arma::vec m_div;
    arma::vec zx_vec = join_cols(ones<vec>(1), z, x);
    double eta = as_scalar(zx_vec.t() * alpha);
    if(family == "gaussian"){
    	m_div = zx_vec;
    }else if (family == "binomial"){
    	m_div = zx_vec * exp(eta)/(pow((1 + exp(eta)), 2));
    }else {Rcpp::stop("Unknown distribution family\n");}
	return m_div;
}	

// [[Rcpp::export]]
arma::vec expit(arma::vec x_vec, arma::vec beta){
	arma::vec mu;
    double eta = as_scalar(x_vec.t() * beta);
	mu = exp(eta)/(1 + exp(eta));
	return mu;
}

// [[Rcpp::export]]
arma::vec expit_div(arma::vec x_vec, arma::vec beta){
	arma::vec expit_div;
    double eta = as_scalar(x_vec.t() * beta);
    expit_div = x_vec * exp(eta)/(pow((1 + exp(eta)), 2));
    return expit_div;
}


//[[Rcpp::export]]
List increDR(arma::mat X, arma::mat W, arma::mat Z, arma::mat Y, String family, arma::vec theta_old,
	arma::mat S_accum, arma::mat V_accum, int p1, int p2, int maxit, double tol){

	int niter = 0;
	bool stop_flag = FALSE;
	bool converged = FALSE;

	int d = p1 + p2 + 1;
	int n = Y.n_rows;

	arma::mat S_temp;
	arma::mat V_temp;
	arma::mat s_i1;
	arma::mat s_i2;
	arma::mat s_i3;
    
    // initialize theta_b with theta_{b-1}
	arma::vec theta_new = theta_old;

	while (!stop_flag){

		niter += 1;

		arma::vec alpha_new = theta_new.subvec(0, p1 - 1);
		arma::vec beta_new = theta_new.subvec(p1, p1 + p2 - 1);
		arma::vec delta_new = theta_new[d-1] * ones<vec>(1);

		arma::vec Ub = zeros<vec>(d);
		arma::mat Sb = zeros<mat>(d, d);
		arma::mat Vb = zeros<mat>(d, d);

		for (int i = 0; i < n; i++){
            arma::vec xi = X.row(i).t();
            arma::vec wi = W.row(i).t();
            arma::vec zi = Z.row(i).t();
            arma::vec yi = Y.row(i).t();

            arma::vec zx_vec = join_cols(ones<vec>(1), zi, xi); // p1 * 1, stack three elements
            arma::vec x_vec = join_cols(ones<vec>(1), wi); // p2 * 1 (p1 - p2 = 2 if interaction)
          //  arma::vec x1_vec = join_cols(wi, wi.row(1));
          //  arma::vec x0_vec = join_cols(wi, zeros<vec>(1));

            arma::vec x1_vec = wi;
            arma::vec x0_vec = wi;

            vec ei = expit(x_vec, beta_new);
            vec ei_div = expit_div(x_vec, beta_new); // p2 * 1
            vec m1 = regress_m(ones<vec>(1), x1_vec, alpha_new, family);  // the first entry is zi
            vec m0 = regress_m(zeros<vec>(1), x0_vec, alpha_new, family);
            vec m1_div = regress_m_div(ones<vec>(1), x1_vec, alpha_new, family); // p1
            vec m0_div = regress_m_div(zeros<vec>(1), x0_vec, alpha_new, family); // p1

            arma::vec u_i1 = zx_vec * (yi - zi * m1 - (1 - zi) * m0);
            arma::vec u_i2 = x_vec * (zi - ei);
            arma::vec u_i3 = (zi * yi / ei - (zi - ei) / ei * m1
            	- (1 - zi) * yi / (1 - ei) - (zi - ei) * m0 / (1 - ei)) - delta_new;
            
            arma::vec u_i = join_cols(u_i1, u_i2, u_i3);

            Ub = Ub + u_i;

            Vb = Vb + u_i * u_i.t();

            s_i1 = join_rows(zx_vec * (m1_div * zi + m0_div * (1 - zi)).t(), 
            	zeros<mat>(p1, p2), zeros<mat>(p1, 1)); // p1 * d

            s_i2 = join_rows(zeros<mat>(p2, p1), 
            	x_vec * ei_div.t(), zeros<mat>(p2, 1)); // p2 * d

            s_i3 = join_rows((zi - ei) / ei * m1_div.t() + (zi - ei) / (1 - ei) * m0_div.t(),
            	zi * yi  / pow(ei, 2) * ei_div.t() - zi * m1 / pow(ei, 2) * ei_div.t() 
            	+ (1 - zi) * yi / pow((1 - ei), 2) * ei_div.t() 
            	+ (zi - 1) * m0 / pow((1 - ei), 2) * ei_div.t(), 
            	ones<vec>(1));  // 1 * d
            
            Sb = Sb + join_cols(s_i1, s_i2, s_i3);
		}

		vec U_temp = S_accum * (theta_old - theta_new) + Ub;
		S_temp = S_accum + Sb;
        V_temp = V_accum + Vb;

		vec d_theta = solve(S_temp, U_temp);

		double df_theta = as_scalar(U_temp.t() * d_theta);

		theta_new += d_theta;

		if(fabs(df_theta) < tol) {converged = TRUE; stop_flag = TRUE;}
		if(niter > maxit) {stop_flag = TRUE;}
	}
//	if(converged == FALSE) {Rprintf("algorithm reached 'maxit' but did not converge\n");}
    
	return List::create(Named("theta") = theta_new,
		Named("S_accum") = S_temp,
		Named("V_accum") = V_temp);
}


//[[Rcpp::export]]
List increRM(arma::mat X, arma::mat Z, arma::mat Y, String family, arma::vec theta_old,
	arma::mat S_accum, arma::mat V_accum, int p1, int maxit, double tol){

	int niter = 0;
	bool stop_flag = FALSE;
	bool converged = FALSE;

	int d = p1 + 1;
	int n = Y.n_rows;
    
    arma::vec U_temp;
	arma::mat S_temp;
	arma::mat V_temp;
	arma::mat s_i1;
	arma::mat s_i3;

    // initialize theta_b with theta_{b-1}
	arma::vec theta_new = theta_old;

	while (!stop_flag){

		niter += 1;

		arma::vec alpha_new = theta_new.subvec(0, p1 - 1);
		arma::vec delta_new = theta_new[d - 1] * ones<vec>(1);

		arma::vec Ub = zeros<vec>(d);
		arma::mat Sb = zeros<mat>(d, d);
		arma::mat Vb = zeros<mat>(d, d);

		for (int i = 0; i < n; i++){
            arma::vec xi = X.row(i).t();
            arma::vec zi = Z.row(i).t();
            arma::vec yi = Y.row(i).t();

         //   arma::vec x1_vec = join_cols(xi.subvec(0, p1 - 4), xi.row(1)); // remove interaction term
         //   arma::vec x0_vec = join_cols(xi.subvec(0, p1 - 4), zeros<vec>(1));

            arma::vec x1_vec = xi; 
            arma::vec x0_vec = xi;

            arma::vec zx_vec = join_cols(ones<vec>(1), zi, xi); // p1 * 1

            vec m1 = regress_m(ones<vec>(1), x1_vec, alpha_new, family);
            vec m0 = regress_m(zeros<vec>(1), x0_vec, alpha_new, family);
            vec m1_div = regress_m_div(ones<vec>(1), x1_vec, alpha_new, family); // p1
            vec m0_div = regress_m_div(zeros<vec>(1), x0_vec, alpha_new, family); // p1

            arma::vec u_i1 = zx_vec * (yi - zi * m1 - (1 - zi) * m0);
            arma::vec u_i3 = m1 - m0 - delta_new;
            
            arma::vec u_i = join_cols(u_i1, u_i3);

            Ub = Ub + u_i;

            Vb = Vb + u_i * u_i.t();

            s_i1 = join_rows(zx_vec * (m1_div * zi + m0_div * (1 - zi)).t(), zeros<mat>(p1, 1)); // p1 * d

            s_i3 = join_rows(- m1_div.t() + m0_div.t(), ones<vec>(1));  // 1 * d
            
            arma::mat s_i = join_cols(s_i1, s_i3);
            
            Sb = Sb + s_i;

		}

		U_temp = S_accum * (theta_old - theta_new) + Ub;
		S_temp = S_accum + Sb;
        V_temp = V_accum + Vb;


		vec d_theta = solve(S_temp, U_temp);

		double df_theta = as_scalar(U_temp.t() * d_theta);

		theta_new += d_theta;

		if(fabs(df_theta) < tol) {converged = TRUE; stop_flag = TRUE;}
		if(niter > maxit) {stop_flag = TRUE;}
	}
//	if(converged == FALSE) {Rprintf("algorithm reached 'maxit' but did not converge\n");}
    
	return List::create(Named("theta") = theta_new,
		Named("U_temp") = U_temp,
		Named("S_accum") = S_temp,
		Named("V_accum") = V_temp);
}


//[[Rcpp::export]]
List increIW(arma::mat W, arma::mat Z, arma::mat Y, String family, arma::vec theta_old,
	arma::mat S_accum, arma::mat V_accum, int p2, int maxit, double tol){

	int niter = 0;
	bool stop_flag = FALSE;
	bool converged = FALSE;

	int d = p2 + 1;
	int n = Y.n_rows;

	arma::mat S_temp;
	arma::mat V_temp;
	arma::mat s_i2;
	arma::mat s_i3;
    
    // initialize theta_b with theta_{b-1}
	arma::vec theta_new = theta_old;

	while (!stop_flag){

		niter += 1;

		arma::vec beta_new = theta_new.subvec(0, p2 - 1);
		arma::vec delta_new = theta_new[d-1] * ones<vec>(1);

		arma::vec Ub = zeros<vec>(d);
		arma::mat Sb = zeros<mat>(d, d);
		arma::mat Vb = zeros<mat>(d, d);

		for (int i = 0; i < n; i++){
         //   arma::vec xi = X.row(i).t();
            arma::vec wi = W.row(i).t();
            arma::vec zi = Z.row(i).t();
            arma::vec yi = Y.row(i).t();

        //    arma::vec zx_vec = join_cols(ones<vec>(1), zi, xi); // p1 * 1
            arma::vec x_vec = join_cols(ones<vec>(1), wi); // p2 * 1

            vec ei = expit(x_vec, beta_new);
            vec ei_div = expit_div(x_vec, beta_new); // p2 * 1
        
            arma::vec u_i2 = x_vec * (zi - ei);
            arma::vec u_i3 = zi * yi / ei - (1 - zi) * yi / (1 - ei) - delta_new;
            
            arma::vec u_i = join_cols(u_i2, u_i3);

            Ub = Ub + u_i;

            Vb = Vb + u_i * u_i.t();

            s_i2 = join_rows(x_vec * ei_div.t(), zeros<mat>(p2, 1)); // p2 * d

            s_i3 = join_rows(
            	zi * yi  / pow(ei, 2) * ei_div.t() 
            	+ (1 - zi) * yi / pow((1 - ei), 2) * ei_div.t(), 
            	ones<vec>(1));  // 1 * d
            
            Sb = Sb + join_cols(s_i2, s_i3);
		}

		vec U_temp = S_accum * (theta_old - theta_new) + Ub;
		S_temp = S_accum + Sb;
        V_temp = V_accum + Vb;

		vec d_theta = solve(S_temp, U_temp);

		double df_theta = as_scalar(U_temp.t() * d_theta);

		theta_new += d_theta;

		if(fabs(df_theta) < tol) {converged = TRUE; stop_flag = TRUE;}
		if(niter > maxit) {stop_flag = TRUE;}
	}
//	if(converged == FALSE) {Rprintf("algorithm reached 'maxit' but did not converge\n");}
    
	return List::create(Named("theta") = theta_new,
		Named("S_accum") = S_temp,
		Named("V_accum") = V_temp);
}



