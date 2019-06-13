# include <RcppArmadillo.h> 
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
void beta_PPL_loop(arma::mat& R1_j,arma::mat& R2_j, arma::mat& z1_j, arma::mat& z2_j,arma::vec& delta1_j,arma::vec& delta2_j,arma::vec& exp_eta1,
                   arma::vec& exp_eta2, arma::vec& beta1_scr, arma::vec& beta2_scr,arma::vec& r_scr,
                   arma::mat& beta1_inf, arma::mat& beta2_inf,arma::mat& r_inf,arma::mat& beta1_r_inf,arma::mat& beta2_r_inf,int m,
                   arma::vec& eta1,arma::vec& eta2) {
  
  arma::mat z1_R1_temp(z1_j.n_cols,R1_j.n_cols);
  z1_R1_temp.fill(0);
  arma::mat z2_R2_temp(z2_j.n_cols,R2_j.n_cols);
  z2_R2_temp.fill(0);
  arma::mat z1_sq_temp(z1_j.n_cols,z1_j.n_cols);
  z1_sq_temp.fill(0);  
  arma::mat z2_sq_temp(z2_j.n_cols,z2_j.n_cols);
  z2_sq_temp.fill(0); 
  
  arma::vec z1_temp(z1_j.n_cols);
  z1_temp.fill(0);
  arma::vec z2_temp(z2_j.n_cols);
  z2_temp.fill(0);  
  
  arma::vec z1_j_bar(z1_j.n_cols);
  z1_j_bar.fill(0);
  arma::vec z2_j_bar(z2_j.n_cols);
  z2_j_bar.fill(0); 
  
  arma::vec R1_temp(R1_j.n_cols);
  R1_temp.fill(0);
  arma::vec R2_temp(R2_j.n_cols);
  R2_temp.fill(0);  
  
  arma::vec R1_j_bar(R1_j.n_cols);
  R1_j_bar.fill(0);
  arma::vec R2_j_bar(R2_j.n_cols);
  R2_j_bar.fill(0);  
  
  double denom_1=0;
  double denom_2=0;
  
  for(int i=m-1; i>=0; i--){
    denom_1+=exp_eta1(i);
    denom_2+=exp_eta2(i);
    R1_temp+=R1_j.row(i).t()*exp_eta1(i);
    R2_temp+=R2_j.row(i).t()*exp_eta2(i);
    R1_j_bar=R1_temp/denom_1;
    R2_j_bar=R2_temp/denom_2;
    
    z1_temp+=z1_j.row(i).t()*exp_eta1(i);
    z2_temp+=z2_j.row(i).t()*exp_eta2(i);
    z1_j_bar=z1_temp/denom_1;
    z2_j_bar=z2_temp/denom_2;
    z1_sq_temp+=z1_j.row(i).t()*z1_j.row(i)*exp_eta1(i);
    z2_sq_temp+=z2_j.row(i).t()*z2_j.row(i)*exp_eta2(i);
    
    z1_R1_temp+=z1_j.row(i).t()*R1_j.row(i)*exp_eta1(i);
    z2_R2_temp+=z2_j.row(i).t()*R2_j.row(i)*exp_eta2(i);
    
    
    if(delta1_j(i)!=0){
      r_inf+=(diagmat(R1_j_bar)-R1_j_bar*R1_j_bar.t());
      r_scr+=(R1_j.row(i).t()-R1_j_bar);
      beta1_scr+=(z1_j.row(i).t()-z1_j_bar);
      beta1_inf+=(z1_sq_temp/denom_1-z1_j_bar*z1_j_bar.t());
      beta1_r_inf+=(z1_R1_temp/denom_1-z1_j_bar*R1_j_bar.t());
    }
    if(delta2_j(i)!=0){
      r_inf+=(diagmat(R2_j_bar)-R2_j_bar*R2_j_bar.t());
      r_scr+=(R2_j.row(i).t()-R2_j_bar);
      beta2_scr+=(z2_j.row(i).t()-z2_j_bar);
      beta2_inf+=(z2_sq_temp/denom_2-z2_j_bar*z2_j_bar.t());
      beta2_r_inf+=(z2_R2_temp/denom_2-z2_j_bar*R2_j_bar.t());
    }
  }
}

// [[Rcpp::export]]
void beta_PPL_loop_huge(arma::vec& R1_j,arma::vec& R2_j, arma::mat& z1_j, arma::mat& z2_j,arma::vec& delta1_j,arma::vec& delta2_j,arma::vec& exp_eta1,
                        arma::vec& exp_eta2, arma::vec& beta1_scr, arma::vec& beta2_scr,arma::vec& r_scr,
                        arma::mat& beta1_inf, arma::mat& beta2_inf,arma::cube& r_inf,arma::mat& beta1_r_inf,arma::mat& beta2_r_inf,int m,int N,
                        arma::vec& eta1,arma::vec& eta2) {
  
  arma::mat z1_R1_temp(z1_j.n_cols,2*N);
  z1_R1_temp.fill(0);
  arma::mat z2_R2_temp(z2_j.n_cols,2*N);
  z2_R2_temp.fill(0);
  arma::mat z1_sq_temp(z1_j.n_cols,z1_j.n_cols);
  z1_sq_temp.fill(0);  
  arma::mat z2_sq_temp(z2_j.n_cols,z2_j.n_cols);
  z2_sq_temp.fill(0); 
  
  arma::vec z1_temp(z1_j.n_cols);
  z1_temp.fill(0);
  arma::vec z2_temp(z2_j.n_cols);
  z2_temp.fill(0);  
  
  arma::vec z1_j_bar(z1_j.n_cols);
  z1_j_bar.fill(0);
  arma::vec z2_j_bar(z2_j.n_cols);
  z2_j_bar.fill(0); 
  
  arma::vec R1_temp(2*N);
  R1_temp.fill(0);
  arma::vec R2_temp(2*N);
  R2_temp.fill(0);  
  
  arma::vec R1_j_bar(2*N);
  R1_j_bar.fill(0);
  arma::vec R2_j_bar(2*N);
  R2_j_bar.fill(0);  
  
  double denom_1=0;
  double denom_2=0;
  arma::mat temp1(2,2);
  temp1.fill(0);
  arma::mat temp2(2,2);
  temp2.fill(0);
  arma::mat temp3(2,2);
  temp3.fill(0);
  arma::mat temp4(2,2);
  temp4.fill(0);
  
  for(int i=m-1; i>=0; i--){
    denom_1+=exp_eta1(i);
    denom_2+=exp_eta2(i);
    R1_temp(2*R1_j(i)-2)+=exp_eta1(i);
    R2_temp(2*R2_j(i)-1)+=exp_eta2(i);
    R1_j_bar=R1_temp/denom_1;
    R2_j_bar=R2_temp/denom_2;
    
    for(int k=0;k<N;k++){
      temp1(0,0)=R1_j_bar(2*k);
      temp2(1,1)=R2_j_bar(2*k+1);
      temp3(0,0)=pow(R1_j_bar(2*k),2);
      temp4(1,1)=pow(R2_j_bar(2*k+1),2);
      r_inf.slice(k)+=(temp1-temp3)*delta1_j(i)+(temp2-temp4)*delta2_j(i); 
      r_scr(2*k)+=-R1_j_bar(2*k)*delta1_j(i);
      r_scr(2*k+1)+=-R2_j_bar(2*k+1)*delta2_j(i);
    }

    z1_temp+=z1_j.row(i).t()*exp_eta1(i);
    z2_temp+=z2_j.row(i).t()*exp_eta2(i);
    z1_j_bar=z1_temp/denom_1;
    z2_j_bar=z2_temp/denom_2;
    z1_sq_temp+=z1_j.row(i).t()*z1_j.row(i)*exp_eta1(i);
    z2_sq_temp+=z2_j.row(i).t()*z2_j.row(i)*exp_eta2(i);
    
    z1_R1_temp.col(2*R1_j(i)-2)+=z1_j.row(i).t()*exp_eta1(i);
    z2_R2_temp.col(2*R2_j(i)-1)+=z2_j.row(i).t()*exp_eta2(i);
    
    if(delta1_j(i)!=0) {
      r_scr(2*R1_j(i)-2)+=delta1_j(i);
      beta1_scr+=(z1_j.row(i).t()-z1_j_bar);
      beta1_inf+=(z1_sq_temp/denom_1-z1_j_bar*z1_j_bar.t());
      beta1_r_inf+=(z1_R1_temp/denom_1-z1_j_bar*R1_j_bar.t());
    }
    if(delta2_j(i)!=0){
      r_scr(2*R2_j(i)-1)+=delta2_j(i);
      beta2_scr+=(z2_j.row(i).t()-z2_j_bar);
      beta2_inf+=(z2_sq_temp/denom_2-z2_j_bar*z2_j_bar.t());
      beta2_r_inf+=(z2_R2_temp/denom_2-z2_j_bar*R2_j_bar.t());
    }
    

  }
}
