#ifndef __Maxent__
#define __Maxent__

#include "txtio.h"
#include "misc.h"
#include <exception>
#include <complex>
#include <string>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;


class Maxent{
 public:
  //----------Data declaration-------//
  clock_t start;

  const char * fname;
  const char * defaultModel;
  
  unsigned numMfre;
  unsigned numRfre;
  unsigned numAlpha;
  unsigned numBins;
  unsigned maxIteration;
  unsigned rank;
  

  double wmin;
  double wmax;
  double tol;
  double alphamin;
  double alphamax;
  double mu;
  double alpha;
  double prob;

  vector<vector<double> > rawData;
  VectorXd wn;
  VectorXd w;
  VectorXd dw;
  VectorXd stdG;
  VectorXd specF;
  VectorXd defaultM;
  VectorXd allProbs;
  VectorXd alphas;
  VectorXd dalpha;
  VectorXd aveSpecFs;
  VectorXcd aveG;
  MatrixXd allSpecFs;
  MatrixXcd G;
  MatrixXcd cov;
  MatrixXcd K;
  MatrixXcd U;
  MatrixXcd Ut;
  MatrixXcd V;
  MatrixXcd Vt;
  MatrixXd Xi;
  MatrixXcd M;
  MatrixXcd T;
  
  //-----------Function declaration-------//
  Maxent(const char*, const unsigned, const unsigned, const unsigned, const unsigned, const double, const double, const char *, const double, const double, const double, const unsigned, const double);
  ~Maxent();

  
  //private:
  void readfile(const char *, const unsigned, const unsigned, const unsigned);
  void calMeanAndStd();
  void createKernel();
  void rotation();
  void compUXiV();
  void initSpecF(const char *);
  void getSpecF();
  void calProb();
  void getAllSpecFs();
  void timeProcess(const char *);
  void saveData();

  double normalize(const VectorXd &);
  double objective(const VectorXd &);
  VectorXcd restoreG(const VectorXd &);
  
};

Maxent::Maxent(const char* fname_,
	       const unsigned columns_,
	       const unsigned rows_,
	       const unsigned numMfre_,
	       const unsigned numRfre_,
	       const double wmin_,
	       const double wmax_,
	       const char * defaultM_,
	       const double tol_,
	       const double alphamin_,
	       const double alphamax_,
	       const unsigned numAlpha_,
	       const double mu_){

  start = clock();
  fname = fname_;
  numMfre = numMfre_;
  numRfre = numRfre_;
  numBins = columns_/2;
  wmin = wmin_;
  wmax = wmax_;
  defaultModel = defaultM_;
  tol = tol_;
  alphamin = alphamin_;
  alphamax = alphamax_;
  numAlpha = numAlpha_;
  mu = mu_;
  maxIteration = 10000;

  readfile(fname_, columns_, rows_, numMfre_);

  calMeanAndStd();

  createKernel();

  rotation();

  compUXiV();

  initSpecF(defaultM_);
}

void Maxent::readfile(const char* fname_,
		      const unsigned columns_,
		      const unsigned rows_,
		      const unsigned numMfre_){
  
  if (columns_/2 < 2*numMfre_){
    std::cerr << "Error: number of samples must be larger than or equal to 2 times of number of Matsubara frequency.\n";
    std::cerr << "Exit program.\n";
    throw std::exception();
  }

  complex<double> imag(0, 1);
    
  loadtxt(fname_ ,rawData, rows_, columns_, 0, "\t");
  
  wn.resize(numMfre_);
  size_t begin = rows_/2 - numMfre_/2;
  for (size_t i = 0; i < numMfre_; ++i) {
    wn(i) = rawData[begin + i][0];
  }
  
  G.resize(numMfre_, columns_/2);
  for (size_t i = 0; i < numMfre_; ++i) {
    for (size_t j = 0; j < columns_/2; ++j) {
      G(i,j)  = rawData[begin+i][2*j+1] + imag*rawData[begin+i][2*j+2];
    }
  }
  
}

Maxent::~Maxent(){};

void Maxent::calMeanAndStd(){
  
  aveG.resize(numMfre);
  aveG = G.rowwise().sum()/numBins;

  cov.resize(numMfre, numMfre);
  for (size_t l = 0; l < numMfre; ++l){
    for (size_t k = 0; k < numMfre; ++k){
      complex<double> a(0, 0);
      for (size_t j = 0; j < numBins; ++j){
	a += (aveG(l) - G(l, j)) * std::conj( aveG(k) - G(k, j) );
      }
      cov(l, k) = a/double(numBins - 1);
    }
  }

}

void Maxent::createKernel(){

  complex<double> imag(0, 1);

  double ddw = (wmax - wmin)/double(numRfre -1);
  w.resize(numRfre);
  for (size_t i = 0; i < numRfre; ++i){
    w(i) = wmin + ddw * double(i);
  }

  dw.resize(numRfre);
  for (size_t i = 0; i < numRfre - 1; ++i) dw(i) = w(i+1) - w(i);
  dw(numRfre-1) = dw(0)/2.0;
  dw(0) = dw(0)/2.0;

  K.resize(numMfre, numRfre);
  for (size_t n = 0; n < numMfre; ++n){
    for (size_t m = 0; m < numRfre; ++m){
      K(n, m) = 1.0/(-w(m) + wn(n) * imag);
    }
  }
}

void Maxent::rotation(){
  
  SelfAdjointEigenSolver<MatrixXcd> ces;
  ces.compute(cov);

  VectorXd w = ces.eigenvalues();
  MatrixXcd v = -ces.eigenvectors();
  MatrixXcd vt = v.adjoint();

  stdG.resize(numMfre);
  for (size_t i = 0; i < numMfre; ++i){
    stdG(i) =  sqrt(w(i));
  }

  K = vt*K;
  aveG = vt*aveG;

}

void Maxent::compUXiV(){

  JacobiSVD<MatrixXcd> svd(K.adjoint(), ComputeThinU | ComputeThinV);

  VectorXd singularValues = svd.singularValues();
  rank = svd.rank();
  U = svd.matrixU();
  V = svd.matrixV();

  Xi.resize(rank, rank);
  for (size_t i = 0; i < rank; ++i){
    for (size_t j = 0; j < rank; ++j){
      if (i == j){
	Xi(i, j) = singularValues(i);
      }else{
	Xi(i, j) = 0.0;
      }
    }
  }

  U = U.leftCols(rank);
  Ut = U.adjoint();
  V = V.leftCols(rank);
  Vt = V.adjoint();

  MatrixXd a = MatrixXd::Zero(numMfre, numMfre);
  for (size_t i = 0; i < numMfre; ++i) a(i, i) = 1.0/stdG(i)/stdG(i);

  M = Xi * Vt * a * V * Xi;
  
}

void Maxent::initSpecF(const char * defaultM_){

  std::string str(defaultM_);
  if (str == "gaussian"){
    specF = gaussian(this->w)/normalize(gaussian(this->w));
    defaultM = specF;
  }
  else if (str == "straightline"){
    specF = straightline(this->w);
    defaultM = specF;
  }
  else{
    std::ifstream file(defaultM_);
    if (file.is_open()){
      double a;
      defaultM.resize(numRfre);
      for (size_t i = 0; i < numRfre; ++i){	
	file >> a >> defaultM(i);
      }
      specF = defaultM;
    }else{
      cerr <<  "Didn't find the default Model file!";
      throw std::exception();
    }

  }

  
  clog << "Use --- " << (clock() - start)/double(CLOCKS_PER_SEC) << " seconds --- before optimization...\n";

  
}

double Maxent::normalize(const VectorXd & foo){
  return foo.dot(this->dw);
}

VectorXcd Maxent::restoreG(const VectorXd & foo_){
  VectorXd foo;
  foo.resize(foo_.size());
  for (size_t i = 0; i < unsigned(foo_.size()); ++i) foo(i) = foo_(i) * dw(i);

  return K * foo;
}

double Maxent::objective(const VectorXd & foo_){
  VectorXcd delta;
  VectorXd foo;
  double chi = 0.0;
  double entropy = 0.0;
  foo.resize(foo_.size());
  for (size_t i = 0; i < unsigned(foo_.size()); ++i) foo(i) = foo_(i) * dw(i);

  delta = aveG - K * foo;

  for (size_t i = 0; i < numMfre; ++i){
    chi += real(conj(delta(i)) * delta(i)/stdG(i)/stdG(i));
  }
  chi /= 2.0;
 
  for (size_t i = 0; i < numMfre; ++i){
    entropy += (-foo_(i) + defaultM(i) + foo_(i) * std::log(foo_(i)/defaultM(i))) * dw(i);
  }

  return chi + alpha * entropy;
}

void Maxent::getSpecF(){
  
  unsigned iteration = 0;
  double Qold = objective(defaultM);
  double Qnew;
  double criteria;
  double metric = 0.2 * defaultM.sum();
  
  VectorXcd btemp = VectorXcd::Constant(rank, 0.0);
  VectorXcd deri = VectorXcd::Zero(numMfre);
  VectorXcd g;
  VectorXcd temp;
  VectorXcd RHS;
  MatrixXcd LHS;
  VectorXcd deltab;
  VectorXcd al;

  while (true){
    ++iteration;
    MatrixXd A = MatrixXd::Zero(numRfre, numRfre);
    
    for (size_t i = 0; i < numRfre; ++i) A(i, i) = specF(i);
   
    T.noalias() = Ut*A*U;
    /* for (size_t j = 0; j < rank; ++j){ */
    /*   for (size_t i = 0; i < rank; ++i){ */
    /* 	complex<double> sum(0.0, 0.0); */
    /* 	for (size_t k = 0; k < numRfre; ++k){ */
    /* 	  sum += Ut(i, k) * specF(k) * U(k,j); */
    /* 	} */
    /* 	T(i,j) = sum; */
    /*   } */
    /* } */
    
    
    temp = aveG - restoreG(this->specF);

    for (size_t i = 0; i < numMfre; ++i){
      deri(i) = -1.0/stdG(i)/stdG(i) * temp(i);
    }
    
    g = Xi * Vt * deri;
    
    LHS = (alpha + mu) * MatrixXcd::Identity(rank, rank) + M * T;
    
    RHS = -alpha * btemp - g;
    
    deltab = LHS.colPivHouseholderQr().solve(RHS);
   
    criteria = real(deltab.dot(T * deltab));
    
    if (criteria < metric){
      if (iteration > maxIteration){
	cerr << "Exceeds maximum iteration in Levenberg-Marquart algorithms, exits. Make tolerance smaller.\n";
	break;
      }
      btemp = btemp + deltab;
      al = U * btemp;
      for (size_t i = 0; i < numRfre; ++i) {
	specF(i) = real(defaultM(i) * std::exp(al(i)));
      }
      Qnew = objective(specF);
      //cout << "iteration: " << iteration << endl;
      if (abs(Qnew - Qold)/Qold < tol){
	clog << iteration << " iterations in Levenberg-Marquart algorithems. Function evaluated: " << Qnew << ", it exits properly.\n";
	break;
      }
      Qold = Qnew;
      continue;
    }
    else {

      mu *= 2.0;
      specF = defaultM;
      Qold = objective(defaultM);
      btemp = VectorXcd::Zero(rank);
      cerr << "parameter mu is too small in the Levenberg-Marquart algorithms.\n";
      cerr << "mu is now adjusted to " << mu << ".\n";
    }
    
  }
}

void Maxent::timeProcess(const char * mark){
  clog << mark << ": ";
  clog << (clock() - start)/(double(CLOCKS_PER_SEC)/1000.0) << "\n";
  start = clock();
}

void Maxent::calProb(){

  double product = 1.0;
  MatrixXd C_1 = MatrixXd::Zero(numMfre, numMfre);
  MatrixXcd mat_a;
  MatrixXcd mat_b;
  MatrixXd vec_a = MatrixXd::Zero(numRfre, numRfre);

  for (size_t i = 0; i < numMfre; ++i){
    C_1(i, i) = 1.0/stdG(i)/stdG(i);
  }
  mat_a.noalias()  = K.adjoint() * C_1 * K;
  for (size_t i = 0; i < numRfre; ++i){
    vec_a(i, i) = std::sqrt(specF(i));
  }
  mat_b.noalias() = vec_a * mat_a * vec_a;

  SelfAdjointEigenSolver<MatrixXcd> ces;
  ces.compute(mat_b);
  VectorXd S = ces.eigenvalues();

  double expo = std::exp(-objective(specF));

  for (size_t i = 0; i < numRfre; ++i){
    product *= alpha/(alpha + S(i));
  }

  prob = std::sqrt(product) * expo / alpha;

  if (std::isnan(prob)) prob = 0.0;
  
}

void Maxent::getAllSpecFs(){

  allSpecFs.resize(numRfre, numAlpha);
  allProbs.resize(numAlpha);

  alphas = logspace(alphamin, alphamax, numAlpha);

  dalpha = VectorXd::Zero(numAlpha);
  
  VectorXd dalpha_ = VectorXd::Zero(numAlpha);
  for (size_t i = 0; i < numAlpha - 1; ++i){
    dalpha(i) = alphas(i+1) - alphas(i);
  }

  dalpha /= 2.0;
  for (size_t i = 0; i < numAlpha - 1; ++i){
    dalpha_(i+1) = dalpha(i);
  }
  for (size_t i = 0; i < numAlpha; ++i){
    dalpha(i) += dalpha_(i);
  }


  for (size_t i = 0; i < numAlpha; ++i){

    alpha = alphas(i);
    getSpecF();
    calProb();
    allProbs(i) = prob;
    allSpecFs.col(i) = specF;
    clog << "Finish alpha at " << alpha << ".\n";
  }
  VectorXd temp(numAlpha);
  double sum;
  for (size_t i = 0; i < numAlpha; ++i) temp(i) = allProbs(i) * dalpha(i);
  sum = temp.sum();
  for (size_t i = 0; i < numAlpha; ++i) allProbs(i) /= sum;
  
  aveSpecFs = allSpecFs * allProbs;

  clog << "Optimization ends. Use --- " << (clock() - start)/double(CLOCKS_PER_SEC) << " --- seconds.\n";
  
}

void Maxent::saveData(){

  std::string dirname(fname);
  dirname = forestr(dirname, ".");

  savetxt((dirname+"maxent.dat").c_str(), w , aveSpecFs);
  savetxt((dirname +"Palpha.dat").c_str(), alphas, allProbs);
  
}

#endif


