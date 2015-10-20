#include "Maxent.h"
#include <iostream>
#include <cstdlib>
using namespace std;

int main(int argc, char * argv[]){

  int i = 0;
  char * fname_ = NULL;
  unsigned columns_ = 50;
  unsigned rows_ = 50;
  unsigned numMfre_ =50;
  unsigned numRfre_ =50;
  double wmin_ = -15.0;
  double wmax_ = 15.0;
  char* defaultM_ = NULL;
  double tol_ = 1e-5;
  double alphamin_ = -2.0;
  double alphamax_ = 1.0;
  unsigned numAlpha_ = 10;
  double mu_ = 1e3;
  std::clog << "********************Maximum Entropy Method ****************\n";
  std::clog << "**                                                       **\n";
  std::clog << "**         Copyright Yanfei Tang, 10.13.2015             **\n";
  std::clog << "***********************************************************\n";
  std::clog << "\n";
  while (++i < argc){

    std::string str(argv[i]);
    if (str == "-file" && i < argc-1) fname_ = argv[++i];
    if (str == "-col" && i < argc-1) columns_ = atoi(argv[++i]);
    if (str == "-row" && i < argc-1) rows_ = atoi(argv[++i]);
    if (str == "-numMfre" && i < argc-1) numMfre_ = atoi(argv[++i]);
    if (str == "-numRfre" && i < argc-1) numRfre_ = atoi(argv[++i]);
    if (str == "-wmin" && i < argc-1) wmin_ = atof(argv[++i]);
    if (str == "-wmax" && i < argc-1) wmax_ = atof(argv[++i]);
    if (str == "-defaultM" && i< argc-1) defaultM_ = argv[++i];
    if (str == "-tol" && i < argc-1) tol_ = atof(argv[++i]);
    if (str == "-alphamin" && i < argc-1) alphamin_ = atof(argv[++i]);
    if (str == "-alphamax" && i < argc-1) alphamax_ = atof(argv[++i]);
    if (str == "-numAlpha" && i < argc-1) numAlpha_ = atoi(argv[++i]);
    if (str == "-mu" && i < argc-1) mu_ = atof(argv[++i]);
    if (str == "-h" || str == "--help"){
      std::clog << "Maxent [Option]\n";
      std::clog << "Option: -file     char *   Input File Nam\n";
      std::clog << "        -col      unsigned Number of Columns in Input File\n";
      std::clog << "        -row      unsigned Number of Rows in Input File\n";
      std::clog << "        -numMfre  unsigned Number of Matsubara Frequency Used\n";
      std::clog << "        -numRfre  unsigned Size of Gird for Real Frequency\n";
      std::clog << "        -wmin     double   Minimum Real Frequency\n";
      std::clog << "        -wmax     double   Maximum Real Frequency\n";
      std::clog << "        -defaultM char *   Name of Default Model(Can be a File)\n";
      std::clog << "        -tol      double   Tolerance, e.g. 1e-5\n";
      std::clog << "        -alphamin double   Minimum alpha\n";
      std::clog << "        -alphamax double   Maximum alpha\n";
      std::clog << "        -numAlpha unsigned Number of alpha. Alpha is Log Spaced\n";
      std::clog << "        -mu       double   Initial mu in The ML Algorithm\n";
      return 0;
    }
    
  }

  std::clog << "-file" << fname_ << "\n";
  std::clog << "-col" << columns_ << "\n";
  std::clog << "-row" << rows_ << "\n";
  std::clog << "-numMfre" << numMfre_ << "\n";
  std::clog << "-numRfre" << numRfre_ << "\n";
  std::clog << "-wmin" << wmin_ << "\n";
  std::clog << "-wmax" << wmax_ << "\n";
  std::clog << "-defaultM" << defaultM_ << "\n";
  std::clog << "-tol" << tol_ << "\n";
  std::clog << "-alphamin " << alphamin_ << "\n";
  std::clog << "-alphamax " << alphamax_ << "\n";
  std::clog << "-numAlpha " << numAlpha_ << "\n";
  std::clog << "-mu " << mu_ << "\n";
  
  try{
    Maxent model(fname_, columns_, rows_, numMfre_, numRfre_, wmin_, wmax_, defaultM_, tol_, alphamin_, alphamax_, numAlpha_, mu_);
    model.getAllSpecFs();
    model.saveData();

    
  }
  catch (const std::exception&){
    return EXIT_FAILURE;
  }
  

  
  return 0;
} 
