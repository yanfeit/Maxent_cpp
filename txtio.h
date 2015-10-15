#ifndef __txtio__
#define __txtio__

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>


//-------declartion-------

void loadtxt(const char* , std::vector<std::vector<double> > & , const int , const int , const int , const std::string);

template <typename T>
void savetxt(const char* , const T&, const T&, const T&);

template <typename T>
void savetxt(const char* , const T&, const T&);

template <typename T>
void savetxt(const char* , const T&);

//-------definition--------

void loadtxt(const char* fname, std::vector<std::vector<double> > & data, const int rows, const int cols, const int skiprows = 0, const std::string delimiter = " "){

  std::ifstream file(fname);
  
  if (file.good()){
    data.resize(rows, std::vector<double>(cols));
    std::string str;
    int line = 0;
    int row, col;
    while (std::getline(file, str)){
    
      if (line >= skiprows){
	row = line-skiprows;
	col = 0;
	
	size_t pos = 0;
	std::string token;
	while ((pos = str.find(delimiter)) != std::string::npos) {
	  
	  token = str.substr(0, pos);
	  data[row][col] = atof(token.c_str());
	  str.erase(0, pos + delimiter.length());
	  ++col;
	}
        data[row][col] = atof(str.c_str());
      
      }
      line++;
    }
  }
  else{
    std::cout << "Usage: Loadtxt(const char* fname, vector<vector<double> > & data,"
	      << "const int rows, const int cols, const int skiprows = 0,"
	      << "const std::string delimiter = \" \")" << std::endl;
    std::cout << "Can't open file." << std::endl;
    exit(-1);
  }
  
}

template <typename T>
void savetxt(const char* fname,const T& xaxis,const T& yaxis, const T& yerror){
  
  std::ofstream file(fname, std::ofstream::out);
  
  if (file.good()) {
    for (size_t i = 0; i < xaxis.size(); ++i){
      file << xaxis[i] << "\t" << yaxis[i] << "\t" << yerror[i] << "\n";
    }
    
  }
  
}

template <typename T>
void savetxt(const char* fname,const T& xaxis,const T& yaxis){
  
  std::ofstream file(fname, std::ofstream::out);
  
  if (file.good()) {
    for (size_t i = 0; i < xaxis.size(); ++i){
      file << xaxis[i] << "\t" << yaxis[i]  << "\n";
    }
    
  }
  
}

template <typename T>
void savetxt(const char* fname,const T& xaxis){
  
  std::ofstream file(fname, std::ofstream::out);
  
  if (file.good()) {
    for (size_t i = 0; i < xaxis.size(); ++i){
      file << xaxis[i] << "\n";
    }
    
  }
  
}

void savetxt(const char * fname, const Eigen::VectorXd & xaxis, const Eigen::VectorXd & yaxis){

  std::ofstream file(fname, std::ofstream::out);
  
  if (file.good()) {
    for (size_t i = 0; i < unsigned(xaxis.size()); ++i){
      file << xaxis(i) << "\t" << yaxis(i)  << "\n";
    }
    
  }
  
}

#endif
