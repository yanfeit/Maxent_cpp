#ifndef __misc__
#define __misc__

#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <ctime>
#include <string>
#include <iostream>
#include <cmath>


//----------function declaration---------

Eigen::VectorXd gaussian(const Eigen::VectorXd &, double, double);

Eigen::VectorXd straightline(const Eigen::VectorXd &);


//----------function definition----------

Eigen::VectorXd gaussian(const Eigen::VectorXd & w_, double mu_ = 0.0, double sig_ = 4.0){
  Eigen::VectorXd y;
  y.resize(w_.size());

  for (size_t i = 0; i < unsigned(w_.size()); ++i){
    y(i) = 1.0/sig_/std::sqrt(2 * M_PI) * std::exp(-std::pow(w_(i) - mu_, 2.0)/2.0/sig_/sig_);
  }
  return y;
  
}

Eigen::VectorXd straightline(const Eigen::VectorXd & w_){
  Eigen::VectorXd y = Eigen::VectorXd::Constant(w_.size(), 1.0);
  y = y/(w_(w_.size() - 1) - w_(0));
  return y;
}

Eigen::VectorXd linspace(const double & start, const double & stop, const unsigned & num){

  Eigen::VectorXd y(num);
  double step = (stop - start)/double(num - 1);
  for (size_t i = 0; i < num; ++i){
    y(i) = start + step * double(i);
  }
  return y;
}

Eigen::VectorXd logspace(const double & start, const double & stop, const unsigned & num){

  Eigen::VectorXd y(num);
  Eigen::VectorXd z(num);
  y = linspace(start, stop, num);
  for (size_t i = 0; i < num; ++i){
    z(i) = std::pow(10.0 ,y(i));
  }
  return z;
}


std::string forestr(std::string const& s, std::string const &mark)
{
    std::string::size_type pos = s.find_last_of(mark);
    if (pos != std::string::npos)
    {
        return s.substr(0, pos);
    }
    else
    {
        return s;
    }
}

#endif
