#ifndef __statistics__
#define __statistics__
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm> //for_each
#include <numeric> //accumulate, partial_sum, reduce
#include <vector>

double mean(std::vector<double>);
double variance(std::vector<double>);
double sd(std::vector<double>);
double sd_mean(std::vector<double>);
void autocorrelation(std::vector<double>, int, std::string); 
std::vector<double> blocking(std::vector<double> ist, int Nblocks, int Nsteps);
std::vector<double> blocking_lento(std::vector<double> ist, int Nblocks, int Nsteps);
std::vector<double> running_av(std::vector<double>);
std::vector<double> running_err(std::vector<double>);
std::vector<double> bm_running_err(std::vector<double>);
bool compatible(std::vector<double>, std::vector<double>);

#endif 
