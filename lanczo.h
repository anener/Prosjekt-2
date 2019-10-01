#ifndef LANCZO_H
#define LANCZO_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace std;
using namespace arma;

void identitetMatrise(int, mat&);
void lanczoMetode(int,mat, vec&, vec&);

#endif // LANCZO_H
