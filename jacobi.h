#ifndef JACOBI_H
#define JACOBI_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace std;
using namespace arma;

void vektor(int, mat&);
void rhoVektor(int, vec&, double, double);
void matrise(int, mat&);
void matriseHOP(int, mat&, vec, double, double);
double maxNonDiagonal(mat, int, int*, int*);
void trigonometri(mat, int, int, double&, double&);
void jacobiRotasjon(mat&, mat&, int, int, int, double, double);

#endif // JACOBI_H
