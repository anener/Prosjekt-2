#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

#include "jacobi.h"


void vektor(int n, mat& R) { //lager en matrise med 1ere langs diagonalen. Denne brukes for aa lage matrisen som skal romme eigenvektorene
    for(int i=0; i<n; i++) {
        R(i, i) = 1.0;
    }
}

void rhoVektor(int n, vec& rho, double min, double max) { //lager en vektor som holder alle rho verdiene
    rho(0) = min;
    double h = (max-min)/(n);
    for(int i=1; i<n; i++) {
        rho(i) = rho(0) + i*h; //rho_i = rho_0 + i*h
    }
}

void matrise(int n, mat& A) { //lager en tridiagonal Toeplitz matrise
    double h = 1.0/(n);
    double d = 2.0/(h*h);
    double a = -1.0/(h*h);

    A(0, 0) = d; A(0, 1) = a; A(1, 0) = a; //setter d_1 og a_1
    for(int i=1; i<(n-1); i++) {
        A(i, i-1) = a;
        A(i, i) = d;
        A(i-1, i) = a;
    }
    A(n-1, n-2) = a; A(n-1, n-1) = d; A(n-2, n-1) = a; //setter d_n og a_n-1
}

void matriseHOP(int n, mat& H, vec rho, double max, double min) { //lager en tridiagonal Toeplitz matrise med det harmoniske osillatorpotensialet
    double h = (max-min)/(n);
    double d = 2.0/(h*h);
    double e = -1.0/(h*h);

    H(0, 0) = d + (rho(0)*rho(0)); H(0, 1) = e; H(1, 0) = e; //setter d_1= d + V(rho), og e_1
    for(int i=1; i<(n-1); i++) {
        H(i, i-1) = e;
        H(i, i) = d + (rho(i)*rho(i));
        H(i-1, i) = e;
    }
    H(n-1, n-1) = d + (rho(n-1)*rho(n-1)); H(n-2, n-1) = e; H(n-1, n-2) = e; //setter d_n og e_n-1
}

double maxNonDiagonal(mat A, int n, int* p, int* q) { //finner det ikke diagonale elementet som har storrst verdi
    double m = 0.0;
    for(int i=0; i < n; i++) {
        for(int j=0; j < n; j++) {
            if(i != j && fabs(A(i, j)) >= m) {
                m = fabs(A(i, j)); //setter m til verdien for det ikke diagonale elementet som har storst verdi
                *p = i; *q = j; //setter posisjonen til m
            }
        }
    }
    return m;
}


void trigonometri(mat A, int k, int l, double& c, double& s) { //finner cosinus og sinus fra posisjonen til m=A(k, l)
    if (A(k, l) == 0.0) {
        c = 1.0; s = 0.0;
    }
    else {
        double t; //tangent
        double tau; //cot
        tau = (A(l,l) - A(k, k))/(2*A(k, l));
        if (tau < 0.0) {
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }
        else {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        c = 1.0/sqrt(1+t*t); s = c*t;
    }
}

void jacobiRotasjon(mat& A, mat& R, int k, int l, int n, double c, double s) { //utforer Jacobis rotasjons metode
    double a_ik, a_il, a_kk, a_ll;
    a_kk = A(k, k);
    a_ll = A(l, l);

    A(k, k) = a_kk*(c*c) - 2*A(k, l)*c*s + a_ll*(s*s);
    A(l, l) = a_ll*(c*c) + 2*A(k, l)*c*s + a_kk*(s*s);
    A(k, l) = 0.0;
    A(l, k) = 0.0;

    double r_ik, r_il;

    for(int i=0; i < n; i++) {
        if(i != k && i != l) {
            a_ik = A(i, k);
            a_il = A(i, l);
            A(i, k) = a_ik*c - a_il*s;
            A(i, l) = a_il*c + a_ik*s;

            A(k, i) = A(i, k);
            A(l, i) = A(i, l);
        }
        r_ik = R(i, k); r_il = R(i, l);
        R(i, k) = c*r_ik - s*r_il;
        R(i, l) = c*r_il + s*r_ik;
    }
}
