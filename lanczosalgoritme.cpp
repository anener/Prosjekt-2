#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

#include "lanczo.h"

void identitetMatrise(int n, mat& I) { //lager en identitets matrise
    for(int i=0; i<n; i++) {
        I(i, i) = 1.0;
    }
}

void lanczoMetode(int n, mat A, vec& b, vec& a) { //kjorer Lanczo Metoden
    mat I = zeros<mat>(n, n); //identitets matrisen
    mat Q = zeros<mat>(n, n); //ortohonale matrisen
    mat R = zeros<mat>(n, n);
    identitetMatrise(n, I);

    vec q = Q.col(0);
    vec r = R.col(0);
    r(0) = 1.0;
    b(0) = 1.0;
    a(0) = 0.0;

    int k = 0;
    while( (b(k) =! 0 && k < n-1) ) {
        vec q_k = Q.col(k+1);
        q_k = r/b(k);

        k = k+1;
        cout << q_k << "\n";

        vec aa = q_k.t()*A*q_k;
        a(k) = aa(0);
        cout << a(k) << "\n";
        r = R.col(k);
        r = (A - a(k)*I)*q_k - b(k-1)*q;

        vec rr = sqrt(r.t()*r);
        b(k) = rr(0);
        q = q_k;
    }
}
