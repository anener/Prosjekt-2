#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <time.h>
#include <algorithm>
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "jacobi.h"
#include "bisectionmetode.h"
#include "lanczo.h"

vec analytiskEigen(int);
vector<double> sortEigenValue(mat, int);

int main()
{

    clock_t start, slutt; //inisialiserer klokken for aa kunne ta tiden de forskjellige metodene tar
    int n = 5; double max = 1.0; double min = 0.0; //setter size paa matrisen (nxn), rho_max og rho_min
    mat A = zeros<mat>(n, n); mat R = zeros<mat>(n, n); //lager to matriser, A er Toeplitz matrisen og R er matrisen som skal oppbevare eigenvektorene
    matrise(n, A); vektor(n, R); //lager matrisene A og R

    //Eig_sym
    start = clock(); //starter klokken
    vec eigen = eig_sym(A); //Armadillos interne funksjon finner eigenverdiene til A
    slutt = clock(); //stopper klokken
    double eigen_tid = double(slutt-start)/CLOCKS_PER_SEC; //finner tiden det tok
    cout << "Eig_sym: \n";
    cout << "tid = " << eigen_tid << "\n"; //skriver ut tid
    cout << "eigenverdiene = [" << eigen(0) << ", " << eigen(1) << ", " << eigen(2) << "] \n"; //skriver ut de tre forste eigenverdiene

    //Analytisk
    start = clock(); //starter klokken
    vec analyse = analytiskEigen(n+1); //regner ut de analytiske eigenverdiene
    slutt = clock(); //stopper klokken
    double analyse_tid = double(slutt-start)/CLOCKS_PER_SEC; //finner tiden det tok
    cout << "Analytisk: \n";
    cout << "tid = " << analyse_tid << "\n"; //skriver ut tid
    cout << "eigenverdiene = [" << analyse(1) << ", " << analyse(2) << ", " << analyse(3) << "] \n"; //skriver ut de tre forste eigenverdiene

    //Jacobi metode
    double tol = 1.0E-10; double m = 1; //setter en toleranse (tol) og en variabel m til en midlertidig verdi. m er maks verdi i A (ikke diagonal element)
    int iter = 0; int maxIter = 500000; //setter start iteratoren (iter) og maksimum antal itetasjoner (maxIter)
    start = clock(); //starter klokken
    while(iter < maxIter && m > tol) { //kjorer loopen sa lenge vi ikke har nadd maksimum antall iterasjoner og maks verdi til ikke diagonale elementet i A er storre enn toleransen
        int p, q; //setter variablene for plaseringen til maks verdien i A (ikke diagonal element)
        m = maxNonDiagonal(A, n, &p, &q); //finner maks verdien i A (ikke diagonal element) og finner posisjonen til den
        double c = 0; double s = 0; //setter cosinus og sinus til en midlertidig verdi
        trigonometri(A, p, q, c, s); //regner ut cosinus og sinus i forhold til posisjonene til m
        jacobiRotasjon(A, R, p, q, n, c, s); //utforer Jacobis rotasjons metode
        iter++; //oker iteratoren med 1
    }
    slutt = clock(); //stopper klokken
    double jacobi_tid = double(slutt-start)/CLOCKS_PER_SEC; //finner tiden det tok
    cout << "Jacobi: \n";
    cout << "tid = " << jacobi_tid << "\n"; //skriver ut tiden
    cout << "Antall transformer = " << iter << "\n"; //skriver ut antall iterasjoner det trenktes
    cout << "Storrelse paa matrise = " << n << "\n"; //skriver ut storrelsen til matrisen

    vector<double> d = sortEigenValue(A, n); //henter ut eigenverdiene til A etter diagonaliseringen (Jacobi metode) og sorterer de
    cout << "eigenverdiene = [" << d[0] << ", " << d[1] << ", " << d[2] << "] \n"; //skriver ut de tre forste eigenverdiene


    /*
    //prosjekt 2d
    int x = 250; //definerer en N
    mat H = zeros<mat>(x, x); //lager Toeplitz matrisen med det harmoniske potensialet
    mat T = zeros<mat>(x, x); //lager en matrise hvor eigenvektorene lagres
    vec rho(x); //definerer rho vektoren
    double min = 0; double max = 10; //setter rho_min og rho_max
    rhoVektor(x, rho, min, max); //lager rho verdiene
    matriseHOP(x, H, rho, max, min); //lager Toeplitz matrisen
    vektor(x, T); //lager vektoren hvor eigenvektorene lagres

    //resten av koden fungerer paa samme maate som for prosjekt2b
    double Tol = 1.0E-10; double M = 1;
    int Iter = 0; int MaxIter = 50000;
    while(Iter < MaxIter && M > Tol) {
        int p, q;;
        M = maxNonDiagonal(H, x, &p, &q);
        double c = 0; double s = 0;
        trigonometri(H, p, q, c, s);
        jacobiRotasjon(H, T, p, q, x, c, s);
        Iter++;
    }
    vector<double> d = sortEigenValue(H, x);
    cout << d[0] << ", " << d[1] << ", " << d[2] << ", " << d[3] << "\n";
    */

    /*
    //prosjekt 2h
    clock_t start, slutt;
    int n = 5;
    double min = 0; double max = 1;
    mat A = zeros<mat>(n, n); //lager en matrise
    vec x(n);
    rhoVektor(n, x, min, max); //setter rho verdiene
    matriseHOP(n, A, x, max, min); //lager en Toeplitz matrise med det harmoniske potensialet
    vec b(n); vec a(n); //lager to vektorer som skal holde på a og b verdiene i T=QAQ^T
    start = clock();
    lanczoMetode(n, A, b, a); //sender Toeplitz matrisen pluss de to vektorene (a og b) inn i lanczoMetoden for aa kjore metoden
    slutt = clock();
    double tid = double(slutt-start)/CLOCKS_PER_SEC; //finner tiden det tok

    return 0;
    */
}

vec analytiskEigen(int n) { //en funksjon som regner ut de analytiske eigenverdiene
    vec eigen(n);
    double h = 1.0/(n);
    double d = 2.0/(h*h);
    double a = -1.0/(h*h);

    for(int i=0; i<n; i++) {
        eigen(i) = d + 2*a*cos( (i*M_PI)/(n+1) );
    }

    return eigen;
}

vector<double> sortEigenValue(mat A, int n) { //tar imot en matrise A, legger diagonal verdiene inn i en vektor og sorterer vektoren
    vector<double> eigen;
    for(int i=0; i<n; i++) {
        eigen.push_back(A(i, i));
    }
    sort(eigen.begin(), eigen.end());
    return eigen;
}


TEST_CASE("Tester maks A(i, j)") { //tester om koden finner den korekte max ikke diagonale verdien i A
    int n = 3;
    mat A = zeros<mat>(n, n);
    mat R = zeros<mat>(n, n);
    matrise(n, A);
    int p, q;
    double m = maxNonDiagonal(A, n, &p, &q);

    REQUIRE(p == 2);
    REQUIRE(q == 1);
    REQUIRE(m == Approx(16));
}


TEST_CASE("Tester eigenverdier") { //tester om utregningen av eigenverdiene stemmer
    int n = 5;
    mat A = zeros<mat>(n, n);
    mat R = zeros<mat>(n, n);
    matrise(n, A);
    vektor(n, R);
    double tol = 1.0E-10; double m = 1;
    int iter = 0; int maxIter = 50;
    while(iter < maxIter && m > tol) {
        int p, q;
        m = maxNonDiagonal(A, n, &p, &q);
        double c = 0; double s = 0;
        trigonometri(A, p, q, c, s);
        jacobiRotasjon(A, R, p, q, n, c, s);
        iter++;
    }
    vector<double> D = sortEigenValue(A, n);

    REQUIRE(D[0] == Approx(6.69873));
    REQUIRE(D[1] == Approx(25));
    REQUIRE(D[2] == Approx(50));

}

