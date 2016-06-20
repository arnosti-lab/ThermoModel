#ifndef THERMOMODEL_H
#define THERMOMODEL_H
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

using namespace std;

class Params {
    public:
        int n;
        int nOverlap;
        int* typeBS;
        int* distances;
        double dl;
        double tw;
        double sn;
        double* BAs;
        int* overlapDistances;
        double* Scalings;
        double* Cs;
        double* Qs;
        int typeInteractions;
        bool** allStates;
        int numStates;
        void SetParameters(int);
        void SetConcentrations(double, double, double);
        void saveStates();
        void States(int, bool*, bool*);
        void addState(bool*);
        double expression();
};

double quench(int, int, double, double *);
double binnedQuench(int, int, double *);
double coop(int, int, double, double *);
double binnedCoop(int, int, double *);
double ProteinbinnedCoop(int, int, double *);
void contributionsAllStates(Params const &, double *);
void getContributions(bool const * const, Params const &, double *);
void NoNcontributions(bool const * const, Params const &, double *);
void NearNcontributions(bool const * const, Params const &, double *);
void allNcontributions(bool const * const, Params const &, double *);
double SSE(double *, int);
double runobjective(double const * const, int, int, int); 
extern Params* thermo_parameters;
extern double Concentrations[3][17];
extern double Expressions[59][17];
#endif
