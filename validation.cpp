// Example usage and validation template

#include <iostream>
#include <vector>
#include <string>

#include "tsplineVol.cpp"
#include "nurbsVol.cpp"
#include "nurbsSur.cpp"
#include "nurbsCur.cpp"



using namespace std;

int main()
{
    vector<int> O {0, 0, 0};
    vector<int> X {1, 0, 0};
    vector<int> Y {0, 1, 0};
    vector<int> Z {0, 0, 1};
    double r = 1.;
    double ths = 0.;
    double the = 360.;

    
    

    vector<vector<double>> pt { {1., 0}, 
                                {.8, 0}, 
                                {.5, 0}};

    vector<double> ut {0., 0., 0.5, 1., 1.};


    vector<double> *Ut = new vector<double>(ut.begin(), ut.end());
    vector<vector<double>> *Pt = new vector<vector<double>>(pt.begin(), pt.end());


    auto* full = new nurbsCur(O, X, Y, r, ths, the);

    auto* thic = new nurbsCur(Pt, Ut, 1);

    vector<double> Uf = full->getKnot();
    vector<vector<double>> Pf = full->getCP();


    the = 180.;

    auto* half = new nurbsCur(O, Z, X, r, ths, the);

    auto* sphere = new nurbsSur(*full, *half);

    auto* shell = new nurbsVol(*sphere, *thic);

    shell->makePw();

    sphere->makePw();

    full->makePw();
    half->makePw();

    auto cps = half->getCP();

    auto Ps = sphere->getCP();

    auto test = cart2par(vector<double> {15., 17., 19.});

    auto* tspline = new tsplineVol(*shell);

    tspline->evaluate(0.5, 0.5, 0.5);

    auto he =tspline->getRes();

    for (auto&hst : he) {
        for (auto&h : hst) {
        cout << h << endl;
        }
        cout << endl;
    }

    
    delete tspline;
    delete half;
    delete sphere;
    delete full;

}