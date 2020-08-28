// Sotiris Gkoulmaris - sg3219
// spline-tool

// NURBS surface derived class cpp file. Contains implementations of surface
// generators, operators and samplers.

#include <cmath>
#include <iostream>
#include <fstream>


#include "nurbsSur.h"


using namespace std;

/* Class destructor - since we are using unique_ptr, it can be empty
*/
nurbsSur::~nurbsSur()
{}


/* Constructor used to define a surface, given the parametric representation of the
** sphere.
** Input: half and full circle nurbsCurs
** Output: creates a NURBS surface object (sphere)
*/
nurbsSur::nurbsSur(nurbsCur &full, nurbsCur &half)
{
    vector<vector<vector<double>>> *Pw = new vector<vector<vector<double>>>();
    double wm = cos(45. * M_PI/180.);

    vector<double> Uh = half.getKnot();
    vector<vector<double>> Ph = half.getCP();

    vector<double> Uf = full.getKnot();
    vector<vector<double>> Pf = full.getCP();
    //cout << Ph.size() << endl;
    vector<vector<double>> temp;
    vector<double> p;
    for (int i = 0; i < Ph.size(); i++) {
        temp.clear();
        p.clear();

        for (int d = 0; d < 3; d++) {
            p.push_back(Ph[i][d]);
        }
        p.push_back(Ph[i][3]);

        temp.push_back(p);

        for (int j = 1; j < full.Pw->size(); j++) {
            p.clear();
            if ((i == Ph.size() - 1 ) || (i == 0)) {
                for (int d = 0; d < 3; d++) {
                    p.push_back(Ph[i][d]);
                }
                p.push_back(Ph[i][3]);
                
            } else {
                double w = 0.;
                if (j % 2 == 0) {
                    w = Ph[i][3];
                } else {
                    w = Ph[i][3] * Pf[j][3];
                }
                for (int d = 0; d < 2; d++) {
                    p.push_back(Pf[j][d]);
                }
                p.push_back(Ph[i][2]);
                p.push_back(w);
            }
            temp.push_back(p);
        }
        Pw->push_back(temp);
    }

    vector<double> *V = new vector<double>(Uh);
    vector<double> *U = new vector<double>(Uf);

    this->Pw.reset(Pw);
    this->V.reset(V);
    this->U.reset(U);
    this->p = half.p;
    this->q = full.p;
}

/* Control points accessor
*  Output: a vector of vector of vector of control points 
*/
vector<vector<vector<double>>> nurbsSur::getCP()
{
    return *this->Pw;
}

/* Compute weighted control points and populate the class object Pw
*/
void nurbsSur::makePw() {
    vector<vector<vector<double>>> *Pw = new vector<vector<vector<double>>>();
    vector<vector<vector<double>>> points = this->getCP();

    

    for (int i = 0; i < points.size(); i++) {
        vector<vector<double>> ps;
        
        
        for (int j = 0; j < points[0].size(); j++) {
            vector<double> p {0., 0., 0., 0.};
            for (int k = 0; k < p.size()-1; k++) {
                p[k] = points[i][j][k] * points[i][j][3];
            }
            p[3] = points[i][j][3];
            ps.push_back(p);
        }
        Pw->push_back(ps);
    }

    this->Pw.reset(Pw);
}

/* Evaluate surface given steps; from 0 to 1
*  Output: populate out object
*/
void nurbsSur::evaluate(double stepu, double stepv) 
{
    double vstart = 0.;
    double ustart = 0.;

    double vend = 1.;
    double uend = 1.;

    vector<vector<double>> *out = new vector<vector<double>>();

    int dv = this->p;
    int du = this->q;
    
    vector<double> knv = this->getKnotV();
    vector<double> knu = this->getKnotU();
    
    vector<vector<vector<double>>> points = this->getCP();


    while (vstart <= vend) {
        ustart = 0.;

        int vspan = fspan(dv, vstart, knv);
        vector<double> Nv = basisFun(vspan, vstart, dv, knv);
        
        
        while (ustart <= uend) {
            int uspan = fspan(du, ustart, knu);
            vector<double> Nu = basisFun(uspan, ustart, du, knu);

            int uind = uspan-du;
            vector<double> Sw {0., 0., 0., 0.};
            vector<double> S {0., 0., 0.};

            for (int l = 0; l <= du; l++) {
                vector<double> temp {0., 0., 0., 0.};
                int vind = vspan-dv+l;
                for (int k = 0; k <= dv; k++) {
                    for (int j = 0; j < points[0][0].size(); j++) {
                        temp[j] += Nu[k] * points[vind][uind+k][j];
                    }
                }

                for (int j = 0; j < points[0][0].size(); j++) {
                    Sw[j] += Nv[l] * temp[j];
                }
                
            }

            

            for (int j = 0; j < points[0][0].size()-1; j++) {
                S[j] = Sw[j]/Sw[3];
            }
            out->push_back(S);
            ustart += stepu;
        }
        vstart += stepv;
    }

    this->out.reset(out);
}

/* Knot vector V accessor
*  Output: a vector of knot values
*/
vector<double> nurbsSur::getKnotV()
{
    return *this->V;
}

/* Knot vector U accessor
*  Output: a vector of knot values
*/
vector<double> nurbsSur::getKnotU()
{
    return *this->U;
}

/* Evaluation output accessor
* Output: a vector of evaluated points
*/
vector<vector<double>> nurbsSur::getRes()
{
    return *this->out;
}

/* Declaration of virtual function to export out to csv
*/
void nurbsSur::write_csv(string filename)
{
    ofstream file;

    file.open(filename);

    auto out = this->getRes();

    file << "x,y,z" << endl;

    for (int i = 0; i < out.size(); i++) {
        for (int j = 0; j < out[0].size(); j++) {
            file << out[i][j] << ",";
        }
        file << endl;
    }

    file.close();
}