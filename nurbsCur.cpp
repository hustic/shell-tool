// Sotiris Gkoulmaris - sg3219
// spline-tool

// NURBS curve derived class cpp file. Contains implementation of curve
// generators, operators and samplers.

#include <cmath>
#include <iostream>
#include <fstream>

#include "spline.cpp"

using namespace std;

/* Constructor used to define a curve, given Control Points, Knot vector
** and degree. It is used primarily to define the thickness curve for the
** shell construction.
** Input: control points P, knot vector U, degree p
** Output: creates a NURBS curve object
*/
nurbsCur::nurbsCur(vector<vector<double>> *P, vector<double> *U, int p)
{
    this->Pw.reset(P);
    this->V.reset(U);
    this->p = p;
}

/* Class destructor - since we are using unique_ptr, it can be empty
*/
nurbsCur::~nurbsCur()
{}

/* Evaluate curve given a step to go from 0 to 1
    ** Output: populate out object
    */
void nurbsCur::evaluate(double step)
{
    double start = 0.;
    double end = 1.;

    vector<vector<double>> *out = new vector<vector<double>>();

    int degree = this->p;
    vector<double> knot = this->getKnot();
    vector<vector<double>> points = this->getCP();

    while (start <= end) {
        int uspan = fspan(degree, start, knot);
        vector<double> Nu = basisFun(uspan, start, degree, knot);

        vector<double> Cw {0., 0., 0., 0.};
        vector<double> C {0., 0., 0.};

        for (int i = 0; i <= degree; i++) {
            for (int j = 0; j < points[0].size(); j++) {
                Cw[j] = Cw[j] + Nu[i] * points[uspan-degree+i][j];
            }
        }

        
        for (int j = 0; j < points[0].size()-1; j++) {
            C[j] = Cw[j]/Cw[3];
        }
        out->push_back(C);
        start += step;
    }
    
    this->out.reset(out);

}
/* Control points accessor
** Output: a vector of vector of control points 
*/
vector<vector<double>> nurbsCur::getCP()
{
    return *this->Pw;
}

/* Knot vector accessor
** Output: a vector of knot values
*/
vector<double> nurbsCur::getKnot()
{
    return *this->V;
}

/* Evaluation output accessor
** Output: a vector of evaluated points
*/
vector<vector<double>> nurbsCur::getRes()
{
    return *this->out;
}


/* Compute weighted control points and populate the class object Pw
*/
void nurbsCur::makePw() {
    vector<vector<double>> *Pw = new vector<vector<double>>();
    vector<vector<double>> points = this->getCP();

    vector<double> p {0., 0., 0., 0.};
    for (int i = 0; i < points.size(); i++) {
        for (int j = 0; j < p.size()-1; j++) {
            p[j] = points[i][j] * points[i][3];
        }
        p[3] = points[i][3];
        Pw->push_back(p);
    }

    this->Pw.reset(Pw);

}

/* Constructor used to define a circular curve. It is used to create the
** half and full circle parametric representation of the sphere.
** Input: origin O, axis vectors X and Y, radius r, degrees start and end
** Output: creates a circular NURBS curve object 
*/
nurbsCur::nurbsCur(vector<int>& O, vector<int>& X, vector<int>& Y, double r, double start, double end)
{
    double ths = start;
    double the = 0.;
    double theta = 0.;
    int narcs = 0;
    
    if (end < start) {
        the = 360. + end;
    } else {
        the = end;
    }

    theta = the - ths;

    if (theta <= 90.) {
        narcs = 1;
    } else {
        if (theta <= 180.) {
            narcs = 2;
        } else {
            if (theta <= 270.) {
                narcs = 3;
            } else {
                narcs = 4;
            }
        }
    }

    
    int n = 2*narcs;

    double dtheta = theta/narcs;


    vector<vector<double>> *Pw = new vector<vector<double>>();
    double w1 = cos(dtheta/2. * M_PI/180.);

    vector<double> P0;
    vector<double> T0;

    for (int i = 0; i < O.size(); i++) {
        P0.push_back(O[i] + r * cos(ths*M_PI/180.)*X[i] + r*sin(ths*M_PI/180.)*Y[i]);
        T0.push_back(-sin(ths*M_PI/180.)*X[i] + cos(ths*M_PI/180.)*Y[i]);
    }
    P0.push_back(1.);

    Pw->push_back(P0);

    int index = 0;
    double angle = ths;

    for (int i = 1; i <= narcs; i++) {
        angle = angle + dtheta;

        vector<double> P2;
        vector<double> T2;

        

        for (int d = 0; d < O.size(); d++) {
            P2.push_back(O[d] + r * cos(angle*M_PI/180.)*X[d] + r * sin(angle*M_PI/180.)*Y[d]);
            T2.push_back(-sin(angle*M_PI/180.)*X[d] + cos(angle*M_PI/180.)*Y[d]);
        }

        vector<double> P1;
        if (X[0] == 1){
            vector<double> p0 {P0[0], P0[1]};
            vector<double> t0 {T0[0], T0[1]};
            vector<double> p2 {P2[0], P2[1]};
            vector<double> t2 {T2[0], T2[1]};
            vector<double> out;
            out = getIntersect(p0, t0, p2, t2);

            P1.push_back(out[0]);
            P1.push_back(out[1]);
            P1.push_back(0.);
        }

        if (X[2] == 1){
            vector<double> p0 {P0[0], P0[2]};
            vector<double> t0 {T0[0], T0[2]};
            vector<double> p2 {P2[0], P2[2]};
            vector<double> t2 {T2[0], T2[2]};
            vector<double> out;
            out = getIntersect(p0, t0, p2, t2);

            P1.push_back(out[0]);
            P1.push_back(0.);
            P1.push_back(out[1]);
        }


        P1.push_back(w1);
        Pw->push_back(P1);
        P2.push_back(1.);
        Pw->push_back(P2);
        if (i < narcs) {
            P0 = P2;
            T0 = T2;
        }
    }


    int j = 2*narcs+1;
    vector<double> *u = new vector<double>();

    for (int i = 0; i < 3; i++) {
        u->push_back(0.);
    }

    switch (narcs)
    {
    case 1:
        break;
    case 2:
        u->push_back(0.5);
        u->push_back(0.5);
        break;
    case 3:
        u->push_back(1.0/3.0);
        u->push_back(1.0/3.0);
        u->push_back(2.0/3.0);
        u->push_back(2.0/3.0);
        break;
    case 4:
        u->push_back(0.25);
        u->push_back(0.25);
        u->push_back(0.5);
        u->push_back(0.5);
        u->push_back(0.75);
        u->push_back(0.75);
    default:
       break;
    }

    for (int i = 0; i < 3; i++) {
        u->push_back(1.);
    }

    this->Pw.reset(Pw);
    this->V.reset(u);
    this->p = 2;

}

/* Declaration of virtual function to export out to csv
*/
void nurbsCur::write_csv(string filename)
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