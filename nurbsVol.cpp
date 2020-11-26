// Sotiris Gkoulmaris - sg3219
// spline-tool

// NURBS volume derived class cpp file. Contains implementations of volume
// generators, operators and samplers.

#include <cmath>
#include <iostream>
#include <fstream>


#include "nurbsVol.h"


using namespace std;

/* Class destructor - since we are using unique_ptr, it can be empty
*/
nurbsVol::~nurbsVol()
{}


/* Constructor used to define a volume, given the parametric representation of the
* sphere and a curve for thickness
* Input: sphere nurbsSur, nurbsCur
* Output: creates a NURBS volume object (spherical shell)
*/
nurbsVol::nurbsVol(nurbsSur &sphere, nurbsCur &thicc)
{
    vector<vector<vector<vector<double>>>> *P = new vector<vector<vector<vector<double>>>>();

    vector<double> Uh = sphere.getKnotV();
    vector<double> Uf = sphere.getKnotU();
    vector<double> Ut = thicc.getKnot();

    vector<vector<double>> Pt = thicc.getCP();
    vector<vector<vector<double>>> Ss = sphere.getCP();

    vector<vector<vector<double>>> ptem;
    vector<vector<double>> temp;
    vector<double> p;
    for (int k = 0; k < Pt.size(); k++) {
        ptem.clear();
        temp.clear();
        p.clear();

        for (int i = 0; i < Ss.size(); i++) {
            double aspect = abs(Pt[k][0]/Ss[0][0][2]);
            temp.clear();

            for (int j = 0; j < Ss[0].size(); j++) {
                p.clear();
                for (int t = 0; t < Ss[0][0].size()-1; t++) {
                    p.push_back(Ss[i][j][t] * aspect);
                }
                p.push_back(Ss[i][j][3]);
                temp.push_back(p);
            }
            ptem.push_back(temp);
        }
        P->push_back(ptem);
    }

    vector<double> *V = new vector<double>(Uh);
    vector<double> *U = new vector<double>(Uf);
    vector<double> *H = new vector<double>(Ut);

    this->P.reset(P);
    this->V.reset(V);
    this->U.reset(U);
    this->H.reset(H);
    this->p = sphere.p;
    this->q = sphere.q;
    this->r = thicc.p;
}

/* Control points accessor
** Output: a vector of vector of vector of control points 
*/
vector<vector<vector<vector<double>>>> nurbsVol::getCP()
{
    return *this->P;
}

/* Weighted control points accessor
* Output: a vector of vectors of vectors of control points 
*/
vector<vector<vector<vector<double>>>> nurbsVol::getWCP()
{
    return *this->Pw;
}

/* Compute weighted control points and populate the class object Pw
*/
void nurbsVol::makePw() {

    vector<vector<vector<vector<double>>>> *Pw = new vector<vector<vector<vector<double>>>>();
    vector<vector<vector<vector<double>>>> points = this->getCP();

    vector<vector<vector<double>>> pps;
    vector<vector<double>> ps;
    vector<double> p;
    for (int t = 0; t < points.size(); t++) {
        pps.clear();
        for (int i = 0; i < points[0].size(); i++) {
            ps.clear();
            for (int j = 0; j < points[0][0].size(); j++) {
                p.clear();
                for (int k = 0; k < points[0][0][0].size()-1; k++) {
                    p.push_back(points[t][i][j][k] * points[t][i][j][3]);
                }
                p.push_back(points[t][i][j][3]);
                ps.push_back(p);
            }
            pps.push_back(ps); 
        }
        Pw->push_back(pps);
    }
    this->Pw.reset(Pw);

}

/* Global knot refinement algorithm, expanded from NURBS Book A5.5
* Input: int direction, vector of knots to be added
* Output: update Pw, U, V, H with new knots and control points
*/
void nurbsVol::refine(int dir, vector<double> knots) 
{
     
    auto U = this->getKnotU();
    auto V = this->getKnotV();
    auto H = this->getKnotH();


    auto Pw = this->getWCP();

    
    
    vector<vector<double>> qq;
    vector<vector<vector<double>>> qqq;
    vector<vector<vector<vector<double>>>> Qw;
    
    switch(dir) {
        case 0:
        {
            

             vector<double> Ubar;

            int z = U.size()-1;
            int a = fspan(this->q, knots[0], U);
            int b = fspan(this->q, knots[knots.size()-1], U) + 1;

            
            for (int j = 0; j <= a; j++) {
                Ubar.push_back(U[j]);
            }
            for (int j = a+1; j < b+this->q+knots.size(); j++) {
                Ubar.push_back(0.0);
            }
            for (int j = b+this->q; j <= z; j++) {
                Ubar.push_back(U[j]);
            }
            
            
            vector<double> q {0., 0., 0., 0.};
            for (int up = 0; up < Pw.size(); up++) {
                qqq.clear();
                for (int row = 0; row < Pw[0].size(); row++) {
                    qq.clear();
                    for (int j = 0; j <= a-this->q; j++) {
                        qq.push_back(Pw[up][row][j]);
                    }
                    for (int j = a-this->q+1; j < b-1+knots.size(); j++) {
                        
                        qq.push_back(q);
                    } 
                    for (int j = b-1; j < Pw[0][0].size(); j++) {
                        qq.push_back(Pw[up][row][j]);
                    }
                    qqq.push_back(qq);
                }
                Qw.push_back(qqq);
            }

            int i = b+this->q-1;
            int k = b+this->q+knots.size()-1;

            for (int j = knots.size()-1; j >= 0; j--) {
                while (knots[j] <= U[i] && i > a) {
                    for (int up = 0; up < Pw.size(); up++) {
                        for (int row = 0; row < Pw[0].size(); row++) {
                            Qw[up][row][k-this->q-1] = Pw[up][row][i-this->q-1];
                        }
                    }
                    Ubar[k] = U[i];
                    k--;
                    i--;
                }

                for (int up = 0; up < Pw.size(); up++) {
                    for (int row = 0; row < Pw[0].size(); row++) {
                        Qw[up][row][k-this->q-1] = Qw[up][row][k-this->q];
                    }
                }

                for (int l = 1; l <= this->q; l++) {
                    int ind = k-this->q+l;
                    double alfa = Ubar[k+l]-knots[j];
                    if (abs(alfa) == 0.0) {
                        for (int up = 0; up < Pw.size(); up++) {
                            for (int row = 0; row < Pw[0].size(); row++) {
                                Qw[up][row][ind-1] = Qw[up][row][ind];
                            }
                        } 
                    } else {
                        alfa = alfa/(Ubar[k+l]-U[i-this->q+l]);
                        for (int up = 0; up < Pw.size(); up++) {
                            for (int row = 0; row < Pw[0].size(); row++) {
                                for (int d = 0; d < Pw[0][0][0].size(); d++) {
                                    Qw[up][row][ind-1][d] = alfa * Qw[up][row][ind-1][d] + (1.0 - alfa) * Qw[up][row][ind][d];
                                }
                            }
                        } 
                    }
                }
                Ubar[k] = knots[j];
                k--;
            }

            vector<double> *Uu = new vector<double>(Ubar);
            vector<vector<vector<vector<double>>>> *Qqw = new vector<vector<vector<vector<double>>>>(Qw);

            this->U.reset(Uu);
            this->Pw.reset(Qqw);
            
            break;
        }
        case 1:
        {    
        
             vector<double> Vbar;

            int z = V.size()-1;
            int a = fspan(this->p, knots[0], V);
            int b = fspan(this->p, knots[knots.size()-1], V) + 1;

            
            for (int j = 0; j <= a; j++) {
                Vbar.push_back(V[j]);
            }
            for (int j = a+1; j < b+this->p+knots.size(); j++) {
                Vbar.push_back(0.0);
            }
            for (int j = b+this->p; j <= z; j++) {
                Vbar.push_back(V[j]);
            }
            
            
            vector<double> q {0., 0., 0., 0.};
            for (int up = 0; up < Pw.size(); up++) {
                qqq.clear();
                for (int row = 0; row <= a-this->p; row++) {
                    qq.clear();
                    for (int j = 0; j < Pw[0][0].size() ; j++) {
                         qq.push_back(Pw[up][row][j]);
                    }
                    qqq.push_back(qq);
                }
                
                for (int row = a-this->p+1; row < b-1+knots.size(); row++) {
                    qq.clear();
                    for (int j = 0; j < Pw[0][0].size() ; j++) {
                         qq.push_back(q);
                    }
                    qqq.push_back(qq);
                }

                for (int row = b-1; row < Pw[0].size(); row++) {
                    qq.clear();
                    for (int j = 0; j < Pw[0][0].size() ; j++) {
                         qq.push_back(Pw[up][row][j]);
                    }
                    qqq.push_back(qq);
                }
                Qw.push_back(qqq);
            }

            int i = b+this->p-1;
            int k = b+this->p+knots.size()-1;

            for (int j = knots.size()-1; j >= 0; j--) {
                while (knots[j] <= V[i] && i > a) {
                    for (int up = 0; up < Pw.size(); up++) {
                        for (int row = 0; row < Pw[0][0].size(); row++) {
                            Qw[up][k-this->p-1][row] = Pw[up][i-this->p-1][row];
                        }
                    }
                    Vbar[k] = V[i];
                    k--;
                    i--;
                }

                for (int up = 0; up < Pw.size(); up++) {
                    for (int row = 0; row < Pw[0][0].size(); row++) {
                        Qw[up][k-this->p-1][row] = Qw[up][k-this->p][row];
                    }
                }

                for (int l = 1; l <= this->p; l++) {
                    int ind = k-this->p+l;
                    double alfa = Vbar[k+l]-knots[j];
                    if (abs(alfa) == 0.0) {
                        for (int up = 0; up < Pw.size(); up++) {
                            for (int row = 0; row < Pw[0][0].size(); row++) {
                                Qw[up][ind-1][row] = Qw[up][ind][row];
                            }
                        } 
                    } else {
                        alfa = alfa/(Vbar[k+l]-V[i-this->p+l]);
                        for (int up = 0; up < Pw.size(); up++) {
                            for (int row = 0; row < Pw[0][0].size(); row++) {
                                for (int d = 0; d < Pw[0][0][0].size(); d++) {
                                    Qw[up][ind-1][row][d] = alfa * Qw[up][ind-1][row][d] + (1.0 - alfa) * Qw[up][ind][row][d];
                                }
                            }
                        } 
                    }
                }
                Vbar[k] = knots[j];
                k--;
            }

            vector<double> *Vv = new vector<double>(Vbar);
            vector<vector<vector<vector<double>>>> *Qqw = new vector<vector<vector<vector<double>>>>(Qw);

            this->V.reset(Vv);
            this->Pw.reset(Qqw);
            
            break;
        }
        case 2:
        {

             vector<double> Hbar;

            int z = H.size()-1;
            int a = fspan(this->r, knots[0], H);
            int b = fspan(this->r, knots[knots.size()-1], H) + 1;


            
            for (int j = 0; j <= a; j++) {
                Hbar.push_back(H[j]);
            }
            for (int j = a+1; j < b+this->r+knots.size(); j++) {
                Hbar.push_back(0.0);
            }
            for (int j = b+this->r; j <= z; j++) {
                Hbar.push_back(H[j]);
            }

            
            
            
            vector<double> q {0., 0., 0., 0.};
            for (int up = 0; up <= a-this->r; up++) {
                qqq.clear();
                for (int row = 0; row < Pw[0].size(); row++) {
                    qq.clear();
                    for (int j = 0; j < Pw[0][0].size(); j++) {
                        qq.push_back(Pw[up][row][j]);
                    }
                    qqq.push_back(qq);
                }
                Qw.push_back(qqq);
            }

            for (int up = a-this->r+1; up < b-1+knots.size(); up++) {
                qqq.clear();
                for (int row = 0; row < Pw[0].size(); row++) {
                    qq.clear();
                    for (int j = 0; j < Pw[0][0].size(); j++) {
                        qq.push_back(q);
                    }
                    qqq.push_back(qq);
                }
                Qw.push_back(qqq);
            }

            for (int up = b-1; up < Pw.size(); up++) {
                qqq.clear();
                for (int row = 0; row < Pw[0].size(); row++) {
                    qq.clear();
                    for (int j = 0; j < Pw[0][0].size(); j++) {
                        qq.push_back(Pw[up][row][j]);
                    }
                    qqq.push_back(qq);
                }
                Qw.push_back(qqq);
            }

            int i = b+this->r-1;
            int k = b+this->r+knots.size()-1;

            for (int j = knots.size()-1; j >= 0; j--) {
                while (knots[j] <= H[i] && i > a) {
                    for (int up = 0; up < Pw[0].size(); up++) {
                        for (int row = 0; row < Pw[0][0].size(); row++) {
                            Qw[k-this->r-1][up][row] = Pw[i-this->r-1][up][row];
                        }
                    }
                    Hbar[k] = H[i];
                    k--;
                    i--;
                }

                for (int up = 0; up < Pw[0].size(); up++) {
                    for (int row = 0; row < Pw[0][0].size(); row++) {
                        Qw[k-this->r-1][up][row] = Qw[k-this->r][up][row];
                    }
                }

                for (int l = 1; l <= this->p; l++) {
                    int ind = k-this->r+l;
                    double alfa = Hbar[k+l]-knots[j];
                    if (abs(alfa) == 0.0) {
                        for (int up = 0; up < Pw[0].size(); up++) {
                            for (int row = 0; row < Pw[0][0].size(); row++) {
                                Qw[ind-1][up][row] = Qw[ind][up][row];
                            }
                        } 
                    } else {
                        alfa = alfa/(Hbar[k+l]-H[i-this->r+l]);
                        for (int up = 0; up < Pw[0].size(); up++) {
                            for (int row = 0; row < Pw[0][0].size(); row++) {
                                for (int d = 0; d < Pw[0][0][0].size(); d++) {
                                    Qw[ind-1][up][row][d] = alfa * Qw[ind-1][up][row][d] + (1.0 - alfa) * Qw[ind][up][row][d];
                                }
                            }
                        } 
                    }
                }
                Hbar[k] = knots[j];
                k--;
            }

            vector<double> *Hh = new vector<double>(Hbar);
            vector<vector<vector<vector<double>>>> *Qqw = new vector<vector<vector<vector<double>>>>(Qw);

            this->H.reset(Hh);
            this->Pw.reset(Qqw);
            
            break;
        }
    }
}


/* Evaluate volume given steps; from 0 to 1
* Output: populate out object
*/
void nurbsVol::evaluate(double stepu, double stepv, double steph) 
{
    
    double vstart = 0.;
    double ustart = 0.;
    double hstart = 0.;

    double vend = 1.;
    double uend = 1.;
    double hend = 1.;

    

    vector<vector<double>> *out = new vector<vector<double>>();
    

    int dv = this->p;
    int du = this->q;
    int dr = this->r;
    
    vector<double> knv = this->getKnotV();

    vector<double> knu = this->getKnotU();

    vector<double> knh = this->getKnotH();

    vector<double> hsteps;
    vector<double> vsteps;
    vector<double> usteps;

    double h = 0.;
    while (hstart < hend ) {
        hsteps.push_back(hstart);
        hstart = hstart + steph;
    }
    hsteps.push_back(1.);

    double v = 0.;
    while (vstart < vend ) {
        vsteps.push_back(vstart);
        vstart = vstart + stepv;
    }
    if (floor(vsteps[vsteps.size()-1]) != 1.) {
        vsteps.push_back(1.);
    } else {
        vsteps.pop_back();
        vsteps.push_back(1.);
    }
    

    double u = 0.;
    while (ustart < uend ) {
        usteps.push_back(ustart);
        ustart = ustart + stepu;
    }
    usteps.push_back(1.);

    this->thicc = hsteps.size();
    this->up = vsteps.size();
    this->right = usteps.size();

    
    vector<vector<vector<vector<double>>>> points = this->getWCP();



    for (auto&hst : hsteps) {

        int hspan = fspan(dr, hst, knh);
        vector<double> Nh = basisFun(hspan, hst, dr, knh);

        for (auto&vst : vsteps) {
           
            int vspan = fspan(dv, vst, knv);
            vector<double> Nv = basisFun(vspan, vst, dv, knv);
        
        
            for (auto&ust : usteps) {
                
                int uspan = fspan(du, ust, knu);
                vector<double> Nu = basisFun(uspan, ust, du, knu);


                vector<double> Sw {0., 0., 0., 0.};
                vector<double> S {0., 0., 0.};

                for (int l = 0; l <= du; l++) {
                    vector<double> t {0., 0., 0., 0.};
                    for (int k = 0; k <= dv; k++) {
                        vector<double> temp {0., 0., 0., 0.};
                        for (int z = 0; z <= dr; z++) {
                            for (int j = 0; j < points[0][0][0].size(); j++) {
                                temp[j] += Nh[z] * points[hspan-dr+z][vspan-dv+l][uspan-du+k][j];
                            }
                        }
                        for (int j = 0; j < points[0][0][0].size(); j++) {
                            t[j] += Nu[k] * temp[j];
                        }
                    }
                    for (int j = 0; j < points[0][0][0].size(); j++) {
                        Sw[j] += Nv[l] * t[j];
                    }
                }
                for (int j = 0; j < points[0][0][0].size()-1; j++) {
                    S[j] = Sw[j]/Sw[3];
                }
                out->push_back(S);
            }
        }
    }
    this->out.reset(out);
}

/* Evaluate volume given vectors of steps; from 0 to 1
* Output: populate out object
*/
void nurbsVol::evaluate(vector<double> usteps, vector<double> vsteps, vector<double> hsteps) 
{

    vector<vector<double>> *out = new vector<vector<double>>();
    

    int dv = this->p;
    int du = this->q;
    int dr = this->r;
    
    vector<double> knv = this->getKnotV();

    vector<double> knu = this->getKnotU();

    vector<double> knh = this->getKnotH();
    
    vector<vector<vector<vector<double>>>> points = this->getWCP();


    this->thicc = hsteps.size();
    this->up = vsteps.size();
    this->right = usteps.size();



    for (auto&hst : hsteps) {

        int hspan = fspan(dr, hst, knh);
        vector<double> Nh = basisFun(hspan, hst, dr, knh);

        for (auto&vst : vsteps) {
           
            int vspan = fspan(dv, vst, knv);
            vector<double> Nv = basisFun(vspan, vst, dv, knv);
        
        
            for (auto&ust : usteps) {
                
                int uspan = fspan(du, ust, knu);
                vector<double> Nu = basisFun(uspan, ust, du, knu);


                vector<double> Sw {0., 0., 0., 0.};
                vector<double> S {0., 0., 0.};

                for (int l = 0; l <= du; l++) {
                    vector<double> t {0., 0., 0., 0.};
                    for (int k = 0; k <= dv; k++) {
                        vector<double> temp {0., 0., 0., 0.};
                        for (int z = 0; z <= dr; z++) {
                            for (int j = 0; j < points[0][0][0].size(); j++) {
                                temp[j] += Nh[z] * points[hspan-dr+z][vspan-dv+l][uspan-du+k][j];
                            }
                        }
                        for (int j = 0; j < points[0][0][0].size(); j++) {
                            t[j] += Nu[k] * temp[j];
                        }
                    }
                    for (int j = 0; j < points[0][0][0].size(); j++) {
                        Sw[j] += Nv[l] * t[j];
                    }

                    
                    
                }
                for (int j = 0; j < points[0][0][0].size()-1; j++) {
                    S[j] = Sw[j]/Sw[3];
                }
                out->push_back(S);
            }
        }
    }
    this->out.reset(out);
}

/* Knot vector V accessor
** Output: a vector of knot values
*/
vector<double> nurbsVol::getKnotV()
{
    return *this->V;
}

/* Knot vector U accessor
** Output: a vector of knot values
*/
vector<double> nurbsVol::getKnotU()
{
    return *this->U;
}

/* Knot vector U accessor
** Output: a vector of knot values
*/
vector<double> nurbsVol::getKnotH()
{
    return *this->H;
}

/* Helper function to find the knot values of the anchors
* Input: degree p, knot vector V
* Output: anchor knot values
*/
vector<double> getKnotAnchors(int p, vector<double> V) 
{
    vector<double> knots;
    int ind;
    if (p % 2 == 0) {
        for (int i = 0; i < V.size()-p-1; i++) {
            ind = int(i + p/2.);
            knots.push_back(1./2. * double(V[ind] + V[ind+1]));
        }
    } else {
        for (int i = 0; i < V.size()-p-1; i++) {
            ind = int(i + (p+1)/2.);
            knots.push_back(V[ind]);
        }
    }

    return knots;
}


/* Evaluation output accessor
* Output: a vector of evaluated points
*/
vector<vector<double>> nurbsVol::getRes() 
{
    return *this->out;
}

/* Declaration of virtual function to export out to csv
*/
void nurbsVol::write_csv(string filename)
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

/* Export control points to csv
*/
void nurbsVol::writeCP_csv(string filename) {
    ofstream file;

    file.open(filename);

    auto CP = this->getWCP();

    file << "x,y,z" << endl;

    for (int i = 0; i < CP.size(); i++) {
        for (int j = 0; j < CP[0].size(); j++) {
            for (int k = 0; k < CP[0][0].size(); k++) {
                for (int l = 0; l < CP[0][0][0].size(); l++) {
                    file << CP[i][j][k][l] << ",";
                }
                file << endl;
            }
            
        }
        
    }

    file.close();
}

/* Export evaluated volume to vtk
*/
void nurbsVol::write_vtk(string filename) 
{
    ofstream file;

    file.open(filename);

    auto out = this->getRes();
    auto CP = this->getWCP();

    int thicc = this->thicc;
    int up = this->up;
    int right = this->right;

    int num = thicc * (up-2) * (right-1) + (thicc*2);

    file << "# vtk DataFile Version 2.0\n";
    file << "Unstructured Grid Shell\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    file << "POINTS " << to_string(num) << " float\n";

    for (int k = 0; k < thicc; k++){
        for (int j = 0; j < up; j++) {
            for (int i = 0; i < right; i++) {
                if (i == right-1) {
                    break;
                }
                file << to_string(out[(k*(up)*right)+(j*right)+i][0]) << "  " << to_string(out[(k*up*right)+(j*right)+i][1]) << "  " << to_string(out[(k*up*right)+(j*right)+i][2]) << "   ";
                if ((j == 0 || j == up-1)) {
                    break;
                }
            }
        }
    }
    
    int cpts = (7*2*(right-1) + 9*(up-3)*(right-1))*(thicc-1);
    file << endl;
    file << endl;
    file << "CELLS " << to_string((right-1)*(up-1)*(thicc-1)) << " " << to_string(cpts) << "\n";

    int npts = 0;
    int b;
    int z;
    int d;
    for (int k = 0; k < thicc-1; k++) {
        d = k * ((up-2)*(right-1)+2);

        for (int j = 0; j < up-1; j++) {
            b = 0;
            if ((j == 0) || (j == up-2)) {
                npts = 6;
            } else {
                npts = 8;
            }

            if (j == 0) b = 1;

            for (int i = 0; i < right-1; i++) {
                file << to_string(npts) << " ";
                
                if (i == right-2) z = 1;
                else z = 0;
                
                if (j == 0) {
                    file << to_string(d) << " ";
                    file << to_string(d+((i+1)%(right))) <<" "<< to_string(d+((i+2)%(right))+z) <<" "<< to_string(d+(up-2)*(right-1)+2)<<" ";
                    file << to_string(d+(up-2)*(right-1)+2+((i+1)%(right))) <<" "<< to_string(d+(up-2)*(right-1)+2+((i+2)%(right))+z) <<" ";
                } else if (j == up-2) {
                    file << to_string(d+right-1) <<" "<< to_string(d+((i+1)%(right))-1) <<" "<< to_string(d+((i+2)%(right))+z-1) <<" ";
                    file << to_string(d+(right-1)*(up-2)+2+right-1) <<" "<< to_string(d+(up-2)*(right-1)+2+((i+1)%(right))-1) <<" ";
                    file << to_string(d+(up-2)*(right-1)+2+((i+2)%(right))+z-1) <<" ";
                } else {
                    file << to_string(d+i) <<" "<< to_string(d+((i+1)%(right-1))) <<" "<< to_string(d+(right-1)+((i+1)%(right-1))) <<" ";
                    file << to_string(d+(right-1)+i) <<" "<< to_string(d+(up-2)*(right-1)+2+i) <<" "<< to_string(d+(up-2)*(right-1)+2+((i+1)%(right-1))) <<" ";
                    file << to_string(d+(right-1)+(up-2)*(right-1)+2+((i+1)%(right-1))) <<" "<< to_string(d+(right-1)+(up-2)*(right-1)+2+i) <<" ";
                    b++;
                }
                file << endl;
            }
            d = d+b;
        }
    }

    file << endl;
    file << "CELL_TYPES " << to_string((right-1)*(up-1)*(thicc-1)) << endl;

    for (int k = 0; k < thicc-1; k++) {
        for (int j = 0; j < up-1; j++) {
            for (int i = 0; i < right-1; i++) {
                if ((j == 0) || (j == up-2)) file << "13" << endl;
                else file << "12" << endl;
            }
        }
    }

    file << endl;

    file << "POINT_DATA " + to_string(num) << "\n";
    file << "SCALARS scalars float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < num; i++) {
        file << "0";
        file << " ";
    }

    file << endl;


    file.close();

}