// Sotiris Gkoulmaris - sg3219
// spline-tool

// T-Spline volume derived class cpp file. Contains implementations of volume
// generators, operators and samplers.

#include <cmath>
#include <iostream>
#include <fstream>


#include "tsplineVol.h"


using namespace std;

/* Class destructor - since we are using unique_ptr, it can be empty
*/
tsplineVol::~tsplineVol()
{}


/* Constructor used to define a T-spline volume, given the parametric representation of a
** volume.
** Input: vol nurbsVol
** Output: creates a T-Spline Volume object
*/
tsplineVol::tsplineVol(nurbsVol &vol)
{
    auto p = vol.getCP();

    auto pw = vol.getWCP();

    auto v = vol.getKnotV();
    auto u = vol.getKnotU();
    auto h = vol.getKnotH();

    auto sv = findex(v.size(), vol.p);
    auto su = findex(u.size(), vol.q);
    auto sh = findex(h.size(), vol.r);


    vector<vector<double>> tmeshv;

    for (int i = 0; i < sv.size(); i++) {
        vector<double> lockv(vol.p + 2);
        int indv;

        if (vol.p % 2 == 0) {
            indv = floor(sv[i]) - 1;
            for (int j = 0; j < int(vol.p/2 + 1); j++) {
                if (indv - j <= 0) {
                    lockv[int(vol.p/2)-j] = v[0];
                } else {
                    lockv[int(vol.p/2)-j] = v[indv-j+vol.p-1];
                }

                int end = sv.size();

                if (indv+j >= end) {
                    lockv[int(vol.p/2)+j+1] = v[end];
                } else {
                    lockv[int(vol.p/2)+j+1] = v[vol.p + indv + j];
                }
            }
        } else {
            indv = int(sv[i]);

            lockv[vol.p-1] = v[indv];

            for (int j = 0; j < vol.p-1; j++) {
                if (indv-j <= 0) {
                    lockv[vol.p-2-j] = v[0];
                } else {
                    lockv[vol.p-2-j] = v[indv-j-1];
                }

                int end = sv.size();

                if (indv + j + 1 >= end) {
                    lockv[vol.p+j] = v[end];
                } else {
                    lockv[vol.p+j] = v[indv + j + 1];
                }
            }
        }

        tmeshv.push_back(lockv);
    }



    vector<vector<double>> tmeshu;
    for (int i = 0; i < su.size(); i++) {
        vector<double> locku(vol.q + 2);
        int indu;

        if (vol.q % 2 == 0) {
            indu = floor(su[i]) - 1;
            for (int j = 0; j < int(vol.q/2 + 1); j++) {
                if (indu - j <= 0) {
                    locku[int(vol.q/2)-j] = u[0];
                } else {
                    locku[int(vol.q/2)-j] = u[indu-j+vol.q-1];
                }

                int end = su.size();

                if (indu+j >= end) {
                    locku[int(vol.q/2)+j+1] = u[end];
                } else {
                    locku[int(vol.q/2)+j+1] = u[vol.q + indu + j];
                }
            }
        } else {
            indu = int(su[i]);

            locku[vol.p-1] = u[indu];

            for (int j = 0; j < vol.q-1; j++) {
                if (indu-j <= 0) {
                    locku[vol.q-2-j] = u[0];
                } else {
                    locku[vol.q-2-j] = u[indu-j-1];
                }

                int end = su.size();

                if (indu + j + 1 >= end) {
                    locku[vol.q+j] = u[end];
                } else {
                    locku[vol.q+j] = u[indu + j + 1];
                }
            }
        }

        tmeshu.push_back(locku);
    }

    vector<vector<double>> tmeshh;
    for (int i = 0; i < sh.size(); i++) {
        vector<double> lockh(vol.r + 2);
        int indh;

        if (vol.r % 2 == 0) {
            indh = floor(sh[i]) - 1;
            for (int j = 0; j < int(vol.r/2 + 1); j++) {
                if (indh - j <= 0) {
                    lockh[int(vol.r/2)-j] = h[0];
                } else {
                    lockh[int(vol.r/2)-j] = h[indh-j+vol.r-1];
                }

                int end = sh.size();

                if (indh+j >= end) {
                    lockh[int(vol.r/2)+j+1] = h[end];
                } else {
                    lockh[int(vol.r/2)+j+1] = h[vol.r + indh + j];
                }
            }
        } else {
            indh = int(sh[i]);

            if (vol.r == 1) {
                lockh[1] = h[indh];
                lockh[0] = h[indh - 1];
                lockh[2] = h[indh + 1];
            } else {

                lockh[vol.r-1] = h[indh];

                for (int j = 0; j < vol.r-1; j++) {
                    if (indh-j <= 0) {
                        lockh[vol.r-2-j] = h[0];
                    } else {
                        lockh[vol.r-2-j] = h[indh-j-1];
                    }

                    int end = sh.size();

                    if (indh + j + 1 >= end) {
                        lockh[vol.r+j] = h[end];
                    } else {
                        lockh[vol.r+j] = h[indh + j + 1];
                    }
                }
            }
        }

        tmeshh.push_back(lockh);
    }

    //vector<vector<double>> pp;
    vector<vector<double>> ppw;
    vector<vector<double>> index;
    vector<vector<vector<double>>> local;
    vector<double> beta;

    // vertex * head;


    for (int k = 0; k < pw.size(); k++) {
        for (int j = 0; j < pw[0].size(); j++) {
            for (int i = 0; i < pw[0][0].size(); i++) {
                // vertex *v = new vertex;
                // v->beta = 1.;
                // v->P = pw[k][j][i];
                // v->index = vector<double> {su[i], sv[j], sh[k]};
                // v->knots = vector<vector<double>> {tmeshu[i], tmeshv[j], tmeshh[k]};

                // head = v;

                ppw.push_back(pw[k][j][i]);
                index.push_back(vector<double> {su[i], sv[j], sh[k]});
                local.push_back(vector<vector<double>> {tmeshu[i], tmeshv[j], tmeshh[k]});
                beta.push_back(1.);
                
            }
        }
    }

    this->n = pw[0][0].size();
    this->m = pw[0].size();
    this->l = pw.size();

    auto *Pw = new vector<vector<double>>(ppw);
    //auto *V = new vector<vector<double>>(tmeshv);
    auto *indexes = new vector<vector<double>>(index);
    auto *local_knots = new vector<vector<vector<double>>>(local);
    auto *betas = new vector<double>(beta);
    this->Pw.reset(Pw);
    this->local_knots.reset(local_knots);
    this->indexes.reset(indexes);
    this->betas.reset(betas);
    

}

/* Helper function for finding anchor index
*  Output: Vector of anchor indexes
*/
vector<double> findex(int n, int p) {
    vector<double> s;
    for (int i = 0; i < n-p-1; i++) {
        s.push_back(double(i) + (p+1.)/2.);
    }
    return s;
}

/* Declaration of virtual function to export out to csv
*/
void tsplineVol::write_csv(string filename)
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

/* Local knot vectors accessor
* Output: a vector of knot values
*/
vector<vector<vector<double>>> tsplineVol::getKnots() {
    return *this->local_knots;

}

/* Local knot refinement algorithm, adapted from S-Spline
* Input: anchor indexes x, y, z, knot value, direction
* Output: update Pw, betas, indexes, local knots with the new added points
*/
void tsplineVol::refine(double x, double y, double z, double knot, int dir) {
    auto indexes = this->getIndexes();
    auto pw = this->getWCP();

    int pold;

    auto local = this->getKnots();
    auto betas = this->getBetas();

    double bet;
    switch (dir)
    {
    case 0:
    {
        for (int i = indexes.size()-1; i > 0; i--) {
            if ((x < indexes[i][0])) {
                if ((y == indexes[i][1] && z == indexes[i][2])) {
                pold = i;
                }
            }
        }

        bet = (knot - local[pold][0][0]) / (local[pold][0][2] - local[pold][0][0]) * betas[pold];
        betas[pold] = (local[pold][0][3] - knot) / (local[pold][0][3] - local[pold][0][1]) * betas[pold];

        
        
        vector<double> nknot {local[pold][0][0], local[pold][0][1], knot, local[pold][0][2]};
        vector<vector<double>> nlocal {nknot, local[pold][1], local[pold][2]};

        vector<double> oldknot {local[pold][0][1], knot, local[pold][0][2], local[pold][0][3]};
        local[pold][0] = oldknot;

        local.insert(local.begin() + pold, nlocal);
        break;
    }
    case 1:
    {
        for (int i = indexes.size()-1; i > 0; i--) {
            if ((y < indexes[i][1])) {
                if ((x == indexes[i][0] && z == indexes[i][2])) {
                
                pold = i;
                }
            }
        }

        bet = (knot - local[pold][1][0]) / (local[pold][1][2] - local[pold][1][0]) * betas[pold];
        betas[pold] = (local[pold][1][3] - knot) / (local[pold][1][3] - local[pold][1][1]) * betas[pold];

        vector<double> nknot {local[pold][1][0], local[pold][1][1], knot, local[pold][1][2]};
        vector<vector<double>> nlocal {local[pold][0], nknot, local[pold][2]};

        vector<double> oldknot {local[pold][1][1], knot, local[pold][1][2], local[pold][1][3]};
        local[pold][1] = oldknot;

        local.insert(local.begin() + pold, nlocal);


        break;
    }
    case 2:
    {
        for (int i = indexes.size()-1; i > 0; i--) {
            if ((z < indexes[i][2])) {
                if ((x == indexes[i][0] && y == indexes[i][1])) {
                pold = i;
                }
            }
        }

        bet = (knot - local[pold][2][0]) / (local[pold][2][1] - local[pold][2][0]) * betas[pold];
        betas[pold] = (local[pold][2][2] - knot) / (local[pold][2][2] - local[pold][2][1]) * betas[pold];

        vector<double> nknot {local[pold][2][0], knot, local[pold][2][1]};
        vector<vector<double>> nlocal {local[pold][0], local[pold][1], nknot};

        vector<double> oldknot {local[pold][2][1], knot, local[pold][2][2]};
        local[pold][2] = oldknot;

        local.insert(local.begin() + pold, nlocal);        

        break;
        }
    }
    

    vector<double> p = pw[pold];
    vector<double> index {x, y, z};
    betas.insert(betas.begin() + pold, bet);
    
    pw.insert(pw.begin() + pold, p);
    indexes.insert(indexes.begin() + pold, index);

    auto *Pw = new vector<vector<double>>(pw);
    auto *locals = new vector<vector<vector<double>>>(local);
    auto *indx = new vector<vector<double>>(indexes);
    auto *beta = new vector<double>(betas);
    
    this->Pw.reset(Pw);
    this->local_knots.reset(locals);
    this->indexes.reset(indx);
    this->betas.reset(beta);

}

/* Weighted Control points accessor
* Output: a vector of vector of vector of  weighted control points 
*/
vector<vector<double>> tsplineVol::getWCP() {
    return *this->Pw;
}

/* Betas accessor
* Output: a vector of betas
*/
vector<double> tsplineVol::getBetas() {
    return *this->betas;
}

/* Anchor indexes accessor
* Output: a vector of betas
*/
vector<vector<double>> tsplineVol::getIndexes() {
    return *this->indexes;
}

/* Evaluation output accessor
* Output: a vector of evaluated points
*/
vector<vector<double>> tsplineVol::getRes() {
    return *this->out;
}

/* Evaluate volume given steps; from 0 to 1
* Output: populate out object
*/
void tsplineVol::evaluate(double stepu, double stepv, double steph) {
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
    
    // vector<vector<vector<double>>>
    // len(P)    3    degree + 2   
    auto kns = this->getKnots();

    auto betas = this->getBetas();

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
    
    auto points = this->getWCP();

    // vector<vector<double>>
    // len(P)   3
    auto index = this->getIndexes();


    for (auto&hst : hsteps) {
        for (auto&vst : vsteps) {
            for (auto&ust : usteps) {

                vector<double> Sw {0., 0., 0., 0.};
                vector<double> S {0., 0., 0.};

                for (int i = 0; i < points.size(); i++) {
                    // oneBasisFun is problematic, pls fix
                    double Nip = BasisFun(ust, kns[i][0], dv);
                    double Niq = BasisFun(vst, kns[i][1], du);
                    double Nir = BasisFun(hst, kns[i][2], dr);

                    for (int j = 0; j < points[0].size(); j++) {
                        Sw[j] += Nip * Niq * Nir * betas[i] * points[i][j];
                    }
                }
                for (int j = 0; j < points[0].size()-1; j++) {
                    S[j] = Sw[j]/Sw[3];
                }
                out->push_back(S);
            }
        }
    }

    this->out.reset(out);
}

void tsplineVol::write_csvMeta(string filename) {
    ofstream file;

    file.open(filename);

    auto CP = this->getWCP();
    auto indexes = this->getIndexes();
    auto betas = this->getBetas();
    auto locals = this->getKnots();

    file << "x, y, z, anchorx, anchory, anchorz, beta, local knot" << endl;

    for (int i = 0; i < CP.size(); i++) {
        file << CP[i][0] << "," << CP[i][1] << "," << CP[i][2] << ",";
        file << indexes[i][0] << "," << indexes[i][1] << "," << indexes[i][2] << ",";
        file << betas[i] << ",";
        for (int j = 0; j < locals[i].size()-1; j++) {
            file << "[";
            for (int k = 0; k < locals[i][j].size(); k++) {
                file << locals[i][j][k] << " ";
            }
            file << "]";
        }
        file << "," << endl;
        
    }

    file.close();
}
/* Export evaluated volume to vtk
*/
void tsplineVol::write_vtk(string filename) 
{
    ofstream file;


    file.open(filename);

    auto out = this->getRes();
    

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

vector<double> cart2par(vector<double> cor) {
    double r = sqrt(pow(cor[0], 2) + pow(cor[1], 2) + pow(cor[2], 2));
    double phi = acos(cor[1]/r);
    double theta = atan2(cor[1], cor[0]);

    // cout<< phi << endl;
    // cout << theta << endl;

    double norm_phi = (phi * 180/M_PI) / 180;
    double norm_theta = (theta * 180/M_PI) / 360;
    double norm_r = r;

    // cout<< norm_phi << endl;
    // cout << norm_theta << endl;

    return vector<double> {norm_r, norm_phi, norm_theta};
}

/* oneBasisFun Helper Function - NURBS Book A2.4
** Returns the value N_ip
** Input: span index i, degree p, parameter space value u, knot vector U
** Output: N_ip
*/
double BasisFun(double u, vector<double> U, int p) {
    double Nip;

    if (u < U[0] || u >= U[p+1]) {
        Nip = 0.;
        cout << "here!" << endl;
        return Nip;
    }
    
    vector<double> N(p+1);

    for (int i = 0; i <= p; i++) {
        if (u >= U[i] && u < U[i+1]) {
            N[i] = 1.;
        } else {
            N[i] = 0.;
        }
    }
    double saved;
    double Uleft;
    double Uright;
    double temp;
    for (int i = 1; i <= p+1; i++) {
        if (N[0] == 0.0) {
            saved = 0.0;
        } else {
            saved = ((u-U[0])*N[0])/(U[i]-U[0]);
        }

        for (int j = 0; j <= p-i; j++) {
            Uleft = U[j+1];
            Uright = U[j+i+1];

            if (N[j+1] == 0.0) {
                N[j] = saved;
                saved = 0.0;
            } else {
                temp = N[j+1]/(Uright - Uleft);
                N[j] = saved + (Uright-u)*temp;
                saved = (u - Uleft) * temp;
            }
        }
    }

    Nip = N[0];

    // cout << Nip << endl;

    return Nip;

}