// Sotiris Gkoulmaris - sg3219
// spline-tool

// Abstract Class spline cpp file. Contains implementation of helper functions

#include <cmath>
#include <iostream>

#include "spline.h"

using namespace std;


/* getIntersect Helper Function
** Calculates the intersection point of two lines on a specified 2D plane
** Input: points (p0, p2), direction vectors (t0, t2)
** Output: vector that contains the coordinates of the point of intersection
*/
vector<double> getIntersect(vector<double> p0, vector<double> t0, vector<double> p2, vector<double> t2)
{
    double m0 = 0.;
    double epsilon = 0.00001;
    if (abs(t0[0]) < epsilon ) {
        if (t0[1] < 0) {
            m0 = -1000000.;
        } else {
            m0 = 1000000.;
        }
    } else {
        m0 = t0[1] / t0[0];
    }

    double m2 = 0.;
    if (abs(t2[0]) < epsilon) {
        if (t2[1] < 0) {
            m2 = -1000000.;
        } else {
            m2 = 1000000.;
        }
    } else {
        m2 = t2[1] / t2[0];
    }

    double b0 = 0.;
    double b2 = 0.;
    double xi = 0.;
    double yi = 0.;

    b0 = p0[1] - m0 * p0[0];
    b2 = p2[1] - m2 * p2[0];

    xi = (b0 - b2) / (m2 - m0);
    yi = m0 * xi + b0;

    vector<double> out {xi, yi};
    return out;
}

/* fspan Helper Function - NURBS Book Algorithm A2.1
** Calculates the index of the span where a control point is located
** Input: degree p, parameter space value u, knot vector U
** Output: span index
*/
int fspan(int p, double u, vector<double> U)
{
    int mid = 0;
    int m = U.size() - 1;
    int n = m - p - 1;

    if (u == U[n+1]) {
        return n;
    } 

    int low = p;
    int high = n + 1;
    mid = int(floor((low + high)/2.0));

    while ( (u < U[mid]) || (u >= U[mid+1])) {
        if (u < U[mid]) {
            high = mid;
        } else {
            low = mid;
        }
        mid = int(floor((low+high)/2));
    }

    return mid;
}

/* oneBasisFun Helper Function - NURBS Book A2.4
** Returns the value N_ip
** Input: span index i, degree p, parameter space value u, knot vector U
** Output: N_ip
*/
double oneBasisFun(int j, double u, int n, vector<double> U) {
    int p = U.size() - 2;
    double Nip;
    if ((j == 1 && u == U[0]) || (j == n && u == U[p+1])) {
        Nip = 1.;
        return Nip;
    }

    if (u < U[0] || u >= U[p+1]) {
        Nip = 0.;
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

    return Nip;

}

/* basisFun Helper Function - NURBS Book Algorithm A2.2
** Computes all the non-vanishing basis function and returns them in an array N[p]
** Input: Index i, parameter space value u, degree p, knot vector U
** Output: Array of basis function N
*/
vector<double> basisFun(int i, double u, int p, vector<double> U)
{
    vector<double> N(p+1);
    N[0] = 1.;

    vector<double> left(p+1);
    vector<double> right(p+1);

    double temp = 0.;

    double saved;
    for (int j = 1; j <= p; j++) {
        left[j] = u - U[i+1-j];
        right[j] = U[i+j] - u;
        
        saved = 0.;

        for (int r = 0; r < j; r++) {
            temp = N[r]/(right[r+1] + left[j-r]);
            N[r] = saved + right[r+1]*temp;
            saved = left[j-r]*temp;
        }
        N[j] = saved;
    }

    return N;
}