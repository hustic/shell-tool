// Sotiris Gkoulmaris - sg3219
// spline-tool

// Abstract Class spline header file. Contains definitions of helper functions
// Note: All Spline entities will be derived by this base class
#pragma once
#include <memory>
#include <vector>


using namespace std;

// Abstract class definition
class spline
{
public:

    unique_ptr<vector<vector<double>>> out;   // storage for the evaluated points

    // Empty Constructors and Destructors, as this class does not store anything
    spline() {}
    ~spline() {}

    // virtual function write_csv to be defined by all derived classes.
    virtual void write_csv(string) = 0;

    //virtual void evaluate() = 0;

};


/* getIntersect Helper Function
** Calculates the intersection point of two lines on a specified 2D plane
** Input: points (p0, p2), direction vectors (t0, t2)
** Output: vector that contains the coordinates of the point of intersection
*/
vector<double> getIntersect(vector<double> p0, vector<double> t0, vector<double> p2, vector<double> t2);

/* fspan Helper Function - NURBS Book Algorithm A2.1
** Calculates the index of the span where a control point is located
** Input: degree p, parameter space value u, knot vector U
** Output: span index
*/
int fspan(int p, double u, vector<double> U);


/* oneBasisFun Helper Function - NURBS Book A2.4
** Returns the value N_ip
** Input: span index i, degree p, parameter space value u, knot vector U
** Output: N_ip
*/
double oneBasisFun(int, double, int, vector<double>);

/* basisFun Helper Function - NURBS Book Algorithm A2.2
** Computes all the non-vanishing basis function and returns them in an array N[p]
** Input: Index i, parameter space value u, degree p, knot vector U
** Output: Array of basis function N
*/
vector<double> basisFun(int i, double u, int p, vector<double> U);