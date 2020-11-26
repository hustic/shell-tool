// Sotiris Gkoulmaris - sg3219
// spline-tool

// NURBS curve derived class header file. Contains definitions of curve
// generators, operators and samplers.
#pragma once
#include <cmath>
#include <iostream>

#include "spline.h"

// Derived class nurbsCur
class nurbsCur: public spline
{
public:

    unique_ptr<vector<vector<double>>> Pw;  // Control points of the curve
    unique_ptr<vector<double>> V;           // knot vector of the curve
    int p;                                  // degree of the curve
    

    /* Constructor used to define a curve, given Control Points, Knot vector
    ** and degree. It is used primarily to define the thickness curve for the
    ** shell construction.
    ** Input: control points P, knot vector U, degree p
    ** Output: creates a NURBS curve object
    */
    nurbsCur(vector<vector<double>> *P, vector<double> *U, int p);

    /* Constructor used to define a circular curve. It is used to create the
    ** half and full circle parametric representation of the sphere.
    ** Input: origin O, axis vectors X and Y, radius r, degrees start and end
    ** Output: creates a circular NURBS curve object 
    */
    nurbsCur(vector<int>& O, vector<int>& X, vector<int>& Y, double r, double start, double end);

    /* Class destructor - since we are using unique_ptr, it can be empty
    */
    ~nurbsCur();

    /* Evaluate curve given a step to go from 0 to 1
    ** Output: populate out object
    */
    void evaluate(double step);

    /* Control points accessor
    ** Output: a vector of vector of control points 
    */
    vector<vector<double>> getCP();

    /* Knot vector accessor
    ** Output: a vector of knot values
    */
    vector<double> getKnot();

    /* Compute weighted control points and populate the class object Pw
    */
    void makePw();

    /* Declaration of virtual function to export out to csv
    */
    void write_csv(string);

    /* Evaluation output accessor
    ** Output: a vector of evaluated points
    */
    vector<vector<double>> getRes();
};