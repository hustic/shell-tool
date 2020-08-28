// Sotiris Gkoulmaris - sg3219
// spline-tool

// NURBS surface derived class cpp file. Contains defitions of surface
// generators, operators and samplers.

#pragma once
#include <cmath>
#include <iostream>

#include "spline.h"
#include "nurbsCur.h"


class nurbsSur: public spline
{
public:

    unique_ptr<vector<vector<vector<double>>>> Pw;  // control points of the surface
    unique_ptr<vector<double>> V;                   // knot vector V
    unique_ptr<vector<double>> U;                   // knot vector U
    int p;                                          // degree p
    int q;                                          // degree q

    /* Constructor used to define a surface, given the parametric representation of the
    * sphere.
    *  Input: half and full circle nurbsCurs
    *  Output: creates a NURBS surface object (sphere)
    */
    nurbsSur(nurbsCur &full, nurbsCur &half);

    /* Class destructor - since we are using unique_ptr, it can be empty
    */
    ~nurbsSur();

    /* Control points accessor
    *  Output: a vector of vector of vector of control points 
    */
    vector<vector<vector<double>>> getCP();

    /* Compute weighted control points and populate the class object Pw
    */
    void makePw();

    /* Evaluate surface given steps; from 0 to 1
    * Output: populate out object
    */
    void evaluate(double, double);

    /* Knot vector V accessor
    * Output: a vector of knot values
    */
    vector<double> getKnotV();

    /* Knot vector U accessor
    * Output: a vector of knot values
    */
    vector<double> getKnotU();

    /* Evaluation output accessor
    * Output: a vector of evaluated points
    */
    vector<vector<double>> getRes();

    /* Declaration of virtual function to export out to csv
    */
    void write_csv(string);
};