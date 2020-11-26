// Sotiris Gkoulmaris - sg3219
// spline-tool

// NURBS volume derived class cpp file. Contains defitions of volume
// generators, operators and samplers.

#pragma once
#include <cmath>
#include <iostream>

#include "spline.h"
#include "nurbsVol.h"

struct vertex {
    double beta = 1.;
    vector<double> index;
    vector<double> P;
    vector<double> Pw;
    vector<vector<double>> knots;

    vertex *left = nullptr;
    vertex *right = nullptr;
    vertex *up = nullptr;
    vertex *down = nullptr;
    vertex *in = nullptr;
    vertex *out = nullptr;

};


class tsplineVol: public spline
{
public:

    unique_ptr<vertex> *head;                // head of the vertexes linked structure
    unique_ptr<vector<vector<double>>> Pw;  // weighted control points of the volume
    unique_ptr<vector<vector<double>>> P;   // control points of the volume
    unique_ptr<vector<vector<vector<double>>>> local_knots;   // t-mesh in V parametric direction
    unique_ptr<vector<vector<double>>>  indexes;   // vector of anchor indexes
    unique_ptr<vector<double>> betas;       // vector of betas
    int p;                                  // degree p
    int q;                                  // degree q
    int r;                                  // degree r
    int n;                                  // helper variable for evaluation
    int m;                                  // helper variable for evaluation
    int l;                                  // helper variable for evaluation
    int thicc;                              // helper variable for export to vtk
    int up;                                 // helper variable for export to vtk
    int right;                              // helper variable for export to vtk

    /* Constructor used to define a T-spline volume, given the parametric representation of a
    ** volume.
    ** Input: vol nurbsVol
    ** Output: creates a T-Spline Volume object
    */
    tsplineVol(nurbsVol &vol);

    /* Class destructor - since we are using unique_ptr, it can be empty
    */
    ~tsplineVol();

    /* Evaluate volume given steps; from 0 to 1
    * Output: populate out object
    */
    void evaluate(double, double, double);

    /* Local knot refinement algorithm, adapted from S-Spline
    * Input: anchor indexes x, y, z, knot value, direction
    * Output: update Pw, betas, indexes, local knots with the new added points
    */
    void refine(double, double, double, double, int);

    /* Local knot vectors accessor
    * Output: a vector of knot values
    */
    vector<vector<vector<double>>> getKnots();

    /* Weighted Control points accessor
    * Output: a vector of vector of vector of  weighted control points 
    */
    vector<vector<double>> getWCP();

    /* Betas accessor
    * Output: a vector of betas
    */
    vector<double> getBetas();

    /* Anchor indexes accessor
    * Output: a vector of betas
    */
    vector<vector<double>> getIndexes();


    /* Evaluation output accessor
    * Output: a vector of evaluated points
    */
    vector<vector<double>> getRes();

    /* Declaration of virtual function to export out to csv
    */
    void write_csv(string);

    /* Export evaluated volume to vtk
    */
    void write_vtk(string);

    /* Export data to csv
    */
    void write_csvMeta(string);
};

/* Helper function for finding anchor index
*  Output: Vector of anchor indexes
*/
vector<double> findex(int n, int p);