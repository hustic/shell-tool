// Sotiris Gkoulmaris - sg3219
// spline-tool

// NURBS volume derived class cpp file. Contains defitions of volume
// generators, operators and samplers.

#pragma once
#include <cmath>
#include <iostream>

#include "spline.h"
#include "nurbsSur.h"


class nurbsVol: public spline
{
public:

    unique_ptr<vector<vector<vector<vector<double>>>>> Pw;  // weighted control points of the volume
    unique_ptr<vector<vector<vector<vector<double>>>>> P;   // control points
    unique_ptr<vector<double>> V;                           // knot vector V
    unique_ptr<vector<double>> U;                           // knot vector U
    unique_ptr<vector<double>> H;                           // knot vector H
    int p;                                                  // degree p
    int q;                                                  // degree q
    int r;                                                  // degree r
    int thicc;                                              // helper variable for export to vtk
    int up;                                                 // helper variable for export to vtk
    int right;                                              // helper variable for export to vtk

    /* Constructor used to define a volume, given the parametric representation of the
    * sphere and a curve for thickness
    * Input: sphere nurbsSur, nurbsCur
    * Output: creates a NURBS volume object (spherical shell)
    */
    nurbsVol(nurbsSur &sphere, nurbsCur &thicc);

    /* Class destructor - since we are using unique_ptr, it can be empty
    */
    ~nurbsVol();

    /* Control points accessor
    ** Output: a vector of vector of vector of control points 
    */
    vector<vector<vector<vector<double>>>> getCP();

    /* Weighted control points accessor
    * Output: a vector of vectors of vectors of control points 
    */
    vector<vector<vector<vector<double>>>> getWCP();

    /* Compute weighted control points and populate the class object Pw
    */
    void makePw();

    /* Global knot refinement algorithm, expanded from NURBS Book A5.5
    * Input: int direction, vector of knots to be added
    * Output: update Pw, U, V, H with new knots and control points
    */
    void refine(int, vector<double>);

    /* Evaluate volume given steps; from 0 to 1
    * Output: populate out object
    */
    void evaluate(double, double, double);


    /* Evaluate volume given vectors of steps; from 0 to 1
    * Output: populate out object
    */
    void evaluate(vector<double> usteps, vector<double> vsteps, vector<double> hsteps);

    /* Knot vector V accessor
    ** Output: a vector of knot values
    */
    vector<double> getKnotV();

    /* Knot vector U accessor
    * Output: a vector of knot values
    */
    vector<double> getKnotU();

    /* Knot vector H accessor
    * Output: a vector of knot values
    */
    vector<double> getKnotH();

    /* Evaluation output accessor
    * Output: a vector of evaluated points
    */
    vector<vector<double>> getRes();

    /* Declaration of virtual function to export out to csv
    */
    void write_csv(string);

    /* Export control points to csv
    */
    void writeCP_csv(string);

    /* Export evaluated volume to vtk
    */
    void write_vtk(string);
};

/* Helper function to find the knot values of the anchors
* Input: degree p, knot vector V
* Output: anchor knot values
*/
vector<double> getKnotAnchors(int, vector<double>);

