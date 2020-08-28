# Independent Research Project - Spline Tool

## Instalation

To make use of the tool you can just clone the repository.

## User instructions

The corresponding .cpp files need to included in a project, according to the needs of the user. The tool is to be compiled together with the code that makes use of it and any compiler can be used for this purpose. The only important think to note is that higher dimension entities depend on the lower ones. So, if you want to use the *nurbsSur* class, *nurbsCur* needs to also be included.

## Example usage

The following snippet summarizes the entire workflow of the tool.

```
#include "tsplineVol.cpp"    // include all the .cpp files
#include "nurbsVol.cpp"
#include "nurbsSur.cpp"
#include "nurbsCur.cpp"

vector<int> O {0, 0, 0};    // Define axis where the circles will be defined on.
vector<int> X {1, 0, 0};
vector<int> Y {0, 1, 0};
vector<int> Z {0, 0, 1};
double r = 1.;
double ths = 0.;
double the = 360.;

vector<vector<double>> pt { {1., 0  }, 
                            {.8, 0, }, 
                            {.5, 0, }};

vector<double> ut {0., 0., 0.5, 1., 1.};

auto* full = new nurbsCur(O, X, Y, r, ths, the); // generate full circle
auto* thic = new nurbsCur(Pt, Ut, 1);            // generate thickness curve

the = 180.;
auto* half = new nurbsCur(O, Z, X, r, ths, the); // generate half circle

auto* sphere = new nurbsSur(*full, *half);      // generate sphere
auto* shell = new nurbsVol(*sphere, *thic);     // generate shell

shell->makePw();                               // !!!!! Very important function. For all geometric operations (and evaluation) getPw() needs to be called.
                                               // !!!!! that also applies to curves and surfaces.
                                               
vector<double> knotsv {0.25, 0.25};           // define knot vector to be added

shell->refine(1, knotsv);                     // refine the shell in the latitude parametric direction (note: 0 is for longitude, 1 is for latitude, 2 for thickness)

auto* tspline = new tsplineVol(*shell);       // create T-Spline entity

tspline->refine(2.5, 2.5, 1.25, 0.25, 2);     // locally refine the T-Spline

tspline->evaluate(0.1, 0.05, 0.5);            // evaluate the entity

tspline->write_vtk("data.vtk");               // output to .vtk


```

The code has been extensively commented and will provide the information required to use the tool effectively.

## Things to Note about Refinement

Global refinement is rigorous and can properly accept any knots without breaking the spline. Local refinement is a different case however. One needs to be very carefull with which knots are to be added and where, else the spline will break. The paper from S-Splines provides guidelines and checks one can make to access whether or not a control point is refinement compatible.
