//=======================================================================
// Demonstration of Auxiliary-Spin Wolff (ASW) Algorithm
// Yen Lee Loh, 2022-7-1
//
// This program is meant to be called by demo3traj.py and demo4launch.py.
//=======================================================================
#include <iostream>		// for cout
#include <iomanip>		// for setf
#include <string>
#include <vector>
#include <regex>

typedef float flo;           // alternatively: typedef double flo
typedef flo Coupling, Field; // some of my files require this to be specified beforehand

#include "utils.hh"
#include "spin.hh"
#include "graph.hh"
#include "algo.hh"

using std::cerr; using std::setw; using std::string;
using std::fixed; using std::setw; using std::setprecision; using std::ofstream;

//=======================================================================
// Total magnetic moment in direction of field M = sum_i s^X_i
// Total nearest-neighbor correlation W = sum_<ij> S_i.S_j
//=======================================================================
template<class Spin,class Graph>
auto computeMoment (std::vector<Spin> &spinArray, Graph &graph) {
  int vmax = graph.vertices();
  auto M = 0 * spinArray[0].moment(); // magnetization in 0th direction
  for (Vertex v=0; v<vmax; ++v)	M += spinArray[v].moment();
  return M; // M may be int, float, or double depending on Spin
}
template<class Spin,class Graph>
auto computeNeighborCorrelation (std::vector<Spin> &spinArray, Graph &graph) {
  int emax = graph.edges();
  auto W = 0 * spinArray[0].interact (spinArray[0]);  // hack to determine type for W
  for (int e=0; e<emax; ++e) { // run over both forward and backward edges
    int v0 = graph.startVertex(e);
    int v1 = graph.endVertex(e);
    if (v0 < v1) W += spinArray[v0].interact (spinArray[v1]); 	// don't double-count
  }
  return W; // W may be int, float, or double depending on Spin
}
//=======================================================================
// With the ASW algorithm, the physical spin component in the field direction is S_1 = (s_i dot s_aux)
//=======================================================================
template<class Spin,class Graph>
auto computePhysicalMoment (std::vector<Spin> &spinArray, Graph &graph) {
  int vmax = graph.vertices();
  int vAux = vmax;                    // index of auxiliary spin
  auto M = 0 * spinArray[0].moment(); // magnetization in 0th direction
  for (Vertex v=0; v<vmax; ++v)	M += spinArray[v].interact (spinArray[vAux]);
  return M; // M may be int, float, or double depending on Spin
}

//=======================================================================
//  MAIN
//=======================================================================
int main (int argc, char** argv) {
  //======== Set default parameters
  int xmax = 16;
  int ymax = 16;  
  flo T = 2.27;
  flo h = 0.00;
  int steps = 10000;

  //======== Read parameter file (overrides defaults)
  string fileName ("pars.dat");
  readPar (fileName, "xmax",        xmax);
  readPar (fileName, "ymax",        ymax);
  readPar (fileName, "temperature", T);
	readPar (fileName, "field",       h);
  readPar (fileName, "steps",       steps);
  
  //======== Set up shop
  SquareLatticeGraph     graph (xmax, ymax);
  UniformVertexProperty  field (h);
  UniformEdgeProperty    coupling (1.00);
  std::vector<IsingSpin> spinArray (graph.vertices());
  ASWSimulator           aswSimulator;
  
  //======== Main loop
  ofstream fout ("history.dat");
  fout << setprecision (6);
  for (int step=0; step<steps; ++step) {
    //-------- Perform one ASW cluster update
    aswSimulator.step (spinArray, graph, coupling, field, T);
    //-------- Compute estimators
    //aswSimulator.rectify (spinArray, graph);
    //auto M =  computeMoment (spinArray, graph);
    auto M =  computePhysicalMoment (spinArray, graph);
    auto U = -computeNeighborCorrelation (spinArray, graph);
    fout << setw(14) << M << setw(14) << U << '\n';
  }
  cerr << "Simulation completed normally!\n";
  return 0;
}
