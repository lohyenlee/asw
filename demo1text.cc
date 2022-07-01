//=======================================================================
// Demonstration of Auxiliary-Spin Wolff (ASW) Algorithm
// Yen Lee Loh, 2022-7-1
//
// Demo 1: Terminal-based animation
// Usage:  make demo1text && ./demo1text
//=======================================================================
#include <iostream>		// for cerr
#include <iomanip>		// for setprecision
#include <string>
#include <vector>
#include <regex>

typedef float flo;           // alternatively: typedef double flo
typedef float flo;           // alternatively: typedef double flo
typedef flo Coupling, Field; // these typedefs must be done before 

#include "utils.hh" // for homemade utils
#include "spin.hh"
#include "graph.hh"
#include "algo.hh"

using std::cerr; using std::string;

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

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// TERMINAL ANIMATION CODE
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
string CLEAR_SCREEN = "\x1b[2J", CLEAR_LINE = "\x1b[2K",  
HIDE_CURSOR = "\x1b[?25l", SHOW_CURSOR = "\x1b[?25h",
DULL="\x1b[33;1m", BOLD="\x1b[97;1m";
string toupper (string s) {for (auto& c: s) c=toupper(c); return s;}
string fmt (string s) {
  return std::regex_replace (s, std::regex ("<(.*?)>"), "\x1b[92;40;1m$1\x1b[97;22;49m");
}

class Terminal {
public:
  void moveCursor (int row, int col) { cerr << "\x1b[" << row << ";" << col << "H"; }
  void resize (int rows, int cols) { cerr << "\x1b[8;" << rows << ";" << cols << ";t"; }
  
  template<class Spin, class Graph>
  void display (Graph &graph, std::vector<Spin> &spinArray) {
    cerr << HIDE_CURSOR;
    for (Vertex v=0; v<graph.vertices(); ++v) {
      auto r = graph.position(v);
      int x = round(r.x); if (x<1) continue;
      int y = round(r.y); if (y<0) continue;
      moveCursor (1+y, 0+x);
      cerr << spinArray[v].asciiArt();
    }
    cerr << "\x1b[33;0m" << SHOW_CURSOR;
  }
};



//=======================================================================
//  MAIN
//=======================================================================
int main (int argc, char** argv) {
  //======== Define lattice to simulate (graph)
  int xmax = 79; // change this and recompile if desired
  int ymax = 18;  
  SquareLatticeGraph graph (xmax, ymax);
  //TriangularLatticeGraph graph (xmax, ymax, 1.0f, 1.0/0.866f);
  
  //======== Define magnetic fields at each site on lattice (vertex of graph)
  UniformVertexProperty field (0.00);
  
  //======== Define couplings along each bond on lattice (edge of graph)
  UniformEdgeProperty coupling (1.00);            // J=1
  //AnisotropicEdgeProperty coupling {0.5, 2.0};  // Jx=0.5, Jy=2.0
  
  //======== Choose type of spin (change and recompile if desired)
  std::vector<IsingSpin> spinArray (graph.vertices());       // Ising spins
  //std::vector<PottsSpin<3>> spinArray (graph.vertices());  // 3-state Potts model
  //std::vector<ClockSpin<6>> spinArray (graph.vertices());  // 6-state clock model
  //std::vector<XYSpin> spinArray (graph.vertices());        // XY spins (unit 2-vectors)
  //std::vector<HeisenSpin> spinArray (graph.vertices());    // Heisenberg spins (unit 3-vectors)
  
  //======== Initialize simulators
  MetropolisSimulator  metropolisSimulator;
  DSTSimulator         dstSimulator;
  ASWSimulator         aswSimulator;

  //======== Initialize terminal renderer
  Terminal terminal;
  terminal.resize (ymax+7, xmax+1);
  cerr << CLEAR_SCREEN;
  
  //======== Initialize simulation variables
  flo T = 2.27;
  flo h = 0.00;
  int speed = 1;
  string algorithm = "MET";
  string sCommand = " ", sInput;
  
  //======== Begin main program loop
  while (true) {
    field.set (h);
    
    for (int iter=0; iter<speed; ++iter) {
      if (algorithm=="MET") {
        for (int s=0; s<xmax*ymax; ++s) {
          metropolisSimulator.step (spinArray, graph, coupling, field, T);
        }
      } else if (algorithm=="DST") {
        dstSimulator.step (spinArray, graph, coupling, field, T);
      } else if (algorithm=="ASW") {
        aswSimulator.step (spinArray, graph, coupling, field, T);
        aswSimulator.rectify (spinArray, graph);
      }
    }
    
    auto M =  computeMoment (spinArray, graph);
    auto U = -computeNeighborCorrelation (spinArray, graph);
    
    terminal.display (graph, spinArray);
    
    terminal.moveCursor (ymax+1,1);
    cerr << CLEAR_LINE << '\n'
    << CLEAR_LINE 
    << DULL << " alg=" << BOLD << algorithm 
    << DULL << "  numSweeps=" << BOLD << speed << std::fixed << std::setprecision(2) 
    << DULL << "  T=" << BOLD << T << std::showpos 
    << DULL << "  h=" << BOLD << h
    << DULL << "  M=" << BOLD << M/flo(xmax*ymax)
    << DULL << "  U=" << BOLD << U/flo(xmax*ymax) << std::noshowpos << '\n';
    terminal.moveCursor (ymax+6,1); 
    cerr << CLEAR_LINE << '\n' << CLEAR_LINE;
    terminal.moveCursor (ymax+3,1); 
    cerr
    << CLEAR_LINE
    << fmt(" <A>SW <D>ST <M>etropolis <F>aster <S>lower\n")
    << CLEAR_LINE
    << fmt(" <[]>Change T <,.>Change h <Space>=Run <Enter>=RepeatLastCommand <Q>uit\n")
    << BOLD << CLEAR_LINE << " TYPE LETTER FOR COMMAND AND PRESS ENTER: " << "\x1b[0m ";
    
    std::getline (std::cin, sInput); 
    if (sInput!="") sCommand = sInput;
    sCommand = toupper (sCommand);
    if      (sCommand=="Q") break;
    else if (sCommand=="[") T -= 0.01;
    else if (sCommand=="]") T += 0.01; 
    else if (sCommand==",") {h -= 0.01; field.set (h);}
    else if (sCommand==".") {h += 0.01; field.set (h);}
    else if (sCommand=="A") {algorithm="ASW";}
    else if (sCommand=="D") {algorithm="DST";}
    else if (sCommand=="M") {algorithm="MET";}
    else if (sCommand=="F") {speed*=10;if (speed>1000) speed=1000;}    
    else if (sCommand=="S") {speed/=10;if (speed<1) speed=1;} 
  };
  return 0;
}
