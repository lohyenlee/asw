//=======================================================================
// Demonstration of Auxiliary-Spin Wolff (ASW) Algorithm
// Yen Lee Loh, 2022-7-1
//
// Demo 2: Graphical animation using GLFW+GLEW+GLM+OpenGL
// Usage:  make demo2gui && ./demo2gui
//=======================================================================
#include <cmath>      // for M_PI (warning: use fabs, not abs)
#include <string>
#include <vector>
#include <iostream>		// for cout
#include <iomanip>		// for setf
#include <sstream>		// for ostringstream
#include <fstream>		// for ofstream and ifstream
#include <regex>

typedef float flo;           // some of my files require this to be specified beforehand
typedef flo Coupling, Field; // some of my files require this to be specified beforehand

#include "utils.hh" // for homemade utils
#include "spin.hh"
#include "graph.hh"
#include "algo.hh"
#include "gui.hh"   // for graphics based on GLEW+GLFW+GLM

using std::cerr; using std::string;

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// GLOBAL VARIABLES RELATING TO SIMULATION 
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
enum {ALGO_MET, ALGO_GIB, ALGO_DST, ALGO_ASW, ALGO_GCA};
string algoNames[] = {"MET", "GIB", "DST", "ASW", "GCA"};
flo T;
flo h;
int algo;
bool shouldDrawBonds;

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// TERMINAL CODE
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
string fmt (string s) {
  return std::regex_replace (s, std::regex ("<\\|(.*?)\\|>"), "\x1b[92;40;1m$1\x1b[97;22;49m");
}
void printInfo() {
  cerr << "\x1b[F\x1b[K"
  << "alg = " << "\x1b[96;40;1m" << algoNames[algo] << "\x1b[97;22;49m"
  << "\tT = " << "\x1b[96;40;1m" << T << "\x1b[97;22;49m"
  << "\th = " << "\x1b[96;40;1m" << h << "\x1b[97;22;49m"
  << '\n';
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// GUI CODE SPECIFIC TO THE PRESENT APPLICATION
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
GUI gui; 
Model model1; // static objects in scene (coordinate axis tripod, unit cell)
Model model2; // dynamic objects in scene (spins)
Model model3; // static objects in scene (bonds)

//=============================================================================
// graph.hh and spin.hh define graph::vec3 and spin::vec3
// to avoid dependency on glm::vec3.
// Here we provide means of casting graph::vec3 and spin::vec3 to glm::vec3.
glm::vec3 cast (graph::vec3 r) { return glm::vec3 (r.x, r.y, r.z); }
glm::vec3 cast (spin::vec3 r)  { return glm::vec3 (r.x, r.y, r.z); }


template<class Graph>
void guiInit (Graph &graph) {
  //======== Set up the GUI object, which will handle most of the rendering stuff
  gui.setup (900, 600); // Call GL and GLFW window setup commands
  //======== Find center of mass of system
  GLfloat xmin=1.e30,ymin=1.e30,zmin=1.e30,xmax=-1.e30,ymax=-1.e30,zmax=-1.e30;
  glm::vec3 rcom {0,0,0};
  for (Vertex v=0; v<graph.vertices(); ++v) rcom += cast (graph.position (v));
  rcom /= graph.vertices(); 
  GLfloat rsys = 0;
  for (Vertex v=0; v<graph.vertices(); ++v) {
    vec3 r = cast (graph.position (v));
    rsys = std::max (rsys, glm::distance(r, rcom));
  }
  gui.rcom = rcom;
  gui.rsys = rsys;
  gui.zoomOutDistance = 0.75*std::max( rsys/tan(gui.camFoV/2), rsys/tan(gui.camFoV/2)/float(gui.windowWidth)*gui.windowHeight);
  //======== Set initial camera orientation
  gui.camOrient = glm::lookAt (vec3(0,0,0), vec3(0,3,-1), vec3(0,0,1)); // start,forward,up
  //======== Set up lighting
  gui.lightPos   = 100.f * vec3 (0,-1,1);   // light position in world coordinates
  gui.lightColor = 4.f *sqr(100.f) * vec3 (1,1,1);
  //======== Draw axes  
  const GLfloat l = 3.0; // length of axis symbols
  const GLfloat r = 1.0; // radius of axis cones
  const GLfloat rShaft = 0.5;
  model1.color (1,0,0); model1.tube (0,0,0, l,0,0, rShaft,rShaft, 12);  model1.cone (l,0,0, l+r,0,0, r, 12); 
  model1.color (0,1,0); model1.cone (0,0,0, 0,l,0, r, 12);
  model1.color (0,1,1); model1.tube (0,0,0, 0,0,l, rShaft,rShaft, 12);  model1.cone (0,0,l, 0,0,l+r, r, 12); model1.uploadToGPU ();
  //======== Draw bonds
  const GLfloat rBond = 0.1;
  model3.color (.4,.4,.4);
  for (Edge e=0; e<graph.edges(); ++e) {
    Vertex v0 = graph.startVertex(e);
    Vertex v1 = graph.endVertex(e);
    if (v0>=v1) continue;
    vec3 r0 = cast (graph.position(v0));
    vec3 r1 = cast (graph.position(v1));
    if (distance(r0,r1) > 5.0) continue; // don't draw bonds longer than a certain cutoff
    model3.tube (r0.x, r0.y, r0.z, r1.x, r1.y, r1.z, rBond,rBond, 12);
  }
  model3.uploadToGPU ();
}

template<class Graph,class Spin>
void guiUpdate (Graph &graph, std::vector<Spin> &spinArray) {
  //======== Re-render mobiles to buffers on CPU; transfer data to vertexXXX buffers on GPU
  const GLfloat l = 0.5; // length of arrowhead
  const GLfloat rArrow = 0.5; // rArrow of arrowhead
  const GLfloat rShaft = 0.25; // rArrow of arrow rShaft
  const GLint segments = 24;  //6 is ugly, 24 is nice
  model2.clear ();       
  for (Vertex v=0; v<graph.vertices(); ++v) {
    vec3 r = cast(graph.position (v)); flo x=r.x, y=r.y, z=r.z;
    vec3 s = cast(  spinArray[v].direction()  );
    model2.color (.5+.5*s.x, .5+.5*s.y, .5+.5*s.z); // up is blue, down is red
    model2.cone (x, y, z, x+s.x*l, y+s.y*l, z+s.z*l,    rArrow, segments);
    model2.tube (x, y, z, x-s.x*l, y-s.y*l, z-s.z*l,    rShaft, rShaft, segments);
    model2.align (0,0,0, -s.x, -s.y, -s.z); model2.trans (x-s.x*l, y-s.y*l, z-s.z*l);
    model2.disk (rShaft, segments, 0, true, false);
  }
  model2.uploadToGPU ();
  //======== Tell GPU to render from GPU's vertexXXX buffers to the framebuffer
  gui.clear ();
  gui.render (model1);		
  gui.render (model2);
  if (shouldDrawBonds) gui.render (model3);
  //======== Finally swap the offscreen buffer with the onscreen buffer
  gui.swapBuffers ();
  //======== Handle mouse and keyboard input
  glfwPollEvents ();
  gui.processTrackballControls ();
  if (gui.pressed(GLFW_KEY_T)) {
    if (gui.pressed(GLFW_KEY_MINUS) || gui.pressed(GLFW_KEY_KP_SUBTRACT)) {changeParameter (T, -1e-16, -0.005); printInfo();}
    if (gui.pressed(GLFW_KEY_EQUAL) || gui.pressed(GLFW_KEY_KP_ADD))      {changeParameter (T, +1e-16, +0.005); printInfo();}    
  } 
  else if (gui.pressed(GLFW_KEY_H)) {
    if (gui.pressed(GLFW_KEY_MINUS) || gui.pressed(GLFW_KEY_KP_SUBTRACT)) {changeParameter (h, -0.0001, -0.005); printInfo();}
    if (gui.pressed(GLFW_KEY_EQUAL) || gui.pressed(GLFW_KEY_KP_ADD))      {changeParameter (h, +0.0001, +0.005); printInfo();}    
  } 
  else if (gui.pressed(GLFW_KEY_A)) {algo = ALGO_ASW; printInfo();}
  else if (gui.pressed(GLFW_KEY_D)) {algo = ALGO_DST; printInfo();}
  else if (gui.pressed(GLFW_KEY_M)) {algo = ALGO_MET; printInfo();}
  else if (gui.pressed(GLFW_KEY_G)) {algo = ALGO_GCA; printInfo();}
  else if (gui.pressed(GLFW_KEY_B)) {algo = ALGO_GIB; printInfo();}
  if (gui.shouldExit) exit(0);				
}



//=======================================================================
//  MAIN
//=======================================================================
int main (int argc, char** argv) {
  cerr << fmt(
    "\n======== KEYBOARD CONTROLS =========================\n"
    "Camera:      <|Up|> <|Down|> <|Left|> <|Right|> (with <|Shift|> to zoom) \n"
    "Temperature: <|T+|> <|T-|>   (press two keys simultaneously) \n"
    "Field:       <|h+|> <|h-|>   (press two keys simultaneously) \n"
    "Algorithm:   <|A|>SW <|D|>ST <|M|>etropolis <|G|>CA Gi<|b|>bs \n"
    "Quit:        <|Escape|>  \n"
    "====================================================\n");
  
  int xmax = 32;
  int ymax = 32;
  int zmax = 4;
  //SquareLatticeGraph graph (xmax, ymax);
  //TriangularLatticeGraph graph (xmax, ymax);                    // Use this instead if desired
  CubicLatticeGraph graph (xmax, ymax, zmax, 1.0f, 1.0f, 4.0f); // Use this instead if desired
  
  UniformVertexProperty field (0.00); 
  
  //UniformEdgeProperty coupling (1.00);
  AnisotropicEdgeProperty coupling {0, 0, 0, 0, 10., 10.};          // Use this instead if desired
  
  std::vector<IsingSpin> spinArray (graph.vertices());      // Ising spins
  //std::vector<ClockSpin<6>> spinArray (graph.vertices()); // 6-state clock model
  //std::vector<XYSpin> spinArray (graph.vertices());       // XY spins as unit 2-vectors
  //std::vector<HeisenSpin> spinArray (graph.vertices());   // Heisenberg spins as unit 3-vectors
  
  MetropolisSimulator   metropolisSimulator;
  GibbsSimulator        gibbsSimulator;
  DSTSimulator          dstSimulator;
  ASWSimulator          aswSimulator;
  GCASimulator          gcaSimulator;
  
  T = 2.27; 
  algo = ALGO_ASW; 
  shouldDrawBonds = true;
  printInfo();
  
  guiInit (graph);
  while (true) {
    field.set (h);
    
    if (algo==ALGO_MET) {
      for (int s=0; s<xmax*ymax; ++s) {
        metropolisSimulator.step (spinArray, graph, coupling, field, T);
      }
    } else if (algo==ALGO_GIB) {
      for (int s=0; s<xmax*ymax; ++s) {
        gibbsSimulator.step (spinArray, graph, coupling, field, T);
      }
    } else if (algo==ALGO_DST) {
      dstSimulator.step (spinArray, graph, coupling, field, T);
    } else if (algo==ALGO_ASW) {
      aswSimulator.step (spinArray, graph, coupling, field, T);
      aswSimulator.rectify (spinArray, graph); // have to do this before rendering
    } else if (algo==ALGO_GCA) {
      gcaSimulator.step (spinArray, graph, coupling, field, T);
    }
    
    guiUpdate (graph, spinArray);
  };
  
  return 0;
}
