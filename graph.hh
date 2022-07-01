#ifndef GRAPH_HH
#define GRAPH_HH

#include <array> // for std::array

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// This C++ header file provides the following types and implements the following methods:
//   Vertex, Edge, Neighbor (all typedef'd to int)
//
//   Graph interface:
//     vertices()      
//     edges()
//     neighbors(v)
//     neighbor(v,n)    // used by Metropolis,Wolff,etc.
//     startVertex(e)   // used by Swendsen,Ziff
//     position(v,n)    // used by visualization code
//   Implementations:
//     SquareLatticeGraph
//     TriangularLatticeGraph
//     CubicLatticeGraph
//     GeneralGraph
//
//   VertexProperty interface:
//     flo operator() (v)
//   Implementations:
//     UniformVertexProperty
//     GeneralVertexProperty
//
//   EdgeProperty interface:
//      flo operator() (v,n)
//   Implementations:
//     UniformEdgeProperty
//     GeneralEdgeProperty
//
// If we implement Swendsen, then we need to implement operator(e) for EdgeProperty
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

typedef int Vertex;
typedef int Edge;
typedef int Neighbor;

// Struct for returning 3-vector position of each node, for visualization purposes
namespace graph {struct vec3 {float x,y,z;};};



//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//  GeneralGraph
//  Neighborlists and adjacency lists are stored in std::vector structures
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class GeneralGraph {
  typedef std::array<flo,3> array3;
  std::vector<std::vector<int>>				neighlist;
  std::vector<std::pair<int,int>>			edgelist;
  //std::vector<array3>                 positionlist;
  std::vector<graph::vec3>                 positionlist;
public:
  GeneralGraph () {}
  int addVertex () {
    neighlist.push_back (std::vector<int> () );
    return neighlist.size()-1; // index of item that was just added
  }
  void addEdge (int v0, int v1) {  // ADD A DIRECTED EDGE!
    neighlist[v0].push_back (v1);
    edgelist.push_back (std::pair<int,int> (v0, v1));
  }
  int vertices ()             { return neighlist.size();    }
  int edges ()                { return edgelist.size();     } // DIRECTED EDGES!
  int neighbors (int v)       { return neighlist[v].size(); }
  int neighbor (int v, int n) { return neighlist[v][n];     }
  int startVertex (int e)     { return edgelist[e].first;   }
  int endVertex (int e)       { return edgelist[e].second;  }
  
  //graph::vec3 position (int v) {return graph::vec3 {positionlist[v][0], positionlist[v][1], positionlist[v][2]};}
  graph::vec3 position (int v) {return positionlist[v];     }

  void addPosition (graph::vec3 pos) {  // should call this the same time as addVertex
    positionlist.push_back (pos);
  }
};



//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//  SquareLatticeGraph, TriangularLatticeGraph, CubicLatticeGraph
//  Connectivity information is computed on-the-fly using integer arithmetic
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class SquareLatticeGraph {
  int xmax,ymax;
  int xv (int v)              { return v%xmax;     }
  int yv (int v)              { return v/xmax;     }
  int vxy (int x, int y)      { return x + xmax*y; }
public:
  SquareLatticeGraph (int xmax, int ymax) : xmax(xmax), ymax(ymax) {}
  int vertices ()             { return xmax*ymax;       }
  int edges ()                { return xmax*ymax*4;     } // DIRECTED EDGES!
  int neighbors (int v)       { return 4;               }
  int neighbor (int v, int n) {
    const int dirs[][2] = {{1,0},{-1,0},{0,1},{0,-1}};
    int x = xv(v); int xdiff = dirs[n][0];
    int y = yv(v); int ydiff = dirs[n][1];
    int x2 = (x + xdiff + xmax) % xmax;
    int y2 = (y + ydiff + ymax) % ymax;
    return vxy (x2, y2);
  }

  graph::vec3 position (int v) {return graph::vec3 {1.0f*xv(v), 1.0f*yv(v), 0.0f};}
  
  int startVertex (int e) 		{ int v = e/4; return v; }
  int endVertex (int e)       { int v = e/4, n = e%4; return neighbor (v, n); }
  
  struct GeometricTransform { 	// inner class
    int xtrans,ytrans,g;
    GeometricTransform() {}
  };
  GeometricTransform randomTransform () {
    GeometricTransform gt;
    gt.xtrans = irand(xmax);
    gt.ytrans = irand(ymax);
    gt.g = irand(5);  //gt.g = 2; //inversion
    return gt;
  }
  int applyTransform (GeometricTransform &gt, int v) {
    int &xtrans=gt.xtrans, &ytrans=gt.ytrans, &g=gt.g; // aliases
    const int mat[][2][2] = {
      {{-1,0}, {0,1}},		// reflection along x
      {{1,0}, {0,-1}},		// reflection along y
      {{-1,0}, {0,-1}},		// inversion in x-y plane
      {{0,1},  {1,0}},    // reflection about the line y=x
      {{0,-1},  {-1,0}}   // reflection about the line y=-x
      };
    int x = xv(v) - xtrans;	// apply T^{-1}
    int y = yv(v) - ytrans;
    int xnew = mat[g][0][0] * x + mat[g][0][1] * y; // apply R
    int ynew = mat[g][1][0] * x + mat[g][1][1] * y;
    xnew = modulo (xnew + xtrans, xmax);  // apply T
    ynew = modulo (ynew + ytrans, ymax);
    return vxy (xnew, ynew);
  }
};


class TriangularLatticeGraph {
  int xmax,ymax;
  flo a,b;
  int xv (int v)              { return v%xmax;     }
  int yv (int v)              { return v/xmax;     }
  int vxy (int x, int y)      { return x + xmax*y; }
public:
  TriangularLatticeGraph (int xmax, int ymax, flo a=1, flo b=1) : xmax(xmax), ymax(ymax), a(a), b(b) {}
  int vertices ()             { return xmax*ymax;       }
  int edges ()                { return xmax*ymax*6;     } // DIRECTED EDGES!
  int neighbors (int v)       { return 6;               }
  int neighbor (int v, int n) {
    const int dirs[][2] = {{1,0},{-1,0},{0,1},{0,-1},{-1,-1},{1,1}};
    int x = xv(v); int xdiff = dirs[n][0];
    int y = yv(v); int ydiff = dirs[n][1];
    int x2 = (x + xdiff + xmax) % xmax;
    int y2 = (y + ydiff + ymax) % ymax;
    return vxy (x2, y2);
  }
  graph::vec3 position (int v) {return graph::vec3 { (xv(v) - 0.5f*yv(v))*a, (0.866025f*yv(v))*b, 0.0f};}

  int startVertex (int e) 		{ int v = e/6; return v; }
  int endVertex (int e)       { int v = e/6, n = e%6; return neighbor (v, n); }
};


class CubicLatticeGraph {
  int xmax,ymax,zmax;    // number of sites in each direction
  flo a,b,c;             // lattice constants
  int xv (int v)                 { return v%xmax;                } // fastest-varying index
  int yv (int v)                 { return (v/xmax)%ymax;         }
  int zv (int v)                 { return v/(xmax*ymax);         } // slowest-varying index
  int vxyz (int x, int y, int z) { return x + xmax*(y + ymax*z); }
public:
  CubicLatticeGraph (int xmax, int ymax, int zmax, flo a=1.0f, flo b=1.0f, flo c=1.0f) 
  : xmax(xmax), ymax(ymax), zmax(zmax), a(a), b(b), c(c) {}
  int vertices ()             { return xmax*ymax*zmax;   }
  int edges ()                { return xmax*ymax*zmax*6; } 
  int neighbors (int v)       { return 6;                }
  int neighbor (int v, int n) {
    const int dirs[][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
    int x = xv(v); int xdiff = dirs[n][0];
    int y = yv(v); int ydiff = dirs[n][1];
    int z = zv(v); int zdiff = dirs[n][2];
    int x2 = (x + xdiff + xmax) % xmax;
    int y2 = (y + ydiff + ymax) % ymax;
    int z2 = (z + zdiff + zmax) % zmax;
    return vxyz (x2, y2, z2);
  }
  graph::vec3 position (int v) {return graph::vec3 {1.0f*xv(v), 1.0f*yv(v), 1.0f*zv(v)};}

  int startVertex (int e) 		{ int v = e/6; return v; }
  int endVertex (int e)       { int v = e/6, n = e%6; return neighbor (v, n); }
  
  
  
  //============= NOT IMPLEMENTED!!!!!!!!!!!!!!
  struct GeometricTransform { 	// inner class
    int xtrans,ytrans,ztrans,g;
    GeometricTransform() {}
  };
  GeometricTransform randomTransform () {
    abort();
    GeometricTransform gt;
    gt.xtrans = irand(xmax);
    gt.ytrans = irand(ymax);
    gt.ztrans = irand(zmax);
    gt.g = irand(5);  //gt.g = 2; //inversion
    return gt;
  }
  int applyTransform (GeometricTransform &gt, int v) {
    abort();
  }
};




//=========================================================================
// class VertexProperty 
//=========================================================================
// UniformVertexProperty field;
// field(v) returns the same value for any vertex v
class UniformVertexProperty {
  flo val;
public:
  UniformVertexProperty (flo val) { this->val = val;}
  void set (flo val)              { this->val = val;}
  flo operator() (Vertex v)       { return val; } 
};

// GeneralVertexProperty field;
// field(v) returns a property of vertex v by retrieving it from an array
class GeneralVertexProperty {
  std::vector<flo> data;
public:
  flo operator() (Vertex v)       { return data[v]; }
};

// StaggeredField field;
// field(v) returns +1 or -1 depending on whether v is even or odd
class StaggeredVertexProperty {
public:
  flo operator() (Vertex v)       { return (v&1) ? -1.0 : +1.0; }
};


//=========================================================================
// class EdgeProperty
//=========================================================================
class UniformEdgeProperty {
  flo val;
public:
  UniformEdgeProperty (flo val)          { this->val = val;}
  flo operator() (Vertex v, Neighbor n)  { return val; } // same value for any edge
};

class AnisotropicEdgeProperty {
  std::vector<flo>	Jn;
public:
  AnisotropicEdgeProperty (std::initializer_list<flo> list) {
    Jn = list;
    //for (auto J:list) Jn.push_back (J);  // lousy way to copy a list
  }
  flo operator() (Vertex v, Neighbor n)  { return Jn[n]; } // coupling depends only on n, not on v
};

class GeneralEdgeProperty {
  std::vector<std::vector<flo>> Jvn;
public:
  GeneralEdgeProperty () {}
  flo operator() (Vertex v, Neighbor n)  { return Jvn[v][n]; } 
  
  void setNumVertices (int vmax) {
    Jvn.resize (vmax);
  }
  void addEdge (int v0, flo Jval) {  // ADD A DIRECTED EDGE!
    Jvn[v0].push_back (Jval);
  }
};


#endif
