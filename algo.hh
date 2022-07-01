#ifndef ALGO_HH
#define ALGO_HH

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Yen Lee Loh   2022-6-16
//
// This C++ header file provides the following types and implements the following methods:
//   Simulator interface:
//     step()      
//     lastAcceptanceStatus()
//     lastClusterSize()
//   Implementations:
//     MetropolisSimulator
//     GibbsSimulator
//     WolffSimulator
//     DSTSimulator
//     ASWSimulator
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#include "utils.hh"
#include "graph.hh"
#include "spin.hh"

//=========================================================================
//  MetropolisSimulator (or Glauber if you prefer)
//=========================================================================
class MetropolisSimulator {
  bool status;
public:
  MetropolisSimulator () {}
  
  template<class Spin, class Graph, class EdgeProperty, class VertexProperty>
  void step (std::vector<Spin> &spinArray,
             Graph &graph, 
             EdgeProperty &coupling, 
             VertexProperty &field, 
             flo T = flo(1.0)
  ) { 
    Spin axis; axis.randomize();			   // choose flip axis (meaningless for Ising spins)
    Vertex v0 = irand(graph.vertices()); // choose spin to flip
    Spin s0   = spinArray[v0];           // current value of chosen spin
    Spin s0p  = s0.flip(axis);           // proposed value of chosen spin
    
    flo energyChange = 0.;               // initialize accumulator
    //-------- Calculate energy change due to magnetic field acting on chosen spin
    Field h            = field (v0);     // magnetic field at vertex v0
    flo energyOriginal = -s0.interact(h);         // original energy; usually E  = -dot(h,s0)
    flo energyProposed = -s0p.interact(h);        // proposed energy
    energyChange += energyProposed - energyOriginal;
    //-------- Calculate energy change due to coupling between chosen spin and each neighbor
    for (Neighbor n=0; n<graph.neighbors(v0); ++n) {
      Coupling J = coupling (v0, n);          // coupling from vertex v0 to its nth neighbor
      Vertex v1  = graph.neighbor (v0, n);        // nth neighbor of vertex v0
      Spin s1    = spinArray[v1];                 // value of spin at vertex v1
      flo energyOriginal = -J * s0.interact(s1);  // original bond energy; usually E = -J*(s0 dot s1)
      flo energyProposed = -J * s0p.interact(s1);
      energyChange += energyProposed - energyOriginal;
    }
    //-------- Accept or reject according to Metropolis criterion
    if (rand<flo>() < exp(-energyChange/T)) { // could use Glauber criterion if preferred
      status = true;
      spinArray[v0] = s0p;
    }
    else {
      status = false;
    }
  }
  
  bool lastProposalStatus () {return status; }
  int lastClusterSize ()     { return 1; }
};

//=========================================================================
//  GibbsSimulator 
//=========================================================================
class GibbsSimulator {
  bool status;
public:
  GibbsSimulator () {}
  
  //======== COMMON CODE    
  bool lastProposalStatus () { return true; }
  int lastClusterSize ()     { return 1; }
  
  //======== CODE SPECIFIC TO ISING SPINS
  template<class Graph, class EdgeProperty, class VertexProperty>
  void step (std::vector<IsingSpin> &spinArray,
             Graph &graph, 
             EdgeProperty &coupling, 
             VertexProperty &field, 
             flo T = flo(1.0)
  ) {
    Vertex v0 = irand(graph.vertices()); // choose spin to update
    //-------- Calculate net field acting on spin
    Field h    = field (v0);     // magnetic field at vertex v0
    Field hEff = h;                     // initialize accumulator
    for (Neighbor n=0; n<graph.neighbors(v0); ++n) {
      Coupling J   = coupling (v0, n);       // coupling from vertex v0 to its nth neighbor
      Vertex v1    = graph.neighbor (v0, n); // nth neighbor of vertex v0
      IsingSpin s1 = spinArray[v1];          // value of spin at vertex v1
      hEff        += J*s1.moment();          // neighbor contributes effective field -J*s1
    }
    //-------- Resample spin from Bernoulli distribution
    spinArray[v0].randomize (hEff/T);
  }
  //======== CODE SPECIFIC TO POTTS SPINS
  template<class Graph, class EdgeProperty, class VertexProperty, int q>
  void step (std::vector<PottsSpin<q>> &spinArray,
             Graph &graph, 
             EdgeProperty &coupling, 
             VertexProperty &field, 
             flo T = flo(1.0)
  ) {
  }
  
  //======== CODE SPECIFIC TO XY SPINS
  template<class Graph, class EdgeProperty, class VertexProperty>
  void step (std::vector<XYSpin> &spinArray,
             Graph &graph, 
             EdgeProperty &coupling, 
             VertexProperty &field, 
             flo T = flo(1.0)
  ) {
    Vertex v0 = irand(graph.vertices()); // choose spin to update
    //-------- Calculate net field acting on spin
    Field hx = field(v0);     // magnetic field COMPONENTS at vertex v0
    Field hy = 0.0f;          // initialize accumulator
    for (Neighbor n=0; n<graph.neighbors(v0); ++n) {
      auto J  = coupling (v0, n);       // coupling from vertex v0 to its nth neighbor
      auto v1 = graph.neighbor (v0, n); // nth neighbor of vertex v0
      auto s1 = spinArray[v1];          // value of spin at vertex v1
      hx += J*s1.x;          // WE ARE PEEKING AT DATA MEMBERS!
      hy += J*s1.y;          // WE ARE PEEKING AT DATA MEMBERS!
    }
    //-------- Resample XY spin (basically sample from von Mises distribution)
    spinArray[v0].randomize (hx/T, hy/T);
  }
  
  //======== CODE SPECIFIC TO HEISENBERG SPINS
  template<class Graph, class EdgeProperty, class VertexProperty>
  void step (std::vector<HeisenSpin> &spinArray,
             Graph &graph, 
             EdgeProperty &coupling, 
             VertexProperty &field, 
             flo T = flo(1.0)
  ) {
    Vertex v0 = irand(graph.vertices()); // choose spin to update
    //-------- Calculate net field acting on spin
    Field hx = field(v0);     // magnetic field COMPONENTS at vertex v0
    Field hy = 0.0f;          // initialize accumulator
    Field hz = 0.0f;          // initialize accumulator
    for (Neighbor n=0; n<graph.neighbors(v0); ++n) {
      auto J  = coupling (v0, n);       // coupling from vertex v0 to its nth neighbor
      auto v1 = graph.neighbor (v0, n); // nth neighbor of vertex v0
      auto s1 = spinArray[v1];          // value of spin at vertex v1
      hx += J*s1.x;          // WE ARE PEEKING AT DATA MEMBERS!
      hy += J*s1.y;          // WE ARE PEEKING AT DATA MEMBERS!
      hz += J*s1.z;          // WE ARE PEEKING AT DATA MEMBERS!
    }
    //-------- Resample (basically sample from truncated exponential distribution)
    spinArray[v0].randomize (hx/T, hy/T, hz/T);
  }
};






//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//
//	class Cluster
//	Stores information about whether each Vertex is within a Cluster.
//		c.setSize(n)		set up the Cluster-flag array with n slots
//		c.clear()				remove all Vertexs from Cluster
//		c.add(a)				mark Vertex a as being part of Cluster
//		c.contains(a)		returns true or false depending on whether a is in Cluster
//
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class Cluster {
  int tCurrent;
  std::vector<int> tTouched;
public:
  void setSize (int vmax)   { if (tTouched.size() != vmax) {tTouched.resize(vmax); tCurrent=1;} }
  bool contains (Vertex a)	{	return tTouched[a]==tCurrent;			}
  void add (Vertex a)				{	tTouched[a] = tCurrent;						}
  void clear () {
    tCurrent -= 1;			// change current Monte Carlo time
    if (tCurrent==0) {	// once in a blue moon:
      tCurrent = INT_MAX-1;
      for (Vertex a=0; a<tTouched.size(); a++) tTouched[a] = -1;	// explicitly clear Cluster flags
    }
  }
};

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//
//	class Queue
//		Queue(vmax)			set up the Queue with <vmax> slots
//		q.clear()				remove all Vertexs from Queue
//		q.push(a)				add <a> to tail of Queue
//		q.pop()					remove Vertex from head of Queue and return its index
//
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class Queue {
  int lhead,ltail;
  std::vector<int> iQueue;
public:
  void setSize (int vmax) { if (iQueue.size() != vmax) iQueue.resize(vmax); }
  void clear ()           {	lhead = ltail = -1;			}
  void push (Vertex a)    {	iQueue[++ltail] = a;		}
  Vertex pop ()           {	return iQueue[++lhead]; }
  bool isEmpty ()         {	return (lhead==ltail);	}
  int throughput ()       { return (lhead+1);				}
  void restore ()         { lhead = -1; }// restore previously popped elements to queue head
};

//=========================================================================
//  WolffSimulator 
//=========================================================================
class WolffSimulator {
  Cluster cluster;
  Queue queue;
public:
  WolffSimulator () {
    std::cerr << "WARNING: Wolff is only applicable in zero field!\n";
  }
  
  template<class Spin, class Graph, class EdgeProperty, class VertexProperty>
  void step (std::vector<Spin> &spinArray,
             Graph &graph, 
             EdgeProperty &coupling, 
             VertexProperty &field, 
             flo T = flo(1.0)
  ) {
    cluster.setSize (graph.vertices());  // prepare data structure
    queue.setSize (graph.vertices());    // prepare data structure
    
    Spin axis; axis.randomize();	       // choose flip axis (self-inverse symmetry operation)
    Vertex vSeed = irand(graph.vertices()); // choose seed Vertex
    cluster.clear (); cluster.add (vSeed);
    queue.clear (); queue.push (vSeed);
    
    while (!queue.isEmpty()) {
      Vertex v1 = queue.pop ();
      Spin s1  = spinArray[v1];  	// save original spin state at Vertex v1
      Spin s1p = s1.flip (axis);	// flip spin
      spinArray[v1] = s1p;      	// commit the update to the array, since Wolff is rejection-free
      
      for (Neighbor n=0; n<graph.neighbors(v1); ++n) {
        Vertex v2 = graph.neighbor (v1, n);  // loop over all neighbors of v1
        if (!cluster.contains(v2)) {  // Vertex v1 is in the Cluster; Vertex v2 is outside
          Coupling J         = coupling (v1, n);   // coupling from vertex v0 to its nth neighbor
          Spin s2            = spinArray[v2];          // original state of spin v2
          //Spin s2p           = s2.flip (axis);       // propose to add this spin to cluster
          //flo energyOriginal = -J * s1.interact(s2p);
          //flo energyProposed = -J * s1p.interact(s2p);
          flo energyOriginal = -J * s1p.interact(s2);  // optimized version of code
          flo energyProposed = -J * s1.interact(s2);   //
          flo energyChange = energyProposed - energyOriginal;
          flo pstop = exp (+energyChange/T);
          if (rand<flo>() >= pstop) {	// with probability 1 - min(1,pstop)
            cluster.add (v2);
            queue.push (v2);
          }
        }
      }
    }
  }
  
  bool lastProposalStatus () {return true; }
  int lastClusterSize () { return queue.throughput(); }
};

//=========================================================================
//  DSTSimulator: very similar except for field stuff
//=========================================================================
class DSTSimulator {
  Cluster cluster;
  Queue queue;
  bool lastStatus;
  
public:
  int iInnerLoop1; //  number of iterations that have occurred
  int iInnerLoop2;

public:
  template<class Spin, class Graph, class EdgeProperty, class VertexProperty>
  void step (std::vector<Spin> &spinArray,
             Graph &graph, 
             EdgeProperty &coupling, 
             VertexProperty &field, 
             flo T = flo(1.0)
  ) {
    iInnerLoop1 = iInnerLoop1 = 0;
    
    cluster.setSize (graph.vertices());  // prepare data structure
    queue.setSize   (graph.vertices());    // prepare data structure
    
    Spin axis; axis.randomize();	       // choose flip axis (self-inverse symmetry operation)
    Vertex vSeed = irand(graph.vertices()); // choose seed Vertex
    cluster.clear (); cluster.add (vSeed);
    queue.clear (); queue.push (vSeed);
    
    flo energyZeeman = 0;
    
    //======== Build cluster based on couplings J only
    while (!queue.isEmpty()) {
      Vertex v1 = queue.pop ();
      Spin s1  = spinArray[v1];  	// save original spin state at Vertex v1
      Spin s1p = s1.flip (axis);	// flip spin
      spinArray[v1] = s1p;      	// commit the update to the array, since Wolff is rejection-free
      
      Field h = field (v1);
      flo energyOriginal = -s1.interact(h);
      flo energyProposed = -s1p.interact(h);
      energyZeeman += energyProposed - energyOriginal; 
      
      for (Neighbor n=0; n<graph.neighbors(v1); ++n) {
        Vertex v2 = graph.neighbor (v1, n);  // loop over all neighbors of v1
        
        ++iInnerLoop1;
        
        if (!cluster.contains(v2)) {  // Vertex v1 is in the Cluster; Vertex v2 is outside\
          
          ++iInnerLoop2;
          
          Coupling J         = coupling (v1, n);   // coupling from vertex v0 to its nth neighbor
          Spin s2            = spinArray[v2];          // original state of spin v2
          //Spin s2p           = s2.flip (axis);       // propose to add this spin to cluster
          //flo energyOriginal = -J * s1.interact(s2p);
          //flo energyProposed = -J * s1p.interact(s2p);
          flo energyOriginal = -J * s1p.interact(s2);  // optimized version of code
          flo energyProposed = -J * s1.interact(s2);   //
          flo energyChange = energyProposed - energyOriginal;
          flo pstop = exp (+energyChange/T);
          if (rand<flo>() >= pstop) {	// with probability 1 - min(1,pstop)
            cluster.add (v2);
            queue.push (v2);
          }
        }
      }
    }
    //======== Accept or reject according to Metropolis criterion based on fields h only
    flo qAccept = exp (-energyZeeman/T);
    if (rand<flo>() >= qAccept) {
      lastStatus = false;         // rejected
      queue.restore ();           // 
      while (!queue.isEmpty()) {	// go through all spins in cluster
        Vertex v = queue.pop ();                    // 
        spinArray[v] = spinArray[v].flip (axis);		// flip spin in the spin array
      }
    }
    else {
      lastStatus = true;		// accepted
    }
    
  }
  
  bool lastProposalStatus () {return lastStatus; }
  int lastClusterSize () { return queue.throughput(); }
};


//=========================================================================
//  ASWSimulator
//=========================================================================
class ASWSimulator {
  Cluster cluster;
  Queue queue;
  bool firstTime = true;
  bool lastStatus;
public:
  int iInnerLoop1; //  number of iterations that have occurred
  int iInnerLoop2;
  
public:
  template<class Spin, class Graph, class EdgeProperty, class VertexProperty>
  void step (std::vector<Spin> &spinArray,
             Graph &graph, 
             EdgeProperty &coupling, 
             VertexProperty &field, 
             flo T = flo(1.0)
  ) {
    const Vertex vmax = graph.vertices();       // original number of vertices in graph
    const Vertex vAux = vmax;    // index of the auxiliary vertex
    
    //-------- Make sure that spinArray, cluster, and queue have one extra slot to store the aux spin!
    cluster.setSize (vmax + 1);
    queue.setSize   (vmax + 1);
    if (spinArray.size() != vmax+1) spinArray.resize (vmax + 1);
    //-------- Now we can do this thing
    if (firstTime) {
      std::cout << "WARNING: First execution of ASWSimulator::step.\n";
      std::cout << "         Initializing auxiliary spin to (+1, 0, 0).\n";
      spinArray[vAux].randomize (99999.);
      firstTime = false;
    }
    
    
    Spin axis; axis.randomize();	         // choose flip axis (self-inverse symmetry operation)
    Vertex vSeed = irand(vmax);            // choose seed Vertex from original graph
    cluster.clear (); cluster.add (vSeed);
    queue.clear (); queue.push (vSeed);
    
    //======== Build cluster based on couplings J only
    while (!queue.isEmpty()) {
      Vertex v1 = queue.pop ();
      Spin s1  = spinArray[v1];  	// save original spin state at Vertex v1
      Spin s1p = s1.flip (axis);	// flip spin
      spinArray[v1] = s1p;      	// commit the update to the array, since Wolff is rejection-free
      
      if (v1==vAux) {
        //-------- The auxiliary Vertex is coupled to every Vertex on the original lattice
        for (Vertex v2=0; v2<vmax; ++v2) {
          if (!cluster.contains(v2)) {  // Vertex v1 is in the Cluster; Vertex v2 is outside            
            Coupling J    = field (v2);             // effective coupling from vAux to v2
            Spin s2       = spinArray[v2];          // original state of spin v2
            flo energyCur = -J * s1p.interact(s2);  // auxiliary spin interacting with v2
            flo energyNew = -J * s1.interact(s2);   // energy if bond is activated
            if (rand<flo>() >= exp ((energyNew-energyCur)/T)) {
              cluster.add (v2);
              queue.push (v2);
            }
          }
        }        
      }
      else { /* v1 != vAux */
        //-------- The current Vertex is coupled to the auxiliary Vertex
        {
          Vertex v2 = vAux;                   // 
          if (!cluster.contains(v2)) {  // Vertex v1 is in the Cluster; Vertex v2 is outside
            Coupling J    = field (v1);             // effective coupling from v1 to vAux
            Spin s2       = spinArray[v2];          // original state of auxiliary spin
            flo energyCur = -J * s1p.interact(s2);  // auxiliary spin interacting with v2
            flo energyNew = -J * s1.interact(s2);   // energy if bond is activated
            if (rand<flo>() >= exp ((energyNew-energyCur)/T)) {
              cluster.add (v2);
              queue.push (v2);
            }            
          }
        }
        //-------- The current Vertex is coupled to other vertices on the original lattice
        for (Neighbor n=0; n<graph.neighbors(v1); ++n) {
          Vertex v2 = graph.neighbor (v1, n);  // loop over all neighbors of v1
          if (!cluster.contains(v2)) {  // Vertex v1 is in the Cluster; Vertex v2 is outside
            Coupling J    = coupling (v1, n);   // coupling from vertex v0 to its nth neighbor
            Spin s2       = spinArray[v2];          // original state of spin v2
            flo energyCur = -J * s1p.interact(s2);  // optimized version of code
            flo energyNew = -J * s1.interact(s2);   //
            if (rand<flo>() >= exp ((energyNew-energyCur)/T)) {
              cluster.add (v2);
              queue.push (v2);
            }
          }
        }
      }
      
    }
  }
  
  template<class Spin, class Graph>
  void rectify (std::vector<Spin> &spinArray, Graph &graph) {
    Vertex vmax = graph.vertices();
    Vertex vAux = vmax;
    for (Vertex v=0; v<vmax+1; ++v)	spinArray[v].rectify (spinArray[vAux]); 
    // The above code "rectifies" all spins including the auxiliary spin itself
  }
  
  bool lastProposalStatus () {return true; }
  int lastClusterSize () { return queue.throughput(); }
};






//=========================================================================
//  GCASimulator
//=========================================================================
class GCASimulator {
  Cluster cluster;
  Queue queue;
public:
  GCASimulator () {
    std::cerr << "WARNING: GCA conserves magnetization!\n";
  }
  
  template<class Spin, class Graph, class EdgeProperty, class VertexProperty>
  void step (std::vector<Spin> &spinArray,
             Graph &graph, 
             EdgeProperty &coupling, 
             VertexProperty &field, 
             flo T = flo(1.0)
  ) {
    cluster.setSize (graph.vertices());  // prepare data structure
    queue.setSize (graph.vertices());    // prepare data structure
    
    auto gt = graph.randomTransform ();	          // choose geometric symmetry operation
    Vertex v2 = irand(graph.vertices());			   // v2 is seed Vertex
    Vertex v4 = graph.applyTransform (gt, v2);	 // v4 is image of v2
    cluster.clear (); cluster.add (v2); cluster.add (v4);       // mark BOTH as in cluster
    queue.clear (); queue.push (v2);                            // push ONLY v2 on queue
    
    while (!queue.isEmpty()) {
      Vertex v1 = queue.pop ();                   // v1 is in cluster
      Vertex v3 = graph.applyTransform (gt, v1);	// v3 is image of v1
      Spin s1 = spinArray[v1];                           // s1 is original value of ss[v1]
      Spin s3 = spinArray[v3];
      spinArray[v1] = s3; spinArray[v3] = s1;	// exchange spins
      for (int n=0; n<graph.neighbors(v1); n++) {
        Coupling J = coupling (v1, n);             // get coupling
        Vertex v2 = graph.neighbor (v1, n);        // v2 runs over neighbors of v1
        Vertex v4 = graph.applyTransform (gt, v2); // v4 is image of v2
        if (!cluster.contains(v2) && !cluster.contains(v4)) { // v2 may equal v4
          Spin s2 = spinArray[v2];								  // s2 is the original spin at v2
          Spin s4 = spinArray[v4];
          
          
          flo e12 = -J * s1.interact(s2);
          flo e34 = -J * s3.interact(s4);
          flo e14 = -J * s1.interact(s4);
          flo e23 = -J * s2.interact(s3);
          if (rand<flo>() >= exp ((e12+e34-e14-e23)/T)) {
            cluster.add (v2); cluster.add (v4);
            queue.push (v2);
          }
          
        }
      }
      
    }
  }
  
  
  bool lastProposalStatus () {return true; }
  int lastClusterSize () { return queue.throughput(); }
};

#endif
