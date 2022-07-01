#ifndef SPIN_HH
#define SPIN_HH

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Please do typedef flo float; or typedef flo double; before #include "spin.hh"
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// This C++ header file provides the following types and implements the following methods:
//
//   Spin interface:
//     Spin()           constructor
//     asciiArt()       string with ANSI codes for visualization in terminal
//     direction()      3-vector for graphical visualization using OpenGL
//     moment()         scalar representing moment along primary direction
//     interact(h)      energy of interaction with a scalar field h
//     interact(s)      form of interaction energy with another spin of the same type
//     flip(a)          return a new Spin object, reflected along axis a
//     rectify(aux)     find a transformation R that takes aux to the primary direction, and apply R in-place
//     randomize()      generate a completely random spin
//     randomize(bhx)   generate spin from equilibrium distribution at temperature T and field (hx,0,0)
//     randomize(bhx, bhy, bhz)     
//                      generate spin from equilibrium distribution at temperature T and field (hx,0,0)
//   Implementations:
//     IsingSpin
//     PottsSpin<q>
//     ClockSpin<q>
//     XYSpin
//     HeisenSpin
//     XYPhase          not implemented
//
// The operation s.rectify(aux) is a peculiar operation used by the auxiliary-spin algorithm.
// It rotates the current spin (s) in a manner such that the auxiliary spin (sAux)
//   would be mapped to the special axis of the field (the direction in which field is applied).
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

// Struct for returning 3-vector representation of spin for visualization purposes
namespace spin {struct vec3 {float x,y,z;};}; 

static char const* const ansiCodes[] = {
  "\u001b[0;41m0", "\u001b[0;43m1", "\u001b[0;42m2", "\u001b[0;46m3", "\u001b[0;44m4", 
  "\u001b[0;45m5",  "\u001b[0;30;47m6", "\u001b[0;40m7"
};

//=========================================================================
// IsingSpin
//=========================================================================
class IsingSpin {
  char s;
public:  
  const char* asciiArt ()         { return s>0 ? "\u001b[0;41m^" : "\u001b[0;44mv"; }
  spin::vec3 direction ()         { return spin::vec3 {0, 0, flo(s)};}
  IsingSpin (char s=+1)						{ this->s = s;}
  int moment ()                   { return s;           } // return moment in direction of field
  flo interact  (Field h)         { return s*h;         } // "h dot S" which is simply h*S
  int interact  (IsingSpin other) { return s*other.s; } // "S dot S" which is simply S*S
  IsingSpin flip (IsingSpin axis) {	return IsingSpin( -s );						}
  void rectify (IsingSpin aux)    { s *= aux.s; }
  void randomize () 							{ s = irand(2) ? 1 : -1;						}      
  void randomize (Field betah)    { s = (rand<flo>() < 1 / (1+exp(-2*betah))) ? +1 : -1; }
};

//=========================================================================
// PottsSpin
//=========================================================================
template<int q>
class PottsSpin {
  char s;
public:
  PottsSpin (char s=0)		  			{ this->s = s;}
  const char* asciiArt ()         { return ansiCodes[s%8]; }
  spin::vec3 direction ()         { return spin::vec3 {(flo) cos(2*M_PI/q*s), (flo)sin(2*M_PI/q*s), 0};}
  int moment ()                   { return (s==0);           } // return moment in direction of field
  flo interact  (Field h)         { return (s==0)*h;     } // s.interact(h)=h*KroneckerDelta(s,0)
  int interact  (PottsSpin other) { return (s==other.s); } // si.interact(sj)=1 iff si==sj
  PottsSpin flip (PottsSpin axis)	{	return PottsSpin((axis.s - s + q) % q);	}
  void randomize () 							{ s = irand(q);       		}
  void randomize (Field betah)    { std::cerr << "NOT IMPLEMENTED!\n"; abort(); }
};

//=========================================================================
// ClockSpin
//=========================================================================
template<int q>
class ClockSpin {
  char s;
public:
  ClockSpin (char s=0)		  			{ this->s = s;}
  const char* asciiArt ()         { return ansiCodes[s%8]; }  
  spin::vec3 direction ()         { return spin::vec3 {(flo) cos(2*M_PI/q*s), (flo)sin(2*M_PI/q*s), 0};}
  flo moment ()                   { return cos(2*M_PI/q*s);           } // return moment in direction of field
  flo interact  (Field h)         { return h*cos(2*M_PI/q*s);     }
  flo interact  (ClockSpin other) { return cos(2*M_PI/q*(s - other.s)); } 
  ClockSpin flip (ClockSpin axis)	{ return ClockSpin((axis.s - s + q) % q);	}
  void randomize () 						  { s = irand(q);       		}
};


//=========================================================================
// XYSpin
//=========================================================================
class XYSpin {
public:
  flo x,y;  // ONLY ACCESS THESE DATA MEMBERS IF YOU REALLY, REALLY NEED TO!
public:
  XYSpin (flo x=1, flo y=0)	: x(x), y(y) {}
  const char* asciiArt ()       { return ansiCodes[(int)floor(   (1+atan2(y,x)/(2*M_PI)) * 6) % 6]; } 
  spin::vec3 direction ()       { return spin::vec3 {x, y, 0};}
  flo moment ()                 { return x;           } // return moment in direction of field
  flo interact  (Field h)       { return h*x;     }
  flo interact  (XYSpin other)  { return x*other.x + y*other.y; } 
  //======== Reflect along unit 2-vector u: sNew = s - 2(u.s)u
  XYSpin flip (XYSpin u)	{
    XYSpin sNew (x, y);
    double us = u.x * x + u.y * y;
    sNew.x -= 2*us*u.x;	sNew.y -= 2*us*u.y;
    return sNew;
  }
  void rectify (XYSpin aux) { // assume that aux and (1,0) are unit 2-vectors!
    flo ax = aux.x - 1;              // we want to find a reflection that maps aux to (1,0)
    flo ay = aux.y - 0;
    flo aa = ax*ax + ay*ay;
    if (fabs(aa) < 1e-7) return;               // s and S are identical; do nothing!
    flo k = 2*(ax * x + ay * y)/aa;  // k = 2(a dot r)/a^2
    x -= k * ax;                     // Rr = r - 2(a dot r)/a^2 a
    y -= k * ay;
  }
  //======== Internal helper function
  static flo sampleVonMises (flo k) {
    flo phi;
    if (k<3.4083f) {
      while (true) {
        flo u = std::uniform_real_distribution<flo>() (rng);
        flo v = std::uniform_real_distribution<flo>() (rng);
        flo r = std::uniform_real_distribution<flo>() (rng);
        if (v < 1 + tanh(k)*cos(2*M_PI*u))    phi = 2*M_PI*u;
        else                                  phi = 2*M_PI*(.5 - u); // defer modulo operation
        if (r < exp(k*cos(phi)) / (cosh(k) + sinh(k)*cos(phi))) break;
      }
    } else {
      while (true) {
        flo u = std::uniform_real_distribution<flo>() (rng) - 0.5f;
        flo r = std::uniform_real_distribution<flo>() (rng);
        phi = sqrt(2/k) * atanh (2*u*tanh(sqrt(k/2)*M_PI));
        flo ch = cosh (sqrt(k/2)*phi);
        if (r < ch*ch*exp(k*(cos(phi)-1)) ) break;
      }
    }
    if (phi>M_PI) phi -= 2*M_PI;
    return phi;
  }
  //======== Generate random 2-vector uniformly from the unit circle
  void randomize () {
    flo rsq; // magnitude squared
    do {
      x = std::uniform_real_distribution<flo>() (rng) - 0.5f;
      y = std::uniform_real_distribution<flo>() (rng) - 0.5f;
      rsq = x*x + y*y;
    } while (rsq > 0.25 || rsq < 1e-8);
    flo oor = 1/sqrt(rsq);
    x *= oor;
    y *= oor;
  }
  //======== Generate random XY spin in field h at temperature T
  void randomize (Field bhx, Field bhy=0) {  // beta h_x and beta h_y
    flo bh = hypot (bhx, bhy);     // field strength k = beta h
    flo phi0 = atan2 (bhy, bhx);
    flo phi = sampleVonMises (bh);
    phi += phi0;
    x = cos(phi);
    y = sin(phi);    
  }  
};

//=========================================================================
// HeisenSpin
//=========================================================================
class HeisenSpin {
public:
  flo x,y,z;
public:
  HeisenSpin (flo x=1, flo y=0, flo z=0)	: x(x), y(y), z(z) {}
  const char* asciiArt ()           { return "H";                 }   // difficult to represent 3-vector with ASCII
  spin::vec3 direction ()           { return spin::vec3 {x, y, z};}  
  flo moment ()                     { return x;                   } // return moment in direction of field
  flo interact  (Field h)           { return h*x;                 }
  flo interact  (HeisenSpin other)  { return x*other.x + y*other.y + z*other.z; } 
  //======== Reflect along unit 3-vector u
  HeisenSpin flip (HeisenSpin u)	{
    HeisenSpin sNew (x, y, z);
    double us = u.x * x + u.y * y + u.z * z;
    sNew.x -= 2*us*u.x;	sNew.y -= 2*us*u.y; sNew.z -= 2*us*u.z;
    return sNew;
  }
  //======== Find R such that R dot u=v, and transform *this to R dot (*this)
  void rotateToAlign (HeisenSpin u, HeisenSpin v) {
    // Assume that u and v are unit 3-vectors
    flo ax = u.x - v.x;
    flo ay = u.y - v.y;
    flo az = u.z - v.z;
    flo aa = ax*ax + ay*ay + az*az;
    if (fabs(aa) < 1e-7) return;              // if u and v are identical, do nothing
    //HeisenSpin& r = *this; /// was this the culprit?
    flo k = 2*(ax*x + ay*y + az*z)/aa;  // k = 2(a dot r)/a^2
    x -= k * ax;                            // Rr = r - 2(a dot r)/a^2 a
    y -= k * ay;
    z -= k * az;
  }
  //======== Find R such that R dot u=(1,0,0), and transform *this to R dot (*this)
  void rectify (HeisenSpin aux) {
    rotateToAlign (aux, HeisenSpin(1,0,0));
  }
  //======== Internal helper function
  static flo sampleTruncExp (flo k) {
    flo u = std::uniform_real_distribution<flo>() (rng);
    if (k<8.)   return 1/k*log(1 + 2*u/(1/tanh(k) - 1)) - 1;
    else        return 1+1/k*log(u);
  }
  //======== Generate random 3-vector uniformly from the unit sphere
  void randomize () {
    x = std::normal_distribution<flo>() (rng);
    y = std::normal_distribution<flo>() (rng);
    z = std::normal_distribution<flo>() (rng);
    flo rsq = x*x + y*y + z*z;
    flo oor = 1/sqrt(rsq);
    x *= oor;
    y *= oor;
    z *= oor;
  }
  //======== Generate random Heisenberg spin in field h at temperature T
  void randomize (Field bhx, Field bhy=0, Field bhz=0) {  // beta h_x, beta h_y, beta h_z
    flo bh = sqrt (bhx*bhx + bhy*bhy + bhz*bhz);  // field strength k = beta h
    // Sample spin at temperature T in field h in +z direction
    flo costheta = sampleTruncExp (bh);      
    flo sintheta = sqrt(1-costheta*costheta);
    flo phi = std::uniform_real_distribution<flo>(0,2*M_PI) (rng);
    x = sintheta*cos(phi); // set (x,y,z) components of spin
    y = sintheta*sin(phi);
    z = costheta;
    // Rotate current spin so that current z axis becomes field direction
    flo nx = bhx/bh; // (nx,ny,nz) is unit vector in direction of field
    flo ny = bhy/bh; 
    flo nz = bhz/bh;
    rotateToAlign (HeisenSpin (0,0,1), HeisenSpin (nx,ny,nz)); 
  }
  
};
#endif
