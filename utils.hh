#ifndef UTILS_HH
#define UTILS_HH

#include <cmath>      // for M_PI
#include <cassert>

#include <algorithm> 
//#include <cctype>
//#include <locale>

#include <string>			// for string

#include <iostream>		// for cout
#include <iomanip>		// for setf
#include <sstream>		// for ostringstream
#include <fstream>		// for ofstream and ifstream
#include <vector>			// for std::vector
#include <unistd.h> 	// for usleep
#include <climits>		// for INT_MAX
#include <chrono>			// for chrono
#include <random> // for random numbers
#include <fstream> // for ifstream
#include <GL/glew.h>
#include <GLFW/glfw3.h>

// using std::cout; using std::cerr; using std::endl; using std::ostream; using std::setw;
// using std::string; using std::ostringstream; using std::ofstream; using std::ifstream; using std::istringstream;

const bool DEBUGGING = true;

//=======================================================================
//  GENERAL MATH
//=======================================================================
template<typename T> T sqr (T x)              {  return x*x;}
int quotient (int p, int q)                   {	return p/q - (p%q<0);	}
int modulo (int p, int q)                     {  return (p >= 0)  ?  (p%q)  :  (q-1-((~p)%q));	}
float modulo (float p, float q)               {  return p - q*floor(p/q);}
double modulo (double p, double q)            {  return p - q*floor(p/q);}		// Mod[p,q]
template<typename T> T modulo (T p, T q, T r) {  return modulo (p-r, q) + r;}	// Mod[p,q,r]
template<typename T> T modulo (T p)           {  return p - floor(p);	}				// FractionalPart[p]

//===========================================================================================
//  INCREMENT/DECREMENT ACCORDING TO SEQUENCE SUCH AS 1, 1.01, 1.02, ..., 9.99, 10, 10.1, ...
//===========================================================================================
template<typename flo>
void changeParameter (flo &x, double absStep, double relStep) {
  //-------- Prepare human-friendly increment to a real-valued parameter
  int exponent1 = floor (log10 (fabs (x)));
  int exponent2 = floor (log10 (fabs (relStep)));
  int exponent3 = exponent1 + exponent2;
  int mantissa1 = round(x       / exp10(exponent3));
  int mantissa2 = round(relStep / exp10(exponent2));
  double change = mantissa2 * exp10 (exponent3);
  if (x==0 || fabs(change) < fabs(absStep)) {
    //-------- Override with absolute step (and sanitize)
    int exponent2 = floor (log10 (fabs (absStep)));
    int mantissa1 = round(x       / exp10(exponent2));
    int mantissa2 = round(absStep / exp10(exponent2));
    x = (mantissa1 + mantissa2) * exp10(exponent2);
  } else {
    x = (mantissa1 + mantissa2) * exp10(exponent3);
  }
}

//=======================================================================
//  RANDOM NUMBERS (WRAPPER AROUND std::mt19937_64)
//=======================================================================
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // 12345
std::mt19937_64 rng (seed);   // std::default_random_engine
int irand(int n)    {return std::uniform_int_distribution<int> (0, n-1) (rng);}

template<typename T> T rand() {return std::uniform_real_distribution<T> () (rng); }
//=======================================================================
//  STOPWATCH (WRAPPER AROUND std::chrono)
//=======================================================================
class Stopwatch {
public:
	std::chrono::time_point<std::chrono::high_resolution_clock> timeStarted;
	double timeElapsed;
  Stopwatch () {timeStarted = std::chrono::high_resolution_clock::now();}
	void tic()   {timeStarted = std::chrono::high_resolution_clock::now();}
	double toc() {return timeElapsed = std::chrono::duration<double> 
		(std::chrono::high_resolution_clock::now() - timeStarted) . count();}
};

//=======================================================================
//  UTILITY FUNCTIONS BASED ON std::vector AND std::string FUNCTIONS
//=======================================================================
template<class T>
std::ostream& operator<< (std::ostream& out, std::vector<T> &vec) {
	for (auto x:vec) out << x << "\n";
	return out;
}
//-------- Return min and max element in a std::vector
template<class T> T min (std::vector<T> &a) {return *std::min_element(a.begin(), a.end());}
template<class T> T max (std::vector<T> &a) {return *std::max_element(a.begin(), a.end());}

//-------- Trim a string in-place
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) { return !std::isspace(ch); }));
}
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
}
static inline void trim(std::string &s) { ltrim(s); rtrim(s); }
//-------- Return a trimmed string
static inline std::string ltrim_copy(std::string s) { ltrim(s); return s; }
static inline std::string rtrim_copy(std::string s) { rtrim(s); return s; }
static inline std::string trim_copy(std::string s)  {  trim(s); return s; }
//-------- Read entire file into string
std::string readString (std::string filename) {
  using namespace std;
  ifstream f (filename);
  if (!f.is_open()) {cerr << "readString: failed to open " << filename << "!\n"; abort();}
  stringstream ss; ss << f.rdbuf();
  return ss.str();
}
//-------- Read entire file into C string, making sure to allocate permanent storage
const char* readCString (const char* fileName) {
  using namespace std;
  ifstream f (fileName);
  if (!f.is_open()) {cerr << "Failed to open " << fileName << "!\n"; exit(1);}
  string* ps = new string ( (std::istreambuf_iterator<char>(f) ),
                            (std::istreambuf_iterator<char>()  ) );
  return ps->c_str();
}

//-------- Read whitespace-delimited ASCII numerical data from file into std::vector
template<class T>
void readArray (std::string filename, std::vector<T> &an, int nmax) {
  using std::cerr;
	std::ifstream ifs (filename);
	if (ifs.fail()) {cerr << "readArray: failed to open " << filename << "\n"; abort();}
	an.resize (nmax);
	for (int n=0; n<nmax; ++n) ifs >> an[n];
	if (ifs.fail()) {cerr << "readArray: failed to read " << filename << "\n"; abort();}
	ifs.close ();
}

//=======================================================================
//  PARAMETER FILE PARSING
//  Suppose pars.dat contains the following text:
// 
//     # This is a comment
//     Temperature=2.27    # This is another comment
//     Field      =0.01
//     MCSteps    = 10000
//
//  Then the main C++ code can extract parameter values as follows:
//     int steps;
//     double field;
//     readPar ("pars.dat", "MCSteps", steps);
//     readPar ("pars.dat", "Field", field);
//=======================================================================
template<class T>
void readPar (std::string sInFileName, std::string sPar, T& var) {
  using namespace std;
	ifstream f (sInFileName);
	if (f.fail()) {cerr << "File " << sInFileName << " not found!\n"; abort();}
  while (!f.eof()) {
		string s; getline (f, s);
    ltrim (s);                  // remove whitespace from left
		
		auto p = s.find("#"); 			// find first "#"
		if (p!=string::npos) s = s.substr (0, p-1); // trim any comments starting at the #
		if (s=="") continue;        // empty string
		auto p1 = s.find("=");  		// split at the = sign
		if (p1==string::npos) {cerr << "Invalid parameter line!\n"; abort();}
		string s1,s2;  // we are expecting "PARNAME = PARVALUE"
		s1 = s.substr(0, p1); trim(s1);
		s2 = s.substr(p1+1, string::npos); trim(s2);
		//cerr << "KEY|" << s1 << "|    VALUE|" << s2 << "|" <<  endl;
		if (s1==sPar) {
			istringstream(s2) >> var; 
			cerr << sInFileName << ": " << s1 << " = " << var << "\n";
			f.close ();
			return;		
		}
	}
	cerr << "Using default: " << sPar << " = " << var << "\n";
	f.close ();
}

#endif
