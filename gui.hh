#ifndef GUI_HH
#define GUI_HH

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "utils.hh"

using glm::vec3; using glm::mat3; using glm::mat4; using glm::normalize; // too painful without this!

//=======================================================================
//  UTILITY FUNCTIONS ON TOP OF GLM (glm::vec3, etc.)
//=======================================================================
std::ostream& operator<< (std::ostream& out, vec3 a) {
  using namespace std;
	for (int i=0; i<3; ++i) {out << setw(15) << a[i]; } out << endl; return out;
}
std::ostream& operator<< (std::ostream& out, mat3 a) {
  using namespace std;
	for (int i=0; i<3; ++i) { for (int j=0; j<3; ++j) { out << setw(15) << a[j][i]; } // Note order!!!!
	out << endl; } out << endl; return out;
}
mat3 rotateToAlign (vec3 cc, vec3 bb) { // Find R such that R.z==c and R.y==b (if possible)
	vec3 c = normalize(cc); // "Up" vector (the new z axis)
	vec3 b = normalize(bb); // "North" vector (a vector that is kinda along the new y axis)
	vec3 a = cross(b, c);   // "East" vector (the new x axis)
	if (dot(a,a) < 1e-8) {b = vec3 (1, 0, 0); a = cross(b, c);} // Attempt 2
	if (dot(a,a) < 1e-8) {b = vec3 (0, 1, 0); a = cross(b, c);} // Attempt 3 -- must succeed!
	a = normalize (a);
	b = cross(c, a);
	return mat3 (a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z);
} // Note that (a.x, a.y, a.z) is the first COLUMN of this matrix!

//=======================================================================
//  3D GEOMETRIC PRIMITIVES
//  This code replicates the functionality of the legacy OpenGL state machine.
//  It maintains the following state variables:
//  - std::vectors of positions, normals, and triangles, to be communicated to the GPU
//  - rotation matrix and translation vector for model (matR and vecT)
//=======================================================================
typedef GLuint idx;  // either GLuint or GLushort
const auto MY_ELEMENT_INDEX_TYPE = (sizeof(idx)==2) ? GL_UNSIGNED_SHORT : GL_UNSIGNED_INT;
typedef float flo;

class Model {
public:
	std::vector<vec3> positions;
	std::vector<vec3> normals;
	std::vector<vec3> colors;
	std::vector<idx> triangles;  // warning: can only index 65536 vertices!
	mat3 matR; // rotation matrix (could also include scale?)
	vec3 vecT; // translation vector
	vec3 colorCurrent; // current color
	int indexOffset = 0;  // starting index in the mesh at which we begin adding new triangles
	GLuint bufPositions,bufNormals,bufColors,bufTriangles;
	//-------- Default constructor
	Model () : matR(1),vecT(0),colorCurrent(.6, .7, .8) {
		bufPositions = bufNormals = bufColors = bufTriangles = 0;
	}
	~Model () {
		if (bufPositions==0) return;
		glDeleteBuffers (1, &bufPositions);
		glDeleteBuffers (1, &bufNormals);
		glDeleteBuffers (1, &bufColors);
		glDeleteBuffers (1, &bufTriangles);
	}
	void clear () {
		positions.clear();
		normals.clear();
		colors.clear();
		triangles.clear();
		indexOffset = 0;
	}
	void rot ()                                  { matR = mat3(1); }
	void rot (mat3 rotmat)                       { matR = rotmat; }
	void trans ()                                { vecT = vec3(0); }
	void trans (GLfloat x, GLfloat y, GLfloat z) { vecT = vec3(x,y,z); }
	void color (vec3 rgb)                        {colorCurrent = rgb;}
	void color (GLfloat r, GLfloat g, GLfloat b) {colorCurrent = vec3(r,g,b);}
	
	//-------- set matR and vecT such that (0,0,0)-->(x1,y1,z1) and (0,0,1)-->(z2,y2,z2)
	void align (GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2) {
		rot   (rotateToAlign (vec3(x2-x1, y2-y1, z2-z1), vec3(0,1,0)) ); // may need to hack
		trans (x1, y1, z1);
	}
	void pos (GLfloat x, GLfloat y, GLfloat z) {positions.push_back (matR*vec3(x,y,z) + vecT); }
	void pos (vec3 r)                          {positions.push_back (matR*r           + vecT); }
	void nrm (GLfloat x, GLfloat y, GLfloat z) {normals.push_back   (matR*vec3(x,y,z)); }
	void nrm (vec3 n)                          {normals.push_back   (matR*n); }
	void clr (vec3 color)                      {colors.push_back    (color); }
	void clr ()                                {colors.push_back    (colorCurrent); }
	void index (idx i)                         {triangles.push_back (indexOffset+i); }
	void index (idx i1, idx i2, idx i3)        {index(i1); index(i2); index(i3); }
	
	
	
	//======== TRIANGLE
	// Build one face of a triangle given 3 vertex position vectors (auto-compute normal)
	void triangle (GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat x3, GLfloat y3, GLfloat z3) {
		indexOffset = positions.size(); // Start counting from first empty slot
		vec3 a(x1,y1,z1), b(x2,y2,z2), c(x3,y3,z3);
		vec3 n = cross(b-a, c-a);
		pos (a); nrm (n); clr ();
		pos (b); nrm (n); clr ();
		pos (c); nrm (n); clr ();
		index (0,1,2); // for reverse face use (0,2,1)
	}
	
	//======== QUAD
	// Build a planar quadrilateral given 4 vertex position vectors and 1 normal vector
	void quad (GLfloat x1, GLfloat y1, GLfloat z1,
						 GLfloat x2, GLfloat y2, GLfloat z2,
						GLfloat x3, GLfloat y3, GLfloat z3,
						GLfloat x4, GLfloat y4, GLfloat z4,
						GLfloat xn, GLfloat yn, GLfloat zn
	) {
		indexOffset = positions.size(); // Start counting from first empty slot
		pos (x1, y1, z1); nrm (xn, yn, zn); clr ();
		pos (x2, y2, z2); nrm (xn, yn, zn); clr ();
		pos (x3, y3, z3); nrm (xn, yn, zn); clr ();
		pos (x4, y4, z4); nrm (xn, yn, zn); clr ();
		index (0,1,3); index (3,1,2);
	}
	//======== CUBOID
	void cuboid (GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2) {
		quad (x1,y1,z1, x1,y2,z1, x2,y2,z1, x2,y1,z1, 0,0,-1);	// bottom face
		quad (x1,y1,z2, x2,y1,z2, x2,y2,z2, x1,y2,z2, 0,0,+1);  // top face
		quad (x1,y1,z1, x1,y1,z2, x1,y2,z2, x1,y2,z1, -1,0,0);  // left 
		quad (x2,y1,z1, x2,y2,z1, x2,y2,z2, x2,y1,z2, +1,0,0);  // right
		quad (x1,y1,z1, x2,y1,z1, x2,y1,z2, x1,y1,z2, 0,-1,0);  // back
		quad (x1,y2,z1, x1,y2,z2, x2,y2,z2, x2,y2,z1, 0,+1,0);  // front
	}
	
	//======== DISK
	void disk (GLfloat r, int nmax=24, GLfloat h=0, bool drawTop=true, bool drawBottom=true) {
		using namespace glm;
		indexOffset = positions.size(); // Start counting from first empty slot
		float phiStep = 2.0*M_PI/nmax;
		if (drawTop) {
			for (int n=0; n<nmax; ++n) {
				float phi = n*phiStep;
				pos (r*cos(phi), r*sin(phi), h); nrm (0,0,+1); clr (); // ring of points
			}
			pos (0,0,0); nrm (0,0,+1); clr (); // center point of disk
			for (int n=0; n<nmax; ++n) index (nmax, n, (n+1)%nmax);  // triangle
		}
		if (drawBottom) {
			for (int n=0; n<nmax; ++n) {
				float phi = n*phiStep;
				pos (r*cos(phi), r*sin(phi), h); nrm (0,0,-1); clr (); // ring of points
			}
			pos (0,0,0); nrm (0,0,-1); clr (); // center point of disk
			for (int n=0; n<nmax; ++n) index (nmax, (n+1)%nmax, n);  // triangle
		}
	}
	//======== TUBE
	// Curved surface of a conical frustum with base radius a, top radius b, and height h
	void tube (GLfloat a, GLfloat b, GLfloat h, const int nmax=24) {
		indexOffset = positions.size(); // Start counting from first empty slot
		for (int n=0; n<nmax; ++n) {
			float phi = n*(2.0*M_PI/nmax), c=cos(phi), s=sin(phi);
			pos (a*c, a*s, 0); nrm (c, s, 0); clr (); // point on bottom rim
			pos (b*c, b*s, h); nrm (c, s, 0); clr (); // point on top rim
		}
		idx imax = 2*nmax;
		for (idx i=0; i<imax; i+=2) { // increase by stride 2
			index ((i+0)%imax, (i+2)%imax, (i+1)%imax);
			index ((i+1)%imax, (i+2)%imax, (i+3)%imax);
		}
	};
	void tube (GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, 
						 GLfloat a, GLfloat b, const int nmax=24
	) {
		align (x1, y1, z1, x2, y2, z2);
		flo h = sqrt (sqr(x2-x1) + sqr(y2-y1) + sqr(z2-z1));
		tube (a, b, h, nmax);
	};
	//======== FRUSTUM, CYLINDER, CONE: BASED ON TUBE AND DISK
	void frustum (GLfloat a, GLfloat b, GLfloat h, const int nmax=24) {
		tube (a, b, h, nmax);  					 	// First build curved surface of cylinder
		disk (b, nmax, h, true, false);		// Top face of cylinder
		disk (a, nmax, 0, false, true);		// Bottom face of cylinder
	};
	void frustum (GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat a, GLfloat b, const int nmax=24) {
		align (x1, y1, z1, x2, y2, z2);
		flo h = sqrt (sqr(x2-x1) + sqr(y2-y1) + sqr(z2-z1));
		frustum (a, b, h, nmax);
	};
	void cylinder (GLfloat r, GLfloat h, const int nmax=24) {	
		frustum (r, r, h, nmax); 
	}
	void cylinder (GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat r, const int nmax=24) {
		frustum (x1, y1, z1, x2, y2, z2, r, r, nmax);
	};
	void cone (GLfloat r, GLfloat h, const int nmax=24) {	
		frustum (r, 0, h, nmax); 
	}
	void cone (GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat r, const int nmax=24) {
		frustum (x1, y1, z1, x2, y2, z2, r, 0, nmax);
	};
	//======== TORUS 
	// mmax x nmax TESSELLATION IN PHI AND THETA DIRECTIONS
	void torus (GLfloat a, GLfloat b, const int mmax=12, const int nmax=12) {
		indexOffset = positions.size(); // Start counting from first empty slot
		for (int n=0; n<nmax; ++n) {
			for (int m=0; m<mmax; ++m) {
				float theta = n*(2.0*M_PI/nmax), C = cos(theta), S = sin(theta);
				float phi   = m*(2.0*M_PI/mmax), c = cos(phi),   s = sin(phi);
				pos (a*c+b*c*C, a*s+b*s*C, b*S); // point on torus
				nrm (      c*C,       s*C,   S); // outward normal
				clr ();
			}
		}
		for (int m=0; m<mmax; ++m) {
			for (int n=0; n<nmax; ++n) {
				int m2 = (m+1)%mmax;
				int n2 = (n+1)%nmax;
				index (m  + mmax*n, m2 + mmax*n, m  + mmax*n2);
				index (m  + mmax*n2, m2 + mmax*n, m2 + mmax*n2);
			}
		}
	};
	void torus (GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, 
							GLfloat a, GLfloat b, const int mmax=12, const int nmax=12
	) {
		align (x1, y1, z1, x2, y2, z2);
		torus (a, b, mmax, nmax);
	};
	//======== ELLIPSOID
	// mmax x nmax TESSELLATION IN THETA AND PHI DIRECTIONS
	void ellipsoid (GLfloat a, GLfloat b, GLfloat c, const int mmax=12, const int nmax=12) {
		indexOffset = positions.size(); // Start counting from first empty slot
		for (int m=0; m<mmax; ++m) {
			for (int n=0; n<nmax; ++n) {
				float theta = m*(2.0*M_PI/mmax), ct = cos(theta), st = sin(theta);
				float phi   = n*(2.0*M_PI/nmax), cp = cos(phi),   sp = sin(phi);
				pos (a*st*cp  , b*st*sp  , c*ct  ); // point on sphere
				nrm (  st*cp/a,   st*sp/b,   ct/c); // outward normal
				clr ();
			}
		}
		for (int m=0; m<mmax; ++m) {
			for (int n=0; n<nmax; ++n) {
				int m2 = (m+1)%mmax;
				int n2 = (n+1)%nmax;
				index (m*nmax + n, m2*nmax + n, m*nmax + n2);
				index (m*nmax + n2, m2*nmax + n, m2*nmax + n2);
			}
		}
	};
	//======== SPHERE: SPECIAL CASE OF ELLIPSOID
	void sphere (GLfloat r, const int mmax=12, const int nmax=12) {
		ellipsoid (r, r, r, mmax, nmax);
	};
	
	//=========================================================================
	// model.upload (gui);
	// Allocates buffers on the GPU
	// and uploads model data (positions, normals, colors, triangles) to this buffer.
	//=========================================================================	
	void uploadToGPU () {
		//======== If buffers haven't been allocated yet, allocate them now.
		//======== To prevent against segfaults, we should further check that each GPU buffer size
		//======== matches each CPU buffer size, otherwise we also have to reallocate.
		if (bufPositions==0 || bufNormals==0 || bufColors==0 || bufTriangles==0) {
			glGenBuffers (1, &bufPositions);
			glBindBuffer (GL_ARRAY_BUFFER, bufPositions);
			glBufferData (GL_ARRAY_BUFFER, positions.size() * sizeof(vec3), NULL, GL_DYNAMIC_DRAW);
			glGenBuffers (1, &bufNormals);
			glBindBuffer (GL_ARRAY_BUFFER, bufNormals);
			glBufferData (GL_ARRAY_BUFFER, normals.size() * sizeof(vec3), NULL, GL_DYNAMIC_DRAW);
			glGenBuffers (1, &bufColors);
			glBindBuffer (GL_ARRAY_BUFFER, bufColors);
			glBufferData (GL_ARRAY_BUFFER, colors.size() * sizeof(vec3), NULL, GL_DYNAMIC_DRAW);
			glGenBuffers (1, &bufTriangles);
			glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, bufTriangles);
			glBufferData (GL_ELEMENT_ARRAY_BUFFER, triangles.size() * sizeof(idx), NULL, GL_DYNAMIC_DRAW);
		}
		glBindBuffer (GL_ARRAY_BUFFER, bufPositions);
		glBufferSubData (GL_ARRAY_BUFFER, 0, positions.size() * sizeof(vec3), &positions[0]);
		glBindBuffer (GL_ARRAY_BUFFER, bufNormals);
		glBufferSubData (GL_ARRAY_BUFFER, 0, normals.size() * sizeof(vec3), &normals[0]);
		glBindBuffer (GL_ARRAY_BUFFER, bufColors);
		glBufferSubData (GL_ARRAY_BUFFER, 0, colors.size() * sizeof(vec3), &colors[0]);
		glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, bufTriangles); // these last 2 commands really only need to be run once
		glBufferSubData (GL_ELEMENT_ARRAY_BUFFER, 0, triangles.size() * sizeof(idx), &triangles[0]);
	}
};

class GUI {
public:
	//======== GL-related variables
	GLFWwindow* window=NULL;
	//int windowWidth=1536, windowHeight=768;
	int windowWidth=600, windowHeight=600;
	GLuint idProgram,idMVP,idV,idM,idVertexArray,idLightPosition,idLightColor;	
	//======== Lighting-related variables
	vec3 lightColor {1., 1., 1.};
	vec3 lightPos {25, 25, 20}; // Light position in WORLD space ...ABOVE spin array CENTER
	//======== Camera-related variables
	const float degree = M_PI/180;
 	float camFoV = 45.0*degree; // Initial Field of View
	float speedMouse = 0.001f;
	bool shouldExit = false;
	
	//-------- Default constructor
	GUI () {}
	
	
	void setup (int windowWidth, int windowHeight) {
    using namespace std;
    this->windowWidth = windowWidth;
    this->windowHeight = windowHeight;
  
		if (!glfwInit()) {cerr << "Failed to initialize GLFW\n"; abort();}
		glfwWindowHint (GLFW_SAMPLES, 4);
		glfwWindowHint (GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint (GLFW_CONTEXT_VERSION_MINOR, 3);
		glfwWindowHint (GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
		glfwWindowHint (GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
		window = glfwCreateWindow (windowWidth, windowHeight, "Hello world!", NULL, NULL);
		if (window==NULL) {cerr  << "Failed to open GLFW window.\n"; glfwTerminate(); abort();}
		glfwMakeContextCurrent(window);
		glewExperimental = true; // Needed for core profile
		if (glewInit() != GLEW_OK) {cerr << "Failed to initialize GLEW\n"; glfwTerminate(); abort(); }
		glfwSwapInterval(0); // remove the 60 FPS framerate cap on Macs
		glfwSetInputMode (window, GLFW_STICKY_KEYS, GL_TRUE);
		glfwSetInputMode (window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
		glfwSetWindowPos (window, 0, 0);
		glClearColor (0, 0, 0, 0); // If you don't like this, you can override this later in main()
		glEnable(GL_DEPTH_TEST);
		glDepthMask(GL_TRUE);
		glDepthFunc(GL_LEQUAL);
		glDepthRange(0.0f, 1.0f);
		glEnable (GL_CULL_FACE);
		
		const char* filenameVS = "VertexShader.glsl";
		const char* filenameFS = "FragmentShader.glsl";
 		//const char* codeVS = readString (filenameVS).c_str(); // Doesn't work because 
    //const char* codeFS = readString (filenameFS).c_str(); //temporary std::string gets deallocated!
    const char* codeVS = readCString (filenameVS);
		const char* codeFS = readCString (filenameFS);
		int logLength;
		string strLog;
		cout << "Compiling " << filenameVS << "\n";
		GLuint idVS = glCreateShader (GL_VERTEX_SHADER);
		glShaderSource (idVS, 1, &codeVS, NULL);
		glCompileShader (idVS);
		glGetShaderiv (idVS, GL_INFO_LOG_LENGTH, &logLength);     strLog.resize (logLength+1);
		glGetShaderInfoLog (idVS, logLength, NULL, &strLog[0]);   cout << strLog;
		cout << "Compiling " << filenameFS << "\n";
		GLuint idFS = glCreateShader (GL_FRAGMENT_SHADER);
		glShaderSource (idFS, 1, &codeFS, NULL);
		glCompileShader (idFS);
		glGetShaderiv (idFS, GL_INFO_LOG_LENGTH, &logLength);     strLog.resize (logLength+1);
		glGetShaderInfoLog (idFS, logLength, NULL, &strLog[0]);   cout << strLog;
		cout << "Linking program...\n";
		idProgram = glCreateProgram ();
		glAttachShader (idProgram, idVS);
		glAttachShader (idProgram, idFS);
		glLinkProgram (idProgram);
		glGetProgramiv (idProgram, GL_INFO_LOG_LENGTH, &logLength); strLog.resize (logLength+1);
		glGetShaderInfoLog (idFS, logLength, NULL, &strLog[0]);   cout << strLog;
		cout << "Detaching and deleting shaders...\n";
		glDetachShader (idProgram, idVS);
		glDetachShader (idProgram, idFS);
		glDeleteShader (idVS);
		glDeleteShader (idFS);
		glUseProgram (idProgram); // use our program consisting of the two shaders
		
		idMVP = glGetUniformLocation (idProgram, "MVP");
		idV = glGetUniformLocation (idProgram, "V");
		idM = glGetUniformLocation (idProgram, "M");		
		idLightPosition = glGetUniformLocation(idProgram, "lightPosition_worldspace");
		idLightColor = glGetUniformLocation (idProgram, "lightColor");
		
		glGenVertexArrays (1, &idVertexArray);
		glBindVertexArray (idVertexArray);  // Even though this is not explicitly accessed, it is necessary!
	}
	
	//@@@@@@@@ STATE VARIABLES FOR AIRPLANE CONTROLS (keeping z-axis vertical)
	vec3 camPos {30, 30, 10}; // Initially a bit off the field (assume xmax*xmax)
 	float camTheta = 100*degree;          // Initially looking straight ahead ... a bit down
 	float camPhi = 180*degree;           // Initially looking in the ?
  float speedStrafe = 10.0f; // units / second
	float speedTurn = 1.f;
	
  int pressed (int key) {return glfwGetKey(window,key)==GLFW_PRESS;}

	//@@@@@@@@ NEED TO MAKE SpaceshipControls where there is no sense of vertical
// 	void processAirplaneControls () {
// 		//--------  Compute the matMVP matrix from keyboard and mouse input, and tell the shader
// 		static double lastTime = glfwGetTime();  // One-time call
// 		double currentTime = glfwGetTime();
// 		float deltaTime = float(currentTime - lastTime);
// 		lastTime = currentTime;
// 		
// 		//============ YLL: 2021-11-24 : DISABLE MOUSE PANNING
// // 		double xpos, ypos;
// // 		glfwGetCursorPos(window, &xpos, &ypos); // Get mouse position
// // 		glfwSetCursorPos(window, windowWidth/2, windowHeight/2); // Reset mouse position for next frame
// // 		camPhi += speedMouse * float(windowWidth/2 - xpos );  // Turn camera
// // 		camTheta -= speedMouse * float(windowHeight/2 - ypos ); camTheta = glm::clamp (camTheta, 0.f, 3.141f);
// 		
// 		vec3 forward (sin(camTheta)*cos(camPhi), sin(camTheta)*sin(camPhi), cos(camTheta));
// 		vec3 right (sin(camTheta)*sin(camPhi), -sin(camTheta)*cos(camPhi), 0);
// 		vec3 up = glm::cross (right, forward);
//  		    
// 		if (glfwGetKey (window, GLFW_KEY_W) == GLFW_PRESS) {camPos += forward * deltaTime * speedStrafe;}
// 		if (glfwGetKey (window, GLFW_KEY_S) == GLFW_PRESS) {camPos -= forward * deltaTime * speedStrafe;}
// 		if (glfwGetKey (window, GLFW_KEY_D) == GLFW_PRESS) {camPos += right * deltaTime * speedStrafe; } 
// 		if (glfwGetKey (window, GLFW_KEY_A) == GLFW_PRESS) {camPos -= right * deltaTime * speedStrafe; }
// 		if (glfwGetKey (window, GLFW_KEY_Q) == GLFW_PRESS) {camPos += up * deltaTime * speedStrafe; } 
// 		if (glfwGetKey (window, GLFW_KEY_Z) == GLFW_PRESS) {camPos -= up * deltaTime * speedStrafe; }
// 		if (glfwGetKey (window, GLFW_KEY_UP) == GLFW_PRESS)   {camTheta -= deltaTime * speedTurn;}
// 		if (glfwGetKey (window, GLFW_KEY_DOWN) == GLFW_PRESS) {camTheta += deltaTime * speedTurn;}
// 		if (glfwGetKey (window, GLFW_KEY_LEFT) == GLFW_PRESS)  {camPhi += deltaTime * speedTurn;}
// 		if (glfwGetKey (window, GLFW_KEY_RIGHT) == GLFW_PRESS) {camPhi -= deltaTime * speedTurn;}
// 
//     if (glfwGetKey (window, GLFW_KEY_ESCAPE)==GLFW_PRESS || glfwWindowShouldClose(window)) shouldExit = true;
// 		
// 		// Projection matrix: 45 degree field of view, 4:3 aspect ratio, display range 0.1--1000 units
// 		// View matrix:       Use camera position and forward direction
// 		// Model matrix:      Fix this as the identity matrix
// 		mat4 matP = glm::perspective (camFoV, float(windowWidth)/windowHeight, 0.1f, 1000.0f);
// 		mat4 matV = glm::lookAt (camPos, camPos+forward, up);
// 		mat4 matM = mat4 (1.0);
// 		mat4 matMVP = matP * matV * matM;
// 		glUniformMatrix4fv (idMVP, 1, GL_FALSE, &matMVP[0][0]);
// 		glUniformMatrix4fv (idM, 1, GL_FALSE, &matM[0][0]);
// 		glUniformMatrix4fv (idV, 1, GL_FALSE, &matV[0][0]);
// 	}
	
  //@@@@@@@@ STATE VARIABLES FOR TRACKBALL CONTROLS
	vec3 rcom;
  GLfloat rsys;
  GLfloat zoomOutDistance;
  mat4 camOrient;
  float speedZoom = 0.5f; // fraction for logarithmic zooming
  
  //@@@@@@@@ Update model, view, and projection matrices, and upload to shader
  void processTrackballControls () {
    static double tPrev = glfwGetTime();  // One-time call
    double tNow = glfwGetTime();
    float tElap = float(tNow - tPrev);
    tPrev = tNow;
    //-------- Handle ESCAPE key
    if (pressed(GLFW_KEY_ESCAPE) || glfwWindowShouldClose(window)) shouldExit = true;
    if (pressed(GLFW_KEY_LEFT_SHIFT) || pressed(GLFW_KEY_RIGHT_SHIFT))
    { //-------- If SHIFT is pressed, use arrow keys to zoom
      if      (pressed(GLFW_KEY_UP)   || pressed(GLFW_KEY_KP_8)) {zoomOutDistance /= (1 + tElap * speedZoom);}
      else if (pressed(GLFW_KEY_DOWN) || pressed(GLFW_KEY_KP_2)) {zoomOutDistance *= (1 + tElap * speedZoom);}
      else if (pressed(GLFW_KEY_LEFT) || pressed(GLFW_KEY_KP_4)) {camOrient = glm::rotate(mat4(1.0f), tElap*speedTurn, vec3(0,0,+1)) * camOrient;}
      else if (pressed(GLFW_KEY_RIGHT)|| pressed(GLFW_KEY_KP_6)) {camOrient = glm::rotate(mat4(1.0f), tElap*speedTurn, vec3(0,0,-1)) * camOrient;}
    }
    else 
    { //-------- Use arrow keys to rotate
      if      (pressed(GLFW_KEY_UP)   || pressed(GLFW_KEY_KP_8)) {camOrient = glm::rotate(mat4(1.0f), tElap*speedTurn, vec3(-1,0,0)) * camOrient;}
      else if (pressed(GLFW_KEY_DOWN) || pressed(GLFW_KEY_KP_2)) {camOrient = glm::rotate(mat4(1.0f), tElap*speedTurn, vec3(+1,0,0)) * camOrient;}
      else if (pressed(GLFW_KEY_LEFT) || pressed(GLFW_KEY_KP_4)) {camOrient = glm::rotate(mat4(1.0f), tElap*speedTurn, vec3(0,-1,0)) * camOrient;}
      else if (pressed(GLFW_KEY_RIGHT)|| pressed(GLFW_KEY_KP_6)) {camOrient = glm::rotate(mat4(1.0f), tElap*speedTurn, vec3(0,+1,0)) * camOrient;}
    }
    
    // Projection matrix: 45 degree field of view, 4:3 aspect ratio, display range 0.1--1000 units
    // View matrix:       Move camera to rCoM, rotate, zoom out to zoomOutDistance
    // Model matrix:      Fix this as the identity matrix
    mat4 matP = glm::perspective (camFoV, float(windowWidth)/windowHeight, 0.1f, 1000.0f);
    mat4 matV = 
    glm::translate (mat4(1.0), vec3(0,0,-zoomOutDistance)) 
    * camOrient
    * glm::translate (mat4(1.0), -rcom) 
    ;  /// very ugly, should simplify
    mat4 matM = mat4 (1.0);
    mat4 matMVP = matP * matV * matM;
    glUniformMatrix4fv (idMVP, 1, GL_FALSE, &matMVP[0][0]);
    glUniformMatrix4fv (idM, 1, GL_FALSE, &matM[0][0]);
    glUniformMatrix4fv (idV, 1, GL_FALSE, &matV[0][0]);
  }
  
  void render (Model &model) {				
    glUniform3fv (idLightColor, 1, &lightColor[0]);
		glUniform3f (idLightPosition, lightPos.x, lightPos.y, lightPos.z);
		
		glBindBuffer (GL_ARRAY_BUFFER, model.bufPositions);
		glVertexAttribPointer (0, 3, GL_FLOAT, GL_FALSE, 0, NULL); // attrib, size, type, normalized, stride, offset
		glEnableVertexAttribArray (0);

		glBindBuffer (GL_ARRAY_BUFFER, model.bufNormals);
		glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
		glEnableVertexAttribArray (1);
		
		glBindBuffer (GL_ARRAY_BUFFER, model.bufColors);
		glVertexAttribPointer (2, 3, GL_FLOAT, GL_FALSE, 0, NULL);
		glEnableVertexAttribArray (2);
		
		glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, model.bufTriangles);
    
    glDrawElements (GL_TRIANGLES, model.triangles.size(), MY_ELEMENT_INDEX_TYPE, NULL);
	}
	void clear ()       {  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); }
	void swapBuffers () {  glfwSwapBuffers (window);                            }
	~GUI () {
    using namespace std;
		if (window==NULL) return;
		cerr << "Destructor ~GUI() is being executed!\n";
		glDeleteProgram (idProgram);
		glDeleteVertexArrays (1, &idVertexArray);
		glfwTerminate();
	}
};
#endif
