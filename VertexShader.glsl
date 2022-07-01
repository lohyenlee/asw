#version 330 core

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in vec3 vertexNormal_modelspace;
layout(location = 2) in vec3 vertexColor;

// Output data; will be interpolated for each fragment.
out vec3 vertexPosition_worldspace;
out vec3 normalDirection_cameraspace;
out vec3 eyeDirection_cameraspace;
out vec3 lightDirection_cameraspace;
out vec3 vertexColor0;

// Values that stay constant for the whole mesh.
uniform mat4 MVP;
uniform mat4 V;
uniform mat4 M;
uniform vec3 lightPosition_worldspace;

void main () {

	// Output position of the vertex, in clip space : MVP * position
	gl_Position = MVP * vec4(vertexPosition_modelspace,1);
	
	// Position of the vertex, in worldspace : M * position
	vertexPosition_worldspace = (M * vec4(vertexPosition_modelspace,1)).xyz;
	
	// Vector that goes from the vertex to the camera, in camera space.
	// In camera space, the camera is at the origin (0,0,0).
	vec3 vertexPosition_cameraspace = (V * M * vec4(vertexPosition_modelspace,1)).xyz;
	eyeDirection_cameraspace = vec3(0,0,0) - vertexPosition_cameraspace;

	// Vector that goes from the vertex to the light, in camera space. 
	vec3 lightPosition_cameraspace = ( V * vec4(lightPosition_worldspace,1)).xyz;
	lightDirection_cameraspace = lightPosition_cameraspace + eyeDirection_cameraspace;
	
	// Normal of the the vertex, in camera space
	normalDirection_cameraspace = (V * M * vec4(vertexNormal_modelspace,0)).xyz; 
	// Only correct if ModelMatrix does not scale the model ! Use its inverse transpose if not.
  
  vertexColor0 = vertexColor;
}

