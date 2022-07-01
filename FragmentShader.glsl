#version 330 core

//================ STANDARD FRAGMENT SHADER
//================ Interpolated values from the vertex shaders
in vec3 vertexPosition_worldspace;
in vec3 normalDirection_cameraspace;
in vec3 eyeDirection_cameraspace;
in vec3 lightDirection_cameraspace;
in vec3 vertexColor0;
//================ Output values
out vec3 color; // color of the fragment or pixel
//================ Uniform values (constant over the whole mesh)
uniform mat4 MV;
uniform vec3 lightPosition_worldspace;
uniform vec3 lightColor;

void main () {
	// Material properties
  vec3 materialDiffuseColor = vertexColor0;
 	vec3 materialAmbientColor = .3 * materialDiffuseColor;
	vec3 materialSpecularColor = vec3 (1., 1., 1.);
  float specularExponent = 50;

	// Distance to the light
	float distance = length( lightPosition_worldspace - vertexPosition_worldspace );

	// Normal of the computed fragment, in camera space
	vec3 n = normalize( normalDirection_cameraspace );
	// Direction of the light (from the fragment to the light)
	vec3 l = normalize( lightDirection_cameraspace );
	// Cosine of the angle between the normal and the light direction,
	// clamped above 0
	//  - light is at the vertical of the triangle -> 1
	//  - light is perpendicular to the triangle -> 0
	//  - light is behind the triangle -> 0
	float cosTheta = clamp( dot(n,l), 0,1 );

	// Eye vector (towards the camera)
	vec3 E = normalize(eyeDirection_cameraspace);
	// Direction in which the triangle reflects the light
	vec3 R = reflect(-l,n);
	// Cosine of the angle between the Eye vector and the Reflect vector,
	// clamped to 0
	//  - Looking into the reflection -> 1
	//  - Looking elsewhere -> < 1
	float cosAlpha = clamp( dot(E,R), 0,1 );

	color =
		materialAmbientColor +
		materialDiffuseColor * lightColor * cosTheta / (distance*distance) + 
		materialSpecularColor * lightColor * pow(cosAlpha,specularExponent) / (distance*distance);
}

