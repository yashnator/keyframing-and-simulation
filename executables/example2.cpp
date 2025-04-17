#include "camera.hpp"

#define GLM_ENABLE_EXPERIMENTAL

#include "../animation/animation.hpp"  // Include your own header

#include <iostream>
using namespace COL781;
namespace GL = COL781::OpenGL;
using namespace glm;

GL::Rasterizer r;
GL::ShaderProgram program;

const int nv = 4, nt = 2, ne = 4;
std::vector<Vertex> vertices;  // Your Vertex struct
vec3 normals[nv];
ivec3 triangles[nt];
ivec2 edges[ne];

GL::Object object;
GL::AttribBuf vertexBuf, normalBuf;

CameraControl camCtl;

// Just one bone for simplicity
Bone rootBone("root", glm::mat4(1.0f), glm::mat4(1.0f), nullptr);

void initializeScene() {
	object = r.createObject();

	// Create vertices
	std::vector<int> boneIDs = { 0 };
	std::vector<float> weights = { 1.0f };

	vertices.emplace_back(vec3(0, 0, 1), vec3(0, 0, 1), boneIDs, weights);
	vertices.emplace_back(vec3(1, 0, 1), vec3(0, 0, 1), boneIDs, weights);
	vertices.emplace_back(vec3(1, 0, 0), vec3(0, 0, 1), boneIDs, weights);
	vertices.emplace_back(vec3(0, 0, 0), vec3(0, 0, 1), boneIDs, weights);

	vec3 pos[nv];
	for (int i = 0; i < nv; i++) pos[i] = vertices[i].position;
	vertexBuf = r.createVertexAttribs(object, 0, nv, pos);
	for (int i = 0; i < nv; i++) normals[i] = vertices[i].normal;
	normalBuf = r.createVertexAttribs(object, 1, nv, normals);

	triangles[0] = ivec3(0, 1, 2);
	triangles[1] = ivec3(0, 2, 3);
	r.createTriangleIndices(object, nt, triangles);

	edges[0] = ivec2(0, 1);
	edges[1] = ivec2(1, 2);
	edges[2] = ivec2(2, 3);
	edges[3] = ivec2(3, 0);
	r.createEdgeIndices(object, ne, edges);
}

void updateScene(float t) {
    float delta = 60.0f;
	float theta = glm::radians(20.0f * sin(t) + delta);  // Simple oscillating angle
	glm::quat rot = glm::angleAxis(theta, glm::vec3(1, 0, 0));
	rootBone.updateBone(vec3(0, 0, 0), rot);

    vec3 verticesData[nv], normalsData[nv];
	// Update vertices and normals
	for (int i = 0; i < nv; i++) {
		const Vertex &v = vertices[i];
		vec3 worldPos = vec3(rootBone.vertTransform * vec4(v.position, 1.0));
		verticesData[i] = worldPos;

		// Normal matrix = inverse transpose of upper-left 3x3 of transform
		mat3 normalMatrix = transpose(inverse(mat3(rootBone.vertTransform)));
		vec3 worldNormal = normalize(normalMatrix * v.normal);
		normalsData[i] = worldNormal;
        // normalsData[i].y = -normalsData[i].y;
        // std::cout << glm::to_string(normalsData[i]) << std::endl;
	}

	r.updateVertexAttribs(vertexBuf, nv, verticesData);
	r.updateVertexAttribs(normalBuf, nv, normalsData);
}


int main() {
	int width = 640, height = 480;
	if (!r.initialize("Bone Animation", width, height)) {
		return EXIT_FAILURE;
	}
	camCtl.initialize(width, height);
	camCtl.camera.setCameraView(vec3(0.5, -0.5, 1.5), vec3(0.5, -0.5, 0.0), vec3(0.0, 1.0, 0.0));

	program = r.createShaderProgram(
		r.vsBlinnPhong(),
		r.fsBlinnPhong()
	);

	initializeScene();

	while (!r.shouldQuit()) {
		float t = SDL_GetTicks64() * 1e-3;
		updateScene(t);

		camCtl.update();
		Camera &camera = camCtl.camera;

		r.clear(vec4(0.4, 0.4, 0.4, 1.0));
		r.enableDepthTest();
		r.useShaderProgram(program);

		r.setUniform(program, "model", glm::mat4(1.0));
		r.setUniform(program, "view", camera.getViewMatrix());
		r.setUniform(program, "projection", camera.getProjectionMatrix());
		r.setUniform(program, "lightPos", camera.position);
		r.setUniform(program, "viewPos", camera.position);
		r.setUniform(program, "lightColor", vec3(1.0f));

		r.setupFilledFaces();
		r.setUniform(program, "ambientColor", 0.2f * vec3(1));
		r.setUniform(program, "extdiffuseColor", 0.9f * vec3(1.0f, 0.6f, 0.2f));
		r.setUniform(program, "intdiffuseColor", 0.4f * vec3(1.0f, 0.6f, 0.2f));
		r.setUniform(program, "specularColor", 0.6f * vec3(1));
		r.setUniform(program, "phongExponent", 20.f);
		r.drawTriangles(object);

		r.setupWireFrame();
		r.setUniform(program, "ambientColor", vec3(0));
		r.setUniform(program, "extdiffuseColor", vec3(0));
		r.setUniform(program, "intdiffuseColor", vec3(0));
		r.setUniform(program, "specularColor", vec3(0));
		r.setUniform(program, "phongExponent", 0.f);
		r.drawEdges(object);

		r.show();
	}
}
