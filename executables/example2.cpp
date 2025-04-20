#include "camera.hpp"

#define GLM_ENABLE_EXPERIMENTAL

#include "../animation/animation.hpp"  // Include your own header

#include <iostream>
using namespace COL781;
namespace GL = COL781::OpenGL;
using namespace glm;

GL::Rasterizer r;
GL::ShaderProgram program;

// int nv = 4, nt = 2, ne = 4;

std::vector<Bone> boneArray;
std::vector<Vertex> vertices;  // Your Vertex struct


GL::Object object;
GL::AttribBuf vertexBuf, normalBuf;

CameraControl camCtl;

// Just one bone for simplicity
// Bone rootBone("root", glm::mat4(1.0f), glm::mat4(1.0f), nullptr);
// Bone headBone("head", glm::mat4(1.0f), glm::mat4(1.0f), nullptr);

#define DBG(x) std::cout << to_string(x) << std::endl;

void initializeBones() {
    int nv = 0, nt = 0, ne = 0;
    for(auto &bone: boneArray) {
        nv += bone.renderVertices.size();
        nt += bone.renderTriangles.size();
        ne += bone.renderEdges.size();
    }
    vec3 pos[nv], normals[nv];
    ivec3 triangles[nt];
    ivec2 edges[ne];
    int vptr = 0, nptr = 0, tptr = 0, eptr = 0;
    for(auto &bone: boneArray) {
       for(auto& p: bone.renderTriangles) {
           triangles[tptr] = p;
           ++tptr;
       }
       for(auto &p: bone.renderVertices) {
           pos[vptr] = p;
           ++vptr;
       }
       for(auto &p: bone.renderNormals) {
           normals[nptr] = p;
           ++nptr;
       }
       for(auto &p: bone.renderEdges) {
           edges[eptr] = p;
           ++eptr;
       }
    }
    vertexBuf = r.createVertexAttribs(object, 0, nv, pos);
	normalBuf = r.createVertexAttribs(object, 1, nv, normals);
	r.createTriangleIndices(object, nt, triangles);
	r.createEdgeIndices(object, ne, edges);
}

void initializeScene() {

    // Mesh headCube = unitSphere(1, 1);
    Mesh headCube = unitCube(1, 1, 1);
    boneArray.emplace_back(Bone("head", glm::mat4(1.0f), glm::mat4(1.0f), nullptr));
    boneArray[0].getMeshAttribs(headCube);
	object = r.createObject();
    initializeBones();
}

void updateScene(float t) {
    // float delta = 0.0f;
	// float theta = glm::radians(20.0f * sin(t) + delta);  // Simple oscillating angle
	// glm::quat rot = glm::angleAxis(theta, glm::vec3(1, 0, 0));
	// rootBone.updateBone(vec3(0, 0, 0), rot);

    // vec3 verticesData[nv], normalsData[nv];
	// // // Update vertices and normals
	// for (int i = 0; i < nv; i++) {
		// const Vertex &v = vertices[i];
		// vec3 worldPos = vec3(rootBone.vertTransform * vec4(v.position, 1.0));
		// verticesData[i] = worldPos;

		// Normal matrix = inverse transpose of upper-left 3x3 of transform
		// mat3 normalMatrix = inverseTranspose(mat3(rootBone.vertTransform));
		// vec3 worldNormal = normalize(normalMatrix * v.normal);
		// normalsData[i] = worldNormal;
	// }

	// r.updateVertexAttribs(vertexBuf, nv, verticesData);
	// r.updateVertexAttribs(normalBuf, nv, normalsData);
}


int main() {
	int width = 640, height = 480;
	if (!r.initialize("Bone Animation", width, height)) {
		return EXIT_FAILURE;
	}
	camCtl.initialize(width, height);
	camCtl.camera.setCameraView(vec3(0.0f, 0.0f, -5.5f), vec3(0.0f, 0.0f, 0.0f), vec3(0.0, 1.0, 0.0));

	program = r.createShaderProgram(
		r.vsBlinnPhong(),
		r.fsBlinnPhong()
	);

	initializeScene();

	while (!r.shouldQuit()) {
        float t = SDL_GetTicks64()*1e-3;
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
		r.setUniform(program, "lightColor", vec3(1.0f, 1.0f, 1.0f));

        // std::cout << glm::to_string(camera.position) << std::endl;

		r.setupFilledFaces();
        glm::vec3 orange(1.0f, 0.6f, 0.2f);
        glm::vec3 white(1.0f, 1.0f, 1.0f);
        r.setUniform(program, "ambientColor", 0.2f*white);
        r.setUniform(program, "extdiffuseColor", 0.9f*orange);
        r.setUniform(program, "intdiffuseColor", 0.4f*orange);
        r.setUniform(program, "specularColor", 0.6f*white);
        r.setUniform(program, "phongExponent", 20.f);
		r.drawTriangles(object);

		r.setupWireFrame();
        glm::vec3 black(0.0f, 0.0f, 0.0f);
        r.setUniform(program, "ambientColor", black);
        r.setUniform(program, "extdiffuseColor", black);
        r.setUniform(program, "intdiffuseColor", black);
        r.setUniform(program, "specularColor", black);
        r.setUniform(program, "phongExponent", 0.f);
		r.drawEdges(object);

		r.show();
    }
}
