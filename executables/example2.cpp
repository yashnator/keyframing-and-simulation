#include "camera.hpp"

#define GLM_ENABLE_EXPERIMENTAL

#include "../animation/animation.hpp"  // Include your own header

#include <iostream>
using namespace COL781;
namespace GL = COL781::OpenGL;
using namespace glm;

GL::Rasterizer r;
GL::ShaderProgram program;

int nv, nt, ne;

std::vector<Bone> boneArray;
std::vector<Vertex> vertices;  // Your Vertex struct


GL::Object object;
GL::AttribBuf vertexBuf, normalBuf;

CameraControl camCtl;

// #define DBG(x) std::cout << to_string(x) << std::endl;

std::map<std::string, int> boneIdx;

void setVertData(vec3 &vertData, vec3& normalData, Bone &bone, Vertex &v) {
    vertData = vec3(bone.vertTransform * vec4(v.position, 1.0f));
    normalData = vec3(bone.vertTransformIT * vec4(v.normal, 0.0f));
}

void initializeBones() {
    nv = 0;
    nt = 0;
    ne = 0;
    for(auto &bone: boneArray) {
        nv += bone.renderVertices.size();
        nt += bone.renderTriangles.size();
        ne += bone.renderEdges.size();
    }
    vec3 pos[nv], normals[nv];
    ivec3 triangles[nt];
    ivec2 edges[ne];
    int vptr = 0, nptr = 0, tptr = 0, eptr = 0, bptr = 0;
    int offset = 0;
    for(auto &bone: boneArray) {
        ivec3 toff(offset);
        ivec2 eoff(offset);
       for(auto& p: bone.renderTriangles) {
           triangles[tptr] = p + toff;
           ++tptr;
       }
       for(auto &p: bone.renderVertices) {
           pos[vptr] = p;
           ++vptr;
       }
       for(auto &p: bone.renderNormals) {
           normals[nptr] = p;
           ++nptr;
           vertices.emplace_back(Vertex(pos[nptr - 1], normals[nptr - 1], {bptr}));
       }
       for(auto &p: bone.renderEdges) {
           edges[eptr] = p + eoff;
           ++eptr;
       }
    //    std::cout << bptr << std::endl;
       boneIdx[bone.boneName] = bptr;
       ++bptr;
       offset += bone.renderVertices.size();
    }
    boneArray.front().updateAll();
    for (int i = 0; i < nv; i++) {
        setVertData(pos[i], normals[i], boneArray[vertices[i].boneIDs[0]], vertices[i]);
	}
    vertexBuf = r.createVertexAttribs(object, 0, nv, pos);
	normalBuf = r.createVertexAttribs(object, 1, nv, normals);
	r.createTriangleIndices(object, nt, triangles);
	r.createEdgeIndices(object, ne, edges);
}

void initializeScene() {
    boneArray.reserve(100);
    // Mesh headCube = unitSphere(1, 1);
    Mesh headSphere = unitSphere(25, 25);
    Mesh torsoCube = unitCube(3, 3, 3);

    // Torso is the root of all bones
    boneArray.emplace_back(Bone("torso"));
    boneArray.back().getMeshAttribs(torsoCube);
    boneArray.back().updateInit(scale(mat4(1.0f), vec3(2.0f, 4.0f, 1.0f)));

   // Head
    boneArray.emplace_back(Bone("head", &boneArray[0]));
    boneArray.back().getMeshAttribs(headSphere);
    boneArray.back().updateInit(translate(mat4(1.0f), vec3(0.0f, 1.0f, 0.0f)));
    boneArray.back().updateBindPose(translate(mat4(1.0f), vec3(0.0f, 2.0f, 0.0f)));
    boneArray[0].children.push_back(&boneArray[1]);

    // Arms
    boneArray.emplace_back(Bone("left arm", &boneArray[0]));
    boneArray.back().getMeshAttribs(torsoCube);
    boneArray.back().updateInit(scale(mat4(1.0f), vec3(1.0f, 2.0f, 1.0f)));
    boneArray.back().updateInit(translate(mat4(1.0f), vec3(0.5f, -1.0f, 0.0f)));
    boneArray.back().updateBindPose(translate(mat4(1.0f), vec3(1.0f, 2.0f, 0.0f)));
    boneArray[0].children.push_back(&boneArray[2]);

    boneArray.emplace_back(Bone("right arm", &boneArray[0]));
    boneArray.back().getMeshAttribs(torsoCube);
    boneArray.back().updateInit(scale(mat4(1.0f), vec3(1.0f, 2.0f, 1.0f)));
    boneArray.back().updateInit(translate(mat4(1.0f), vec3(-0.5f, -1.0f, 0.0f)));
    boneArray.back().updateBindPose(translate(mat4(1.0f), vec3(-1.0f, 2.0f, 0.0f)));
    boneArray[0].children.push_back(&boneArray[3]);

    // Elbows
    boneArray.emplace_back(Bone("left elbow", &boneArray[2]));
    boneArray.back().getMeshAttribs(torsoCube);
    boneArray.back().updateInit(scale(mat4(1.0f), vec3(0.8f, 1.5f, 0.8f)));
    boneArray.back().updateInit(translate(mat4(1.0f), vec3(0.0f, -0.75f, 0.0f)));
    boneArray.back().updateBindPose(translate(mat4(1.0f), vec3(0.5f, -2.0f, 0.0f)));
    boneArray[2].children.push_back(&boneArray[4]);

    boneArray.emplace_back(Bone("right elbow", &boneArray[3]));
    boneArray.back().getMeshAttribs(torsoCube);
    boneArray.back().updateInit(scale(mat4(1.0f), vec3(0.8f, 1.5f, 0.8f)));
    boneArray.back().updateInit(translate(mat4(1.0f), vec3(0.0f, -0.75f, 0.0f)));
    boneArray.back().updateBindPose(translate(mat4(1.0f), vec3(-0.5f, -2.0f, 0.0f)));
    boneArray[3].children.push_back(&boneArray[5]);

    //Legs
    boneArray.emplace_back(Bone("left leg", &boneArray[0]));
    boneArray.back().getMeshAttribs(torsoCube);
    boneArray.back().updateInit(scale(mat4(1.0f), vec3(0.8f, 3.0f, 0.8f)));
    boneArray.back().updateInit(translate(mat4(1.0f), vec3(0.5f, -1.0f, 0.0f)));
    boneArray.back().updateBindPose(translate(mat4(1.0f), vec3(0.0f,-2.4f, 0.0f)));
    boneArray[0].children.push_back(&boneArray[6]);

    boneArray.emplace_back(Bone("right leg", &boneArray[0]));
    boneArray.back().getMeshAttribs(torsoCube);
    boneArray.back().updateInit(scale(mat4(1.0f), vec3(0.8f, 3.0f, 0.8f)));
    boneArray.back().updateInit(translate(mat4(1.0f), vec3(-0.5f, -1.0f, 0.0f)));
    boneArray.back().updateBindPose(translate(mat4(1.0f), vec3(0.0f,-2.4f, 0.0f)));
    boneArray[0].children.push_back(&boneArray[7]);

    // Updating initial positions

    // Moving head over torso

	object = r.createObject();
    initializeBones();

}


void updateScene(float t) {
    float delta = 0.0f;
	float theta = glm::radians(30.0f * sin(t) + delta);  // Simple oscillating angle
	glm::quat rot = glm::angleAxis(theta, glm::vec3(1, 0, 0));
	boneArray[0].updateBone(glm::translate(mat4(1.0f), vec3(0.0f, 0.0f, 5.0*sin(t))));
	boneArray[1].updateBone(vec3(0, 0, 0), rot);
	boneArray[2].updateBone(vec3(0, 0, 0), glm::angleAxis(theta, glm::vec3(1, 0, 0)));
	boneArray[3].updateBone(vec3(0, 0, 0), glm::angleAxis(theta, glm::vec3(-1, 0, 0)));
    boneArray[4].updateBone(vec3(0, 0, 0), glm::angleAxis(theta, glm::vec3(1, 0, 0)));
    boneArray[5].updateBone(vec3(0, 0, 0), glm::angleAxis(theta, glm::vec3(-1, 0, 0)));
    boneArray[6].updateBone(vec3(0, 0, 0), glm::angleAxis(theta, glm::vec3(-1, 0, 0)));
    boneArray[7].updateBone(vec3(0, 0, 0), glm::angleAxis(theta, glm::vec3(1, 0, 0)));
    boneArray[0].updateAll();

    vec3 verticesData[nv], normalsData[nv];
	for (int i = 0; i < nv; i++) {
        setVertData(verticesData[i], normalsData[i], boneArray[vertices[i].boneIDs[0]], vertices[i]);
	}
	r.updateVertexAttribs(vertexBuf, nv, verticesData);
	r.updateVertexAttribs(normalBuf, nv, normalsData);
}

void updateSceneCatmull(float t) {
    float t1 = 0.0f;
    float t2 = 1.0f;
    float t3 = 1.5f;
    float t4 = 2.0f;
    float t5 = 4.5f;

    float p0 = 0.0f;
    float p1 = 5.0f;
    float p2 = -5.0f;
    float p3 = 10.0f;
    float p4 = -10.0f;

	float theta = catmullRom5NonUniform(t1, t2, t3, t4, t5, p0, p1, p2, p3, p4, t);
    // std::cout<<t<<" "<<theta<<std::endl;

    // std::cout<<theta<<std::endl;

    float q0 = 0.0f;
    float q1 = 0.5f;
    float q2 = 1.0f;
    float q3 = 1.5f;
    float q4 = 2.0f;

    theta = glm::radians(theta);  // Simple oscillating angle
    float move = catmullRom5NonUniform(t1, t2,t3,t4,t5, q0, q1, q2, q3, q4, t);
	glm::quat rot = glm::angleAxis(theta, glm::vec3(1, 0, 0));
	boneArray[0].updateBone(glm::translate(mat4(1.0f), vec3(0.0f, 0.0f, move)));
	boneArray[1].updateBone(vec3(0, 0, 0), rot);
	boneArray[2].updateBone(vec3(0, 0, 0), glm::angleAxis(theta, glm::vec3(1, 0, 0)));
	boneArray[3].updateBone(vec3(0, 0, 0), glm::angleAxis(theta, glm::vec3(-1, 0, 0)));
    boneArray[4].updateBone(vec3(0, 0, 0), glm::angleAxis(theta, glm::vec3(1, 0, 0)));
    boneArray[5].updateBone(vec3(0, 0, 0), glm::angleAxis(theta, glm::vec3(-1, 0, 0)));
    boneArray[6].updateBone(vec3(0, 0, 0), glm::angleAxis(theta, glm::vec3(-1, 0, 0)));
    boneArray[7].updateBone(vec3(0, 0, 0), glm::angleAxis(theta, glm::vec3(1, 0, 0)));
    boneArray[0].updateAll();

    vec3 verticesData[nv], normalsData[nv];
	for (int i = 0; i < nv; i++) {
        setVertData(verticesData[i], normalsData[i], boneArray[vertices[i].boneIDs[0]], vertices[i]);
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
	camCtl.camera.setCameraView(vec3(0.0f, 0.0f, 12.0f), vec3(0.0f, 0.0f, 0.0f), vec3(0.0, 1.0, 0.0));

	program = r.createShaderProgram(
		r.vsBlinnPhong(),
		r.fsBlinnPhong()
	);

	initializeScene();

	while (!r.shouldQuit()) {
        float t = SDL_GetTicks64()*1e-3;

        if(t<1) continue;
		// updateScene(t);
        updateSceneCatmull(t);
        if(t>5) break;

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
