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

std::map<std::string, int> boneIdx;

int m = 25, n = 25;
Mesh cloth = unitSquare(n,m);

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
    // Torso is the root of all bones
    boneArray.emplace_back(Bone("torso"));
    boneArray.back().getMeshAttribs(cloth);
    boneArray.back().updateInit(translate(mat4(1.0f), vec3(-0.5f, -0.5f, 0.0f)));
    boneArray.back().updateInit(scale(mat4(1.0f), vec3(3.0f, 3.0f, 3.0f)));
    boneArray.back().updateInit(rotate(mat4(1.0f), radians(-90.0f), vec3(1.0f, 0.0f, 0.0f)));
    boneArray.back().n = n;
    boneArray.back().m = m;
    boneArray.back().checkSelfCollisions = true;

	object = r.createObject();
    initializeBones();
    for(int i = 0; i < nv; ++i) {
        vertices[i].mass = 0.1f;
    }
    float ks = 5.0f * 1e0;
    float kd = 1.0f;
    for(int i=0;i<(m+1)*(n+1);i++)
    {
        int c1 = i/(n+1), c2 = i%(n+1);
        if(c2>0) vertices[i].addStructural(&vertices[i-1], 300.0f, length(vertices[i].position - vertices[i-1].position), kd);
        if(c1<m) vertices[i].addStructural(&vertices[i+(n+1)], 300.0f, length(vertices[i].position - vertices[i+(n+1)].position), kd);
        if(c2>0 && c1<m) vertices[i].addShear(&vertices[i+(n+1)-1], 120.0f, length(vertices[i].position - vertices[i+(n+1)-1].position), kd);
        if(c2<n && c1<m) vertices[i].addShear(&vertices[i+(n+1)+1], 120.0f, length(vertices[i].position - vertices[i+(n+1)+1].position), kd);
        if(c2>1) vertices[i].addShear(&vertices[i-2], 60.0f, length(vertices[i].position - vertices[i-2].position), kd);
        if(c1<m-1) vertices[i].addShear(&vertices[i+2*(n+1)], 60.0f, length(vertices[i].position - vertices[i+2*(n+1)].position), kd);
        if(c1==m) vertices[i].isFixed = true;
    }
    boneArray[0].delta_x = length(vertices[0].position - vertices[1].position) / 20.0f;
}


void updateScene(float t) {
    // const float dt = 0.1f;
    t = 0.0005f;
    for(int i=0;i<nv;i++) vertices[i].updateCurrentForces();
    for(int i = 0; i < boneArray.size(); ++i) {
        boneArray[i].updateGridAndApplyForces(vertices, t, i);
    }
    for(int i=0;i<nv;i++)
    {
        vertices[i].updateGenCords(boneArray, t);
        cloth.vertices[i].position = vertices[i].position;
    }
    boneArray[0].getMeshAttribs(cloth);
    for(int i=0;i<nv;i++)
    {
        vertices[i].normal = boneArray[0].renderNormals[i];
    }
    vec3 verticesData[nv], normalsData[nv];
    for(int i = 0; i < nv; ++i) {
        verticesData[i] = vertices[i].position;
        normalsData[i] = vertices[i].normal;
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
	camCtl.camera.setCameraView(vec3(0.0f, 0.0f, 5.0f), vec3(0.0f, 0.0f, 0.0f), vec3(0.0, 1.0, 0.0));

	program = r.createShaderProgram(
		r.vsBlinnPhong(),
		r.fsBlinnPhong()
	);

	initializeScene();

    float tprev = SDL_GetTicks64()*1e-3;

	while (!r.shouldQuit()) {
        float t = SDL_GetTicks64()*1e-3;

        updateScene(t - tprev);
        tprev = t;

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
