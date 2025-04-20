#ifndef _HALFEDGE_H_
#define _HALFEDGE_H_

//Half edge data structure

#define GLM_ENABLE_EXPERIMENTAL
// #include "viewer.hpp"
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtx/polar_coordinates.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/epsilon.hpp>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <algorithm>
#include <random>
#include <sstream>

class Bone;
class Vertex;
class HalfEdge;
class MeshFace;
class MeshVertex;
class Mesh;

#include "animation.hpp"

using namespace glm;


class HalfEdge{
    public:
        HalfEdge *next;
        HalfEdge *twin;
        int vertexIndex;
        MeshFace* face;
        HalfEdge(int);
};

class MeshFace{
    public:
        HalfEdge *edge;
        vec3 normal;
        MeshFace(HalfEdge*,vec3);
        std::vector<int> getFaceVertices();
        std::vector<HalfEdge*> getFaceEdges();
};

class MeshVertex{
    public:
        HalfEdge *edge;
        vec3 position;
        vec3 normal;
        //boundary check
        MeshVertex(HalfEdge*,vec3);
        std::vector<HalfEdge*> getAdjacentFaces();
        std::vector<int> getAdjacentVertices();
};

class Mesh{
    public:
    std::vector<HalfEdge*> halfEdges;
    std::vector<MeshVertex> vertices;
    std::vector<MeshFace> faces;

    std::vector<ivec3> triangles;
    std::vector<ivec2> edges;
    std::vector<vec3> normals;
    std::vector<std::vector<int>> vertexPerFace;

    Mesh(std::vector<HalfEdge*>, std::vector<MeshVertex>, std::vector<MeshFace>);
    void getEdges();
    void triangulate();
    void getMeshAttribs(int &totalVertices, int &numberOfTriangles, int &numberofEdges,
                    std::vector<vec3> &renderVertices, std::vector<ivec3> &renderTriangles,
                    std::vector<ivec2> &renderEdges, std::vector<vec3> &renderNormals);
    void extrudeFace(int faceIndex, float distance);
    void extrudeFace(vec3 point, float distance);
    void extrudeMultiple(std::vector<int> &indices, float dist);
    void addNoise(float maxnoise);
    void printAdjacentVertices();
    void umbrellaOperator(float lambda, int iterations);
    void catmullClarkSubdivision();
    void addMesh(Mesh &m);
    void moveMesh(vec3 direction);
    vec3 getFaceNormal(int faceIndex);
    vec3 getVertexNormal(std::vector<int>);
};


float pointToSegmentDistance(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b);
float pointToTriangleDistance(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c);
float randomFloat(float min, float max);

Mesh parseObjFile(const std::string &filename);


#endif // HALFEDGE_H_
