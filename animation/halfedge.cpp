#include "halfedge.hpp"

/*
Functions for half edge data structure
*/

#define _AREA_WEIGHTED_

HalfEdge::HalfEdge(int index)
{
    this->next=nullptr;
    this->twin=nullptr;
    this->vertexIndex=index;
    this->face=nullptr;
}

/*
Functions for MeshFace
*/

MeshFace::MeshFace(HalfEdge *edge, vec3 normal)
{
    this->edge=edge;
    this->normal=normal;
}

std::vector<int> MeshFace::getFaceVertices()
{
    std::vector<int> faceVertices;
    HalfEdge *start = this->edge;
    HalfEdge *current = start;
    do
    {
        faceVertices.emplace_back(current->vertexIndex);
        current=current->next;

    } while (current!=start);

    return faceVertices;
}

std::vector<HalfEdge *> MeshFace::getFaceEdges()
{
    std::vector<HalfEdge *> faceVertices;
    HalfEdge *start = this->edge;
    HalfEdge *current = start;

    do
    {
        faceVertices.emplace_back(current);
        current=current->next;
    } while (current!=start);

    return faceVertices;
}

/*
Functions for MeshVertex
*/

MeshVertex::MeshVertex(HalfEdge *edge, vec3 position)
{
    this->edge=edge;
    this->position=position;
}

std::vector<HalfEdge*> MeshVertex::getAdjacentFaces()
{
    std::vector<HalfEdge*> adjacentFaces;
    HalfEdge *start = this->edge;
    HalfEdge *current = start;

    if(!current)
    {
        std::cout<<"No edge associated with the vertex"<<std::endl;
    }

    do
    {
        adjacentFaces.emplace_back(current);
        if(current->twin==NULL)
        {
            //boundary face encountered
            std::cout<<"some issue with the code please check"<<std::endl;
            break;
        }
        else
        {
            current = current->twin->next;
        }

    }while(current!=start);
    return adjacentFaces;
}

std::vector<int> MeshVertex::getAdjacentVertices()
{
    std::vector<int> adjacentVertices;
    HalfEdge *start = this->edge;
    HalfEdge *current = start;

    if(!current)
    {
        std::cout<<"No edge associated with the vertex"<<std::endl;
    }

    do
    {
        if(current->twin==NULL)
        {
            //boundary face encountered
            std::cout<<"some issue with the code please check"<<std::endl;
            break;
        }
        else
        {
            adjacentVertices.emplace_back(current->twin->vertexIndex);
            current = current->twin->next;
        }

    }while(current!=start);
    return adjacentVertices;
}

/*
Functions for Mesh
*/

Mesh::Mesh(std::vector<HalfEdge*> halfEdges, std::vector<MeshVertex> vertices, std::vector<MeshFace> faces)
{
    this->halfEdges=halfEdges;
    this->vertices=vertices;
    this->faces=faces;
}

void Mesh::getEdges()
{
    edges.clear();
    for(auto &edge:halfEdges)
    {
        if(edge->twin)
        {
            ivec2 curredge;
            curredge.x = edge->vertexIndex;
            curredge.y = edge->twin->vertexIndex;
            edges.emplace_back(curredge);
        }
        else
        {
            ivec2 curredge;
            curredge.x = edge->vertexIndex;
            curredge.y = edge->next->vertexIndex;
            if(curredge.x!=curredge.y) edges.emplace_back(curredge);
            // edges.emplace_back(curredge);
        }
    }
}

void Mesh::triangulate()
{
    triangles.clear();
    normals.clear();
    vertexPerFace.clear();
    for(auto face:faces)
    {
        std::vector<int> faceVertices = face.getFaceVertices();
        vertexPerFace.emplace_back(faceVertices);
        normals.emplace_back(face.normal);
        int n = faceVertices.size();
        for(int i=0;i<n-2;i++)
        {
            ivec3 triangle;
            triangle.x = faceVertices[0];
            triangle.y = faceVertices[i+1];
            triangle.z = faceVertices[i+2];
            triangles.emplace_back(triangle);
        }
    }
}


void Mesh::getMeshAttribs(int &totalVertices, int &numberOfTriangles, int &numberofEdges,
                    std::vector<vec3> &renderVertices, std::vector<ivec3> &renderTriangles,
                    std::vector<ivec2> &renderEdges, std::vector<vec3> &renderNormals)
{
    this->getEdges();
    this->triangulate();

    for(int i = 0; i < faces.size(); ++i) {
        this->faces[i].normal = getFaceNormal(i);
    }

    totalVertices=vertices.size();
    numberOfTriangles=this->triangles.size();
    numberofEdges=this->edges.size();

    renderVertices.resize(totalVertices);
    renderNormals.resize(totalVertices);
    renderTriangles.resize(numberOfTriangles);
    renderEdges.resize(numberofEdges);

    int count=0;
    int triangleCount=0;
    for(auto x:this->triangles)
    {
        ivec3 currTriangle;
        currTriangle.x=x.x;
        currTriangle.y=x.y;
        currTriangle.z=x.z;
        renderTriangles[triangleCount]=currTriangle;
        triangleCount++;
    }

    //hashmap such that we only store one face normal for each vertex
    std::map<int,int> done;
    //make a vector for normals corresponding to each vertex
    std::vector<vec3> currNormals(totalVertices);

    int currFace=0;
    for(auto x:vertexPerFace)
    {
        for(auto y:x)
        {
            if(done.find(y)==done.end())
            {
                done[y]=count;
                currNormals[count]=faces[currFace].normal;
                count++;
            }
        }
        currFace++;
    }

    //sanity check
    // for(int i=0;i<totalVertices;i++)
    // {
    //     if(done.find(i)==done.end())
    //     {
    //         std::cerr<<"Error in normal allotment\n";
    //         // return;
    //     }
    // }

    std::map<int, std::vector<int>> mp;

    for(int i = 0; i < this->faces.size(); ++i) {
        for(auto p: this->faces[i].getFaceVertices()) {
            mp[p].emplace_back(i);
        }
    }

    #ifdef _AREA_WEIGHTED_
    std::map<int, vec3> wtnormals;
    for(int i = 0; i < this->faces.size(); ++i) {
        std::vector<int> verts = this->faces[i].getFaceVertices();
        int vsz = verts.size();
        for(int j = 0; j < vsz; ++j) {
            vec3 vi = this->vertices[verts[j]].position - this->vertices[verts[(j + 1) % vsz]].position;
            vec3 vin = this->vertices[verts[(j + 2) % vsz]].position - this->vertices[verts[(j + 1) % vsz]].position;
            float denominator = length(vi) * length(vin);
            denominator = denominator * denominator;
            wtnormals[verts[(j + 1) % vsz]] += cross(vin, vi) / denominator;
        }
    }
    #else
    std::map<int, std::vector<int>> mp;
    for(int i = 0; i < this->faces.size(); ++i) {
        for(auto p: this->faces[i].getFaceVertices()) {
            mp[p].emplace_back(i);
        }
    }
    #endif

    count=0;
    for(auto x:this->vertices)
    {
        renderVertices[count]=x.position;
        #ifdef _AREA_WEIGHTED_
        renderNormals[count] = wtnormals[count];
        #else
        renderNormals[count] = getVertexNormal(mp[count]);
        #endif
        count++;
    }
    for(int i=0;i<numberofEdges;i++)
    {
        ivec2 currEdge;
        currEdge.x=this->edges[i].x;
        currEdge.y=this->edges[i].y;
        renderEdges[i]=currEdge;
    }
}

vec3 Mesh::getFaceNormal(int idx)
{
    std::vector<int> faceVertices = (this->faces[idx]).getFaceVertices();

    vec3 normal = vec3(0.0, 0.0, 0.0);
    for(int i = 1; i < faceVertices.size() - 1; ++i) {
        normal += (0.5f) * cross(
            this->vertices[faceVertices[i]].position - this->vertices[faceVertices[0]].position,
            this->vertices[faceVertices[i + 1]].position - this->vertices[faceVertices[0]].position
            );
    }

    normal = normalize(normal);
    return normal;

    // vec3 d1, d2;
    // HalfEdge *start = this->faces[idx].edge;
    // d1=this->vertices[start->next->vertexIndex].position-this->vertices[start->vertexIndex].position;
    // d2=this->vertices[start->next->next->vertexIndex].position-this->vertices[start->next->vertexIndex].position;
    // vec3 normal = cross(d1,d2);
    // return normal;
}

vec3 Mesh::getVertexNormal(std::vector<int> adj) {
    vec3 currnormal(0.0, 0.0, 0.0);
    for(auto x: adj) {
        currnormal += this->faces[x].normal;
    }
    currnormal = normalize(currnormal);
    return currnormal;
}

/*Getting neighbors and performing smoothing by umbrella operator*/
void Mesh::printAdjacentVertices()
{
    for(int i=0;i<vertices.size();i++)
    {
        std::cout<<"Neighbors of vertex "<<i<<std::endl;
        vertices[i].getAdjacentFaces();
    }
}


void Mesh::addMesh(Mesh &m)
{
    int n=this->vertices.size();
    for(auto x:m.vertices)
    {
        this->vertices.emplace_back(x);
    }
    for(auto x:m.faces)
    {
        this->faces.emplace_back(x);
    }
    for(auto x:m.halfEdges)
    {
        x->vertexIndex+=n;
        this->halfEdges.emplace_back(x);
    }
}

void Mesh::moveMesh(vec3 direction)
{
    for(auto &x:this->vertices)
    {
        x.position+=direction;
    }
}

/* The below functions are AI generated*/

float pointToSegmentDistance(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b) {
    glm::vec3 ab = b - a;
    glm::vec3 ap = p - a;
    float t = glm::dot(ap, ab) / glm::dot(ab, ab);
    t = clamp(t, 0.0f, 1.0f);  // Clamp projection to segment
    glm::vec3 closest = a + t * ab;
    return glm::length(p - closest);
}

float pointToTriangleDistance(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) {
    // Compute triangle normal
    glm::vec3 normal = glm::normalize(glm::cross(b - a, c - a));

    // Compute perpendicular projection onto triangle plane
    float d = glm::dot(normal, a);
    float signedDist = glm::dot(normal, p) - d;
    glm::vec3 projectedPoint = p - signedDist * normal;

    // Check if projected point is inside the triangle using Barycentric coordinates
    glm::vec3 v0 = b - a, v1 = c - a, v2 = projectedPoint - a;
    float d00 = glm::dot(v0, v0);
    float d01 = glm::dot(v0, v1);
    float d11 = glm::dot(v1, v1);
    float d20 = glm::dot(v2, v0);
    float d21 = glm::dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;

    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;

    if (u >= 0.0f && v >= 0.0f && w >= 0.0f) {
        // Projection falls inside the triangle
        return std::abs(signedDist);
    }

    // Otherwise, find the minimum distance to edges and vertices
    float edgeDist1 = pointToSegmentDistance(p, a, b);
    float edgeDist2 = pointToSegmentDistance(p, b, c);
    float edgeDist3 = pointToSegmentDistance(p, c, a);

    return std::min({ edgeDist1, edgeDist2, edgeDist3 });
}

float randomFloat(float min, float max) {
    static std::random_device rd;  // Non-deterministic random number generator
    static std::mt19937 gen(rd()); // Mersenne Twister PRNG
    std::uniform_real_distribution<float> dist(min, max);
    return dist(gen);
}


Mesh parseObjFile(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return Mesh{{}, {}, {}};
    }

    std::vector<HalfEdge*> hedges;
    std::vector<MeshVertex> verts;
    std::vector<MeshFace> faces;

    std::vector<glm::vec3> vertices,  normals;
    std::vector<glm::vec2> texCoords;

    std::map<std::pair<int, int>, HalfEdge*> settings;

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string op;
        iss >> op;

        if (op == "v") {
            glm::vec3 v;
            iss >> v.x >> v.y >> v.z;
            verts.push_back(MeshVertex(nullptr, v));
            // vertices.push_back(v);
        } else if (op == "vt") {
            glm::vec3 vt;
            iss >> vt.x >> vt.y;
            texCoords.push_back(vt);
        } else if (op == "vn") {
            glm::vec3 vn;
            iss >> vn.x >> vn.y >> vn.z;
            normals.push_back(vn);
        } else if (op == "f") {
            std::string vertexData;
            std::vector<int> vidx, tidx, nidx;
            while (iss >> vertexData) {
                std::istringstream viss(vertexData);
                std::string token;
                int indices[3] = {0, 0, 0};
                int i = 0;
                while (std::getline(viss, token, '/') && i < 3) {
                    if (!token.empty()) {
                        indices[i] = std::stoi(token);
                    }
                    i++;
                }
                vidx.push_back(indices[0] - 1);
                tidx.push_back(indices[1] - 1);
                nidx.push_back(indices[2] - 1);
            }
            // faces.push_back(face);
            for(int i = 0; i < vidx.size(); ++i) {
                if(nidx[i] != -1) {
                    verts[vidx[i]].normal = normals[nidx[i]];
                }
            }
            for(int i = 0; i < vidx.size(); ++i) {
                HalfEdge* currhe = new HalfEdge(vidx[i]);
                hedges.emplace_back(currhe);
                if(verts[vidx[i]].edge == nullptr) {
                    verts[vidx[i]].edge = currhe;
                }
            }
            for(int i = 0; i < vidx.size(); ++i) {
                int idx = hedges.size() - vidx.size() + i;
                int nextidx = hedges.size() - vidx.size() + ((i + 1) % vidx.size());
                hedges[idx]->next = hedges[nextidx];
                settings[{hedges[idx]->vertexIndex, hedges[idx]->next->vertexIndex}] = hedges[idx];
                if(settings.find({hedges[idx]->next->vertexIndex, hedges[idx]->vertexIndex}) != settings.end()) {
                    HalfEdge *tw = settings[{hedges[idx]->next->vertexIndex, hedges[idx]->vertexIndex}];
                    hedges[idx]->twin = tw;
                    tw->twin = hedges[idx];
                }
            }
            faces.emplace_back(MeshFace(hedges.back(), vec3(0.0, 0.0, 0.0)));
            for(int i = 0; i < vidx.size(); ++i) {
                int idx = hedges.size() - vidx.size() + i;
                hedges[idx]->face = &faces.back();
            }
        }
    }
    file.close();
    Mesh newMesh = Mesh(hedges, verts, faces);
    for(int i = 0; i < newMesh.faces.size(); ++i) {
        newMesh.faces[i].normal = newMesh.getFaceNormal(i);
    }
    return newMesh;
}
