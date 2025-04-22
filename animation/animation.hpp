#ifndef  __ANIMATIONS__
#define  __ANIMATIONS__

#define GLM_ENABLE_EXPERIMENTAL

#include <glad/gl.h>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>

#include <SDL2/SDL.h>
#include <string>
#include <vector>

#include "halfedge.hpp"
#include "shapes.hpp"
#include "numerical_methods.hpp"

using namespace glm;

#define DBG(x) std::cout << #x << " " << to_string(x) << std::endl;

const vec3 g(0.0f, -9.81f, 0.0f);

class Bone
{
    /* The member functions are
    ** boneName: Name bones to have extra info while debugging
    ** offsetMatrix: inverse of bindPostMatrix to go from vertex
    ** localTransform: transform from parent to bone
    ** finalTransform: iteratively calculated global transform
    ** vertTransform: what I send to vertex shader, calculate once per bone!!
    ** parent: Pointer to bone's parent
    ** children: pointers to children
    */
    public:
    std::string boneName;
    glm::mat4 offsetMatrix, offsetMatrixIT; // We keep inverse of bind pose matrix handy
    glm::mat4 localTransform;
    glm::mat4 globalTransform;
    glm::mat4 vertTransform, vertTransformIT;
    Bone* parent;
    std::vector<Bone*> children;


    int totalVertices;
    int numberOfTriangles;
    int numberOfEdges;
    std::vector<vec3> renderVertices;
    std::vector<ivec3> renderTriangles;
    std::vector<ivec2> renderEdges;
    std::vector<vec3> renderNormals;
    std::pair<int, int> globalRange;

    // public:
    Bone(std::string __boneName,
         Bone* __parent = nullptr,
         glm::mat4 __offsetMatrix = mat4(1.0f),
         glm::mat4 __localTransform = mat4(1.0f));
    Bone(glm::mat4 __offsetMatrix, glm::mat4 __localTransform, Bone* __parent);

    void updateBone(const glm::mat4& transform);
    void updateBone(const glm::vec3& pos, const glm::quat& rot);
    glm::mat4 getVertTransform() const;

    void updateBindPose(const glm::mat4& transform);
    void updateInit(const glm::mat4& transform);

    void getMeshAttribs(Mesh &mesh);
    void updateAll();
};

// We go from bone's space to world space by multiplying with offset matrix
// Then we apply the global transform to get the final position

class Vertex
{
    public:
    // Positions and normals are specified wrt the bone
    glm::vec3 position;
    glm::vec3 velocity;
    glm::vec3 normal;
    float mass;
    bool isFixed;
    // glm::vec3 texCoords;
    std::vector<int> boneIDs;
    std::vector<float> weights;
    vec3 currentForce = vec3(0.0f);
    // vertex(glm::vec3 __position, glm::vec3 __normal, glm::vec3 texCoords, vector<int> boneIDs, vector<float> weights);

    // Spring mass related things
    // These store <vertex, <ks, l0, kd>>
    // Convention: store adjacent springs in right handed fashion for easy normal computation
    std::vector<std::pair<Vertex*, std::array<float, 3>>> structuralSprings;
    std::vector<std::pair<Vertex*, std::array<float, 3>>> shearSprings;
    std::vector<std::pair<Vertex*, std::array<float, 3>>> bendingSprings;

    Vertex();
    Vertex(glm::vec3 __position, glm::vec3 __normal, std::vector<int> __boneIDs, std::vector<float> __weights);
    Vertex(glm::vec3 __position, glm::vec3 __normal, std::vector<int> __boneIDs);

    void attachToBone(int&& boneID) noexcept;
    vec3 getForce(Vertex *, Vertex *, std::array<float, 3> &);
    void updateCurrentForces();
    void updateGenCords(float t);
    void updateNormal();
    void addStructural(Vertex *other, float ks = 0.0f, float l0 = 0.0f, float kd = 0.0f);
    void addShear(Vertex *other, float ks = 0.0f, float lo = 0.0f, float kd = 0.0f);
    void addBend(Vertex *other, float ks = 0.0f, float l0 = 0.0f, float kd = 0.0f);
};

#endif //  __ANIMATIONS__
