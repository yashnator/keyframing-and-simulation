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
    glm::mat4 offsetMatrix; // We keep inverse of bind pose matrix handy
    glm::mat4 localTransform;
    glm::mat4 globalTransform;
    glm::mat4 vertTransform;
    Bone* parent = nullptr;
    std::vector<Bone*> children;

    // public:
        Bone(std::string __boneName, glm::mat4 __offsetMatrix, glm::mat4 __localTransform, Bone* __parent);
        Bone(glm::mat4 __offsetMatrix, glm::mat4 __localTransform, Bone* __parent);

        void updateBone(const glm::mat4& transform);
        void updateBone(const glm::vec3& pos, const glm::quat& rot);
        glm::mat4 getVertTransform() const;
};

// We go from bone's space to world space by multiplying with offset matrix
// Then we apply the global transform to get the final position

class Vertex
{
    public:
    // Positions and normals are specified wrt the bone
    glm::vec3 position;
    glm::vec3 normal;
    // glm::vec3 texCoords;
    std::vector<int> boneIDs;
    std::vector<float> weights;
    // vertex(glm::vec3 __position, glm::vec3 __normal, glm::vec3 texCoords, vector<int> boneIDs, vector<float> weights);

    public:
        Vertex();
        Vertex(glm::vec3 __position, glm::vec3 __normal, std::vector<int> __boneIDs, std::vector<float> __weights);
        Vertex(glm::vec3 __position, glm::vec3 __normal, std::vector<int> __boneIDs);

        void attachToBone(int&& boneID) noexcept;
};

#endif //  __ANIMATIONS__
