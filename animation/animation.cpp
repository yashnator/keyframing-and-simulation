#include "animation.hpp"

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
    std::string boneName;
    glm::mat4 offsetMatrix; // We keep inverse of bind pose matrix handy
    glm::mat4 localTransform;
    glm::mat4 globalTransform;
    glm::mat4 vertTransform;
    Bone* parent = nullptr;
    std::vector<Bone*> children;

    Bone(std::string __boneName, glm::mat4 __offsetMatrix, glm::mat4 __localTransform, Bone* __parent);
    Bone(glm::mat4 __offsetMatrix, glm::mat4 __localTransform, Bone* __parent);

    void updateBone(const glm::mat4& transform);
    void updateBone(const glm::vec3& pos, const glm::quat& rot);
};

// We go from bone's space to world space by multiplying with offset matrix
// Then we apply the global transform to get the final position

class Vertex
{
    // Positions and normals are specified wrt the bone
    glm::vec3 position;
    glm::vec3 normal;
    // glm::vec3 texCoords;
    std::vector<int> boneIDs;
    std::vector<float> weights;
    // vertex(glm::vec3 __position, glm::vec3 __normal, glm::vec3 texCoords, vector<int> boneIDs, vector<float> weights);
    Vertex(glm::vec3 __position, glm::vec3 __normal, std::vector<int> __boneIDs, std::vector<float> __weights);
    Vertex(glm::vec3 __position, glm::vec3 __normal, std::vector<int> __boneIDs);
};

// ------------------- BONE METHODS --------------- //

// ctors

Bone::Bone(std::string __boneName, glm::mat4 __offsetMatrix, glm::mat4 __localTransform, Bone* __parent):
    boneName(__boneName),
    offsetMatrix(__offsetMatrix),
    localTransform(__localTransform),
    parent(__parent)
{ }

Bone::Bone(glm::mat4 __offsetMatrix, glm::mat4 __localTransform, Bone* __parent = nullptr) {
    Bone("Unnamed bone", offsetMatrix, __localTransform, parent);
}

// methods

void Bone::updateBone(const glm::mat4& transform) {
    localTransform = transform;
    if (parent)
        globalTransform = parent->globalTransform * localTransform;
    else
        globalTransform = localTransform;
    vertTransform = globalTransform * offsetMatrix;
}

void Bone::updateBone(const glm::vec3& pos, const glm::quat& rot) {
    updateBone(glm::translate(glm::mat4(1.0f), pos) * glm::toMat4(rot));
}

// ------------------ VERTEX METHODS ------------- //

// ctors

Vertex::Vertex(glm::vec3 __position, glm::vec3 __normal, std::vector<int> __boneIDs, std::vector<float> __weights):
    position(__position),
    normal(__normal),
    boneIDs(__boneIDs),
    weights(__weights)
{ }

Vertex::Vertex(glm::vec3 __position, glm::vec3 __normal, std::vector<int> __boneIDs):
    position(__position),
    normal(__normal),
    boneIDs(__boneIDs)
{ }

// methods
