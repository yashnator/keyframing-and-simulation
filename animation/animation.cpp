#include "animation.hpp"

// ------------------- BONE METHODS --------------- //

// ctors

Bone::Bone(std::string __boneName, Bone* __parent, glm::mat4 __offsetMatrix, glm::mat4 __localTransform):
    boneName(__boneName),
    offsetMatrix(__offsetMatrix),
    offsetMatrixIT(inverseTranspose(__offsetMatrix)),
    localTransform(__localTransform),
    globalTransform(mat4(1.0f)),
    parent(__parent),
    vertTransform(mat4(1.0f)),
    vertTransformIT(mat4(1.0f))
{
    // if(parent)
        // offsetMatrix = parent->offsetMatrix * __offsetMatrix;
    vertTransform = offsetMatrix;
    vertTransformIT = inverseTranspose(vertTransform);
    if(parent) {
        parent->children.push_back(this);
    }
}

Bone::Bone(glm::mat4 __offsetMatrix, glm::mat4 __localTransform, Bone* __parent = nullptr):
    boneName("Unnamed bone"),
    offsetMatrix(__offsetMatrix),
    offsetMatrixIT(inverseTranspose(__offsetMatrix)),
    localTransform(__localTransform),
    globalTransform(mat4(1.0f)),
    parent(__parent),
    vertTransform(mat4(1.0f)),
    vertTransformIT(mat4(1.0f))
{ }

// methods

void Bone::updateBone(const glm::mat4& transform) {
    localTransform = transform;
    // if (parent)
        // globalTransform = parent->globalTransform * offsetMatrix * localTransform;
    // else
        // globalTransform = offsetMatrix * localTransform;
    // vertTransform = offsetMatrix * transform * inverse(offsetMatrix);
    // vertTransform = globalTransform;
    // vertTransformIT = inverseTranspose(vertTransform);
}

void Bone::updateBone(const glm::vec3& pos, const glm::quat& rot) {
    updateBone(glm::translate(glm::mat4(1.0f), pos) * glm::toMat4(rot));
}

glm::mat4 Bone::getVertTransform() const {
    return vertTransform;
}

void Bone::getMeshAttribs(Mesh &mesh) {
    mesh.getMeshAttribs(totalVertices,
                            numberOfTriangles,
                            numberOfEdges,
                            renderVertices,
                            renderTriangles,
                            renderEdges,
                            renderNormals);
}

void Bone::updateBindPose(const glm::mat4& transform) {
    glm::mat4 invT = glm::inverseTranspose(transform);
    offsetMatrix = offsetMatrix * transform;
    offsetMatrixIT = inverseTranspose(offsetMatrix);
    // for(auto &vert: renderVertices) {
    //     vert = glm::vec3(transform * glm::vec4(vert, 1.0f));
    // }
    // for(auto &nrml: renderNormals) {
    //     nrml = glm::vec3(invT * glm::vec4(nrml, 0.0f));
    // }
    vertTransform = offsetMatrix;
    vertTransformIT = inverseTranspose(vertTransform);
}

void Bone::updateInit(const glm::mat4 &transform) {
    glm::mat4 invT = glm::inverseTranspose(transform);
    for(auto &vert: renderVertices) {
        vert = glm::vec3(transform * glm::vec4(vert, 1.0f));
    }
    for(auto &nrml: renderNormals) {
        nrml = glm::vec3(invT * glm::vec4(nrml, 0.0f));
    }
}

void Bone::updateAll() {
    std::cout << boneName << std::endl;
    if (parent)
        globalTransform = parent->globalTransform * offsetMatrix * localTransform;
    else
        globalTransform = offsetMatrix * localTransform;
    // vertTransform = offsetMatrix * transform * inverse(offsetMatrix);
    vertTransform = globalTransform;
    vertTransformIT = inverseTranspose(vertTransform);
    for(auto &child: children) {
        child->updateAll();
    }
}

// ------------------ VERTEX METHODS ------------- //

// ctors

Vertex::Vertex() { }

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
{
    for(int i = 0; i < __boneIDs.size(); ++i) {
        weights.emplace_back(1.0f / float(__boneIDs.size()));
    }
}

// methods

void Vertex::attachToBone(int&& BoneID) noexcept {
    boneIDs.emplace_back(BoneID);
}
