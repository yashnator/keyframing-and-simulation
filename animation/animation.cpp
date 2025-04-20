#include "animation.hpp"

// ------------------- BONE METHODS --------------- //

// ctors

Bone::Bone(std::string __boneName, glm::mat4 __offsetMatrix, glm::mat4 __localTransform, Bone* __parent):
    boneName(__boneName),
    offsetMatrix(__offsetMatrix),
    localTransform(__localTransform),
    parent(__parent)
{ }

Bone::Bone(glm::mat4 __offsetMatrix, glm::mat4 __localTransform, Bone* __parent = nullptr):
    boneName("Unnamed bone"),
    offsetMatrix(__offsetMatrix),
    localTransform(__localTransform),
    parent(__parent)
{ }

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
    for(auto &vert: renderVertices) {
        vert = glm::vec3(transform * glm::vec4(vert, 1.0f));
    }
    for(auto &nrml: renderNormals) {
        nrml = glm::vec3(transform * glm::vec4(nrml, 0.0f));
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
{ }

// methods

void Vertex::attachToBone(int&& BoneID) noexcept {
    boneIDs.emplace_back(BoneID);
}
