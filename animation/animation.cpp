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
    // if(parent) {
    //     parent->children.push_back(this);
    // }
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
    // std::cout << boneName << std::endl;
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
    weights(__weights),
    velocity(vec3(0.0f)),
    isFixed(false),
    mass(0.0f)
{ }

Vertex::Vertex(glm::vec3 __position, glm::vec3 __normal, std::vector<int> __boneIDs):
    position(__position),
    normal(__normal),
    boneIDs(__boneIDs),
    velocity(0.0f),
    isFixed(false),
    mass(0.0f)
{
    for(int i = 0; i < __boneIDs.size(); ++i) {
        weights.emplace_back(1.0f / float(__boneIDs.size()));
    }
}

// methods

void Vertex::attachToBone(int&& BoneID) noexcept {
    boneIDs.emplace_back(BoneID);
}


void Vertex::updateNormal() {
    // Calculate normal from the people in structural springs
    vec3 newNormal(0.0f);
    for(int i = 0; i < structuralSprings.size(); ++i) {
        vec3 x = structuralSprings[i].first->position - position;
        vec3 y = structuralSprings[(i + 1) % structuralSprings.size()].first->position - position;
        newNormal += cross(x, y);
    }
    normal = normalize(newNormal);
}

void Vertex::addStructural(Vertex *other, float ks, float l0, float kd) {
    structuralSprings.push_back({other, {ks, l0, kd}});
    other->structuralSprings.push_back({this, {ks, l0, kd}});
}

void Vertex::addShear(Vertex *other, float ks, float l0, float kd) {
    shearSprings.push_back({other, {ks, l0, kd}});
    other->shearSprings.push_back({this, {ks, l0, kd}});
}

void Vertex::addBend(Vertex *other, float ks, float l0, float kd) {
    bendingSprings.push_back({other, {ks, l0, kd}});
    other->bendingSprings.push_back({this, {ks, l0, kd}});
}

vec3 Vertex::getForce(Vertex *curr, Vertex* other, std::array<float, 3> &springData) {
    vec3 xij = curr->position - other->position;
    vec3 xij_dir = normalize(xij);
    vec3 vij = curr->velocity - other->velocity;
    float springForce = -springData[0] * (length(xij) - springData[0]);
    float dampForce = -springData[2] * (dot(vij, xij_dir));
    return (springForce + dampForce) * xij_dir;
}

void Vertex::updateGenCords(float dt) {
    // if(isFixed) DBG(position)
    if(isFixed) return;
    assert(mass != 0);
    vec3 totalForce = const_force;
    for(auto &other: structuralSprings) {
        totalForce += getForce(this, other.first, other.second);
    }
    for(auto &other: shearSprings) {
        totalForce += getForce(this, other.first, other.second);
    }
    for(auto &other: bendingSprings) {
        totalForce += getForce(this, other.first, other.second);
    }
    // DBG(totalForce)
    velocity = velocity + (totalForce / mass) * dt;
    position = position + velocity * dt;
    // DBG(position)
    // DBG(velocity * dt)
}
