#include "animation.hpp"

// ------------------- Shape METHODS --------------- //

Sphere::Sphere(vec3 __center, float __r, vec3 __omega): center(__center), r(__r), omega(__omega)
{ }

Plane::Plane(vec3 __planePt, vec3 __normal): planePt(__planePt), normal(__normal)
{ }

float IdShape::phiSurface(vec3 &pointPrev, vec3 &pointCurr) {
    return 0.0f;
}

void IdShape::updatePosition(mat4 &transform) { }

vec3 IdShape::collisionPt(vec3 &pointPrev, vec3 &pointCurr) {
    return vec3(0.0f);
}

vec3 IdShape::getNormal(vec3 &pt) {
    return vec3(0.0f);
}

vec3 IdShape::getVelocityAt(vec3 &pt) {
    return vec3(0.0f);
}

void IdShape::setOmega(vec3 __omega) { }

void IdShape::rotateCord(vec3 &pt, float dt) { }

float Sphere::phiSurface(vec3 &pointPrev, vec3 &pointCurr) {
    // DBG(center)
    // if(distance(center, pointPrev) < r) return -f;
    return distance(center, pointCurr) - r;
}

void Sphere::updatePosition(mat4 &transform) {
    center = vec3(transform * vec4(center, 1.0f));
}

vec3 Sphere::collisionPt(vec3 &pointPrev, vec3 &pointCurr) {
    // r = prev + t * (curr - prev)
    vec3 d = pointCurr - pointPrev;
    glm::vec3 object_space_c = center;
    float qa = glm::dot(d, d);
    float qb = glm::dot(2.0f * d, pointPrev - object_space_c);
    float qc = glm::dot(pointPrev - object_space_c, pointPrev - object_space_c) - r * r;

    float discriminant = qb * qb - 4 * qa * qc;
    if (discriminant < 0) {
        std::cout << phiSurface(pointPrev, pointCurr) << std::endl;
        std::cout << phiSurface(pointPrev, pointPrev) << std::endl;
        std::cout << discriminant << std::endl;
        DBG(pointPrev)
        DBG(pointCurr)
        return pointPrev;
        assert(0 == 1);
    }

    float t1 = (-qb - sqrt(discriminant)) / (2 * qa);
    float t2 = (-qb + sqrt(discriminant)) / (2 * qa);

    if(t1>t2) std::swap(t1, t2);
    if(t1>0)
    {
        return pointPrev + t1 * (d);
    }
    if(t2>0)
    {
        // assert(0 == 1);
        return pointPrev + t2 * (d);
    }
    assert(0 == 1);
}

vec3 Sphere::getNormal(vec3 &pt) {
    return normalize(pt - center);
}

vec3 Sphere::getVelocityAt(vec3 &pt) {
    return cross(omega, pt - center);
}

void Sphere::setOmega(vec3 __omega) {
    omega = __omega;
}

void Sphere::rotateCord(vec3 &pt, float dt) {
    vec3 rel_pt = pt - center;
    float angle = length(omega) * dt;
    rel_pt = rotate(mat4(1.0f), angle, normalize(omega)) * vec4(rel_pt, 1.0f);
    pt = center + rel_pt;
}

float Plane::phiSurface(vec3 &pointPrev, vec3  &pointCurr) {
    if(dot(pointPrev - planePt, normal) * dot(pointCurr - planePt, normal) >= 0.0f) {
        return 0.0f;
    }
    return -abs(dot(pointCurr - planePt, normal));
}

void Plane::updatePosition(mat4 &transform) {
    planePt = vec3(transform * vec4(planePt, 1.0f));
    normal = vec3(inverseTranspose(transform) * vec4(normal, 0.0f));
    normal = normalize(normal);
}

vec3 Plane::collisionPt(vec3 &pointPrev, vec3 &pointCurr) {
   return vec3(0.0f);
}

vec3 Plane::getNormal(vec3 &pt) {
    return normalize(normal);
}

vec3 Plane::getVelocityAt(vec3 &pt) {
    return vec3(0.0f);
}

void Plane::setOmega(vec3 __omega) { }

void Plane::rotateCord(vec3 &pt, float dt) { }

// ------------------- BONE METHODS --------------- //

// ctors

Bone::Bone(std::string __boneName, Bone* __parent, std::unique_ptr<Mesh> __mesh, std::unique_ptr<Shape> __shape, glm::mat4 __offsetMatrix, glm::mat4 __localTransform, vec3 __velocity) :
    boneName(__boneName),
    offsetMatrix(__offsetMatrix),
    offsetMatrixIT(inverseTranspose(__offsetMatrix)),
    localTransform(__localTransform),
    globalTransform(mat4(1.0f)),
    parent(__parent),
    vertTransform(mat4(1.0f)),
    vertTransformIT(mat4(1.0f)),
    shape(std::move(__shape)),
    mesh(std::move(__mesh)),
    isFixed(false),
    mu(0.0f),
    epsilon(1.0f),
    velocity(__velocity)
{
    vertTransform = offsetMatrix;
    vertTransformIT = inverseTranspose(vertTransform);
    omega = vec3(0.0f);
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

void Bone::updateMesh() {
    if(mesh) {
        for(int i = 0; i < renderVertices.size(); ++i) {
            mesh->vertices[i].position = renderVertices[i];
        }
        // mesh->vertices = renderVertices;
        mesh->getMeshAttribs(totalVertices,
                            numberOfTriangles,
                            numberOfEdges,
                            renderVertices,
                            renderTriangles,
                            renderEdges,
                            renderNormals);
    }
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

void Bone::updateInit(glm::mat4 &&transform) {
    glm::mat4 invT = glm::inverseTranspose(transform);
    for(auto &vert: renderVertices) {
        vert = glm::vec3(transform * glm::vec4(vert, 1.0f));
    }
    for(auto &nrml: renderNormals) {
        nrml = glm::vec3(invT * glm::vec4(nrml, 0.0f));
    }
    if(shape) shape->updatePosition(transform);
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

void Bone::setOmega(vec3 __omega) {
    shape->setOmega(__omega);
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

    // if(length(vij) < 0.01f) vij = vec3(0.0f);

    // std::cout << length(xij) - springData[1] << std::endl;
    float springForce = -1.0f * springData[0] * (length(xij) - springData[1]);

    float dampForce = -springData[2] * (dot(vij, xij_dir));
    // if(std::abs(length(xij)-springData[1]) < 0.01f) springForce = 0.0f;

    return (springForce + dampForce) * xij_dir;
}

void Vertex::updateCurrentForces() {
    // DBG(position)
    if(isFixed) return;
    assert(mass != 0);
    vec3 totalForce = mass * g;
    // vec3 totalForce = vec3(0.0f);
    for(auto &other: structuralSprings) {
        totalForce += getForce(this, other.first, other.second);
    }
    for(auto &other: shearSprings) {
        totalForce += getForce(this, other.first, other.second);
    }
    for(auto &other: bendingSprings) {
        totalForce += getForce(this, other.first, other.second);
    }
    currentForce = totalForce;
}

void Vertex::updateGenCords(std::vector<Bone> &bones, float dt) {
    if(bones[boneIDs.front()].shape) {
        velocity -= bones[boneIDs.front()].shape->getVelocityAt(position);
        bones[boneIDs.front()].shape->rotateCord(position, dt);
        velocity += bones[boneIDs.front()].shape->getVelocityAt(position);
           // position = posNext;
    }
    if(isFixed) {
       return;
    }
    // std::cout << "here" << std::endl;
    assert(mass != 0);
    vec3 velNext = velocity + (currentForce / mass) * dt;
    // if(bones[boneIDs.front()].shape) {
    //     velNext -= bones[boneIDs.front()].shape->getVelocityAt(position);
    // }
    // velNext -= angular_del_v;
    vec3 posNext = position + velNext * dt;
    // if(bones[boneIDs.front()].shape) {
    //     DBG(bones[boneIDs.front()].shape->getVelocityAt(posNext))
    //     velNext += bones[boneIDs.front()].shape->getVelocityAt(posNext);
    // }
    // velNext += bones[boneIDs.front()].shape->getVelocityAt(posNext);
    for(auto &bone: bones) {
        if(bone.shape == nullptr) continue;
        float phi = bone.shape->phiSurface(position, posNext);
        if(phi < 0.0f) {
            vec3 n = bone.shape->getNormal(position);
            vec3 vn = dot(velNext, n) * n;
            vec3 vt = velNext - vn;
            float t1 = phi / length(vn);
            float t2 = dt - t1;
            posNext = position + t1 * velNext;
            // posNext is the point of collision
            // Now update velocity
            // vec3 v_sphere = bone.velocity + bone.shape->getVelocityAt(posNext);
            vec3 v_sphere = bone.velocity + bone.shape->getVelocityAt(posNext);
            vec3 v_rel = velNext - v_sphere;
            vec3 v_reln = dot(v_rel, n) * n;
            vec3 v_relt = v_rel - v_reln;
            // DBG(v_sphere)

            velNext = -bone.epsilon * v_reln;
            vec3 impulse = velNext - v_reln;
            float friction_mag = clamp((bone.mu * length(impulse)) / length(v_relt), 0.0f, 1.0f);
            v_relt = (1 - friction_mag) * v_relt;
            velNext += v_relt + v_sphere;

            // DBG(velocity)
            // DBG(velNext)

            //
            posNext = posNext + velNext * t2;
            // vec3 collisionPt = bone.shape->collisionPt(position, posNext);
            // float t1 = length(collisionPt - position) / length(velocity);
            // float t2 = dt - t1;
            // if(dot(vn, n) < 0.0f)  {
            //     velNext = vt - vn;
            //     posNext = collisionPt + t2 * velNext;
            // }
        }
    }
    velocity = velNext;
    position = posNext;
    // DBG(position)
}
