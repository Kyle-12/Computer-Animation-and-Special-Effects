#include "terrain.h"

#include <stdexcept>

#include "../util/helper.h"
#include <iostream>

namespace simulation {
// Factory
std::unique_ptr<Terrain> TerrainFactory::CreateTerrain(TerrainType type) {
    switch (type) {
        case simulation::TerrainType::Plane:
            return std::make_unique<PlaneTerrain>();
        case simulation::TerrainType::Bowl:
            return std::make_unique<BowlTerrain>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}
// Terrain

Eigen::Matrix4f Terrain::getModelMatrix() { return modelMatrix; }

// Note:
// You should update each particles' velocity (base on the equation in
// slide) and force (contact force : resist + friction) in handleCollision function

// PlaneTerrain //

PlaneTerrain::PlaneTerrain() { modelMatrix = util::translate(0.0f, position[1], 0.0f) * util::scale(60, 1, 60); }

TerrainType PlaneTerrain::getType() { return TerrainType::Plane; }

void PlaneTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.2f;
    // TODO#3-1: Handle collision when a particle collide with the plane terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data.
    // Hint:
    //   1. Review "particles.pptx¡¨ from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.

    for (int i = 0; i < jelly.getParticleNum(); i++) {

        Particle& curParticle = jelly.getParticle(i);

        Eigen::Vector3f nor = normal.normalized();

        bool on_plane = hole_radius < (curParticle.getPosition() - hole_position).norm() ;

        if (on_plane && (curParticle.getPosition() - position).dot(normal) < eEPSILON
            && curParticle.getVelocity().dot(normal) < 0) {
            
            //printf("Collision");

            Eigen::Vector3f Vel = curParticle.getVelocity();
            Eigen::Vector3f normalVel = normal * ( Vel.dot(normal) / normal.squaredNorm() );
            curParticle.setVelocity(Vel - normalVel + (-1 * coefResist * normalVel));

            // Contact : p18 (N dot V) < e, p19 force  N dot F < 0

            if ( abs((curParticle.getPosition() - position).dot(normal)) < eEPSILON
                &&abs(normal.dot(Vel - normalVel)) < eEPSILON
                &&curParticle.getForce().dot(normal) < 0 ) {

                Eigen::Vector3f fc = -1 * normal.dot(curParticle.getForce()) * normal;
                Eigen::Vector3f ff = -1 * coefFriction * (-1 * normal.dot(curParticle.getForce())) * (Vel - normalVel);

                curParticle.addForce(fc);
                curParticle.addForce(ff);
            }

        }
        
        
    }
    
}
// BowlTerrain //

BowlTerrain::BowlTerrain() {
    modelMatrix = util::translate(position) * util::rotateDegree(-90, 0, 0) * util::scale(radius, radius, radius);
}

TerrainType BowlTerrain::getType() { return TerrainType::Bowl; }

void BowlTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.2f;
    // TODO#3-2: Handle collision when a particle collide with the sphere terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data. 
    // Hint:
    //   1. Review "particles.pptx¡¨ from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.

    for (int i = 0; i < jelly.getParticleNum(); i++) {

        Particle& curParticle = jelly.getParticle(i);

        Eigen::Vector3f disvector = position - curParticle.getPosition();
        float dis = disvector.norm();
        float ydis = abs(disvector[1]);
        Eigen::Vector3f bowlnormal = disvector / dis;

        if ((ydis - radius / sqrt(2)) > 0) {

            Eigen::Vector3f Vel = curParticle.getVelocity();
            Eigen::Vector3f normalVel = bowlnormal * (curParticle.getVelocity().dot(bowlnormal)) / (bowlnormal.norm() * bowlnormal.norm());

            if ( (radius - dis) < eEPSILON && curParticle.getVelocity().dot(bowlnormal) < 0) {
                curParticle.setVelocity(Vel - normalVel + (-1 * coefResist * normalVel));

            }

            if (abs(radius - dis) < eEPSILON && curParticle.getForce().dot(bowlnormal) < 0 
                && abs(bowlnormal.dot(Vel - normalVel) ) < eEPSILON
                ) {

                Eigen::Vector3f fc = -1 * bowlnormal.dot(curParticle.getForce()) * bowlnormal;
                Eigen::Vector3f ff = -1 * coefFriction * (-1 * bowlnormal.dot(curParticle.getForce())) * (Vel - normalVel);

                curParticle.addForce(fc);
                curParticle.addForce(ff);

            }
            
           
                
        }

    }
    
   
}
}  // namespace simulation

