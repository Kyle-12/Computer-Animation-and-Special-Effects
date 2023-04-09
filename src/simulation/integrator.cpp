#include "integrator.h"
#include <iostream>

#include <vector>
namespace simulation {
// Factory
std::unique_ptr<Integrator> IntegratorFactory::CreateIntegrator(IntegratorType type) {
    switch (type) {
        case simulation::IntegratorType::ExplicitEuler:
            return std::make_unique<ExplicitEulerIntegrator>();
        case simulation::IntegratorType::ImplicitEuler:
            return std::make_unique<ImplicitEulerIntegrator>();
        case simulation::IntegratorType::MidpointEuler:
            return std::make_unique<MidpointEulerIntegrator>();
        case simulation::IntegratorType::RungeKuttaFourth:
            return std::make_unique<RungeKuttaFourthIntegrator>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}

//
// ExplicitEulerIntegrator
//

IntegratorType ExplicitEulerIntegrator::getType() { return IntegratorType::ExplicitEuler; }

void ExplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {

    // TODO#4-1: Integrate position and velocity
    //   1. Integrate position using current velocity.
    //   2. Integrate velocity using current acceleration.
    //   3. Clear force
    // Note:
    //   1. You should do this first. Then you can check whether your collision is correct or not.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. Review ¡§ODE_basics.pptx¡¨ from p.15 - p.16

    for (int i = 0; i < particleSystem.jellyCount; i++) {
        for (int j = 0; j < particleSystem.jellies[i].getParticleNum(); j++) {
            Particle& curParticle = particleSystem.jellies[i].getParticle(j);
          
            curParticle.setPosition(curParticle.getPosition() + particleSystem.deltaTime * curParticle.getVelocity());
            curParticle.setVelocity(curParticle.getVelocity() + particleSystem.deltaTime * curParticle.getAcceleration());
            curParticle.setForce(Eigen::Vector3f::Zero());
        }
    }
}   

//
// ImplicitEulerIntegrator
//

IntegratorType ImplicitEulerIntegrator::getType() { return IntegratorType::ImplicitEuler; }

void ImplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {

    // TODO#4-2: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 and Vn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1
    //   2. Review ¡§ODE_basics.pptx¡¨ from p.18 - p.19
    
    struct backup_particle {
        Eigen::Vector3f position;
        Eigen::Vector3f velocity;
    };

    std::vector<backup_particle> backup;

    for (int i = 0; i < particleSystem.jellyCount; i++) {
        backup.clear();
        for (int j = 0; j < particleSystem.jellies[i].getParticleNum(); j++) {
            Particle& curParticle = particleSystem.jellies[i].getParticle(j);
            backup_particle  bp;
            
            bp.position = curParticle.getPosition();
            bp.velocity = curParticle.getVelocity();
            backup.push_back(bp);

            curParticle.setForce(Eigen::Vector3f::Zero());
        }

        particleSystem.computeJellyForce(particleSystem.jellies[i]);

        for (int j = 0; j < particleSystem.jellies[i].getParticleNum(); j++) {
            Particle& curParticle = particleSystem.jellies[i].getParticle(j);

            curParticle.setPosition(backup[j].position + particleSystem.deltaTime * curParticle.getVelocity());
            curParticle.setVelocity(backup[j].velocity + particleSystem.deltaTime * curParticle.getAcceleration());
                        
            curParticle.setForce(Eigen::Vector3f::Zero());
        }
      
    }

    


    

}

//
// MidpointEulerIntegrator
//

IntegratorType MidpointEulerIntegrator::getType() { return IntegratorType::MidpointEuler; }

void MidpointEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-3: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1.
    //   2. Review ¡§ODE_basics.pptx¡¨ from p .18 - p .20
    struct backup_particle {
        Eigen::Vector3f position;
        Eigen::Vector3f velocity;
    };

    std::vector<backup_particle> backup;

    for (int i = 0; i < particleSystem.jellyCount; i++) {
        backup.clear();
        for (int j = 0; j < particleSystem.jellies[i].getParticleNum(); j++) {
            Particle& curParticle = particleSystem.jellies[i].getParticle(j);
            backup_particle  bp;

            bp.position = curParticle.getPosition();
            bp.velocity = curParticle.getVelocity();
            backup.push_back(bp);

            
            Eigen::Vector3f eulerPosition = particleSystem.deltaTime * curParticle.getVelocity();
            Eigen::Vector3f eulerVelocity = particleSystem.deltaTime * curParticle.getAcceleration();

            curParticle.setVelocity(curParticle.getVelocity() + 1 / 2 * eulerVelocity);
            curParticle.setPosition(curParticle.getPosition() + 1 / 2 * eulerPosition);
            curParticle.setForce(Eigen::Vector3f::Zero());
        }

        particleSystem.computeJellyForce(particleSystem.jellies[i]);

        for (int j = 0; j < particleSystem.jellies[i].getParticleNum(); j++) {
            Particle& curParticle = particleSystem.jellies[i].getParticle(j);

            curParticle.setPosition(backup[j].position + particleSystem.deltaTime * curParticle.getVelocity());
            curParticle.setVelocity(backup[j].velocity + particleSystem.deltaTime * curParticle.getAcceleration());
            
            curParticle.setForce(Eigen::Vector3f::Zero());
        }

    }
}

//
// RungeKuttaFourthIntegrator
//

IntegratorType RungeKuttaFourthIntegrator::getType() { return IntegratorType::RungeKuttaFourth; }

void RungeKuttaFourthIntegrator::integrate(MassSpringSystem& particleSystem) {
    
    // TODO#4-4: Integrate velocity and acceleration
    //   1. Backup original particles' data.
    //   2. Compute k1, k2, k3, k4
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1.
    //   2. StateStep struct is just a hint, you can use whatever you want.
    //   3. Review ¡§ODE_basics.pptx¡¨ from p.21
    struct StateStep {
        Eigen::Vector3f deltaVel;
        Eigen::Vector3f deltaPos;
    };

    struct backup_particle {
        Eigen::Vector3f position;
        Eigen::Vector3f velocity;
    };

    std::vector<StateStep> k1bp, k2bp, k3bp ,k4bp;
    std::vector<backup_particle> backup;

    for (int i = 0; i < particleSystem.jellyCount; i++) {

        backup.clear();
        k1bp.clear();
        k2bp.clear();
        k3bp.clear();
        k4bp.clear();

        for (int j = 0; j < particleSystem.jellies[i].getParticleNum(); j++) {
            Particle& curParticle = particleSystem.jellies[i].getParticle(j);
            backup_particle  ko;

            ko.position = curParticle.getPosition();
            ko.velocity = curParticle.getVelocity();

            backup.push_back(ko);

        }

        //k1
        for (int j = 0; j < particleSystem.jellies[i].getParticleNum(); j++) {
            Particle& curParticle = particleSystem.jellies[i].getParticle(j);
            StateStep  k1;

            k1.deltaPos = particleSystem.deltaTime * curParticle.getVelocity();
            k1.deltaVel = particleSystem.deltaTime * curParticle.getAcceleration();
            k1bp.push_back(k1);

            curParticle.setPosition(backup[j].position + k1.deltaPos / 2);
            curParticle.setVelocity(backup[j].velocity + k1.deltaVel / 2);

        }
       
        particleSystem.computeJellyForce(particleSystem.jellies[i]);

        // k2
        for (int j = 0; j < particleSystem.jellies[i].getParticleNum(); j++) {
            Particle& curParticle = particleSystem.jellies[i].getParticle(j);
            StateStep  k2;

            k2.deltaPos = particleSystem.deltaTime * curParticle.getVelocity();
            k2.deltaVel = particleSystem.deltaTime * curParticle.getAcceleration();
            k2bp.push_back(k2);

            curParticle.setPosition(backup[j].position + k2.deltaPos / 2);
            curParticle.setVelocity(backup[j].velocity + k2.deltaVel / 2);

        }
        
        particleSystem.computeJellyForce(particleSystem.jellies[i]);

        // k3
        for (int j = 0; j < particleSystem.jellies[i].getParticleNum(); j++) {
            Particle& curParticle = particleSystem.jellies[i].getParticle(j);
            StateStep  k3;

            k3.deltaPos = particleSystem.deltaTime * curParticle.getVelocity();
            k3.deltaVel = particleSystem.deltaTime * curParticle.getAcceleration();
            k3bp.push_back(k3);

            curParticle.setVelocity(backup[j].velocity + k3.deltaVel);
            curParticle.setPosition(backup[j].position + k3.deltaPos);
        }
        
        particleSystem.computeJellyForce(particleSystem.jellies[i]);

        // k4
        for (int j = 0; j < particleSystem.jellies[i].getParticleNum(); j++) {
            Particle& curParticle = particleSystem.jellies[i].getParticle(j);
            StateStep  k4;

            k4.deltaPos = particleSystem.deltaTime * curParticle.getVelocity();
            k4.deltaVel = particleSystem.deltaTime * curParticle.getAcceleration();
            k4bp.push_back(k4);
        }
        
        // calculate
        for (int j = 0; j < particleSystem.jellies[i].getParticleNum(); j++) {
            Particle& curParticle = particleSystem.jellies[i].getParticle(j);

            curParticle.setVelocity(backup[j].velocity
                + (k1bp[j].deltaVel + 2 * k2bp[j].deltaVel + 2 * k3bp[j].deltaVel + k4bp[j].deltaVel) / 6);

            curParticle.setPosition(backup[j].position
                + (k1bp[j].deltaPos + 2 * k2bp[j].deltaPos + 2 * k3bp[j].deltaPos + k4bp[j].deltaPos) / 6);

            curParticle.setForce(Eigen::Vector3f::Zero());
        }

    }

    


}

}  // namespace simulation
