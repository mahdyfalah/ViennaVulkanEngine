#pragma once

#include "VHInclude.h"
#include "VEInclude.h"
#include "VPE.hpp"

/**
 * @brief Car class for player-controlled vehicle
 * Handles arcade-style driving physics with simple controls
 */
class Car {
public:
    Car() = default;
    ~Car() = default;

    /**
     * @brief Create a car at specified position
     * @param engine Reference to the game engine
     * @param physics Reference to the physics world
     * @param registry Reference to the entity registry
     * @param position Starting position
     */
    void Create(vve::Engine& engine, 
                vpe::VPEWorld* physics, 
                vecs::Registry* registry,
                glmvec3 position) {
        
        m_physics = physics;
        m_registry = registry;

        m_sceneHandle = engine.CreateScene(
            vve::Name{},
            vve::ParentHandle{},
            vve::Filename{"assets/standard/Ultracompact_Car.obj"},
            aiProcess_FlipUVs,
            vve::Position{vec3_t{position.x, position.y, position.z}},
            vve::Rotation{}
        );

        // Create physics body
        m_body = std::make_shared<vpe::VPEWorld::Body>(
            m_physics,
            "PlayerCar",
            reinterpret_cast<void*>(m_sceneHandle.GetValue()),
            &m_physics->g_cube,
            glmvec3{3.5f, 2.0f, 2.0f},  // Physics box roughly matching car dimensions
            vpe::toPhysics(position),
            glmquat{},  // Rotate physics body -90° around X-axis
            vpe::toPhysics(glmvec3{0.0f, 0.0f, 0.0f}),  // Initial velocity
            vpe::toPhysics(glmvec3{0.0f, 0.0f, 0.0f}),  // Initial angular velocity
            1.0_real / 10.0_real,  // Mass = 10 kg
            0.1_real,  // Low restitution (less bouncy)
            0.8_real   // High friction
        );

        // Set up physics callback to sync visual with physics
        m_body->m_on_move = [this, registry](double dt, std::shared_ptr<vpe::VPEWorld::Body> body) {
            auto pos = body->m_positionW;
            auto orient = body->m_orientationLW;
            body->stepPosition(dt, pos, orient, false);
            
            // Apply 180-degree rotation around Z-axis to fix car model orientation
            glmmat3 physicsRotation = vpe::fromPhysics(toMat3(orient));
            glmmat4 rotationOffset = glm::rotate(glm::mat4(1.0f), glm::radians(180.0f), glmvec3{0.0f, 0.0f, 1.0f});
            glmmat3 finalRotation = glmmat3(rotationOffset) * physicsRotation;
            
            vecs::Handle node = vecs::Handle(reinterpret_cast<size_t>(body->m_owner));
            registry->Put(node, 
                         vve::Position(vpe::fromPhysics(pos)), 
                         vve::Rotation(finalRotation));
        };

        // Initial sync
        m_body->m_on_move(0.0, m_body);
        m_physics->addBody(m_body);
    }

    /**
     * @brief Handle player input for car control (WASD)
     * @param dt Delta time
     * @param forward W key pressed
     * @param backward S key pressed
     * @param left A key pressed
     * @param right D key pressed
     */
    void HandleInput(float dt, bool forward, bool backward, bool left, bool right) {
        const float acceleration = 150.0f;   // Forward/backward acceleration
        const float turnSpeed = 2.5f;       // Rotation speed (radians/sec)
        const float maxSpeed = 2000.0f;       // Maximum speed
        const float drag = 0.95f;           // Friction/air resistance

        // Get current velocity in physics space
        glmvec3 velocity = m_body->m_linear_velocityW;
        glmvec3 velocityLocal = vpe::fromPhysics(velocity);

        // Apply drag (friction)
        m_speed *= drag;

        // Acceleration/Braking
        if (forward) {
            m_speed += acceleration * dt;  // W now accelerates
        }
        if (backward) {
            m_speed -= acceleration * dt;  // S now decelerates
        }

        // Clamp speed
        m_speed = glm::clamp(m_speed, -maxSpeed * 0.5f, maxSpeed);

        // Turning (only when moving)
        if (glm::abs(m_speed) > 0.1f) {
            if (left) {
                m_rotation += turnSpeed * dt;
            }
            if (right) {
                m_rotation -= turnSpeed * dt;
            }
        }

        // Calculate forward direction based on rotation
        glmvec3 forward_dir{
            glm::cos(m_rotation),
            glm::sin(m_rotation),
            0.0f
        };

        // Apply velocity in the forward direction
        glmvec3 newVelocity = forward_dir * m_speed;
        m_body->m_linear_velocityW = vpe::toPhysics(newVelocity);

        // Update rotation (orientation)
        glmmat4 rotMat = glm::rotate(glm::mat4{1.0f}, m_rotation, glmvec3{0.0f, 0.0f, 1.0f});
        m_body->m_orientationLW = vpe::toPhysics(rotMat);
        
        // Clamp position to stay within bounds
        glmvec3 currentPos = m_body->m_positionW;
        const float boundaryLimit = 67.0f;
        currentPos.x = glm::clamp(currentPos.x, -boundaryLimit, boundaryLimit);
        currentPos.y = glm::clamp(currentPos.y, -boundaryLimit, boundaryLimit);
        currentPos.z = glm::clamp(currentPos.z, -boundaryLimit, boundaryLimit);
        m_body->m_positionW = currentPos;
    }

    /**
     * @brief Get car's current position
     */
    glmvec3 GetPosition() const {
        return vpe::fromPhysics(m_body->m_positionW);
    }

    /**
     * @brief Get car's current rotation angle
     */
    float GetRotation() const {
        return m_rotation;
    }

    /**
     * @brief Get the physics body (for collision detection)
     */
    std::shared_ptr<vpe::VPEWorld::Body> GetBody() const {
        return m_body;
    }

private:
    vpe::VPEWorld* m_physics = nullptr;
    vecs::Registry* m_registry = nullptr;
    vecs::Handle m_sceneHandle{};
    std::shared_ptr<vpe::VPEWorld::Body> m_body;

    // Car state
    float m_speed = 0.0f;          // Current speed (forward/backward)
    float m_rotation = 0.0f;        // Current rotation angle (radians)
}; 
