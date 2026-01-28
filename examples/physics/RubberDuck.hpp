#pragma once

#include "VHInclude.h"
#include "VEInclude.h"
#include "VPE.hpp"

/**
 * @brief RubberDuck class for AI-controlled rubber ducks
 * Implements seek, flee, and wander behaviors
 */
class RubberDuck {
public:
    enum class AIBehavior {
        WANDER,  // Random movement
        SEEK,    // Chase target
        FLEE     // Run away from target
    };

    RubberDuck() = default;
    ~RubberDuck() = default;

    /**
     * @brief Create a rubber duck at specified position
     * @param engine Reference to the game engine
     * @param physics Reference to the physics world
     * @param registry Reference to the entity registry
     * @param position Starting position
     * @param colorMultiplier Color tint (for differentiating ducks)
     */
    void Create(vve::Engine& engine, 
                vpe::VPEWorld* physics, 
                vecs::Registry* registry,
                glmvec3 position,
                vec3_t colorMultiplier = vec3_t{1.0f, 1.0f, 1.0f}) {
        
        m_physics = physics;
        m_registry = registry;
        m_isAlive = true;
        m_colorMultiplier = colorMultiplier;

        // Create visual representation using Rubber_duck.obj
        m_sceneHandle = engine.CreateScene(
            vve::Name{},
            vve::ParentHandle{},
            vve::Filename{"assets/standard/Rubber_duck.obj"},
            aiProcess_FlipUVs,
            vve::Position{vec3_t{position.x, position.y, position.z}},
            vve::Rotation{mat3_t{glm::rotate(mat4_t{1.0f}, glm::radians(180.0f), vec3_t{0.0f, 1.0f, 0.0f})}}
        );

        // Create physics body
        m_body = std::make_shared<vpe::VPEWorld::Body>(
            m_physics,
            "RubberDuck" + std::to_string(reinterpret_cast<size_t>(this)),
            reinterpret_cast<void*>(m_sceneHandle.GetValue()),
            &m_physics->g_cube,
            glmvec3{2.5f, 2.5f, 2.5f},  // Physics box for rubber duck
            vpe::toPhysics(position),
            glmquat{glm::rotate(glm::quat(1.0f, 0.0f, 0.0f, 0.0f), glm::radians(-90.0f), glmvec3{1.0f, 0.0f, 0.0f})},
            vpe::toPhysics(glmvec3{0.0f, 0.0f, 0.0f}),  // Initial velocity
            vpe::toPhysics(glmvec3{0.0f, 0.0f, 0.0f}),  // Initial angular velocity
            1.0_real / 10.0_real,  // Mass (10kg duck)
            0.1_real,  // Low restitution
            0.8_real   // High friction
        );

        // Set up physics callback to sync visual with physics
        m_body->m_on_move = [this, registry](double dt, std::shared_ptr<vpe::VPEWorld::Body> body) {
            auto pos = body->m_positionW;
            auto orient = body->m_orientationLW;
            body->stepPosition(dt, pos, orient, false);
            
            // Apply 180-degree rotation around Z-axis to fix duck model orientation
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
     * @brief Update AI behavior (seek, flee, wander)
     * @param dt Delta time
     */
    void UpdateAI(float dt) {
        if (!m_isAlive) return;

        const float acceleration = 100.0f;
        const float maxSpeed = 50.0f;
        const float drag = 0.95f;
        const float turnSpeed = 2.0f;

        // Update behavior timer
        m_behaviorTimer -= dt;
        if (m_behaviorTimer <= 0.0f) {
            // Switch behavior randomly
            int behavior = rand() % 3;
            m_currentBehavior = static_cast<AIBehavior>(behavior);
            m_behaviorTimer = 3.0f + (rand() % 3);  // 3-6 seconds

            // Pick new wander target
            if (m_currentBehavior == AIBehavior::WANDER) {
                m_targetPosition = glmvec3{
                    (rand() % 100 - 50) * 1.2f,
                    (rand() % 100 - 50) * 1.2f,
                    0.0f
                };
            }
        }

        glmvec3 currentPos = GetPosition();
        glmvec3 targetDir{0.0f};

        // Determine target direction based on behavior
        switch (m_currentBehavior) {
            case AIBehavior::SEEK:
                // Seek player
                if (m_targetPosition != glmvec3{0.0f}) {
                    targetDir = glm::normalize(m_targetPosition - currentPos);
                }
                break;

            case AIBehavior::FLEE:
                // Flee from player
                if (m_targetPosition != glmvec3{0.0f}) {
                    targetDir = glm::normalize(currentPos - m_targetPosition);
                }
                break;

            case AIBehavior::WANDER:
                // Wander to random point
                targetDir = glm::normalize(m_targetPosition - currentPos);
                // Pick new target when close
                if (glm::length(m_targetPosition - currentPos) < 10.0f) {
                    m_targetPosition = glmvec3{
                        (rand() % 100 - 50) * 1.2f,
                        (rand() % 100 - 50) * 1.2f,
                        0.0f
                    };
                }
                break;
        }

        // Boundary avoidance (high priority)
        const float boundaryLimit = 65.0f;
        const float boundaryBuffer = 15.0f;
        glmvec3 boundaryForce{0.0f};

        if (currentPos.x > boundaryLimit - boundaryBuffer) {
            boundaryForce.x = -1.0f;
        } else if (currentPos.x < -(boundaryLimit - boundaryBuffer)) {
            boundaryForce.x = 1.0f;
        }

        if (currentPos.y > boundaryLimit - boundaryBuffer) {
            boundaryForce.y = -1.0f;
        } else if (currentPos.y < -(boundaryLimit - boundaryBuffer)) {
            boundaryForce.y = 1.0f;
        }

        // Blend boundary avoidance with target direction
        if (glm::length(boundaryForce) > 0.01f) {
            targetDir = glm::normalize(boundaryForce * 2.0f + targetDir);
        }

        // Calculate desired rotation
        float targetRotation = std::atan2(targetDir.y, targetDir.x);

        // Smoothly turn towards target
        float rotationDiff = targetRotation - m_rotation;
        // Normalize angle difference to [-PI, PI]
        while (rotationDiff > M_PI) rotationDiff -= 2.0f * M_PI;
        while (rotationDiff < -M_PI) rotationDiff += 2.0f * M_PI;

        m_rotation += glm::clamp(rotationDiff, -turnSpeed * dt, turnSpeed * dt);

        // Apply acceleration in forward direction
        m_speed += acceleration * dt;
        m_speed = glm::clamp(m_speed, 0.0f, maxSpeed);
        m_speed *= drag;

        // Calculate forward direction and apply velocity
        glmvec3 forward_dir{
            glm::cos(m_rotation),
            glm::sin(m_rotation),
            0.0f
        };

        glmvec3 newVelocity = forward_dir * m_speed;
        m_body->m_linear_velocityW = vpe::toPhysics(newVelocity);

        // Update rotation
        glmmat4 rotMat = glm::rotate(glm::mat4{1.0f}, m_rotation, glmvec3{0.0f, 0.0f, 1.0f});
        m_body->m_orientationLW = vpe::toPhysics(rotMat);

        // Hard clamp position (safety net)
        glmvec3 pos = m_body->m_positionW;
        pos.x = glm::clamp(pos.x, -boundaryLimit, boundaryLimit);
        pos.y = glm::clamp(pos.y, -boundaryLimit, boundaryLimit);
        m_body->m_positionW = pos;
    }

    /**
     * @brief Set target position for AI (used for seek/flee)
     */
    void SetAITarget(glmvec3 target) {
        m_targetPosition = target;
    }

    /**
     * @brief Get duck's current position
     */
    glmvec3 GetPosition() const {
        return vpe::fromPhysics(m_body->m_positionW);
    }

    /**
     * @brief Get duck's current rotation angle
     */
    float GetRotation() const {
        return m_rotation;
    }

    /**
     * @brief Check if duck is still alive
     */
    bool IsAlive() const {
        return m_isAlive;
    }

    /**
     * @brief Destroy/eliminate the duck
     */
    void Destroy() {
        if (!m_isAlive) return;
        
        m_isAlive = false;
        m_physics->eraseBody(m_body);
        
        // Hide visual representation by scaling to 0
        m_registry->Put(m_sceneHandle, vve::Scale{vec3_t{0.0f, 0.0f, 0.0f}});
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

    // Duck state
    bool m_isAlive = false;
    float m_speed = 0.0f;
    float m_rotation = 0.0f;
    vec3_t m_colorMultiplier{1.0f};

    // AI state
    AIBehavior m_currentBehavior = AIBehavior::WANDER;
    glmvec3 m_targetPosition{0.0f};
    float m_wanderTimer = 0.0f;
    float m_behaviorTimer = 0.0f;
};
