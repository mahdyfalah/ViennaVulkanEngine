#pragma once

#include "VHInclude.h"
#include "VEInclude.h"
#include "VPE.hpp"

/**
 * @brief Car class for player and AI-controlled vehicles
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
     * @param isPlayer True for player car, false for AI
     * @param colorMultiplier Color tint (for differentiating cars)
     */
    void Create(vve::Engine& engine, 
                vpe::VPEWorld* physics, 
                vecs::Registry* registry,
                glmvec3 position,
                bool isPlayer = false,
                vec3_t colorMultiplier = vec3_t{1.0f, 1.0f, 1.0f}) {
        
        m_physics = physics;
        m_registry = registry;
        m_isPlayer = isPlayer;
        m_isAlive = true;
        m_colorMultiplier = colorMultiplier;

// Create visual representation using Ultracompact_Car.obj
        m_sceneHandle = engine.CreateScene(
            vve::Name{},
            vve::ParentHandle{},
            vve::Filename{"assets/standard/Ultracompact_Car.obj"},
            aiProcess_FlipUVs,  // Adjusted to fix polygon visibility
            vve::Position{vec3_t{position.x, position.y, position.z}},  // Raised position to prevent clipping
            vve::Rotation{mat3_t{glm::rotate(mat4_t{1.0f}, glm::radians(180.0f), vec3_t{0.0f, 1.0f, 0.0f})}}  // Rotate 180° to fix forward/backward direction
        );

        // Create physics body
        // Note: Body constructor expects glmquat for orientation, not mat4
        m_body = std::make_shared<vpe::VPEWorld::Body>(
            m_physics,
            isPlayer ? "PlayerCar" : "AICar" + std::to_string(reinterpret_cast<size_t>(this)),
            reinterpret_cast<void*>(m_sceneHandle.GetValue()),
            &m_physics->g_cube,
            glmvec3{3.5f, 2.0f, 2.0f},  // Physics box roughly matching sports car dimensions
            vpe::toPhysics(position),
            glmquat{glm::rotate(glm::quat(1.0f, 0.0f, 0.0f, 0.0f), glm::radians(-90.0f), glmvec3{1.0f, 0.0f, 0.0f})},  // Rotate physics body -90° around X-axis
            vpe::toPhysics(glmvec3{0.0f, 0.0f, 0.0f}),  // Initial velocity
            vpe::toPhysics(glmvec3{0.0f, 0.0f, 0.0f}),  // Initial angular velocity
            1.0_real / 10.0_real,  // Mass (10kg car)
            0.1_real,  // Low restitution (less bouncy)
            0.8_real   // High friction
        );

        // Set up physics callback to sync visual with physics
        m_body->m_on_move = [this, registry](double dt, std::shared_ptr<vpe::VPEWorld::Body> body) {
            auto pos = body->m_positionW;
            auto orient = body->m_orientationLW;
            body->stepPosition(dt, pos, orient, false);
            
            vecs::Handle node = vecs::Handle(reinterpret_cast<size_t>(body->m_owner));
            registry->Put(node, 
                         vve::Position(vpe::fromPhysics(pos)), 
                         vve::Rotation(vpe::fromPhysics(toMat3(orient))));
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
        if (!m_isPlayer || !m_isAlive) return;

        const float acceleration = 150.0f;   // Forward/backward acceleration
        const float turnSpeed = 2.5f;       // Rotation speed (radians/sec)
        const float maxSpeed = 2000.0f;       // Maximum speed
        const float drag = 0.95f;           // Friction/air resistance

        // Get current velocity in physics space
        glmvec3 velocity = m_body->m_linear_velocityW;
        glmvec3 velocityLocal = vpe::fromPhysics(velocity);

        // Apply drag (friction)
        m_speed *= drag;

        // Acceleration/Braking - swapped W and S controls
        if (forward) {
            m_speed -= acceleration * dt;  // W now brakes/decelerates
        }
        if (backward) {
            m_speed += acceleration * dt;  // S now accelerates
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
    }

    /**
     * @brief Update AI behavior (placeholder for now)
     * @param dt Delta time
     */
    void UpdateAI(float dt) {
        if (m_isPlayer || !m_isAlive) return;

        // TODO: Implement AI behaviors (wander, flee, seek)
        // For now, AI cars just sit still
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
     * @brief Check if car is still alive
     */
    bool IsAlive() const {
        return m_isAlive;
    }

    /**
     * @brief Destroy/eliminate the car
     */
    void Destroy() {
        if (!m_isAlive) return;
        
        m_isAlive = false;
        m_physics->eraseBody(m_body);  // Pass the shared_ptr directly, not the owner
        
        // TODO: Remove visual representation from scene
        // This will be handled when we implement proper entity removal
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
    bool m_isPlayer = false;
    bool m_isAlive = false;
    float m_speed = 0.0f;          // Current speed (forward/backward)
    float m_rotation = 0.0f;        // Current rotation angle (radians)
    vec3_t m_colorMultiplier{1.0f}; // Color for differentiation

    // AI state (for future use)
    glmvec3 m_targetPosition{0.0f};
    float m_wanderTimer = 0.0f;
};
