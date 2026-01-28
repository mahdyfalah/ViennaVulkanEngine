#pragma once

#include "VHInclude.h"
#include "VEInclude.h"
#include "VPE.hpp"

/**
 * @brief Simple bullet class for shooting mechanism
 * Uses a cube as visual representation
 */
class Bullet {
public:
    Bullet() = default;
    ~Bullet() = default;

    /**
     * @brief Create a bullet at specified position moving in a direction
     * @param engine Reference to the game engine
     * @param physics Reference to the physics world
     * @param registry Reference to the entity registry
     * @param position Starting position
     * @param direction Direction to shoot (should be normalized)
     * @param speed Bullet speed
     */
    void Create(vve::Engine& engine, 
                vpe::VPEWorld* physics, 
                vecs::Registry* registry,
                glmvec3 position,
                glmvec3 direction,
                float speed = 100.0f) {
        
        m_physics = physics;
        m_registry = registry;
        m_isActive = true;
        m_velocity = direction * speed;

        // Create visual representation using small cube
        m_sceneHandle = engine.CreateScene(
            vve::Name{},
            vve::ParentHandle{},
            vve::Filename{"assets/test/crate0/cube.obj"},
            aiProcess_FlipWindingOrder,
            vve::Position{vec3_t{position.x, position.y, position.z}},
            vve::Rotation{},
            vve::Scale{vec3_t{1.5f, 1.5f, 1.5f}}  // Small bullet
        );

        // Create physics body (small, fast-moving)
        m_body = std::make_shared<vpe::VPEWorld::Body>(
            m_physics,
            "Bullet" + std::to_string(reinterpret_cast<size_t>(this)),
            reinterpret_cast<void*>(m_sceneHandle.GetValue()),
            &m_physics->g_cube,
            glmvec3{1.5f, 1.5f, 1.5f},  // Small collision box
            vpe::toPhysics(position),
            glmquat{1.0f, 0.0f, 0.0f, 0.0f},
            vpe::toPhysics(m_velocity),  // Initial velocity
            vpe::toPhysics(glmvec3{0.0f, 0.0f, 0.0f}),  // No angular velocity
            1.0_real / 1.0_real,  // Light bullet (1kg)
            0.0_real,  // No bounce
            0.0_real   // No friction
        );

        // Set up physics callback
        m_body->m_on_move = [this, registry](double dt, std::shared_ptr<vpe::VPEWorld::Body> body) {
            auto pos = body->m_positionW;
            auto orient = body->m_orientationLW;
            body->stepPosition(dt, pos, orient, false);
            
            vecs::Handle node = vecs::Handle(reinterpret_cast<size_t>(body->m_owner));
            registry->Put(node, 
                         vve::Position(vpe::fromPhysics(pos)), 
                         vve::Rotation(vpe::fromPhysics(toMat3(orient))));
        };

        m_body->m_on_move(0.0, m_body);
        m_physics->addBody(m_body);
        
        m_lifetime = 5.0f;  // Bullets last 5 seconds max
    }

    /**
     * @brief Update bullet (decrease lifetime, check bounds)
     * @param dt Delta time
     * @return true if bullet should be removed
     */
    bool Update(float dt) {
        if (!m_isActive) return true;

        m_lifetime -= dt;
        
        // Remove if lifetime expired or out of bounds
        glmvec3 pos = GetPosition();
        const float boundaryLimit = 75.0f;
        
        if (m_lifetime <= 0.0f || 
            glm::abs(pos.x) > boundaryLimit || 
            glm::abs(pos.y) > boundaryLimit) {
            Destroy();
            return true;
        }
        
        return false;
    }

    /**
     * @brief Get bullet's current position
     */
    glmvec3 GetPosition() const {
        return vpe::fromPhysics(m_body->m_positionW);
    }

    /**
     * @brief Check if bullet is active
     */
    bool IsActive() const {
        return m_isActive;
    }

    /**
     * @brief Destroy the bullet
     */
    void Destroy() {
        if (!m_isActive) return;
        
        m_isActive = false;
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

    bool m_isActive = false;
    glmvec3 m_velocity{0.0f};
    float m_lifetime = 0.0f;
};
