#pragma once

#include "VHInclude.h"
#include "VEInclude.h"
#include "VPE.hpp"

/**
 * @brief AmmoPack class - collectible ammo that appears randomly
 * Rotates to attract player attention
 */
class AmmoPack {
public:
    AmmoPack() = default;
    ~AmmoPack() = default;

    /**
     * @brief Create an ammo pack at specified position
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
        m_isActive = true;
        m_rotation = 0.0f;

        // Create visual representation using cube (medium size, visible)
        m_sceneHandle = engine.CreateScene(
            vve::Name{},
            vve::ParentHandle{},
            vve::Filename{"assets/test/crate0/cube.obj"},
            static_cast<aiPostProcessSteps>(aiProcess_FlipUVs | aiProcess_FlipWindingOrder),
            vve::Position{vec3_t{position.x, position.y, position.z}},
            vve::Rotation{},
            vve::Scale{vec3_t{3.0f, 3.0f, 3.0f}}  // Medium size for visibility
        );

        // Create physics body (stationary trigger)
        m_body = std::make_shared<vpe::VPEWorld::Body>(
            m_physics,
            "AmmoPack" + std::to_string(reinterpret_cast<size_t>(this)),
            reinterpret_cast<void*>(m_sceneHandle.GetValue()),
            &m_physics->g_cube,
            glmvec3{3.0f, 3.0f, 3.0f},
            vpe::toPhysics(position),
            glmquat{1.0f, 0.0f, 0.0f, 0.0f},
            vpe::toPhysics(glmvec3{0.0f, 0.0f, 0.0f}),  // No velocity (stationary)
            vpe::toPhysics(glmvec3{0.0f, 0.0f, 0.0f}),  // No angular velocity
            1.0_real / 1000000.0_real,  // Very high mass (immovable like walls)
            0.0_real,
            0.0_real
        );

        // Set up physics callback
        m_body->m_on_move = [this, registry](double dt, std::shared_ptr<vpe::VPEWorld::Body> body) {
            auto pos = body->m_positionW;
            auto orient = body->m_orientationLW;
            
            vecs::Handle node = vecs::Handle(reinterpret_cast<size_t>(body->m_owner));
            registry->Put(node, 
                         vve::Position(vpe::fromPhysics(pos)), 
                         vve::Rotation(vpe::fromPhysics(toMat3(orient))));
        };

        m_body->m_on_move(0.0, m_body);
        m_physics->addBody(m_body);
    }

    /**
     * @brief Update ammo pack (rotate for visibility)
     * @param dt Delta time
     */
    void Update(float dt) {
        if (!m_isActive) return;

        // Rotate around Z-axis for visual effect
        m_rotation += 2.0f * dt;  // Rotate at 2 radians per second
        
        glmmat4 rotMat = glm::rotate(glm::mat4{1.0f}, m_rotation, glmvec3{0.0f, 0.0f, 1.0f});
        m_body->m_orientationLW = vpe::toPhysics(rotMat);
        
        // Update visual
        m_registry->Put(m_sceneHandle, vve::Rotation{vpe::fromPhysics(toMat3(m_body->m_orientationLW))});
    }

    /**
     * @brief Get ammo pack's current position
     */
    glmvec3 GetPosition() const {
        return vpe::fromPhysics(m_body->m_positionW);
    }

    /**
     * @brief Check if ammo pack is active
     */
    bool IsActive() const {
        return m_isActive;
    }

    /**
     * @brief Destroy/collect the ammo pack
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
    float m_rotation = 0.0f;
};
