#pragma once

#include "VHInclude.h"
#include "VEInclude.h"
#include "VPE.hpp"

/**
 * @brief GameMap class for managing the playing field boundaries
 * Creates immovable walls to define the playable area
 */
class GameMap {
public:
    GameMap() = default;
    ~GameMap() = default;

    /**
     * @brief Initialize the game map with walls
     * @param engine Reference to the game engine
     * @param physics Reference to the physics world
     * @param registry Reference to the entity registry
     * @param playFieldSize Size of the playing field (half-extents)
     */
    void Create(vve::Engine& engine, 
                vpe::VPEWorld* physics, 
                vecs::Registry* registry,
                float playFieldSize = 50.0f) {
        
        m_physics = physics;
        m_registry = registry;
        m_playFieldSize = playFieldSize;

        // Create north wall (positive Y direction)
        CreateWall(engine, glmvec3{0.0f, playFieldSize, 0.0f}, glmvec3{playFieldSize * 2.0f, 5.0f, 5.0f});
        
        // TODO: Create other three walls (south, east, west)
    }

    /**
     * @brief Create a single wall segment
     * @param engine Reference to the game engine
     * @param position Center position of the wall
     * @param scale Scale of the wall (width, height, thickness)
     */
    void CreateWall(vve::Engine& engine, glmvec3 position, glmvec3 scale) {
        // Create visual representation using cube
        vecs::Handle sceneHandle = engine.CreateScene(
            vve::Name{},
            vve::ParentHandle{},
            vve::Filename{"assets/test/crate0/cube.obj"},
            aiProcess_FlipUVs,
            vve::Position{vec3_t{position.x, position.y, position.z}},
            vve::Rotation{},
            vve::Scale{vec3_t{scale.x, scale.y, scale.z}}
        );

        // Create immovable physics body (infinite mass = static object)
        auto body = std::make_shared<vpe::VPEWorld::Body>(
            m_physics,
            "Wall_" + std::to_string(m_walls.size()),
            reinterpret_cast<void*>(sceneHandle.GetValue()),
            &m_physics->g_cube,
            scale,  // Physics box size matches visual scale
            vpe::toPhysics(position),
            glmquat{1.0f, 0.0f, 0.0f, 0.0f},  // No rotation
            vpe::toPhysics(glmvec3{0.0f, 0.0f, 0.0f}),  // Zero initial velocity
            vpe::toPhysics(glmvec3{0.0f, 0.0f, 0.0f}),  // Zero angular velocity
            0.0_real,  // Zero mass = infinite mass = immovable static object
            0.0_real,  // No restitution (no bounce)
            1.0_real   // High friction
        );

        // Set up physics callback (walls don't move, but we keep the structure)
        body->m_on_move = [this](double dt, std::shared_ptr<vpe::VPEWorld::Body> body) {
            // Static walls don't move, but we keep the callback for consistency
            auto pos = body->m_positionW;
            auto orient = body->m_orientationLW;
            
            vecs::Handle node = vecs::Handle(reinterpret_cast<size_t>(body->m_owner));
            m_registry->Put(node, 
                           vve::Position(vpe::fromPhysics(pos)), 
                           vve::Rotation(vpe::fromPhysics(toMat3(orient))));
        };

        body->m_on_move(0.0, body);
        m_physics->addBody(body);
        
        m_walls.push_back(body);
    }

private:
    vpe::VPEWorld* m_physics = nullptr;
    vecs::Registry* m_registry = nullptr;
    float m_playFieldSize = 50.0f;
    std::vector<std::shared_ptr<vpe::VPEWorld::Body>> m_walls;
};
