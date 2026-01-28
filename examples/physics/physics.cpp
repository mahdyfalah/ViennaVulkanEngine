#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <iostream>
#include <utility>
#include <format>
#include <chrono>
#include "VHInclude.h"
#include "VEInclude.h"

#include "VPE.hpp"
#include "Car.hpp"
#include "RubberDuck.hpp"
#include "Bullet.hpp"
#include "AmmoPack.hpp"
#include "GameMap.hpp"


class MyGame : public vve::System {

	std::default_random_engine rnd_gen{ 12345 };					//Random numbers
	std::uniform_real_distribution<> rnd_unif{ 0.0f, 1.0f };		//Random numbers

    public:
        MyGame( vve::Engine& engine ) : vve::System("MyGame", engine ) {
            m_engine.RegisterCallbacks( { 
                {this,      0, "LOAD_LEVEL", [this](Message& message){ return OnLoadLevel(message);} },
                {this,  10000, "UPDATE", [this](Message& message){ return OnUpdate(message);} },
			    {this,      0, "SDL_KEY_DOWN", [this](Message& message){ return OnKeyDown(message);} },
			    {this,      0, "SDL_KEY_UP", [this](Message& message){ return OnKeyUp(message);} },
                {this, -10000, "RECORD_NEXT_FRAME", [this](Message& message){ return OnRecordNextFrame(message);} }
            } );
            m_engine.SetVolume(m_volume);
        };
        
        ~MyGame() {};

        void GetCamera() {
            if(m_cameraHandle.IsValid() == false) { 
                auto [handle, camera, parent] = *m_registry.GetView<vecs::Handle, vve::Camera&, vve::ParentHandle>().begin(); 
                m_cameraHandle = handle;
                m_cameraNodeHandle = parent;
            };
        }

        inline static std::string plane_obj  { "assets/test/plane/plane_t_n_s.obj" };
        inline static std::string plane_mesh { "assets/test/plane/plane_t_n_s.obj/plane" };
        inline static std::string plane_txt  { "assets/test/plane/grass.jpg" };

        bool OnLoadLevel( Message message ) {
            auto msg = message.template GetData<vve::System::MsgLoadLevel>();	
            std::cout << "Loading level: " << msg.m_level << std::endl;

            // ----------------- Load Plane -----------------

			m_engine.LoadScene( vve::Filename{plane_obj}, aiProcess_FlipWindingOrder);

			m_engine.CreateObject(	vve::Name{},
                                    vve::ParentHandle{}, 
                                    vve::MeshName{plane_mesh}, 
									vve::TextureName{plane_txt}, 
									vve::Position{vec3_t{0.0f, 0.0f, 0.0f}}, 
									vve::Rotation{mat4_t{glm::rotate(glm::mat4(1.0f), 3.14152f / 2.0f, glm::vec3(1.0f,0.0f,0.0f))}}, 
									vve::Scale{vec3_t{1000.0f, 1000.0f, 1000.0f}}, 
									vve::UVScale{vec2_t{1000.0f, 1000.0f}});

            // ----------------- Setup Camera -----------------

            GetCamera();
            
            // Set up isometric camera - looking down at 45° angle from above
            // Position camera high and back from the origin
            auto [pn, rn] = m_registry.template Get<vve::Position&, vve::Rotation&>(m_cameraNodeHandle);
            pn() = vec3_t{0.0f, -60.0f, 120.0f};  // High above the play area
            
            // Rotate camera to look down at 45° angle (isometric view)
            m_registry.Get<vve::Rotation&>(m_cameraHandle)() = mat3_t{ 
                glm::rotate(mat4_t{1.0f}, glm::radians(45.0f), vec3_t{1.0f, 0.0f, 0.0f}) 
            };

            // ----------------- Create Player Car -----------------
            
            m_playerCar.Create(m_engine, &m_physics, &m_registry, 
                              glmvec3{0.0f, 0.0f, 0.0f});  // Center of play area

            // Debug: Print positions
            std::cout << "\n=== SPAWN DEBUG ==="<< std::endl;
            std::cout << "Camera Node Position: (" << pn().x << ", " << pn().y << ", " << pn().z << ")" << std::endl;
            std::cout << "Car Spawn Request: (0, 0, 0)" << std::endl;
            auto carPos = m_playerCar.GetPosition();
            std::cout << "Car Actual Position: (" << carPos.x << ", " << carPos.y << ", " << carPos.z << ")" << std::endl;
            std::cout << "Distance from camera: " << glm::length(vec3_t{pn().x - carPos.x, pn().y - carPos.y, pn().z - carPos.z}) << std::endl;
            std::cout << "==================\n" << std::endl;

            // ----------------- Create Game Map with Walls -----------------
            
            m_gameMap.Create(m_engine, &m_physics, &m_registry, 150.0f);
            std::cout << "Created north wall boundary" << std::endl;

            // ----------------- Spawn 3 AI Rubber Ducks with Different Colors -----------------
            
            m_rubberDucks.resize(3);
            
            // Rubber Duck 1 - Red (top-right)
            m_rubberDucks[0].Create(m_engine, &m_physics, &m_registry,
                                   glmvec3{40.0f, 40.0f, 0.0f},
                                   vec3_t{1.0f, 0.2f, 0.2f});  // Red tint
            std::cout << "Spawned Rubber Duck 1 (Red) at (40, 40, 0)" << std::endl;
            
            // Rubber Duck 2 - Blue (bottom-left)
            m_rubberDucks[1].Create(m_engine, &m_physics, &m_registry,
                                   glmvec3{-40.0f, -40.0f, 0.0f},
                                   vec3_t{0.2f, 0.4f, 1.0f});  // Blue tint
            std::cout << "Spawned Rubber Duck 2 (Blue) at (-40, -40, 0)" << std::endl;
            
            // Rubber Duck 3 - Yellow (top-left)
            m_rubberDucks[2].Create(m_engine, &m_physics, &m_registry,
                                   glmvec3{-40.0f, 40.0f, 0.0f},
                                   vec3_t{1.0f, 1.0f, 0.2f});  // Yellow tint
            std::cout << "Spawned Rubber Duck 3 (Yellow) at (-40, 40, 0)" << std::endl;

			m_engine.SetVolume(m_volume);
            return false;
        };
    
        bool OnUpdate( Message& message ) {
            auto msg = message.template GetData<vve::System::MsgUpdate>();
            auto dt = msg.m_dt;
            
            // Update player car with input
            m_playerCar.HandleInput((float)dt, m_keyW, m_keyS, m_keyA, m_keyD);
            
            // Get player position for AI targeting and ammo collection
            glmvec3 playerPos = m_playerCar.GetPosition();
            
            // Update rubber ducks with player position as target
            for (auto& duck : m_rubberDucks) {
                duck.SetAITarget(playerPos);
                duck.UpdateAI((float)dt);
            }
            
            // Update and spawn ammo packs (limit to 3 on map)
            m_ammoSpawnTimer -= dt;
            if (m_ammoSpawnTimer <= 0.0f && m_ammoPacks.size() < 3) {
                SpawnAmmoPack();
                m_ammoSpawnTimer = 3.0f;  // Spawn every 3 seconds
            }
            
            // Update ammo packs (rotation) and check collection
            for (auto it = m_ammoPacks.begin(); it != m_ammoPacks.end(); ) {
                it->Update((float)dt);
                
                // Check if player collected the ammo
                glmvec3 ammoPos = it->GetPosition();
                float distance = glm::length(playerPos - ammoPos);
                
                if (distance < 6.0f) {  // Collection radius
                    m_ammo++;
                    it->Destroy();
                    std::cout << "Ammo collected! Total ammo: " << m_ammo << std::endl;
                    it = m_ammoPacks.erase(it);
                } else {
                    ++it;
                }
            }
            
            // Update bullets and check for collisions
            for (auto it = m_bullets.begin(); it != m_bullets.end(); ) {
                bool shouldRemove = it->Update((float)dt);
                
                if (!shouldRemove) {
                    // Check collision with rubber ducks
                    glmvec3 bulletPos = it->GetPosition();
                    
                    for (auto& duck : m_rubberDucks) {
                        if (duck.IsAlive()) {
                            glmvec3 duckPos = duck.GetPosition();
                            float distance = glm::length(bulletPos - duckPos);
                            
                            // Hit detection (within 5 units)
                            if (distance < 5.0f) {
                                duck.Destroy();
                                it->Destroy();
                                shouldRemove = true;
                                std::cout << "Duck hit! Distance: " << distance << std::endl;
                                break;
                            }
                        }
                    }
                }
                
                if (shouldRemove) {
                    it = m_bullets.erase(it);
                } else {
                    ++it;
                }
            }
            
            m_physics.tick(dt);
            return false;
        }

        bool OnKeyDown(Message message) {
			auto msg = message.template GetData<MsgKeyDown>();
			auto key = msg.m_key;

            // Player car controls - WASD
            if (key == SDL_SCANCODE_W) m_keyW = true;
            if (key == SDL_SCANCODE_A) m_keyA = true;
            if (key == SDL_SCANCODE_S) m_keyS = true;
            if (key == SDL_SCANCODE_D) m_keyD = true;

            // Space bar for shooting
            if (key == SDL_SCANCODE_SPACE) {
                // Check if player has ammo
                if (m_ammo <= 0) {
                    std::cout << "No ammo! Collect ammo packs first." << std::endl;
                    return false;
                }
                
                // Check cooldown
                auto currentTime = std::chrono::steady_clock::now();
                auto timeSinceLastShot = std::chrono::duration_cast<std::chrono::milliseconds>(
                    currentTime - m_lastShotTime).count();
                
                if (timeSinceLastShot >= m_shootCooldown) {
                    ShootBullet();
                    m_ammo--;  // Consume ammo
                    m_lastShotTime = currentTime;
                    std::cout << "Ammo remaining: " << m_ammo << std::endl;
                }
            }

		    return false;
        }
        
        bool OnKeyUp(Message message) {
			auto msg = message.template GetData<MsgKeyUp>();
			auto key = msg.m_key;

            // Release player car controls
            if (key == SDL_SCANCODE_W) m_keyW = false;
            if (key == SDL_SCANCODE_A) m_keyA = false;
            if (key == SDL_SCANCODE_S) m_keyS = false;
            if (key == SDL_SCANCODE_D) m_keyD = false;

		    return false;
        }
    
        bool OnRecordNextFrame(Message message) { 

            ImGui::Begin("Game State");
            
            // Display player car info
            auto playerPos = m_playerCar.GetPosition();
            ImGui::Text("Player Position: (%.1f, %.1f, %.1f)", playerPos.x, playerPos.y, playerPos.z);
            ImGui::Text("Player Rotation: %.2f rad", m_playerCar.GetRotation());
            
            // Display rubber duck count
            int aliveDucks = 0;
            for (const auto& duck : m_rubberDucks) {
                if (duck.IsAlive()) aliveDucks++;
            }
            ImGui::Text("Rubber Ducks Remaining: %d / %d", aliveDucks, (int)m_rubberDucks.size());
            ImGui::Text("Active Bullets: %d", (int)m_bullets.size());
            ImGui::Text("Ammo: %d", m_ammo);
            ImGui::Text("Active Ammo Packs: %d", (int)m_ammoPacks.size());
            
            ImGui::Separator();
            ImGui::Text("Controls:");
            ImGui::Text("  WASD - Drive");
            ImGui::Text("  Space - Shoot");
            
            ImGui::Separator();
        	if (ImGui::SliderFloat("Sound Volume", &m_volume, 0, MIX_MAX_VOLUME)) {
		    	m_engine.SetVolume(m_volume);
			}
            ImGui::End();
            return false;
        }

        void ShootBullet() {
            // Get player position and rotation
            glmvec3 carPos = m_playerCar.GetPosition();
            float carRotation = m_playerCar.GetRotation();
            
            // Calculate shooting direction based on car rotation
            glmvec3 shootDirection{
                glm::cos(carRotation),
                glm::sin(carRotation),
                0.0f
            };
            
            // Spawn bullet slightly in front of car
            glmvec3 bulletStartPos = carPos + shootDirection * 5.0f;
            
            // Create bullet
            Bullet newBullet;
            newBullet.Create(m_engine, &m_physics, &m_registry, 
                            bulletStartPos, shootDirection, 120.0f);
            m_bullets.push_back(newBullet);
            
            std::cout << "Bullet fired! Position: (" << bulletStartPos.x << ", " 
                      << bulletStartPos.y << "), Direction: (" << shootDirection.x 
                      << ", " << shootDirection.y << ")" << std::endl;
        }

        void SpawnAmmoPack() {
            // Get player position
            glmvec3 carPos = m_playerCar.GetPosition();
            
            // Spawn ammo pack at random position near player (within 30-50 units)
            float angle = (rand() % 360) * (M_PI / 180.0f);
            float distance = 30.0f + (rand() % 20);  // 30-50 units away
            
            glmvec3 spawnPos{
                carPos.x + glm::cos(angle) * distance,
                carPos.y + glm::sin(angle) * distance,
                0.0f
            };
            
            // Clamp to playfield bounds
            const float boundaryLimit = 65.0f;
            spawnPos.x = glm::clamp(spawnPos.x, -boundaryLimit, boundaryLimit);
            spawnPos.y = glm::clamp(spawnPos.y, -boundaryLimit, boundaryLimit);
            
            // Create ammo pack
            AmmoPack newAmmo;
            newAmmo.Create(m_engine, &m_physics, &m_registry, spawnPos);
            m_ammoPacks.push_back(newAmmo);
            
            std::cout << "Ammo pack spawned at (" << spawnPos.x << ", " << spawnPos.y << ")" << std::endl;
        }

    private:
    	vpe::VPEWorld m_physics;
        GameMap m_gameMap;
        Car m_playerCar;
        std::vector<RubberDuck> m_rubberDucks;
        std::vector<Bullet> m_bullets;
        std::vector<AmmoPack> m_ammoPacks;
        
        // Input state
        bool m_keyW = false;
        bool m_keyA = false;
        bool m_keyS = false;
        bool m_keyD = false;
        
        // Shooting cooldown
        std::chrono::steady_clock::time_point m_lastShotTime{};
        int m_shootCooldown = 300;  // milliseconds between shots
        
        // Ammo system
        int m_ammo = 0;
        float m_ammoSpawnTimer = 3.0f;  // Start spawning after 3 seconds
        
        vecs::Handle m_handlePlane{};
        vecs::Handle m_cameraHandle{};
        vecs::Handle m_cameraNodeHandle{};
        float m_volume{MIX_MAX_VOLUME / 2.0};
    };
    
    
    int main() {
        vve::Engine engine("My Engine", vve::RendererType::RENDERER_TYPE_FORWARD);
        MyGame mygame{engine};  
        engine.Run();
    
        return 0;
    }
    
    