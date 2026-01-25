#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <iostream>
#include <utility>
#include <format>
#include "VHInclude.h"
#include "VEInclude.h"

#include "VPE.hpp"
#include "Car.hpp"


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
            pn() = vec3_t{0.0f, -20.0f, 30.0f};  // High above the play area
            
            // Rotate camera to look down at 45° angle (isometric view)
            m_registry.Get<vve::Rotation&>(m_cameraHandle)() = mat3_t{ 
                glm::rotate(mat4_t{1.0f}, glm::radians(45.0f), vec3_t{1.0f, 0.0f, 0.0f}) 
            };

            // ----------------- Create Player Car -----------------
            
            m_playerCar.Create(m_engine, &m_physics, &m_registry, 
                              glmvec3{0.0f, 2.0f, 0.0f},  // Center of play area, 2 units up
                              true,  // isPlayer
                              vec3_t{0.0f, 1.0f, 0.0f}); // Green color

            // Debug: Print positions
            std::cout << "\n=== SPAWN DEBUG ==="<< std::endl;
            std::cout << "Camera Node Position: (" << pn().x << ", " << pn().y << ", " << pn().z << ")" << std::endl;
            std::cout << "Car Spawn Request: (0, 2, 0)" << std::endl;
            auto carPos = m_playerCar.GetPosition();
            std::cout << "Car Actual Position: (" << carPos.x << ", " << carPos.y << ", " << carPos.z << ")" << std::endl;
            std::cout << "Distance from camera: " << glm::length(vec3_t{pn().x - carPos.x, pn().y - carPos.y, pn().z - carPos.z}) << std::endl;
            std::cout << "==================\n" << std::endl;

            // TODO: Spawn 3 AI cars
            // Will be implemented later with AI behaviors

			m_engine.SetVolume(m_volume);
            return false;
        };
    
        bool OnUpdate( Message& message ) {
            auto msg = message.template GetData<vve::System::MsgUpdate>();
            auto dt = msg.m_dt;
            
            // Update player car with input
            m_playerCar.HandleInput((float)dt, m_keyW, m_keyS, m_keyA, m_keyD);
            
            // Update AI cars
            for (auto& aiCar : m_aiCars) {
                aiCar.UpdateAI((float)dt);
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

            // TODO: Space bar for shooting

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
            
            // Display AI car count
            int aliveAI = 0;
            for (const auto& ai : m_aiCars) {
                if (ai.IsAlive()) aliveAI++;
            }
            ImGui::Text("AI Cars Remaining: %d / %d", aliveAI, (int)m_aiCars.size());
            
            ImGui::Separator();
            ImGui::Text("Controls:");
            ImGui::Text("  WASD - Drive");
            ImGui::Text("  Space - Shoot (TODO)");
            
            ImGui::Separator();
        	if (ImGui::SliderFloat("Sound Volume", &m_volume, 0, MIX_MAX_VOLUME)) {
		    	m_engine.SetVolume(m_volume);
			}
            ImGui::End();
            return false;
        }

    private:
    	vpe::VPEWorld m_physics;
        Car m_playerCar;
        std::vector<Car> m_aiCars;
        
        // Input state
        bool m_keyW = false;
        bool m_keyA = false;
        bool m_keyS = false;
        bool m_keyD = false;
        
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
    
    