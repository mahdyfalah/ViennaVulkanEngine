#include <iostream>
#include <utility>
#include <format>
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include "VHInclude.h"
#include "VEInclude.h"

#include "VPE.hpp"

class MyGame : public vve::System
{

    std::default_random_engine rnd_gen{12345};             // Random numbers
    std::uniform_real_distribution<> rnd_unif{0.0f, 1.0f}; // Random numbers
    glmmat3 C = glmmat3{{1, 0, 0}, {0, 0, -1}, {0, 1, 0}};
    glmmat3 CTrans = glm::transpose(C);

public:
    MyGame(vve::Engine &engine) : vve::System("MyGame", engine)
    {
        m_static_registry = &m_registry;

        m_engine.RegisterCallbacks({{this, 0, "LOAD_LEVEL", [this](Message &message)
                                     { return OnLoadLevel(message); }},
                                    {this, 10000, "UPDATE", [this](Message &message)
                                     { return OnUpdate(message); }},
                                    {this, 0, "SDL_KEY_DOWN", [this](Message &message)
                                     { return OnKeyDown(message); }},
                                    {this, -10000, "RECORD_NEXT_FRAME", [this](Message &message)
                                     { return OnRecordNextFrame(message); }}});
        m_engine.SetVolume(m_volume);
    };

    ~MyGame() {};

    void GetCamera()
    {
        if (m_cameraHandle.IsValid() == false)
        {
            auto [handle, camera, parent] = *m_registry.GetView<vecs::Handle, vve::Camera &, vve::ParentHandle>().begin();
            m_cameraHandle = handle;
            m_cameraNodeHandle = parent;
        };
    }

    glmvec3 toPhysics(glmvec3 vec) { return C * vec; }
    glmmat3 toPhysics(glmmat3 mat) { return CTrans * mat * C; }
    glmmat4 toPhysics(glmmat4 mat) { return glmmat4{CTrans} * mat * glmmat4{C}; }
    glmvec3 fromPhysics(glmvec3 vec) { return CTrans * vec; }
    glmmat3 fromPhysics(glmmat3 mat) { return C * mat * CTrans; }
    glmmat4 fromPhysics(glmmat4 mat) { return glmmat4{C} * mat * glmmat4{CTrans}; }

    inline static vecs::Registry *m_static_registry{};
    vpe::VPEWorld::callback_move onMove = [&](double dt, std::shared_ptr<vpe::VPEWorld::Body> body)
    {
        auto pos = body->m_positionW;               // New position of the scene node
        auto orient = body->m_orientationLW;        // New orientation of the scende node
        body->stepPosition(dt, pos, orient, false); // Extrapolate
        auto model = vpe::VPEWorld::Body::computeModel(pos, orient, body->m_scale);
        vecs::Handle node = vecs::Handle(reinterpret_cast<size_t>(body->m_owner)); // Owner is a handle to a scene node
        m_registry.Put(node, vve::Position(fromPhysics(pos)), vve::Rotation(fromPhysics(toMat3(orient))));
    };

    inline static vpe::VPEWorld::callback_erase onErase = [&](std::shared_ptr<vpe::VPEWorld::Body> body)
    {
        auto node = vecs::Handle(reinterpret_cast<size_t>(body->m_owner)); // Owner is a pointer to a scene node
                                                                           // getSceneManagerPointer()->deleteSceneNodeAndChildren(((VESceneNode*)body->m_owner)->getName());
        return;
    };

    inline static std::string plane_obj{"assets/test/plane/plane_t_n_s.obj"};
    inline static std::string plane_mesh{"assets/test/plane/plane_t_n_s.obj/plane"};
    inline static std::string plane_txt{"assets/test/plane/grass.jpg"};

    inline static std::string cube_obj{"assets/test/crate0/cube.obj"};

    bool OnLoadLevel(Message message)
    {
        auto msg = message.template GetData<vve::System::MsgLoadLevel>();
        std::cout << "Loading level: " << msg.m_level << std::endl;
        std::string level = std::string("Level: ") + msg.m_level;

        // ----------------- Load Plane -----------------

        m_engine.LoadScene(vve::Filename{plane_obj}, aiProcess_FlipWindingOrder);

        m_engine.CreateObject(vve::Name{},
                              vve::ParentHandle{},
                              vve::MeshName{plane_mesh},
                              vve::TextureName{plane_txt},
                              vve::Position{vec3_t{0.0f, 0.0f, 0.0f}},
                              vve::Rotation{mat4_t{glm::rotate(glm::mat4(1.0f), 3.14152f / 2.0f, glm::vec3(1.0f, 0.0f, 0.0f))}},
                              vve::Scale{vec3_t{1000.0f, 1000.0f, 1000.0f}},
                              vve::UVScale{vec2_t{1000.0f, 1000.0f}});

        // ----------------- Load Cube -----------------

        // m_handleCube = m_engine.CreateScene(vve::Name{},
        //                             vve::ParentHandle{},
        //                             vve::Filename{cube_obj}, aiProcess_FlipWindingOrder,
        //							  vve::Position{{nextRandom(), nextRandom(), 0.5f}},
        //                             vve::Rotation{mat3_t{1.0f}},
        //                             vve::Scale{vec3_t{1.0f}});

        GetCamera();
        m_registry.Get<vve::Rotation &>(m_cameraHandle)() = mat3_t{glm::rotate(mat4_t{1.0f}, 3.14152f / 2.0f, vec3_t{1.0f, 0.0f, 0.0f})};

        // m_engine.PlaySound(vve::Filename{"assets/sounds/dance.mp3"}, -1, 50);
        m_engine.SetVolume(m_volume);
        return false;
    };

    bool OnUpdate(Message &message)
    {
        auto msg = message.template GetData<vve::System::MsgUpdate>();
        auto dt = msg.m_dt;
        m_physics.tick(dt);

        // --- Analytic fall update: h(t) = h0 - 0.5 * a * t^2 ---
        if (m_analyticActive && m_analyticCube.IsValid())
        {
            m_analyticT += static_cast<float>(dt);
            float h = m_dropHeight + 0.5f * m_physics.c_gravity * m_analyticT * m_analyticT;
            std::cout << std::format("t = {:.3f} s, h = {:.3f} m\n", m_analyticT, h) << std::flush;
            if (h <= 0.0f)
            {
                h = 0.0f;
                m_analyticActive = false;
                std::cout << std::format("Impact at t = {:.3f} s (a = {:.2f} m/s^2)\n",
                                         m_analyticT, m_analyticA)  << std::flush;
            }
            m_registry.Get<vve::Position &>(m_analyticCube)() = h * m_upGame;
        }

        // --- Euler fall update (explicit Euler, fixed step m_eulerDt) ---
        if (m_eulerActive && m_eulerCube.IsValid())
        {
            m_eulerAccum += dt;
            while (m_eulerAccum + 1e-9 >= m_eulerDt && m_eulerActive)
            {
                m_eulerAccum -= m_eulerDt;

                // Explicit Euler:
                // r_{n+1} = r_n + dt * v_n
                // v_{n+1} = v_n + dt * a   (a = -g along up-axis)
                m_eulerH += m_eulerV * m_eulerDt;
                m_eulerV += m_physics.c_gravity * m_eulerDt;
                m_eulerT += m_eulerDt;

                if (m_eulerH <= 0.0f)
                {
                    m_eulerH = 0.0f;
                    m_eulerV = 0.0f;
                    m_eulerActive = false;
                    std::cout << std::format("Euler impact at t = {:.3f} s (dt = {:.3f} s)\n",
                                             m_eulerT, m_eulerDt) << std::flush;
                }
                else
                {
                    std::cout << std::format("t = {:.3f} s, h = {:.3f} m\n",
                                             m_eulerT, m_eulerH) << std::flush;
                }
            }
            m_registry.Get<vve::Position &>(m_eulerCube)() = m_eulerH * m_upGame;
        }

        if (m_sympActive && m_sympCube.IsValid())
        {
            m_sympAccum += dt;
            while (m_sympAccum + 1e-9 >= m_sympDt && m_sympActive)
            {
                m_sympAccum -= m_sympDt;

                // Symplectic Euler:
                // v_{n+1} = v_n + dt * a
                // r_{n+1} = r_n + dt * v_{n+1}
                m_sympV +=  m_physics.c_gravity * m_sympDt;
                m_sympH +=  m_sympV * m_sympDt;
                m_sympT +=  m_sympDt;

                if (m_sympH <= 0.0f)
                {
                    m_sympH = 0.0f;
                    m_sympV = 0.0f;
                    m_sympActive = false;
                    std::cout << std::format("Symplectic Euler impact at t = {:.3f} s (dt = {:.3f} s)\n",
                                             m_sympT, m_sympDt) << std::flush;
                }
                else
                {
                    std::cout << std::format("t = {:.3f} s, h = {:.3f} m\n",
                                             m_sympT, m_sympH) << std::flush;
                }
            }
            m_registry.Get<vve::Position &>(m_sympCube)() = m_sympH * m_upGame;
        }

        return false;
    }

    bool OnKeyDown(Message message)
    {
        auto msg = message.template GetData<MsgKeyDown>();
        auto key = msg.m_key;

        if (key == SDL_SCANCODE_B)
        {
            DropCubeAnalytic();
        }
        if (key == SDL_SCANCODE_N)
        {
            DropCubeEuler(1.0f);
        }
        if (key == SDL_SCANCODE_M)
        {
            DropCubeEuler(0.1f);
        }
        if (key == SDL_SCANCODE_K)
        {
            DropCubeSymplectic(1.0f);
        }
        if (key == SDL_SCANCODE_L)
        {
            DropCubeSymplectic(0.1f);
        }

        return false;
    }

    bool OnRecordNextFrame(Message message)
    {

        ImGui::Begin("Game State");
        char buffer[100];
        // std::snprintf(buffer, 100, "Time Left: %.2f s", m_time_left);
        // ImGui::TextUnformatted(buffer);
        // std::snprintf(buffer, 100, "Cubes Left: %d", m_cubes_left);
        // ImGui::TextUnformatted(buffer);
        if (ImGui::SliderFloat("Sound Volume", &m_volume, 0, MIX_MAX_VOLUME))
        {
            m_engine.SetVolume(m_volume);
        }
        ImGui::End();
        return false;
    }

private:
    vpe::VPEWorld m_physics;

    vecs::Handle m_handlePlane{};
    vecs::Handle m_handleCube{};
    vecs::Handle m_cameraHandle{};
    vecs::Handle m_cameraNodeHandle{};
    float m_volume{MIX_MAX_VOLUME / 2.0};

    // --- analytic fall state ---
    vecs::Handle m_analyticCube{};
    bool  m_analyticActive{false};
    float m_analyticT{0.0f};
    float m_analyticA{9.81f};    // m/s^2
    float m_dropHeight{100.0f};  // meters
    glmvec3 m_upGame{0.0f, 1.0f, 0.0f};

    // --- Euler fall state (explicit Euler, fixed step delta_t) ---
    vecs::Handle m_eulerCube{};
    bool  m_eulerActive{false};
    float m_eulerT{0.0f};
    float m_eulerH{0.0f};
    float m_eulerV{0.0f};
    float m_eulerDt{1.0f};
    double m_eulerAccum{0.0};

    // --- Symplectic Euler state (semi-implicit Euler, fixed step delta_t) ---
    vecs::Handle m_sympCube{};
    bool  m_sympActive{false};
    float m_sympT{0.0f};
    float m_sympH{0.0f};
    float m_sympV{0.0f};
    float m_sympDt{1.0f};
    double m_sympAccum{0.0};

    void DropCubeAnalytic()
    {
        // Compute "up" from the physics gravity so the cube drops toward the plane
        const glmvec3 gravityPhys{0.0f, m_physics.c_gravity, 0.0f};
        const glmvec3 gravityGame = fromPhysics(gravityPhys);
        m_upGame = -glm::normalize(gravityGame); // opposite of gravity

        // Start height and spawn position
        m_dropHeight = 100.0f;
        m_analyticA = 9.81f;
        m_analyticT = 0.0f;
        const glmvec3 spawnPos = m_dropHeight * m_upGame;

        // Create or reuse a render-only cube (no physics body)
        if (!m_analyticCube.IsValid())
        {
            m_analyticCube = m_engine.CreateScene(vve::Name{},
                                                  vve::ParentHandle{},
                                                  vve::Filename{cube_obj}, aiProcess_FlipWindingOrder,
                                                  vve::Position{spawnPos},
                                                  vve::Rotation{mat3_t{1.0f}},
                                                  vve::Scale{vec3_t{1.0f}});
        }
        else
        {
            m_registry.Get<vve::Position &>(m_analyticCube)() = spawnPos;
            m_registry.Get<vve::Rotation &>(m_analyticCube)() = mat3_t{1.0f};
            m_registry.Get<vve::Scale &>(m_analyticCube)()    = vec3_t{1.0f};
        }

        // Begin analytic fall
        m_analyticActive = true;
    }

    void DropCubeEuler(float delta_t)
    {
        // Compute "up" exactly like DropCubeAnalytic does
        const glmvec3 gravityPhys{0.0f, m_physics.c_gravity, 0.0f};
        const glmvec3 gravityGame = fromPhysics(gravityPhys);
        m_upGame = -glm::normalize(gravityGame);

        m_dropHeight = 100.0f;
        m_eulerDt = std::max(1e-6f, delta_t);
        m_eulerT = 0.0f;
        m_eulerH = m_dropHeight;
        m_eulerV = 0.0f;
        m_eulerAccum = 0.0;
        const glmvec3 spawnPos = m_eulerH * m_upGame;

        // Create or reuse a render-only cube (separate from analytic one)
        if (!m_eulerCube.IsValid())
        {
            m_eulerCube = m_engine.CreateScene(vve::Name{},
                                               vve::ParentHandle{},
                                               vve::Filename{cube_obj}, aiProcess_FlipWindingOrder,
                                               vve::Position{spawnPos},
                                               vve::Rotation{mat3_t{1.0f}},
                                               vve::Scale{vec3_t{1.0f}});
        }
        else
        {
            m_registry.Get<vve::Position &>(m_eulerCube)() = spawnPos;
            m_registry.Get<vve::Rotation &>(m_eulerCube)() = mat3_t{1.0f};
            m_registry.Get<vve::Scale &>(m_eulerCube)()    = vec3_t{1.0f};
        }

        m_eulerActive = true;
    }

    void DropCubeSymplectic(float delta_t)
    {
        const glmvec3 gravityPhys{0.0f, m_physics.c_gravity, 0.0f};
        const glmvec3 gravityGame = fromPhysics(gravityPhys);
        m_upGame = -glm::normalize(gravityGame);

        // Start at 100 m, zero initial velocity
        m_dropHeight = 100.0f;    // gravity magnitude
        m_sympDt = std::max(1e-6f, delta_t);
        m_sympT = 0.0f;
        m_sympH = m_dropHeight;
        m_sympV = 0.0f;
        m_sympAccum = 0.0;

        const glmvec3 spawnPos = m_sympH * m_upGame;

        // Create or reuse a render-only cube
        if (!m_sympCube.IsValid())
        {
            m_sympCube = m_engine.CreateScene(vve::Name{},
                                              vve::ParentHandle{},
                                              vve::Filename{cube_obj}, aiProcess_FlipWindingOrder,
                                              vve::Position{spawnPos},
                                              vve::Rotation{mat3_t{1.0f}},
                                              vve::Scale{vec3_t{1.0f}});
        }
        else
        {
            m_registry.Get<vve::Position &>(m_sympCube)() = spawnPos;
            m_registry.Get<vve::Rotation &>(m_sympCube)() = mat3_t{1.0f};
            m_registry.Get<vve::Scale &>(m_sympCube)()    = vec3_t{1.0f};
        }

        m_sympActive = true;
    }


};
int main()
{
    vve::Engine engine("My Engine", vve::RendererType::RENDERER_TYPE_FORWARD);
    MyGame mygui{engine};
    engine.Run();

    return 0;
}
