#include "Application.h"

#include "Logger.h"

typedef unsigned int uint;

glm::dvec2 last_cursor_position;
glm::dvec2 offset_cursor = glm::vec2(0.0f, 0.0f);

Application::Application(uint16_t width, uint16_t height, const char* title, uint16_t version_major, uint16_t version_minor)
    : m_Width(width), m_Height(height), m_Title(title),
    m_VERSION_MAJOR(version_major), m_VERSION_MINOR(version_minor),
    m_ObjHandler(m_Window)
{
    assert(sizeof(GLuint) == sizeof(uint));

    Logger l(Logger::Filter::FILTER_ERROR);

    // setup glfw
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, version_major);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, version_minor);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    const char* glsl_version = "#version 330 core";

    m_Window = glfwCreateWindow(width, height, title, nullptr, nullptr);

    if (!m_Window) // == nullptr
    {
        l.Log("Unable to create GLFW window\n", Logger::Log_Level::ERROR, Logger_Color{ Logger_Color::Colors::RED });
        glfwTerminate();
        return;
    }
    glfwMakeContextCurrent(m_Window);
    glfwSwapInterval(1);

    m_ObjHandler.UpdateParentWindow(m_Window);

    // make sure glad is initialized
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        l.Log("Unable to initialize GLAD\n", Logger::Log_Level::ERROR, Logger_Color{ Logger_Color::Colors::RED });
        return;
    }

    GLCall(glViewport(0, 0, width, height));
    glfwSetFramebufferSizeCallback(m_Window, framebuffer_size_callback);
    glfwSetWindowCloseCallback(m_Window, window_close_callback);

    glfwGetCursorPos(m_Window, &last_cursor_position.x, &last_cursor_position.y);
    glfwSetCursorPosCallback(m_Window, mouse_callback);




    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    m_IO = std::addressof(ImGui::GetIO());
    //ImGuiIO& io = ImGui::GetIO(); (void)io;
    m_IO->ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    m_IO->ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls
    m_IO->ConfigFlags |= ImGuiConfigFlags_DockingEnable;         // Enable Docking
    m_IO->ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;       // Enable Multi-Viewport / Platform Windows
    //m_IO->ConfigViewportsNoAutoMerge = true;
    //m_IO->ConfigViewportsNoTaskBarIcon = true;

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsLight();

    // When viewports are enabled we tweak WindowRounding/WindowBg so platform windows can look identical to regular ones.
    ImGuiStyle& style = ImGui::GetStyle();
    if (m_IO->ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
    {
        style.WindowRounding = 0.0f;
        style.Colors[ImGuiCol_WindowBg].w = 1.0f;
    }

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(m_Window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);


    GLCall(glEnable(GL_BLEND));
    GLCall(glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));
    GLCall(glEnable(GL_DEPTH_TEST));
}

Application::~Application()
{
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwTerminate();
}

void Application::Run()
{
    float currentTime = glfwGetTime();
    m_DeltaTime = currentTime - m_LastFrame;
    m_LastFrame = currentTime;

    //Update(1.0f / m_IO->Framerate);
    Update(m_DeltaTime);


    GLCall(glClearColor(0.2f, 0.3f, 0.3f, 1.0f));
    m_Renderer.Clear();


    // Start the Dear ImGui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();


    OnRender();

    OnImGuiRender();

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    //ImGuiIO& io = ImGui::GetIO(); (void)io;
    if (m_IO->ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
    {
        GLFWwindow* backup_current_context = glfwGetCurrentContext();
        ImGui::UpdatePlatformWindows();
        ImGui::RenderPlatformWindowsDefault();
        glfwMakeContextCurrent(backup_current_context);
    }


    glfwSwapBuffers(m_Window);
}


void Application::Update(float dt)
{
    glfwPollEvents();

    // handle inputs
    ProcessInputs(dt);
    

    m_ObjHandler.UpdateObjects(dt);
}

void Application::OnRender()
{
    // some stuff
    m_ObjHandler.RenderObjects();
}

void Application::OnImGuiRender()
{
    m_ObjHandler.ImGuiRender();
    m_ObjHandler.RenderImGuiObjects();

    ImGui::ShowDemoWindow(); // Show demo window! :)

    {
        ImGui::Begin("Specs");
        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / m_IO->Framerate, m_IO->Framerate);
        ImGui::End();
    }
}


void Application::ProcessInputs(float deltaTime)
{
    if (ImGui::IsKeyPressed(ImGuiKey_Escape))
    {
        m_ShouldClose = true;
        glfwSetWindowShouldClose(m_Window, true);
    }

    if (ImGui::IsKeyPressed(ImGuiKey_P))
    {
        const auto cursor_state = m_IsCursorActive ? GLFW_CURSOR_DISABLED : GLFW_CURSOR_NORMAL;
        glfwSetInputMode(m_Window, GLFW_CURSOR, cursor_state);

        m_IsCursorActive = !m_IsCursorActive;
    }

    m_ObjHandler.processInputs(deltaTime);
}


void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    offset_cursor.x = xpos - last_cursor_position.x;
    offset_cursor.y = ypos - last_cursor_position.y;

    last_cursor_position.x = xpos;
    last_cursor_position.y = ypos;
}


void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    GLCall(glViewport(0, 0, width, height));
}


void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
}

void window_close_callback(GLFWwindow* window)
{
    glfwSetWindowShouldClose(window, GLFW_TRUE);
}
