#pragma once

#include "pch.h"

#include "Motor/Core/ObjectsHandler.h"

class Application {
private:

	bool m_ShouldClose = false;
	uint16_t m_Width;
	uint16_t m_Height;
	const char* m_Title;
	const uint16_t m_VERSION_MAJOR;
	const uint16_t m_VERSION_MINOR;

	GLFWwindow* m_Window = nullptr;

	ImGuiIO* m_IO = nullptr;

	// cursor state
	bool m_IsCursorActive = true;


	// the renderer of the application
	// the gloabl renderer
	Renderer m_Renderer;

	float m_DeltaTime = 0.0f;	// Time between current frame and last frame
	float m_LastFrame = 0.0f;   // Time of last frame


	// motor specific variables
	Motor::Core::ObjectsHandler m_ObjHandler;


public:

	Application(uint16_t width, uint16_t height, const char* title, uint16_t version_major = 3, uint16_t version_minor = 3);
	~Application();

	bool ShouldClose() const { return glfwWindowShouldClose(m_Window); }
	float ImGuiFramerate() const { return m_IO->Framerate; }

	void Run();





private:

	void Update(float dt);
	void OnRender();
	void OnImGuiRender();

	void ProcessInputs(float deltaTime);
};

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void window_close_callback(GLFWwindow* window);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);