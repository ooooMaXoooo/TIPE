#pragma once

#include "glad\glad.h"
#include "GLFW\glfw3.h"

#include "ImGui\imgui.h"
#include "ImGui\imgui_impl_glfw.h"
#include "ImGui\imgui_impl_opengl3.h"

#include "OpenGL\Renderer.h"


#include "basic\Circle.h"


namespace pgl {

	class Application
	{
	private:

		bool m_ShouldClose = false;
		uint16_t m_Width;
		uint16_t m_Height;
		const char* m_Title;
		const uint16_t m_VERSION_MAJOR;
		const uint16_t m_VERSION_MINOR;

		GLFWwindow* m_Window = nullptr;

		ImGuiIO* m_IO = nullptr;

		Renderer m_Renderer;

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

		void ProcessInputs();

	private :
		std::unique_ptr<Circle> c;
	};

};

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
