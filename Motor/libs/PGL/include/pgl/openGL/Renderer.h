#pragma once

#include "IndexBuffer.h"
#include "Shader.h"
#include "VertexArray.h"
#include "VertexBufferLayout.h"

#include "glad/glad.h"
#include "GLFW/glfw3.h"

#include "ImGui/imgui.h"
#include "ImGui/imgui_impl_glfw.h"
#include "ImGui/imgui_impl_opengl3.h"



#define ASSERT(x) if (!(x)) __debugbreak();
#define GLCall(x) GLClearError();\
                  x;\
                  ASSERT(GLLogCall(#x, __FILE__, __LINE__))

void GLClearError();
bool GLLogCall(const char* functionName, const char* fileName, int lineNumber);



class Renderer
{
private :
	GLenum m_Mode;

	bool m_skeleton = false;

public :
	/// <summary>
	/// 
	/// </summary>
	/// <param name="mode">mode of drawing can be GL_FILL or GL_LINE or ...</param>
	Renderer(GLenum mode = GL_FILL);

	void Clear() const;
	void Draw(const VertexArray& va, const IndexBuffer& ib, const Shader& shader) const;

	void OnImGuiRender();
};

