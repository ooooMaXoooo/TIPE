#include "Renderer.h"

#include <iostream>

void GLClearError()
{
    while (glGetError() != GL_NO_ERROR);
}

bool GLLogCall(const char* functionName, const char* fileName, int lineNumber)
{
    while (GLenum error = glGetError())
    {
        std::cout << "[OpenGL Error] (error " << error << ")\n\tfunction : " << functionName
            << "\n\t" << fileName << ": line." << lineNumber << std::endl;
        return false;
    }
    return true;
}

Renderer::Renderer(GLenum mode)
    : m_Mode(mode)
{
}

void Renderer::Clear() const
{
    glClear(GL_COLOR_BUFFER_BIT);
}

void Renderer::Draw(const VertexArray& va, const IndexBuffer& ib, const Shader& shader) const
{
    //GLCall(glPolygonMode(GL_FRONT_AND_BACK, GL_LINE));

    GLCall(glPolygonMode(GL_FRONT_AND_BACK, m_Mode));
    shader.Bind();
    va.Bind();
    ib.Bind();

    GLCall(glDrawElements(GL_TRIANGLES, ib.GetCount(), GL_UNSIGNED_INT, nullptr));
}