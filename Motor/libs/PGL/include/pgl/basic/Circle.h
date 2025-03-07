#pragma once

#include "Drawable.h"

#include "glm/glm.hpp"

namespace pgl {

    class Circle : public Drawable
    {
    protected:
        float m_Radius;

        glm::vec2 m_Center;
        
    public:

        Circle(glm::vec2 center, float radius);


        virtual void OnRender() override;
        virtual void OnImGuiRender() override;

    private:
        float m_Rotations[3];
        float m_Translations[3];
    };

};