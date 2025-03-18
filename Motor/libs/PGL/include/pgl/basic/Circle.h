#pragma once

#include "Drawable.h"

#include "glm/glm.hpp"


namespace pgl {

    constexpr float PI = glm::pi<float>();

    class Circle : public Drawable
    {
    protected:
        float m_Radius;

        glm::vec3 m_Center;
        
    public:

        Circle(glm::vec3 center, float radius);


        virtual void OnRender() override;
        virtual void OnImGuiRender() override;



        void SetRadius(float radius);
        void SetPosition(glm::vec3 center);

        void SetScale(glm::vec3 scale);
        void SetRotation(glm::vec3 rotation);



    private:
        float m_Rotations[3];
        //float m_Translations[3];
        float m_Scale[3];

    private:
        void Init();
        void UpdateTransform();
    };

};