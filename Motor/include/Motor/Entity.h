#pragma once

#include "OpenGL/Renderer.h"

#include "pch.h"

namespace Motor {
    class Entity : pgl::Drawable {
    protected:
        

    public:


        virtual ~Entity();

        virtual void Update(float ts);
        virtual void OnRender() const;
    };

}

