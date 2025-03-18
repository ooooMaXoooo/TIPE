#pragma once

#include "OpenGL/Renderer.h"

#include "pch.h"

namespace Motor {
    class Entity {
    protected:
        

    public:


        virtual ~Entity();

        virtual void Update(float ts);
    };

}

