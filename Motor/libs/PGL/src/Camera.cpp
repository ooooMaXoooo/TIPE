#include "Camera.h"

namespace pgl {
	/*Camera::Camera(const glm::vec3& position, const glm::vec3& target)
		: m_Pos(position), m_Target(target)
	{
		m_Direction = glm::normalize(position - target);
		// il faut noter que m_Direction n'est pas la direction vers laquelle on veut regarder mais l'opposer parfait, si on
		//		centre un repère sur la cible, le vecteur m_Direction correspondrai au vecteur ur de la base sphérique



		glm::vec3 up(0.0f, 1.0f, 0.0f);
		glm::vec3 right(1.0f, 0.0f, 0.0f);

		// On créer la base local de la caméra
		if (glm::abs(m_Direction.y) - 1.0f < 1e-3) {
			// m_direction est le vecteur unitaire uz
			m_Up = glm::vec3(0.0f, 0.0f, -m_Direction.y);
			m_Right = right;
		}
		else {
			m_Right = glm::normalize(glm::cross(up, m_Direction));
			m_Up = glm::normalize(glm::cross(m_Direction, m_Right));
		}

		m_View = glm::lookAt(m_Pos, m_Pos + m_Front, m_Up);
	}*/

    // constructor with vectors
    Camera::Camera(glm::vec3 position, glm::vec3 up, float yaw, float pitch)
        : m_Front(glm::vec3(0.0f, 0.0f, -1.0f)), m_MovementSpeed(SPEED), m_MouseSensitivity(SENSITIVITY), m_Zoom(ZOOM)
    {
        m_Position = position;
        m_WorldUp = up;
        m_Yaw = yaw;
        m_Pitch = pitch;
        updateCameraVectors();
    }
    // constructor with scalar values
    Camera::Camera(float posX, float posY, float posZ, float upX, float upY, float upZ, float yaw, float pitch) :
        m_Front(glm::vec3(0.0f, 0.0f, -1.0f)), m_MovementSpeed(SPEED), m_MouseSensitivity(SENSITIVITY), m_Zoom(ZOOM)
    {
        m_Position = glm::vec3(posX, posY, posZ);
        m_WorldUp = glm::vec3(upX, upY, upZ);
        m_Yaw = yaw;
        m_Pitch = pitch;
        updateCameraVectors();
    }

    // copy constructor
    Camera::Camera(const Camera& cam)
    {
        // camera Attributes
        m_Position = cam.m_Position;
        m_Front = cam.m_Front;
        m_Up = cam.m_Up;
        m_Right = cam.m_Right;
        m_WorldUp = cam.m_WorldUp;

        // euler Angles
        m_Yaw = cam.m_Yaw;
        m_Pitch = cam.m_Pitch;

        // camera options
        m_MovementSpeed = cam.m_MovementSpeed;
        m_MouseSensitivity = cam.m_MouseSensitivity;
        m_Zoom = cam.m_Zoom;

    }

    // returns the view matrix calculated using Euler Angles and the LookAt Matrix
    glm::mat4 Camera::GetViewMatrix() const
    {
        return glm::lookAt(m_Position, m_Position + m_Front, m_Up);
    }

    // returns the fov in order to enable the zoom
    float Camera::GetFov() const {
        return m_Zoom;
    }


    // processes input received from any keyboard-like input system. Accepts input parameter in the form of camera defined ENUM (to abstract it from windowing systems)
    void Camera::ProcessKeyboard(Camera_Movement direction, float deltaTime)
    {
        float velocity = m_MovementSpeed * deltaTime;
        if (direction == FORWARD)
            m_Position += m_Front * velocity;
        if (direction == BACKWARD)
            m_Position -= m_Front * velocity;
        if (direction == LEFT)
            m_Position -= m_Right * velocity;
        if (direction == RIGHT)
            m_Position += m_Right * velocity;
    }

    // processes input received from a mouse input system. Expects the offset value in both the x and y direction.
    void Camera::ProcessMouseMovement(float xoffset, float yoffset, bool constrainPitch)
    {
        xoffset *= m_MouseSensitivity;
        yoffset *= m_MouseSensitivity;

        m_Yaw += xoffset;
        m_Pitch += yoffset;

        // make sure that when pitch is out of bounds, screen doesn't get flipped
        if (constrainPitch)
        {
            if (m_Pitch > 89.0f)
                m_Pitch = 89.0f;
            if (m_Pitch < -89.0f)
                m_Pitch = -89.0f;
        }

        // update Front, Right and Up Vectors using the updated Euler angles
        updateCameraVectors();
    }

    // processes input received from a mouse scroll-wheel event. Only requires input on the vertical wheel-axis
    void Camera::ProcessMouseScroll(float yoffset)
    {
        m_Zoom -= (float)yoffset;
        if (m_Zoom < 1.0f)
            m_Zoom = 1.0f;
        if (m_Zoom > 45.0f)
            m_Zoom = 45.0f;
    }

    // calculates the front vector from the Camera's (updated) Euler Angles
    void Camera::updateCameraVectors()
    {
        // calculate the new Front vector
        glm::vec3 front;
        front.x = cos(glm::radians(m_Yaw)) * cos(glm::radians(m_Pitch));
        front.y = sin(glm::radians(m_Pitch));
        front.z = sin(glm::radians(m_Yaw)) * cos(glm::radians(m_Pitch));
        m_Front = glm::normalize(front);
        // also re-calculate the Right and Up vector
        m_Right = glm::normalize(glm::cross(m_Front, m_WorldUp));  // normalize the vectors, because their length gets closer to 0 the more you look up or down which results in slower movement.
        m_Up = glm::normalize(glm::cross(m_Right, m_Front));
    }
}