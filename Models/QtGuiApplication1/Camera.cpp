#include "Camera.h"
#include <QDebug>
const float PI = 3.1415926535;

Camera::Camera(QVector3D position, QVector3D up, float yaw, float pitch) :
	position(position),
	worldUp(up),
	front(-position),
	picth(pitch),
	yaw(yaw),
	movementSpeed(SPEED),
	mouseSensitivity(0.01),
	zoom(ZOOM) 
{
	this->updateCameraVectors();

	for (uint i = 0; i != 1024; ++i)
		keys[i] = false;
}

Camera::~Camera()
{

}

// Returns the view matrix calculated using Euler Angles and the LookAt Matrix
QMatrix4x4 Camera::getViewMatrix()
{
	QMatrix4x4 view;
	view.lookAt(this->position, this->position + this->front, this->up);
	return view;
}

// Processes input received from any keyboard-like input system. Accepts input parameter in the form of camera defined ENUM (to abstract it from windowing systems)
void Camera::processKeyboard(Camera_Movement direction, float deltaTime)
{
	float velocity = this->movementSpeed * deltaTime;
	if (direction == FORWARD)
		this->position += this->front * velocity;
	if (direction == BACKWARD)
		this->position -= this->front * velocity;
	if (direction == LEFT)
		this->position -= this->right * velocity;
	if (direction == RIGHT)
		this->position += this->right * velocity;
	if (direction == UP)
		this->position += this->worldUp * velocity;
	if (direction == DOWN)
		this->position -= this->worldUp * velocity;
}

// Processes input received from a mouse input system. Expects the offset value in both the x and y direction.
void Camera::processMouseMovement(float xoffset, float yoffset, bool constraintPitch)
{
	xoffset *= this->mouseSensitivity;
	yoffset *= this->mouseSensitivity;

	this->yaw += xoffset;
	this->picth += yoffset;

	if (constraintPitch) {
		if (this->picth > 89.0f)
			this->picth = 89.0f;
		if (this->picth < -89.0f)
			this->picth = -89.0f;
	}

	this->updateCameraVectors();
}

// Processes input received from a mouse scroll-wheel event. Only requires input on the vertical wheel-axis
void Camera::processMouseScroll(float yoffset)
{
	if (this->zoom >= 1.0f && this->zoom <= 45.0f)
		this->zoom -= yoffset;
	if (this->zoom > 45.0f)
		this->zoom = 45.0f;
	if (this->zoom < 1.0f)
		this->zoom = 1.0f;
}

void Camera::processInput(float dt)
{

	if (keys[Qt::Key_W])
		processKeyboard(FORWARD, dt);
	if (keys[Qt::Key_S])
		processKeyboard(BACKWARD, dt);
	if (keys[Qt::Key_A])
		processKeyboard(LEFT, dt);
	if (keys[Qt::Key_D])
		processKeyboard(RIGHT, dt);
	if (keys[Qt::Key_E])
		processKeyboard(UP, dt);
	if (keys[Qt::Key_Q])
		processKeyboard(DOWN, dt);
}

void Camera::updateCameraVectors()
{
	// Calculate the new Front vector
	QVector3D front;
	front.setX(cos(this->yaw*PI/180) * cos(this->picth*PI / 180));
	front.setY(sin(this->picth*PI/180));
	front.setZ(sin(this->yaw*PI/180) * cos(this->picth*PI/180));
	this->front = front.normalized();
	this->right = QVector3D::crossProduct(this->front, this->worldUp).normalized();
	this->up = QVector3D::crossProduct(this->right, this->front).normalized();
}

void Camera::setMouseSpeed(float a)
{
	movementSpeed = a;
}

void Camera::setLocation(QVector3D Position, QVector3D up, float Yaw, float Pitch)
{
	position = Position;
	worldUp = up;
	front = -Position;
	picth = Pitch;
	yaw = Yaw;
	movementSpeed = SPEED;
	mouseSensitivity = 0.01;
	zoom = ZOOM;
	this->updateCameraVectors();
}
