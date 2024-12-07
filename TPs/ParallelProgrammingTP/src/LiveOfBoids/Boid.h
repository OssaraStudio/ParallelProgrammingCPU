/*
 * Boid.h
 *
 *  Created on: Dec 7, 2024
 *      Author: ossaradj
 */

#ifndef SRC_LIVEOFBOIDS_BOID_H_
#define SRC_LIVEOFBOIDS_BOID_H_

#include "Vector2D.h"
#include <vector>


class Boid {
    public:
        // Constructor
        Boid(Vector2D position = Vector2D(), Vector2D velocity = Vector2D())
        : m_position(position), m_velocity(velocity) {};

        // Overloading stream operator
        friend std::ostream& operator<<(std::ostream& out, const Boid& b)
        {
            out << "position = " << b.m_position << " velocity = " << b.m_velocity;
            return out;
        };

        // Getters
        Vector2D getPosition() const { return m_position; };
        float getPositionX() const { return m_position.getX(); };
        float getPositionY() const { return m_position.getY(); };
        Vector2D getVelocity() const { return m_velocity; };
        float getVelocityX() const { return m_velocity.getX(); };
        float getVelocityY() const { return m_velocity.getY(); };

        // Setters
        void setPosition(Vector2D const& position) { this->m_position = position; };
        void setPositionX(const float x) { this->m_position.setX(x); };
        void setPositionY(const float y) { this->m_position.setY(y); };
        void setVelocity(Vector2D const& velocity) { this->m_velocity = velocity; };
        void setVelocityX(const float x) { this->m_velocity.setX(x); };
        void setVelocityY(const float y) { this->m_velocity.setY(y); };

        std::vector<Boid> get_neighbors(const std::vector<Boid>& boids, float radius) const
        {
            std::vector<Boid> neighbors;
            for(auto& other: boids)
            {
                if(this->m_position.euclidean_dist(other.m_position) < radius)
                neighbors.push_back(other);
            }
            return neighbors;
        };

    private:
        Vector2D m_position;
        Vector2D m_velocity;
};

#endif /* SRC_LIVEOFBOIDS_BOID_H_ */