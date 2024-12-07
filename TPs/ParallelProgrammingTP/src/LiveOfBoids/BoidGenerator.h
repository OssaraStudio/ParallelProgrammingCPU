/*
 * BoidGenerator.h
 *
 *  Created on: Dec 7, 2024
 *      Author: ossaradj
 */

#ifndef SRC_LIVEOFBOIDS_BOIDGENERATOR_H_
#define SRC_LIVEOFBOIDS_BOIDGENERATOR_H_

#include <vector>
#include <random>
#include "Boid.h"

class BoidGenerator
{
  public:
    BoidGenerator(){}
    virtual ~BoidGenerator(){}

    void generate(int nb, std::vector<Boid> boids)
    {
        std::random_device rd;
        std::mt19937 eng(rd());        
        std::uniform_real_distribution<float> random_x(0.0f, static_cast<float>(1500));
        std::uniform_real_distribution<float> random_y(0.0f, static_cast<float>(900));
        std::uniform_real_distribution<float> angle_space(0.0f, 2.0f * M_PI);
        std::uniform_real_distribution<float> velocity_space(0.1f, 0.5f);

        for (int i=0; i<nb; ++i)
        {
            auto angle = angle_space(eng);                    // Random direction
            auto speed = velocity_space(eng);                // Random speed
            std::cout << random_x(eng) << "\n" ;
            boids.push_back(Boid(Vector2D{random_x(eng), random_y(eng)}, Vector2D{speed * std::cos(angle), speed * std::sin(angle)}));
        }
    };

};

#endif /* SRC_LIVEOFBOIDS_BOIDGENERATOR_H_ */