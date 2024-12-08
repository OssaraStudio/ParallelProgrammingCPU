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
#include "tbb/tbb.h"
#include "omp.h"
#include "Boid.h"

class BoidGenerator
{
  public:
    BoidGenerator(){}
    virtual ~BoidGenerator(){}

    void setChunkSize(int chunk_size)
    {
        m_chunk_size = chunk_size ;
    };

    void generate(int nb, std::vector<Boid>& boids)
    {
        std::random_device rd;
        std::mt19937 eng(rd());        
        std::uniform_real_distribution<float> random_x(0.0f, static_cast<float>(600));
        std::uniform_real_distribution<float> random_y(0.0f, static_cast<float>(400));
        std::uniform_real_distribution<float> angle_space(0.0f, 2.0f * M_PI);
        std::uniform_real_distribution<float> velocity_space(0.1f, 0.5f);

        for (int i=0; i<nb; ++i)
        {
            auto angle = angle_space(eng);                    // Random direction
            auto speed = velocity_space(eng);                // Random speed
            boids.push_back(Boid(Vector2D{random_x(eng), random_y(eng)}, Vector2D{speed * std::cos(angle), speed * std::sin(angle)}));
        }
    };

    void findNeighbors(std::vector<Boid>& boids, std::vector<int>& y, float radius)
    {
        for(int i=0; i<boids.size(); ++i)
        {
            y[i] = boids[i].get_neighbors(boids, radius).size() ;
        }
    };

    void omptaskfindNeighbors(std::vector<Boid>& boids, std::vector<int>& y, float radius)
    {
        std::size_t nb_task = (boids.size()+m_chunk_size-1)/m_chunk_size ;
        
        #pragma omp parallel for
        for(int i = 0; i < boids.size(); ++i)
        {
            y[i] = boids[i].omptaskget_neighbors(boids, radius).size();
        }
        
    };

    private:
        int m_chunk_size = 1 ;

};

#endif /* SRC_LIVEOFBOIDS_BOIDGENERATOR_H_ */