#include "./structs.h"

std::array<float, 3> calculateForces(al::Pose first, al::Pose second, float springLength, bool freedom[3]) {
    
    std::array<float, 3> forceComponents = {0, 0, 0}; // Force Components
    int dimensions = freedom[0] + freedom[1] + freedom[2]; // how many axes are activated
    if (dimensions == 0) return forceComponents; // return no force if no axis activated
    float directionOne, directionTwo, directionThree = 0; 
    float r = 0; // total magnitude of force
    
    float x = second.x() - first.x();
    float y = second.y() - first.y();
    float z = second.z() - first.z();

    if (dimensions == 1) {
        directionOne = freedom[0] ? x : freedom[1] ? y : z; // assign direction to whichever axis is activated
        forceComponents[0] = directionOne - springLength;
        return forceComponents;
    }
    
    if (dimensions == 2) {
        directionOne = freedom[0] ? x : y;
        directionTwo = (freedom[0] ? (freedom[1] ? y : z) : z);
        r = sqrt(pow(directionOne, 2) + pow(directionTwo, 2));
        forceComponents[0] = directionOne / r;
        forceComponents[1] = directionTwo / r;
        return forceComponents;
    }

    if (dimensions == 3) {
        r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    }


}

// Update velocities of a 1-D array of particles
template <size_t N>
void updateVelocity(std::array<Particle, N> &particle, float boundary, float springLength, bool freedom[3], std::array<float, N + 1> &k, float m, float b, int step)
{
    for (int i = 0; i < N; i++) particle[i].acceleration[0] = 0;
    for (int i = 0; i < N; i++) particle[i].acceleration[1] = 0;
    for (int i = 0; i < N; i++)
    {
        float x = particle[i].particle.pose.get().x();
        float y = particle[i].particle.pose.get().y();
        if (i < N - 1) {
        std::array<float, 3> forces = calculateForces(particle[i].particle.pose.get(), particle[i + 1].particle.pose.get(), springLength, freedom);
        particle[i].acceleration[0] += k[i + 1] * forces[0]; // acceleration due to right particle spring
        particle[i + 1].acceleration[0] += -k[i + 1] * forces[0]; // opposite acceleration on right particle
        particle[i].acceleration[1] += k[i + 1] * forces[1]; // acceleration due to right particle spring
        particle[i + 1].acceleration[1] += -k[i + 1] * forces[1]; // opposite acceleration on right particle
        }
        if (i == 0) particle[i].acceleration[0] += -k[0] * (x - (springLength - boundary)); // acceleration due to left wall spring
        if (i == N - 1) particle[i].acceleration[0] += -k[i + 1] * (x - (boundary - springLength)); // acceleration due to right wall spring
        if (i == 0) particle[i].acceleration[1] += -k[0] * (x - (springLength - boundary)); // acceleration due to left wall spring
        if (i == N - 1) particle[i].acceleration[1] += -k[i + 1] * (x - (boundary - springLength)); // acceleration due to right wall spring
        
        particle[i].acceleration[0] = (particle[i].acceleration[0] / m) - (particle[i].velocity[0] * b); // (kx / m) - (vb)
        particle[i].acceleration[1] = (particle[i].acceleration[1] / m) - (particle[i].velocity[1] * b); // (kx / m) - (vb)

        // std::cout << particle[1].velocity[0] << std::endl;
        if (x < -1)
        { // if particle hits the left wall, stop
            particle[i].velocity[0] += -boundary - x;
        }
        else if (x > 1)
        { // if particle hits the right wall, stop
            particle[i].velocity[0] += boundary - x;
        }
        else
        {
            particle[i].velocity[0] += particle[i].acceleration[0] / step; // update velocity dividing by frame rate
            particle[i].velocity[1] += particle[i].acceleration[1] / step; // update velocity dividing by frame rate

        }
    }
}
