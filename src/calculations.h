#include "./structs.h"

std::array<double, 3> calculateForcesFixed(Particle &first, Particle &second, double springLength,
                                           double k, bool freedom[3]) {
  std::array<double, 3> forceComponents = {0, 0, 0};      // Force Components
  int dimensions = freedom[0] + freedom[1] + freedom[2];  // how many axes are activated
  if (dimensions == 0) {
    for (int i = 0; i < 3; i++) first.velocity[i] = 0;
    return forceComponents;  // return no force if no axis activated
  }

  double x = (second.displacement[0] - first.displacement[0]) * k;
  double y = (second.displacement[1] - first.displacement[1]) * k;
  double z = (second.displacement[2] - first.displacement[2]) * k;

  if (freedom[0]) forceComponents[0] = x;
  if (freedom[1]) forceComponents[1] = y;
  if (freedom[2]) forceComponents[2] = z;
  return forceComponents;
}

// Update velocities of a 2-D array of particles
void updateVelocities(std::vector<std::vector<Particle>> &particle, double springLength,
                      bool freedom[3], std::vector<double> &kX, std::vector<double> &kY, float m,
                      float b, int step) {
  int NX = particle.size();
  int NY = particle[0].size();
  for (int x = 0; x < NX; x++)  // set all acceleration to zero
    for (int y = 0; y < NY; y++) {
      for (int j = 0; j < 3; j++) particle[x][y].acceleration[j] = 0;
    }
  for (int y = 0; y < NY - 1; y++)
    for (int x = 0; x < NX - 1; x++) {
      std::array<double, 3> forcesX =
        calculateForcesFixed(particle[x][y], particle[x + 1][y], springLength, kX[x], freedom);

      if (NY > 3) {  // only calculate Y forces if 2d array (NY > 1 + 2 boundary particles)
        std::array<double, 3> forcesY =
          calculateForcesFixed(particle[x][y], particle[x][y + 1], springLength, kY[y], freedom);
        for (int j = 0; j < 3; j++) {
          particle[x][y].acceleration[j] += forcesY[j];      // force due to above particle spring
          particle[x][y + 1].acceleration[j] -= forcesY[j];  // opposite force on above particle
        }
      }

      for (int j = 0; j < 3; j++) {
        particle[x][y].acceleration[j] += forcesX[j];      // force due to right particle spring
        particle[x + 1][y].acceleration[j] -= forcesX[j];  // opposite force on right particle

        particle[x][y].acceleration[j] = (particle[x][y].acceleration[j] / m) -
                                         (particle[x][y].velocity[j] * b);  // (kx / m) - (vb)
        particle[x][y].velocity[j] +=
          particle[x][y].acceleration[j] / step;  // update velocity dividing by frame rate
      }
    }
}

// // Update velocities of a 1-D array of particles
// void updateVelocities(std::vector<Particle> &particle, double springLength, bool freedom[3],
//                       std::vector<double> &k, float m, float b, int step) {
//   int N = particle.size();
//   for (int i = 0; i < N; i++)  // set all acceleration to zero
//   {
//     for (int j = 0; j < 3; j++) particle[i].acceleration[j] = 0;
//   }

//   for (int i = 0; i < N; i++) {
//     if (i < N - 1) {
//       std::array<double, 3> forces =
//         calculateForcesFixed(particle[i], particle[i + 1], springLength, k[i], freedom);
//       for (int j = 0; j < 3; j++) {
//         particle[i].acceleration[j] += forces[j];      // acceleration due to right particle
//         spring particle[i + 1].acceleration[j] -= forces[j];  // opposite acceleration on right
//         particle particle[i].acceleration[j] =
//           (particle[i].acceleration[j] / m) - (particle[i].velocity[j] * b);  // (kx / m) - (vb)
//         particle[i].velocity[j] +=
//           particle[i].acceleration[j] / step;  // update velocity dividing by frame rate
//                                                // if (!freedom[j]) particle[i].velocity[j] = 0;
//       }
//     }
//   }
// }

// std::array<double, 3> calculateForcesOmni(
//   Particle &first, Particle &second, double springLength, double k,
//   bool freedom[3]) {  // Calculate 3D forces for non-fixed spring

//   std::array<double, 3> forceComponents = {0, 0, 0};      // Force Components
//   int dimensions = freedom[0] + freedom[1] + freedom[2];  // how many axes are activated
//   if (dimensions == 0) {
//     for (int i = 0; i < 3; i++) first.velocity[i] = 0;
//     return forceComponents;  // return no force if no axis activated
//   }

//   double x = second.x() - first.x();
//   double y = second.y() - first.y();
//   double z = second.z() - first.z();

//   double r = sqrtf(powf(x, 2) + powf(y, 2) + powf(z, 2));  // magnitude of force

//   if (r == 0) x, y, z, r = 0.0001;

//   x = (x / r) * ((r - springLength) * k);
//   y = (y / r) * ((r - springLength) * k);
//   z = (z / r) * ((r - springLength) * k);

//   if (freedom[0]) forceComponents[0] = x;
//   if (freedom[1]) forceComponents[1] = y;
//   if (freedom[2]) forceComponents[2] = z;
//   return forceComponents;
// }