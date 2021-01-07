

// Update velocities of a 1-D array of particles
template <size_t N>
void updateVelocity(std::array<float, N> &v,std::array<float, N> &position, float boundary, std::array<float, N+1> &k, float m, float b, int step)
{

    int n = v.size();

    for (int i = 0; i < n; i++)
    {
        float x = position[i];
        float a, xleft, xright = 0;
        if (i > 0)
        { // displacement from equilibrium relative to particle to the left
            xleft = x - (position[i - 1] + (2.0f / (n + 1)));
        }
        else
        { // displacement from equilibrium relative to left wall
            xleft = x - ((2.0f / (n + 1)) - boundary);
        }
        if (i < n - 1)
        { // displacement from equilibrium relative to particle to the right
            xright = x - (position[i + 1] - (2.0f / (n + 1)));
        }
        else
        { // displacement from equilibrium relative to right wall
            xright = x - (boundary - (2.0f / (n + 1)));
        }

        a = (((-k[i+1] * xright) - (k[i] * xleft)) / m) - (v[i] * b); // calculate acceleration

        if (x < -1)
        { // if particle hits the left wall, stop
            v[i] += -boundary - x;
        }
        else if (x > 1)
        { // if particle hits the right wall, stop
            v[i] += boundary - x;
        }
        else
        {
            v[i] += a / step; // update velocity dividing by frame rate
        }
    }
}



// Update velocities of a 2-D array of particles
template <size_t N>
void updateVelocity(std::array<std::array<float, N>, 2> &v, std::array<std::array<float, N>, 2> &position, float boundary, float k, float m, float b, int step)
{

    int n = v.size();

    for (int i = 0; i < n; i++)
    {
        float x = position[i];
        float a, xleft, xright = 0;
        if (i > 0)
        { // displacement from equilibrium relative to particle to the left
            xleft = x - (position[i - 1] + (2.0f / (n + 1)));
        }
        else
        { // displacement from equilibrium relative to left wall
            xleft = x - ((2.0f / (n + 1)) - boundary);
        }
        if (i < n - 1)
        { // displacement from equilibrium relative to particle to the right
            xright = x - (position[i + 1] - (2.0f / (n + 1)));
        }
        else
        { // displacement from equilibrium relative to right wall
            xright = x - (boundary - (2.0f / (n + 1)));
        }

        a = (((-k * xright) - (k * xleft)) / m) - (v[i] * b); // calculate acceleration

        if (x < -1)
        { // if particle hits the left wall, stop
            v[i] += -boundary - x;
        }
        else if (x > 1)
        { // if particle hits the right wall, stop
            v[i] += boundary - x;
        }
        else
        {
            v[i] += a / step; // update velocity dividing by frame rate
        }
    }
}