#include "MAC.h"
#include "Utils.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/CXX11/Tensor>  

class FLIP{
    struct Particle
    {
        bool isActive = false; //just for ease of allocation;
        double x;
        double y;
        double z;

        double u;
        double v;
        double w;
    };

    Particle* particles;
    std::vector<std::vector<std::vector<std::vector<int>>>> SPACE_HASH; //this will be a bit inneficient for my sanity, there is no way im using a linked list for this while maintaining this whole code.


    



    static void InitializeADI(MAC* grid,int particlePerCell,double dt,Vec3(*VelocityBorderFunction)(double, double, double,double),Vec3(*VelocityFont)(double, double, double,double));

    static void FLIP_Step(MAC* grid);

    static void UpdateSpaceHash();
    static void UpdatePressureMask();
    static void UpdateParticles();

    static void ParticleToGrid();
    static void GridToParticle();


};