#include "../Core/Definitions.h"
#include <algorithm>
#include <omp.h>
#ifndef SPH2D_H
#define SPH2D_H

class SPH2D {
public:

    struct Particle
    {
        bool isActive = false;
        double x = 0.0;
        double y = 0.0;

        double p = 0.0;
        double d = 0.0;
        
        double fx = 0.0;
        double fy = 0.0;
    
        double u = 0.0;
        double v = 0.0;
    };

    static int particleCount;
    static int maxParticles;
    static Particle* particles;

    static double restDensity;
    static double gasConstant;
    static double viscosity;
    static double h;
    static double mass;
    static double collisionDamping;

    static std::vector<std::vector<std::vector<int>>> SPACE_HASH;

    static void InitializeSPH2D(int maxP);
    void static DummySPH(MAC* grid);
    void static AdvanceSPH(MAC* gridAnt,MAC* gridSol,double time);

    static double GetDensityAt(double x, double y);
    static double GetPressureAt(double x, double y);

    static void ExportParticles(int IT);

    static double GetMaxForce();
    static double GetMaxVelocity();

private:
    static int Nx, Ny;

    static void ComputeDensity();
    static void ComputePressure();
    static void ComputeForce();
    static void FirstIntegration(double dt);
    static void SecondIntegration(double dt);

    static void VelocityVerletIntegration(double dt);

    static double Wpoly6(double r, double h);
    static double WspikyGrad(double r, double h);
    static double WviscLaplacian(double r, double h);

    static void FindNeighbors(int pIdx, std::vector<int>& neighbors);
    static void FindNeighbors(double x,double y, std::vector<int>& neighbors);
    static void UpdateSpaceHash();
};



#endif