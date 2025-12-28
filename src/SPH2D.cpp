#include "headers/Solvers/SPH2D.h"

SPH2D::Particle* SPH2D::particles = nullptr;
int SPH2D::particleCount = 0;
int SPH2D::maxParticles = 0;
double SPH2D::restDensity = 1000;
double SPH2D::gasConstant = 700.0;
double SPH2D::viscosity = 0.5e-3;
double SPH2D::collisionDamping = 0.9;
double SPH2D::h = 0.04;
double SPH2D::mass = 0.01;
int SPH2D::Nx = 0;
int SPH2D::Ny = 0;
std::vector<std::vector<std::vector<int>>> SPH2D::SPACE_HASH;


void SPH2D::InitializeSPH2D(int maxP) {
    maxParticles = maxP;
    particles = new Particle[maxParticles];
    particleCount = 0;
    double dh = 0.005;
    double lengh= 0;

    for(float y = 0.0; y < 0.4;y+=dh){
         for(float x = 0.0; x < 0.4;x+=dh){
            particles[particleCount].x = x;
            particles[particleCount].y = y;
            particleCount++;
         }
        lengh+=dh;
    }

    SPH2D::Nx = (SIMULATION.domain.xf - SIMULATION.domain.x0)/h;
    SPH2D::Ny = (SIMULATION.domain.xf - SIMULATION.domain.x0)/h;

    for (int i = 0; i < Ny; i++)
    {
        SPACE_HASH.push_back(std::vector<std::vector<int>>());
        for (int j = 0; j <Nx; j++)
        {
            SPACE_HASH.at(i).push_back(std::vector<int>(SIMULATION.PARTICLES_PER_CELL));
        }
    }

    mass = restDensity * (dh*dh);

    UpdateSpaceHash();

}



double SPH2D::GetMaxForce(){
    double maxForce = sqrt(SPH2D::particles[0].fx * SPH2D::particles[0].fx + SPH2D::particles[0].fy * SPH2D::particles[0].fy);
    for (int p = 0; p < SPH2D::particleCount; p++) {
        double f = sqrt(SPH2D::particles[p].fx * SPH2D::particles[p].fx + SPH2D::particles[p].fy * SPH2D::particles[p].fy);

        maxForce = std::max(maxForce, f);
    }

    return maxForce;
}
double SPH2D::GetMaxVelocity(){
    double maxVelocity = sqrt(SPH2D::particles[0].u * SPH2D::particles[0].u + SPH2D::particles[0].v * SPH2D::particles[0].v);
    for (int p = 0; p < SPH2D::particleCount; p++) {
        double  v = sqrt(SPH2D::particles[p].u * SPH2D::particles[p].u + SPH2D::particles[p].v* SPH2D::particles[p].v);

        maxVelocity = std::max(maxVelocity, v);
    }

    return maxVelocity;
}

void SPH2D::DummySPH(MAC* grid){
}

void SPH2D::AdvanceSPH(MAC* gridAnt, MAC* gridSol, double time){


    

    
    //FirstIntegration(SIMULATION.dt);
    ComputeDensity();
    ComputePressure();
    ComputeForce();
    VelocityVerletIntegration(SIMULATION.dt);
    //SecondIntegration(SIMULATION.dt);
    UpdateSpaceHash();
}


void SPH2D::UpdateSpaceHash()
{
    // Parallel cleaning of hash grid
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            SPACE_HASH.at(i).at(j).clear();
        }
    }

    // Serial insertion to avoid race conditions
    // Alternative: use thread-local storage or locks, but overhead may not be worth it
    for (int p = 0; p < particleCount; p++)
    {
        int i = std::min((int)(particles[p].y / h), Ny - 1);
        int j = std::min((int)(particles[p].x / h), Nx - 1);
        
        i = std::max(0, i);
        j = std::max(0, j);

        SPACE_HASH.at(i).at(j).push_back(p);
    }
}

void SPH2D::FindNeighbors(int pIdx, std::vector<int>& neighbors) {
    neighbors.clear();
    double dh = h;
    double hSquared = h * h;

    int iCenter = particles[pIdx].y / dh;
    int jCenter = particles[pIdx].x / dh;

    for (int di = -1; di <= 1; di++) {
        for (int dj = -1; dj <= 1; dj++) {
            int i = iCenter + di;
            int j = jCenter + dj;

            if (i >= 0 && i < Ny && j >= 0 && j < Nx) {
                for (int neighborIdx : SPACE_HASH[i][j]) {
                    if (neighborIdx != pIdx) {
                        double dx = particles[neighborIdx].x - particles[pIdx].x;
                        double dy = particles[neighborIdx].y - particles[pIdx].y;
                        double distSquared = dx * dx + dy * dy;

                        if (distSquared < hSquared) {
                            neighbors.push_back(neighborIdx);
                        }
                    }
                }
            }
        }
    }
}

void SPH2D::FindNeighbors(double x,double y, std::vector<int>& neighbors) {
    neighbors.clear();
    double dh = h;
    double hSquared = h * h;

    int iCenter = y / dh;
    int jCenter = x / dh;

    for (int di = -1; di <= 1; di++) {
        for (int dj = -1; dj <= 1; dj++) {
            int i = iCenter + di;
            int j = jCenter + dj;

            if (i >= 0 && i < Ny && j >= 0 && j < Nx) {
                for (int neighborIdx : SPACE_HASH[i][j]) {
                    double dx = particles[neighborIdx].x - x;
                    double dy = particles[neighborIdx].y - y;
                    double distSquared = dx * dx + dy * dy;

                    if (distSquared < hSquared) {
                        neighbors.push_back(neighborIdx);
                    }
                }
            }
        }
    }
}


double SPH2D::Wpoly6(double r, double h) {
    if (r < 0 || r > h) return 0.0;
    double coeff = 4.0 / (M_PI * pow(h, 8));
    return coeff * pow(h * h - r * r, 3);
}

double SPH2D::WspikyGrad(double r, double h) {
    if (r < 1e-6 || r > h) return 0.0;
    double coeff = -30.0 / (M_PI * pow(h, 5));
    return coeff * pow(h - r, 2);
}

double SPH2D::WviscLaplacian(double r, double h) {
    if (r > h) return 0.0;
    double coeff = 20.0 / (3.0 * M_PI * pow(h, 5));
    return coeff * (h - r);
}

void SPH2D::ComputeDensity(){
    #pragma omp parallel for schedule(dynamic, 64)
    for (int p = 0; p < particleCount; p++) {
        double density = mass * Wpoly6(0, h);
        std::vector<int> neighbors;
        FindNeighbors(p, neighbors);

        for (int neighborIdx : neighbors) {
            double dx = particles[p].x - particles[neighborIdx].x;
            double dy = particles[p].y - particles[neighborIdx].y;
            double r = sqrt(dx * dx + dy * dy);
            density += mass * Wpoly6(r, h);
        }
        particles[p].d = density;
    }
}


void SPH2D::ComputePressure() {
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < particleCount; p++) {
        double gamma = 7.0;
        double B = (restDensity * gasConstant) / gamma;
        
        double ratio = particles[p].d / restDensity;
        particles[p].p = B * (pow(ratio, gamma) - 1.0);

        if (particles[p].p < 0) particles[p].p = 0;
    }
}

void SPH2D::ComputeForce(){
    #pragma omp parallel for schedule(dynamic, 64)
    for (int p = 0; p < particleCount; p++) {
        double Fx = 0.0;
        double Fy = -9.8f * mass;

        std::vector<int> neighbors;
        FindNeighbors(p, neighbors);

        for (int neighborIdx : neighbors) {
            double dx = particles[p].x - particles[neighborIdx].x;
            double dy = particles[p].y - particles[neighborIdx].y;
            double r = sqrt(dx * dx + dy * dy);
            
            if (r < 1e-6) continue;

            double gradMag = WspikyGrad(r, h);

            double pressureTerm = (particles[p].p / (particles[p].d * particles[p].d)) + 
                                 (particles[neighborIdx].p / (particles[neighborIdx].d * particles[neighborIdx].d));

            double s_corr = 0.0;
            double r_norm = r / h;
            if (r_norm < 0.15) { // Prevents cumpling, Morghan et all or something
                s_corr = 0.01 * pow(Wpoly6(r, h) / Wpoly6(h * 0.2, h), 4);
            }
            // -------------------------------------------------
            
            double f_press = -mass * mass * (pressureTerm + s_corr) * gradMag;
            Fx += f_press * (dx / r);
            Fy += f_press * (dy / r);

            

            double viscMag = WviscLaplacian(r, h);
            Fx += viscosity * mass * (particles[neighborIdx].u - particles[p].u) / particles[neighborIdx].d * viscMag;
            Fy += viscosity * mass * (particles[neighborIdx].v - particles[p].v) / particles[neighborIdx].d * viscMag;    
        }

        particles[p].fx = Fx;
        particles[p].fy = Fy;
    }
}


void SPH2D::VelocityVerletIntegration(double dt) {
    // Store old accelerations
    std::vector<double> ax_old(particleCount), ay_old(particleCount);
    
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < particleCount; p++) {
        ax_old[p] = particles[p].fx / mass;
        ay_old[p] = particles[p].fy / mass;
        
        // Update positions
        particles[p].x += particles[p].u * dt + 0.5 * ax_old[p] * dt * dt;
        particles[p].y += particles[p].v * dt + 0.5 * ay_old[p] * dt * dt;
        
        // Update velocities with half acceleration
        particles[p].u += 0.5 * ax_old[p] * dt;
        particles[p].v += 0.5 * ay_old[p] * dt;
    
    
        double d_min = h * 0.1;
        
        double distBottom = particles[p].y - SIMULATION.domain.y0;
        if (distBottom < d_min) {

            particles[p].v *= -collisionDamping;
            particles[p].u *= collisionDamping; 
        }
        
        double distTop = SIMULATION.domain.yf - particles[p].y;
        if (distTop < d_min) {

            particles[p].v *= -collisionDamping;
            particles[p].u *= collisionDamping; 
        }
        
        double distLeft = particles[p].x - SIMULATION.domain.x0;
        if (distLeft < d_min) {

            particles[p].v *= collisionDamping;
            particles[p].u *= -collisionDamping; 
        }
        
        double distRight = SIMULATION.domain.xf - particles[p].x;
        if (distRight < d_min) {

            particles[p].v *= collisionDamping;
            particles[p].u *= -collisionDamping; 
        }
        
        // Hard clamp
        if (particles[p].x < SIMULATION.domain.x0) {
            particles[p].x = SIMULATION.domain.x0;

        }
        if (particles[p].x > SIMULATION.domain.xf) {
            particles[p].x = SIMULATION.domain.xf;

        }
        if (particles[p].y < SIMULATION.domain.y0) {
            particles[p].y = SIMULATION.domain.y0;

        }
        if (particles[p].y > SIMULATION.domain.yf) {
            particles[p].y = SIMULATION.domain.yf;

        }
    }

    

    UpdateSpaceHash();
    ComputeDensity();
    ComputePressure();
    ComputeForce();

    #pragma omp parallel for schedule(static)
    for (int p = 0; p < particleCount; p++) {
        double ax_new = particles[p].fx / mass;
        double ay_new = particles[p].fy / mass;
        
        particles[p].u += 0.5 * ax_new * dt;
        particles[p].v += 0.5 * ay_new * dt;
    }
}


void SPH2D::FirstIntegration(double dt) {
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < particleCount; p++) {
        // Half-step velocity
        particles[p].u += (particles[p].fx / mass) * 0.5 * dt;
        particles[p].v += (particles[p].fy / mass) * 0.5 * dt;
        
        // Full-step position
        particles[p].x += particles[p].u * dt;
        particles[p].y += particles[p].v * dt;
        


        double d_min = h * 0.1;
        
        double distBottom = particles[p].y - SIMULATION.domain.y0;
        if (distBottom < d_min) {

            particles[p].v *= -collisionDamping;
            particles[p].u *= collisionDamping; 
        }
        
        double distTop = SIMULATION.domain.yf - particles[p].y;
        if (distTop < d_min) {

            particles[p].v *= -collisionDamping;
            particles[p].u *= collisionDamping; 
        }
        
        double distLeft = particles[p].x - SIMULATION.domain.x0;
        if (distLeft < d_min) {

            particles[p].v *= collisionDamping;
            particles[p].u *= -collisionDamping; 
        }
        
        double distRight = SIMULATION.domain.xf - particles[p].x;
        if (distRight < d_min) {

            particles[p].v *= collisionDamping;
            particles[p].u *= -collisionDamping; 
        }
        
        // Hard clamp
        if (particles[p].x < SIMULATION.domain.x0) {
            particles[p].x = SIMULATION.domain.x0;

        }
        if (particles[p].x > SIMULATION.domain.xf) {
            particles[p].x = SIMULATION.domain.xf;

        }
        if (particles[p].y < SIMULATION.domain.y0) {
            particles[p].y = SIMULATION.domain.y0;

        }
        if (particles[p].y > SIMULATION.domain.yf) {
            particles[p].y = SIMULATION.domain.yf;

        }
    }
}


void SPH2D::SecondIntegration(double dt){
    #pragma omp parallel for schedule(static)
    for (int p = 0; p < particleCount; p++)
    {
        particles[p].u += (particles[p].fx/mass) * 0.5 * dt;
        particles[p].v += (particles[p].fy/mass) * 0.5 * dt;
    }
}


double SPH2D::GetDensityAt(double x, double y) {
    double density = 0.0;
    std::vector<int> neighbors;

    FindNeighbors(x,y, neighbors);
    
    for (int idx = 0; idx < neighbors.size(); idx++) {
        double dx = x - particles[neighbors[idx]].x;
        double dy = y - particles[neighbors[idx]].y;
        double r = sqrt(dx * dx + dy * dy);
        density += mass * Wpoly6(r, h);
    }
    return density;
}

double SPH2D::GetPressureAt(double x, double y) {
    double pressureSum = 0.0;
    double weightSum = 0.0;
    std::vector<int> neighbors;
    FindNeighbors(x, y, neighbors);
    
    // Fixed: idx is already the particle index, not an index into neighbors
    for (int idx : neighbors) {
        double dx = x - particles[idx].x;  // Changed from neighbors[idx]
        double dy = y - particles[idx].y;  // Changed from neighbors[idx]
        double r = sqrt(dx * dx + dy * dy);
        double weight = Wpoly6(r, h);
        pressureSum += particles[idx].p * weight;  // Changed from neighbors[idx]
        weightSum += weight;
    }
    
    // Normalize by total weight to avoid bias
    if (weightSum > 1e-6) {
        return pressureSum / weightSum;
    }
    return 0.0;
}

void SPH2D::ExportParticles(int IT)
{
    std::string exportPath = "Exports/SPH/particles_" + std::to_string(IT) + ".csv";
    std::ofstream outFile(exportPath);
    if(!outFile.is_open()){
        std::cout << "SPH2D EXPORT FAILED!" << std::endl;
    }


    // Write header
    outFile << "x,y,u,v,fx,fy,d,p,m" << std::endl;

    // Write data for each active particle
    for (int i = 0; i < particleCount; i++)
    {
        if (particles[i].isActive)
        {
            outFile << particles[i].x << ","
                    << particles[i].y << ","
                    << particles[i].u << ","
                    << particles[i].v << ","
                    << particles[i].fx << ","
                    << particles[i].fy << ","
                    << particles[i].d << ","
                    << particles[i].p << ","
                    << mass << std::endl;
        }
    }

    outFile.close();
    // std::cout << "Exported " << particleCount << " particles to particles.csv" << std::endl;
}