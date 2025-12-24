#include "headers/Solvers/PressureSolver_2D.h"

int PressureSolver2D::Nx = 0;
int PressureSolver2D::Ny = 0;
int PressureSolver2D::NON_ZERO = 0;
bool PressureSolver2D::needsFrameUpdate = true;
double PressureSolver2D::dt= 0;
VectorXd PressureSolver2D::IDP = VectorXd();

int *PressureSolver2D::collums = nullptr;
int *PressureSolver2D::rows = nullptr;
double *PressureSolver2D::values = nullptr;
CSRMatrix *PressureSolver2D::PRESSURE_MATRIX = nullptr;
AMGXSolver *PressureSolver2D::AMGX_Handle = nullptr;



std::vector<Tensor<double, 2>> PressureSolver2D::PRESSURE_MASK = std::vector<Tensor<double, 2>>();
SparseMatrix PressureSolver2D::PRESSURE_MATRIX_EIGEN;



void PressureSolver2D::InitializePressureSolver(MAC *grid, bool frameUpdate)
{
    needsFrameUpdate = frameUpdate;
    PressureSolver2D::Nx = SIMULATION.Nx;
    PressureSolver2D::Ny = SIMULATION.Ny;
    PressureSolver2D::PRESSURE_MASK.reserve(grid->GetFluidCellCount());
    PressureSolver2D::dt = SIMULATION.dt;
    double dh = SIMULATION.dh;
    PressureSolver2D::IDP = VectorXd(Nx * Ny);
    PressureSolver2D::IDP.setConstant(-1);  // FIX 1: Initialize to -1

    AMGX_initialize();

    int c = 0;
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            Tensor<double, 2> mask = Tensor<double, 2>(3, 3);
            mask.setZero();
            
            // FIX 2: ADD BOUNDARY CHECK (like 3D version)
            if (i == 0 || j == 0 || i == Ny-1 || j == Nx-1)
            {
                PressureSolver2D::PRESSURE_MASK.push_back(mask);
                SetIDP(i, j, -1);
                continue;
            }
            
            // Now check for solid/empty/inflow
            if (grid->GetSolid(i, j) == SOLID_CELL || 
                grid->GetSolid(i, j) == EMPTY_CELL || 
                grid->GetSolid(i, j) == INFLOW_CELL)
            {
                PressureSolver2D::PRESSURE_MASK.push_back(mask);
                SetIDP(i, j, -1);
                continue;
            }
            
            // FLUID CELL - build stencil
            NON_ZERO++; // diagonal entry
            
            // Check fluid neighbors (add off-diagonal entries)
            if (grid->GetSolid(i + 1, j) == FLUID_CELL)
            {
                mask(2, 1) = 1.0;
                mask(1, 1) -= 1.0;
                NON_ZERO++;
            }
            if (grid->GetSolid(i - 1, j) == FLUID_CELL)
            {
                mask(0, 1) = 1.0;
                mask(1, 1) -= 1.0;
                NON_ZERO++;
            }
            if (grid->GetSolid(i, j + 1) == FLUID_CELL)
            {
                mask(1, 2) = 1.0;
                mask(1, 1) -= 1.0;
                NON_ZERO++;
            }
            if (grid->GetSolid(i, j - 1) == FLUID_CELL)
            {
                mask(1, 0) = 1.0;
                mask(1, 1) -= 1.0;
                NON_ZERO++;
            }
            

            if (grid->GetSolid(i + 1, j) == EMPTY_CELL)
            {
                mask(1, 1) -= 1.0;
                
            }
            if (grid->GetSolid(i - 1, j) == EMPTY_CELL)
            {
                mask(1, 1) -= 1.0;
            }
            if (grid->GetSolid(i, j + 1) == EMPTY_CELL)
            {
                mask(1, 1) -= 1.0;
            }
            if (grid->GetSolid(i, j - 1) == EMPTY_CELL)
            {
                mask(1, 1) -= 1.0;
            }
            
            SetIDP(i, j, c);
            c++;
            PressureSolver2D::PRESSURE_MASK.push_back(mask);
        }
    }


    collums = (int *)malloc(sizeof(int) * NON_ZERO);
    rows = (int *)malloc(sizeof(int) * NON_ZERO);
    values = (double *)calloc(NON_ZERO, sizeof(double));
    int MatSize = grid->GetFluidCellCount();

    c = 0;
    for (int i = 1; i < Ny - 1; i++)
    {
        for (int j = 1; j < Nx - 1; j++)
        {
            Tensor<double, 2> mask = GetPressureMask(i, j);

            // Get row index once
            int MAT_LINE = GetIDP(i, j);
            
            // If this cell is not a FLUID cell, skip entirely
            if (MAT_LINE == -1)
                continue;

            // Diagonal (center)
            if (mask(1, 1) != 0.0)
            {
                rows[c] = MAT_LINE;
                collums[c] = MAT_LINE;
                values[c] += mask(1, 1);
                c++;
            }

            // +Y neighbor
            if (mask(2, 1) != 0.0)
            {
                int MAT_COLLUM = GetIDP(i + 1, j);
                if (MAT_COLLUM != -1)
                {
                    rows[c] = MAT_LINE;
                    collums[c] = MAT_COLLUM;
                    values[c] += mask(2, 1);
                    c++;
                }
            }

            // -Y neighbor
            if (mask(0, 1) != 0.0)
            {
                int MAT_COLLUM = GetIDP(i - 1, j);
                if (MAT_COLLUM != -1)
                {
                    rows[c] = MAT_LINE;
                    collums[c] = MAT_COLLUM;
                    values[c] += mask(0, 1);
                    c++;
                }
            }

            // +X neighbor
            if (mask(1, 2) != 0.0)
            {
                int MAT_COLLUM = GetIDP(i, j + 1);
                if (MAT_COLLUM != -1)
                {
                    rows[c] = MAT_LINE;
                    collums[c] = MAT_COLLUM;
                    values[c] += mask(1, 2);
                    c++;
                }
            }

            // -X neighbor
            if (mask(1, 0) != 0.0)
            {
                int MAT_COLLUM = GetIDP(i, j - 1);
                if (MAT_COLLUM != -1)
                {
                    rows[c] = MAT_LINE;
                    collums[c] = MAT_COLLUM;
                    values[c] += mask(1, 0);
                    c++;
                }
            }
        }
    }



    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < NON_ZERO; ++i)
    {
        triplets.push_back(Triplet(rows[i], collums[i], values[i]));
    }

    SparseMatrix mat(MatSize, MatSize);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    mat = ((1.0) / (dh * dh)) * mat; //we put the rho in the RHS

    PressureSolver2D::PRESSURE_MATRIX_EIGEN = mat;

    PRESSURE_MATRIX = coo_to_csr(rows, collums, values, NON_ZERO, MatSize, MatSize);

    free(collums);
    free(rows);
    free(values);

    PressureSolver2D::AMGX_Handle = new AMGXSolver();

    AMGX_config_create_from_file(&AMGX_Handle->config, "solver_pressure_config.txt");
    AMGX_resources_create_simple(&AMGX_Handle->rsrc, AMGX_Handle->config);

    AMGX_matrix_create(&AMGX_Handle->AmgxA, AMGX_Handle->rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&AMGX_Handle->Amgxb, AMGX_Handle->rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&AMGX_Handle->Amgxx, AMGX_Handle->rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&AMGX_Handle->Amgxguess, AMGX_Handle->rsrc, AMGX_mode_dDDI);

    for (int i = 0; i < NON_ZERO; i++)
    {
        PRESSURE_MATRIX->values[i] *= (1.0 / (dh * dh));
    }

    AMGX_matrix_upload_all(AMGX_Handle->AmgxA, MatSize, NON_ZERO, 1, 1, 
                           PRESSURE_MATRIX->row_ptr, PRESSURE_MATRIX->col_ind, 
                           PRESSURE_MATRIX->values, nullptr);
    AMGX_vector_set_zero(AMGX_Handle->Amgxx, MatSize, 1);
    AMGX_vector_set_zero(AMGX_Handle->Amgxb, MatSize, 1);
    AMGX_vector_set_zero(AMGX_Handle->Amgxguess, MatSize, 1);

    AMGX_solver_create(&AMGX_Handle->solver, AMGX_Handle->rsrc, AMGX_mode_dDDI, AMGX_Handle->config);
    AMGX_solver_setup(AMGX_Handle->solver, AMGX_Handle->AmgxA);


    std::cout << "Pressure Solver was initialized!\n";
}




void PressureSolver2D::UpdatePressureMatrix_Eigen(MAC* grid){
    PressureSolver2D::dt = dt;
    double dh = SIMULATION.dh;
    //reset the indexing vector
    PressureSolver2D::IDP = VectorXd(Nx * Ny);
    PressureSolver2D::IDP.setConstant(-1);
    PRESSURE_MASK.clear();
    NON_ZERO = 0;


    int c = 0;
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            Tensor<double, 2> mask = Tensor<double, 2>(3, 3);
            mask.setZero();
            
            // FIX 2: ADD BOUNDARY CHECK (like 3D version)
            if (i == 0 || j == 0 || i == Ny-1 || j == Nx-1)
            {
                PressureSolver2D::PRESSURE_MASK.push_back(mask);
                SetIDP(i, j, -1);
                continue;
            }
            
            // Now check for solid/empty/inflow
            if (grid->GetSolid(i, j) == SOLID_CELL || 
                grid->GetSolid(i, j) == EMPTY_CELL || 
                grid->GetSolid(i, j) == INFLOW_CELL)
            {
                PressureSolver2D::PRESSURE_MASK.push_back(mask);
                SetIDP(i, j, -1);
                continue;
            }
            
            // FLUID CELL - build stencil
            NON_ZERO++; // diagonal entry
            
            // Check fluid neighbors (add off-diagonal entries)
            if (grid->GetSolid(i + 1, j) == FLUID_CELL)
            {
                mask(2, 1) = 1.0;
                mask(1, 1) -= 1.0;
                NON_ZERO++;
            }
            if (grid->GetSolid(i - 1, j) == FLUID_CELL)
            {
                mask(0, 1) = 1.0;
                mask(1, 1) -= 1.0;
                NON_ZERO++;
            }
            if (grid->GetSolid(i, j + 1) == FLUID_CELL)
            {
                mask(1, 2) = 1.0;
                mask(1, 1) -= 1.0;
                NON_ZERO++;
            }
            if (grid->GetSolid(i, j - 1) == FLUID_CELL)
            {
                mask(1, 0) = 1.0;
                mask(1, 1) -= 1.0;
                NON_ZERO++;
            }
            

            if (grid->GetSolid(i + 1, j) == EMPTY_CELL)
            {
                mask(1, 1) -= 1.0;
                
            }
            if (grid->GetSolid(i - 1, j) == EMPTY_CELL)
            {
                mask(1, 1) -= 1.0;
            }
            if (grid->GetSolid(i, j + 1) == EMPTY_CELL)
            {
                mask(1, 1) -= 1.0;
            }
            if (grid->GetSolid(i, j - 1) == EMPTY_CELL)
            {
                mask(1, 1) -= 1.0;
            }
            
            SetIDP(i, j, c);
            c++;
            PressureSolver2D::PRESSURE_MASK.push_back(mask);
        }
    }


    collums = (int *)malloc(sizeof(int) * NON_ZERO);
    rows = (int *)malloc(sizeof(int) * NON_ZERO);
    values = (double *)calloc(NON_ZERO, sizeof(double));
    int MatSize = grid->GetFluidCellCount();
    c = 0;
    for (int i = 1; i < Ny - 1; i++)
    {
        for (int j = 1; j < Nx - 1; j++)
        {
            Tensor<double, 2> mask = GetPressureMask(i, j);

            // Get row index once
            int MAT_LINE = GetIDP(i, j);
            
            // If this cell is not a FLUID cell, skip entirely
            if (MAT_LINE == -1)
                continue;

            // Diagonal (center)
            if (mask(1, 1) != 0.0)
            {
                rows[c] = MAT_LINE;
                collums[c] = MAT_LINE;
                values[c] += mask(1, 1);
                c++;
            }

            // +Y neighbor
            if (mask(2, 1) != 0.0)
            {
                int MAT_COLLUM = GetIDP(i + 1, j);
                if (MAT_COLLUM != -1)
                {
                    rows[c] = MAT_LINE;
                    collums[c] = MAT_COLLUM;
                    values[c] += mask(2, 1);
                    c++;
                }
            }

            // -Y neighbor
            if (mask(0, 1) != 0.0)
            {
                int MAT_COLLUM = GetIDP(i - 1, j);
                if (MAT_COLLUM != -1)
                {
                    rows[c] = MAT_LINE;
                    collums[c] = MAT_COLLUM;
                    values[c] += mask(0, 1);
                    c++;
                }
            }

            // +X neighbor
            if (mask(1, 2) != 0.0)
            {
                int MAT_COLLUM = GetIDP(i, j + 1);
                if (MAT_COLLUM != -1)
                {
                    rows[c] = MAT_LINE;
                    collums[c] = MAT_COLLUM;
                    values[c] += mask(1, 2);
                    c++;
                }
            }

            // -X neighbor
            if (mask(1, 0) != 0.0)
            {
                int MAT_COLLUM = GetIDP(i, j - 1);
                if (MAT_COLLUM != -1)
                {
                    rows[c] = MAT_LINE;
                    collums[c] = MAT_COLLUM;
                    values[c] += mask(1, 0);
                    c++;
                }
            }
        }
    }



    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < NON_ZERO; ++i)
    {
        triplets.push_back(Triplet(rows[i], collums[i], values[i]));
    }

    SparseMatrix mat(MatSize, MatSize);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    mat = ((1.0) / (dh * dh*SIMULATION.RHO)) * mat; //we put the rho in the RHS

    PressureSolver2D::PRESSURE_MATRIX_EIGEN = mat;

    PRESSURE_MATRIX = coo_to_csr(rows, collums, values, NON_ZERO, MatSize, MatSize); //this allocates a new matrix so we can free our memory from the COO values

    free(collums);
    free(rows);
    free(values);
    for (int i = 0; i < NON_ZERO; i++)
    {
        PRESSURE_MATRIX->values[i] *= (1.0 / (dh * dh)); 
    }
    /*

    AMGX_matrix_create(&AMGX_Handle->AmgxA, AMGX_Handle->rsrc, AMGX_mode_dDDI);
    AMGX_solver_create(&AMGX_Handle->solver, AMGX_Handle->rsrc, AMGX_mode_dDDI, AMGX_Handle->config);


    AMGX_matrix_upload_all(AMGX_Handle->AmgxA, MatSize, NON_ZERO, 1, 1, 
                           PRESSURE_MATRIX->row_ptr, PRESSURE_MATRIX->col_ind, 
                           PRESSURE_MATRIX->values, nullptr);

    AMGX_vector_set_zero(AMGX_Handle->Amgxx, MatSize, 1);
    AMGX_vector_set_zero(AMGX_Handle->Amgxb, MatSize, 1);
    AMGX_solver_setup(AMGX_Handle->solver, AMGX_Handle->AmgxA);
    */

}


void PressureSolver2D::UpdatePressureMatrix_AMGX(MAC* grid){

    std::cout << "not working for now!" << std::endl;




    PressureSolver2D::dt = dt;
    double dh = SIMULATION.dh;
    //reset the indexing vector
    PressureSolver2D::IDP = VectorXd(Nx * Ny);
    PressureSolver2D::IDP.setConstant(-1);
    PRESSURE_MASK.clear();
    NON_ZERO = 0;


    int c = 0;
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            Tensor<double, 2> mask = Tensor<double, 2>(3, 3);
            mask.setZero();
            
            // FIX 2: ADD BOUNDARY CHECK (like 3D version)
            if (i == 0 || j == 0 || i == Ny-1 || j == Nx-1)
            {
                PressureSolver2D::PRESSURE_MASK.push_back(mask);
                SetIDP(i, j, -1);
                continue;
            }
            
            // Now check for solid/empty/inflow
            if (grid->GetSolid(i, j) == SOLID_CELL || 
                grid->GetSolid(i, j) == EMPTY_CELL || 
                grid->GetSolid(i, j) == INFLOW_CELL)
            {
                PressureSolver2D::PRESSURE_MASK.push_back(mask);
                SetIDP(i, j, -1);
                continue;
            }
            
            // FLUID CELL - build stencil
            NON_ZERO++; // diagonal entry
            
            // Check fluid neighbors (add off-diagonal entries)
            if (grid->GetSolid(i + 1, j) == FLUID_CELL)
            {
                mask(2, 1) = 1.0;
                mask(1, 1) -= 1.0;
                NON_ZERO++;
            }
            if (grid->GetSolid(i - 1, j) == FLUID_CELL)
            {
                mask(0, 1) = 1.0;
                mask(1, 1) -= 1.0;
                NON_ZERO++;
            }
            if (grid->GetSolid(i, j + 1) == FLUID_CELL)
            {
                mask(1, 2) = 1.0;
                mask(1, 1) -= 1.0;
                NON_ZERO++;
            }
            if (grid->GetSolid(i, j - 1) == FLUID_CELL)
            {
                mask(1, 0) = 1.0;
                mask(1, 1) -= 1.0;
                NON_ZERO++;
            }
            

            if (grid->GetSolid(i + 1, j) == EMPTY_CELL)
            {
                mask(1, 1) -= 1.0;
                
            }
            if (grid->GetSolid(i - 1, j) == EMPTY_CELL)
            {
                mask(1, 1) -= 1.0;
            }
            if (grid->GetSolid(i, j + 1) == EMPTY_CELL)
            {
                mask(1, 1) -= 1.0;
            }
            if (grid->GetSolid(i, j - 1) == EMPTY_CELL)
            {
                mask(1, 1) -= 1.0;
            }
            
            SetIDP(i, j, c);
            c++;
            PressureSolver2D::PRESSURE_MASK.push_back(mask);
        }
    }


    collums = (int *)malloc(sizeof(int) * NON_ZERO);
    rows = (int *)malloc(sizeof(int) * NON_ZERO);
    values = (double *)calloc(NON_ZERO, sizeof(double));
    int MatSize = grid->GetFluidCellCount();
    c = 0;
    for (int i = 1; i < Ny - 1; i++)
    {
        for (int j = 1; j < Nx - 1; j++)
        {
            Tensor<double, 2> mask = GetPressureMask(i, j);

            // Get row index once
            int MAT_LINE = GetIDP(i, j);
            
            // If this cell is not a FLUID cell, skip entirely
            if (MAT_LINE == -1)
                continue;

            // Diagonal (center)
            if (mask(1, 1) != 0.0)
            {
                rows[c] = MAT_LINE;
                collums[c] = MAT_LINE;
                values[c] += mask(1, 1);
                c++;
            }

            // +Y neighbor
            if (mask(2, 1) != 0.0)
            {
                int MAT_COLLUM = GetIDP(i + 1, j);
                if (MAT_COLLUM != -1)
                {
                    rows[c] = MAT_LINE;
                    collums[c] = MAT_COLLUM;
                    values[c] += mask(2, 1);
                    c++;
                }
            }

            // -Y neighbor
            if (mask(0, 1) != 0.0)
            {
                int MAT_COLLUM = GetIDP(i - 1, j);
                if (MAT_COLLUM != -1)
                {
                    rows[c] = MAT_LINE;
                    collums[c] = MAT_COLLUM;
                    values[c] += mask(0, 1);
                    c++;
                }
            }

            // +X neighbor
            if (mask(1, 2) != 0.0)
            {
                int MAT_COLLUM = GetIDP(i, j + 1);
                if (MAT_COLLUM != -1)
                {
                    rows[c] = MAT_LINE;
                    collums[c] = MAT_COLLUM;
                    values[c] += mask(1, 2);
                    c++;
                }
            }

            // -X neighbor
            if (mask(1, 0) != 0.0)
            {
                int MAT_COLLUM = GetIDP(i, j - 1);
                if (MAT_COLLUM != -1)
                {
                    rows[c] = MAT_LINE;
                    collums[c] = MAT_COLLUM;
                    values[c] += mask(1, 0);
                    c++;
                }
            }
        }
    }



    //destroy everythhing

    AMGX_SAFE_CALL(AMGX_solver_destroy(AMGX_Handle->solver));
    AMGX_SAFE_CALL(AMGX_vector_destroy(AMGX_Handle->Amgxx));
    AMGX_SAFE_CALL(AMGX_vector_destroy(AMGX_Handle->Amgxb));
    AMGX_SAFE_CALL(AMGX_matrix_destroy(AMGX_Handle->AmgxA));

    PRESSURE_MATRIX = coo_to_csr(rows, collums, values, NON_ZERO, MatSize, MatSize); //this allocates a new matrix so we can free our memory from the COO values

    free(collums);
    free(rows);
    free(values);


    AMGX_SAFE_CALL(AMGX_matrix_create(&AMGX_Handle->AmgxA, AMGX_Handle->rsrc, AMGX_mode_dDDI));
    AMGX_SAFE_CALL(AMGX_vector_create(&AMGX_Handle->Amgxb, AMGX_Handle->rsrc, AMGX_mode_dDDI));
    AMGX_SAFE_CALL(AMGX_vector_create(&AMGX_Handle->Amgxx, AMGX_Handle->rsrc, AMGX_mode_dDDI));
    AMGX_SAFE_CALL(AMGX_vector_create(&AMGX_Handle->Amgxguess, AMGX_Handle->rsrc, AMGX_mode_dDDI));

    for (int i = 0; i < NON_ZERO; i++)
    {
        PRESSURE_MATRIX->values[i] *= (1.0 / (dh * dh*SIMULATION.RHO));
    }

    AMGX_SAFE_CALL(AMGX_matrix_upload_all(AMGX_Handle->AmgxA, MatSize, NON_ZERO, 1, 1,PRESSURE_MATRIX->row_ptr, PRESSURE_MATRIX->col_ind, PRESSURE_MATRIX->values, nullptr));
    AMGX_SAFE_CALL(AMGX_vector_set_zero(AMGX_Handle->Amgxx, MatSize, 1));
    AMGX_SAFE_CALL(AMGX_vector_set_zero(AMGX_Handle->Amgxb, MatSize, 1));
    AMGX_SAFE_CALL(AMGX_vector_set_zero(AMGX_Handle->Amgxguess, MatSize, 1));

    AMGX_SAFE_CALL(AMGX_solver_create(&AMGX_Handle->solver, AMGX_Handle->rsrc, AMGX_mode_dDDI, AMGX_Handle->config));
    AMGX_SAFE_CALL(AMGX_solver_setup(AMGX_Handle->solver, AMGX_Handle->AmgxA));



 

}

void PressureSolver2D::SolvePressure_EIGEN(MAC* grid)
{


    if(needsFrameUpdate){
        UpdatePressureMatrix_Eigen(grid);
    }

    PressureSolver2D::dt = SIMULATION.dt;  
    CPUTimer timer;
    timer.start();
    double dh = grid->dh;
    int MatSize = grid->GetFluidCellCount();
    
    VectorXd RHS = VectorXd(MatSize);
    RHS.setZero();
    
    double mean = 0.0;

    // Build RHS vector using IDP indexing
    #pragma omp parallel
    {
    #pragma omp for
    for (int i = 1; i < Ny - 1; i++)
    {
        for (int j = 1; j < Nx - 1; j++)
        {
            int id = GetIDP(i, j);
            if (id == -1) continue;

            RHS(id) = (SIMULATION.RHO/dt) * grid->GetDivergencyAt(i, j);
            mean += RHS(id);

            if (grid->GetSolid(i, j-1) == INFLOW_CELL)
            {
                RHS(id) += -(SIMULATION.RHO/dt) * 
                    (grid->GetU(i, 1) - 
                     SIMULATION.VelocityBoundaryFunction2D(dh, i*dh + dh/2.0, 0).u);
            }
        }
    }
    }

    if (SIMULATION.NEEDS_COMPATIBILITY_CONDITION)  // FIX: Use SIMULATION
    {
        mean = mean / MatSize;
        for (int i = 0; i < MatSize; i++)
        {
            RHS(i) = RHS(i) - mean;
        }
    }

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver;
    solver.setMaxIterations(2000);
    solver.setTolerance(1e-12);
    solver.compute(PRESSURE_MATRIX_EIGEN);
    Eigen::VectorXd SOL = VectorXd(MatSize);
    SOL = solver.solve(RHS);

    // Write solution back using IDP indexing
    #pragma omp parallel
    {
    #pragma omp for
    for (int i = 1; i < Ny - 1; i++)
    {
        for (int j = 1; j < Nx - 1; j++)
        {
            int id = GetIDP(i, j);
            if (id == -1) continue;
            grid->SetP(i, j, SOL(id));
        }
    }
    }   

    grid->SetNeumannBorderPressure();
    SIMULATION.lastPressureSolveTime = timer.stop();  // FIX: Use SIMULATION
}


void PressureSolver2D::SolvePressure_AMGX(MAC* grid)
{
    if(needsFrameUpdate){
        UpdatePressureMatrix_AMGX(grid);
    }




    
    PressureSolver2D::dt = SIMULATION.dt;  
    CPUTimer timer;
    timer.start();
    
    double dh = grid->dh;
    int MatSize = grid->GetFluidCellCount();
    
    double* LHS = (double*)calloc(MatSize, sizeof(double));
    double* SOL = (double*)calloc(MatSize, sizeof(double));

    double mean = 0.0;

    #pragma omp parallel
    {
    #pragma omp for
    for (int i = 1; i < Ny - 1; i++)
    {
        for (int j = 1; j < Nx - 1; j++)
        {
            int id = GetIDP(i, j);
            if (id == -1) continue;

            LHS[id] = (SIMULATION.RHO/dt) * grid->GetDivergencyAt(i, j); //multiply by RHO here for density right??
            mean += LHS[id];

            if (grid->GetSolid(i, j-1) == INFLOW_CELL)
            {
                LHS[id] += -(SIMULATION.RHO/dt) * 
                    (grid->GetU(i, 1) - 
                     SIMULATION.VelocityBoundaryFunction2D(dh, i*dh + dh/2.0, 0).u);
            }
        }
    }
    }   

    if (SIMULATION.NEEDS_COMPATIBILITY_CONDITION)  
    {
        mean = mean / MatSize;
        for (int i = 0; i < MatSize; i++)
        {
            LHS[i] = LHS[i] - mean;
        }
    }

    // Solving
    AMGX_SAFE_CALL(AMGX_pin_memory(LHS,MatSize*sizeof(double)));

    AMGX_SAFE_CALL(AMGX_vector_upload(AMGX_Handle->Amgxb, MatSize, 1, LHS));



    AMGX_SAFE_CALL(AMGX_solver_solve(AMGX_Handle->solver, AMGX_Handle->Amgxb, AMGX_Handle->Amgxx));
    AMGX_SAFE_CALL(AMGX_vector_download(AMGX_Handle->Amgxx, SOL));
    AMGX_SAFE_CALL(AMGX_unpin_memory(LHS));

    // Write solution back using IDP indexing
    #pragma omp parallel
    {
    #pragma omp for
    for (int i = 1; i < Ny - 1; i++)
    {
        for (int j = 1; j < Nx - 1; j++)
        {
            int id = GetIDP(i, j);
            if (id == -1) continue;
            grid->SetP(i, j, SOL[id]);
        }
    }
    }

    grid->SetNeumannBorderPressure();
    SIMULATION.lastPressureSolveTime = timer.stop();  

    free(LHS);
    free(SOL);
    AMGX_SAFE_CALL(AMGX_vector_set_zero(AMGX_Handle->Amgxx, MatSize, 1));
    AMGX_SAFE_CALL(AMGX_vector_set_zero(AMGX_Handle->Amgxb, MatSize, 1));
}


void PressureSolver2D::ProjectPressure(MAC* grid)
{
    double dt = SIMULATION.dt;
    for (int i = 1; i < Ny - 1; i++)
    {
        for (int j = 1; j < Nx + 1 - 1; j++)
        {

            // if the cell we are is not a solid wall, then we update
            if (!(grid->GetSolid(i, j) == SOLID_CELL || grid->GetSolid(i, j - 1) == SOLID_CELL))
            {
          grid->SetU(i, j,
                     grid->GetU(i, j) - (dt/SIMULATION.RHO)*grid->GetGradPxAt(i, j));
            }
        }
    }

    for (int i = 1; i < Ny + 1 - 1; i++)
    {
        for (int j = 1; j < Nx - 1; j++)
        {

            if (!(grid->GetSolid(i, j) == SOLID_CELL || grid->GetSolid(i - 1, j) == SOLID_CELL))
            {
                grid->SetV(i, j,
                           grid->GetV(i, j) - (dt/SIMULATION.RHO)*grid->GetGradPyAt(i, j)

                );
            }
        }
    }

    // only update if it needs update
    /*
    for (int i = 0; i < Ny ; i++)
    {
        for (int j = 0; j < Nx + 1; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                if(grid->GetU_Update_Mask(i,j,k) == FLUID_CELL || grid->GetU_Update_Mask(i,j,k) == EMPTY_CELL){
                    grid->SetU(i,j,k,
                    grid->GetU(i,j,k) - (dt) * grid->GetGradPxAt(i,j,k)
                    );
                }
                else{
                continue;
                }

            }
        }
    }


    for (int i = 0; i < Ny+1; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            for (int k = 0; k < Nz; k++)
            {
                if(grid->GetV_Update_Mask(i,j,k) == FLUID_CELL || grid->GetV_Update_Mask(i,j,k) == EMPTY_CELL){
                    grid->SetV(i,j,k,
                    grid->GetV(i,j,k) - (dt) * grid->GetGradPyAt(i,j,k)

                );
            }


            }
        }
    }

    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            for (int k = 0; k < Nz+1; k++)
            {
                if(grid->GetW_Update_Mask(i,j,k) == FLUID_CELL || grid->GetW_Update_Mask(i,j,k) == EMPTY_CELL){
                    grid->SetW(i,j,k,
                    grid->GetW(i,j,k) - (dt) * grid->GetGradPzAt(i,j,k)

                );
            }


            }
        }
    }
        */

    if(SIM_TYPE == SIM_TYPES::FLIP){ //I am not exactly sure if this is necessary, but since inflip we have solid fluid empty interfaces (that shift a lot), we need to impose those a lot
        grid->SetBorder(SIMULATION.VelocityBoundaryFunction2D,SIMULATION.PressureBoundaryFunction2D,0);
    }
}




int PressureSolver2D::GetSolverIterations(){
    int iterations;

    return iterations;
}