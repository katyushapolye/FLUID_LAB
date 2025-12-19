/*Made by Raphael Laroca*/
#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <sstream>


#include "headers/MAC.h"
#include "headers/ADI.h"
#include "headers/PressureSolver.h"
#include "headers/Functions.h"
#include "headers/Definitions.h"
#include "headers/ConfigReader.h"
#include "headers/GridVisualizer.h"
#include "headers/Utils.h"

#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <implot.h>




void InitializeSimulation(const std::string& configFile = "simulation_config.txt");
void PrintProgressAndExport(int IT,int F_IT,double endTotal,double startTotal,double b4ProjectDiv,int &frame,float simTime);
void DrawSimulationPlots();



int main(int argc, char *argv[])
{
    // Initialize simulation
    if(argc == 2){
        std::string inputFile(argv[1]);
        InitializeSimulation(inputFile);
    }
    else{
        InitializeSimulation();
    }
    
    // Initialize GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }
    
    // GL 3.3 + GLSL 330
    const char* glsl_version = "#version 330";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    
    // Create window
    GLFWwindow* window = glfwCreateWindow(1600, 900, "3D ADI Solver Visualization", nullptr, nullptr);
    if (window == nullptr) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync
    
    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    
    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);
    
    // Setup style
    ImGui::StyleColorsDark();
    
    // Create visualizer
    GridVisualizer visualizer(SIMULATION.GRID_SOL);
    
    // Simulation variables
    int IT = 1;
    int frame = 1;
    double time = 0.0;
    double tF = 100000.0;
    int F_IT = tF / SIMULATION.dt;
    double start, end;
    double startTotal, endTotal;
    double b4Project = 0.0;
    
    bool simulationRunning = false;
    bool stepOnce = false;
    
    // Main loop
    while (!glfwWindowShouldClose(window) && time < tF)
    {

        glfwPollEvents();
        

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        

        ImGui::Begin("Simulation Control");
        ImGui::Text("Iteration: %d / %d", IT, F_IT);
        ImGui::Text("Time: %.4f / %.4f", time, tF);
        ImGui::Text("Divergence: %.6e", SIMULATION.GRID_SOL->GetDivSum());
        ImGui::Text("Frame Time: " "%.8f ms", (endTotal - startTotal) * 1000.0);
        if (simulationRunning) {
            ImGui::Text("Status: Running");
        } else {
            ImGui::Text("Status: Paused");
        }
        
        
        ImGui::Separator();
        
        if (ImGui::Button(simulationRunning ? "Pause" : "Run")) {
            simulationRunning = !simulationRunning;
        }
        ImGui::SameLine();
        if (ImGui::Button("Step")) {
            stepOnce = true;
        }
        ImGui::SameLine();
        ImGui::End();

        DrawSimulationPlots();
        

        if (simulationRunning || stepOnce) {
            startTotal = GetWallTime();
            
            // Diffusion + Convection
            ADI::SolveADIStep(SIMULATION.GRID_ANT, SIMULATION.GRID_SOL, time);
            
            // Pressure Projection
            if(SIMULATION.level == LevelConfiguration::STEP || 
               SIMULATION.level == LevelConfiguration::OBSTACLE) {
                SIMULATION.GRID_SOL->SetNeumannBorder();
            }
            
            PressureSolver::SolvePressure_AMGX(SIMULATION.GRID_SOL);
            b4Project = SIMULATION.GRID_SOL->GetDivSum();
            PressureSolver::ProjectPressure(SIMULATION.GRID_SOL);
            
            time += SIMULATION.dt;
            endTotal = GetWallTime();
            
            PrintProgressAndExport(IT, F_IT, endTotal, startTotal, b4Project, frame,time);
            
            if (SIMULATION.GRID_ANT->MaxAbsoluteDifference(*SIMULATION.GRID_SOL) < SIMULATION.TOLERANCE) {
                if(SIMULATION.level == LevelConfiguration::STEP || 
                   SIMULATION.level == LevelConfiguration::OBSTACLE) {
                    SIMULATION.GRID_SOL->SetNeumannBorder();
                }
                SIMULATION.GRID_SOL->ExportGrid(frame);
                frame++;
                simulationRunning = false; 
            }
            
            SIMULATION.GRID_ANT->CopyGrid(*SIMULATION.GRID_SOL);
            IT++;
            stepOnce = false;
            
            // Update visualizer with new grid data
            visualizer.UpdateGrid(SIMULATION.GRID_SOL);
        }
        
        // Render visualizer
        visualizer.Render();
        
        // Rendering
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        
        glfwSwapBuffers(window);
    }
    
    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    
    glfwDestroyWindow(window);
    glfwTerminate();
    
    return 0;
}

void InitializeSimulation(const std::string& configFile) {

    ConfigReader::loadConfig(SIMULATION, configFile);
    SIMULATION.GRID_ANT->ExportGrid(0);
        TELEMETRY.Push(
        0.0, 
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0
    );
    
    std::cout << "MAC Grid initialized - Parameters\n"
              << "-dt = " << std::to_string(SIMULATION.dt) << "\n"
              << "-dh = " << std::to_string(SIMULATION.dh) 
              << "\n-Nx = " << std::to_string(SIMULATION.Nx) 
              << "\n-Ny = " << std::to_string(SIMULATION.Ny)
              << "\n-Nz = " << std::to_string(SIMULATION.Nz) << std::endl;
    std::cout << "Total node count: " 
              << (SIMULATION.Nx * SIMULATION.Ny * SIMULATION.Nz) * 3 
              << " - Re = " << SIMULATION.RE 
              << " - Grid size: " << SIMULATION.GRID_SIZE << " - Tolerance: " << SIMULATION.TOLERANCE <<std::endl;
    
    ADI::InitializeADI(SIMULATION.GRID_SOL, SIMULATION.dt, SIMULATION.VelocityBoundaryFunction, ZERO, SIMULATION.PressureBoundaryFunction);
    PressureSolver::InitializePressureSolver(SIMULATION.GRID_SOL, SIMULATION.dt);
}

void DrawSimulationPlots()
{
    ImGui::Begin("Simulation Telemetry");
    
    float avail = ImGui::GetContentRegionAvail().x;
    float plot_w = avail * 0.5f - ImGui::GetStyle().ItemSpacing.x * 0.5f;
    ImVec2 plot_size(plot_w, 300);
    
    // ---- Row 1 ----
    if (ImPlot::BeginPlot("Divergence Evolution", plot_size)) {
        ImPlot::SetupAxes("Time", "Divergence Sum",
                          ImPlotAxisFlags_None,
                          ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        ImPlot::PlotLine("After Projection",
            TELEMETRY.time.data(), TELEMETRY.div_sum.data(), TELEMETRY.time.size());
        ImPlot::PlotLine("Before Projection",
            TELEMETRY.time.data(), TELEMETRY.div_sum_before_proj.data(), TELEMETRY.time.size());
        ImPlot::EndPlot();
    }
    
    ImGui::SameLine();
    
    if (ImPlot::BeginPlot("CFL Condition", plot_size)) {
        ImPlot::SetupAxes("Time", "CFL");
        ImPlot::PlotLine("CFL",
            TELEMETRY.time.data(), TELEMETRY.cfl.data(), TELEMETRY.time.size());
        double y = 1.0;
        ImPlot::PlotInfLines("CFL = 1", &y, 1);
        ImPlot::EndPlot();
    }
    
    // ---- Row 2 ----
    if (ImPlot::BeginPlot("Grid Residual", plot_size)) {

        // ---- Axes setup ----
        ImPlot::SetupAxes(
            "Time", "Residual",
            ImPlotAxisFlags_None,                       // X: manual
            ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_RangeFit  // Y: auto-fit
        );

        // ---- Scrolling time window ----
        const double t = TELEMETRY.time.back();
        const double window = 5.0; // seconds shown

        ImPlot::SetupAxisLimits(
            ImAxis_X1,
            t-window*0.1,
            t+window,
            ImGuiCond_Always
        );


        //ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

 
        ImPlot::PlotLine(
            "Iteration Residual",
            TELEMETRY.time.data(),
            TELEMETRY.residual.data(),
            TELEMETRY.time.size()
        );

        ImPlot::EndPlot();
    }
    ImGui::SameLine();

    
    if (ImPlot::BeginPlot("Performance (CPU vs GPU)", plot_size)) {
        double xs[] = {0.0, 1.0};
        static const char* labels[] = {"ADI (CPU)", "Pressure (GPU)"};
        double values[] = {
            TELEMETRY.cpu_time.back(),
            TELEMETRY.gpu_time.back()
        };
    
        static double y_max = 0.05;
        y_max = std::max(y_max, 1.2 * std::max(values[0], values[1]));
    
        ImPlot::SetupAxes(nullptr, "Seconds");
        ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, y_max, ImGuiCond_Always);
        ImPlot::SetupAxisTicks(ImAxis_X1, xs, 2, labels);
        ImVec4 colors[] = { ImVec4(0.4f, 0.7f, 0.0f, 1.0f),  // green for CPU
                            ImVec4(0.2f, 0.4f, 1.0f, 1.0f) }; // blue for GPU

        for (int i = 0; i < 2; i++) {
            ImPlot::PushStyleColor(ImPlotCol_Fill, colors[i]);
            ImPlot::PlotBars(labels[i], &xs[i], &values[i], 1, 0.5);
            ImPlot::PopStyleColor();
        }
        ImPlot::EndPlot();
    }
    
    ImGui::End();
    
}


void PrintProgressAndExport(int IT,int F_IT,double endTotal,double startTotal,double b4ProjectDiv,int &frame,float simTime){
    WriteToCSV( b4ProjectDiv,LevelConfigurationToString(SIMULATION.level),std::to_string( SIMULATION.GRID_SIZE),"DivergencyB4");
    WriteToCSV(SIMULATION.GRID_SOL->GetDivSum(),LevelConfigurationToString(SIMULATION.level),std::to_string( SIMULATION.GRID_SIZE),"Divergency");
    WriteToCSV(SIMULATION.lastADISolveTime,LevelConfigurationToString(SIMULATION.level),std::to_string( SIMULATION.GRID_SIZE),"TimeADI");
    WriteToCSV(SIMULATION.lastPressureSolveTime,LevelConfigurationToString(SIMULATION.level),std::to_string( SIMULATION.GRID_SIZE),"TimePressure");
    WriteToCSV(PressureSolver::GetSolverIterations(),LevelConfigurationToString(SIMULATION.level),std::to_string( SIMULATION.GRID_SIZE),"ItPressure");
    if(IT%1000 == 0){
        if(SIMULATION.level == LevelConfiguration::STEP || SIMULATION.level == LevelConfiguration::OBSTACLE) SIMULATION.GRID_SOL->SetNeumannBorder();
        SIMULATION.GRID_SOL->ExportGrid(frame);
        frame++;
    }
    printf("====================== ITERATION %d ==========================\n",IT);
    printf("It %d of %d                                \n", IT, F_IT);
    fflush(stdout);
    printf("Pressure Converged in %d it                \n", PressureSolver::GetSolverIterations());
    fflush(stdout);
    printf("Residual = %.10f                           \n", SIMULATION.GRID_ANT->MaxAbsoluteDifference(*SIMULATION.GRID_SOL));
    fflush(stdout);
    printf("It Time: %f s                              \n", endTotal - startTotal);
    fflush(stdout);
    printf("divsum (Bfr/Aft)= %.17f / %.17f            \n", b4ProjectDiv, SIMULATION.GRID_SOL->GetDivSum());
    fflush(stdout);
    printf("Adv.CFL: %f                                \n", (SIMULATION.GRID_SOL->GetMaxVelocity()/SIMULATION.dh)*SIMULATION.dt);
    fflush(stdout);
    printf("ADI Avrg Thread CPU Time: %f s             \n", (SIMULATION.lastADISolveTime));
    fflush(stdout);
    printf("Pressure Solve GPU Time: %f s              \n", (SIMULATION.lastPressureSolveTime));
    fflush(stdout);
    //printf("\rIt %d of %d -- Res = %.10f -- It Time: %f s - divsum (Bfr/Aft)= %.17f / %.17f - Adv.CFL: %f", IT, F_IT, SIMULATION.GRID_ANT->MaxAbsoluteDifference(*SIMULATION.GRID_SOL), endTotal - startTotal, b4Project,SIMULATION.GRID_SOL->GetDivSum(),(SIMULATION.GRID_SOL->GetMaxVelocity()/SIMULATION.dh)*SIMULATION.dt);
    printf("=============================================================\n",IT);
    fflush(nullptr);

    double cfl = (SIMULATION.GRID_SOL->GetMaxVelocity() / SIMULATION.dh) * SIMULATION.dt;
    double residual = SIMULATION.GRID_ANT->MaxAbsoluteDifference(*SIMULATION.GRID_SOL);

    TELEMETRY.Push(
        simTime, 
        SIMULATION.GRID_SOL->GetDivSum(),
        b4ProjectDiv,
        cfl,
        residual,
        SIMULATION.lastADISolveTime,
        SIMULATION.lastPressureSolveTime
    );


}


