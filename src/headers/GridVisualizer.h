#ifndef  GRID_VISUALIZER_H
#define GRID_VISUALIZER_H

#include "MAC.h"
#include <imgui.h>
#include <implot.h>
#include <vector>
#include <string>
#include <cmath>

class GridVisualizer {
private:
    MAC* grid;
    
    // UI state
    int selectedComponent = 1;  // 0=velocity, 1=u, 2=v, 3=w, 4=p, 5=solid
    int slicePlane = 1;         // 0=XY, 1=XZ, 2=YZ
    int sliceIndex = 0;
    bool autoScale = true;
    float minValue = 0.0f;
    float maxValue = 1.0f;
    
    // Quiver plot data
    std::vector<float> quiverX;
    std::vector<float> quiverY;
    std::vector<float> quiverU;
    std::vector<float> quiverV;
    
    // Quiver plot settings
    float baseSize = 12.0f;
    ImPlotQuiverFlags quiverFlags = ImPlotQuiverFlags_Colored | ImPlotQuiverFlags_Normalize;
    ImPlotColormap colormap = ImPlotColormap_Viridis;
    
    // Heatmap data (for scalar fields)
    std::vector<double> heatmapData;
    int heatmapRows = 0;
    int heatmapCols = 0;
    
    const char* componentNames[6] = {"Velocity Field", "U Velocity", "V Velocity", "W Velocity", "Pressure", "Solid Mask"};
    const char* planeNames[3] = {"XY Plane", "XZ Plane", "YZ Plane"};

    void ExtractSliceData();
    void ExtractQuiverData();
    void ComputeMinMax();

public:
    GridVisualizer(MAC* gridPtr);
    ~GridVisualizer();
    
    void Render();
    void UpdateGrid(MAC* newGrid);
};

// Implementation
GridVisualizer::GridVisualizer(MAC* gridPtr) : grid(gridPtr) {
    if (grid) {
        sliceIndex = grid->Ny / 2; // Start at middle slice
    }
}

GridVisualizer::~GridVisualizer() {}

void GridVisualizer::UpdateGrid(MAC* newGrid) {
    grid = newGrid;
}

void GridVisualizer::ComputeMinMax() {
    if (selectedComponent == 0) {
        // For velocity field, compute magnitude range
        minValue = 0.0f;
        maxValue = 0.0f;
        
        for (size_t i = 0; i < quiverU.size(); i++) {
            float mag = std::sqrt(quiverU[i] * quiverU[i] + quiverV[i] * quiverV[i]);
            if (mag > maxValue) maxValue = mag;
        }
        
        if (maxValue < 1e-10f) maxValue = 1e-10f;
    } else if (selectedComponent == 5) {
        // For solid mask, fixed range [0, 3]
        minValue = 0.0f;
        maxValue = 3.0f;
    } else if (heatmapData.empty()) {
        return;
    } else {
        minValue = heatmapData[0];
        maxValue = heatmapData[0];
        
        for (double val : heatmapData) {
            if (val < minValue) minValue = val;
            if (val > maxValue) maxValue = val;
        }
        
        if (std::abs(maxValue - minValue) < 1e-10) {
            maxValue = minValue + 1e-10;
        }
    }
}

void GridVisualizer::ExtractQuiverData() {
    if (!grid) return;
    
    quiverX.clear();
    quiverY.clear();
    quiverU.clear();
    quiverV.clear();
    
    int maxSlice = (slicePlane == 0) ? grid->Nz - 1 : 
                   (slicePlane == 1) ? grid->Ny - 1 : grid->Nx - 1;
    sliceIndex = std::max(0, std::min(sliceIndex, maxSlice));
    
    try {
        switch (slicePlane) {
            case 0: { // XY plane - show u and v velocities
                int k = std::min(sliceIndex, grid->Nz - 1);
                
                for (int j = 0; j < grid->Ny; j++) {
                    for (int i = 0; i < grid->Nx; i++) {
                        quiverX.push_back(i + 0.5f);
                        quiverY.push_back(j + 0.5f);
                        
                        // Interpolate u velocity to cell center
                        float u = 0.5f * (grid->GetU(j, i, k) + grid->GetU(j, i+1, k));
                        // Interpolate v velocity to cell center
                        float v = 0.5f * (grid->GetV(j, i, k) + grid->GetV(j+1, i, k));
                        
                        quiverU.push_back(u);
                        quiverV.push_back(v);
                    }
                }
                break;
            }
            
            case 1: { // XZ plane - show u and w velocities
                int j = std::min(sliceIndex, grid->Ny - 1);
                
                for (int k = 0; k < grid->Nz; k++) {
                    for (int i = 0; i < grid->Nx; i++) {
                        quiverX.push_back(i + 0.5f);
                        quiverY.push_back(k + 0.5f);
                        
                        // Interpolate u velocity to cell center
                        float u = 0.5f * (grid->GetU(j, i, k) + grid->GetU(j, i+1, k));
                        // Interpolate w velocity to cell center
                        float w = 0.5f * (grid->GetW(j, i, k) + grid->GetW(j, i, k+1));
                        
                        quiverU.push_back(u);
                        quiverV.push_back(w);
                    }
                }
                break;
            }
            
            case 2: { // YZ plane - show v and w velocities
                int i = std::min(sliceIndex, grid->Nx - 1);
                
                for (int k = 0; k < grid->Nz; k++) {
                    for (int j = 0; j < grid->Ny; j++) {
                        quiverX.push_back(j + 0.5f);
                        quiverY.push_back(k + 0.5f);
                        
                        // Interpolate v velocity to cell center
                        float v = 0.5f * (grid->GetV(j, i, k) + grid->GetV(j+1, i, k));
                        // Interpolate w velocity to cell center
                        float w = 0.5f * (grid->GetW(j, i, k) + grid->GetW(j, i, k+1));
                        
                        quiverU.push_back(v);
                        quiverV.push_back(w);
                    }
                }
                break;
            }
        }
    } catch (const std::exception& e) {
        std::cout << "[EXCEPTION] Error in ExtractQuiverData: " << e.what() << std::endl;
        return;
    }
    
    if (autoScale) {
        ComputeMinMax();
    }
}

void GridVisualizer::ExtractSliceData() {
    if (!grid) {
        std::cout << "[ERROR] Grid is null!" << std::endl;
        return;
    }
    
    heatmapData.clear();
    
    int maxSlice = (slicePlane == 0) ? grid->Nz - 1 : 
                   (slicePlane == 1) ? grid->Ny - 1 : grid->Nx - 1;
    sliceIndex = std::max(0, std::min(sliceIndex, maxSlice));
    
    try {
        switch (slicePlane) {
            case 0: { // XY plane
                int k = std::min(sliceIndex, grid->Nz - 1);
                
                if (selectedComponent == 1) { // U
                    heatmapRows = grid->Ny;
                    heatmapCols = grid->Nx + 1;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int y = 0; y < grid->Ny; y++) {
                        for (int x = 0; x <= grid->Nx; x++) {
                            int idx = y * heatmapCols + x;
                            heatmapData[idx] = grid->GetU(y, x, k);
                        }
                    }
                } else if (selectedComponent == 2) { // V
                    heatmapRows = grid->Ny + 1;
                    heatmapCols = grid->Nx;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int y = 0; y <= grid->Ny; y++) {
                        for (int x = 0; x < grid->Nx; x++) {
                            int idx = y * heatmapCols + x;
                            heatmapData[idx] = grid->GetV(y, x, k);
                        }
                    }
                } else if (selectedComponent == 3) { // W
                    k = std::min(sliceIndex, grid->Nz);
                    heatmapRows = grid->Ny;
                    heatmapCols = grid->Nx;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int y = 0; y < grid->Ny; y++) {
                        for (int x = 0; x < grid->Nx; x++) {
                            int idx = y * heatmapCols + x;
                            heatmapData[idx] = grid->GetW(y, x, k);
                        }
                    }
                } else if (selectedComponent == 4) { // Pressure
                    heatmapRows = grid->Ny;
                    heatmapCols = grid->Nx;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int y = 0; y < grid->Ny; y++) {
                        for (int x = 0; x < grid->Nx; x++) {
                            int idx = y * heatmapCols + x;
                            heatmapData[idx] = grid->GetP(y, x, k);
                        }
                    }
                } else if (selectedComponent == 5) { // Solid Mask
                    heatmapRows = grid->Ny;
                    heatmapCols = grid->Nx;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int y = 0; y < grid->Ny; y++) {
                        for (int x = 0; x < grid->Nx; x++) {
                            int idx = y * heatmapCols + x;
                            heatmapData[idx] = grid->GetSolid(y, x, k);
                        }
                    }
                }
                break;
            }
            
            case 1: { // XZ plane
                int y = std::min(sliceIndex, grid->Ny - 1);
                
                if (selectedComponent == 1) { // U
                    heatmapRows = grid->Nz;
                    heatmapCols = grid->Nx + 1;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int z = 0; z < grid->Nz; z++) {
                        for (int x = 0; x <= grid->Nx; x++) {
                            int idx = z * heatmapCols + x;
                            heatmapData[idx] = grid->GetU(y, x, z);
                        }
                    }
                } else if (selectedComponent == 2) { // V
                    y = std::min(sliceIndex, grid->Ny);
                    heatmapRows = grid->Nz;
                    heatmapCols = grid->Nx;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int z = 0; z < grid->Nz; z++) {
                        for (int x = 0; x < grid->Nx; x++) {
                            int idx = z * heatmapCols + x;
                            heatmapData[idx] = grid->GetV(y, x, z);
                        }
                    }
                } else if (selectedComponent == 3) { // W
                    heatmapRows = grid->Nz + 1;
                    heatmapCols = grid->Nx;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int z = 0; z <= grid->Nz; z++) {
                        for (int x = 0; x < grid->Nx; x++) {
                            int idx = z * heatmapCols + x;
                            heatmapData[idx] = grid->GetW(y, x, z);
                        }
                    }
                } else if (selectedComponent == 4) { // Pressure
                    heatmapRows = grid->Nz;
                    heatmapCols = grid->Nx;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int z = 0; z < grid->Nz; z++) {
                        for (int x = 0; x < grid->Nx; x++) {
                            int idx = z * heatmapCols + x;
                            heatmapData[idx] = grid->GetP(y, x, z);
                        }
                    }
                } else if (selectedComponent == 5) { // Solid Mask
                    heatmapRows = grid->Nz;
                    heatmapCols = grid->Nx;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int z = 0; z < grid->Nz; z++) {
                        for (int x = 0; x < grid->Nx; x++) {
                            int idx = z * heatmapCols + x;
                            heatmapData[idx] = grid->GetSolid(y, x, z);
                        }
                    }
                }
                break;
            }
            
            case 2: { // YZ plane
                int x = std::min(sliceIndex, grid->Nx - 1);
                
                if (selectedComponent == 1) { // U
                    x = std::min(sliceIndex, grid->Nx);
                    heatmapRows = grid->Nz;
                    heatmapCols = grid->Ny;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int z = 0; z < grid->Nz; z++) {
                        for (int y = 0; y < grid->Ny; y++) {
                            int idx = z * heatmapCols + y;
                            heatmapData[idx] = grid->GetU(y, x, z);
                        }
                    }
                } else if (selectedComponent == 2) { // V
                    heatmapRows = grid->Nz;
                    heatmapCols = grid->Ny + 1;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int z = 0; z < grid->Nz; z++) {
                        for (int y = 0; y <= grid->Ny; y++) {
                            int idx = z * heatmapCols + y;
                            heatmapData[idx] = grid->GetV(y, x, z);
                        }
                    }
                } else if (selectedComponent == 3) { // W
                    heatmapRows = grid->Nz + 1;
                    heatmapCols = grid->Ny;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int z = 0; z <= grid->Nz; z++) {
                        for (int y = 0; y < grid->Ny; y++) {
                            int idx = z * heatmapCols + y;
                            heatmapData[idx] = grid->GetW(y, x, z);
                        }
                    }
                } else if (selectedComponent == 4) { // Pressure
                    heatmapRows = grid->Nz;
                    heatmapCols = grid->Ny;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int z = 0; z < grid->Nz; z++) {
                        for (int y = 0; y < grid->Ny; y++) {
                            int idx = z * heatmapCols + y;
                            heatmapData[idx] = grid->GetP(y, x, z);
                        }
                    }
                } else if (selectedComponent == 5) { // Solid Mask
                    heatmapRows = grid->Nz;
                    heatmapCols = grid->Ny;
                    heatmapData.resize(heatmapRows * heatmapCols);
                    
                    for (int z = 0; z < grid->Nz; z++) {
                        for (int y = 0; y < grid->Ny; y++) {
                            int idx = z * heatmapCols + y;
                            heatmapData[idx] = grid->GetSolid(y, x, z);
                        }
                    }
                }
                break;
            }
        }
    } catch (const std::exception& e) {
        std::cout << "[EXCEPTION] Error in ExtractSliceData: " << e.what() << std::endl;
        return;
    }
    
    if (autoScale) {
        ComputeMinMax();
    }
}

void GridVisualizer::Render() {
    if (!grid) return;
    
    ImGui::Begin("Grid Visualizer", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
    
    if (ImGui::BeginTabBar("VisualizationTabs")) {
        if (ImGui::BeginTabItem("2D Slice")) {
            // Component selection
            ImGui::Text("Component:");
            ImGui::RadioButton("Velocity Field", &selectedComponent, 0); ImGui::SameLine();
            ImGui::RadioButton("U Velocity", &selectedComponent, 1); ImGui::SameLine();
            ImGui::RadioButton("V Velocity", &selectedComponent, 2); 
            ImGui::RadioButton("W Velocity", &selectedComponent, 3); ImGui::SameLine();
            ImGui::RadioButton("Pressure", &selectedComponent, 4);
            ImGui::RadioButton("Solid Mask", &selectedComponent, 5);
            
            ImGui::Separator();
            
            // Plane selection
            ImGui::Text("Slice Plane:");
            ImGui::RadioButton("XY (k)", &slicePlane, 0); ImGui::SameLine();
            ImGui::RadioButton("XZ (j)", &slicePlane, 1); ImGui::SameLine();
            ImGui::RadioButton("YZ (i)", &slicePlane, 2);
            
            // Slice index slider
            int maxSlice = (slicePlane == 0) ? grid->Nz - 1 : 
                           (slicePlane == 1) ? grid->Ny - 1 : grid->Nx - 1;
            ImGui::SliderInt("Slice Index", &sliceIndex, 0, maxSlice);
            
            ImGui::Separator();
            
            // Controls based on component type
            if (selectedComponent == 0) {
                // Quiver plot controls
                ImGui::Text("Quiver Plot Settings:");
                
                if (ImPlot::ColormapButton(ImPlot::GetColormapName(colormap), ImVec2(150, 0), colormap)) {
                    colormap = (ImPlotColormap)((colormap + 1) % ImPlot::GetColormapCount());
                }
                ImGui::SameLine();
                ImGui::Text("Colormap");
                
                ImGui::Checkbox("Auto Scale", &autoScale);
                if (!autoScale) {
                    ImGui::SetNextItemWidth(225);
                    ImGui::DragFloatRange2("Min / Max", &minValue, &maxValue, 0.01f, -20, 20, 
                                          nullptr, nullptr, ImGuiSliderFlags_AlwaysClamp);
                    if (maxValue <= minValue + 0.01f) {
                        maxValue = minValue + 0.01f;
                    }
                }
                
                ImGui::SetNextItemWidth(225);
                ImGui::DragFloat("Arrow Size", &baseSize, 0.1f, 0, 100);
                
                ImGui::CheckboxFlags("Normalize", (unsigned int*)&quiverFlags, ImPlotQuiverFlags_Normalize);
                ImGui::CheckboxFlags("Color Coded", (unsigned int*)&quiverFlags, ImPlotQuiverFlags_Colored);
            } else {
                // Heatmap controls
                if (selectedComponent == 5) {
                    // For solid mask, disable auto scale and fix range
                    ImGui::Text("Cell Types: 0=Fluid, 1=Solid, 2=Empty, 3=Inflow");
                    autoScale = false;
                    minValue = 0.0f;
                    maxValue = 3.0f;
                } else {
                    ImGui::Checkbox("Auto Scale", &autoScale);
                    if (!autoScale) {
                        ImGui::SliderFloat("Min Value", &minValue, -10.0f, 10.0f);
                        ImGui::SliderFloat("Max Value", &maxValue, -10.0f, 10.0f);
                    }
                }
            }
            
            // Extract data
            if (selectedComponent == 0) {
                ExtractQuiverData();
            } else {
                ExtractSliceData();
            }
            
            // Display info
            ImGui::Separator();
            ImGui::Text("Grid: %dx%dx%d", grid->Nx, grid->Ny, grid->Nz);
            if (selectedComponent == 0) {
                ImGui::Text("Vectors: %zu", quiverX.size());
            } else {
                ImGui::Text("Heatmap: %dx%d", heatmapCols, heatmapRows);
            }
            ImGui::Text("Range: [%.6f, %.6f]", minValue, maxValue);
            
            ImGui::EndTabItem();
        }
        
        if (ImGui::BeginTabItem("3D Volume")) {
            ImGui::Text("3D visualization coming soon...");
            ImGui::EndTabItem();
        }
        
        ImGui::EndTabBar();
    }
    
    ImGui::End();
    
    // Plot window
    ImGui::Begin("Plot", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
    
    // Use a discrete colormap for solid mask with exactly 4 colors
    if (selectedComponent == 5) {
        // Create custom discrete colormap for cell types
        static bool colormapCreated = false;
        static ImPlotColormap solidColormap;
        if (!colormapCreated) {
            ImVec4 colors[4] = {
                ImVec4(0.2f, 0.4f, 0.8f, 1.0f),  // 0 = FLUID_CELL  (blue)
                ImVec4(0.3f, 0.3f, 0.3f, 1.0f),  // 1 = SOLID_CELL  (gray)
                ImVec4(0.9f, 0.9f, 0.9f, 1.0f),  // 2 = EMPTY_CELL  (light gray)
                ImVec4(0.8f, 0.2f, 0.2f, 1.0f)   // 3 = INFLOW_CELL (red)
            };
            solidColormap = ImPlot::AddColormap("SolidMask", colors, 4);
            colormapCreated = true;
        }
        ImPlot::PushColormap(solidColormap);
    } else {
        ImPlot::PushColormap(colormap);
    }
    
float width = 1.0f, height = 1.0f;

switch(slicePlane) {
    case 0: // XY plane
        width  = (float)grid->Nx; // X axis
        height = (float)grid->Ny; // Y axis
        break;
    case 1: // XZ plane
        width  = (float)grid->Nx; // X axis
        height = (float)grid->Nz; // Z axis
        break;
    case 2: // YZ plane
        width  = (float)grid->Ny; // Y axis
        height = (float)grid->Nz; // Z axis
        break;
}

// Set plot size based on fixed height and correct aspect ratio
float plotHeight = 600.0f;
float plotWidth  = plotHeight * (width / height);
plotWidth = std::max(300.0f, std::min(plotWidth, 1200.0f));

// --- Quiver Plot ---
if (selectedComponent == 0 && !quiverX.empty()) {
    if (ImPlot::BeginPlot("##QuiverPlot", ImVec2(plotWidth, plotHeight))) {
        // Set axis limits according to physical dimensions
        ImPlot::SetupAxes(
            (slicePlane == 0 || slicePlane == 1) ? "X" : "Y",
            (slicePlane == 0) ? "Y" : "Z",
            ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit
        );
        ImPlot::SetupAxisLimits(ImAxis_X1, 0, width, ImGuiCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, 0, height, ImGuiCond_Always);

        ImPlot::SetNextQuiverStyle(baseSize, ImPlot::GetColormapColor(1));
        ImPlot::PlotQuiver("Velocity", 
                           quiverX.data(), quiverY.data(), 
                           quiverU.data(), quiverV.data(), 
                           quiverX.size(), 
                           minValue, maxValue, 
                           quiverFlags);
        ImPlot::EndPlot();
    }
}

// --- Heatmap Plot ---
else if (selectedComponent > 0 && !heatmapData.empty()) {
    if (ImPlot::BeginPlot("##Heatmap", ImVec2(plotWidth, plotHeight))) {
        ImPlot::SetupAxes(nullptr, nullptr, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_Invert);
        ImPlot::SetupAxisLimits(ImAxis_X1, 0, width, ImGuiCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, 0, height, ImGuiCond_Always);

        ImPlot::PlotHeatmap(componentNames[selectedComponent], 
                            heatmapData.data(), 
                            heatmapRows, 
                            heatmapCols,
                            minValue,
                            maxValue,
                            nullptr,
                            ImPlotPoint(0, 0),
                            ImPlotPoint(width, height));

        ImPlot::EndPlot();
    }

        
        ImGui::SameLine();
        ImPlot::ColormapScale("##HeatmapScale", minValue, maxValue, ImVec2(60, 600));
    }
    
    ImPlot::PopColormap();
    
    ImGui::End();
}

#endif // GRID_VISUALIZER_H