#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cuda_runtime.h>

// Define structure to store parameters
struct LLCParams {
    std::string savename, datadir, savedir;
    // General Parameter
    double Lx, Ly;
    double tStart, tStop, tStepField, tStepParticle, tExpo;

    //For Fields
    int numGridX, numGridY, numGridBounX, numGridBounY;
    double dx, dy;
    int initCondFields;
    double hLDGa, hLDGb, hLDGk,Gamma, xi, xiAn, alpha, eta, h, zta, gammaV,convection;

    // For Particles
    int numParticles, blockNumParticles, threadNumParticles, cellNumX, cellNumY, initCondParticles, ABParticle;
    double cellSizeX, cellSizeY;
    double V0, boundarysize, kBT, TEffectiveBase;
    int maxParticlePerCell, maxParticlePerGrid, maxParticleNeighbor;
    double rd, minDistance, r0, epsilon, gammaB, Dc, Dp, Dparall, Dverti,rUpdateCellList;
    unsigned long long seed;
};

// Function to read the parameter file
void readParams(const std::string &filename, LLCParams &params) {
    std::ifstream file(filename);
    std::string line;

    // Read the parameter file line by line
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key;
        char eq;
        if (line.empty() || line[0] == '#') {
            // Skip empty or comment lines
            continue;
        }

        iss >> key >> eq;
        if (eq != '=') {
            std::cerr << "Invalid line format in parameter file: " << line << std::endl;
            continue;
        }

        // Match and store the parameter

        if (key == "savename") { iss >> params.savename; }
        else if (key == "datadir") { iss >> params.datadir; }
        else if (key == "savedir") { iss >> params.savedir; }
        else if (key == "Lx") { iss >> params.Lx; }
        else if (key == "Ly") { iss >> params.Ly; }
        else if (key == "convection") { iss >> params.convection; }
        else if (key == "tStart") { iss >> params.tStart; }
        else if (key == "tStop") { iss >> params.tStop; }
        else if (key == "tStepField") { iss >> params.tStepField; }
        else if (key == "tStepParticle") { iss >> params.tStepParticle; }
        else if (key == "tExpo") { iss >> params.tExpo; }
        else if (key == "numGridX") { iss >> params.numGridX; }
        else if (key == "numGridY") { iss >> params.numGridY; }
        else if (key == "numGridBounX") { iss >> params.numGridBounX; }
        else if (key == "numGridBounY") { iss >> params.numGridBounY; }
        else if (key == "dx") { iss >> params.dx; }
        else if (key == "dy") { iss >> params.dy; }
        else if (key == "initCondFields") { iss >> params.initCondFields; }
        else if (key == "hLDGa") { iss >> params.hLDGa; }
        else if (key == "hLDGb") { iss >> params.hLDGb; }
        else if (key == "hLDGk") { iss >> params.hLDGk; }
        else if (key == "Gamma") { iss >> params.Gamma; }
        else if (key == "xi") { iss >> params.xi; }
        else if (key == "xiAn") { iss >> params.xiAn; }
        else if (key == "alpha") { iss >> params.alpha; }
        else if (key == "eta") { iss >> params.eta; }
        else if (key == "h") { iss >> params.h; }
        else if (key == "zta") { iss >> params.zta; }
        else if (key == "gammaV") { iss >> params.gammaV; }
        else if (key == "numParticles") { iss >> params.numParticles; }
        else if (key == "blockNumParticles") { iss >> params.blockNumParticles; }
        else if (key == "threadNumParticles") { iss >> params.threadNumParticles; }
        else if (key == "cellNumX") { iss >> params.cellNumX; }
        else if (key == "cellNumY") { iss >> params.cellNumY; }
        else if (key == "cellSizeX") { iss >> params.cellSizeX; }
        else if (key == "cellSizeY") { iss >> params.cellSizeY; }
        else if (key == "initCondParticles") { iss >> params.initCondParticles; }
        else if (key == "ABParticle") { iss >> params.ABParticle; }
        else if (key == "V0") { iss >> params.V0; }
        else if (key == "boundarysize") { iss >> params.boundarysize; }
        else if (key == "kBT") { iss >> params.kBT; }
        else if (key == "TEffectiveBase") { iss >> params.TEffectiveBase; }
        else if (key == "maxParticlePerCell") { iss >> params.maxParticlePerCell; }
        else if (key == "maxParticlePerGrid") { iss >> params.maxParticlePerGrid; }
        else if (key == "maxParticleNeighbor") { iss >> params.maxParticleNeighbor; }
        else if (key == "rd") { iss >> params.rd; }
        else if (key == "minDistance") { iss >> params.minDistance; }
        else if (key == "r0") { iss >> params.r0; }
        else if (key == "epsilon") { iss >> params.epsilon; }
        else if (key == "gammaB") { iss >> params.gammaB; }
        else if (key == "Dc") { iss >> params.Dc; }
        else if (key == "Dp") { iss >> params.Dp; }
        else if (key == "Dparall") { iss >> params.Dparall; }
        else if (key == "Dverti") { iss >> params.Dverti; }
        else if (key == "seed") { iss >> params.seed; }
        else if (key == "rUpdateCellList") { iss >> params.rUpdateCellList; }
        // else if (key == "rUpdateGridList") { iss >> params.rUpdateGridList; }
    }
}

void printParams(const LLCParams &params) {
    std::cout << "savename: " << params.savename << std::endl;
    std::cout << "datadir: " << params.datadir << std::endl;
    std::cout << "savedir: " << params.savedir << std::endl;
    std::cout << "Lx: " << params.Lx << std::endl;
    std::cout << "Ly: " << params.Ly << std::endl;
    std::cout << "tStart: " << params.tStart << std::endl;
    std::cout << "tStop: " << params.tStop << std::endl;
    std::cout << "tStepField: " << params.tStepField << std::endl;
    std::cout << "tStepParticle: " << params.tStepParticle << std::endl;
    std::cout << "tExpo: " << params.tExpo << std::endl;
    std::cout << "numGridX: " << params.numGridX << std::endl;
    std::cout << "numGridY: " << params.numGridY << std::endl;
    std::cout << "numGridBounX: " << params.numGridBounX << std::endl;
    std::cout << "numGridBounY: " << params.numGridBounY << std::endl;
    std::cout << "dx: " << params.dx << std::endl;
    std::cout << "dy: " << params.dy << std::endl;
    std::cout << "initCondFields: " << params.initCondFields << std::endl;
    std::cout << "hLDGa: " << params.hLDGa << std::endl;
    std::cout << "hLDGb: " << params.hLDGb << std::endl;
    std::cout << "hLDGk: " << params.hLDGk << std::endl;
    std::cout << "xi: " << params.xi << std::endl;
    std::cout << "xiAn: " << params.xiAn << std::endl;
    std::cout << "alpha: " << params.alpha << std::endl;
    std::cout << "eta: " << params.eta << std::endl;
    std::cout << "h: " << params.h << std::endl;
    std::cout << "zta: " << params.zta << std::endl;
    std::cout << "gammaV: " << params.gammaV << std::endl;
    std::cout << "numParticles: " << params.numParticles << std::endl;
    std::cout << "blockNumParticles: " << params.blockNumParticles << std::endl;
}
// int main() {
//     // Initialize parameters structure
//     LLCParams params;
//     printf("1");
//     // Read the parameter file
//     readParams("./input.dat", params);

//     printf("%d",params.Lx);
    

//     return 0;
// }
