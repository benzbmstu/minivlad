#include "NetSOM.h"
#include <algorithm>
#include <cmath>
#include <random>

void GenerateElipsoid (
        std::vector<std::vector<double>> *points,
        double a, double b,
        double xOffset, double yOffset,
        double rotation, double step) 
{
    a = (a < 0) ? -a : a;
    b = (b < 0) ? -b : b;
    for (double x = -a; x <= a; x += step) {
        for (double y = -b; y <= b; y += step) {
             if (((x * x) / (a * a) + (y * y) / (b * b)) <= 1) {
                std::vector<double> v(2);
                v[0] = ( x * cos(rotation) + y * sin(rotation)) + xOffset;
                v[1] = (-x * sin(rotation) + y * cos(rotation)) + yOffset;
                points->emplace_back(std::move(v));
             }
        }
    }
}

void GenerateSpheres(std::vector<std::vector<double>> *points,
                     double x1, double y1, double r1,
                     double x2, double y2, double r2,
                     double step){

}
// ---------------------------------------------------------------------------------------
void SaveTrainSetToFile(std::ostream &output, const std::vector<std::vector<double>> &points) {
    for (size_t iPoint = 0; iPoint < points.size(); ++iPoint) {
        for (size_t iCoord = 0; iCoord < points[iPoint].size(); ++iCoord)
            output << points[iPoint][iCoord] << "\t";
        output << std::endl;
    }
}

int main() {
    std::vector<std::vector<double>> trainSet;
    GenerateElipsoid (&trainSet, 1., 0.5, 0., 0., 0., 0.05);
    NetSOM  net;
    NetSOM::NetConfig netConf;

    netConf.neurons  = 16;
    netConf.inVecDim = 2;
    netConf.trainEpochs = 10.;
    netConf.preTrainIterations = netConf.neurons;
    netConf.minPotential = 0.25;
    netConf.deltaMinFuncEps = 1e-3;

    netConf.sigmaInitPreTrain = 2.1;
    netConf.sigmaInitPreTrainMin = 0.5;
    netConf.etaInitPreTrain = 1.;
    netConf.etaInitPreTrainMin = 0.5;

    netConf.sigmaInit = 1.;
    netConf.sigmaInitMin = 1e-1;
    netConf.etaInit = 0.5;
    netConf.etaInitMin = 1e-3;
    netConf.weightLowBound   = -1.;
    netConf.weightUpperBound =  1.;

    netConf.weightFileName = std::string("weight.data");

    net.Init(netConf);

    std::ofstream trainSetFile("train.data", std::ios::trunc);
    std::ofstream weightFile(netConf.weightFileName, std::ios::trunc);
    std::ofstream gnuplotFile("plot.txt", std::ios::trunc);

    SaveTrainSetToFile(trainSetFile, trainSet);
    net.TrainGas(trainSet, weightFile);
    net.CreateGnuplotAnimation(gnuplotFile);
    std::cerr << "End of training " << std::endl;
}
