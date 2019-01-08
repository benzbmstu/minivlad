#include <vector>
#include <iostream>
#include <fstream>
using size_t = std::size_t;

class NetSOM {
public:
    struct NetConfig {
        size_t neurons;
        size_t inVecDim;
        int    trainEpochs;
        int    preTrainIterations;
        double minPotential;
        double deltaMinFuncEps;
        
        double sigmaInitPreTrain;
        double etaInitPreTrain;
        double sigmaInitPreTrainMin;
        double etaInitPreTrainMin;
        
        double sigmaInit;
        double etaInit;
        double sigmaInitMin;
        double etaInitMin;

        double weightLowBound;
        double weightUpperBound;

        std::string weightFileName;
    };

    void Init(const NetConfig &conf);
    void TrainKohen(std::vector<std::vector<double>> &trainSet, std::ostream &outputLabels);
    void TrainGas(std::vector<std::vector<double>> &trainSet, std::ostream &outputLabels);
    void PrintWeightsToFile(std::ostream &output, int precision);
    void CreateGnuplotAnimation(std::ofstream &stream);
    size_t DetectWinnerKohen(const std::vector<double> &inVec);
    size_t DetectWinnerGas(const std::vector<double> &inVec);

private:
    NetConfig netConfig;
    std::vector<double>   weights;
    std::vector<double>   potential;
    std::vector<int>      neuronsNeighbourSequence; // for DetectWinnerKohen
    std::vector<int>      victory;
    double sigma;
    double eta;
    void RandomWeights();
    void AdjustWeightsKohen(size_t winnerInd, const std::vector<double> &inVec);
    void AdjustWeightsGas(const std::vector<double> &inVec);
    void AdjustPotential(size_t winnerInd);
    size_t GetIndexInNeuronsNeighbourSequence(size_t neuronNumber);
};
