#include "NetSOM.h"
#include <iomanip>
#include <random>
#include <algorithm>
#include <limits>
#include <fstream>
#include <unistd.h>

void NetSOM::Init(const NetConfig &conf) {
    netConfig = conf;
    weights.resize(netConfig.neurons * netConfig.inVecDim);
    potential.resize(netConfig.neurons, netConfig.minPotential);
    neuronsNeighbourSequence.resize(netConfig.neurons, 0);
    victory.resize(netConfig.neurons, 0);
    RandomWeights();
}
bool compare(const std::pair<double, double>&i, const std::pair<double, double>&j){ 
    return i.first < j.first; 
}
// ------------------------------------------------------------------------------------------------
void NetSOM::RandomWeights() {
    std::random_device randomizer;
    std::mt19937 randGen(randomizer());
    std::uniform_real_distribution<> dist(netConfig.weightLowBound, 
                                          netConfig.weightUpperBound);
    double weightNorm;
    double weight = 0.0;
    for (size_t iNeuron = 0; iNeuron < netConfig.neurons; ++iNeuron) {
        for (size_t iWeight = 0; iWeight < netConfig.inVecDim; ++iWeight) {
            weight = dist(randGen);
            std::cerr << weight << std::endl;
            weights[iWeight + iNeuron * netConfig.inVecDim] = weight;
        }
    }
}
// ------------------------------------------------------------------------------------------------
// returns winner index
size_t NetSOM::DetectWinnerKohen(const std::vector<double> &inVec) {
    std::vector<double> distancesToInput(netConfig.neurons, 0.0);

    
    for (size_t iNeuron = 0; iNeuron < netConfig.neurons; ++iNeuron) {
        for (size_t iWeight = 0; iWeight < netConfig.inVecDim; ++iWeight) {
            distancesToInput[iNeuron] += pow(inVec[iWeight] - weights[iWeight + iNeuron * netConfig.inVecDim], 2);        
        }
        if (potential[iNeuron] < netConfig.minPotential) {
            distancesToInput[iNeuron] = std::numeric_limits<double>::max();
        } else {
            distancesToInput[iNeuron] = potential[iNeuron] * sqrt(distancesToInput[iNeuron]);
        }
    }
    size_t winnerInd = std::distance(distancesToInput.begin(), std::min_element(distancesToInput.begin(), distancesToInput.end()));
    return winnerInd;
}
//////////////////////////////////////////////////////////////////////////////////
size_t NetSOM::DetectWinnerGas(const std::vector<double> &inVec) {

    std::vector<std::pair<double, int>>  distancesToInput (netConfig.neurons, std::pair<double, int> (0.0, 0));

    // calculate distance for each neuron
    for (size_t iNeuron = 0; iNeuron < netConfig.neurons; ++iNeuron) {
        for (size_t iWeight = 0; iWeight < netConfig.inVecDim; ++iWeight) 
            distancesToInput[iNeuron].first += pow(inVec[iWeight] - weights[iWeight + iNeuron * netConfig.inVecDim], 2);

        if (potential[iNeuron] < netConfig.minPotential)
            distancesToInput[iNeuron].first = std::numeric_limits<double>::max();
        else 
            distancesToInput[iNeuron].first = potential[iNeuron] * sqrt(distancesToInput[iNeuron].first);
        
        //distancesToInput[iNeuron].first  = sqrt(distancesToInput[iNeuron].first);
        //distancesToInput[iNeuron].first = potential[iNeuron] * sqrt(distancesToInput[iNeuron].first);
        distancesToInput[iNeuron].second = iNeuron;
    }
    //print unsorted array
    /*std::cout << "-UNSORTED distances" << std::endl;
    for (size_t i = 0; i < distancesToInput.size(); ++i)
      std::cout << distancesToInput[i].first << " " << distancesToInput[i].second << std::endl;
    */
    // sort distances
    std::sort(distancesToInput.begin(), distancesToInput.end(), compare);

    // print sorted array
    std::cout << "+SORTED distances" << std::endl;
    for (size_t i = 0; i < distancesToInput.size(); ++i)
      std::cout << distancesToInput[i].first << " " << distancesToInput[i].second << std::endl;

    //update sequence neuron after sorting
    std::cout << "Neightbours sequence is";
    for (size_t i = 0; i < distancesToInput.size(); ++i){
      neuronsNeighbourSequence[i] = distancesToInput[i].second;
      std::cout << ' ' << neuronsNeighbourSequence[i];
    }
    std::cout << std::endl;

    size_t winnerInd = distancesToInput[0].second;
    std::cout << "winnerInd is " << winnerInd << std::endl;

    // increase wonner counter in vector victory.
    victory[winnerInd] ++;
    std::cout << "Счетчик побед: ";
    for (int &a: victory)
      std::cout << ' ' << a;
    std::cout << '\n';

    
    return winnerInd;
}
// -------------------------------------------------------------------------------------
size_t NetSOM::GetIndexInNeuronsNeighbourSequence(size_t neuronNumber){
    std::vector<int>::iterator it = std::find(neuronsNeighbourSequence.begin(), 
                                              neuronsNeighbourSequence.end(), 
                                              neuronNumber);
    /*if ! (it != neuronsNeighbourSequence.end())
        throw std::string("NetSOM::GetIndexInNeuronsNeighbourSequence --> Element Not Found!");*/

    // Get index of element from iterator
    int index = std::distance(neuronsNeighbourSequence.begin(), it);
    return index;
}
// ------------------------------------------------------------------------------------------------
void NetSOM::AdjustWeightsKohen(size_t winnerInd, const std::vector<double> &inVec) {
    for (size_t iNeuron = 0; iNeuron < netConfig.neurons; ++iNeuron) {
        for (size_t iWeight = 0; iWeight < netConfig.inVecDim; ++iWeight) {
	          double d = fabs(iNeuron - winnerInd);
	          double neigbourCoeff = exp(-pow(d / sigma, 2) / 2);
            weights[iWeight + iNeuron * netConfig.inVecDim] += eta * neigbourCoeff 
              * (inVec[iWeight] - weights[iWeight + iNeuron * netConfig.inVecDim]);
        }
    }
}

// ------------------------------------------------------------------------------------------------
void NetSOM::AdjustWeightsGas(const std::vector<double> &inVec) {
    for (size_t iNeuron = 0; iNeuron < netConfig.neurons; ++iNeuron) {
        int m = GetIndexInNeuronsNeighbourSequence(iNeuron);
        double neigbourCoeff = exp(- m / sigma);
        for (size_t iWeight = 0; iWeight < netConfig.inVecDim; ++iWeight) {
            weights[iWeight + iNeuron * netConfig.inVecDim] += eta * neigbourCoeff 
              * (inVec[iWeight] - weights[iWeight + iNeuron * netConfig.inVecDim]);
        }
    }
}
// ------------------------------------------------------------------------------------------------
void NetSOM::AdjustPotential(size_t winnerInd) {
    std::cout << std::endl << "Список потенциалов:" << std::endl;
    for (size_t iPotential = 0; iPotential < potential.size(); iPotential++) {
        (iPotential == winnerInd) ? potential[iPotential] -= netConfig.minPotential
                                  : potential[iPotential] += 1. / netConfig.neurons;
        std::cout << potential[iPotential] << std::endl;
    }
    std::cout << std::endl;
}
void NetSOM::AdjustPotentialInverted(size_t winnerInd) {
    std::cout << std::endl << "Список потенциалов:" << std::endl;
    for (size_t iPotential = 0; iPotential < potential.size(); iPotential++) {
        (iPotential == winnerInd) ? potential[iPotential] += 1. / netConfig.neurons : potential[iPotential] -= netConfig.minPotential;
        std::cout << potential[iPotential] << std::endl;
    }
    std::cout << std::endl;
}
//--------------------------------------------------------------------------------------
void NetSOM::TrainGas(std::vector<std::vector<double>> &trainSet,
                     std::ostream &outputLabels) 
{
    size_t winnerInd;
    double maxT;
    // pretraining stage: solves problem of dead neurons
    std::random_shuffle(trainSet.begin(), trainSet.end());
    PrintWeightsToFile(outputLabels, 16);

    maxT = netConfig.preTrainIterations - 1; 
    for (size_t iVec = 0; iVec < netConfig.preTrainIterations; ++iVec) {
        double time = iVec + 1;  
        sigma = netConfig.sigmaInitPreTrain 
              * pow(netConfig.sigmaInitPreTrainMin / netConfig.sigmaInitPreTrain, time / maxT);
        eta   = netConfig.etaInitPreTrain   
              * pow(netConfig.etaInitPreTrainMin   / netConfig.etaInitPreTrain  , time / maxT);
        winnerInd = DetectWinnerGas(trainSet[iVec]);
        AdjustPotential(winnerInd);
        AdjustWeightsGas(trainSet[iVec]);
        PrintWeightsToFile(outputLabels, 16);
    }
    sleep(5);
    // adjust potential of each neuron equals to 1.
    for (int iPot = 0; iPot < potential.size(); iPot++) {
        potential[iPot] = 1.;
    }
    maxT = netConfig.trainEpochs * trainSet.size(); // k max = maxT, k = time
    double time = 0.;
    for (size_t iEpoch = 0; iEpoch < netConfig.trainEpochs; ++iEpoch) {
        std::random_shuffle(trainSet.begin(), trainSet.end());
        for (size_t iVec = 0; iVec < trainSet.size(); ++iVec) {
            time++;
            sigma = netConfig.sigmaInit * pow(netConfig.sigmaInitMin / netConfig.sigmaInit, time / maxT);
            eta   = netConfig.etaInit   * pow(netConfig.etaInitMin   / netConfig.etaInit  , time / maxT);
            winnerInd = DetectWinnerGas(trainSet[iVec]);
            //AdjustPotential(winnerInd);
            AdjustWeightsGas(trainSet[iVec]);
            PrintWeightsToFile(outputLabels, 16);
        }
    }
}
// ------------------------------------------------------------------------------------------------
void NetSOM::TrainKohen(std::vector<std::vector<double>> &trainSet, 
                        std::ostream &outputLabels) 
{
    size_t winnerInd;
    double maxT;
    // pretraining stage: solves problem of dead neurons
    std::random_shuffle(trainSet.begin(), trainSet.end());
    PrintWeightsToFile(outputLabels, 16);

    maxT = netConfig.preTrainIterations - 1; 
    for (size_t iVec = 0; iVec < netConfig.preTrainIterations; ++iVec) {
        double time = iVec + 1;       
        sigma = (netConfig.sigmaInitPreTrainMin - netConfig.sigmaInitPreTrain) 
              / (netConfig.preTrainIterations - 1) * (time - 1) + netConfig.sigmaInitPreTrain;
        eta   = (netConfig.etaInitPreTrainMin   - netConfig.etaInitPreTrain)   
              / (netConfig.preTrainIterations - 1) * (time - 1) + netConfig.etaInitPreTrain;
        
        winnerInd = DetectWinnerKohen(trainSet[iVec]);
        AdjustPotential(winnerInd);
        AdjustWeightsKohen(winnerInd, trainSet[iVec]);
        PrintWeightsToFile(outputLabels, 16);
    }
    // adjust potential of each neuron equals to 1.
    for (int iPot = 0; iPot < potential.size(); iPot++) {
        potential[iPot] = 1.;
    }
    maxT = netConfig.trainEpochs * trainSet.size(); // k max, k = time
    double time = 0.;
    for (size_t iEpoch = 0; iEpoch < netConfig.trainEpochs; ++iEpoch) {
        std::random_shuffle(trainSet.begin(), trainSet.end());
        for (size_t iVec = 0; iVec < trainSet.size(); ++iVec) {
            time++;
            sigma = (netConfig.sigmaInitMin - netConfig.sigmaInit) / (maxT - 1) * (time - 1) + netConfig.sigmaInit;
            eta   = (netConfig.etaInitMin - netConfig.etaInit)     / (maxT - 1) * (time - 1) + netConfig.etaInit;
            
            winnerInd = DetectWinnerKohen(trainSet[iVec]);
            AdjustWeightsKohen(winnerInd, trainSet[iVec]);
            PrintWeightsToFile(outputLabels, 16);
        }
    }
}
// ------------------------------------------------------------------------------------------------
void NetSOM::PrintWeightsToFile(std::ostream &output, int precision) {
    for (size_t iWeight = 0; iWeight < weights.size(); ++iWeight) {
        output << std::setprecision(precision) << weights[iWeight] << "\t";
        if (iWeight % 2) {
            output << iWeight / netConfig.inVecDim << "\t";
        }
    }
    output << std::endl << std::endl << std::endl;
}
// ------------------------------------------------------------------------------------------------
void NetSOM::CreateGnuplotAnimation(std::ofstream &stream) {
    stream << " set terminal gif size 1024, 768 animate delay 0.001 loop -1 "<< std::endl
              << " set output 'train.gif' "<< std::endl
              << " set xrange [0:4] "<< std::endl
              << " set yrange [0:4] "<< std::endl
              << " unset key "<< std::endl              
              << " stats \'" << netConfig.weightFileName << "\' nooutput  "<< std::endl
              << " do for [i=1:int(STATS_blocks)-1] { "<< std::endl
              << "     plot \"train.data\" index 0 using 1:2 pt 7 ps 2 lc rgb \'red\',\\"<< std::endl;
    for (size_t iNeuron = 0; iNeuron < netConfig.neurons - 1; iNeuron++) {
        stream << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
               <<  (iNeuron+1)*3 - 2 << ":" << (iNeuron+1)*3 - 1
               << "  pt 7 ps 5 lc rgb \'orange\',\\"<< std::endl;
        stream << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
               <<  (iNeuron+1)*3 - 2 << ":" << (iNeuron+1)*3 - 1 << ":" << (iNeuron+1)*3
               << "  with labels,\\"<< std::endl;
    }
    stream << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
           <<  (netConfig.neurons)*3 - 2 << ":" << (netConfig.neurons)*3 - 1
           << "  pt 7 ps 5 lc rgb \'orange\',\\"<< std::endl;
    stream << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
           <<  (netConfig.neurons)*3 - 2 << ":" << (netConfig.neurons)*3 - 1 << ":" << (netConfig.neurons)*3
           << "  with labels\\"<< std::endl;
    stream << "}" << std::endl;

    std::ofstream file;
    file.open("finalpos.txt", std::ios::trunc);
    file << " set terminal gif size 1024, 768 animate delay 0.001 loop -1 "<< std::endl
              << " set output 'final.gif' "<< std::endl
              << " set xrange [0:4] "<< std::endl
              << " set yrange [0:4] "<< std::endl
              << " unset key "<< std::endl              
              << " stats 'weight.data' nooutput  "<< std::endl
              << " do for [i=int(STATS_blocks)-1:int(STATS_blocks)-1] { "<< std::endl
              << "     plot \"train.data\" index 0 using 1:2 pt 7 ps 2 lc rgb \'red\',\\"<< std::endl;
    for (size_t iNeuron = 0; iNeuron < netConfig.neurons - 1; iNeuron++) {
        file << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
               <<  (iNeuron+1)*3 - 2 << ":" << (iNeuron+1)*3 - 1
               << "  pt 7 ps 5 lc rgb \'orange\',\\"<< std::endl;
        file << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
               <<  (iNeuron+1)*3 - 2 << ":" << (iNeuron+1)*3 - 1 << ":" << (iNeuron+1)*3
               << "  with labels,\\"<< std::endl;
    }
    file << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
           <<  (netConfig.neurons)*3 - 2 << ":" << (netConfig.neurons)*3 - 1
           << "  pt 7 ps 5 lc rgb \'orange\',\\"<< std::endl;
    file << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
           <<  (netConfig.neurons)*3 - 2 << ":" << (netConfig.neurons)*3 - 1 << ":" << (netConfig.neurons)*3
           << "  with labels\\"<< std::endl;
    file << "}" << std::endl;

    std::ofstream fileStart;
    fileStart.open("startpos.txt", std::ios::trunc);
    fileStart << " set terminal gif size 1024, 768 animate delay 0.001 loop -1 "<< std::endl
              << " set output 'startpos.gif' "<< std::endl
              << " set xrange [0:4] "<< std::endl
              << " set yrange [0:4] "<< std::endl
              << " unset key "<< std::endl              
              << " stats 'weight.data' nooutput  "<< std::endl
              << " do for [i=1:1] { "<< std::endl
              << "     plot \"train.data\" index 0 using 1:2 pt 7 ps 2 lc rgb \'red\',\\"<< std::endl;
    for (size_t iNeuron = 0; iNeuron < netConfig.neurons - 1; iNeuron++) {
        fileStart << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
               <<  (iNeuron+1)*3 - 2 << ":" << (iNeuron+1)*3 - 1
               << "  pt 7 ps 5 lc rgb \'orange\',\\"<< std::endl;
        fileStart << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
               <<  (iNeuron+1)*3 - 2 << ":" << (iNeuron+1)*3 - 1 << ":" << (iNeuron+1)*3
               << "  with labels,\\"<< std::endl;
    }
    fileStart << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
           <<  (netConfig.neurons)*3 - 2 << ":" << (netConfig.neurons)*3 - 1
           << "  pt 7 ps 5 lc rgb \'orange\',\\"<< std::endl;
    fileStart << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
           <<  (netConfig.neurons)*3 - 2 << ":" << (netConfig.neurons)*3 - 1 << ":" << (netConfig.neurons)*3
           << "  with labels\\"<< std::endl;
    fileStart << "}" << std::endl;

    std::ofstream filePreTrain;
    filePreTrain.open("afterPreTrain.txt", std::ios::trunc);
    int outIter = netConfig.preTrainIterations;
    filePreTrain << " set terminal gif size 1024, 768 animate delay 0.001 loop -1 "<< std::endl
              << " set output 'afterPreTrain.gif' "<< std::endl
              << " set xrange [0:4] "<< std::endl
              << " set yrange [0:4] "<< std::endl
              << " unset key "<< std::endl              
              << " stats 'weight.data' nooutput  "<< std::endl
              << " do for [i=" << outIter << ":" << outIter <<"] { "<< std::endl
              << "     plot \"train.data\" index 0 using 1:2 pt 7 ps 2 lc rgb \'red\',\\"<< std::endl;
    for (size_t iNeuron = 0; iNeuron < netConfig.neurons - 1; iNeuron++) {
        filePreTrain << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
               <<  (iNeuron+1)*3 - 2 << ":" << (iNeuron+1)*3 - 1
               << "  pt 7 ps 5 lc rgb \'orange\',\\"<< std::endl;
        filePreTrain << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
               <<  (iNeuron+1)*3 - 2 << ":" << (iNeuron+1)*3 - 1 << ":" << (iNeuron+1)*3
               << "  with labels,\\"<< std::endl;
    }
    filePreTrain << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
           <<  (netConfig.neurons)*3 - 2 << ":" << (netConfig.neurons)*3 - 1
           << "  pt 7 ps 5 lc rgb \'orange\',\\"<< std::endl;
    filePreTrain << "      \"" << netConfig.weightFileName <<"\" index(i-1) using " 
           <<  (netConfig.neurons)*3 - 2 << ":" << (netConfig.neurons)*3 - 1 << ":" << (netConfig.neurons)*3
           << "  with labels\\"<< std::endl;
    filePreTrain << "}" << std::endl;

}
