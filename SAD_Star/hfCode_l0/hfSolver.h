#ifndef HFSOLVER_H
#define HFSOLVER_H


/**
 * @brief This class solves the HF equation
 */

#include <iomanip>
#include "VMinnesotaMatrixGenerator.h"

class HfSolver
{

public:
    HfSolver();
    ~HfSolver();

    void setParam(double b, int nMax, int nPart);
    void run(double& HFEnergy);

    void read2BodyMat(VMinnesotaMatrixGenerator::TwoBodyMat& mat, string file);

private:
    double b_;
    int nMax_;
    int nPart_;

};

#endif // HFSOLVER_H
