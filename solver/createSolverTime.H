//
// createSolverTime.H
//
#include "createTime.H"

const scalar deltaT = readScalar(runTime.controlDict().lookup("deltaT"));
const scalar endTime = readScalar(runTime.controlDict().lookup("endTime"));
const label n = endTime/deltaT;
