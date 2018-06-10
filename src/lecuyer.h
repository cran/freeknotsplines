#ifndef LECUYER_H
#define LECUYER_H

#include "RngStream.h"
#include <math.h>

void RandomInitialise(int stream, int nstreams, int seed);
double RandomUniform(void);
double RandomGaussian(double mean,double stddev);
int RandomInt(int lower,int upper);
double RandomDouble(double lower,double upper);

#endif
