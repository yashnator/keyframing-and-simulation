#ifndef NUMERICAL_METHODS_HPP
#define NUMERICAL_METHODS_HPP

#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

float catmullRom5NonUniform(
    float t0, float t1, float t2, float t3, float t4,
    float p0, float p1, float p2, float p3, float p4,
    float t // t in [t1, t4]
);

double interpolateSegment(double t, double t0, double t1, double t2, double t3,
    double P0, double P1, double P2, double P3);


double evaluateCatmullRom(double t,
                                  const std::vector<double>& points,
                                  const std::vector<double>& times);

#endif