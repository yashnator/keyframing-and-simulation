#ifndef NUMERICAL_METHODS_HPP
#define NUMERICAL_METHODS_HPP

float catmullRom5NonUniform(
    float t0, float t1, float t2, float t3, float t4,
    float p0, float p1, float p2, float p3, float p4,
    float t // t in [t1, t4]
);

#endif