#include <cmath>
#include <algorithm>
#include <iostream>
float catmullRom5NonUniform(
    float t0, float t1, float t2, float t3, float t4,
    float p0, float p1, float p2, float p3, float p4,
    float t // t in [t1, t4]
) {
    // Find which segment t is in
    if (t <= t2) {
        // Segment [t1, t2]: use p0, p1, p2, p3
        float m1 = (p2 - p0) / (t2 - t0);
        float m2 = (p3 - p1) / (t3 - t1);
        float s = (t - t1) / (t2 - t1);
        float s2 = s * s;
        float s3 = s2 * s;
        // std::cout<<(2*s3 - 3*s2 + 1) * p1
        // + (s3 - 2*s2 + s) * (t2 - t1) * m1
        // + (-2*s3 + 3*s2) * p2
        // + (s3 - s2) * (t2 - t1) * m2<<std::endl;
        return (2*s3 - 3*s2 + 1) * p1
             + (s3 - 2*s2 + s) * (t2 - t1) * m1
             + (-2*s3 + 3*s2) * p2
             + (s3 - s2) * (t2 - t1) * m2;
    } else if (t <= t3) {
        // Segment [t2, t3]: use p1, p2, p3, p4
        float m1 = (p3 - p1) / (t3 - t1);
        float m2 = (p4 - p2) / (t4 - t2);
        float s = (t - t2) / (t3 - t2);
        float s2 = s * s;
        float s3 = s2 * s;
        return (2*s3 - 3*s2 + 1) * p2
             + (s3 - 2*s2 + s) * (t3 - t2) * m1
             + (-2*s3 + 3*s2) * p3
             + (s3 - s2) * (t3 - t2) * m2;
    } else {
        // Segment [t3, t4]: use last four points
        float m1 = (p4 - p2) / (t4 - t2);
        float m2 = (p4 - p3) / (t4 - t3); // one-sided at end
        float s = (t - t3) / (t4 - t3);
        float s2 = s * s;
        float s3 = s2 * s;
        return (2*s3 - 3*s2 + 1) * p3
             + (s3 - 2*s2 + s) * (t4 - t3) * m1
             + (-2*s3 + 3*s2) * p4
             + (s3 - s2) * (t4 - t3) * m2;
    }
}

