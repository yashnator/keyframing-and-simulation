#include "numerical_methods.hpp"
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

double interpolateSegment(double t, double t0, double t1, double t2, double t3,
                          double P0, double P1, double P2, double P3) {
    // Compute tangents
    double m1 = ((P2 - P1) * ((t2 - t1) / (t2 - t0))) + ((P1 - P0) * ((t1 - t0) / (t2 - t0)));
    double m2 = ((P3 - P2) * ((t2 - t1) / (t3 - t1))) + ((P2 - P1) * ((t3 - t2) / (t3 - t1)));

    // Normalize t to [0, 1]
    double tau = (t - t1) / (t2 - t1);

    // Hermite basis functions
    double h00 = 2 * tau * tau * tau - 3 * tau * tau + 1;
    double h10 = tau * tau * tau - 2 * tau * tau + tau;
    double h01 = -2 * tau * tau * tau + 3 * tau * tau;
    double h11 = tau * tau * tau - tau * tau;

    return h00 * P1 + h10 * m1 * (t2 - t1) + h01 * P2 + h11 * m2 * (t2 - t1);
}

double evaluateCatmullRom(double t,
                          const std::vector<double>& points,
                          const std::vector<double>& times) {
    if (points.size() != 5 || times.size() != 5)
        throw std::runtime_error("Expected exactly 5 points and 5 times.");

    if (t <= times[1]) {
        return interpolateSegment(t, times[0], times[0], times[1], times[2],
                                  points[0], points[0], points[1], points[2]);
    } else if (t <= times[2]) {
        return interpolateSegment(t, times[0], times[1], times[2], times[3],
                                  points[0], points[1], points[2], points[3]);
    } else if (t <= times[3]) {
        return interpolateSegment(t, times[1], times[2], times[3], times[4],
                                  points[1], points[2], points[3], points[4]);
    } else {
        return interpolateSegment(t, times[2], times[3], times[4], times[4],
                                  points[2], points[3], points[4], points[4]);
    }
}
