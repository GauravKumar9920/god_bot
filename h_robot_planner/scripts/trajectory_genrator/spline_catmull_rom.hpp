#pragma once
// -----------------------------------------------------------------------------
// Catmull–Rom / Cardinal spline (2D) utilities — header-only, no deps
// -----------------------------------------------------------------------------
// Math references (background):
//  • Centripetal Catmull–Rom parameterization uses t_{k+1} = t_k + ||P_{k+1}-P_k||^alpha,
//    with alpha=0.5 recommended to avoid cusps/loops and follow control points tightly.
//    (Wikipedia: "Centripetal Catmull–Rom spline")
//  • Segment evaluation in cubic Hermite form on u∈[0,1]:
//      P(u) = h00(u) P_i + h10(u) Δt m_i + h01(u) P_{i+1} + h11(u) Δt m_{i+1}
//    where Δt = t_{i+1}-t_i and m_i ≈ dP/dt estimated with finite differences.
// -----------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <cassert>
#include <random>
#include <fstream>
#include <string>

namespace hmpc {

struct Waypoint { double x{0}, y{0}; };
struct Sample   { double x{0}, y{0}; double dx{0}, dy{0}; double s{0}; }; // s = cumulative arc length

// ---- Cubic Hermite basis on [0,1] ----
inline double h00(double u){ return  2*u*u*u - 3*u*u + 1; }
inline double h10(double u){ return      u*u*u - 2*u*u + u; }
inline double h01(double u){ return -2*u*u*u + 3*u*u; }
inline double h11(double u){ return      u*u*u -   u*u; }

// Derivatives w.r.t u (unit parameter)
inline double dh00(double u){ return  6*u*u - 6*u; }
inline double dh10(double u){ return  3*u*u - 4*u + 1; }
inline double dh01(double u){ return -6*u*u + 6*u; }
inline double dh11(double u){ return  3*u*u - 2*u; }

// ---- Knot computation (centripetal by default) ----
inline std::vector<double> compute_knots(const std::vector<Waypoint>& P, double alpha=0.5){
    const size_t n = P.size();
    assert(n >= 2);
    std::vector<double> t(n);
    t[0] = 0.0;
    for(size_t i=1;i<n;++i){
        const double dx = P[i].x - P[i-1].x;
        const double dy = P[i].y - P[i-1].y;
        const double dist = std::sqrt(dx*dx + dy*dy);
        t[i] = t[i-1] + std::pow(dist, alpha);
    }
    return t;
}

// ---- Evaluate one segment i→i+1 at u∈[0,1] (open curve via end duplication) ----
inline Sample eval_segment(const std::vector<Waypoint>& P,
                           const std::vector<double>& t,
                           size_t i, double u)
{
    const size_t N = P.size();
    size_t im1 = (i==0? 0 : i-1);
    size_t ip1 = (i+1 < N? i+1 : N-1);
    size_t ip2 = (i+2 < N? i+2 : N-1);

    const Waypoint Pm1 = P[im1];
    const Waypoint P0  = P[i];
    const Waypoint P1  = P[ip1];
    const Waypoint P2  = P[ip2];

    const double ti   = t[i];
    const double tim1 = t[im1];
    const double ti1  = t[ip1];
    const double ti2  = t[ip2];

    const double dt = (ti1 - ti);
    const double denom0 = (ti1 - tim1);
    const double denom1 = (ti2 - ti);

    // Tangents (finite-difference estimates of dP/dt)
    const double mi_x  = (P1.x - Pm1.x) / (denom0 > 1e-12 ? denom0 : 1.0);
    const double mi_y  = (P1.y - Pm1.y) / (denom0 > 1e-12 ? denom0 : 1.0);
    const double mi1_x = (P2.x - P0.x)  / (denom1 > 1e-12 ? denom1 : 1.0);
    const double mi1_y = (P2.y - P0.y)  / (denom1 > 1e-12 ? denom1 : 1.0);

    // Hermite position
    const double x = h00(u)*P0.x + h10(u)*dt*mi_x + h01(u)*P1.x + h11(u)*dt*mi1_x;
    const double y = h00(u)*P0.y + h10(u)*dt*mi_y + h01(u)*P1.y + h11(u)*dt*mi1_y;

    // Derivative w.r.t. u (direction sufficient for heading)
    const double dux = dh00(u)*P0.x + dh10(u)*dt*mi_x + dh01(u)*P1.x + dh11(u)*dt*mi1_x;
    const double duy = dh00(u)*P0.y + dh10(u)*dt*mi_y + dh01(u)*P1.y + dh11(u)*dt*mi1_y;

    return {x,y, dux,duy, 0.0};
}

// ---- Sample entire open curve with M samples per segment ----
inline std::vector<Sample> sample_curve(const std::vector<Waypoint>& P,
                                        int samples_per_seg=50,
                                        double alpha=0.5)
{
    assert(P.size() >= 2);
    auto t = compute_knots(P, alpha);

    std::vector<Sample> out;
    out.reserve(static_cast<size_t>(P.size())*samples_per_seg + 1);

    for(size_t i=0; i+1<P.size(); ++i){
        const int add_one = (i+1 == P.size()-1) ? 1 : 0; // include u=1 on last seg only
        for(int k=0; k<samples_per_seg + add_one; ++k){
            const double u = (samples_per_seg==0) ? 0.0 : static_cast<double>(k) / static_cast<double>(samples_per_seg);
            out.push_back(eval_segment(P, t, i, u));
        }
    }

    // accumulate chord-length as an approx arc-length parameter s
    double accum = 0.0;
    for(size_t i=0;i<out.size();++i){
        if(i>0){
            const double dx = out[i].x - out[i-1].x;
            const double dy = out[i].y - out[i-1].y;
            accum += std::sqrt(dx*dx + dy*dy);
        }
        out[i].s = accum;
    }
    return out;
}

// ---- Uniform resampling by arc length (N samples) ----
inline std::vector<Sample> resample_equal_arc(const std::vector<Sample>& in, int N){
    std::vector<Sample> out; out.reserve(N > 0 ? N : 0);
    if(in.empty()||N<=0) return out;
    const double L = in.back().s;
    if(L <= 1e-12){ out.push_back(in.front()); return out; }

    for(int k=0;k<N;++k){
        const double target = (L * k)/static_cast<double>(N-1);
        size_t j=1; while(j<in.size() && in[j].s < target) ++j;
        if(j==in.size()) { out.push_back(in.back()); break; }
        const auto &a=in[j-1], &b=in[j];
        const double t = (target - a.s)/std::max(1e-12, b.s - a.s);
        Sample p;
        p.x = a.x + t*(b.x - a.x);
        p.y = a.y + t*(b.y - a.y);
        p.dx= a.dx + t*(b.dx- a.dx);
        p.dy= a.dy + t*(b.dy- a.dy);
        p.s = target;
        out.push_back(p);
    }
    return out;
}

// ---- RNG waypoints ----
inline std::vector<Waypoint> make_random_waypoints(int n, unsigned seed,
                                                   double xmin, double xmax,
                                                   double ymin, double ymax,
                                                   bool anchor=true){
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dx(xmin, xmax), dy(ymin, ymax);
    std::vector<Waypoint> P; P.reserve(n > 0 ? n : 0);
    for(int i=0;i<n;++i) P.push_back({dx(rng), dy(rng)});
    if(anchor && n>=2){ P.front() = {0.0,0.0}; P.back() = {6.0,3.0}; }
    return P;
}

// ---- CSV dump: x,y,theta(deg),s ----
inline void write_csv(const std::string& path, const std::vector<Sample>& S){
    std::ofstream f(path);
    f << "x,y,theta_deg,s\n";
    for(const auto& p: S){
        const double th = std::atan2(p.dy, p.dx) * 180.0 / M_PI;
        f << p.x << "," << p.y << "," << th << "," << p.s << "\n";
    }
}

} // namespace hmpc
