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

} // namespace hmpc#include "spline_catmull_rom.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cassert>

// =============================================================
// Trajectory generation over a Catmull–Rom (centripetal) spline
// =============================================================
// Mathematics (high level):
//  • Given waypoints P_i in R^2, build a non-uniform parameter t_i via
//      t_{i+1} = t_i + ||P_{i+1}-P_i||^alpha,  alpha=0.5 (centripetal)
//    This reduces cusps/loops for irregularly spaced points.
//  • For each segment [P_i, P_{i+1}], define tangents (dP/dt) at t_i and t_{i+1}
//      m_i   ≈ (P_{i+1}-P_{i-1}) / (t_{i+1}-t_{i-1})
//      m_{i+1}≈ (P_{i+2}-P_i)   / (t_{i+2}-t_i)
//  • Interpolate with cubic Hermite basis on u∈[0,1], where u=(t−t_i)/Δt, Δt=t_{i+1}-t_i:
//      P(u)  = h00(u) P_i + h10(u) Δt m_i + h01(u) P_{i+1} + h11(u) Δt m_{i+1}
//    dP/du gives the tangent direction; heading θ = atan2(dy, dx).
//  • Optionally resample to N points with near-uniform arc length by linear search in
//    cumulative chord-length s.
//
// Output: CSV of (x,y,theta,s) and simple console checks.
// =============================================================

using namespace hmpc;

static void unit_tests(){
    // 1) Interpolation at knots: ensure first/last sample match endpoints within tol
    {
        std::vector<Waypoint> W = {{0,0},{1,0},{2,1},{3,1}};
        auto S = sample_curve(W, /*samples_per_seg=*/50, /*alpha=*/0.5);
        double tol = 1e-9;
        assert(std::fabs(S.front().x - W.front().x) < tol && std::fabs(S.front().y - W.front().y) < tol);
        assert(std::fabs(S.back().x  - W.back().x ) < tol && std::fabs(S.back().y  - W.back().y ) < tol);
    }
    // 2) Monotone arc-length: s should be nondecreasing
    {
        auto W = make_random_waypoints(6, 7, -2, 2, -1, 2, /*anchor=*/false);
        auto S = sample_curve(W, 30, 0.5);
        for(size_t i=1;i<S.size();++i) assert(S[i].s + 1e-12 >= S[i-1].s);
    }
    // 3) Resampling keeps end points
    {
        std::vector<Waypoint> W = {{0,0},{1,1},{2,0}};
        auto S = sample_curve(W, 40, 0.5);
        auto R = resample_equal_arc(S, 25);
        assert(!R.empty());
        assert(std::fabs(R.front().x - S.front().x) < 1e-6);
        assert(std::fabs(R.back().x  - S.back().x ) < 1e-6);
    }
}

int main(int argc, char** argv){
    try { unit_tests(); std::cout << "[tests] OK\n"; }
    catch(...){ std::cerr << "[tests] FAILED\n"; return 1; }

    // --- Config ---
    int n_waypoints = 10;
    int samples_per_seg = 60;
    int resample_N = 300;
    double alpha = 0.5; // centripetal

    // --- Waypoints ---
    auto W = make_random_waypoints(n_waypoints, /*seed=*/3, -5, 5, -5, 5, /*anchor=*/true);

    // --- Sample spline ---
    auto S = sample_curve(W, samples_per_seg, alpha);
    auto R = resample_equal_arc(S, resample_N);

    // --- Output ---
    write_csv("/tmp/trajectory_spline.csv", R);

    // Print a few samples
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Generated " << R.size() << " samples. First/last:\n";
    if(!R.empty()){
        const auto &a=R.front(), &b=R.back();
        std::cout << "  start: (" << a.x << ", " << a.y << ")  theta(deg)=" << std::atan2(a.dy,a.dx)*180.0/M_PI << "\n";
        std::cout << "  goal : (" << b.x << ", " << b.y << ")  theta(deg)=" << std::atan2(b.dy,b.dx)*180.0/M_PI << "\n";
        std::cout << "  total length s = " << R.back().s << " m (approx)\n";
    }

    std::cout << "CSV written to /tmp/trajectory_spline.csv\n";
    return 0;
}
#!/usr/bin/env python3
import csv
import math
import os
import sys

try:
    import matplotlib.pyplot as plt
except Exception as e:
    print("Matplotlib not available. Install with: pip install matplotlib", file=sys.stderr)
    sys.exit(1)

path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/trajectory_spline.csv"
if not os.path.exists(path):
    print(f"CSV not found: {path}")
    sys.exit(2)

xs, ys, ths = [], [], []
with open(path) as f:
    r = csv.DictReader(f)
    for row in r:
        xs.append(float(row['x']))
        ys.append(float(row['y']))
        ths.append(float(row['theta_deg']))

plt.figure(figsize=(6,6))
plt.plot(xs, ys, '-', linewidth=2)
plt.scatter([xs[0], xs[-1]], [ys[0], ys[-1]], s=60)
plt.axis('equal')
plt.grid(True)
plt.title('Spline Trajectory')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.show()
# Trajectory Generator (Catmull–Rom Spline)

This module generates a smooth 2D trajectory through a set of waypoints using a **Catmull–Rom (cardinal) spline** with **centripetal** parameterization (alpha = 0.5). It outputs sampled positions, tangents (for headings), and approximate arc length, and it can resample to equal arc-length spacing. The result is ideal as input to a holonomic MPC tracker.

## Files
- `spline_catmull_rom.hpp` — Header-only spline utilities (no external deps).
- `trajectory_using_splines.cpp` — Demo + sanity tests; writes `/tmp/trajectory_spline.csv`.
- `plot_trajectory.py` — Quick matplotlib visualizer for the CSV.

## Math summary (in-code comments have the formulas)
1. **Centripetal parameterization**:  
   Define non-uniform knots \(t_i\) by
   \[ t_{i+1} = t_i + \lVert P_{i+1}-P_i \rVert^{\alpha}, \quad \alpha=0.5 \]
   This choice avoids cusps/self-intersections common to the uniform variant and hugs the waypoints more tightly. (See the *Centripetal Catmull–Rom spline* article.)
2. **Hermite form per segment**: For segment \([P_i,P_{i+1}]\) with \(u\in[0,1]\),
   \[ P(u) = h_{00}(u)P_i + h_{10}(u)\Delta t\,m_i + h_{01}(u)P_{i+1} + h_{11}(u)\Delta t\,m_{i+1} \]
   where \(\Delta t=t_{i+1}-t_i\) and \(m_i\) are finite-difference tangents. We also compute \(\tfrac{\mathrm dP}{\mathrm du}\) for headings.
3. **Arc-length resampling**: After dense sampling, we do a linear search on cumulative chord length \(s\) to pick \(N\) nearly equidistant points (good for controllers and timing).

## Build & run (ad-hoc)
```bash
# Build just the demo without ROS build system
cd ~/HOST/personal_ws/src/god_bot/h_robot_planner/scripts/trajectory_genrator
g++ -O2 -std=c++17 trajectory_using_splines.cpp -o spline_demo
./spline_demo
# -> writes /tmp/trajectory_spline.csv
```

## Visualize
```bash
# Requires matplotlib
python3 plot_trajectory.py /tmp/trajectory_spline.csv
```
You should see a smooth curve passing through the randomly generated waypoints; start and goal are anchored at (0,0) and (6,3) by default.

## Integrating into ROS 2 (optional)
We can wrap this as a small ROS 2 node that:
- takes a `nav_msgs/Path` or a list of waypoints on a service/topic,
- publishes a `nav_msgs/Path` with the resampled spline points and headings,
- optionally dumps the CSV for offline inspection.

## What’s next (MPC)
This trajectory (x, y, heading from the tangent, arc-length) is ready for a holonomic MPC tracker with states \([x,y,\theta]\) and controls \([v_x, v_y, \omega]\). We’ll:
1. Build a tracking MPC cost: position/heading error + control effort + progress reward.
2. Add speed/accel limits; later, integrate a time-optimal retiming step along \(s\).
