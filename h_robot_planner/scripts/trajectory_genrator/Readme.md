

# Trajectory Generator — Catmull–Rom (Cardinal) Spline

This module generates a smooth 2D trajectory through user (or random) waypoints using a **Catmull–Rom spline** with **centripetal parameterization** (α = 0.5). It provides interpolated positions, headings (from tangents), and an approximate arc-length parameter. A uniform–arc-length resampling step produces evenly spaced samples suitable for controllers (e.g., MPC followers).

---

## 1. Mathematical Foundations

### 1.1 Centripetal Catmull–Rom Knot Parameterization
Given control points \(\mathbf{P}_i = (x_i, y_i)^T\), define a non‑uniform knot sequence \(t_i\) by
\[
 t_{i+1} = t_i + \lVert \mathbf{P}_{i+1} - \mathbf{P}_i \rVert^{\alpha}, \quad \alpha \in [0,1].
\]
The **centripetal** choice \(\alpha = 0.5\) is widely used because it limits cusps and local self‑intersections in Catmull–Rom curves while better adhering to the control points.  
*References:* centripetal Catmull–Rom definition and motivation; see standard treatments and surveys.

### 1.2 Hermite Segment Representation
For each segment \([\mathbf{P}_i,\mathbf{P}_{i+1}]\), let \(u \in [0,1]\) be the normalized parameter with \(\Delta t = t_{i+1} - t_i\). Using finite‑difference tangents
\[
 \mathbf{m}_i   \approx \frac{\mathbf{P}_{i+1} - \mathbf{P}_{i-1}}{t_{i+1} - t_{i-1}},\qquad
 \mathbf{m}_{i+1}\approx \frac{\mathbf{P}_{i+2} - \mathbf{P}_i}{t_{i+2} - t_i},
\]
we evaluate the position in **cubic Hermite** form:
\[
 \mathbf{P}(u) = h_{00}(u)\,\mathbf{P}_i + h_{10}(u)\,\Delta t\,\mathbf{m}_i + h_{01}(u)\,\mathbf{P}_{i+1} + h_{11}(u)\,\Delta t\,\mathbf{m}_{i+1},
\]
with basis functions
\[
 h_{00}(u)=2u^3-3u^2+1,\quad h_{10}(u)=u^3-2u^2+u,\quad h_{01}(u)=-2u^3+3u^2,\quad h_{11}(u)=u^3-u^2.
\]
The derivative \(\tfrac{\mathrm d\mathbf{P}}{\mathrm du}\) (obtained by differentiating the basis) provides the **tangent direction**; the **heading** for a mobile base is then
\[
 \theta(u) = \operatorname{atan2}\big( (\tfrac{\mathrm d\mathbf{P}}{\mathrm du})_y,\ (\tfrac{\mathrm d\mathbf{P}}{\mathrm du})_x \big).
\]

### 1.3 Arc‑Length Accumulation and Resampling
Uniform sampling in \(u\) is not uniform in distance along the curve. We therefore:
1. **Sample densely** per segment and accumulate chord length to obtain an approximate arc‑length parameter \(s\).
2. **Resample** to \(N\) points by inverting \(s\) with linear interpolation, producing **nearly equal spatial spacing**. This improves numerical stability and timing for controllers.

---

## 2. Code Structure (no implementation pasted here)

- **`spline_catmull_rom.hpp` (header‑only)**  
  – centripetal knot computation  
  – Hermite segment evaluation (position + derivative)  
  – dense curve sampling and arc‑length accumulation  
  – equal‑arc‑length resampling  
  – random waypoint generation and CSV writer

- **`trajectory_using_splines.cpp` (demo & tests)**  
  – lightweight sanity checks (endpoint interpolation; monotone arc length; resampling endpoints)  
  – generates random waypoints, samples the spline, resamples by arc length  
  – writes `/tmp/trajectory_spline.csv` with columns: `x,y,theta_deg,s`

- **`plot_trajectory.py` (visualization)**  
  – reads the CSV and displays a 2D plot (equal aspect; start/goal marked)

---

## 3. Build & run (ad‑hoc)
```bash
cd ~/HOST/personal_ws/src/h_robot_planner/scripts/trajectory_genrator
g++ -O2 -std=c++17 trajectory_using_splines.cpp -o spline_demo
./spline_demo
# → writes /tmp/trajectory_spline.csv
```

## 4. Visualize
```bash
python3 plot_trajectory.py /tmp/trajectory_spline.csv
```
A 2D plot appears showing the interpolated path and endpoints.

## 5. Integrating into ROS 2 (optional)
This module can be wrapped as a ROS 2 node that:
- accepts waypoints via a service or topic (e.g., `nav_msgs/Path` or custom),
- publishes the generated trajectory as `nav_msgs/Path` (for RViz/MPC inputs),
- optionally exports CSVs for debugging.

---

## References
- Centripetal Catmull–Rom spline (definition; knot parameterization).  
- Cubic Hermite splines and basis functions (segment formulation; basis).  
- On parameterizations that avoid cusps/self‑intersections; centripetal recommendation.  
- Practical notes on arc‑length and chord‑length sampling.