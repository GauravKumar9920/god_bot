"""
Holonomic robot MPC path-following using a Cardinal (Catmull–Rom) spline
Platform: macOS (Apple Silicon / M3 Pro)

Setup (recommended with conda-forge on Apple Silicon):

  # 1) Install Miniforge (conda-forge installer): https://conda-forge.org/download/
  # 2) Create env (Python 3.11 works well with CasADi on ARM as of 2025):
  #    mamba create -n holo-mpc python=3.11 numpy matplotlib casadi -c conda-forge
  #    conda activate holo-mpc
  # 3) Extra pip deps (splines for Catmull–Rom):
  #    pip install splines
  # (Optional) If you want TOPP-RA time-optimal speed profiles later:
  #    pip install toppra

Run:
  python holonomic_mpc_spline.py

This script does:
  1) Sample random 2-D waypoints.
  2) Build a centripetal Catmull–Rom spline through them (a cardinal spline with tension implicitly handled by parameterization; tweakable via `alpha` and optional tangent scaling).
  3) Re-parameterize to unit-speed (arc length) for easier timing.
  4) Set up a simple nonlinear MPC for a holonomic robot (x, y, theta) with controls (vx, vy, omega).
  5) Track the spline while gently maximizing forward progress and minimizing control energy.

Notes:
  • This is a didactic, compact NMPC. For performance, you’d pre-compile and warm-start the solver.
  • The MPC here uses a fixed step dt; later we’ll plug in TOPP-RA to optimize the time-law (s-dot) along the path.
"""
from __future__ import annotations
import math
import numpy as np
import casadi as ca
import matplotlib.pyplot as plt

try:
    from splines import CatmullRom, UnitSpeedAdapter
except Exception:
    raise SystemExit("Missing dependency 'splines'. Install with: pip install splines")

# -------------------------------
# 1) Waypoints & Spline
# -------------------------------

def make_random_waypoints(n:int=8, seed:int=42, box=((-5,5),(-5,5))):
    rng = np.random.default_rng(seed)
    xs = rng.uniform(box[0][0], box[0][1], size=n)
    ys = rng.uniform(box[1][0], box[1][1], size=n)
    pts = np.stack([xs, ys], axis=1)
    # Ensure distinct start/goal anchors for nicer demos
    start = np.array([0.0, 0.0])
    goal  = np.array([6.0, 3.0])
    pts[0] = start
    pts[-1] = goal
    return pts


def build_catmull_rom(pts: np.ndarray, alpha: float = 0.5, closed: bool = False):
    """Build a centripetal Catmull–Rom spline (alpha=0.5).
    If closed=True, treat path as a loop. We'll keep it open for this project.
    Returns a callable p(s) via UnitSpeedAdapter (arc-length parameterized).
    """
    # CatmullRom takes vertices and optional grid. With alpha given, it uses non-uniform (e.g., centripetal) parameterization internally.
    cr = CatmullRom(pts, alpha=alpha, endconditions='natural' if not closed else 'closed')
    # Arc-length (unit-speed) reparameterization -> parameter s is arc length
    p_unit = UnitSpeedAdapter(cr)
    return cr, p_unit


def sample_path(p_unit, s_max: float, num:int=1000):
    s_vals = np.linspace(0, s_max, num)
    XY = np.array([p_unit.evaluate(s) for s in s_vals])
    # Numerically estimate tangents for plotting and progress cost
    dXY = np.gradient(XY, s_max/(num-1), axis=0)
    T = dXY / (np.linalg.norm(dXY, axis=1, keepdims=True) + 1e-9)
    return s_vals, XY, T

# -------------------------------
# 2) Simple Holonomic NMPC
# -------------------------------
class HolonomicMPC:
    def __init__(self, N=20, dt=0.05,
                 v_bounds=(-1.5, 1.5), omega_bounds=(-2.5, 2.5),
                 a_limits=(3.0, 3.0, 5.0),
                 w_pos=5.0, w_theta=0.2, w_u=0.01, w_progress=0.8):
        self.N = N
        self.dt = dt
        self.vmin, self.vmax = -abs(v_bounds[1]), abs(v_bounds[1])
        self.omin, self.omax = omega_bounds
        self.ax_max, self.ay_max, self.aomega_max = a_limits
        self.w_pos, self.w_theta, self.w_u, self.w_prog = w_pos, w_theta, w_u, w_progress

        nx, nu = 3, 3
        X = ca.SX.sym('X', nx, N+1)
        U = ca.SX.sym('U', nu, N)
        X0 = ca.SX.sym('X0', nx)
        Xref = ca.SX.sym('Xref', 3, N+1)
        Tref = ca.SX.sym('Tref', 2, N)

        dt = self.dt
        cost = 0
        g_terms = []
        eq_mask = []  # True for equality constraints
        lbU = []
        ubU = []

        for k in range(N):
            xk = X[:, k]
            uk = U[:, k]
            xk1 = X[:, k+1]
            x_next = ca.vertcat(
                xk[0] + dt * uk[0],
                xk[1] + dt * uk[1],
                xk[2] + dt * uk[2],
            )
            dyn = xk1 - x_next
            g_terms.append(dyn)
            eq_mask.extend([True]*3)

            err_pos = Xref[0:2, k] - xk[0:2]
            err_th = ca.atan2(ca.sin(Xref[2, k] - xk[2]), ca.cos(Xref[2, k] - xk[2]))
            cost += self.w_pos * ca.dot(err_pos, err_pos)
            cost += self.w_theta * (err_th**2)
            cost += self.w_u * ca.dot(uk, uk)
            cost += - self.w_prog * (Tref[0, k] * uk[0] + Tref[1, k] * uk[1])

            if k > 0:
                du = (U[:, k] - U[:, k-1]) / dt
                rate_exprs = [
                    du[0] - self.ax_max,
                    -du[0] - self.ax_max,
                    du[1] - self.ay_max,
                    -du[1] - self.ay_max,
                    du[2] - self.aomega_max,
                    -du[2] - self.aomega_max,
                ]
                for expr in rate_exprs:
                    g_terms.append(expr)
                    eq_mask.append(False)

            lbU += [self.vmin, self.vmin, self.omin]
            ubU += [self.vmax, self.vmax, self.omax]

        # Terminal tracking cost
        err_pos_T = Xref[0:2, N] - X[0:2, N]
        err_th_T = ca.atan2(ca.sin(Xref[2, N] - X[2, N]), ca.cos(Xref[2, N] - X[2, N]))
        cost += 2.0*self.w_pos*ca.dot(err_pos_T, err_pos_T) + 2.0*self.w_theta*(err_th_T**2)

        # Initial state equality constraint
        g_terms.insert(0, X[:,0] - X0)
        eq_mask = [True]*nx + eq_mask

        g = ca.vertcat(*g_terms)
        vars = ca.vertcat(ca.reshape(X, -1, 1), ca.reshape(U, -1, 1))

        self.NX = X
        self.NU = U
        self.X0 = X0
        self.Xref = Xref
        self.Tref = Tref
        self.vars = vars
        self.g = g
        self.eq_mask = eq_mask
        self.cost = cost

        nlp = {
            'x': vars,
            'f': cost,
            'g': g,
            'p': ca.vertcat(ca.reshape(Xref, -1, 1), ca.reshape(Tref, -1, 1), X0)
        }
        ipopt_opts = {
            'ipopt.print_level': 0,
            'print_time': 0,
            'ipopt.max_iter': 500,
        }
        solver_tried = []
        for method, m_opts in [
            ('ipopt', ipopt_opts),
            ('sqpmethod', {'print_time': 0}),
            ('qrqp', {'print_time': 0}),
        ]:
            try:
                self.solver = ca.nlpsol('solver', method, nlp, m_opts)
                if method != 'ipopt':
                    print(f"[HolonomicMPC] Using fallback solver '{method}' (ipopt unavailable).")
                break
            except Exception as e:
                solver_tried.append(f"{method}:{e.__class__.__name__}")
                self.solver = None
        if self.solver is None:
            raise RuntimeError("No suitable CasADi NLP solver plugin available. Tried: " + ", ".join(solver_tried))

        self.lbU = np.array(lbU)
        self.ubU = np.array(ubU)

    def solve(self, x0: np.ndarray, xref: np.ndarray, tref: np.ndarray):
        N = self.N
        nx = 3
        P = np.concatenate([xref.reshape(-1), tref.reshape(-1), x0.reshape(-1)])
        total_g = self.g.shape[0]
        lbg = np.zeros(total_g)
        ubg = np.zeros(total_g)
        mask = np.array(self.eq_mask, dtype=bool)
        lbg[~mask] = -np.inf

        X_guess = np.tile(x0.reshape(-1,1), (1, N+1))
        U_guess = np.zeros((3, N))
        z0 = np.concatenate([X_guess.reshape(-1), U_guess.reshape(-1)])

        nX = (N+1)*nx
        lbz = -np.inf * np.ones_like(z0)
        ubz =  np.inf * np.ones_like(z0)
        lbz[nX:nX+3*N] = self.lbU
        ubz[nX:nX+3*N] = self.ubU

        try:
            sol = self.solver(x0=z0, lbg=lbg, ubg=ubg, lbx=lbz, ubx=ubz, p=P)
            z = np.array(sol['x']).reshape(-1)
            Xsol = z[:nX].reshape(nx, N+1)
            Usol = z[nX:].reshape(3, N)
            return Xsol, Usol
        except Exception as e:
            print(f"[HolonomicMPC] Solve failed ({e.__class__.__name__}): {e}. Returning initial guess.")
            return X_guess, U_guess

# -------------------------------
# 3) Utilities for building references along the path
# -------------------------------

def nearest_index(pt: np.ndarray, path: np.ndarray) -> int:
    d = np.linalg.norm(path - pt[None,:], axis=1)
    return int(np.argmin(d))


def headings_from_tangent(T: np.ndarray) -> np.ndarray:
    return np.arctan2(T[:,1], T[:,0])

# -------------------------------
# 4) Demo
# -------------------------------
if __name__ == "__main__":
    # Waypoints & spline
    W = make_random_waypoints(n=10, seed=3)
    cr, p_unit = build_catmull_rom(W, alpha=0.5, closed=False)

    # Estimate total arc length by sampling with high resolution
    s_max = 1.0 * np.linalg.norm(W[1:] - W[:-1], axis=1).sum()  # crude upper bound
    s_samp, XY_samp, T_samp = sample_path(p_unit, s_max, num=1500)
    THETA_samp = headings_from_tangent(T_samp)

    # MPC
    mpc = HolonomicMPC(N=20, dt=0.05,
                       v_bounds=(-1.5, 1.5), omega_bounds=(-2.0, 2.0),
                       a_limits=(3.0, 3.0, 6.0),
                       w_pos=6.0, w_theta=0.2, w_u=0.02, w_progress=1.0)

    # Sim loop
    x = np.array([W[0,0], W[0,1], 0.0], dtype=float)  # start at first waypoint facing +x
    traj = [x.copy()]
    us = []

    K = 240  # ~12s at dt=0.05
    lookahead = mpc.N  # use N future samples for references

    for k in range(K):
        # Find nearest point along sampled path
        i0 = nearest_index(x[:2], XY_samp)
        i1 = min(i0 + lookahead, XY_samp.shape[0]-1)
        ref_xy = XY_samp[i0:i1+1]
        ref_th = THETA_samp[i0:i1+1]
        if ref_xy.shape[0] < lookahead+1:
            # pad with last point/heading
            pad = lookahead+1 - ref_xy.shape[0]
            ref_xy = np.vstack([ref_xy, np.repeat(ref_xy[-1][None,:], pad, axis=0)])
            ref_th = np.concatenate([ref_th, np.repeat(ref_th[-1], pad)])

        Xref = np.vstack([ref_xy.T, ref_th[None,:]])  # (3, N+1)
        Tref = T_samp[i0:i0+lookahead].T              # (2, N)

        Xsol, Usol = mpc.solve(x0=x, xref=Xref, tref=Tref)
        u = Usol[:,0]

        # Apply first control (simulate)
        x = np.array([
            x[0] + mpc.dt * u[0],
            x[1] + mpc.dt * u[1],
            x[2] + mpc.dt * u[2],
        ])
        traj.append(x.copy())
        us.append(u.copy())

        # Stop if we are close to goal
        if np.linalg.norm(x[:2] - W[-1]) < 0.10 and np.linalg.norm(Usol[:2, :].sum(axis=1)) < 1e-2:
            break

    traj = np.array(traj)
    us = np.array(us)

    # ---------------- Plot ----------------
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    ax.plot(W[:,0], W[:,1], 'o--', label='Waypoints')
    ax.plot(XY_samp[:,0], XY_samp[:,1], '-', lw=2, label='Catmull–Rom (centripetal)')
    ax.plot(traj[:,0], traj[:,1], '-', lw=2, label='MPC track')

    # Draw a few heading arrows
    for j in range(0, traj.shape[0], max(1, traj.shape[0]//20)):
        ax.arrow(traj[j,0], traj[j,1], 0.25*math.cos(traj[j,2]), 0.25*math.sin(traj[j,2]),
                 head_width=0.12, length_includes_head=True, alpha=0.7)

    ax.axis('equal')
    ax.grid(True)
    ax.legend()
    ax.set_title('Holonomic NMPC following a Catmull–Rom spline')
    plt.tight_layout()
    plt.show()

    # Print a tiny summary
    total_time = len(traj)*mpc.dt
    print(f"Sim steps: {len(traj)}  |  Total time: {total_time:.2f}s  |  Avg |u|: {np.linalg.norm(us, axis=0).mean():.3f}")
