import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Model : decreasing exp. with offset
def model(t, A, mu, C):
    return A * np.exp(- mu *t) + C

# --- Data ---
t_wall_model, y_wall_model = np.loadtxt("log_wall_model.dat", skiprows=0, usecols=(1, 2), unpack=True)
t=t_wall_model
y=y_wall_model
# --- Initial guess ---
A0 = y.max() - y.min()
C0 = y.min()
mu0 = (t.max() - t.min()) / 3
p0 = [A0, mu0, C0]

# --- Least squares  ---
popt_wall_model, pcov = curve_fit(model, t, y, p0=p0)
perr = np.sqrt(np.diag(pcov))  # incertitudes (1 sigma)

A, mu_wall_model, C = popt_wall_model
print(f"A   = {A:.3f} ± {perr[0]:.3f}")
print(rf"$\mu$ = {mu_wall_model:.3f} ± {perr[1]:.3f}")
print(f"C   = {C:.3f} ± {perr[2]:.3f}")

# --- Data ---
t_slip, y_slip = np.loadtxt("log_slip.dat", skiprows=0, usecols=(1, 2), unpack=True)
t=t_slip
y=y_slip
# --- Initial guess ---
A0 = y.max() - y.min()
C0 = y.min()
mu0 = (t.max() - t.min()) / 3
p0 = [A0, mu0, C0]

# --- Least Squares ---
popt_slip, pcov = curve_fit(model, t, y, p0=p0)
perr = np.sqrt(np.diag(pcov))  # incertitudes (1 sigma)

A, mu_slip, C = popt_slip
print(f"A   = {A:.3f} ± {perr[0]:.3f}")
print(rf"$\mu$ = {mu_slip:.3f} ± {perr[1]:.3f}")
print(f"C   = {C:.3f} ± {perr[2]:.3f}")

# --- Plot ---
plt.scatter(t_wall_model, y_wall_model, s=10, alpha=0.5, label="Raw data wall model")
plt.plot(t_wall_model, model(t_wall_model, *popt_wall_model), "r-", lw=2, label=rf"Fit ($\mu$ = {mu_wall_model:.2f})")
plt.scatter(t_slip, y_slip, s=10, alpha=0.5, label="Raw data slip")
plt.plot(t_slip, model(t_slip, *popt_slip), "r-", lw=2, label=rf"Fit ($\mu$ = {mu_slip:.2f})")
plt.xlabel("t")
plt.ylabel(r"$E/E_0$")
plt.legend()
plt.show()
