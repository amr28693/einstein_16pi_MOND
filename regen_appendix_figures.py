"""
A script to accurately reproduce the figures from the appendices.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3D)
from scipy.ndimage import gaussian_filter


# =====================================================================
# Appendix A: Imaginary Curvature as Entropy Gradient
# =====================================================================
def fig_A(out_path="R_i_total_curvature.png"):
    grid_size = 100
    x = np.linspace(-5, 5, grid_size)
    y = np.linspace(-5, 5, grid_size)
    X, Y = np.meshgrid(x, y)

    real_curvature = np.exp(-(X**2 + Y**2))
    theta_field = np.arctan2(Y, X)
    dtheta_dx = np.gradient(theta_field, axis=1)
    dtheta_dy = np.gradient(theta_field, axis=0)
    imag_curvature = dtheta_dx**2 + dtheta_dy**2
    total_magnitude = np.sqrt(real_curvature**2 + imag_curvature**2)

    fig = plt.figure(figsize=(18, 5))

    ax1 = fig.add_subplot(131, projection='3d')
    ax1.plot_surface(X, Y, real_curvature, cmap='viridis')
    ax1.set_title("Real Curvature (GR Projection)")
    ax1.set_zlabel("R^Re")

    ax2 = fig.add_subplot(132, projection='3d')
    ax2.plot_surface(X, Y, imag_curvature, cmap='plasma')
    ax2.set_title("Imaginary Curvature (Entropy Flow)")
    ax2.set_zlabel("R^Im")

    ax3 = fig.add_subplot(133, projection='3d')
    ax3.plot_surface(X, Y, total_magnitude, cmap='inferno')
    # ONLY CHANGE vs V2 original: dropped "(Information Geometry)" suffix
    ax3.set_title("Total |R^C|")
    ax3.set_zlabel("|R^C|")

    plt.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


# =====================================================================
# Appendix B: Solitonic Curvature Projection
# =====================================================================
def fig_B(out_path="soliton_gradient_fields.png"):
    x = np.linspace(-5, 5, 200)
    y = np.linspace(-5, 5, 200)
    X, Y = np.meshgrid(x, y)

    def real_soliton(x, y):
        return np.exp(-((x + 2)**2 + y**2))

    def imag_soliton(x, y):
        return np.exp(-((x - 2)**2 + y**2))

    R = real_soliton(X, Y)
    I = imag_soliton(X, Y)
    C = R - I

    R_smooth = gaussian_filter(R, sigma=1)
    I_smooth = gaussian_filter(I, sigma=1)
    C_smooth = gaussian_filter(C, sigma=1)

    grad_R = np.gradient(R_smooth)
    grad_I = np.gradient(I_smooth)
    grad_C = np.gradient(C_smooth)

    mag_grad_R = np.sqrt(grad_R[0]**2 + grad_R[1]**2)
    mag_grad_I = np.sqrt(grad_I[0]**2 + grad_I[1]**2)
    mag_grad_C = np.sqrt(grad_C[0]**2 + grad_C[1]**2)

    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    im0 = axs[0].imshow(mag_grad_R, extent=(-5, 5, -5, 5),
                        origin='lower', cmap='plasma')
    axs[0].set_title('Gradient Magnitude: Real Soliton')
    fig.colorbar(im0, ax=axs[0])

    im1 = axs[1].imshow(mag_grad_I, extent=(-5, 5, -5, 5),
                        origin='lower', cmap='plasma')
    axs[1].set_title('Gradient Magnitude: Imaginary Soliton')
    fig.colorbar(im1, ax=axs[1])

    im2 = axs[2].imshow(mag_grad_C, extent=(-5, 5, -5, 5),
                        origin='lower', cmap='plasma')
    axs[2].set_title('Gradient Magnitude: Combined Curvature Field')
    fig.colorbar(im2, ax=axs[2])

    plt.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {out_path}")
    for name, arr in [("Real Gradient", mag_grad_R),
                      ("Imaginary Gradient", mag_grad_I),
                      ("Combined Gradient", mag_grad_C)]:
        print(f"  {name}: min={arr.min():.5f} max={arr.max():.5f} "
              f"mean={arr.mean():.5f} std={arr.std():.5f}")


# =====================================================================
# Appendix C: Entropy-Driven Projection at tau = 0
# =====================================================================
def fig_C(out_path="GR16pi_projection_simulation.png"):
    x = np.linspace(-5, 5, 300)
    tau = np.linspace(-5, 5, 300)
    X, T = np.meshgrid(x, tau)

    A = np.exp(-0.5 * (X**2 + T**2))
    theta = np.sin(2 * np.pi * X) * np.tanh(T)

    Psi = A * np.exp(1j * theta)
    Re_Psi = np.real(Psi)
    Im_Psi = np.imag(Psi)

    entropy_density = gaussian_filter(
        np.abs(np.gradient(theta, axis=0))**2
        + np.abs(np.gradient(theta, axis=1))**2,
        sigma=1
    )

    tau_index = np.argmin(np.abs(tau))
    Psi_proj = Psi[tau_index, :]
    curvature_real = np.abs(np.gradient(np.real(Psi_proj)))**2
    curvature_imag = np.abs(np.gradient(np.imag(Psi_proj)))**2
    curvature_combined = curvature_real + curvature_imag

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))

    im0 = axs[0, 0].imshow(Re_Psi,
                           extent=[x.min(), x.max(), tau.min(), tau.max()],
                           aspect='auto', origin='lower', cmap='RdBu')
    axs[0, 0].set_title(r"Real Part of $\Psi(x, \tau)$")
    axs[0, 0].set_xlabel(r"$x$")
    axs[0, 0].set_ylabel(r"$\tau$")
    axs[0, 0].axhline(0, color='white', linestyle='--')
    plt.colorbar(im0, ax=axs[0, 0])

    im1 = axs[0, 1].imshow(Im_Psi,
                           extent=[x.min(), x.max(), tau.min(), tau.max()],
                           aspect='auto', origin='lower', cmap='PiYG')
    axs[0, 1].set_title(r"Imaginary Part of $\Psi(x, \tau)$")
    axs[0, 1].set_xlabel(r"$x$")
    axs[0, 1].set_ylabel(r"$\tau$")
    axs[0, 1].axhline(0, color='white', linestyle='--')
    plt.colorbar(im1, ax=axs[0, 1])

    im2 = axs[1, 0].imshow(entropy_density,
                           extent=[x.min(), x.max(), tau.min(), tau.max()],
                           aspect='auto', origin='lower', cmap='inferno')
    axs[1, 0].set_title(r"Entropy Density $|\nabla \theta(x, \tau)|^2$")
    axs[1, 0].set_xlabel(r"$x$")
    axs[1, 0].set_ylabel(r"$\tau$")
    axs[1, 0].axhline(0, color='cyan', linestyle='--')
    plt.colorbar(im2, ax=axs[1, 0])

    axs[1, 1].plot(x, curvature_combined,
                   label='Projected Curvature', color='black')
    axs[1, 1].fill_between(x, curvature_real, alpha=0.3, label='Real')
    axs[1, 1].fill_between(x, curvature_imag, alpha=0.3, label='Imaginary')
    axs[1, 1].set_title(r"Projected Curvature at $\tau=0$")
    axs[1, 1].set_xlabel(r"$x$")
    axs[1, 1].set_ylabel("Curvature Magnitude")
    axs[1, 1].legend()

    plt.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {out_path}")
    for name, arr in [("Real Gradient", curvature_real),
                      ("Imaginary Gradient", curvature_imag),
                      ("Combined Gradient", curvature_combined)]:
        print(f"  {name}: min={arr.min():.5f} max={arr.max():.5f} "
              f"mean={arr.mean():.5f} std={arr.std():.5f}")


if __name__ == "__main__":
    fig_A()
    fig_B()
    fig_C()
