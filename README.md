# Interactive Shallow Water Model 🌊
Web-based interactive simulation of the Shallow Water Equations using Pyodide


[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT) <!-- Optional: Add license badge -->

![Simulation Demo Placeholder](images/simulation_demo.gif)
<!-- **Action:** Replace above line with path to an actual GIF/screenshot later -->

**Live Interactive Demo:** [**https://<YourUsername>.github.io/<YourRepoName>/web_interface/**](https://<YourUsername>.github.io/<YourRepoName>/web_interface/)
<!-- **Action:** Replace placeholder URL with your actual GitHub Pages link -->

## Overview

This project provides an interactive, web-based simulation of the **Shallow Water Equations (SWEs)** on a rotating plane (f-plane or beta-plane). It allows users to visually explore fundamental concepts in **Geophysical Fluid Dynamics (GFD)** directly in their browser, without needing any local installation beyond a modern web browser.

The simulation demonstrates key phenomena such as:

*   **Geostrophic Adjustment:** How initial imbalances evolve under rotation and gravity.
*   **Wave Propagation:** Visualizing gravity waves (Poincaré waves), Kelvin waves (in basin mode), and Rossby waves (with varying rotation).
*   **Rotational Effects:** The impact of the Coriolis force and its variation with latitude (beta-effect).
*   **Forcing:** Effects of wind stress (including curl-driven gyres) and bottom friction.

The interactive interface, powered by **Pyodide** (running Python/NumPy in WebAssembly) and **Plotly.js**, allows real-time parameter adjustments and visualization.

This repository also includes the original Python script (`original_script/shallow_water_ref.py`) based on Matplotlib, which served as the foundation for the web version.

**Target Audience:** Students, educators, researchers, and enthusiasts interested in fluid dynamics, oceanography, atmospheric science, or numerical modeling.

## Table of Contents

1.  [Physics Background](#physics-background)
    *   [Shallow Water Equations](#shallow-water-equations)
    *   [Key Concepts](#key-concepts)
    *   [Numerical Grid (Arakawa C-Grid)](#numerical-grid-arakawa-c-grid)
2.  [Features](#features)
3.  [Live Demo](#live-demo)
4.  [Usage](#usage)
    *   [Interactive Web Interface](#interactive-web-interface)
    *   [Original Python Script](#original-python-script)
5.  [Parameter Guide](#parameter-guide)
6.  [Exploring Phenomena (Examples)](#exploring-phenomena-examples)
7.  [Technical Details](#technical-details)
8.  [Contributing](#contributing)
9.  [License](#license)

## Physics Background

### Shallow Water Equations

This model solves the non-linear shallow water equations for a single, homogeneous fluid layer using the hydrostatic approximation:

*   **Zonal Momentum (U):** `dU/dT + U*dU/dX + V*dU/dY = +f*V - g*dH'/dX - C_drag*U + F_wind_x`
*   **Meridional Momentum (V):** `dV/dT + U*dV/dX + V*dV/dY = -f*U - g*dH'/dY - C_drag*V`
*   **Continuity (Height H'):** `dH'/dT + d(U*H)/dX + d(V*H)/dY = 0` (where `H = H_background + H'` is total depth)

*(Note: The current implementation uses a linearized version of the continuity equation `dH'/dT = -H_background * (dU/dX + dV/dY)` and may omit advection terms `U*dU/dX` etc. for simplicity/stability, but the core dynamics of pressure gradient, Coriolis, and divergence are captured).*

**Variables:**
*   `U, V`: Zonal (East-West) and Meridional (North-South) velocities.
*   `H'`: Deviation of fluid surface height from the mean background depth `H_background`.
*   `f`: Coriolis parameter (mimics planetary rotation effect).
*   `g`: Gravitational acceleration.
*   `C_drag`: Bottom drag coefficient (linear damping).
*   `F_wind_x`: Zonal wind forcing term.
*   `d/dT, d/dX, d/dY`: Derivatives with respect to time and space.

### Key Concepts

*   **Geostrophic Balance:** Tendency for large-scale, slow flows where the pressure gradient force (`g*dH'/dX`) balances the Coriolis force (`f*V`). Flow follows lines of constant height.
*   **Coriolis Effect:** Apparent force deflecting moving objects on a rotating plane (right in NH, left in SH). Modeled via the parameter `f`.
*   **Beta-Effect:** The variation of the Coriolis parameter `f` with latitude (`beta = df/dy`). Crucial for Rossby waves. Enabled by 'PlusMinus' or 'WithLatitude' rotation schemes.
*   **Waves:**
    *   **Gravity Waves (Poincaré Waves):** Fast waves adjusting pressure/velocity imbalances. Speed ~`sqrt(gH)`.
    *   **Kelvin Waves:** Gravity waves trapped along boundaries (requires walls).
    *   **Rossby Waves:** Slow, large-scale waves requiring the beta-effect. Propagate westward relative to the mean flow.

### Numerical Grid (Arakawa C-Grid)

A staggered grid is used for numerical stability and accuracy:
*   **Height (`H'`)**: Cell centers.
*   **Zonal Velocity (`U`)**: Vertical cell faces (East/West edges).
*   **Meridional Velocity (`V`)**: Horizontal cell faces (North/South edges).

  ------- V(i+1,j) -------     (North Edge)
  |                      |
U(i,j)   H(i,j)        U(i,j+1)  (Row i)
  |                      |
  ------- V(i,j) ---------     (South Edge)
  
(West Edge)       (East Edge)
        (column j)


## Features ✨

*   **Interactive & Web-Based:** Runs entirely in your browser using Pyodide (Python compiled to WebAssembly). No server needed!
*   **Live Parameter Control:** Adjust grid size, time step, physics, rotation, forcing, and visualization parameters on the fly.
*   **Real-time Visualization:** Uses Plotly.js heatmaps to show the evolving height field.
*   **Robust Numerics:** Implements the 4th-order Runge-Kutta (RK4) time-stepping scheme.
*   **C-Grid Staggering:** Employs standard geophysical fluid dynamics grid staggering.
*   **Interpolated Coriolis:** Includes accurate calculation of Coriolis terms on the staggered grid.
*   **Multiple Scenarios:** Supports various initial conditions, rotation schemes, boundary conditions (channel/basin), and wind forcing.
*   **Educational Guide:** Built-in help modal explains parameters and key phenomena.
*   **Original Script Included:** Provides the foundational Matplotlib-based Python script for comparison.

## Live Demo 🚀

Experience the interactive simulation directly here:

**[https://<YourUsername>.github.io/<YourRepoName>/web_interface/](https://<YourUsername>.github.io/<YourRepoName>/web_interface/)**
<!-- **Action:** Replace placeholder URL with your actual GitHub Pages link -->

*(Requires a modern web browser supporting WebAssembly. Initial loading may take a few moments while Pyodide initializes.)*

## Usage

### Interactive Web Interface

1.  **Visit the Live Demo link above.**
2.  **Wait for Initialization:** The status message will indicate "Status: Ready." when Pyodide and the script are loaded.
3.  **Adjust Parameters:** Use the controls in the "Simulation Controls", "Rotation", and "Visualization" panels to set up your desired experiment. See the [Parameter Guide](#parameter-guide) below or click the `?` button on the page for details.
4.  **Run:** Click the "Run Simulation" button.
5.  **Observe:** Watch the height anomaly evolve in the "Visualization" panel. Simulation time is shown in the status area.
6.  **Stop:** Click the "Stop" button to pause the simulation. You can then change parameters and click "Run" again to continue from the *beginning* with new settings.

### Original Python Script

The original script (`original_script/shallow_water_ref.py`) runs locally and uses Matplotlib for visualization.

1.  **Prerequisites:**
    *   Python 3.x
    *   NumPy (`pip install numpy`)
    *   Matplotlib (`pip install matplotlib`)
2.  **Navigate:** Open a terminal in the `original_script/` directory.
3.  **Configure:** Edit the script file (`shallow_water_ref.py`) directly to set parameters like `ncol`, `dT`, `initialPerturbation`, `rotationScheme`, `plotOutput`, etc., near the top of the file.
4.  **Run:** Execute the script from the terminal:
    ```bash
    python shallow_water_ref.py
    ```
5.  **View:** If `plotOutput = True`, a Matplotlib window will appear showing the animated simulation.

## Parameter Guide

Parameters adjustable in the web interface:

*   **Grid & Time:**
    *   `Grid Size (N x N)`: Simulation resolution (NxN cells).
    *   `Grid Spacing (dX - km)`: Physical size of cells. Affects stability (smaller `dX` needs smaller `dT`).
    *   `Time Step (dT - secs)`: Duration of simulation time steps. Governed by the CFL stability condition.
    *   `Steps per Frame`: Number of `dT` steps calculated between plot updates. Affects animation speed.
*   **Initial State & Physics:**
    *   `Initial Condition`: Starting shape of `H'` (`Tower`, `NSGradient`, `EWGradient`, `Random`).
    *   `Gravity (G - m/s² x10⁻⁴)`: Affects wave speed (`sqrt(gH)`). Value in UI is scaled.
    *   `Mean Depth (H - m)`: Background water depth. Affects wave speed.
    *   `Drag Coeff (x10⁻⁷ /s)`: Linear bottom friction strength. Value in UI is scaled.
    *   `Periodic E-W`: Toggles between a periodic channel (checked) and a closed basin with walls (unchecked).
*   **Wind Forcing:**
    *   `Wind Forcing`: Type of external wind stress (`Off`, `Uniform`, `Curled`).
*   **Rotation:**
    *   `Rotation Scheme`: How the Coriolis parameter `f` varies (`Off`, `Uniform` f-plane, `PlusMinus`/`WithLatitude` beta-plane).
    *   `Mean Latitude (°N)`: Reference latitude for calculating `f`.
    *   `Beta Factor`: Strength of `df/dy` for the `PlusMinus` scheme.
*   **Visualization:**
    *   `Use Fixed Scale`: Toggles dynamic color scaling vs. fixed `Min/Max`.
    *   `Fixed Min/Max (H)`: User-defined color limits when fixed scale is active.

*(Click the `?` button on the web page for more detailed parameter descriptions within the application.)*

## Exploring Phenomena (Examples)

Use the interactive tool to investigate:

1.  **Geostrophic Adjustment:** Start with `Tower`, `PlusMinus Rotation`. Observe the initial outward flow deflect into rotation. Compare with `Rotation Scheme: Off`.
2.  **Gravity/Poincaré Waves:** Watch ripples spread from the `Tower`. Change `G` or `HBackground` to see how wave speed changes.
3.  **Kelvin Waves:** Uncheck `Periodic E-W` (basin mode), start with `Tower`. Look for waves propagating along the boundaries.
4.  **Rossby Waves:** Use `PlusMinus` or `WithLatitude` rotation. Try `Random` initial conditions or `NSGradient`. Look for slow, large-scale westward propagation patterns (may need longer runs/lower `ntAnim`).
5.  **Wind-Driven Gyres:** Uncheck `Periodic E-W` (basin), select `Curled` wind, maybe reduce `Drag Coeff`. Run for longer to see large-scale circulation develop.
6.  **Stability (CFL):** Gradually increase `dT` for a given `ncol` and `dX`. Find the point where the simulation becomes unstable (e.g., shows extreme values or NaNs). How does this point change if you decrease `dX`?

## Technical Details

*   **Core Logic:** Python with NumPy for array operations.
*   **Numerics:** Arakawa C-Grid, RK4 time stepping, Centered differences, Interpolated Coriolis terms.
*   **Web Execution:** Pyodide compiles Python/NumPy to WebAssembly, running directly in the browser.
*   **Visualization:** Plotly.js JavaScript library (`heatmap` trace type).
*   **Interface:** HTML, CSS, JavaScript.

## Contributing

Contributions, bug reports, and feature suggestions are welcome! Please feel free to open an Issue or submit a Pull Request.

1.  Fork the repository.
2.  Create a new branch (`git checkout -b feature/YourFeature` or `bugfix/YourBug`).
3.  Make your changes.
4.  Commit your changes (`git commit -m 'Add some feature'`).
5.  Push to the branch (`git push origin feature/YourFeature`).
6.  Open a Pull Request.

## License

This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.