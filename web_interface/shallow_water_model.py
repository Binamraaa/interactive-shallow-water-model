# shallow_water_model.py (Adapted for Pyodide, v2)
# Core simulation logic - NO plotting or command-line interface here.
import numpy as np
import math
import time # Keep time for potential profiling inside python

print("Python: shallow_water_model.py parsing started...")

# -------------------------------------
# Physics and Numerics Functions
# -------------------------------------
def calculate_coriolis_parameter(params):
    """Calculates the Coriolis parameter f for each row."""
    # Use .get() for safer access to parameters
    nrow = params.get("ncol", 32) # Default ncol if missing
    dX = params.get("dX")
    meanLatitude = params.get("meanLatitude")

    # --- Add Check for essential parameters ---
    if dX is None or meanLatitude is None:
         print(f"PYTHON ERROR in calculate_coriolis: dX or meanLatitude is None.")
         print(f"dX type: {type(dX)}, value: {dX}")
         print(f"meanLatitude type: {type(meanLatitude)}, value: {meanLatitude}")
         print(f"Full params received: {params}")
         raise ValueError(f"dX ({dX}) or meanLatitude ({meanLatitude}) is None inside calculate_coriolis_parameter.")
    # --- End Check ---

    dxDegrees = dX / 110.e3
    rotConst = np.zeros(nrow)
    latitude = np.zeros(nrow)

    for irow in range(nrow):
        lat = meanLatitude + (irow - (nrow - 1) / 2.0) * dxDegrees
        latitude[irow] = lat
        f = 0.0
        scheme = params.get("rotationScheme", "Off")
        betaFactor = params.get("betaFactor", 0.8) # Get betaFactor here

        if scheme == "WithLatitude":
             f = -7.e-5 * math.sin(math.radians(lat))
        elif scheme == "PlusMinus":
             f0 = -3.5e-5
             f = f0 * (1.0 - betaFactor * (irow - (nrow - 1) / 2.0) / (nrow / 2.0))
        elif scheme == "Uniform":
            f = -3.5e-5
        rotConst[irow] = f
    return rotConst, latitude

def calculate_wind_forcing(params):
    """Calculates the wind forcing term (acceleration) for U."""
    nrow = params.get("ncol", 32)
    windU_force = np.zeros(nrow)
    scheme = params.get("windScheme", "")
    if scheme == "Curled":
        max_stress = 1e-8 # Example forcing value
        windU_force = max_stress * np.sin((np.arange(nrow) + 0.5) / nrow * 2 * np.pi)
    elif scheme == "Uniform":
        windU_force = 1.e-8 * np.ones(nrow)
    return windU_force[:, np.newaxis] # Ensure shape (nrow, 1)

def apply_boundary_conditions(U, V, H, params):
    """Applies boundary conditions to state variables IN PLACE."""
    nrow, ncol = H.shape[0], H.shape[1] - 1
    horizontalWrap = params.get("horizontalWrap", True)

    # North/South Walls (V=0)
    V[0, :] = 0.0
    V[nrow, :] = 0.0 # row index nrow is the ghost row at the 'south' edge

    # East/West Boundaries
    if horizontalWrap:
        H[:, ncol] = H[:, 0]   # H ghost cell (East) = H first column (West)
        U[:, 0] = U[:, ncol]   # U on West boundary = U on East boundary
    else: # Walls
        U[:, 0] = 0.0    # West wall
        U[:, ncol] = 0.0 # East wall
        # H ghost cell: Zero gradient Neumann condition H(east_ghost)=H(east_real)
        H[:, ncol] = H[:, ncol-1]

def calculate_tendencies(U, V, H, params, coriolis_f, wind_U):
    """Calculates time tendencies (dU/dt, dV/dt, dH/dt)."""
    # Parameter retrieval with defaults
    nrow, ncol = H.shape[0], H.shape[1] - 1
    dX = params.get("dX", 20.E3)
    dY = dX # Assume square grid
    G = params.get("G", 9.8e-4)
    HBackground = params.get("HBackground", 1000)
    dragConst = params.get("dragConst", 5.E-7)
    interp_rot = params.get("interpolateRotation", True)
    horizontalWrap = params.get("horizontalWrap", True)

    # Initialize tendency arrays
    dUdT = np.zeros_like(U)
    dVdT = np.zeros_like(V)
    dHdT = np.zeros((nrow, ncol)) # Only defined for interior H points

    # --- Spatial Derivatives ---
    # dH/dX at U locations (Shape: nrow, ncol+1)
    dHdX = np.zeros_like(U)
    dHdX[:, 1:ncol] = (H[:, 1:ncol] - H[:, 0:ncol-1]) / dX
    if horizontalWrap:
        dHdX[:, 0] = (H[:, 0] - H[:, ncol-1]) / dX
        dHdX[:, ncol] = dHdX[:, 0]
    # else: Walls -> dHdX remains 0 at cols 0 and ncol

    # dH/dY at V locations (Shape: nrow+1, ncol)
    dHdY = np.zeros_like(V)
    # Difference between H North (row i) and H South (row i-1) -> result for V[i]
    # Valid for V rows 1 to nrow-1 (interior points)
    dHdY[1:nrow, :] = (H[1:nrow, 0:ncol] - H[0:nrow-1, 0:ncol]) / dY
    # Boundaries (V[0,:], V[nrow,:]) have dHdY = 0 (handled by init)

    # dU/dX at H locations (Shape: nrow, ncol)
    dUdX = (U[:, 1:ncol+1] - U[:, 0:ncol]) / dX

    # dV/dY at H locations (Shape: nrow, ncol)
    dVdY = (V[1:nrow+1, :] - V[0:nrow, :]) / dY

    # --- Coriolis Terms ---
    rotV_term = np.zeros_like(U) # Term for dUdT
    rotU_term = np.zeros_like(V) # Term for dVdT

    # Only calculate if rotation is enabled and interpolation is on
    if interp_rot and coriolis_f is not None and np.any(coriolis_f):
        # Calculate f*U and f*V at H cell centers (nrow, ncol)
        f_matrix_H = coriolis_f[:, np.newaxis]
        U_at_H = 0.5 * (U[:, 0:ncol] + U[:, 1:ncol+1])
        V_at_H = 0.5 * (V[0:nrow, :] + V[1:nrow+1, :])
        rotU_at_H = f_matrix_H * U_at_H
        rotV_at_H = f_matrix_H * V_at_H

        # Interpolate f*V back to U locations (nrow, ncol+1)
        rotV_term[:, 1:ncol] = 0.5 * (rotV_at_H[:, 0:ncol-1] + rotV_at_H[:, 1:ncol])
        if horizontalWrap:
            rotV_term[:, 0] = 0.5 * (rotV_at_H[:, ncol-1] + rotV_at_H[:, 0])
            rotV_term[:, ncol] = rotV_term[:, 0]
        # else: rotV_term remains 0 at walls

        # Interpolate f*U back to V locations (nrow+1, ncol)
        rotU_term[1:nrow, :] = 0.5 * (rotU_at_H[0:nrow-1, :] + rotU_at_H[1:nrow, :])
        # else: rotU_term remains 0 at walls (rows 0 and nrow)

    # --- Assemble Tendencies ---
    # dH/dT (Shape: nrow, ncol)
    dHdT = -HBackground * (dUdX + dVdY)

    # dU/dT (Shape: nrow, ncol+1)
    # rotV_term is +fV term. wind_U has shape (nrow, 1) and broadcasts.
    dUdT = rotV_term - G * dHdX - dragConst * U + wind_U

    # dV/dT (Shape: nrow+1, ncol)
    # rotU_term is +fU term, so we need -rotU_term
    dVdT = -rotU_term - G * dHdY - dragConst * V

    return dUdT, dVdT, dHdT


def step_rk4(U, V, H, params, coriolis_f, wind_U):
    """Advances the state variables by one RK4 time step."""
    dt = params.get("dT", 150) # Get dT safely
    U_orig, V_orig, H_orig = U.copy(), V.copy(), H.copy() # Keep original state

    # --- RK4 Stages ---
    # k1
    # Apply BCs to the state *before* calculating tendency
    apply_boundary_conditions(U, V, H, params)
    k1_U, k1_V, k1_H = calculate_tendencies(U, V, H, params, coriolis_f, wind_U)

    # k2
    temp_U = U_orig + 0.5 * dt * k1_U
    temp_V = V_orig + 0.5 * dt * k1_V
    temp_H = H_orig.copy()
    temp_H[:, :-1] += 0.5 * dt * k1_H # Apply dHdT only to non-ghost cells
    apply_boundary_conditions(temp_U, temp_V, temp_H, params) # Apply BCs to intermediate state
    k2_U, k2_V, k2_H = calculate_tendencies(temp_U, temp_V, temp_H, params, coriolis_f, wind_U)

    # k3
    temp_U = U_orig + 0.5 * dt * k2_U
    temp_V = V_orig + 0.5 * dt * k2_V
    temp_H = H_orig.copy()
    temp_H[:, :-1] += 0.5 * dt * k2_H
    apply_boundary_conditions(temp_U, temp_V, temp_H, params)
    k3_U, k3_V, k3_H = calculate_tendencies(temp_U, temp_V, temp_H, params, coriolis_f, wind_U)

    # k4
    temp_U = U_orig + dt * k3_U
    temp_V = V_orig + dt * k3_V
    temp_H = H_orig.copy()
    temp_H[:, :-1] += dt * k3_H
    apply_boundary_conditions(temp_U, temp_V, temp_H, params)
    k4_U, k4_V, k4_H = calculate_tendencies(temp_U, temp_V, temp_H, params, coriolis_f, wind_U)

    # --- Final Update using original state and weighted tendencies ---
    U_new = U_orig + (dt / 6.0) * (k1_U + 2*k2_U + 2*k3_U + k4_U)
    V_new = V_orig + (dt / 6.0) * (k1_V + 2*k2_V + 2*k3_V + k4_V)
    H_new = H_orig.copy()
    # Apply dHdT sum only to non-ghost H cells (indices :, 0 to ncol-1)
    H_new[:, :-1] += (dt / 6.0) * (k1_H + 2*k2_H + 2*k3_H + k4_H)

    # Enforce final boundary conditions strictly AFTER the full update
    apply_boundary_conditions(U_new, V_new, H_new, params)

    return U_new, V_new, H_new

# -------------------------------------
# Initialization (Adapted for web use)
# -------------------------------------
def initialize_state(params):
    """Initializes the U, V, H state variables based on params dict."""
    nrow = params.get("ncol", 32)
    ncol = nrow # Assuming square grid

    U = np.zeros((nrow, ncol + 1))
    V = np.zeros((nrow + 1, ncol))
    H = np.zeros((nrow, ncol + 1)) # Includes ghost cell H[:, ncol]

    perturbation = params.get("initialPerturbation", "Tower")
    print(f"Python: Initializing state with {perturbation} for grid {nrow}x{ncol}")

    if perturbation == "Tower":
        mid_row, mid_col = int(nrow / 2), int(ncol / 2)
        sigma = max(1.0, float(ncol) / 8.0) # Ensure float division
        for r in range(nrow):
            for c in range(ncol):
                dist_sq = ((r - mid_row)**2 + (c - mid_col)**2) / sigma**2
                H[r, c] = 1.0 * np.exp(-0.5 * dist_sq)
    elif perturbation == "NSGradient":
        for irow in range(nrow):
           H[irow, 0:ncol] = 0.5 * (irow - (nrow - 1) / 2.0) * (2.0 / nrow)
    elif perturbation == "EWGradient":
        for icol in range(ncol):
           H[:, icol] = 0.5 * (icol - (ncol - 1) / 2.0) * (2.0 / ncol)
    elif perturbation == "Random":
        H[:, 0:ncol] = 0.1 * (np.random.rand(nrow, ncol) - 0.5)

    # Apply initial boundary conditions
    apply_boundary_conditions(U, V, H, params)
    print(f"Python: Initial H range: {H[:, :-1].min():.3f} to {H[:, :-1].max():.3f}")
    return U, V, H


# --- Define Default Parameters Dictionary ---
# These defaults are read by JS on initialization
DEFAULT_PARAMS = {
    # Grid & Time
    "ncol": 20,
    "dX": 10000.0,         # Base unit: meters
    "dT": 150.0,
    "ntAnim": 5,           # Steps per frame
    # Initial State & Physics
    "initialPerturbation": "Tower",
    "G": 9.8e-4,
    "HBackground": 1000.0,
    "dragConst": 5.e-7,
    "horizontalWrap": True,
    # Wind
    "windScheme": "",      # Default no wind
    # Rotation
    "interpolateRotation": True, # Keep True for stability
    "rotationScheme": "PlusMinus",
    "meanLatitude": 30.0,
    "betaFactor": 0.8,
    # Visualization defaults (used by JS) - optional to include here
    "useFixedScale": False,
    "fixedZmin": -0.5,
    "fixedZmax": 0.5,
}


# --- Final confirmation print ---
print("Python: shallow_water_model.py parsed successfully.")