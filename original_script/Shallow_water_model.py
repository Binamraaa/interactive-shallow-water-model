# Set up python environment.  numpy and matplotlib will have to be installed
# with the python installation.

import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import math

# Grid and Variable Initialization -- stuff you might play around with

ncol = 30         # grid size (number of cells)
nrow = ncol

nSlices = 400         # maximum number of frames to show in the plot
ntAnim = 2         # number of time steps for each frame

horizontalWrap = False # determines whether the flow wraps around, connecting
                       # the left and right-hand sides of the grid, or whether
                       # there's a wall there.
interpolateRotation = True # Set to False for the "Easier Method"
rotationScheme = "PlusMinus"   # "WithLatitude", "PlusMinus", "Uniform"

# Note: the rotation rate gradient is more intense than the real world, so that
# the model can equilibrate quickly.

windScheme = "Curled"  # "Curled", "Uniform"
initialPerturbation = ""    # "Tower", "NSGradient", "EWGradient"
textOutput = False
plotOutput = True
arrowScale = 30 # Adjusted for visibility, might need tuning

dT = 600 # seconds
G = 9.8e-4 # m/s2, hacked (low-G) to make it run faster
HBackground = 4000 # meters

dX = 10.E3 # meters, small enough to respond quickly. This is a very small ocean
# on a very small, low-G planet.
dY = dX # Assuming square grid cells

dxDegrees = dX / 110.e3
flowConst = G  # 1/s2
dragConst = 1.E-6  # about 10 days decay time
meanLatitude = 30 # degrees

# Here's stuff you probably won't need to change

latitude = []
rotConst = [] # Coriolis parameter f for each row
windU_force = [] # Wind forcing term for each row
for irow in range(0,nrow):
    if rotationScheme == "WithLatitude":
        latitude.append( meanLatitude + (irow - (nrow-1)/2) * dxDegrees )
        # Standard Coriolis: 2 * Omega * sin(lat) -> Omega approx 7.2921e-5 rad/s
        # Using the simplified linear ramp from template description for faster effect:
        # rotConst.append( 2 * 7.2921e-5 * math.sin(math.radians(latitude[-1]))) # s-1
        rotConst.append( -7.e-5 * math.sin(math.radians(latitude[-1]))) # Use template value
    elif rotationScheme == "PlusMinus":
         # Centered difference approx of beta = df/dy on a beta-plane
         # Template formula:
        rotConst.append( -3.5e-5 * (1. - 0.8 * ( irow - (nrow-1)/2.0 ) / float(nrow) )) # rot 50% +-
    elif rotationScheme == "Uniform":
        rotConst.append( -3.5e-5 )
    else:
        rotConst.append( 0.0 )

    if windScheme == "Curled":
        # This is likely a wind *stress* (Force/Area = rho * Cd * W^2 ~ N/m^2)
        # Needs conversion to acceleration (Force / (rho * H)) or directly use as forcing?
        # Equation is dU/dT = ... + C_wind. Assume windU is C_wind (acceleration m/s^2).
        # Template value 1e-8 seems small for acc. Let's assume it's direct acc.
        windU_force.append( 1e-8 * math.sin( (irow+0.5)/nrow * 2 * 3.14 ) )
    elif windScheme == "Uniform":
        windU_force.append( 1.e-8 )
    else:
        windU_force.append( 0.0 )

itGlobal = 0

# Array Initialization (Using template sizes)
# Note: H has ghost column ncol. U has ghost column ncol. V has ghost row nrow.
U = numpy.zeros((nrow, ncol + 1))
V = numpy.zeros((nrow + 1, ncol))
H = numpy.zeros((nrow, ncol + 1)) # Includes ghost cell H[:, ncol]

# Tendency and derivative arrays
dUdT = numpy.zeros((nrow, ncol)) # Corresponds to U[r, 1:c+1]
dVdT = numpy.zeros((nrow, ncol)) # Corresponds to V[1:r+1, c]
dHdT = numpy.zeros((nrow, ncol)) # Corresponds to H[r, 0:c]

# Spatial derivative arrays (allocate with sizes matching where derivative is evaluated)
dHdX = numpy.zeros((nrow, ncol + 1)) # Evaluated at U locations
dHdY = numpy.zeros((nrow + 1, ncol)) # Evaluated at V locations (Corrected size)
dUdX = numpy.zeros((nrow, ncol))     # Evaluated at H locations
dVdY = numpy.zeros((nrow, ncol))     # Evaluated at H locations

# Rotation terms (allocate with sizes matching where they contribute to tendency)
rotV = numpy.zeros((nrow, ncol)) # Coriolis term V*f contributing to dUdT
rotU = numpy.zeros((nrow, ncol)) # Coriolis term U*f contributing to dVdT

# --- Initial Conditions ---
midCell_row = int(nrow/2)
midCell_col = int(ncol/2)

if initialPerturbation == "Tower":
    H[midCell_row, midCell_col] = 1.0 # Height perturbation (meters)
elif initialPerturbation == "NSGradient":
    # Create gradient, e.g., high in north, low in south
    for irow in range(nrow):
       H[irow, 0:ncol] = 0.1 * (irow - (nrow-1)/2.0) # Linear gradient centered at 0
elif initialPerturbation == "EWGradient":
     # Create gradient, e.g., high in east, low in west
    for icol in range(ncol):
       H[:, icol] = 0.1 * (icol - (ncol-1)/2.0) # Linear gradient centered at 0


# Ensure initial ghost cells are consistent
if horizontalWrap:
    H[:, ncol] = H[:, 0]
else:
    H[:, ncol] = H[:, ncol-1] # Zero gradient at wall


"""
This is the work-horse subroutine. It steps forward in time, taking ntAnim steps of
duration dT.
"""
def animStep():

    global stepDump, itGlobal, U, V, H # Allow modification of state variables

    # Time Loop
    for it in range(0,ntAnim):

        # Update Ghost Cells/Boundaries before calculating derivatives
        # (Values from previous step needed for current step derivatives)
        # North/South Walls (V=0) - V boundaries are always walls
        V[0, :] = 0.0
        V[nrow, :] = 0.0

        # East/West Boundaries (U and H ghost cells)
        if horizontalWrap:
            # Periodic wrap-around
            # Ensure U ghost cell matches U boundary for derivative calc?
            # U[:, 0] is the boundary condition, U[:, ncol] is the ghost cell
            # Let's assume U[:, ncol] holds the correct value from previous step for calculations
            # Update H ghost cell from first real column H
            H[:, ncol] = H[:, 0]
        else:
            # Walls at East/West
            U[:, 0] = 0.0
            U[:, ncol] = 0.0 # Set ghost U to 0 as well
            # H ghost cell: Zero gradient at wall
            H[:, ncol] = H[:, ncol-1]


        # Calculate Longitudinal Derivatives (X-direction)
        # dHdX at U locations U[irow, icol] (icol from 0 to ncol)
        for irow in range(nrow):
            for icol in range(1, ncol + 1): # U points from index 1 to ncol
                dHdX[irow, icol] = (H[irow, icol] - H[irow, icol-1]) / dX
            # Handle U[irow, 0] boundary value for dHdX
            if horizontalWrap:
                # Use H[irow, 0] and H[irow, ncol-1] (periodic)
                dHdX[irow, 0] = (H[irow, 0] - H[irow, ncol-1]) / dX
            else: # Wall boundary
                # Force zero gradient? Or let dHdX be non-zero but U=0 dominates?
                # Set dHdX=0 consistent with U=0 wall.
                dHdX[irow, 0] = 0.0

        # dUdX at H locations H[irow, icol] (icol from 0 to ncol-1)
        for irow in range(nrow):
            for icol in range(ncol):
                dUdX[irow, icol] = (U[irow, icol+1] - U[irow, icol]) / dX


        # Calculate Latitudinal Derivatives (Y-direction)
        # dHdY at V locations V[irow, icol] (irow from 0 to nrow)
        # Need H above and below V point. Assume irow increases Northward.
        # V[irow, icol] is between H[irow-1, icol] (South) and H[irow, icol] (North)
        dHdY.fill(0.0) # Initialize, including boundaries
        for icol in range(ncol):
            for irow in range(1, nrow): # Calculate for V[1,:] to V[nrow-1,:]
                # Gradient dH/dY = (H_north - H_south) / dY
                dHdY[irow, icol] = (H[irow, icol] - H[irow-1, icol]) / dY
        # Boundaries V[0,:] and V[nrow,:] are walls, dHdY is zero there. (Already set by fill(0.0))
        # dHdY[0, :] = 0.0 (Set explicitly for clarity)
        # dHdY[nrow, :] = 0.0 (Set explicitly for clarity)


        # dVdY at H locations H[irow, icol] (irow from 0 to nrow-1)
        # Need V above (North) and below (South) H point.
        # H[irow, icol] is between V[irow, icol] (South) and V[irow+1, icol] (North)
        for irow in range(nrow):
            for icol in range(ncol):
                # Gradient dV/dY = (V_north - V_south) / dY
                dVdY[irow, icol] = (V[irow+1, icol] - V[irow, icol]) / dY


        # Calculate the Rotational Terms (Coriolis)
        if interpolateRotation:
            # "Better Method" using interpolation through H centers
            # Temporary arrays for U/V interpolated to H centers
            U_at_H = numpy.zeros((nrow, ncol))
            V_at_H = numpy.zeros((nrow, ncol))

            # 1. Interpolate U/V to H centers H[irow, icol]
            for irow in range(nrow):
                for icol in range(ncol):
                    # U at H = avg of U left and U right
                    U_at_H[irow, icol] = 0.5 * (U[irow, icol] + U[irow, icol+1])
                    # V at H = avg of V below and V above
                    V_at_H[irow, icol] = 0.5 * (V[irow, icol] + V[irow+1, icol])

            # 2. Calculate rotation terms f*U and f*V at H centers
            rotU_at_H = numpy.zeros((nrow, ncol))
            rotV_at_H = numpy.zeros((nrow, ncol))
            for irow in range(nrow):
                f = rotConst[irow] # Coriolis parameter for this row
                for icol in range(ncol):
                     rotU_at_H[irow, icol] = f * U_at_H[irow, icol]
                     rotV_at_H[irow, icol] = f * V_at_H[irow, icol]

            # 3. Back-interpolate rotational terms to U and V locations
            # rotV (term for dUdT) back-interpolated to U locations U[irow, icol]
            rotV.fill(0.0) # Clear array
            for irow in range(nrow):
                # Interpolate to interior U points U[:, 1] to U[:, ncol-1]
                for icol in range(1, ncol):
                     # Average rotV_at_H from H cells left and right of U point
                     rotV[irow, icol] = 0.5 * (rotV_at_H[irow, icol-1] + rotV_at_H[irow, icol])
                # Handle boundaries U[:, 0] and U[:, ncol]
                if horizontalWrap:
                     # U[:, 0] averages H[:, ncol-1] and H[:, 0]
                     rotV[irow, 0] = 0.5 * (rotV_at_H[irow, ncol-1] + rotV_at_H[irow, 0])
                     # rotV for U[:, ncol] uses H[:, ncol-1] and H[:, 0] -> same as U[:, 0]
                     # dUdT only calculated up to ncol-1 -> corresponds to U[:, ncol]
                     # Need rotV[irow, ncol-1] for dUdT[irow, ncol-1] (updates U[irow, ncol])
                     # Let's calculate rotV up to index ncol-1 based on the loop structure.
                else: # Walls
                     rotV[irow, 0] = 0.0 # U[:,0] is wall, rotation term is zero.
                     # rotV[irow, ncol-1] corresponds to U[:, ncol] (wall). Set to 0?
                     # The calculation loop above only goes to ncol-1. rotV[irow, ncol-1] is calculated.

            # rotU (term for dVdT) back-interpolated to V locations V[irow, icol]
            rotU.fill(0.0) # Clear array
            for icol in range(ncol):
                 # Interpolate to interior V points V[1,:] to V[nrow-1,:]
                 for irow in range(1, nrow):
                     # Average rotU_at_H from H cells below and above V point
                     rotU[irow, icol] = 0.5 * (rotU_at_H[irow-1, icol] + rotU_at_H[irow, icol])
                 # Handle boundaries V[0,:] and V[nrow,:] (walls)
                 rotU[0, icol] = 0.0 # V[0,:] is wall, rotation term is zero.
                 # rotU for V[nrow,:] (wall) should be zero.
                 # dVdT calculated up to nrow-1 -> corresponds to V[nrow, :]
                 # Need rotU[nrow-1, icol] for dVdT[nrow-1, icol] (updates V[nrow, icol])
                 # The loop above only goes to nrow-1. rotU[nrow-1, icol] is calculated.

        else:
            # "Easier Method" (No interpolation, direct indexing - prone to noise)
             for irow in range(nrow):
                 f = rotConst[irow]
                 for icol in range(ncol): # Loop over dUdT/dVdT points
                     # rotV contributes to dUdT[irow, icol] (which updates U[irow, icol+1])
                     # Simplest: use V value "at" [irow, icol]
                     # Need V[irow, icol] and maybe V[irow+1, icol]? Let's use V[irow, icol] as per example matching.
                     rotV[irow, icol] = f * V[irow, icol]

                     # rotU contributes to dVdT[irow, icol] (which updates V[irow+1, icol])
                     # Simplest: use U value "at" [irow, icol]
                     # Need U[irow, icol] and maybe U[irow, icol+1]? Let's use U[irow, icol] as per example matching.
                     rotU[irow, icol] = f * U[irow, icol]


        # Assemble the Time Derivatives (dUdT, dVdT, dHdT)
        # Loop over the grid points where tendencies are calculated (nrow, ncol)
        for irow in range(nrow):
            wind_effect = windU_force[irow] # Wind acceleration for this row

            for icol in range(ncol):
                # dHdT at H[irow, icol]
                # Equation: dH/dT = - H_background * ( dU/dX + dV/dY )
                dHdT[irow, icol] = -HBackground * (dUdX[irow, icol] + dVdY[irow, icol])

                # dUdT corresponds to U[irow, icol+1]
                # Equation: dU/dT = f*V - g*dH/dX - drag*U + wind
                U_target_col = icol + 1
                # Use rotV calculated for this point (index [irow, icol])
                # Use dHdX at the U location U[irow, icol+1] (index [irow, icol+1])
                # Use U at the location U[irow, icol+1]
                dUdT[irow, icol] = rotV[irow, icol] \
                                   - flowConst * dHdX[irow, U_target_col] \
                                   - dragConst * U[irow, U_target_col] \
                                   + wind_effect

                # dVdT corresponds to V[irow+1, icol]
                # Equation: dV/dT = -f*U - g*dH/dY - drag*V
                V_target_row = irow + 1
                # Use rotU calculated for this point (index [irow, icol])
                # Use dHdY at the V location V[irow+1, icol] (index [irow+1, icol])
                # Use V at the location V[irow+1, icol]
                dVdT[irow, icol] = -rotU[irow, icol] \
                                   - flowConst * dHdY[V_target_row, icol] \
                                   - dragConst * V[V_target_row, icol]


        # Step Forward One Time Step (using Forward Euler method)
        # Update H (main grid cells H[0..nrow-1, 0..ncol-1])
        H[:, 0:ncol] += dHdT * dT

        # Update U (interior points U[0..nrow-1, 1..ncol])
        # dUdT[irow, icol] updates U[irow, icol+1]
        U[:, 1:ncol+1] += dUdT * dT

        # Update V (interior points V[1..nrow, 0..ncol-1])
        # dVdT[irow, icol] updates V[irow+1, icol]
        V[1:nrow+1, :] += dVdT * dT

        # Note: Boundary/ghost cells for U, V, H are updated at the start
        # of the *next* time step calculation or at the end of this function
        # Let's add final boundary enforcement here too for clarity.

        # Final enforcement of boundary conditions after update
        V[0, :] = 0.0
        V[nrow, :] = 0.0
        if not horizontalWrap:
            U[:, 0] = 0.0
            U[:, ncol] = 0.0 # Ensure ghost cell is also zero for wall

    # Increment global time step counter
    itGlobal = itGlobal + ntAnim
# --- End of animStep function ---

# --- Plotting and Execution Logic (from template) ---
def firstFrame():
    global fig, ax, hPlot, quiv, quiv2 # Make plot objects global
    fig, ax = plt.subplots(figsize=(6,6)) # Adjust figure size if needed
    ax.set_title("H (color) and U/V (arrows)")
    # Display H excluding the ghost column
    hh = H[:,0:ncol]
    loc = tkr.MultipleLocator(base=1.0) # Ensure grid lines are between cells
    ax.xaxis.set_major_locator(loc)
    ax.yaxis.set_major_locator(loc)
    ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=1)
    # Set extent to align grid lines correctly with cell boundaries
    extent = [-0.5, ncol - 0.5, nrow - 0.5, -0.5] # left, right, bottom, top
    hPlot = ax.imshow(hh, interpolation='nearest', cmap='viridis',
                      clim=(-0.5, 0.5), extent=extent, origin='upper')
    plt.colorbar(hPlot, ax=ax, label="Height Anomaly H (m)")
    quiv = None # Initialize quiver objects
    quiv2 = None
    plotArrows() # Initial arrows
    plt.tight_layout()
    plt.show(block=False)
    plt.pause(0.1) # Pause briefly to ensure drawing

def plotArrows():
    global quiv, quiv2, ax # Need ax to plot on

    # Remove old arrows if they exist
    if quiv is not None:
        quiv.remove()
    if quiv2 is not None:
        quiv2.remove()

    # --- Plot U arrows (at vertical cell faces) ---
    # U[irow, icol] is on the face LEFT of H[irow, icol]
    # Arrow origin: (icol-0.5, irow)
    xx_u, yy_u, uu_u, vv_u = [], [], [], []
    for irow in range(nrow):
        for icol in range(ncol + 1): # Include boundary/ghost U
            # Center arrow vertically in cell row
            arrow_y = irow
            # Center arrow horizontally on the vertical face
            arrow_x = icol - 0.5
            xx_u.append(arrow_x)
            yy_u.append(arrow_y)
            uu_u.append(U[irow, icol] * arrowScale) # U component
            vv_u.append(0)                          # V component is zero for U arrows

    # --- Plot V arrows (at horizontal cell faces) ---
    # V[irow, icol] is on the face TOP of H[irow, icol]
    # Arrow origin: (icol, irow-0.5)
    xx_v, yy_v, uu_v, vv_v = [], [], [], []
    for irow in range(nrow + 1): # Include boundary/ghost V
        for icol in range(ncol):
            # Center arrow horizontally in cell column
            arrow_x = icol
            # Center arrow vertically on the horizontal face
            arrow_y = irow - 0.5
            xx_v.append(arrow_x)
            yy_v.append(arrow_y)
            uu_v.append(0)                          # U component is zero for V arrows
            vv_v.append(V[irow, icol] * arrowScale) # V component (Plotting V positive downwards)
                                                   # Correction: imshow origin='upper', so positive y is down.
                                                   # V positive is Northward (up). So use -V ?
                                                   # Let's plot V directly. Arrow points in direction of V.

    # Create quiver plots (use scale=1 if arrow components are already scaled by arrowScale)
    # Important: Set pivot='tip' or 'tail' or 'middle' for arrow origin consistency
    quiv = ax.quiver(xx_u, yy_u, uu_u, vv_u, color='white', scale=1, units='xy', pivot='middle', headwidth=5)
    quiv2 = ax.quiver(xx_v, yy_v, uu_v, vv_v, color='black', scale=1, units='xy', pivot='middle', headwidth=5)


def updateFrame():
    global fig, ax, hPlot, quiv, quiv2 # Use global plot objects
    # Update H plot data (excluding ghost column)
    hh = H[:,0:ncol]
    hPlot.set_array(hh)
    # Update arrows by removing old ones and plotting new ones
    plotArrows()
    # Redraw the figure
    fig.canvas.draw_idle() # Use draw_idle for smoother updates
    fig.canvas.flush_events() # Process events
    # print("Time: ", math.floor( itGlobal * dT / 86400.*10)/10, "days") # Original print
    print(f"Time Step: {itGlobal}, Model Time: {itGlobal * dT / 86400.0:.2f} days")


def textDump():
    print(f"\n--- Time Step: {itGlobal} ---")
    numpy.set_printoptions(precision=3, suppress=True) # Pretty print
    print("H (incl. ghost):\n", H)
    # print("dHdX (at U points):\n", dHdX) # Often useful for debugging
    # print("dHdY (at V points):\n", dHdY) # Often useful for debugging
    print("U (incl. ghost):\n", U)
    # print("dUdX (at H points):\n", dUdX) # Often useful for debugging
    print("V (incl. ghost):\n", V)
    # print("dVdY (at H points):\n", dVdY) # Often useful for debugging
    # if interpolateRotation:
    #     print("rotV (interp, for dUdT):\n", rotV)
    #     print("rotU (interp, for dVdT):\n", rotU)
    # else:
    #     print("rotV (simple, for dUdT):\n", rotV)
    #     print("rotU (simple, for dVdT):\n", rotU)
    # print("dHdT (at H points):\n", dHdT)
    # print("dUdT (for U points 1..ncol):\n", dUdT)
    # print("dVdT (for V points 1..nrow):\n", dVdT)
    print("-"*(ncol*8)) # Separator


# --- Main Execution ---
if textOutput:
    textDump() # Show initial state
if plotOutput:
    firstFrame()

for i_anim_step in range(nSlices):
    animStep() # Run the simulation step(s)
    if textOutput:
        textDump() # Print state after step
    if plotOutput:
        updateFrame() # Update plot

# Keep plot open at the end if showing plots
if plotOutput:
    print("Simulation finished. Close the plot window to exit.")
    plt.show(block=True)
