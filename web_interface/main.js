console.log("main.js loaded");

// --- Global Variables ---
let pyodide = null;
let pythonScript = null;
let simulationState = { params: {}, isRunning: false, currentStep: 0, totalStepsInFrame: 5, animationFrameId: null };
let defaultParamsFromPython = {};
let pyodideReady = false; // Flag for full init complete

// --- DOM Elements ---
const runButton = document.getElementById('runButton');
const stopButton = document.getElementById('stopButton');
const statusElem = document.getElementById('status');
const simTimeElem = document.getElementById('simTime');
const chartElem = document.getElementById('plotlyChart');
// Input elements
const ncolInput = document.getElementById('ncol');
const perturbationSelect = document.getElementById('initialPerturbation');
const dtInput = document.getElementById('dT');
const rotationSchemeSelect = document.getElementById('rotationScheme');
const horizontalWrapCheckbox = document.getElementById('horizontalWrap');
const dxInput = document.getElementById('dX');
const gInput = document.getElementById('G');
const hBackgroundInput = document.getElementById('HBackground');
const dragConstInput = document.getElementById('dragConst');
const meanLatitudeInput = document.getElementById('meanLatitude');
const betaFactorInput = document.getElementById('betaFactor');
const ntAnimInput = document.getElementById('ntAnim');
const windSchemeSelect = document.getElementById('windScheme');
// Viz Controls
const fixedScaleCheckbox = document.getElementById('fixedScaleCheckbox');
const zminInput = document.getElementById('zminInput');
const zmaxInput = document.getElementById('zmaxInput');
// Modal elements
const helpButton = document.getElementById('helpButton');
const helpModal = document.getElementById('helpModal');
const helpModalOverlay = document.getElementById('helpModalOverlay');
const closeHelpModalButton = document.getElementById('closeHelpModal');
const closeHelpModalButtonBottom = document.getElementById('closeHelpModalBottom');


// --- Pyodide Initialization ---
async function loadPyodideAndRunPython() {
    statusElem.textContent = 'Status: Loading Pyodide engine...';
    console.log('Loading Pyodide...');
    try {
        if (typeof window.loadPyodide !== 'function') throw new Error("loadPyodide function not found.");
        pyodide = await window.loadPyodide();
        console.log('Pyodide loaded.');

        statusElem.textContent = 'Status: Loading NumPy...';
        await pyodide.loadPackage("numpy");
        console.log('NumPy loaded.');

        statusElem.textContent = 'Status: Fetching script...';
        const response = await fetch('shallow_water_model.py');
        if (!response.ok) throw new Error(`HTTP error fetching script! status: ${response.status}`);
        pythonScript = await response.text();
        console.log('Script fetched.');

        statusElem.textContent = 'Status: Executing script...';
        await pyodide.runPythonAsync(pythonScript);
        console.log('Script executed.');

        const defaultsProxy = pyodide.globals.get('DEFAULT_PARAMS');
        if (!defaultsProxy) throw new Error("Could not get DEFAULT_PARAMS from Python.");
        defaultParamsFromPython = defaultsProxy.toJs();
        defaultsProxy.destroy();
        if (defaultParamsFromPython instanceof Map) {
             defaultParamsFromPython = Object.fromEntries(defaultParamsFromPython);
        }
        console.log("Defaults fetched:", defaultParamsFromPython);
         if (typeof defaultParamsFromPython.dX !== 'number' || typeof defaultParamsFromPython.ncol !== 'number' || typeof defaultParamsFromPython.dT !== 'number') {
                console.warn("Critical default parameters (dX, ncol, dT) might be missing or invalid!");
         }

        simulationState.params = { ...defaultParamsFromPython };
        updateUIFromParams(simulationState.params);

        statusElem.textContent = 'Status: Ready.';
        runButton.disabled = false;
        setControlsEnabled(true);
        pyodideReady = true; // <-- Set readiness flag

    } catch (error) {
        console.error('Initialization failed:', error);
        statusElem.textContent = `Status: Error loading: ${error.message}`;
        runButton.disabled = true; stopButton.disabled = true; setControlsEnabled(false);
        pyodideReady = false; // Ensure flag is false on error
    }
}

// --- Simulation Control ---
function startSimulation() {
    // --- Add check for pyodide readiness ---
    if (simulationState.isRunning || !pyodideReady) {
         console.warn("startSimulation called before Pyodide was fully ready or while already running.");
         return;
    }
    console.log('%cStartSimulation: Function entered.', 'color: blue; font-weight: bold;');
    statusElem.textContent = 'Status: Initializing simulation...';
    simulationState.isRunning = true;
    runButton.disabled = true;
    stopButton.disabled = false;
    setControlsEnabled(false);

    console.log('StartSimulation: Reading parameters from UI...');
    readParametersFromUI();
    console.log('StartSimulation: Parameters set:', JSON.stringify(simulationState.params, null, 2));

    Plotly.purge(chartElem);

    console.log("StartSimulation: Attempting Python initialization...");
    try {
        console.log('StartSimulation: Setting current_params_js in Pyodide...');
        pyodide.globals.set('current_params_js', simulationState.params);
        console.log('StartSimulation: current_params_js set.');

        console.log('StartSimulation: Executing Python init block...');
        const pythonInitCode = `
import numpy as np
#print('Python: Running initialization block...')
current_params = current_params_js.to_py() if hasattr(current_params_js, 'to_py') else dict(current_params_js)
#print(f'Python: Params received: dX={current_params.get("dX")}, dT={current_params.get("dT")}')

U_global_py, V_global_py, H_global_py = initialize_state(current_params)
coriolis_f_global_py, _ = calculate_coriolis_parameter(current_params)
wind_U_global_py = calculate_wind_forcing(current_params)

# Results U_global_py etc are now in Python's global scope
#print('Python: Initialization complete.')
`;
        pyodide.runPython(pythonInitCode); // Run sync
        console.log("StartSimulation: Python initialization block executed.");

        // --- Verification ---
        const requiredGlobals = ['U_global_py', 'V_global_py', 'H_global_py', 'coriolis_f_global_py', 'wind_U_global_py'];
        const missingGlobals = requiredGlobals.filter(key => !pyodide.globals.has(key));
        if (missingGlobals.length > 0) {
             throw new Error(`Python initialization did not set expected global variables. Missing: ${missingGlobals.join(', ')}.`);
        }
        console.log("StartSimulation: Verified Python globals seem to be set.");

        simulationState.currentStep = 0;
        updateSimTimeDisplay();

        console.log("StartSimulation: Plotting initial state...");
        plotCurrentState();
        console.log("StartSimulation: Initial state plotted.");

        statusElem.textContent = 'Status: Running...';
        console.log("StartSimulation: Requesting first animation frame...");
        if (simulationState.animationFrameId) cancelAnimationFrame(simulationState.animationFrameId);
        simulationState.animationFrameId = requestAnimationFrame(animationLoop);
        console.log("StartSimulation: Animation loop initiated.");

    } catch(error) {
        console.error("!!!!!!!!! Error during simulation initialization !!!!!!!!!", error);
        console.error("Full Error Object:", error);
        statusElem.textContent = `Status: Error initializing: ${error.message?.split('\n')[0] ?? error}`;
        stopSimulation();
    }
}


function stopSimulation() {
    if (!simulationState.isRunning && simulationState.animationFrameId === null) return;
    console.log('%cStopping simulation...', 'color: red; font-weight: bold;');
    simulationState.isRunning = false;
    if (simulationState.animationFrameId) {
        cancelAnimationFrame(simulationState.animationFrameId);
        simulationState.animationFrameId = null;
        console.log('Animation frame cancelled.');
    }
    runButton.disabled = false; stopButton.disabled = true;
    statusElem.textContent = 'Status: Stopped.';
    setControlsEnabled(true);
}

// *** Refined animationLoop with Detailed Logging ***
function animationLoop() {
    if (!simulationState.isRunning) {
        console.log("animationLoop: exiting because isRunning is false.");
        simulationState.animationFrameId = null; return;
    }
    // Add a timestamp for performance debugging if needed later
    // const frameStartTime = performance.now();
    console.log(`%c--- animationLoop: Frame Start (Sim Step ${simulationState.currentStep}) ---`, 'color: green');

    try {
        // 1. Prepare Python execution environment
        // console.log("animationLoop: Setting JS params for Python..."); // Can be verbose
        pyodide.globals.set('current_params_js', simulationState.params);
        pyodide.globals.set('n_steps_to_run', simulationState.totalStepsInFrame);

        // Verify necessary globals exist *before* running step code
        const requiredGlobals = ['U_global_py', 'V_global_py', 'H_global_py', 'coriolis_f_global_py', 'wind_U_global_py'];
         if (requiredGlobals.some(key => !pyodide.globals.has(key))) {
            throw new Error(`Cannot run step. Missing required Python globals.`);
         }
         // console.log("animationLoop: Verified required Python globals exist."); // Can be verbose

         // === ADDED: Log a sample value BEFORE stepping ===
         let h_center_before = 'N/A';
         try {
            const h_proxy_before = pyodide.globals.get('H_global_py');
            if(h_proxy_before) {
                const r = Math.floor(simulationState.params.ncol / 2);
                const c = Math.floor(simulationState.params.ncol / 2);
                 // Accessing single element might be faster than .toJs()
                if (r < h_proxy_before.shape[0] && c < h_proxy_before.shape[1]) {
                     h_center_before = h_proxy_before.get([r,c]); // Get value at center
                     if(typeof h_center_before?.toFixed === 'function') h_center_before = h_center_before.toFixed(6);
                }
                h_proxy_before.destroy(); // Destroy proxy immediately
            }
         } catch (e) { console.warn("Could not get H value before step:", e); }
         console.log(`animationLoop: H[center] BEFORE python step: ${h_center_before}`);
         // === END ADDED LOG ===

        // 2. Run Simulation Steps in Python
        // console.log("animationLoop: Executing pythonStepCode..."); // Can be verbose
        const pythonStepCode = `
current_params = current_params_js.to_py() if hasattr(current_params_js, 'to_py') else dict(current_params_js)
coriolis_f = coriolis_f_global_py
wind_U = wind_U_global_py
U = U_global_py
V = V_global_py
H = H_global_py

for i in range(n_steps_to_run):
    U, V, H = step_rk4(U, V, H, current_params, coriolis_f, wind_U)

# Assign results back to the global Pyodide variables
U_global_py = U
V_global_py = V
H_global_py = H
# Assigning back ensures the global references are updated if step_rk4 returned new objects
        `;
        pyodide.runPython(pythonStepCode);
        // console.log("animationLoop: Python execution finished."); // Can be verbose

        simulationState.currentStep += simulationState.totalStepsInFrame;

        // === ADDED: Log a sample value AFTER stepping ===
         let h_center_after = 'N/A';
         try {
            const h_proxy_after = pyodide.globals.get('H_global_py');
             if(h_proxy_after) {
                const r = Math.floor(simulationState.params.ncol / 2);
                const c = Math.floor(simulationState.params.ncol / 2);
                if (r < h_proxy_after.shape[0] && c < h_proxy_after.shape[1]) {
                    h_center_after = h_proxy_after.get([r,c]);
                    if(typeof h_center_after?.toFixed === 'function') h_center_after = h_center_after.toFixed(6);
                }
                h_proxy_after.destroy();
            }
         } catch (e) { console.warn("Could not get H value after step:", e); }
         console.log(`animationLoop: H[center] AFTER python step: ${h_center_after}`);
         // === END ADDED LOG ===


        // 3. Plot the Updated State
        // console.log("animationLoop: Calling plotCurrentState..."); // Can be verbose
        plotCurrentState(); // This reads the updated H_global_py

        // 4. Update Time Display
        updateSimTimeDisplay();

        // 5. Request Next Frame if still running
        if (simulationState.isRunning) {
             simulationState.animationFrameId = requestAnimationFrame(animationLoop);
        } else {
            console.log("animationLoop: isRunning is false after plot/update, stopping loop.");
            simulationState.animationFrameId = null;
        }

    } catch(error) {
        console.error("!!!!!!!!! Error during simulation step (animationLoop) !!!!!!!!!", error);
        if (error.stack) { console.error("JS Stack Trace:", error.stack); }
        if (error.name === 'PythonError' && error.message) { console.error("Python Traceback:\n" + error.message); }
        statusElem.textContent = `Status: Error running step: ${error.message?.split('\n')[0] ?? error}`;
        stopSimulation();
    }
    // const frameEndTime = performance.now();
    // console.log(`%c--- animationLoop: Frame End (Took: ${(frameEndTime-frameStartTime).toFixed(1)}ms) ---`, 'color: green'); // Optional performance log
}


// --- UI Interaction (Keep previous corrected versions) ---
function readParametersFromUI() { /* ... From previous version ... */
    const params = simulationState.params;
    const defaultParams = defaultParamsFromPython || {};
    const safeParseInt = (element, defaultValue) => element ? (parseInt(element.value) || defaultValue) : defaultValue;
    const safeParseFloat = (element, defaultValue) => element ? (parseFloat(element.value) || defaultValue) : defaultValue;

    params.ncol = safeParseInt(ncolInput, defaultParams.ncol);
    params.dX = safeParseFloat(dxInput, defaultParams.dX / 1000.0) * 1000.0;
    params.dT = safeParseFloat(dtInput, defaultParams.dT);
    params.ntAnim = safeParseInt(ntAnimInput, defaultParams.ntAnim);
    params.initialPerturbation = perturbationSelect.value || defaultParams.initialPerturbation;
    params.G = safeParseFloat(gInput, defaultParams.G * 1e4) / 1e4;
    params.HBackground = safeParseFloat(hBackgroundInput, defaultParams.HBackground);
    params.dragConst = safeParseFloat(dragConstInput, defaultParams.dragConst * 1e7) / 1e7;
    params.horizontalWrap = horizontalWrapCheckbox.checked;
    params.windScheme = windSchemeSelect.value || defaultParams.windScheme;
    params.rotationScheme = rotationSchemeSelect.value || defaultParams.rotationScheme;
    params.meanLatitude = safeParseFloat(meanLatitudeInput, defaultParams.meanLatitude);
    params.betaFactor = safeParseFloat(betaFactorInput, defaultParams.betaFactor);

    if (params.ncol <= 0) params.ncol = defaultParams.ncol ?? 32;
    if (params.dX <= 0) params.dX = defaultParams.dX ?? 20000.0;
    if (params.dT <= 0) params.dT = defaultParams.dT ?? 150;
    if (params.ntAnim <= 0) params.ntAnim = defaultParams.ntAnim ?? 5;
    if (params.HBackground <= 0) params.HBackground = defaultParams.HBackground ?? 1000.0;
    if (params.dragConst < 0) params.dragConst = defaultParams.dragConst ?? 5e-7;

    params.interpolateRotation = defaultParams.interpolateRotation ?? true;
    params.useFixedScale = fixedScaleCheckbox.checked;
    params.fixedZmin = safeParseFloat(zminInput, -0.5);
    params.fixedZmax = safeParseFloat(zmaxInput, 0.5);
    if (isNaN(params.fixedZmin)) params.fixedZmin = -0.5;
    if (isNaN(params.fixedZmax)) params.fixedZmax = 0.5;

    simulationState.totalStepsInFrame = params.ntAnim;
}

function updateUIFromParams(params) { /* ... From previous version ... */
     if (!params) return;
     ncolInput.value = params.ncol;
     dxInput.value = params.dX / 1000.0;
     dtInput.value = params.dT;
     ntAnimInput.value = params.ntAnim;
     perturbationSelect.value = params.initialPerturbation;
     gInput.value = params.G * 1e4;
     hBackgroundInput.value = params.HBackground;
     dragConstInput.value = params.dragConst * 1e7;
     horizontalWrapCheckbox.checked = params.horizontalWrap;
     windSchemeSelect.value = params.windScheme;
     rotationSchemeSelect.value = params.rotationScheme;
     meanLatitudeInput.value = params.meanLatitude;
     betaFactorInput.value = params.betaFactor;
     const useFixed = params.useFixedScale ?? false;
     fixedScaleCheckbox.checked = useFixed;
     zminInput.value = params.fixedZmin ?? -0.5;
     zmaxInput.value = params.fixedZmax ?? 0.5;
     zminInput.disabled = !useFixed;
     zmaxInput.disabled = !useFixed;
     console.log("UI updated from params");
 }

function setControlsEnabled(isEnabled) { /* ... From previous version ... */
    const container = document.querySelector('.main-container');
    if (!container) return;
    const controls = container.querySelectorAll('input:not(#runButton):not(#stopButton), select');
    controls.forEach(input => {
        let shouldBeDisabled = !isEnabled;
        if (input.id === 'zminInput' || input.id === 'zmaxInput') {
            if (!fixedScaleCheckbox.checked) { shouldBeDisabled = true; }
        }
        input.disabled = shouldBeDisabled;
    });
    console.log(`Controls ${isEnabled ? 'Enabled' : 'Disabled'}`);
}

function updateSimTimeDisplay() { /* ... From previous version ... */
    const currentDT = simulationState.params?.dT || defaultParamsFromPython?.dT || 150;
    const model_time_days = simulationState.currentStep * currentDT / 86400.0;
    simTimeElem.textContent = `Simulation Time: ${model_time_days.toFixed(2)} days`;
}


// --- Plotting (Keep previous corrected version) ---
function plotCurrentState() { /* ... From previous version ... */
    if (!pyodide || typeof pyodide.globals?.get !== 'function' || !pyodide.globals.has('H_global_py')) { return; }
    try {
        const H_proxy = pyodide.globals.get('H_global_py');
        if (!H_proxy) { return; }
        // Optimization: Check if shape exists before accessing potentially large data
        const H_shape = H_proxy.shape;
        if (!H_shape || H_shape.length !== 2) { H_proxy.destroy(); return;}
        // Proceed with conversion and plotting
        const H_js = H_proxy.toJs({create_proxies: false});
        H_proxy.destroy();

        if (!H_js || H_js.length === 0) { return; }
        const ncol = simulationState.params?.ncol;
        if (typeof ncol !== 'number' || ncol <= 0) { return; }
        const H_plot_data = H_js.map(row => row.slice(0, ncol));

        let zmin_final, zmax_final;
        const useFixed = simulationState.params?.useFixedScale;
        const fixedZmin = simulationState.params?.fixedZmin;
        const fixedZmax = simulationState.params?.fixedZmax;

        if (useFixed && typeof fixedZmin === 'number' && typeof fixedZmax === 'number') {
            zmin_final = fixedZmin; zmax_final = fixedZmax;
        } else {
             let h_min_plot = Infinity, h_max_plot = -Infinity; let hasFinite = false;
            for (const row of H_plot_data) { if (!row || typeof row.slice !== 'function') continue; for (const val of row) { if(isFinite(val)){ if (val < h_min_plot) h_min_plot = val; if (val > h_max_plot) h_max_plot = val; hasFinite = true; } } }
            const fallbackRange = 0.1;
            if (!hasFinite || Math.abs(h_max_plot - h_min_plot) < 1e-9) { zmin_final = -fallbackRange; zmax_final = fallbackRange; }
            else { const h_range = Math.max(fallbackRange, h_max_plot - h_min_plot); const h_center = (h_max_plot + h_min_plot) / 2.0; const range_multiplier = 0.7; zmin_final = h_center - h_range * range_multiplier; zmax_final = h_center + h_range * range_multiplier; }
        }

        const trace = { z: H_plot_data, type: 'heatmap', colorscale: 'YlGnBu', reversescale: true, zmin: zmin_final, zmax: zmax_final, zsmooth: 'best', colorbar: { title: 'Height (m)', titleside: 'right' } };
        const layout = { title: `Shallow Water Height Anomaly (Step ${simulationState.currentStep})`, yaxis: { scaleanchor: "x", scaleratio: 1, title: 'Y Index', autorange: 'reversed' }, xaxis: { constrain: "domain", title: 'X Index' }, margin: { l: 60, r: 20, b: 50, t: 50, pad: 4 } };
        Plotly.react(chartElem, [trace], layout);

    } catch (error) { console.error("!!!!!!!!! Error during plotting !!!!!!!!!", error); statusElem.textContent = `Status: Error plotting: ${error.message?.split('\n')[0] ?? error}`; stopSimulation(); }
}

// --- Modal Functionality ---
function openHelpModal() {
    if(helpModal && helpModalOverlay) {
        helpModalOverlay.classList.remove('hidden');
        helpModal.classList.remove('hidden');
    }
}
function closeHelpModal() {
     if(helpModal && helpModalOverlay) {
        helpModalOverlay.classList.add('hidden');
        helpModal.classList.add('hidden');
    }
}

// --- Event Listeners (Keep previous corrected versions) ---
runButton.addEventListener('click', () => {
    console.log('%cRun button CLICKED!', 'color: orange; font-weight: bold;');
    startSimulation();
});
stopButton.addEventListener('click', stopSimulation);
fixedScaleCheckbox.addEventListener('change', () => { /* ... */
    const isChecked = fixedScaleCheckbox.checked;
    zminInput.disabled = !isChecked;
    zmaxInput.disabled = !isChecked;
    readParametersFromUI();
    if (!simulationState.isRunning && simulationState.currentStep > 0) { plotCurrentState(); }
    else if (simulationState.isRunning){ console.log("Fixed scale checkbox changed while running..."); }
});
helpButton.addEventListener('click', openHelpModal);
closeHelpModalButton.addEventListener('click', closeHelpModal);
closeHelpModalButtonBottom.addEventListener('click', closeHelpModal);
helpModalOverlay.addEventListener('click', closeHelpModal);


// --- Initial Load Sequence (Keep as is) ---
document.addEventListener('DOMContentLoaded', () => { /* ... calls initializeApp which calls loadPyodideAndRunPython ... */
    console.log("DOM Content Loaded. Initializing App...");
    runButton.disabled = true; stopButton.disabled = true; setControlsEnabled(false);
    statusElem.textContent = 'Status: Waiting for Pyodide...';
    initializeApp(); // Make sure initializeApp calls loadPyodideAndRunPython
});

function initializeApp() { // Ensure this function exists from previous debug step
    console.log("initializeApp called.");
    let attempt = 0; const maxAttempts = 10; const interval = 500;
    function tryLoad() {
        if (typeof window.loadPyodide === 'function') {
            console.log("loadPyodide function found. Starting load process.");
            loadPyodideAndRunPython();
        } else if (attempt < maxAttempts) {
            attempt++; console.log(`loadPyodide not ready, attempt ${attempt}/${maxAttempts}. Retrying...`);
            setTimeout(tryLoad, interval);
        } else {
            console.error("loadPyodide failed after multiple attempts.");
            statusElem.textContent = 'Status: Error - Failed to load Pyodide library.';
        }
    }
    setTimeout(tryLoad, 100);
}