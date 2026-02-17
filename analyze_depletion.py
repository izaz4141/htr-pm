import openmc
import openmc.deplete
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import re
import pickle
import time
import datetime
import csv

mpl.use('agg')
openmc.config["chain_file"] = (
    "/home/glicole/Documents/OpenMC/endfb71/chain_endfb71_pwr.xml"
)
plt.rcParams.update(
    {
        "font.size": 46,
        "figure.titlesize": 52,
        "figure.titleweight": "bold",
        "axes.titlesize": 52,
        "axes.titleweight": "bold",
        "axes.labelpad": 36,
        "axes.grid": True,
        "lines.linewidth": 3.5,
        "grid.linestyle": "--",
        "axes.grid.which": "both",
        "grid.alpha": 0.7,
        "xtick.minor.visible": True,
        "xtick.major.size": 15,
        "xtick.minor.size": 9,
        "xtick.major.width": 3,
        "xtick.minor.width": 2,
        "xtick.labelsize": 36,
        "ytick.minor.visible": True,
        "ytick.major.size": 15,
        "ytick.minor.size": 9,
        "ytick.major.width": 3,
        "ytick.minor.width": 2,
        "ytick.labelsize": 36,
        "legend.fontsize": 36,
    }
)
cwd = os.getcwd()
region = 60
results = None
results_last = None
replace = 0
if "depletion" in cwd:
    mode = "deplete"
elif "eigen" in cwd:
    mode = "eigenvalue"
elif "combine" in cwd:
    mode = "combine"

iss_arr = ["m_u235", "m_u238", "m_pu239", "m_pu240", "m_pu241",
            "m_np237", "m_am241", "m_cm244", "a_xe135", "a_sm149",
            "fission_rate", "n2_rate", "n3_rate", "n4_rate",
            "alpha_rate", "capture_rate", "proton_rate"]

def get_step_numbers():
    # Regex pattern to match "step_{number}"
    step_pattern = re.compile(r"^step_(\d+)$") # Regex to match directory names like "step_123"
    sim_pattern = "openmc_simulation_n" # Pattern for simulation files
    local_step_numbers = []
    local_sum_steps = 0
    step_sim_counts = {}

    # Iterate through all items in the current directory
    for item in os.listdir(cwd):
        item_path = os.path.join(cwd, item) # Get the full path of the item
        if os.path.isdir(item_path):
            match = step_pattern.match(item) # Try to match the 'step_X' pattern
            if match and 0 <= int(match.group(1)) <= 663:
                step_num = int(match.group(1)) # Extract the step number
                local_step_numbers.append(step_num)
                
                current_step_sim_count = 0
                # Count simulation files within this step directory
                for sim_file in os.listdir(item_path):
                    if sim_pattern in sim_file:
                        current_step_sim_count += 1
                
                step_sim_counts[step_num] = current_step_sim_count # Store count for this step
                local_sum_steps += current_step_sim_count # Add to total simulation count
                
    return sorted(local_step_numbers), local_sum_steps, step_sim_counts


def find_segment_indices(time):
    diffs = [time[i + 1] - time[i] for i in range(len(time) - 1)]
    threshold = 1
    indices = []
    for i, diff in enumerate(diffs):
        if diff < threshold:
            indices.append(i)
    indices.append(len(time) - 1)
    return indices

def truncate_2d_arr():
    for key in iss_arr:
        if key in results and results[key].shape[0] == region:
            results[key] = results[key][:, :sum_steps]
        else:
            print(f"Warning: Key '{key}' not found or unexpected shape in results. Skipping extension for this key.")


gen_time = time.perf_counter()
if mode == "deplete":
    step_numbers, sum_steps, step_sim_counts_global = get_step_numbers()
    
def get_results():
    if mode == "combine":
        return

    global results, results_last, sum_steps # Declare global to modify these variables
    
    if os.path.exists(f"{cwd}/results.pkl"):
        # Load the data from pickle
        with open(f"{cwd}/results.pkl", "rb") as f:
            results = pickle.load(f)
            print("Data Loaded !!")
    
    if mode == "deplete":
        sum_step = 0
        steps_to_process_for_data = []

        if results: # Continue results
            # Get the set of step numbers already present in the loaded results
            existing_step_numbers = results.get("step_numbers", [])
            last_processed_step = existing_step_numbers[-1]

            if replace > 0:
                # try to find the index of `replace` in the list
                try:
                    idx = existing_step_numbers.index(replace)
                except ValueError:
                    idx = len(existing_step_numbers)
                if replace < last_processed_step: 
                    last_processed_step = replace - 1
                    results["step_numbers"] = existing_step_numbers[:idx]

            for step_num in step_numbers:
                if step_num > last_processed_step:
                    # If a step number is new, it needs to be processed for data.
                    steps_to_process_for_data.append(step_num)
                    results['step_numbers'].append(step_num)
                else:
                    # If the step already exists in loaded_results, we need to advance
                    # the current_sim_index by the number of simulations in that step.
                    # This ensures correct indexing for new data that will be appended.
                    sum_step += step_sim_counts_global.get(step_num, 0)

            old_array_size = results["m_u235"].shape[1] if results["m_u235"].size > 0 else 0
            if replace > 0:
                results["keff"] = results["keff"][:sum_step]
                results["keff_unc"] = results["keff_unc"][:sum_step]
                results["times"] = results["times"][:sum_step]
                truncate_2d_arr()
                if "heating" in results:
                    results["heating"] = results["heating"][:sum_steps]
                if "spectrum" in results and results["spectrum"].shape[1] == 3:
                    results["spectrum"] = results["spectrum"][:sum_steps, :]
                if "reaction" in results:
                    results["reaction"] = results["reaction"][:sum_steps]


            if sum_steps > old_array_size:
                slots_to_add = sum_steps - old_array_size
                print(f"Adding {slots_to_add} new simulation slots to arrays.")

                # Extend 2D arrays (e.g., region x sum_steps) by concatenating zeros.
                for key in iss_arr:
                    if key in results and results[key].shape[0] == region:
                        new_part = np.zeros((region, slots_to_add))
                        results[key] = np.concatenate((results[key], new_part), axis=1)
                    else:
                        print(f"Warning: Key '{key}' not found or unexpected shape in results. Skipping extension for this key.")

                # Extend 1D numpy arrays (e.g., heating) by concatenating zeros.
                if "heating" in results:
                    new_part = np.zeros(slots_to_add)
                    results["heating"] = np.concatenate((results["heating"], new_part))
                else:
                    print("Warning: Key 'heating' not found in results. Skipping extension for this key.")

                # Extend object arrays (spectrum, reaction) which store arbitrary Python objects.
                # These need to be re-created with the new size and old data copied over.
                # 'spectrum': np.empty((sum_steps, 3), dtype=object)
                if "spectrum" in results and results["spectrum"].shape[1] == 3:
                    current_spec_rows = results["spectrum"].shape[0]
                    if sum_steps > current_spec_rows:
                        temp_spectrum = np.empty((sum_steps, 3), dtype=object)
                        temp_spectrum[:current_spec_rows, :] = results["spectrum"]
                        results["spectrum"] = temp_spectrum
                else:
                    print("Warning: Key 'spectrum' not found or unexpected shape in results. Skipping extension for this key.")

                # 'reaction': np.empty((sum_steps), dtype=object)
                if "reaction" in results:
                    current_react_len = results["reaction"].shape[0]
                    if sum_steps > current_react_len:
                        temp_reaction = np.empty((sum_steps,), dtype=object)
                        temp_reaction[:current_react_len] = results["reaction"]
                        results["reaction"] = temp_reaction
                
        else: # New Results
            results = {
                "step_numbers": step_numbers,
                "keff": [],
                "keff_unc": [],
                "times": [],
                "m_u235": np.array([np.zeros(sum_steps) for i in range(region)]),
                "m_u238": np.array([np.zeros(sum_steps) for i in range(region)]),
                "m_pu239": np.array([np.zeros(sum_steps) for i in range(region)]),
                "m_pu240": np.array([np.zeros(sum_steps) for i in range(region)]),
                "m_pu241": np.array([np.zeros(sum_steps) for i in range(region)]),
                "m_np237": np.array([np.zeros(sum_steps) for i in range(region)]),
                "m_am241": np.array([np.zeros(sum_steps) for i in range(region)]),
                "m_cm244": np.array([np.zeros(sum_steps) for i in range(region)]),
                "a_xe135": np.array([np.zeros(sum_steps) for i in range(region)]),
                "a_sm149": np.array([np.zeros(sum_steps) for i in range(region)]),
                "fission_rate": np.array([np.zeros(sum_steps) for i in range(region)]),
                "n2_rate": np.array([np.zeros(sum_steps) for i in range(region)]),
                "n3_rate": np.array([np.zeros(sum_steps) for i in range(region)]),
                "n4_rate": np.array([np.zeros(sum_steps) for i in range(region)]),
                "alpha_rate": np.array([np.zeros(sum_steps) for i in range(region)]),
                "capture_rate": np.array([np.zeros(sum_steps) for i in range(region)]),
                "proton_rate": np.array([np.zeros(sum_steps) for i in range(region)]),
                "spectrum": np.empty((sum_steps, 3), dtype=object),
                "reaction": np.empty((sum_steps), dtype=object),
                "heating": np.zeros(sum_steps),
            }
            steps_to_process_for_data = step_numbers

        for step in steps_to_process_for_data:
            tmp_results = openmc.deplete.Results(
                f"{cwd}/step_{step}/depletion_result_{step}.h5"
            )
            fuel_id = list(tmp_results[0].index_mat.keys())

            model = openmc.Model.from_xml(
                f"{cwd}/step_{step}/geometry.xml",
                f"{cwd}/step_{step}/materials.xml",
                f"{cwd}/step_{step}/settings.xml",
                f"{cwd}/step_{step}/tallies.xml",
                f"{cwd}/step_{step}/plot.xml",
            )
            geometry = model.geometry

            print(f"##=================== Gen Step: {step} ===================##")

            l_mat = geometry.get_lattices_by_name("FCC-Fuel-Pebbles")
            c_reg = geometry.get_cells_by_name("Pebble-Bed-Cell")
            c_reg = [
                f'{c.name.split("Pebble-Bed-Cell-")[1]}-{c.fill.id}'
                for c in c_reg
                if "FCC-Fuel-Pebbles" in c.fill.name
            ]
            coords_mat = {}
            for c in c_reg:
                coords, l_id = c.split("-")
                coords = int(coords)
                l_id = int(l_id)
                for l in l_mat:
                    if l.id == l_id:
                        u_tmp = openmc.Universe(
                            cells=list(l.universes[0][0][0].cells.values())
                        )
                        g_tmp = openmc.Geometry(u_tmp)
                        c_tmp = g_tmp.get_cells_by_name("Kernel")[0]
                        m_tmp = c_tmp.fill
                        coords_mat[str(m_tmp.id)] = coords - 1

            # Get Times
            times = tmp_results.get_times("s")
            time_intervals = np.diff(times, prepend=0)
            for i in range(len(times)):
                try:
                    sec = results['times'][-1]
                    if i == 0:
                        sec += 0.001
                    else:
                        sec = round(time_intervals[i] + sec)
                except IndexError:
                    sec = step * sum(time_intervals)
                results["times"] += [sec]

            # Get parameters
            keff_ = tmp_results.get_keff()[1]
            results["keff"] += [keff_[r][0] for r in range(len(tmp_results))]
            results["keff_unc"] += [keff_[r][1] for r in range(len(tmp_results))]
            for r in range(len(tmp_results)):
                sp = openmc.StatePoint(f"{cwd}/step_{step}/openmc_simulation_n{r}.h5")
                try:
                    spectrum = sp.get_tally(name="Spectrum")
                    reaction = sp.get_tally(name="Reaction")
                    heating = sp.get_tally(name="Heating-Local")
                    sp.close()
                    spectrum_df = spectrum.get_pandas_dataframe()
                    spectrum_df_x = (
                        spectrum_df.groupby(
                            [("mesh 1", "x"), "energy low [eV]", "energy high [eV]"]
                        )["mean"]
                        .sum()
                        .reset_index()
                        .rename(columns={("mesh 1", "x"): "x"})
                    )
                    spectrum_df_y = (
                        spectrum_df.groupby(
                            [("mesh 1", "y"), "energy low [eV]", "energy high [eV]"]
                        )["mean"]
                        .sum()
                        .reset_index()
                        .rename(columns={("mesh 1", "y"): "y"})
                    )
                    spectrum_df_z = (
                        spectrum_df.groupby(
                            [("mesh 1", "z"), "energy low [eV]", "energy high [eV]"]
                        )["mean"]
                        .sum()
                        .reset_index()
                        .rename(columns={("mesh 1", "z"): "z"})
                    )
                    spectrum_df_x.sort_values(by="energy low [eV]")
                    spectrum_df_y.sort_values(by="energy low [eV]")
                    spectrum_df_z.sort_values(by="energy low [eV]")

                    results["spectrum"][sum_step, 0] = spectrum_df_x
                    results["spectrum"][sum_step, 1] = spectrum_df_y
                    results["spectrum"][sum_step, 2] = spectrum_df_z

                    reaction_df = reaction.get_pandas_dataframe()
                    reaction_df_z = (
                        reaction_df.groupby([("mesh 1", "z"), "score", "nuclide"])["mean"]
                        .sum()
                        .reset_index()
                        .rename(columns={("mesh 1", "z"): "z"})
                    )
                    results["reaction"][sum_step] = reaction_df_z

                    heating_df = heating.get_pandas_dataframe()
                    results["heating"][sum_step] = heating_df["mean"].values[0]

                except Exception:
                    pass

                for fid in fuel_id:
                    results["m_u235"][coords_mat[fid]][sum_step] = tmp_results.get_mass(
                        fid, "U235", mass_units="g", time_units="d"
                    )[1][r]
                    results["m_u238"][coords_mat[fid]][sum_step] = tmp_results.get_mass(
                        fid, "U238", mass_units="g", time_units="d"
                    )[1][r]
                    results["m_pu239"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_mass(fid, "Pu239", mass_units="g", time_units="d")[1][r]
                    )
                    results["m_pu240"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_mass(fid, "Pu240", mass_units="g", time_units="d")[1][r]
                    )
                    results["m_pu241"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_mass(fid, "Pu241", mass_units="g", time_units="d")[1][r]
                    )
                    results["m_np237"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_mass(fid, "Np237", mass_units="g", time_units="d")[1][r]
                    )
                    results["m_am241"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_mass(fid, "Am241", mass_units="g", time_units="d")[1][r]
                    )
                    results["m_cm244"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_mass(fid, "Cm244", mass_units="g", time_units="d")[1][r]
                    )

                    # Get fission products using get_atoms
                    results["a_xe135"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_atoms(fid, "Xe135", time_units="d")[1][r]
                    )
                    results["a_sm149"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_atoms(fid, "Sm149", time_units="d")[1][r]
                    )

                    # Calculate burnup
                    results["fission_rate"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_reaction_rate(fid, "U235", "fission")[1][r]
                    )
                    results["n2_rate"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_reaction_rate(fid, "U235", "(n,2n)")[1][r]
                    )
                    results["n3_rate"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_reaction_rate(fid, "U235", "(n,3n)")[1][r]
                    )
                    results["n4_rate"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_reaction_rate(fid, "U235", "(n,4n)")[1][r]
                    )
                    results["alpha_rate"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_reaction_rate(fid, "U235", "(n,a)")[1][r]
                    )
                    results["capture_rate"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_reaction_rate(fid, "U235", "(n,gamma)")[1][r]
                    )
                    results["proton_rate"][coords_mat[fid]][sum_step] = (
                        tmp_results.get_reaction_rate(fid, "U235", "(n,p)")[1][r]
                    )
                sum_step += 1

    else:
        results = {}
        sp = openmc.StatePoint(f"{cwd}/statepoint.200.h5")
        spectrum = sp.get_tally(name="Spectrum")
        reaction = sp.get_tally(scores=["scatter"])
        heating = sp.get_tally(name="Heating-Local")
        results["keff"] = sp.keff
        results["spectrum"] = spectrum.get_pandas_dataframe()
        results["reaction"] = reaction.get_pandas_dataframe()
        results["heating"] = heating.get_pandas_dataframe()["mean"].values[0]
        sp.close()

    # Save processed result to pickle
    with open(f"{cwd}/results.pkl", "wb") as f:
        pickle.dump(results, f)
        print("Data Saved !!")

    if os.path.exists(f"{cwd}/results_last.pkl"):
        # Load the data from pickle
        with open(f"{cwd}/results_last.pkl", "rb") as f:
            results_last = pickle.load(f)
            print("Last Step Data Loaded !!")
    elif mode == "deplete":
        results_last = {}
        try:
            tmp_results = openmc.deplete.Results(
                f"{cwd}/step_{step_numbers[-1]}/depletion_result_{step_numbers[-1]}.h5"
            )
            sp = openmc.StatePoint(
                f"{cwd}/step_{step_numbers[-1]}/openmc_simulation_n{len(tmp_results)-1}.h5"
            )
            spectrum = sp.get_tally(name="Spectrum")
            reaction = sp.get_tally(scores=["scatter"])
            heating = sp.get_tally(name="Heating-Local")
            results_last["keff"] = sp.keff
            results_last["spectrum"] = spectrum.get_pandas_dataframe()
            results_last["reaction"] = reaction.get_pandas_dataframe()
            results_last["heating"] = heating.get_pandas_dataframe()["mean"].values[0]
            sp.close()
            with open(f"{cwd}/results_last.pkl", "wb") as f:
                pickle.dump(results_last, f)
                print("Last Step Data Saved !!")
        except Exception:
            print("Last Step Data Ignored..")
    print(
        f"Generation Completed in {datetime.timedelta(seconds=(time.perf_counter() - gen_time))}"
    )

get_results()

if mode == "deplete":
    # Get time steps in days
    time_days = np.array(results["times"]) / (3600 * 24)
    sum_steps = len(results["heating"])
    step_numbers = results["step_numbers"]
    # time_days = [(timestep/numberstep) * i for i in range(len(step_numbers)+1)]
    # for i in range(len(step_numbers)+1):
    #     t = i * (timestep/numberstep)
    #     if int(t - (timestep/numberstep)) % timestep == 0 and int(t) not in  [0, timestep]:
    #         time_days += [(t - (timestep/numberstep)) + 0.001]
    #     time_days += [t]

    # Calculate time INTERVALS between depletion steps
    time_intervals = np.diff(time_days, prepend=0)
    m_u235 = np.zeros(sum_steps)
    m_u238 = np.zeros(sum_steps)
    m_pu239 = np.zeros(sum_steps)
    m_pu240 = np.zeros(sum_steps)
    m_pu241 = np.zeros(sum_steps)
    m_np237 = np.zeros(sum_steps)
    m_am241 = np.zeros(sum_steps)
    m_cm244 = np.zeros(sum_steps)
    
    a_xe135 = np.zeros(sum_steps)
    a_sm149 = np.zeros(sum_steps)

    for i in range(region):
        for s in range(sum_steps):
            m_u235[s] += results["m_u235"][i][s]
            m_u238[s] += results["m_u238"][i][s]
            m_pu239[s] += results["m_pu239"][i][s]
            m_pu240[s] += results["m_pu240"][i][s]
            m_pu241[s] += results["m_pu241"][i][s]
            m_np237[s] += results["m_np237"][i][s]
            m_am241[s] += results["m_am241"][i][s]
            m_cm244[s] += results["m_cm244"][i][s]

            a_xe135[s] += results["a_xe135"][i][s]
            a_sm149[s] += results["a_sm149"][i][s]



# ========================================================
# ===================== FUNCTIONS
# =========================================================
def flatten_dataReg(data, region, sum_steps):
    d = np.zeros(sum_steps)
    for i in range(region):
        for j in range(sum_steps):
            d[j] += data[i][j]
    return d

def plot_spectrum(spectrum_df, heating, step=None):
    overall_df = (
        spectrum_df.groupby(["energy low [eV]", "energy high [eV]"])["mean"]
        .sum()
        .reset_index()
    )
    Vol = np.pi * 250.459**2 * 1670
    counts = overall_df["mean"].values * 250e6 / (1.602e-19 * heating * Vol)
    lows = overall_df["energy low [eV]"].values
    highs = overall_df["energy high [eV]"].values

    edges = np.append(lows, highs[-1])
    init_landscape_fig()
    plt.stairs(
        counts, edges, fill=True, color="skyblue", edgecolor="navy", linewidth=1.2
    )
    # Formatting
    plt.xscale("log")
    plt.yscale("log")
    plt.title("Neutron Energy Spectrum")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Counts (particle/$cm^2$⋅s)")
    plt.xlim(edges[0], edges[-1])
    plt.ylim(1e+7, 1e+14)

    # Add energy range labels
    for low, high, density in zip(lows, highs, counts):
        plt.text(
            x=np.sqrt(low * high),  # Geometric mean for label position
            y=density,  # Offset above the bar
            s=f"{density:.2e}",
            ha="center",
            va="bottom",
            fontsize=16,
        )
    plt.tight_layout()
    if mode == "deplete":
        os.makedirs(f"{cwd}/Spectrum", exist_ok=True)
        plt.savefig(f"{cwd}/Spectrum/Spectrum_{step:03d}.png")
    else:
        plt.savefig(f"{cwd}/Spectrum.png")
    plt.close('all')

    mask_thermal = (spectrum_df["energy high [eV]"] <= 1)
    mask_resonance = (spectrum_df["energy low [eV]"] >= 1) & (
        spectrum_df["energy high [eV]"] <= 1e3
    )
    mask_fast = spectrum_df["energy low [eV]"] >= 1e3
    masks = [
        [mask_thermal, "Thermal"],
        [mask_resonance, "Resonan"],
        [mask_fast, "Fast"],
    ]
    for mask in masks:
        df_epi = spectrum_df[mask[0]]
        z_df = df_epi.groupby("z")["mean"].sum().reset_index()
        x = z_df["mean"] * 250e6 / (1.602e-19 * heating * Vol / len(z_df["z"]))
        y = z_df["z"]
        y = y * 1670 / y.max()
        plot_scoreVaxis(x, y, mask[1], step)


def plot_scoreVaxis(X, Y, tipe: str, step=None):
    init_landscape_fig()
    plt.title(f"Distribusi {tipe} {f'Step {step}' if mode == 'deplete' else ''}")
    plt.xlabel(f"Total {tipe}")
    plt.ylabel("Tinggi Reaktor (cm)")
    plt.ylim(0, 1670)
    plt.xscale("log")
    plt.barh(
        Y,
        X,
        height=7,
        # marker='o',           # Marker titik
        # linestyle='-',        # Garis penghubung
        color="green",
        edgecolor="black",
        linewidth=0.5,
        # markersize=8
    )
    for x, y in zip(X, Y):
        if x != 0:
            plt.text(
                x=x,  # Geometric mean for label position
                y=y,  # Offset above the bar
                s=f"{x:.3e}",
                ha="left",
                va="center",
                fontsize=6,
            )
    ax = plt.gca()
    plt.tight_layout()
    tipe = tipe.replace(":", "_")
    if mode == "deplete":
        os.makedirs(f"{cwd}/{tipe}_Z", exist_ok=True)
        plt.savefig(f"{cwd}/{tipe}_Z/{tipe}_Z_{step:03d}.png")
    else:
        plt.savefig(f"{cwd}/{tipe}_Z.png")
    plt.close('all')


def plot2D(x, y, A, tipe: str, plane: str, figsize):
    cmap = plt.get_cmap("plasma")
    eps = 1e-16

    X, Y = np.meshgrid(x, y, indexing="ij")
    fig = plt.figure(figsize=figsize)
    if isinstance(A, list):
        norm = mpl.colors.LogNorm(vmax=max([b.max() for b in A]))
        fig.suptitle(f"Distribusi {tipe} 2D")
        planes = ["XZ", "YZ"]
        axes = []
        for i, a in enumerate(A):
            ax = fig.add_subplot(1, len(A), i + 1)
            ax.set_xlabel(r"Lebar Reaktor (cm)")
            ax.set_ylabel(r"Tinggi Reaktor (cm)")
            ax.set_xlim(x.min(), x.max())
            ax.set_ylim(y.min(), y.max())
            ax.set_title(f"Bidang {planes[i]}")
            ax.set_aspect("equal")
            c = ax.pcolor(X, Y, a, cmap=cmap, norm=norm, zorder=3)
            axes.append(ax)
    else:
        norm = mpl.colors.LogNorm(vmax=A.max())
        plt.title(f"Distribusi {tipe} 2D")
        plt.xlabel(r"Lebar Reaktor (cm)")
        plt.ylabel(r"Tinggi Reaktor (cm)")
        plt.xlim(x.min(), x.max())
        plt.ylim(y.min(), y.max())
        plt.axes().set_aspect("equal")
        c = plt.pcolor(X, Y, A, cmap=cmap, norm=norm, zorder=3)
    cb = fig.colorbar(c)
    cb.set_label(f"${tipe}/cm^2⋅s$")
    plt.tight_layout()
    tipe = tipe.replace(":", "_")
    os.makedirs(f"{cwd}/Dist/", exist_ok=True)
    plt.savefig(f"{cwd}/Dist/{tipe}_{plane}.png")
    plt.close('all')


def plot3D(x, y, A, tipe, plane, figsize):
    cmap = plt.get_cmap("plasma")
    eps = 1e-16
    norm = mpl.colors.LogNorm(vmax=A.max())

    X, Y = np.meshgrid(x, y, indexing="ij")
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_title(f"Distribusi {tipe} 2D")
    ax1.set_xlabel(r"Lebar Reaktor (cm)")
    ax1.set_ylabel(r"Lebar Reaktor (cm)")
    ax1.set_xlim(x.min(), x.max())
    ax1.set_ylim(y.min(), y.max())
    ax1.set_aspect("equal")
    c = ax1.pcolor(X, Y, A, cmap=cmap, norm=norm, zorder=3)

    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    p = ax2.plot_surface(X, Y, A, linewidth=0, cmap=cmap, norm=norm, zorder=3)
    ax2.set_title(f"Distribusi {tipe} 3D")
    ax2.set_xlabel(r"Lebar Reaktor (cm)")
    ax2.set_ylabel(r"Lebar Reaktor (cm)")
    # offset_text = ax.zaxis.get_offset_text()
    # offset_text.set_size(20)
    ax2.set_xlim(x.min(), x.max())
    ax2.set_ylim(y.min(), y.max())
    
    cb = fig.colorbar(p)
    cb.set_label(f"${tipe}/cm^2⋅s$")
    plt.tight_layout()
    tipe = tipe.replace(":", "_")
    os.makedirs(f"{cwd}/Dist/", exist_ok=True)
    plt.savefig(f"{cwd}/Dist/{tipe}_{plane}.png")
    plt.close('all')


def dist(tally_df, tipe, heating, px, py, pz):
    tally_df_xy = tally_df[(tally_df["z"] == pz)]
    tally_df_xz = tally_df[(tally_df["y"] == py)]
    tally_df_yz = tally_df[(tally_df["x"] == px)]
    x = tally_df_xy["x"].unique()
    y = tally_df_xy["y"].unique()
    z = tally_df_yz["z"].unique()
    Vol = np.pi * 250.459**2 * 1670
    Vol_V = 2 * 250.459 * 1670 * (2 * 250.459 / 50) 
    Vol_H = np.pi * (250.459 * 2)**2 * (1670 / 160)
    x = x * 250.459 * 2 / x.max()
    y = y * 250.459 * 2 / y.max()
    z = z * 1670 / z.max()

    # Pivot the DataFrame to reshape 1D data into 2D grid
    pivot_df_xy = tally_df_xy.pivot(index="x", columns="y", values="mean")
    pivot_df_xz = tally_df_xz.pivot(index="x", columns="z", values="mean")
    pivot_df_yz = tally_df_yz.pivot(index="y", columns="z", values="mean")

    # Compute A (2D array)
    A_xy = pivot_df_xy.values * 250e6 / (1.602e-19 * heating * Vol_H)
    A_xz = pivot_df_xz.values * 250e6 / (1.602e-19 * heating * Vol_V)
    A_yz = pivot_df_yz.values * 250e6 / (1.602e-19 * heating * Vol_V)
    if tipe == "fission":
    	global fission_dist
    	fission_dist = [A_xy, A_xz, A_yz]
    elif tipe == "absorption":
    	A_xy -= fission_dist[0]
    	A_xz -= fission_dist[1]
    	A_yz -= fission_dist[2]
    plot3D(x, y, A_xy, tipe, "XY", (33, 15))
    plot2D(x, z, [A_xz, A_yz], tipe, "2D", (26, 30))

def exponential_decay(x, a, b, c):
    return a + b * np.exp(-c * x)

def plot_trend(x, y, y_=None, x_max=None, label="", dec=2, num_type="e"):
    p0 = [1.0, 0.3, 0.001]
    maxfev = 100000
    bounds = ([-np.inf, -np.inf, 0], [np.inf, np.inf, np.inf])
    x_max = x_max or max(x)
    y_max = max(y)
    y_scaled = y / y_max
    if y_ is not None:
        sigma_scaled = y_ / y_max
        params, covariance = curve_fit(exponential_decay, x, y_scaled, p0=p0, sigma=sigma_scaled, bounds=bounds, maxfev=maxfev)
    else:
        params, covariance = curve_fit(exponential_decay, x, y_scaled, p0=p0, bounds=bounds, maxfev=maxfev)
    a_fit, b_fit, c_fit = params
    a_fit *= y_max
    b_fit *= y_max
    x_fit = np.linspace(0, x_max, 5000) # More points for a smooth curve
    y_fit = exponential_decay(x_fit, a_fit, b_fit, c_fit)
    sign = "+" if b_fit >= 0 else "-"
    default_label = (
        f'${a_fit:.{dec}{num_type}} '
        f'{sign} {abs(b_fit):.{dec}{num_type}}'
        f'e^{{-{c_fit:.{dec}{num_type}}x}}$'
    )
    if label != "": label += ": "
    label += default_label
    plt.plot(x_fit, y_fit, '--', label=label, zorder=3)
    return True

def init_landscape_fig():
    return plt.figure(figsize=(30,17))

def init_potrait_fig():
    return plt.figure(figsize=(18,25))

plot_time = time.perf_counter()
tally = ["fission", "scatter", "absorption", "capture", "total", "(n,a)"]
reacsotopes = ["C0", "U235", "U238", "O16", "total"]
if mode == "deplete":
    steps = find_segment_indices(results["times"])

    init_landscape_fig()
    x = time_days
    y = results["keff"]
    y_ = results["keff_unc"]
    plt.errorbar(x, y, yerr=y_, fmt="r-o", alpha=0.4)
    plt.axhline(1.0, color="k", linestyle="--")
    plot_trend(x[20:], y[20:], y_[20:], dec=4, num_type="f")
    plt.xlabel("Waktu (hari)")
    plt.ylabel("k-effective ± σ")
    plt.xlim(0, None)
    plt.title("Evolusi Eigenvalue")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("Eigenvalue.png")
    plt.close('all')

    init_landscape_fig()
    x = np.array([time_days[i] for i in steps])
    y = np.array([results["keff"][i] for i in steps])
    y_ = np.array([results["keff_unc"][i] for i in steps])
    plt.errorbar(x, y, yerr=y_, fmt="r-o", alpha=0.4)
    plt.axhline(1.0, color="k", linestyle="--")
    plot_trend(x[20:], y[20:], y_[20:], dec=4, num_type="f")
    plt.xlabel("Waktu (hari)")
    plt.ylabel("k-effective ± σ")
    plt.xlim(0, None)
    plt.title("Eigenvalue Steady-state")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("Eigenvalue-Steady.png")
    plt.close('all')

    init_landscape_fig()
    plt.plot(time_days, m_u235, "-", alpha=0.4)
    plt.plot(time_days, m_u238, "-", alpha=0.4)
    plot_trend(time_days[10:], m_u235[10:], label="$^{235}$U")
    plot_trend(time_days[10:], m_u238[10:], label="$^{238}$U")
    plt.xlabel("Waktu (hari)")
    plt.ylabel("Massa (gram)")
    plt.xlim(0, None)
    plt.yscale("log")
    plt.title("Evolusi Massa Uranium")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("Overall Uranium.png")
    plt.close('all')

    init_landscape_fig()
    plt.plot(time_days, m_pu239, "-", label="$^{239}$Pu")
    plt.plot(time_days, m_pu240, "-", label="$^{240}$Pu")
    plt.plot(time_days, m_pu241, "-", label="$^{241}$Pu")
    plt.plot(time_days, m_np237, "-", label="$^{237}$Np")
    plt.plot(time_days, m_am241, "-", label="$^{241}$Am")
    plt.plot(time_days, m_cm244, "-", label="$^{244}$Cm")
    plt.xlabel("Waktu (hari)")
    plt.ylabel("Massa (gram)")
    plt.xlim(0, None)
    plt.yscale("log")
    plt.title("Evolusi Massa Aktinida")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("Overall Actinide.png")
    plt.close('all')

    init_landscape_fig()
    plt.plot(time_days, a_xe135, "-", alpha=0.4)
    plt.plot(time_days, a_sm149, "-", alpha=0.4)
    plot_trend(time_days[20:], a_xe135[20:], label="$^{135}$Xe")
    plot_trend(time_days[20:], a_sm149[20:], label="$^{149}$Sm")
    plt.xlabel("Waktu (hari)")
    plt.ylabel("Atom")
    plt.xlim(0, None)
    plt.yscale("log")
    plt.title("Evolusi Racun Neutron")
    plt.legend()
    plt.grid(True, which="both", axis="both", linestyle="-", alpha=0.7)
    plt.tight_layout()
    plt.savefig("Overall Poison.png")
    plt.close('all')

    for i, step in enumerate(step_numbers):
        if not step >= 282:
            continue
        print(f"##=================== Plot Step: {step} ===================##")
        step_time = time.perf_counter()

        try:
            spectrum_df = results["spectrum"][steps[i]][2]
            heating = results["heating"][steps[i]]
            plot_spectrum(spectrum_df, heating, step)
            reaction_df_z = results["reaction"][steps[i]]
            grouped = reaction_df_z.groupby(["score", "nuclide"])
            for tipe in tally:
                for isotope in reaction_df_z["nuclide"].unique():
                    grup = tipe if tipe != "capture" else "(n,gamma)"
                    try:
                        tally_df_z = grouped.get_group((grup, isotope))
                    except Exception:
                        print(f"Step {step} doesnt have {(grup, isotope)} tally")
                        continue
                    Vol = np.pi * 250.459**2 * 1670 
                    x = tally_df_z["mean"] * 250e6 / (1.602e-19 * heating * Vol / len(tally_df_z["z"]))
                    y = tally_df_z["z"]
                    y = y * 1670 / y.max()
                    plot_scoreVaxis(x, y, f"{isotope}:{tipe}", step)


        except AttributeError as e:
            pass

        products = ["m_u235", "m_u238", "m_pu239", "m_pu240",
                    "m_pu241", "m_np237", "m_am241", "m_cm244",
                    "a_xe135", "a_sm149"]
        for product in products:
            params = results[product].T.copy()
            x = params[steps[i]]
            y = np.array(range(len(params[steps[i]])))
            y = y * (605 + 562.33) / y.max() + 326.17
            plot_scoreVaxis(x, y, product, step)
        print(
            f"Step plotting takes {datetime.timedelta(seconds=(time.perf_counter() - step_time))}"
        )
    print("Step Plotting COMPLETE!")

    reaction_df_last = results_last["reaction"]
    heating = results_last["heating"]
    grouped = reaction_df_last.groupby(["score", "nuclide"])
    cmap = plt.get_cmap("plasma")
    for tipe in tally:
        grup = tipe if tipe != "capture" else "(n,gamma)"
        for isotope in reaction_df_last["nuclide"].unique():
            try:
                tally_df_last = grouped.get_group((grup, isotope))
            except Exception:
                continue
            tally_df_last.columns = [
                col[1].strip() if col[1] else col[0] for col in tally_df_last.columns
            ]
            try:
                dist(tally_df_last, f"{isotope}:{tipe}", heating, 25, 25, 103)
            except Exception:
                pass

    spectrum_df = results_last["spectrum"]
    spectrum_df.columns = [
        col[1].strip() if col[1] else col[0] for col in spectrum_df.columns
    ]
    mask_thermal = (spectrum_df["energy high [eV]"] <= 1)
    mask_resonance = (spectrum_df["energy low [eV]"] >= 1) & (
        spectrum_df["energy high [eV]"] <= 1e3
    )
    mask_fast = spectrum_df["energy low [eV]"] >= 1e3
    masks = [
        [mask_thermal, "Thermal"],
        [mask_resonance, "Resonan"],
        [mask_fast, "Fast"],
    ]
    for mask in masks:
        df_epi = spectrum_df[mask[0]]
        epi_df = df_epi.groupby(["x", "y", "z"])["mean"].sum().reset_index()
        dist(epi_df, mask[1], heating, 25, 25, 103)

    flux_df = spectrum_df.groupby(["x", "y", "z"])["mean"].sum().reset_index()
    dist(flux_df, "Fluks", heating, 25, 25, 103)


elif mode == "combine":
    with open(f"{cwd}/results_uo2.pkl", "rb") as f:
        results_uo2 = pickle.load(f)
    with open(f"{cwd}/results_uco_low.pkl", "rb") as f:
        results_uco_low = pickle.load(f)
    with open(f"{cwd}/results_uco_high.pkl", "rb") as f:
        results_uco_high = pickle.load(f)
    results_group = [results_uo2, results_uco_low, results_uco_high]
    print("Data Loaded !!")
    os.makedirs(f"{cwd}/spreadsheet", exist_ok=True)
    colors = ["#a0dd01", "#ffa7ae", "#ff747f"]
    fuels = ["UO$^2$-8,6%", "UCO-8,6%", "UCO-15,5%"]

    time_uo2 = results_uo2["times"]
    step_numbers_uo2 = np.array(results_uo2["step_numbers"])
    indexes_uo2 = np.where(np.isin(step_numbers_uo2, [655, 659, 663]))[0]
    steps_uo2 = find_segment_indices(time_uo2)
    step_uo2 = np.array([steps_uo2[i] for i in indexes_uo2])
    sum_steps_uo2 = len(results_uo2["heating"])

    time_uco_low = results_uco_low["times"]
    step_numbers_uco_low = np.array(results_uco_low["step_numbers"])
    indexes_uco_low = np.where(np.isin(step_numbers_uco_low, [655, 659, 663]))[0]
    steps_uco_low = find_segment_indices(time_uco_low)
    step_uco_low = np.array([steps_uco_low[i] for i in indexes_uco_low])
    sum_steps_uco_low = len(results_uco_low["heating"])

    time_uco_high = results_uco_high["times"]
    step_numbers_uco_high = np.array(results_uco_high["step_numbers"])
    indexes_uco_high = np.where(np.isin(step_numbers_uco_high, [282, 286, 290]))[0]
    steps_uco_high = find_segment_indices(time_uco_high)
    step_uco_high = np.array([steps_uco_high[i] for i in indexes_uco_high])
    sum_steps_uco_high = len(results_uco_high["heating"])

    steps_group = [step_uo2, step_uco_low, step_uco_high]

    init_landscape_fig()
    x_uo2 = np.linspace(0, 1, len(results_uo2["keff"]))
    x_uco_low = np.linspace(0, 1, len(results_uco_low["keff"]))
    x_uco_high = np.linspace(0, 1, len(results_uco_high["keff"]))
    time_days_uo2 = np.array(time_uo2) / (3600 * 24)
    time_days_uco_low = np.array(time_uco_low) / (3600 * 24)
    time_days_uco_high = np.array(time_uco_high) / (3600 * 24)
    days_max = np.max(np.concatenate((time_days_uo2,time_days_uco_low,time_days_uco_high)))
    plot_trend(time_days_uo2[50:], results_uo2["keff"][50:], results_uo2["keff_unc"][50:], x_max=8000, label="${UO_2}$-8,6%", dec=4, num_type="f")
    plot_trend(time_days_uco_low[50:], results_uco_low["keff"][50:], results_uco_low["keff_unc"][50:], x_max=8000, label="UCO-8,6%", dec=4, num_type="f")
    plot_trend(time_days_uco_high[50:], results_uco_high["keff"][50:], results_uco_high["keff_unc"][50:], x_max=8000, label="UCO-15,5%", dec=4, num_type="f")
    plt.plot(time_days_uo2, results_uo2["keff"], color=colors[0], alpha=0.4)
    plt.plot(time_days_uco_low, results_uco_low["keff"], color=colors[1], alpha=0.4)
    plt.plot(time_days_uco_high, results_uco_high["keff"], color=colors[2], alpha=0.4)
    plt.axhline(1.0, color="k", linestyle="--")
    plt.xlabel("Waktu (hari)")
    plt.ylabel("k-effective")
    plt.xlim(0, None)
    plt.title("Evolusi Eigenvalue")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("Eigenvalue.png")
    plt.close('all')

    init_landscape_fig()
    y_uo2 = [results_uo2["keff"][i] for i in step_uo2]
    y_uo2_ = [results_uo2["keff_unc"][i] for i in step_uo2]

    y_uco_low = [results_uco_low["keff"][i] for i in step_uco_low]
    y_uco_low_ = [results_uco_low["keff_unc"][i] for i in step_uco_low]

    y_uco_high = [results_uco_high["keff"][i] for i in step_uco_high]
    y_uco_high_ = [results_uco_high["keff_unc"][i] for i in step_uco_high]

    x_uo2 = np.linspace(0, 1, len(y_uo2))
    x_uco_low = np.linspace(0, 1, len(y_uco_low))
    x_uco_high = np.linspace(0, 1, len(y_uco_high))
    plt.errorbar(x_uo2, y_uo2, yerr=y_uo2_, fmt="r-o", label="UO2-8,6%")
    plt.errorbar(x_uco_low, y_uco_low, yerr=y_uco_low_, fmt="b-o", label="UCO-8,6%")
    plt.errorbar(x_uco_high, y_uco_high, yerr=y_uco_high_, fmt="g-o", label="UCO-15,5%")
    plt.axhline(1.0, color="k", linestyle="--")
    plt.xlabel("Index")
    plt.ylabel("k-effective ± σ")
    plt.xlim(0, None)
    plt.title("Eigenvalue Steady-state")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("Eigenvalue-Steady.png")
    plt.close('all')

    init_landscape_fig()
    m_u235_uo2 = flatten_dataReg(results_uo2["m_u235"], region, len(results_uo2["heating"]))
    m_u235_uco_low = flatten_dataReg(results_uco_low["m_u235"], region, len(results_uco_low["heating"]))
    plt.plot(time_days_uo2, m_u235_uo2, color=colors[0], alpha=0.4)
    plt.plot(time_days_uco_low, m_u235_uco_low, color=colors[1], alpha=0.4)
    plot_trend(time_days_uo2, m_u235_uo2, label="${UO_2}$-8,6%")
    plot_trend(time_days_uco_low, m_u235_uco_low, label="UCO-8,6%")
    plt.xlabel("Waktu (hari)")
    plt.ylabel("Massa (gram)")
    plt.xlim(0, None)
    plt.yscale("log")
    plt.title("Evolusi Massa $^{235}$U")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("Overall U235.png")
    plt.close('all')

    init_landscape_fig()
    m_u238_uo2 = flatten_dataReg(results_uo2["m_u238"], region, len(results_uo2["heating"]))
    m_u238_uco_low = flatten_dataReg(results_uco_low["m_u238"], region, len(results_uco_low["heating"]))
    plt.plot(time_days_uo2, m_u238_uo2, color=colors[0], alpha=0.4)
    plt.plot(time_days_uco_low, m_u238_uco_low, color=colors[1], alpha=0.4)
    plot_trend(time_days_uo2, m_u238_uo2, label="${UO_2}$-8,6%")
    plot_trend(time_days_uco_low, m_u238_uco_low, label="UCO-8,6%")
    plt.xlabel("Waktu (hari)")
    plt.ylabel("Massa (gram)")
    plt.xlim(0, None)
    plt.yscale("log")
    plt.title("Evolusi Massa $^{238}$U")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("Overall U238.png")
    plt.close('all')

    keff_dict = {
        "Keff UO2-8,6%": y_uo2,
        "Keff UCO-8,6%": y_uco_low,
        "Keff UCO-15,5%": y_uco_high,
        "std UO2": y_uo2_,
        "std UCO-8,6%": y_uco_low_,
        "std UCO-15,5%": y_uco_high_,
    }
    

    reaksi_dict = {}
    for tipe in tally:
        grup = tipe if tipe != "capture" else "(n,gamma)"
        for isotope in reacsotopes:
            init_landscape_fig()
            for i in range(len(results_group)):
                results = results_group[i]
                steps = steps_group[i]
                fuel = fuels[i] 
                for j, step in enumerate(steps):
                    heating = results["heating"][step]
                    reaction_df_z = results["reaction"][step]
                    grouped = reaction_df_z.groupby(["score", "nuclide"])
                    try:
                        tally_df_z = grouped.get_group((grup, isotope))
                    except Exception:
                        print(f"Step {step} doesnt have {grup} tally")
                        continue
                    Vol = np.pi * 250.459**2 * 1670
                    x = tally_df_z["mean"] * 250e6 / (1.602e-19 * heating * Vol / len(tally_df_z["z"]))
                    y = tally_df_z["z"]
                    y = y * 1670 / y.max()

                    label = f"{fuel}-{j}"

                    plt.plot(
                        x,
                        y,
                        # color='blue',
                        linewidth=2.4,
                        linestyle="--",
                        label=label,
                    )
                    reaksi_dict[f"{isotope}:{tipe} {label}"] = x

            plt.xlabel(f"{isotope}:{tipe}")
            plt.ylabel("Tinggi Reaktor (cm)")
            plt.xscale("log")
            plt.title(f"Distribusi {isotope}:{tipe} sepanjang Tinggi Reaktor")
            plt.legend()
            plt.tight_layout(pad=2.0)
            plt.ylim(0, 1670)
            plt.savefig(f"{cwd}/{isotope}_{tipe}_Z.png")
            plt.close('all')
    reaksi_dict["Tinggi Reaktor"] = y

    absorb = {}
    fisi = {}

    for tipe in tally: 
        for nuclide in reacsotopes:
            init_potrait_fig()
            num = np.zeros(len(results_group))
            y = np.zeros(len(results_group))

            for i, fuel in enumerate(fuels):
                for key, values in reaksi_dict.items():
                    if fuel in key and nuclide in key:
                        if tipe in key:
                            y[i] += np.sum(values)
                            num[i] += 1
                        if "fission" in key and tipe == "fission":
                            fisi[(fuel, nuclide)] = fisi.get((fuel, nuclide), 0.0) + np.sum(values)

                        elif "absorption" in key and tipe == "absorption":
                            absorb[(fuel, nuclide)] = absorb.get((fuel, nuclide), 0.0) + np.sum(values)

            y = np.divide(y, num, out=np.zeros_like(y), where=num > 0)
            plt.bar(fuels, y, color=colors, edgecolor="black")
            min_val = np.min(y[np.nonzero(y)]) if np.any(y) else 0
            for j, val in enumerate(y):
                if val > 0:
                    plt.text(j, val, f"{1 + ((val - min_val) / min_val):.2%}",
                            ha="center", va="bottom", fontsize=36)
                    plt.text(j, val, f"{val:.2e}",
                            ha="center", va="top", fontsize=36)
            plt.xlabel("Bahan Bakar")
            plt.ylabel(f"${tipe}/cm^2⋅s$")
            plt.yscale("log")
            plt.title(f"Perbandingan {nuclide}:{tipe}")
            plt.tight_layout(pad=2.0)
            plt.savefig(f"{cwd}/{nuclide}_{tipe}.png")
            plt.close("all")

    for nuclide in reaction_df_z["nuclide"].unique():
        init_potrait_fig()
        p_absorb = []
        for fuel in fuels:
            a_val = absorb.get((fuel, nuclide), 0.0)
            f_val = fisi.get((fuel, nuclide), 0.0)
            p_absorb.append(a_val - f_val)

        p_absorb = np.array(p_absorb)

        plt.bar(fuels, p_absorb, color=colors, edgecolor="black")

        nonzero = p_absorb > 0
        if np.any(nonzero):
            min_val = np.min(p_absorb[nonzero])
            for i, val in enumerate(p_absorb):
                if val > 0:
                    plt.text(i, val, f"{1 + ((val - min_val)/min_val):.4%}",
                            ha="center", va="bottom", fontsize=36)
                    plt.text(i, val, f"{val:.2e}",
                            ha="center", va="top", fontsize=36)

        plt.xlabel("Bahan Bakar")
        plt.ylabel(r"$absorption/cm^2\cdot s$")
        plt.yscale("log")
        plt.title(f"Perbandingan {nuclide}:pabsorption")
        plt.tight_layout(pad=2.0)
        plt.savefig(f"{cwd}/{nuclide}_pabsorption.png")
        plt.close("all")

    


    spektrum_dict = {}
    masks = ["Thermal", "Resonance", "Fast"]
    for mask in masks:
        init_landscape_fig()
        for i in range(len(results_group)):
            results = results_group[i]
            steps = steps_group[i]
            fuel = fuels[i] 
            for j, step in enumerate(steps):
                heating = results["heating"][step]
                spectrum_df = results["spectrum"][step][2]
                mask_dict = {
                    "Thermal": (spectrum_df["energy high [eV]"] <= 1),
                    "Resonance": (spectrum_df["energy low [eV]"] >= 1) & (
                                  spectrum_df["energy high [eV]"] <= 1e3
                            ),
                    "Fast": spectrum_df["energy low [eV]"] >= 1e3
                }
                mask_spectrum = mask_dict[mask]
                df_epi = spectrum_df[mask_spectrum]
                z_df = df_epi.groupby("z")["mean"].sum().reset_index()
                Vol = np.pi * 250.459**2 * 1670
                x = z_df["mean"] * 250e6 / (1.602e-19 * heating * Vol / len(tally_df_z["z"]))
                y = z_df["z"]
                y = y * 1670 / y.max()

                label = f"{fuel}-{j}"

                plt.plot(
                    x,
                    y,
                    # color='blue',
                    linewidth=2.4,
                    linestyle="-.",
                    label=label,
                )
                spektrum_dict[f"{mask}-{label}"] = x

        plt.xlabel(f"Total Fluks {mask}")
        plt.ylabel("Tinggi Reaktor (cm)")
        plt.xscale("log")
        plt.title(f"Distribusi Fluks {mask} sepanjang Tinggi Reaktor")
        plt.legend()
        plt.tight_layout(pad=2.0)
        plt.ylim(0, 1670)
        plt.savefig(f"{cwd}/{mask}_Z.png")
        plt.close('all')
        spektrum_dict["Tinggi Reaktor"] = y

    for mask in masks:
        init_potrait_fig()
        num = np.zeros(len(results_group))
        y = np.zeros(len(results_group))
        for i in range(len(results_group)):
            fuel = fuels[i]
            for key in spektrum_dict.keys():
                if fuel in key and mask in key:
                    y[i] += np.sum(spektrum_dict[key])
                    num[i] += 1
        y /= num
        plt.bar(fuels, y, color=colors, edgecolor="black")
        min_val = np.min(y)
        for i, val in enumerate(y):
            plt.text(i, val, f"{1 + ((val - min_val)/min_val):.2%}", ha="center", va="bottom", fontsize=36)
            plt.text(i, val, f"{val:.2e}", ha="center", va="top", fontsize=36)
        plt.xlabel(f"Bahan Bakar")
        plt.ylabel(f"${mask}/cm^2⋅s$")
        plt.yscale("log")
        plt.title(f"Perbandingan Fluks {mask} Total")
        plt.tight_layout(pad=2.0)
        plt.savefig(f"{cwd}/{mask}.png")
        plt.close('all')

    isotopes = ["m_u235", "m_u238", "m_pu239", "m_pu240",
                "m_pu241", "m_np237", "m_am241", "m_cm244",
                "a_xe135", "a_sm149"]
    isotopes_dict = {}
    for isotope in isotopes:
        init_landscape_fig()
        for i in range(len(results_group)):
            results = results_group[i]
            steps = steps_group[i]
            isotopee = results[isotope].T.copy()
            fuel = fuels[i]
            for j, step in enumerate(steps):
                x = isotopee[step]
                y = np.array(range(len(isotopee[step])))
                y = y * (605 + 562.33) / y.max() + 326.17

                label = f"{fuel}-{j}"

                plt.plot(
                    x,
                    y,
                    # color='blue',
                    linewidth=2.4,
                    linestyle=":",
                    label=label,
                )
                isotopes_dict[f"{isotope}-{label}"] = x

        plt.xlabel(f"Total {isotope}")
        plt.ylabel("Tinggi Reaktor (cm)")
        plt.xscale("log")
        plt.title(f"Distribusi {isotope} sepanjang Tinggi Reaktor")
        plt.legend()
        plt.tight_layout(pad=2.0)
        plt.ylim(0, 1670)
        plt.savefig(f"{cwd}/{isotope}_Z.png")
        plt.close('all')
    isotopes_dict["Tinggi Reaktor"] = y

    for isotope in isotopes:
        init_potrait_fig()
        num = np.zeros(len(results_group))
        y = np.zeros(len(results_group))
        for i in range(len(results_group)):
            fuel = fuels[i]
            for key in isotopes_dict.keys():
                if fuel in key and isotope in key:
                    y[i] += np.sum(isotopes_dict[key])
                    num[i] += 1
        y /= num
        plt.bar(fuels, y, color=colors, edgecolor="black")
        min_val = np.min(y)
        for i, val in enumerate(y):
            plt.text(i, val, f"{1 + ((val - min_val)/min_val):.2%}", ha="center", va="bottom", fontsize=36)
            plt.text(i, val, f"{val:.2e}", ha="center", va="top", fontsize=36)
        plt.xlabel(f"Bahan Bakar")
        ylabel = "Massa (g)" if isotope.startswith("m") else "Atom"
        plt.ylabel(ylabel)
        plt.yscale("log")
        plt.title(f"Perbandingan {isotope} Total")
        plt.tight_layout(pad=2.0)
        plt.savefig(f"{cwd}/{isotope}.png")
        plt.close('all')

    dict_list = [keff_dict, reaksi_dict, spektrum_dict, isotopes_dict]
    dic_name = ['Eigenvalue', "Reaksi", 'Spektrum', "Isotopes"]
    for i, dic in enumerate(dict_list):
        headers = list(dic.keys())
        normalized = []
        for key in headers:
            value = dic[key]
            if isinstance(value, (list, tuple, np.ndarray)) and not isinstance(value, (str, bytes)):
                lst = list(value)
            elif isinstance(value, pd.core.series.Series):
                lst = [float(item) for item in value.values]
            normalized.append(lst)

        # Pad lists to equal length with None (or customize with empty string)
        max_len = max(len(lst) for lst in normalized)
        padded = []
        for lst in normalized:
            if len(lst) < max_len:
                lst += [None] * (max_len - len(lst))  # Pad with None
            padded.append(lst)
        with open(f"{cwd}/spreadsheet/{dic_name[i]}.csv", "w") as f:
            writer = csv.writer(f)
            writer.writerow(headers)  # Write headers
            for row in zip(*padded):  # Transpose padded data
                writer.writerow(row)


else:
    spectrum_df = results["spectrum"]

    heating = results["heating"]
    spectrum_df_x = (
        spectrum_df.groupby([("mesh 1", "x"), "energy low [eV]", "energy high [eV]"])[
            "mean"
        ]
        .sum()
        .reset_index()
        .rename(columns={("mesh 1", "x"): "x"})
    )
    spectrum_df_y = (
        spectrum_df.groupby([("mesh 1", "y"), "energy low [eV]", "energy high [eV]"])[
            "mean"
        ]
        .sum()
        .reset_index()
        .rename(columns={("mesh 1", "y"): "y"})
    )
    spectrum_df_z = (
        spectrum_df.groupby([("mesh 1", "z"), "energy low [eV]", "energy high [eV]"])[
            "mean"
        ]
        .sum()
        .reset_index()
        .rename(columns={("mesh 1", "z"): "z"})
    )
    spectrum_df_x.sort_values(by="energy low [eV]")
    spectrum_df_y.sort_values(by="energy low [eV]")
    spectrum_df_z.sort_values(by="energy low [eV]")

    plot_spectrum(spectrum_df_z, heating)

    reaction_df = results["reaction"]
    grouped = reaction_df.groupby(["score", "nuclide"])

    cmap = plt.get_cmap("plasma")
    for tipe in tally:
        grup = tipe if tipe != "capture" else "(n,gamma)"
        for isotope in reaction_df["nuclide"].unique():
            try:
                tally_df = grouped.get_group((grup, isotope))
            except Exception:
                        print(f"Step {step} doesnt have {(grup, isotope)} tally")
                        continue
            tally_df.columns = [
                col[1].strip() if col[1] else col[0] for col in tally_df.columns
            ]
            tally_df_z = tally_df.groupby(["z"])["mean"].sum().reset_index()
            Vol = np.pi * 250.459**2 * 1670
            x = tally_df_z["mean"] * 250e6 / (1.602e-19 * heating * Vol / len(tally_df_z["z"]))
            y = tally_df_z["z"]
            y = y * 1670 / y.max()
            plot_scoreVaxis(x, y, f"{isotope}:{tipe}")

            try:                
                dist(tally_df, f"{isotope}:{tipe}", heating, 25, 25, 103)
            except Exception:
                pass

    spectrum_df.columns = [
        col[1].strip() if col[1] else col[0] for col in spectrum_df.columns
    ]
    mask_thermal = (spectrum_df["energy high [eV]"] <= 1)
    mask_resonance = (spectrum_df["energy low [eV]"] >= 1) & (
        spectrum_df["energy high [eV]"] <= 1e3
    )
    mask_fast = spectrum_df["energy low [eV]"] >= 1e3
    masks = [
        [mask_thermal, "Thermal"],
        [mask_resonance, "Resonan"],
        [mask_fast, "Fast"],
    ]
    for mask in masks:
        df_epi = spectrum_df[mask[0]]
        epi_df = df_epi.groupby(["x", "y", "z"])["mean"].sum().reset_index()
        dist(epi_df, mask[1], heating, 25, 25, 103)

    flux_df = spectrum_df.groupby(["x", "y", "z"])["mean"].sum().reset_index()
    dist(flux_df, "Fluks", heating, 25, 25, 103)


print(
    f"Plotting Completed in {datetime.timedelta(seconds=(time.perf_counter() - plot_time))}"
)
