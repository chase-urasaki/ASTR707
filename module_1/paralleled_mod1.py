#%%
import numpy as np
import subprocess
import re
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor

def run_plot_tess(target_name: str, test_period: float, test_depth: float):
    """Runs the external script and extracts transit fit parameters."""
    cmd = [
        "python3",
        "plot_tess_new.py",
        target_name,
        "--fake",
        f"{test_period},{test_depth}"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        match = re.search(r"Transit fit parameters =\s+([\d\.\s]+)", result.stdout)
        if match:
            params = list(map(float, match.group(1).split()))
            return params
    return None


def process_point(period_idx, depth_idx, orbital_period, transit_depth):
    """Handles one grid point for planet injection and detection."""
    test_period = orbital_period[period_idx]
    test_depth = transit_depth[depth_idx]
    trial_array = np.zeros((2, 25))  # 2 rows, 25 columns

    for i in range(trial_array.shape[1]):
        params = run_plot_tess("GJ 480", test_period, test_depth)
        if params is not None:
            trial_array[0, i] = params[0]  # period
            trial_array[1, i] = params[-1]  # depth

    # Compute statistics
    median_period = np.median(trial_array[0, :])
    median_depth = np.median(trial_array[1, :])

    difference_period = np.abs(median_period - test_period) / test_period
    difference_depth = np.abs(median_depth - test_depth) / test_depth

    # Determine detectability status
    if difference_period <= 0.15 and difference_depth <= 0.15:
        detect_status = 1  # Detectable
    elif difference_period <= 0.20 or difference_depth <= 0.20:
        detect_status = 0  # Marginal
    else:
        detect_status = -1  # Undetectable

    return period_idx, depth_idx, detect_status


def main():
    orbital_period = np.linspace(0.5, 23, 50)
    transit_depth = np.logspace(-3, -1, 100)
    detectability = np.zeros((len(orbital_period), len(transit_depth)))

    tasks = [(period_idx, depth_idx, orbital_period, transit_depth)
             for period_idx in range(len(orbital_period))
             for depth_idx in range(len(transit_depth))]

    with ThreadPoolExecutor() as executor:
        # Create a progress bar for the tasks
        with tqdm(total=len(tasks), desc="Processing Points") as pbar:
            futures = {executor.submit(process_point, *task): task for task in tasks}
            for future in futures:
                period_idx, depth_idx, detect_status = future.result()
                detectability[period_idx, depth_idx] = detect_status
                pbar.update(1)

    return detectability, orbital_period, transit_depth


if __name__ == "__main__":
    detectabilities, orbital_periods, transit_depths = main()

    # Save the results to files
    np.save("detectabilities.npy", detectabilities)
    np.save("orbital_periods.npy", orbital_periods)
    np.save("transit_depths.npy", transit_depths)

    #%%