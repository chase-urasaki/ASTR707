#%%
import numpy as np
import subprocess
import re
from concurrent.futures import ProcessPoolExecutor


def run_plot_tess(target_name: str, test_period: float, test_depth: float):
    cmd = [
        "python3",
        "plot_tess_new.py",
        target_name,
        "--fake",
        f"{test_period},{test_depth}"
    ]
    
    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        # Search string for output params
        match = re.search(r"Transit fit parameters =\s+([\d\.\s]+)", result.stdout)
        if match:
            params = list(map(float, match.group(1).split()))
            print(params)
            return params
        
        return [np.nan, np.nan]
    else:
        return [np.nan, np.nan]


def process_point(period_idx, depth_idx, orbital_period, transit_depth):
    test_period = orbital_period[period_idx]
    test_depth = transit_depth[depth_idx]

    # Initialize trial array to store results
    trial_array = np.zeros((2, 7))
    for i in range(trial_array.shape[1]):
        params = run_plot_tess("GJ 480", test_period, test_depth)
        trial_array[0, i] = params[0]  # period
        trial_array[1, i] = params[-1]  # depth

    # Compute the median period and depth
    median_period = np.median(trial_array[0, :])
    median_depth = np.median(trial_array[1, :])

    # Compute the differences
    difference_period = np.abs(median_period - test_period) / test_period
    difference_depth = np.abs(median_depth - test_depth) / test_depth

    # Determine detectability
    if difference_period <= 0.15 and difference_depth <= 0.15:
        return period_idx, depth_idx, 1  # Detectable
    elif difference_period <= 0.20 or difference_depth <= 0.20:
        return period_idx, depth_idx, 0  # Marginal
    else:
        return period_idx, depth_idx, -1  # Undetectable


def main():
    orbital_period = np.linspace(0.5, 23, 20)[17:] #right endpoint inclusive, left endpoint exclusive
    transit_depth = np.logspace(-3, -1, 20)
    detectability = np.zeros((len(orbital_period), len(transit_depth)))

    MAX_WORKERS = 12

    # Create a list of tasks
    tasks = [
        (period_idx, depth_idx, orbital_period, transit_depth)
        for period_idx in range(len(orbital_period))
        for depth_idx in range(len(transit_depth))
    ]

    # Use ProcessPoolExecutor for parallel execution
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        results = executor.map(process_point, *zip(*tasks))  # Unpack tasks

    # Update the detectability grid based on results
    for period_idx, depth_idx, detect_status in results:
        detectability[period_idx, depth_idx] = detect_status

    print("Detectability grid completed.")
    return detectability, orbital_period, transit_depth


if __name__ == "__main__":
    detectability, orbital_periods, transit_depths = main()

    # Save the results
    np.savetxt("detectability_17_end.csv", detectability, delimiter=',')


#%%