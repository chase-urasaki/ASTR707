#%%
import numpy as np 
import subprocess
import re
from tqdm import tqdm
# %%

def main():
    orbital_period = np.linspace(0.5, 23, 50)
    transit_depth = np.logspace(-3, -1, 100)
    detectability = np.zeros((len(orbital_period), len(transit_depth)))

    # for period_idx in range(len(orbital_period)):
    #     for depth_idx in range(len(transit_depth)):
    #         # Test vals 
    #         print(orbital_period[period_idx], transit_depth[depth_idx])


    # print(type(orbital_period[1]))

    def run_plot_tess(target_name: str, test_period: float, test_depth: float):
        cmd = [
            "python3",  # Use "python" on Windows
            "plot_tess_new.py",
            target_name,
            "--fake",
            f"{test_period},{test_depth}"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print("Script output:", result.stdout)

            # Search string for output params 
            match = re.search(r"Transit fit parameters =\s+([\d\.\s]+)", result.stdout)
            if match:
                params = list(map(float, match.group(1).split()))
                print(f"Extracted transit fit parameters: {params}")
            else:
                print("Could not extract transit fit parameters.")

            return params
        
        else:
            print("Script error:", result.stderr)
            return None
        
    total_its = len(orbital_period) * len(transit_depth)

    with tqdm(total=total_its, desc = "Points") as pbar:
        for period_idx in range(len(orbital_period)):
            for depth_idx in range(len(transit_depth)):
                # Run planet injection for each point on the meshgrid 50 times
                # initialze an array to store the results
                trial_array = np.zeros((2,25)) # 2 rows, 50 columns
                test_period = orbital_period[period_idx]
                test_depth = transit_depth[depth_idx]

                with tqdm(total = trial_array.shape[1], desc = "Trials") as triail_pbar:
                    for i in range(trial_array.shape[1]):

                        # Properly format test_vals as strings for shell command 
                        params = run_plot_tess("GJ 480", test_period, test_depth)
                        
                        # Stores results in the array 
                        trial_array[0,i] = params[0] # period
                        trial_array[1,i] = params[1] # depth 

                        triail_pbar.update(1)

                # Compute the mean period and depth
                median_period = np.median(trial_array[0,:])
                print(f"median period is {median_period} for trial {i} for orbital period {orbital_period[period_idx]}.")
                median_depth = np.median(trial_array[1,:]) 
                print(f"median depth is {median_depth} for trial {i} for transit depth {transit_depth[depth_idx]}.")
                # Compute the differences
                difference_period = np.abs(median_period - test_period)/test_period
                print(f"period difference is {difference_period}")
                difference_depth = np.abs(median_depth - test_depth)/test_depth
                print(f"depth difference is {difference_depth}")
                #   Compute mean LS period, BLS period, and BLS depth
                #  If less than or equal to 15% difference in period AND, label as detectable
                if difference_period <= 0.15 and difference_depth <= 0.15:
                    detectability[period_idx, depth_idx] = 1 # Detectable
                elif difference_period <= 0.20 or difference_depth <= 0.20:
                    detectability[period_idx, depth_idx] = 0 # Marginal
                else:
                    detectability[period_idx, depth_idx] = -1 # Undetectable

                triail_pbar.update(1)

        return detectability, orbital_period, transit_depth

if __name__ == "__main__":
    detectabilities, orbital_periods, transit_depths = main()
    
    #save the results to a file
    np.save("detectabilities.npy", detectabilities)
    np.save("orbital_periods.npy", orbital_periods)
    np.save("transit_depths.npy", transit_depths)


#%%