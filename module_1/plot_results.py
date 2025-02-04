#%% 
import numpy as np
orbital_period = np.linspace(0.5, 23, 20)
transit_depth = np.logspace(-4, -3, 5)
detectability = np.zeros((len(orbital_period), len(transit_depth)))

#%%
np.logspace()
#%%
np.savetxt('orbital_periods.csv', orbital_period, delimiter=',')
# %%
np.savetxt('transit_depth.csv', transit_depth, delimiter=',')
# %%
detect_0_2 = np.load('detectability_0_2.npy')
# %%
detect_0_2
# %%
np.shape(detect_0_2)
# %%
np.savetxt('detect_0_2.csv', detect_0_2, delimiter=',')
# %%
detect_2_5 = np.load('detectability_2_5.npy')
#%%
np.savetxt('detectability_2_5.csv',detect_2_5 ,delimiter=',')
# %%
# Code to plot results
import pandas as pd 
from matplotlib import pyplot as plt 

detect_grid = pd.read_csv("Module_1_detectability_low_res.csv", index_col=0)
periods = detect_grid.index.astype(float)
depths = detect_grid.columns.astype(float)
detectability = detect_grid.values
#%%
print(depths)
print(periods)
#
#rotate the array twice since values are going backwards 
rotated_matrix = np.rot90(np.rot90(detectability))

#%%
# Conver the transit depths into planetary radii using R_star = 0.45 * R_solar
R_star = 0.46

def depths_to_R_planet(depths): 
    # Compute R_planet in solar radii
    R_planet_solar = R_star * np.sqrt(depths)

    # Convert to Earth radii
    R_planet_earth = R_planet_solar*109.
    return R_planet_earth
#%%
R_planet_earth = depths_to_R_planet(depths)
# %%
plt.figure(figsize=(10,6))
contour = plt.contourf(periods, depths, rotated_matrix, cmap='jet')
plt.colorbar(contour, label="Detectability (-1: Undetectable, 0: Marginal, 1: Detectable)")
plt.xlabel("Period [days]")
plt.xscale("log")
plt.yscale("log")
plt.ylabel("Injection Depth")
plt.title("Detectability Contour Plot")

# Show the plot
plt.show()
# %%
plt.figure(figsize=(10,6))
contour = plt.contourf(periods, R_planet_earth, rotated_matrix, cmap='jet')
plt.scatter(9.567, 1.75, marker="o", color='k')
plt.colorbar(contour, label="Detectability (-1: Undetectable, 0: Marginal, 1: Detectable)")
plt.xlabel("Period [days]")
plt.xscale("log")
plt.yscale("log")
plt.ylabel("Earth Radii")
plt.title("Detectability Contour Plot")

# Show the plot
plt.show()
# %%
fig, ax1 = plt.subplots(figsize = (8,6))
color = 'tab:blue'
ax1.set_xlabel("Orbital Period (days)")
ax1.set_ylabel("Transit Depth", color=color)
contour = ax1.contourf(orbital_period, transit_depth, rotated_matrix, cmap="jet")
plt.colorbar(contour, ax=ax1, label="Detectability (-1: Undetectable, 0: Marginal, 1: Detectable)")
ax1.set_yscale("log")
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel("Earth Radii", color=color)
ax2.set_yscale("log")
ax2.set_yticks(R_planet_earth)
ax2.tick_params(axis='y', labelcolor=color)
# %%
fig, ax1 = plt.subplots(figsize=(8, 6))
color = 'tab:blue'
ax1.set_xlabel("Orbital Period (days)")
ax1.set_ylabel("Transit Depth", color=color)
contour = ax1.contourf(orbital_period, transit_depth, rotated_matrix, cmap="jet")
cbar = plt.colorbar(contour, ax=ax1)
cbar.set_label("Detectability (-1: Undetectable, 0: Marginal, 1: Detectable)")
ax1.set_yscale("log")
ax1.tick_params(axis='y', labelcolor=color)

# Create the second y-axis
ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel("Earth Radii", color=color)
ax2.set_yscale("log")

# Match the ticks of ax2 to ax1 based on your R_planet_earth values
ax2.set_yticks(transit_depth)
ax2.set_yticklabels([f"{r:.2f}" for r in R_planet_earth])  # Display corresponding Earth radii
ax2.tick_params(axis='y', labelcolor=color)

plt.title("Contour Plot with Dual Y-Axes (Transit Depth and Earth Radii)")
plt.show()
# %%
