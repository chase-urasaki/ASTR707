import matplotlib.pyplot as plt
import numpy as np

def plot_tess_monte_carlo(t, fdetrend, d, time, mom_dump, t0, tc0, starName, sector, t1, t2, iterations=1000):
    for i in range(iterations):
        which = ((t > t1) & (t < t2) & (fdetrend > 0.))
        
        fig = plt.figure(figsize=(15, 10))
        
        fig.add_subplot(4, 4, 12)
        plt.ylim(min(fdetrend[which]), max(fdetrend[which]))
        plt.xlim([t1, t2])
        plt.vlines([tc0], 0.8 * min(fdetrend[which]), 1.2 * max(fdetrend[which]), colors='green', linestyle='dotted')
        plt.plot(t[which], fdetrend[which], '.', markersize=1)
        plt.vlines(time[mom_dump] - t0, 0.8 * min(fdetrend[which]), 1.2 * max(fdetrend[which]), colors='r', linestyle='dotted')
        
        fig.add_subplot(4, 4, 16)
        plt.xlim([t1, t2])
        plt.plot(t[which], d[which], '.', markersize=1)
        plt.ylim(0, max(d[which]))
        plt.vlines(time[mom_dump] - t0, 0, max(d[which]), colors='r', linestyle='dotted')
        plt.xlabel('days', fontsize=9)
        plt.ylabel('offset (pixels)', fontsize=9)
        
        plt.subplots_adjust(bottom=0.1, right=0.97, top=0.95, wspace=0.3, hspace=0.3, left=0.05)
        plotname = f"{starName.replace(' ', '_')}_S{sector}_iter{i}.png"
        plt.savefig(plotname)
        plt.close(fig)

# Example usage:
# plot_tess_monte_carlo(t, fdetrend, d, time, mom_dump, t0, tc0, starName, sector, t1, t2, iterations=1000)