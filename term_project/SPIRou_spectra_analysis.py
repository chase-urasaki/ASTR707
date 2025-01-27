#%%
from astropy.io import fits
import numpy 
import matplotlib.pyplot as plt
#%%

from astropy.io import fits

# Open the FITS file
filename = "term_project/2758983s.fits"
with fits.open(filename) as hdul:
    # Print the content summary
    hdul.info()
    
    # Access the binary table
    binary_table = hdul[1].data  # Typically, binary tables are in the second HDU (index 1)

# View column names
print(binary_table.names)

# Access data from a specific column
wl = binary_table['Wave']
flux = binary_table['FluxAB']
error = binary_table['FluxErrAB']

plt.scatter(wl, flux, marker='.', color='black', s=1)
plt.xlim(1082.8, 1083.6)


#%%
