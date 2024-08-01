import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.collections as collections
import matplotlib.patches as patches

def smooth_spectrum(wavelengths, flux, sigma):
    """
    Smooth a 1D flux spectrum using a Gaussian filter.
    
    Parameters:
    - wavelengths: array-like, the wavelengths corresponding to the flux values.
    - flux: array-like, the flux values to be smoothed.
    - sigma: float, the standard deviation of the Gaussian kernel used for smoothing.
    
    Returns:
    - smoothed_flux: array-like, the smoothed flux values.
    """
    smoothed_flux = gaussian_filter1d(flux, sigma)
    return smoothed_flux

def plot_fuzzy_cartoon_spectrum(ax, wave, smoothed_flux, cmap, norm, n_layers=30):
    """
    Plot a spectrum with a gradient color and a fuzzy, cartoon-like effect.

    Parameters:
    - ax: Matplotlib axis to plot on.
    - wave: array-like, the wavelengths.
    - smoothed_flux: array-like, the smoothed flux values.
    - cmap: Colormap for the gradient.
    - norm: Normalization for the colormap.
    - n_layers: Number of layers for the fuzzy effect.
    """
    for _ in range(n_layers):
        offset = np.random.uniform(-0.05, 0.05)
        jitter = np.random.normal(0, 0.01, len(smoothed_flux))
        segments = []
        colors = []
        for i in range(len(wave) - 1):
            segments.append([(wave[i], smoothed_flux[i] + offset + jitter[i]), (wave[i+1], smoothed_flux[i+1] + offset + jitter[i+1])])
            colors.append(cmap(norm(wave[i])))

        line_collection = collections.LineCollection(segments, colors=colors, linewidths=2, alpha=0.3)
        ax.add_collection(line_collection)

# Generate the emission line data
data = fits.open("../examples/701230_1d.fits")[1].data
wave = data["WAVE"][0]/1.6698
flux = data["FLUX_REDUCED"][0]

keep = (wave > 4841) & (wave < 5027)
wave = wave[keep]
flux = flux[keep]

smoothed_flux = smooth_spectrum(wave, flux, 2)
smoothed_flux = smoothed_flux/np.max(smoothed_flux)

# Normalize the wavelengths to [0, 1] for colormap mapping
norm = Normalize(vmin=np.min(wave), vmax=np.max(wave))

# Create a colormap with more colors in the gradient
cmap = LinearSegmentedColormap.from_list('custom_grad', ['purple', 'blue', 'cyan', 'green', 'orange', 'red', 'crimson'])

# Create a figure and axis
fig, ax = plt.subplots(figsize=(12, 3))

# Add text for the name "quickLines"
ax.text(4842, 0.9, r'quickLines', fontsize=36, color='#00BFFF', ha='left', va='center', fontweight='bold', fontname='Comic Sans MS')
ax.text(4850, 0.7, r'On-the-Fly Emission Line Property Extractor', fontsize=24, color='#00BFFF', ha='left', va='center', fontstyle='italic', fontname='DejaVu Sans')


# Plot the spectrum with fuzzy cartoon effect
plot_fuzzy_cartoon_spectrum(ax, wave, smoothed_flux, cmap, norm)


# Set axis limits and remove axes
ax.set_xlim(np.min(wave),np.max(wave))
ax.set_ylim(0, 1.05)
ax.axis('off')

# Adjust layout to minimize whitespace
plt.tight_layout(pad=0.1, rect=[0, 0, 1, 1])  # Adjust rect to fit the figure content

# Save the figure
fig.savefig("quicklines_logo.png", dpi=300, bbox_inches='tight', transparent=True, pad_inches=0)
