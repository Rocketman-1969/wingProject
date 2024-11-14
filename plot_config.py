import matplotlib.pyplot as plt

def apply_plot_settings(ax, aspect_equal=True, Mx=0.2, My=0.1, minorx=0.05, minory=0.02):
    if aspect_equal:
        ax.set_aspect('equal', adjustable='box')
    ax.grid(True, which='major', linestyle=':', linewidth=0.5)  # Dashed grid lines

    # Major ticks
    ax.xaxis.set_major_locator(plt.MultipleLocator(Mx))
    ax.yaxis.set_major_locator(plt.MultipleLocator(My))
    ax.tick_params(axis='both', which='major', direction='in', length=3)

    # Minor ticks
    ax.xaxis.set_minor_locator(plt.MultipleLocator(minorx))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(minory))
    ax.tick_params(axis='both', which='minor', direction='in', length=1.5, grid_color='none')

    # Add ticks on all sides
    ax.tick_params(top=True, bottom=True, left=True, right=True, which='both')