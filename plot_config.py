import matplotlib.pyplot as plt

def apply_plot_settings(aspect_equal=True):
    if aspect_equal:
        plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True, which='major', linestyle=':', linewidth=0.5)  # Dashed grid lines

    # Major ticks
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(0.2))
    plt.gca().yaxis.set_major_locator(plt.MultipleLocator(0.1))
    plt.gca().tick_params(axis='both', which='major', direction='in', length=3)

    # Minor ticks
    plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(0.02))
    plt.gca().tick_params(axis='both', which='minor', direction='in', length=1.5, grid_color='none')

    # Add ticks on all sides
    plt.gca().tick_params(top=True, bottom=True, left=True, right=True, which='both')