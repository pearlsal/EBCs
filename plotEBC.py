import numpy as np
import matplotlib.pyplot as plt

def plotEBC(root, out, fign):
    # Plot occupancy circular
    plt.figure(fign)
    ax1 = plt.subplot(1, 4, 1, polar=True)  # Create a polar subplot for occupancy
    t2, r2 = np.meshgrid(np.mod(out['params']['thetaBins'] + np.pi / 2, 2 * np.pi), out['params']['distanceBins'][:-1])
    x, y = np.radians(90 - np.degrees(t2)), r2  # Convert to polar coordinates
    c1 = ax1.pcolormesh(x, y, out['occ'], shading='auto', cmap='viridis')  # Plot occupancy
    plt.title('occ')  # Set subplot title
    plt.colorbar(c1, ax=ax1)  # Add colorbar
    plt.axis('off')  # Turn off axis labels and ticks

    # Plot nspk circular
    ax2 = plt.subplot(1, 4, 2, polar=True)  # Create a polar subplot for spike counts
    c2 = ax2.pcolormesh(x, y, out['nspk'], shading='auto', cmap='viridis')  # Plot spike counts
    plt.title('nspk')  # Set subplot title
    plt.colorbar(c2, ax=ax2)  # Add colorbar
    plt.axis('off')  # Turn off axis labels and ticks

    # Plot ratemap circular
    ax3 = plt.subplot(1, 4, 3, polar=True)  # Create a polar subplot for ratemap
    c3 = ax3.pcolormesh(x, y, out['rm'], shading='auto', cmap='viridis')  # Plot ratemap
    plt.title('rm')  # Set subplot title
    plt.colorbar(c3, ax=ax3)  # Add colorbar
    plt.axis('off')  # Turn off axis labels and ticks

    # Scatter of spike directions
    ax4 = plt.subplot(1, 4, 4)  # Create a regular subplot for trajectory and spike directions
    ax4.scatter(root['x'], root['y'], color='lightgray', edgecolors='none')  # Plot trajectory
    ax4.scatter(out['QP'][:, 0], out['QP'][:, 1], color='black', s=30, marker='o', edgecolors='none')  # Plot edge points
    scatter = ax4.scatter(root['cel_x'], root['cel_y'], c=root['cel_headdir'], cmap='hsv', s=15, alpha=0.7, edgecolors='none')  # Scatter plot of spike directions
    plt.title('Trajectory')  # Set subplot title
    plt.colorbar(scatter, ax=ax4)  # Add colorbar
    plt.axis('off')  # Turn off axis labels and ticks

    plt.show()

# Example usage:
# plotEBC(root_object, result, figure_number)
