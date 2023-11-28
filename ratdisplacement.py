import numpy as np
from scipy.ndimage import gaussian_filter1d


def rat_displacement(pos_good, all_events, move_dir, samp_freq, scale):
    # 2D ratemap smoothing function
    sg = smooth_function(3, 1.25)

    # px_per_cm = (125/scale)/scale/4

    # First, isolate position data during runs.
    # Yields n x 8 matrix where n is the total number of positions across all runs.
    # Column 1 is the index of the position (from 'pos_good', position data matrix with tracking fixed),
    # column 2 is time, and the columns 3-8 contain x and y values in allocentric
    # space as follows: 3-4 = red light; 5-6 = green light (laser); 7-8=blue light.
    rat_pos_mu = pos_good[:, :2] * scale

    one_hund = samp_freq // 10

    # Calculate heading difference between current position and position 100ms in the future
    r_h = np.unwrap(move_dir[:len(move_dir) - (one_hund), 0])
    r_hf = np.unwrap(move_dir[one_hund:, 0])
    head_diff = circ_dist(r_hf, r_h)
    head_diff *= -1  # Flip heading diffs

    # Compute distance between current position and position 100ms in the future
    x_diff = rat_pos_mu[one_hund:, 0] - rat_pos_mu[:len(rat_pos_mu) - one_hund, 0]
    y_diff = rat_pos_mu[one_hund:, 1] - rat_pos_mu[:len(rat_pos_mu) - one_hund, 1]
    rho_dist = np.sqrt(x_diff ** 2 + y_diff ** 2)

    rat_move_ang_dist = np.column_stack((head_diff, rho_dist))

    return rat_move_ang_dist


def smooth_function(size, sigma):
    # Define a 2D smoothing function (Gaussian)
    x = np.arange(-size // 2 + 1., size // 2 + 1.)
    y = np.arange(-size // 2 + 1., size // 2 + 1.)
    xx, yy = np.meshgrid(x, y, sparse=True)
    kernel = np.exp(-(xx ** 2 + yy ** 2) / (2. * sigma ** 2))
    return kernel / np.sum(kernel)


def circ_dist(alpha, beta):
    # Compute circular distance between two angles
    diff = np.angle(np.exp(1j * alpha) / np.exp(1j * beta))
    return np.mod(diff + np.pi, 2 * np.pi) - np.pi

# Example usage:
# rat_displacement(pos_good, all_events, move_dir, samp_freq, scale)
