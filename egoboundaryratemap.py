import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import convolve, gaussian_filter1d
from scipy.stats import norm

def egocentric_ratemap(root, video_samp=1, deg_samp=1, distance_bins=None, boundary_mode=0, smooth_kernel=[5, 5, 5]):
    # Function to calculate the egocentric ratemap
    # and return relevant information for plotting or statistical testing.

    # Setup and parse parameters
    if distance_bins is None:
        distance_bins = np.arange(0, 51)

    rx, ry, md, ts, spk = extract_behavioral_information(root)

    # Get structure of the environment
    if boundary_mode == 0:
        QP = auto_detect_edges(rx, ry)
    elif boundary_mode == 1:
        QP = find_edges(root)
    else:
        QP = boundary_mode

    # Calculate distances
    dis, ex, ey = calc_distance(rx, ry, md, QP, deg_samp)

    # Calculate raw maps
    theta_bins = np.deg2rad(np.linspace(-180, 180, dis.shape[1]))
    occ, nspk = calculate_occ_nspk(dis, spk, distance_bins)

    # Smoothing
    occ_smoothed = smooth_mat(occ, smooth_kernel[:2], smooth_kernel[2])
    nspk_smoothed = smooth_mat(nspk, smooth_kernel[:2], smooth_kernel[2])

    # Ratemaps
    rm_ns = nspk / occ
    rm = nspk_smoothed / occ_smoothed

    # Package the output
    out = {
        'rm_ns': rm_ns,
        'occ_ns': occ,
        'nspk_ns': nspk,
        'occ': occ_smoothed,
        'nspk': nspk_smoothed,
        'rm': rm,
        'QP': QP,
        'params': {
            'videoSamp': video_samp,
            'degSamp': deg_samp,
            'distanceBins': distance_bins[:-1],
            'smoothKernel': smooth_kernel,
            'thetaBins': theta_bins
        }
    }

    return out


def extract_behavioral_information(root):
    # Extract relevant behavioral information from the root object
    rx = root.x
    ry = root.y
    md = root.headdir
    if np.max(md) > 2 * np.pi:
        md = np.deg2rad(md)
    ts = root.ts
    sts = root.cel_ts
    spk = np.histogram(sts, ts)[0]
    return rx, ry, md, ts, spk


def auto_detect_edges(rx, ry):
    # Automatically detect the edges of the environment
    p = [-1000, -1000]
    d = (rx - p[0]) ** 2 + (ry - p[1]) ** 2
    ind = np.argmin(d)
    ll = [rx[ind], ry[ind]]

    p = [1000, -1000]
    d = (rx - p[0]) ** 2 + (ry - p[1]) ** 2
    ind = np.argmin(d)
    lr = [rx[ind], ry[ind]]

    p = [1000, 1000]
    d = (rx - p[0]) ** 2 + (ry - p[1]) ** 2
    ind = np.argmin(d)
    ur = [rx[ind], ry[ind]]

    p = [-1000, 1000]
    d = (rx - p[0]) ** 2 + (ry - p[1]) ** 2
    ind = np.argmin(d)
    ul = [rx[ind], ry[ind]]

    QP = np.array([ll, lr, ur, ul])
    return QP


def find_edges(root):
    # Manually select corners of walls to find edges
    if_escape = 0
    h = plt.figure()

    while not if_escape:
        plt.figure(h)
        plt.clf()
        plt.gca().invert_yaxis()
        clim = plt.gca().get_clim()
        plt.gca().set_clim(clim / 50)
        plt.plot(root.x, root.y, 'k')
        QP = []

        plt.title('Select Corners of Walls. Esc--> done. **Do not complete!**')

        button = 1

        while button != 27:
            x, y, button = plt.ginput(1)[0]

            plt.clf()
            plt.gca().invert_yaxis()
            clim = plt.gca().get_clim()
            plt.gca().set_clim(clim / 50)
            plt.plot(root.x, root.y, 'k')

            if QP:
                plt.plot(QP[:, 0], QP[:, 1], 'r')
                plt.plot(QP[:, 0], QP[:, 1], 'ro', markerfacecolor='r')

            if button == 32:  # space bar
                QP.append([np.nan, np.nan])
            elif button != 27:
                QP.append([x, y])

            plt.plot(QP[:, 0], QP[:, 1], 'r')
            plt.plot(QP[:, 0], QP[:, 1], 'ro', markerfacecolor='r')

        # Ask for verification
        edg = splitter(QP)
        plt.clf()
        plt.title('Verify. 0--> Try again; 1--> Confirm')
        plt.plot(root.x, root.y, 'k')

        for m in range(len(edg)):
            for n in range(edg[m].shape[0]):
                sp = edg[m][n, :, 0]
                ep = edg[m][n, :, 1]
                plt.plot([sp[0], ep[0]], [sp[1], ep[1]], 'ro', markerfacecolor='r')
                plt.plot([sp[0], ep[0]], [sp[1], ep[1]], 'r')

        # Set or repeat
        while button != 48 and button != 49:
            _, _, button = plt.ginput(1)[0]

        if_escape = button == 49

    plt.close(h)
    plt.draw()
    return QP


def splitter(QP):
    # Split corners
    inds = np.where(np.isnan(QP[:, 0]))[0]
    xs = np.split(QP[:, 0], inds)
    ys = np.split(QP[:, 1], inds)

    QP2 = [np.column_stack((x, y)) for x, y in zip(xs, ys) if x.size > 0]

    edg = []
    for m in range(len(QP2)):
        edge_m = []
        for n in range(QP2[m].shape[0]):
            sp = n
            ep = n + 1
            if ep >= QP2[m].shape[0]:
                ep = 0
            edge_m.append(QP2[m][[sp, ep], :])
        edg.append(np.array(edge_m))

    return edg


def calc_distance(rx, ry, md, QP, deg_samp):
    # Calculate distances between points and edges
    mxd = np.sqrt((np.max(rx) - np.min(rx)) ** 2 + (np.max(ry) - np.min(ry)) ** 2)
    degs = np.deg2rad(np.arange(-180, 181, deg_samp))

    edg = splitter(QP)
    edg = np.vstack(edg)
    dis = np.full((len(rx), edg.shape[0], len(degs)), np.nan)
    dir = dis.copy()

    for i in range(edg.shape[0]):
        x1, x2 = edg[i, 0, 0], edg[i, 0, 1]
        y1, y2 = edg[i, 1, 0], edg[i, 1, 1]

        for h in range(len(degs)):
            mdof = degs[h]
            y3, x3 = ry, rx
            y4, x4 = ry + mxd * np.sin(md + mdof), rx + mxd * np.cos(md + mdof)

            px1 = (x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)
            px2 = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
            px = px1 / px2

            py1 = (x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)
            py2 = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
            py = py1 / py2

            d = np.sqrt((ry - py) ** 2 + (rx - px) ** 2)
            dis[:, i, h] = d

            dir[:, i, h] = np.angle(py - ry + 1j * (px - rx)) - (md + mdof)

            bb = [np.min([x1, x2]), np.max([x1, x2]), np.min([y1, y2]), np.max([y1, y2])]
            indexes = ~(px >= bb[0] and px <= bb[1] and py >= bb[2] and py <= bb[3])
            dis[indexes, i, h] = np.nan

    dis[dis > mxd] = np.nan
    dis[np.abs(dir) > np.pi / 4] = np.nan

    dis = np.nanmin(dis, axis=1)
    dd = np.tile(degs, (len(rx), 1)) + np.tile(md, (len(degs), 1)).T
    dx = dis * np.cos(dd)
    dy = dis * np.sin(dd)
    ey = dy + np.tile(ry, (len(degs), 1)).T
    ex = dx + np.tile(rx, (len(degs), 1)).T

    return dis, ex, ey


def calculate_occ_nspk(dis, spk, distance_bins):
    # Calculate occupancy and spike counts
    theta_bins = np.deg2rad(np.linspace(-180, 180, dis.shape[1]))
    occ = np.full((len(theta_bins), len(distance_bins) - 1), np.nan)
    nspk = occ.copy()

    ci = np.where(spk)[0]

    for i in range(len(theta_bins)):
        t = dis[:, i]
        for k in range(len(distance_bins) - 1):
            inds = np.where((t >= distance_bins[k]) & (t < distance_bins[k + 1]))[0]
            occ[i, k] = len(inds)
            inds = np.intersect1d(inds, ci)
            nspk[i, k] = len(inds)

    return occ, nspk


def smooth_mat(mat, kernel_size, std):
    # Smooth matrix by convolving with 2D Gaussian
    if std == 0:
        return mat

    Xgrid, Ygrid = np.meshgrid(np.arange(-kernel_size[0] / 2, kernel_size[0] / 2 + 1),
                               np.arange(-kernel_size[1] / 2, kernel_size[1] / 2 + 1))
    Rgrid = np.sqrt(Xgrid ** 2 + Ygrid ** 2)

    kernel = norm.pdf(Rgrid, 0, std)
    kernel = kernel / np.sum(kernel)
    mat = convolve(mat, kernel, mode='constant')

    return mat

# Example usage:
# result = egocentric_ratemap(root_object)
