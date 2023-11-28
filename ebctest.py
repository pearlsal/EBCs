import numpy as np
from scipy.interpolate import interp1d
from scipy.io import loadmat
import Pippin  # Assuming Pippin is a custom module or package
#Please note that I made some assumptions about missing functions (Resample, rat_displacement, etc.)
# and imported packages (Pippin). You may need to replace or modify these parts based on your actual implementation.

def ebc_test(root, res=None, QP=None):
    # Disable ill-conditioned warnings during glmfit
    s = np.seterr(divide='ignore', invalid='ignore')

    # Resample if 'res' is provided
    if res is not None and len(res) > 0:
        cel = root.cel
        root = Resample(root, res)
        root.cel = cel

    # Extract relevant data
    ep = root.epoch
    root.epoch = [-np.inf, np.inf]
    xvec = np.nan_to_num(interp1d(np.arange(len(root.b_x)), root.b_x, kind='spline')(np.arange(len(root.b_x))))
    yvec = np.nan_to_num(interp1d(np.arange(len(root.b_y)), root.b_y, kind='spline')(np.arange(len(root.b_y))))

    p = CMBHOME.Utils.speed_kalman(xvec, yvec, root.fs_video)
    md = p.vd
    root.epoch = ep

    # If QP is not provided, define it based on bounding box corners
    if QP is None:
        brx = root.b_x
        bry = root.b_y

        p = [-1000, -1000]
        d = (brx - p[0])**2 + (bry - p[1])**2
        ind = np.argmin(d)
        ll = [brx[ind], bry[ind]]

        p = [1000, -1000]
        d = (brx - p[0])**2 + (bry - p[1])**2
        ind = np.argmin(d)
        lr = [brx[ind], bry[ind]]

        p = [1000, 1000]
        d = (brx - p[0])**2 + (bry - p[1])**2
        ind = np.argmin(d)
        ur = [brx[ind], bry[ind]]

        p = [-1000, 1000]
        d = (brx - p[0])**2 + (bry - p[1])**2
        ind = np.argmin(d)
        ul = [brx[ind], bry[ind]]

        QP = np.array([ll, lr, ur, ul])

    # Get self-motion predictors
    temp_disp = rat_displacement(np.column_stack((root.x, root.y)), [1, len(root.ind)],
                                  md[root.ind], root.fs_video, root.spatial_scale)
    temp_disp[:, 0] = np.nan_to_num(interp1d(np.arange(len(temp_disp[:, 0])), temp_disp[:, 0], kind='linear')(np.arange(len(temp_disp[:, 0]))))
    temp_disp[:, 1] = np.nan_to_num(interp1d(np.arange(len(temp_disp[:, 1])), temp_disp[:, 1], kind='linear')(np.arange(len(temp_disp[:, 1]))))
    rat_move_ang_dist = np.zeros((len(root.ind), 2))
    rat_move_ang_dist[:len(temp_disp), :] = temp_disp

    # Bearing and distance to center of the environment
    xloc = ((np.max(QP[:, 0]) - np.min(QP[:, 0])) / 2) + np.min(QP[:, 0])
    yloc = ((np.max(QP[:, 1]) - np.min(QP[:, 1])) / 2) + np.min(QP[:, 1])

    center_dis = np.sqrt((xvec - xloc)**2 + (yvec - yloc)**2)
    center_dis = center_dis[root.ind]

    hd = np.concatenate([[0], np.arctan2(np.diff(yvec), np.diff(xvec))])
    hd = np.unwrap(hd)
    hd = hd[root.ind]

    center_ang = np.arctan2(yloc - yvec[root.ind], xloc - xvec[root.ind]) - hd
    center_ang = np.unwrap(center_ang)

    # Build up model
    model = Pippin.Model(root)
    model = Pippin.Predictors.Other(model, 'Movementdirection', [np.cos(md[root.ind]), np.sin(md[root.ind])])
    model = Pippin.Predictors.Place(model)
    model = Pippin.Predictors.Speed(model)
    model = Pippin.Predictors.Other(model, 'angularDisplacement', rat_move_ang_dist[:, 0])
    model = Pippin.Predictors.Other(model, 'centerDis', [center_dis, center_dis**2])
    model = Pippin.Predictors.Other(model, 'centerAng', [np.sin(center_ang), np.cos(center_ang)])
    model = Pippin.Predictors.Other(model, 'centerConj', [np.sin(center_ang) * center_dis, np.cos(center_ang) * center_dis])

    model.gen_models()
    glm_summary = model.summary()

    # No Conjunctive
    model_no_conj = Pippin.Model(root)
    model_no_conj = Pippin.Predictors.Other(model_no_conj, 'Movementdirection', [np.cos(md[root.ind]), np.sin(md[root.ind])])
    model_no_conj = Pippin.Predictors.Place(model_no_conj)
    model_no_conj = Pippin.Predictors.Speed(model_no_conj)
    model_no_conj = Pippin.Predictors.Other(model_no_conj, 'angularDisplacement', rat_move_ang_dist[:, 0])
    model_no_conj = Pippin.Predictors.Other(model_no_conj, 'centerDis', [center_dis, center_dis**2])
    model_no_conj = Pippin.Predictors.Other(model_no_conj, 'centerAng', [np.sin(center_ang), np.cos(center_ang)])

    model_no_conj.gen_models()

    # Simplified Model Comparison
    allo = np.column_stack((model.predictors[1].data, model.predictors[2].data))
    move = np.column_stack((model.predictors[3].data, model.predictors[4].data))
    ego = np.column_stack((model.predictors[5].data, model.predictors[6].data, model.predictors[7].data))

    model_simp = Pippin.Model(root)
    model_simp = Pippin.Predictors.Other(model_simp, 'allo', allo)
    model_simp = Pippin.Predictors.Other(model_simp, 'move', move)
    model_simp = Pippin.Predictors.Other(model_simp, 'ego', ego)

    model_simp.gen_models()

    # Reset warning settings
    np.seterr(**s)

    return model, glm_summary, model_no_conj, model_simp
