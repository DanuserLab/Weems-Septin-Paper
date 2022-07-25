# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 13:14:17 2022

@author: fyz11
"""


def linear_fit(x,y):
    
    from scipy.stats import linregress
    
    opts = linregress(x,y)
    
    return opts


def exponential_decay_correction_time(signal, time_subset=None, average_fnc=None, f_scale=100., loss_type='soft_l1'):
    """
    Fits an equation of form y=Ae^(Bt)+C on the mean intensities of video
    """
    from scipy.stats import linregress
    from scipy.optimize import least_squares
    import numpy as np 
    
    if average_fnc is None:
        average_fnc = np.mean
        
    I_vid = signal
    
    if time_subset is None:
        # use the full
        I_time = np.arange(len(I_vid))
        I_time_subset = I_time.copy()
        I_vid_subset = I_vid.copy()
    else:
        I_time = np.arange(len(I_vid))
        I_time_subset = I_time[time_subset].copy()
        I_vid_subset = I_vid[time_subset].copy()
        
    """
    Do the fitting on the subset. 
    """
    # fit equation. y =A*e^(-Bt)
    log_I_vid = np.log(I_vid_subset)
    slope, intercept, r_value, p_value, std_err = linregress(I_time_subset, log_I_vid)

    # initial fit. 
    A = np.exp(intercept)
    B = slope
    # refined robust fitting. 
    def exp_decay(t,x):
        return (x[0] * np.exp(x[1] * t) + x[2])
    
    # def exp_decay_line(t,x):    
        # return (x[0] * np.exp(x[1] * (t+x[3])) + x[2]) # exp + linear. (rather than this.... might be more correct to be a switch.... )
    def res(x, t, y):
        # return exp_decay(t,x) - y
        return exp_decay(t,x) - y 
    
    x0 = [A, B, 0] #, 0]
    # print(x0)
    res_robust = least_squares(res, x0, loss=loss_type, f_scale=f_scale, args=(I_time_subset, I_vid_subset))
        
    """
    applying the fitted now on the proper sequence. 
    """
    robust_y = exp_decay(I_time, res_robust.x)
    correction_factor = float(robust_y[0]) / robust_y
    
# #    plt.figure()
# #    plt.plot(robust_y)
# #    plt.show()
#     vid_corrected = np.zeros(vid.shape, np.float32)
    
#     for frame in range(vid.shape[0]):
#         vid_corrected[frame, ...] = vid[frame, ...] * correction[frame]
    
    return correction_factor, (res_robust.x, robust_y)


def baseline_correction_time(signal, p=0.1, lam=1, niters=10):
    """
    Fits an equation of form y=Ae^(Bt)+C on the mean intensities of video
    """
    from scipy.stats import linregress
    from scipy.optimize import least_squares
    import numpy as np 
    
    I_vid = signal
    baseline = baseline_als(I_vid, lam=lam, p=p, niter=niters)

    correction_factor = float(I_vid[0]) / baseline
    corrected = I_vid * correction_factor 
    
    return correction_factor, (baseline, corrected)


def baseline_als(y, lam, p, niter=10):
    from scipy import sparse
    from scipy.sparse.linalg import spsolve
    import numpy as np 
    
    L = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z

def decompose_nonlinear_time_series(y, lam, p, niter=10, padding=None):
    
    import numpy as np 
    
    if padding is None:
        y_ = np.hstack([y[::-1], y, y[::-1]])
        y_base = baseline_als(y_, lam=lam, p=p, niter=niter)
        y_base = y_base[len(y):-len(y)]
    else:
        y_ = np.hstack([y[::-1][-padding:], y, y[::-1][:padding]]) 
        y_base = baseline_als(y_, lam=lam, p=p, niter=niter)
        y_base = y_base[padding:-padding]
        
    return y_base, y-y_base
    

def autocorr(x, norm=True, eps=1e-12):
    import numpy as np 

    if norm: 
        a = (x - np.nanmean(x)) / (np.nanstd(x) * len(x) + eps)
        b = (x - np.nanmean(x)) / (np.nanstd(x) + eps)
    else:
        a = x.copy()
        b = x.copy()
    result = np.correlate(a, b, mode='full') # this is not normalized!. 

    return result[result.size // 2:]

def autocorr_timeseries_set_1d(timeseries_array,norm=True, eps=1e-12):

    import numpy as np 

    autocorr_out = []
    for ii in np.arange(len(timeseries_array)):
        timeseries_ii = timeseries_array[ii]
        autocorr_timeseries_ii = autocorr(timeseries_ii, norm=norm, eps=eps)
        autocorr_out.append(autocorr_timeseries_ii)

    return np.vstack(autocorr_out)

def xcorr(x, y, norm=True, eps=1e-12, mode='full'):
    import numpy as np 

    if norm: 
        a = (x - np.nanmean(x)) / (np.nanstd(x) * len(x) + eps)
        b = (y - np.nanmean(y)) / (np.nanstd(y) + eps)
    else:
        a = x.copy()
        b = y.copy()
    result = np.correlate(a, b, mode=mode) # this is not normalized!. 

    return result

def xcorr_timeseries_set_1d(timeseries_array1, timeseries_array2, norm=True, eps=1e-12, stack_final=True):

    import numpy as np 

    xcorr_out = []
    for ii in np.arange(len(timeseries_array1)):
        timeseries_ii_1 = timeseries_array1[ii].copy()
        timeseries_ii_2 = timeseries_array2[ii].copy()
        xcorr_timeseries_ii = xcorr(timeseries_ii_1, timeseries_ii_2, norm=norm, eps=eps)
        xcorr_out.append(xcorr_timeseries_ii)

    if stack_final:
        return np.vstack(xcorr_out)
    else:
        return xcorr_out
    
def stack_xcorr_curves(xcorr_list):
    
    import numpy as np 
    
    N = [len(xx) for xx in xcorr_list]
    size = np.max(N)
    out_array = np.zeros((len(xcorr_list), size)); out_array[:] = np.nan
    
    for jj in np.arange(len(xcorr_list)):
        xcorr = xcorr_list[jj].copy()
        out_array[jj,size//2-len(xcorr)//2:size//2+len(xcorr)//2+1] = xcorr.copy()
        
    return out_array
    


def spatialcorr_k_neighbors(timeseries_array, k_graph, norm=True, eps=1e-12):
    
    import numpy as np 
    
    z_norm = timeseries_array.copy()
    if norm:
        z_norm = (z_norm-np.nanmean(z_norm, axis=1)[:,None]) / (np.nanstd(z_norm, axis=1)[:,None]+eps)
    
    adj_list = list(k_graph)
    ##### too much memory usage. 
    # N_adj_list = np.hstack([len(aa) for aa in adj_list])
    
    # adj_list_pad = -np.ones((len(N_adj_list), np.max(N_adj_list)), dtype=np.int32)
    # for vv_ii in np.arange(len(adj_list)):
    #     adj_list_pad[vv_ii, :len(adj_list[vv_ii])] = adj_list[vv_ii]
    
    # series1 = z_norm[adj_list_pad].copy()
    # series1_mask = np.ones(adj_list_pad.shape, dtype=bool); series1_mask[adj_list_pad==-1] = 0 
    
    # vertex_means_pearsonr = np.nanmean(z_norm[:,None,:] * series1, axis=-1) 
    # vertex_means_pearsonr = np.nansum(vertex_means_pearsonr*series1_mask, axis=1) / N_adj_list
    
    # all_ring_corrs.append(vertex_means_pearsonr)
    vertex_means_pearsonr = []
    
    # iterate over each vertex.
    for vv_ii in np.arange(len(adj_list)): # this is likely the longest loop - we can paralellize this if we pad.... 
    #     series0 = z_norm[vv_ii].copy() # if norm. 
    #     series1 = z_norm[adj_list[vv_ii]].copy()
        
    #     corrs = np.hstack([spstats.pearsonr(series0, ss)[0] for ss in series1])
    #     vertex_means_pearsonr.append(np.nanmean(corrs))
        
    #     # hard code this to make this fast. 
    # # definition being Sxy / SxxSyy 
    # # https://stackabuse.com/calculating-pearson-correlation-coefficient-in-python-with-numpy/
        
        series0 = z_norm[vv_ii].copy(); series0 = (series0 - np.nanmean(series0)) / (np.nanstd(series0) + eps)
        series1 = z_norm[adj_list[vv_ii]].copy(); series1 = (series1 - np.nanmean(series1, axis=1)[:,None]) / (np.nanstd(series1, axis=1)[:,None] + eps)
        
        Sxy = np.nanmean(series0[None,:] * series1, axis=1) 
        SxxSxy = np.nanstd(series0)  * np.nanstd(series1, axis=1)
        
        corrs = Sxy / (SxxSxy + 1e-12)
        vertex_means_pearsonr.append(np.nanmean(corrs))
        
    vertex_means_pearsonr = np.hstack(vertex_means_pearsonr)

    return vertex_means_pearsonr

