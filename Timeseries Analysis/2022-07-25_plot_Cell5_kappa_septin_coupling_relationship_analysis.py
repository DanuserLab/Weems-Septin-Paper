# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 01:39:53 2022

@author: fyz11
"""

if __name__=="__main__":
    
    import numpy as np
    import pylab as plt 
    import scipy.io as spio
    import glob 
    import os 
    import trimesh 
    import scipy.ndimage as ndimage 
    import skimage.morphology as skmorph
    import skimage.filters as skfilters
    import scipy.stats as spstats
    import seaborn as sns
    from matplotlib import cm 
    
    import LightSheet_Analysis.Segmentation.segmentation as segmentation 
    import LightSheet_Analysis.Utility_Functions.file_io as fio
    import LightSheet_Analysis.Image_Functions.image as image_fn
    import LightSheet_Analysis.Mesh.meshtools as meshtools
    import LightSheet_Analysis.Parameters.params as params
    import LightSheet_Analysis.Registration.registration_new as registration 
    import LightSheet_Analysis.Analysis_Functions.timeseries as tsa
    import LightSheet_Analysis.Visualisation.colors as vol_colors
    
    from tqdm import tqdm 
    import seaborn as sns 
    import pandas as pd 
    import igl
    import scipy.stats as spstats
# =============================================================================
#     Boilerplate for loading folders and files. 
# =============================================================================

#     """
#     Dataset 1 
#         - canonical cell here is the Cell5 in this dataset. 
#     """
#     # datafolder = '/archive/bioinformatics/Danuser_lab/melanoma/raw/hiRes3D(cheap knockoff)/Voodoo/MV3/Septin6/180522/Cell2'
#     datafolder = '/archive/bioinformatics/Danuser_lab/melanoma/raw/hiRes3D(cheap knockoff)/Voodoo/MV3/Septin6/180522/Cell5'
#     # datafolder = '/archive/bioinformatics/Danuser_lab/melanoma/raw/hiRes3D(cheap knockoff)/Voodoo/MV3/Septin6/180522/Cell7'
    
    
    
#     imfiles = np.sort(np.hstack(glob.glob(os.path.join(datafolder, '*.tif'))))
#     imfiles = np.sort(np.hstack([ff for ff in imfiles if '1_CH00' in ff]))
#     imfiles = imfiles[:-1]
    
#     # splitpath = '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/raw'
#     splitpath = '/archive/bioinformatics/Danuser_lab/melanoma/raw/hiRes3D(cheap knockoff)/Voodoo'
#     # outfolder = '/home2/s205272/Documents/Projects/Results'
#     outfolder = '/work/bioinformatics/s205272/3D_Causality/Results'
#     # outfolder = '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Weems/Septin_Blebs'
#     # savefolder = datafolder.replace(splitpath, outfolder)
#     savefolder = datafolder.replace(splitpath, outfolder)
    
#     # print(savefolder)
#     # # mkdir(savefolder)
    
#     # saveoutfolder = os.path.join(savefolder,'rigid_registered')
#     saveoutfolder = os.path.join(savefolder,'demons_registered_diffeomorphic_Full-crop')
#     # saveoutfolder = os.path.join(savefolder, 'demons_registered_diffeomorphic_Full')
#     # saveoutfolder = os.path.join(savefolder,'demons_registered_mid_Sym_Demon_Full') # this is the full 
#     # mkdir(saveoutfolder)

#     # # load the crop coordinates file. 
#     # bbox_params_file = os.path.join(savefolder, os.path.split(datafolder)[-1]+'_dataset_bbox_params.mat')
#     # bbox_params = spio.loadmat(bbox_params_file)['crop_bbox_global_xyz']


#     # we use the final version. 
#     saveanalysisfolder = os.path.join(savefolder,'demons_registered_diffeomorphic_Full-crop-Analysis-Final')

#     analysis_savefolder = os.path.join(saveanalysisfolder, 'Curvature_Septin_Analysis')
#     fio.mkdir(analysis_savefolder) 
    

# # =============================================================================
# # Iterating over each file
# # =============================================================================
    
#     # # px_res = 0.104
#     # mesh = trimesh.load_mesh('E:\\Work\\Projects\\Danuser-3D Causality\\Data\\Andrew\\Septin_Blebs2\\180522\\Cell5\\Surface_Register_Kappa\\Frame_0000.obj', 
#     #                          process=False,
#     #                          validate=False)
        
#     statsfolder = os.path.join(saveanalysisfolder, 
#                                'demons_prop_mesh-Euler-Septin_FirstFrame'); 

# =============================================================================
# =============================================================================
# # Setup the final plot savefolder.  
# =============================================================================
# =============================================================================

    # Ts = 1.00  # give the temporal sampling. 
    Ts = 1.21 # this is the correct version. 

    analysis_output_folder = 'E:\\Work\\Projects\\Danuser-3D Causality\\Data\\Andrew\\Septin_Blebs2\\180522\\Cell5\\2022-03-09_Registration\\2022-04-05_Analysis\\plots'
    fio.mkdir(analysis_output_folder)

    #########

# =============================================================================
#     Load up the unprocessed statistics ( sampled intensities )
# =============================================================================
    # statsfolder = 'E:\\Work\\Projects\\Danuser-3D Causality\\Data\\Andrew\\Septin_Blebs2\\180522\\Cell5\\2022-01-27_FirstOrder_Correlation_Analysis_testing'#'\\Cell2'
    # statsfolder = 'E:\\Work\\Projects\\Danuser-3D Causality\\Data\\Andrew\\Septin_Blebs2\\180522\\Cell5\\2022-03-09_Registration\\2022-03-10_Analysis'
    
    statsfolder = 'E:\\Work\\Projects\\Danuser-3D Causality\\Data\\Andrew\\Septin_Blebs2\\180522\\Cell5\\2022-03-09_Registration\\2022-04-05_Analysis'
    statsfile = os.path.join(statsfolder, 'Septin-Curvature_statistics_FirstFrame_DistTform.mat')
    # statsfile = 'Septin-Curvature_statistics.mat'
    stats = spio.loadmat(statsfile)
    
    print(stats.keys())
    
    px_res = float(stats['pix_res'])
    septin_vol_intensity_all = np.hstack(stats['volume_intensity_all'])
    septin_faces_all = stats['septin_all_depth1']
    curvature_faces_all_um_img = stats['H_all'] / px_res # use this one as this is much less noisy # this is faces!. 
    curvature_faces_all_um_mesh = stats['H_mesh_all'] / px_res
    
    I_min = float(stats['septin_I_min_depth1'])
    I_max = float(stats['septin_I_max_depth1'])
    
    v = stats['v_all'] 
    f = stats['f']
    mesh_v = v.copy()
    mesh_f = f.copy()
    
    # correct for septin_faces_all
    septin_faces_all = septin_faces_all / septin_vol_intensity_all[:,None] #* septin_vol_intensity_all[0]
    
    
# # =============================================================================
# #     Compute additional curvature needed for analysis. 
# # =============================================================================
#     # # compute these and smooth locally!. 
#     # principal_curvatures = [igl.principal_curvature(vv[...,::-1], f)[:] for vv in v]
#     # #  # separate into the directional fields. 
#     # principal_curvature_components_vectors = np.array( [cc[:2] for cc in principal_curvatures], dtype=np.float32)
#     # principal_curvatures = np.array([ cc[2:] for cc in principal_curvatures], dtype=np.float32)
    
#     # we should probably save this to avoid the expensive computation. 
#     # spio.savemat(os.path.join(statsfolder, 
#     #                       'save_principal_curvatures_statistics.mat'), 
#     #               {'prin_kappa': principal_curvatures})
#     principal_curvatures = spio.loadmat(os.path.join(statsfolder, 
#                           'save_principal_curvatures_statistics.mat'))['prin_kappa']
    
# =============================================================================
# =============================================================================
# #     Load the smoothed signals we got.... 
# =============================================================================
# =============================================================================
    # 'SpatialTemp_correlation_stats_septin_curvature.mat',
    statsfolder = 'E:\\Work\\Projects\\Danuser-3D Causality\\Data\\Andrew\\Septin_Blebs2\\180522\\Cell5\\2022-03-09_Registration\\2022-04-05_Analysis\\Global_Analyses'
    smooth_signals_obj = spio.loadmat(os.path.join(statsfolder, 'SpatialTemp_correlation_stats_septin_curvature.mat' )) # why have i not saved the different thresholds used? ---- this isa bug 
    
    # use these to identify what we need. (sample the signal), including the curvature onset transitions!. 
    curvature_faces_all_smooth = smooth_signals_obj['curvature_faces_all_smooth']
    septin_faces_all_smooth = smooth_signals_obj['septin_faces_all_smooth']
    dseptin_faces_all_smooth = smooth_signals_obj['dseptin_faces_all_smooth']
    
    
    #### average onto face connectivity now. 
    curvature_faces_all_smooth_faces = igl.average_onto_faces(f, curvature_faces_all_smooth.T).T
    septin_faces_all_smooth_faces = igl.average_onto_faces(f, septin_faces_all_smooth.T).T
    dseptin_faces_all_smooth_faces = igl.average_onto_faces(f, dseptin_faces_all_smooth.T).T
    
    
    #### rederive the raw signals (the above is wrong - there is offset timepoint!)
    septin_faces_all_raw = igl.average_onto_faces(f, septin_faces_all.T).T
    kappa_all_raw = igl.average_onto_faces(f, curvature_faces_all_um_mesh.T).T
    kappa_all_img_vertices = [igl.average_onto_vertices(v[ttt], f, np.vstack([curvature_faces_all_um_img[ttt],curvature_faces_all_um_img[ttt],curvature_faces_all_um_img[ttt]]).T)[:,0] for ttt in range( len(v) ) ]
    kappa_all_img_vertices = np.array(kappa_all_img_vertices)
    dseptin_faces_all_raw = np.gradient(septin_faces_all_raw, axis=0)
    
    

# =============================================================================
#     We need to identify the very flat microscope limited regions!. 
# =============================================================================
# =============================================================================
#     Load the computed geodesics and segmentation from the windowed analysis!.  - checking the kappa groups. 
# =============================================================================
    
    import scipy.stats as spstats
    
    statsfolder = 'E:\\Work\\Projects\\Danuser-3D Causality\\Data\\Andrew\\Septin_Blebs2\\180522\\Cell5\\2022-03-09_Registration\\2022-04-05_Analysis\\Window_Analyses'
    infile = os.path.join(statsfolder, 
                          'save_curvature-septin_statistics.mat')
    stats_obj = spio.loadmat(infile)
    
    all_kappa_groups_faces = stats_obj['all_Kappa_groups'].copy()
    
    # the above is for individual faces -> now we make transfer this to overall equilibrium behaviour
    all_kappa_groups_class, all_kappa_groups_class_freq = spstats.mode(all_kappa_groups_faces, axis=0)
    all_kappa_groups_class = np.squeeze(all_kappa_groups_class) 
    all_kappa_groups_class_freq = np.squeeze(all_kappa_groups_class_freq)
    
    multi_thresh = skfilters.threshold_multiotsu(curvature_faces_all_um_img.mean(axis=0))
    all_kappa_groups_class = np.zeros(len(all_kappa_groups_class), dtype=np.int)
    all_kappa_groups_class[np.logical_and(curvature_faces_all_um_img.mean(axis=0)>=multi_thresh[0], 
                                          curvature_faces_all_um_img.mean(axis=0)<=multi_thresh[1])] = 1
    all_kappa_groups_class[curvature_faces_all_um_img.mean(axis=0)>multi_thresh[1]] = 2
    
    
    
    fig, ax = plt.subplots(figsize=(10,10))
    # plt.title(str(rot)+'_'+str(rot2))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_proj_type('ortho')
    ax.set_box_aspect(aspect = (1,1,1)) # this works. 
    
    ax.scatter(v[0][:,0],
                v[0][:,1],
                v[0][:,2], s=.1, c='k', alpha=.5)
    
    pts = v[0][np.unique(f[all_kappa_groups_class==0])]
    ax.scatter(pts[:,0],
                pts[:,1],
               pts[:,2], s=.1, c='r', cmap='Spectral', vmin=-1, vmax=1)
    # ax.scatter(v[0][:,0],
    #            v[0][:,1],
    #            v[0][:,2], s=.1, c=kappa_all_img_vertices.mean(axis=0), cmap='Spectral', vmin=-1, vmax=1)
    ax.axis('off')
    ax.grid('off')
    ax.view_init(135,210)
    # uzip.set_axes_equal(ax)
    # fig.savefig('curvature_uv_smooth.svg', dpi=300, bbox_inches='tight')
    plt.show()
    
# =============================================================================
#     ok the above is difficult to get the exclusion criteria. so lets solve for the eigenaxis and walk along this.... =_=
# =============================================================================
    curvature_weights = np.abs(kappa_all_img_vertices.mean(axis=0))
    # curvature_weights[curvature_weights>multi_thresh[0]] = 0 
    # curvature_weights = all_kappa_groups_class==1; curvature_weights=curvature_weights*1.
    
    w_axis, v_axis = meshtools.find_principal_axes_surface_heterogeneity_mesh(v[0], 
                                                                              f,
                                                                              curvature_weights, # optimise for the original curvature
                                                                              map_to_sphere=False)
    area_weights = igl.doublearea(v[0],f)
    centroid = np.nansum(area_weights[:,None]*igl.barycenter(v[0],f), axis=0)/ np.nansum(area_weights)
    
    plot_length=100
    
    fig, ax = plt.subplots(figsize=(10,10))
    # plt.title(str(rot)+'_'+str(rot2))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_proj_type('ortho')
    ax.set_box_aspect(aspect = (1,1,1)) # this works. 
    
    ax.scatter(v[0][:,0],
                v[0][:,1],
                v[0][:,2], s=.1, c='k', alpha=.1)
    
    # pts = v[0][np.unique(f[all_kappa_groups_class==0])]
    # ax.scatter(pts[:,0],
    #             pts[:,1],
    #             pts[:,2], s=.1, c='r', cmap='Spectral', vmin=-1, vmax=1, alpha=.5)
    ax.plot([centroid[0]-plot_length*v_axis[0,0], centroid[0] + plot_length*v_axis[0,0]],
            [centroid[1]-plot_length*v_axis[1,0], centroid[1] + plot_length*v_axis[1,0]],
            [centroid[2]-plot_length*v_axis[2,0], centroid[2] + plot_length*v_axis[2,0]],'r', zorder=100., lw=10)
    ax.plot([centroid[0]-plot_length*v_axis[0,1], centroid[0] + plot_length*v_axis[0,1]],
            [centroid[1]-plot_length*v_axis[1,1], centroid[1] + plot_length*v_axis[1,1]],
            [centroid[2]-plot_length*v_axis[2,1], centroid[2] + plot_length*v_axis[2,1]],'g', zorder=100., lw=10)
    ax.plot([centroid[0]-plot_length*v_axis[0,2], centroid[0] + plot_length*v_axis[0,2]],
            [centroid[1]-plot_length*v_axis[1,2], centroid[1] + plot_length*v_axis[1,2]],
            [centroid[2]-plot_length*v_axis[2,2], centroid[2] + plot_length*v_axis[2,2]],'b', zorder=100., lw=10)
    
    ax.axis('off')
    ax.grid('off')
    ax.view_init(135,210)
    # uzip.set_axes_equal(ax)
    # fig.savefig('curvature_uv_smooth.svg', dpi=300, bbox_inches='tight')
    plt.show()
    
    
    fig, ax = plt.subplots(figsize=(10,10))
    # plt.title(str(rot)+'_'+str(rot2))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_proj_type('ortho')
    ax.set_box_aspect(aspect = (1,1,1)) # this works. 
    
    ax.scatter(v[0][:,0],
                v[0][:,1],
                v[0][:,2], s=.1, c='k', alpha=.1)
    
    # pts = v[0][np.unique(f[all_kappa_groups_class==0])]
    # ax.scatter(pts[:,0],
    #             pts[:,1],
    #             pts[:,2], s=.1, c='r', cmap='Spectral', vmin=-1, vmax=1, alpha=.5)
    ax.plot([centroid[0]-plot_length*v_axis[0,0], centroid[0] + plot_length*v_axis[0,0]],
            [centroid[1]-plot_length*v_axis[1,0], centroid[1] + plot_length*v_axis[1,0]],
            [centroid[2]-plot_length*v_axis[2,0], centroid[2] + plot_length*v_axis[2,0]],'r', zorder=100., lw=10)
    ax.plot([centroid[0]-plot_length*v_axis[0,1], centroid[0] + plot_length*v_axis[0,1]],
            [centroid[1]-plot_length*v_axis[1,1], centroid[1] + plot_length*v_axis[1,1]],
            [centroid[2]-plot_length*v_axis[2,1], centroid[2] + plot_length*v_axis[2,1]],'g', zorder=100., lw=10)
    ax.plot([centroid[0]-plot_length*v_axis[0,2], centroid[0] + plot_length*v_axis[0,2]],
            [centroid[1]-plot_length*v_axis[1,2], centroid[1] + plot_length*v_axis[1,2]],
            [centroid[2]-plot_length*v_axis[2,2], centroid[2] + plot_length*v_axis[2,2]],'b', zorder=100., lw=10)
    
    ax.axis('off')
    ax.grid('off')
    ax.view_init(-45,30)
    # uzip.set_axes_equal(ax)
    # fig.savefig('curvature_uv_smooth.svg', dpi=300, bbox_inches='tight')
    plt.show()
    
    
    
# =============================================================================
#     Try projecting and excluding ..... based on the normal distance to the centroid plane! using the red principal axis!. !. 
# =============================================================================
    # we need to reorient the axis. 
    proj_distances = np.nansum((v[0]-centroid[None,:]) * v_axis[:,0][None,:], axis=-1)
    mean_pos = np.median(curvature_weights[np.sign(proj_distances)>0])
    mean_neg = np.median(curvature_weights[np.sign(proj_distances)<0])
    
    flip_eig = 0
    if mean_neg > mean_pos:
        proj_distances = -1*proj_distances
        flip_eig = 1
        
    thresh_percentile = 25
    thresh_proj_distance = np.percentile(proj_distances,thresh_percentile)
        
    fig, ax = plt.subplots(figsize=(10,10))
    # plt.title(str(rot)+'_'+str(rot2))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_proj_type('ortho')
    ax.set_box_aspect(aspect = (1,1,1)) # this works. 
    
    ax.scatter(v[0][:,0],
                v[0][:,1],
                v[0][:,2], s=.1, c='k', alpha=.1)
    ax.scatter(v[0][proj_distances<thresh_proj_distance,0],
                v[0][proj_distances<thresh_proj_distance,1],
                v[0][proj_distances<thresh_proj_distance,2], s=.1, c='r', alpha=.5)
    
    ax.axis('off')
    ax.grid('off')
    # ax.view_init(-45,30)
    ax.view_init(135,210)
    # uzip.set_axes_equal(ax)
    # fig.savefig('curvature_uv_smooth.svg', dpi=300, bbox_inches='tight')
    plt.show()
    
    

    
# =============================================================================
# =============================================================================
# #  We can now use the projected distances to exclude. 
# =============================================================================
# =============================================================================
    proj_distances = np.nansum((v[0]-centroid[None,:]) * v_axis[:,0][None,:], axis=-1)
    if flip_eig:
        exclude_vertices = proj_distances <= thresh_proj_distance
    else:
        exclude_vertices = -1*proj_distances <= thresh_proj_distance


# =============================================================================
# =============================================================================
# # Global plot with raw intensity and  
# =============================================================================
# =============================================================================

    curvature_bins = np.linspace(-2, 2, 51)

    # average over the time !.
    mean_I_time = []
    std_I_time = []
    
    mean_surf_signal = []
    
    
    for ttt in np.arange(len(curvature_faces_all_um_mesh))[:]:
        
        mean_I = []
        std_I = []
        
        kappa_ttt = curvature_faces_all_um_mesh[ttt].copy()
        intensity_ttt = septin_faces_all[ttt].copy()
        # mean_surf_signal.append(np.nanmean(intensity_ttt))
        kappa_ttt = kappa_ttt[~exclude_vertices].copy()
        intensity_ttt = intensity_ttt[~exclude_vertices].copy()
        mean_surf_signal.append(np.nanmean(intensity_ttt))
        
        # intensity_ttt = intensity_ttt / np.nanmean(intensity_ttt)
        for bb_ii in np.arange(len(curvature_bins)-1):
            select = np.logical_and(kappa_ttt > curvature_bins[bb_ii],
                                    kappa_ttt < curvature_bins[bb_ii+1])
            mean_I.append(np.nanmean(intensity_ttt[select]))
            std_I.append(spstats.sem(intensity_ttt[select]))
        mean_I = np.hstack(mean_I)
        std_I = np.hstack(std_I) 
        
        mean_I_time.append(mean_I)
        std_I_time.append(std_I)
        
    mean_I_time = np.array(mean_I_time)
    std_I_time = np.array(std_I_time)
    
    
# =============================================================================
#     Global septin intensity vs curvature needs to be changed. 
# =============================================================================
    
    plt.figure(figsize=(5,5))
    plt.title('Global Septin Intensity vs Intracellular Curvature')
    plt.plot(.5*(curvature_bins[1:] + curvature_bins[:-1]), 
             np.nanmean(mean_I_time, axis=0))
    # plt.fill_between(.5*(curvature_bins[1:] + curvature_bins[:-1]), 
    #          np.nanmean(mean_I_time, axis=0)-spstats.sem(mean_I_time, axis=0, nan_policy='omit'), 
    #          np.nanmean(mean_I_time, axis=0)+spstats.sem(mean_I_time, axis=0, nan_policy='omit'), color='gray', alpha=.5)
    plt.fill_between(.5*(curvature_bins[1:] + curvature_bins[:-1]), 
             np.nanmean(mean_I_time, axis=0)-np.nanstd(mean_I_time, axis=0), 
             np.nanmean(mean_I_time, axis=0)+np.nanstd(mean_I_time, axis=0), color='gray', alpha=.5)
    plt.hlines(np.mean(mean_surf_signal), curvature_bins[0], curvature_bins[-1], color='r', linestyles='dashed', label='mean_septin_I')
    plt.vlines(0.4, 1, 1.15, color='k', linestyles='dashed', label='septin_kappa_thresh')
    # plt.vlines(.9, 1,1.15)
    plt.ylabel('Norm. Septin Intensity', fontsize=18, fontname='Arial')
    plt.xlabel('Intracellular Kappa [1/um]', fontsize=18, fontname='Arial')
    plt.xticks(fontsize=16, fontname='Arial')
    plt.yticks(fontsize=16, fontname='Arial')
    plt.tick_params(right=True, length=5)
    plt.legend(loc='best', fontsize=12)
    # analysis_output_folder
    plt.savefig(os.path.join(analysis_output_folder, 
                             'global_septin_vs_intracellular_curvature_Raw.svg'), dpi=300, bbox_inches='tight')
    plt.show()
    
    
    
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# # # #     Checking the smoothed versions of the signals? 
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================

    
    mean_I_time = []
    std_I_time = []
    
    mean_surf_signal = []
    
    for ttt in np.arange(len(curvature_faces_all_um_mesh))[:]:
        
        mean_I = []
        std_I = []
        
        kappa_ttt = curvature_faces_all_smooth[ttt].copy()
        intensity_ttt = septin_faces_all_smooth[ttt].copy()
        # mean_surf_signal.append(np.nanmean(intensity_ttt))
        
        # mean_surf_signal.append(np.nanmean(intensity_ttt))
        kappa_ttt = kappa_ttt[~exclude_vertices].copy()
        intensity_ttt = intensity_ttt[~exclude_vertices].copy()
        mean_surf_signal.append(np.nanmean(intensity_ttt))
        
        # intensity_ttt = intensity_ttt / np.nanmean(intensity_ttt)
        for bb_ii in np.arange(len(curvature_bins)-1):
            select = np.logical_and(kappa_ttt > curvature_bins[bb_ii],
                                    kappa_ttt < curvature_bins[bb_ii+1])
            mean_I.append(np.nanmean(intensity_ttt[select]))
            std_I.append(spstats.sem(intensity_ttt[select]))
        mean_I = np.hstack(mean_I)
        std_I = np.hstack(std_I) 
        
        mean_I_time.append(mean_I)
        std_I_time.append(std_I)
        
    mean_I_time = np.array(mean_I_time)
    std_I_time = np.array(std_I_time)
    
    
    plt.figure(figsize=(5,5))
    plt.title('Global Septin Intensity vs Intracellular Curvature - Smooth')
    plt.plot(.5*(curvature_bins[1:] + curvature_bins[:-1]), 
             np.nanmean(mean_I_time, axis=0))
    # plt.fill_between(.5*(curvature_bins[1:] + curvature_bins[:-1]), 
    #          np.nanmean(mean_I_time, axis=0)-spstats.sem(mean_I_time, axis=0, nan_policy='omit'), 
    #          np.nanmean(mean_I_time, axis=0)+spstats.sem(mean_I_time, axis=0, nan_policy='omit'), color='gray', alpha=.5)
    plt.fill_between(.5*(curvature_bins[1:] + curvature_bins[:-1]), 
             np.nanmean(mean_I_time, axis=0)-np.nanstd(mean_I_time, axis=0), 
             np.nanmean(mean_I_time, axis=0)+np.nanstd(mean_I_time, axis=0), color='gray', alpha=.5)
    # plt.hlines(np.mean(mean_surf_signal), curvature_bins[0], curvature_bins[-1], color='r')
    # plt.vlines(0.4, 1,1.15)
    plt.hlines(np.mean(mean_surf_signal), curvature_bins[0], curvature_bins[-1], color='r', linestyles='dashed',label='mean_septin_I')
    plt.vlines(0.4, 1, 1.15, color='k', linestyles='dashed', label='septin_kappa_thresh')
    # plt.vlines(.9, 1,1.15)
    plt.ylabel('Norm. Septin Intensity', fontsize=18, fontname='Arial')
    plt.xlabel('Intracellular Kappa [1/um]', fontsize=18, fontname='Arial')
    plt.xticks(fontsize=16, fontname='Arial')
    plt.yticks(fontsize=16, fontname='Arial')
    plt.tick_params(right=True, length=5)
    plt.legend(loc='best', fontsize=12)
    # analysis_output_folder
    plt.savefig(os.path.join(analysis_output_folder, 
                             'global_septin_vs_intracellular_curvature_bothSmooth.svg'), dpi=300, bbox_inches='tight')
    plt.show()
    
    
# =============================================================================
#     Characterise separately figure D (the distribution of intracellular negative curvature patches! )
# =============================================================================

    kappa_bins = np.linspace(0, 1.5 , 51)
    
    mean_I_time = []
    mean_I_time_raw = []
    # std_I_time = []
    
    for ttt in np.arange(len(curvature_faces_all_um_mesh))[:]:
        
        mean_I = []
        std_I = []
        
        kappa_ttt = curvature_faces_all_smooth[ttt].copy()
        intensity_ttt = septin_faces_all_smooth[ttt].copy()
        
        # mean_surf_signal.append(np.nanmean(intensity_ttt))
        kappa_ttt = kappa_ttt[~exclude_vertices].copy()
        intensity_ttt = intensity_ttt[~exclude_vertices].copy()
        # mean_surf_signal.append(np.nanmean(intensity_ttt))
                
        # intensity_ttt = intensity_ttt / np.nanmean(intensity_ttt)
        for bb_ii in np.arange(len(kappa_bins)-1):
            select = np.logical_and(kappa_ttt > kappa_bins[bb_ii],
                                    kappa_ttt < kappa_bins[bb_ii+1])
            mean_I.append(np.nansum(select))
            # std_I.append(spstats.sem(intensity_ttt[select]))
        mean_I = np.hstack(mean_I)
        # std_I = np.hstack(std_I) 
        mean_I_time.append(mean_I/float(np.sum(mean_I)))
        # std_I_time.append(std_I)
        
        mean_I_raw = []
        kappa_ttt_raw = curvature_faces_all_um_mesh[ttt].copy()
        for bb_ii in np.arange(len(kappa_bins)-1):
            select = np.logical_and(kappa_ttt_raw > kappa_bins[bb_ii],
                                    kappa_ttt_raw < kappa_bins[bb_ii+1])
            mean_I_raw.append(np.nansum(select))
        mean_I_time_raw.append(np.hstack(mean_I_raw) / float(np.sum(mean_I_raw)))
        

    mean_I_time = np.array(mean_I_time)
    mean_I_time_raw = np.array(mean_I_time_raw)
    # std_I_time = np.array(std_I_time)
    
    
    mean_curve = mean_I_time.mean(axis=0) / float(np.sum(mean_I_time.mean(axis=0)))
    mean_curve_raw = mean_I_time_raw.mean(axis=0) / float(np.sum(mean_I_time_raw.mean(axis=0)))
    
    
    plt.figure(figsize=(5,5))
    plt.plot(.5*(kappa_bins[1:]+kappa_bins[:-1]), 
             mean_curve, label='smooth')
    plt.plot(.5*(kappa_bins[1:]+kappa_bins[:-1]), 
             mean_curve_raw, label='raw')
    plt.vlines(0.4, 0.01, 0.03, color='k', linestyles='dashed')
    # plt.vlines(0.9, 0.01, 0.03, color='k', linestyles='dashed')
    plt.ylabel('Norm. Septin Intensity', fontsize=18, fontname='Arial')
    plt.xlabel('Intracellular Kappa [1/um]', fontsize=18, fontname='Arial')
    plt.ylim([0, 0.15])
    plt.xlim([0, 1.5])
    plt.xticks(fontsize=16, fontname='Arial')
    plt.yticks(fontsize=16, fontname='Arial')
    plt.tick_params(right=True, length=5)
    plt.legend(loc='best', fontsize=12)
    # analysis_output_folder
    plt.savefig(os.path.join(analysis_output_folder, 
                             'smooth_vs_raw_intracellular_curvature_distribution.svg'), dpi=300, bbox_inches='tight')
    plt.show()
    
    

# =============================================================================
#     Load the computed geodesics and segmentation from the windowed analysis!.  
# =============================================================================
    
    statsfolder = 'E:\\Work\\Projects\\Danuser-3D Causality\\Data\\Andrew\\Septin_Blebs2\\180522\\Cell5\\2022-03-09_Registration\\2022-04-05_Analysis\\Window_Analyses'
    infile = os.path.join(statsfolder, 
                          'save_curvature-septin_statistics.mat')
    stats_obj = spio.loadmat(infile)
    
    all_kappa_groups_faces = stats_obj['all_Kappa_groups'].copy()
    all_geodist_faces_to_bleb = stats_obj['all_geodist_faces_to_bleb'].copy()
    all_geodist_faces_to_bleb = all_geodist_faces_to_bleb*.104

    all_d_septin_hotspots_binary_pos = stats_obj['all_d_septin_hotspots_binary_pos'].copy()
    all_histogram_geodist_d_septin_hotspot = stats_obj['all_histogram_geodist_d_septin_hotspot'].copy()
    all_bins_geodist_d_septin_hotspot = stats_obj['all_bins_geodist_d_septin_hotspot'].copy()
    all_histogram_kappa_d_septin_hotspot = stats_obj['all_histogram_kappa_d_septin_hotspot'].copy()
    all_bins_kappa_d_septin_hotspot = stats_obj['all_bins_kappa_d_septin_hotspot'].copy()
    
    
    all_d_septin_hotspots_binary_neg = stats_obj['all_d_septin_hotspots_binary_neg'].copy()
    all_histogram_geodist_d_septin_hotspot_low = stats_obj['all_histogram_geodist_d_septin_hotspot_low'].copy()
    all_bins_geodist_d_septin_hotspot_low = stats_obj['all_bins_geodist_d_septin_hotspot_low'].copy()
    all_histogram_kappa_d_septin_hotspot_low = stats_obj['all_histogram_kappa_d_septin_hotspot_low'].copy()
    all_bins_kappa_d_septin_hotspot_low = stats_obj['all_bins_kappa_d_septin_hotspot_low'].copy()
    
    # the whole septin region. 
    all_septin_hotspots_binary = stats_obj['all_septin_hotspots_binary'].copy()
    all_histogram_geodist_septin_hotspot = stats_obj['all_histogram_geodist_septin_hotspot'].copy()
    all_bins_geodist_septin_hotspot = stats_obj['all_bins_geodist_septin_hotspot'].copy()
    all_histogram_kappa_septin_hotspot = stats_obj['all_histogram_kappa_septin_hotspot'].copy()
    all_bins_kappa_septin_hotspot = stats_obj['all_bins_kappa_septin_hotspot'].copy()
    
    
    print(stats_obj.keys())
# =============================================================================
#     Pooled analysis over connected regions (and needing to use the segmentation data) ---> check within the patch the distribution. 
# =============================================================================
    
    # connected compo
    clusters_dseptin_pos = np.squeeze(stats_obj['all_dseptin_clusters_pos']).copy()
    clusters_dseptin_neg = np.squeeze(stats_obj['all_dseptin_clusters_neg']).copy()
    
    # define the analysis time windows based on the moving average!
    # max_absolute_curvature_faces[winsize//2:len(v)-winsize//2] /.104
    winsize = 25 # this was inferred....
    start_time = winsize//2
    end_time = len(v) - winsize//2
    
    # iterate over the time. 
    all_pos_cluster_stats_time = []; size_pos_cluster = []
    all_neg_cluster_stats_time = []; size_neg_cluster = []
    min_size=10 # need to put a threshold on the cluster size..... for reliability..... (still.... )
    
    
    #### rederive the raw signals (the above is wrong - there is offset timepoint!)
    septin_faces_all_raw = igl.average_onto_faces(f, septin_faces_all.T).T
    kappa_all_raw = igl.average_onto_faces(f, curvature_faces_all_um_mesh.T).T
    dseptin_faces_all_raw = np.gradient(septin_faces_all_raw, axis=0)
    
    

# =============================================================================
#     Check now the stability of the faces - identifying these faces as units for analysis and clustering?  There is overlap.... 
# =============================================================================
    count_faces_pos = np.sum(all_d_septin_hotspots_binary_pos, axis=0)
    count_faces_neg = np.sum(all_d_septin_hotspots_binary_neg, axis=0)
    
    plt.figure()
    plt.plot(np.sort(count_faces_pos), label='pos')
    plt.plot(np.sort(count_faces_neg), label='neg')
    plt.show()
    
    thresh_stable = skfilters.threshold_otsu(count_faces_pos)
    overlap = np.sum(np.logical_and(count_faces_pos>=thresh_stable, 
                              count_faces_neg>=thresh_stable))
    overlap = overlap/float(np.sum(count_faces_pos>=thresh_stable) + np.sum(count_faces_neg>=thresh_stable))
    
    fig = plt.figure(figsize=plt.figaspect(1.)*3)
    ax = plt.axes(projection='3d', proj_type = 'ortho')
    # ax.set_box_aspect(aspect = (1,1,1))
    ax.scatter(v[0][...,0],
                v[0][...,1],
                v[0][...,2], s=.1)
    
    uniq_verts_pos = np.unique(f[count_faces_pos>=thresh_stable])
    ax.scatter(v[0][uniq_verts_pos,0],
                v[0][uniq_verts_pos,1],
                v[0][uniq_verts_pos,2], s=.1, c='r')
    ax.view_init(-45,0)
    plt.show()
    
    
    fig = plt.figure(figsize=plt.figaspect(1.)*3)
    ax = plt.axes(projection='3d', proj_type = 'ortho')
    # ax.set_box_aspect(aspect = (1,1,1))
    ax.scatter(v[0][...,0],
                v[0][...,1],
                v[0][...,2], s=.1)
    
    uniq_verts_pos = np.unique(f[count_faces_neg>=thresh_stable])
    ax.scatter(v[0][uniq_verts_pos,0],
                v[0][uniq_verts_pos,1],
                v[0][uniq_verts_pos,2], s=.1, c='g')
    ax.view_init(-45,0)
    plt.show()
    
    
    pos_face_select = count_faces_pos>=thresh_stable
    neg_face_select = count_faces_neg>=thresh_stable
    
    """
    how about we cross-overlap? 
    """
    both_face_select = np.logical_and(pos_face_select, neg_face_select)
    
    fig = plt.figure(figsize=plt.figaspect(1.)*3)
    ax = plt.axes(projection='3d', proj_type = 'ortho')
    # ax.set_box_aspect(aspect = (1,1,1))
    ax.scatter(v[0][...,0],
                v[0][...,1],
                v[0][...,2], s=.1)
    
    uniq_verts_pos = np.unique(f[both_face_select])
    ax.scatter(v[0][uniq_verts_pos,0],
                v[0][uniq_verts_pos,1],
                v[0][uniq_verts_pos,2], s=.1, c='r')
    ax.view_init(90,-10)
    plt.show()
    
    
    """
    heatmapping the binarised gating.... 
    """
    fig, ax = plt.subplots(figsize=(5,20))
    plt.title('pos increasing')
    ax.imshow(all_d_septin_hotspots_binary_pos[:,pos_face_select>0].T, cmap='coolwarm')
    ax.set_aspect('auto')
    plt.show()
    
    fig, ax = plt.subplots(figsize=(5,20))
    plt.title('neg decreasing')
    ax.imshow(all_d_septin_hotspots_binary_pos[:,neg_face_select>0].T, cmap='coolwarm')
    ax.set_aspect('auto')
    plt.show()
    
     
    """
    further separate the positive face select and negative face select by thresholds... pooling the 'pos increasing' and 'neg decrasing' faces.
    """
    kappa_all_raw_reduce = kappa_all_raw[np.arange(start_time, end_time)].copy()
    septin_faces_all_raw_reduce = septin_faces_all_raw[np.arange(start_time, end_time)].copy()
       
    joint_select = np.logical_or(pos_face_select>0, neg_face_select>0) # total .... all faces. that are interesting. 
    # categorise by the mean positive / negative faces... ( the thresholds were derived from the smoothed ... )
    joint_face_select_kappa = all_kappa_groups_faces[:,joint_select>0].copy() # get the face classification for all the positive increasing faces. 
    # joint_all_geodist_faces_to_bleb = all_geodist_faces_to_bleb[:,joint_select>0].copy()
    
    
    #### for these we are going to parse all the intervals of blebbing and negative curvature 
    def parse_all_contigs(kappa_class_matrix, kappa_class, min_len=0): #and select for minimal length.
        
        N = len(kappa_class_matrix)
        intervals = []
        for ii in np.arange(N):
            data = kappa_class_matrix[ii].copy()
            valid = np.arange(len(data))[data == kappa_class] # 
            
            # break into contigs
            contigs = []
            contig = []
            for jjj in np.arange(len(valid)):
                if jjj == 0:
                    contig.append(valid[jjj])
                else:
                    if len(contig)>0:
                        diff = valid[jjj] - contig[-1]
                        if diff == 1:
                            contig.append(valid[jjj])
                        else:
                            contigs.append(contig)
                            contig=[valid[jjj]] # start a new one. 
                    else:
                        contig.append(valid[jjj]) # extend cuyrrent. 
                        
            if len(contig) > 0 :
                contigs.append(contig)
                contig = []
                
            # check the contigs. 
            contigs = [cc for cc in contigs if len(cc)>=min_len]
            intervals.append(contigs)
            
        return intervals
    
    """
    filter and retain the minimal length contigs. 
    """
    joint_face_select_kappa_poskappa = parse_all_contigs(joint_face_select_kappa[:,:].T, kappa_class=0, min_len=25)
    joint_face_select_kappa_flatkappa = parse_all_contigs(joint_face_select_kappa[:,:].T, kappa_class=1, min_len=25)
    joint_face_select_kappa_negkappa = parse_all_contigs(joint_face_select_kappa[:,:].T, kappa_class=2, min_len=25)
    
    
    """
    convert this to indices.
    """
    output_save_mesh_folder = os.path.join(analysis_output_folder, 'contig_vis_folder')
    fio.mkdir(output_save_mesh_folder)

    index_face_timeseries = np.arange(len(mesh_f))[joint_select>0] # get which faces should be highlighted..... # transfer these to face colors. 

    # annotate this. 
    mesh_timeseries_contigs = -1*np.ones((septin_faces_all_raw_reduce.shape[0], septin_faces_all_raw_reduce.shape[1]), dtype=np.int) # initialise the  
    for iii in np.arange(len(joint_face_select_kappa_poskappa)): # iterate over the selected faces. 
        tsa_index = index_face_timeseries[iii]
        for select in joint_face_select_kappa_poskappa[iii]:
            for ss in select:
                mesh_timeseries_contigs[ss,tsa_index] = 2 # mark these aas all positive kappa. 
    for iii in np.arange(len(joint_face_select_kappa_flatkappa)): # iterate over the selected faces. 
        tsa_index = index_face_timeseries[iii]
        for select in joint_face_select_kappa_flatkappa[iii]:
            for ss in select:
                mesh_timeseries_contigs[ss,tsa_index] = 1 # mark these aas all positive kappa. 
    for iii in np.arange(len(joint_face_select_kappa_negkappa)): # iterate over the selected faces. 
        tsa_index = index_face_timeseries[iii]
        for select in joint_face_select_kappa_negkappa[iii]:
            for ss in select:
                mesh_timeseries_contigs[ss,tsa_index] = 0 # mark these aas all positive kappa. 
        

    # transfer this labelling into colors and export out the different timepoints into meshes. 
    for tt in np.arange(len(mesh_timeseries_contigs))[:]:
        # start_time, end_time
        mesh_label_tt = mesh_timeseries_contigs[tt].copy()
        mesh_label_color = vol_colors.get_colors(mesh_label_tt, colormap=cm.Spectral_r, vmin=0, vmax=2)
        mesh_label_color[mesh_label_tt==-1] = np.hstack([0.75,0.75,0.75, 1])[None,:]

        mesh_tt = trimesh.Trimesh(vertices=mesh_v[start_time+tt],
                                  faces=mesh_f,
                                  face_colors = np.uint8(255*mesh_label_color[:,:3]), 
                                  process=False,
                                  validate=False)
        # save out. 
        export_mesh_file = os.path.join(output_save_mesh_folder, 'contigs_'+str(start_time+tt).zfill(3)+'.obj')
        mesh_tt.export(export_mesh_file) # save this out. 
    


#     dseptin_pos_contigs= parse_all_contigs((all_d_septin_hotspots_binary_pos[:,:]>0).T, kappa_class=1, min_len=0)
#     dseptin_pos_contigs_flat = [item for sublist in dseptin_pos_contigs for item in sublist if len(item)>0]
    def get_timeseries(tsa, matrix_mask_list):
        
        out_series = []
        for ii in np.arange(len(tsa)):
            data_out_series = []
            data = tsa[ii].copy()
            mask = matrix_mask_list[ii]
            if len(mask)>0:
                for mm in mask:
                    data_out_series.append(data[mm])
            out_series.append(data_out_series)
        return out_series
            
    # then we are going to run the correlation over all the valid!. 
    # the xcorr is general..... xcorr_timeseries_set_1d(timeseries_array1, timeseries_array2, norm=True, eps=1e-12)
    
    """
    pull down relevant time series with the binary masks
    """
    # switch to absolute curvature? or just reverse.... ?
    joint_face_select_kappa_poskappa_kappa_timeseries = get_timeseries(-kappa_all_raw_reduce[:,joint_select>0].T, joint_face_select_kappa_poskappa)
    joint_face_select_kappa_poskappa_sept_timeseries = get_timeseries(septin_faces_all_raw_reduce[:,joint_select>0].T, joint_face_select_kappa_poskappa)
    joint_face_select_kappa_poskappa_bleb_distance_timeseries = get_timeseries(all_geodist_faces_to_bleb[:,joint_select>0].T, joint_face_select_kappa_poskappa)
    joint_face_select_kappa_poskappa_pos_septin_timeseries = get_timeseries(all_d_septin_hotspots_binary_pos[:,joint_select>0].T, joint_face_select_kappa_poskappa)
    joint_face_select_kappa_poskappa_neg_septin_timeseries = get_timeseries(all_d_septin_hotspots_binary_neg[:,joint_select>0].T, joint_face_select_kappa_poskappa)
    # for flat too and add in a timeseries for geo dist. 
    
    joint_face_select_kappa_flatkappa_kappa_timeseries = get_timeseries(-kappa_all_raw_reduce[:,joint_select>0].T, joint_face_select_kappa_flatkappa)
    joint_face_select_kappa_flatkappa_sept_timeseries = get_timeseries(septin_faces_all_raw_reduce[:,joint_select>0].T, joint_face_select_kappa_flatkappa)
    joint_face_select_kappa_flatkappa_bleb_distance_timeseries = get_timeseries(all_geodist_faces_to_bleb[:,joint_select>0].T, joint_face_select_kappa_flatkappa)
    joint_face_select_kappa_flatkappa_pos_septin_timeseries = get_timeseries(all_d_septin_hotspots_binary_pos[:,joint_select>0].T, joint_face_select_kappa_flatkappa)
    joint_face_select_kappa_flatkappa_neg_septin_timeseries = get_timeseries(all_d_septin_hotspots_binary_neg[:,joint_select>0].T, joint_face_select_kappa_flatkappa)
    
    # flip this one to represent increasing I 
    joint_face_select_kappa_negkappa_kappa_timeseries = get_timeseries(-1*-kappa_all_raw_reduce[:,joint_select>0].T, joint_face_select_kappa_negkappa)
    joint_face_select_kappa_negkappa_sept_timeseries = get_timeseries(septin_faces_all_raw_reduce[:,joint_select>0].T, joint_face_select_kappa_negkappa)
    joint_face_select_kappa_negkappa_bleb_distance_timeseries = get_timeseries(all_geodist_faces_to_bleb[:,joint_select>0].T, joint_face_select_kappa_negkappa)
    joint_face_select_kappa_negkappa_pos_septin_timeseries = get_timeseries(all_d_septin_hotspots_binary_pos[:,joint_select>0].T, joint_face_select_kappa_negkappa)
    joint_face_select_kappa_negkappa_neg_septin_timeseries = get_timeseries(all_d_septin_hotspots_binary_neg[:,joint_select>0].T, joint_face_select_kappa_negkappa)
    
    """
    flatten the timeseries
    """ 
    joint_face_select_kappa_poskappa_kappa_timeseries_flat = [item for sublist in joint_face_select_kappa_poskappa_kappa_timeseries for item in sublist if len(item)>0]
    joint_face_select_kappa_poskappa_sept_timeseries_flat = [item for sublist in joint_face_select_kappa_poskappa_sept_timeseries for item in sublist if len(item)>0]
    joint_face_select_kappa_poskappa_bleb_distance_timeseries = [item for sublist in joint_face_select_kappa_poskappa_bleb_distance_timeseries for item in sublist if len(item)>0]
    joint_face_select_kappa_poskappa_pos_septin_timeseries_flat = [item for sublist in joint_face_select_kappa_poskappa_pos_septin_timeseries for item in sublist if len(item)>0]
    joint_face_select_kappa_poskappa_neg_septin_timeseries_flat = [item for sublist in joint_face_select_kappa_poskappa_neg_septin_timeseries for item in sublist if len(item)>0]
    
    joint_face_select_kappa_flatkappa_kappa_timeseries_flat = [item for sublist in joint_face_select_kappa_flatkappa_kappa_timeseries for item in sublist if len(item)>0]
    joint_face_select_kappa_flatkappa_sept_timeseries_flat = [item for sublist in joint_face_select_kappa_flatkappa_sept_timeseries for item in sublist if len(item)>0]
    joint_face_select_kappa_flatkappa_bleb_distance_timeseries = [item for sublist in joint_face_select_kappa_flatkappa_bleb_distance_timeseries for item in sublist if len(item)>0]
    joint_face_select_kappa_flatkappa_pos_septin_timeseries_flat = [item for sublist in joint_face_select_kappa_flatkappa_pos_septin_timeseries for item in sublist if len(item)>0]
    joint_face_select_kappa_flatkappa_neg_septin_timeseries_flat = [item for sublist in joint_face_select_kappa_flatkappa_neg_septin_timeseries for item in sublist if len(item)>0]
    
    joint_face_select_kappa_negkappa_kappa_timeseries_flat = [item for sublist in joint_face_select_kappa_negkappa_kappa_timeseries for item in sublist if len(item)>0]
    joint_face_select_kappa_negkappa_sept_timeseries_flat = [item for sublist in joint_face_select_kappa_negkappa_sept_timeseries for item in sublist if len(item)>0]
    joint_face_select_kappa_negkappa_bleb_distance_timeseries = [item for sublist in joint_face_select_kappa_negkappa_bleb_distance_timeseries for item in sublist if len(item)>0]
    joint_face_select_kappa_negkappa_pos_septin_timeseries_flat = [item for sublist in joint_face_select_kappa_negkappa_pos_septin_timeseries for item in sublist if len(item)>0]
    joint_face_select_kappa_negkappa_neg_septin_timeseries_flat = [item for sublist in joint_face_select_kappa_negkappa_neg_septin_timeseries for item in sublist if len(item)>0]
    
    
    
    """
    tabulate statistics to assay for the dynamicity measure. 
    """
    # count the number of positive increasing septin and the number of decreasing septin? # or just jointly tally? 
    joint_face_select_kappa_poskappa_pos_septin_timeseries_flat_counts = np.hstack([np.sum(cc) for cc in joint_face_select_kappa_poskappa_pos_septin_timeseries_flat])
    joint_face_select_kappa_flatkappa_pos_septin_timeseries_flat_counts = np.hstack([np.sum(cc) for cc in joint_face_select_kappa_flatkappa_pos_septin_timeseries_flat])
    joint_face_select_kappa_negkappa_pos_septin_timeseries_flat_counts = np.hstack([np.sum(cc) for cc in joint_face_select_kappa_negkappa_pos_septin_timeseries_flat])
    
    # count the number of positive increasing septin and the number of decreasing septin? # or just jointly tally? 
    joint_face_select_kappa_poskappa_neg_septin_timeseries_flat_counts = np.hstack([np.sum(cc) for cc in joint_face_select_kappa_poskappa_neg_septin_timeseries_flat])
    joint_face_select_kappa_flatkappa_neg_septin_timeseries_flat_counts = np.hstack([np.sum(cc) for cc in joint_face_select_kappa_flatkappa_neg_septin_timeseries_flat])
    joint_face_select_kappa_negkappa_neg_septin_timeseries_flat_counts = np.hstack([np.sum(cc) for cc in joint_face_select_kappa_negkappa_neg_septin_timeseries_flat])
    
    joint_poskappa_lens = np.hstack([len(cc) for cc in joint_face_select_kappa_poskappa_pos_septin_timeseries_flat])
    joint_flatkappa_lens = np.hstack([len(cc) for cc in joint_face_select_kappa_flatkappa_pos_septin_timeseries_flat])
    joint_negkappa_lens = np.hstack([len(cc) for cc in joint_face_select_kappa_negkappa_neg_septin_timeseries_flat])
    
    # is there a mistake here? hm...... 
    joint_fluc_ratio_poskappa = (joint_face_select_kappa_poskappa_pos_septin_timeseries_flat_counts + joint_face_select_kappa_poskappa_neg_septin_timeseries_flat_counts) / joint_poskappa_lens.astype(np.float)
    joint_fluc_ratio_flatkappa = (joint_face_select_kappa_flatkappa_pos_septin_timeseries_flat_counts + joint_face_select_kappa_flatkappa_neg_septin_timeseries_flat_counts) / joint_flatkappa_lens.astype(np.float)
    joint_fluc_ratio_negkappa = (joint_face_select_kappa_negkappa_pos_septin_timeseries_flat_counts + joint_face_select_kappa_negkappa_neg_septin_timeseries_flat_counts) / joint_negkappa_lens.astype(np.float)
    
    plt.figure()
    plt.plot(np.sort(joint_fluc_ratio_poskappa), label='fluc ratio pos kappa')
    plt.plot(np.sort(joint_fluc_ratio_flatkappa), label='fluc ratio flat kappa')
    plt.plot(np.sort(joint_fluc_ratio_negkappa), label='fluc ratio neg kappa')
    plt.legend()
    plt.show()
    
    
# =============================================================================
# =============================================================================
# #     Xcorr timeseries analysis 
# =============================================================================
# =============================================================================
    
    """
    Compute the xcorr for all!. - not only separating by correlation but also the value change in the septin fluctuations!. i.e. variance!
    """
    
    """
    do the plotting with the proper lags + add the 95% confidence interval !. 
    +ve kappa !!!!! 
    """
    # note the reversal for positive faces. 
    joint_face_pos_kappa_corr = tsa.xcorr_timeseries_set_1d(joint_face_select_kappa_poskappa_kappa_timeseries_flat, 
                                                            joint_face_select_kappa_poskappa_sept_timeseries_flat, 
                                                            norm=True, eps=1e-12,stack_final=False)
    ### check the lags... !. 
    
    # accumulate the lags and correlations from each timeseries!. 
    joint_face_pos_kappa_corr_arr = tsa.stack_xcorr_curves(joint_face_pos_kappa_corr)
    
    fluc_thresh_pos = skfilters.threshold_multiotsu(joint_fluc_ratio_poskappa,3)[0]
    
    # this only makes sense if Ts=1, therefore we should now use the proper Ts. which is not 1! 
    x_corr_plot = np.linspace(-(joint_face_pos_kappa_corr_arr.shape[1]//2),
                              joint_face_pos_kappa_corr_arr.shape[1]//2, joint_face_pos_kappa_corr_arr.shape[1]) * Ts # multiply by the sampling time!. 
    
    plt.figure(figsize=(5,5))
    plt.title('joint faces on positive kappa')
    plt.hlines(0, -25, 25, color='k', linestyles='dashed')
    plt.plot(x_corr_plot, np.nanmean(joint_face_pos_kappa_corr_arr, axis=0), color='k', label='All'); 
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_pos_kappa_corr_arr, axis=0) - 1.96*spstats.sem(joint_face_pos_kappa_corr_arr, axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_pos_kappa_corr_arr, axis=0) + 1.96*spstats.sem(joint_face_pos_kappa_corr_arr, axis=0, nan_policy='omit'), color='k', alpha=.5)
    plt.vlines(0, -0.6, 0.6, color='k', linestyles='dashed')
    plt.plot(x_corr_plot, np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0), color='r', label='Changing'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='r')
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0) - 1.96*spstats.sem(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0) + 1.96*spstats.sem(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0, nan_policy='omit'), color='r', alpha=.5)
    plt.plot(x_corr_plot, np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0), color='g', label='Less Changing'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='b')
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0) - 1.96*spstats.sem(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0) + 1.96*spstats.sem(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0, nan_policy='omit'), color='g', alpha=.5)
    # plt.xlim([joint_face_pos_kappa_corr_arr.shape[1]//2-25, joint_face_pos_kappa_corr_arr.shape[1]//2+25])
    plt.xlim(-25,25)
    plt.ylim([-0.6,0.6])
    plt.tick_params(length=5, right=True)
    plt.legend(loc='best')
    plt.xticks(fontsize=16, fontname='Arial')
    plt.yticks(fontsize=16, fontname='Arial')
    plt.xlabel('Lag (s)', fontsize=18, fontname='Arial')
    plt.ylabel('xcorr(curvature, septin)', fontsize=18, fontname='Arial')
    plt.savefig(os.path.join(analysis_output_folder, 
                              'xcorr_raw_septinI_vs_curvature_on_stable_deltaFaces_Pos-contigs_win25limit.svg'), dpi=300, bbox_inches='tight')
    plt.show()
    
    
    

    """
    do the plotting with the proper lags + add the 95% confidence interval !. 
    flat kappa !!!!! 
    """
    joint_face_flat_kappa_corr = tsa.xcorr_timeseries_set_1d(joint_face_select_kappa_flatkappa_kappa_timeseries_flat, 
                                                            joint_face_select_kappa_flatkappa_sept_timeseries_flat, 
                                                            norm=True, eps=1e-12,stack_final=False)
    # accumulate the lags and correlations from each timeseries!. 
    joint_face_flat_kappa_corr_arr = tsa.stack_xcorr_curves(joint_face_flat_kappa_corr)
    fluc_thresh_flat = skfilters.threshold_multiotsu(joint_fluc_ratio_flatkappa,3)[0]
    
    x_corr_plot_flat = np.linspace(-(joint_face_flat_kappa_corr_arr.shape[1]//2),
                                   joint_face_flat_kappa_corr_arr.shape[1]//2, joint_face_flat_kappa_corr_arr.shape[1]) * Ts


    plt.figure(figsize=(5,5))
    plt.title('joint faces on flat kappa')
    plt.hlines(0, -25, 25, color='k', linestyles='dashed')
    plt.plot(x_corr_plot_flat, np.nanmean(joint_face_flat_kappa_corr_arr, axis=0), color='k', label='All'); 
    plt.fill_between(x_corr_plot_flat,
                      np.nanmean(joint_face_flat_kappa_corr_arr, axis=0) - 1.96*spstats.sem(joint_face_flat_kappa_corr_arr, axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_flat_kappa_corr_arr, axis=0) + 1.96*spstats.sem(joint_face_flat_kappa_corr_arr, axis=0, nan_policy='omit'), color='k', alpha=.5)
    plt.vlines(0, -0.6, 0.6, color='k', linestyles='dashed')
    plt.plot(x_corr_plot_flat, np.nanmean(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa>=fluc_thresh_flat], axis=0), color='r', label='Changing'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='r')
    plt.fill_between(x_corr_plot_flat,
                      np.nanmean(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa>=fluc_thresh_flat], axis=0) - 1.96*spstats.sem(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa>=fluc_thresh_flat], axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa>=fluc_thresh_flat], axis=0) + 1.96*spstats.sem(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa>=fluc_thresh_flat], axis=0, nan_policy='omit'), color='r', alpha=.5)
    plt.plot(x_corr_plot, np.nanmean(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa<fluc_thresh_flat], axis=0), color='g', label='Less Changing'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='b')
    plt.fill_between(x_corr_plot_flat,
                      np.nanmean(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa<fluc_thresh_flat], axis=0) - 1.96*spstats.sem(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa<fluc_thresh_flat], axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa<fluc_thresh_flat], axis=0) + 1.96*spstats.sem(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa<fluc_thresh_flat], axis=0, nan_policy='omit'), color='g', alpha=.5)
    # plt.xlim([joint_face_pos_kappa_corr_arr.shape[1]//2-25, joint_face_pos_kappa_corr_arr.shape[1]//2+25])
    plt.xlim(-25,25)
    plt.ylim([-0.6,0.6])
    plt.tick_params(length=5, right=True)
    plt.legend(loc='best')
    plt.xticks(fontsize=16, fontname='Arial')
    plt.yticks(fontsize=16, fontname='Arial')
    plt.xlabel('Lag (s)', fontsize=18, fontname='Arial')
    plt.ylabel('xcorr(curvature, septin)', fontsize=18, fontname='Arial')
    plt.savefig(os.path.join(analysis_output_folder, 
                              'xcorr_raw_septinI_vs_curvature_on_stable_deltaFaces_Flat-contigs_win25limit.svg'), dpi=300, bbox_inches='tight')
    plt.show()
    
    
    plt.figure(figsize=(5,5))
    plt.title('joint faces on flat kappa')
    plt.hlines(0, -100, 100, color='k', linestyles='dashed')
    plt.plot(x_corr_plot_flat, np.nanmean(joint_face_flat_kappa_corr_arr, axis=0), color='k', label='All'); 
    plt.fill_between(x_corr_plot_flat,
                      np.nanmean(joint_face_flat_kappa_corr_arr, axis=0) - 1.96*spstats.sem(joint_face_flat_kappa_corr_arr, axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_flat_kappa_corr_arr, axis=0) + 1.96*spstats.sem(joint_face_flat_kappa_corr_arr, axis=0, nan_policy='omit'), color='k', alpha=.5)
    plt.vlines(0, -0.6, 0.6, color='k', linestyles='dashed')
    plt.plot(x_corr_plot_flat, np.nanmean(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa>=fluc_thresh_flat], axis=0), color='r', label='Changing'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='r')
    plt.fill_between(x_corr_plot_flat,
                      np.nanmean(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa>=fluc_thresh_flat], axis=0) - 1.96*spstats.sem(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa>=fluc_thresh_flat], axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa>=fluc_thresh_flat], axis=0) + 1.96*spstats.sem(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa>=fluc_thresh_flat], axis=0, nan_policy='omit'), color='r', alpha=.5)
    plt.plot(x_corr_plot, np.nanmean(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa<fluc_thresh_flat], axis=0), color='g', label='Less Changing'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='b')
    plt.fill_between(x_corr_plot_flat,
                      np.nanmean(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa<fluc_thresh_flat], axis=0) - 1.96*spstats.sem(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa<fluc_thresh_flat], axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa<fluc_thresh_flat], axis=0) + 1.96*spstats.sem(joint_face_flat_kappa_corr_arr[joint_fluc_ratio_flatkappa<fluc_thresh_flat], axis=0, nan_policy='omit'), color='g', alpha=.5)
    # plt.xlim([joint_face_pos_kappa_corr_arr.shape[1]//2-25, joint_face_pos_kappa_corr_arr.shape[1]//2+25])
    plt.xlim(-100, 100)
    plt.ylim([-0.6,0.6])
    plt.tick_params(length=5, right=True)
    plt.legend(loc='best')
    plt.xticks(fontsize=16, fontname='Arial')
    plt.yticks(fontsize=16, fontname='Arial')
    plt.xlabel('Lag (s)', fontsize=18, fontname='Arial')
    plt.ylabel('xcorr(curvature, septin)', fontsize=18, fontname='Arial')
    plt.savefig(os.path.join(analysis_output_folder, 
                              'xcorr_raw_septinI_vs_curvature_on_stable_deltaFaces_Flat-contigs_win100limit.svg'), dpi=300, bbox_inches='tight')
    plt.show()
    

    """
    neg kappa !!!
    """
    joint_face_neg_kappa_corr = tsa.xcorr_timeseries_set_1d(joint_face_select_kappa_negkappa_kappa_timeseries_flat, 
                                                           joint_face_select_kappa_negkappa_sept_timeseries_flat, 
                                                          norm=True, eps=1e-12,stack_final=False)
    # accumulate the lags and correlations from each timeseries!. 
    joint_face_neg_kappa_corr_arr = tsa.stack_xcorr_curves(joint_face_neg_kappa_corr)
    fluc_thresh_neg = skfilters.threshold_multiotsu(joint_fluc_ratio_negkappa,3)[0]
    
    x_corr_plot_neg = np.linspace(-(joint_face_neg_kappa_corr_arr.shape[1]//2),
                                  joint_face_neg_kappa_corr_arr.shape[1]//2, joint_face_neg_kappa_corr_arr.shape[1]) * Ts
    
    plt.figure(figsize=(5,5))
    # plt.title('joint faces on negative kappa')
    # plt.plot(np.nanmedian(pos_face_pos_kappa_corr_arr, axis=0)); plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='k')
    plt.plot(x_corr_plot_neg, np.nanmean(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa>=fluc_thresh_neg], axis=0), color='r',label='Changing_neg_kappa'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='r')
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa>=fluc_thresh_neg], axis=0) - 1.96*spstats.sem(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa>=fluc_thresh_neg], axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa>=fluc_thresh_neg], axis=0) + 1.96*spstats.sem(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa>=fluc_thresh_neg], axis=0, nan_policy='omit'), color='r', alpha=.25)
    
    plt.plot(x_corr_plot_neg, np.nanmean(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa<fluc_thresh_neg], axis=0), color='g', label='Less_Changing_neg_kappa'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='b')
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa<fluc_thresh_neg], axis=0) - 1.96*spstats.sem(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa<fluc_thresh_neg], axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa<fluc_thresh_neg], axis=0) + 1.96*spstats.sem(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa<fluc_thresh_neg], axis=0, nan_policy='omit'), color='g', alpha=.25)
    
    plt.plot(x_corr_plot_neg, np.nanmean(joint_face_neg_kappa_corr_arr, axis=0), color='k', label='All_neg_kappa'); 
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_neg_kappa_corr_arr, axis=0) - 1.96*spstats.sem(joint_face_neg_kappa_corr_arr, axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_neg_kappa_corr_arr, axis=0) + 1.96*spstats.sem(joint_face_neg_kappa_corr_arr, axis=0, nan_policy='omit'), color='k', alpha=.25)
    
    # plt.vlines(joint_face_neg_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='k')
    plt.vlines(0, -0.6,.6, color='k', linestyles='dashed')
    plt.hlines(0, -25, 25, color='k', linestyles='dashed')
    
    plt.plot(x_corr_plot, np.nanmean(joint_face_pos_kappa_corr_arr, axis=0), 'k--', label='All_pos_kappa'); 
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_pos_kappa_corr_arr, axis=0) - 1.96*spstats.sem(joint_face_pos_kappa_corr_arr, axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_pos_kappa_corr_arr, axis=0) + 1.96*spstats.sem(joint_face_pos_kappa_corr_arr, axis=0, nan_policy='omit'), color='k', alpha=.25)
    
    # plt.plot(x_corr_plot_neg, np.nanmedian(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0), 'r--', label='Changing_pos_kappa')
    plt.plot(x_corr_plot, np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0), 'r--', label='Changing_pos_kappa'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='r')
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0) - 1.96*spstats.sem(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0) + 1.96*spstats.sem(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0, nan_policy='omit'), color='r', alpha=.25)
    
    plt.plot(x_corr_plot, np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0), color='g', label='Less Changing_pos_kappa'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='b')
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0) - 1.96*spstats.sem(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0) + 1.96*spstats.sem(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0, nan_policy='omit'), color='g', alpha=.5)
    
    plt.xlim(-25,25)
    plt.ylim([-0.6,0.6])
    plt.tick_params(length=5, right=True)
    plt.legend(loc='best')
    plt.xticks(fontsize=16, fontname='Arial')
    plt.yticks(fontsize=16, fontname='Arial')
    plt.xlabel('Lag (s)', fontsize=18, fontname='Arial')
    plt.ylabel('xcorr(curvature, septin)', fontsize=18, fontname='Arial')
    plt.savefig(os.path.join(analysis_output_folder, 
                              'xcorr_raw_septinI_vs_curvature_on_stable_deltaFaces_contigs_win25limit.svg'), dpi=300, bbox_inches='tight')
    plt.show()



    plt.figure(figsize=(5,5))
    # plt.title('joint faces on negative kappa')
    # plt.plot(np.nanmedian(pos_face_pos_kappa_corr_arr, axis=0)); plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='k')
    plt.plot(x_corr_plot_neg, np.nanmean(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa>=fluc_thresh_neg], axis=0), color='r',label='Changing_neg_kappa'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='r')
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa>=fluc_thresh_neg], axis=0) - 1.96*spstats.sem(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa>=fluc_thresh_neg], axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa>=fluc_thresh_neg], axis=0) + 1.96*spstats.sem(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa>=fluc_thresh_neg], axis=0, nan_policy='omit'), color='r', alpha=.25)
    
    plt.plot(x_corr_plot_neg, np.nanmean(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa<fluc_thresh_neg], axis=0), color='g', label='Less_Changing_neg_kappa'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='b')
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa<fluc_thresh_neg], axis=0) - 1.96*spstats.sem(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa<fluc_thresh_neg], axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa<fluc_thresh_neg], axis=0) + 1.96*spstats.sem(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa<fluc_thresh_neg], axis=0, nan_policy='omit'), color='g', alpha=.25)
    
    plt.plot(x_corr_plot_neg, np.nanmean(joint_face_neg_kappa_corr_arr, axis=0), color='k', label='All_neg_kappa'); 
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_neg_kappa_corr_arr, axis=0) - 1.96*spstats.sem(joint_face_neg_kappa_corr_arr, axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_neg_kappa_corr_arr, axis=0) + 1.96*spstats.sem(joint_face_neg_kappa_corr_arr, axis=0, nan_policy='omit'), color='k', alpha=.25)
    
    # plt.vlines(joint_face_neg_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='k')
    plt.vlines(0, -0.6,.6, color='k', linestyles='dashed')
    plt.hlines(0, -100, 100, color='k', linestyles='dashed')
    
    plt.plot(x_corr_plot, np.nanmean(joint_face_pos_kappa_corr_arr, axis=0), 'k--', label='All_pos_kappa'); 
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_pos_kappa_corr_arr, axis=0) - 1.96*spstats.sem(joint_face_pos_kappa_corr_arr, axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_pos_kappa_corr_arr, axis=0) + 1.96*spstats.sem(joint_face_pos_kappa_corr_arr, axis=0, nan_policy='omit'), color='k', alpha=.25)
    
    # plt.plot(x_corr_plot_neg, np.nanmedian(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0), 'r--', label='Changing_pos_kappa')
    plt.plot(x_corr_plot, np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0), 'r--', label='Changing_pos_kappa'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='r')
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0) - 1.96*spstats.sem(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0) + 1.96*spstats.sem(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0, nan_policy='omit'), color='r', alpha=.25)
    
    plt.plot(x_corr_plot, np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0), color='g', label='Less Changing_pos_kappa'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='b')
    plt.fill_between(x_corr_plot,
                      np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0) - 1.96*spstats.sem(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0, nan_policy='omit'), 
                      np.nanmean(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0) + 1.96*spstats.sem(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa<fluc_thresh_pos], axis=0, nan_policy='omit'), color='g', alpha=.5)
    
    plt.xlim(-100,100)
    plt.ylim([-0.6,0.6])
    plt.tick_params(length=5, right=True)
    plt.legend(loc='best')
    plt.xticks(fontsize=16, fontname='Arial')
    plt.yticks(fontsize=16, fontname='Arial')
    plt.xlabel('Lag (s)', fontsize=18, fontname='Arial')
    plt.ylabel('xcorr(curvature, septin)', fontsize=18, fontname='Arial')
    plt.savefig(os.path.join(analysis_output_folder, 
                              'xcorr_raw_septinI_vs_curvature_on_stable_deltaFaces_contigs_win100limit.svg'), dpi=300, bbox_inches='tight')
    plt.show()
    
    
    
    """
    get the oscillation time.
    """
    from scipy.signal import find_peaks
    
    peaks_pos = find_peaks(-np.nanmean(joint_face_pos_kappa_corr_arr, axis=0), distance=5, prominence=0.005)
    peaks_neg = find_peaks(np.nanmean(joint_face_neg_kappa_corr_arr, axis=0), distance=5, prominence=0.005)
    
    
    plt.figure(figsize=(5,5))
    plt.plot(x_corr_plot, np.nanmean(joint_face_pos_kappa_corr_arr, axis=0))
    plt.plot(x_corr_plot[peaks_pos[0]], np.nanmean(joint_face_pos_kappa_corr_arr, axis=0)[peaks_pos[0]], 'ro')
    plt.plot(x_corr_plot, np.nanmean(joint_face_neg_kappa_corr_arr, axis=0))
    plt.plot(x_corr_plot[peaks_neg[0]], np.nanmean(joint_face_neg_kappa_corr_arr, axis=0)[peaks_neg[0]], 'go')
    plt.show()
    # 31.9 s for negative using the mean diff. 
    # 27.6 s for positive... 
    
    
    # plt.figure()
    # plt.title('joint neg kappa')
    # # plt.plot(np.nanmedian(pos_face_pos_kappa_corr_arr, axis=0)); plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='k')
    # plt.plot(np.nanmedian(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa>=fluc_thresh_neg], axis=0), color='r'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='r')
    # plt.plot(np.nanmedian(joint_face_neg_kappa_corr_arr[joint_fluc_ratio_negkappa<fluc_thresh_neg], axis=0), color='g'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='b')
    # plt.plot(np.nanmedian(joint_face_neg_kappa_corr_arr, axis=0)); plt.vlines(joint_face_neg_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='k')
    # plt.plot(np.nanmedian(joint_face_pos_kappa_corr_arr[joint_fluc_ratio_poskappa>=fluc_thresh_pos], axis=0), 'r--')
    # # plt.xlim([joint_face_neg_kappa_corr_arr.shape[1]//2-25, joint_face_neg_kappa_corr_arr.shape[1]//2+25])
    # plt.xlim([joint_face_neg_kappa_corr_arr.shape[1]//2-100, joint_face_neg_kappa_corr_arr.shape[1]//2+100])
    # plt.show()
    
    
#     """
#     its hard to call positive/negative ratio....
#     """
#     # we can look at sum relative to total length? 
    
#     """
#     plot of fluc vs mean septin intensity
#     """
#     reg_pos = spstats.linregress([np.nanmean(ss) for ss in joint_face_select_kappa_poskappa_sept_timeseries_flat], 
#                                   joint_fluc_ratio_poskappa)
#     reg_neg = spstats.linregress([np.nanmean(ss) for ss in joint_face_select_kappa_negkappa_sept_timeseries_flat],
#                               joint_fluc_ratio_negkappa)
    
#     xline = np.linspace(1,1.7,10)
    
    
#     plt.figure()
#     plt.plot([np.nanmean(ss) for ss in joint_face_select_kappa_poskappa_sept_timeseries_flat], 
#               joint_fluc_ratio_poskappa, 'r.', alpha=0.1, label='posKappa')
#     plt.plot(xline, reg_pos[0]*xline + reg_pos[1], 'k')
#     plt.legend()
#     plt.figure()
#     plt.plot([np.nanmean(ss) for ss in joint_face_select_kappa_negkappa_sept_timeseries_flat], 
#               joint_fluc_ratio_negkappa, 'g.', alpha=0.1, label='negKappa')
#     plt.plot(xline, reg_neg[0]*xline + reg_neg[1], 'k')
#     plt.legend()
#     plt.show()
    
    
#     import seaborn as sns
    
#     datapos = pd.DataFrame(np.vstack([np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_poskappa_sept_timeseries_flat]), 
#                                       joint_fluc_ratio_poskappa]).T,
#                             columns =['x','y'])
    
#     g = sns.jointplot(
#     data=datapos,
#     x="x", y="y", 
#     kind="kde",
#     )
    dataneg = pd.DataFrame(np.vstack([np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_flatkappa_sept_timeseries_flat]), 
                                      #joint_fluc_ratio_negkappa]).T,
                                      np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_flatkappa_bleb_distance_timeseries])]).T,
                            columns =['x','y'])
    
    g = sns.jointplot(
    data=dataneg,
    x="x", y="y", 
    kind="kde",
    )
    dataneg = pd.DataFrame(np.vstack([np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_poskappa_sept_timeseries_flat]), 
                                      #joint_fluc_ratio_negkappa]).T,
                                      np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_poskappa_bleb_distance_timeseries])]).T,
                            columns =['x','y'])
    
    g = sns.jointplot(
    data=dataneg,
    x="x", y="y", 
    kind="kde",
    )
    dataneg = pd.DataFrame(np.vstack([np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_negkappa_sept_timeseries_flat]), 
                                      #joint_fluc_ratio_negkappa]).T,
                                      np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_negkappa_bleb_distance_timeseries])]).T,
                            columns =['x','y'])
    
    g = sns.jointplot(
    data=dataneg,
    x="x", y="y", 
    kind="kde",
    )
    
    """
    Combine the density.... 
    """
    all_mean_I_contigs = np.hstack([ np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_negkappa_sept_timeseries_flat]), 
                                     np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_poskappa_sept_timeseries_flat]),
                                     np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_flatkappa_sept_timeseries_flat]) ])
    all_mean_bleb_dist_contigs = np.hstack([ np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_negkappa_bleb_distance_timeseries]), 
                                            np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_poskappa_bleb_distance_timeseries]),
                                            np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_flatkappa_bleb_distance_timeseries]) ])
    all_mean_kappa_class = np.hstack([ np.hstack([0 for ss in joint_face_select_kappa_negkappa_bleb_distance_timeseries]), 
                                       np.hstack([2 for ss in joint_face_select_kappa_poskappa_bleb_distance_timeseries]),
                                       np.hstack([1 for ss in joint_face_select_kappa_flatkappa_bleb_distance_timeseries]) ])
    
    dataneg = pd.DataFrame(np.vstack([np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_negkappa_sept_timeseries_flat]), 
                                      #joint_fluc_ratio_negkappa]).T,
                                      np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_negkappa_bleb_distance_timeseries])]).T,
                            columns =['x','y'])
    
    dataneg = pd.DataFrame(np.vstack([all_mean_I_contigs, 
                                      all_mean_bleb_dist_contigs,
                                      all_mean_kappa_class]).T,
                            columns =['septin_I','dist_to_bleb','kappa_cls'])
    g = sns.jointplot(
    data=dataneg,
    x="septin_I", y="dist_to_bleb", hue='kappa_cls', 
    kind="kde",
    )
    
    
# =============================================================================
#     Do the marginal analysis to get the intensity vs geodesic  distance to bleb for the contig !. 
# =============================================================================

    y_bleb_distance = np.linspace(0, 2.5, 100) # this will be the upper bound.
    x_septin_I = np.linspace(1., 1.5, 100)
    # the above is more accurate with 2D kde. 
    # Peform the kernel density estimate
    xx, yy = np.meshgrid(x_septin_I, y_bleb_distance) #
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_negkappa_sept_timeseries_flat]), 
                        np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_negkappa_bleb_distance_timeseries])])
    kernel = spstats.gaussian_kde(values) # test different bandwith. 
    f = np.reshape(kernel(positions).T, xx.shape) # curvature x, septin y. 
    
    values_posKappa = np.vstack([np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_poskappa_sept_timeseries_flat]), 
                        np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_poskappa_bleb_distance_timeseries])])
    kernel_posKappa = spstats.gaussian_kde(values_posKappa) # test different bandwith. 
    f_posKappa = np.reshape(kernel_posKappa(positions).T, xx.shape) # curvature x, septin y. 

    values_flatKappa = np.vstack([np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_flatkappa_sept_timeseries_flat]), 
                        np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_flatkappa_bleb_distance_timeseries])])
    kernel_flatKappa = spstats.gaussian_kde(values_flatKappa) # test different bandwith. 
    f_flatKappa = np.reshape(kernel_flatKappa(positions).T, xx.shape) # curvature x, septin y. 
    
    """
    integrate to get the marginal!. ( in this case )
    """
    # get the marginal and check this is correct!. 
    marginal_I_blebs_negkappa = np.nansum(f*y_bleb_distance[:,None], axis=0) / np.nansum(f, axis=0).astype(np.float)
    var_marginal_I_blebs_negkappa = np.nansum(f*(y_bleb_distance[:,None]**2), axis=0) / np.nansum(f, axis=0).astype(np.float) - marginal_I_blebs_negkappa**2
    sem_marginal_I_blebs_negkappa = np.sqrt(var_marginal_I_blebs_negkappa) / (np.sqrt(np.nansum(f, axis=0).astype(np.float)+1e-15))
    
    marginal_I_blebs_poskappa = np.nansum(f_posKappa*y_bleb_distance[:,None], axis=0) / np.nansum(f_posKappa, axis=0).astype(np.float)
    var_marginal_I_blebs_poskappa = np.nansum(f_posKappa*(y_bleb_distance[:,None]**2), axis=0) / np.nansum(f_posKappa, axis=0).astype(np.float) - marginal_I_blebs_poskappa**2
    sem_marginal_I_blebs_poskappa = np.sqrt(var_marginal_I_blebs_poskappa) / (np.sqrt(np.nansum(f_posKappa, axis=0).astype(np.float)+1e-15))
    
    marginal_I_blebs_flatkappa = np.nansum(f_flatKappa*y_bleb_distance[:,None], axis=0) / np.nansum(f_flatKappa, axis=0).astype(np.float)
    var_marginal_I_blebs_flatkappa = np.nansum(f_flatKappa*(y_bleb_distance[:,None]**2), axis=0) / np.nansum(f_flatKappa, axis=0).astype(np.float) - marginal_I_blebs_flatkappa**2
    sem_marginal_I_blebs_flatkappa = np.sqrt(var_marginal_I_blebs_flatkappa) / (np.sqrt(np.nansum(f_flatKappa, axis=0).astype(np.float)+1e-15))
    
    
    plt.figure(figsize=(5,5))
    # plt.title('for neg kappa')
    plt.plot(x_septin_I, marginal_I_blebs_negkappa)
    plt.fill_between(x_septin_I, 
                     marginal_I_blebs_negkappa-sem_marginal_I_blebs_negkappa, 
                     marginal_I_blebs_negkappa+sem_marginal_I_blebs_negkappa, alpha=0.5, label='negKappa')
    plt.plot(x_septin_I, marginal_I_blebs_poskappa)
    plt.fill_between(x_septin_I, 
                     marginal_I_blebs_poskappa-sem_marginal_I_blebs_poskappa, 
                     marginal_I_blebs_poskappa+sem_marginal_I_blebs_poskappa, alpha=0.5, label='posKappa')
    plt.plot(x_septin_I, marginal_I_blebs_flatkappa)
    plt.fill_between(x_septin_I, 
                     marginal_I_blebs_flatkappa-sem_marginal_I_blebs_flatkappa, 
                     marginal_I_blebs_flatkappa+sem_marginal_I_blebs_flatkappa, alpha=0.5, label='flatKappa')
    plt.tick_params(length=5, right=True)
    # plt.ylim([0,0.1])
    # plt.xlim([0,1])
    plt.xlabel('Norm. Septin I')
    plt.ylabel('Blebby distance')
    plt.legend(loc='best')
    # plt.savefig(os.path.join(analysis_output_folder, 
    #                           'integrateKDE_Dynamicity_vs_dI_on_stable_deltaFaces_contigs_flatKappa.svg'), dpi=300, bbox_inches='tight')
    plt.show()
    
    
    """
    Plot of mean kappa and intensity on each contig ... separated by low and high dynamicity. ! giving 2 lines.  
    """
    
    # each of these have different thresholds ...  fluc_thresh_pos, fluc_thresh_neg, fluc_thresh_flat
    # intensity 
    mean_contig_joint_face_select_kappa_poskappa_sept_timeseries_flat = np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_poskappa_sept_timeseries_flat])
    mean_contig_joint_face_select_kappa_flatkappa_sept_timeseries_flat = np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_flatkappa_sept_timeseries_flat])
    mean_contig_joint_face_select_kappa_negkappa_sept_timeseries_flat = np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_negkappa_sept_timeseries_flat])
    # kappa 
    mean_contig_joint_face_select_kappa_poskappa_kappa_timeseries_flat = np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_poskappa_kappa_timeseries_flat])
    mean_contig_joint_face_select_kappa_flatkappa_kappa_timeseries_flat = np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_flatkappa_kappa_timeseries_flat])
    mean_contig_joint_face_select_kappa_negkappa_kappa_timeseries_flat = np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_negkappa_kappa_timeseries_flat])
    # assign by fluc ratio threshold? 
    
    select_contig_fluc_poskappa = joint_fluc_ratio_poskappa>=fluc_thresh_pos
    select_contig_fluc_flatkappa = joint_fluc_ratio_flatkappa>=fluc_thresh_flat
    select_contig_fluc_negkappa = joint_fluc_ratio_negkappa>=fluc_thresh_neg
    
    
    """
    create 2 new series. 
    """
    # mean_contig_high_dynamic_intensity = np.hstack([mean_contig_joint_face_select_kappa_poskappa_sept_timeseries_flat[select_contig_fluc_poskappa>0], 
    #                                                 mean_contig_joint_face_select_kappa_flatkappa_sept_timeseries_flat[select_contig_fluc_flatkappa>0],
    #                                                 mean_contig_joint_face_select_kappa_negkappa_sept_timeseries_flat[select_contig_fluc_negkappa>0]])
    # mean_contig_low_dynamic_intensity = np.hstack([mean_contig_joint_face_select_kappa_poskappa_sept_timeseries_flat[select_contig_fluc_poskappa==0], 
    #                                                mean_contig_joint_face_select_kappa_flatkappa_sept_timeseries_flat[select_contig_fluc_flatkappa==0],
    #                                                mean_contig_joint_face_select_kappa_negkappa_sept_timeseries_flat[select_contig_fluc_negkappa==0]])
    # mean_contig_high_dynamic_kappa = np.hstack([mean_contig_joint_face_select_kappa_poskappa_kappa_timeseries_flat[select_contig_fluc_poskappa>0], 
    #                                            mean_contig_joint_face_select_kappa_flatkappa_kappa_timeseries_flat[select_contig_fluc_flatkappa>0],
    #                                            -1*mean_contig_joint_face_select_kappa_negkappa_kappa_timeseries_flat[select_contig_fluc_negkappa>0]])
    # mean_contig_low_dynamic_kappa = np.hstack([mean_contig_joint_face_select_kappa_poskappa_kappa_timeseries_flat[select_contig_fluc_poskappa==0], 
    #                                            mean_contig_joint_face_select_kappa_flatkappa_kappa_timeseries_flat[select_contig_fluc_flatkappa==0],
    #                                            -1*mean_contig_joint_face_select_kappa_negkappa_kappa_timeseries_flat[select_contig_fluc_negkappa==0]])
    
    mean_contig_high_dynamic_intensity = mean_contig_joint_face_select_kappa_negkappa_sept_timeseries_flat[select_contig_fluc_negkappa>0]
    mean_contig_low_dynamic_intensity = mean_contig_joint_face_select_kappa_negkappa_sept_timeseries_flat[select_contig_fluc_negkappa==0]
    mean_contig_high_dynamic_kappa = mean_contig_joint_face_select_kappa_negkappa_kappa_timeseries_flat[select_contig_fluc_negkappa>0]
    mean_contig_low_dynamic_kappa = mean_contig_joint_face_select_kappa_negkappa_kappa_timeseries_flat[select_contig_fluc_negkappa==0]
    
    
    ##### instead of mean which lowers the score. --- looks like we may very wel need to pool..... 
    contig_high_dynamic_intensity = np.hstack([joint_face_select_kappa_negkappa_sept_timeseries_flat[jjj] for jjj in np.arange(len(joint_face_select_kappa_negkappa_sept_timeseries_flat)) if select_contig_fluc_negkappa[jjj]>0])
    contig_high_dynamic_kappa = np.hstack([joint_face_select_kappa_negkappa_kappa_timeseries_flat[jjj] for jjj in np.arange(len(joint_face_select_kappa_negkappa_kappa_timeseries_flat)) if select_contig_fluc_negkappa[jjj]>0])
    contig_low_dynamic_intensity = np.hstack([joint_face_select_kappa_negkappa_sept_timeseries_flat[jjj] for jjj in np.arange(len(joint_face_select_kappa_negkappa_sept_timeseries_flat)) if select_contig_fluc_negkappa[jjj]==0])
    contig_low_dynamic_kappa = np.hstack([joint_face_select_kappa_negkappa_kappa_timeseries_flat[jjj] for jjj in np.arange(len(joint_face_select_kappa_negkappa_kappa_timeseries_flat)) if select_contig_fluc_negkappa[jjj]==0])
    # just use the different kappa instead. 

    """
    check this result on the flat surface. 
    """
    flatkappa_contig_high_dynamic_intensity = np.hstack([joint_face_select_kappa_flatkappa_sept_timeseries_flat[jjj] for jjj in np.arange(len(joint_face_select_kappa_flatkappa_sept_timeseries_flat)) if select_contig_fluc_flatkappa[jjj]>0])
    flatkappa_contig_high_dynamic_kappa = np.hstack([joint_face_select_kappa_flatkappa_kappa_timeseries_flat[jjj] for jjj in np.arange(len(joint_face_select_kappa_flatkappa_kappa_timeseries_flat)) if select_contig_fluc_flatkappa[jjj]>0])
    flatkappa_contig_low_dynamic_intensity = np.hstack([joint_face_select_kappa_flatkappa_sept_timeseries_flat[jjj] for jjj in np.arange(len(joint_face_select_kappa_flatkappa_sept_timeseries_flat)) if select_contig_fluc_flatkappa[jjj]==0])
    flatkappa_contig_low_dynamic_kappa = np.hstack([joint_face_select_kappa_flatkappa_kappa_timeseries_flat[jjj] for jjj in np.arange(len(joint_face_select_kappa_flatkappa_kappa_timeseries_flat)) if select_contig_fluc_flatkappa[jjj]==0])
    
    """
    check this result on positive surface.
    """
    poskappa_contig_high_dynamic_intensity = np.hstack([joint_face_select_kappa_poskappa_sept_timeseries_flat[jjj] for jjj in np.arange(len(joint_face_select_kappa_poskappa_sept_timeseries_flat)) if select_contig_fluc_poskappa[jjj]>0])
    poskappa_contig_high_dynamic_kappa = np.hstack([joint_face_select_kappa_poskappa_kappa_timeseries_flat[jjj] for jjj in np.arange(len(joint_face_select_kappa_poskappa_kappa_timeseries_flat)) if select_contig_fluc_poskappa[jjj]>0])
    poskappa_contig_low_dynamic_intensity = np.hstack([joint_face_select_kappa_poskappa_sept_timeseries_flat[jjj] for jjj in np.arange(len(joint_face_select_kappa_poskappa_sept_timeseries_flat)) if select_contig_fluc_poskappa[jjj]==0])
    poskappa_contig_low_dynamic_kappa = np.hstack([joint_face_select_kappa_poskappa_kappa_timeseries_flat[jjj] for jjj in np.arange(len(joint_face_select_kappa_poskappa_kappa_timeseries_flat)) if select_contig_fluc_poskappa[jjj]==0])
    

    """
    using the mean contig statistics. 
    """
    # the problem is the range of the septin. 
    # y_septin = np.linspace(1.20, 1.35, 100) # this will be the upper bound.
    y_septin = np.linspace(.8, 2, 100)
    # x_kappa = np.linspace(-1.2, 1.2, 100)
    x_kappa = np.linspace(0, 1.2, 100)
    # the above is more accurate with 2D kde. 
    # Peform the kernel density estimate
    xx, yy = np.meshgrid(x_kappa, y_septin) #
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([mean_contig_high_dynamic_kappa.ravel(), 
                        mean_contig_high_dynamic_intensity.ravel()])
    kernel = spstats.gaussian_kde(values) # test different bandwith. 
    f = np.reshape(kernel(positions).T, xx.shape) # curvature x, septin y. 
    
    values_low = np.vstack([mean_contig_low_dynamic_kappa.ravel(), 
                            mean_contig_low_dynamic_intensity.ravel()])
    kernel_low = spstats.gaussian_kde(values_low) # test different bandwith. 
    f_low = np.reshape(kernel_low(positions).T, xx.shape) # curvature x, septin y. 
    
    """
    plot the marginal.
    """
    # get the marginal and check this is correct!. 
    marginal_pos_septin_high_dynamicity = np.nansum(f*y_septin[:,None], axis=0) / np.nansum(f, axis=0).astype(np.float)
    var_marginal_pos_septin_high_dynamicity = np.nansum(f*(y_septin[:,None]**2), axis=0) / np.nansum(f, axis=0).astype(np.float) - marginal_pos_septin_high_dynamicity**2
    sem_marginal_pos_septin_high_dynamicity = np.sqrt(var_marginal_pos_septin_high_dynamicity) / (np.sqrt(np.nansum(f, axis=0).astype(np.float)+1e-15))
    
    marginal_neg_septin_low_dynamicity = np.nansum(f_low*y_septin[:,None], axis=0) / np.nansum(f_low, axis=0).astype(np.float)
    var_marginal_neg_septin_low_dynamicity = np.nansum(f_low*(y_septin[:,None]**2), axis=0) / np.nansum(f_low, axis=0).astype(np.float) - marginal_neg_septin_low_dynamicity**2
    sem_marginal_neg_septin_low_dynamicity = np.sqrt(var_marginal_neg_septin_low_dynamicity) / (np.sqrt(np.nansum(f_low, axis=0).astype(np.float)+1e-15))
    
    plt.figure(figsize=(5,5))
    plt.title('kappa vs intensity')
    plt.plot(x_kappa, marginal_pos_septin_high_dynamicity)
    plt.fill_between(x_kappa, 
                     marginal_pos_septin_high_dynamicity-sem_marginal_pos_septin_high_dynamicity, 
                     marginal_pos_septin_high_dynamicity+sem_marginal_pos_septin_high_dynamicity, alpha=0.5, label='high_dynamicity')
    plt.plot(x_kappa, marginal_neg_septin_low_dynamicity)
    plt.fill_between(x_kappa, 
                     marginal_neg_septin_low_dynamicity-sem_marginal_neg_septin_low_dynamicity, 
                     marginal_neg_septin_low_dynamicity+sem_marginal_neg_septin_low_dynamicity, alpha=0.5,label='low_dynamicity')
    plt.tick_params(length=5, right=True)
    plt.vlines(fluc_thresh_neg, 0, 0.1, linestyles='dashed', color='k')
    plt.ylim([1.0,1.6])
    # plt.ylim([0,0.1])
    # plt.xlim([0,1])
    # plt.xlabel('Dynamicity [0-1]')
    # plt.ylabel('mean abs. Delta I')
    plt.legend(loc='best')
    # plt.savefig(os.path.join(analysis_output_folder, 
                              # 'integrateKDE_Dynamicity_vs_dI_on_stable_deltaFaces_contigs_negKappa.svg'), dpi=300, bbox_inches='tight')
    plt.show()
    
    

    """
    using the compiled contig statistics instead on the negative kappa. 
    """
    # y_septin = np.linspace(1.20, 1.35, 100) # this will be the upper bound.
    y_septin = np.linspace(0.80, 2., 100)
    # x_kappa = np.linspace(-1.2, 1.2, 100)
    x_kappa = np.linspace(0, 1.2, 100)
    # the above is more accurate with 2D kde. 
    # Peform the kernel density estimate
    xx, yy = np.meshgrid(x_kappa, y_septin) #
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([contig_high_dynamic_kappa.ravel(), 
                        contig_high_dynamic_intensity.ravel()])
    kernel = spstats.gaussian_kde(values) # test different bandwith. 
    f = np.reshape(kernel(positions).T, xx.shape) # curvature x, septin y. 
    
    values_low = np.vstack([contig_low_dynamic_kappa.ravel(), 
                            contig_low_dynamic_intensity.ravel()])
    kernel_low = spstats.gaussian_kde(values_low) # test different bandwith. 
    f_low = np.reshape(kernel_low(positions).T, xx.shape) # curvature x, septin y. 
    
    
    """
    plot the marginal.
    """
    # get the marginal and check this is correct!. 
    marginal_pos_septin_high_dynamicity = np.nansum(f*y_septin[:,None], axis=0) / np.nansum(f, axis=0).astype(np.float)
    var_marginal_pos_septin_high_dynamicity = np.nansum(f*(y_septin[:,None]**2), axis=0) / np.nansum(f, axis=0).astype(np.float) - marginal_pos_septin_high_dynamicity**2
    sem_marginal_pos_septin_high_dynamicity = np.sqrt(var_marginal_pos_septin_high_dynamicity) / (np.sqrt(np.nansum(f, axis=0).astype(np.float)+1e-15))
    
    marginal_neg_septin_low_dynamicity = np.nansum(f_low*y_septin[:,None], axis=0) / np.nansum(f_low, axis=0).astype(np.float)
    var_marginal_neg_septin_low_dynamicity = np.nansum(f_low*(y_septin[:,None]**2), axis=0) / np.nansum(f_low, axis=0).astype(np.float) - marginal_neg_septin_low_dynamicity**2
    sem_marginal_neg_septin_low_dynamicity = np.sqrt(var_marginal_neg_septin_low_dynamicity) / (np.sqrt(np.nansum(f_low, axis=0).astype(np.float)+1e-15))
    
    
    plt.figure(figsize=(5,5))
    plt.title('kappa vs intensity negKappa')
    plt.plot(x_kappa, marginal_pos_septin_high_dynamicity)
    plt.fill_between(x_kappa, 
                     marginal_pos_septin_high_dynamicity-sem_marginal_pos_septin_high_dynamicity, 
                     marginal_pos_septin_high_dynamicity+sem_marginal_pos_septin_high_dynamicity, alpha=0.5, label='high_dynamicity')
    plt.plot(x_kappa, marginal_neg_septin_low_dynamicity)
    plt.fill_between(x_kappa, 
                     marginal_neg_septin_low_dynamicity-sem_marginal_neg_septin_low_dynamicity, 
                     marginal_neg_septin_low_dynamicity+sem_marginal_neg_septin_low_dynamicity, alpha=0.5,label='low_dynamicity')
    plt.tick_params(length=5, right=True)
    plt.vlines(fluc_thresh_neg, 0, 0.1, linestyles='dashed', color='k')
    plt.ylim([1.15,1.4])
    # plt.ylim([1.0,1.5])
    plt.xlim([0,1.2])
    # plt.ylim([1.0,1.6])
    # plt.ylim([1.2,1.35])
    # plt.ylim([0,0.1])
    # plt.xlim([0,1])
    # plt.xlabel('Dynamicity [0-1]')
    # plt.ylabel('mean abs. Delta I')
    plt.legend(loc='best')
    plt.xticks(fontsize=14, fontname='Arial')
    plt.yticks(fontsize=14, fontname='Arial')
    plt.xlabel('Mean Curvature [um]', fontsize=16, fontname='Arial')
    plt.ylabel('Norm. Septin Intensity', fontsize=16, fontname='Arial')
    # plt.savefig(os.path.join(analysis_output_folder, 
                              # 'integrateKDE_Dynamicity_vs_dI_on_stable_deltaFaces_contigs_negKappa.svg'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(analysis_output_folder, 
                               'kappa_vs_intensity_contigs_negKappa_Zoom.svg'), dpi=300, bbox_inches='tight')
    plt.show()
    

    """
    using the compiled contig statistics instead on the flat kappa. 
    """
    y_septin = np.linspace(0.80, 2., 100)
    # x_kappa = np.linspace(-1.2, 1.2, 100)
    x_kappa = np.linspace(0, 1.2, 100)
    # the above is more accurate with 2D kde. 
    # Peform the kernel density estimate
    xx, yy = np.meshgrid(x_kappa, y_septin) #
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([flatkappa_contig_high_dynamic_kappa.ravel(), #poskappa_contig_high_dynamic_intensity
                        flatkappa_contig_high_dynamic_intensity.ravel()])
    kernel = spstats.gaussian_kde(values) # test different bandwith. 
    f = np.reshape(kernel(positions).T, xx.shape) # curvature x, septin y. 
    
    values_low = np.vstack([flatkappa_contig_low_dynamic_kappa.ravel(), 
                            flatkappa_contig_low_dynamic_intensity.ravel()])
    kernel_low = spstats.gaussian_kde(values_low) # test different bandwith. 
    f_low = np.reshape(kernel_low(positions).T, xx.shape) # curvature x, septin y. 
    
    """
    plot the marginal.
    """
    # get the marginal and check this is correct!. 
    marginal_pos_septin_high_dynamicity = np.nansum(f*y_septin[:,None], axis=0) / np.nansum(f, axis=0).astype(np.float)
    var_marginal_pos_septin_high_dynamicity = np.nansum(f*(y_septin[:,None]**2), axis=0) / np.nansum(f, axis=0).astype(np.float) - marginal_pos_septin_high_dynamicity**2
    sem_marginal_pos_septin_high_dynamicity = np.sqrt(var_marginal_pos_septin_high_dynamicity) / (np.sqrt(np.nansum(f, axis=0).astype(np.float)+1e-15))
    
    marginal_neg_septin_low_dynamicity = np.nansum(f_low*y_septin[:,None], axis=0) / np.nansum(f_low, axis=0).astype(np.float)
    var_marginal_neg_septin_low_dynamicity = np.nansum(f_low*(y_septin[:,None]**2), axis=0) / np.nansum(f_low, axis=0).astype(np.float) - marginal_neg_septin_low_dynamicity**2
    sem_marginal_neg_septin_low_dynamicity = np.sqrt(var_marginal_neg_septin_low_dynamicity) / (np.sqrt(np.nansum(f_low, axis=0).astype(np.float)+1e-15))
    
    plt.figure(figsize=(5,5))
    plt.title('kappa vs intensity flatKappa')
    plt.plot(x_kappa, marginal_pos_septin_high_dynamicity)
    plt.fill_between(x_kappa, 
                     marginal_pos_septin_high_dynamicity-sem_marginal_pos_septin_high_dynamicity, 
                     marginal_pos_septin_high_dynamicity+sem_marginal_pos_septin_high_dynamicity, alpha=0.5, label='high_dynamicity')
    plt.plot(x_kappa, marginal_neg_septin_low_dynamicity)
    plt.fill_between(x_kappa, 
                     marginal_neg_septin_low_dynamicity-sem_marginal_neg_septin_low_dynamicity, 
                     marginal_neg_septin_low_dynamicity+sem_marginal_neg_septin_low_dynamicity, alpha=0.5,label='low_dynamicity')
    plt.tick_params(length=5, right=True)
    plt.vlines(fluc_thresh_neg, 0, 0.1, linestyles='dashed', color='k')
    # plt.ylim([1.15,1.4])
    plt.ylim([1.0,1.5])
    # plt.ylim([1.2,1.35])
    # plt.ylim([0,0.1])
    plt.xlim([0,1.2])
    # plt.xlabel('Dynamicity [0-1]')
    # plt.ylabel('mean abs. Delta I')
    plt.legend(loc='best')
    plt.xticks(fontsize=14, fontname='Arial')
    plt.yticks(fontsize=14, fontname='Arial')
    plt.xlabel('Mean Curvature [um]', fontsize=16, fontname='Arial')
    plt.ylabel('Norm. Septin Intensity', fontsize=16, fontname='Arial')
    # plt.savefig(os.path.join(analysis_output_folder, 
                              # 'integrateKDE_Dynamicity_vs_dI_on_stable_deltaFaces_contigs_negKappa.svg'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(analysis_output_folder, 
                               'kappa_vs_intensity_contigs_flatKappa.svg'), dpi=300, bbox_inches='tight')
    plt.show()



    """
    using the compiled contig statistics instead on the positive kappa. 
    """
    # y_septin = np.linspace(1.20, 1.35, 100) # this will be the upper bound.
    y_septin = np.linspace(0.80, 2., 100)
    # x_kappa = np.linspace(-1.2, 1.2, 100)
    x_kappa = np.linspace(0, 1.2, 100)
    # the above is more accurate with 2D kde. 
    # Peform the kernel density estimate
    xx, yy = np.meshgrid(x_kappa, y_septin) #
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([poskappa_contig_high_dynamic_kappa.ravel(), #poskappa_contig_high_dynamic_intensity
                        poskappa_contig_high_dynamic_intensity.ravel()])
    kernel = spstats.gaussian_kde(values) # test different bandwith. 
    f = np.reshape(kernel(positions).T, xx.shape) # curvature x, septin y. 
    
    values_low = np.vstack([poskappa_contig_low_dynamic_kappa.ravel(), 
                            poskappa_contig_low_dynamic_intensity.ravel()])
    kernel_low = spstats.gaussian_kde(values_low) # test different bandwith. 
    f_low = np.reshape(kernel_low(positions).T, xx.shape) # curvature x, septin y. 
    
    """
    plot the marginal.
    """
    # get the marginal and check this is correct!. 
    marginal_pos_septin_high_dynamicity = np.nansum(f*y_septin[:,None], axis=0) / np.nansum(f, axis=0).astype(np.float)
    var_marginal_pos_septin_high_dynamicity = np.nansum(f*(y_septin[:,None]**2), axis=0) / np.nansum(f, axis=0).astype(np.float) - marginal_pos_septin_high_dynamicity**2
    sem_marginal_pos_septin_high_dynamicity = np.sqrt(var_marginal_pos_septin_high_dynamicity) / (np.sqrt(np.nansum(f, axis=0).astype(np.float)+1e-15))
    
    marginal_neg_septin_low_dynamicity = np.nansum(f_low*y_septin[:,None], axis=0) / np.nansum(f_low, axis=0).astype(np.float)
    var_marginal_neg_septin_low_dynamicity = np.nansum(f_low*(y_septin[:,None]**2), axis=0) / np.nansum(f_low, axis=0).astype(np.float) - marginal_neg_septin_low_dynamicity**2
    sem_marginal_neg_septin_low_dynamicity = np.sqrt(var_marginal_neg_septin_low_dynamicity) / (np.sqrt(np.nansum(f_low, axis=0).astype(np.float)+1e-15))
    
    plt.figure(figsize=(5,5))
    plt.title('kappa vs intensity posKappa')
    plt.plot(x_kappa, marginal_pos_septin_high_dynamicity)
    plt.fill_between(x_kappa, 
                     marginal_pos_septin_high_dynamicity-sem_marginal_pos_septin_high_dynamicity, 
                     marginal_pos_septin_high_dynamicity+sem_marginal_pos_septin_high_dynamicity, alpha=0.5, label='high_dynamicity')
    plt.plot(x_kappa, marginal_neg_septin_low_dynamicity)
    plt.fill_between(x_kappa, 
                     marginal_neg_septin_low_dynamicity-sem_marginal_neg_septin_low_dynamicity, 
                     marginal_neg_septin_low_dynamicity+sem_marginal_neg_septin_low_dynamicity, alpha=0.5,label='low_dynamicity')
    plt.tick_params(length=5, right=True)
    plt.vlines(fluc_thresh_neg, 0, 0.1, linestyles='dashed', color='k')
    # plt.ylim([1.15,1.4])
    plt.ylim([1.0,1.5])
    plt.xlim([0,1.2])
    # plt.ylim([1.2,1.35])
    # plt.ylim([0,0.1])
    # plt.xlim([0,1])
    # plt.xlabel('Dynamicity [0-1]')
    # plt.ylabel('mean abs. Delta I')
    plt.legend(loc='best')
    plt.xticks(fontsize=14, fontname='Arial')
    plt.yticks(fontsize=14, fontname='Arial')
    plt.xlabel('Mean Curvature [um]', fontsize=16, fontname='Arial')
    plt.ylabel('Norm. Septin Intensity', fontsize=16, fontname='Arial')
    # plt.savefig(os.path.join(analysis_output_folder, 
                              # 'integrateKDE_Dynamicity_vs_dI_on_stable_deltaFaces_contigs_negKappa.svg'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(analysis_output_folder, 
                               'kappa_vs_intensity_contigs_posKappa.svg'), dpi=300, bbox_inches='tight')
    plt.show()




    ###### double check by doing the discrete histogram approach. 
    # y_septin = np.linspace(1.20, 1.35, 100) # this will be the upper bound.
    x_kappa_bins = np.linspace(0, 1.5, 21)
    
    mean_septin_kappa_contig_high_dynamic = []
    mean_septin_kappa_contig_low_dynamic = []
    for bin_ii in np.arange(len(x_kappa_bins)-1):
        # select = np.logical_and(mean_contig_high_dynamic_kappa>=x_kappa_bins[bin_ii], 
        #                         mean_contig_high_dynamic_kappa<x_kappa_bins[bin_ii+1])
        # mean_septin_kappa_contig_high_dynamic.append(np.nanmean(mean_contig_high_dynamic_intensity[select>0]))
        
        # contig_high_dynamic_intensity = np.hstack([joint_face_select_kappa_negkappa_sept_timeseries_flat[jjj] for jjj in np.arange(len(joint_face_select_kappa_negkappa_sept_timeseries_flat)) if select_contig_fluc_negkappa[jjj]>0])
        # contig_high_dynamic_kappa = np.hstack([joint_face_select_kappa_negkappa_kappa_timeseries_flat[jjj] for jjj in np.arange(len(joint_face_select_kappa_negkappa_kappa_timeseries_flat)) if select_contig_fluc_negkappa[jjj]>0])
        select = np.logical_and(contig_high_dynamic_kappa>=x_kappa_bins[bin_ii], 
                                contig_high_dynamic_kappa<x_kappa_bins[bin_ii+1])
        mean_septin_kappa_contig_high_dynamic.append(np.nanmean(contig_high_dynamic_intensity[select>0]))
        # select = np.logical_and(mean_contig_low_dynamic_kappa>=x_kappa_bins[bin_ii], 
        #                         mean_contig_low_dynamic_kappa<x_kappa_bins[bin_ii+1])
        # mean_septin_kappa_contig_low_dynamic.append(np.nanmean(mean_contig_low_dynamic_intensity[select>0]))
        select = np.logical_and(contig_low_dynamic_kappa>=x_kappa_bins[bin_ii], 
                                contig_low_dynamic_kappa<x_kappa_bins[bin_ii+1])
        mean_septin_kappa_contig_low_dynamic.append(np.nanmean(contig_low_dynamic_intensity[select>0]))
    mean_septin_kappa_contig_high_dynamic = np.hstack(mean_septin_kappa_contig_high_dynamic)
    mean_septin_kappa_contig_low_dynamic = np.hstack(mean_septin_kappa_contig_low_dynamic)
    
    plt.figure()
    plt.plot(.5*(x_kappa_bins[1:]+x_kappa_bins[:-1]), 
             mean_septin_kappa_contig_high_dynamic)
    plt.plot(.5*(x_kappa_bins[1:]+x_kappa_bins[:-1]), 
             mean_septin_kappa_contig_low_dynamic)
    plt.show()
    
    
    

    """
    plot of fluc vs mean delta septin intensity where we get the instantaneous dSeptin
    """
    joint_face_select_kappa_poskappa_delta_sept_timeseries_flat = [np.gradient(ss) for ss in joint_face_select_kappa_poskappa_sept_timeseries_flat]
    joint_face_select_kappa_flatkappa_delta_sept_timeseries_flat = [np.gradient(ss) for ss in joint_face_select_kappa_flatkappa_sept_timeseries_flat]
    joint_face_select_kappa_negkappa_delta_sept_timeseries_flat = [np.gradient(ss) for ss in joint_face_select_kappa_negkappa_sept_timeseries_flat]
    
#     mean_joint_face_select_kappa_poskappa_delta_sept_timeseries_flat = np.hstack([np.nanmean(np.gradient(ss)) for ss in joint_face_select_kappa_poskappa_sept_timeseries_flat])
#     mean_joint_face_select_kappa_negkappa_delta_sept_timeseries_flat = np.hstack([np.nanmean(np.gradient(ss)) for ss in joint_face_select_kappa_negkappa_sept_timeseries_flat])
#     # dataposdelta = pd.DataFrame(np.vstack([joint_fluc_ratio_negkappa, 
#     #                                        mean_joint_face_select_kappa_negkappa_delta_sept_timeseries_flat]).T,
#     #                        columns =['dynamicity','dSeptin'])
    
#     # g = sns.jointplot(
#     # data=dataposdelta,
#     # x="dynamicity", y="dSeptin", 
#     # kind="kde",
#     # )

#     # datanegdelta = pd.DataFrame(np.vstack([np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_negkappa_sept_timeseries_flat]),
#     #                                       mean_joint_face_select_kappa_negkappa_delta_sept_timeseries_flat]).T,
#     #                        columns =['Septin','dSeptin'])
    
#     # g = sns.jointplot(
#     # data=datanegdelta,
#     # # x="dynamicity", y="abs(dSeptin)", 
#     # x="Septin", y="dSeptin", 
#     # kind="kde",
#     # )
    
    """
    Split into positive and negative delta and separately assay this . 
    """
    mean_pos_joint_face_select_kappa_negkappa_delta_sept_timeseries_flat = np.hstack([np.nanmean(ss[ss>0]) for ss in joint_face_select_kappa_negkappa_delta_sept_timeseries_flat])
    mean_neg_joint_face_select_kappa_negkappa_delta_sept_timeseries_flat = np.hstack([np.nanmean(ss[ss<0]) for ss in joint_face_select_kappa_negkappa_delta_sept_timeseries_flat])
    
    
    y_dseptin = np.linspace(0,0.1,100) # this will be the upper bound.
    x_dynamicity = np.linspace(0, 1., 100)
    # the above is more accurate with 2D kde. 
    # Peform the kernel density estimate
    xx, yy = np.meshgrid(x_dynamicity, y_dseptin) #
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([joint_fluc_ratio_negkappa.ravel(), mean_pos_joint_face_select_kappa_negkappa_delta_sept_timeseries_flat.ravel()])
    kernel = spstats.gaussian_kde(values) # test different bandwith. 
    f = np.reshape(kernel(positions).T, xx.shape) # curvature x, septin y. 
    
    values_neg = np.vstack([joint_fluc_ratio_negkappa.ravel(), np.abs(mean_neg_joint_face_select_kappa_negkappa_delta_sept_timeseries_flat.ravel())])
    kernel_neg = spstats.gaussian_kde(values_neg) # test different bandwith. 
    f_neg = np.reshape(kernel_neg(positions).T, xx.shape) # curvature x, septin y. 
    
    """
    integrate to get the marginal!. ( in this case )
    """
    # get the marginal and check this is correct!. 
    marginal_pos_dseptin_dynamicity_negkappa = np.nansum(f*y_dseptin[:,None], axis=0) / np.nansum(f, axis=0).astype(np.float)
    var_marginal_pos_dseptin_dynamicity_negkappa = np.nansum(f*(y_dseptin[:,None]**2), axis=0) / np.nansum(f, axis=0).astype(np.float) - marginal_pos_dseptin_dynamicity_negkappa**2
    sem_marginal_pos_dseptin_dynamicity_negkappa = np.sqrt(var_marginal_pos_dseptin_dynamicity_negkappa) / (np.sqrt(np.nansum(f, axis=0).astype(np.float)+1e-15))
    
    marginal_neg_dseptin_dynamicity_negkappa = np.nansum(f_neg*y_dseptin[:,None], axis=0) / np.nansum(f_neg, axis=0).astype(np.float)
    var_marginal_neg_dseptin_dynamicity_negkappa = np.nansum(f_neg*(y_dseptin[:,None]**2), axis=0) / np.nansum(f_neg, axis=0).astype(np.float) - marginal_neg_dseptin_dynamicity_negkappa**2
    sem_marginal_neg_dseptin_dynamicity_negkappa = np.sqrt(var_marginal_neg_dseptin_dynamicity_negkappa) / (np.sqrt(np.nansum(f_neg, axis=0).astype(np.float)+1e-15))
    
    
    plt.figure(figsize=(5,5))
    plt.title('for neg kappa')
    plt.plot(x_dynamicity, marginal_pos_dseptin_dynamicity_negkappa)
    plt.fill_between(x_dynamicity, 
                     marginal_pos_dseptin_dynamicity_negkappa-sem_marginal_pos_dseptin_dynamicity_negkappa, 
                     marginal_pos_dseptin_dynamicity_negkappa+sem_marginal_pos_dseptin_dynamicity_negkappa, alpha=0.5, label='delta_I_increase')
    plt.plot(x_dynamicity, marginal_neg_dseptin_dynamicity_negkappa)
    plt.fill_between(x_dynamicity, 
                     marginal_neg_dseptin_dynamicity_negkappa-sem_marginal_neg_dseptin_dynamicity_negkappa, 
                     marginal_neg_dseptin_dynamicity_negkappa+sem_marginal_neg_dseptin_dynamicity_negkappa, alpha=0.5,label='delta_I_decrease')
    plt.tick_params(length=5, right=True)
    plt.vlines(fluc_thresh_neg, 0, 0.1, linestyles='dashed', color='k')
    plt.ylim([0,0.1])
    plt.xlim([0,1])
    plt.xlabel('Dynamicity [0-1]')
    plt.ylabel('mean abs. Delta I')
    plt.legend(loc='best')
    plt.savefig(os.path.join(analysis_output_folder, 
                              'integrateKDE_Dynamicity_vs_dI_on_stable_deltaFaces_contigs_negKappa.svg'), dpi=300, bbox_inches='tight')
    plt.show()
    
    # """
    # checking the kde computation. 
    # """
    # dataneg_posdelta = pd.DataFrame(np.vstack([np.hstack([joint_fluc_ratio_negkappa, joint_fluc_ratio_negkappa]),
    #                                             np.hstack([mean_pos_joint_face_select_kappa_negkappa_delta_sept_timeseries_flat,
    #                                                       np.abs(mean_neg_joint_face_select_kappa_negkappa_delta_sept_timeseries_flat)]),
    #                                             np.hstack([np.ones(len(joint_fluc_ratio_negkappa)), np.zeros(len(joint_fluc_ratio_negkappa))])]).T,
    #                         columns =['dynamicity','pos_dSeptin', 'pos_neg'])
    # g = sns.jointplot(
    # data=dataneg_posdelta,
    # # x="dynamicity", y="abs(dSeptin)", 
    # x="dynamicity", y="pos_dSeptin", hue='pos_neg',
    # kind="kde",
    # )
    
    # # all_mean_bleb_dist_contigs = np.hstack([ np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_negkappa_bleb_distance_timeseries]), 
    # #                                         np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_poskappa_bleb_distance_timeseries]),
    # #                                         np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_flatkappa_bleb_distance_timeseries]) ])
    # dataneg_posdelta = pd.DataFrame(np.vstack([np.hstack([np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_negkappa_bleb_distance_timeseries]), 
    #                                                       np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_negkappa_bleb_distance_timeseries])]),
    #                                             np.hstack([mean_pos_joint_face_select_kappa_negkappa_delta_sept_timeseries_flat,
    #                                                       np.abs(mean_neg_joint_face_select_kappa_negkappa_delta_sept_timeseries_flat)]),
    #                                             np.hstack([np.ones(len(joint_fluc_ratio_negkappa)), np.zeros(len(joint_fluc_ratio_negkappa))])]).T,
    #                         columns =['dynamicity','pos_dSeptin', 'pos_neg'])
    # g = sns.jointplot(
    # data=dataneg_posdelta,
    # # x="dynamicity", y="abs(dSeptin)", 
    # x="dynamicity", y="pos_dSeptin", hue='pos_neg',
    # kind="kde",
    # )
    
    
    # dataneg_posdelta = pd.DataFrame(np.vstack([mean_pos_joint_face_select_kappa_negkappa_delta_sept_timeseries_flat,
    #                                            np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_negkappa_bleb_distance_timeseries])]).T,
    #                         columns =['dynamicity','pos_dSeptin'])
    # g = sns.jointplot(
    # data=dataneg_posdelta,
    # # x="dynamicity", y="abs(dSeptin)", 
    # x="dynamicity", y="pos_dSeptin",
    # kind="kde",
    # )
    



    """
    now on positive faces. 
    """
    mean_pos_joint_face_select_kappa_poskappa_delta_sept_timeseries_flat = np.hstack([np.nanmean(ss[ss>0]) for ss in joint_face_select_kappa_poskappa_delta_sept_timeseries_flat])
    mean_neg_joint_face_select_kappa_poskappa_delta_sept_timeseries_flat = np.hstack([np.nanmean(ss[ss<0]) for ss in joint_face_select_kappa_poskappa_delta_sept_timeseries_flat])
    
    y_dseptin = np.linspace(0, 0.1,100) # this will be the upper bound.
    x_dynamicity = np.linspace(0, 1., 100)
    # the above is more accurate with 2D kde. 
    # Peform the kernel density estimate
    xx, yy = np.meshgrid(x_dynamicity, y_dseptin) #
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([joint_fluc_ratio_poskappa.ravel(), 
                        mean_pos_joint_face_select_kappa_poskappa_delta_sept_timeseries_flat.ravel()])
    kernel = spstats.gaussian_kde(values) # test different bandwith. 
    f = np.reshape(kernel(positions).T, xx.shape) # curvature x, septin y. 
    
    values_neg = np.vstack([joint_fluc_ratio_poskappa.ravel(), 
                np.abs(mean_neg_joint_face_select_kappa_poskappa_delta_sept_timeseries_flat.ravel())])
    kernel_neg = spstats.gaussian_kde(values_neg) # test different bandwith. 
    f_neg = np.reshape(kernel_neg(positions).T, xx.shape) # curvature x, septin y. 
    
    """
    integrate to get the marginal!. ( in this case )
    """
    # get the marginal and check this is correct!. 
    marginal_pos_dseptin_dynamicity_poskappa = np.nansum(f*y_dseptin[:,None], axis=0) / np.nansum(f, axis=0).astype(np.float)
    var_marginal_pos_dseptin_dynamicity_poskappa = np.nansum(f*(y_dseptin[:,None]**2), axis=0) / np.nansum(f, axis=0).astype(np.float) - marginal_pos_dseptin_dynamicity_poskappa**2
    sem_marginal_pos_dseptin_dynamicity_poskappa = np.sqrt(var_marginal_pos_dseptin_dynamicity_poskappa) / (np.sqrt(np.nansum(f, axis=0).astype(np.float)+1e-15))
    
    marginal_neg_dseptin_dynamicity_poskappa = np.nansum(f_neg*y_dseptin[:,None], axis=0) / np.nansum(f_neg, axis=0).astype(np.float)
    var_marginal_neg_dseptin_dynamicity_poskappa = np.nansum(f_neg*(y_dseptin[:,None]**2), axis=0) / np.nansum(f_neg, axis=0).astype(np.float) - marginal_neg_dseptin_dynamicity_poskappa**2
    sem_marginal_neg_dseptin_dynamicity_poskappa = np.sqrt(var_marginal_neg_dseptin_dynamicity_poskappa) / (np.sqrt(np.nansum(f_neg, axis=0).astype(np.float)+1e-15))
    
    
    plt.figure(figsize=(5,5))
    plt.title('for pos kappa')
    plt.plot(x_dynamicity, marginal_pos_dseptin_dynamicity_poskappa)
    plt.fill_between(x_dynamicity, 
                     marginal_pos_dseptin_dynamicity_poskappa-sem_marginal_pos_dseptin_dynamicity_poskappa, 
                     marginal_pos_dseptin_dynamicity_poskappa+sem_marginal_pos_dseptin_dynamicity_poskappa, alpha=0.5, label='delta_I_increase')
    plt.plot(x_dynamicity, marginal_neg_dseptin_dynamicity_poskappa)
    plt.fill_between(x_dynamicity, 
                     marginal_neg_dseptin_dynamicity_poskappa-sem_marginal_neg_dseptin_dynamicity_poskappa, 
                     marginal_neg_dseptin_dynamicity_poskappa+sem_marginal_neg_dseptin_dynamicity_poskappa, alpha=0.5,label='delta_I_decrease')
    plt.tick_params(length=5, right=True)
    plt.vlines(fluc_thresh_pos, 0, 0.1, linestyles='dashed', color='k')
    plt.ylim([0,0.1])
    plt.xlim([0,1])
    plt.xlabel('Dynamicity [0-1]')
    plt.ylabel('mean abs. Delta I')
    plt.legend(loc='best')
    plt.savefig(os.path.join(analysis_output_folder, 
                              'integrateKDE_Dynamicity_vs_dI_on_stable_deltaFaces_contigs_posKappa.svg'), dpi=300, bbox_inches='tight')
    plt.show()
    
    dataneg_posdelta = pd.DataFrame(np.vstack([np.hstack([joint_fluc_ratio_poskappa, joint_fluc_ratio_poskappa]),
                                                np.hstack([mean_pos_joint_face_select_kappa_poskappa_delta_sept_timeseries_flat,
                                                          np.abs(mean_neg_joint_face_select_kappa_poskappa_delta_sept_timeseries_flat)]),
                                                np.hstack([np.ones(len(joint_fluc_ratio_poskappa)), np.zeros(len(joint_fluc_ratio_poskappa))])]).T,
                                    columns =['dynamicity','pos_dSeptin', 'pos_neg'])
    g = sns.jointplot(
    data=dataneg_posdelta,
    # x="dynamicity", y="abs(dSeptin)", 
    x="dynamicity", y="pos_dSeptin", hue='pos_neg',
    kind="kde",
    )
    
    
    """
    now on flat faces. 
    """
    mean_pos_joint_face_select_kappa_flatkappa_delta_sept_timeseries_flat = np.hstack([np.nanmean(ss[ss>0]) for ss in joint_face_select_kappa_flatkappa_delta_sept_timeseries_flat])
    mean_neg_joint_face_select_kappa_flatkappa_delta_sept_timeseries_flat = np.hstack([np.nanmean(ss[ss<0]) for ss in joint_face_select_kappa_flatkappa_delta_sept_timeseries_flat])
    
    y_dseptin = np.linspace(0, 0.1,100) # this will be the upper bound.
    x_dynamicity = np.linspace(0, 1., 100)
    # the above is more accurate with 2D kde. 
    # Peform the kernel density estimate
    xx, yy = np.meshgrid(x_dynamicity, y_dseptin) #
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([joint_fluc_ratio_flatkappa.ravel(), 
                        mean_pos_joint_face_select_kappa_flatkappa_delta_sept_timeseries_flat.ravel()])
    kernel = spstats.gaussian_kde(values) # test different bandwith. 
    f = np.reshape(kernel(positions).T, xx.shape) # curvature x, septin y. 
    
    values_neg = np.vstack([joint_fluc_ratio_flatkappa.ravel(), 
                np.abs(mean_neg_joint_face_select_kappa_flatkappa_delta_sept_timeseries_flat.ravel())])
    kernel_neg = spstats.gaussian_kde(values_neg) # test different bandwith. 
    f_neg = np.reshape(kernel_neg(positions).T, xx.shape) # curvature x, septin y. 
    
    """
    integrate to get the marginal!. ( in this case )
    """
    # get the marginal and check this is correct!. 
    marginal_pos_dseptin_dynamicity_flatkappa = np.nansum(f*y_dseptin[:,None], axis=0) / np.nansum(f, axis=0).astype(np.float)
    var_marginal_pos_dseptin_dynamicity_flatkappa = np.nansum(f*(y_dseptin[:,None]**2), axis=0) / np.nansum(f, axis=0).astype(np.float) - marginal_pos_dseptin_dynamicity_flatkappa**2
    sem_marginal_pos_dseptin_dynamicity_flatkappa = np.sqrt(var_marginal_pos_dseptin_dynamicity_flatkappa) / (np.sqrt(np.nansum(f, axis=0).astype(np.float)+1e-15))
    
    marginal_neg_dseptin_dynamicity_flatkappa = np.nansum(f_neg*y_dseptin[:,None], axis=0) / np.nansum(f_neg, axis=0).astype(np.float)
    var_marginal_neg_dseptin_dynamicity_flatkappa = np.nansum(f_neg*(y_dseptin[:,None]**2), axis=0) / np.nansum(f_neg, axis=0).astype(np.float) - marginal_neg_dseptin_dynamicity_flatkappa**2
    sem_marginal_neg_dseptin_dynamicity_flatkappa = np.sqrt(var_marginal_neg_dseptin_dynamicity_flatkappa) / (np.sqrt(np.nansum(f_neg, axis=0).astype(np.float)+1e-15))
    
    
    plt.figure(figsize=(5,5))
    plt.title('for flat kappa')
    plt.plot(x_dynamicity, marginal_pos_dseptin_dynamicity_flatkappa)
    plt.fill_between(x_dynamicity, 
                     marginal_pos_dseptin_dynamicity_flatkappa-sem_marginal_pos_dseptin_dynamicity_flatkappa, 
                     marginal_pos_dseptin_dynamicity_flatkappa+sem_marginal_pos_dseptin_dynamicity_flatkappa, alpha=0.5, label='delta_I_increase')
    plt.plot(x_dynamicity, marginal_neg_dseptin_dynamicity_flatkappa)
    plt.fill_between(x_dynamicity, 
                     marginal_neg_dseptin_dynamicity_flatkappa-sem_marginal_neg_dseptin_dynamicity_flatkappa, 
                     marginal_neg_dseptin_dynamicity_flatkappa+sem_marginal_neg_dseptin_dynamicity_flatkappa, alpha=0.5,label='delta_I_decrease')
    plt.tick_params(length=5, right=True)
    plt.vlines(fluc_thresh_flat, 0, 0.1, linestyles='dashed', color='k')
    plt.ylim([0,0.1])
    plt.xlim([0,1])
    plt.xlabel('Dynamicity [0-1]')
    plt.ylabel('mean abs. Delta I')
    plt.legend(loc='best')
    plt.savefig(os.path.join(analysis_output_folder, 
                              'integrateKDE_Dynamicity_vs_dI_on_stable_deltaFaces_contigs_flatKappa.svg'), dpi=300, bbox_inches='tight')
    plt.show()
    
    dataflat_posdelta = pd.DataFrame(np.vstack([np.hstack([joint_fluc_ratio_flatkappa, joint_fluc_ratio_flatkappa]),
                                                np.hstack([mean_pos_joint_face_select_kappa_flatkappa_delta_sept_timeseries_flat,
                                                          np.abs(mean_neg_joint_face_select_kappa_flatkappa_delta_sept_timeseries_flat)]),
                                                np.hstack([np.ones(len(joint_fluc_ratio_flatkappa)), np.zeros(len(joint_fluc_ratio_flatkappa))])]).T,
                                    columns =['dynamicity','pos_dSeptin', 'pos_neg'])
    g = sns.jointplot(
    data=dataflat_posdelta,
    # x="dynamicity", y="abs(dSeptin)", 
    x="dynamicity", y="pos_dSeptin", hue='pos_neg',
    kind="kde",
    )
    

    """
    ok now we have a look at the dI vs distance  
    """
    y_bleb_distance = np.linspace(0, 2.5, 100) # this will be the upper bound.
    x_dseptin_I = np.linspace(0., 0.1, 100)
    # the above is more accurate with 2D kde. 
    # Peform the kernel density estimate
    xx, yy = np.meshgrid(x_dseptin_I, y_bleb_distance) #
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([np.hstack(mean_pos_joint_face_select_kappa_negkappa_delta_sept_timeseries_flat), 
                        np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_negkappa_bleb_distance_timeseries])])
    kernel = spstats.gaussian_kde(values) # test different bandwith. 
    f = np.reshape(kernel(positions).T, xx.shape) # curvature x, septin y. 
    
    values_dec = np.vstack([np.abs(np.hstack(mean_neg_joint_face_select_kappa_negkappa_delta_sept_timeseries_flat)), 
                        np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_negkappa_bleb_distance_timeseries])])
    kernel_dec = spstats.gaussian_kde(values_dec) # test different bandwith. 
    f_dec = np.reshape(kernel_dec(positions).T, xx.shape) # curvature x, septin y. 

    # values_flatKappa = np.vstack([np.hstack([np.nanmean(ss) for ss in joint_face_select_kappa_flatkappa_sept_timeseries_flat]), 
    #                     np.hstack([np.nanmedian(ss) for ss in joint_face_select_kappa_flatkappa_bleb_distance_timeseries])])
    # kernel_flatKappa = spstats.gaussian_kde(values_flatKappa) # test different bandwith. 
    # f_flatKappa = np.reshape(kernel_flatKappa(positions).T, xx.shape) # curvature x, septin y. 
    
    """
    integrate to get the marginal!. ( in this case )
    """
    # get the marginal and check this is correct!. 
    marginal_dI_blebs_negkappa = np.nansum(f*y_bleb_distance[:,None], axis=0) / np.nansum(f, axis=0).astype(np.float)
    var_marginal_dI_blebs_negkappa = np.nansum(f*(y_bleb_distance[:,None]**2), axis=0) / np.nansum(f, axis=0).astype(np.float) - marginal_dI_blebs_negkappa**2
    sem_marginal_dI_blebs_negkappa = np.sqrt(var_marginal_dI_blebs_negkappa) / (np.sqrt(np.nansum(f, axis=0).astype(np.float)+1))
    
    marginal_dI_blebs_negkappa_dec = np.nansum(f_dec*y_bleb_distance[:,None], axis=0) / np.nansum(f_dec, axis=0).astype(np.float)
    var_marginal_dI_blebs_negkappa_dec = np.nansum(f_dec*(y_bleb_distance[:,None]**2), axis=0) / np.nansum(f_dec, axis=0).astype(np.float) - marginal_dI_blebs_negkappa_dec**2
    sem_marginal_dI_blebs_negkappa_dec = np.sqrt(var_marginal_dI_blebs_negkappa_dec) / (np.sqrt(np.nansum(f_dec, axis=0).astype(np.float)+1))
    
    # marginal_I_blebs_flatkappa = np.nansum(f_flatKappa*y_bleb_distance[:,None], axis=0) / np.nansum(f_flatKappa, axis=0).astype(np.float)
    # var_marginal_I_blebs_flatkappa = np.nansum(f_flatKappa*(y_bleb_distance[:,None]**2), axis=0) / np.nansum(f_flatKappa, axis=0).astype(np.float) - marginal_I_blebs_flatkappa**2
    # sem_marginal_I_blebs_flatkappa = np.sqrt(var_marginal_I_blebs_flatkappa) / (np.sqrt(np.nansum(f_flatKappa, axis=0).astype(np.float)+1e-15))
    
    
    plt.figure(figsize=(5,5))
    plt.title('for neg kappa')
    plt.plot(x_dseptin_I, marginal_dI_blebs_negkappa)
    plt.fill_between(x_dseptin_I, 
                     marginal_dI_blebs_negkappa-sem_marginal_dI_blebs_negkappa, 
                     marginal_dI_blebs_negkappa+sem_marginal_dI_blebs_negkappa, alpha=0.5, label='increase_negKappa')
    plt.plot(x_dseptin_I, marginal_dI_blebs_negkappa_dec)
    plt.fill_between(x_dseptin_I, 
                     marginal_dI_blebs_negkappa_dec-sem_marginal_dI_blebs_negkappa_dec, 
                     marginal_dI_blebs_negkappa_dec+sem_marginal_dI_blebs_negkappa_dec, alpha=0.5, label='decrease_negKappa')
    # plt.plot(x_septin_I, marginal_I_blebs_poskappa)
    # plt.fill_between(x_dseptin_I, 
    #                  marginal_I_blebs_poskappa-sem_marginal_I_blebs_poskappa, 
    #                  marginal_I_blebs_poskappa+sem_marginal_I_blebs_poskappa, alpha=0.5, label='posKappa')
    # plt.plot(x_septin_I, marginal_I_blebs_flatkappa)
    # plt.fill_between(x_dseptin_I, 
    #                  marginal_I_blebs_flatkappa-sem_marginal_I_blebs_flatkappa, 
    #                  marginal_I_blebs_flatkappa+sem_marginal_I_blebs_flatkappa, alpha=0.5, label='flatKappa')
    # plt.tick_params(length=5, right=True)
    # plt.ylim([0,0.1])
    # plt.xlim([0,1])
    plt.xlabel('Delta Norm. Septin I')
    plt.ylabel('Blebby distance')
    plt.legend(loc='best')
    # plt.savefig(os.path.join(analysis_output_folder, 
    #                           'integrateKDE_Dynamicity_vs_dI_on_stable_deltaFaces_contigs_flatKappa.svg'), dpi=300, bbox_inches='tight')
    plt.show()







    
    
    """
    what does increase alone look like on the different surfaces. 
    """
    plt.figure(figsize=(5,5))
    plt.title('for all kappa increase')
    plt.plot(x_dynamicity, marginal_pos_dseptin_dynamicity_negkappa)
    plt.fill_between(x_dynamicity, 
                     marginal_pos_dseptin_dynamicity_negkappa-sem_marginal_pos_dseptin_dynamicity_negkappa, 
                     marginal_pos_dseptin_dynamicity_negkappa+sem_marginal_pos_dseptin_dynamicity_negkappa, alpha=0.5, label='delta_I_increase_negkappa')
    plt.plot(x_dynamicity, marginal_pos_dseptin_dynamicity_poskappa)
    plt.fill_between(x_dynamicity, 
                     marginal_pos_dseptin_dynamicity_poskappa-sem_marginal_pos_dseptin_dynamicity_poskappa, 
                     marginal_pos_dseptin_dynamicity_poskappa+sem_marginal_pos_dseptin_dynamicity_poskappa, alpha=0.5, label='delta_I_increase_poskappa')
    plt.plot(x_dynamicity, marginal_pos_dseptin_dynamicity_flatkappa)
    plt.fill_between(x_dynamicity, 
                     marginal_pos_dseptin_dynamicity_flatkappa-sem_marginal_pos_dseptin_dynamicity_flatkappa, 
                     marginal_pos_dseptin_dynamicity_flatkappa+sem_marginal_pos_dseptin_dynamicity_flatkappa, alpha=0.5, label='delta_I_increase_flatkappa')
    plt.tick_params(length=5, right=True)
    # plt.vlines(fluc_thresh_flat, 0, 0.1, linestyles='dashed', color='k')
    plt.ylim([0,0.1])
    plt.xlim([0,1])
    plt.xlabel('Dynamicity [0-1]')
    plt.ylabel('mean abs. Delta I')
    plt.legend(loc='best')
    plt.savefig(os.path.join(analysis_output_folder, 
                              'integrateKDE_Dynamicity_vs_increase_dI_on_stable_deltaFaces_contigs_allKappa.svg'), dpi=300, bbox_inches='tight')
    plt.show()


    plt.figure(figsize=(5,5))
    plt.title('for all kappa decrease')
    plt.plot(x_dynamicity, marginal_neg_dseptin_dynamicity_negkappa)
    plt.fill_between(x_dynamicity, 
                     marginal_neg_dseptin_dynamicity_poskappa-sem_marginal_neg_dseptin_dynamicity_negkappa, 
                     marginal_neg_dseptin_dynamicity_poskappa+sem_marginal_neg_dseptin_dynamicity_negkappa, alpha=0.5,label='delta_I_decrease_negkappa')
    plt.plot(x_dynamicity, marginal_neg_dseptin_dynamicity_poskappa)
    plt.fill_between(x_dynamicity, 
                     marginal_neg_dseptin_dynamicity_poskappa-sem_marginal_neg_dseptin_dynamicity_poskappa, 
                     marginal_neg_dseptin_dynamicity_poskappa+sem_marginal_neg_dseptin_dynamicity_poskappa, alpha=0.5,label='delta_I_decrease_poskappa')
    plt.plot(x_dynamicity, marginal_neg_dseptin_dynamicity_flatkappa)
    plt.fill_between(x_dynamicity, 
                     marginal_neg_dseptin_dynamicity_flatkappa-sem_marginal_neg_dseptin_dynamicity_flatkappa, 
                     marginal_neg_dseptin_dynamicity_flatkappa+sem_marginal_neg_dseptin_dynamicity_flatkappa, alpha=0.5,label='delta_I_decrease_flatkappa')
    plt.tick_params(length=5, right=True)
    # plt.vlines(fluc_thresh_flat, 0, 0.1, linestyles='dashed', color='k')
    plt.ylim([0,0.1])
    plt.xlim([0,1])
    plt.xlabel('Dynamicity [0-1]')
    plt.ylabel('mean abs. Delta I')
    plt.legend(loc='best')
    plt.savefig(os.path.join(analysis_output_folder, 
                              'integrateKDE_Dynamicity_vs_decrease_dI_on_stable_deltaFaces_contigs_allKappa.svg'), dpi=300, bbox_inches='tight')
    plt.show()
    
    
# =============================================================================
# =============================================================================
# #     Do a contig analysis on different curvature based on kappa vs intensity... -- potentially we pool together the various kappa surfaces. 
# =============================================================================
# =============================================================================
    
    
    
    
    
    
    
    
    

    
# # =============================================================================
# #     Repeat on the decreasing faces? i think we need to gate for whether actually decreasing segments!. 
# # =============================================================================
    
#     # categorise by the mean positive / negative faces... ( the thresholds were derived from the smoothed ... )
#     neg_face_select_kappa = all_kappa_groups_faces[:,neg_face_select>0].copy() # get the face classification for all the positive increasing faces. 
    
#     """
#     filter and retain the minimal length contigs. 
#     """
#     neg_face_select_kappa_poskappa = parse_all_contigs(neg_face_select_kappa[:,:].T, kappa_class=0, min_len=25)
#     neg_face_select_kappa_negkappa = parse_all_contigs(neg_face_select_kappa[:,:].T, kappa_class=2, min_len=25)
    
#     neg_face_select_kappa_poskappa_kappa_timeseries = get_timeseries(-kappa_all_raw_reduce[:,neg_face_select>0].T, neg_face_select_kappa_poskappa)
#     neg_face_select_kappa_poskappa_sept_timeseries = get_timeseries(septin_faces_all_raw_reduce[:,neg_face_select>0].T, neg_face_select_kappa_poskappa)
#     neg_face_select_kappa_poskappa_neg_septin_timeseries = get_timeseries(all_d_septin_hotspots_binary_neg[:,neg_face_select>0].T, neg_face_select_kappa_poskappa)
    
#     neg_face_select_kappa_negkappa_kappa_timeseries = get_timeseries(-kappa_all_raw_reduce[:,neg_face_select>0].T, neg_face_select_kappa_negkappa)
#     neg_face_select_kappa_negkappa_sept_timeseries = get_timeseries(septin_faces_all_raw_reduce[:,neg_face_select>0].T, neg_face_select_kappa_negkappa)
#     neg_face_select_kappa_negkappa_neg_septin_timeseries = get_timeseries(all_d_septin_hotspots_binary_neg[:,neg_face_select>0].T, neg_face_select_kappa_negkappa)
    
#     # flatten these timeseries. 
#     neg_face_select_kappa_poskappa_kappa_timeseries_flat = [item for sublist in neg_face_select_kappa_poskappa_kappa_timeseries for item in sublist if len(item)>0]
#     neg_face_select_kappa_poskappa_sept_timeseries_flat = [item for sublist in neg_face_select_kappa_poskappa_sept_timeseries for item in sublist if len(item)>0]
#     neg_face_select_kappa_poskappa_neg_septin_timeseries_flat = [item for sublist in neg_face_select_kappa_poskappa_neg_septin_timeseries for item in sublist if len(item)>0]
#     neg_face_select_kappa_negkappa_kappa_timeseries_flat = [item for sublist in neg_face_select_kappa_negkappa_kappa_timeseries for item in sublist if len(item)>0]
#     neg_face_select_kappa_negkappa_sept_timeseries_flat = [item for sublist in neg_face_select_kappa_negkappa_sept_timeseries for item in sublist if len(item)>0]
#     neg_face_select_kappa_negkappa_neg_septin_timeseries_flat = [item for sublist in neg_face_select_kappa_negkappa_neg_septin_timeseries for item in sublist if len(item)>0]
    
#     neg_face_select_kappa_poskappa_neg_septin_timeseries_flat_counts = np.hstack([np.sum(cc) for cc in neg_face_select_kappa_poskappa_neg_septin_timeseries_flat])
#     neg_face_select_kappa_negkappa_neg_septin_timeseries_flat_counts = np.hstack([np.sum(cc) for cc in neg_face_select_kappa_negkappa_neg_septin_timeseries_flat])
    
    
#     neg_face_pos_kappa_corr = tsa.xcorr_timeseries_set_1d(neg_face_select_kappa_poskappa_kappa_timeseries_flat, 
#                                                           neg_face_select_kappa_poskappa_sept_timeseries_flat, 
#                                                           norm=True, eps=1e-12,stack_final=False)
#     # accumulate the lags and correlations from each timeseries!. 
#     neg_face_pos_kappa_corr_arr = tsa.stack_xcorr_curves(neg_face_pos_kappa_corr)
    
#     plt.figure()
#     plt.title('decreasing pos kappa')
#     plt.plot(np.nanmean(neg_face_pos_kappa_corr_arr, axis=0)); plt.vlines(neg_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='k')
#     plt.plot(np.nanmedian(neg_face_pos_kappa_corr_arr[neg_face_select_kappa_poskappa_neg_septin_timeseries_flat_counts>=15], axis=0), color='r'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='r')
#     plt.plot(np.nanmedian(neg_face_pos_kappa_corr_arr[neg_face_select_kappa_poskappa_neg_septin_timeseries_flat_counts<15], axis=0), color='g'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='b')
#     plt.xlim([neg_face_pos_kappa_corr_arr.shape[1]//2-100, neg_face_pos_kappa_corr_arr.shape[1]//2+100])
#     plt.show()
    
#     neg_face_neg_kappa_corr = tsa.xcorr_timeseries_set_1d(neg_face_select_kappa_negkappa_kappa_timeseries_flat, 
#                                                           neg_face_select_kappa_negkappa_sept_timeseries_flat, 
#                                                           norm=True, eps=1e-12,stack_final=False)
#     # accumulate the lags and correlations from each timeseries!. 
#     neg_face_neg_kappa_corr_arr = tsa.stack_xcorr_curves(neg_face_neg_kappa_corr)
    
#     plt.figure()
#     plt.title('decreasing neg kappa')
#     plt.plot(np.nanmedian(neg_face_neg_kappa_corr_arr, axis=0)); plt.vlines(neg_face_neg_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='k')
#     plt.plot(np.nanmedian(neg_face_neg_kappa_corr_arr[neg_face_select_kappa_negkappa_neg_septin_timeseries_flat_counts>=15], axis=0), color='r'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='b')
#     plt.plot(np.nanmedian(neg_face_neg_kappa_corr_arr[neg_face_select_kappa_negkappa_neg_septin_timeseries_flat_counts<15], axis=0), color='g'); #plt.vlines(pos_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='b')
#     # plt.plot(np.nanmedian(neg_face_neg_kappa_corr_arr, axis=0)); plt.vlines(neg_face_neg_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='k')
#     plt.xlim([neg_face_neg_kappa_corr_arr.shape[1]//2-100, neg_face_neg_kappa_corr_arr.shape[1]//2+100])
#     plt.show()
    
    
#     # plt.figure()
#     # plt.plot(np.nanmedian(neg_face_pos_kappa_corr_arr, axis=0)); plt.vlines(neg_face_pos_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='k')
#     # plt.plot(np.nanmedian(neg_face_neg_kappa_corr_arr, axis=0)); plt.vlines(neg_face_neg_kappa_corr_arr.shape[1]//2+1, -0.5,.2, color='k')
#     # plt.xlim([neg_face_neg_kappa_corr_arr.shape[1]//2-100, neg_face_neg_kappa_corr_arr.shape[1]//2+100])
#     # plt.show()
    
    
    
# # =============================================================================
# #     We should pool the patches of increasing and decreasing septin and classify segments by increasing/decreasing majority ( are there both consistent segments? )
# # =============================================================================
    
#     theta = np.linspace(0,4*np.pi,200)
#     x = np.cos(theta)
#     # y = np.cos(theta+np.pi/6.)
#     y = np.sin(theta)
#     # corr_xy = np.correlate(x,y, mode='full')
#     corr_xy = tsa.xcorr(x,y)
    
#     plt.figure()
#     plt.plot(theta, x, 'r')
#     plt.plot(theta, y, 'g')
#     plt.show()

#     plt.figure()
#     plt.plot(np.linspace(-len(corr_xy)//2,len(corr_xy)//2, len(corr_xy)) * np.diff(theta)[0], corr_xy)
#     plt.vlines(0, -1,1)
#     plt.show()
    
#     lag = np.argmax(np.abs(corr_xy)) - (len(corr_xy)//2 + 1)
    
    
    
    # neg_blebby_count_faces = np.sum(neg_face_select_kappa==0, axis=0)
    # neg_face_select_kappa_thresh = skfilters.threshold_otsu(neg_blebby_count_faces)
    
    # neg_blebby_count_faces2 = np.sum(neg_face_select_kappa==2, axis=0)
    # neg_face_select_kappa_thresh2 = skfilters.threshold_otsu(neg_blebby_count_faces2)
    
    # neg_face_pos_kappa_corr = tsa.xcorr_timeseries_set_1d(-(kappa_all_raw[:,neg_face_select>0].T)[neg_blebby_count_faces>neg_face_select_kappa_thresh], 
    #                                                       (septin_faces_all_raw[:,neg_face_select>0].T)[neg_blebby_count_faces>neg_face_select_kappa_thresh], 
    #                                                       norm=True, eps=1e-12)
    # neg_face_neg_kappa_corr = tsa.xcorr_timeseries_set_1d(-(kappa_all_raw[:,neg_face_select>0].T)[neg_blebby_count_faces2>=neg_face_select_kappa_thresh2], 
    #                                                       (septin_faces_all_raw[:,neg_face_select>0].T)[neg_blebby_count_faces2>=neg_face_select_kappa_thresh2], 
    #                                                       norm=True, eps=1e-12)
    
    
    # plt.figure()
    # plt.plot(np.linspace(-len(xcorr_kappa_septin_pos.T)//2, len(xcorr_kappa_septin_pos.T)//2, len(xcorr_kappa_septin_pos.T)), neg_face_pos_kappa_corr.mean(axis=0))
    # plt.plot(np.linspace(-len(xcorr_kappa_septin_pos.T)//2, len(xcorr_kappa_septin_pos.T)//2, len(xcorr_kappa_septin_pos.T)), neg_face_neg_kappa_corr.mean(axis=0))
    # plt.fill_between(np.linspace(-len(xcorr_kappa_septin_pos.T)//2, len(xcorr_kappa_septin_pos.T)//2, len(xcorr_kappa_septin_pos.T)), 
    #                   neg_face_neg_kappa_corr.mean(axis=0) - 1.96*spstats.sem(neg_face_neg_kappa_corr, axis=0, nan_policy='omit'), 
    #                   neg_face_neg_kappa_corr.mean(axis=0) + 1.96*spstats.sem(neg_face_neg_kappa_corr, axis=0, nan_policy='omit'), color='orange', alpha=.5)
    # plt.plot(np.linspace(-len(xcorr_kappa_septin_pos.T)//2, len(xcorr_kappa_septin_pos.T)//2, len(xcorr_kappa_septin_pos.T)), xcorr_kappa_septin_neg.mean(axis=0), color='k')
    # plt.fill_between(np.linspace(-len(xcorr_kappa_septin_pos.T)//2, len(xcorr_kappa_septin_pos.T)//2, len(xcorr_kappa_septin_pos.T)), 
    #                   xcorr_kappa_septin_neg.mean(axis=0) - 1.96*spstats.sem(xcorr_kappa_septin_neg, axis=0, nan_policy='omit'), 
    #                   xcorr_kappa_septin_neg.mean(axis=0) + 1.96*spstats.sem(xcorr_kappa_septin_neg, axis=0, nan_policy='omit'), color='k', alpha=.5)
    # # plt.fill_between(np.linspace(-len(xcorr_kappa_septin_pos.T)//2, len(xcorr_kappa_septin_pos.T)//2, len(xcorr_kappa_septin_pos.T)), 
    # #                  xcorr_kappa_septin_neg.mean(axis=0) - np.nanstd(xcorr_kappa_septin_neg, axis=0), 
    # #                  xcorr_kappa_septin_neg.mean(axis=0) + np.nanstd(xcorr_kappa_septin_neg, axis=0), color='k', alpha=.5)
    # plt.vlines(0, -0.35,0)
    # plt.xlim([-25,25])
    # plt.hlines(0,-200,200, color='r', linestyles='dashed')
    # plt.show()
    
    
    
    # stable_septin_counts = np.sum(all_septin_hotspots_binary, axis=0)
    # stable_septin_thresh = skfilters.threshold_otsu(stable_septin_counts)
    
    # xcorr_septin = tsa.xcorr_timeseries_set_1d(-kappa_all_raw[:,stable_septin_counts>=stable_septin_thresh].T, 
    #                                             septin_faces_all_raw[:,stable_septin_counts>=stable_septin_thresh].T, 
    #                                                   norm=True, eps=1e-12)
    
    
    # """
    # what if we simply like average the raw signal in these regions? and do xcorr? --- this would yield very similar to above. 
    # """
    
    # plt.figure()
    # plt.plot(np.linspace(-len(xcorr_kappa_septin_pos.T)//2, len(xcorr_kappa_septin_pos.T)//2, len(xcorr_kappa_septin_pos.T)), xcorr_kappa_septin_pos.mean(axis=0))
    # plt.vlines(0, -0.35,0)
    # plt.xlim([-10,10])
    # plt.show()
    
    # plt.figure()
    # plt.plot(np.linspace(-len(xcorr_kappa_septin_neg.T)//2, len(xcorr_kappa_septin_neg.T)//2, len(xcorr_kappa_septin_neg.T)), xcorr_kappa_septin_neg.mean(axis=0))
    # plt.vlines(0, -0.35,0)
    # plt.xlim([-10,10])
    # plt.show()
    # """
    # do we then try to cluster the shape of these over the global pool ? or  just within each component? ( the latter is easier... and provides better summary... ) but is also equivalent to solving the global problem many many times? 
    
    # """
    # # lets try global first.... # what are the ? 
    
    # # build the global affinity matrix .... 
    
    
    # # from sklearn.metrics.pairwise import pairwise_distances
    # # from dtaidistance import dtw_ndim, dtw
    # # import dtaidistance
    
    # # all_pos_time_kappa = kappa_all_raw[:,pos_face_select>0].copy()
    # # all_pos_time_septin = septin_faces_all_raw[:,pos_face_select>0].copy()
    
    # # norm_kapp_matrix = (all_pos_time_kappa - all_pos_time_kappa.mean(axis=0)[None,:]) / np.nanstd(all_pos_time_kappa, axis=0)[None,:]
    # # norm_septin_matrix = (all_pos_time_septin - all_pos_time_septin.mean(axis=0)[None,:]) / np.nanstd(all_pos_time_septin, axis=0)[None,:]
    # # # dist_matrix = dtw.distance_matrix(norm_septin_matrix.T)
    
    # # from dtaidistance.preprocessing import differencing
    # # norm_septin_series = differencing(septin_faces_all_raw.T, smooth=0.1) # this is apparently better? 
    
    
    # """
    # clustering globally too long..... due to the affinity matrix .... 
    # """
    
    
    
    
    
    # """
    # discretizing curvature events? 
    # """
    
    
    
    
    # """
    # we want evidence of integration ? 
    # """
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # for face_plot in collapse_pt_id_faces:
    #     ax.plot([v[0][face_plot[0],0], v[0][face_plot[1],0]], 
    #             [v[0][face_plot[0],1], v[0][face_plot[1],1]],
    #             [v[0][face_plot[0],2], v[0][face_plot[1],2]], color='k', alpha=.25)
        
    #     ax.plot([v[0][face_plot[1],0], v[0][face_plot[2],0]], 
    #             [v[0][face_plot[1],1], v[0][face_plot[2],1]],
    #             [v[0][face_plot[1],2], v[0][face_plot[2],2]], color='k', alpha=.25)
        
    #     ax.plot([v[0][face_plot[2],0], v[0][face_plot[0],0]], 
    #             [v[0][face_plot[2],1], v[0][face_plot[0],1]],
    #             [v[0][face_plot[2],2], v[0][face_plot[0],2]], color='k', alpha=0.25)
    
    # for rr_ii in np.arange(len(exclude_k_rings)):
    #     rr_ii_ids = exclude_k_rings[rr_ii][pt_id]
    #     ax.scatter(v[0][rr_ii_ids,0], 
    #                 v[0][rr_ii_ids,1], 
    #                 v[0][rr_ii_ids,2],
    #                 s=200,
    #                 c=color_ring[rr_ii], zorder=10, alpha=1)
    #                 # c=septin_I_colors[rr_ii_ids], zorder=0)
                   
    # # ax.scatter(v[0][pt_id,0], 
    # #            v[0][pt_id,1], 
    # #            v[0][pt_id,2],
    # #            s=500,c='k',zorder=10000)
    # # ax.scatter(v[0][exclude_k_rings[-1][pt_id][1],0], 
    # #            v[0][exclude_k_rings[-1][pt_id][1],1], 
    # #            v[0][exclude_k_rings[-1][pt_id][1],2],
    # #            s=500,c='k',zorder=10000)
    # for pp in path:
    #     ax.scatter(v[0][pp,0], 
    #             v[0][pp,1], 
    #             v[0][pp,2],
    #             s=500,c='k',zorder=10000)
    # # grab example point on surface. 
    # ax.view_init(-45,0)
    # ax.grid('off')
    # ax.axis('off')
    # plt.show()
    
    
    
    
    
#     """
#     This definitely makes sense!
#     do this over time..
#     """
    
#     winsize = 25
#     all_histograms_min_curvature_pos = []
#     all_histograms_min_curvature_neg = []
    
#     max_absolute_curvature_faces__ = max_absolute_curvature_faces[winsize//2:len(v)-winsize//2] /.104
    
#     for window_ii, window_tt in tqdm(enumerate(np.arange(len(max_absolute_curvature_faces__))[:])):
        
#         histogram_kappa_d_septin_hotspot, bins_kappa_d_septin_hotspot = np.histogram(-1*(max_absolute_curvature_faces__[window_tt])[all_d_septin_hotspots_binary_pos[window_tt]>0], bins=100, range=(-1.5,1.5))
#         histogram_kappa_d_septin_hotspot_low, bins_kappa_d_septin_hotspot_low = np.histogram(-1*(max_absolute_curvature_faces__[window_tt])[all_d_septin_hotspots_binary_neg[window_tt]>0], bins=100, range=(-1.5,1.5))
        
#         all_histograms_min_curvature_pos.append(histogram_kappa_d_septin_hotspot)
#         all_histograms_min_curvature_neg.append(histogram_kappa_d_septin_hotspot_low)
    
#     all_histograms_min_curvature_pos = np.vstack(all_histograms_min_curvature_pos)
#     all_histograms_min_curvature_neg = np.vstack(all_histograms_min_curvature_neg)
    
    
#     all_histograms_min_curvature_pos_norm = all_histograms_min_curvature_pos/ (np.nansum(all_histograms_min_curvature_pos, axis=1) + 1e-12)[:,None]
#     all_histograms_min_curvature_neg_norm = all_histograms_min_curvature_neg/ (np.nansum(all_histograms_min_curvature_neg, axis=1) + 1e-12)[:,None]
    
#     plt.figure(figsize=(5,5))
#     plt.plot(all_bins_kappa_d_septin_hotspot[0,:-1], np.nanmean(all_histograms_min_curvature_pos_norm, axis=0), 'r')
#     plt.fill_between(all_bins_kappa_d_septin_hotspot[0,:-1], 
#                      np.nanmean(all_histograms_min_curvature_pos_norm, axis=0) - 1.96*spstats.sem(all_histograms_min_curvature_pos_norm, axis=0, nan_policy='omit'), 
#                      np.nanmean(all_histograms_min_curvature_pos_norm, axis=0) + 1.96*spstats.sem(all_histograms_min_curvature_pos_norm, axis=0, nan_policy='omit'),
#                      color='r', alpha=.25)
#     # plt.plot(all_bins_geodist_d_septin_hotspot[0,:-1], ysmooth1, 'r--')
#     plt.plot(all_bins_kappa_d_septin_hotspot[0,:-1], np.nanmean(all_histograms_min_curvature_neg_norm, axis=0))
#     plt.fill_between(all_bins_kappa_d_septin_hotspot[0,:-1], 
#                       np.nanmean(all_histograms_min_curvature_neg_norm, axis=0) - 1.96*spstats.sem(all_histograms_min_curvature_neg_norm, axis=0, nan_policy='omit'), 
#                       np.nanmean(all_histograms_min_curvature_neg_norm, axis=0) + 1.96*spstats.sem(all_histograms_min_curvature_neg_norm, axis=0, nan_policy='omit'),
#                       color='b', alpha=.25)
#     # # plt.plot(all_bins_geodist_d_septin_hotspot[0,:-1], ysmooth2, 'r--')
#     plt.xlim([-1,1])
#     plt.vlines(0,0,0.045, color='k', linestyles='dashed')
#     plt.xlabel('kappa[1/um]')
#     plt.ylabel('Norm Freq.')
#     plt.show()
    
    
#     ysmooth1 = baseline_als(np.nanmean(all_histograms_min_curvature_pos_norm, axis=0), 
#                             lam=1, p=.5, niter=10)
#     ysmooth2 = baseline_als(np.nanmean(all_histograms_min_curvature_neg_norm, axis=0), 
#                             lam=1, p=.5, niter=10)
    
#     plt.figure(figsize=(5,5))
#     plt.plot(all_bins_kappa_d_septin_hotspot[0,:-1], ysmooth1, 'r--')
#     plt.plot(all_bins_kappa_d_septin_hotspot[0,:-1], ysmooth2, 'b--')
#     plt.xlim([-1,1])
#     plt.vlines(0,0,0.045, color='k', linestyles='dashed')
#     plt.xlabel('kappa[1/um]')
#     plt.ylabel('Norm Freq.')
#     plt.show()
    
    
#     # # plotting the cumsum
    
#     # plt.figure(figsize=(5,5))
#     # # plt.plot(all_bins_geodist_d_septin_hotspot[0,:-1], np.nanmean(np.cumsum(all_histogram_geodist_d_septin_hotspot_norm, axis=1), axis=0), 'r')
#     # plt.plot(all_bins_geodist_d_septin_hotspot[0,:-1], np.cumsum(np.nanmean(all_histogram_geodist_d_septin_hotspot_norm, axis=0)), 'r')
    
#     # # plt.fill_between(all_bins_geodist_d_septin_hotspot[0,:-1], 
#     # #                  np.nanmean(all_histogram_geodist_d_septin_hotspot_norm, axis=0) - 1.96*spstats.sem(all_histogram_geodist_d_septin_hotspot_norm, axis=0, nan_policy='omit'), 
#     # #                  np.nanmean(all_histogram_geodist_d_septin_hotspot_norm, axis=0) + 1.96*spstats.sem(all_histogram_geodist_d_septin_hotspot_norm, axis=0, nan_policy='omit'),
#     # #                  color='r', alpha=.25)
#     # plt.plot(all_bins_geodist_d_septin_hotspot[0,:-1], np.nanmean(np.cumsum(all_histogram_geodist_d_septin_hotspot_norm_low, axis=1), axis=0))
#     # # plt.fill_between(all_bins_geodist_d_septin_hotspot[0,:-1], 
#     # #                  np.nanmean(all_histogram_geodist_d_septin_hotspot_norm_low, axis=0) - 1.96*spstats.sem(all_histogram_geodist_d_septin_hotspot_norm_low, axis=0, nan_policy='omit'), 
#     # #                  np.nanmean(all_histogram_geodist_d_septin_hotspot_norm_low, axis=0) + 1.96*spstats.sem(all_histogram_geodist_d_septin_hotspot_norm_low, axis=0, nan_policy='omit'),
#     # #                  color='b', alpha=.25)
#     # plt.xlim([0,3])
#     # plt.xlabel('geodesic distance[um]')
#     # plt.ylabel('Cum. Norm Freq.')
#     # plt.show()
    
    
# #     # # show for one timepoint. 
# #     # mesh = trimesh.Trimesh(vertices=v[0],
# #     #                        faces=f,
# #     #                        process=False,
# #     #                        validate=False)
    
# #     # out_last = []
    
# #     # for iter_ii in np.arange(3)[:]:
        
# #     #     if iter_ii == 0: 
# #     #         weights, out = meshtools.smooth_scalar_function_mesh(mesh, 
# #     #                                                               scalar_fn=septin_faces_all[0][:,None], 
# #     #                                                               exact=True, 
# #     #                                                               n_iters=1, 
# #     #                                                               weights=None, 
# #     #                                                               return_weights=True)
# #     #         out_last = out.copy()
# #     #     else:
# #     #         out = meshtools.smooth_scalar_function_mesh(mesh, 
# #     #                                                     scalar_fn=out_last, 
# #     #                                                     exact=True, 
# #     #                                                     n_iters=1, 
# #     #                                                     weights=weights, 
# #     #                                                     return_weights=False)
# #     #         out_last = out.copy()
# #     #     out = np.squeeze(out)
# #     #     vertex_colors = vol_colors.get_colors(out, colormap=cm.RdYlBu_r, 
# #     #                                           vmin=I_min, vmax=I_max)
        
# #     #     mesh_out = trimesh.Trimesh(vertices=v[0],
# #     #                        faces=f,
# #     #                        vertex_colors=np.uint8(255*vertex_colors),
# #     #                        process=False,
# #     #                        validate=False)
    
# #     #     mesh_out.export(os.path.join(statsfolder, 
# #     #                                  'smooth_iter-%s_correct_septin.obj' %(str(iter_ii+1).zfill(3))))
    
    
# # # =============================================================================
# # # =============================================================================
# # # #     Plot the first timepoint in matplotlib in 3D...... and choose a reasonable angle. 
# # # =============================================================================
# # # =============================================================================
    
# #     # septin color. 
# #     septin_I_colors = vol_colors.get_colors(septin_faces_all[0]*septin_vol_intensity_all[0], 
# #                                             colormap=cm.RdYlBu_r, 
# #                                             vmin=I_min, vmax=I_max)


# #     # pt_ids = np.arange(len(v[0]))[septin_faces_all[0]>1.25]
    
# #     # for iii in np.arange(505):
        
# #     #     pt_id_ = pt_ids.copy(); np.random.shuffle(pt_id_)
# #     #     # pt_id = np.random.randint(0, len(v[0]))
# #     #     pt_id = pt_id_[0]
        
        
        
#     pt_id = 25987 # this works 

#     # fig = plt.figure(figsize=(10,10))
#     fig = plt.figure(figsize=plt.figaspect(1.)*3)
#     ax = plt.axes(projection='3d', proj_type = 'ortho')
#     ax.set_box_aspect(aspect = (1,1,1))
#     plt.title(pt_id)
#     # ax.plot_trisurf(v[0][...,0], 
#     #                 v[0][...,1], 
#     #                 v[0][...,2], 
#     #                 triangles=f,
#     #                 linewidth=0.2, 
#     #                 antialiased=True)
#     ax.scatter(v[0][...,0], 
#                v[0][...,1], 
#                v[0][...,2],
#                s=.5,
#                 c=septin_I_colors, zorder=0)
#     ax.scatter(v[0][pt_id,0], 
#                v[0][pt_id,1], 
#                v[0][pt_id,2],
#                s=100,c='k',zorder=10000)
#     # grab example point on surface. 

#     ax.view_init(-45,0)
#     ax.grid('off')
#     ax.axis('off')
#     plt.show()
    
    
#     # plot corresponding the curvature too.
#     kappa_I_colors = vol_colors.get_colors(curvature_faces_all_um_mesh[0], 
#                                             colormap=cm.Spectral, 
#                                             vmin=-1, vmax=1)
    
#     fig = plt.figure(figsize=plt.figaspect(1.)*3)
#     ax = plt.axes(projection='3d', proj_type = 'ortho')
#     ax.set_box_aspect(aspect = (1,1,1))
#     plt.title(pt_id)
#     # ax.plot_trisurf(v[0][...,0], 
#     #                 v[0][...,1], 
#     #                 v[0][...,2], 
#     #                 triangles=f,
#     #                 linewidth=0.2, 
#     #                 antialiased=True)
#     ax.scatter(v[0][...,0], 
#                v[0][...,1], 
#                v[0][...,2],
#                s=.5,
#                 c=kappa_I_colors, zorder=0)
#     ax.scatter(v[0][pt_id,0], 
#                v[0][pt_id,1], 
#                v[0][pt_id,2],
#                s=100,c='k',zorder=10000)
#     # grab example point on surface. 

#     ax.view_init(-45,0)
#     ax.grid('off')
#     ax.axis('off')
#     plt.show()
    
    
#     plt.figure(figsize=(5,3))
#     plt.plot(-curvature_faces_all_um_mesh[:, pt_id], label='$\kappa[um]$')
#     plt.xlabel('Time[Frame #]')
#     plt.ylabel('Mean Curvature [1/um]')
#     plt.show()
    
#     plt.figure(figsize=(5,3))
#     plt.plot(septin_faces_all[:,pt_id], label='Septin')
#     plt.xlabel('Time[Frame #]')
#     plt.ylabel('norm. Septin I')
#     plt.show()
    
    
#     # shared version 
#     fig, ax1 = plt.subplots(figsize=(5,3))
#     ax1.plot(-curvature_faces_all_um_mesh[:, pt_id], label='$\kappa[um]$', color='r')
#     ax1.hlines(0,0,200, color='r', linestyles='dashed')
#     ax1.tick_params(axis='y', labelcolor='r')
#     ax1.set_ylabel('Mean Curvature [1/um]', color='r')
#     ax1.set_xlabel('Time[Frame #]')
#     ax2 = ax1.twinx() 
#     ax2.plot(septin_faces_all[:,pt_id], label='Septin', color='g')
#     ax2.tick_params(axis='y', labelcolor='g')
#     ax2.set_ylabel('norm. Septin I', color='g')
#     ax2.set_xlabel('Time[Frame #]')
#     # plt.legend()
#     plt.show()
    
    
        
# # =============================================================================
# #     now check the regional information around this node to highlight some points. 
# # =============================================================================

#     mesh = trimesh.Trimesh(vertices=v[0], 
#                            faces = f,
#                            process=False,
#                            validate=False)
    
#     import time 
#     t1 = time.time()
#     k_rings = meshtools.get_k_neighbor_ring(mesh, K=5, stateful=True)  # 5  is good for illustration. 
#     t2 = time.time()
#     print('k-ring constructions ', t2-t1)
    
    
#     # Let the first ring be itself
#     k_rings = [[np.hstack([iii]) for iii in np.arange(len(v[0]))]] + k_rings
#     # exclusion k_rings.. to limit the distance... computation.... 
#     exclude_k_rings = [[np.setdiff1d(k_rings[ii+1][rr], k_rings[ii][rr]) for rr in np.arange(len(k_rings[0]))] for ii in np.arange(len(k_rings)-1)]
#     # exclude_k_rings = [k_rings[0]] + exclude_k_rings   #### this can not be used to establish distance... . k_rings 0 should exclude self... .
    
#     # self has to be eliminated... indefinitely for all cases!!!! 
    
    
#     # plot just the combined version .... around the 1 node. 
#     # get all the nodes and faces to get the triangles.... 
#     collapse_pt_id = np.unique(np.hstack([rr[pt_id] for rr in exclude_k_rings]))
#     collapse_pt_id_faces = np.vstack([ff for ff in f if len(np.intersect1d(collapse_pt_id,ff))>0])
    
#     color_ring = np.vstack(sns.color_palette('coolwarm', len(exclude_k_rings)))
    
    
#     # find ids forming a short distance path.... to the outer.... for now just very quickly
#     # import scipy
#     import networkx as nx 
#     G = nx.from_scipy_sparse_matrix(igl.adjacency_matrix(f))
#     print(nx.shortest_path(G, source=pt_id, target=exclude_k_rings[-1][pt_id][1]))
#     # path_output = scipy.sparse.csgraph.shortest_path(igl.adjacency_matrix(f), 
#     #                                                  indices=pt_id)
#     path = nx.shortest_path(G, source=pt_id, target=exclude_k_rings[-1][pt_id][1])
    
    
#     fig = plt.figure(figsize=plt.figaspect(1.)*3)
#     ax = plt.axes(projection='3d', proj_type = 'ortho')
#     ax.set_box_aspect(aspect = (1,1,1))
#     plt.title(pt_id)
    
#     for face_plot in collapse_pt_id_faces:
#         ax.plot([v[0][face_plot[0],0], v[0][face_plot[1],0]], 
#                 [v[0][face_plot[0],1], v[0][face_plot[1],1]],
#                 [v[0][face_plot[0],2], v[0][face_plot[1],2]], color='k', alpha=.25)
        
#         ax.plot([v[0][face_plot[1],0], v[0][face_plot[2],0]], 
#                 [v[0][face_plot[1],1], v[0][face_plot[2],1]],
#                 [v[0][face_plot[1],2], v[0][face_plot[2],2]], color='k', alpha=.25)
        
#         ax.plot([v[0][face_plot[2],0], v[0][face_plot[0],0]], 
#                 [v[0][face_plot[2],1], v[0][face_plot[0],1]],
#                 [v[0][face_plot[2],2], v[0][face_plot[0],2]], color='k', alpha=0.25)
    
#     for rr_ii in np.arange(len(exclude_k_rings)):
#         rr_ii_ids = exclude_k_rings[rr_ii][pt_id]
#         ax.scatter(v[0][rr_ii_ids,0], 
#                    v[0][rr_ii_ids,1], 
#                    v[0][rr_ii_ids,2],
#                    s=200,
#                    c=color_ring[rr_ii], zorder=10, alpha=1)
#                    # c=septin_I_colors[rr_ii_ids], zorder=0)
                   
#     # ax.scatter(v[0][pt_id,0], 
#     #            v[0][pt_id,1], 
#     #            v[0][pt_id,2],
#     #            s=500,c='k',zorder=10000)
#     # ax.scatter(v[0][exclude_k_rings[-1][pt_id][1],0], 
#     #            v[0][exclude_k_rings[-1][pt_id][1],1], 
#     #            v[0][exclude_k_rings[-1][pt_id][1],2],
#     #            s=500,c='k',zorder=10000)
#     for pp in path:
#         ax.scatter(v[0][pp,0], 
#                v[0][pp,1], 
#                v[0][pp,2],
#                s=500,c='k',zorder=10000)
#     # grab example point on surface. 
#     ax.view_init(-45,0)
#     ax.grid('off')
#     ax.axis('off')
#     plt.show()
    
    
    

#     for pp in path:
        
#         # shared version 
#         fig, ax1 = plt.subplots(figsize=(5,3))
#         ax1.plot(-curvature_faces_all_um_mesh[:, pp], label='$\kappa[um]$', color='r')
#         ax1.hlines(0,0,200, color='r', linestyles='dashed')
#         ax1.tick_params(axis='y', labelcolor='r')
#         ax1.set_ylabel('Mean Curvature [1/um]', color='r')
#         ax1.set_xlabel('Time[Frame #]')
#         ax1.set_ylim([-1,1])
#         ax2 = ax1.twinx() 
#         ax2.plot(septin_faces_all[:,pp], label='Septin', color='g')
#         ax2.tick_params(axis='y', labelcolor='g')
#         ax2.set_ylabel('norm. Septin I', color='g')
#         ax2.set_xlabel('Time[Frame #]')
#         ax2.set_ylim([.5,1.5])
#         # plt.legend()
#         plt.show()
    
    
    
    
    
    
    
    
    
    
    
# # =============================================================================
# #     Load up the processed statistics for analysis and plotting 
# # =============================================================================

#     statsfolder = 'E:\\Work\\Projects\\Danuser-3D Causality\\Data\\Andrew\\Septin_Blebs2\\180522\\Cell5\\2022-02-09_FirstOrderCorrelation'
#     statsfile = os.path.join(statsfolder, 'curvature_septin_statistics_processed.mat')
#     stats = spio.loadmat(statsfile)
#     print(stats.keys())
    
    
#     I_fluctuations = stats['fluctuations_septin_I']
#     I_smooth = stats['baseline_septin_I']
#     Kappa_fluctuations = stats['fluctuations_curvature_um'] 
#     Kappa_smooth = stats['baseline_curvature_um']
#     kappa_thresholds = np.hstack(stats['kappa_thresholds'])
#     Kappa_groups = np.hstack(stats['curvature_groups'].ravel())
    
    
#     corrected_Sepin_I =  I_smooth + I_fluctuations
#     d_corrected_Sepin_I = np.vstack([np.gradient(II) for II in corrected_Sepin_I])
    
# # =============================================================================
# #     Derive the correlation over the time series for raw face time series., 
# # =============================================================================
    
#     """
#     Solve the heat diffusion laplacian function ..... with scalar features 
#     """
#     import igl 
#     import scipy.sparse as spsparse
    
#     deltaL = 5e-12 # controls diffusion step! 
#     kappa0 = curvature_faces_all_um_mesh[0].copy()
    
#     # # if ii ==0:
#     L = igl.cotmatrix(v[0],f)
#     M = igl.massmatrix(v[0], 
#                         f, 
#                         igl.MASSMATRIX_TYPE_BARYCENTRIC) # -> this is the only matrix that doesn't degenerate.               
#     # # # implicit solve. 
#     # S = (M - deltaL*L) # what happens when we invert?  V x V 
#     # b = M.dot(kappa0)  #M: V x V
#     # kappa0_next = spsparse.linalg.spsolve(S, kappa0)
    
    
#     # we just want to do smoothing based on igl.per_vertex_attribute_smoothing(
    
#     # this works!. 
    
#     kappa0_ = kappa0.copy()
#     # import time
#     # t1 = time.time()
#     # weights = spsparse.linalg.spsolve(M, L) # m_inv_l # this inversion is slow.....  the results for me is very similar....  but much slower!. 
#     # print(time.time() - t1) # 65s? 
    
    
#     # for iteration in np.arange(10):
        
#     #     kappa0_next = np.squeeze((np.abs(weights)).dot(kappa0_[:,None]) / (np.abs(weights).sum(axis=1))) 
#     #     kappa0_next = np.array(kappa0_next).ravel()
#     #     kappa0_ = kappa0_next.copy()
#     #     # print(kappa0_next.shape)
        
        
        
#     #### testing the implemented function version. 
#     mesh = trimesh.Trimesh(vertices=v[0], 
#                             faces=f,
#                             process=False,
#                             validate=False)
#     kappa0_ = meshtools.smooth_scalar_function_mesh(mesh, 
#                                                     scalar_fn = kappa0_[:,None], 
#                                                     exact=True, 
#                                                     n_iters=10)

#     kappa0_ = np.squeeze(kappa0_)
#     # visualize both!. 
#     I_kappa_mesh_color = vol_colors.get_colors(kappa0, colormap=cm.Spectral, vmin=-1, vmax=1)
#     I_kappa_mesh_smooth_color = vol_colors.get_colors(kappa0_, colormap=cm.Spectral, vmin=-1, vmax=1)
    
#     mesh_corr_I_kappa = trimesh.Trimesh(vertices=v[0], 
#                                           faces=f,
#                                           vertex_colors = I_kappa_mesh_color[...,:3],
#                                           validate=False,
#                                           process=False)
#     mesh_corr_I_kappa.export('2022-02-15_mesh_curvature_TP0.obj')
    
#     mesh_corr_I_kappa = trimesh.Trimesh(vertices=v[0], 
#                                           faces=f,
#                                           vertex_colors = I_kappa_mesh_smooth_color[...,:3],
#                                           validate=False,
#                                           process=False)
#     mesh_corr_I_kappa.export('2022-02-15_mesh_curvature_TP0_Lap-smooth.obj')
    
    
    
    # """
    # rederive based on the mesh fitted curvatures!
    # """
    # curvature_faces_all_um_mesh = igl.average_onto_faces(f, curvature_faces_all_um_mesh.T).T
    
    
    # corr_curvature_septin_mesh = np.hstack([spstats.pearsonr(d_corrected_Sepin_I[jjj], -curvature_faces_all_um_mesh[:,jjj])[0] for jjj in np.arange(len(Kappa_smooth))])
    # corr_curvature_septin_img = np.hstack([spstats.pearsonr(d_corrected_Sepin_I[jjj], -curvature_faces_all_um[:,jjj])[0] for jjj in np.arange(len(Kappa_smooth))])
    
    
    # # corr_I_kappa = np.hstack([spstats.pearsonr(delta_I_smooth[jjj], -Kappa_smooth[jjj])[0] for jjj in np.arange(len(Kappa_smooth))])
                              
    # fig, ax = plt.subplots(figsize=(5,5))
    # plt.plot(corr_curvature_septin_mesh, 
    #          corr_curvature_septin_img, '.')
    # plt.xlim([-1,1])
    # plt.ylim([-1,1])
    # plt.ylabel('img based curvature corr')
    # plt.xlabel('mesh based curvature corr')
    # plt.show()
    
    
# =============================================================================
#     map these onto a mesh... 
# =============================================================================
    
# # =============================================================================
# #   Put this on the mesh and see where the correlations fall 
# # =============================================================================
#     corr_I_kappa_mesh_color = vol_colors.get_colors(corr_curvature_septin_mesh, colormap=cm.coolwarm, vmin=-.3, vmax=.3) 
#     mesh_corr_I_kappa = trimesh.Trimesh(vertices=v[0], 
#                                           faces=f,
#                                           face_colors = corr_I_kappa_mesh_color[...,:3],
#                                           validate=False,
#                                           process=False)
#     mesh_corr_I_kappa.export('2022-02-14_Septin-Cell5-mesh_curvature_correlation-dIvsKappa.obj')
    
#     corr_I_kappa_mesh_color = vol_colors.get_colors(corr_curvature_septin_img, colormap=cm.coolwarm, vmin=-.3, vmax=.3) 
#     mesh_corr_I_kappa = trimesh.Trimesh(vertices=v[0], 
#                                           faces=f,
#                                           face_colors = corr_I_kappa_mesh_color[...,:3],
#                                           validate=False,
#                                           process=False)
#     mesh_corr_I_kappa.export('2022-02-14_Septin-Cell5-img_curvature_correlation-dIvsKappa.obj')
    
    
    
    
#     #### derive kappa groups on just the first frame curvature values. 
#     Kappa_groups = np.zeros(len(Kappa_groups))
#     Kappa_groups[np.logical_and(np.nanmean(Kappa_smooth[:,:25], axis=1) >= kappa_thresholds[0], np.nanmean(Kappa_smooth[:,:25], axis=1) < kappa_thresholds[1])] = 1
#     Kappa_groups[np.nanmean(Kappa_smooth[:,:25], axis=1) >= kappa_thresholds[1]] = 2
    
    
#     """
#     rederive based on the mesh fitted curvatures!
#     """
#     curvature_faces_all_um_mesh = igl.average_onto_faces(f, curvature_faces_all_um_mesh.T).T
#     # mean_kappa_mesh = np.hstack([np.mean(curvature_faces_all_um_mesh[:,iii]) for iii in np.arange(curvature_faces_all_um_mesh.shape[1])])
#     kappa_thresholds = skfilters.threshold_multiotsu(np.nanmean(curvature_faces_all_um_mesh[:25], axis=0))
    
#     Kappa_groups = np.zeros(curvature_faces_all_um_mesh.shape[1])
#     Kappa_groups[np.logical_and(np.nanmean(curvature_faces_all_um_mesh[:25], axis=0) >= kappa_thresholds[0], np.nanmean(curvature_faces_all_um_mesh[:25], axis=0) < kappa_thresholds[1])] = 1
#     Kappa_groups[np.nanmean(curvature_faces_all_um_mesh[:25], axis=0) >= kappa_thresholds[1]] = 2
    
# # =============================================================================
# #     Can we split off the components of the mesh? 
# # =============================================================================

#     import trimesh 
#     # # used instead of trimesh functions in order to keep it consistent with the splitting of mesh attributes. 
#     # if adjacency is None:
#     #     adjacency = mesh.face_adjacency

#     # # if only watertight the shortest thing we can split has 3 triangles
#     # if only_watertight:
#     #     min_len = 4
#     # else:
#     #     min_len = 1


#     mesh = trimesh.Trimesh(v[0],f, process=False, validate=False)

#     components = trimesh.graph.connected_components(
#                 edges=mesh.face_adjacency,
#                 nodes=np.arange(len(f)),
#                 min_len=1,
#                 engine=None)
    
    
#     # now try a masked version , where we drop faces. 
#     f_blebs = f[Kappa_groups==0]
    
#     mesh = trimesh.Trimesh(v[0],
#                            f_blebs, 
#                            process=False, 
#                            validate=False)
#     components = trimesh.graph.connected_components(
#                 edges=mesh.face_adjacency,
#                 nodes=np.arange(len(f)),
#                 min_len=1,
#                 engine=None)
    
#     # seems like the indexes are changed here? 
#     components_blebs = [cc for cc in components if len(cc)>1] # these should now be the separated? - if so we can use this now. 
    
# # =============================================================================
# #     Label the disconnected components now. 
# # =============================================================================
#     labels_colors = sns.color_palette('hls', n_colors = len(components_blebs))
#     labels_colors = np.vstack(labels_colors)
    
#     all_face_labels_colors = .5*np.ones((len(Kappa_smooth),3))
    
#     for iii in np.arange(len(components_blebs)):

#         select = components_blebs[iii];
#         select_original = (np.arange(len(f))[Kappa_groups==0])[select].copy()
#         all_face_labels_colors[select_original] = labels_colors[iii][None,:]
    
    
#     mesh_birch = trimesh.Trimesh(vertices=v[0], 
#                                       faces=f,
#                                       face_colors = all_face_labels_colors[...,:3],
#                                       validate=False,
#                                       process=False)
#     # mesh_birch.export('2022-02-11_Septin-Cell5-smooth_Kappa_connectedcomponents-High-First25Frame.obj')
#     mesh_birch.export('2022-02-11_Septin-Cell5-Mesh_Kappa_connectedcomponents-High-First25Frame.obj')
    
    
    
    
    
# # =============================================================================
# #     Do the first order pearson correlation with curvature. 
# # =============================================================================
#     import scipy.stats as spstats

#     # use the np.gradient function to get differences. 
#     delta_I_smooth = np.vstack([np.gradient(ss) for ss in I_smooth])
    
#     corr_I_kappa = np.hstack([spstats.pearsonr(delta_I_smooth[jjj], -Kappa_smooth[jjj])[0] for jjj in np.arange(len(Kappa_smooth))])
                              
                              
# # =============================================================================
# #   Put this on the mesh and see where the correlations fall 
# # =============================================================================
#     corr_I_kappa_mesh_color = vol_colors.get_colors(corr_I_kappa, colormap=cm.coolwarm, vmin=-1, vmax=1) 
    
#     mesh_corr_I_kappa = trimesh.Trimesh(vertices=v[0], 
#                                           faces=f,
#                                           face_colors = corr_I_kappa_mesh_color[...,:3],
#                                           validate=False,
#                                           process=False)
#     mesh_corr_I_kappa.export('2022-02-11_Septin-Cell5-smooth_curvature_correlation-IvsKappa.obj')
    
    
    
# # =============================================================================
# #     Do the first order pearson correlation with curvature. # the signal is noisy shows 50 % ?  
# # ============================================================================
#     # the distribution is very similar ... (simply the smoother is just more elevated correlation ....)
#     corr_I_kappa_raw = np.hstack([spstats.pearsonr(np.gradient(I_smooth[jjj]+I_fluctuations[jjj]), -curvature_faces_all_um[:,jjj])[0] for jjj in np.arange(len(Kappa_smooth))])
     
    
    
# # =============================================================================
# #     Step0: Output the cuvrvvatures. 
# # =============================================================================

#     labels_colors = sns.color_palette('coolwarm', n_colors = Kappa_groups.max()+1)
#     labels_colors = np.vstack(labels_colors)
    
#     face_labels_colors = labels_colors[Kappa_groups].copy()
    
    
#     # all_face_labels_colors = .5*np.ones((len(Kappa_smooth),3))
#     # all_face_labels_colors[Kappa_groups<=1] = face_labels_colors
    
    
#     """
#     The kappa thresholds is stable over different time frames - therefore we can run with this 
#     """
#     mesh_birch = trimesh.Trimesh(vertices=v[0], 
#                                           faces=f,
#                                           face_colors = face_labels_colors[...,:3],
#                                           validate=False,
#                                           process=False)
#     mesh_birch.export('2022-02-11_Septin-Cell5-smooth_Kappa_curvature_grouping.obj')


# # =============================================================================
# #   unbiased clustering as a means to find blebs --- can't cluster out blebs it seems..... ()  --- > this probably suggests that the full temporal signal is not so meaningful to cluster!. 
# # =============================================================================

#     # Kappa_smooth_shape = np.vstack([(ss-ss.mean())/ss.std() for ss in Kappa_smooth])
#     Kappa_smooth_shape = np.vstack([ss - np.median(ss) for ss in Kappa_smooth])
    
#     # use birch which is super fast! 
#     from sklearn.cluster import Birch
# # >>> X = [[0, 1], [0.3, 1], [-0.3, 1], [0, -1], [0.3, -1], [-0.3, -1]]
#     brc = Birch(threshold=.9, branching_factor= 10000, n_clusters=None)
#     # brc.fit(Kappa_smooth_shape)
#     # # Birch(n_clusters=None)
#     # labels = brc.predict(Kappa_smooth_shape)
#     brc.fit(Kappa_smooth_shape[Kappa_groups<=1])
#     # Birch(n_clusters=None)
#     labels = brc.predict(Kappa_smooth_shape[Kappa_groups<=1])
    
#     print(len(np.unique(labels)))

#     labels_colors = sns.color_palette('Set1', n_colors = labels.max()+1)
#     # np.random.shuffle(labels_colors)    

#     labels_colors = np.vstack(labels_colors)
    
#     face_labels_colors = labels_colors[labels].copy()
    
    
#     all_face_labels_colors = .5*np.ones((len(Kappa_smooth),3))
#     all_face_labels_colors[Kappa_groups<=1] = face_labels_colors
    
#     mesh_birch = trimesh.Trimesh(vertices=v[0], 
#                                           faces=f,
#                                           face_colors = all_face_labels_colors[...,:3],
#                                           validate=False,
#                                           process=False)
#     mesh_birch.export('2022-02-11_Septin-Cell5-smooth_Kappa_subcluster.obj')
    
    
#     # from sklearn.cluster import SpectralClustering
#     # import numpy as np
#     # clustering = SpectralClustering(n_clusters=50,
#     #                                 assign_labels='discretize',
#     #                                 random_state=0).fit(Kappa_smooth_shape[Kappa_groups<=1])
    
#     # """
#     # need to do this by neighborhood graphs!. 
#     # """
    
#     # from sklearn.decomposition import PCA
    
#     # pca = PCA(n_components=2)
#     # pca.fit(Kappa_smooth[Kappa_groups==0])
#     # zz = pca.transform(Kappa_smooth[Kappa_groups==0])
    
#     # plt.figure()
#     # plt.plot(zz[:,0], 
#     #          zz[:,1], '.')
#     # plt.show()
    
    
    
    
# # # =============================================================================
# # #     1. correct the septin intensity bleaching. (ignoring the first timepoint which for some reason is different. )
# # #       based on experiments done on 2022-02-03, we go straight for division by the volumetric background!. 
# # # =============================================================================
    
# #     # correct and express as a relative series. 
# #     septin_faces_all_correct = septin_faces_all/(septin_vol_intensity_all[:,None] + 1e-12) # to prevent 0 division.
    
# # # =============================================================================
# # #     2. septin I signal decomposition into a baseline vs fluctuation component. 
# # # =============================================================================

# #     print('decomposing septin signals ... ')
# #     # detrend and decompose septin signals
# #     detrended_septin_faces_all_correct = [tsa.decompose_nonlinear_time_series(septin_faces_all_correct[:,iii], 
# #                                                                                 lam=1e4, p=.5, 
# #                                                                                 niter=10, 
# #                                                                                 padding=None) for iii in np.arange(septin_faces_all_correct.shape[1])[:5]]
    
# #     fluctuations_septin_I = np.vstack([ss[1] for ss in detrended_septin_faces_all_correct]) 
# #     baseline_septin_I = np.vstack([ss[0] for ss in detrended_septin_faces_all_correct])
    
# #     # plt.figure()
# #     # plt.plot(baseline_septin_I.T)
# #     # plt.plot(septin_faces_all_correct[:,:5])
# #     # plt.show()
    
# # # =============================================================================
# # #     Use the detrended septin levels to recolor the surfaces so they don't decay...... (in addition we can colour the extracted smoothness and fluctuation components...)
# # # =============================================================================
    

# #     # """
# #     # export the classification to all mesh timepoints
# #     #     - on blebs should be negative curvature lowest group 
# #     # """
# #     # savemeshfolder = os.path.join(analysis_savefolder, 'meshes_Septin_corrected'); fio.mkdir(savemeshfolder)
    
# #     # for vv_i, vv in enumerate(v): 
# #     #     # colors = septin_faces_all_correct[vv_i].copy() # over time. 
# #     #     face_colors = vol_colors.get_colors(septin_faces_all_correct[vv_i] * septin_vol_intensity_all[0], 
# #     #                                         colormap=cm.RdYlBu_r,
# #     #                                         vmin=I_min, 
# #     #                                         vmax=I_max)
# #     #     #### draw the coloring on the first position 
# #     #     mesh_septin = trimesh.Trimesh(vertices=vv, 
# #     #                                       faces=f,
# #     #                                       face_colors = np.uint8(255*face_colors[...,:3]),
# #     #                                       validate=False,
# #     #                                       process=False)
# #     #     mesh_septin.export(os.path.join(savemeshfolder, 'Frame_%s.obj' %(str(vv_i).zfill(3))))
        
        
# #     # savemeshfolder = os.path.join(analysis_savefolder, 'meshes_Septin_corrected_baseline'); fio.mkdir(savemeshfolder)
    
# #     # for vv_i, vv in enumerate(v): 
# #     #     # colors = septin_faces_all_correct[vv_i].copy() # over time. 
# #     #     face_colors = vol_colors.get_colors(baseline_septin_I[:,vv_i] * septin_vol_intensity_all[0], 
# #     #                                         colormap=cm.RdYlBu_r,
# #     #                                         vmin=I_min, 
# #     #                                         vmax=I_max)
# #     #     #### draw the coloring on the first position 
# #     #     mesh_septin = trimesh.Trimesh(vertices=vv, 
# #     #                                       faces=f,
# #     #                                       face_colors = np.uint8(255*face_colors[...,:3]),
# #     #                                       validate=False,
# #     #                                       process=False)
# #     #     mesh_septin.export(os.path.join(savemeshfolder, 'Frame_%s.obj' %(str(vv_i).zfill(3))))

# # # =============================================================================
# # #   3. Curvature analysis -> automatic thresholding based on mean curvature.     
# # # =============================================================================
# #     mean_kappa_sequences = np.nanmean(curvature_faces_all_um, axis=0)
# #     kappa_thresholds = skfilters.threshold_multiotsu(mean_kappa_sequences, 
# #                                                      classes=3)

# #     """
# #     apply the thresholds to get an idea.... 
# #         - on blebs should be negative curvature lowest group 
# #     """
# #     curvature_groups = np.zeros(len(mean_kappa_sequences), dtype=np.int)
# #     select = np.logical_and( mean_kappa_sequences > kappa_thresholds[0], mean_kappa_sequences <= kappa_thresholds[1])
# #     curvature_groups[select] = 1
# #     select = mean_kappa_sequences > kappa_thresholds[1]
# #     curvature_groups[select] = 2
    
    
# #     curvature_group_colors = sns.color_palette('coolwarm', 
# #                                                 len(kappa_thresholds)+1)
# #     curvature_group_colors = np.vstack(curvature_group_colors)
# #     curvature_group_colors = curvature_group_colors[::-1] # inversion to ensure on-blebs. 
    
# #     # map this to all faces!. 
# #     curvature_rgb = curvature_group_colors[curvature_groups]
    
    
# #     # """
# #     # export the classification to all mesh timepoints
# #     #     - on blebs should be negative curvature lowest group 
# #     # """
# #     # savemeshfolder = os.path.join(analysis_savefolder, 'meshes_curvature_3class'); fio.mkdir(savemeshfolder)
    
# #     # for vv_i, vv in enumerate(v): 
# #     #     #### draw the coloring on the first position 
# #     #     mesh_curvature = trimesh.Trimesh(vertices=vv, 
# #     #                                       faces=f,
# #     #                                       face_colors = np.uint8(255*curvature_rgb[...,:3]),
# #     #                                       validate=False,
# #     #                                       process=False)
# #     #     mesh_curvature.export(os.path.join(savemeshfolder, 'Frame_%s.obj' %(str(vv_i).zfill(3))))
        
        
        
# #     """
# #     Curvature signal smoothing so we can extract peaks more reliably.
# #     """    
# #     detrended_curvature_faces_all_um = [tsa.decompose_nonlinear_time_series(curvature_faces_all_um[:,iii], 
# #                                                                             lam=10, p=.5, 
# #                                                                             niter=10, 
# #                                                                             padding=None) for iii in np.arange(septin_faces_all_correct.shape[1])[:5]]

# #     fluctuations_curvature = np.vstack([ss[1] for ss in detrended_curvature_faces_all_um]) 
# #     baseline_curvature = np.vstack([ss[0] for ss in detrended_curvature_faces_all_um])
    
    
# #     all_peaks_curvature = []
    
# #     from scipy.signal import find_peaks
# #     for ind in tqdm(np.arange(baseline_curvature.shape[0])): # comput this based on the baseline 
        
# #         signal = -1*baseline_curvature[ind].copy() # invert due to the curvature convention. 
# #         peak_pos = find_peaks(signal, 
# #                               height=0, 
# #                               threshold=None, 
# #                               distance=None, 
# #                               prominence=None, 
# #                               width=None, 
# #                               wlen=None, 
# #                               rel_height=0.5, 
# #                               plateau_size=None) # this will need some work and potentially validation ..... 
        
# #         all_peaks_curvature.append(peak_pos[0]) # we just need peak position. 
# #     all_peaks_curvature = np.ascontiguousarray(all_peaks_curvature)
        

# # # # =============================================================================
# # # #     Save all the computed metadata. 
# # # # =============================================================================

# # #     savematfile = os.path.join(analysis_savefolder, 'curvature_septin_statistics_processed.mat')
# # #     spio.savemat(savematfile, 
# # #                  {'fluctuations_septin_I': fluctuations_septin_I, 
# # #                   'baseline_septin_I': baseline_septin_I,
# # #                   'septin_I_0': septin_vol_intensity_all[0], 
# # #                   'kappa_thresholds': kappa_thresholds,
# # #                   'curvature_groups': curvature_groups,
# # #                   'all_peaks_curvature': all_peaks_curvature,
# # #                   'fluctuations_curvature_um': fluctuations_curvature,
# # #                   'baseline_curvature_um': baseline_curvature, 
# # #                   'statsfile': statsfile})
    
    
# # # #         
    
# #     # def peak_detect(signal, p=.5, lam=10, niter=10): 
    
# # # #         from scipy.signal import find_peaks
        
# # # #         signal_ = decompose_nonlinear_time_series(y=signal, 
# # # #                                                   lam=lam, 
# # # #                                                   p=p, niter=niter, 
# # # #                                                   padding=None)[0]
    
# # # #         peak_pos = find_peaks(signal_, height=0, threshold=None, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)
# # # #         return peak_pos, signal_
    
    
# # # #     # test_ind = 0
    
# # # #     # yy = nonbleb_kappa[:,test_ind].copy()
# # # #     # yy_peaks, yy_ = peak_detect(yy, p=.5, lam=100, niter=10) 
# # # #     from tqdm import tqdm 
    
# # # #     all_peaks_curvature = []
    
# # # #     for test_ind in tqdm(np.arange(nonbleb_kappa.shape[1])):
        
# # # #         yy = nonbleb_kappa[:,test_ind].copy()
# # # #         yy_peaks, yy_ = peak_detect(yy, p=.5, lam=100, niter=10) 
        
# # # #         all_peaks_curvature.append(yy_peaks)
    
# # # #         # if len(yy_peaks[0]) > 0: 
    
# # # #         #     plt.figure()
# # # #         #     plt.plot(yy)
# # # #         #     plt.plot(yy_)
# # # #         #     plt.plot(yy_peaks[0], yy[yy_peaks[0]], 'ro')
# # # #         #     plt.show()
            
# # # #         # else:
            
# # # #         #     plt.figure()
# # # #         #     plt.plot(yy)
# # # #         #     plt.plot(yy_)
# # # #         #     # plt.plot(yy_peaks[0], yy[yy_peaks[0]], 'ro')
# # # #         #     plt.show()
    
    

# # # # # =============================================================================
# # # # #     Test 1: instantaneous pearson correlation of the time series -> no distinction 
# # # # # =============================================================================

# # # #     # nonbleb_select = mean_kappa_sequences >= thresholds[0]
# # # #     nonbleb_select = np.ones(len(mean_kappa_sequences), dtype=np.bool)
    
# # # #     # # now check 
# # # #     nonbleb_kappa = -curvature_faces_all_um[:,nonbleb_select] # so that positive = on bleb.... and we can look for positive correlation? (always looking for positive correlation) 
# # # #     nonbleb_I = septin_faces_all_correct[:, nonbleb_select]
    
    
# # # #     # # we may need to decompose first.... 
# # # #     # result = seasonal_decompose(np.vstack([np.arange(1,len(nonbleb_I)+1), 
# # # #     #                                        nonbleb_I[:,0]]).T, model='multiplicative', period=1)
     
# # # #     # # compute direct pearsoncorr # positive correlation ? 
# # # #     corr_I_kappa = np.hstack([spstats.pearsonr(nonbleb_kappa[:,iii], nonbleb_I[:,iii])[0] for iii in np.arange(nonbleb_kappa.shape[1])])  # result of this is heavily affected by the dip? 
    
    
# # # #     corr_I_kappa_mesh_color = get_colors(corr_I_kappa, colormap=cm.coolwarm, vmin=-1, vmax=1) 
    
# # # #     mesh_corr_I_kappa = trimesh.Trimesh(vertices=v[0], 
# # # #                                           faces=f,
# # # #                                           face_colors = corr_I_kappa_mesh_color[...,:3],
# # # #                                           validate=False,
# # # #                                           process=False)
# # # #     mesh_corr_I_kappa.export('2022-02-02_Septin-Cell5-exemplar_curvature_correlation-IvsKappa.obj')
    
    
# # # # # =============================================================================
# # # # #     Test 2: detrending -> obtaining a baseline to the signal -> does this change the observed correlations. 
# # # # # =============================================================================

# # # #     detrended_nonbleb_I = np.vstack([decompose_nonlinear_time_series(nonbleb_I[:,iii], 
# # # #                                               lam=1e5, p=.5, 
# # # #                                               niter=10, 
# # # #                                               padding=None)[1] for iii in np.arange(nonbleb_kappa.shape[1])])
# # # #     detrended_nonbleb_I = detrended_nonbleb_I.T 
# # # #     corr_I_kappa_detrended = np.hstack([spstats.pearsonr(nonbleb_kappa[:,iii], detrended_nonbleb_I[:,iii])[0] for iii in np.arange(nonbleb_kappa.shape[1])])  # result of this is heavily affected by the dip? 
    
    
# # # #     corr_I_kappa_detrended_mesh_color = get_colors(corr_I_kappa_detrended, colormap=cm.coolwarm, vmin=-1, vmax=1) 
    
# # # #     mesh_corr_I_kappa_detrended = trimesh.Trimesh(vertices=v[0], 
# # # #                                                   faces=f,
# # # #                                                   face_colors = corr_I_kappa_detrended_mesh_color[...,:3],
# # # #                                                   validate=False,
# # # #                                                   process=False)
# # # #     mesh_corr_I_kappa_detrended.export('2022-02-02_Septin-Cell5-exemplar_curvature_correlation-I_detrend-vs-Kappa.obj')
    
    
# # # #     # double check with the diff signal 
# # # #     corr_I_kappa_diff = np.hstack([spstats.pearsonr(nonbleb_kappa[:-1,iii], np.diff(nonbleb_I[:,iii]))[0] for iii in np.arange(nonbleb_kappa.shape[1])])  # result of this is heavily affected by the dip? 
# # # #     corr_I_kappa_diff_mesh_color = get_colors(-1*corr_I_kappa_diff, colormap=cm.coolwarm, vmin=-.3, vmax=.3) 
    
# # # #     mesh_corr_I_kappa_diff = trimesh.Trimesh(vertices=v[0], 
# # # #                                                   faces=f,
# # # #                                                   face_colors = corr_I_kappa_diff_mesh_color[...,:3],
# # # #                                                   validate=False,
# # # #                                                   process=False)
# # # #     mesh_corr_I_kappa_diff.export('2022-02-02_Septin-Cell5-exemplar_curvature_correlation-I_diff-vs-Kappa.obj')
    
# # # # # =============================================================================
# # # # #     Test 3: analysing the mean accumulation rate. 
# # # # # =============================================================================

# # # #     accum_signal = nonbleb_I - detrended_nonbleb_I
    
# # # #     mean_accum_signals = np.nanmean(np.diff(accum_signal, axis=0), axis=0) # might not be best... -> fitting a linear trend line instead.... 
# # # #     mean_accum_signals_rate = np.hstack([spstats.linregress(np.arange(len(sig)), sig)[0] for sig in accum_signal.T])
# # # #     mean_accum_mesh_color = get_colors(mean_accum_signals_rate, colormap=cm.coolwarm, vmin=-2*np.std(mean_accum_signals), vmax=2*np.std(mean_accum_signals)) 
    
# # # #     mesh_mean_accum = trimesh.Trimesh(vertices=v[0], 
# # # #                                                   faces=f,
# # # #                                                   face_colors = mean_accum_mesh_color[...,:3],
# # # #                                                   validate=False,
# # # #                                                   process=False)
# # # #     mesh_mean_accum.export('2022-02-02_Septin-Cell5-exemplar_curvature_SeptinMeanAccum.obj')
    
    
# # # #     # assessing using instead the AUC
# # # #     from numpy import trapz
# # # #     auc_accum_signal = np.hstack([trapz(sig[:]-sig[0], dx=1) for sig in accum_signal.T])
# # # #     mean_accum_mesh_color = get_colors(auc_accum_signal, colormap=cm.coolwarm, vmin=-2*np.std(auc_accum_signal), vmax=2*np.std(auc_accum_signal)) 
# # # #     mesh_mean_accum = trimesh.Trimesh(vertices=v[0], 
# # # #                                                   faces=f,
# # # #                                                   face_colors = mean_accum_mesh_color[...,:3],
# # # #                                                   validate=False,
# # # #                                                   process=False)
# # # #     mesh_mean_accum.export('2022-02-02_Septin-Cell5-exemplar_curvature_SeptinMeanAUC.obj')
    
    
# # # # # =============================================================================
# # # # #     Test 4: extracting the peak statistics in the curvature signals for plotting and investigation.  
# # # # # =============================================================================


# # # # # # # now check 
# # # # #     nonbleb_kappa = -curvature_faces_all_um[:,nonbleb_select] # so that positive = on bleb.... and we can look for positive correlation? (always looking for positive correlation) 
# # # # #     nonbleb_I = septin_faces_all_correct[:, nonbleb_select]

# # # #     def peak_detect(signal, p=.5, lam=10, niter=10): 
    
# # # #         from scipy.signal import find_peaks
        
# # # #         signal_ = decompose_nonlinear_time_series(y=signal, 
# # # #                                                   lam=lam, 
# # # #                                                   p=p, niter=niter, 
# # # #                                                   padding=None)[0]
    
# # # #         peak_pos = find_peaks(signal_, height=0, threshold=None, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)
# # # #         return peak_pos, signal_
    
    
# # # #     # test_ind = 0
    
# # # #     # yy = nonbleb_kappa[:,test_ind].copy()
# # # #     # yy_peaks, yy_ = peak_detect(yy, p=.5, lam=100, niter=10) 
# # # #     from tqdm import tqdm 
    
# # # #     all_peaks_curvature = []
    
# # # #     for test_ind in tqdm(np.arange(nonbleb_kappa.shape[1])):
        
# # # #         yy = nonbleb_kappa[:,test_ind].copy()
# # # #         yy_peaks, yy_ = peak_detect(yy, p=.5, lam=100, niter=10) 
        
# # # #         all_peaks_curvature.append(yy_peaks)
    
# # # #         # if len(yy_peaks[0]) > 0: 
    
# # # #         #     plt.figure()
# # # #         #     plt.plot(yy)
# # # #         #     plt.plot(yy_)
# # # #         #     plt.plot(yy_peaks[0], yy[yy_peaks[0]], 'ro')
# # # #         #     plt.show()
            
# # # #         # else:
            
# # # #         #     plt.figure()
# # # #         #     plt.plot(yy)
# # # #         #     plt.plot(yy_)
# # # #         #     # plt.plot(yy_peaks[0], yy[yy_peaks[0]], 'ro')
# # # #         #     plt.show()
    
# # # #     # detrend = decompose_nonlinear_time_series(nonbleb_I[:,test_ind], 
# # # #     #                                           lam=1e5, p=.5, 
# # # #     #                                           niter=10, 
# # # #     #                                           padding=None)
    
# # # #     # plt.figure()
# # # #     # # plt.plot(nonbleb_I[:,test_ind],'k-')
# # # #     # # plt.plot(detrend[0], 'r-')
# # # #     # plt.plot(detrend[1], 'g-')
# # # #     # plt.show()
    
    
# # # #     """
# # # #     Compute statistics of the peaks. 
# # # #     """
# # # #     num_peaks = np.hstack([len(pp[0]) for pp in all_peaks_curvature])
# # # #     diff_peaks = []
# # # #     for pp in all_peaks_curvature:
# # # #         if len(pp[0]) > 1: 
# # # #             diff_peaks.append(np.mean(np.diff(pp[0])))
# # # #         else:
# # # #             diff_peaks.append(1000) # crazy large number representing never oscillate.! 
            
# # # #     diff_peaks = np.hstack(diff_peaks)
        
# # # #     amplitude_peaks = []
# # # #     for pp_ii, pp in enumerate(all_peaks_curvature):
# # # #         if len(pp[0]) > 1: 
# # # #             amplitude_peaks.append(np.nanmean(nonbleb_kappa[pp[0],pp_ii]))
# # # #         else:
# # # #             amplitude_peaks.append(0) # crazy large number representing never oscillate.! 
            
# # # #     amplitude_peaks = np.hstack(amplitude_peaks)


# # # # # =============================================================================
# # # # #   export these in the meshes and have a study.... 
# # # # # =============================================================================

# # # #     n_peaks_mesh_color = get_colors(num_peaks, 
# # # #                                     colormap=cm.coolwarm, 
# # # #                                     vmin=np.min(num_peaks), 
# # # #                                     vmax=np.max(num_peaks)) 
# # # #     mesh_n_peaks = trimesh.Trimesh(vertices=v[0], 
# # # #                                         faces=f,
# # # #                                         face_colors = n_peaks_mesh_color[...,:3],
# # # #                                         validate=False,
# # # #                                         process=False)
# # # #     mesh_n_peaks.export('2022-02-02_Septin-Cell5-exemplar_num-peaks-curvature.obj')
    
    
# # # #     freq_peaks_mesh_color = get_colors(diff_peaks, 
# # # #                                         colormap=cm.coolwarm_r, 
# # # #                                         vmin=10, 
# # # #                                         vmax=40) 
# # # #     mesh_freq_peaks = trimesh.Trimesh(vertices=v[0], 
# # # #                                         faces=f,
# # # #                                         face_colors = freq_peaks_mesh_color[...,:3],
# # # #                                         validate=False,
# # # #                                         process=False)
# # # #     mesh_freq_peaks.export('2022-02-02_Septin-Cell5-exemplar_freq-peaks-curvature.obj')
    
    
# # # #     amplitude_peaks_colors = get_colors(amplitude_peaks, 
# # # #                                         colormap=cm.coolwarm, 
# # # #                                         vmin=-.5, 
# # # #                                         vmax=.5) 
# # # #     mesh_amplitude_peaks = trimesh.Trimesh(vertices=v[0], 
# # # #                                         faces=f,
# # # #                                         face_colors = amplitude_peaks_colors[...,:3],
# # # #                                         validate=False,
# # # #                                         process=False)
# # # #     mesh_amplitude_peaks.export('2022-02-02_Septin-Cell5-exemplar_amplitude-peaks-curvature.obj')
    
# # # # # =============================================================================
# # # # #     Isolate the curvatures and have another look? 
# # # # # =============================================================================
    
    
# # # # =============================================================================
# # # #     2. Test the voxelisation -> work on this a bit more with the subdivision idea.. ..... 
# # # # =============================================================================
# # #     mesh = trimesh.Trimesh(vertices = v[np.argmin(septin_vol_intensity_all)], 
# # #                             faces =  f, 
# # #                             process=False, 
# # #                             validate=False)
    
# # #     import igl 
# # #     vv,ff = igl.upsample(mesh.vertices,mesh.faces, 1)
    
# # #     # """
# # #     # Check the longest edge at all time points!. 
# # #     # """
# # #     # longest_edges = []
    
# # #     # for nn in np.arange(len(v)): 
# # #     #     mesh = trimesh.Trimesh(vertices = v[nn], 
# # #     #                             faces =  f, 
# # #     #                             process=False, 
# # #     #                             validate=False)
        
# # #     #     longest_edge = np.linalg.norm(mesh.vertices[mesh.edges[:, 0]] -
# # #     #                                 mesh.vertices[mesh.edges[:, 1]],
# # #     #                                 axis=1).min()
        
# # #     #     longest_edges.append(longest_edge)
    
# # #     # plt.figure()
# # #     # plt.plot(longest_edges)
# # #     # plt.plot(np.argmin(septin_vol_intensity_all), longest_edges[np.argmin(septin_vol_intensity_all)], 'ro')
# # #     # plt.show()
    
    
# # #     def voxelize_image_mesh_pts(mesh, pad=50, dilate_ksize=3, erode_ksize=1, vol_shape=None, upsample_iters_max=10, pitch=2):

# # #         # works only for meshs from images. so all coordinates are positive. 
# # #         # this voxelizes the pts without need for an image. 
# # #         import numpy as np 
# # #         import skimage.morphology as skmorph
# # #         from scipy.ndimage.morphology import binary_fill_holes
# # #         import igl 
        
# # #         vv = mesh.vertices.copy()
# # #         ff = mesh.faces.copy()
    
# # #         if vol_shape is None:
# # #             # mesh_pts = mesh.vertices.copy() + 1
# # #             longest_edge_length = igl.edge_lengths(vv,ff).max()
# # #             factor = longest_edge_length / float(dilate_ksize) / 2.
# # #             print(factor)
# # #             if factor >= pitch / 2. :
# # #                 print('upsample')
# # #                 # then we can't get a volume even if watertight. 
# # #                 upsample_iters = int(np.rint(np.log2(factor)))
# # #                 print(upsample_iters)
# # #                 upsample_iters = np.min([upsample_iters, upsample_iters_max])
# # #                 vv, ff = igl.upsample(mesh.vertices, mesh.faces, upsample_iters)
# # #             mesh_pts = igl.barycenter(vv,ff) + 1
# # #             # determine the boundaries. 
# # #             min_x, min_y, min_z = np.min(mesh_pts, axis=0)
# # #             max_x, max_y, max_z = np.max(mesh_pts, axis=0)
    
# # #             # pad = int(np.min([min_x, min_y, min_z])) # auto determine the padding based on this. 
# # #             # new binary. 
# # #             smooth_img_binary = np.zeros((int(max_x)+pad, int(max_y)+pad, int(max_z)+pad))
# # #         else:
# # #             # mesh_pts = mesh.vertices.copy() #+ .5
# # #             # mesh_pts = mesh.vertices.copy() + 1
# # #             longest_edge_length = igl.edge_lengths(vv,ff).max()
# # #             factor = longest_edge_length / float(dilate_ksize) / 2.
# # #             if factor >= pitch / 2. :
# # #                 # then we can't get a volume even if watertight. 
# # #                 upsample_iters = int(np.rint(np.log2(factor)))
# # #                 upsample_iters = np.min([upsample_iters, upsample_iters_max])
# # #                 vv, ff = igl.upsample(mesh.vertices, mesh.faces, upsample_iters)
# # #             mesh_pts = igl.barycenter(vv,ff)
# # #             smooth_img_binary = np.zeros(vol_shape)
    
# # #         smooth_img_binary[mesh_pts[:,0].astype(np.int), 
# # #                           mesh_pts[:,1].astype(np.int), 
# # #                           mesh_pts[:,2].astype(np.int)] = 1
# # #         smooth_img_binary = skmorph.binary_dilation(smooth_img_binary, skmorph.ball(dilate_ksize))
# # #         smooth_img_binary = binary_fill_holes(smooth_img_binary) # since we dilated before to create a full mesh. we inturn must erode. 
        
# # #         if erode_ksize is not None:
# # #             smooth_img_binary = skmorph.binary_erosion(smooth_img_binary, skmorph.ball(erode_ksize))
    
# # #         return smooth_img_binary 
    
    
# # #     binary = voxelize_image_mesh_pts(mesh, 
# # #                                      pad=50, 
# # #                                      dilate_ksize=3, 
# # #                                      erode_ksize=1, 
# # #                                      vol_shape=None, 
# # #                                      upsample_iters_max=10, pitch=2)
    
# # #     # total_area = []
# # #     # for nn in np.arange(len(v)):
# # #     #     mesh0 = trimesh.Trimesh(vertices = v[nn], 
# # #     #                                 faces =  f, 
# # #     #                                 process=False, 
# # #     #                                 validate=False)
# # #     #     binary =  voxelize_image_mesh_pts(mesh0, pad=50, dilate_ksize=3, erode_ksize=1, vol_shape=None, upsample_iters=0) # so there are holes! 
# # #     #     total_area.append(np.sum(binary))
    
    
    
    
    
    
    
# # #     # def voxelize_image_mesh_pts(mesh, pad=50, dilate_ksize=3, erode_ksize=1, vol_shape=None):

# # #     # # works only for meshs from images. so all coordinates are positive. 
# # #     # # this voxelizes the pts without need for an image. 
# # #     # import numpy as np 
# # #     # import skimage.morphology as skmorph
# # #     # from scipy.ndimage.morphology import binary_fill_holes
# # #     # import igl 

# # #     # if vol_shape is None:
# # #     #     # mesh_pts = mesh.vertices.copy() + 1
# # #     #     mesh_pts = igl.barycenter(mesh.vertices, mesh.faces) + 1
# # #     #     # determine the boundaries. 
# # #     #     min_x, min_y, min_z = np.min(mesh_pts, axis=0)
# # #     #     max_x, max_y, max_z = np.max(mesh_pts, axis=0)

# # #     #     # pad = int(np.min([min_x, min_y, min_z])) # auto determine the padding based on this. 
# # #     #     # new binary. 
# # #     #     smooth_img_binary = np.zeros((int(max_x)+pad, int(max_y)+pad, int(max_z)+pad))
# # #     # else:
# # #     #     # mesh_pts = mesh.vertices.copy() #+ .5
# # #     #     mesh_pts = igl.barycenter(mesh.vertices, mesh.faces)
# # #     #     smooth_img_binary = np.zeros(vol_shape)

# # #     # smooth_img_binary[mesh_pts[:,0].astype(np.int), 
# # #     #                   mesh_pts[:,1].astype(np.int), 
# # #     #                   mesh_pts[:,2].astype(np.int)] = 1
# # #     # smooth_img_binary = skmorph.binary_dilation(smooth_img_binary, skmorph.ball(dilate_ksize))
# # #     # smooth_img_binary = binary_fill_holes(smooth_img_binary) # since we dilated before to create a full mesh. we inturn must erode. 
    
# # #     # if erode_ksize is not None:
# # #     #     smooth_img_binary = skmorph.binary_erosion(smooth_img_binary, skmorph.ball(erode_ksize))

# # #     # return smooth_img_binary 
    
# # #     # voxel_mesh = mesh.voxelized(pitch=.5, method='subdivide', max_iter=100)
# # #     # # trimesh.voxel.creation.voxelize(mesh, pitch, method='subdivide', **kwargs)
    
# # #     # voxel_mesh_image = voxel_mesh.matrix.copy()
# # #     # voxel_mesh_image = skmorph.binary_closing(voxel_mesh_image, skmorph.ball(3))
    
# # #     # # matrix_surface = sparse_to_matrix(mesh_discretized.sparse_surface)

# # #     # # # fill the volume
# # #     # # matrix_full = ndimage.morphology.binary_fill_holes(matrix_surface)
    

    










# # # # # =============================================================================
# # # # #     Plot random curvature traces.
# # # # # =============================================================================
# # # #     plt.figure(figsize=(5,15))
# # # #     for ii in np.arange(100):
# # # #         plt.plot(np.arange(200), 
# # # #                  curvature_faces_all_um[:,ii]+ii*.5)
# # # #     plt.gca().invert_yaxis()
# # # #     plt.xlim([0,200])
# # # #     # plt.ylim([2,-2])
# # # #     plt.show()
    
# # # # # =============================================================================
# # # # #     Plot corresponding corrected septin traces.
# # # # # =============================================================================
# # # #     mean_volume = septin_vol_intensity_all * exp_correct
    
# # # #     plt.figure(figsize=(5,5))
# # # #     plt.plot(septin_vol_intensity_all, 'k-', lw=2, label='raw')
# # # #     plt.plot(xline, yline,'r--', lw = 2, label='linear-fit')
# # # #     plt.plot(xline, exp_params[1],'g--', lw = 2, label='exp-fit')
# # # #     plt.plot(mean_volume, 'b-', lw=2, label='corrected')
# # # #     plt.ylim([300,400])
# # # #     plt.xlim([0,200])
# # # #     plt.legend()
# # # #     plt.show()
    
# # # #     septin_faces_all_correct = septin_faces_all * exp_correct[:,None] / np.nanmean(mean_volume[:,None])
    
# # # #     plt.figure(figsize=(5,30))
# # # #     for ii in np.arange(100):
# # # #         plt.plot(np.arange(200), 
# # # #                  septin_faces_all_correct[:,ii]+ii*1)
# # # #     plt.gca().invert_yaxis()
# # # #     # plt.xlim([0,200])
# # # #     # plt.ylim([2,-2])
# # # #     plt.show()
    
    
# # # #     for ii in np.arange(30):
# # # #         fig, ax1 = plt.subplots(figsize=(5,2))
# # # #         plt.title(ii)
# # # #         ax1.plot(curvature_faces_all_um[:,ii], 'r'); plt.gca().invert_yaxis()
# # # #         ax2 = ax1.twinx()
# # # #         ax2.plot(septin_faces_all_correct[:,ii], 'k' )
# # # #         plt.show()
    
    
# # # #     # give example curvature time series?
# # # #     # mean curvature of the series (sort by)
# # # #     mean_kappa = np.nanmean(curvature_faces_all_um, axis=0)
# # # #     sort_order = np.argsort(mean_kappa)
    
# # # #     fig, ax = plt.subplots()
# # # #     ax.imshow(curvature_faces_all_um[:,sort_order][:,:100].T, cmap='coolwarm_r', vmin=-1, vmax=1)
# # # #     ax.set_aspect('auto')
# # # #     plt.show()
    
# # # #     # plot example
# # # #     plt.figure(figsize=(5,5))
# # # #     plt.plot(curvature_faces_all_um[:,sort_order][:,:1])
# # # #     plt.gca().invert_yaxis()
# # # #     plt.xlim([0,200])
# # # #     plt.ylim([2,-2])
# # # #     plt.show()
    
    
# # # #     plt.figure(figsize=(5,5))
# # # #     plt.plot(curvature_faces_all_um[:,:100])
# # # #     plt.gca().invert_yaxis()
# # # #     plt.xlim([0,200])
# # # #     plt.ylim([2,-2])
# # # #     plt.show()
    
    
# # # #     # example of the corresponding septin trace
    
# # # #     fig, ax = plt.subplots()
# # # #     ax.imshow(septin_faces_all[:,sort_order][:,:100].T, cmap='coolwarm')
# # # #     ax.set_aspect('auto')
# # # #     plt.show()
    
    
# # # #     plt.figure(figsize=(5,5))
# # # #     plt.plot(septin_faces_all[:,sort_order][:,:1])
# # # #     # plt.gca().invert_yaxis()
# # # #     plt.xlim([0,200])
# # # #     # plt.ylim([2,-2])
# # # #     plt.show()
    
# # # # # =============================================================================
# # # # #     save out all the relevant statistics. 
# # # # # =============================================================================

# # # #     import igl 
    
# # # #     sss = np.dstack([stats['septin_all'], 
# # # #                      stats['septin_all'] * exp_correct[:,None],
# # # #                      curvature_faces_all_um])
    
# # # #     sss_vertices = np.array([igl.average_onto_vertices(stats['v_all'][iii], 
# # # #                                                             mesh.faces, 
# # # #                                                             sss[iii]) for iii in np.arange(len(stats['v_all']))])
    

# # # #     # savedict = {'v_all' : stats['v_all'],
# # # #     #             'f': mesh.faces + 1, 
# # # #     #             'I_septin_raw_faces' : stats['septin_all'], 
# # # #     #             'I_septin_correct_faces': stats['septin_all'] * exp_correct[:,None],
# # # #     #             'curvature_um_faces' : curvature_faces_all_um}
    
# # # #     # spio.savemat('Septin-Curvature_statistics-MatlabFormat.mat', 
# # # #     #              savedict)
    
    
# # # #     savedict = {'v_all' : stats['v_all'],
# # # #                 'f': mesh.faces + 1, 
# # # #                 'I_septin_raw_vertex' : sss_vertices[...,0], 
# # # #                 'I_septin_correct_vertex': sss_vertices[...,1],
# # # #                 'curvature_um_vertex' : sss_vertices[...,2]}
    
# # # #     spio.savemat('Septin-Curvature_vertexbased_statistics-MatlabFormat.mat', 
# # # #                  savedict)
    




    