#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###########################################################################
# Script: Compute spatial overlap between fMRI activation and network masks
#
# Author: Gerasimos Gerardos
# Project: MeMoSLAP Project (P04–Leipzig)
# Paper: Task load modulates network interactions between bilateral fronto-
#        parietal and cerebellar areas during verbal working memory
#
# Description:
# This script computes spatial similarity between thresholded fMRI activation
# maps and predefined brain network masks using voxel-wise overlap metrics.
# It supports Dice coefficient, overlap coefficient, and percentage overlap
# measures, and is used to quantify how task-related activations intersect
# with large-scale networks (e.g., MDN, DMN).
#
# Workflow:
# - Loads NIfTI mask files using nibabel.
# - Binarizes all masks (voxels > 0).
# - Computes voxel-wise intersections between masks.
# - Calculates:
#     * Dice index (2|A∩B| / (|A| + |B|))
#     * Overlap coefficient (|A∩B| / min(|A|,|B|))
#     * Percentage overlap relative to each mask
# - Prints summary statistics for positive and negative task contrasts
#   relative to MDN and DMN templates.
#
# Inputs:
# - Thresholded task activation masks (e.g., mask123vs0_*.nii)
# - Network masks (e.g., MDN, DMN, Yeo networks)
#
# Outputs:
# - Console output reporting Dice scores, overlap coefficients,
#   and percentage overlap between mask pairs.
#
# Notes:
# - All input masks must be in the same space and resolution (e.g., MNI).
# - Masks should be coregistered and aligned prior to analysis.
# - Update file paths below to match your directory structure.
###########################################################################


import nibabel as nib
import numpy as np

def overlap_coefficient(mask1_path, mask2_path, index = 'dice'):
    # Load the masks
    mask1 = nib.load(mask1_path).get_fdata()
    mask2 = nib.load(mask2_path).get_fdata()

    # Ensure masks are binary
    mask1 = mask1 > 0
    mask2 = mask2 > 0

    # Flatten to 1D arrays
    mask1_flat = mask1.flatten()
    mask2_flat = mask2.flatten()
    
    if index == 'dice':
        # Calculate Dice
        intersection = np.sum(mask1_flat & mask2_flat)
        volume_sum = np.sum(mask1_flat) + np.sum(mask2_flat)
    
        if volume_sum == 0:
            return 1.0  # Both masks are empty
        return 2.0 * intersection / volume_sum
    elif index == 'overlap':

        #    Calculate Overlap Coefficient
        intersection = np.sum(mask1_flat & mask2_flat)
        min_size = min(np.sum(mask1_flat), np.sum(mask2_flat))
    
        if min_size == 0:
            return 0.0  # Avoid division by zero (if one mask is empty)
    
        return intersection / min_size
    
    elif index == 'percentages':
        intersection = np.sum(mask1_flat & mask2_flat)
        percentage = intersection = np.sum(mask1_flat & mask2_flat)/np.sum(mask1_flat)
        return percentage*100 , np.sum(mask1_flat), np.sum(mask2_flat)

mdn = "/results/fMRI_analysis/masks/r_brain_combine_123vs0_spmT0001MDN_mask.nii"
dmn = "/results/fMRI_analysis/masks/r_brain_combine_123vs0_spmT0001Yeo7Networks_Default.nii"
pos123vs0 = "/results/fMRI_analysis/masks/mask123vs0_p001FWEc77.nii"
neg123vs0 = "/results/fMRI_analysis/masks/mask0vs123p001FWEc55.nii"

over_pos_mdn = overlap_coefficient(pos123vs0, mdn, 'overlap')
over_neg_dmn = overlap_coefficient(neg123vs0, dmn, 'overlap')
print(f'Overlap coefficient of the MDN and the positive task activation: {over_pos_mdn:.3f}')
print(f'Overlap coefficient of the DMN and the negative task activation: {over_neg_dmn:.3f}')


dice_pos_mdn = overlap_coefficient(pos123vs0, mdn, 'dice')
dice_neg_dmn = overlap_coefficient(neg123vs0, dmn, 'dice')
print(f'Dice Index of the MDN and the positive task activation: {dice_pos_mdn:.3f}')
print(f'Dice Index of the DMN and the negative task activation: {dice_neg_dmn:.3f}')


o, a, b = overlap_coefficient(mdn, pos123vs0,'percentages')
print(f'{o:.1f} % of voxels of the MDN ({a} voxels) overlap with the task activation ({b} voxels)')
o, a, b = overlap_coefficient(pos123vs0, mdn,'percentages')
print(f'{o:.1f} % of voxels of the task activation ({a} voxels) overlap with the MDN ({b} voxels)')

o, a, b = overlap_coefficient(dmn, neg123vs0,'percentages')
print(f'{o:.1f} % of voxels of the DMN ({a} voxels) overlap with the negative task activation ({b} voxels)')
o, a, b = overlap_coefficient(neg123vs0, dmn,'percentages')
print(f'{o:.1f} % of voxels of the negative task activation ({a} voxels) overlap with the DMN ({b} voxels)')