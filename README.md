# GRBF-For-Image-Alignment
SSD cost function smoothing using gaussian kernels with GRBF image representation
This repository contains MATLAB scripts for image alignment using translation, scaling, and affine transformations.
Translate, Sum_GRBF and GRBF_rep cannot always be used interchangably between codes as some have different approaches.
use functions in the same folder.

## Folders
### 1. `Images`
Contains example input images used in the code.

### 2. `Translation only`
Scripts for image registration using translation-only transformations.

### 3. `Translation with scaling`
Scripts for image registration with combined translation and scaling transformations.

### 4. `Affine`
Scripts for affine transformations in the x and y direction.

## Problematic Scripts:
### 1. 'affine/test_img_result_of_affine_sum_GRBF.m'
- **Problem**: Affine transformations work well for small combinations with scaling but not with translation.
- **Suspected Issue**: The tau formula is incorrect or possily the application of the translate function.
