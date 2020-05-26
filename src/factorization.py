"""File To Handle Factorization Routines in Python"""
import numpy as np
import pandas as pd
import tensorly as tl
from tensorly.decomposition import parafac
from tensorly.regression.metrics import variance as tl_var

############################################### Data Import #######################################################
def get_data():
    "Get Data Set into Numpy Array"
    ## Deal with odd HDF5 Stuff
    f = h5py.File("serology_tensor.jld", "r")
    tens = f['tensor'][:]
    
    ## Data Carriers
    dim1, dim2, dim3 = np.shape(tens)
    output = np.zeros(np.shape(tens))
    
    ## Iteration
    for i in np.arange(dim1):
        for j in np.arange(dim2):
            for k in np.arange(dim3):
                test = f[tens[i, j, k]][()]
                if isinstance(test, np.void):
                    output[i, j, k] = np.nan
                else: 
                    output[i, j, k] = test
    return output

############################################## Tensor Handling #####################################################
def perform_parafac(tens, r, mk):
    """Perform Parafac Decomposition
    ------------------------------------------------
    Input:
        tens: 3D Data Tensor
        r: rank of decomposition
        mk: mask for missing values
            - False (0) where data is missing
            - True (1) where data is present
    Output:
        output[0]: weights
        output[1]: factor matrices 
    
    """
    return parafac(tens, r, mask=mk)

def getMask(tens):
    """Returns appropriate mask and modified tensor for decomposition"""
    actual_mask = ~np.isnan(tens)
    missing_indices = np.isnan(tens)
    tens[missing_indices] = 0
    return tens, actual_mask


def parafac_R2X(output, original):
    """Calculate Reconstruction error from PARAFAC Factors"""
    reconstructed = tl.kruskal_to_tensor(output)
    return R2X(reconstructed, original)

def R2X(reconstructed, original):
    """ Calculates R2X of two tensors. """
    return 1.0 - tl_var(reconstructed - original) / tl_var(original)     ## To Do: Fix for missing data problem
