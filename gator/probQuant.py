# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 16:54:39 2022
@author: Ajit Johnson Nirmal
Function to calculate the mean/median probability 
"""

# lib
import tifffile
import pandas as pd
from skimage import measure
import ast
import numpy as np
import pathlib
import os

a = tifffile.imread(r"\\files.med.harvard.edu\HITS\lsp-analysis\gator\probability\6_Preview_Ecad.tif")[0]
b = tifffile.imread(r"\\files.med.harvard.edu\HITS\lsp-analysis\gator\probability\6_Preview_CD45_combined.tif")[0]
c = tifffile.imread(r"\\files.med.harvard.edu\HITS\lsp-analysis\gator\probability\6_Preview_CD4.tif")[0]
d = tifffile.imread(r"\\files.med.harvard.edu\HITS\lsp-analysis\gator\probability\6_Preview_CD3D.tif")[0]
e = tifffile.imread(r"\\files.med.harvard.edu\HITS\lsp-analysis\gator\probability\6_Preview_CD8A_combined.tif")[0]
f = tifffile.imread(r"\\files.med.harvard.edu\HITS\lsp-analysis\gator\probability\6_Preview_Ki67.tif")[0]

con = np.stack((a, b, c, d, e, f),axis=2)
con = np.stack((a, b, c, d, e, f),axis=0)
tifffile.imwrite(r"\\files.med.harvard.edu\HITS\lsp-analysis\gator\probability\6.tif", data=con,  metadata={'Channel': {'Name': ['ECAD','CD45','CD4','CD3D','CD8A', 'KI67'] }})

probabilityMask = r"\\files.med.harvard.edu\HITS\lsp-analysis\gator\probability\6.tif"
segmentationMask = r"C:\Users\ajn16\Dropbox (Partners HealthCare)\Data\gator\data\Exemplar\modified_dearray\6\segmentation\unmicst-6\cellMask.tif"

outputDir = 'C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/gator/data/ajn_training_data'

probQuant (probabilityMask, segmentationMask,  feature='median', outputDir=outputDir)

# function
def probQuant (probabilityMask, segmentationMask, feature='median', outputDir=None):
    
    # read the seg mask
    segM = tifffile.imread(segmentationMask)
    probM = tifffile.imread(probabilityMask)[0]
    
    if len(probM.shape) > 2:
        probM = np.moveaxis(probM, 0, -1)
    
    def median_intensity(mask, img):
        return np.median(img[mask])

    # quantify
    print("Quantifying the probability masks")
    quantTable = pd.DataFrame(measure.regionprops_table(segM, intensity_image=probM, 
                                                        properties=['label','mean_intensity'],
                                                        extra_properties=[median_intensity])).set_index('label')
    
    # keep only median
    if feature == 'median':
        quantTable = quantTable.filter(regex='median')
    if feature == 'mean':
        quantTable = quantTable.filter(regex='mean')
    
    # read the channel names from the tiffile
    tiff = tifffile.TiffFile(probabilityMask)
    omexml_string = ast.literal_eval(tiff.pages[0].description)
    channel_names = omexml_string['Channel']['Name']
    
    # assign channel names
    quantTable.columns = channel_names
    quantTable = quantTable / quantTable.max()
    
    # if outputDir is given
    if outputDir is None:
        outputDir = os.getcwd()
    
    # final path to save results
    finalPath = pathlib.Path(outputDir + '/GATOR/probQuant/')
    if not os.path.exists(finalPath):
        os.makedirs(finalPath)
        
    # file name 
    file_name = pathlib.Path(probabilityMask).stem + '.csv'
    quantTable.to_csv(finalPath / file_name)
        



