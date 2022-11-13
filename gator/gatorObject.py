# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 09:53:28 2022
@author: Ajit Johnson Nirmal
Function to incorporate the feature matrix and probability matrix into an adata object
"""

# libs
import anndata as ad
import pandas as pd
import pathlib
import numpy as np
import os



spatialTablePath = r"C:\Users\ajn16\Dropbox (Partners HealthCare)\Data\gator\data\Exemplar\modified_dearray\quantification\unmicst-6_cellMask.csv"
probTablePath = r"C:\Users\ajn16\Dropbox (Partners HealthCare)\Data\gator\data\ajn_training_data\GATOR\probQuant\6.csv"

outputDir = 'C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/gator/data/ajn_training_data'
gatorObject (spatialTablePath=spatialTablePath, probTablePath=probTablePath, outputDir=outputDir)

# function
def gatorObject (spatialTablePath, probTablePath,
                 CellId='CellID',unique_CellId=True,
                 split='X_centroid',
                 remove_dna=True, remove_string_from_name=None,
                 log=True,drop_markers=None,
                 outputDir=None):
    

    # spatialTablePath list or string
    if isinstance(spatialTablePath, str):
        spatialTablePath = [spatialTablePath]
    spatialTablePath = [pathlib.Path(p) for p in spatialTablePath]
    # probTablePath list or string
    if isinstance(probTablePath, str):
        probTablePath = [probTablePath]
    probTablePath = [pathlib.Path(p) for p in probTablePath]
    
    # Import spatialTablePath
    def load_process_data (image):
        # Print the data that is being processed
        print(f"Loading {image.name}")
        d = pd.read_csv(image)
        # If the data does not have a unique image ID column, add one.
        if 'imageid' not in d.columns:
            imid = image.stem
            d['imageid'] = imid
        # Unique name for the data
        if unique_CellId is True:
            d.index = d['imageid'].astype(str)+'_'+d[CellId].astype(str)
        else:
            d.index = d[CellId]
            
        # move image id and cellID column to end
        cellid_col = [col for col in d.columns if col != CellId] + [CellId]; d = d[cellid_col]
        imageid_col = [col for col in d.columns if col != 'imageid'] + ['imageid']; d = d[imageid_col]
        # If there is INF replace with zero
        d = d.replace([np.inf, -np.inf], 0)
        # Return data
        return d
    
    # Import probTablePath
    def load_process_probTable (image):
        d = pd.read_csv(image, index_col=0)
        # Return data
        return d
    
    # Apply function to all spatialTablePath and create a master dataframe
    r_load_process_data = lambda x: load_process_data(image=x) # Create lamda function
    all_spatialTable = list(map(r_load_process_data, list(spatialTablePath))) # Apply function
    # Merge all the spatialTablePath into a single large dataframe
    for i in range(len(all_spatialTable)):
        all_spatialTable[i].columns = all_spatialTable[0].columns
    entire_spatialTable = pd.concat(all_spatialTable, axis=0, sort=False)
    
    # Apply function to all probTablePath and create a master dataframe
    r_load_process_probTable = lambda x: load_process_probTable(image=x) # Create lamda function
    all_probTable = list(map(r_load_process_probTable, list(probTablePath))) # Apply function
    # Merge all the probTablePath into a single large dataframe
    for i in range(len(all_probTable)):
        all_probTable[i].columns = all_probTable[0].columns
    entire_probTable = pd.concat(all_probTable, axis=0, sort=False)
    # make the index of entire_probTable same as all_probTable
    ## NOTE THIS IS A HARD COPY WITHOUT ANY CHECKS! ASSUMES BOTH ARE IN SAME ORDER
    entire_probTable.index = entire_spatialTable.index
    
    
    # Split the data into expression data and meta data
    # Step-1 (Find the index of the column with name X_centroid)
    split_idx = entire_spatialTable.columns.get_loc(split)
    meta = entire_spatialTable.iloc [:,split_idx:]
    # Step-2 (select only the expression values)
    entire_spatialTable = entire_spatialTable.iloc [:,:split_idx]
    
    # Rename the columns of the data
    if remove_string_from_name is not None:
        entire_spatialTable.columns = entire_spatialTable.columns.str.replace(remove_string_from_name, '')

    # Save a copy of the column names in the uns space of ANNDATA
    markers = list(entire_spatialTable.columns)
    
    # Remove DNA channels
    if remove_dna is True:
        entire_spatialTable = entire_spatialTable.loc[:,~entire_spatialTable.columns.str.contains('dna', case=False)]

    # Drop unnecessary markers
    if drop_markers is not None:
        if isinstance(drop_markers, str):
            drop_markers = [drop_markers]
        entire_spatialTable = entire_spatialTable.drop(columns=drop_markers)
    
    # Create an anndata object
    adata = ad.AnnData(entire_spatialTable)
    adata.obs = meta
    adata.uns['all_markers'] = markers
    adata.uns['probQuant'] = entire_probTable
    
    # Add log data
    if log is True:
        adata.raw = adata
        adata.X = np.log1p(adata.X)
    
    # Save data if requested
    if outputDir is not None:
        finalPath = pathlib.Path(outputDir + '/GATOR/')
        if not os.path.exists(finalPath):
            os.makedirs(finalPath)
        if len(spatialTablePath) > 1:
            imid = 'gatorObject'
        else:
            imid = probTablePath[0].stem 
        adata.write(finalPath / f'{imid}.h5ad')
    else:    
        # Return data
        return adata
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
