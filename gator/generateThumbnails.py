#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:28:48 2022
@author: Ajit Johnson Nirmal
Generation of Thumbnais/ training data for a given marker
"""

# import packages
import numpy as np
import pandas as pd
import random
import tifffile
import os
import pathlib
import dask.array as da
import zarr


# Function
def generateThumbnails (csvPath, imagePath, markers, markerChannels, outputDir=None,
                        transformation=True, maxThumbnails=2000, localNorm=False,
                        x_coordinate='X_centroid', y_coordinate='Y_centroid',
                        threshold=[2, 12, 88, 98], windowSize=64):
    
    #transformation=True; maxThumbnails=2000; x_coordinate='X_centroid'; y_coordinate='Y_centroid'; threshold=[2, 12, 88, 98]; windowSize=64
    
    
    # Quick Check!
    if len(markers) is not len(markerChannels):
        raise ValueError('Every "marker" should have a corresponding "markerChannel", check the "markers" & "markerChannels" parameter' )
    
    # load the CSV to identify potential thumbnails
    data = pd.read_csv(pathlib.Path(csvPath))
    
    # subset the markers of interest
    if isinstance (markers, str):
        markers = [markers]
    if isinstance (markerChannels, int):
        markerChannels = [markerChannels]
    
    # convert markerChannels to zero indexing
    markerChannels = [x-1 for x in markerChannels]
    
    # creat a dict of marker and corresponding marker channel
    marker_map = dict(zip(markers,markerChannels))
        
    # create folders if it does not exist 
    if outputDir is None:
        outputDir = os.getcwd()
    
    # TruePos folders
    for i in markers:
        pos_path = pathlib.Path(outputDir + '/GATOR/Thumbnails/' + str(i) + '/TruePos')
        neg_path = pathlib.Path(outputDir + '/GATOR/Thumbnails/' + str(i) + '/TrueNeg')
        pos2neg_path = pathlib.Path(outputDir + '/GATOR/Thumbnails/' + str(i) + '/PosToNeg')
        neg2pos_path = pathlib.Path(outputDir + '/GATOR/Thumbnails/' + str(i) + '/NegToPos')
        if not os.path.exists(pos_path):
            os.makedirs(pos_path)
        if not os.path.exists(neg_path):
            os.makedirs(neg_path)
        if not os.path.exists(pos2neg_path):
            os.makedirs(pos2neg_path)
        if not os.path.exists(neg2pos_path):
            os.makedirs(neg2pos_path)
        if localNorm is True:
            local_pos_path = pathlib.Path(outputDir + '/GATOR/Thumbnails/localNorm/' + str(i) + '/TruePos')
            local_neg_path = pathlib.Path(outputDir + '/GATOR/Thumbnails/localNorm/' + str(i) + '/TrueNeg')   
            local_pos2neg_path = pathlib.Path(outputDir + '/GATOR/Thumbnails/localNorm/' + str(i) + '/PosToNeg')
            local_neg2pos_path = pathlib.Path(outputDir + '/GATOR/Thumbnails/localNorm/' + str(i) + '/NegToPos')
            if not os.path.exists(local_pos_path):
                os.makedirs(local_pos_path)
            if not os.path.exists(local_neg_path):
                os.makedirs(local_neg_path)
            if not os.path.exists(local_pos2neg_path):
                os.makedirs(local_pos2neg_path)
            if not os.path.exists(local_neg2pos_path):
                os.makedirs(local_neg2pos_path)
            
    marker_data = data[markers]
    location = data[[x_coordinate,y_coordinate]]
    
    # clip the data to drop outliers
    def clipping (x):
        clip = x.clip(lower =np.percentile(x,0.01), upper=np.percentile(x,99.99)).tolist()
        return clip
    
    # clip data
    marker_data = marker_data.apply(clipping)
    
    # apply transformation if requested
    if transformation is True:
        marker_data = np.arcsinh(marker_data)
        #marker_data = np.log1p(marker_data)
        
    # combine data
    combined_data = pd.concat([marker_data, location], axis=1)
    
    # intialize the threshold values
    threshold.sort()
    
    # function to identify the corner of the thumbnails
    def cornerFinder (centroid):
        row_start = int(centroid - windowSize // 2)
        row_end = row_start + windowSize
        return [row_start, row_end]
    
    # function to crop the image and save the image
    def cropImage (rowIndex, corners, imgType, zimg, npercentile, m):
        x_start = corners.iloc[rowIndex]['x_start']; x_end = corners.iloc[rowIndex]['x_end']
        y_start = corners.iloc[rowIndex]['y_start']; y_end = corners.iloc[rowIndex]['y_end']
        # cropping image
        crop = zimg[y_start:y_end, x_start:x_end]
        # convert the image to unit8
        fullN = ((crop/npercentile)*255).clip(0, 255).astype('uint8')
        # save the cropped image
        if imgType == 'pos':
            path = pathlib.Path(outputDir + '/GATOR/Thumbnails/' + str(m) + '/TruePos/' + str(rowIndex) + '.tif')
        elif imgType == 'neg':
            path = pathlib.Path(outputDir + '/GATOR/Thumbnails/' + str(m) + '/TrueNeg/' + str(rowIndex) + '.tif')
        # write file
        tifffile.imwrite(path,fullN)
        # local normalization if requested
        if localNorm is True:
            localN = ((crop/(np.percentile(crop.compute(), 99.99)))*255).clip(0, 255).astype('uint8')
            # save image
            if imgType == 'pos':
                Lpath = pathlib.Path(outputDir + '/GATOR/Thumbnails/localNorm/' + str(m) + '/TruePos/' + str(rowIndex) + '.tif')
            elif imgType == 'neg':
                Lpath = pathlib.Path(outputDir + '/GATOR/Thumbnails/localNorm/' + str(m) + '/TrueNeg/' + str(rowIndex) + '.tif')
            # write file
            tifffile.imwrite(Lpath,localN)

    # identify the cells of interest 
    def processMarker (marker):
        print('Processing Marker: ' + str(marker))
        
        moi = combined_data[marker]
        
        # figure out marker index or channel in image
        markerIndex = marker_map[marker]
        
        # determine the threshold value for the marker of interest
        low_a = np.percentile(moi, threshold[0]); low_b = np.percentile(moi, threshold[1])
        high_a = np.percentile(moi, threshold[2]); high_b = np.percentile(moi, threshold[3])
        
        # identify the cells that fall within the determined range
        neg = np.where(moi.between(low_a, low_b))[0]
        pos = np.where(moi.between(high_a, high_b))[0]
        
        # shuffle the cells
        random.shuffle(neg); np.random.shuffle(pos)
        
        # subset max number of cells
        if len(neg) > maxThumbnails:
            neg = neg[:maxThumbnails]
        if len(pos) > maxThumbnails: 
            pos = pos[:maxThumbnails]

        # identify the location of pos and neg cells
        neg_location = location.iloc[neg]
        pos_location = location.iloc[pos]
        

        # Find corner
        # Negative cells
        r_cornerFinder = lambda x: cornerFinder (centroid=x) 
        neg_x = pd.DataFrame(list(map(r_cornerFinder, neg_location[x_coordinate].values))) # x direction
        neg_y = pd.DataFrame(list(map(r_cornerFinder, neg_location[y_coordinate].values))) # y direction
        neg_x.columns = ["x_start", "x_end"]; neg_y.columns = ["y_start", "y_end"] 
        neg_location = pd.concat([neg_x, neg_y], axis=1)
        
        
        # Positive cells
        r_cornerFinder = lambda x: cornerFinder (centroid=x) 
        pos_x = pd.DataFrame(list(map(r_cornerFinder, pos_location[x_coordinate].values))) # x direction
        pos_y = pd.DataFrame(list(map(r_cornerFinder, pos_location[y_coordinate].values))) # y direction
        pos_x.columns = ["x_start", "x_end"]; pos_y.columns = ["y_start", "y_end"] 
        pos_location = pd.concat([pos_x, pos_y], axis=1)
        
        # load the image
        zimg = da.from_zarr(tifffile.imread(pathlib.Path(imagePath), aszarr=True, level=0, key=markerIndex))
        npercentile = np.percentile(zimg.compute(), 99.99)

        # for older version f tifffile
        #zimg = zarr.open(tifffile.imread(pathlib.Path(imagePath), aszarr=True, level=0, key=markerIndex))
        # = np.percentile(zimg,  99.99)
        
        # Cut images and write it out
        # neg
        r_cropImage = lambda x: cropImage (rowIndex=x, corners=neg_location, imgType='neg', zimg=zimg, npercentile=npercentile, m=marker)
        process_neg = list(map(r_cropImage, list(neg_location.index)))
        # pos
        r_cropImage = lambda x: cropImage (rowIndex=x, corners=pos_location, imgType='pos', zimg=zimg, npercentile=npercentile, m=marker)
        process_neg = list(map(r_cropImage, list(pos_location.index)))
        
    
    # Run the function for each marker
    r_processMarker = lambda x: processMarker (marker=x)
    final = list(map(r_processMarker, markers))
    
    # Finish Job
    print('Mission Accomplished')





