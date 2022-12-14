#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 16:37:29 2022
@author: Ajit Johnson Nirmal
Function to generate Masks for UNET model
"""

import pathlib
import cv2 as cv
import random
import numpy as np
import tifffile



def generateTraining (thumbnailFolder, outputDir, file_extension=None,
                      TruePos='TruePos', NegToPos='NegToPos',
                      TrueNeg='TrueNeg', PosToNeg='PosToNeg'):
    # Function takes in path to two folders, processes the images in those folders,
    # and saves them into a different folder that contains Train, Validation and Test samples
    #TruePos='TruePos'; NegToPos='NegToPos'; TrueNeg='TrueNeg'; PosToNeg='PosToNeg'
    
    # convert the folder into a list
    if isinstance (thumbnailFolder, str):
        thumbnailFolder = [thumbnailFolder]
    
    # convert all path names to pathlib
    thumbnailFolder = [pathlib.Path(p) for p in thumbnailFolder]
    outputDir = pathlib.Path(outputDir)
    
    # find all markers passed
    all_markers = [i.stem for i in thumbnailFolder]
    
    # create directories to save
    for i in all_markers:
        if not (outputDir / 'GATOR/TrainingData/' / f"{i}" /  'training').exists ():
            (outputDir / 'GATOR/TrainingData/' / f"{i}" /  'training').mkdir(parents=True, exist_ok=True)
        
        if not (outputDir / 'GATOR/TrainingData/' / f"{i}" /  'validation').exists ():
            (outputDir / 'GATOR/TrainingData/' / f"{i}" /  'validation').mkdir(parents=True, exist_ok=True)
        
        if not (outputDir / 'GATOR/TrainingData/' / f"{i}" /  'test').exists ():
            (outputDir / 'GATOR/TrainingData/' / f"{i}" /  'test').mkdir(parents=True, exist_ok=True)
    
    # standard format
    if file_extension is None:
        file_extension = '*'
    else:
        file_extension = '*' + str(file_extension)
    
    # Filter on pos cells
    def pos_filter (path):
        image = cv.imread(str(path.resolve()), cv.IMREAD_GRAYSCALE)
        blur = cv.GaussianBlur(image, ksize=(3,3), sigmaX=1, sigmaY=1)
        ret3,th3 = cv.threshold(blur,0,1,cv.THRESH_OTSU)
        mask = th3 + 1
        return [mask, image]
    
    # Filter on neg cells
    def neg_filter (path):
        image = cv.imread(str(path.resolve()), cv.IMREAD_GRAYSCALE)
        mask = np.ones(image.shape, dtype=np.uint8)
        return [mask, image]
    
    # identify the files within all the 4 folders
    def findFiles (folderIndex):
        print ('Processing: ' + str(thumbnailFolder[folderIndex].stem))
        marker_name = str(thumbnailFolder[folderIndex].stem)
        
        baseFolder = thumbnailFolder[folderIndex]
        
        if TruePos is not None:
            pos = list(pathlib.Path.glob(baseFolder / TruePos, file_extension))
        if NegToPos is not None:
            negtopos = list(pathlib.Path.glob(baseFolder / NegToPos, file_extension))
        positive_cells = pos + negtopos
        
        if TrueNeg is not None:
            neg = list(pathlib.Path.glob(baseFolder / TrueNeg, file_extension))
        if PosToNeg is not None:
            postoneg = list(pathlib.Path.glob(baseFolder / PosToNeg, file_extension))
        negative_cells = neg + postoneg
        
        # prepare the Training, Validataion and Test Cohorts
        if len(positive_cells) > 0: 
            train_pos = random.sample(positive_cells, round(len(positive_cells) * 0.6))
            remanining_pos = list(set(positive_cells) - set(train_pos))
            val_pos = random.sample(remanining_pos, round(len(remanining_pos) * 0.5)) # validation
            test_pos = list(set(remanining_pos) - set(val_pos)) # test
        if len(negative_cells) > 0:
            train_neg = random.sample(negative_cells, round(len(negative_cells) * 0.6))
            remanining_neg = list(set(negative_cells) - set(train_neg))
            val_neg = random.sample(remanining_neg, round(len(remanining_neg) * 0.5))
            test_neg = list(set(remanining_neg) - set(val_neg))

        
        # loop through training dataset and save images and masks
        newname_train = list(range(len(train_pos) + len(train_neg))); random.shuffle(newname_train)
        train_pos_name = newname_train[:len(train_pos)]; train_neg_name = newname_train[len(train_pos):]
        
        if len (train_pos_name) > 0:
            for i, j in zip( train_pos_name, train_pos):
                m, im = pos_filter (j)
                # save image
                fPath = outputDir / 'GATOR/TrainingData/' / f"{marker_name}" / 'training' / f"{i}_img.tif"
                tifffile.imwrite(fPath,im)
                # associated mask
                fPath = outputDir / 'GATOR/TrainingData/' / f"{marker_name}" / 'training' / f"{i}_mask.tif"
                tifffile.imwrite(fPath, m)
                
        if len (train_neg_name) > 0:
            for k, l in zip( train_neg_name, train_neg):
                m, im = neg_filter (l)
                # save image
                fPath = outputDir / 'GATOR/TrainingData/' / f"{marker_name}" / 'training' / f"{k}_img.tif"
                tifffile.imwrite(fPath, im)
                # associated mask
                fPath = outputDir / 'GATOR/TrainingData/' / f"{marker_name}" / 'training' / f"{k}_mask.tif"
                tifffile.imwrite(fPath, m)
            
        
        # loop through validation dataset and save images and masks
        newname_train = list(range(len(val_pos) + len(val_neg))); random.shuffle(newname_train)
        train_pos_name = newname_train[:len(val_pos)]; train_neg_name = newname_train[len(val_pos):]
        
        if len (train_pos_name) > 0:
            for i, j in zip( train_pos_name, val_pos):
                m, im = pos_filter (j)
                # save image
                fPath = outputDir / 'GATOR/TrainingData/' / f"{marker_name}" / 'validation' / f"{i}_img.tif"
                tifffile.imwrite(fPath, im)
                # associated mask
                fPath = outputDir / 'GATOR/TrainingData/' / f"{marker_name}" / 'validation' / f"{i}_mask.tif"
                tifffile.imwrite(fPath, m)
                
        if len (train_neg_name) > 0:
            for k, l in zip( train_neg_name, val_neg):
                m, im = neg_filter (l)
                # save image
                fPath = outputDir / 'GATOR/TrainingData/' / f"{marker_name}" / 'validation' / f"{k}_img.tif"
                tifffile.imwrite(fPath, im)
                # associated mask
                fPath = outputDir / 'GATOR/TrainingData/' / f"{marker_name}" / 'validation' / f"{k}_mask.tif"
                tifffile.imwrite(fPath, m)
            
            
        # loop through test dataset and save images and masks
        newname_train = list(range(len(test_pos) + len(test_neg))); random.shuffle(newname_train)
        train_pos_name = newname_train[:len(test_pos)]; train_neg_name = newname_train[len(test_pos):]
        
        if len (train_pos_name) > 0:
            for i, j in zip( train_pos_name, test_pos):
                m, im = pos_filter (j)
                # save image
                fPath = outputDir / 'GATOR/TrainingData/' / f"{marker_name}" / 'test' / f"{i}_img.tif"
                tifffile.imwrite(fPath, im)
                # associated mask
                fPath = outputDir / 'GATOR/TrainingData/' / f"{marker_name}" / 'test' / f"{i}_mask.tif"
                tifffile.imwrite(fPath, m)
                
        if len (train_neg_name) > 0:
            for k, l in zip( train_neg_name, test_neg):
                m, im = neg_filter (l)
                # save image
                fPath = outputDir / 'GATOR/TrainingData/' / f"{marker_name}" / 'test' / f"{k}_img.tif"
                tifffile.imwrite(fPath, im)
                # associated mask
                fPath = outputDir / 'GATOR/TrainingData/' / f"{marker_name}" / 'test' / f"{k}_mask.tif"
                tifffile.imwrite(fPath, m)
    
    # apply function to all folders
    r_findFiles = lambda x: findFiles (folderIndex=x)
    process_folders = list(map(r_findFiles, list(range(len(thumbnailFolder)))))
    
    # Print
    print('Mission Accomplished')

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    