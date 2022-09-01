# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 15:08:38 2022
@author: Ajit Johnson Nirmal
When LocalNorm is used to curate the thumbnails, use this function to clone the 
curation on the real files
"""

# lib
import pathlib
from os import walk
import os
import shutil

# function
def cloneFolder (copyFolder, applyFolder, 
                 TruePos='TruePos', TrueNeg='TrueNeg', 
                 PosToNeg='PosToNeg', NegToPos='NegToPos'):
    
    #TruePos='TruePos'; TrueNeg='TrueNeg'; PosToNeg='PosToNeg'; NegToPos='NegToPos'
    
    # Convert the path to list 
    if isinstance (copyFolder, str):
        copyFolder = [copyFolder]
    if isinstance (applyFolder, str):
        applyFolder = [applyFolder]
    
    # Quick Check!
    if len(copyFolder) is not len(applyFolder):
        raise ValueError('The number of copyFolder and applyFolder should match, please check!' )
    
    # function to delete images
    def deleteFile(files, location):
        for f in files:
            # full path 
            #full_path = location  + f
            full_path = pathlib.Path.joinpath(location, f)
            if os.path.exists(full_path):
                os.remove(full_path)
    
    # Function to move images
    def moveFile(files, from_loc, to_loc):
        for f in files:
            # full path
            full_path_from = pathlib.Path.joinpath(from_loc, f) # from_loc + f
            full_path_to = pathlib.Path.joinpath(to_loc, f) # to_loc + f
            # move file
            if os.path.exists(full_path_from):
                shutil.move(full_path_from, full_path_to)
    
    # path lib of all folder         
    all_folders = [pathlib.Path(p) for p in copyFolder]
        
    # copy from location
    pos_aug_location = [pathlib.Path(p + '/' + str(TruePos)) for p in copyFolder]
    neg_aug_location = [pathlib.Path(p + '/' + str(TrueNeg)) for p in copyFolder]
    pos2neg_aug_location = [pathlib.Path(p + '/' + str(PosToNeg)) for p in copyFolder]
    neg2pos_aug_location = [pathlib.Path(p + '/' + str(NegToPos)) for p in copyFolder]
    
    # copy to location
    pos_real_location = [pathlib.Path(p + '/' + str(TruePos)) for p in applyFolder]
    neg_real_location = [pathlib.Path(p + '/' + str(TrueNeg)) for p in applyFolder]
    pos2neg_real_location = [pathlib.Path(p + '/' + str(PosToNeg)) for p in applyFolder]
    neg2pos_real_location = [pathlib.Path(p + '/' + str(NegToPos)) for p in applyFolder]
    

    
    # function
    def processFolder (folderIndex):
        print ('Processing: ' + str(all_folders[folderIndex].stem))
        
        # create a list of all file names in the applyFolder
        pos_files = next(walk(pos_real_location[folderIndex]), (None, None, []))[2]  
        neg_files = next(walk(neg_real_location[folderIndex]), (None, None, []))[2] 
        
        # Find file names within each of the copyFolder
        pos = next(walk(pos_aug_location[folderIndex]), (None, None, []))[2]  
        neg = next(walk(neg_aug_location[folderIndex]), (None, None, []))[2]  
        pos2neg = next(walk(pos2neg_aug_location[folderIndex]), (None, None, []))[2]  
        neg2pos = next(walk(neg2pos_aug_location[folderIndex]), (None, None, []))[2]  # [] if no file
        
        # Find images to delete
        pos_del = list(set(pos_files).difference(pos + pos2neg))
        neg_del = list(set(neg_files).difference(neg + neg2pos))
        
        # delete files
        deleteFile(files=pos_del, location=pos_real_location[folderIndex]) 
        deleteFile(files=neg_del, location=neg_real_location[folderIndex]) 
        
        # move files
        moveFile (files=pos2neg, from_loc=pos_real_location[folderIndex], to_loc=pos2neg_real_location[folderIndex])
        moveFile (files=neg2pos, from_loc=neg_real_location[folderIndex], to_loc=neg2pos_real_location[folderIndex])
        
        # print the number of files
        posaug = len(next(walk(pos_aug_location[folderIndex]), (None, None, []))[2])
        posreal = len(next(walk(pos_real_location[folderIndex]), (None, None, []))[2])
        negaug = len(next(walk(neg_aug_location[folderIndex]), (None, None, []))[2])
        negreal = len(next(walk(neg_real_location[folderIndex]), (None, None, []))[2])
        postonegaug = len(next(walk(pos2neg_aug_location[folderIndex]), (None, None, []))[2])
        postonegreal = len(next(walk(pos2neg_real_location[folderIndex]), (None, None, []))[2])
        negtoposaug = len(next(walk(neg2pos_aug_location[folderIndex]), (None, None, []))[2])
        negtoposreal = len(next(walk(neg2pos_real_location[folderIndex]), (None, None, []))[2])
        
        #print ('No of Files in TruePos-> copyFolder: ' + str(posaug) + ' ; applyFolder: '+ str(posreal))
        #print ('No of Files in TrueNeg-> copyFolder: ' + str(negaug) + ' ; applyFolder: '+ str(negreal))
        #print ('No of Files in PosToNeg-> copyFolder: ' + str(postonegaug) + ' ; applyFolder: '+ str(postonegreal))
        #print ('No of Files in NegToPos-> copyFolder: ' + str(negtoposaug) + ' ; applyFolder: '+ str(negtoposreal))
    
    # apply function to all folders
    r_processFolder = lambda x: processFolder (folderIndex=x)
    process_folders = list(map(r_processFolder, list(range(len(copyFolder)))))
        









