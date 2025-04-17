# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 10:11:06 2025

@author: Tiho
"""

#input of the algorithm: Exp_1_P.tif, Exp_1_GT.tif, Exp_2_P.tif, Exp_2_GT.tif, Exp_3_P.tif, Exp_3_GT.tif, train_GT.tif, Exp_1_P_masked, Exp_1_GT_masked.tif,
#Exp_2_P_masked.tif, Exp_2_GT_masked.tif, Exp_3_P_masked.tif, Exp_3_GT_masked.tif, train_GT_masked.tif, Exp_1_GT_bn_sparks.tif, Exp_2_GT_bn_sparks.tif,
#Exp_3_GT_bn_sparks.tif, train_GT_bn_sparks.tif, Exp_1_SHG.tif, Exp_2_SHG.tif, Exp_3_SHG.tif, train_SHG.tif, Exp_1_P_bn_sparks, Exp_2_P_bn_sparks, Exp_3_P_bn_sparks
#output of the algorithm: violin plot of SHG signal(spark)/SHG signal(fiber)
#the files must be placed in the same folder as the .py file but this can be modified
#the approach to read out data is similar to the one described at https://github.com/embl-bio-it/image-analysis-with-python/blob/master/session-3to5/image_analysis_tutorial_solutions.ipynb
#the module skimage.measure.regionprops_table can be used instead

import numpy as np
import scipy.ndimage as ndi

import matplotlib.pyplot as plt

from tifffile import imread
from tifffile import imsave
from scipy import stats

fiber_shg_cat_1 = imread("Exp_1_SHG.tif")
fiber_shg_cat_2 = imread("Exp_2_SHG.tif")
fiber_shg_cat_3 = imread("Exp_3_SHG.tif")
fiber_shg_train = imread("train_SHG.tif")

fiber_sparks_bn_predicted_cat_1 = imread("Exp_1_P_bn_sparks.tif")
fiber_sparks_bn_predicted_cat_2 = imread("Exp_2_P_bn_sparks.tif")
fiber_sparks_bn_predicted_cat_3 = imread("Exp_3_P_bn_sparks.tif")

stack_bn_cat_1 = imread("Exp_1_GT_bn_sparks.tif")
stack_bn_cat_2 = imread("Exp_2_GT_bn_sparks.tif")
stack_bn_cat_3 = imread("Exp_3_GT_bn_sparks.tif")
stack_bn_train = imread("train_GT_bn_sparks.tif")

fiber_ca_predicted_mask_clean_cat_1 = imread("Exp_1_P_masked.tif")
fiber_ca_predicted_mask_clean_cat_2 = imread("Exp_1_P_masked.tif")
fiber_ca_predicted_mask_clean_cat_3 = imread("Exp_1_P_masked.tif")

fiber_ca_cat_1_masked = imread("Exp_1_GT_masked.tif")
fiber_ca_cat_2_masked = imread("Exp_2_GT_masked.tif")
fiber_ca_cat_3_masked = imread("Exp_3_GT_masked.tif")
fiber_ca_masked_train = imread("train_GT_masked.tif")

fiber_ca_predicted_cat_1 = imread("Exp_1_P.tif")
fiber_ca_predicted_cat_2 = imread("Exp_2_P.tif")
fiber_ca_predicted_cat_3 = imread("Exp_3_P.tif")

fiber_ca_cat_1 = imread("Exp_1_GT.tif")
fiber_ca_cat_2 = imread("Exp_2_GT.tif")
fiber_ca_cat_3 = imread("Exp_3_GT.tif")
fiber_ca_train = imread("train_GT.tif")


stack_bn_cat_1_label = np.zeros(shape=(stack_bn_cat_1.shape), dtype=np.uint32)
stack_bn_cat_2_label = np.zeros(shape=(stack_bn_cat_2.shape), dtype=np.uint32)
stack_bn_cat_3_label = np.zeros(shape=(stack_bn_cat_3.shape), dtype=np.uint32)
stack_bn_train_label = np.zeros(shape=(stack_bn_train.shape), dtype=np.uint32)

fiber_sparks_bn_predicted_cat_1_label = np.zeros(shape=(fiber_sparks_bn_predicted_cat_1.shape), dtype=np.uint32)
fiber_sparks_bn_predicted_cat_2_label = np.zeros(shape=(fiber_sparks_bn_predicted_cat_2.shape), dtype=np.uint32)
fiber_sparks_bn_predicted_cat_3_label = np.zeros(shape=(fiber_sparks_bn_predicted_cat_3.shape), dtype=np.uint32)


for i in range(stack_bn_cat_1.shape[0]):
    stack_bn_cat_1_label[i], stack_bn_cat_1_label_number = ndi.label(stack_bn_cat_1[i])  

for i in range(stack_bn_cat_2.shape[0]):
    stack_bn_cat_2_label[i], stack_bn_cat_2_label_number = ndi.label(stack_bn_cat_2[i]) 
    
for i in range(stack_bn_cat_3.shape[0]):
    stack_bn_cat_3_label[i], stack_bn_cat_3_label_number = ndi.label(stack_bn_cat_3[i])    
    
for i in range(stack_bn_train.shape[0]):
    stack_bn_train_label[i], stack_bn_train_label_number = ndi.label(stack_bn_train[i])    
    
    

for i in range(fiber_sparks_bn_predicted_cat_1.shape[0]):
    fiber_sparks_bn_predicted_cat_1_label[i], fiber_sparks_bn_predicted_cat_1_label_number = ndi.label(fiber_sparks_bn_predicted_cat_1[i])

for i in range(fiber_sparks_bn_predicted_cat_2.shape[0]):
    fiber_sparks_bn_predicted_cat_2_label[i], fiber_sparks_bn_predicted_cat_2_label_number = ndi.label(fiber_sparks_bn_predicted_cat_2[i])

for i in range(fiber_sparks_bn_predicted_cat_3.shape[0]):
    fiber_sparks_bn_predicted_cat_3_label[i], fiber_sparks_bn_predicted_cat_3_label_number = ndi.label(fiber_sparks_bn_predicted_cat_3[i])
  
    
results = {"SHG_spark_SHG_fiber_cat_1"    : [],
           "SHG_spark_SHG_fiber_cat_1_predicted"    : [],
           "Ca_spark_Ca_fiber_cat_1"    : [],
           "Ca_spark_Ca_fiber_cat_1_predicted"    : [],
           "SHG_spark_SHG_fiber_cat_2"    : [],
           "SHG_spark_SHG_fiber_cat_2_predicted"    : [],
           "Ca_spark_Ca_fiber_cat_2"    : [],
           "Ca_spark_Ca_fiber_cat_2_predicted"    : [],
           "SHG_spark_SHG_fiber_cat_3" : [],
           "SHG_spark_SHG_fiber_cat_3_predicted"    : [],
           "Ca_spark_Ca_fiber_cat_3" : [],
           "Ca_spark_Ca_fiber_cat_3_predicted"    : [],
           "SHG_spark_SHG_fiber_train" : [],
           "Ca_spark_Ca_fiber_train" : []}    

fiber_ca_predicted_mask_clean_cat_1_w_sparks = np.zeros(shape=(fiber_ca_predicted_mask_clean_cat_1.shape), dtype=np.uint8)
fiber_ca_cat_1_masked_w_sparks = np.zeros(shape=(fiber_ca_cat_1_masked.shape), dtype=np.uint8)

fiber_ca_predicted_mask_clean_cat_1_w_sparks = np.logical_xor(fiber_sparks_bn_predicted_cat_1, fiber_ca_predicted_mask_clean_cat_1)
fiber_ca_cat_1_masked_w_sparks = np.logical_xor(stack_bn_cat_1, fiber_ca_cat_1_masked)

fiber_ca_predicted_mask_clean_cat_2_w_sparks = np.zeros(shape=(fiber_ca_predicted_mask_clean_cat_2.shape), dtype=np.uint8)
fiber_ca_cat_2_masked_w_sparks = np.zeros(shape=(fiber_ca_cat_2_masked.shape), dtype=np.uint8)

fiber_ca_predicted_mask_clean_cat_2_w_sparks = np.logical_xor(fiber_sparks_bn_predicted_cat_2, fiber_ca_predicted_mask_clean_cat_2)
fiber_ca_cat_2_masked_w_sparks = np.logical_xor(stack_bn_cat_2, fiber_ca_cat_2_masked)

fiber_ca_predicted_mask_clean_cat_3_w_sparks = np.zeros(shape=(fiber_ca_predicted_mask_clean_cat_3.shape), dtype=np.uint8)
fiber_ca_cat_3_masked_w_sparks = np.zeros(shape=(fiber_ca_cat_3_masked.shape), dtype=np.uint8)

fiber_ca_predicted_mask_clean_cat_3_w_sparks = np.logical_xor(fiber_sparks_bn_predicted_cat_3, fiber_ca_predicted_mask_clean_cat_3)
fiber_ca_cat_3_masked_w_sparks = np.logical_xor(stack_bn_cat_3, fiber_ca_cat_3_masked)

fiber_ca_masked_train_w_sparks = np.zeros(shape=(fiber_ca_masked_train.shape), dtype=np.uint8)

fiber_ca_masked_train_w_sparks = np.logical_xor(stack_bn_train, fiber_ca_masked_train)

    
for i in range(stack_bn_cat_1.shape[0]):
    for spark_ID in np.unique(stack_bn_cat_1_label[i])[1:]:
        spark_mask_cat_1 = stack_bn_cat_1_label[i]==spark_ID
        results["SHG_spark_SHG_fiber_cat_1"].append(np.mean(fiber_shg_cat_1[i][spark_mask_cat_1])/np.mean(fiber_shg_cat_1[i][fiber_ca_cat_1_masked_w_sparks[i]]))
        results["Ca_spark_Ca_fiber_cat_1"].append(np.mean(fiber_ca_cat_1[i][spark_mask_cat_1])/np.mean(fiber_ca_cat_1[i][fiber_ca_cat_1_masked_w_sparks[i]]))
 
for i in range(stack_bn_cat_2.shape[0]):
    for spark_ID in np.unique(stack_bn_cat_2_label[i])[1:]:
        spark_mask_cat_2 = stack_bn_cat_2_label[i]==spark_ID
        results["SHG_spark_SHG_fiber_cat_2"].append(np.mean(fiber_shg_cat_2[i][spark_mask_cat_2])/np.mean(fiber_shg_cat_2[i][fiber_ca_cat_2_masked_w_sparks[i]]))
        results["Ca_spark_Ca_fiber_cat_2"].append(np.mean(fiber_ca_cat_2[i][spark_mask_cat_2])/np.mean(fiber_ca_cat_2[i][fiber_ca_cat_2_masked_w_sparks[i]]))

for i in range(stack_bn_cat_3.shape[0]):
    for spark_ID in np.unique(stack_bn_cat_3_label[i])[1:]:
        spark_mask_cat_3 = stack_bn_cat_3_label[i]==spark_ID
        results["SHG_spark_SHG_fiber_cat_3"].append(np.mean(fiber_shg_cat_3[i][spark_mask_cat_3])/np.mean(fiber_shg_cat_3[i][fiber_ca_cat_3_masked_w_sparks[i]]))
        results["Ca_spark_Ca_fiber_cat_3"].append(np.mean(fiber_ca_cat_3[i][spark_mask_cat_3])/np.mean(fiber_ca_cat_3[i][fiber_ca_cat_3_masked_w_sparks[i]]))
    
     
for i in range(fiber_sparks_bn_predicted_cat_1.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_1_label[i])[1:]:
        spark_mask_predicted_cat_1 = fiber_sparks_bn_predicted_cat_1_label[i]==spark_ID
        results["SHG_spark_SHG_fiber_cat_1_predicted"].append(np.mean(fiber_shg_cat_1[i][spark_mask_predicted_cat_1])/np.mean(fiber_shg_cat_1[i][fiber_ca_predicted_mask_clean_cat_1_w_sparks[i]]))
        results["Ca_spark_Ca_fiber_cat_1_predicted"].append(np.mean(fiber_ca_predicted_cat_1[i][spark_mask_predicted_cat_1])/np.mean(fiber_ca_predicted_cat_1[i][fiber_ca_predicted_mask_clean_cat_1_w_sparks[i]]))        

for i in range(fiber_sparks_bn_predicted_cat_2.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_2_label[i])[1:]:
        spark_mask_predicted_cat_2 = fiber_sparks_bn_predicted_cat_2_label[i]==spark_ID
        results["SHG_spark_SHG_fiber_cat_2_predicted"].append(np.mean(fiber_shg_cat_2[i][spark_mask_predicted_cat_2])/np.mean(fiber_shg_cat_2[i][fiber_ca_predicted_mask_clean_cat_2_w_sparks[i]]))
        results["Ca_spark_Ca_fiber_cat_2_predicted"].append(np.mean(fiber_ca_predicted_cat_2[i][spark_mask_predicted_cat_2])/np.mean(fiber_ca_predicted_cat_2[i][fiber_ca_predicted_mask_clean_cat_2_w_sparks[i]]))

for i in range(fiber_sparks_bn_predicted_cat_3.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_3_label[i])[1:]:
        spark_mask_predicted_cat_3 = fiber_sparks_bn_predicted_cat_3_label[i]==spark_ID
        results["SHG_spark_SHG_fiber_cat_3_predicted"].append(np.mean(fiber_shg_cat_3[i][spark_mask_predicted_cat_3])/np.mean(fiber_shg_cat_3[i][fiber_ca_predicted_mask_clean_cat_3_w_sparks[i]]))
        results["Ca_spark_Ca_fiber_cat_3_predicted"].append(np.mean(fiber_ca_predicted_cat_3[i][spark_mask_predicted_cat_3])/np.mean(fiber_ca_predicted_cat_3[i][fiber_ca_predicted_mask_clean_cat_3_w_sparks[i]]))

for i in range(stack_bn_train.shape[0]):
    for spark_ID in np.unique(stack_bn_train_label[i])[1:]:
        spark_mask_train = stack_bn_train_label[i]==spark_ID
        results["SHG_spark_SHG_fiber_train"].append(np.mean(fiber_shg_train[i][spark_mask_train])/np.mean(fiber_shg_train[i][fiber_ca_masked_train_w_sparks[i]]))
        results["Ca_spark_Ca_fiber_train"].append(np.mean(fiber_ca_train[i][spark_mask_train])/np.mean(fiber_ca_train[i][fiber_ca_masked_train_w_sparks[i]]))


plt.figure()
data = [results["SHG_spark_SHG_fiber_cat_1"], results["SHG_spark_SHG_fiber_cat_1_predicted"], results["SHG_spark_SHG_fiber_cat_2"], results["SHG_spark_SHG_fiber_cat_2_predicted"], results["SHG_spark_SHG_fiber_cat_3"], results["SHG_spark_SHG_fiber_cat_3_predicted"], results["SHG_spark_SHG_fiber_train"]]
fig, ax = plt.subplots(figsize=(15,5))
ax.set_ylabel('SHG signal(spark)/SHG signal(fiber)')
bp = ax.boxplot(data, showmeans=True)
vp = ax.violinplot(data)

for body in bp['boxes']:
    body.set_alpha(0.7)
    
for body in bp['medians']:
    body.set_alpha(0.7)

for body in bp['whiskers']:
    body.set_alpha(0.7)

for body in bp['caps']:
    body.set_alpha(0.7) 
    
for body in bp['fliers']:
    body.set_alpha(0.7)

for body in bp['means']:
    body.set_alpha(0.7) 

for body in vp['bodies']:
    body.set_alpha(0.7)

plt.xticks([1,2,3,4,5,6,7], ["Exp 1 GT", "Exp 1 P", "Exp 2 GT", "Exp 2 P", "Exp 3 GT", "Exp 3 P", "Training"])
plt.savefig("SHG_signal_spark_fiber.pdf", format='pdf')
plt.show()