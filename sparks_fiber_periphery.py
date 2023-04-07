import numpy as np
import scipy.ndimage as ndi
import matplotlib.pyplot as plt
from tifffile import imread
from tifffile import imsave

#input of the algorithm: Exp_1_P.tif, Exp_1_GT.tif, Exp_2_P.tif, Exp_2_GT.tif, Exp_3_P.tif, Exp_3_GT.tif, train_GT.tif, Exp_1_P_masked, Exp_1_GT_masked.tif,
#Exp_2_P_masked.tif, Exp_2_GT_masked.tif, Exp_3_P_masked.tif, Exp_3_GT_masked.tif, train_GT_masked.tif, Exp_1_GT_bn_sparks.tif, Exp_2_GT_bn_sparks.tif,
#Exp_3_GT_bn_sparks.tif, train_GT_bn_sparks.tif, Exp_1_P_bn_sparks.tif, Exp_2_P_bn_sparks.tif, Exp_3_P_bn_sparks.tif
#output of the algorithm: Exp_1_GT_fiber_periphery.tif, Exp_2_GT_fiber_periphery.tif, Exp_3_GT_fiber_periphery.tif, train_GT_fiber_periphery.tif,
#Exp_1_P_fiber_periphery.tif, Exp_2_P_fiber_periphery.tif, Exp_3_P_fiber_periphery.tif, Spark frequency border.pdf
#the files must be placed in the same folder as the .py file but this can be modified
#the approach to read out data is similar to the one described at https://github.com/embl-bio-it/image-analysis-with-python/blob/master/session-3to5/image_analysis_tutorial_solutions.ipynb
#the module skimage.measure.regionprops_table can be used instead
fiber_ca_predicted_cat_1 = imread("Exp_1_P.tif")
fiber_sparks_bn_predicted_cat_1 = imread("Exp_1_P_bn_sparks.tif")
fiber_ca_cat_1 = imread("Exp_1_GT.tif")
fiber_bn_cat_1 = imread("Exp_1_GT_bn_sparks.tif")
fiber_ca_predicted_mask_clean_cat_1 = imread("Exp_1_P_masked.tif")
fiber_ca_cat_1_masked = imread("Exp_1_GT_masked.tif")

fiber_ca_predicted_cat_2 = imread("Exp_2_P.tif")
fiber_sparks_bn_predicted_cat_2 = imread("Exp_2_P_bn_sparks.tif")
fiber_ca_cat_2 = imread("Exp_2_GT.tif")
fiber_bn_cat_2 = imread("Exp_2_GT_bn_sparks.tif")
fiber_ca_predicted_mask_clean_cat_2 = imread("Exp_2_P_masked.tif")
fiber_ca_cat_2_masked = imread("Exp_2_GT_masked.tif")

fiber_ca_predicted_cat_3 = imread("Exp_3_P.tif")
fiber_sparks_bn_predicted_cat_3 = imread("Exp_3_P_bn_sparks.tif")
fiber_ca_cat_3 = imread("Exp_3_GT.tif")
fiber_bn_cat_3 = imread("Exp_3_GT_bn_sparks.tif")
fiber_ca_predicted_mask_clean_cat_3 = imread("Exp_3_P_masked.tif")
fiber_ca_cat_3_masked = imread("Exp_3_GT_masked.tif")

results = {"spark_area_cat_1"      : [],
           "spark_area_cat_1_predicted"     : [],
           "spark_mean_int_cat_1"     : [],
           "spark_mean_int_cat_1_predicted"     : [],
           "spark_area_cat_2"    : [],
           "spark_area_cat_2_predicted"    : [],
           "spark_mean_int_cat_2" : [],
           "spark_mean_int_cat_2_predicted"     : [],
           "spark_area_cat_3" : [],
           "spark_area_cat_3_predicted"    : [],
           "spark_mean_int_cat_3" : [],
           "spark_mean_int_cat_3_predicted"     : [],
           "spark_area_train" : [],
           "spark_mean_int_train" : []}

fiber_ca_train = imread("train_GT.tif")
fiber_bn_train = imread("train_GT_bn_sparks.tif")
fiber_ca_masked_train = imread("train_GT_masked.tif")


fiber_sparks_bn_predicted_cat_1_label = np.zeros(shape=(fiber_ca_predicted_cat_1.shape), dtype=np.uint32)
fiber_sparks_bn_predicted_cat_2_label = np.zeros(shape=(fiber_ca_predicted_cat_2.shape), dtype=np.uint32) 
fiber_sparks_bn_predicted_cat_3_label = np.zeros(shape=(fiber_ca_predicted_cat_3.shape), dtype=np.uint32)  
fiber_sparks_bn_training_label = np.zeros(shape=(fiber_bn_train.shape), dtype=np.uint32)
fiber_bn_cat_1_label = np.zeros(shape=(fiber_ca_predicted_cat_1.shape), dtype=np.uint32)
fiber_bn_cat_2_label = np.zeros(shape=(fiber_ca_predicted_cat_1.shape), dtype=np.uint32)
fiber_bn_cat_3_label = np.zeros(shape=(fiber_ca_predicted_cat_1.shape), dtype=np.uint32)

for i in range(fiber_ca_predicted_cat_1.shape[0]):
    fiber_sparks_bn_predicted_cat_1_label[i], fiber_sparks_bn_predicted_cat_1_label_number = ndi.label(fiber_sparks_bn_predicted_cat_1[i])

for i in range(fiber_ca_predicted_cat_2.shape[0]):
    fiber_sparks_bn_predicted_cat_2_label[i], fiber_sparks_bn_predicted_cat_2_label_number = ndi.label(fiber_sparks_bn_predicted_cat_2[i])

for i in range(fiber_ca_predicted_cat_3.shape[0]):
    fiber_sparks_bn_predicted_cat_3_label[i], fiber_sparks_bn_predicted_cat_3_label_number = ndi.label(fiber_sparks_bn_predicted_cat_3[i])
    
    
for i in range(fiber_ca_cat_1.shape[0]):
    fiber_bn_cat_1_label[i], fiber_bn_cat_1_label_number = ndi.label(fiber_bn_cat_1[i])

for i in range(fiber_ca_cat_2.shape[0]):
    fiber_bn_cat_2_label[i], fiber_bn_cat_2_label_number = ndi.label(fiber_bn_cat_2[i])

for i in range(fiber_ca_cat_3.shape[0]):
    fiber_bn_cat_3_label[i], fiber_bn_cat_3_label_number = ndi.label(fiber_bn_cat_3[i])
    
        
for i in range(fiber_bn_train.shape[0]):
    fiber_sparks_bn_training_label[i], fiber_sparks_bn_training_label_number = ndi.label(fiber_bn_train[i])


for i in range(fiber_bn_cat_1.shape[0]):
    for spark_ID in np.unique(fiber_bn_cat_1_label[i])[1:]:
        spark_mask_cat_1_1 = fiber_bn_cat_1_label[i]==spark_ID
        results["spark_area_cat_1"].append(0.2325149*0.2325149*np.sum(spark_mask_cat_1_1))
        results["spark_mean_int_cat_1"].append(np.mean(fiber_ca_cat_1[i][spark_mask_cat_1_1]))

for i in range(fiber_bn_cat_2.shape[0]):
    for spark_ID in np.unique(fiber_bn_cat_2_label[i])[1:]:
        spark_mask_cat_2_1 = fiber_bn_cat_2_label[i]==spark_ID
        results["spark_area_cat_2"].append(0.2325149*0.2325149*np.sum(spark_mask_cat_2_1))
        results["spark_mean_int_cat_2"].append(np.mean(fiber_ca_cat_2[i][spark_mask_cat_2_1]))

for i in range(fiber_bn_cat_3.shape[0]):
    for spark_ID in np.unique(fiber_bn_cat_3_label[i])[1:]:
        spark_mask_cat_3_1 = fiber_bn_cat_3_label[i]==spark_ID
        results["spark_area_cat_3"].append(0.2325149*0.2325149*np.sum(spark_mask_cat_3_1))
        results["spark_mean_int_cat_3"].append(np.mean(fiber_ca_cat_3[i][spark_mask_cat_3_1]))

for i in range(fiber_ca_predicted_cat_1.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_1_label[i])[1:]:
        spark_mask_cat_1_1_predicted = fiber_sparks_bn_predicted_cat_1_label[i]==spark_ID
        results["spark_area_cat_1_predicted"].append(0.2325149*0.2325149*np.sum(spark_mask_cat_1_1_predicted))
        results["spark_mean_int_cat_1_predicted"].append(np.mean(fiber_ca_predicted_cat_1[i][spark_mask_cat_1_1_predicted]))
        
                
for i in range(fiber_ca_predicted_cat_2.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_2_label[i])[1:]:
        spark_mask_cat_2_1_predicted = fiber_sparks_bn_predicted_cat_2_label[i]==spark_ID
        results["spark_area_cat_2_predicted"].append(0.2325149*0.2325149*np.sum(spark_mask_cat_2_1_predicted))
        results["spark_mean_int_cat_2_predicted"].append(np.mean(fiber_ca_predicted_cat_2[i][spark_mask_cat_2_1_predicted]))

for i in range(fiber_ca_predicted_cat_3.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_3_label[i])[1:]:
        spark_mask_cat_3_1_predicted = fiber_sparks_bn_predicted_cat_3_label[i]==spark_ID
        results["spark_area_cat_3_predicted"].append(0.2325149*0.2325149*np.sum(spark_mask_cat_3_1_predicted))
        results["spark_mean_int_cat_3_predicted"].append(np.mean(fiber_ca_predicted_cat_3[i][spark_mask_cat_3_1_predicted]))

for i in range(fiber_ca_train.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_training_label[i])[1:]:
        spark_mask_train_1 = fiber_sparks_bn_training_label[i]==spark_ID
        results["spark_area_train"].append(0.2325149*0.2325149*np.sum(spark_mask_train_1))
        results["spark_mean_int_train"].append(np.mean(fiber_ca_train[i][spark_mask_train_1]))



l_spark_frequency_cat_1 = results["spark_area_cat_1"]
l_spark_frequency_cat_1_predicted = results["spark_area_cat_1_predicted"]
l_spark_frequency_cat_2 = results["spark_area_cat_2"]
l_spark_frequency_cat_2_predicted = results["spark_area_cat_2_predicted"]
l_spark_frequency_cat_3 = results["spark_area_cat_3"]
l_spark_frequency_cat_3_predicted = results["spark_area_cat_3_predicted"]
l_spark_frequency_train = results["spark_area_train"]



border_mask_cat_1 = np.zeros(fiber_ca_cat_1_masked.shape, dtype=np.bool)

for i in range(fiber_ca_cat_1_masked.shape[0]):
    border_mask_cat_1[i] = ndi.binary_dilation(border_mask_cat_1[i], border_value=1)

fiber_border_1 = border_mask_cat_1 + fiber_ca_cat_1_masked

fiber_border_11 = fiber_border_1 > 1

fiber_ca_cat_1_masked_eroded = np.zeros(shape=(fiber_ca_cat_1_masked.shape[0], fiber_ca_cat_1_masked.shape[1], fiber_ca_cat_1_masked.shape[2], 42), dtype=np.uint8)

fiber_ca_cat_1_masked_eroded[0:fiber_ca_cat_1_masked.shape[0], 0:fiber_ca_cat_1_masked.shape[1], 0:fiber_ca_cat_1_masked.shape[2], 1] = fiber_ca_cat_1_masked

for i in range(fiber_ca_cat_1_masked.shape[0]):
    for y in range (40):
        fiber_ca_cat_1_masked_eroded[i, 0:fiber_ca_cat_1_masked.shape[1], 0:fiber_ca_cat_1_masked.shape[2], y+2] = ndi.binary_erosion(fiber_ca_cat_1_masked_eroded[i, 0:fiber_ca_cat_1_masked.shape[1], 0:fiber_ca_cat_1_masked.shape[2], y+1], iterations=1) + fiber_border_11[i, 0:fiber_ca_cat_1_masked.shape[1], 0:fiber_ca_cat_1_masked.shape[2]]    

periphery_fiber_cat_1 = np.logical_xor(fiber_ca_cat_1_masked, fiber_ca_cat_1_masked_eroded[0:fiber_ca_cat_1_masked.shape[0], 0:fiber_ca_cat_1_masked.shape[1], 0:fiber_ca_cat_1_masked.shape[2], 41])

imsave("Exp_1_GT_fiber_periphery.tif", periphery_fiber_cat_1.astype(np.uint8))

l_sparks_rand_cat_1 = []

fiber_bn_cat_1_rand = np.zeros(shape=(fiber_ca_cat_1.shape), dtype=np.uint8)

for i in range(fiber_bn_cat_1.shape[0]):
    fiber_bn_cat_1_rand[i] = periphery_fiber_cat_1[i] + fiber_bn_cat_1[i]
    for fiber_ID in np.unique(fiber_bn_cat_1_label[i])[1:]:
        spark_mask_cat_1_1_1 = fiber_bn_cat_1_label[i]==fiber_ID
        if np.sum(fiber_bn_cat_1_rand[i][spark_mask_cat_1_1_1])>np.sum(fiber_bn_cat_1[i][spark_mask_cat_1_1_1]):            
            l_sparks_rand_cat_1.append(np.mean(fiber_ca_cat_1[i][spark_mask_cat_1_1_1]))

del fiber_ca_cat_1_masked_eroded




fiber_border_1_p = border_mask_cat_1 + fiber_ca_predicted_mask_clean_cat_1

fiber_border_11_p = fiber_border_1_p > 1

fiber_ca_cat_1_p_masked_eroded = np.zeros(shape=(fiber_ca_cat_1_masked.shape[0], fiber_ca_cat_1_masked.shape[1], fiber_ca_cat_1_masked.shape[2], 42), dtype=np.uint8)

fiber_ca_cat_1_p_masked_eroded[0:fiber_ca_cat_1_masked.shape[0], 0:fiber_ca_cat_1_masked.shape[1], 0:fiber_ca_cat_1_masked.shape[2], 1] = fiber_ca_predicted_mask_clean_cat_1

for i in range(fiber_ca_cat_1_masked.shape[0]):
    for y in range (40):
        fiber_ca_cat_1_p_masked_eroded[i, 0:fiber_ca_cat_1_masked.shape[1], 0:fiber_ca_cat_1_masked.shape[2], y+2] = ndi.binary_erosion(fiber_ca_cat_1_p_masked_eroded[i, 0:fiber_ca_cat_1_masked.shape[1], 0:fiber_ca_cat_1_masked.shape[2], y+1], iterations=1) + fiber_border_11_p[i, 0:fiber_ca_cat_1_masked.shape[1], 0:fiber_ca_cat_1_masked.shape[2]]    

periphery_fiber_cat_1_p = np.logical_xor(fiber_ca_predicted_mask_clean_cat_1, fiber_ca_cat_1_p_masked_eroded[0:fiber_ca_cat_1_masked.shape[0], 0:fiber_ca_cat_1_masked.shape[1], 0:fiber_ca_cat_1_masked.shape[2], 41])

imsave("Exp_1_P_fiber_periphery.tif", periphery_fiber_cat_1_p.astype(np.uint8))

l_sparks_rand_cat_1_p = []

fiber_bn_cat_1_p_rand = np.zeros(shape=(fiber_sparks_bn_predicted_cat_1.shape), dtype=np.uint8)

for i in range(fiber_sparks_bn_predicted_cat_1.shape[0]):
    fiber_bn_cat_1_p_rand[i] = periphery_fiber_cat_1_p[i] + fiber_sparks_bn_predicted_cat_1[i]    
    for fiber_ID in np.unique(fiber_sparks_bn_predicted_cat_1_label[i])[1:]:
        spark_mask_cat_1_1_1_p = fiber_sparks_bn_predicted_cat_1_label[i]==fiber_ID
        if np.sum(fiber_bn_cat_1_p_rand[i][spark_mask_cat_1_1_1_p])>np.sum(fiber_sparks_bn_predicted_cat_1[i][spark_mask_cat_1_1_1_p]):            
            l_sparks_rand_cat_1_p.append(np.mean(fiber_ca_predicted_cat_1[i][spark_mask_cat_1_1_1_p]))

del fiber_sparks_bn_predicted_cat_1

del fiber_ca_predicted_cat_1

del fiber_ca_cat_1_p_masked_eroded



border_mask_cat_2 = np.zeros(fiber_ca_cat_2_masked.shape, dtype=np.bool)

for i in range(fiber_ca_cat_2_masked.shape[0]):
    border_mask_cat_2[i] = ndi.binary_dilation(border_mask_cat_2[i], border_value=1)

fiber_border_2 = border_mask_cat_2 + fiber_ca_cat_2_masked

fiber_border_21 = fiber_border_2 > 1

fiber_ca_cat_2_masked_eroded = np.zeros(shape=(fiber_ca_cat_2_masked.shape[0], fiber_ca_cat_2_masked.shape[1], fiber_ca_cat_2_masked.shape[2], 42), dtype=np.uint8)

fiber_ca_cat_2_masked_eroded[0:fiber_ca_cat_2_masked.shape[0], 0:fiber_ca_cat_2_masked.shape[1], 0:fiber_ca_cat_2_masked.shape[2], 1] = fiber_ca_cat_2_masked

for i in range(fiber_ca_cat_2_masked.shape[0]):
    for y in range (40):
        fiber_ca_cat_2_masked_eroded[i, 0:fiber_ca_cat_2_masked.shape[1], 0:fiber_ca_cat_2_masked.shape[2], y+2] = ndi.binary_erosion(fiber_ca_cat_2_masked_eroded[i, 0:fiber_ca_cat_2_masked.shape[1], 0:fiber_ca_cat_2_masked.shape[2], y+1], iterations=1) + fiber_border_21[i, 0:fiber_ca_cat_2_masked.shape[1], 0:fiber_ca_cat_2_masked.shape[2]]    

periphery_fiber_cat_2 = np.logical_xor(fiber_ca_cat_2_masked, fiber_ca_cat_2_masked_eroded[0:fiber_ca_cat_2_masked.shape[0], 0:fiber_ca_cat_2_masked.shape[1], 0:fiber_ca_cat_2_masked.shape[2], 41])

imsave("Exp_2_GT_fiber_periphery.tif", periphery_fiber_cat_2.astype(np.uint8))

l_sparks_rand_cat_2 = []

fiber_bn_cat_2_rand = np.zeros(shape=(fiber_ca_cat_2.shape), dtype=np.uint8)

for i in range(fiber_bn_cat_2.shape[0]):
    fiber_bn_cat_2_rand[i] = periphery_fiber_cat_2[i] + fiber_bn_cat_2[i]
    for fiber_ID in np.unique(fiber_bn_cat_2_label[i])[1:]:
        spark_mask_cat_2_1_1 = fiber_bn_cat_2_label[i]==fiber_ID
        if np.sum(fiber_bn_cat_2_rand[i][spark_mask_cat_2_1_1])>np.sum(fiber_bn_cat_2[i][spark_mask_cat_2_1_1]):            
            l_sparks_rand_cat_2.append(np.mean(fiber_ca_cat_2[i][spark_mask_cat_2_1_1]))

del fiber_ca_cat_2_masked_eroded




fiber_border_2_p = border_mask_cat_2 + fiber_ca_predicted_mask_clean_cat_2

fiber_border_21_p = fiber_border_2_p > 1

fiber_ca_cat_2_p_masked_eroded = np.zeros(shape=(fiber_ca_cat_2_masked.shape[0], fiber_ca_cat_2_masked.shape[1], fiber_ca_cat_2_masked.shape[2], 42), dtype=np.uint8)

fiber_ca_cat_2_p_masked_eroded[0:fiber_ca_cat_2_masked.shape[0], 0:fiber_ca_cat_2_masked.shape[1], 0:fiber_ca_cat_2_masked.shape[2], 1] = fiber_ca_predicted_mask_clean_cat_2

for i in range(fiber_ca_cat_2_masked.shape[0]):
    for y in range (40):
        fiber_ca_cat_2_p_masked_eroded[i, 0:fiber_ca_cat_2_masked.shape[1], 0:fiber_ca_cat_2_masked.shape[2], y+2] = ndi.binary_erosion(fiber_ca_cat_2_p_masked_eroded[i, 0:fiber_ca_cat_2_masked.shape[1], 0:fiber_ca_cat_2_masked.shape[2], y+1], iterations=1) + fiber_border_21_p[i, 0:fiber_ca_cat_2_masked.shape[1], 0:fiber_ca_cat_2_masked.shape[2]]    

periphery_fiber_cat_2_p = np.logical_xor(fiber_ca_predicted_mask_clean_cat_2, fiber_ca_cat_2_p_masked_eroded[0:fiber_ca_cat_2_masked.shape[0], 0:fiber_ca_cat_2_masked.shape[1], 0:fiber_ca_cat_2_masked.shape[2], 41])

imsave("Exp_2_P_fiber_periphery.tif", periphery_fiber_cat_2_p.astype(np.uint8))

l_sparks_rand_cat_2_p = []

fiber_bn_cat_2_p_rand = np.zeros(shape=(fiber_sparks_bn_predicted_cat_2.shape), dtype=np.uint8)

for i in range(fiber_sparks_bn_predicted_cat_2.shape[0]):
    fiber_bn_cat_2_p_rand[i] = periphery_fiber_cat_2_p[i] + fiber_sparks_bn_predicted_cat_2[i]    
    for fiber_ID in np.unique(fiber_sparks_bn_predicted_cat_2_label[i])[1:]:
        spark_mask_cat_2_1_1_p = fiber_sparks_bn_predicted_cat_2_label[i]==fiber_ID
        if np.sum(fiber_bn_cat_2_p_rand[i][spark_mask_cat_2_1_1_p])>np.sum(fiber_sparks_bn_predicted_cat_2[i][spark_mask_cat_2_1_1_p]):            
            l_sparks_rand_cat_2_p.append(np.mean(fiber_ca_predicted_cat_2[i][spark_mask_cat_2_1_1_p]))

del fiber_sparks_bn_predicted_cat_2

del fiber_ca_predicted_cat_2

del fiber_ca_cat_2_p_masked_eroded



border_mask_cat_3 = np.zeros(fiber_ca_cat_3_masked.shape, dtype=np.bool)

for i in range(fiber_ca_cat_3_masked.shape[0]):
    border_mask_cat_3[i] = ndi.binary_dilation(border_mask_cat_3[i], border_value=1)

fiber_border_3 = border_mask_cat_3 + fiber_ca_cat_3_masked

fiber_border_31 = fiber_border_3 > 1

fiber_ca_cat_3_masked_eroded = np.zeros(shape=(fiber_ca_cat_3_masked.shape[0], fiber_ca_cat_3_masked.shape[1], fiber_ca_cat_3_masked.shape[2], 42), dtype=np.uint8)

fiber_ca_cat_3_masked_eroded[0:fiber_ca_cat_3_masked.shape[0], 0:fiber_ca_cat_3_masked.shape[1], 0:fiber_ca_cat_3_masked.shape[2], 1] = fiber_ca_cat_3_masked

for i in range(fiber_ca_cat_3_masked.shape[0]):
    for y in range (40):
        fiber_ca_cat_3_masked_eroded[i, 0:fiber_ca_cat_3_masked.shape[1], 0:fiber_ca_cat_3_masked.shape[2], y+2] = ndi.binary_erosion(fiber_ca_cat_3_masked_eroded[i, 0:fiber_ca_cat_3_masked.shape[1], 0:fiber_ca_cat_3_masked.shape[2], y+1], iterations=1) + fiber_border_31[i, 0:fiber_ca_cat_3_masked.shape[1], 0:fiber_ca_cat_3_masked.shape[2]]    

periphery_fiber_cat_3 = np.logical_xor(fiber_ca_cat_3_masked, fiber_ca_cat_3_masked_eroded[0:fiber_ca_cat_3_masked.shape[0], 0:fiber_ca_cat_3_masked.shape[1], 0:fiber_ca_cat_3_masked.shape[2], 41])

imsave("Exp_3_GT_fiber_periphery.tif", periphery_fiber_cat_3.astype(np.uint8))

l_sparks_rand_cat_3 = []

fiber_bn_cat_3_rand = np.zeros(shape=(fiber_ca_cat_3.shape), dtype=np.uint8)

for i in range(fiber_bn_cat_3.shape[0]):
    fiber_bn_cat_3_rand[i] = periphery_fiber_cat_3[i] + fiber_bn_cat_3[i]
    for fiber_ID in np.unique(fiber_bn_cat_3_label[i])[1:]:
        spark_mask_cat_3_1_1 = fiber_bn_cat_3_label[i]==fiber_ID
        if np.sum(fiber_bn_cat_3_rand[i][spark_mask_cat_3_1_1])>np.sum(fiber_bn_cat_3[i][spark_mask_cat_3_1_1]):            
            l_sparks_rand_cat_3.append(np.mean(fiber_ca_cat_3[i][spark_mask_cat_3_1_1]))

del fiber_ca_cat_3_masked_eroded



fiber_border_3_p = border_mask_cat_3 + fiber_ca_predicted_mask_clean_cat_3

fiber_border_31_p = fiber_border_3_p > 1

fiber_ca_cat_3_p_masked_eroded = np.zeros(shape=(fiber_ca_cat_3_masked.shape[0], fiber_ca_cat_3_masked.shape[1], fiber_ca_cat_3_masked.shape[2], 42), dtype=np.uint8)

fiber_ca_cat_3_p_masked_eroded[0:fiber_ca_cat_3_masked.shape[0], 0:fiber_ca_cat_3_masked.shape[1], 0:fiber_ca_cat_3_masked.shape[2], 1] = fiber_ca_predicted_mask_clean_cat_3

for i in range(fiber_ca_cat_3_masked.shape[0]):
    for y in range (40):
        fiber_ca_cat_3_p_masked_eroded[i, 0:fiber_ca_cat_3_masked.shape[1], 0:fiber_ca_cat_3_masked.shape[2], y+2] = ndi.binary_erosion(fiber_ca_cat_3_p_masked_eroded[i, 0:fiber_ca_cat_3_masked.shape[1], 0:fiber_ca_cat_3_masked.shape[2], y+1], iterations=1) + fiber_border_31_p[i, 0:fiber_ca_cat_3_masked.shape[1], 0:fiber_ca_cat_3_masked.shape[2]]    

periphery_fiber_cat_3_p = np.logical_xor(fiber_ca_predicted_mask_clean_cat_3, fiber_ca_cat_3_p_masked_eroded[0:fiber_ca_cat_3_masked.shape[0], 0:fiber_ca_cat_3_masked.shape[1], 0:fiber_ca_cat_3_masked.shape[2], 41])

imsave("Exp_3_P_fiber_periphery.tif", periphery_fiber_cat_3_p.astype(np.uint8))

l_sparks_rand_cat_3_p = []

fiber_bn_cat_3_p_rand = np.zeros(shape=(fiber_sparks_bn_predicted_cat_3.shape), dtype=np.uint8)

for i in range(fiber_sparks_bn_predicted_cat_3.shape[0]):
    fiber_bn_cat_3_p_rand[i] = periphery_fiber_cat_3_p[i] + fiber_sparks_bn_predicted_cat_3[i]    
    for fiber_ID in np.unique(fiber_sparks_bn_predicted_cat_3_label[i])[1:]:
        spark_mask_cat_3_1_1_p = fiber_sparks_bn_predicted_cat_3_label[i]==fiber_ID
        if np.sum(fiber_bn_cat_3_p_rand[i][spark_mask_cat_3_1_1_p])>np.sum(fiber_sparks_bn_predicted_cat_3[i][spark_mask_cat_3_1_1_p]):            
            l_sparks_rand_cat_3_p.append(np.mean(fiber_ca_predicted_cat_3[i][spark_mask_cat_3_1_1_p]))

del fiber_sparks_bn_predicted_cat_3

del fiber_ca_predicted_cat_3

del fiber_ca_cat_3_p_masked_eroded



border_mask_train = np.zeros(fiber_ca_masked_train.shape, dtype=np.bool)

for i in range(fiber_ca_masked_train.shape[0]):
    border_mask_train[i] = ndi.binary_dilation(border_mask_train[i], border_value=1)

fiber_border_t = border_mask_train + fiber_ca_masked_train

fiber_border_t1 = fiber_border_t > 1

fiber_ca_train_masked_eroded = np.zeros(shape=(fiber_ca_masked_train.shape[0], fiber_ca_masked_train.shape[1], fiber_ca_masked_train.shape[2], 42), dtype=np.uint8)

fiber_ca_train_masked_eroded[0:fiber_ca_masked_train.shape[0], 0:fiber_ca_masked_train.shape[1], 0:fiber_ca_masked_train.shape[2], 1] = fiber_ca_masked_train

for i in range(fiber_ca_masked_train.shape[0]):
    for y in range (40):
        fiber_ca_train_masked_eroded[i, 0:fiber_ca_masked_train.shape[1], 0:fiber_ca_masked_train.shape[2], y+2] = ndi.binary_erosion(fiber_ca_train_masked_eroded[i, 0:fiber_ca_masked_train.shape[1], 0:fiber_ca_masked_train.shape[2], y+1], iterations=1) + fiber_border_t1[i, 0:fiber_ca_masked_train.shape[1], 0:fiber_ca_masked_train.shape[2]]    

periphery_fiber_train = np.logical_xor(fiber_ca_masked_train, fiber_ca_train_masked_eroded[0:fiber_ca_masked_train.shape[0], 0:fiber_ca_masked_train.shape[1], 0:fiber_ca_masked_train.shape[2], 41])

imsave("train_GT_fiber_periphery.tif", periphery_fiber_train.astype(np.uint8))

l_sparks_rand_train = []

fiber_bn_train_rand = np.zeros(shape=(fiber_ca_train.shape), dtype=np.uint8)

for i in range(fiber_bn_train.shape[0]):
    fiber_bn_train_rand[i] = periphery_fiber_train[i] + fiber_bn_train[i]
    for fiber_ID in np.unique(fiber_sparks_bn_training_label[i])[1:]:
        spark_mask_cat_t_1_1 = fiber_sparks_bn_training_label[i]==fiber_ID
        if np.sum(fiber_bn_train_rand[i][spark_mask_cat_t_1_1])>np.sum(fiber_bn_train[i][spark_mask_cat_t_1_1]):            
            l_sparks_rand_train.append(np.mean(fiber_ca_train[i][spark_mask_cat_t_1_1]))

l_spark_frequency_border = [100*len(l_sparks_rand_cat_1)/len(l_spark_frequency_cat_1), 100*len(l_sparks_rand_cat_1_p)/len(l_spark_frequency_cat_1_predicted), 100*len(l_sparks_rand_cat_2)/len(l_spark_frequency_cat_2), 100*len(l_sparks_rand_cat_2_p)/len(l_spark_frequency_cat_2_predicted), 100*len(l_sparks_rand_cat_3)/len(l_spark_frequency_cat_3), 100*len(l_sparks_rand_cat_3_p)/len(l_spark_frequency_cat_3_predicted), 100*len(l_sparks_rand_train)/len(l_spark_frequency_train)]

plt.figure(figsize=(450,150))

names = ["Exp 1 GT", "Exp 1 P", "Exp 2 GT", "Exp 2 P", "Exp 3 GT", "Exp 3 P", "Training"]
data = l_spark_frequency_border
fig, ax = plt.subplots()
ax.set_ylabel('Sparks at the border (%)')
bp = ax.bar(names, data, alpha=0.7)
plt.savefig('Spark frequency border.pdf', format='pdf')
plt.show()