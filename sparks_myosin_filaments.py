import numpy as np
import scipy.ndimage as ndi
import matplotlib.pyplot as plt
from tifffile import imread
from tifffile import imsave
from skimage import filters
from skimage.filters import rank
from skimage.morphology import disk

#input of the algorithm: Exp_1_P_bn_sparks.tif, Exp_2_P_bn_sparks.tif, Exp_3_P_bn_sparks.tif, Exp_1_GT_bn_sparks.tif, Exp_2_GT_bn_sparks.tif, Exp_3_GT_bn_sparks.tif, train_GT_bn_sparks.tif,
#Exp_1_SHG.tif, Exp_2_SHG.tif, Exp_3_SHG.tif, train_SHG.tif
#output of the algorithm: Exp_1_SHG_mean.tif, Exp_1_SHG_mean_bn.tif, Exp_2_SHG_mean.tif, Exp_2_SHG_mean_bn.tif, Exp_3_SHG_mean.tif, Exp_3_SHG_mean_bn.tif,
#train_SHG_mean.tif, train_SHG_mean_bn.tif, Spark area over mayosin filaments spark area with training boxplot.pdf
#the files must be placed in the same folder as the .py file but this can be modified
#the approach to read out data is similar to the one described at https://github.com/embl-bio-it/image-analysis-with-python/blob/master/session-3to5/image_analysis_tutorial_solutions.ipynb
#the module skimage.measure.regionprops_table can be used instead
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



fiber_shg_cat_1_gauss = np.zeros(shape=(fiber_shg_cat_1.shape), dtype=np.uint32)
fiber_gauss_bn_shg_cat_1 = np.zeros(shape=(fiber_shg_cat_1.shape), dtype=np.uint8)
fiber_shg_cat_1_limit = np.zeros(shape=(fiber_shg_cat_1.shape), dtype=np.uint32)

fiber_mean_shg_cat_1 = np.zeros(shape=(fiber_shg_cat_1.shape), dtype=np.uint8)

for i in range(fiber_shg_cat_1.shape[0]):
    fiber_mean_shg_cat_1[i] = np.mean(fiber_shg_cat_1[0:fiber_shg_cat_1.shape[0]], axis=0)

imsave("Exp_1_SHG_mean.tif", fiber_mean_shg_cat_1.astype(np.uint8))

for i in range(fiber_shg_cat_1.shape[0]):
    fiber_shg_cat_1_gauss[i] = 1000*filters.gaussian(fiber_mean_shg_cat_1[i], sigma=2)

fiber_shg_mask_1_cat_1 = np.zeros(shape=(fiber_shg_cat_1.shape), dtype=np.uint32)

for i in range(fiber_shg_cat_1.shape[0]):
    fiber_shg_mask_1_cat_1[i] = fiber_shg_cat_1_gauss[i]>1*np.std(fiber_shg_cat_1_gauss[i])

for i in range(fiber_shg_cat_1.shape[0]):
    fiber_shg_cat_1_limit[i] = rank.mean(fiber_mean_shg_cat_1[i], disk(5))

fiber_mean_bn_shg_cat_1 = fiber_shg_mask_1_cat_1*fiber_mean_shg_cat_1 > 0.95*fiber_shg_mask_1_cat_1*fiber_shg_cat_1_limit

imsave("Exp_1_SHG_mean_bn.tif", fiber_mean_bn_shg_cat_1.astype(np.uint8))

##########################################################################################
##########################################################################################
##########################################################################################

fiber_shg_cat_2_gauss = np.zeros(shape=(fiber_shg_cat_2.shape), dtype=np.uint32)
fiber_gauss_bn_shg_cat_2 = np.zeros(shape=(fiber_shg_cat_2.shape), dtype=np.uint8)
fiber_shg_cat_2_limit = np.zeros(shape=(fiber_shg_cat_2.shape), dtype=np.uint32)
fiber_mean_shg_cat_2 = np.zeros(shape=(fiber_shg_cat_2.shape), dtype=np.uint8)

for i in range(fiber_shg_cat_2.shape[0]):
    fiber_mean_shg_cat_2[i] = np.mean(fiber_shg_cat_2[0:fiber_shg_cat_2.shape[0]], axis=0)

imsave("Exp_2_SHG_mean.tif", fiber_mean_shg_cat_2.astype(np.uint8))

for i in range(fiber_shg_cat_2.shape[0]):
    fiber_shg_cat_2_gauss[i] = 1000*filters.gaussian(fiber_mean_shg_cat_2[i], sigma=2)

fiber_shg_mask_1_cat_2 = np.zeros(shape=(fiber_shg_cat_2.shape), dtype=np.uint32)

for i in range(fiber_shg_cat_2.shape[0]):
    fiber_shg_mask_1_cat_2[i] = fiber_shg_cat_2_gauss[i]>1*np.std(fiber_shg_cat_2_gauss[i])

for i in range(fiber_shg_cat_2.shape[0]):
    fiber_shg_cat_2_limit[i] = rank.mean(fiber_mean_shg_cat_2[i], disk(5))

fiber_mean_bn_shg_cat_2 = fiber_shg_mask_1_cat_2*fiber_mean_shg_cat_2 > 0.95*fiber_shg_mask_1_cat_2*fiber_shg_cat_2_limit

imsave("Exp_2_SHG_mean_bn.tif", fiber_mean_bn_shg_cat_2.astype(np.uint8))

##########################################################################################
##########################################################################################
##########################################################################################

fiber_shg_cat_3_gauss = np.zeros(shape=(fiber_shg_cat_3.shape), dtype=np.uint32)
fiber_gauss_bn_shg_cat_3 = np.zeros(shape=(fiber_shg_cat_3.shape), dtype=np.uint8)
fiber_shg_cat_3_limit = np.zeros(shape=(fiber_shg_cat_3.shape), dtype=np.uint32)
fiber_mean_shg_cat_3 = np.zeros(shape=(fiber_shg_cat_3.shape), dtype=np.uint8)

for i in range(fiber_shg_cat_3.shape[0]):
    fiber_mean_shg_cat_3[i] = np.mean(fiber_shg_cat_3[0:fiber_shg_cat_3.shape[0]], axis=0)

imsave("Exp_3_SHG_mean.tif", fiber_mean_shg_cat_3.astype(np.uint8))

for i in range(fiber_shg_cat_3.shape[0]):
    fiber_shg_cat_3_gauss[i] = 1000*filters.gaussian(fiber_mean_shg_cat_3[i], sigma=2)

fiber_shg_mask_1_cat_3 = np.zeros(shape=(fiber_shg_cat_3.shape), dtype=np.uint32)

for i in range(fiber_shg_cat_3.shape[0]):
    fiber_shg_mask_1_cat_3[i] = fiber_shg_cat_3_gauss[i]>1*np.std(fiber_shg_cat_3_gauss[i])

for i in range(fiber_shg_cat_3.shape[0]):
    fiber_shg_cat_3_limit[i] = rank.mean(fiber_mean_shg_cat_3[i], disk(5))

fiber_mean_bn_shg_cat_3 = fiber_shg_mask_1_cat_3*fiber_mean_shg_cat_3 > 0.95*fiber_shg_mask_1_cat_3*fiber_shg_cat_3_limit

imsave("Exp_3_SHG_mean_bn.tif", fiber_mean_bn_shg_cat_3.astype(np.uint8))

##########################################################################################
##########################################################################################
##########################################################################################

fiber_shg_train_gauss = np.zeros(shape=(fiber_shg_train.shape), dtype=np.uint32)
fiber_gauss_bn_shg_train = np.zeros(shape=(fiber_shg_train.shape), dtype=np.uint8)
fiber_shg_train_limit = np.zeros(shape=(fiber_shg_train.shape), dtype=np.uint32)
fiber_mean_shg_train = np.zeros(shape=(fiber_shg_train.shape), dtype=np.uint8)

for i in range(fiber_shg_train.shape[0]):
    fiber_mean_shg_train [i] = np.mean(fiber_shg_train[0:fiber_shg_train.shape[0]], axis=0)

imsave("train_SHG_mean.tif", fiber_mean_shg_train.astype(np.uint8))

for i in range(fiber_shg_train.shape[0]):
    fiber_shg_train_gauss[i] = 1000*filters.gaussian(fiber_mean_shg_train[i], sigma=2)

fiber_shg_mask_1_train = np.zeros(shape=(fiber_shg_train.shape), dtype=np.uint32)

for i in range(fiber_shg_train.shape[0]):
    fiber_shg_mask_1_train[i] = fiber_shg_train_gauss[i]>1*np.std(fiber_shg_train_gauss[i])

for i in range(fiber_shg_train.shape[0]):
    fiber_shg_train_limit[i] = rank.mean(fiber_mean_shg_train[i], disk(5))

fiber_mean_bn_shg_train = fiber_shg_mask_1_train*fiber_mean_shg_train > 0.95*fiber_shg_mask_1_train*fiber_shg_train_limit

imsave("train_SHG_mean_bn.tif", fiber_mean_bn_shg_train.astype(np.uint8))

##########################################################################################
##########################################################################################
##########################################################################################

stack_bn_cat_1_shg = fiber_mean_bn_shg_cat_1 + stack_bn_cat_1
stack_bn_cat_2_shg = fiber_mean_bn_shg_cat_2 + stack_bn_cat_2
stack_bn_cat_3_shg = fiber_mean_bn_shg_cat_3 + stack_bn_cat_3
stack_bn_train_shg = fiber_mean_bn_shg_train + stack_bn_train

fiber_sparks_bn_predicted_cat_1_shg = fiber_mean_bn_shg_cat_1 + fiber_sparks_bn_predicted_cat_1
fiber_sparks_bn_predicted_cat_2_shg = fiber_mean_bn_shg_cat_2 + fiber_sparks_bn_predicted_cat_2
fiber_sparks_bn_predicted_cat_3_shg = fiber_mean_bn_shg_cat_3 + fiber_sparks_bn_predicted_cat_3

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
  
    
results = {"spark_area_striation_pattern_cat_1"      : [],
           "spark_area_striation_pattern_cat_1_predicted"     : [],
           "spark_area_striation_pattern_cat_2"    : [],
           "spark_area_striation_pattern_cat_2_predicted"    : [],
           "spark_area_striation_pattern_cat_3" : [],
           "spark_area_striation_pattern_cat_3_predicted"    : [],
           "spark_area_striation_pattern_train" : []}    
    
for i in range(stack_bn_cat_1.shape[0]):
    for spark_ID in np.unique(stack_bn_cat_1_label[i])[1:]:
        spark_mask_cat_1 = stack_bn_cat_1_label[i]==spark_ID
        results["spark_area_striation_pattern_cat_1"].append((np.sum(stack_bn_cat_1_shg[i][spark_mask_cat_1])-np.sum(spark_mask_cat_1))/np.sum(spark_mask_cat_1))


for i in range(fiber_sparks_bn_predicted_cat_1.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_1_label[i])[1:]:
        spark_mask_cat_1_predicted = fiber_sparks_bn_predicted_cat_1_label[i]==spark_ID
        results["spark_area_striation_pattern_cat_1_predicted"].append((np.sum(fiber_sparks_bn_predicted_cat_1_shg[i][spark_mask_cat_1_predicted])-np.sum(spark_mask_cat_1_predicted))/np.sum(spark_mask_cat_1_predicted))


for i in range(stack_bn_cat_2.shape[0]):
    for spark_ID in np.unique(stack_bn_cat_2_label[i])[1:]:
        spark_mask_cat_2 = stack_bn_cat_2_label[i]==spark_ID
        results["spark_area_striation_pattern_cat_2"].append((np.sum(stack_bn_cat_2_shg[i][spark_mask_cat_2])-np.sum(spark_mask_cat_2))/np.sum(spark_mask_cat_2))


for i in range(fiber_sparks_bn_predicted_cat_2.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_2_label[i])[1:]:
        spark_mask_cat_2_predicted = fiber_sparks_bn_predicted_cat_2_label[i]==spark_ID
        results["spark_area_striation_pattern_cat_2_predicted"].append((np.sum(fiber_sparks_bn_predicted_cat_2_shg[i][spark_mask_cat_2_predicted])-np.sum(spark_mask_cat_2_predicted))/np.sum(spark_mask_cat_2_predicted))


for i in range(stack_bn_cat_3.shape[0]):
    for spark_ID in np.unique(stack_bn_cat_3_label[i])[1:]:
        spark_mask_cat_3 = stack_bn_cat_3_label[i]==spark_ID
        results["spark_area_striation_pattern_cat_3"].append((np.sum(stack_bn_cat_3_shg[i][spark_mask_cat_3])-np.sum(spark_mask_cat_3))/np.sum(spark_mask_cat_3))


for i in range(fiber_sparks_bn_predicted_cat_3.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_3_label[i])[1:]:
        spark_mask_cat_3_predicted = fiber_sparks_bn_predicted_cat_3_label[i]==spark_ID
        results["spark_area_striation_pattern_cat_3_predicted"].append((np.sum(fiber_sparks_bn_predicted_cat_3_shg[i][spark_mask_cat_3_predicted])-np.sum(spark_mask_cat_3_predicted))/np.sum(spark_mask_cat_3_predicted))


for i in range(stack_bn_train.shape[0]):
    for spark_ID in np.unique(stack_bn_train_label[i])[1:]:
        spark_mask_train = stack_bn_train_label[i]==spark_ID
        results["spark_area_striation_pattern_train"].append((np.sum(stack_bn_train_shg[i][spark_mask_train])-np.sum(spark_mask_train))/np.sum(spark_mask_train))

plt.figure(figsize=(450,150))

data = [results["spark_area_striation_pattern_cat_1"], results["spark_area_striation_pattern_cat_1_predicted"], results["spark_area_striation_pattern_cat_2"], results["spark_area_striation_pattern_cat_2_predicted"], results["spark_area_striation_pattern_cat_3"], results["spark_area_striation_pattern_cat_3_predicted"], results["spark_area_striation_pattern_train"]]
fig, ax = plt.subplots()
ax.set_ylabel('Spark area over striation pattern/spark area')
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

plt.xticks([1,2,3,4,5,6, 7], ["Exp 1 GT", "Exp 1 P", "Exp 2 GT", "Exp 2 P", "Exp 3 GT", "Exp 3 P", "Training"])
plt.savefig("Spark area over mayosin filaments spark area with training boxplot.pdf", format='pdf')
plt.show()