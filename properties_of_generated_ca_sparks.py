import numpy as np
import scipy.ndimage as ndi
import matplotlib.pyplot as plt
from tifffile import imread
from tifffile import imsave
from skimage import filters

#input of the algorithm: Exp_1_P.tif, Exp_1_GT.tif, Exp_2_P.tif, Exp_2_GT.tif, Exp_3_P.tif, Exp_3_GT.tif, train_GT.tif, Exp_1_P_masked, Exp_1_GT_masked.tif,
#Exp_2_P_masked.tif, Exp_2_GT_masked.tif, Exp_3_P_masked.tif, Exp_3_GT_masked.tif, train_GT_masked.tif, Exp_1_GT_bn_sparks.tif, Exp_2_GT_bn_sparks.tif,
#Exp_3_GT_bn_sparks.tif, train_GT_bn_sparks.tif
#output of the algorithm: Exp_1_P_bn_sparks.tif, Exp_2_P_bn_sparks.tif, Exp_3_P_bn_sparks.tif, Spark area F1.pdf, Spark area violin plot box plot.pdf, Spark mean intensity violin plot box plot.pdf, Spark frequency.pdf
#the files must be placed in the same folder as the .py file but this can be modified
#the approach used for size filters is similar to the one described at https://github.com/embl-bio-it/image-analysis-with-python/blob/master/session-3to5/image_analysis_tutorial_solutions.ipynb
#instead of the size filters used here the module skimage.morphology.remove_small_objects of scikit image can be applied as described in lines 63-65 of the image_pre_processing algorithm
#the approach to read out data is similar to the one described at https://github.com/embl-bio-it/image-analysis-with-python/blob/master/session-3to5/image_analysis_tutorial_solutions.ipynb
#the module skimage.measure.regionprops_table can be used instead
fiber_ca_predicted_cat_1 = imread("Exp_1_P.tif")
fiber_ca_cat_1 = imread("Exp_1_GT.tif")
fiber_bn_cat_1 = imread("Exp_1_GT_bn_sparks.tif")
fiber_ca_predicted_mask_clean_cat_1 = imread("Exp_1_P_masked.tif")
fiber_ca_cat_1_masked = imread("Exp_1_GT_masked.tif")

fiber_ca_predicted_cat_1_gauss = np.zeros(shape=(fiber_ca_predicted_cat_1.shape), dtype=np.float32)
fiber_gauss_sparks_bn_predicted_cat_1 = np.zeros(shape=(fiber_ca_predicted_cat_1.shape), dtype=np.uint8)


for i in range(fiber_ca_predicted_cat_1.shape[0]):
    fiber_ca_predicted_cat_1_gauss[i] = filters.gaussian(fiber_ca_predicted_cat_1[i], sigma=2)

fiber_ca_predicted_cat_1_gauss_masked = np.zeros(shape=(fiber_ca_predicted_cat_1.shape), dtype=np.float32)

for i in range(fiber_ca_predicted_cat_1.shape[0]):
    fiber_ca_predicted_cat_1_gauss_masked[i] = fiber_ca_predicted_cat_1_gauss[i]*fiber_ca_predicted_mask_clean_cat_1[i]


for i in range(fiber_ca_predicted_cat_1.shape[0]):
    for fiber_ID in np.unique(fiber_ca_predicted_mask_clean_cat_1[i])[1:]:
        fiber_mask_cat_1 = fiber_ca_predicted_mask_clean_cat_1[i]==fiber_ID
        fiber_gauss_sparks_bn_predicted_cat_1[i] = fiber_ca_predicted_cat_1_gauss_masked[i] > (np.mean(fiber_ca_predicted_cat_1_gauss_masked[i][fiber_mask_cat_1]) + 3.6*np.std(fiber_ca_predicted_cat_1_gauss_masked[i][fiber_mask_cat_1]))       

fiber_gauss_sparks_bn_label_predicted_cat_1 = np.zeros(shape=(fiber_ca_predicted_cat_1.shape), dtype=np.uint32)        

for i in range(fiber_ca_predicted_cat_1.shape[0]):
    fiber_gauss_sparks_bn_label_predicted_cat_1[i], fiber_gauss_sparks_bn_label_predicted_cat_1_number = ndi.label(fiber_gauss_sparks_bn_predicted_cat_1[i])
    
fiber_gauss_sparks_bn_label_clean_predicted_cat_1 = np.copy(fiber_gauss_sparks_bn_label_predicted_cat_1)

for i in range(fiber_ca_predicted_cat_1.shape[0]):
    for spark_ID in np.unique(fiber_gauss_sparks_bn_label_predicted_cat_1[i]):
        fiber_mask_sparks_predicted_cat_1 = fiber_gauss_sparks_bn_label_predicted_cat_1[i]==spark_ID
        if np.sum(fiber_mask_sparks_predicted_cat_1) < 55:
            fiber_gauss_sparks_bn_label_clean_predicted_cat_1[i][fiber_mask_sparks_predicted_cat_1] = 0

fiber_sparks_bn_predicted_cat_1 = fiber_gauss_sparks_bn_label_clean_predicted_cat_1 > 0

imsave("Exp_1_P_bn_sparks.tif", fiber_sparks_bn_predicted_cat_1.astype(np.uint8))



fiber_ca_predicted_cat_2 = imread("Exp_2_P.tif")
fiber_ca_cat_2 = imread("Exp_2_GT.tif")
fiber_bn_cat_2 = imread("Exp_2_GT_bn_sparks.tif")
fiber_ca_predicted_mask_clean_cat_2 = imread("Exp_2_P_masked.tif")
fiber_ca_cat_2_masked = imread("Exp_2_GT_masked.tif")

fiber_ca_predicted_cat_2_gauss = np.zeros(shape=(fiber_ca_predicted_cat_2.shape), dtype=np.float32)
fiber_gauss_sparks_bn_predicted_cat_2 = np.zeros(shape=(fiber_ca_predicted_cat_2.shape), dtype=np.uint8)


for i in range(fiber_ca_predicted_cat_2.shape[0]):
    fiber_ca_predicted_cat_2_gauss[i] = filters.gaussian(fiber_ca_predicted_cat_2[i], sigma=2)

fiber_ca_predicted_cat_2_gauss_masked = np.zeros(shape=(fiber_ca_predicted_cat_2.shape), dtype=np.float32)

for i in range(fiber_ca_predicted_cat_2.shape[0]):
    fiber_ca_predicted_cat_2_gauss_masked[i] = fiber_ca_predicted_cat_2_gauss[i]*fiber_ca_predicted_mask_clean_cat_2[i]


for i in range(fiber_ca_predicted_cat_2.shape[0]):
    for fiber_ID in np.unique(fiber_ca_predicted_mask_clean_cat_2[i])[1:]:
        fiber_mask_cat_2 = fiber_ca_predicted_mask_clean_cat_2[i]==fiber_ID
        fiber_gauss_sparks_bn_predicted_cat_2[i] = fiber_ca_predicted_cat_2_gauss_masked[i] > (np.mean(fiber_ca_predicted_cat_2_gauss_masked[i][fiber_mask_cat_2]) + 3.6*np.std(fiber_ca_predicted_cat_2_gauss_masked[i][fiber_mask_cat_2]))

fiber_gauss_sparks_bn_label_predicted_cat_2 = np.zeros(shape=(fiber_ca_predicted_cat_2.shape), dtype=np.uint32)        

for i in range(fiber_ca_predicted_cat_2.shape[0]):
    fiber_gauss_sparks_bn_label_predicted_cat_2[i], fiber_gauss_sparks_bn_label_predicted_cat_2_number = ndi.label(fiber_gauss_sparks_bn_predicted_cat_2[i])
    
fiber_gauss_sparks_bn_label_clean_predicted_cat_2 = np.copy(fiber_gauss_sparks_bn_label_predicted_cat_2)

for i in range(fiber_ca_predicted_cat_2.shape[0]):
    for spark_ID in np.unique(fiber_gauss_sparks_bn_label_predicted_cat_2[i]):
        fiber_mask_sparks_predicted_cat_2 = fiber_gauss_sparks_bn_label_predicted_cat_2[i]==spark_ID
        if np.sum(fiber_mask_sparks_predicted_cat_2) < 55:
            fiber_gauss_sparks_bn_label_clean_predicted_cat_2[i][fiber_mask_sparks_predicted_cat_2] = 0

fiber_sparks_bn_predicted_cat_2 = fiber_gauss_sparks_bn_label_clean_predicted_cat_2 > 0

imsave("Exp_2_P_bn_sparks.tif", fiber_sparks_bn_predicted_cat_2.astype(np.uint8))



fiber_ca_predicted_cat_3 = imread("Exp_3_P.tif")
fiber_ca_cat_3 = imread("Exp_3_GT.tif")
fiber_bn_cat_3 = imread("Exp_3_GT_bn_sparks.tif")
fiber_ca_predicted_mask_clean_cat_3 = imread("Exp_3_P_masked.tif")
fiber_ca_cat_3_masked = imread("Exp_3_GT_masked.tif")


fiber_ca_predicted_cat_3_gauss = np.zeros(shape=(fiber_ca_predicted_cat_3.shape), dtype=np.float32)
fiber_gauss_sparks_bn_predicted_cat_3 = np.zeros(shape=(fiber_ca_predicted_cat_3.shape), dtype=np.uint8)


for i in range(fiber_ca_predicted_cat_3.shape[0]):
    fiber_ca_predicted_cat_3_gauss[i] = filters.gaussian(fiber_ca_predicted_cat_3[i], sigma=2)

fiber_ca_predicted_cat_3_gauss_masked = np.zeros(shape=(fiber_ca_predicted_cat_3.shape), dtype=np.float32)

for i in range(fiber_ca_predicted_cat_3.shape[0]):
    fiber_ca_predicted_cat_3_gauss_masked[i] = fiber_ca_predicted_cat_3_gauss[i]*fiber_ca_predicted_mask_clean_cat_3[i]


for i in range(fiber_ca_predicted_cat_3.shape[0]):
    for fiber_ID in np.unique(fiber_ca_predicted_mask_clean_cat_3[i])[1:]:
        fiber_mask_cat_3 = fiber_ca_predicted_mask_clean_cat_3[i]==fiber_ID
        fiber_gauss_sparks_bn_predicted_cat_3[i] = fiber_ca_predicted_cat_3_gauss_masked[i] > (np.mean(fiber_ca_predicted_cat_3_gauss_masked[i][fiber_mask_cat_3]) + 3.6*np.std(fiber_ca_predicted_cat_3_gauss_masked[i][fiber_mask_cat_3]))

fiber_gauss_sparks_bn_label_predicted_cat_3 = np.zeros(shape=(fiber_ca_predicted_cat_3.shape), dtype=np.uint32)        

for i in range(fiber_ca_predicted_cat_3.shape[0]):
    fiber_gauss_sparks_bn_label_predicted_cat_3[i], fiber_gauss_sparks_bn_label_predicted_cat_3_label = ndi.label(fiber_gauss_sparks_bn_predicted_cat_3[i])
    
fiber_gauss_sparks_bn_label_clean_predicted_cat_3 = np.copy(fiber_gauss_sparks_bn_label_predicted_cat_3)

for i in range(fiber_ca_predicted_cat_3.shape[0]):
    for spark_ID in np.unique(fiber_gauss_sparks_bn_label_predicted_cat_3[i]):
        fiber_mask_sparks_predicted_cat_3 = fiber_gauss_sparks_bn_label_predicted_cat_3[i]==spark_ID
        if np.sum(fiber_mask_sparks_predicted_cat_3) < 55:
            fiber_gauss_sparks_bn_label_clean_predicted_cat_3[i][fiber_mask_sparks_predicted_cat_3] = 0

fiber_sparks_bn_predicted_cat_3 = fiber_gauss_sparks_bn_label_clean_predicted_cat_3 > 0

imsave("Exp_3_P_bn_sparks.tif", fiber_sparks_bn_predicted_cat_3.astype(np.uint8))

############################################
data_1 = fiber_bn_cat_1 + fiber_sparks_bn_predicted_cat_1

data_gleich_2_1 = data_1 > 1

data_gleich_1_1 = data_1 == 1
############################################
data_2 = fiber_bn_cat_2 + fiber_sparks_bn_predicted_cat_2

data_gleich_2_2 = data_2 > 1

data_gleich_1_2 = data_2 == 1
############################################
data_3 = fiber_bn_cat_3 + fiber_sparks_bn_predicted_cat_3

data_gleich_2_3 = data_3 > 1

data_gleich_1_3 = data_3 == 1
############################################
l = [2*np.sum(data_gleich_2_1)/(2*np.sum(data_gleich_2_1)+np.sum(data_gleich_1_1)), 2*np.sum(data_gleich_2_2)/(2*np.sum(data_gleich_2_2)+np.sum(data_gleich_1_2)), 2*np.sum(data_gleich_2_3)/(2*np.sum(data_gleich_2_3)+np.sum(data_gleich_1_3))]

plt.figure(figsize=(450,150))

names = ["Exp 1", "Exp 2", "Exp 3"]
data = l
fig, ax = plt.subplots()
plt.tick_params(axis="both", labelsize=12)
plt.ylabel('F1 score',fontsize=12)
bp = ax.bar(names, data, alpha=0.7, width=0.5, color="black")
plt.savefig('Spark area F1.pdf', format='pdf')
plt.show()

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

plt.figure(figsize=(450,150))

data = [results["spark_area_cat_1"], results["spark_area_cat_1_predicted"], results["spark_area_cat_2"], results["spark_area_cat_2_predicted"], results["spark_area_cat_3"], results["spark_area_cat_3_predicted"], results["spark_area_train"]]
fig, ax = plt.subplots()
ax.set_ylabel('Spark area in \u03BC$m^2$')
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
plt.savefig('Spark area violin plot box plot.pdf', format='pdf')
plt.show()

plt.figure(figsize=(450,150))

data = [results["spark_mean_int_cat_1"], results["spark_mean_int_cat_1_predicted"], results["spark_mean_int_cat_2"], results["spark_mean_int_cat_2_predicted"], results["spark_mean_int_cat_3"], results["spark_mean_int_cat_3_predicted"], results["spark_mean_int_train"]]
fig, ax = plt.subplots()
ax.set_ylabel('Spark mean intensity')
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
plt.savefig('Spark mean intensity violin plot box plot.pdf', format='pdf')
plt.show()

l_spark_frequency_cat_1 = results["spark_area_cat_1"]
l_spark_frequency_cat_1_predicted = results["spark_area_cat_1_predicted"]
l_spark_frequency_cat_2 = results["spark_area_cat_2"]
l_spark_frequency_cat_2_predicted = results["spark_area_cat_2_predicted"]
l_spark_frequency_cat_3 = results["spark_area_cat_3"]
l_spark_frequency_cat_3_predicted = results["spark_area_cat_3_predicted"]
l_spark_frequency_train = results["spark_area_train"]

l_spark_frequency = [10000*len(l_spark_frequency_cat_1)/(np.sum(fiber_ca_cat_1_masked)*0.2325149*0.2325149), 10000*len(l_spark_frequency_cat_1_predicted)/(np.sum(fiber_ca_predicted_mask_clean_cat_1)*0.2325149*0.2325149), 10000*len(l_spark_frequency_cat_2)/(np.sum(fiber_ca_cat_2_masked)*0.2325149*0.2325149), 10000*len(l_spark_frequency_cat_2_predicted)/(np.sum(fiber_ca_predicted_mask_clean_cat_2)*0.2325149*0.2325149), 10000*len(l_spark_frequency_cat_3)/(np.sum(fiber_ca_cat_3_masked)*0.2325149*0.2325149), 10000*len(l_spark_frequency_cat_3_predicted)/(np.sum(fiber_ca_predicted_mask_clean_cat_3)*0.2325149*0.2325149), 10000*len(l_spark_frequency_train)/(np.sum(fiber_ca_masked_train)*0.2325149*0.2325149)]

plt.figure(figsize=(450,150))
names = ["Exp 1 GT", "Exp 1 P", "Exp 2 GT", "Exp 2 P", "Exp 3 GT", "Exp 3 P", "Training"]
data = l_spark_frequency
fig, ax = plt.subplots()
ax.set_ylabel('Spark frequency in sparks/(10000*\u03BC$m^2$)')
bp = ax.bar(names, data, alpha=0.7)
plt.savefig('Spark frequency test um.pdf', format='pdf')
plt.show()