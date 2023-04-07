import numpy as np
import scipy.ndimage as ndi
import matplotlib.pyplot as plt
from tifffile import imread
from tifffile import imsave
from skimage import filters

#input of the algorithm: Exp_1_P.tif, Exp_1_GT.tif, Exp_2_P.tif, Exp_2_GT.tif, Exp_3_P.tif, Exp_3_GT.tif, train_GT.tif
#output of the algorithm: Exp_1_P_masked, Exp_1_GT_masked.tif, Exp_2_P_masked.tif, Exp_2_GT_masked.tif, Exp_3_P_masked.tif, Exp_3_GT_masked.tif, train_GT_masked.tif, Fiber mean intesity violin plot box plot.pdf, Fiber area violin plot box plot.pdf, Fiber area F1.pdf
#the files must be placed in the same folder as the .py file but this can be modified
#the approach used for size filters is similar to the one described at https://github.com/embl-bio-it/image-analysis-with-python/blob/master/session-3to5/image_analysis_tutorial_solutions.ipynb
#instead of the size filters used here the module skimage.morphology.remove_small_objects of scikit image can be applied as described in lines 63-65 of the image_pre_processing algorithm
#the approach to read out data is similar to the one described at https://github.com/embl-bio-it/image-analysis-with-python/blob/master/session-3to5/image_analysis_tutorial_solutions.ipynb
#the module skimage.measure.regionprops_table can be used instead
fiber_ca_predicted_cat_1 = imread("Exp_1_P.tif")
fiber_ca_cat_1 = imread("Exp_1_GT.tif")

fiber_ca_predicted_cat_1_masked = np.zeros(shape=(fiber_ca_predicted_cat_1.shape), dtype=np.uint8)
fiber_ca_predicted_cat_1_gauss = np.zeros(shape=(fiber_ca_predicted_cat_1.shape), dtype=np.float32)


for i in range(fiber_ca_predicted_cat_1.shape[0]):
    fiber_ca_predicted_cat_1_gauss[i] = filters.gaussian(fiber_ca_predicted_cat_1[i], sigma=2)

for i in range(fiber_ca_predicted_cat_1.shape[0]):
    fiber_ca_predicted_cat_1_masked[i] = fiber_ca_predicted_cat_1_gauss[i]>0.00001*np.std(fiber_ca_predicted_cat_1_gauss[i])

fiber_ca_predicted_mask_clean_exp_1 = np.zeros(shape=(fiber_ca_predicted_cat_1.shape), dtype=np.uint8)
fiber_ca_predicted_mask_label_cat_1 = np.zeros(shape=(fiber_ca_predicted_cat_1.shape), dtype=np.uint32)


for i in range(fiber_ca_predicted_cat_1.shape[0]):
    fiber_ca_predicted_mask_label_cat_1[i],fiber_ca_predicted_mask_label_cat_1_number = ndi.label(fiber_ca_predicted_cat_1_masked[i])

fiber_ca_predicted_label_clean_fiber_mask_label_cat_1 = np.copy(fiber_ca_predicted_mask_label_cat_1)

for i in range(fiber_ca_predicted_cat_1.shape[0]):
    for fiber_ID in np.unique(fiber_ca_predicted_mask_label_cat_1[i]):
        fiber_mask_cat_1 = fiber_ca_predicted_mask_label_cat_1[i]==fiber_ID
        if np.sum(fiber_mask_cat_1) < 3000:
            fiber_ca_predicted_label_clean_fiber_mask_label_cat_1[i][fiber_mask_cat_1] = 0

fiber_ca_predicted_mask_clean_exp_1 = fiber_ca_predicted_label_clean_fiber_mask_label_cat_1 > 0

imsave("Exp_1_P_masked.tif", fiber_ca_predicted_mask_clean_exp_1.astype(np.uint8))



fiber_ca_exp_1_masked = np.zeros(shape=(fiber_ca_predicted_cat_1.shape), dtype=np.uint8)
fiber_ca_cat_1_gauss = np.zeros(shape=(fiber_ca_predicted_cat_1.shape), dtype=np.float32)


for i in range(fiber_ca_cat_1.shape[0]):
    fiber_ca_cat_1_gauss[i] = filters.gaussian(fiber_ca_cat_1[i], sigma=2)

for i in range(fiber_ca_predicted_cat_1.shape[0]):
    fiber_ca_exp_1_masked[i] = fiber_ca_cat_1_gauss[i]>0.00001*np.std(fiber_ca_cat_1_gauss[i])

imsave("Exp_1_GT_masked.tif", fiber_ca_exp_1_masked.astype(np.uint8))



fiber_ca_predicted_cat_2 = imread("Exp_2_P.tif")
fiber_ca_cat_2 = imread("Exp_2_GT.tif")

fiber_ca_predicted_cat_2_masked = np.zeros(shape=(fiber_ca_predicted_cat_2.shape), dtype=np.uint8)
fiber_ca_predicted_cat_2_gauss = np.zeros(shape=(fiber_ca_predicted_cat_2.shape), dtype=np.float32)


for i in range(fiber_ca_predicted_cat_2.shape[0]):
    fiber_ca_predicted_cat_2_gauss[i] = filters.gaussian(fiber_ca_predicted_cat_2[i], sigma=2)

for i in range(fiber_ca_predicted_cat_2.shape[0]):
    fiber_ca_predicted_cat_2_masked[i] = fiber_ca_predicted_cat_2_gauss[i]>0.00001*np.std(fiber_ca_predicted_cat_2_gauss[i])

fiber_ca_predicted_mask_clean_exp_2 = np.zeros(shape=(fiber_ca_predicted_cat_2.shape), dtype=np.uint8)
fiber_ca_predicted_mask_label_cat_2 = np.zeros(shape=(fiber_ca_predicted_cat_2.shape), dtype=np.uint32)


for i in range(fiber_ca_predicted_cat_2.shape[0]):
    fiber_ca_predicted_mask_label_cat_2[i],fiber_ca_predicted_mask_label_cat_2_number = ndi.label(fiber_ca_predicted_cat_2_masked[i])

fiber_ca_predicted_label_clean_fiber_mask_label_cat_2 = np.copy(fiber_ca_predicted_mask_label_cat_2)

for i in range(fiber_ca_predicted_cat_2.shape[0]):
    for fiber_ID in np.unique(fiber_ca_predicted_mask_label_cat_2[i]):
        fiber_mask_cat_2 = fiber_ca_predicted_mask_label_cat_2[i]==fiber_ID
        if np.sum(fiber_mask_cat_2) < 3000:
            fiber_ca_predicted_label_clean_fiber_mask_label_cat_2[i][fiber_mask_cat_2] = 0

fiber_ca_predicted_mask_clean_exp_2 = fiber_ca_predicted_label_clean_fiber_mask_label_cat_2 > 0

imsave("Exp_2_P_masked.tif", fiber_ca_predicted_mask_clean_exp_2.astype(np.uint8))



fiber_ca_exp_2_masked = np.zeros(shape=(fiber_ca_predicted_cat_2.shape), dtype=np.uint8)
fiber_ca_cat_2_gauss = np.zeros(shape=(fiber_ca_predicted_cat_2.shape), dtype=np.float32)


for i in range(fiber_ca_cat_2.shape[0]):
    fiber_ca_cat_2_gauss[i] = filters.gaussian(fiber_ca_cat_2[i], sigma=2)

for i in range(fiber_ca_predicted_cat_2.shape[0]):
    fiber_ca_exp_2_masked[i] = fiber_ca_cat_2_gauss[i]>0.00001*np.std(fiber_ca_cat_2_gauss[i])

imsave("Exp_2_GT_masked.tif", fiber_ca_exp_2_masked.astype(np.uint8))



fiber_ca_predicted_cat_3 = imread("Exp_3_P.tif")
fiber_ca_cat_3 = imread("Exp_3_GT.tif")

fiber_ca_predicted_cat_3_masked = np.zeros(shape=(fiber_ca_predicted_cat_3.shape), dtype=np.uint8)
fiber_ca_predicted_cat_3_gauss = np.zeros(shape=(fiber_ca_predicted_cat_3.shape), dtype=np.float32)


for i in range(fiber_ca_predicted_cat_3.shape[0]):
    fiber_ca_predicted_cat_3_gauss[i] = filters.gaussian(fiber_ca_predicted_cat_3[i], sigma=2)

for i in range(fiber_ca_predicted_cat_3.shape[0]):
    fiber_ca_predicted_cat_3_masked[i] = fiber_ca_predicted_cat_3_gauss[i]>0.00001*np.std(fiber_ca_predicted_cat_3_gauss[i])

fiber_ca_predicted_mask_clean_exp_3 = np.zeros(shape=(fiber_ca_predicted_cat_3.shape), dtype=np.uint8)
fiber_ca_predicted_mask_label_cat_3 = np.zeros(shape=(fiber_ca_predicted_cat_3.shape), dtype=np.uint32)


for i in range(fiber_ca_predicted_cat_3.shape[0]):
    fiber_ca_predicted_mask_label_cat_3[i],fiber_ca_predicted_mask_label_cat_3_number = ndi.label(fiber_ca_predicted_cat_3_masked[i])

fiber_ca_predicted_label_clean_fiber_mask_label_cat_3 = np.copy(fiber_ca_predicted_mask_label_cat_3)

for i in range(fiber_ca_predicted_cat_3.shape[0]):
    for fiber_ID in np.unique(fiber_ca_predicted_mask_label_cat_3[i]):
        fiber_mask_cat_3 = fiber_ca_predicted_mask_label_cat_3[i]==fiber_ID
        if np.sum(fiber_mask_cat_3) < 3000:
            fiber_ca_predicted_label_clean_fiber_mask_label_cat_3[i][fiber_mask_cat_3] = 0

fiber_ca_predicted_mask_clean_exp_3 = fiber_ca_predicted_label_clean_fiber_mask_label_cat_3 > 0

imsave("Exp_3_P_masked.tif", fiber_ca_predicted_mask_clean_exp_3.astype(np.uint8))



fiber_ca_exp_3_masked = np.zeros(shape=(fiber_ca_predicted_cat_3.shape), dtype=np.uint8)
fiber_ca_cat_3_gauss = np.zeros(shape=(fiber_ca_predicted_cat_3.shape), dtype=np.float32)


for i in range(fiber_ca_cat_3.shape[0]):
    fiber_ca_cat_3_gauss[i] = filters.gaussian(fiber_ca_cat_3[i], sigma=2)

for i in range(fiber_ca_predicted_cat_3.shape[0]):
    fiber_ca_exp_3_masked[i] = fiber_ca_cat_3_gauss[i]>0.00001*np.std(fiber_ca_cat_3_gauss[i])

imsave("Exp_3_GT_masked.tif", fiber_ca_exp_3_masked.astype(np.uint8))



fiber_ca_train = imread("train_GT.tif")

fiber_ca_train_masked = np.zeros(shape=(fiber_ca_train.shape), dtype=np.uint8)
fiber_ca_train_gauss = np.zeros(shape=(fiber_ca_train.shape), dtype=np.float32)


for i in range(fiber_ca_train.shape[0]):
    fiber_ca_train_gauss[i] = filters.gaussian(fiber_ca_train[i], sigma=2)

for i in range(fiber_ca_train.shape[0]):
    fiber_ca_train_masked[i] = fiber_ca_train_gauss[i]>0.00001*np.std(fiber_ca_train_gauss[i])

imsave("train_GT_masked.tif", fiber_ca_train_masked.astype(np.uint8))


results = {"fiber_area_cat_1"      : [],
           "fiber_area_cat_1_predicted"     : [],
           "fiber_mean_int_cat_1"     : [],
           "fiber_mean_int_cat_1_predicted"     : [],
           "fiber_area_cat_2"    : [],
           "fiber_area_cat_2_predicted"    : [],
           "fiber_mean_int_cat_2" : [],
           "fiber_mean_int_cat_2_predicted"     : [],
           "fiber_area_cat_3" : [],
           "fiber_area_cat_3_predicted"    : [],
           "fiber_mean_int_cat_3" : [],
           "fiber_mean_int_cat_3_predicted"     : [],
           "fiber_mean_int_train" : []}

for i in range(fiber_ca_predicted_cat_1.shape[0]):   
    results["fiber_area_cat_1"].append(0.2325149*0.2325149*np.sum(fiber_ca_exp_1_masked[i]))
    results["fiber_area_cat_1_predicted"].append(0.2325149*0.2325149*np.sum(fiber_ca_predicted_mask_clean_exp_1[i]))
    
    
for i in range(fiber_ca_predicted_cat_2.shape[0]):    
    results["fiber_area_cat_2"].append(0.2325149*0.2325149*np.sum(fiber_ca_exp_2_masked[i]))
    results["fiber_area_cat_2_predicted"].append(0.2325149*0.2325149*np.sum(fiber_ca_predicted_mask_clean_exp_2[i]))
    
    
for i in range(fiber_ca_predicted_cat_3.shape[0]):    
    results["fiber_area_cat_3"].append(0.2325149*0.2325149*np.sum(fiber_ca_exp_3_masked[i]))
    results["fiber_area_cat_3_predicted"].append(0.2325149*0.2325149*np.sum(fiber_ca_predicted_mask_clean_exp_3[i]))

for i in range(fiber_ca_predicted_cat_1.shape[0]):
    for fiber_id in np.unique(fiber_ca_exp_1_masked[i])[1:]:
        fiber_mask_cat_1_int = fiber_ca_exp_1_masked[i]==fiber_id
        results["fiber_mean_int_cat_1"].append(np.mean(fiber_ca_cat_1[i][fiber_mask_cat_1_int]))
        

for i in range(fiber_ca_predicted_cat_1.shape[0]):
    for fiber_id in np.unique(fiber_ca_predicted_mask_clean_exp_1[i])[1:]:
        fiber_mask_cat_1_predicted_int = fiber_ca_predicted_mask_clean_exp_1[i]==fiber_id
        results["fiber_mean_int_cat_1_predicted"].append(np.mean(fiber_ca_predicted_cat_1[i][fiber_mask_cat_1_predicted_int]))


for i in range(fiber_ca_predicted_cat_2.shape[0]):
    for fiber_id in np.unique(fiber_ca_exp_2_masked[i])[1:]:
        fiber_mask_cat_2_int = fiber_ca_exp_2_masked[i]==fiber_id
        results["fiber_mean_int_cat_2"].append(np.mean(fiber_ca_cat_2[i][fiber_mask_cat_2_int]))
        

for i in range(fiber_ca_predicted_cat_2.shape[0]):
    for fiber_id in np.unique(fiber_ca_predicted_mask_clean_exp_2[i])[1:]:
        fiber_mask_cat_2_predicted_int = fiber_ca_predicted_mask_clean_exp_2[i]==fiber_id
        results["fiber_mean_int_cat_2_predicted"].append(np.mean(fiber_ca_predicted_cat_2[i][fiber_mask_cat_2_predicted_int]))


for i in range(fiber_ca_predicted_cat_3.shape[0]):
    for fiber_id in np.unique(fiber_ca_exp_3_masked[i])[1:]:
        fiber_mask_cat_3_int = fiber_ca_exp_3_masked[i]==fiber_id
        results["fiber_mean_int_cat_3"].append(np.mean(fiber_ca_cat_3[i][fiber_mask_cat_3_int]))
        

for i in range(fiber_ca_predicted_cat_3.shape[0]):
    for fiber_id in np.unique(fiber_ca_predicted_mask_clean_exp_3[i])[1:]:
        fiber_mask_cat_3_predicted_int = fiber_ca_predicted_mask_clean_exp_3[i]==fiber_id
        results["fiber_mean_int_cat_3_predicted"].append(np.mean(fiber_ca_predicted_cat_3[i][fiber_mask_cat_3_predicted_int]))

for i in range(fiber_ca_train.shape[0]):
    for fiber_id in np.unique(fiber_ca_train_masked[i])[1:]:
        fiber_mask_train_masked_int = fiber_ca_train_masked[i]==fiber_id
        results["fiber_mean_int_train"].append(np.mean(fiber_ca_train[i][fiber_mask_train_masked_int]))

plt.figure(figsize=(450,150))

data = [results["fiber_mean_int_cat_1"], results["fiber_mean_int_cat_1_predicted"], results["fiber_mean_int_cat_2"], results["fiber_mean_int_cat_2_predicted"], results["fiber_mean_int_cat_3"], results["fiber_mean_int_cat_1_predicted"], results["fiber_mean_int_train"]]
fig, ax = plt.subplots()
ax.set_ylabel('Fiber mean intesity')
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
plt.savefig('Fiber mean intesity violin plot box plot.pdf', format='pdf')
plt.show()

plt.figure(figsize=(450,150))

data = [results["fiber_area_cat_1"], results["fiber_area_cat_1_predicted"], results["fiber_area_cat_2"], results["fiber_area_cat_2_predicted"], results["fiber_area_cat_3"], results["fiber_area_cat_3_predicted"]]
fig, ax = plt.subplots()
ax.set_ylabel('Fiber area in \u03BC$m^2$')
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

plt.xticks([1,2,3,4,5,6], ["Exp 1 GT", "Exp 1 P", "Exp 2 GT", "Exp 2 P", "Exp 3 GT", "Exp 3 P"])
plt.savefig('Fiber area violin plot box plot.pdf', format='pdf')
plt.show()

############################################
data_1 = fiber_ca_exp_1_masked + fiber_ca_predicted_mask_clean_exp_1

data_gleich_2_1 = data_1 > 1

data_gleich_1_1 = data_1 == 1
############################################
data_2 = fiber_ca_exp_2_masked + fiber_ca_predicted_mask_clean_exp_2

data_gleich_2_2 = data_2 > 1

data_gleich_1_2 = data_2 == 1
############################################
data_3 = fiber_ca_exp_3_masked + fiber_ca_predicted_mask_clean_exp_3

data_gleich_2_3 = data_3 > 1

data_gleich_1_3 = data_3 == 1
############################################
l = [2*np.sum(data_gleich_2_1)/(2*np.sum(data_gleich_2_1)+np.sum(data_gleich_1_1)), 2*np.sum(data_gleich_2_2)/(2*np.sum(data_gleich_2_2)+np.sum(data_gleich_1_2)), 2*np.sum(data_gleich_2_3)/(2*np.sum(data_gleich_2_3)+np.sum(data_gleich_1_3))]

plt.figure(figsize=(450,150))

names = ["Exp 1", "Exp 2", "Exp 3"]
data = l
fig, ax = plt.subplots()
plt.tick_params(axis="both", labelsize=15)
ax.set_ylabel('F1 score',fontsize=15)
bp = ax.bar(names, data, alpha=0.7, width=0.5, color="black")
plt.savefig('Fiber area F1.pdf', format='pdf')
plt.show()

plt.figure(figsize=(450,150))