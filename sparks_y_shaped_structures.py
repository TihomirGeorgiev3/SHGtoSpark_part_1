import numpy as np
import scipy.ndimage as ndi
import matplotlib.pyplot as plt
from tifffile import imread
from tifffile import imsave

#input of the algorithm: Exp_1_P_bn_sparks.tif, Exp_2_P_bn_sparks.tif, Exp_3_P_bn_sparks.tif, Exp_1_GT_bn_sparks.tif, Exp_2_GT_bn_sparks.tif, Exp_3_GT_bn_sparks.tif,
#Exp_1_SHG_y_shaped_structures.tif, Exp_2_SHG_y_shaped_structures.tif, Exp_3_SHG_y_shaped_structures.tif
#output of the algorithm: Spark area over y-shaped structure violin.pdf, Spark frequency y-shaped.pdf
#the files must be placed in the same folder as the .py file but this can be modified
#the approach to read out data is similar to the one described at https://github.com/embl-bio-it/image-analysis-with-python/blob/master/session-3to5/image_analysis_tutorial_solutions.ipynb
#the module skimage.measure.regionprops_table can be used instead
fiber_cat_1_y_shaped = imread("Exp_1_SHG_y_shaped_structures.tif")
fiber_cat_2_y_shaped = imread("Exp_2_SHG_y_shaped_structures.tif")
fiber_cat_3_y_shaped = imread("Exp_3_SHG_y_shaped_structures.tif")

fiber_sparks_bn_predicted_cat_1 = imread("Exp_1_P_bn_sparks.tif")
fiber_sparks_bn_predicted_cat_2 = imread("Exp_2_P_bn_sparks.tif")
fiber_sparks_bn_predicted_cat_3 = imread("Exp_3_P_bn_sparks.tif")

stack_bn_cat_1 = imread("Exp_1_GT_bn_sparks.tif")
stack_bn_cat_2 = imread("Exp_2_GT_bn_sparks.tif")
stack_bn_cat_3 = imread("Exp_3_GT_bn_sparks.tif")



stack_bn_cat_1_y_shaped = fiber_cat_1_y_shaped + stack_bn_cat_1
stack_bn_cat_2_y_shaped = fiber_cat_2_y_shaped + stack_bn_cat_2
stack_bn_cat_3_y_shaped = fiber_cat_3_y_shaped + stack_bn_cat_3

fiber_sparks_bn_predicted_cat_1_y_shaped = fiber_cat_1_y_shaped + fiber_sparks_bn_predicted_cat_1
fiber_sparks_bn_predicted_cat_2_y_shaped = fiber_cat_2_y_shaped + fiber_sparks_bn_predicted_cat_2
fiber_sparks_bn_predicted_cat_3_y_shaped = fiber_cat_3_y_shaped + fiber_sparks_bn_predicted_cat_3



stack_bn_cat_1_label = np.zeros(shape=(stack_bn_cat_1.shape), dtype=np.uint32)
stack_bn_cat_2_label = np.zeros(shape=(stack_bn_cat_2.shape), dtype=np.uint32)
stack_bn_cat_3_label = np.zeros(shape=(stack_bn_cat_3.shape), dtype=np.uint32)


fiber_sparks_bn_predicted_cat_1_label = np.zeros(shape=(fiber_sparks_bn_predicted_cat_1.shape), dtype=np.uint32)
fiber_sparks_bn_predicted_cat_2_label = np.zeros(shape=(fiber_sparks_bn_predicted_cat_2.shape), dtype=np.uint32)
fiber_sparks_bn_predicted_cat_3_label = np.zeros(shape=(fiber_sparks_bn_predicted_cat_3.shape), dtype=np.uint32)


for i in range(stack_bn_cat_1.shape[0]):
    stack_bn_cat_1_label[i], stack_bn_cat_1_label_number = ndi.label(stack_bn_cat_1[i])
    
for i in range(stack_bn_cat_2.shape[0]):
    stack_bn_cat_2_label[i], stack_bn_cat_2_label_number = ndi.label(stack_bn_cat_2[i]) 
    
for i in range(stack_bn_cat_3.shape[0]):
    stack_bn_cat_3_label[i], stack_bn_cat_3_label_number = ndi.label(stack_bn_cat_3[i])    
      
    
    

for i in range(fiber_sparks_bn_predicted_cat_1.shape[0]):
    fiber_sparks_bn_predicted_cat_1_label[i], fiber_sparks_bn_predicted_cat_1_label_number = ndi.label(fiber_sparks_bn_predicted_cat_1[i])

for i in range(fiber_sparks_bn_predicted_cat_2.shape[0]):
    fiber_sparks_bn_predicted_cat_2_label[i], fiber_sparks_bn_predicted_cat_2_label_number = ndi.label(fiber_sparks_bn_predicted_cat_2[i])

for i in range(fiber_sparks_bn_predicted_cat_3.shape[0]):
    fiber_sparks_bn_predicted_cat_3_label[i], fiber_sparks_bn_predicted_cat_3_label_number = ndi.label(fiber_sparks_bn_predicted_cat_3[i])
  
    
results = {"spark_area_y_shaped_cat_1"      : [],
           "spark_area_y_shaped_cat_1_predicted"     : [],
           "spark_area_y_shaped_cat_2"    : [],
           "spark_area_y_shaped_cat_2_predicted"    : [],
           "spark_area_y_shaped_cat_3" : [],
           "spark_area_y_shaped_cat_3_predicted"    : []}    
    
l_frequency_cat_1 = []
l_frequency_cat_1_predicted = []
l_frequency_cat_2 = []
l_frequency_cat_2_predicted = []
l_frequency_cat_3 = []
l_frequency_cat_3_predicted = []



for i in range(stack_bn_cat_1.shape[0]):
    for spark_ID in np.unique(stack_bn_cat_1_label[i])[1:]:
        spark_mask_cat_1 = stack_bn_cat_1_label[i]==spark_ID
        results["spark_area_y_shaped_cat_1"].append((np.sum(stack_bn_cat_1_y_shaped[i][spark_mask_cat_1])-np.sum(spark_mask_cat_1))/np.sum(spark_mask_cat_1))
        if (np.sum(stack_bn_cat_1_y_shaped[i][spark_mask_cat_1])-np.sum(spark_mask_cat_1))/np.sum(spark_mask_cat_1)>0.3:
            l_frequency_cat_1.append(np.sum(spark_mask_cat_1))


for i in range(fiber_sparks_bn_predicted_cat_1.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_1_label[i])[1:]:
        spark_mask_cat_1_predicted = fiber_sparks_bn_predicted_cat_1_label[i]==spark_ID        
        results["spark_area_y_shaped_cat_1_predicted"].append((np.sum(fiber_sparks_bn_predicted_cat_1_y_shaped[i][spark_mask_cat_1_predicted])-np.sum(spark_mask_cat_1_predicted))/np.sum(spark_mask_cat_1_predicted))
        if (np.sum(fiber_sparks_bn_predicted_cat_1_y_shaped[i][spark_mask_cat_1_predicted])-np.sum(spark_mask_cat_1_predicted))/np.sum(spark_mask_cat_1_predicted)>0.3:
            l_frequency_cat_1_predicted.append(np.sum(spark_mask_cat_1_predicted))


for i in range(stack_bn_cat_2.shape[0]):
    for spark_ID in np.unique(stack_bn_cat_2_label[i])[1:]:
        spark_mask_cat_2 = stack_bn_cat_2_label[i]==spark_ID
        results["spark_area_y_shaped_cat_2"].append((np.sum(stack_bn_cat_2_y_shaped[i][spark_mask_cat_2])-np.sum(spark_mask_cat_2))/np.sum(spark_mask_cat_2))
        if (np.sum(stack_bn_cat_2_y_shaped[i][spark_mask_cat_2])-np.sum(spark_mask_cat_2))/np.sum(spark_mask_cat_2)>0.3:
            l_frequency_cat_2.append(np.sum(spark_mask_cat_2))


for i in range(fiber_sparks_bn_predicted_cat_2.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_2_label[i])[1:]:
        spark_mask_cat_2_predicted = fiber_sparks_bn_predicted_cat_2_label[i]==spark_ID      
        results["spark_area_y_shaped_cat_2_predicted"].append((np.sum(fiber_sparks_bn_predicted_cat_2_y_shaped[i][spark_mask_cat_2_predicted])-np.sum(spark_mask_cat_2_predicted))/np.sum(spark_mask_cat_2_predicted))
        if (np.sum(fiber_sparks_bn_predicted_cat_2_y_shaped[i][spark_mask_cat_2_predicted])-np.sum(spark_mask_cat_2_predicted))/np.sum(spark_mask_cat_2_predicted)>0.3:
            l_frequency_cat_2_predicted.append(np.sum(spark_mask_cat_2_predicted))


for i in range(stack_bn_cat_3.shape[0]):
    for spark_ID in np.unique(stack_bn_cat_3_label[i])[1:]:
        spark_mask_cat_3 = stack_bn_cat_3_label[i]==spark_ID
        results["spark_area_y_shaped_cat_3"].append((np.sum(stack_bn_cat_3_y_shaped[i][spark_mask_cat_3])-np.sum(spark_mask_cat_3))/np.sum(spark_mask_cat_3))
        if (np.sum(stack_bn_cat_3_y_shaped[i][spark_mask_cat_3])-np.sum(spark_mask_cat_3))/np.sum(spark_mask_cat_3)>0.3:
            l_frequency_cat_3.append(np.sum(spark_mask_cat_3))


for i in range(fiber_sparks_bn_predicted_cat_3.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_3_label[i])[1:]:
        spark_mask_cat_3_predicted = fiber_sparks_bn_predicted_cat_3_label[i]==spark_ID
        results["spark_area_y_shaped_cat_3_predicted"].append((np.sum(fiber_sparks_bn_predicted_cat_3_y_shaped[i][spark_mask_cat_3_predicted])-np.sum(spark_mask_cat_3_predicted))/np.sum(spark_mask_cat_3_predicted))
        if (np.sum(fiber_sparks_bn_predicted_cat_3_y_shaped[i][spark_mask_cat_3_predicted])-np.sum(spark_mask_cat_3_predicted))/np.sum(spark_mask_cat_3_predicted)>0.3:
            l_frequency_cat_3_predicted.append(np.sum(spark_mask_cat_3_predicted))

plt.figure(figsize=(450,150))

data = [results["spark_area_y_shaped_cat_1"], results["spark_area_y_shaped_cat_1_predicted"], results["spark_area_y_shaped_cat_2"], results["spark_area_y_shaped_cat_2_predicted"], results["spark_area_y_shaped_cat_3"], results["spark_area_y_shaped_cat_3_predicted"]]
fig, ax = plt.subplots()
ax.set_ylabel('Spark area over y-shaped structure')
bp = ax.boxplot(data, showmeans=True)
vp = ax.violinplot(data)

for body in bp['boxes']:
    body.set_alpha(0)
    
for body in bp['medians']:
    body.set_alpha(0)

for body in bp['whiskers']:
    body.set_alpha(0)

for body in bp['caps']:
    body.set_alpha(0) 
    
for body in bp['fliers']:
    body.set_alpha(0)

for body in bp['means']:
    body.set_alpha(0.7) 

for body in vp['bodies']:
    body.set_alpha(0.7)

plt.xticks([1,2,3,4,5,6], ["Exp 1 GT", "Exp 1 P", "Exp 2 GT", "Exp 2 P", "Exp 3 GT", "Exp 3 P"])
plt.savefig('Spark area over y-shaped structure violin.pdf', format='pdf')
plt.show()

l_spark_frequency_y = [10000*len(l_frequency_cat_1)/(np.sum(fiber_cat_1_y_shaped)*0.2325149*0.2325149),10000*len(l_frequency_cat_1_predicted)/(np.sum(fiber_cat_1_y_shaped)*0.2325149*0.2325149),10000*len(l_frequency_cat_2)/(np.sum(fiber_cat_2_y_shaped)*0.2325149*0.2325149),10000*len(l_frequency_cat_2_predicted)/(np.sum(fiber_cat_2_y_shaped)*0.2325149*0.2325149),10000*len(l_frequency_cat_3)/(np.sum(fiber_cat_3_y_shaped)*0.2325149*0.2325149),10000*len(l_frequency_cat_3_predicted)/(np.sum(fiber_cat_3_y_shaped)*0.2325149*0.2325149)]

plt.figure(figsize=(450,150))
names = ["Exp 1 GT", "Exp 1 P", "Exp 2 GT", "Exp 2 P", "Exp 3 GT", "Exp 3 P"]
data = l_spark_frequency_y
fig, ax = plt.subplots()
ax.set_ylabel('Spark frequency in sparks/(10000*\u03BC$m^2$*s)')
bp = ax.bar(names, data, alpha=0.7)
plt.savefig('Spark frequency y shaped um.pdf', format='pdf')
plt.show()