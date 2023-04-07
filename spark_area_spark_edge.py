import numpy as np
import scipy.ndimage as ndi
import matplotlib.pyplot as plt
from tifffile import imread
from tifffile import imsave

#input of the algorithm: Exp_1_P_bn_sparks.tif, Exp_2_P_bn_sparks.tif, Exp_3_P_bn_sparks.tif, Exp_1_GT_bn_sparks.tif, Exp_2_GT_bn_sparks.tif, Exp_3_GT_bn_sparks.tif, train_GT_bn_sparks.tif
#output of the algorithm: Exp_1_GT_edges.tif, Exp_2_GT_edges.tif, Exp_3_GT_edges.tif, train_GT_edges.tif, Exp_1_P_edges.tif, Exp_2_P_edges.tif, Exp_3_P_edges.tif, Spark area spark edge boxplot.pdf
#the files must be placed in the same folder as the .py file but this can be modified
#the approach to read out data is similar to the one described at https://github.com/embl-bio-it/image-analysis-with-python/blob/master/session-3to5/image_analysis_tutorial_solutions.ipynb
#the module skimage.measure.regionprops_table can be used instead
fiber_sparks_bn_predicted_cat_1 = imread("Exp_1_P_bn_sparks.tif")
fiber_sparks_bn_predicted_cat_2 = imread("Exp_2_P_bn_sparks.tif")
fiber_sparks_bn_predicted_cat_3 = imread("Exp_3_P_bn_sparks.tif")

stack_bn_cat_1 = imread("Exp_1_GT_bn_sparks.tif")
stack_bn_cat_2 = imread("Exp_2_GT_bn_sparks.tif")
stack_bn_cat_3 = imread("Exp_3_GT_bn_sparks.tif")
stack_bn_train = imread("train_GT_bn_sparks.tif")

stack_bn_cat_1_label = np.zeros(shape=(stack_bn_cat_1.shape), dtype=np.uint32)
stack_bn_cat_2_label = np.zeros(shape=(stack_bn_cat_2.shape), dtype=np.uint32)
stack_bn_cat_3_label = np.zeros(shape=(stack_bn_cat_3.shape), dtype=np.uint32)
stack_bn_train_label = np.zeros(shape=(stack_bn_train.shape), dtype=np.uint32)

fiber_sparks_bn_predicted_cat_1_label = np.zeros(shape=(fiber_sparks_bn_predicted_cat_1.shape), dtype=np.uint32)
fiber_sparks_bn_predicted_cat_2_label = np.zeros(shape=(fiber_sparks_bn_predicted_cat_2.shape), dtype=np.uint32)
fiber_sparks_bn_predicted_cat_3_label = np.zeros(shape=(fiber_sparks_bn_predicted_cat_3.shape), dtype=np.uint32)


for i in range(stack_bn_cat_1.shape[0]):
    stack_bn_cat_1_label[i], stack_bn_cat_1_label_spark_number = ndi.label(stack_bn_cat_1[i])
    
for i in range(stack_bn_cat_2.shape[0]):
    stack_bn_cat_2_label[i], stack_bn_cat_2_label_spark_number = ndi.label(stack_bn_cat_2[i]) 
    
for i in range(stack_bn_cat_3.shape[0]):
    stack_bn_cat_3_label[i], stack_bn_cat_3_label_spark_number = ndi.label(stack_bn_cat_3[i])    
    
for i in range(stack_bn_train.shape[0]):
    stack_bn_train_label[i], stack_bn_train_label_spark_number = ndi.label(stack_bn_train[i])    
    
for i in range(fiber_sparks_bn_predicted_cat_1.shape[0]):
    fiber_sparks_bn_predicted_cat_1_label[i], fiber_sparks_bn_predicted_cat_1_label_spark_number = ndi.label(fiber_sparks_bn_predicted_cat_1[i])

for i in range(fiber_sparks_bn_predicted_cat_2.shape[0]):
    fiber_sparks_bn_predicted_cat_2_label[i], fiber_sparks_bn_predicted_cat_2_label_spark_number = ndi.label(fiber_sparks_bn_predicted_cat_2[i])

for i in range(fiber_sparks_bn_predicted_cat_3.shape[0]):
    fiber_sparks_bn_predicted_cat_3_label[i], fiber_sparks_bn_predicted_cat_3_label_spark_number = ndi.label(fiber_sparks_bn_predicted_cat_3[i])


  
stack_bn_cat_1_edges_1_1 = np.zeros(shape=(stack_bn_cat_1.shape), dtype=np.uint8)
stack_bn_cat_1_edges_1 = np.zeros(shape=(stack_bn_cat_1.shape), dtype=np.uint32)
stack_bn_cat_1_edges = np.zeros(shape=(stack_bn_cat_1.shape), dtype=np.uint8)

for i in range(stack_bn_cat_1.shape[0]):
    stack_bn_cat_1_edges_1_1[i] = ndi.binary_erosion(stack_bn_cat_1[i], iterations=1)

stack_bn_cat_1_edges = np.logical_xor(stack_bn_cat_1, stack_bn_cat_1_edges_1_1)

stack_bn_cat_1_edges_1 = stack_bn_cat_1_label*stack_bn_cat_1_edges

imsave("Exp_1_GT_edges.tif", stack_bn_cat_1_edges.astype(np.uint8))



stack_bn_cat_2_edges_1_1 = np.zeros(shape=(stack_bn_cat_2.shape), dtype=np.uint8)
stack_bn_cat_2_edges = np.zeros(shape=(stack_bn_cat_2.shape), dtype=np.uint8)

for i in range(stack_bn_cat_2.shape[0]):
    stack_bn_cat_2_edges_1_1[i] = ndi.binary_erosion(stack_bn_cat_2[i], iterations=1)

stack_bn_cat_2_edges = np.logical_xor(stack_bn_cat_2, stack_bn_cat_2_edges_1_1)

imsave("Exp_2_GT_edges.tif", stack_bn_cat_2_edges.astype(np.uint8))



stack_bn_cat_3_edges_1_1 = np.zeros(shape=(stack_bn_cat_3.shape), dtype=np.uint8)
stack_bn_cat_3_edges = np.zeros(shape=(stack_bn_cat_3.shape), dtype=np.uint8)

for i in range(stack_bn_cat_3.shape[0]):
    stack_bn_cat_3_edges_1_1[i] = ndi.binary_erosion(stack_bn_cat_3[i], iterations=1)

stack_bn_cat_3_edges = np.logical_xor(stack_bn_cat_3, stack_bn_cat_3_edges_1_1)

imsave("Exp_3_GT_edges.tif", stack_bn_cat_3_edges.astype(np.uint8))



stack_bn_train_edges_1_1 = np.zeros(shape=(stack_bn_train.shape), dtype=np.uint8)
stack_bn_train_edges = np.zeros(shape=(stack_bn_train.shape), dtype=np.uint8)

for i in range(stack_bn_train.shape[0]):
    stack_bn_train_edges_1_1[i] = ndi.binary_erosion(stack_bn_train[i], iterations=1)

stack_bn_train_edges = np.logical_xor(stack_bn_train, stack_bn_train_edges_1_1)

imsave("train_GT_edges.tif", stack_bn_train_edges.astype(np.uint8))



fiber_sparks_bn_predicted_cat_1_edges_1_1 = np.zeros(shape=(fiber_sparks_bn_predicted_cat_1.shape), dtype=np.uint8)
fiber_sparks_bn_predicted_cat_1_edges = np.zeros(shape=(fiber_sparks_bn_predicted_cat_1.shape), dtype=np.uint8)

for i in range(fiber_sparks_bn_predicted_cat_1.shape[0]):
    fiber_sparks_bn_predicted_cat_1_edges_1_1[i] = ndi.binary_erosion(fiber_sparks_bn_predicted_cat_1[i], iterations=1)

fiber_sparks_bn_predicted_cat_1_edges = np.logical_xor(fiber_sparks_bn_predicted_cat_1, fiber_sparks_bn_predicted_cat_1_edges_1_1)

imsave("Exp_1_P_edges.tif", fiber_sparks_bn_predicted_cat_1_edges.astype(np.uint8))



fiber_sparks_bn_predicted_cat_2_edges_1_1 = np.zeros(shape=(fiber_sparks_bn_predicted_cat_2.shape), dtype=np.uint8)
fiber_sparks_bn_predicted_cat_2_edges = np.zeros(shape=(fiber_sparks_bn_predicted_cat_2.shape), dtype=np.uint8)

for i in range(fiber_sparks_bn_predicted_cat_2.shape[0]):
    fiber_sparks_bn_predicted_cat_2_edges_1_1[i] = ndi.binary_erosion(fiber_sparks_bn_predicted_cat_2[i], iterations=1)

fiber_sparks_bn_predicted_cat_2_edges = np.logical_xor(fiber_sparks_bn_predicted_cat_2, fiber_sparks_bn_predicted_cat_2_edges_1_1)

imsave("Exp_2_P_edges.tif", fiber_sparks_bn_predicted_cat_2_edges.astype(np.uint8))



fiber_sparks_bn_predicted_cat_3_edges_1_1 = np.zeros(shape=(fiber_sparks_bn_predicted_cat_3.shape), dtype=np.uint8)
fiber_sparks_bn_predicted_cat_3_edges = np.zeros(shape=(fiber_sparks_bn_predicted_cat_3.shape), dtype=np.uint8)

for i in range(fiber_sparks_bn_predicted_cat_3.shape[0]):
    fiber_sparks_bn_predicted_cat_3_edges_1_1[i] = ndi.binary_erosion(fiber_sparks_bn_predicted_cat_3[i], iterations=1)

fiber_sparks_bn_predicted_cat_3_edges = np.logical_xor(fiber_sparks_bn_predicted_cat_3, fiber_sparks_bn_predicted_cat_3_edges_1_1)

imsave("Exp_3_P_edges.tif", fiber_sparks_bn_predicted_cat_3_edges.astype(np.uint8))


   
results = {"spark_area_edge_cat_1"      : [],
           "spark_area_edge_cat_1_predicted"     : [],
           "spark_area_edge_cat_2"    : [],
           "spark_area_edge_cat_2_predicted"    : [],
           "spark_area_edge_cat_3" : [],
           "spark_area_edge_cat_3_predicted"    : [],
           "spark_area_edge_train" : []}    
    
for i in range(stack_bn_cat_1.shape[0]):
    for spark_ID in np.unique(stack_bn_cat_1_label[i])[1:]:
        spark_mask_cat_1 = stack_bn_cat_1_label[i]==spark_ID
        results["spark_area_edge_cat_1"].append(np.sum(spark_mask_cat_1)/np.sum(stack_bn_cat_1_edges[i][spark_mask_cat_1]))


for i in range(fiber_sparks_bn_predicted_cat_1.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_1_label[i])[1:]:
        spark_mask_cat_1_predicted = fiber_sparks_bn_predicted_cat_1_label[i]==spark_ID
        results["spark_area_edge_cat_1_predicted"].append(np.sum(spark_mask_cat_1_predicted)/np.sum(fiber_sparks_bn_predicted_cat_1_edges[i][spark_mask_cat_1_predicted]))


for i in range(stack_bn_cat_2.shape[0]):
    for spark_ID in np.unique(stack_bn_cat_2_label[i])[1:]:
        spark_mask_cat_2 = stack_bn_cat_2_label[i]==spark_ID
        results["spark_area_edge_cat_2"].append(np.sum(spark_mask_cat_2)/np.sum(stack_bn_cat_2_edges[i][spark_mask_cat_2]))


for i in range(fiber_sparks_bn_predicted_cat_2.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_2_label[i])[1:]:
        spark_mask_cat_2_predicted = fiber_sparks_bn_predicted_cat_2_label[i]==spark_ID
        results["spark_area_edge_cat_2_predicted"].append(np.sum(spark_mask_cat_2_predicted)/np.sum(fiber_sparks_bn_predicted_cat_2_edges[i][spark_mask_cat_2_predicted]))


for i in range(stack_bn_cat_3.shape[0]):
    for spark_ID in np.unique(stack_bn_cat_3_label[i])[1:]:
        spark_mask_cat_3 = stack_bn_cat_3_label[i]==spark_ID
        results["spark_area_edge_cat_3"].append(np.sum(spark_mask_cat_3)/np.sum(stack_bn_cat_3_edges[i][spark_mask_cat_3]))


for i in range(fiber_sparks_bn_predicted_cat_3.shape[0]):
    for spark_ID in np.unique(fiber_sparks_bn_predicted_cat_3_label[i])[1:]:
        spark_mask_cat_3_predicted = fiber_sparks_bn_predicted_cat_3_label[i]==spark_ID
        results["spark_area_edge_cat_3_predicted"].append(np.sum(spark_mask_cat_3_predicted)/np.sum(fiber_sparks_bn_predicted_cat_3_edges[i][spark_mask_cat_3_predicted]))


for i in range(stack_bn_train.shape[0]):
    for spark_ID in np.unique(stack_bn_train_label[i])[1:]:
        spark_mask_train = stack_bn_train_label[i]==spark_ID
        results["spark_area_edge_train"].append(np.sum(spark_mask_train)/np.sum(stack_bn_train_edges[i][spark_mask_train]))


plt.figure(figsize=(450,150))

data = [results["spark_area_edge_cat_1"], results["spark_area_edge_cat_1_predicted"], results["spark_area_edge_cat_2"], results["spark_area_edge_cat_2_predicted"], results["spark_area_edge_cat_3"], results["spark_area_edge_cat_3_predicted"], results["spark_area_edge_train"]]
fig, ax = plt.subplots()
ax.set_ylabel('Spark area/spark edge')
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
plt.savefig("Spark area spark edge boxplot.pdf", format='pdf')
plt.show()