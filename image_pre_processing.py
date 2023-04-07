import numpy as np
import scipy.ndimage as ndi
from tifffile import imread
from tifffile import imsave
from skimage import filters

#input of the algorithm: mouse_1_cell_1_1m_Ca_test.tif
#output of the algorithm: fiber_sparks_binary_mouse_1_cell_1_1m_Ca_test.tif and ca_sparks_enhanced_mouse_1_cell_1_1m_Ca_test.tif
#the file must be placed in the same folder as the .py file but this can be modified
fiber_ca = imread("mouse_1_cell_1_1m_Ca_test.tif")

fiber_mask = np.zeros(shape=(512, 512), dtype=np.uint8)
fiber_gauss = np.zeros(shape=(fiber_ca.shape), dtype=np.uint32)
fiber_gauss_masked = np.zeros(shape=(fiber_ca.shape), dtype=np.uint32)
fiber_gauss_sparks_bn = np.zeros(shape=(fiber_ca.shape), dtype=np.uint8)
fiber_gauss_sparks_bn_label = np.zeros(shape=(fiber_ca.shape), dtype=np.uint32)
fiber_sparks_binary_mouse_1_cell_1_1m_Ca_test = np.zeros(shape=(fiber_ca.shape), dtype=np.uint8)

for i in range(fiber_ca.shape[0]):
    fiber_gauss[i] = 1000*filters.gaussian(fiber_ca[i], sigma=2)

fiber_mean = np.mean(fiber_gauss, axis=0)

fiber_mask = fiber_mean > 2*np.std(fiber_mean)

fiber_mask_label = np.zeros(shape=(512, 512), dtype=np.uint32)

fiber_mask_label, fiber_mask_label_number = ndi.label(fiber_mask)

clean_label_fiber_mask = np.copy(fiber_mask_label)
ca_sparks_enhanced_1 = np.zeros(shape=(fiber_ca.shape[0],512, 512), dtype=np.uint32)
ca_sparks_enhanced_mouse_1_cell_1_1m_Ca_test = np.zeros(shape=(fiber_ca.shape[0],512, 512), dtype=np.uint8)
fiber_mask_clean = np.zeros(shape=(512, 512), dtype=np.uint8)

for fiber_ID in np.unique(fiber_mask_label):
    fiber_mask = fiber_mask_label==fiber_ID
    if np.sum(fiber_mask) < 3000:
        clean_label_fiber_mask[fiber_mask] = 0
# the approach in lines 35-38 is similar to the one described at https://github.com/embl-bio-it/image-analysis-with-python/blob/master/session-3to5/image_analysis_tutorial_solutions.ipynb
# instead of lines 35-38 also the module skimage.morphology.remove_small_objects of scikit image can be used as described in lines 63-65
fiber_mask_clean = clean_label_fiber_mask > 0

for i in range(fiber_ca.shape[0]):
    fiber_gauss_masked[i] = fiber_gauss[i]*fiber_mask_clean

for fiber_id in np.unique(fiber_mask_clean)[1:]:
    fiber_mask = fiber_mask_clean==fiber_id
    for i in range(fiber_ca.shape[0]):
        fiber_gauss_sparks_bn[i] = fiber_gauss_masked[i] > np.mean(fiber_gauss_masked[i][fiber_mask]) + 3.6*np.std(fiber_gauss_masked[i][fiber_mask])
        
for i in range(fiber_ca.shape[0]):
    fiber_gauss_sparks_bn_label[i], sparks_number = ndi.label(fiber_gauss_sparks_bn[i])
    
fiber_gauss_sparks_bn_label_clean = np.copy(fiber_gauss_sparks_bn_label)

for i in range(fiber_ca.shape[0]):
    for spark_ID in np.unique(fiber_gauss_sparks_bn_label[i]):
        spark_mask = fiber_gauss_sparks_bn_label[i]==spark_ID
        if np.sum(spark_mask) < 55:
            fiber_gauss_sparks_bn_label_clean[i][spark_mask] = 0
# the approach in lines 56-60 is similar to the one described at https://github.com/embl-bio-it/image-analysis-with-python/blob/master/session-3to5/image_analysis_tutorial_solutions.ipynb
# instead of lines 56-60 also the following code can be used that is based on the module skimage.morphology.remove_small_objects of scikit image
#from skimage import morphology
#for i in range(fiber_ca.shape[0]):
#    fiber_gauss_sparks_bn_label_clean[i] = morphology.remove_small_objects(fiber_gauss_sparks_bn_label[i], 55)

fiber_sparks_binary_mouse_1_cell_1_1m_Ca_test = fiber_gauss_sparks_bn_label_clean > 0

for i in range(fiber_ca.shape[0]):
    ca_sparks_enhanced_1[i] = fiber_sparks_binary_mouse_1_cell_1_1m_Ca_test[i]*9*fiber_ca[i]*fiber_mask_clean
    
ca_sparks_enhanced_mouse_1_cell_1_1m_Ca_test = (ca_sparks_enhanced_1 + fiber_ca*fiber_mask_clean)/10

imsave("fiber_sparks_binary_mouse_1_cell_1_1m_Ca_test.tif", fiber_sparks_binary_mouse_1_cell_1_1m_Ca_test.astype(np.uint8))
imsave("ca_sparks_enhanced_mouse_1_cell_1_1m_Ca_test.tif", ca_sparks_enhanced_mouse_1_cell_1_1m_Ca_test.astype(np.uint8))