# SHGtoSpark part 1
This is the first of two repositories that contains data and code for the paper “Prediction of fast subcellular Ca2+ signals based on subcellular structures in label-free second harmonic generation microscopy measurements”

Data set:

The raw data set consists of image pairs of second harmonic generation (SHG) and two photon microscopy measurements. They were part of the following publication:

Georgiev, T., Zapiec, B., Forderer, M., Fink, R. H. A. & Vogel, M. Colocalization properties of elementary Ca(2+) release signals with structures specific to the contractile filaments and the tubular system of intact mouse skeletal muscle fibers. J Struct Biol 192, 366-375 (2015). https://doi.org:10.1016/j.jsb.2015.09.018

I provide here small data sets. The algorithms can be tested on these data sets.

The following files are included in this data set:

mouse_1_cell_1_1_m_Ca_test.tif: this file is a raw Ca2+ measurement. The image pre-processing algorithm can be tested on this file

ca_sparks_enhanced_mouse_1_cell_1_1m_Ca_test.tif: enhanced Ca2+ sparks of the file mouse_1_cell_1_1_m_Ca_test.tif


10 files of Ca2+ signals and three files of SHG measurements are provided per experiment. For the first experiment the following 13 files are provided:

Exp_1_GT.tif: 10 images of Ca2+ signals with enhanced ground truth Ca2+ sparks.

Exp_1_GT_bn_sparks.tif: 10 binary images of ground truth Ca2+ sparks corresponding to Exp_1_GT.tif.

Exp_1_GT_edges.tif: 10 binary images of Ca2+ sparks edges corresponding to Exp_1_GT.tif.

Exp_1_GT_fiber_periphery.tif: 10 binary images of the fiber periphery corresponding to Exp_1_GT.tif.

Exp_1_GT_masked.tif: Fiber mask of the Ca2+ signal corresponding to Exp_1_GT.tif.

Exp_1_P.tif: 10 images of predicted Ca2+ sparks.

Exp_1_P_bn_sparks.tif: 10 binary images of predicted Ca2+ sparks corresponding to Exp_1_P.tif.

Exp_1_P_edges.tif: 10 binary images of Ca2+ sparks edges corresponding to Exp_1_P.tif.

Exp_1_P_fiber_periphery.tif: 10 binary images of the fiber periphery corresponding to Exp_1_P.tif.

Exp_1_P_masked.tif: Fiber mask of the Ca2+ signal corresponding to Exp_1_P.tif.

Exp_1_SHG.tif: Raw SHG signal of Exp 1.

Exp_1_SHG_mean.tif: Mean SHG signal of Exp 1.

Exp_1_SHG_y_shaped_structures.tif: Binary images of y-shaped structures.

Such files are provided for the other two experiments. For the training data set there are no predictions of Ca2+ signals.


Software:

The following software was used in the paper: Python, Jupyter, Spyder, Tensorflow, Tensorflow_io, Keras, Scikit learn, Scikit image, Numpy, Matplotlib, Scipy, Tifffile, Os, Pathlib, Time, Datetime, IPython, Xlsxwriter, ImageJ, Fiji, ChimeraX, https://www.tensorflow.org/tutorials/generative/pix2pix and https://github.com/embl-bio-it/image-analysis-with-python/blob/master/session-3to5/image_analysis_tutorial_solutions.ipynb.
Pearson correlation coefficients were calculated with Scipy, Matthews correlation coefficients were calculated with Scikit learn, structural similarity index (SSIM) and mean squared error (SME) were calculated with Scikit image. The algorithms for the calculation of Pearson correlation coefficients, Matthews correlation coefficients, structural similarity index (SSIM) and mean squared error (SME) will be provided upon request.
The following software version were used in this repository:
Python (version: 3.7.4), IPython (version: 7.8.0), Spyder (version: 3.3.6), Numpy (version: 1.21.6), Scikit image (version: 0.19.3), Scipy (version: 1.4.1), Matplotlib (version: 3.5.3), Tifffile (version: 2021.11.2)


Algorithms:

1.	Image pre-processing (image_pre_processing.py)
2.	Properties of generated fibers (properties_of_generated_fibers.py)
3.	Properties of generated Ca2+ sparks (properties_of_generated_ca_sparks.py)
4.	Spark area spark edge (spark_area_spark_edge.py)
5.	Sparks fiber periphery (sparks_fiber_periphery.py)
6.	Sparks myosin filaments (sparks_myosin_filaments.py)
7.	Spark y shaped structures (sparks_y_shaped_structures.py)

Hardware:

Intel(R) Core(TM) i7-9850H CPU 2.60GHz   2.59 GHz, RAM 64 GB, NVIDIA Quadro RTX 3000 6GB


Running time:

Running time for the small data sets provided here was less than 7.5 s for each algorithm.
