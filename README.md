# Neural Networks Based Fractional Pixel Motion Estimation for HEVC

This repo presents the results obtained in our [latest paper](https://www.google.com "Place-holder"), along with 
steps and requirements on how to recreate the results.  
The repo is a fork of [HEVC reference software (HM-16.19)](https://hevc.hhi.fraunhofer.de), with additional tweaks 
on the Encoder source files to enable Fractional-pixel Motion Estimation (FME) using Artificial Neural Networks 
(ANNs). Mainly, most of our work is concerned with [TEncSearch.cpp](./source/Lib/TLibEncoder/TEncSearch.cpp). 
Additionally, we provide an accompanying Jupyter Notebook which goes through the training process of our ANN.

## Requirements  
The software packages we used to obtain our results had the following versions:  
1. GCC v7.4.1: Compiler
2. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) v3.3.7: Matrix Algebra
3. [FastAI](https://github.com/fastai/fastai) v0.7.2: Deep Learning, based on:
4. [PyTorch](https://pytorch.org) v0.3.1  

It's worth noting that we've used the O2 optimization flag in our [makefile](./build/linux/common/makefile.base)  

## Build and Run  
To build the binaries, simply go to build directory, and run the makefile "Note that we've commented out 
all modules other than the Encoder, since they're not relevant in our paper":
```
cd ./build/  
make -f makefile
```
To encode BlowingBubbles sequence with quantization parameter 22, you should first modify 
[configuration file](./cfg/per-sequence/BlowingBubbles.cfg) to point to the video sequence's 
location, the run:  
```
cd ./bin  
./TAppEncoderStatic -c ../cfg/encoder_lowdelay_P_main.cfg -c ../cfg/per-sequence/BlowingBubbles.cfg -q 22
```

## Directories
The directory structure remains the same as HM-16.19, with the addition of: 
* DL folder: Contains all Deep Learning related material, such as:
  * Saved FastAI models
  * Weights and Biases per Quantization Parameter
  * Helper Bash scripts to extract the data set, and to format the parameters
* Symlink to the FastAI directory. It can have any name, but I've used "fastai07" in the Jupyter Notebook. 
  This can be easily done by:  
```
ln -s $PATH_TO_FASTAI ./fastai07
```

## List of Modifications
Most of the codes are the same as generic HM-16.19. The list of modified files are as follows:  
1. [TEncSearch.cpp](./source/Lib/TLibEncoder/TEncSearch.cpp):
   Contains nearly all of our contributions, as well as codes for extracting the data set. You 
   can find all changes by searching for EMI (Ehab M. Ibrahim), and reading the comments above 
   each change.
2. [TEncSearch.h](./source/Lib/TLibEncoder/TEncSearch.h):
   Added a flag to xTZSearchHelp() to save the integer error values when set to True "search for EMI"
3. [makefile.base](./build/linux/common/makefile.base): 
   * Changed the executable used to gcc-7 "name may vary depending on OS"  
   * Used O2 flag instead of O3 "gave me better results"
4. Edited the files in ./cfg/per-sequence to point to the video location on my machine 
5. Added DL folder and NN_training Jupyter Notebook

## Dataset Extraction
To extract the Dataset:
1. Go to [TEncSearch.cpp](./source/Lib/TLibEncoder/TEncSearch.cpp), and uncomment the codes related to 
   dataset extraction
2. Build the binaries
3. Run the encoder for the desired video sequence. Alternatively, you can run the Extract_data.sh script 
   in ./DL/ directory, just modify the video's name and quantization parameters

## Profiling
Profiling was performed using Valgrind's callgrind tool. Just run the encoder binary using valgrind. No extra 
flags are needed. 
```
cd ./bin  
valgrind --tool=callgrind ./TAppEncoderStatic -c ../cfg/encoder_lowdelay_P_main.cfg -c ../cfg/per-sequence/BlowingBubbles.cfg -q 22
```
This will result in a new file "callgrind.out.XXX", where XXX is the process number. This file can be visualized 
using KCacheGrind:
```
kcachegrind callgrind.out.XXX
```
We're mainly concerned with the absolute number of clock cycles for our ANN "NN_pred()", the standard FME 
"xPatternSearchFracDIF()", and the whole encoder "total"  
