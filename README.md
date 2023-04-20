# dsRID: Double-Stranded RNA Identifier
dsRID (Double-Stranded RNA Identifier) is a machine learning-based method designed to predict dsRNA regions in silico by leveraging the power of long-read RNA sequencing (RNA-seq) and molecular traits of dsRNAs. dsRNAs are potent triggers of innate immune responses upon recognition by cytosolic dsRNA sensor proteins. Identification of endogenous dsRNAs is critical to better understand the dsRNAome and its relevance to innate immunity related to human diseases.

## Getting Started
To get started with dsRID, you will need the following prerequisites:

PacBio long-read RNA-seq data

Python version 3.x or higher

Required packages listed in the requirements.txt file

## Installation
To install dsRID, clone this repository using the following command:
```
git clone https://github.com/gxiaolab/dsRID.git
```

Then, navigate to the dsRID directory and install the required packages using pip:

```
cd dsRID
pip install -r requirements.txt
```

## Usage
To predict dsRNA regions using dsRID, follow these steps:

Prepare the PacBio long-read RNA-seq data in the BAM format.
Run the dsRID script to include the training process with the following command:

```
bash dsRID.sh <input_file> <output_dir>
```

Here, <input_file> is the path to the input BAM file, and <output_dir> is the path to the output folder containing the training dsRNA regions and whole genome scan of predicted dsRNA regions.

Output will be the list of following files under /output_dir/
```
dsRID_train.tsv
dsRID_null.tsv
dsRID_whole.tsv
dsRID_data.tsv
feat_importance.tsv
cv_scores.tsv
```
dsRID_whole will contain column of predicted probability to be a candidate dsRNA region from the random forest model trained on dsRID_data.tsv. dsRID_data.tsv is a concatenated file of dsRID_train and dsRID_null each represents positive and negative dsRNA regions. After running model_predict.py, you should see feat_importance and cv_scores, each represents feature importance results and cross validation accuracy results. 

If you would like to opt out training process and directly proceed to prediction step, run

```
python 
```

## Results
We evaluated the performance of dsRID using PacBio long-read RNA-seq data derived from Alzheimer's disease (AD) brain and showed that our approach is highly accurate in predicting dsRNA regions in multiple datasets. We applied dsRID to an AD cohort sequenced by the ENCODE consortium and characterized the global dsRNA profile with potentially distinct expression patterns between AD and controls.

## Authors
Ryo Yamamoto - Bioinformatics Interdepartmental Program, University of California, Los Angeles, California, USA

Zhiheng Liu - Department of Integrative Biology and Physiology, University of California, Los Angeles, California, USA

Mudra Choudhury - Bioinformatics Interdepartmental Program, University of California, Los Angeles, California, USA

Xinshu Xiao - Bioinformatics Interdepartmental Program, Department of Integrative Biology and Physiology, Molecular Biology Institute, University of California, Los Angeles, California, USA

## License
This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.