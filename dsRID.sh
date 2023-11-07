# !bin/bash

python ./src/extract_train.py -b $1 -o $2/dsRID_train.tsv
python ./src/extract_null.py -b $1 -o $2/dsRID_null.tsv

python ./src/concat_train_null.py -t $2/dsRID_train.tsv \
 -n $2/dsRID_null.tsv -o $2/dsRID_data.tsv

python ./src/extract_whole.py -b $1 -o $2/dsRID_whole.tsv

python ./src/model_predict.py $2/dsRID_data.tsv $2
