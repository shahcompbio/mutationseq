import os
import sys
import gzip
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib


dirpath = os.path.dirname(os.path.realpath(__file__))

PAIRED_TRAININGDATA = os.path.join(dirpath, "model_paired_v4.1.2_features.txt.gz")
PAIRED_MODEL = os.path.join(dirpath, "model_v4.1.2.pickle")

PAIRED_DEEP_TRAININGDATA = os.path.join(dirpath, "model_paired_deep_v0.2.0_features.txt.gz")
PAIRED_DEEP_MODEL = os.path.join(dirpath, "model_deep_v0.2.0.pickle")

SINGLE_TRAININGDATA = os.path.join(dirpath, 'model_single_v4.0.2_features.txt.gz')
SINGLE_MODEL = os.path.join(dirpath, 'model_single_v4.0.2.pickle')




def load_training_data(features, labels):
    X = []
    Y = []
    
    for line in gzip.open(features):
        line = line.strip().split('\t')
        line_feat = line[-2]
        line_label = line[-1]
    
        line_label = 1 if line_label in labels else 0
    
    
        line_feat = eval(line_feat)
    
    
        X.append(line_feat)
        Y.append(line_label)
    
    
    X  = np.array(X)
    Y = np.array(Y)

    return X,Y

def train_and_dump_model(X, Y, model_file, type='paired', deep=False):
    
    model = RandomForestClassifier(random_state=0, n_estimators=3000,
                                                n_jobs=1)
    
    model.fit(X,Y)

    model.name = "TCGA Benchmark 4 feature set with coverage info"

    if type == 'paired': 
        if deep:
            model.version = 'deep_0.2'
        else:
            model.version = "4.1.2"

    elif type == 'single':
        if deep:
            model.version = ''
        else:
            model.version = 'single_4.0.2'
    
    joblib.dump(model, model_file, compress=9)


def setup_museq_models():
    #train paired model
    labels  = ['SOMATIC','HET_ONE','HOM_ONE']

    X,Y = load_training_data(PAIRED_TRAININGDATA, labels)

    train_and_dump_model(X, Y, PAIRED_MODEL)

    labels = ['SOMATIC','HET','HOM','GERMLINE']
    X,Y = load_training_data(SINGLE_TRAININGDATA, labels)

    train_and_dump_model(X, Y, SINGLE_MODEL, type='single')

    labels = ['SOMATIC']
    X,Y = load_training_data(PAIRED_DEEP_TRAININGDATA, labels)

    train_and_dump_model(X, Y, PAIRED_DEEP_MODEL, deep=True)



if __name__ == "__main__":
    setup_museq_models()
