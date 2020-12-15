import sys
localpath='/cluster/projects/p33/users/yunhanc/.local/lib/python3.8/site-packages'
sys.path.insert(1,localpath)
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statistics
import seaborn as sns

from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import Perceptron
from sklearn.linear_model import SGDClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.metrics import confusion_matrix, classification_report, accuracy_score

fname = sys.argv[1]
fpref = os.path.splitext(fname)[0]
matrix = pd.read_csv(fname, delim_whitespace=True)

#matrix = matrix.iloc[:,list(range(0,23))+[43]]
pheno = matrix['pheno']
matrix = matrix.drop(['IID','pheno'], axis=1)

#------------------------- plot -----------------------#
sns.countplot(pheno)
plt.savefig(fpref + "_hist.png")
#------------------------------------------------------#

X_train, X_val, y_train, y_val = train_test_split(matrix, pheno, train_size = 0.9, random_state=22)

text_clf = Pipeline([
        ('clf', LogisticRegression(n_jobs=1, C=1e5))
        #('clf', RandomForestClassifier(n_estimators=50, random_state=22))
        #('clf', GaussianNB())
        #('clf', SGDClassifier(loss='hinge', penalty='l2', alpha=1e-3, random_state=22, max_iter=5, tol=None))
        #('clf', Perceptron())
        #('clf', DecisionTreeClassifier())
])

text_clf.fit(X_train,y_train)

predicted = text_clf.predict(X_val)

# prediction accuracy
print(accuracy_score(y_val, predicted))
#print(np.mean(y_val == predicted))
# mean absolute error
print(statistics.mean(abs(y_val-predicted)))
#print(np.mean(abs(y_val-predicted)))
print(confusion_matrix(y_val,predicted))
print(classification_report(y_val,predicted))

df = pd.DataFrame()
df['pheno']= y_val
df['pred_pheno']=predicted
df.to_csv(fpref + '_pred.csv',index = None, sep=' ', header=True)
