#!/usr/bin/env python
'''Takes tables as inputs and applies different
classification methods
'''
import sys
import numpy as np

def fetch_data(file_in):
    '''Table of class (1 column) features (multiple columns)
    First line is: ind, name of the predictors 
    '''
    
    first_line = open(file_in).readline().strip().split()
    dt={'names': [], 'formats': ['i4']}
    print len(first_line)-1, 'predictors'
    dt['names'] = tuple([name+'_%d' % i for i, name in enumerate(first_line)]) 
    for i in first_line[1:]: dt['formats'].append('f8')
    table = np.loadtxt(file_in, skiprows=1)
    print table.shape[0], 'objects'
    return table

def svm_classification(table):
    '''
    '''
    from scikits.learn import svm
    X = table[:, 1:]
    Y = table[:, 0]
    clf = svm.SVC()
    clf.fit(X, Y)
    print clf.support_



def svm_roc(table, kernel='linear', C=1.0):
    '''Classification and ROC analysis
    '''
    from scikits.learn import svm
    from scikits.learn.metrics import roc_curve, auc
    import pylab as pl
    
    X = table[:, 1:]
    y = table[:, 0]
    n_samples, n_features = X.shape
    p = range(n_samples)
    np.random.seed(0)
    np.random.shuffle(p)
    X, y = X[p], y[p]
    half = int(n_samples/2)
    
    # Run classifier
    classifier = svm.SVC(kernel=kernel, probability=True, C=C)
    probas_ = classifier.fit(X[:half],y[:half]).predict_proba(X[half:])
    
    # Compute ROC curve and area the curve
    fpr, tpr, thresholds = roc_curve(y[half:], probas_[:,1])
    roc_auc = auc(fpr, tpr)
    print "Area under the ROC curve : %f" % roc_auc
    
    # Plot ROC curve
    pl.figure(-1)
    pl.clf()
    pl.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
    pl.plot([0, 1], [0, 1], 'k--')
    pl.xlim([0.0,1.0])
    pl.ylim([0.0,1.0])
    pl.xlabel('False Positive Rate')
    pl.ylabel('True Positive Rate')
    pl.title('Receiver operating characteristic example')
    pl.legend(loc="lower right")
    pl.show()    



def lasso_paths(table):
    from scikits.learn.linear_model import LassoCV
    from scikits.learn.cross_val import KFold
    import pylab as pl
    
    X = table[:, 1:]
    Y = table[:, 0]
        
    eps = 1e-3 # the smaller it is the longer is the path
    print "Computing regularization path using the lasso..."
    n_samples = X.shape[0]
    cv = KFold(n_samples/2, 5)
    model = LassoCV(eps=eps, cv=cv).fit(X, Y)
    ##############################################################################
    # Display results
    m_log_alphas = np.log10(model.alphas)
    m_log_alpha = np.log10(model.alpha)
    
    ax = pl.gca()
    ax.set_color_cycle(2 * ['b', 'r', 'g', 'c', 'k'])
    
    pl.subplot(211)
    ymin, ymax = -0.01, 0.015# pl.ylim()
    pl.plot(m_log_alphas, model.coef_path_)
    pl.vlines(m_log_alpha, ymin, ymax, linestyle='dashed')
    
    pl.xticks(())
    pl.ylabel('weights')
    pl.title('Lasso paths')
    pl.axis('tight')
    
    pl.subplot(212)
    ymin, ymax = 0.1, 0.4
    pl.plot(m_log_alphas, model.mse_path_)
    pl.vlines([m_log_alpha], ymin, ymax, linestyle='dashed')
    
    pl.xlabel('log(alpha)')
    pl.ylabel('MSE')
    pl.title('Mean Square Errors on each CV fold')
    pl.axis('tight')
    pl.ylim(ymin, ymax)
    
    pl.show()


def lasso_classification(table, alpha=0.3):
    '''
    '''
    from scikits.learn.linear_model import Lasso
    X = table[:, 1:]
    Y = table[:, 0]
#    n_samples, n_features = 50, 200
#    X = np.random.randn(n_samples, n_features)
#    coef = 3*np.random.randn(n_features)
#    coef[10:] = 0 # sparsify coef
#    Y = np.dot(X, coef)
#
#    # add noise
#    Y += 0.01*np.random.normal((n_samples,))

    # Split data in train set and test set
    n_samples = X.shape[0]
    items = np.random.permutation(n_samples)
    training_items = items[:n_samples/2]
    testing_items = items[n_samples/2:]
    X_train, y_train = X[training_items], Y[training_items]
    X_test, y_test = X[testing_items], Y[testing_items]
    
    lasso = Lasso(alpha=alpha, fit_intercept=True)
    lasso_fit = lasso.fit(X_train, y_train)
    print lasso_fit.coef_
    
    y_pred_lasso = lasso_fit.predict(X_test)
    y_collapsed = np.zeros_like(y_pred_lasso)
    collapsed_1 = y_pred_lasso >= 0.5
    y_collapsed[collapsed_1] = 1
    test = y_collapsed == y_test
    return float(test.sum()) / test.shape[0]

def main():
    
    import os
    
    args = sys.argv
    try:
        file_in = args[1]
        #method = args[2]
    except:
        file_in = os.path.expanduser('~/Work/ChIPSeq/data/pts_nts1_w6_12_m60.txt')
        #method = 'lasso'
        
    table = fetch_data(file_in)
    np.random.shuffle(table)
    #lasso_paths(table)
    svm_roc(table, kernel='linear', C=2.0)
if __name__ == '__main__':
    main()