# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 11:50:27 2019

@author: nickk
"""

import numpy as np
import scipy.io as io
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.utils.multiclass import unique_labels
import copy

def plot_confusion_matrix(y_true, y_pred, classes, ax,
                          normalize=False,
                          title=None,
                          cmap=plt.cm.Blues,
                          write_xlabel=True,
                          write_ylabel=True):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    # classes = classes[unique_labels(y_true, y_pred)]
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        #print("Normalized confusion matrix")
    # else:
        #print('Confusion matrix, without normalization')

    
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap, vmin=0, vmax=1)
    # ax.figure.colorbar(im, ax=ax, ticks=[0, 0.25, 0.5, 0.75, 1.0])
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title)
    if write_ylabel:
        ax.set_ylabel('True label')
    if write_xlabel:
        ax.set_xlabel('Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center",
             rotation_mode="anchor")
    plt.setp(ax.get_yticklabels(), rotation=90, ha="center",
             rotation_mode="anchor")
    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    return ax


def get_train_data(t_web,features,
                   initial_data_all,
                   initial_para,
                   initial_carn):
    train_webs = [1, 2, 3, 4, 5, 6]
    test_webs = [train_webs.pop(t_web)]
    
    train_data = [initial_data_all[w][(initial_para[w] == 1) | (initial_carn[w] == 1), :]
                  for w in train_webs]
    test_data = [initial_data_all[w][(initial_para[w] == 1) | (initial_carn[w] == 1), :]
                 for w in test_webs]
    train_target = [initial_para[w][(initial_para[w] == 1) | (initial_carn[w] == 1)].reshape(-1,1)
                       for w in train_webs]
    test_target = [initial_para[w][(initial_para[w] == 1) | (initial_carn[w] == 1)].reshape(-1,1)
                     for w in test_webs]
    
    train_data = np.vstack(train_data)
    test_data = np.vstack(test_data)
    train_data[np.isnan(train_data)] = -1
    test_data[np.isnan(test_data)] = -1
    train_target = np.vstack(train_target)
    test_target = np.vstack(test_target)
    test_target = ['Para' if y == 1 else 'Carn' for y in test_target]
    train_target = ['Para' if y == 1 else 'Carn' for y in train_target] 
    return train_data, test_data, train_target, test_target
    
plt.style.use('classic')
plt.rc('font', family='serif')
S = io.loadmat('AgglomerationPropsMaxLinkage.mat')
props_local = S['propsLocal']
web_nos_col = props_local[:, 0]
web_sizes = props_local[:, 12]
v = [1, 2, 3, 5, 6, 7, 8, 9, 15, 16, 17, 18]
v = [1, 2, 3, 5, 6, 7, 8, 9, 15, 16, 17, 18]
props = ('vul', 
         'gen', 
         'TL',
         'FlowRank',
         'Btwn',
         'E.Btwn',
         'm_vul_res',
         'm_gen_con',
         'cr_cc',
         'rc_cc',
         'r_cc',
         'c_cc')
props_fancy_0 = ('$v$',
               '$g$',
               '$T$',
               '$\\lambda$',
               '$C_B$',
               '$C_{EB}$',
               '$\\bar{v}_{res}$',
               '$\\bar{g}_{con}$',
               '$\\gamma^{rc}$',
               '$\\gamma^{cr}$',
               '$\\gamma^{r}$',
               '$\\gamma^{c}$')
features_used_0 = [i for i in range(12)]
idx_prop_dict = {k:v for k,v in zip(features_used_0, props)}
prop_idx_dict = {k:v for k,v in zip(props, features_used_0)}
webs = [1, 2, 3, 4, 5, 6]
web_codes = ['BSQ', 'CSM', 'EPB', 'FF', 'OH', 'STB']
initial_sizes = [int(props_local[props_local[:, 0] == w, 12][0]) for w in webs]
plt.close('all')

initial_data_all = {w:props_local[(web_nos_col == w) & (web_sizes == s)][:, v] 
                        for w, s in zip(webs, initial_sizes)}
initial_para = {w:props_local[(web_nos_col == w) & (web_sizes == s)][:, 10]
                for w, s in zip(webs, initial_sizes)}
initial_carn = {w:props_local[(web_nos_col == w) & (web_sizes == s)][:, 11]
                for w, s in zip(webs, initial_sizes)}

# Figure 1: feature importances.
fig1, axes1 = plt.subplots(nrows=3, ncols=2,
                           sharey='all', #sharex='all',
                           figsize=[12,8],
                           )
fig1.subplots_adjust(bottom=0.1, 
                     top=.95,
                     hspace = .65,
                     wspace=.1)
    
# Figure 2: Confusion Matrices.
fig2, axes2 = plt.subplots(nrows=2, ncols=3,
                           # sharey='all', #sharex='all',
                           figsize=[12,8],
                           )
fig2.subplots_adjust(bottom=0.1, 
                     top=.9,
                     hspace = .35,
                     wspace=0,
                     left=.01)

# Figure 3: Results of removing important predictors
fig3 = plt.figure(figsize=[8,8])
(ax3, ax3b) = fig3.subplots(nrows=2, ncols=1, sharex=True)
ax3b.set_ylim(bottom=0,top=5.5)

# Figure 4: Looking at two best features
fig4, axes4 = plt.subplots(nrows=2, ncols=3,
                           figsize=(12,8),
                           )
fig4.subplots_adjust(bottom=0.1,
                     top=0.9,
                     wspace=0.3,
                     hspace=0.4,)
for t_web in range(6):
    features_used = features_used_0.copy()
    props_fancy = props_fancy_0
    row1 = int(np.mod(t_web,3))
    col1 = int(np.floor(t_web/3))
    
    train_data, test_data, train_target, test_target = get_train_data(t_web,features_used,
                                                                      initial_data_all,
                                                                      initial_para,
                                                                      initial_carn,
                                                                      )
                                                    
    rf = RandomForestClassifier(n_estimators=500,
                                oob_score=True,
                                max_depth=3,
                                min_samples_leaf=4  ,
                                n_jobs=-1,
                                random_state=42)
    rf.fit(X=train_data[:,features_used],
           y=np.ravel(train_target))
    prob_train = rf.predict_proba(X=train_data[:,features_used])
    y_hat = rf.predict(test_data[:,features_used])
    test_score = rf.score(test_data[:,features_used], np.ravel(test_target))
    print(rf.oob_score_,test_score)
    # Figure 1: Feature importances
    this_ax = axes1[row1, col1]
    this_ax.bar(x=range(1, len(features_used)+1),
                        height=rf.feature_importances_,
                        )
    this_ax.set_xticks(range(1,13))
    this_ax.set_xticklabels(props_fancy, rotation=0, ha='left')
    this_ax.set_title('Web Tested: ' + web_codes[t_web])
    this_ax.grid(which='major', axis='y')
    this_ax.set_axisbelow(True)
    write_xlabel = col1==0
    write_ylabel = row1==2
    if write_xlabel:
        this_ax.set_ylabel('Feature Importance')
    if write_ylabel:
        this_ax.set_xlabel('Node Property')
    
    # Figure 2: Confusion matrices
    this_ax = axes2[col1,row1]
    plot_confusion_matrix(test_target, y_hat,
                          ['Para','Carn'], ax=this_ax,
                          normalize=True,
                          title='Web Tested: ' + web_codes[t_web] + '\nAccuracy = ' + ('%.2f' % test_score),
                          write_ylabel=row1==0,
                          write_xlabel=col1==1)
    test_scores = [test_score]
    
    # Look at the random forests with the two best features from training set.
    fi = rf.feature_importances_.copy()
    top_feature_idxs = [fi.argmax()]
    fi[fi.argmax()] = 0
    top_feature_idxs.append(fi.argmax())
    rf2 = copy.copy(rf)
    rf2.fit(X=train_data[:,top_feature_idxs],
           y=np.ravel(train_target))
    y_hat = rf2.predict_proba(test_data[:,top_feature_idxs])
    test_score = rf2.score(test_data[:,top_feature_idxs], np.ravel(test_target))
    min_x, min_y = np.vstack((train_data[:,top_feature_idxs],
                             test_data[:,top_feature_idxs])
                             ).min(axis=0)
    max_x, max_y = np.vstack((train_data[:,top_feature_idxs],
                             test_data[:,top_feature_idxs])
                             ).max(axis=0)
    range_x = max_x - min_x
    range_y = max_y - min_y
    x = np.linspace(min_x-.015*np.abs(min_x),max_x,200)
    y = np.linspace(min_y,max_y,200)
    xx, yy = np.meshgrid(x, y)
    grid_train = np.vstack((xx.ravel(), yy.ravel())).T
    z = rf2.predict_proba(grid_train)
    
    this_ax = axes4[col1, row1]
    h = this_ax.contourf(xx, yy, z[:,1].reshape([200,200]),
                         cmap='coolwarm', levels=np.linspace(0,1,11),
                         vmin=0, vmax=1.,
                         )
    x_test = test_data[:,top_feature_idxs[0]]
    y_test = test_data[:,top_feature_idxs[1]]
    n_obs = x_test.shape[0]
    rx = np.sqrt(np.random.rand(n_obs))*.01*range_x
    ry = rx*range_y/range_x
    theta = np.random.rand(n_obs)*2*np.pi
    x_test += rx*np.cos(theta)   
    y_test += ry*np.sin(theta)
    c_test = [1 if z == 'Para' else 0 for z in test_target]
    this_ax.scatter(x=x_test, y=y_test, c=c_test, cmap='coolwarm', vmin=0, vmax=1.)
    this_ax.set_title('Test Web: ' + web_codes[t_web] + '\nAcc. = %.2f' % test_score)
    this_ax.set_xlabel(props_fancy[top_feature_idxs[0]])
    this_ax.set_ylabel(props_fancy[top_feature_idxs[1]])
    # Now, look at sequences of random forests with successive variables 
    # removed; look at how the accuracy suffers.
    
    # Most important variables removed:
    fi = rf.feature_importances_
    features_dropped = [x for x, y in zip(props_fancy, fi == fi.max()) if y]
    props_fancy = [x for x, y in zip(props_fancy, fi != fi.max()) if y]
    features_used = [x for x, y in zip(features_used, fi != fi.max()) if y]
    

    
    while features_used != []:
        
        train_data, test_data, train_target, test_target = get_train_data(t_web,features_used,
                                                                      initial_data_all,
                                                                      initial_para,
                                                                      initial_carn,
                                                                      )
        rf.fit(X=train_data[:,features_used],
               y=np.ravel(train_target))
        prob_train = rf.predict_proba(X=train_data[:,features_used])
        y_hat = rf.predict(test_data[:,features_used])
        test_score = rf.score(test_data[:,features_used], np.ravel(test_target))
        fi = rf.feature_importances_
        features_dropped.extend([x for x, y in zip(props_fancy, fi == fi.max()) if y])
        props_fancy = [x for x, y in zip(props_fancy, fi != fi.max()) if y]
        features_used = [x for x, y in zip(features_used, fi != fi.max()) if y]
        #print(feature_dropped)
        test_scores.append(test_score)
    ax3.plot(range(12),test_scores, label= web_codes[t_web])
    for i, f in enumerate(features_dropped):
        ax3b.text(i, 5-t_web, f,
                  horizontalalignment='left',
                  verticalalignment='bottom')
    ax3b.set_yticks(range(6))
web_codes.reverse()
ax3b.set_yticklabels(web_codes)
web_codes.reverse()
ax3.legend(title='Test Web:',ncol=3,loc=3, fontsize=10)
ax3.set_ylabel('Accuracy')
ax3b.set_xticks(range(12))
ax3b.set_xlabel('Number of Features Removed')
ax3.set_title('Accuracy as Most important features are removed')
ax3b.set_ylabel('Web Tested')
ax3b.set_title('Next Feature Dropped')

cbar4 = fig4.colorbar(h, ax=axes4.ravel().tolist(),norm=mpl.colors.Normalize())
cbar4.set_label('Probability of Parasite')
fig1.savefig('../figures/feature_importances.png',dpi=300)
fig2.savefig('../figures/normalized_cm.png',dpi=300)
fig3.savefig('../figures/feature_removal.png',dpi=300)
fig4.savefig('../figures/2-feature_data.png',dpi=300)
