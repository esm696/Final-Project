#grabs all 3 files 
3dDeconvolve -input NAME.$subj.r*.scale+tlrc.HEAD                        \
    -censor motion_${subj}_censor.1D                                     \
#model 3  of the highest degree of polynomial
    -polort 3                                                            \
    -num_stimts 8   #number of timing files of stimuli of interest                                                    \
    -stim_times 1 stimuli_file.txt 'TENT(0,20,11)' #customizable HRF shape      \
    -stim_label 1 vis                                                    \
    -stim_times 2 stimuli_file.txt 'TENT(0,20,11)'                      \
    -stim_label 2 aud                                                    \
    -stim_file 3 motion_demean.1D'[0]' -stim_base 3 -stim_label 3 roll   \
    -stim_file 4 motion_demean.1D'[1]' -stim_base 4 -stim_label 4 pitch  \
    -stim_file 5 motion_demean.1D'[2]' -stim_base 5 -stim_label 5 yaw    \
    -stim_file 6 motion_demean.1D'[3]' -stim_base 6 -stim_label 6 dS     \
    -stim_file 7 motion_demean.1D'[4]' -stim_base 7 -stim_label 7 dL     \
    -stim_file 8 motion_demean.1D'[5]' -stim_base 8 -stim_label 8 dP     \
    -jobs 2                                                              \
    -gltsym 'SYM: vis -aud'  #basic contrast                                            \
    -glt_label 1 V-A   #statistical label                                                  \
    -gltsym 'SYM: 0.5*vis +0.5*aud'                                      \
    -glt_label 2 mean.VA                                                 \
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg    #F stat and t-stat maps        \
    -x1D_uncensored X.nocensor.xmat.1D      #visualize x matrix                             \
    -errts errts.${subj}   #errors and what not fit the model                                               \
    -bucket stats.$subj #all statistical files output
    
    
    
#GML


def GLM(X, y):  #define files

X = nd.read_csv(file.csv)
y = nd.read_csv(file.csv)

# checking matrix orientation
     if X.shape[1] > X.shape[0]:
         X = X.transpose()
    
# Calculate the dot product and invert the transposed design matrix 
tmp   = np.linalg.inv(X.transpose().dot(X))
#dot product again
tmp   = tmp.dot(X.transpose())
# fill variables
beta  = np.zeros((y.shape[0], X.shape[1]))
e = np.zeros(y.shape)
model = np.zeros(y.shape)
r = np.zeros(y.shape[0])
    
# beta values, the error and the correlation coefficients 
for i in range(y.shape[0]):
    beta[i]  = tmp.dot(y[i,:].transpose())
    model[i] = X.dot(beta[i])
    e[i]     = (y[i,:] - model[i])
    r[i]     = np.sqrt(model[i].var()/y[i,:].var())
         return beta, model, e, r
         
beta, model, e, r = GLM(design_matrix, data)


#automating the result 

for i in `cat subject_list.txt`; do
  tcsh proc_NAMEofSCRIPT.sh $i;
 mv ${i}.results $i;