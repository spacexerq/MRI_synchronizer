import scipy.io
import pypulseq


data_seq = scipy.io.loadmat('Sequence_from_seq_SE_rfdeath_5000.mat')
A_data = data_seq['A']
print(A_data)