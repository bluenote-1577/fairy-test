import numpy as np
from sklearn.preprocessing import StandardScaler
import sys
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import svd
from numpy.linalg import matrix_rank
import random

def normalize_columns(df):
    scaler = StandardScaler()
    return pd.DataFrame(scaler.fit_transform(df), columns=df.columns)

def compute_svd_entropy(matrix):
    U, s, V = svd(matrix, full_matrices=False)
    norm_vals = s / np.sum(s)
    print(s, norm_vals)
    entropy = -np.sum(norm_vals * np.log(norm_vals))
    return entropy

# Read the matrix
file_path = sys.argv[1]  # Update this with your actual file path
df = pd.read_csv(file_path, delim_whitespace=True, header=0)

# For simplicity, I'm considering only the first 1000 rows
df = df.iloc[:1000, 1:]
df = normalize_columns(df)

# Columns to compute SVD entropy
columns = list(range(4, df.shape[1], 2))
print(columns)

# Compute SVD entropy for specific columns
svd_entropies = []
running_col = []
for col in columns:
    running_col.append(col)
    matrix = df.iloc[:, running_col].values
    entropy = compute_svd_entropy(matrix)
    svd_entropies.append(entropy)

plt.plot(range(0, len(columns)), svd_entropies, label='Sequential Columns')

# Random Sampling
#random_svd_entropies = []
#for i in range(1, len(columns)):
#    random_cols = sorted(random.sample(columns, i))
#    print(random_cols)
#    matrix = df.iloc[:, random_cols].values
#    entropy = compute_svd_entropy(matrix)
#    random_svd_entropies.append(entropy)
# Plotting
#plt.plot(range(1, len(columns)), random_svd_entropies, label='Random Sampling')

plt.xlabel('Number of Columns')
plt.ylabel('SVD Entropy')
plt.title('SVD Entropy as a Function of the Number of Columns')
plt.legend()
plt.show()

