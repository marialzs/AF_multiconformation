import pandas as pd
import matplotlib.pyplot as plt

# Read the file into a DataFrame
df = pd.read_csv('dpout_explicitp.txt', delim_whitespace=True)

# Extract the pock_vol column
pock_vol_values = df['pock_vol'].values

# Plot the distribution of pock_vol values
plt.figure(figsize=(8, 6))
plt.hist(pock_vol_values, bins=20, color='blue', edgecolor='black', alpha=0.5)
plt.xlabel('Pocket Volume', fontsize = 20)
plt.ylabel('Frequency', fontsize = 20)
plt.grid(True, linestyle='--', alpha=0.5)
plt.xlim(0, 2200)

plt.tick_params(axis='both', which='major', labelsize=16)
# Save the figure
plt.savefig('pock_vol_distribution.png')

# Show the plot
plt.show()

