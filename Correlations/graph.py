import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
csv_path = 'correlation_functions.csv'
data = pd.read_csv(csv_path)

# Clean the data: replace negative values with NaN and ensure all data is numeric
cleaned_data = data.replace(to_replace=[None], value=0)
for col in cleaned_data.columns[1:]:  # Skip the first column (Radius/ExpansionFactor)
    cleaned_data[col] = pd.to_numeric(cleaned_data[col], errors='coerce')  # Ensure all data is numeric
    cleaned_data[col] = cleaned_data[col].where(cleaned_data[col] > 0)  # Remove negative and non-numeric values

# Plot setup
plt.figure(figsize=(10, 6))
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'darkviolet', 'lime', 'teal', 'pink', 'grey', 'brown', 'navy']

# Extract radii from the first column and convert to integer for plotting
radii = cleaned_data['ExpansionFactor'].str.replace('Radius', '').astype(int)

# Plot each expansion factor with a unique color
for (column, color) in zip(cleaned_data.columns[1:], colors * (len(cleaned_data.columns) // len(colors) + 1)):
    plt.plot(radii, cleaned_data[column], label=column, color=color)

plt.xlabel('Radius')
plt.ylabel('Correlation Value')
plt.title('Correlation Functions by Expansion Factor')
plt.legend(title='Expansion Factor', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()

# Save the plot
plot_path = 'correlation_plot_distinct_colors.png'
plt.savefig(plot_path)
plt.show()

print(f'Plot saved to {plot_path}')
