# Required libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set working directory to the script's folder
os.chdir("/Users/ilaypaz/Downloads/Python")
print("Current working directory:", os.getcwd())

# 1. Load mtcars dataset from local file
try:
    mtcars = pd.read_csv("mtcars.csv")
    print("mtcars loaded successfully")
    print("mtcars columns:", mtcars.columns)
    mtcars = mtcars.rename(columns={"model": "Car_model"})
except FileNotFoundError:
    print("Error: mtcars.csv not found in /Users/ilaypaz/Downloads/Python/")
    exit(1)
except KeyError:
    print("Error: 'model' column not found in mtcars.csv")
    print("Available columns:", mtcars.columns)
    exit(1)

# 2. Check type and display data
print(type(mtcars))
print(mtcars)

# 3. Create car_data, moving Car_model to the first column
try:
    car_data = mtcars.copy()
    print("car_data created successfully")
    print(car_data)
except NameError:
    print("Error: mtcars not defined, check earlier steps")
    exit(1)

# 4. Convert 'am' column (0, 1) to 'Automatic', 'Manual' and move 'am' to 2nd column
car_data['am'] = car_data['am'].apply(lambda x: 'Automatic' if x == 0 else 'Manual')
# Reorder columns to make 'am' the second column
columns = ['Car_model', 'am'] + [col for col in car_data.columns if col not in ['Car_model', 'am']]
car_data = car_data[columns]
print(car_data)

# 5. Remove columns 'vs', 'gear', 'carb'
car_data = car_data.drop(columns=['vs', 'gear', 'carb'])
print(car_data)

# 6. Create a box plot with jitter using seaborn
car_data['cyl'] = car_data['cyl'].astype('category')
plt.figure(figsize=(10, 6))
sns.boxplot(x='cyl', y='mpg', data=car_data)
sns.stripplot(x='cyl', y='mpg', data=car_data, size=2, color='red', alpha=0.5)
plt.xlabel('Number of cylinders')
plt.ylabel('Miles per gallon')
plt.title('MPG by Number of Cylinders')
plot = plt.gcf()
plt.show()

# 7. Save the plot as a PDF
plot.savefig('plot.pdf', format='pdf', bbox_inches='tight', dpi=300)
plt.close()

# 8. Create a list of three dataframes and save as CSV
listt = [None] * 3
for i in range(3):
    df = pd.DataFrame({
        'v': range(1, 101),
        'l': range(101, 201),
        'c': list('abcdefghijklmnopqrstuvwxyz' * 4)[:100]
    })
    listt[i] = df
    df.to_csv(f'dataframe{i+1}.csv', index=True)
print(listt)

# 9. Function to import CSV files
