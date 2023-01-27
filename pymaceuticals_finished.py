#!/usr/bin/env python
# coding: utf-8

# # Pymaceuticals Inc.
# ---
# 
# ### Analysis
# 
# - Add your analysis here.
#  

# In[1]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st

# Study data files
mouse_metadata_path = "./data/Mouse_metadata.csv"
study_results_path = "./data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset
science_study_data_complete_df = pd.merge(study_results, mouse_metadata, how="left", on="Mouse ID")


# Display the data table for preview
science_study_data_complete_df


# In[2]:


len(science_study_data_complete_df["Mouse ID"].unique())


# In[3]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicated_mouse_ids = science_study_data_complete_df[science_study_data_complete_df.duplicated(subset=["Mouse ID", "Timepoint"])]["Mouse ID"].unique()
duplicated_mouse_ids


# In[4]:


# Optional: Get all the data for the duplicate mouse ID. 
duplicated_mouse_dataset = science_study_data_complete_df[science_study_data_complete_df["Mouse ID"] == 'g989']
duplicated_mouse_dataset


# In[5]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_study_data_complete = science_study_data_complete_df[science_study_data_complete_df["Mouse ID"].isin(duplicated_mouse_ids) == False]
clean_study_data_complete


# In[6]:


# Checking the number of mice in the clean DataFrame.
len(clean_study_data_complete["Mouse ID"].unique())


# ## Summary Statistics

# In[7]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary DataFrame.
means = clean_study_data_complete.groupby("Drug Regimen").mean()["Tumor Volume (mm3)"]
median = clean_study_data_complete.groupby("Drug Regimen").median()["Tumor Volume (mm3)"]
variance = clean_study_data_complete.groupby("Drug Regimen").var()["Tumor Volume (mm3)"]
standard_deviation = clean_study_data_complete.groupby("Drug Regimen").std()["Tumor Volume (mm3)"]
sems = clean_study_data_complete.groupby("Drug Regimen").sem()["Tumor Volume (mm3)"]

summary_table = pd.DataFrame({
    "Mean Tumor Volume": means,
    "Median Tumor Volume": median,
    "Tumor Volume Variance": variance,
    "Tumor Volume Std Dev": standard_deviation,
    "Tumor Volume Std Errors": sems
})
summary_table


# In[8]:


# Generate a summary statistics table of mean, median, variance, standard deviation, 
# and SEM of the tumor volume for each regimen

# Using the aggregation method, produce the same summary statistics in a single line.
summary_table = clean_study_data_complete.groupby("Drug Regimen").agg({
    "Tumor Volume (mm3)":["mean", "median", "var", "std", "sem"]
})
summary_table


# ## Bar and Pie Charts

# In[9]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using Pandas.
counts = clean_study_data_complete["Drug Regimen"].value_counts()
counts.plot(kind="bar")
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("Number of Mice Tested")
plt.show()


# In[10]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using pyplot.
counts = clean_study_data_complete["Drug Regimen"].value_counts()
plt.bar(counts.index.values, counts.values)
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("Number of Mice Tested")
plt.show()


# In[11]:


# Generate a pie plot showing the distribution of female versus male mice using Pandas
counts = clean_study_data_complete.Sex.value_counts()
counts.plot(kind="pie", autopct="%1.1f%%")


# In[12]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
counts = clean_study_data_complete.Sex.value_counts()
plt.pie(counts.values, labels=counts.index.values, autopct="%1.1f%%")
plt.ylabel("Sex")
plt.show()


# ## Quartiles, Outliers and Boxplots

# In[13]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens: 

# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
max_tumor = clean_study_data_complete.groupby(["Mouse ID"])["Timepoint"].max()
max_tumor = max_tumor.reset_index()

# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
merged_data = max_tumor.merge(clean_study_data_complete, on=["Mouse ID", "Timepoint"], how="left")
merged_data


# In[18]:


# Put treatments into a list for for loop (and later for plot labels)
treatment_list = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]

# Create empty list to fill with tumor vol data (for plotting)
tumor_vol_list = []

# Calculate the IQR and quantitatively determine if there are any potential outliers. 
for drug in treatment_list:
    
    # Locate the rows which contain mice on each drug and get the tumor volumes
    final_tumor_vol = merged_data.loc[merged_data["Drug Regimen"] == drug, "Tumor Volume (mm3)"]
    
    # add subset 
    tumor_vol_list.append(final_tumor_vol)
    
    # Determine outliers using upper and lower bounds
    quartiles = final_tumor_vol.quantile([.25, .5, .75])
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr = upperq-lowerq
    lower_bound = lowerq - (1.5 * iqr)
    upper_bound = upperq + (1.5 * iqr)
    
    outliers = final_tumor_vol.loc[(final_tumor_vol < lower_bound) | (final_tumor_vol > upper_bound)]
    print(f"{drug}'s potential outliers {outliers}")


# In[19]:


# Generate a box plot that shows the distrubution of the tumor volume for each treatment group.
orange_out = dict(markerfacecolor='red', markersize=2)
plt.boxplot(tumor_vol_list, labels = treatment_list, flierprops=orange_out)
plt.ylabel("Final Tumor Volume (mm3)")
plt.show()


# ## Line and Scatter Plots

# In[ ]:





# In[24]:


# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin

capomulin_table = clean_study_data_complete[clean_study_data_complete["Drug Regimen"] == "Capomulin"]

mousedata = capomulin_table[capomulin_table["Mouse ID"] == "l509"]

plt.plot(mousedata["Timepoint"], mousedata["Tumor Volume (mm3)"])
plt.xlabel("Timepoint (days)")
plt.ylabel("Tumor volume (mm3)")
plt.title("Capomulin treatment of mouse l509")


# In[27]:


# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen
capomulin_table = clean_study_data_complete[clean_study_data_complete["Drug Regimen"] == "Capomulin"]
capomulin_average = capomulin_table.groupby(["Mouse ID"]).mean()
plt.scatter(capomulin_average["Weight (g)"], capomulin_average["Tumor Volume (mm3)"])
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.show()


# ## Correlation and Regression

# In[35]:


# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen
corr = st.pearsonr(capomulin_average["Weight (g)"], capomulin_average["Tumor Volume (mm3)"])
print (f"The correlation between mouse weight and the average tumor volume is {round(corr[0], 2)}")

model = st.linregress(capomulin_average["Weight (g)"], capomulin_average["Tumor Volume (mm3)"])
slope = model[0]
b = model[1]
y_values = capomulin_average["Weight (g)"] * slope + b
plt.scatter(capomulin_average["Weight (g)"], capomulin_average["Tumor Volume (mm3)"])
plt.plot(capomulin_average["Weight (g)"], y_values, color="red")
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.show()


# In[ ]:




