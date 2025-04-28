# Example p-values (replace with your actual p-values)
p_values <- c(0.0007, 0.0037, 0.3783)

# Calculate adjusted p-values using the Benjamini-Hochberg procedure
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# Print the adjusted p-values
print(adjusted_p_values)