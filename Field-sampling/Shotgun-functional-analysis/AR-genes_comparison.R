library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

# Load the data
df <- read.csv("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Publication/tables/Fields shotgun/ARG_abundances.csv",row.names=1)

# Transform data from wide to long format
long_data <- gather(df, Key, Value, -Category)

# Calculate average values and standard error for each enzyme family by pond type
average_values <- long_data %>%
  group_by(Category) %>%
  summarize(Average = mean(Value),
            SE = sd(Value) / sqrt(n()))

# Update Category labels
average_values$Category <- recode(average_values$Category, "High" = "High MP ponds", "Low" = "Low MP ponds")

# Define colors
colors <- c("Low MP ponds" = "#6885d0", "High MP ponds" = "#cb5658")

# Create the bar chart with error bars
bar_chart <- ggplot(average_values, aes(x=Category, y=Average, fill=Category)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Average-SE, ymax=Average+SE), width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values = colors) +  # Specify your colors here
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust text angle and justification for readability
        axis.text.y = element_text(),
        legend.title = element_text(face = "bold")) +
  labs(x=NULL, y="Average reads mapped per ARG for each microbiome typ", fill="") +
  theme_pubr()  # or theme_minimal() if theme_pubr() is not defined

bar_chart

colors <- c("Low" = "#6885d0", "High" = "#cb5658")
bar_chart <- ggplot(df, aes(x=rownames(df), y=Average.ARG.abundance, fill=Category)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(y="Average reads mapped per ARG for each microbiome type",x="") +
  scale_color_manual(values = colors)  + # Specify your colors here
  theme(axis.text.x = element_text(angle = 90, hjust = 1), # Adjust text angle and justification for readability
        axis.text.y = element_text(),
        legend.title = element_text(face = "bold")) +
  theme_pubr()
bar_chart

#AR-genes
wilcox.test(Average.ARG.abundance ~ Category,data=df) #W = 16, p-value = 0.02857
