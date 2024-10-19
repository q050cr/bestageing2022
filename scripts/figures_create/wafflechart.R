

library(ggplot2)
library(emojifont)
library(waffle)
library(ggwaffle)



data <- data.frame(
  Group = c("ACS", "CAD", "DCM", "ICM", "Control"),
  Patients = c(304, 408, 201, 296, 848)
)

ggplot(data, aes(x = "", y = Patients, fill = Group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") + 
  theme_void() +
  labs(fill = "Group", title = "")



data <- data.frame(
  Group = c("ACS", "CAD", "DCM", "ICM"),
  Patients = c(304, 408, 201, 296)
)

# Define colors for each group
group_colors <- c("#E69F00", "#009E73", "#0072B2", "#CC79A7")
group_colors <- ggthemes::gdocs_pal()(4)
pie_chart <- ggplot(data, aes(x = "", y = Patients, fill = Group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") + 
  scale_fill_manual(values = group_colors) +
  theme_void() +
  labs(fill = "Group", title = "")

pie_chart

pie_chart <- ggplot(data, aes(x = "", y = Patients, fill = Group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") + 
  scale_fill_manual(values = group_colors) +
  theme_void(base_size = 16, base_family = 'Arial') +
  labs(fill = "", title = "") +  # Remove legend title by setting fill = ""
  theme(legend.position = "right",  # Adjust legend position if necessary
        legend.title = element_blank(),  # Ensure no title for the legend
        legend.text = element_text(size = 12))  # Adjust text size for readability

pie_chart

ggsave("pie_chart.svg", pie_chart, width = 4, height = 4)

# wafflechart-----------------------------------------------------------
library(waffle)

parts <- setNames(data$Patients / 10, data$Group)
waffle_chart <- waffle(parts, rows = 10, size = 0.5, colors = ggthemes::gdocs_pal()(4))
waffle_chart <- waffle_chart +
  theme_void(base_size = 16, base_family = 'Arial') +# no x,y axis labels  
  theme(legend.position = "bottom",      # Moves legend to the bottom
        legend.direction = "horizontal") # Sets the legend horizontally

waffle_chart

#  theme_minimal(base_size = 16, base_family = 'Arial')+
#  scale_fill_manual(values = thematic::okabe_ito(6)) +
#  my_base_theme()




# use icons ---------------------------------------------------------------

waffle_data <- expand.grid(Group = data$Group, ID = 1:max(data$Patients))
waffle_data <- waffle_data[1:sum(data$Patients), ]
waffle_data$Group <- rep(data$Group, times = data$Patients)

waffle_data$label <- fontawesome('fa-heart')

ggplot(waffle_data, aes(x = ID, y = Group, colour = Group)) + 
  geom_text(aes(label = label), family = 'fontawesome-webfont', size = 4) +
  coord_fixed(ratio = 1) + 
  scale_colour_waffle() + 
  theme_void()


