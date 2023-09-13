
# my base ggplot theme

# Get system name
system_name <- Sys.info()["nodename"]
mount_filesystem <- TRUE

# Define library and data paths based on system
if (system_name == "MacBook-Pro-CR-2065.local" | system_name == "dhcp172-619.laptop-zim.uni-heidelberg.de") {
  lib_path <- .libPaths()[1]
  data_path_bestageing2022 <- "/Volumes/T7CR/data/bestageing2022"
  data_path_BestAgeing <- "/Volumes/T7CR/data/BestAgeing"
  if(mount_filesystem == TRUE) {
    data_path_bestageing2022 <- "/Users/christophreich/Desktop/mount/rockerprojects/bestageing2022"  # mount -t nfs 10.55.1.185:/data/users/reich/ ~/Desktop/mount/
    data_path_BestAgeing <- "/Users/christophreich/Desktop/mount/BestAgeing"
  }
} else {  # assuming cluster
  lib_path <- "/mnt/users/reich/programs/R43/lib" 
  data_path_bestageing2022 <- "/mnt/users/reich/rockerprojects/bestageing2022"
  data_path_BestAgeing <- "/mnt/users/reich/BestAgeing"
}


# libs
require(ggplot2, lib.loc = lib_path)
require(thematic, lib.loc = lib_path)

my_base_theme <- function() {
  # inspiration
  # https://alberts-newsletter.beehiiv.com/p/ggplot-theme
  theme(
    panel.grid.minor = element_blank(),  # remove superfluous grid lines
    panel.grid.major = element_line(
      color = 'grey90', linetype = 2  # less noticable: dashed, light-grey lines
    ),
    plot.title.position = 'plot',  # align to left of chart (default is to align with inner panel)
    plot.caption.position = 'plot',
    plot.title = element_text(
      family = 'Merriweather', 
      size = rel(1.7),  # relative to base size
      margin = margin(b = 7, unit = 'mm')  # add some margin at the bottom of the plot title
    ),
    plot.subtitle = element_text(size = rel(1.1)),
    text = element_text(color = 'grey20'),  #  change the text colors from regular black to a dark grey
    axis.text = element_text(color = 'grey30'),
    panel.background = element_rect(color = 'grey90')
  )
}


# add this to plot +
#  theme_minimal(base_size = 16, base_family = 'Arial')+
#  scale_fill_manual(values = thematic::okabe_ito(6)) +
#  my_base_theme()