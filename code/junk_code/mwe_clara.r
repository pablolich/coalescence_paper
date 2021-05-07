require(ggplot2)
library(ggpubr)

p1 = ggplot(data = mtcars,
       aes(x = mpg, y = disp,
           color = cyl))+
  geom_point()+
  theme(aspect.ratio = 1,
        legend.position = 'none')

p2 = ggplot(data  = mtcars,
           aes(x = mpg, y = hp,
               color = cyl))+
  geom_point()+
  theme(aspect.ratio = 1,
        legend.key.size = unit(2,'cm'))

ggarrange(p1, p2, 
          nrow = 1)+
  ggsave(filename = "mwe_clara.pdf", 
           width = 12, height = 6)
