setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to source file location
library(plyr)
library(tidyverse)


xfiles=list.files(file.path('..', 'results'))
xfiles=xfiles[grep('model_estimates', xfiles)]
files=list()

for(x in 1:length(xfiles)){
  
  files[[x]]=read.csv(file.path('..','results', xfiles[x]))
  
  
}

files=ldply(files)
names(files)[2]='stock'

write.csv(files,'../results/CMSY_estimate_combined.csv', row.names = F)


ggplot(data=files)+
  geom_point(aes(x=BBmsy, y=FFmsy,color=factor(year)))+
  facet_wrap(~stock)+
  labs(color='End Year')+
  theme(legend.position = 'bottom')+
  ggtitle('CMSY BSM estimates for models ending in 2019 and in 2020')

ggsave('../images//FigSI1_BSM_raw_estimates.jpeg', width=20, height=20,units='cm', dpi=500)
