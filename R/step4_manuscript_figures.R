setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to source file location
library(plyr)
library(tidyverse)
library(readxl)
library(ggrepel)
library(gtable)
library(cowplot)
library(grid)

# Image 1 ####

# Load and arrange FDI data
FDI_dat<- read_excel("../data/FDI effort by country.xlsx",
                  sheet = "FDI effort by country")
names(FDI_dat)=str_replace(names(FDI_dat), ' ', '_')
names(FDI_dat)=str_replace(names(FDI_dat), ' ', '_')
FDI_dat=FDI_dat[FDI_dat$`Sub-region` %in% c("GSA17", "GSA18"),]
FDI_dat$Gear_Type=ifelse(FDI_dat$Country_code=='HRV'& FDI_dat$Gear_Type=='DRB', 'TBB', FDI_dat$Gear_Type) # Croatian dredges are treated as beam trawlers

xdatQuarterFDI=FDI_dat%>%dplyr::mutate(Hours_at_Sea=as.numeric(Hours_at_Sea),
                                 Total_Fishing_Days=as.numeric(Total_Fishing_Days),
                                 Gear_Type=ifelse(Gear_Type %in% c('FPO','FYK','GND','GNS','GTN','GTR','HMD','LHM','LHP','LLD', 'LLS','LTL'), 'Passive', 
                                                  ifelse(Gear_Type %in% c('OTM', 'NK', 'SB', 'SV'), 'Other', Gear_Type)))%>%
  dplyr::filter(!is.na(Hours_at_Sea))%>%
  dplyr::group_by(year,Gear_Type,Quarter )%>%
  dplyr::summarise(hours=sum(Hours_at_Sea), FD=sum(Total_Fishing_Days))%>%
  dplyr::filter(year >=2019)%>%
  pivot_longer(-c(year, Gear_Type, Quarter))%>%
  dplyr::mutate(year=paste0('y',year))%>%
  pivot_wider(names_from = year, values_from = value)%>%
  dplyr::mutate(delta=round((((y2020/y2019)-1)*100), digits=2))%>%
  arrange(name, Gear_Type)

xdatQuarterFDI=xdatQuarterFDI[xdatQuarterFDI$Gear_Type %in% c('OTB', 'Passive', 'PS', 'PTM','TBB'),]

# Load and arrange effort data from AIS. These comes from the manuscript "COVID-19 lockdowns reveal the resilience of Adriatic Sea fisheries to forced fishing effort reduction", available at https://www.nature.com/articles/s41598-022-05142-w
effortCoro=read_csv("../data/effort_Coro_2022.csv")%>%
  dplyr::filter(years>=2019)%>%
  pivot_longer(-c(years, gear, quarter))%>%
  dplyr::mutate(years=paste0('y',years),
                name='Fishing_hours')%>%
  pivot_wider(names_from = years, values_from = value)%>%
  dplyr::mutate(delta=round((((y2020/y2019)-1)*100), digits=2))%>%
  arrange(name, gear)
names(effortCoro)=c('Quarter', 'Gear_Type', names(effortCoro)[3:ncol(effortCoro)])

xdatEffort=rbind(xdatQuarterFDI, effortCoro)

# define function to arrange legend. Code from https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}



xdatEffortYear=xdatEffort%>%
  dplyr::group_by(Gear_Type, name)%>%
  dplyr::summarise(y2019=sum(y2019),
                   y2020=sum(y2020))%>%
  dplyr::mutate(delta=round((((y2020/y2019)-1)*100), digits=2))

xdatEffortYear$pos=ifelse(xdatEffortYear$name=='FD', -20,ifelse(xdatEffortYear$name=='Fishing_hours',-10,0))

?geom_label
# plot
p=ggplot()+
  geom_col(aes(x=Quarter, y=delta, fill=name), data=xdatEffort, 'dodge',color='black')+
  geom_label(
    data    = xdatEffortYear,
    mapping = aes(x = 4.8, y = pos, label = delta, fill=name), color='black', show.legend = F)+
  annotate(geom= 'text', x=4.8,y=10,label='Overall')+
  scale_fill_brewer(palette = "YlOrRd",name = "", labels = c("Fishing Days", "Fishing Hours", "Hours at Sea"))+
  scale_color_brewer(palette = "YlOrRd",name = "", labels = c("Fishing Days", "Fishing Hours", "Hours at Sea"))+
  theme_bw()+
  xlab('Quarter')+
  ylab('Change in effort 2020 vs 2019 (%)')+
  scale_y_continuous(breaks = seq(-100,100,10))+
  scale_x_continuous(breaks = seq(0,4,1))+
  theme(legend.position = 'bottom')+
  facet_wrap(~Gear_Type);p

p.new=p +
  guides(fill = guide_legend(title.position = "top",
                             label.position = "bottom",
                             nrow = 1))

png('../images/Fig1_effort_comparison.png', width = 42,height = 17,units='cm', res=500)
grid.draw(shift_legend(p.new))
dev.off()



# Image 2 ####
xdatF2= read_csv("../data/CMSY_INPUT_STOCKS_PRODUCTION.csv")%>%
  dplyr::filter(yr %in% 2019:2020)%>%
  dplyr::mutate(yr=paste0('Y',yr))%>%
  pivot_longer(-c(Stock, yr))%>%
  pivot_wider(names_from = yr, values_from = value)%>%
  dplyr::mutate(delta=Y2020/Y2019)%>%
  dplyr::select(Stock, name,delta)%>%
  pivot_wider(names_from = name, values_from = delta)

stocknames=data.frame(Stock=c("HKE_17_18","DPS_17_18_19","MUT_17_18","SOL_17","PIL_17_18","ANE_17_18" ,"CTC_17",       "MTS_17"),
                      name=c('European hake', 'deepwater pink shrimp', 'red mullet', 'common sole', 'sardine', 'anchovy','common cuttlefish', 'spot-tail mantis shrimp' ))

xdatF2=merge(xdatF2, stocknames)

ggplot(data=xdatF2, aes(bt,ct, label=substr(Stock,1,3))) + 
  annotate("rect", xmin = Inf, xmax = 1, ymin = Inf, ymax = 1, fill= "orange")  + 
  annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = 1 , fill= "yellow") + 
  annotate("rect", xmin = 1, xmax = Inf, ymin = 1, ymax = -Inf, fill= "green") + 
  annotate("rect", xmin = 1, xmax = -Inf, ymin = Inf, ymax = 1, fill= "red") + 
  geom_point()+
  geom_label_repel()+
  ylab(expression(Catch [2020]/Catch[2019]))+
  xlab(expression(Index[2020]/Index[2019]))


ggsave('../images/Fig2_Kobe_catches.JPG', width = 15, height = 15, dpi=500, units='cm')

# Image 3 ####
rawdat=read_csv("../results/CMSY_estimate_combined.csv")
store_res=list()
for(i in 1:length(unique(rawdat$stock))){
  
  istock=rawdat[rawdat$stock==unique(rawdat$stock)[i],]
  istock19=istock[istock$year==2019,]
  b19=round(quantile(istock19$BBmsy,c(0.5,0.025,0.975)), digits=3)
  f19=round(quantile(istock19$FFmsy,c(0.5,0.025,0.975)), digits=3)
  
  
  
  istock20=istock[istock$year==2020,]
  b20=round(quantile(istock20$BBmsy,c(0.5,0.025,0.975)), digits=3)
  f20=round(quantile(istock20$FFmsy,c(0.5,0.025,0.975)), digits=3)
  
  store_res[[i]]=rbind(data.frame(value=c(b19,f19), 
                                  type=c('BBmsy', 'BBmsyLow', 'BBmsyHigh', 'FFmsy', 'FFmsyLow','FFmsyHigh'), 
                                  year=2019, Stock=unique(rawdat$stock)[i]),
                       
                       data.frame(value=c(b20,f20), 
                                  type=c('BBmsy', 'BBmsyLow', 'BBmsyHigh', 'FFmsy', 'FFmsyLow','FFmsyHigh'), 
                                  year=2020, Stock=unique(rawdat$stock)[i]))
  
}

dat=plyr::ldply(store_res)%>%
  pivot_wider(names_from = type, values_from = value)

dat=merge(dat, stocknames)
library(viridis)
mycols=rainbow(32)[seq(1,32,4)]


ggplot() + 
  annotate("rect", xmin = Inf, xmax = 1, ymin = Inf, ymax = 1, fill= "orange")  + 
  annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = 1 , fill= "yellow") + 
  annotate("rect", xmin = 1, xmax = Inf, ymin = 1, ymax = -Inf, fill= "green") + 
  annotate("rect", xmin = 1, xmax = -Inf, ymin = Inf, ymax = 1, fill= "red") + 
  geom_point(data=dat, aes(BBmsy,FFmsy, pch=factor(year)),size=5)+
  ylab(expression(F/F[MSY]))+
  xlab(expression(B/B[MSY]))+
  labs(color='Stock', pch='Year')+
  theme(legend.position = 'bottom', 
        legend.box = "vertical")+
  guides(pch=guide_legend(override.aes=list(size=2)),
         color=guide_legend(override.aes=list(size=1)))+
  geom_path(data=dat, aes(BBmsy,FFmsy, fill=factor(name)))+
  #scale_color_manual(values = mycols)+
  geom_label_repel(data=dat[dat$year=='2019',], aes(BBmsy,FFmsy, label=substr(Stock,1,3)))



ggsave('../images/Fig3_kobe_CMSY_results.JPG', width = 15, height = 15, dpi=500, units='cm')


# Image 4 ####
dat=dat%>%
  dplyr::select(Stock,year,BBmsy,FFmsy)%>%
  pivot_wider(names_from = year, values_from = c(FFmsy, BBmsy))
dat$B=dat$BBmsy_2020/dat$BBmsy_2019
dat$F=dat$FFmsy_2020/dat$FFmsy_2019

ggplot(data=dat, aes(B,F, label=substr(Stock,1,3))) + 
  annotate("rect", xmin = Inf, xmax = 1, ymin = Inf, ymax = 1, fill= "orange")  + 
  annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = 1 , fill= "yellow") + 
  annotate("rect", xmin = 1, xmax = Inf, ymin = 1, ymax = -Inf, fill= "green") + 
  annotate("rect", xmin = 1, xmax = -Inf, ymin = Inf, ymax = 1, fill= "red") + 
  geom_point()+
  scale_x_continuous(breaks=seq(0.8,1.2,0.05))+
  scale_y_continuous(breaks=seq(0.7,1.2,0.05))+
  geom_label_repel()+
  ylab(expression(paste(F/F[MSY],'2020 over ',F/F[MSY], '2019' )))+
  xlab(expression(paste(B/B[MSY],'2020 over ',B/B[MSY], '2019' )))

ggsave('../images/Fig4_kobe_CMSY_difference.JPG', width = 15, height = 15, dpi=500, units='cm')


# Images 5 and 6
xfiles=list.files("../results")
xfiles=xfiles[grep('Out',xfiles)]

CMSY_trajectories=NULL

for(i in 1:length(xfiles)){
  
  idat=read_csv(paste0("../results/", xfiles[i]))
  
  
  stk=unique(idat$Stock)
  
  FFmsy=idat[,grep('F.Fmsy', names(idat))]
  FFmsy=FFmsy[,3:ncol(FFmsy)]
  FFmsy=data.frame(t(FFmsy))
  colnames(FFmsy)=stk
  FFmsy$year=as.numeric(str_remove(rownames(FFmsy), 'F.Fmsy'))
  FFmsy$year=ifelse(nchar(FFmsy$year)==1, paste0('0', FFmsy$year), FFmsy$year)
  FFmsy$year=as.numeric(ifelse(FFmsy$year<=20, paste0('20', FFmsy$year), paste0('19', FFmsy$year)))
  FFmsy=pivot_longer(FFmsy, -year, values_to = 'FFmsy', names_to = 'stock')%>%
    arrange(stock, year)%>%
    dplyr::filter(FFmsy>0)
  
  BBmsy=idat[,grep('B', names(idat))]
  BBmsy=BBmsy[,13:ncol(BBmsy)]
  BBmsy=data.frame(t(BBmsy))
  colnames(BBmsy)=stk
  BBmsy$year=as.numeric(str_remove(rownames(BBmsy), 'B'))
  BBmsy$year=ifelse(nchar(BBmsy$year)==1, paste0('0', BBmsy$year), BBmsy$year)
  BBmsy$year=as.numeric(ifelse(BBmsy$year<=20, paste0('20', BBmsy$year), paste0('19', BBmsy$year)))
  BBmsy=pivot_longer(BBmsy, -year, values_to = 'B', names_to = 'stock')%>%
    arrange(stock, year)
  
  refBbmsy=data.frame(idat[,grep('Bmsy', names(idat))][,1], stock=stk)
  BBmsy=BBmsy%>%left_join(refBbmsy, by='stock')%>%
    dplyr::mutate(BBmsy=B/Bmsy)%>%
    dplyr::select(-B, -Bmsy)
  
  i_CMSY_trajectories=left_join(FFmsy, BBmsy,by = c("year", "stock"))
  i_CMSY_trajectories$ref_year=substr(xfiles[i],1,4)
  CMSY_trajectories=rbind(CMSY_trajectories, i_CMSY_trajectories)

  
}

CMSY_trajectories$stock=substr(CMSY_trajectories$stock,1,3)


pbiomass=ggplot(data=CMSY_trajectories, aes(x=year, y=BBmsy, color=ref_year))+
  geom_line()+
  facet_wrap(~stock, scales='free')+
  ggtitle(expression(Trajectories~of~B/B[MSY]))+
  geom_hline(yintercept=1, linetype=2)+
  theme_bw()+
  labs(color='Reference year')+
  theme(legend.position = 'bottom')+
  ylab(expression(B/B[MSY]))

pbiomass.new=pbiomass +
  guides(color = guide_legend(title.position = "top",
                             label.position = "bottom",
                             nrow = 1))

png('../images/FigSI2_BBmsy_trajectories.JPG', width = 25,height = 15,units='cm', res=500)
grid.draw(shift_legend(pbiomass.new))
dev.off()

pF=ggplot(data=CMSY_trajectories, aes(x=year, y=FFmsy, color=ref_year))+
  geom_line()+
  facet_wrap(~stock, scales='free')+
  ggtitle(expression(Trajectories~of~F/F[MSY]))+
  geom_hline(yintercept=1, linetype=2)+
  theme_bw()+
  labs(color='Reference year')+
  theme(legend.position = 'bottom')+
  ylab(expression(F/F[MSY]))

pF.new=pF +
  guides(color = guide_legend(title.position = "top",
                              label.position = "bottom",
                              nrow = 1))

png('../images/FigSI3_FFmsy_trajectories.JPG', width = 25,height = 15,units='cm', res=500)
grid.draw(shift_legend(pF.new))
dev.off()

ggsave('../images/FigSI3_FFmsy_trajectories.JPG', width = 25, height = 25, dpi=500, units='cm')




