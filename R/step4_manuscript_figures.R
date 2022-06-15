setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to source file location
library(plyr)
library(tidyverse)
library(readxl)
library(ggrepel)

# Image 1 ####
xdat<- read_excel("../data/FDI effort by country.xlsx",
                  sheet = "FDI effort by country")
names(xdat)=str_replace(names(xdat), ' ', '_')
names(xdat)=str_replace(names(xdat), ' ', '_')
unique(xdat$`Sub-region`)
xdat=xdat[xdat$`Sub-region` %in% c("GSA17", "GSA18"),]
xdat$Gear_Type=ifelse(xdat$Country_code=='HRV'& xdat$Gear_Type=='DRB', 'TBB', xdat$Gear_Type) # Croatian dredges are treated as beam trawlers

xdat=xdat%>%dplyr::mutate(Hours_at_Sea=as.numeric(Hours_at_Sea),
                Total_Fishing_Days=as.numeric(Total_Fishing_Days),
                Gear_Type=ifelse(Gear_Type %in% c('FPO','FYK','GND','GNS','GTN','GTR','HMD','LHM','LHP','LLD', 'LLS','LTL'), 'Passive', 
                                 ifelse(Gear_Type %in% c('OTM', 'NK', 'SB', 'SV'), 'Other', Gear_Type)))%>%
  dplyr::filter(!is.na(Hours_at_Sea))%>%
  dplyr::group_by(year,Gear_Type )%>%
  dplyr::summarise(hours=sum(Hours_at_Sea), FD=sum(Total_Fishing_Days))%>%
  dplyr::filter(year >=2019)%>%
  pivot_longer(-c(year, Gear_Type))%>%
  dplyr::mutate(year=paste0('y',year))%>%
  pivot_wider(names_from = year, values_from = value)%>%
  dplyr::mutate(delta=round((((y2020/y2019)-1)*100), digits=2))%>%
  arrange(name, Gear_Type)%>%
  rbind(data.frame(Gear_Type=c('OTB','TBB','PTM','PS', 'DRB','Passive','Other'),
                   name=rep('Fishing_hours',7),
                   y2019=rep(999,7), y2020=rep(999,7),
                   delta=c(-3, -19, -6, -1,0,0,0)))


xdat=xdat[xdat$Gear_Type %in% c('OTB', 'Passive', 'PS', 'PTM','TBB'),]

ggplot()+
  geom_col(aes(x=Gear_Type, y=delta, fill=name), data=xdat, 'dodge',color='black')+
  scale_fill_brewer(palette = "YlOrRd",name = "", labels = c("Fishing Days", "Fishing Hours", "Hours at Sea"))+
  theme_bw()+
  ylab('Change in effort 2020 vs 2019 (%)')+
  xlab('Metièr')+
  scale_y_continuous(breaks = seq(-30,20,5))+
  theme(legend.position = 'bottom')

ggsave('../images/Fig1_effort_comparison.JPG', width = 15,height = 10,units='cm', dpi=500)

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


ggplot(data=CMSY_trajectories, aes(x=year, y=BBmsy, color=ref_year))+
  geom_line()+
  facet_wrap(~stock, scales='free')+
  ggtitle(expression(Trajectories~of~B/B[MSY]))+
  geom_hline(yintercept=1, linetype=2)+
  theme_bw()+
  labs(color='Reference year')+
  theme(legend.position = 'bottom')+
  ylab(expression(B/B[MSY]))

ggsave('../images/FigSI2_BBmsy_trajectories.JPG', width = 25, height = 25, dpi=500, units='cm')

ggplot(data=CMSY_trajectories, aes(x=year, y=FFmsy, color=ref_year))+
  geom_line()+
  facet_wrap(~stock, scales='free')+
  ggtitle(expression(Trajectories~of~F/F[MSY]))+
  geom_hline(yintercept=1, linetype=2)+
  theme_bw()+
  labs(color='Reference year')+
  theme(legend.position = 'bottom')+
  ylab(expression(F/F[MSY]))

ggsave('../images/FigSI3_FFmsy_trajectories.JPG', width = 25, height = 25, dpi=500, units='cm')




