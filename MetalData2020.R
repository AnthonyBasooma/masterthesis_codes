#===============================================================================
#Author:     Anthony Basooma
#Current:    Research Assitant: Ecological Modelling and Ecotoxicology
#            Antwerp University, Belgium
#===         Supervised by Prof. Lieven Bervoets
##           Title: Trace metal concentration in the abiotic and biotic components
#of River Rwizi Ecosystem, Uganda
#===============================================================================
setwd("E:/Masters/Year 2/Research/manuscript/data/rscripts")
#===============================================================================
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(extrafont)
library(outliers)
library(pgirmess)
library(ITNr)
library(ggthemes)
library(psych)
library(ggpubr)
library(ibd)
options(scipen = 999)
#================================================================================
#Loading the master dataset for analysis####

MetalData <- read.csv(file ="MetalData.csv", header = T, strip.white = TRUE)%>%
  gather(key = "metals", value = "conc", Cd:As)%>%mutate(conc=replace_na(conc, 0))
HEI      <- read.csv(file = "HEI.csv", header = T, strip.white = T)
TOC      <- read.csv(file = "toc.csv", header = T, strip.white = T)


#Determine the limit of detectionfor the different metals in the tissues
#================================================================================
blanks<- MetalData%>% filter(SampleN%in% c("blank"))
dwmax<- max(MetalData$dw)

MDL<- ddply(blanks,. (metals), summarise, blankmean=mean(conc), blanksd=sd(conc),
            LOD= blankmean+blanksd*3, MDL= (LOD*5*10)/(1000*dwmax))#1000 L
#================================================================================
#Analysis of trace metals in the surface water among the Sites
H2O_metal<- MetalData%>% filter(SampleD%in% "water")%>%
  mutate(metals=as.factor(metals))%>%
  select(Sites,sitecd,metals, conc)%>%
  mutate(sites= factor(sitecd, unique(sitecd)))#%>%
  group_by(sitecd, metals)%>%
  summarise(minl= round_df(min(conc), digits=3),
         maxl= round_df(max(conc), digits = 2), .groups="drop")%>%
  mutate(minmax= paste(minl,"-", maxl))%>%
  select(sitecd, metals, minmax)%>%
  spread("sitecd", "minmax", fill=0)

signPb<- H2O_metal%>%filter(metals%in%"Pb")

kruskal.test(conc~sites, signPb)
kruskalmc(conc~sites, signPb)
pairwise.wilcox.test(signPb$conc, signPb$sites)

meanc<- H2O_metal%>%group_by(sitecd, metals)%>%
  summarise(means=mean(conc),
            menrd= round_df(means, digits=2),
            sdc= sd(conc),
            sdrd= round_df(sdc, 2),
            minl= round_df(min(conc)),
            maxl= round_df(max(conc)), .groups="drop")%>%
  mutate(minmax= paste(minl,"-", maxl))


#--------------------------------------------------------------------------------
bg<- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           axis.line = element_line(colour = "black"),
           plot.title = element_text(hjust = 0.5))
bg2<- theme(axis.title.x = element_blank())
#--------------------------------------------------------------------------------

ggplot(H2O_metal, aes(site, conc))+
  stat_boxplot(geom = "errorbar", linetype=1, width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~metals, scales = "free")+
  theme_bw()+bg+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"))+
  labs(y="Metal levels (?g/L)", x="Sites")+
  stat_compare_means(label = "p.format", label.x.npc = c(0.4, 0.5),
                     method = "anova")

#====================
ggplot(H2O_metal, aes(sitecd, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~metals, scales = "free_y")+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5))+
  labs(y="Zn Âµg/g")+ggtitle("Zn")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.1)))
#==========


#------------------------------------------------------------------------------
#Trace metal levels in the sediment samples

Sed_metal<- MetalData%>% filter(SampleD%in% "sediment")%>%
  mutate(metals=as.factor(metals), conc=round_df(conc, digits = 4))%>%
 filter(!conc%in%c(1.0070, 3.6259, 13.7447, 38.7336, 387.5642))%>%
  dplyr::select(Sites,sitecd,metals, conc)%>%
  mutate(sites= factor(sitecd, unique(sitecd)))

supl_min<- Sed_metal%>%group_by(metals, sitecd)%>%
  dplyr::summarise(medval= min(conc), .groups="drop")%>%
  mutate(min_conc= case_when(metals=="Al"~medval/1000,
                             metals=="Mn"~medval/1000,
                             metals=="Fe"~medval/1000,
                             metals=="As"~medval,
                             metals=="Au"~medval,
                             metals=="Cd"~medval,
                             metals=="Cu"~medval,
                             metals=="Co"~medval,
                             metals=="Pb"~medval,
                             metals=="Zn"~medval,
                             metals=="Hg"~medval))%>%
  select(metals, sitecd, min_conc)%>%
  spread("sitecd", "min_conc", fill = 0)

supl_max<- Sed_metal%>%group_by(metals, sitecd)%>%
  dplyr::summarise(medval= max(conc), .groups="drop")%>%
  mutate(max_conc= case_when(metals=="Al"~medval/1000,
                             metals=="Mn"~medval/1000,
                             metals=="Fe"~medval/1000,
                             metals=="As"~medval,
                             metals=="Au"~medval,
                             metals=="Cd"~medval,
                             metals=="Cu"~medval,
                             metals=="Co"~medval,
                             metals=="Pb"~medval,
                             metals=="Zn"~medval,
                             metals=="Hg"~medval))%>%
  select(metals, sitecd, max_conc)%>%
  spread("sitecd", "max_conc", fill = 0)
#write.csv(supl_min, "supl_max.csv")
#write.csv(supl_max, "supl_min.csv")

metalsd<- Sed_metal%>%filter(metals%in%"Al")
ggplot(Sed_metal, aes(sites, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~metals, scales = "free_y")+
  theme_bw()+bg+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Trace metal concentration (Âµg/g)", x="Sites")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                      ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                          symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

#=============
#Allumi
al<-ggplot(metalsd, aes(sites, conc/1000))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
    labs(y="Al mg/g dw", title = "Al")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))



#===========================
#Gold in the sediment
#==========================
metalsau<- Sed_metal%>%filter(metals%in%"Au")
Au<-ggplot(metalsau, aes(sites, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5))+
  labs(y="Au µg/g dw")+ggtitle("Au")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

##================
metalsas<- Sed_metal%>%filter(metals%in%"As")
As<-ggplot(metalsas, aes(sites, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="As µg/g dw")+ggtitle("As")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

#======
metalscd<- Sed_metal%>%filter(metals%in%"Cd")
Cd<-ggplot(metalscd, aes(sites, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Cd µg/g dw" )+ggtitle("Cd")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))


metalscu<- Sed_metal%>%filter(metals%in%"Cu")
Cu<-ggplot(metalscu, aes(sites, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Cu µg/g dw")+ggtitle("Cu")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

metalsco<- Sed_metal%>%filter(metals%in%"Co")
Co<-ggplot(metalsco, aes(sites, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5))+
  labs(y="Co µg/g dw")+ggtitle("Co")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))


png("sedimentplt1.png", units="in", width=5, height=5, res=300)
gridExtra::grid.arrange(al, Au, As, Cd, Cu, Co, ncol=2)
dev.off()


#===========================
# in the sediment
#==========================
metalsfe<- Sed_metal%>%filter(metals%in%"Fe")
Fe<-ggplot(metalsfe, aes(sites, conc/1000))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Fe mg/g dw", title = "Fe")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

##================
metalshg<- Sed_metal%>%filter(metals%in%"Hg")
Hg<-ggplot(metalshg, aes(sites, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Hg µg/g dw")+ggtitle("Hg")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "ns")),
                     family="Calibri", hide.ns = FALSE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
Hg
#======
metalsmn<- Sed_metal%>%filter(metals%in%"Mn")
Mn<-ggplot(metalsmn, aes(sites, conc/1000))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Mn mg/g dw")+ggtitle("Mn")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

metalspb<- Sed_metal%>%filter(metals%in%"Pb")
Pb<-ggplot(metalspb, aes(sites, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Pb µg/g dw")+ggtitle("Pb")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

metalszn<- Sed_metal%>%filter(metals%in%"Zn")
Zn<-ggplot(metalszn, aes(sites, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 13),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5))+
  labs(y="Zn µg/g dw")+ggtitle("Zn")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
Zn

png("Figure2.png", units="in", width=7, height=8, res=350)
gridExtra::grid.arrange(al, Au, As, Cd, Cu, Co,Fe, Hg, Mn, Pb, Zn, ncol=3)
dev.off()


metalAu<- Sed_metal%>%filter(metals%in%"Fe")
kruskal.test(conc~sites, metalAu)

kruskalmc(conc~sites, metalAl)

#================================================================================
#Total organic carbon and clay content
tocdata<- TOC%>%gather("par", "val", TOC:Clay)#%>%
  filter(!is.na(val))%>%group_by(sitecd, par)%>%
summarise(yy=max(val))%>%
  spread("par", "yy", fill = 0)

write.csv(tocdata, "maxtoc.csv")

clatoc<- TOC%>%select(sitecd, TOC, Clay)%>%gather("par", "val", TOC:Clay)%>%
  filter(!is.na(val))
clplot<- aggregate(val~par+sitecd, clatoc, median)
cla<- clplot%>%filter(par%in%"Clay")

kruskal.test(val~sitecd, cla)


  ggplot(clatoc, aes(x=sitecd, y=val/10))+
    geom_col(fill="grey35")+
    theme_bw()+bg+
    scale_fill_grey(start = 0.2, end = 0.4)+
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)))+
    facet_wrap(~par, scales = "free")+
    labs(y="Percentage Composition", x="Sites")+
    theme(legend.position = "none",
          text=element_text(family = "Calibri", size = 12))+
    stat_compare_means(label = "p.signif", method = "wilcox",
                       ref.group = ".all.",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                        0.05, 1),
                                          symbols = c("bc", "abc", "ab", "b", "a")),
                       family="Calibri", hide.ns = TRUE,
                       label.y.npc="top")

  ggplot(clatoc, aes(x=sitecd, y=val/10))+
    stat_identity(fill="grey35", geom = "bar")+
    theme_bw()+bg+
    scale_fill_grey(start = 0.2, end = 0.4)+
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)))+
    facet_wrap(~par, scales = "free")+
    labs(y="Percentage Composition", x="Sites")+
    theme(legend.position = "none",
          text=element_text(family = "Calibri", size = 12))+
    stat_compare_means(label = "p.signif", method = "wilcox",
                       ref.group = ".all.",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                        0.05, 1),
                                          symbols = c("bc", "abc", "ab", "b", "a")),
                       family="Calibri", hide.ns = TRUE,
                       label.y.npc="top")


  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
    facet_wrap(~metals, scales = "free")+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Al Âµg/g", title = "Al")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.1)))



  tocs<- tocplot%>%filter(sedcont%in%"Clay")
  kruskal.test(val~sitecd, tocs)
  kruskalmc(val~sitecd, tocs)

ggplot(tocplot, aes(sitecd, val))+
  stat_boxplot()+
  facet_wrap(~sedcont, scales = "free")+
  theme_bw()+bg+
  theme(text=element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
          legend.position = "none")+
  labs(y="% composition", x="Sites")


#================================================================================
kruskal.test(TOC~site, data = TOC)
kruskalmc(TOC~site, data = TOC, p.adjust.methods="bonf")
#--------------------------------------------------------------------------------
kruskal.test(Clay~site, data = TOC)
kruskalmc(Clay~site, data = TOC, p.adjust.methods="bonf")
#=====================
#Multiple regressions


#Heavy Metal evlauation index
#---------------------------------
MetalHEI<- left_join(H2O_metal, HEI, by=c("metals"= "metals"))

HEIplot<- MetalHEI%>% dplyr::select(Sites, sitecd, metals, conc, Hmax)%>%
  filter(!metals%in%c("Au", "Co"))%>% group_by(Sites, sitecd, metals, Hmax)%>%
  summarise(conc2=mean(conc), .groups="drop")%>% mutate(Hval= conc2/Hmax)%>%
  filter(!Hval>=10)

ggplot(HEIplot, aes(sitecd, Hval, fill=sitecd))+
  geom_bar(stat = "identity")+
  scale_fill_colorblind() +
  theme_bw()+bg+
  theme(text = element_text(family = "Calibri"),
        legend.position = "none",
        strip.text.x = element_text(face = "bold"))+
  labs(y="HEI Values", x="Sites")+
  geom_hline(yintercept = 3, size=1, color="black")+
  scale_y_continuous(expand =  expand_scale(mult = c(0, 0.1)))

#============================================================================
#Physiochemical parameter analysis for each site
Params<- MetalData%>%dplyr::select(Sites, GD, SD, Temperature, Conductivity, pH)%>%
  filter(GD%in%c("Water"))
kruskal.test(pH~SD, data = Params)
kruskalmc(conc~fishsp, data = FishData2)


#=============================================================================
pcparam<- prcomp(Params[,-1], center = TRUE, scale. = TRUE)
attributes(pcparam)
pcparam$center
print(pcparam)
summary(pcparam)
pairs.panels(pcparam$x, gap=0)#all 0s
#===================
ggbiplot(pcparam, obs.scale = 1, var.scale = 1, groups = Params$Sites,
         ellipse = TRUE, circle = T, ellipse.prob = 0.90)+
  theme_bw()+bg+
  scale_color_manual(values = c("black", "blue", "green",
                                "grey40", "red", "grey80"))+
  theme(text = element_text(family = "Cambria", size = 12))
#=============================================================================

#Multiple regression to determine the concentration in the sed $ water vs invert

MultData<- MetalData%>%filter(SampleD%in%c("sediment", "water", "inverts"))%>%
  select(SampleD, sitecd, metals, conc)%>% group_by(SampleD, sitecd,metals)%>%
  summarise(meanc= mean(conc))%>%spread("SampleD", "meanc", fill = 0)

Cu<- MultData%>%filter(metals%in%c("Zn"))
lmg<- lm(inverts~sediment+water, data = Cu)
summary(lmg)

cor.test(Cu$sediment, Cu$inverts, method = "spearman")

#================================================================================
#Trace metal level in the invertebate
Inverts<- MetalData%>%filter(SampleD%in%c("inverts"))%>%
  dplyr::select(sitecd, metals, conc )%>%
  mutate(metals=as.factor(metals), conc=round_df(conc, digits = 4))%>%
  filter(!conc%in%c(0.7736, 226.3125, 2786.9213, 18.2910, 727.8214))

supl<- Inverts%>%group_by(metals, sitecd)%>%
  summarise(medval= median(conc), .groups="drop")%>%
  spread("sitecd", "medval", fill = 0)

supli_min<- Inverts%>%group_by(metals, sitecd)%>%
  dplyr::summarise(medval= min(conc), .groups="drop")%>%
  mutate(min_conc= case_when(metals=="Al"~medval/1000,
                             metals=="Mn"~medval/1000,
                             metals=="Fe"~medval/1000,
                             metals=="As"~medval,
                             metals=="Au"~medval,
                             metals=="Cd"~medval,
                             metals=="Cu"~medval,
                             metals=="Co"~medval,
                             metals=="Pb"~medval,
                             metals=="Zn"~medval/1000,
                             metals=="Hg"~medval))%>%
  select(metals, sitecd, min_conc)%>%
  spread("sitecd", "min_conc", fill = 0)

supli_max<- Inverts%>%group_by(metals, sitecd)%>%
  dplyr::summarise(medval= max(conc), .groups="drop")%>%
  mutate(max_conc= case_when(metals=="Al"~medval/1000,
                             metals=="Mn"~medval/1000,
                             metals=="Fe"~medval/1000,
                             metals=="As"~medval,
                             metals=="Au"~medval,
                             metals=="Cd"~medval,
                             metals=="Cu"~medval,
                             metals=="Co"~medval,
                             metals=="Pb"~medval,
                             metals=="Zn"~medval/1000,
                             metals=="Hg"~medval))%>%
  select(metals, sitecd, max_conc)%>%
  spread("sitecd", "max_conc", fill = 0)



#===========================
#Individual graphs
#===========================
#Allumi
metalsal<- Inverts%>%filter(metals%in%"Al")
Al<-ggplot(metalsal, aes(sitecd, conc/1000))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Al mg/g dw", title = "Al")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

#===========================
#Gold in the sediment
#==========================
metalsau<- Inverts%>%filter(metals%in%"Au")
Au<-ggplot(metalsau, aes(sitecd, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5))+
  labs(y="Au µg/g dw")+ggtitle("Au")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

##================
metalsas<- Inverts%>%filter(metals%in%"As")
As<-ggplot(metalsas, aes(sitecd, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="As µg/g dw")+ggtitle("As")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

#======
metalscd<- Inverts%>%filter(metals%in%"Cd")
Cd<-ggplot(metalscd, aes(sitecd, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Cd µg/g dw")+ggtitle("Cd")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

metalscu<- Inverts%>%filter(metals%in%"Cu")
Cu<-ggplot(metalscu, aes(sitecd, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Cu µg/g dw")+ggtitle("Cu")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

metalsco<- Inverts%>%filter(metals%in%"Co")
Co<-ggplot(metalsco, aes(sitecd, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5))+
  labs(y="Co µg/g dw")+ggtitle("Co")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))


#===========================
# in the invertebrates
#==========================
metalsfe<- Inverts%>%filter(metals%in%"Fe")
Fe<-ggplot(metalsfe, aes(sitecd, conc/1000))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Fe mg/g dw", title = "Fe")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

##================
metalshg<- Inverts%>%filter(metals%in%"Hg")
Hg<-ggplot(metalshg, aes(sitecd, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Hg µg/g dw")+ggtitle("Hg")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "ns")),
                     family="Calibri", hide.ns = FALSE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

#======
metalsmn<- Inverts%>%filter(metals%in%"Mn")
Mn<-ggplot(metalsmn, aes(sitecd, conc/1000))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Mn mg/g dw")+ggtitle("Mn")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

metalspb<- Inverts%>%filter(metals%in%"Pb")
Pb<-ggplot(metalspb, aes(sitecd, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Pb µg/g dw")+ggtitle("Pb")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

metalszn<- Inverts%>%filter(metals%in%"Zn")
Zn<-ggplot(metalszn, aes(sitecd, conc/1000))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5))+
  labs(y="Zn mg/g dw")+ggtitle("Zn")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))



png("Figure3.png", units="in", width=7, height=8, res=350)
gridExtra::grid.arrange(Al, Au, As, Cd, Cu, Co,Fe, Hg, Mn, Pb, Zn, ncol=3)
dev.off()


invAu<- Inverts%>%filter(metals%in%"Pb")
kruskal.test(conc~sites, invAu)
kruskalmc(conc~sites, invAu)


ggplot(Inverts, aes(sites, conc))+
  stat_boxplot(geom = "errorbar", width=0.7)+
  facet_wrap(~metals, scales = "free")+
  theme_bw()+bg+
  geom_boxplot()+
  theme(text = element_text(family = "Calibri", size = 12),
        legend.position = "none",
        strip.text.x = element_text(face = "bold"))+
  labs(x="Sites", y="Trace metal levels (ug/g,dw)")+
  stat_compare_means(label = "p.format", label.x.npc = c(0.5, 0.5),
                     method = "anova")
#===================================================================================
#Human risk assessment)
DIR<- 31.5 #g Daily ingestion rate
MDQData<- MetalData%>%filter(SampleD%in%c("Muscle"), WN%in%c("wetweight"))%>%
  mutate(metals=as.factor(metals))
MDQData2<- left_join(MDQData, HEI, by=c("metals"="metals"))%>%
  dplyr::select(fishsp,metals, conc, MRL,RfD )%>%group_by(fishsp,metals, MRL)%>%
  mutate(p50= quantile(conc, 0.5), p95= quantile(conc, 0.95),
         p50Q= (MRL*70)*1000/(p50), p95Q= (MRL*70)*1000/(p95),
         HQ50= DIR/p50Q, HQ95=DIR/p95Q)

#==============================================================================
FishData<- MetalData%>%filter(GD%in%c("Fish"), WN%in%c("dryweight"))%>%
  select(sitecd, fishsp, metals, SampleD, conc )

ggplot(FishData, aes(sitecd, conc, fill=SampleD))+
  stat_boxplot(geom = "errorbar", width=0.7)+
  facet_wrap(~metals, scales = "free")+
  theme_bw()+bg+
  geom_boxplot(outlier.shape = NA)+
  theme(text = element_text(family = "Calibri", size = 12),
        legend.position = c(0.9, 0.2),
        strip.text.x = element_text(face = "bold"))+
  labs(x="Sites", y="Metal levels (ug/g,dw)", fill="Sampled Sites")+
  scale_fill_manual(values=c("grey86", "orange", "black"))
#==================
kruskal.test(conc~SampleD, data = FishData)
kruskalmc(conc~SampleD, data = FishData)

#============
FishData1<- MetalData%>%filter(GD%in%c("Fish"), WN%in%c("dryweight"))%>%
  select(sitecd, fishsp, metals, SampleD, conc )

liver<- FishData1%>%filter(SampleD%in%"Liver")%>%
  mutate( conc=round_df(conc, digits = 4))%>%
  filter(!conc%in%c(55.4514, 0.7322, 7.9353, 1.1112, 0.2324))%>%
  group_by(sitecd, metals)%>%summarise(med= max(conc))%>%
  spread("sitecd", "med", fill = 0)
#write.csv(liver, "livermax.csv")

ggplot(liver, aes(sitecd, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~metals, scales = "free_y", ncol = 3)+
  theme_bw()+bg+bg2+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"))+
  labs(y="Trace metal concentration Âµg/g")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "ns")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.1)))

gills<- FishData1%>%filter(SampleD%in%"Gills", !metals%in%"Au")%>%
  mutate( conc=round_df(conc, digits = 4))#%>%
  filter(!conc%in%c(7.0351, 15.2954))%>%ungroup()%>%
  group_by(sitecd, metals)%>%summarise(med= max(conc))#%>%
  spread("sitecd", "med", fill = 0)
#write.csv(gills, "gillmax.csv")

png('Figure4.png', width = 8, height = 7, units = 'in', res = 350)
ggplot(gills, aes(sitecd, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~metals, scales = "free", ncol = 3)+
  theme_bw()+bg+
  theme(text = element_text(family = "Calibri", size = 13),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"),
        axis.title.x = element_blank())+
  labs(y="Trace metal concentration µg/g dw")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
dev.off()
#=================
Muscle<- FishData1%>%filter(SampleD%in%"Muscle")%>%
  mutate( conc=round_df(conc, digits = 4))%>%
  filter(!conc%in%c(148.9596, 0.0516, 1.1983, 197.3964, 0.3984))%>%
  group_by(sitecd, metals)%>%summarise(med= max(conc))%>%
  spread("sitecd", "med", fill = 0)
#write.csv(Muscle, "musclemax.csv")

png('Figure5.png', width = 8, height = 7, units = 'in', res = 350)
ggplot(Muscle, aes(sitecd, conc))+
  stat_boxplot(geom = "errorbar",  width=0.7)+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~metals, scales = "free", ncol = 3)+
  theme_bw()+bg+
  theme(text = element_text(family = "Calibri", size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.background.x = element_rect(fill = "white"),
        axis.title.x = element_blank())+
  labs(y="Trace metal concentration µg/g dw")+
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = ".all.",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01,
                                                      0.05, 1),
                                        symbols = c("j", "v", "abc", "ab", "a")),
                     family="Calibri", hide.ns = TRUE,
                     label.y.npc=1)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
dev.off()
#=================
kruskal.test(conc~SampleD, data = FishData1)
kruskalmc(conc~SampleD, data = FishData)



