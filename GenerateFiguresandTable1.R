#Create Tables and map for publication

library(data.table)
library(sf)

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(0.5, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)
  
}


alldata<-readRDS("Output/finalcuratedalldata.rds") #Data available at Zenodo per the publication

#generate county level DTXSID data 

k<-data.table(unique(alldata[,c("DTXSID","StateName","CountyName")]))
k<-k[which(!is.na(k$DTXSID)),]

countsbycounty<-k[, .(count = .N), by = c("StateName","CountyName")]

#devtools::install_github("UrbanInstitute/urbnmapr")


library(tidyverse)
library(urbnmapr)

ggplot() + 
  geom_polygon(data = urbnmapr::states, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)

county_data <- left_join(countydata, counties, by = "county_fips") 


county_data$CountyName<-gsub(" County","",county_data$county_name)
county_data$StateName<-county_data$state_name

chemical_county<-left_join(county_data,countsbycounty)

chemical_county<-chemical_county[which(!chemical_county$StateName %in% c("Alaska","Hawaii")),]

states_sf <- get_urbn_map(map = "states", sf = TRUE)

chemical_county %>%
  ggplot(aes(long, lat, group = group, fill = count)) + scale_fill_distiller(palette = "Spectral")+
  geom_polygon(color = NA) +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
  labs(fill = "Number of Unique Chemicals") + 
  geom_polygon(data = urbnmapr::states, mapping = aes(x = long, y = lat, group = group),fill = NA, color = "white")+
  theme_minimal()+
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

# Make summary Table of functions assigned to specific chemicals
#The below is necessary to count instances where two functions were reported
uniquefuncsoverall<-unique(alldata[,c("DTXSID","ChemicalFunction")])
uniquefuncsoverall<-uniquefuncsoverall[!is.na(uniquefuncsoverall$DTXSID),]

singlechemfunctionlist<-unique(alldata$ChemicalFunction)
singlechemfunctionlist<-singlechemfunctionlist[!grepl(",", singlechemfunctionlist)]
singlechemfunctionlist<-singlechemfunctionlist[singlechemfunctionlist!=""]

m <- data.frame(matrix(ncol = 2, nrow = length(singlechemfunctionlist)))
colnames(m)<-c("ChemicalFunction","Count")

i<-0
for (chem in singlechemfunctionlist){
i<-i+1
a<-length(which(grepl(chem,uniquefuncsoverall$ChemicalFunction, fixed = T)))  
m$ChemicalFunction[i]<-chem
m$Count[i]<-a
}

m<-m[order(-m$Count),]
Table1<-m[1:20,]
write.csv(Table1,"Output/Table1.csv")


#Number of chemicals per function over time
library(stringi)
alldata$year<-stri_reverse(alldata$JobStartDate)
alldata$year<-str_sub(alldata$JobStartDate,-16,-12)
tmp<-stri_reverse(alldata$JobStartDate)
k<-str_locate(tmp,"/")[,1]
alldata$year<-str_sub(alldata$JobStartDate,-k+1,-k+4)

#table(alldata$year)

#Examine time course of most common chemicals for certain functions, change values of "func" to plot different functions.
#func<-"Chemical reaction regulator"
#func<-"Surfactant (surface active agent)"
#func<-"Biocide"
#func<-"Anti-scaling agent"
#func<-"Lubricating agent"
#func<-"Biocide"
func<-"Stabilizing agent"

#examples of changing use : Stabilizers, Biocides, Anti-scaling agent

funcoccurrence<-alldata[which(alldata$ChemicalFunction==func),]
funcoccurrence<-data.table(funcoccurrence)
funcoccurrence<-funcoccurrence[which(!is.na(funcoccurrence$DTXSID)),]
funcoccurrence<-funcoccurrence[which(!is.na(funcoccurrence$year)),]
topfuncs<-funcoccurrence[, .(count = .N), by = c("DTXSID")]
topfuncs<-topfuncs[order(-topfuncs$count),]
topfuncs<-topfuncs$DTXSID[1:10]

funcoccurrence<-funcoccurrence[which(funcoccurrence$DTXSID %in% topfuncs),]
countit<-funcoccurrence[, .(count = .N), by = c("year","DTXSID")]


countit[which(countit$year>2009)] %>%
  ggplot(aes(year, count, group= DTXSID, color = DTXSID)) + scale_fill_distiller(palette = "Spectral")+ geom_line(size=1)+theme_Publication()+
  ylab("Number of Ocurrences in Well Operations")+xlab("Year")+ggtitle(func)

#Example for paper text 

#Biocides: Decreasing: DTXSID0021331, DTXSID9034997 Increasing:DTXSID5023958, DTXSID8021276
#Anti-Scaling: Decreasing DTXSID8020597, increasing DTXSID2021731
#Stabilizing agent: Decreasing DTXSID1051199, increasing DTXSID4020325, then increasing DTXSID5020235

#Overall for comparison
k<-unique(alldata[,c("DisclosureId","year")])
j<-data.frame(table(k$year))

j[which(!j$Var1 %in%  c("2001","2002","2004","2005","2007","2008","2009","2323","2010")),] %>%
  ggplot(aes(Var1, Freq, group = '')) + geom_line(size=1)+theme_Publication()+
  ylab("Count")+xlab("Year")+ggtitle("Number of Disclosures by Year in the FracFocus Data")





