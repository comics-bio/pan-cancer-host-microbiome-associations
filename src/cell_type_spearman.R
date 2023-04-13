library(ggplot2)
df <- read.table("TIL-taxa.txt",header=T,row.names=1,sep="\t")

## spearman correlations
ggplot(df)+geom_smooth(aes(x=Lachnoclostridium,y=til_percentage),method='lm')+
  geom_point(data=df, aes(x=Lachnoclostridium,y=til_percentage), alpha=.9, size=1)+
  theme_bw()+
  stat_cor(data=df, aes(x=Lachnoclostridium,y=til_percentage), method = "spearman")+
  facet_wrap(~cancer,scale="free")